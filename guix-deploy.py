#!/usr/bin/env python2.7

#=======================================================================
# Author(s): Ben Woodcroft
#
# Copyright 2016
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import argparse
import subprocess
import logging
import re
import os, os.path

from subprocess import CalledProcessError

class GuixInfo:
    def __init__(self, guix_package_name, guix_exe):
        cmd = "%s package --show=%s" % (guix_exe, guix_package_name)
        logging.info("Running cmd: %s" % cmd)
        output = subprocess.check_output(cmd, shell=True)

        if output=='':
            raise Exception("Unfortunately it appears there is no '%s' package in guix" % guix_package_name)

        infos = {}
        r = re.compile("(.+?): (.*)")
        last_key = None
        for l in output.split("\n"):
            if len(l) == 0: continue
            if l[0] == '+':
                if not last_key:
                    raise Exception("info started with a + ?!!?!: %s" % l)
                infos[last_key] += ' '
                infos[last_key]+= l[1:]
            else:
                o = r.search(l)
                if o:
                    last_key = o.groups(0)[0].strip()
                    value = o.groups(0)[1].strip()
                    if last_key in infos:
                        raise Exception("There appears to be multiple packages available, use "
                                        "\"%s package --show=%s\" and guix_module -i %s --guix_package_name @.. with a specific version" % (guix_exe, guix_package_name, guix_package_name))
                    infos[last_key] = value
                else:
                    raise Exception("Badly parsed info line: %s" % l)

        # these throw errors if they don't exist, so don't need to check
        self.name = infos['name']
        self.version = infos['version']
        self.synopsis = infos['synopsis']
        self.description = infos['description']

class ProfileInfo:
    def __init__(self, base, name, guix_info):
        self.base = base
        self.profile_name = name
        self.guix_info = guix_info

    def profile_directory(self):
        return os.path.join(self.base, self.profile_name)

    def profile_path(self):
        return os.path.join(self.profile_directory(), self.guix_info.version)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--install', metavar='package_name', help='Name of the profile to be installed')
    parser.add_argument('--guix_package_name', metavar='name', help='the name of the package in guix (default: just the one specified with --install)')
    parser.add_argument('--no_rsync', action='store_true',
                        help='skip rsync step (build step only)', default=False)
    parser.add_argument('--rsync_only', action='store_true',
                        help='skip build step (rsync step only)', default=False)
    args = parser.parse_args()
    if not args.install and not args.rsync_only:
        raise Exception("Either --install or --rsync_only must be specified.")

    logging.basicConfig(level=logging.INFO)

    # Constant paths
    euramoo = "uqbwoodc@ssh1.qriscloud.org.au"
    awoonga = "uqbwoodc@awoonga.qriscloud.org.au"
    base_path = '/RDS/Q0227/profiles/base/'
    guix_git = '/home/ben/git/guix-euramoo'
    # There is also direct references to /RDS/Q0227 here and in the guix git configure.

    os.chdir(guix_git)
    logging.info("Running make on the guix directory ..")
    #subprocess.check_call("guix environment guix -- make", shell=True)


    if not args.rsync_only:
        guix_exe = "%s/pre-inst-env guix" % guix_git
        logging.info("Using guix binary: '%s'" % guix_exe)

        guix_package_name = args.install
        if args.guix_package_name:
            guix_package_name = args.guix_package_name

        # First check if all packages are known
        # make sure that the guix package exists, and gather info on it at the same time
        # using guix package --show=guix_package_name
        logging.info("Testing if %s is a known Guix package" % guix_package_name)
        guix_info = GuixInfo(guix_package_name, guix_exe)
        profile = ProfileInfo(base_path, args.install, guix_info)

        # create or make sure that the profile directory exists
        profile_dir = profile.profile_directory()
        profile_path = profile.profile_path()
        logging.info("Profile to be created at %s" % profile_path)

        # ensure profile path doesn't exist already, we don't want to update packages
        if os.path.exists(profile_path):
            raise Exception("A profile already exists at %s" % profile_path)
        if not os.path.exists(profile_dir):
            os.makedirs(profile_dir)

        # Build the packages. Do this separately so if the build fails none of
        # the packages get installed into the profile directory.
        logging.info("Building packages before installing into the profile")
        cmd = "%s build %s" % (guix_exe, guix_package_name)
        logging.info("Running cmd: %s" % cmd)
        subprocess.check_call(cmd, shell=True)

        # run guix package -i to do the actual installation to the profile
        logging.info("Installing packages into the profile")
        cmd = "%s package -p '%s' -i %s" % (guix_exe, profile_path, guix_package_name)
        logging.info("Running cmd: %s" % cmd)
        subprocess.check_call(cmd, shell=True)

    if not args.no_rsync:
        for folder in ['var','profiles','store']:
            # use rsync '-az' i.e. '-rlptgoD' except don't preserve owner, group, or specials or devices.
            cmd = "rsync -rlptz --delete --exclude=.links --exclude=var/guix/userpool --exclude=var/guix/db/big-lock --exclude=var/guix/db/reserved --exclude=var/guix/gc.lock --exclude=var/guix/daemon-socket/socket /RDS/Q0227/%s %s:/data/Q0227/" % (folder, euramoo)
            logging.info("Running %s" % cmd)
            try:
                subprocess.check_call(cmd, shell=True)
            except CalledProcessError, e:
                if e.returncode == 23 and folder == 'store':
                    logging.info("Ignoring returncode 23 as expected when rsync'ing the store")
                else:
                    raise e

        cmd = "ssh %s setfacl --recursive -m m::rwx /data/Q0227/store" % euramoo
        logging.info("Running %s" % cmd)
        subprocess.check_call(cmd, shell=True)

        # It is not entirely clear why this is necessary, why the permissions
        # get changed. But no time to figure it out. Run on awoonga because the
        # links actually resolve, so dead link errors are not thrown.
        cmd = "ssh %s chmod +x '/RDS/Q0227/store/*/bin/* /RDS/Q0227/store/*/bin/.*'" % awoonga
        logging.info("Running %s" % cmd)
        subprocess.check_call(cmd, shell=True)

    logging.info("Finished, it seems")
