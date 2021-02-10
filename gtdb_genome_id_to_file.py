#!/usr/bin/env python3
import logging
import argparse
import csv
import re
import os

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information',
                               action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors',
                               action="store_true")

    parent_parser.add_argument('-i','--ids', help='input file of NCBI/GTDB/UBA IDs. UBA IDs can be in U_123 or UBA123 form', required=True)
    parent_parser.add_argument('-r','--release',help='GTDB release [default 95]',default='95')
    parent_parser.add_argument('--file-type',help='which file from each genome [default fna]', default='fna')
    args = parent_parser.parse_args()

    input_accessions = args.ids

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    if args.release=='95':
        genome_dirs = '/srv/db/gtdb/genomes/ncbi/release95/genome_dirs.tsv'
    elif args.release=='89':
        genome_dirs = '/srv/db/gtdb/genomes/ncbi/release89/genome_dirs.tsv'
    else:
        raise Exception("Unexpected gtdb release specified")
    release = args.release
    uba_mapping_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), 'uba_ncbi_accessions.tsv' )

    uba_genomes_base_dir = '/srv/db/gtdb/genomes/user/uqdparks/'

    id_to_folder = {}
    uba_accession_to_u = {}
    gtdb_accession_to_u = {}

    logging.info("Reading genome_dirs file ..")
    with open(genome_dirs) as f:
        csv_reader = csv.reader(f, delimiter="\t")
        for row in csv_reader:
            if (release=='89' and len(row) != 2) or (release=='95' and len(row) != 3):
                raise Exception("Unexpected genome_dirs row found: {}".format(row))

            id_to_folder[row[0]] = row[1].strip()

    logging.info("Reading UBA accession mapping file ..")
    with open(uba_mapping_file) as f:
        csv_reader = csv.reader(f, delimiter="\t")
        for row in csv_reader:
            if len(row) != 3:
                raise Exception("Unexpected UBA mapping file row found: {}".format(row))

            if row[2].strip() != '':
                gtdb_accession_to_u[row[2].strip()] = row[0].strip()
            uba_accession_to_u[row[1].strip()] = row[0].strip()

    logging.info("Converting IDs ..")
    with open(input_accessions) as f:
        for genome_id in f:
            genome_id = genome_id.strip()
            q = genome_id
            if q in id_to_folder:
                folder = id_to_folder[q]
            else:
                q2 = re.sub('^GB_','',re.sub("^RS_",'',q))
                if q2 in id_to_folder:
                    folder = id_to_folder[q2]
                    q = q2 # required for prodigal_faa type
                elif q in gtdb_accession_to_u:
                    folder = os.path.join(uba_genomes_base_dir, gtdb_accession_to_u[q])
                elif q in uba_accession_to_u:
                    folder = os.path.join(uba_genomes_base_dir, uba_accession_to_u[q])
                else:
                    raise Exception("Unable to find accession {}".format(genome_id))
            stub = folder.split('/')[-1]
            if args.file_type =='fna':
                print("\t".join([genome_id, "{}/{}_genomic.fna".format(folder,stub)]))
            elif args.file_type == 'prodigal_faa':
                print("\t".join([genome_id, "{}/prodigal/{}_protein.faa".format(folder,q)]))
            else:
                raise Exception("bad file type")
