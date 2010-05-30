# To change this template, choose Tools | Templates
# and open the template in the editor.

$:.unshift File.join(File.dirname(__FILE__),'..','lib')

require 'test/unit'
require 'rpbs'

class RpbsTest < Test::Unit::TestCase
  def test_default
    expected = <<DATA
!/bin/bash
#PBS -S /bin/bash
#PBS -M b.woodcroft@pgrad.unimelb.edu.au
#PBS -n MyJob
#PBS -l nodes=1
#PBS -l walltime=24:00:00
#PBS -m ae
# Changes directory to your execution directory (Leave as is)
cd $PBS_O_WORKDIR
# Actual command to run is below
one
DATA
    assert_equal expected.split("\n").sort.join("\n"), Pbs.get_pbs_script('one').split("\n").sort.join("\n")
  end

  def test_change_email
    expected = <<DATA
!/bin/bash
#PBS -S /bin/bash
#PBS -M yeh right
#PBS -n MyJob
#PBS -l nodes=1
#PBS -l walltime=24:00:00
#PBS -m ae
# Changes directory to your execution directory (Leave as is)
cd $PBS_O_WORKDIR
# Actual command to run is below
one
DATA
    puts Pbs.get_pbs_script('one', {:email => 'yeh right'})
    assert_equal expected.split("\n").sort.join("\n"),
      Pbs.get_pbs_script('one', {:email => 'yeh right'}).split("\n").sort.join("\n")
  end
end
