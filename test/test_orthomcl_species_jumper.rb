
#$:.unshift File.join(File.dirname(__FILE__),'..')

require 'test/unit'

class OrthoMCLSpeciesJumperTest < Test::Unit::TestCase
  # Path to the script
  def jumper
    "#{File.join('.','..','orthomcl_species_jumper.rb')}"    
  end
  
  def test_groups_input
    # no idea if this stuff works on platforms other than linux
    assert_equal "OG4_10036\n", `echo 'OG4_10036' |#{jumper}`
    assert_equal "OG4_10036\tpfal|PFF0510w,pfal|PFF0865w\n", `echo 'OG4_10036' |#{jumper} -o pfal`
  end
end