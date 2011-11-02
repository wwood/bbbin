
#$:.unshift File.join(File.dirname(__FILE__),'..')

require 'test/unit'

class OrthoMCLFilerTest < Test::Unit::TestCase
  # Path to the script
  def filter
    "#{File.join('.','..','orthomcl_filter.rb')}"    
  end
  
  def test_easy
    # no idea if this stuff works on platforms other than linux
    assert_equal "OG2\n", `echo 'OG4_10036' |#{filter} -g '#{File.join(File.dirname(__FILE__),'data','orthomcl.txt')}' -i pfal,pcha`
  end
end