# 
# To change this template, choose Tools | Templates
# and open the template in the editor.
 

$:.unshift File.join(File.dirname(__FILE__),'..','lib')

require 'test/unit'
require 'uniq'

class UniqTest < Test::Unit::TestCase
  def test_obvious
    u = Uniq.new
    assert_equal 'two', u.make_unique('two')
    assert_equal 'two-1', u.make_unique('two')
    assert_equal 'three', u.make_unique('three')
  end
end
