require 'test/unit'
require 'enzyme_matcher'

class EnzymeNameMatcherTest < Test::Unit::TestCase
	def test_simple
		e = EnzymeMatcher
		assert_equal true, e.match?('EcoRI','EcoRI'), "EcoRI vs. EcoRI"
		assert_equal true, e.match?('EcoRI','EcoRI,'), "EcoRI vs. EcoRI,"
		assert_equal false, e.match?('EcoRI','GrEcoRI'), "EcoRI vs. GrEcoRI"
	end
end