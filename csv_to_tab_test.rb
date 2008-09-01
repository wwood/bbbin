

require 'test/unit'
require 'csv_to_tab'

class CsvToTabTest <Test::Unit::TestCase
def test_normal
c = CsvToTab.new
assert_equal 'one	two	three', c.parse_line('one,two,three',3,2)
assert_equal 'one	two, yah	blah', c.parse_line('one,two, yah,blah',3,2)

assert_equal 'one	three	two, yah	blah', c.parse_line('one,three,two, yah,blah',4,3)
end
end
