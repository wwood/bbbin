
#$:.unshift File.join(File.dirname(__FILE__),'..')

require 'test/unit'
require 'tempfile'

class DNDS2SeqsTest < Test::Unit::TestCase
  # Path to the script
  def jumper
    "#{File.join('.','..','dNdS_2seqs.rb')}"    
  end
  
  def test_1
    # no idea if this stuff works on platforms other than linux
Tempfile.open('testDNDS') do |tempfile|
tempfile.puts <<EOF
>one
atgatgatg
>two
atgatgatg
EOF
tempfile.close
`cat #{tempfile.path}`
answer = `#{jumper} #{tempfile.path}`
    assert_equal "dN\tdS\tdN/dS\n0\t0\tNaN\n", answer
end


Tempfile.open('testDNDS') do |tempfile|
tempfile.puts <<EOF
>one
atgatgatt
>two
atgatgatg
EOF
tempfile.close
`cat #{tempfile.path}`
answer = `#{jumper} #{tempfile.path}`
    assert_equal "dN\tdS\tdN/dS\n1\t0\tInfinity\n", answer
end
end

def test_2

Tempfile.open('testDNDS') do |tempfile|
tempfile.puts <<EOF
>one
atgatgattgaa
>two
atgatgatggag
EOF
tempfile.close
`cat #{tempfile.path}`
answer = `#{jumper} #{tempfile.path}`
puts answer
    assert_equal "dN\tdS\tdN/dS\n1\t1\t0.5\n", answer
end

  end
end
