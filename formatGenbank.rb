#!/usr/bin/ruby

f = File.open('unformatted.gb')

f.readlines.each do |line| 
if !line.strip.match(/\:/)

i=0

while i+59>line.length
puts line[i..i+59]
i += 60
end
puts line[i..line.length-1]

else
puts line

end
end

