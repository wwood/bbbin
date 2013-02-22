#!/usr/bin/env ruby

require 'open3'

output_file = File.join(ENV['HOME'],"#{`hostname`.strip}.memory_watcher.log")
out = File.open(output_file,'a')
out.sync = true #disable buffering
user = `whoami`.strip
$stderr.puts "Logging memory usage to #{output_file}"
$stderr.puts "Watching user '#{user}'"
out.puts "Initialised logging of user #{user} at #{`date`.strip}"
out.close

last_memory_usage = nil
while true
  cmd = "ps -U '#{user}' -o pid,%mem,cmd"
  Open3.popen3(cmd) do |stdin, stdout, stderr|
    err = stderr.read
    raise err unless err==''
    velvets = stdout.readlines.collect{|line| line.strip.split(/\s+/).reject{|s| s==""}}
    memories = velvets.collect{|v| v[1].to_f}
    total_memory = memories.inject{|sum,m|sum+=m}
    out = File.open(output_file,'a')
    if last_memory_usage != total_memory
      out.puts `date`; out.puts total_memory
      out.flush
      last_memory_usage = total_memory
    end
    if !memories.empty? and total_memory > 75.0
      out.puts "total memory #{total_memory}"
      out.flush
      
      victim = velvets.max{|v1, v2| v1[1].to_f<=>v2[1].to_f}
      to_kill = victim[0]
      out.puts "Killing the most memory intensive process, pid #{to_kill}, which took #{victim[1]}% of memory: #{victim[2...victim.length]}"
      out.puts `kill -9 #{to_kill}`
      out.flush
    end
    out.close
  end

  sleep 2
end
