#!/usr/bin/env ruby

if __FILE__ == $0
  # E.g. 
  # V,K,W,L,A,M,Y,N,C,D,P,E,Q,F,R,G,S,H,T,I,U
  # 96912,320297,11384,204768,45673,60002,158248,412149,46587,175730,47622,184987,72767,118590,69392,67773,171322,66871,104377,255389,0

  lines = ARGF.read.split("\n")
  aa_names = lines[0].split(",")
  numbers = lines[1].split(',').collect{|a| a.to_f}
  raise unless aa_names.length == numbers.length
  
  total = 0.0
  numbers.each {|n| total += n}
  percentages = numbers.collect{|n| n/total}
  
  steps = []
  percentages.each_with_index do |n,i|
    if steps.empty?
      steps.push n
    else
      steps.push n+steps[i-1]
    end
  end
  
  (1..10000).each do |iter|
    length = 100
    seq = 'M' #always start with a methionine
    while seq.length < length
      r = rand
      steps.each_with_index do |step, i|
        # puts [seq, r, step, i].join("\t") #debug
        if r < step
          seq = seq+aa_names[i]
          break
        end
      end
    end
    puts ">#{iter}"
    puts seq
  end
end