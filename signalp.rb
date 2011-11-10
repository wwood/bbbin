#!/usr/bin/env ruby

require 'rubygems'
require 'bio'
require 'bio-signalp'

# if this was not called as a module, run as a script.
if $0 == __FILE__
  require 'bio'
  require 'optparse'
  
  runner = Bio::SignalP::Wrapper.new
  
  options = ARGV.getopts("sShvfF") #s for summary, no args required
  if options['h']
    $stderr.puts "Usage: signalp.rb [-svf] <my.fasta>"
    $stderr.puts "Where my.fasta is the name of the fasta file you want to analyse. Default output is all the sequences with their signal sequences cleaved."
    $stderr.puts "-s: summary: print a tab separated table indicating if the sequence had a signal peptide according to the HMM and NN results, respectively."
    $stderr.puts "-S: bigger_summary: like -s, except also includes where the cleavage site is predicted"
    $stderr.puts "-v: verbose summary: much like -s except more details of the prediction are predicted."
    $stderr.puts "-f: filter in: print those sequences that have a signal peptide"
    $stderr.puts "-F: filter out: print those sequences that don't have a signal peptide"
    exit
  end
  
  # Print headers if required
  if options['s']
    puts [
      'Name',
      'NN Prediction',
      'HMM Prediction'
    ].join("\t")
  elsif options['S']
    puts [
      'Name',
      'NN Prediction',
      'HMM Prediction',
      'Predicted?',
      'Cleavege site (if predicted)'
    ].join("\t")

  elsif options['v']
    #       [:nn_Cmax, :nn_Cmax_position, :nn_Cmax_prediction, 
    #      :nn_Ymax, :nn_Ymax_position, :nn_Ymax_prediction, 
    #      :nn_Smax, :nn_Smax_position, :nn_Smax_prediction, 
    #      :nn_Smean, :nn_Smean_prediction,
    #      :nn_D, :nn_D_prediction]
    #    @@hmm_results = [
    #      :hmm_result, :hmm_Cmax, :hmm_Cmax_position, :hmm_Cmax_prediction, :hmm_Sprob, :hmm_Sprob_prediction]
    puts [
      'Name',
      'NN Cmax',
      'NN Cmax position',
      'NN Cmax prediction',
      'NN Ymax',
      'NN Ymax position',
      'NN Ymax prediction',
      'NN Smax',
      'NN Smax position',
      'NN Smax prediction',
      'NN Smean',
      'NN Smean prediction',
      'NN D',
      'NN D prediction',
      'HMM result',
      'HMM Cmax',
      'HMM Cmax position',
      'HMM Cmax prediction',
      'HMM Sprob',
      'HMM Sprob prediction',
    ].join("\t")
  end
  
  Bio::FlatFile.open(ARGV[0]).each do |seq|
    result = runner.calculate(seq.seq)
    if options['s']
      puts [
      seq.entry_id,
      result.nn_D_prediction ? 'T' : 'F',
      result.hmm_Sprob_prediction ? 'T' : 'F'
      ].join("\t")
    elsif options['S']
      puts [
      seq.entry_id,
      result.nn_D_prediction ? 'T' : 'F',
      result.hmm_Sprob_prediction ? 'T' : 'F',
      result.signal? ? 'T' : 'F',
      result.signal? ? result.cleavage_site : 0,
      ].join("\t")
    elsif options['v']
      taputs = [seq.definition]
      [:nn_Cmax, :nn_Cmax_position, :nn_Cmax_prediction, 
      :nn_Ymax, :nn_Ymax_position, :nn_Ymax_prediction, 
      :nn_Smax, :nn_Smax_position, :nn_Smax_prediction, 
      :nn_Smean, :nn_Smean_prediction,
      :nn_D, :nn_D_prediction,
      :hmm_result, :hmm_Cmax, :hmm_Cmax_position, :hmm_Cmax_prediction, 
      :hmm_Sprob, :hmm_Sprob_prediction].each do |meth|
        taputs.push result.send(meth)
      end
      puts taputs.join("\t")
    elsif options['f']
      if result.signal?
        puts seq
      end
    elsif options['F']
      if !result.signal?
        puts seq
      end
    else
      puts ">#{seq.entry_id}\n#{result.cleave(seq.seq)}"
    end
  end
end
