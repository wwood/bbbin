# A loose collection of classes for primer design

class EnzymeCut
  attr_accessor :start, :stop, :enzyme
  def to_s
    [@enzyme, @start, @stop].join("\t")
  end
end

# For design of oligo sequences.

class OligoDesigner
  # Given a sequence, find the subsequence that starts at the 5' end (ie the
  # start of the string), and ends when the melting temperature is maximal but
  # below the max_temperature requires the oligotm program to be available on the
  # cmd line.
  #
  # * nucleotide_string: the full sequence that we are choosing oligos from
  # * max_temperature: the maximal temperature to start things off at
  # * gc_clamp: require this many G or C residues at the 3' end of the oligo.
  def just_below(nucleotide_string, max_temperature, gc_clamp=0)
    # initial conditions
    guess = 0
    guess_temp = 0

    # loop around
    while guess_temp < max_temperature
      guess += 1
      guess_temp = melting_temperature nucleotide_string[0..guess-1]

      # break if there's we've reached the end of the line
      return nucleotide_string if guess > nucleotide_string.length
    end
    return nucleotide_string[0..guess-2]
  end

  # Rank oligomers within some constraints.
  def possible_oligos_ordered_by_temperature_difference(nucleotide_string, min_temperature, best_temperature, max_temperature, gc_clamp)
    default_distance = lambda do |seq, tm|
    # fails constraints if not enough GC clamp
      if seq[seq.length-gc_clamp..seq.length-1].gsub(/[gc]/i,'').length > 0
      false
      else
      # the sequence is within contraints. The melting temperature closest to the best wins.
      tm_diff = (best_temperature-tm).abs
      tm_diff
      end
    end

    # initial conditions
    guess = 0
    guess_temp = 0
    # arrays to fill with possible possibles
    oligos = []

    # loop around, until max temperature is reached
    while guess_temp < max_temperature
      guess += 1
      seq = nucleotide_string[0..guess-1]
      guess_temp = melting_temperature seq

      # Add it to the list if there is enough temperature
      if guess_temp > min_temperature and guess_temp < max_temperature
      o = Oligo.new
      o.sequence = seq
      o.tm = guess_temp
      oligos.push o
      end

      # break if there's we've reached the end of the line
      break if guess > nucleotide_string.length-1
    end

    # Convert sequences into distances
    oligos.each do |oligo|
      oligo.distance = default_distance.call(oligo.sequence, oligo.tm)
    end

    # remove sequences that don't meet the constraints, and sort the rest with
    # smallest distance first
    return oligos.reject{|o| o.distance == false}.sort{|a,b|
      a.distance<=>b.distance
    }.collect{|o| o.sequence}
  end
  alias_method :order, :possible_oligos_ordered_by_temperature_difference

  # A simple method to return the melting temperature of a particular nucleotide
  # string. Uses the command line
  #   oligotm -tp 1 -sc 1 '#{nucleotide_string}'
  def melting_temperature(nucleotide_string)
    `oligotm -tp 1 -sc 1 -n 0.8 -d 500 -mv 0 -dv 50 '#{nucleotide_string}'`.to_f
  end

  private

  class Oligo
    attr_accessor :sequence, :tm, :distance
  end
end

# A class for input and output of Boulder I/O files, for instance as used in
# primer3. This (I don't think) is a full implementation of the Boulder I/O
# format, but serves for my purposes.

class BoulderIO
  include Enumerable
  class Record
    include Enumerable
    attr_accessor :attributes # a hash of key-value pairs
    # Initialise, setting the hash optionally. If no hash is specified, then the
    # attributes will be an empty hash
    def initialize(hash=nil)
      @attributes = hash
      @attributes ||= {}
    end

    def to_s
      ats = @attributes.collect do |key,value|
        "#{key}=#{value}\n"
      end
      "#{ats}=\n"
    end

    # Given a Boulder I/O formatted string, parse into a record
    def self.create(boulder_io_string)
      baby = self.new

      boulder_io_string.split("\n").each do |line|
        line.strip!
        splits = line.split('=')

        # error checking
        if splits.length != 2
          raise ParseException, "Could not parse Boulder I/O line: `#{line}', quitting"
        end

        baby[splits[0]] = splits[1]
      end
      baby #return the new'un for convenience
    end

    # for Enumerable-compatibility
    def each
      @attributes.each do |k,v|
        yield k,v
      end
    end

    # Equivalent to attributes[] - this is a convenient method, isn't it?
    def [](key)
      @attributes[key]
    end

    # Equivalent to attributes[]= - this is a convenient method, isn't it?
    def []=(key,value)
      @attributes[key] = value
    end
  end

  # an array of records
  attr_accessor :records

  # Open a Boulder I/O file and parse it. It is possible to go through all the
  # records:
  #
  #    BoulderIO.open('/path/to/boulderio_file').each {|record| record}
  def self.open(filename)
    record_strings = File.open(filename).read.split("\n=\n")
    @records = record_strings.collect do |s|
      Record.create(s)
    end
    @records # return the records for convenience
  end

  # for Enumerable-compatibility
  def each
    @records.each do |r|
      yield r
    end
  end
end

# Methods for interacting with Primer3 output. Probably overly specific to my
# problems

class Primer3Result
  attr_accessor :output_hash
  def initialize
    @output_hash ||= {}
  end

  def self.create_from_primer3_output_filename(filename)
    unparsed_result = BoulderIO.open(filename)[0]
    baby = self.new
    unparsed_result.each do |key,value|
      baby[key] = value
    end
    baby
  end

  # Was there any primers found? Assumes you were looking for a left primer, with
  # a right primer, which won't always be te case.
  def primer_found?
    @output_hash.keys.include?('PRIMER_LEFT_SEQUENCE')
  end
  alias_method :yeh?, :primer_found?

  # Return the requested part of the result
  def [](key)
    @output_hash[key]
  end

  # set the requested part of the result
  def []=(key,value)
    @output_hash[key] = value
  end
end