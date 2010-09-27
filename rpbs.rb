#!/usr/bin/env ruby

# A script for automatic creation and setup of pbs scripts, so
# that a script doesn't have to be set up for each job that is to be submitted

require 'tempfile'

class Pbs
  
  DEFAULT_PARAMETERS = {
    :email => 'b.woodcroft@pgrad.unimelb.edu.au',
    :job_name => 'MyJob',
    :nodes => 'nodes=1', #default to serial jobs
    :walltime => '24:00:00',
    :receive_email => 'ae',
    :shell => '/bin/bash'
  }
  
  PARAMETER_CONVERSIONS = {
    :email => 'PBS -M ',
    :job_name => 'PBS -N ',
    :nodes => 'PBS -l ',
    :walltime => 'PBS -l walltime=',
    :receive_email => 'PBS -m ',
    :shell => 'PBS -S '
  }
  
  def self.qsub(command, pbs_parameters = {}, verbose = false)
    # Create a tempfile to create the script in
    Tempfile.open('rpbs_qsub_script') do |tempfile|
      
      # write the command to the tempfile
      contents = get_pbs_script(command, pbs_parameters)
      tempfile.puts contents
      
      puts contents if verbose
      
      tempfile.flush
      
      # qsub in reality
      output = system "qsub #{tempfile.path}"
      unless output
        $stderr.puts "Rpbs: Failed to submit to qsub: #{command}"
        $stderr.puts "Error: #{$?.inspect}"
      end
    end
  end
  
  # Return the string that is to be written to the tempfile that is to be qsub'd
  def self.get_pbs_script(command, pbs_parameters={})
    to_return = "#!/bin/bash\n"
    
    # write each of the parameters to the tempfile
    final_hash = DEFAULT_PARAMETERS.merge pbs_parameters
    final_hash.each do |parameter, endpoint|
      conv = PARAMETER_CONVERSIONS[parameter]
      raise Exception, "RPBS: Don't know how to handle parameter: #{parameter}" if conv.nil?
      
      to_return << "##{conv}#{endpoint}\n"
    end
    
    # some other stuff
    to_return << "# Changes directory to your execution directory (Leave as is)\n"
    to_return << "cd $PBS_O_WORKDIR\n"
    to_return << "# Actual command to run is below\n"
    
    # write the command to the tempfile
    to_return << "#{command}\n"
    
    return to_return
  end
end


# If run as a script, this allows arguments to the script, compared to qsub
# which does not, as far as I know
if $0 == __FILE__
  require 'optparse'
  
  options = {}
  
  USAGE = "rpbs.rb [OPTIONS] [COMMAND]"
  
  OptionParser.new do |opts|
    opts.banner = USAGE
    
    opts.on("-e", "--email EMAIL_ADDRESS", String, "email address to be notified of job updates e.g. 'someone@somewhere.com'. Default: #{Pbs::DEFAULT_PARAMETERS[:email]}") do |v|
      options[:email] = v
    end
 
    opts.on("-j", "--job_name JOB_NAME", String, "name of the job e.g. 'myJob'. Default: #{Pbs::DEFAULT_PARAMETERS[:job_name]}") do |v|
      options[:job_name] = v
    end
    
    opts.on("-n", "--nodes NODES_STRING", String, "nodes string to be used in the pbs script e.g. 'nodes=1:ppn=4'. Default: #{Pbs::DEFAULT_PARAMETERS[:nodes]}") do |v|
      options[:nodes] = v
    end
    
    opts.on("-w", "--walltime WALLTIME", String, "maximum time to be used by the job e.g. '24:00:00' for 1 day. Default: #{Pbs::DEFAULT_PARAMETERS[:walltime]}") do |v|
      options[:walltime] = v
    end
    
    opts.on("-r", "--receive_email EMAIL_NOTIFICATION_CODE", String, "when should the PBS system notify you by email? e.g. 'ae' (send email when job aborts or terminates). Default: #{Pbs::DEFAULT_PARAMETERS[:receive_email]}") do |v|
      options[:receive_email] = v
    end
    
    opts.on("-s", "--shell SHELL", String, "shell to run the script e.g. '/bin/bash'. Default: #{Pbs::DEFAULT_PARAMETERS[:shell]}") do |v|
      options[:shell] = v
    end
  end.parse!
  
  command = ARGV.join ' '
  
  # Be verbose about things
  $stderr.puts "Running command: #{command}"
  $stderr.puts "With options: #{options.inspect}"
  $stderr.puts
  
  Pbs.qsub ARGV.join(' '), options, true
end
