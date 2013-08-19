require 'rake'  


taskset = Dir.entries(pwd) & ['genomes', 'probesets', 'predictionsets']

taskset.each do |t|  
  # prevent the task from being run multiple times.
  unless Rake::Task.task_defined? t
    # Load the rakefile so users of the gem get the default metric_fu task
    load File.join(TFBS::LIB_ROOT, 'tasks', "#{t}.rake")
  end
end

task :default => :list_tasks

#p display_tasks_and_comments

 

task :list_tasks do |t|
  puts 'Hi'
end