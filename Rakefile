

require 'rubygems'
require 'tfbs'
require 'erb'
require 'rake/testtask'
require 'rake/packagetask'
require 'rake/gempackagetask'
require 'rake/rdoctask'
#require 'code_statistics'
#require 'metric_fu'


task :default => "test"



Rake::TestTask.new do |t|
  t.libs << 'preload'
  t.test_files = FileList["test/{unit,functional}/**/test_*.rb"]
  t.verbose = true
end



Rake::RDocTask.new do |r|
  r.rdoc_dir = "rdoc"
  main = 'README.rdoc'
  r.main = main
  r.options = ['-N', '-x old_stuff']
  r.options << "--diagram"
end

desc 'print citations'
task :citations do
  require 'tfbs'
  TFBS::Cite.print
  
end  

task :citations => :libs

task :libs do
  require 'tfbs'
end  

