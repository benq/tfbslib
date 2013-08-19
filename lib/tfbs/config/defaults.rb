module TFBS
  TEST_ROOT=File.expand_path(File.join(File.dirname(__FILE__),'..','..','..','test'))
  METHODS_ROOT=File.expand_path('~/devel/tfbsdig/methods')
  DATA_ROOT=File.expand_path('~/FH/tfbs-data')
  
  #   require 'neo4j-rails'
  # Neo4j::Config[:storage_path] = File.expand_path(File.join(DATA_ROOT,'db'))
  # Lucene::Config[:store_on_file] = true
  # Lucene::Config[:storage_path] = File.expand_path(File.join(DATA_ROOT,'db','lucene'))
  # require 'rails'
  # #Rails.application='tfbs'
  # require "action_controller/railtie"
  # 
  # class Application < Rails::Application
  # 
  # end

  # TFBS::Application.configure {|app|
  #   app.config.neo4j.storage_path = File.expand_path(File.join(DATA_ROOT,'db'))
  #   app.config.lucene.storage_path = File.expand_path(File.join(DATA_ROOT,'db','lucene'))
  #   app.config.lucene.store_on_file = true
  # }

  #Neo4j.start
  

  
    
  def self.log(msg)
    File.open(File.join(DATA_ROOT, "tfbs.log"),"a"){|f| f.puts(msg)}
  end  
  
  # this has to change thus we wrap them at the moment
  def self.lib_root
    LIB_ROOT
  end  
  
  def self.data_root
    DATA_ROOT
  end
  
end #module TFBS  