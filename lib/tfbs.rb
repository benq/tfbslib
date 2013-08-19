
module TFBS
  # make sure we can use this
  LIB_ROOT = File.dirname(__FILE__)
  require 'rubygems'
  require 'bio'
  require "rinruby"
  #require 'neo4j'
 
  # require 'neo4j/extensions/reindexer'
  require 'tfbs/config/defaults'
  
  autoload :Cite,           'tfbs/config/citations'
  autoload :HasDataDir,     'tfbs/config/has_data_dir'
  autoload :EachTx,     'tfbs/config/each_tx'
  
  #require 'tfbs/config/citations'
  
  # common for background creation and analysis tools

  autoload :Genome,         'tfbs/genome' 
  
  # probeset
  autoload :ProbeSet,       'tfbs/probe_set'
    autoload :Taxon,          'tfbs/probe_set/taxon' 
  autoload :Promoters,      'tfbs/probe_set/promoters'
  
  # results
  autoload :Motif,          'tfbs/motif'
 
  autoload :Motiflist,      'tfbs/motiflist'
  
  # prediction tools
  autoload :Weeder,         'tfbs/appl/weeder'
  autoload :Bioprospector,  'tfbs/appl/bioprospector'
  autoload :Meme,           'tfbs/appl/meme'
  autoload :Motifsampler,   'tfbs/appl/motifsampler'
  

  
  # run predictions 
  autoload :PredictionSet,    'tfbs/prediction_set'  
  autoload :Run,            'tfbs/prediction_set/run' 

  
  # analysis tools
  autoload :Fragrep,        'tfbs/appl/fragrep' 
  autoload :Stamp,          'tfbs/appl/stamp'
  
  # helpers 
  autoload :MotiflistHelper,          'tfbs/motiflist/motiflist_helper' 
  
  attr_accessor :motifs
  
  # constants for installed components
  PREDICTORS=[ :Weeder,:Bioprospector,:Meme,:Motifsampler]
  PROGS=[:bioprospector, :meme, :motifsampler, :weeder]
  INPUT_FORMATS=[:tamo]
  OUTPUT_METHODS=[:transfac,:jasper,:fragrep]
  
  #defaults
  GENOME="/srv/tfbsdig/data/organisms/1/1-upload.gff"
  ALPHA=[:A,:C,:G,:T]
  BACKGROUND_FREQ={:A=>0.293976, :C=>0.205515, :G=>0.20598, :T=>0.294529}      
  DATAFILES={:weeder=>"out.wee", :meme=>"meme.txt", :bioprospector=>"out.txt", :motifsampler=>"out.mtrx"}   
  WIDTHS = [6,8,10,12]
  # motifs to discover per run
  MPR=10
  

 
  
end # module TFBS