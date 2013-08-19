module TFBS
  
  # this class is a meta/abstraction for a prediction run, since we mostly do multiple runs
  # it allow the various objects to track folders and options
  # typically one would create an array of runs...
  class Run
  
    def initialize(promoters,dir,predictor,widths)
      @promoters = promoters
      @widths    = widths
      @dir       = dir  
      @predictor = eval("TFBS::#{predictor}.new(promoters,dir,widths)")
      @status    = nil
    end  
  
    # promoter file for the prediction
    attr_accessor :promoters
  
    # array of widths for this run e.g. [6,8,10,12] (default)
    attr_accessor :widths
  
    # base dir usually, there will be a predictor and width subdir
    # e.g. .../predictions/:id/motifsearch/:predictor/:width
    attr_accessor :dir
  
    # currently one of TFBS::PREDICTORS
    attr_reader :predictor
  
    # done, running, failed
    attr_accessor :status
  
    def motifs
      @predictor.parse(:foo=>'bar')
    end
  
    def perform
      @predictor.perform
    
    
      # prediction
      #    parse
      #    weblogo
      #    matches
    end  
  
  end # class Run
end # module TFBS  