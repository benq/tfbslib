module TFBS

  # this class is a setup for a run, since we mostly do multiple runs
  # it allow the various objects to track folders and options
  # typically one would create an array of runs...
  class PredictionSet
  
    # creates a new flow that does predictions and then clustering and reports
    # requires
    #   promoter:   path to the file
    # optional
    #   dir:        default = promoterdir/motifsearch
    #               base directory where the prediction data will be placed     
    #   predictors: default = all (see TFBS::PREDICTORS)
    #   widths:     default = [6,8,10,12]
    def initialize(promoters,dir=nil,predictors=nil,widths=[6,8,10,12])
      @promoters  = promoters
      @widths     = widths 
      @dir        = dir || File.join(File.dirname(promoter),'motifsearch')  
      @predictors = predictors || TFBS::PREDICTORS
      @status     = nil

    end  
  
    # rails compatible?
    attr_accessor :id
  

    # Promoters for the prediction
    attr_accessor :promoters
  
    # array of widths for this run e.g. [6,8,10,12] (default)
    attr_accessor :widths
  
    # base dir usually, there will be a predictor and width subdir
    # e.g. .../predictions/:id/motifsearch/:predictor/:width
    attr_accessor :dir
  
    # currently one of TFBS::PREDICTORS
    attr_accessor :predictors
  
    # integrety status
    # final state: complete
    # .*ed (like requested, indexed) is a intermediate state
    # .*ing (like indexing, pareing,..) marks a locked state 
    attr_accessor :istatus
  
    attr_accessor :motifs
  
    attr_accessor :runs
  
    def run
      @runs=predictors.collect {|p| 
        TFBS::Run.new(@promoters, @dir, p, @widths)
      }
    
      @runs.each {|run| run.perform } # prediction, parse, weblogo, matches
      @motifs=TFBS::Motiflist.new()
      @runs.each {|run| 
        run.motifs.each{|m|
          if m != nil
            @motifs << m
          else
            puts "found a nil motif, check your code"    
          end  
        }
      }

      #clusters=Cluster.new(@motifs)
      #@motifs << clusters.collect {|cluster| cluster.motifs}
      #rank
      #prune
      #report
    end  
  
    # just for compatibility of the tests
    def perform
      run
    end
  end # class PredictionSet
end # module TFBS  