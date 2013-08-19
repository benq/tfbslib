module TFBS
  require 'tmpdir'
  require 'rubygems'
  require 'rinruby'    
  
  # while a motiflist is derived from array it is assumed that it operates on motifs
  # thus all elements have to be motifs
  class Motiflist < Array
    include TFBS::MotiflistHelper

    # Use Motiflists to load and use reference motifs (e.g. jaspar, Harbinger)
    # as well a prediction, the Report will expect a Motiflist
    #
    # also use to convert into a different format
    # e.g. Motiflist.new(Jasper.parse(filename)).output(:transfac,filename)
    #
    # allow to inherit from Array
    #def initialize(motifs, source, description)
    #  @motifs=motifs
    #  @source= source
    #  @description= description
    #  
    #end  

    # name for the list
    attr_accessor :name
    
    
    # set an output dir
    def dir=(outputdir)
        @dir = outputdir
        Dir.mkdir(@dir) unless File.exists?(@dir)
    end


    # use a tempdir if none is specified
    # gets cleaned up ;-)
    def dir
      if @dir
        @dir
      elsif @origin
        @dir=File.join(@origin.dir,name)
        Dir.mkdir(@dir) unless File.exists?(@dir)
      else    
        @dir ||=  Dir.mktmpdir
      end
      @dir  
    end

    def sort_by_property(attr)
      self.sort { |a,b|  b.send(attr) <=> a.send(attr) }
    end 
     
    def clip_by_auc(min)
      self.collect { |m| m if m.roc_auc > min }
    end  
     
    def sensible_sort
        #sort_by {|k| k.split(/([-+]?\d+(?:\.\d+)?(?:[-+]?[eE]\d+)?)/).map {|v| Float(v) rescue v}}

    end   
      


    # old        
    def parse(format,in_fn,tag)
      if PROGS.include?(format.to_sym)
        self.send("parse_#{format}", in_fn,tag)
      elsif INPUT_FORMATS.include?(format.to_sym)
        self.send("parse_#{format}", in_fn,tag)
      else
        raise("called unknown input format use one of #{PROGS}")
      end
    end  
  


    # after running the predictions this procedure reads the results to a motiflist
    #
    def parse_motifsearch_dir(p,dataset,i)
      PROGS.each{|m|
        unless m==:weeder 
        WIDTHS.each{ |w|    
          tag={:dataset=>dataset, :dataidx=>i, :method=>m, :width=>w}
          parse(m.to_s,"#{p}/motifsearch/#{m}/#{w}/#{DATAFILES[m]}",tag) 
        }
        else
          tag=tag={:dataset=>dataset, :dataidx=>i, :method=>m}
          parse(m.to_s,"#{p}/motifsearch/#{m}/#{DATAFILES[m]}",tag)
        end    
      }
    end    
  end  # end class Motiflist
  

  module Motiflist::MotiflistHelper
    def names
      n=[]
      self.each{|e| n << e.name}
      n
    end
  end 
end # module TFBS