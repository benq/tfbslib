module TFBS
  class Promoters
  
    # wraps up all methods todo with promoters
    #   lookup of upstream seq or genelist if seq is provided
    #   lookup of homologues and creation of a sequenceset for weederH and Phylocon
    def initialize(name, file, taxon, genelist=nil)
      @name = name
      @file = file
      @taxon= taxon
      @genelist= genelist
    end
  
    attr_accessor :name
    attr_accessor :file  
    attr_accessor :taxon  
    attr_accessor :genelist
  
    def clean_file
      @clean_file=File.join(dir,"promoters_of_#{ename}_cleaned.fa")
      unless File.exists?(@clean_file)
        # clean upstream.fa 
        # TODO this should be improved removing sequences shorter than 8 bp
        `grep -v '#' "#{@file}" |grep -v '^$' > "#{@clean_file}"`
      end
      @clean_file  
    end  
  
    # allows us to put the motif files there as a default (jasper, weblogos,..)
    def dir
      @dir ||= File.dirname(@file)
    end  
  
    # escaped name, we dont need spaces,... meme cant handle them for one
    def ename
      name.gsub(" ","_")
    end  
  
    # return 3 randomized promoter files
    # create if neccessarz
    def randomized_promoters
    
    
    end
    
    def genelist
      pr=Bio::FlatFile.auto(clean_file)
      dl=pr.collect {|e| e.definition.split('|')}
    end  
  
  end # class Promoters
end # module TFBS    