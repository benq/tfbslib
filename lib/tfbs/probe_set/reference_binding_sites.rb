module TFBS
  class ReferenceBindingSites
  
    # currently TOMPA only
    def initialize(name, filename)
      @name = name
      @filename = filename
      @datasets=parse
    end
  
    attr_accessor :name
    attr_accessor :file  
  
    def parse
      fh = File.open(filename)
      fh.each |line| do
        if line =~ /^/
      end   
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
  
  end # class ReferenceBindingSites
end # module TFBS  