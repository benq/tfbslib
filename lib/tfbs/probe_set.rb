module TFBS
  
  # Container for Prediction and Testing
  # requires:
  # * name        without blanks
  # * condition   condition which suggests coregulation
  # * genes/seqs  a list of genes or and sequences
  class ProbeSet
    include HasDataDir
    
    def initialize(name,condition,inputdata=nil,type=:fasta)
      @name=name#.gsub(' ','_')
      @condition=condition
      @inputtype=type
      Dir.mkdir(self.dir)
      # case type
      #  when :fasta
      #    File.new(File.join(self.dir,"input.fasta")) {|f| f.puts inputdata}
      #  when :genelist  
      #    File.new(File.join(self.dir,"input.genelist")) {|f| f.puts inputdata}
      #  when :seqlist
      #    File.new(File.join(self.dir,"input.seqlist")) {|f| f.puts inputdata}
      #  end    
    end
    
    attr_reader  :name
    
    
    # parse input
    # input has to be input.fasta (may implement also input.xls)
    def parse_input
      File.join(self.dir, 'input.fasta')
    end
    

  end # class ProbeSet
  
end # module TFBS 