module TFBS
  
  # This is both Taxon and Genome in one class
  # assuming that most data is downloaded from KEGG/NCBI it defaults to KEGG taxon ids
  # Taxon: entry_id, definition
  # Genome: name
  # if you want to add a genome use the name attribute
  # istatus: requested, ..., complete => so we can use is_complete? 
  class Genome 
    include HasDataDir


    def self.kegg
      @@kegg ||= Bio::KEGG::API.new
    end   
    
    def self.organisms
      @@org ||= kegg.list_organisms.collect {|e| [e.entry_id.to_s,e.definition] }
    end  
    
    def definition
      #d=self.class.organisms.collect {|k,v| v if k==@name}.compact.to_s
      #d.gsub!(" \(.*\)","")
      self.class.kegg.binfo(@name).first[17..-1].gsub(" KEGG Genes Database\n",'').gsub(/ \(.*\)/,"")
    end

    # return a hash of kegg taxa
    # use as type ahead for webform
    # def self.update_taxa
    #   
    #     kegg=Bio::KEGG::API.new
    #     kegg.list_organisms.each do |e|  
    #       n = find(:entry_id => e.entry_id.to_s )
    #       if n.size >= 1  
    #         n[0].definition = e.definition
    #         n[0].name       = e.entry_id unless n[0].name 
    #       else
    #         n = Genome.new
    #         n.entry_id   = e.entry_id
    #         n.definition = e.definition
    #         n.name       = e.entry_id unless n.name 
    #         n.save
    #       end    
    #     end
    #   
    # end
    
    # old method to stay compatible is used for lookup and to store backgrounds
    def initialize(name,short='',genome_path='')
      @name=name
    end
    
    attr_reader :name
  
    # # new mothod, check if it is a kegg org
    # # if so grab the genome from ncbi
    # # else raise an name error
    # def initialize(name)  
    #   taxa=self.class.taxa
    #   if taxa.has_key?(name)
    #     @istatus='requested'
    #   else
    #     raise "Name ERROR: there is no kegg taxon with the entry_id #{name} or kegg api not available"
    #   end
    # end  
    #   

    # backward compat
    def short; name end
  
    # path to a genome.gff3 file
    def genome_path
      File.join(dir,"genome.gff3")
    end    
  
    def genome_gb
      File.join(dir,"genome.gb")
    end  
  
    # return a hash of kegg taxa
    # use as type ahead for webform
    def self.taxa
      unless @@taxa
        @@taxa={}
        # TODO KEGG handle to a class method and handle disconnects if this can happen
        srv=Bio::KEGG::API.new
        srv.list_organisms.each{|e| @@taxa[e.entry_id]=e.definition}
      end
      @@taxa  
    end
    
    # download from ncbi and create everything we need to fetch upstream regions
    # http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html
    # dependencies might give this a pretty look
    # we rely on kegg ids and kegg knowing all other ids
    def fetch_genome
        # fetch and save genbank file
        ncbi_es = Bio::NCBI::REST::ESearch.new
        ncbi_ef = Bio::NCBI::REST::EFetch.new
        chromosomes=ncbi_es.genome("#{definition}[Organism]")
        chromosomes.each do |c|
          fh=File.new(File.join(dir,"#{c}.gb"),"w") 
          c_gb = ncbi_ef.sequence(c,'gb',{'Db'=>'genome','retmode'=>'text'})
          puts "writting genome sequence #{c} of #{name}"
          fh.puts c_gb
          fh.close 
        end 
    end
    
    def genes
      Dir["#{dir}/*.gb"].each do |c|
        ff=Bio::FlatFile.auto(c)
        ff.entries
        ff.entry.each_gene{|g| p g['locus_tag']}
      end  
    end 
    
    def save_fasta
      fsa=File.open("#{dir}/genome.fasta","w")
      Dir["#{dir}/*.gb"].each do |c|
        ff=Bio::FlatFile.auto(c)
        ff.entries
        fsa.puts ff.entry.seq.to_fasta
      end
      fsa.close
    end  
    
    # search for upstream region return a Bio:Sequence object with the gene name prefixed by i
    def upstream_of(gene)
      seq='not found'
      Dir["#{dir}/*.gb"].each do |c|
        ff=Bio::FlatFile.auto(c)
        ff.entries
        ff.entry.each_gene do |g| 
          if g['locus_tag'].include?(gene) 
            
            puts pos=g.position
            p m=/(|complement)\(?\<?(\d+)..\>?(\d+)\)?/.match(pos)
            if m[1]=="complement"
              puts "#{m[3].to_i+1},#{m[3].to_i+1001}"
              seq=ff.entry.seq.subseq(m[3].to_i+1,m[3].to_i+1001).complement
            else  
              puts "#{m[2].to_i-1000}..#{m[2].to_i-1}"
              seq=ff.entry.seq.subseq(m[2].to_i-1001,m[2].to_i-1) 
            end  
          end  
        end
      end
      return seq
    end   
      
    
  end # class Genome
end # module TFBS