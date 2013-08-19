module TFBS
  
  # the initial idea was to use BioSQL but I had to split the UI (rails based) from the backend so I could run on different boxes
  # so I extracted what I need to reference the right genome
  class Taxon
  
    # Taxon is used for lookup and to store backgrounds
    def initialize(name,short,genome_path)
      @name=name
      @short=short
      @genome_path=genome_path
    end
  
    # taxon name  
    attr_accessor :name
  
    # short name e.g. scer
    attr_accessor :short
  
    # path to a genome.gff file
    attr_accessor :genome_path  
    
    # fragrep does not like gff
    def genome_fasta
      File.join(File.dirname(genome_path),"#{File.basename(genome_path,".*")}.fasta")
    end  
    
    # obviously, the ppas.gff is a fasta file
    def genome_gff3
      File.join(File.dirname(genome_path),"#{File.basename(genome_path,".*")}.gff3")
    end

    #format {:gene_id=>[seq_id,rec.start,rec.end,rec.strand,rec.feature],..}
    def genes
      g={}
      gff=Bio::GFF::GFF3.new(File.read(genome_gff3))
        
      gff.records.each do |rec|
        
        if rec.attributes.length > 0
          rec_id = nil
          rec_name = nil
          begin
            (k, rec_id) = rec.attributes.assoc("ID")
            (k, rec_name) = rec.attributes.assoc("Name")

            # as GFF3 is not everywhere available I will make a bold assumtion
            # blame me if I am wrong, I will assume the first attribute record holds the ID
            # whether it is called ID, Name, gene_id,... in gff2/gtf
            #
            (k, rec_id) = rec.attributes.first unless rec_id
            (k, rec_name) = rec.attributes.first unless rec_name

          rescue
            errors.add(:Organism_parse_feature_records, "attributes is not an array!!")
            (k, rec_id) = rec.attributes.first unless rec_id
            (k, rec_name) = rec.attributes.first unless rec_name
          end

          #another fix for gff2/gtf format
          # drop anythig after a space as some annotations use seqname like
          # supercontig of Candida lustinae
          seq_id = rec.seqname.scan(/[\.\w]+/)[0] if rec.methods.include?("seqname")
          seq_id = rec.seqid.scan(/[\.\w]+/)[0] if rec.methods.include?("seqid")
          if rec.feature == "gene" or rec.feature == "CDS"
            g[rec_id]=[seq_id,rec.start,rec.end,rec.strand,rec.feature] 
          end  
        end
      end  
      return g
    end
    
    #use like {"ROT2" => "RPPA08267","SEC31"=> "RPPA06211","UBC1"=>"RPPA09581","ARG1"=>"orf_11738","ACS1"=>"RPPA07570","GCN4"=>"RPPA07905","TRR1"=>"RPPA06201","CUP5"=>"RPPA09067","GCR1"=>"RPPA06205","BMH2"=>"RPPA07190"}.each{|k,x| puts ">#{x}|ppas|#{k}\n#{u[x][6]}"}
    
    def upstream
      g={}
      gff=Bio::GFF::GFF3.new(File.read(genome_gff3))
      q={}
      gff.sequences.each{|s| s.na; q[s.entry_id] = s }
      gff.records.each do |rec|
        
        if rec.attributes.length > 0
          rec_id = nil
          rec_name = nil
          #begin
            (k, rec_id) = rec.attributes.assoc("ID")
            (k, rec_name) = rec.attributes.assoc("Name")

            # as GFF3 is not everywhere available I will make a bold assumtion
            # blame me if I am wrong, I will assume the first attribute record holds the ID
            # whether it is called ID, Name, gene_id,... in gff2/gtf
            #
            #(k, rec_id) = rec.attributes.first unless rec_id
            #(k, rec_name) = rec.attributes.first unless rec_name

          #rescue
          #  errors.add(:Organism_parse_feature_records, "attributes is not an array!!")
          #  (k, rec_id) = rec.attributes.first unless rec_id
          #  (k, rec_name) = rec.attributes.first unless rec_name
          #end

          
          #another fix for gff2/gtf format
          # drop anythig after a space as some annotations use seqname like
          # supercontig of Candida lustinae
          seq_id = rec.seqname.scan(/[\.\w]+/)[0] if rec.methods.include?("seqname")
          seq_id = rec.seqid.scan(/[\.\w]+/)[0] if rec.methods.include?("seqid")
          if rec.feature == "gene" or rec.feature == "CDS"
            if rec.strand =="+"
              upstream = q[seq_id][rec.start-1000,1000]
            else
              upstream = q[seq_id][rec.end,1000].complement
            end  
            g[rec_id]=[seq_id,rec.start,rec.end,rec.strand,rec.feature,rec_name,upstream] 
          end  
        end
        
      end
      return g
    end   
  
  end # class Taxon
end # module TFBS  