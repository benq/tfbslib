module TFBS

CITATIONS={ 'Pavesi et al., 2004'=> 15215380,
            'Bailey and Elkan, 1994'=> 7584402
  }

  class Cite

    def self.citations
      CITATIONS
    end

    def self.print(list=TFBS::CITATIONS)
      puts
      list.each do |k,v|
        #puts "#{k}: #{v}"
        puts "[#{k}]: #{Bio::MEDLINE.new(Bio::PubMed.query(v)).reference.general}"
        puts 
      end
    end  

    def self.endnote(list=TFBS::CITATIONS)
      list.each do |k,v|
        puts "#{k}: #{Bio::MEDLINE.new(Bio::PubMed.query(v)).reference.endnote}"
      end   
    end
          
  end
  
end