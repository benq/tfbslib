require 'rubygems'
require 'erb'
require 'rake'
require 'bio'
require 'tfbs'

namespace :genomes do 
  
  desc "dummy task"
  task :list do
    puts 'complete genomes:'
    puts 
    puts complete
    
    puts 'pending'
    puts fresh
    
  end
  
  def fresh
    FileList['genomes/*']
  end
  
  def complete 
    list.collect{|e| e if File.exists?(File.join(e,'genome.fasta'))}
  end
  
  def list
    FileList['genomes/*']
  end  
    
end



namespace :taxa do 
  
  desc "fetch taxa ncbi ids from ncbi"
  file "genomes/taxdump/names.dmp" => ["genomes/taxdump","genomes/taxdump.tar.gz"] do
    `tar xzf genomes/taxdump.tar.gz -C genomes/taxdump`
  end
  
  file "genomes/taxdump" do 
    mkdir "genomes/taxdump"
  end
  
  file "genomes/taxdump.tar.gz" do
    `wget --directory-prefix=genomes ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`
  end

  desc "fetch taxa id from kegg"
  task "kegg_ids" do
    kegg=Bio::KEGG::API.new
    kegg.list_organisms.each do |e|
        n = find(:entry_id => e.entry_id.to_s )
        if n.size >= 1  
          n[0].definition = e.definition
        else
          n = Taxon.new
          n.entry_id   = e.entry_id
          n.definition = e.definition
        end  
    end      
  end
  
  desc "fetch taxa from kegg and update db"
  task :update do
    TFBS::Taxon.update_taxa
  end

  desc "list taxa in db"
  task :list do
   taxa=TFBS::Taxon.all.nodes
   taxa.each {|t| p t }
  end  

  task :update_index do
    Neo4j::Transaction.run do
      TFBS::Taxon.update_index
    end
  end
end



 
desc "list taxa"
task :list_taxa do
  FileList['genomes/*'].each {|taxon| puts taxon }
end



task :list_taxa => :genomes

task :load_tfbs do
  
end  