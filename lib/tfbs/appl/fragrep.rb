require 'ftools'
require 'tfbs/motif'

class TFBS::Fragrep
  
  
  
    # Create TFBS::Fragrep factory 

    def initialize
      @pattern2eps = "/srv/tfbsdig/db/vendor/fragrep-2/Darwin-i386/pattern2eps"
      @fragrep     = "/srv/tfbsdig/db/vendor/fragrep-2/Darwin-i386/fragrep"
    end
    #   
  
    def pattern_to_logo(pattern_filename)
      dir=File.dirname(pattern_filename)
      name=File.basename(pattern_filename,".*")
      unless File.exists?("#{dir}/#{name}.png")
        `cd "#{dir}"; #{@pattern2eps} "#{pattern_filename}"  2>pattern2eps.err`
        File.rename("#{dir}/pattern2eps.eps","#{dir}/#{name}.eps")
        `cd "#{dir}"; convert "#{name}.eps" #{name}.png 2>convert.err`
      end
      "#{dir}/#{name}.png"
    end
  
    # Options 
    attr_accessor :options
    
    
    # run fragrep
    # fragrep [options] <motif-file> [<genome-file>],
    # 
    # where <motif-file> is a filename as well as <genome-file>.
    # if <genome-file> is not specified, stdin is read instead, which is expected to be in fasta format.
    # valid options are:
    # -u (unweighted) to switch off weighting and optimize for intersection cardinality;
    # -w select weighting based on based on p-values
    # -v verbous output
    # -c scaffold cutoff scoring
    # -r include reverse complemented matches (default for fragrep-0.4 and later)
    # -f do not include reverse complemented matches
    # -s <prefixlength> <suffixlength> to extend extraction of result sequences by a specified number of nucleotides
    # -S use simplified output (leave out aligned query sequences and print reverse sequences for reverse matches)
    #    output produced with -S can be plugged into alignment programs such as clustalw immediately, but is less legible.
    def run_fragrep(infile,pattern_file,outfile)
      unless File.exists?(outfile)
        dir=File.dirname(infile)
        fn= File.basename(infile)
        `cd "#{dir}"; #{@fragrep} -c -r "#{pattern_file}"  "#{fn}"> "#{outfile}"`
      end
    end
    
    # return position of hits
    def parse_fragrep_S(fn)
      #runner.puts "#{fragrep} -c -r #{pp}/#{ntag} #{opposite_promoter} > #{pp}/opposite_promoter_matches/#{ntag}.matches &"
      # count them
      #puts "echo #{ntag}, \"genome: \", `grep matchseq  #{pp}/genome_matches/#{ntag}.matches| wc -l`, \"promoter: \", `grep matchseq  #{pp}/promoter_matches/#{ntag}.matches| wc -l`, \"opposite_promoter: \",`grep matchseq  #{pp}/opposite_promoter_matches/#{ntag}.matches| wc -l` >>count.log"
      # count them
      data=[]
      rec=[]
      fh=File.open(fn)
      fh.each {|line|
        if line =~ /^>(.*)([+-])(\d*)/
           rec=[$1,$2,$3]
        else
           rec << line.strip
           data << rec 
        end  
       }
       data
    end  
    # return position of hits
    def parse_fragrep(fn)
      #runner.puts "#{fragrep} -c -r #{pp}/#{ntag} #{opposite_promoter} > #{pp}/opposite_promoter_matches/#{ntag}.matches &"
      # count them
      #puts "echo #{ntag}, \"genome: \", `grep matchseq  #{pp}/genome_matches/#{ntag}.matches| wc -l`, \"promoter: \", `grep matchseq  #{pp}/promoter_matches/#{ntag}.matches| wc -l`, \"opposite_promoter: \",`grep matchseq  #{pp}/opposite_promoter_matches/#{ntag}.matches| wc -l` >>count.log"
      # count them
      # returns [[chromosome or seq, strand, pos, weight, p-value ]]
      data=[]
      rec=[]
      fh=File.open(fn)
      nx=false
      fh.each {|line|
        if line =~ /^>(.*)-(.*):pos(\d*)\sweight=(.*)\sp-value=(.*)/
          $2=="fwd" ? strand="+" : strand="-"
           rec=[$1,strand,$3,nil,$4,$5]
           nx=true
        else
          if nx
            rec[3]=line.strip
            data << rec 
            nx=false
          end  
        end  
       }
       data
    end
 
end