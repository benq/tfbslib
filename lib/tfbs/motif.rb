 #
 # = tfbs/motif.rb - Motif classes
 #
 # Copyright::  Copyright (C) 2009-2010 Benjamin Almeida <benjamin.q.almeida@gmail.com>,
 #
 # License::    The Ruby License
 #
 # == Description
 #
 # This file contains classes that implement a TFBS Motif as created by various prediction tools
 #
 # == References
 #
 # * Hoge F. et al., The Hoge database, Nucleic. Acid. Res. 123:100--123 (2030)
 # * http://hoge.db/
 #

module TFBS
  
  
  
  class Motif


    autoload :Jasper,          'motif/jasper'
   # autoload :Fragrep,         'motif/fragrep'
    autoload :Transfac,        'motif/transfac'
    autoload :Tamo,            'motif/tamo'
    autoload :Parser,          'motif/parser'
    autoload :Metrics,          'tfbs/motif/metrics'

    # Motif class
    # this should relieve the pain of all the different output formats that prediction tools create
    # inputs:
    # * matrix = takes PWM/PFM as input (if sum over ACGT > 1 it is a PFM otherwise PWM)
    # * tag = just a list to keep track of where data is from, what prediction tool,...
    def initialize(matrix,tag=[])
      @tag=tag
      @iupac_consensus=nil
      if matrix[0].values.inject{|a,b| a.to_f+b.to_f} > 1.05 # be round error friendly a freq matr is much bigger
        @matrix=matrix
      else  
        @prob_matrix=matrix
      end
      
    end


    ### matrix data stuff

    attr_accessor :matrix
    attr_accessor :tag
    attr_accessor :rank
    
    # return simplified consensus 
    def consensus;  @consensus ||= "no_consensus" end
    def rank;       @rank ||= "no_rank" end
    
    def n=(i);      @tag[:n]=i end




    ### linking with other objecs..

    # parents
    attr_accessor :promoters


    # parent
    attr_accessor :motiflist

    # children
    attr_accessor :motifs

    attr_accessor :predictor

    def taxon
      if @taxon
        @taxon
      elsif promoters
        promoters.taxon
      else
        raise 'a motif needs eighter a promoter or taxon set to reference use some functions'
      end
    end        

    def taxon=(arg)
      @taxon=arg
    end  

    def width
      matrix.size
    end  
    
    # set an output dir
    def dir=(arg)
      @dir = arg
    end

    # replace the use of ntag
    def name
      if @name
        @name
      else  
        if motifs != nil
          @name="cluster_#{iupac_consensus}_#{predictor.appl_name}"
        elsif promoters
          @name="#{promoters.ename}_#{iupac_consensus}_#{predictor.appl_name}_#{width}_#{rank}"
        else
          @name="unknown_motif_type_#{iupac_consensus}"   
        end
      end      
    end  

    #allow to set name if we import a list (e.g. harbison, jasper2010)
    def name=(n)
      @name=n
    end  

    # use a tempdir if none is specified
    # gets cleaned up ;-)
    def dir
      if @dir 
        @dir
      else
        if motiflist
          @dir=motiflist.dir
        elsif promoters
          #not the best option as different prediction runs would overwrite
          # thus recent_motifs is a hint, that this is not to keep
          @dir=File.join(promoters.dir, 'motifs', name )
        elsif @origin
          @dir=@origin.dir
        else
          @dir=Dir.mktmpdir
        end
      end  
      unless File.exist?(@dir)
        Dir.mkdir(@dir)
      end
      return @dir  
    end


    ### queries for more data
    ### many do external calls

    # run fragrep and make a gff of all hits
    #
    def gff
      "TODO"
    end  

    def rand_hits
      "TODO"
    end   
  
  
  
    def hitlist
      data=""
      #hl =~ /^>(\w+)|ppas|(\w+)-\w+:pos\d+ .*/
      hl.each{|line| line =~ /^>(\w+)\|ppas\|(.+)/ ; data<< "#{$1.rjust(10)}(#{$2.center(10)})\t" }
      return data
    end  

    # accessor handles missing matrix by creating it from the opposite matrix
    # position probability matrix
    # corrected probability as in WW Waserman, A Sandelin 2004
    # p(b,i)=(f(b,i)+s(b))/N+sum(s(b'))
    #         f(b,i) .. counts of base b in position i
    #         b' = [A,C,G,T]
    #         N .. number of sites
    def prob_matrix
      unless @prob_matrix
        @prob_matrix= Marshal.load(Marshal.dump(@matrix))
        @prob_matrix.each{|pos|
          sum=pos.values.inject{|a,b| a.to_f+b.to_f} 
          [:A,:C,:G,:T].each{|nuc| pos[nuc]=pos[nuc].to_f+pseudo([nuc])/sum+pseudo([:A,:C,:G,:T])} 
        }
      end
      return @prob_matrix
    end
    
    # position specific score matrix .. PSSW or PWM
    # as in WW Waserman, A Sandelin 2004
    # W(b,i) = log2(p(b,i)/p(b))
    # TODO change from hash to array
    def pssm_matrix
      unless @pssm_matrix
        @pssm=Marshal.load(Marshal.dump(prob_matrix))
        @pssm.each{|w|
          [:A,:C,:G,:T].each{|nuc| w[nuc]= Math.log(w[nuc]/genome_p(nuc))/Math.log(2)}
        }
      end      
    end
    
    # genome background probability for base
    # TODO this should use the background of the genome
    # currently yeast
    def genome_p(arg)
      if TFBS::ALPHA.include?(arg.upcase)
        TFBS::BACKGROUND_FREQ[arg.upcase]
      else
        raise "ERROR: wrong key need ACGT or 0..3"
      end    
    end  
    
    # pseudocount function
    # currently 0
    def pseudo(arg)
      0
    end  

    # accessor handles missing matrix by creating it from the opposite matrix
    # position frequency matrix ... PFM
    def matrix
      unless @matrix
        @matrix= Marshal.load(Marshal.dump(@prob_matrix))
        min=@matrix.collect{ |pos| pos.values.min.to_f }.min
        min=0.1 if min == 0.0 # if the smalest is 0 we use 0.1 
        @matrix.each{|pos| [:A,:C,:G,:T].each{|nuc| pos[nuc]=pos[nuc].to_f/min} 
        }
        max=@matrix.collect{ |pos| pos.values.max.to_f }.max
        if max > 100 then 
          # reduce fragrep load
          @matrix.each{|pos| [:A,:C,:G,:T].each{|nuc| pos[nuc]=pos[nuc]-1} }
        end  
      end
      return @matrix
    end

    # create iupac_consensus 
    # TODO use bioruby to create these
    #
    # M = A or C (i.e. A=0.5 C=0.5 G=0 T=0)
    # R = A or G
    # W = A or T
    # S = C or G
    # Y = C or T
    # K = G or T
    # V = not T (i.e. A=1/3 C=1/3 G=1/3 T=0)
    # H = not G
    # D = not C
    # B = not A
    # N = A=C=G=T=0.25
    # use similar rules as used with stamp
    # probability thresholds:
    # A/C/G/T is used if the appropriate single base frequency is >0.6
    # M/R/W/S/Y/K is used if the sum of the appropriate two bases is >0.8
    #
    def iupac_consensus 
      unless @iupac_consensus
        cons=prob_matrix.collect{|pos|
          p=pos.sort{|a,b| b[1].to_f<=>a[1].to_f} # sort desending by value
          c="_"
            if p[0][1].to_f > 0.6
              c=p[0][0]
            elsif p[0][1].to_f+p[1][1].to_f > 0.8
              k=[p[0][0].to_s,p[1][0].to_s].sort.join
              c=case k
              when "AC" then "M"
              when "AG" then "R"
              when "AT" then "W"
              when "CG" then "S"
              when "CT" then "Y"
              when "GT" then "K"
              end
            elsif p[3][1].to_f==0.0
              k=p[3][0]
              c=case k
              when "A" then "B"
              when "C" then "D"
              when "G" then "H"
              when "T" then "V"
              end
            else
              c="N"                            
            end
            c   
          }
          @iupac_consensus=cons.join  
    
      end
      return @iupac_consensus  
    end  


      ### output stuff
    def to_jasper
      render_to_text(:jasper)
    end  

    def to_transfac
      render_to_text(:transfac)
    end


    # make a record
    def render_to_text(format)
      data=''
      case format
      when :list
        #data << "#{tag[:n].to_s.ljust(4)} #{name.ljust(22)} #{consensus.ljust(12)} #{iupac_consensus.ljust(12)} #{"%8.3g" % bla_quot.to_s} | #{"%8.3g" % a} #{"%8.3g" % b} #{"%8.3g" % c} #{"%8.3g" % d} #{"%8.3g" % e} #{"%8.3g" % f} #{"%8.3g" % g} #{"%8.3g" % h.to_f} #{"%8.3g" % l.to_f} | #{hl.scan(/\n/).length} of #{DATALENGTH[tag[:dataidx]].to_i} #{tag[:dataset].ljust(18)} #{hits.inspect}"
        data << "#{tag[:n].to_s.ljust(4)} #{name.ljust(22)} #{consensus.ljust(12)} #{iupac_consensus.ljust(12)} "
      when :hitlist
        data << hitlist
      when :jasper
        data << "> #{iupac_consensus}_#{name}-#{consensus}\n"
        data << [:A,:C,:G,:T].collect{|nuc| "#{nuc} [#{matrix.collect{|pos| pos[nuc]}.join(' ')}]"}.join("\n")+"\n"
      when :transfac
          #data << "NA\t#{name}\n"
          #data << "XX\n"
          data << "DE\t#{name}\n"
          #data << "XX\n"
          data << "P0\tA\tC\tG\tT\n"
          i=0
          data << matrix.collect{|v|
            ((i+=1)<11) ?  frow="0".to_s+(i-1).to_s : frow="".to_s+(i-1).to_s
            [frow,[:A,:C,:G,:T].collect{|nuc| sprintf("%0.3g",matrix[i-1][nuc])} ].join"\t"}.join("\n") 
          data << "\n"
          data << "XX\n"        
      when :tree_tf
              #data << "NA\t#{iupac_consensus}\n"
              #data << "XX\n"
              data << "DE\t#{iupac_consensus}\n"
              #data << "XX\n"
              data << "P0\tA\tC\tG\tT\n"
              i=0
              data << matrix.collect{|v|
                ((i+=1)<11) ?  frow="0".to_s+(i-1).to_s : frow="".to_s+(i-1).to_s
                [ frow,[:A,:C,:G,:T].collect{|nuc| "#{matrix[i-1][nuc]}"}].join("\t") }.join("\n") 
              data << "\n"
              data << "XX\n"
      when :fragrep
        #1 matrices
        #    0     0 M0:UCGCGC 0.83 0
        #M0:UCGCGC
        #   1    1    1    1    1    1 
        #   1   16    1   16    1   16 
        #   1    1   16    1   16    1 
        #  16    1    1    1    1    1 
        ##  U    C    G    C    G    C   
        data << "1 matrices\n"
        data << "    0     0 M0:#{iupac_consensus} 0.9 0\n"
        data << "M0:#{iupac_consensus}\n"
        data << [:A,:C,:G,:T].collect{|nuc| "   #{matrix.collect{|pos| pos[nuc].to_i}.join("\t")}"}.join("\n")+"\n"
        data << "#  " + iupac_consensus.split('').join("\t") + "\n"
      end
      return data
    end

    # create fragrep file and return file name
    def render_fragrep
      format=:fragrep
        Dir.mkdir(dir) unless File.exists?(dir)
        fn="#{dir}/#{name}.pat"
        fh=File.new(fn,"w+")
        fh.puts render_to_text(:fragrep)
        fh.close
      fn  
    end

    # create weblogo and return the path
    def weblogo
      Dir.mkdir(dir) unless File.exists?(dir)
      f=TFBS::Fragrep.new
      f.pattern_to_logo(render_fragrep)
    end

    # run fragrep and parse output 
    def fragrep
      Dir.mkdir("#{dir}/matches/") unless File.exists?("#{dir}/matches/")
      f=TFBS::Fragrep.new
      fra="#{dir}/matches/#{name}.prom.fragrep"
      f.run_fragrep(promoters.clean_file, render_fragrep, fra )
      f.parse_fragrep(fra)
    end  

    # run fragrep to index all motif hits in the genome
    def full_genome_fragrep
      unless name =~ /TOMPA/
        Dir.mkdir("#{dir}/matches/") unless File.exists?("#{dir}/matches/")
        f=TFBS::Fragrep.new
        fra="#{dir}/matches/#{name}.full_genome.fragrep"
        f.run_fragrep(promoters.taxon.genome_fasta, render_fragrep, fra ) unless File.exists?(fra)
        f.parse_fragrep(fra)   
      end  
    end
  


    
    
    
    def metrics
      fn=File.join(dir,'matches',"#{name}_metrics.yaml")
      TFBS.log(fn)
      if @metrics !=nil
        return @metrics
      #elsif File.exist?(fn)
      #  @metrics = YAML::load_file(fn)
      else
        @metrics = metrics! 
      end
      return @metrics
    end  
    
    # recalculate metrics
    def metrics!
      metrics= Metrics.load_or_create(self)
      metrics.perform
      return metrics
    end  
    

    
    # number of promoters in this probeset
    def targets
      promoters.genelist.count
    end
    
 
    
    def roc_img_path
      File.join(dir,'matches',"#{name}_roc_auc.png").gsub('/Users/benjamin/','/images/')
    end  
  
  end #class Motif
end # module TFBS  