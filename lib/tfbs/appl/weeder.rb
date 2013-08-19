module TFBS

  require 'ftools'
  require 'tfbs/motif'
  require 'bio'


  class Weeder
  
    
      
    # Create TFBS::Weeder factory 
    # 
    # 
    #   weeder_factory = TFBS::Weeder.new('PPAS',8, 2, options)
    #
    #   to run "weederTFBS.out -f inputfile -R 50 -O PPAS -W width -e number_of_allowed_mutations -N -N -T #{mpr}\n"
    # 
    #  Usage : weederTFBS.out -f inputfilename -O speciescode -W motif_width -e error -R repeatspercentage <options>
    # 
    # speciescode: two-letter species code (e.g. HS, MM, RN, SC, and so on)
    # -W : motif length
    # -e : number of mutations allowed
    # -R : sequences percentage on which motif has to appear
    # 
    # options (not required):
    # -S: process both strands of input sequences (default: single strand)
    # -M: multiple motif occurrences in each sequence (default: expect zero or one occurrences per sequence)
    # -T <number>: report top <number> motifs
    # -V: verbose mode
    # 
    def initialize( promoters, basedir = nil, options = nil, widths = [6,8,10,12] )
      @promoters   = promoters
      @options    = options
      @weeder     = "#{METHODS_ROOT}/prediction/Weeder1.4.2/bin/Darwin-i386/weederTFBS.out"
      @adviser    = "#{METHODS_ROOT}/prediction/Weeder1.4.2/bin/Darwin-i386/adviser.out"
      @widths     = widths
      @mpr        = 10
    end
    # promoter path
    attr_accessor :promoters   
  
    # Program path
    attr_accessor :weeder
    
    # Program path adviser weeder uses adviser to interpret raw results
    attr_accessor :adviser
    
    # species we want to run the prediction for   
    attr_accessor :taxon

    # Options for weeder
    attr_accessor :options
    
    # runs as the adviser cannot handle a single run, it does not make sense to do this seperatly
    # format: [6,8,10,12]
    attr_accessor :widths
    
    # predictor name, use this method whenever refering to the prediction tool Method 
    def appl_name
      'weeder' 
    end
    
    def appl_citation
      ['Pavesi 2004']
    end  
    
    def appl_manual
      'TODO: ... http://159.149.109.9/modtools/downloads/weedermanual.pdf '
    end  
    
    def basedir
      @basedir ||= File.join(File.dirname(promoters.file), appl_name)
    end
    
    # this should be configurable for now assume defaults
    def var(width)
      case 
        when width > 12 : 4
        when width > 8 : 3
        when width > 6 : 2
        else 1
      end      
    end  

    # weeder stores background in the FREQ directory which has to underneath the bin dir
    # create a frequency file as described in http://159.149.109.9/modtools/downloads/weedermanual.pdf section 4.1
    # to use with the weeder TFBS prediction software
    def create_background_from_genome
      mere6file=File.join(File.dirname(@weeder),'FreqFiles',"#{promoters.taxon.short}.6.freq")
      mere8file=File.join(File.dirname(@weeder),'FreqFiles',"#{promoters.taxon.short}.8.freq")
      unless File.exists?(mere6file) && File.exists?(mere8file)
        @mere6 = {}
        @mere8 = {}

        def addtofreqfiles(naseq)
          (0..naseq.to_s.length-6).each {|i|
            m6=naseq[i..i+5]
            m8=naseq[i..i+7]
            @mere6[m6].nil? ? @mere6[m6]=1 : @mere6[m6]+=1  
            if m8.length==8 then
              @mere8[m8].nil? ? @mere8[m8]=1 : @mere8[m8]+=1
            end
          } 
        end
    
        gff=Bio::GFF::GFF3.new(File.open(promoters.taxon.genome_path))
        gff.sequences.each do |e|
           addtofreqfiles(e.to_s)
        end

        f6=File.open(mere6file,"w")
        @mere6.each {|k,v|
          f6.puts "#{k}\t#{v}"
        }
        f6.close

        f8=File.open(mere8file,"w")
        @mere8.each {|k,v|
          f8.puts "#{k}\t#{v}"
        }
        f8.close
      end
      true  
    end  
    
    def background?
      @background ||= create_background_from_genome
    end
    
    # create a subdir weeder next to the inpu file and
    # run the prediction commands
    def prediction
      background?
      Dir.mkdir(basedir)
      weederin = File.join(basedir,'out')
      logfile = File.join(basedir,'log')
      File.copy(promoters.clean_file,weederin) # as weeder names all files liek the infile 
      @widths.each do |width|
         `cd "#{File.dirname(@weeder)}"; #{@weeder} -f "#{weederin}" -R 50 -O #{promoters.taxon.short} -W #{width} -e #{var(width)} -N -N -T #{@mpr}\n >> "#{logfile}"`
         #puts "rm #{p}/#{i}/#{m}/#{w}/out"
      end
      `#{@adviser} "#{weederin}" >>#{logfile}` 
    end
    
    # create a subdir <appl_name> next to the input file and
    # run the prediction commands
    def perform
      unless File.exists?(basedir)
        prediction
      else
        "dir exists, prediction skipped for #{appl_name} in #{promoters.file}"   
      end  
    end
    
    # parse and return an arry of TFBS::Motifs
    # should be a report which contains motifs, but this will have to do for now  
    def parse(gtag)
      in_fn=File.join(basedir, 'out.wee')
      motifs=[]
        section=:undef
        matrix=[]
        tag=Marshal.load(Marshal.dump(gtag))
        prog=appl_name
        ranklist={}
        File.open(in_fn,"r").each {|line|

          s = case line
          when /^\*\*\*\* MY ADVICE \*\*\*\*/ then :new_motif
          when /^Best occurences/             then :occurrence_list
          when /Frequency Matrix/             then :matrix
          when /==========================================/ then :save_motif
          else :do  
          end
          #p s, section
          if s!=:do
            section=s  
          else
            if line =~ /^(\d+)\)\s(\w+)\s(\d\.\d\d)\s+$/ then ranklist[$2]=[$1,$3] end # ranklist[consensus]=[rank, precedence] 
            if (section==:new_motif and line =~ /^(\w+)$/) then tag[:consensus]=$1 end
            if (section==:new_motif and line =~ /^(\d) redundant motifs found/) then tag[:redundant_motifs]=$1 end
            if (section==:matrix    and line =~ /(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ )
               matrix[$1.to_i-1]={:A=>$2,:C=>$3,:G=>$4,:T=>$5}
            end   

            unless (section != :save_motif) then
              tag[:rank],tag[:precedence]=ranklist[tag[:consensus]]
              m=TFBS::Motif.new(matrix, tag)
              m.promoters=promoters
              m.predictor=self
              m.rank=tag[:rank]
              motifs.push(m)
              tag=Marshal.load(Marshal.dump(gtag))
              matrix=[]
              section = :new_motif
            else   
            end
          end      
        }
      return motifs
    end
  end # class Weeder
end # module TFBS  