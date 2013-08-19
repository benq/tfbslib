
require 'ftools'
require 'tfbs/motif'

class TFBS::Motifsampler

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
        @promoters     = promoters
        @options      = options
        @motifsampler = "#{TFBS::METHODS_ROOT}/prediction/MotifSampler/motifsampler_ppc"
        @bg_bin       = "#{TFBS::METHODS_ROOT}/prediction/MotifSampler/CreateBackgroundModel.dms"
        @widths       = widths
        @mpr          = 10
      end 
      
      # promoter path
      attr_accessor :promoters

      # Program path
      attr_accessor :motifsampler



      # Options for weeder
      attr_accessor :options

      # runs as the adviser cannot handle a single run, it does not make sense to do this seperatly
      # format: [6,8,10,12]
      attr_accessor :widths
      
      # predictor name, use this method whenever refering to the prediction tool Method 
      def appl_name
        'motifsampler' 
      end

      def appl_citation
        'TODO: ... http://... '
      end
      
      
      def basedir
        @basedir ||= File.join(File.dirname(promoters.file), appl_name)
      end

      # create a background file next to the genome e.g. scer.motifsampler.bg and return path
      def create_background_from_genome
        bg=File.join(File.dirname(promoters.taxon.genome_path),"#{promoters.taxon.short}.#{appl_name}.bg")
        unless File.exists?(bg)
          `#{@bg_bin} -f "#{promoters.taxon.genome_path}" -b "#{bg}"`
        end
        return bg
      end  

      # return background file or create it
      def background
        @background ||= create_background_from_genome
      end  

      def background=(path)
        @background = path
      end


      # create a subdir weeder next to the input file and
      # run the prediction commands
      def prediction
        Dir.mkdir(basedir)
        logfile = File.join(basedir,'log')
        @widths.each do |width|
          Dir.mkdir("#{basedir}/#{width}")
          outfile = "#{basedir}/#{width}/out.txt"
          outfile2 = "#{basedir}/#{width}/out.mtrx"
          `#{@motifsampler} -f "#{promoters.clean_file}" -s 0 -b "#{background}" -w #{width} -o "#{outfile}" -m "#{outfile2}" -n #{@mpr}\n`
        end
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

      # parse and return an array of TFBS::Motifs
      # should be a report which contains motifs, but this will have to do for now  
    def parse(gtag)
      motifs=[]
      @widths.each do |width|
        in_fn=File.join(basedir,width.to_s, 'out.mtrx')
         section=:undef
         matrix=[]
         tag=Marshal.load(Marshal.dump(gtag))
         prog=appl_name
         rank=0
         File.open(in_fn,"r").each {|line|
           s = case line
             when /^#ID =/           then :new_motif
             when /^#Score/          then :occurrence_list
             when /^#Consensus/      then :matrix
             when /^$/               then :save_motif
             else 
               :do  
           end
           #p s, section
           if s!=:do
             section=s
             if (section==:new_motif and line =~ /^#ID = (\w+)$/)        then tag[:name]=$1; rank+=1 end
             if (section==:matrix    and line =~ /^#Consensus = (\w+)$/) then tag[:consensus]=$1 end
             if (section==:save_motif) then
               tag[:rank]=rank
               unless matrix.empty?
                 m=TFBS::Motif.new(matrix, tag)
                 m.promoters=promoters
                 m.predictor=self
                 m.rank=rank
                 motifs.push(m)
               end
               
               matrix=[]
               tag=Marshal.load(Marshal.dump(gtag))
               section = :new_motif
             end  
           else
             if (section==:matrix    and line =~ /(\d\.[\d|e|\-]+)\s+(\d\.[\d|e|\-]+)\s+(\d\.[\d|e|\-]+)\s+(\d\.[\d|e|\-]+)/ )
               #p $1,$2,$3,$4
               matrix.push({:A=>$1.to_f,:C=>$2.to_f,:G=>$3.to_f,:T=>$4.to_f})
             end   
           end    
        }
      end  
      return motifs 
    end   
end
 