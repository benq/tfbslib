require 'ftools'
require 'tfbs/motif'

# Wrap Meme
class TFBS::Meme
  
  CITATIONS=['Timothy L. Bailey and Charles Elkan, "Fitting a mixture model by expectation maximization to discover motifs in biopolymers", Proceedings of the Second International Conference on Intelligent Systems for Molecular Biology, pp. 28-36, AAAI Press, Menlo Park, California, 1994.']
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
        @meme       = "#{TFBS::METHODS_ROOT}/prediction/meme4/bin/meme"
        @widths     = widths
        @mpr        = 10
        @background = ""

      end 

      # promoter path
      attr_accessor :promoters

      # Program path
      attr_accessor :motifsampler

      # species we want to run the prediction for   
      attr_accessor :taxon

      # Options for weeder
      attr_accessor :options

      # runs as the adviser cannot handle a single run, it does not make sense to do this seperatly
      # format: [[6,1],[8,2],[10,3],[12,3]]
      attr_accessor :widths
      
      # predictor name, use this method whenever refering to the prediction tool Method 
      def appl_name
        'meme' 
      end

      def appl_citation
        'TODO: ... http://... '
      end
      
      def basedir
        @basedir ||= File.join(File.dirname(promoters.file), appl_name)
      end  
      
      # TODO add background (meme creates one from the promoter when none is supplied)
      # weeder stores background in the FREQ directory which has to underneath the bin dir
      def create_background_from_genome
        # ...
      end  

      # return background file or create it
      def background
         @background ||= create_background_from_genome
      end  

      # create a subdir weeder next to the inpu file and
      # run the prediction commands
      def prediction
        Dir.mkdir(basedir)
        logfile = File.join(basedir,'log')
        @widths.each do |width|
          # careful meme cannot handle spaces in file or directory names
           `#{@meme} "#{promoters.file.gsub("/Volumes/Data HD",'')}" -oc "#{basedir.gsub("/Volumes/Data HD",'')}/#{width}" -maxsize 1000000 -nmotifs #{@mpr} -w #{width} -dna >>"#{logfile}" `
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

      # parse and return an arry of TFBS::Motifs
      # should be a report which contains motifs, but this will have to do for now  
      def parse(gtag)
        motifs=[]
        @widths.each do |width|
          in_fn=File.join(basedir, width.to_s,"meme.txt")
          section=:undef
          matrix=[]
          tag=Marshal.load(Marshal.dump(gtag))
          name=""
          prog=appl_name
          File.open(in_fn,"r").each {|line|
            s = case line
              when /^MOTIF/                                           then :new_motif
              when /Motif (\d+) Description/                          then :description
              when /Motif (\d+) position-specific probability matrix/ then :matrix
              when /^Time/                                            then :save_motif
              else 
                :do  
            end
            #p s, section
            if s!=:do
              section=s
              if line =~ /^MOTIF\s+(\d+)/ then tag[:name]="Motif #{$1}"; tag[:rank]=$1  end  
            else
              if (section==:matrix    and line =~ /(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)/ )
                 matrix.push({:A=>$1,:C=>$2,:G=>$3,:T=>$4})
              end   

              unless (section != :save_motif) then
                m=TFBS::Motif.new(matrix, tag)
                m.promoters=promoters
                m.predictor=self
                m.rank=tag[:rank]
                motifs.push(m)
                 matrix=[]
                 tag=Marshal.load(Marshal.dump(gtag))
                 section = :new_motif
              else   
              end
            end    
          }
        end  
        return motifs
      end #parse
end
