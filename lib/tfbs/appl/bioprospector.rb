
require 'ftools'
require 'tfbs/motif'
require 'bio'

class TFBS::Bioprospector

      # Create TFBS::Bioprospector factory 
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
        @promoters      = promoters
        @options       = options
        @bioprospector = "#{TFBS::METHODS_ROOT}/prediction/BioProspector/BioProspector.mac"
        @widths        = widths
        @mpr           = 10
        @bg_bin        = "#{TFBS::METHODS_ROOT}/prediction/BioProspector/genomebg.mac"
      end 

      # promoter path
      attr_accessor :promoters
      
      # Program path
      attr_accessor :bioprospector


      # Options for weeder
      attr_accessor :options

      # runs as the adviser cannot handle a single run, it does not make sense to do this seperatly
      # format: [[6,1],[8,2],[10,3],[12,3]]
      attr_accessor :widths
      
      # predictor name, use this method whenever refering to the prediction tool Method 
      def appl_name
        'bioprospector' 
      end

      def appl_citation
        'TODO: ... http://... '
      end
      
      def basedir
        @basedir ||= File.join(File.dirname(promoters.file), appl_name)
      end

      # create a background file next to the genome e.g. scer.motifsampler.bg and return path
      # the background utility only accepts up to 200kb filesize this puts some bias on telomere regions :-(
      def create_background_from_genome
        bg=File.join(File.dirname(promoters.taxon.genome_path),"#{promoters.taxon.short}.#{appl_name}.bg")
        genome_fasta=File.join(File.dirname(promoters.taxon.genome_path),"#{promoters.taxon.short}.#{appl_name}.fasta")
        unless File.exists?(bg)
          gff=Bio::GFF::GFF3.new(File.open(promoters.taxon.genome_path))
          ff=File.new(genome_fasta,"w")
          gff.sequences.each{|s| 
            ff.puts ">#{s.seq.id}"
            ff.puts s.seq.slice(1,1000) }
          ff.close
          sleep 1
          `#{@bg_bin} -i "#{genome_fasta}" -o "#{bg}"`
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
      
      def prediction
        Dir.mkdir(basedir)              
        logfile = File.join(@basedir,'log')
        @widths.each do |width|
          Dir.mkdir("#{basedir}/#{width}")
          outfile = "#{basedir}/#{width}/out.txt"
          `#{@bioprospector} -i "#{promoters.clean_file}" -f "#{background}" -W  #{width} -o "#{outfile}" -r #{@mpr}\n`
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
      in_fn=File.join(basedir, width.to_s,'out.txt')
       section=:undef
       matrix=[]
       tag=Marshal.load(Marshal.dump(gtag))
       prog=appl_name
       File.open(in_fn,"r").each {|line|

         s = case line
           when /^Motif/               then :new_motif
           when /^Width/               then :detail
           when /^Blk/                 then :matrix
           when /^>/                   then :occurrence_list  
           when /\**\*/                then :save_motif
           else 
             :do  
         end
         #p s, section
         if s!=:do
           section=s
           if (section==:new_motif and line =~ /^Motif #(\d+): \((\w+)\/(\w+)\)$/) then 
             tag[:rank]=$1; tag[:consensus]=$2; tag[:reverse_complement]=$3 
           end  
         else

           if (section==:matrix    and line =~ /(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\w)\s+(\w)\s+(\w)\s+(\w)/ )
              matrix[$1.to_i-1]={:A=>$2,:C=>$3,:G=>$4,:T=>$5}
           end   

           unless (section != :save_motif) then
             unless matrix.empty?
               m=TFBS::Motif.new(matrix, tag)
               m.promoters=promoters
               m.predictor=self
               m.rank=tag[:rank]
               motifs.push(m)
               matrix=[]
               tag=Marshal.load(Marshal.dump(gtag))
               section = :new_motif
             end  
           else   
           end
         end    
       }
      end 
    return motifs 
  end   
end
