
module TFBS::MotiflistHelper
  def names
    n=[]
    self.each{|e| n << e.name}
    n
  end
  
  
    
    # eigther promoter or a motiflist
    attr_accessor :origin
    

    
    def name
      @name ||= @origin.name ||= default_list
    end  
    
    # source of this data
    attr_accessor :source
    
    # more description 
    attr_accessor :description
    

    
    def origin
      @origin ||=@promoter
    end
    

    
    # make life simple?
    def cluster_dir
      File.join(dir,'clusters')
    end 
    
    # once a list is created, we should have an index to uniquly address each
    def create_index
      
    end  
    
    def find_by_cluster
      
    end  
    
    def find_by_n(n)
      handle=nil
      self.each{|m| handle=m  if n.to_i==m.tag[:n].to_i}
      handle
    end
    
    # get a motiflist
    def clusters
      @clusters ||= create_clusters
    end  
    
    #run STAMP
    def create_clusters
       f=TFBS::Stamp.new
       cl= f.cluster(render_to_file(:transfac,cluster_dir), {:cluster => true})
       sublist=parse_transfac(cl, :foo=>:bar)
       sublist.origin=self
       sublist.name= "cl_" + name
       sublist.each {|c| c.dir=File.join(dir,"cl_" + name)}
       sublist
    end  
    

    

    
    def to_jasper
      self.collect{ |motif| motif.output(:jasper) }.join("\n")
    end  
    
    # may be a hack with large lists but currently more than 200 is unusual
    def jasper_file
      fn="#{dir}/#{name}.motifs.jasper"
      fh=File.new(fn,"w")
      fh.puts to_jasper
      fh.close
      fn
    end  
    
    def to_transfac
      self.collect{ |motif| motif.render_to_text(:transfac) }.join("\n")
    end  
    
    # may be a hack with large lists but currently more than 200 is unusual
    def transfac_file
      fn="#{dir}/#{name}.motifs.transfac"
      fh=File.new(fn,"w")
      fh.puts to_transfac
      fh.close
      fn
    end
    
    def fragrep_files
      fragrep_dir="#{dir}/motifs.#{:fragrep}"
      self.each do |motif|
        motif.output_fragrep_pattern(fragrep_dir)
        motif.output_weblogo("#{dir}/weblogo")
      end
      fragrep_dir
    end  
    
    # creates file and returns the filename
    # e.g. render_to_file(:jasper, 'some_motifs')
    def render_to_file(format,outdir)
      Dir.mkdir(outdir) unless File.exists?(outdir)
      p fn="#{outdir}/motifs.#{format}"
      data =case format
      when :jasper
        "#jasper format\n" +
        self.collect{ |motif| motif.render_to_text(format) }.join("\n")
      when :transfac
        #data="#transfac format\n"
        self.collect{ |motif| motif.render_to_text(format) }.join("\n")
      when :tree_tf
        #data="#transfac format\n"
        self.collect{ |motif| motif.render_to_text(format) }.join("\n")  
      when :fragrep
        self.each do |motif|
          motif.output_fragrep_pattern("#{outdir}/#{format}")
          motif.output_weblogo("#{outdir}/weblogo")
        end
      end
      if data.kind_of?(String) then
        fh=File.new(fn,"w+")
        fh.puts data
        fh.close
      end
      fn unless format==:fragrep
    end
  
  
  #run R to create a tree
  # what:
  # * plain
  # * clusters (default)
  # * both (crazy)
  # compare_lists
  # TFBS::HARBINGER
  # TFBS::MCISAACK
  def tree(options={})
    unless self.size < 2
      f=TFBS::Stamp.new
      tf=render_to_file(:tree_tf,dir)
        if options.has_key?(:match)
          f.cluster(tf,{:match=> '/Users/benjamin/devel/tfbs/test/data/harbison/motifs.transfac'})
        else  
          tree= f.cluster(tf)
          draw_tree(tree)
        end  
      else
      "not_enough_input"
    end   
  end
  
  # append png to name and draw tree
  def draw_tree(infile)
    # run R to create images
    puts "R infile: #{infile}"
    unless File.exists?("#{infile}.png")
      #r=RSRuby.instance
      R.eval('require(ape)') 
      # render tree 
      puts R.eval("png(\"#{infile}.png\",1280,1024)")
      puts R.eval("t=read.tree(\"#{infile}\")")
      puts R.eval('plot(t,"u",cex=0.3,lab4ut="axial")')
      puts R.eval('dev.off()')
    else
      puts "xxx.tree.png found"  
    end  
    "#{infile}.png"
  end
  
  
  # TODO test if this works
  # put this in its own class or module
  # children:
  # XX	Cluster_Members:	tompa_2_CCGGCGSWCT_meme	tompa_2_ACCGGCGS_meme	tompa_2_CCGGCG_meme	
  # has to be adapted to the transfac rules mazbe use the bio::transfac...
  # currently it works for stamp output
  def parse_transfac(in_fn,gtag)
    puts 'parsing the clusters transfac'
    section=:undef
    matrix=[]
    motifs=TFBS::Motiflist.new()
    tag=Marshal.load(Marshal.dump(gtag))
    name=""
    prog='transfac'
    mtrx=[]
    cluster_members=[]
    File.open(in_fn,"r").each {|line|
      case line
        when /^DE\s(.*)/        
          section = :new_motif
          name=$1
        when /^(\d+)\s(\d\.\d*)\s(\d\.\d*)\s(\d\.\d*)\s(\d\.\d*)\s(\w)$/  
          section = :matrix
                        matrix.push({:A=>$2,:C=>$3,:G=>$4,:T=>$5})
        when /^XX\s.*Cluster_Members:(.*)/
          section = :cluster_members  
          cluster_members=$1.split
        when /^$/         
          section = :save_motif
          if matrix.length > 0
            m=TFBS::Motif.new(matrix, tag)
            # this may not work for a mixed promoter
            m.promoters=self.first.promoters
            m.predictor=TFBS::Stamp
            m.motifs=self.find_by_name(cluster_members)
            m.extend(TFBS::MotiflistHelper)
            m.name=name
            m.origin=self
            motifs.push(m)
          end
          p matrix
          # reset
          matrix=[]
          mtrx=[]
          tag=Marshal.load(Marshal.dump(gtag))  
          section = :new_motif
        else
          section = :do
          # nice file format nothing to be done albeit the above
      end
    }
    motifs
  end
  
  # parse the TAMO format
  # to use macIsaak & Harbinger
  #
  def parse_tamo(in_fn,gtag)
    gtag.has_key?(:prefix) ? prefix=gtag[:prefix] : prefix='motif'
    motifs=TFBS::Motiflist.new()
    section=:undef
    matrix=[]
    tag=Marshal.load(Marshal.dump(gtag))
    name="";consensus='';source=[]
    prog='tamo'
    mtrx=[]
    File.open(in_fn,"r").each {|line|
      s = case line
        when /^Log-odds matrix for Motif\s+(\d+)\s(.+)\s\(1\)$/                  
          section=:new_motif
          p name=[prefix,$1].compact.join('_')
          consensus=$2
        when /^Probability matrix for Motif\s+(\d+)\s(.+)\s\(1\)$/ 
          section=:matrix
          # see do in section below
        when /^Sequence Logo$/
          section=:sequence_logo
          # skip
        when /^Source:\s(.*)$/
          section=:source
          source=$1.split  
        when /^MAP Score:/  
          section=:save_motif
          # flip matrix
          #p mtrx
          mtrx.transpose.each{|pos|
            matrix.push({:A=>pos[0],:C=>pos[1],:G=>pos[2],:T=>pos[3]})
          }
          p matrix
          if matrix.length > 0
            m=TFBS::Motif.new(matrix, tag) 
            m.name=[source[0],name].compact.join('_')
            m.motiflist=motifs
            motifs.push(m)
          else
            puts "empty matrix!"  
          end
          # reset
          matrix=[]
          mtrx=[]
          tag=Marshal.load(Marshal.dump(gtag))
          section = :new_motif
        else
          # do: for lines within a section
          if (section==:matrix    and line =~ /^#([ACGT])\s(.*)/ )
             #p $1
             #p $2
             mtrx.push($2.to_s.split)
             p mtrx
          end
      end
      #p s, section
      # if s!=:do
      #   section=s
      #   if line =~ /^Log-odds matrix for Motif\s(\d+)\s(.+)\s\(1\)$/ then tag[:name]="Motif #{$1}"; tag[:consensus]=$2  end
      # else
      #   if (section==:matrix    and line =~ /^#([ACGT])\s(.*)/ )
      #      #p $1
      #      #p $2
      #      mtrx.push($2.to_s.split)
      #   end
      # 
      #   unless (section != :save_motif) then
      #     # flip matrix
      #     #p mtrx
      #     mtrx.transpose.each{|pos|
      #       matrix.push({:A=>pos[0],:C=>pos[1],:G=>pos[2],:T=>pos[3]})
      #     }
      #     #p matrix
      #     m=TFBS::Motif.new(matrix, tag) if matrix.length > 0
      #     m.name=
      #     motifs.push(m)
      #      
      #      matrix=[]
      #      mtrx=[]
      #      tag=Marshal.load(Marshal.dump(gtag))
      #      section = :new_motif
      #   else
      #   end
      # end

    }
    motifs
  end
  
  
  # motif have names and these names are used in the transfac file for clustering
  # we want to search for taht name and return a subselection of self
  # TODO accept string or array
  def find_by_name(query)
    p query
    subset=TFBS::Motiflist.new()
    query.each{|q|
      self.each{|m|
        if "#{m.name}"=="#{q}"
          subset.push(m)
        end  
      }  
    }
    #p subset.names
    subset if subset.size > 0
  end
  
end