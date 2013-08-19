class TFBS::Motif::Tamo
  
  # parse the TAMO format
  # to use macIsaak & Harbinger
  #
  def self.parse(in_fn,gtag)
    motifs=[]
    section=:undef
    matrix=[]
    tag=Marshal.load(Marshal.dump(gtag))
    name=""
    prog='tamo'
    mtrx=[]
    File.open(in_fn,"r").each {|line|
      s = case line
        when /^Log-odds matrix for Motif\s\d/                   then :new_motif
        when /Probability matrix for Motif   (\d+) (\w+) \(1\)/ then :matrix
        when /^Sequence Logo$/                                  then :save_motif
        else
          :do
      end
      #p s, section
      if s!=:do
        section=s
        if line =~ /^Log-odds matrix for Motif\s(\d+)\s(.+)\s\(1\)$/ then tag[:name]="Motif #{$1}"; tag[:consensus]=$2  end
      else
        if (section==:matrix    and line =~ /^#([ACGT])\s(.*)/ )
           #p $1
           #p $2
           mtrx.push($2.to_s.split)
        end

        unless (section != :save_motif) then
          # flip matrix
          #p mtrx
          mtrx.transpose.each{|pos|
            matrix.push({:A=>pos[0],:C=>pos[1],:G=>pos[2],:T=>pos[3]})
          }
          #p matrix
          motifs.push(Motif.new(matrix, tag)) if matrix.length > 0
           
           matrix=[]
           mtrx=[]
           tag=Marshal.load(Marshal.dump(gtag))
           section = :new_motif
        else
        end
      end

    }
    motifs
  end
  
end  