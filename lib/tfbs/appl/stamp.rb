# require 'rubygems'
# require 'rsruby'

class TFBS::Stamp
  
  
  # stamp factory
  # -tf [input file]
  #     Input dataset of motifs in TRANSFAC format [required!]
  # -sd [score file]
  #     Input file with random score distributions [required!]
  # -chp
	#	    Calculate Calinski & Harabasz statistic and output the resulting clusters to a file.
  # -match [TF matrix file]
  #     Match the input PSSMs against this dataset.
  #
  def initialize
    @stamp= "/srv/tfbsdig/methods/clustering/STAMP.v1.1/STAMP"
    @scoretable = '/srv/tfbsdig/methods/clustering/STAMP.v1.1/ScoreDists/JaspRand_PCC_SWU.scores'
  end
  
  def self.appl_name
    'stamp'
  end
  
  # cluster our motifs 
  # assume we use a dir ../clusters relative to the motifs dir (for now)
  # this also works for a second run, but not a third unless you rename the output
  def cluster(motifs_tf, options={})
    #dir           = File.join(File.dirname(motifs_tf),'..','clusters')
    dir           = File.dirname(motifs_tf)
    Dir.mkdir(dir) unless File.exists?(dir)
    name          = File.basename(motifs_tf,'.*')
    unless File.exists?("#{dir}/#{name}.tree")  
      opt = {}
      opt['-sd']    = "#{@scoretable}"
      opt['-tf']    = motifs_tf
      opt['-match']   = options[:match] if options.has_key?(:match)
      opt['-chp'] = nil if options.has_key?(:cluster)
      o             = opt.collect{|k,v| "#{k} \"#{v}\""}.join(' ')
      puts `cd "#{dir}"; #{@stamp} #{o}`
      puts "cd \"#{dir}\"; #{@stamp} #{o}"
      File.rename("#{dir}/out.tree", "#{dir}/#{name}.tree")
      File.rename("#{dir}/out_tree_clusters.txt", "#{dir}/#{name}.clusters.transfac") if options.has_key?(:cluster)
    end
    if options.has_key?(:cluster)
      "#{dir}/#{name}.clusters.transfac"
    else
      "#{dir}/#{name}.tree"
    end      
  end  



  
  
  
end  