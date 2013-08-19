#require 'rsruby'
module TFBS::Motiflist::TreeHelper
  
  def tree_of(what)
    case what
    when :motifs
      tree_to_file("#{self.name}.motifs.tree")
    when :clusters
      tree_to_file("#{self.name}.motifs.clusters.tree")
    end  
  end
  
  # append png to name and draw tree
  def tree_to_file(infile)
    "#{self.name}"
    
    # run R to create images
    
    r=RSRuby.instance
    r.eval_R('require(ape)') 

    # render tree 
    r.png("#{infile}.png")
    r.eval_R("t=read.tree(\"#{infile}\")")
    r.eval_R('plot(t,"u",cex=0.3,lab4ut="axial")')
    r.eval_R('dev.off()')
    r.pdf("#{infile}.pdf")
    r.eval_R('plot(t,"u",cex=0.3,lab4ut="axial")')
    r.eval_R('dev.off()')

  end  
  
end  



