class TFBS::REPORT
  
  
  # Report Builder...
  # eigentlich wollt ich das so machen aber wurscht
  # r = TFBS::REPORT.new(title)
  # r.desc "..."
  # r.section {|s| s.img prediction.tree
  # prediction.clusters.each{|c|
  #   r.section {|s|
  #    s.title "custer1 ..."
  #    s.img c.tree 
  #    s.sidetext "custer1 ..."
  #    s.motif c.motif
  #    c.each {|m| s.section {|sub| sub.motif m.motif} }  
  #   }
  # }
  def initialize(title)
    @title = title
    @sections=[]
  end
  
  def add_section=(subtitle, block)
    @sections << TFBS.Report.new(subtitle, block)
  end
  
  def render
    puts title
    @sections.each {|s| s.render }
    puts footer
  end  
  
end  