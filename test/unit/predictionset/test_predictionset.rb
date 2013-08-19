# libraries needed for the tests
require 'test/unit'
require 'tfbs'
require 'tfbs/motif'


scer=TFBS::Genome.new('Saccharomyces Cerevisiae', 'scer', "#{TFBS::TEST_ROOT}/data/genome/scer.gff")



#[1,2,3,4,5,6,8,
  [2].each {|p|
  t=TFBS::Promoters.new("tompa #{p}","#{TFBS::TEST_ROOT}/data/tompa/#{p}/upstream.fa",scer)
  f=TFBS::PredictionSet.new(t,'motifsearch')
  f.perform

  l=f.motifs
  l.origin=f.promoters
  l.name="pred_for_#{f.promoters.ename}"
  #l.dir="#{TFBS_ROOT}/test/data/tompa/#{p}/default_list"
  l.each{|m| puts 
    puts "failed: render_to_text" unless m.render_to_text(:list)
    puts "failed: weblogo" unless m.weblogo
    puts "failed: fragrep" unless m.fragrep
  }
  l.clusters.each{|c|
    puts c.name
    puts "failed: render_to_text" unless c.render_to_text(:list)
    puts "failed: weblogo" unless c.weblogo
    puts "failed: fragrep" unless c.fragrep
    # p c
    c.motifs.dir=c.dir
    c.motifs.tree
    c.motifs.each{|m|
        puts "    #{m.name}"
      }
  }  
  puts l.dir
  l.tree
  l.clusters.tree
  
  
  p l.clusters.names
  # l.tree_to_file('/Volumes/Data\ HD/Users/almeida/devel/tfbs/test/data/tompa/1/clusters/motifs.tree')
   # l.tree_to_file('/Volumes/Data\ HD/Users/almeida/devel/tfbs/test/data/tompa/1/clusters/motifs.clusters.tree')

  #p l.clusters
}
