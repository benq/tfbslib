# libraries needed for the tests
require 'test/unit'
require 'tfbs'
require 'tfbs/motif'


scer=TFBS::Taxon.new('Saccharomyces Cerevisiae', 'scer', "#{TFBS_ROOT}/test/data/genome/scer.gff")
any=TFBS::Motiflist.new()


harb_dir="/Volumes/Data\ HD/Users/almeida/devel/tfbs/test/data/harbison/"

harb=any.parse_tamo("#{harb_dir}/Final_InTableS2_v24.motifs",:prefix=>'harbison')
puts 'parsed harbison'
harb.dir=harb_dir
 harb.each{|e| 
   puts "failed: weblogo" unless e.weblogo
 }

 harb.render_to_file(:transfac,harb_dir)
 puts harb.names
