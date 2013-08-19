# libraries needed for the tests
require 'test/unit'
require 'tfbs'
require 'tfbs/motif'


scer=TFBS::Genome.new('Saccharomyces Cerevisiae', 'scer', "#{TFBS::TEST_ROOT}/data/genome/scer.gff")
t1=TFBS::Promoters.new("tompa 1","#{TFBS::TEST_ROOT}/data/tompa/1/upstream.fa",scer)

TFBS::PREDICTORS.each {|p|
  
  r=TFBS::Run.new(t1,'motifsearch',p,[6,8,10,12])
  r.perform
}


