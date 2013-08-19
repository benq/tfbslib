# libraries needed for the tests
require 'test/unit'
require 'tfbs'
require 'tfbs/motif'



scer=TFBS::Genome.new('Saccharomyces Cerevisiae', 'scer', "#{TFBS::TEST_ROOT}/data/genomes/scer/scer.gff")
t1=TFBS::Promoters.new("tompa 1","#{TFBS::TEST_ROOT}/data/tompa/1/upstream.fa",scer)

f=TFBS::Bioprospector.new(t1)
f.perform
f.parse(:foo=>'bar')
