module TFBS
  class Motif
    #calculate Quality metrics on a motif
    class Metrics
      
      attr_accessor :motif
      
      # wrap the load if file available process 
      def self.load_or_create(m)
        fn=File.join(m.dir,'matches',"#{m.name}_metrics.yaml")
        if File.exist?(fn)
          metrics=YAML::load_file(fn)
          unless metrics==false
            metrics.motif=m
          else # just if a file is broken
            metrics    =TFBS::Motif::Metrics.new(m)
          end
        else
          metrics=TFBS::Motif::Metrics.new(m)
        end    
        return metrics
      end
        
      def initialize(m)
        @motif=m
        begin
          @church_on_rankedr.delete!
        rescue
          
        end    
      end  
      
      def update
        fn=File.join(@motif.dir,'matches',"#{@motif.name}_metrics.yaml")
        m=@motif
        @motif=nil
        File.open(fn, 'w') { |out| YAML.dump(self, out) } 
        @motif=m
      end  
      
      # calculate all metrics
      def perform
        ranked_genes
        ranked_by_class
        sortfield
        phit
        church
        roc_auc
        church_on_ranked
      end
      
      #
      def hits
        @hits ||= @motif.fragrep.size
      end

      def full_genome_hits
        @full_genome_hits ||= @motif.full_genome_fragrep.size
      end
      
      # percentage of promoters in probeset with a hit
      def phit
        @phit ||= phit!
      end
      
      def phit!
        puts "calc phits"
        numhits=0
        ranked_by_class.each { |e| numhits+=1 if e[2]>0 and e[0]==1    }
        @phit= numhits.to_f/@motif.targets
        update
        return @phit
      end  
        
      def sortfield
        @sortfield ||= sprintf('%0.3g %0.3g %0.3g',church.to_f,1-weighted_auc.to_f, church_on_ranked.to_f)
      end
        
      #return group specific score
      def church
        @church ||= church!
      end  
      
      # calculate group specific score Huges et al (alignace)
      # min(s1,s2)
      # sum dhypergeom i,hits,all_genes-hits, probe
      # i=intersect. 
      # in R  phypergeom(x,m,n,k)
      def church!
        puts "calc curch"    
        rk=scores_genes
        t=rk.size
        t>300 ? s1=300 : s1=rk.size
        s2=@motif.promoters.genelist.size      
        intersect=0
        rk[0..s1].each do |e| 
          intersect+=1 if e[2].to_f > 0.0 and e[0]==1    
        end  
        R.m=s1
        R.i=intersect
        R.n=t-s1
        R.k=s2
        R.eval("x=i:min(m,k)")
        R.eval("p=sum(dhyper(x,m,n,k))")
        p=R.p
        TFBS.log("church: i:#{intersect} m:#{s1} n:#{t-s1} k:#{s2} p:#{p} ")
        @church=p
        update
        return @church
      end
      
      # return saved rgss
      def church_on_ranked
        @church_on_ranked ||= church_on_ranked!
      end
      
        
      
      # calculate group specific score Huges et al (alignace)
      # min(s1,s2)
      # sum dhypergeom i,hits,all_genes-hits, probe
      # i=intersect. 
      # in R  phypergeom(x,m,n,k)
      def church_on_ranked!      
        puts "calc curchr"
        rk=ranked_by_class
        t=rk.size
        t>200 ? s1=200 : s1=rk.size
        s2=@motif.promoters.genelist.size      
        intersect=0
        rk[0..s1].each do |e| 
          intersect+=1 if e[2]>0 and e[0]==1    
        end  
        R.m=s1
        R.i=intersect
        R.n=t-s1
        R.k=s2
        R.eval("x=i:min(m,k)")
        R.eval("p=sum(dhyper(x,m,n,k))")
        p=R.p
        TFBS.log("church: i:#{intersect} m:#{s1} n:#{t-s1} k:#{s2} p:#{p} ")
        @church_on_ranked=p
        update
        return @church_on_ranked
      end
        
      # to sort we degrade any prediction with phit>2
      def weighted_auc
        phit>2 ? roc_auc/(phit-1) : roc_auc
      end  
      
      
      #instead of a full genome scan we check TP,FP,TN,FN
      # ref is format [seq_num, back from TSS, seq, length]
      #        eg.     8,-423,acacc,26
      # fragrep is format [seq_id,strand,forward count ]
      #        eg.     seq_8, +, 586, ACCCAATC
      # the example is a tp for each nucletide 1000-423>586, 586+8<1000-423+26
      def tompa_tp_fp
        ref_fn="/Users/benjamin/devel/tfbs/test/data/tompa/source/Tompa-dataset/transfac.txt"
        ref_fh=File.open(ref_fn,"r")
        ref_all={}
        section=""
        ds=""
        ref_fh.each do |line|
          if line =~ /^>(.*)/
            section='data set' if line =~ />data set/ 
            section='instances' if line =~ />instances/
          else
            ds=line.chomp if section=='data set'
            if section =='instances'
              ref_all.has_key?(ds) ? ref_all[ds]<<line.chomp.split(",") : ref_all[ds]=[line.chomp.split(",")]
            end
          end     
        end
        tompa_id=promoters.ename.split("_")[1]
        ref=ref_all["yst0#{tompa_id}"]
        TFBS.log("fragrep: #{fragrep}")
        res=[]
        fragrep.each do|hit|


          seq_num=hit[0].gsub("seq_","")
          ref.each do |r|
            if r[0]==seq_num
              #puts  "    #{1000+r[1].to_i} > #{hit[2]} #{hit[2].to_i+hit[3].size-1} > #{1000+r[1].to_i+26}   #{r[2]} ~ #{hit[3]}"
              tp=0;fp=0
              (hit[2].to_i..(hit[2].to_i+hit[3].size-1)).each do |n|  
                if (((1000+r[1].to_i)<=n) and (n<= (1000+r[1].to_i+26))) 
                   tp+=1 
                   #print "#{n} "
                else
                  fp+=1
                  #print "#{n} "
                end   
              end
              #puts " tp: #{tp} fp: #{fp}"
              res<<[tp,fp]
            end
          end

        end
        TFBS.log("tp,fp: #{res}")
        return res
      end  

      # this is a nucletide level tp/fp/tn/fn analysis
      def simple_tp_fp
        ref_fn="/Users/benjamin/devel/tfbs/test/data/tompa/source/Tompa-dataset/transfac.txt"
        ref_fh=File.open(ref_fn,"r")
        ref_all={}
        section=""
        ds=""
        ref_fh.each do |line|
          if line =~ /^>(.*)/
            section='data set' if line =~ />data set/ 
            section='instances' if line =~ />instances/
          else
            ds=line.chomp if section=='data set'
            if section =='instances'
              ref_all.has_key?(ds) ? ref_all[ds]<<line.chomp.split(",") : ref_all[ds]=[line.chomp.split(",")]
            end
          end     
        end
        tompa_id=promoters.ename.split("_")[1]
        ref=ref_all["yst0#{tompa_id}"]
        TFBS.log("fragrep: #{fragrep}")



        # mark true
        s={} 
        ref.each do |rs|
          seq_num="seq_#{rs[0]}"
          if s[seq_num].nil? 
            s[seq_num]="x".*1000 # x = TN, acgt = FN, ACGT = TP, X = FP
          end
          start=rs[1].to_i
          stop=rs[1].to_i+rs[3].to_i-1
          s[seq_num][start..stop]=rs[2].downcase 
        end

        # upcase all hits 
        fragrep.each do |hit|
          seq_num=hit[0]
          if s[seq_num].nil? 
             s[seq_num]="x".*1000 # x = TN, acgt = FN, ACGT = TP, X = FP
           end

          start=hit[2].to_i
          stop=hit[2].to_i+hit[3].size-1
          s[seq_num][start..stop]=s[seq_num][start..stop].upcase
        end

        res={}
        tn=0;fn=0;tp=0;fp=0
        # count
        s.each do |k,v|
          res[k]={}
          tn+=res[k][:tn]=s[k].scan(/x/).size
          fn+=res[k][:fn]=s[k].scan(/[a|c|g|t|n]/).size
          tp+=res[k][:tp]=s[k].scan(/[A|C|G|T|N]/).size
          fp+=res[k][:fp]=s[k].scan(/X/).size
          res[k][:fpr]=fp.to_f/(fp+tn)
          res[k][:tpr]=tp.to_f/(tp+fn)
        end    
        tpr=tp.to_f/(tp+fn)
        fpr=fp.to_f/(fp+tn)
        acc=(tp+tn).to_f/(fp+fn+tp+tn)
        spc= 1 - fpr
        ppv = tp.to_f / (tp + fp)
        fdr = fp.to_f/ (fp + tp)
        mmc= (tp*tn-fp*fn)/Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))  
        return {:tn=>tn,:fn=>fn,:tp=>tp,:fp=>fp,:tpr=>tpr,:fpr=>fpr,:acc=>acc,:spc=>spc,:ppv=>ppv,:fdr=>fdr,:mmc=>mmc,:list=>res}
      end


      def tompa_validation
        t=0;f=0
        tompa_tp_fp.each{|e| e[0]>0 ? t+=1 : f+=1}
        return t,f
      end  

      # compare positions of hits with positions of genes and keep trac of the number of motifs in each genes promoter
      # then rank genes by the most occurences 
      # calculate as r[gene_id]=hitcount
      # but return as [[gene_id,hitcount],...]
      def ranked_genes!
        #rs={}
        r={}
        # t [chr] = start, stop, strand, type, gene, count
        t={}
        @motif.taxon.genes.each do |k,g|
          t[g[0]] ||=[]
          t[g[0]] << [g[1],g[2],g[3],g[4],k,0] 
        end

        # r [chr] = pos
        f={}
        @motif.full_genome_fragrep.each do |e|
          f[e[0]] ||=[]
          f[e[0]] << e[2]
        end

        # 
        t.each_key do |chr|
          t[chr].each do |feature|
            if f.has_key?(chr)
              f[chr].each do |hit|
                if feature[2] == "+" # pos strand feature
                  if (feature[0]-1000 < hit.to_i) and (hit.to_i < feature[0].to_i)
                     q=feature[4]
                     r[q]=r[q].to_i+1
                  end    
                else # neg strand feature
                  if (feature[1].to_i < hit.to_i) and (hit.to_i < feature[1].to_i + 1000)
                    q=feature[4]
                    r[q]=r[q].to_i+1
                  end  
                end
              end
            else
              r[feature[4]]=0
            end    
          end
        end
        @ranked_genes=r.sort {|a,b| b[1] <=> a[1]}      
        update
        return @ranked_genes
      end

      def ranked_genes
        @ranked_genes ||= ranked_genes!
      end  

      # r[gene_id] = best score
      def scores_genes
        fn=File.join(@motif.dir,'matches',"#{@motif.name}_scores_genes.yaml")
        rs={}
        r={}
        cr=[]
        if File.exist?(fn)
          cr= YAML::load_file(fn)
        end
        # recover from broken file 
        if cr==nil or cr==false or cr.size==0
           cr=[] 
          # t [chr] = start, stop, strand, type, gene, count
          t={}
          @motif.taxon.genes.each {|k,g|
            t[g[0]] ||=[]
            t[g[0]] << [g[1],g[2],g[3],g[4],k,0] 
          }

          f={}
          @motif.full_genome_fragrep.each {|e|
            f[e[0]] ||=[]
            f[e[0]] << [e[2],e[4]]
          }
          gl=@motif.promoters.genelist.collect{|e| e[0] }
          #TFBS.log(gl.inspect)
          #TFBS.log(ranked_genes.inspect)

          t.each_key do |gene_id|
            t[gene_id].each do |feature|
              if f.has_key?(gene_id)
                f[gene_id].each do |hit,score|
                  if feature[2] == "+" # pos strand feature
                    if (feature[0]-1000 < hit.to_i) and (hit.to_i < feature[0].to_i)
                       r[feature[4]]=score.to_f if r[feature[4]].to_f < score.to_f
                    end    
                  else # neg strand feature
                    if (feature[1].to_i < hit.to_i) and (hit.to_i < feature[1].to_i + 1000)
                      r[feature[4]]=score.to_f if r[feature[4]].to_f < score.to_f
                    end  
                  end
                end
              else
                r[feature[4]]=0
              end    
            end
          end

          rr=r.sort {|a,b| b[1] <=> a[1]}
          rr.each do |k,v|
            if gl.include?(k)
              cr << [1,k,v]
            else
              if k !~ /orf/            
                cr << [0,k,v]
              end  
            end    
          end
          File.open(fn, 'w') { |out| YAML.dump(cr, out) }          
        end 
        return cr   
      end


      # return a genelist
      # [ probe? as 0/1, gene_id, hitcount]
      def ranked_by_class

        fn=File.join(@motif.dir,'matches',"#{@motif.name}_ranked_by_class.yaml")
        TFBS.log(fn)
        cr=[]
        if File.exist?(fn)
          cr= YAML::load_file(fn)
        end
        # handle empty files
        cr=[] if cr==nil  or cr==false or true
        if  cr.size==0
          gl=@motif.promoters.genelist.collect{|e| e[0] }
          #TFBS.log(gl.inspect)
          #TFBS.log(ranked_genes.inspect)
          ranked_genes.each do |k,v|
            if gl.include?(k)
              cr << [1,k,v]
            else
              cr << [0,k,v] if k !~ /orf/
            end    
          end
          File.open(fn, 'w') { |out| YAML.dump(cr, out) }   
        end
        TFBS.log "cr size: #{cr.size}"
        cr=cr.sort{|a,b| b[2] <=> a[2]} if cr.class==Array
        return cr
      end  

      # eighter a ranked AUC or based on annotated tp fp check
      def roc_auc
        @motif.name =~ /tompa/i ? tompa_auc : full_genome_auc
      end  


      # calculated on per nucleotide evaluation against the reference
      def tompa_auc
        labels=[]
        val=[]
        simple_tp_fp[:list].each{|k,v| labels<<0; val<<v[:tp]; labels<<1; val<<v[:fp] }
        auc=0
        fn=File.join(@motif.dir,'matches',"#{@motif.name}_roc_auc.png")
        fn2=File.join(@motif.dir,'matches',"#{@motif.name}_roc_auc.yaml")
        labels
        val
        # if  File.exists?(fn) and File.exists?(fn2)
        #   auc= YAML::load_file(fn2)
        # else
          R.eval('require(caTools)') 

          R.assign("labels",labels)
          R.assign("val",val)
          R.eval("png(\"#{fn}\",240,240)")

          R.eval("auc=colAUC(val, labels, plotROC=TRUE)")
          auc=R.pull( "auc[1]")
          R.eval('dev.off()')
          File.open(fn2, 'w') { |out| YAML.dump(auc, out) } 
        # end  
        return auc
      end  
      
      # return full genome auc 
      def full_genome_auc
        @full_genome_auc||=full_genome_auc!
      end

      # calculate and save ROC and the AUC
      def full_genome_auc!
        #if @full_genome_auc == nil
          puts "calc full genome auc"
          auc=0
          fn=File.join(@motif.dir,'matches',"#{@motif.name}_roc_auc.png")
          fn2=File.join(@motif.dir,'matches',"#{@motif.name}_roc_auc.yaml")
        # if  File.exists?(fn) and File.exists?(fn2) and false
        #  @full_genome_auc= YAML::load_file(fn2)
        # else  
            # # run R to create images
            # puts "R infile: #{infile}"
            # unless File.exists?("#{infile}.png")
            #   #r=RSRuby.instance
            rank_data={}
            rreg=[]
            rgenes=[]
            rcount=[]
            ranked_by_class.each{|row| rreg<<row[0].to_i; rgenes<<row[1].to_s; rcount<<row[2].to_i }

            R.eval('require(caTools)') 

            R.assign("rreg",rreg)
            R.assign("rcount",rcount)
            R.eval("png(\"#{fn}\",240,240)")

            R.eval("auc=colAUC(rcount, rreg, plotROC=TRUE)")
            auc=R.pull( "auc[1]")
            R.eval('dev.off()')
              #   puts R.eval("t=read.tree(\"#{infile}\")")
              #   puts R.eval('plot(t,"u",cex=0.3,lab4ut="axial")')
              #   puts R.eval('dev.off()')
            # else
            #   puts "xxx.tree.png found"  
            # end  
            # "#{infile}.png"
            @full_genome_auc=auc
            #File.open(fn2, 'w') { |out| YAML.dump(@full_genome_auc, out) } 
          #end
        #end
        update
        return @full_genome_auc  
      end  

      
        
        
    end #end class Metrics
  end # end class Motif
end # end module TFBS