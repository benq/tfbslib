module TFBS
      require 'active_support/inflector'
  # to ease use from Rails
  module HasDataDir

    def self.included(base)
      base.extend ClassMethods
      base.class_eval do
        # check if parent Dir is there # only reasnoable for instances
      end  
    end  
     
    # def save
    #   File.new(File.join(dir,'cache.yaml'))self.to_yaml
    # end  
    # 
    # # TODO this is still unsafe!
    # def self.find(query)
    #   if query == :all
    #     return Dir.glob("#{self.dir}/**").collect {|d| File.basename(d) }
    #   else 
    #     return Dir.glob("#{self.dir}/#{query}").collect {|d| File.basename(d) }
    #   end    
    # end  
    # 
    # def self.load(name)
    #   return YAML::load(File.open(File.join(TFBS::ProbeSet.dir,"cache.yaml"),'r').read)
    # end  
    
    module ClassMethods
      def dir
        File.join(TFBS.data_root,self.to_s.gsub('TFBS::','').downcase.pluralize)
        #File.Join(TFBS.data_root, self.class.gsub('TFBS::','').pluralize)
      end  
    end
    
    # each probeset gets a directory in the TFBS_DATA/probesets folder
    public 
    def dir
      if self.methods.include?('parent')
        d=File.join(self.parent.dir,self.class.to_s.gsub('TFBS::','').downcase.pluralize, self.name)
      else
        d=File.join(self.class.dir,self.name)
      end
      Dir.mkdir(d) unless File.exist?(d)
      return d    
    end
    
  end
end    