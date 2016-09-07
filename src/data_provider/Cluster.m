classdef Cluster
    
    properties
        features
        centroid
        num_points
        variance
    end
    
    methods
        
        function [this] = Cluster(f,data)
            this.features = f;
            x = vertcat(data.flw_out_vph);
            this.centroid = Utils.meanwithnan(x,1);
            this.num_points = numel(data);
            this.variance = sqrt(Utils.varwithnan(x,1)); 
        end
        
    end
    
end

