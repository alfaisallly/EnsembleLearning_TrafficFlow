classdef ClusterManager
    
    properties
        time
        clusters
    end
    
    methods( Access = public )
        
        function [this] = ClusterManager(ids,days,time,data,method)
            this.time = time;
            switch method
                case {'','day_of_week'}
                    this.clusters = ClusterManager.cluster_by_day_of_week(ids,days,data);
                case 'k-means'
                    this.clusters = ClusterManager.cluster_by_k_means(ids,days,data);
                otherwise
                    error('unknown method')
            end
                   
        end
        
        function [x] = get_keys(this)
            keys = this.clusters.keys;
            for i=1:length(keys)
                x(i) = this.dehash_features(keys{i});
            end
        end
        
        function [this] = add_cluster(this,c)
            this.clusters(ClusterManager.hash_features(c.features)) = c;
        end
        
        function [c] = get_cluster(this,f)
            c = this.clusters(ClusterManager.hash_features(f));
        end
        
    end
    
    methods( Static, Access = private )
       
        function [h] = hash_features(f)
            %h = DataHash(orderfields(f));
            h = native2unicode(getByteStreamFromArray(orderfields(f)));
        end
        
        function [h] = dehash_features(f)
            %h = DataHash(orderfields(f));
            h = getArrayFromByteStream(unicode2native(f));
        end
        
        function [cluster_map] = cluster_by_day_of_week(ids,days,data)
            if any([numel(days) numel(ids)] ~= size(data))
                error('any([numel(days) numel(ids)] ~= size(data))')
            end
            day_nums = weekday(datestr(days));
            cluster_map = containers.Map;
            for i = 1:7
                day_ind = day_nums==i;
                if ~any(day_ind)
                    continue
                end
                for j = 1:numel(ids)
                    h = ClusterManager.hash_features(struct('day_of_week',i,'id',ids(j)));
                    cluster_map(h) = Cluster(struct('day_of_week',i,'id',ids(j)),data(day_ind,j));
                end
            end
        end
        
        function [cluster_map] = cluster_by_k_means(ids,days,data)
            
            if any([numel(days) numel(ids)] ~= size(data))
                error('any([numel(days) numel(ids)] ~= size(data))')
            end
            
%             
%             % HONGYUAN'S CLUSTERING METHOD
% %             day_nums = weekday(datestr(days));
            cluster_map = containers.Map;
%             
%             flow = vertcat(data.flw_out_vph);
%             
%             ValidCluster = 0 ;
%             trial = 3;
%             while ValidCluster == 0
%                 
%                 if trial == 0
%                     K=K-1;
%                     trial = 3;
%                     %warningmessage = sprintf('warning: no valid clustering found, automatic decrease K to %d',K);
%                     %display(warningmessage);
%                 end
%                 
%                 [clusterID,clusterCentre]=kmeans(flow,K,'EmptyAction','singleton');
%                 
%                 CheckCluster = zeros(K,1);
%                 for k = 1:K
%                     r = find(clusterID == k);
%                     if length(r) > 2
%                         CheckCluster(k)=1;
%                     end
%                 end
%                 
%                 if sum(CheckCluster == K)
%                     ValidCluster = 1 ;
%                 end
%                 
%                 trial = trial - 1;
%             end
%             
%             
% 
% 
%             % COPY OF DOW
%             for i = 1:7
%                 day_ind = day_nums==i;
%                 if ~any(day_ind)
%                     continue
%                 end
%                 for j = 1:numel(ids)
%                     h = ClusterManager.hash_features(struct('day_of_week',i,'id',ids(j)));
%                     cluster_map(h) = Cluster(struct('day_of_week',i,'id',ids(j)),data(day_ind,j));
%                 end
%             end
            
        end
            
   
        
        
    end
    
end

