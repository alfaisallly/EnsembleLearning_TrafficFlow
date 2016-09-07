classdef FreewayInformation < NetworkInformation
    
    properties
        ordered_vds_table           % array of structures with mainline and ramp stations
        ordered_ml_vds              % list of mainline vds
        ordered_ml_link_ids
        ordered_ml_link_lengths_km
        ordered_vds2link_ind        % into ordered_ml_link_ids and ordered_ml_link_lengths_km
    end
    
    methods
        
        function [this] = FreewayInformation(configfile)
            this = this@NetworkInformation(configfile);
            this.ordered_vds_table = this.scenario.get_ordered_vds;
            
            all_vds = [this.ordered_vds_table.sensor_vds];
            ml_ind = strcmpi({this.ordered_vds_table.link_type},'freeway');
            this.ordered_ml_vds = all_vds(ml_ind);
            
            T = this.scenario.get_freeway_structure();
            this.ordered_ml_link_ids = T.linear_ml_link_ids';

            map = [this.ordered_vds_table.link_id;this.ordered_vds_table.sensor_vds]';
            vds_link_ids = map(Utils.index_into(this.ordered_ml_vds,map(:,2)),1);
            this.ordered_vds2link_ind = Utils.index_into(vds_link_ids,this.ordered_ml_link_ids)';
            
            if any(this.ordered_vds2link_ind==0)
                error('any(this.ordered_vds2link_ind==0)')
            end
            
            all_lengths = this.scenario.get_link_lengths('si')/1000;
            this.ordered_ml_link_lengths_km = all_lengths(T.linear_ml_link_ind);
            
        end
        
        function [x] = get_ordered_ml_links(this)
            x = this.ordered_ml_link_ids;
        end
                
        function [x]=get_ordered_mainline_vds(this)
            x = this.ordered_ml_vds;
        end
        
        function [x] = get_ordered_mainline_vds_lengths_km(this)
            x = this.ordered_ml_link_lengths_km(this.ordered_vds2link_ind);
        end
        
    end
    
end

