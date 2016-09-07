classdef (Abstract) ModelRunner < handle
    
    properties
        sim_dt
        out_dt
        network_information @NetworkInformation
    end
    
    methods(Access=public)
        
        function [this] = ModelRunner(network_information,sim_dt,out_dt)
            this.network_information = network_information;
            this.sim_dt = sim_dt;
            this.out_dt = out_dt;
        end
        
        function [X] = get_outputs(this)
            link_ids = this.network_information.link_ids;
            n = this.get_num_time;
            X.flw_vph = nan(n-1,length(link_ids));
            X.dty_vpm = nan(n,length(link_ids));
            for i=1:length(link_ids)
                p = this.get_data_for_link_id(link_ids(i));
                X.flw_vph(:,i) = p.flw_out_vph;
                X.dty_vpm(:,i) = p.dty_vpm;
            end
        end
             
        function [dt] = get_output_dt(this)
            dt = this.out_dt;
        end
                
    end
    
    methods( Abstract )        
        
        % size(linkid_dtyvpk) numlinks,2
        % first column is link ids, second column is density in vpk
        [this] = set_state_veh(this,given_linkids,dty_veh)      
        [linkid_dtyveh] = get_state_veh(this)     
        
        % size(flw_vph) numval,1
        % start_time, dt in seconds
        [this] = set_demand(this,link_id,dt,flw_vph)
        [start_time,dt,flw_vph] = get_demand(this,link_id)
        
        % sr_profiles is an array of structures, each one with fields
        % link_in ... input link
        % link_out ... output link
        % profile ... split ratio profile
        [this] = set_split_ratio(this,node_id,start_time,dt,sr_profiles)
        [start_time,dt,sr_profiles] = get_split_ratio(this,node_id)
        
        [this] = reset(this)
        [this] = run_simulation(this,duration)
        [that] = clone(this)
        [link_ids] = get_link_ids(this)
        [data] = get_data_for_link_id(this,link_id)  % returns a DataProfile
        [n] = get_num_time(this)    % number of time points currently stored
        
        [this] = set_output_dt(this,dt)
        [vds2linkid] = get_vds2linkid(this);
        
    end
    
end
