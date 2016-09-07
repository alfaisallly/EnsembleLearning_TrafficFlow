classdef SimVDSDataProvider < DataProvider
    
    properties
        sim_dt
        configfile
        end_time
    end
    
    methods
        
        function [this] = SimVDSDataProvider(mr,end_time)
                    
            % collect data
            vds2linkid = mr.get_vds2linkid;
            vdss = vds2linkid(:,1)';
            link_ids = vds2linkid(:,2)';
            link_lengths = mr.network_information.get_lengths_for_linkids_in_km(link_ids);

            % super constructor
            this = this@DataProvider(nan,vdss,link_ids,link_lengths,mr.out_dt,'');
            this.configfile = mr.network_information.configfile;
            this.end_time = end_time;

            % this model should export every time step, so that the 
            % density estimate is correct (not the avg over the output dt)
            mr.set_output_dt(mr.sim_dt);            
            
            % run until end_time
            mr.reset();
            mr.run_simulation(end_time);
                        
            % store sampled simulation result in DataProvider.data
            data_time = 0:this.measurement_dt:end_time;
                        
            sim_per_data = this.measurement_dt/mr.sim_dt;
            this.time = data_time;
            time_ind = reshape(1:mr.get_num_time-1,sim_per_data,length(data_time)-1)';
            for i=1:length(link_ids)
                X = mr.get_data_for_link_id(link_ids(i));
                this.data(i).flw_in_vph = mean(X.flw_in_vph(time_ind),2);
                this.data(i).flw_out_vph = mean(X.flw_out_vph(time_ind),2);
                this.data(i).dty_vpk = X.dty_vpm([time_ind(:,1);end]) / 1.60934;
            end
            
        end
        
        function [b] = check_have_day(this,day)
            b = true;
        end
        
        function [x] = get_health(~,getdays,getvds)
            if ischar(getdays)
                getdays = 1;
            end
            x = true(length(getdays),length(getvds));
        end
        
        function [x] = uses_template(~,getdays,getvds)
            x = false(length(getdays),length(getvds));
        end
        
    end
    
end

