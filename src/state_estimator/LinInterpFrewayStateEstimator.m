classdef LinInterpFrewayStateEstimator < FreewayStateEstimator

    properties
       state_link_ind
       state_link_lengths
       
       internal_dp_ids
       link_midpoint_positions % in km, position of the mid point of each mainline link
    end
    
    methods(Access = public)
        
        function [this] = LinInterpFrewayStateEstimator(data_provider,freeway_information)
            this = this@FreewayStateEstimator([],data_provider,freeway_information); 

            % use only data sources that are internal to the network
            internal_link_ids = freeway_information.get_internal_link_ids;
            this.internal_dp_ids = data_provider.get_ids_for_linkids(internal_link_ids); 
            this.internal_dp_ids = [this.internal_dp_ids{:}];
        end

        function [this] = set_day(this,this_day)
            
            set_day@StateEstimator(this,this_day);
                        
            % remove dp_ids that are not good today
            is_good = this.data_provider.get_health(this_day,this.internal_dp_ids);
            acc_dp_ids = this.internal_dp_ids(is_good);
            
            % make them accessible
            this.set_accessible_dp_ids(acc_dp_ids);
            
            % calculate midpoints for each mainline link
            x = cumsum(this.fwy.ordered_ml_link_lengths_km);
            this.link_midpoint_positions = mean([[0 x(1:end-1)];x]);

        end        
        
        function [this] = set_force_bad_sensors(this,this_day,bad_vdss)
            this.set_force_bad_sensors@StateEstimator(this_day,bad_vdss);
            this.set_day(this_day);
        end
        
        function update_estimate(this,new_time)
            if isempty(this.day)
                error('isempty(this.day)')
            end
            if new_time<this.xhat_time
                error('new_time<this.xhat_time')
            end
            if new_time==this.xhat_time
                return
            end
            
            % get the current measured density
            A = this.data_provider.get_data( ...
                                this.day , ...
                                this.accessible_dp_ids,new_time , ...
                                new_time);
            meas_vpm = horzcat(A.dty_vpk)*1.609;
            clear A
                     
            % data positions
            dp_midpoint_positions = this.link_midpoint_positions(this.accessible_link_ind);
            
            % interpolate, translate to veh units
            x_vpm = interp1(dp_midpoint_positions,meas_vpm,this.link_midpoint_positions,'linear','extrap');
            this.xhat = x_vpm.*this.link_lengths_mile;
            this.xhat_time = new_time;
        end
        
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,configfile,sim_dt,output_dt,end_time,update_dt}

            config  = Config.get(params.config);
            fi = ObjectFactory.freeway_information(config.xml_file);
			pems_dp = Utils.get_pems_dp(fi,params.config);                  
            se = ObjectFactory.lininterp_state_estimator(pems_dp,fi);
            se.set_day(config.model_day);
            
            X = run_and_report_error@StateEstimator( struct( ...
                    'ppt_file'          , params.ppt_file , ...
                    'xls_file'          , params.ppt_file , ...
                    'network_info'      , fi , ...
                    'state_estimator'   , se , ...
                    'update_dt'         , params.update_dt ));
        end
        
    end
    
end

