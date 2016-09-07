classdef ZOHFlowPredictor < FlowPredictor

    methods( Access = public )
        
        function [this]=ZOHFlowPredictor(data_provider)
            this = this@FlowPredictor(data_provider);  
        end
        
        function [x] = predict(this,ids,days,from,to,dt)
            
            if numel(days)>1
                error('numel(days)>1')
            end
            
            if isnan(dt)
                dt = this.data_provider.measurement_dt;
            end
            
            % get latest data, exception if from<dt.
            if from>=dt
                start_time = from-dt;
                end_time = from;
            else
                start_time = from;
                end_time = from+dt;
            end
            
            % get detector health
            vds_is_good = this.data_provider.get_health(days,ids);
            
            time = from:dt:to;
            n = length(time)-1;
            x = repmat(DataProfile,1,length(ids));
            for i=1:length(ids)
                
                % run prediction
                if vds_is_good(i)
                    % use today's data
                    xc = this.data_provider.get_data(days,ids(i),start_time,end_time);
                    x(i) = DataProfile(ids(i),nan,time, ...
                        xc.flw_in_vph(1)*ones(1,n) , ...
                        xc.flw_out_vph(1)*ones(1,n) , ...
                        xc.dty_vpk(end)*ones(1,n+1) );
                else
                    % use average historical
                    xc = this.data_provider.get_data('all',ids(i),start_time,end_time);
                    dty = vertcat(xc.dty_vpk);
                    x(i) = DataProfile(ids(i),nan,time, ...
                        mean([xc.flw_in_vph])*ones(1,n) , ...
                        mean([xc.flw_out_vph])*ones(1,n) , ...
                        mean(dty(:,end))*ones(1,n+1) );
                end
                
            end
        end
        
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            
            config = Config.get(params.config);

            ni      = ObjectFactory.network_information(config.xml_file);  
            mr      = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            sim_dp  = ObjectFactory.sim_data_provider(mr,params.end_time);
            fp      = ZOHFlowPredictor(sim_dp);
            
            X = run_and_report_error@FlowPredictor( struct( ...
                    'ppt_file'          , params.ppt_file , ...
                    'xls_file'          , params.ppt_file , ...
                    'flow_predictor'    , fp , ...
                    'day'               , config.model_day , ...
                    'update_dt'         , params.update_dt, ...
                    'horizon'           , params.horizon, ...
                    'sim_data_provider' , sim_dp ));
                
        end
        
    end
        
end

