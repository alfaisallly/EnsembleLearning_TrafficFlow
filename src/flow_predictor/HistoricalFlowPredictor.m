classdef HistoricalFlowPredictor < FlowPredictor
    
    methods
        
        function [this]=HistoricalFlowPredictor(data_provider)
            this = this@FlowPredictor(data_provider);  
        end
        
        function [y] = predict(this,ids,day,from,to,dt)
                       
            if numel(day)~=1
                error('numel(day)~=1')
            end
            
            % get the data
            days_and_ids = this.data_provider.get_days_and_ids;
            x = this.data_provider.get_data('all',ids,from,to,dt);
            vds_is_good = this.data_provider.get_health('all',ids);
            
            % remove day from data
            if(length(days_and_ids.days)>1) %more than one day
                ind = days_and_ids.days==day;
                if any(ind)
                    x(ind,:) = [];
                    vds_is_good(ind,:) = [];
                end
                clear days_and_ids
            end
            
            % TEMPORARY HACK!!!
            if numel(ids)==1 && all(~vds_is_good)
%                 warning(sprintf('VDS %d is not good on %d. Setting to good.\n',ids,days))
                vds_is_good = true;
            end
            
            % average over good days
            time = x(1,1).time;            
            y = repmat(DataProfile(),1,length(ids));
            for i=1:length(ids)
                is_good = vds_is_good(:,i); 
                if ~any(is_good)
                    y(i) = DataProfile( ...
                        ids(i) , nan , time, ...
                        nan(1,length(time)-1) , ...
                        nan(1,length(time)-1) , ...
                        nan(1,length(time)) );
                else
                    y(i) = DataProfile( ...
                        ids(i) , nan , time, ...
                        mean(vertcat(x(is_good,i).flw_in_vph),1) , ...
                        mean(vertcat(x(is_good,i).flw_out_vph),1) , ...
                        mean(vertcat(x(is_good,i).dty_vpk),1) );
                end
                
            end
            
        end
        
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,configfile,sim_dt,output_dt,end_time,update_dt,horizon}
                    
            config = Config.get(params.config);

            ni = ObjectFactory.network_information(config.xml_file);  
            pems_dp = Utils.get_pems_dp(ni,params.config);
            fp      = ObjectFactory.historical_predictor(pems_dp);
            mr      = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            sim_dp  = ObjectFactory.sim_data_provider(mr,params.end_time);
            
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

