classdef ScaledHistoricalFlowPredictor < HistoricalFlowPredictor

    properties
        simulation_data_provider @SimDataProvider
    end
    
    methods( Access = public )
        
        function [this]=ScaledHistoricalFlowPredictor(simulation_data_provider,pems_data_provider)
            this = this@HistoricalFlowPredictor(pems_data_provider);
            this.simulation_data_provider = simulation_data_provider;
        end
        
        function [y_scaled] = predict(this,ids,days,from,to,dt)            

            if isnan(dt)
                dt = this.data_provider.measurement_dt;
            end
            
            % get historical data
            y_pred = predict@HistoricalFlowPredictor(this,ids,days,from,to,dt);
            
            % get link ids on simulation network
            link_ids = this.data_provider.link_ids(Utils.index_into(ids,this.data_provider.ids));    
            
            % get flows over previous time interval
            if from>=dt
                y_curr = this.simulation_data_provider.get_data(days,link_ids,from-dt,from);
            else
                y_curr = this.simulation_data_provider.get_data(days,link_ids,from,from+dt);
            end
                        
            % scal historical to meet historical at
            y_scaled = repmat(DataProfile,1,length(ids));
            for i=1:length(ids)
                alpha = y_curr(i).flw_out_vph(1)/y_pred(i).flw_out_vph(1);
                if isnan(alpha) || isinf(alpha)
                    alpha = 1;
                end
                y_scaled(i) = DataProfile( y_pred(i).id , nan , ...
                    y_pred(i).time, [] , alpha*y_pred(i).flw_out_vph , [] );   
            end
            
        end
        
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            
            config = Config.get(params.config);

            ni = ObjectFactory.network_information(config.xml_file);  
			pems_dp = Utils.get_pems_dp(ni,params.config);
            mr      = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            sim_dp  = ObjectFactory.sim_data_provider(mr,params.end_time);
            fp      = ObjectFactory.scaled_historical_predictor(sim_dp,pems_dp);
            
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

