classdef SimulationFlowPredictor < FlowPredictor

    methods( Access = public )
        
        function [this]=SimulationFlowPredictor(data_provider)
            this = this@FlowPredictor(data_provider);
            
            %check type of data_provider
            if ~isa(data_provider,'SimDataProvider')
                error('~isa(data_provider,''SimDataProvider'')')
            end
                
        end
        
        function [x] = predict(this,ids,days,from,to,dt)
            x = this.data_provider.get_data(days,ids,from,to,dt);
        end
        
    end
        
    methods(Static)
        
        function [X] = run_and_report_error(params)
            config = Config.get(params.config);
            ni      = ObjectFactory.network_information(config.xml_file);  
            mr      = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            sim_dp  = ObjectFactory.sim_data_provider(mr,params.end_time);
            fp      = SimulationFlowPredictor(sim_dp);
            
            X = run_and_report_error@FlowPredictor( struct( ...
                    'ppt_file'          , params.ppt_file , ...
                    'xls_file'          , params.xls_file , ...
                    'flow_predictor'    , fp , ...
                    'day'               , nan , ...
                    'update_dt'         , params.update_dt, ...
                    'horizon'           , params.horizon, ...
                    'sim_data_provider' , sim_dp ));
                
        end
        
    end
    
end

