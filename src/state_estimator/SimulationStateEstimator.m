classdef SimulationStateEstimator < StateEstimator
    
    methods( Access = public )
        
        function [this] = SimulationStateEstimator(model_runner,data_provider)
            this = this@StateEstimator(model_runner,data_provider,model_runner.network_information.link_ids);
        end
        
        function [this] = set_day(this,this_day)
            this.set_day@StateEstimator(nan)
        end
        
        function update_estimate(this,new_time)
            if new_time<this.xhat_time
                error('new_time<this.xhat_time')
            end
            if new_time==this.xhat_time
                return
            end
            state = this.data_provider.get_data(nan,this.link_ids,new_time,new_time);
            this.xhat = this.data_provider.extract_density_veh(state);
            this.xhat_time = new_time;
        end
        
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,configfile,sim_dt,output_dt,end_time,update_dt}

            config  = Config.get(params.config);
            ni      = ObjectFactory.network_information(config.xml_file);  
            mr      = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            sim_dp  = ObjectFactory.sim_data_provider(mr,params.end_time);
            se      = ObjectFactory.simulation_state_estimator(mr,sim_dp);
                                      
            X = run_and_report_error@StateEstimator( struct( ...
                    'ppt_file'          , params.ppt_file , ...
                    'xls_file'          , params.ppt_file , ...
                    'network_info'      , ni, ...
                    'state_estimator'   , se , ...
                    'update_dt'         , params.update_dt ));
        end
        
    end
    
end