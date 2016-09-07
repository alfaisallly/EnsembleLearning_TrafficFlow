classdef HistoricalOrZOHFlowPredictor < FlowPredictor
    
    properties
        zoh_predictor         @ZOHFlowPredictor
        historical_predictor  @HistoricalFlowPredictor
    end

    methods
        
        function [this]=HistoricalOrZOHFlowPredictor(pems_dp,sim_dp)
            this = this@FlowPredictor(pems_dp);
            this.historical_predictor = ObjectFactory.historical_predictor(pems_dp);
            
            if isa(sim_dp,'SimDataProvider')
                this.zoh_predictor = ObjectFactory.zoh_predictor_sim(sim_dp);
            else
                this.zoh_predictor = ObjectFactory.zoh_predictor_pems(pems_dp);
            end
                
        end
        
        function [y] = predict(this,ids,day,from,to,dt)
            % assume that the ids passed here are for the data_provider,
            % which is the same as historical_predictor.data_provider (and
            % also is a PeMS data provider). ie ids are vdss. b
            % This may be different from zoh_predictor.data_provider
            % which could be either PeMS or simulation. 

            if numel(day)~=1
                error('numel(day)~=1')
            end
            
            % get health and template use from pems
            ids_health = this.data_provider.get_health(day,ids);
            uses_template = this.data_provider.uses_template(day,ids);
                        
            % get the ids for the zoh
            switch class(this.zoh_predictor.data_provider)
                case 'PeMSDataProvider'
                    % use vdss
                    zoh_ids = ids;
                case 'SimDataProvider'
                    % use linke ids
                    zoh_ids = this.data_provider.get_linkids_for_ids(ids);
            end
            
            % loop through ids, if it is good use zoh, otherwise use historical 
            y = repmat(DataProfile,1,length(ids));
            for i=1:length(ids)
                if ids_health(i) & ~uses_template(i)
                    y(i) = this.zoh_predictor.predict(zoh_ids(i),day,from,to,dt);
                else
                    y(i) = this.historical_predictor.predict(ids(i),day,from,to,dt);
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

