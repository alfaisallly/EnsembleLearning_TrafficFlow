classdef ModelFullPredictor < Predictor
    
    properties
        model_runner            @ModelRunner
        demand_predictor        @DemandPredictor
        split_predictor         @SplitPredictor
        state_estimator         @StateEstimator
        sim_dt
        
        flow_predictor          @FlowPredictor
       
        estimation_error
        demand_error
        split_error
    end
    
    methods(Access=public)
        
        function [this] = ModelFullPredictor(params)
            this = this@Predictor(params);
        end
        
        function [this] = reset(this)
            this.reset@Predictor();
            this.estimation_error = [];
            this.demand_error = [];
            this.split_error = [];
        end
        
        function [this] = set_day(this,this_day)
            this = this.set_day@Predictor(this_day);
            this.demand_predictor.set_day(this_day);
            this.split_predictor.set_day(this.ni,this_day);
            this.state_estimator.set_day(this_day);
        end
        
        function [this] = set_force_bad_sensors(this,day,ids)
            this = this.set_force_bad_sensors@Predictor(day,ids);
            this.demand_predictor.set_force_bad_sensors(day,ids);
            this.split_predictor.set_force_bad_sensors(day,ids);
            this.state_estimator.set_force_bad_sensors(day,ids);
        end
        
        function [this]=run(this)
            
            if isempty(this.day) || isnan(this.day)
                error('isempty(this.day) || isnan(this.day)')
            end
            
            % run
            current_time = this.start_time;
            
            % get unique dp_ids for flow prediction          
            dpids_unique = unique([this.split_predictor.all_dps horzcat(this.demand_predictor.dp_ids{:})]);
            
            for counter=1:this.num_predictions
                
                fprintf('%d : %d\n',current_time,this.end_time)
                
                % reset the model runner
                this.model_runner.reset;
                
                % update the state estimate (current state in vehicles per
                % link)
                this.state_estimator.update_estimate(current_time);
                current_state = this.state_estimator.get_state_for_link_ids( ...
                    this.ni.link_ids ); % veh/link
                
                % package and send to runner
                this.model_runner.set_state_veh(this.ni.link_ids,current_state);

                %% revised sequence of function calls: flow-prediction -> demand/split-ratio prediction
                % flow prediction       
                flow_pred=repmat(struct(...
                    'dp_ids', nan,...
                    'flow', nan,...
                    'time', nan),1,length(dpids_unique));
                for i=1:length(dpids_unique);
                    flow_pred(i).dp_ids=dpids_unique(i);
                    [flow_pred(i).flow,flow_pred(i).time] = ...
                        this.flow_predictor.predict_flw_for_dpids( ...
                            dpids_unique(i) , ...
                            this.day , ...
                            current_time , ...
                            current_time + this.prediction_horizon , ...
                            nan );
                end

                % predict demand
                demand_prediction = this.demand_predictor.predict_given_flows(flow_pred);
                  
                % package and send to runner
                for i=1:length(demand_prediction)
                    demand = demand_prediction(i);
                    this.model_runner.set_demand( demand.link_id , ...
                        demand.time(2)-demand.time(1), ...
                        demand.flw_vph );
                end
                
                % predict splits
                split_prediction = this.split_predictor.predict_given_flows(flow_pred);
        
                % package and send to runner
                for i=1:length(split_prediction)
                    split = split_prediction(i);
                    
                    sr_profiles = repmat(struct('link_in',nan,'link_out',nan,'profile',[]),1,numel(split.split_ratios));
                    for j=1:length(sr_profiles)
                        sr_profiles(j).link_in = split.links_in;
                        sr_profiles(j).link_out = split.links_out(j);
                        sr_profiles(j).profile = split.split_ratios{j};
                    end
                    
                    this.model_runner.set_split_ratio( ...
                        split.node_id, ...
                        0, ...
                        split.time(2)-split.time(1), ...
                        sr_profiles);
                end

                % run the model
                this.model_runner.run_simulation(this.prediction_horizon);
                
                % get future states
                state_prediction = this.model_runner.get_outputs();
                
                % log
                this.log(counter).start_time = current_time;
                this.log(counter).predicted_state = state_prediction;
                this.log(counter).initial_state = Utils.row_vector(current_state);
                this.log(counter).predicted_demands = demand_prediction;
                this.log(counter).predicted_splits = split_prediction;
                
                %  advance the clock
                current_time = current_time + this.update_dt;
                
            end
            
        end
        
    end
    
    methods(Access=protected)
        
        function [this]=set_values(this,params)
            
            this.set_values@Predictor(params);
            
            % load configuration information
            config = Config.get(params.configname);
            this.sim_dt = config.sim_dt;
            
            %             switch params.error_data_provider
            %                 case 'simulation'
            %                     this.measurement_dt = config.sim_dt;
            %                 case 'pems'
            %                     this.measurement_dt = PeMSDataProvider.MEASUREMENT_DT;
            %                 otherwise
            %                     error('bad data provider')
            %             end
            this.num_pred_meas = this.prediction_horizon / this.measurement_dt;
            
        end
        
        function [] = check_values(this)
            
            this.check_values@Predictor();
            
            if mod(this.measurement_dt,this.sim_dt)~=0
                error('mod(this.measurement_dt,this.sim_dt)~=0')
            end
            
        end
        
        function [this]=instantiate_objects(this,params)
            
            this.instantiate_objects@Predictor(params);
            
            config = Config.get(params.configname);
            
            % who needs a historical pems data provider?
            if  strcmp(params.input_flow_predictor,'historical') || ...
                    strcmp(params.input_flow_predictor,'historical_or_zoh_pems') || ...
                    strcmp(params.input_flow_predictor,'scaled_historical') || ...
                    strcmp(params.input_flow_predictor,'historical_or_zoh_sim') || ...
                    strcmp(params.input_flow_predictor,'zoh_pems') || ...
                    strcmp(params.state_estimator,'linear_interpolation') || ...
                    strcmp(params.input_flow_predictor,'recursiveARMAX')
                
                pems_dp = Utils.get_pems_dp(this.ni,this.configname);
            end
            
            % who needs a simulation data provider?
            if  strcmp(params.input_flow_predictor,'simulation') || ...
                    strcmp(params.input_flow_predictor,'zoh_sim') || ...
                    strcmp(params.input_flow_predictor,'historical_or_zoh_sim') || ...
                    strcmp(params.input_flow_predictor,'scaled_historical') || ...
                    strcmp(params.state_estimator,'simulation')
                mr = ObjectFactory.beats_model_runner( ...
                    this.ni , ...
                    this.sim_dt , ...
                    this.measurement_dt );
                sim_dp  = ObjectFactory.sim_data_provider(mr,params.end_time);
                clear mr
            end
            
            % who needs a simulation vds data provider?
            if strcmp(params.input_flow_predictor,'zoh_sim_vds') || ...
                    strcmp(params.state_estimator,'linear_interpolation_sim')
                mr = ObjectFactory.beats_model_runner( ...
                    this.ni , ...
                    this.sim_dt , ...
                    this.measurement_dt );
                sim_vds_dp  = ObjectFactory.sim_vds_data_provider(mr,params.end_time);
                clear mr
            end
            
            % make a flow predictor to pass to demand and split prediction
            % This uses the data provider created above.
            switch params.input_flow_predictor
                case 'simulation'
                    input_flow_predictor = ObjectFactory.simulation_flow_predictor(sim_dp);
                case 'zoh_sim'
                    input_flow_predictor = ObjectFactory.zoh_predictor_sim(sim_dp);
                case 'zoh_sim_vds'
                    input_flow_predictor = ObjectFactory.zoh_predictor_sim(sim_vds_dp);
                case 'zoh_pems'
                    input_flow_predictor = ObjectFactory.zoh_predictor_pems(pems_dp);
                case 'historical'
                    input_flow_predictor = ObjectFactory.historical_predictor(pems_dp);
                case 'historical_or_zoh_pems'
                    input_flow_predictor = ObjectFactory.historical_or_zoh_predictor(pems_dp);
                case 'historical_or_zoh_sim'
                    input_flow_predictor = ObjectFactory.historical_or_zoh_predictor(pems_dp,sim_dp);
                case 'scaled_historical'
                    input_flow_predictor = ObjectFactory.scaled_historical_predictor(sim_dp,pems_dp);
                case 'recursiveARMAX'
                    input_flow_predictor = ObjectFactory.recursiveARMAX_predictor(pems_dp,params.armax_params);
                otherwise
                    error('unknown input_flow_predictor')
            end
            % instantiate flow predictor
            this.flow_predictor = input_flow_predictor;
            
            % instantiate demand predictor
            this.demand_predictor = DemandPredictor(this.ni,input_flow_predictor);
            
            % instantiate split predictor
            this.split_predictor = SplitPredictor(this.ni,input_flow_predictor);
            
            % instantiate model runner
            switch params.model_runner
                case 'beats'
                    this.model_runner = ObjectFactory.beats_model_runner( ...
                        this.ni, ...
                        this.sim_dt , ...
                        this.measurement_dt );
                otherwise
                    error('unknown model_runner')
            end
            
            % instantiate state estimator
            switch params.state_estimator
                case 'simulation'
                    this.state_estimator = ObjectFactory.simulation_state_estimator(this.model_runner,sim_dp);
                case 'linear_interpolation'
                    fwy_info = ObjectFactory.freeway_information(config.xml_file);
                    this.state_estimator = ObjectFactory.lininterp_state_estimator(pems_dp,fwy_info);
                case 'linear_interpolation_sim'
                    fwy_info = ObjectFactory.freeway_information(config.xml_file);
                    this.state_estimator = ObjectFactory.lininterp_state_estimator(sim_vds_dp,fwy_info);
                otherwise
                    error('unknown state estimator')
            end
            
        end
        
    end
    
    methods(Access=protected,Static)
        
        function []=check_inputs(params)
            
            check_inputs@Predictor(params)
            
            if ~isfield(params,'model_runner')
                error('~isfield(params,''model_runner'')')
            end
            
            if ~isfield(params,'input_flow_predictor')
                error('~isfield(params,''input_flow_predictor'')')
            end
            
            if ~isfield(params,'state_estimator')
                error('~isfield(params,''state_estimator'')')
            end
            
        end
    end
end
