classdef (Abstract) PredictionJobManager < handle
    
    properties( Access = public )
        
        inputfolder
        outputfolder
        configname
        predictor_params
        
        config
        ni
        days
        
        all_ml_vds
        all_rp_vds
        
        task_filters   % array of function handles used to decide whether or not to execute the task
                       % The task executes if any of the filters evaluates
                       % to true. 
                       % Default: task_filters={}, all tasks execute. 
        
    end
    
    properties( Access = protected )
        tasks_file
        task_filters_locked
    end
    
    methods( Access = public )
        
        function [this] = PredictionJobManager(confname,infolder,outfolder)
            
            this.task_filters = {};
            this.task_filters_locked = false;
            this.inputfolder = infolder;
            this.outputfolder = outfolder;
            this.configname = confname;
            this.config = Config.get(confname);
            this.ni = ObjectFactory.network_information(this.config.xml_file);
            this.days = this.config.good_days;
            
            % all mainline and ramp vdss
            T = this.ni.scenario.get_sensor_table;
            T = T(T.is_attached==1,:);
            this.all_ml_vds = T.vds( strcmp(T.sensor_link_type,'ML') )';
            this.all_rp_vds = T.vds( strcmp(T.sensor_link_type,'OR') | strcmp(T.sensor_link_type,'FR') )';
            
            % initialize
            this.initialize;
            
            % build task list
            this.build_job;
        end
        
        function [predictor,predictor_file] = build_predictor(this,params)
            
            if ~isfield(params,'predictor_class')
                error('~isfield(params,''predictor_class'')')
            end
            
            if ~isfield(params,'measurement_dt')
                params.measurement_dt = 300;
            end
            
            if ~isfield(params,'model_runner')
                params.model_runner = 'beats';
            end
            
            if ~isfield(params,'input_flow_predictor')
                params.input_flow_predictor = 'historical_or_zoh_pems';
            end
            
            if ~isfield(params,'state_estimator')
                params.state_estimator = 'linear_interpolation';
            end
            
            if ~isfield(params,'prediction_horizon')
                params.prediction_horizon = 3600;
            end
            
            if ~isfield(params,'update_dt')
                params.update_dt = 1800;
            end
            
            if ~isfield(params,'start_time')
                params.start_time = 0;
            end
            
            if ~isfield(params,'end_time')
                params.end_time = 86400;
            end
            
            params.configname = this.configname;
            
            % create the predictor
            switch params.predictor_class
                case 'ModelFullPredictor'
                    predictor = ModelFullPredictor(params);
                case 'ModelLessPredictor'
                    predictor = ModelLessPredictor(params);
            end
            
            this.predictor_params = params;
            
            % serialize only if not returning
            if nargout<3
                [~,name] = fileparts(tempname);
                predictor_file = fullfile(this.inputfolder,name);
                save(predictor_file,'predictor')
            else
                predictor_file = [];
            end
            
        end
        
        function [this] = add_filters(this,f)
            if this.task_filters_locked
                error('this.task_filters_locked')
            end
            
            if iscell(f)
                if ~all(cellfun(@(x)isa(x,'function_handle'),f))
                    error('non function handle')
                end
                n = numel(this.task_filters);
                this.task_filters(n+1:n+numel(f)) = f;
            else
                if ~isa(f,'function_handle')
                    error('non function handle')
                end
                this.task_filters{end+1} = f;
            end
        end
        
        function [this] = clear_filters(this)
            if this.task_filters_locked
                error('this.task_filters_locked')
            end
            this.task_filters = {};
        end
        
        function [tasks] = get_filtered_tasks(this)
            
            tasks = [];
            
            % load tasks
            load(this.tasks_file);
            
            % apply filters
            if ~isempty(this.task_filters)
                filter_result = nan( numel(tasks) , length(this.task_filters) );
                for i=1:length(this.task_filters)
                    filter_result(:,i) = arrayfun( @(t) this.task_filters{i}(t) , tasks );
                end
                filter_result = any(filter_result,2);
                tasks = tasks(filter_result);
            end
            
        end
        
        function [n] = get_number_filtered_tasks(this)
            n = numel(this.get_filtered_tasks);
        end
        
        function [results_files,tasks] = execute_job_parallel_build(this,pred_params)
            
            fprintf('Executing tasks in parallel, building predictors\n')
            
            % load tasks
            tasks = this.get_filtered_tasks;
            
            % lock task list
            this.task_filters_locked = true;

            build = @this.build_predictor;
            myparsave = @this.parsave;
            results_files = cell(1,length(tasks));
            
            % build one predictor so that all mat files get generated
            build(pred_params);
            
            parfor i=1:length(tasks)
                
                fprintf('\ttask %d of %d\n',i,length(tasks))
                
                % build the predictor
                predictor = build(pred_params);
                
                % set and run
                day = tasks(i).day;
                remove_vdss = [tasks(i).remove_ml_vdss tasks(i).remove_rp_vdss];
                predictor.set_day(day);
                predictor.set_force_bad_sensors(day,remove_vdss);
                predictor.run;
                
                % save to file
                results_files{i} = myparsave(struct( ...
                    'log',predictor.log, ...
                    'day',day, ...
                    'remove_vdss',remove_vdss ) );
                
            end
                        
        end
        
        function [results_files,tasks] = execute_job_parallel(this,predictor_file)
            
            fprintf('Executing tasks in parallel, loading predictors\n')
            
            if ~exist([predictor_file '.mat'],'file')
                error('predictor file not found')
            end
                        
            % load tasks
            tasks = this.get_filtered_tasks;
            results_files = cell(1,length(tasks));

            % lock task list
            this.task_filters_locked = true;
            
            parfor i=1:length(tasks)
                
                fprintf('\ttask %d of %d\n',i,length(tasks))
                
                % load the predictor
                fprintf('Loading predictor for task %d\n',i);
                P = load(predictor_file);
                
                % set and run
                day = tasks(i).day;
                remove_vdss = [tasks(i).remove_ml_vdss tasks(i).remove_rp_vdss];
                P.predictor.set_day(day);
                P.predictor.set_force_bad_sensors(day,remove_vdss);
                P.predictor.run;
                
                % save to file
                results_files{i} = PredictionJobManager.parsave(struct( ...
                    'log',P.predictor.log, ...
                    'day',day, ...
                    'remove_vdss',remove_vdss ) );
                
            end
            
        end
        
        function [results_files,tasks] = execute_job_series(this,params)
            
            fprintf('Executing tasks in series\n')
            
            % load tasks
            tasks = this.get_filtered_tasks;
            
            % lock task list
            this.task_filters_locked = true;
            
            results_files = cell(1,length(tasks));
            
            % create the predictor
            predictor = this.build_predictor(params);
            
            for i=1:length(tasks)
                
                fprintf('\ttask %d of %d\n',i,length(tasks))
                
                % reset predictor
                predictor.reset;
                
                % set and run
                day = tasks(i).day;
                remove_vdss = [tasks(i).remove_ml_vdss tasks(i).remove_rp_vdss];
                predictor.set_day(day);
                predictor.set_force_bad_sensors(day,remove_vdss);
                predictor.run;
                
                % save to file
                results_files{i} = this.parsave(struct( ...
                    'log',predictor.log, ...
                    'day',day, ...
                    'remove_vdss',remove_vdss ) );
                
            end
            
        end
        
    end
    
    methods( Access = protected )
        
        function [tasks_file]=save_job(this,tasks) %#ok<INUSD>
            [~,name] = fileparts(tempname);
            tasks_file = fullfile(this.inputfolder,name);
            save(tasks_file,'tasks')
        end
        
        function [name]=parsave(this,X) %#ok<INUSD>
            name = sprintf('out_%d',randi(1e10,1));
            save(fullfile(this.outputfolder,name),'X')
        end
        
        function [dty_vpm]=process_log_estimation(this,log)
            L = repmat(this.ni.link_lengths_meters/1609.34,48,1);
            dty_vpm = vertcat(log.initial_state)./L;
        end
        
        function [prediction_data,prediction_error,one_hour_flow_error,vds,prediction_demand,prediction_sr,prediction_boundaryflow] = process_log_prediction(this,error_data_provider,day,log)
            % [prediction_data,prediction_error,one_hour_flow_error] = process_log(this,error_data_provider,day,log)
            % prediction_data(i) ... data for predictions of length prediction_data(i).horizon
            % prediction_data(i).time ... array of time instances marking the end of a prediction interval
            % prediction_data(i).flw_vph(v,k) ... predicted flows for time prediction_data(i).time(k). This is the flow predicted for time inteval [time(k)-5min,time(k)]
            % prediction_data(i).dty_vpm(v,k) ... predicted density for time
            % prediction_data(i).time(k). This is the density predicted for time instant time(k)
            % one_hour_flow_error
            % preliminary
            good_ml_vds = this.good_ml_vds_for_day{day==this.days};
            good_ml_links = this.ni.get_link_ids_for_vds(good_ml_vds);
            measurement_dt = this.predictor_params.measurement_dt;
            start_time = this.predictor_params.start_time;
            end_time = this.predictor_params.end_time;
            prediction_horizon = this.predictor_params.prediction_horizon;
            num_pred_meas = prediction_horizon/measurement_dt;
            num_predictions = (end_time-start_time)/this.predictor_params.update_dt;
            
            % get measured data from data provider
            measured_data = error_data_provider.get_data( ...
                day, ...
                good_ml_vds, ...
                start_time , ...
                end_time , ...
                this.predictor_params.measurement_dt, 'snap');
            measured_time = measured_data(1).time(2:end);
            measured_dty_vpm = vertcat(measured_data.dty_vpk)*1.60934;
            measured_flw_vph = vertcat(measured_data.flw_out_vph);
            measured_health = error_data_provider.get_health(day,good_ml_vds);
            clear measured_data eval_meas_ids
            
            % remove bad detectors
            num_sensors = sum(measured_health);
            good_ml_vds(~measured_health) = [];
            good_ml_links(~measured_health) = [];
            measured_dty_vpm(~measured_health,:) = [];
            measured_flw_vph(~measured_health,:) = [];
            clear measured_health
            
            vds.ids = good_ml_vds;
            vds.link_ids = good_ml_links;
            
            % organize into matrices per prediction time interval
            prediction_data = repmat(struct( ...
                'horizon',nan, ...
                'time', [] , ...
                'flw_vph', nan(num_sensors,num_predictions) , ...
                'dty_vpm', nan(num_sensors,num_predictions) ) , ...
                1,num_pred_meas);
            
            prediction_times = [log.start_time];
            for j=1:num_pred_meas
                prediction_data(j).horizon = j*measurement_dt;
                prediction_data(j).time = prediction_times + j*measurement_dt;
            end
            
            link_ind = Utils.index_into(good_ml_links,this.ni.link_ids);
            for i=1:num_predictions
                for j=1:num_pred_meas
                    % prediction of state at time i, made at time i-j.
                    prediction_data(j).flw_vph(:,i) = log(i).predicted_state.flw_vph(j,link_ind)';
                    prediction_data(j).dty_vpm(:,i) = log(i).predicted_state.dty_vpm(j+1,link_ind)';
                end
            end
            
            %save predicted boundary flows
            link_ids = this.ni.scenario.get_link_ids;
            is_source = this.ni.scenario.is_source_link;
            is_sink = this.ni.scenario.is_sink_link;
            source_link_ids = link_ids(is_source);
            sink_link_ids = link_ids(is_sink);
            source_link_types = this.ni.scenario.get_link_types(is_source);
            sink_link_types = this.ni.scenario.get_link_types(is_sink);
            
            rp_links=[source_link_ids,sink_link_ids];
            rp_link_type=[source_link_types,sink_link_types];
            num_rp=length(rp_links);

            prediction_boundaryflow = repmat(struct( ...
                'link_id', nan,...
                'link_type', nan,...
                'horizon', nan, ...
                'time', [] , ...
                'flw_vph', nan(num_pred_meas,num_predictions) , ...
                'dty_vpm', nan(num_pred_meas,num_predictions) ) , ...
                1,num_rp);
            
            rp_ind = Utils.index_into(rp_links,this.ni.link_ids);
            for i=1:num_rp
                prediction_boundaryflow(i).link_id=rp_links(i);
                prediction_boundaryflow(i).link_type=rp_link_type(i);
                prediction_boundaryflow(i).horizon=measurement_dt*(1:1:num_pred_meas);
                prediction_boundaryflow(i).time=prediction_times;
                
                for j=1:num_predictions
                    prediction_boundaryflow(i).flw_vph(:,j) = log(j).predicted_state.flw_vph(:,rp_ind(i));
                    prediction_boundaryflow(i).dty_vpm(:,j) = log(j).predicted_state.dty_vpm(2:end,rp_ind(i));
                end
            end
            
            %save predicted demands
            prediction_demand=[];
            if isfield(log,'predicted_demands')              
                num_ramp=size(log(1).predicted_demands,2);
                prediction_demand = repmat(struct( ...
                    'link_id',nan, ...
                    'horizon',nan, ...
                    'time', [] , ...
                    'flw_vph', nan(num_pred_meas,num_predictions)),...
                    1,num_ramp);
                for i=1:num_ramp
                    prediction_demand(i).link_id=log(1).predicted_demands(i).link_id;
                    prediction_demand(i).horizon=measurement_dt*(1:1:num_pred_meas);
                    prediction_demand(i).time=prediction_times;
                    for j=1:num_predictions
                        prediction_demand(i).flw_vph(:,j)=log(j).predicted_demands(i).flw_vph';
                    end
                end
            end
            
            %save predicted split ratios, only consider split ratios for
            %the first out-link
            prediction_sr=[];
            if isfield(log,'predicted_splits')
                num_diverge=size(log(1).predicted_splits,2);
                prediction_sr = repmat(struct( ...
                    'node_id',nan, ...
                    'links_in',nan, ...
                    'links_out',nan,...
                    'horizon',nan, ...
                    'time', [] , ...
                    'split_ratios', nan(num_pred_meas,num_predictions)),...
                    1,num_diverge);
                for i=1:num_diverge
                    prediction_sr(i).node_id=log(1).predicted_splits(i).node_id;
                    prediction_sr(i).links_in=log(1).predicted_splits(i).links_in;
                    prediction_sr(i).links_out=log(1).predicted_splits(i).links_out;
                    prediction_sr(i).horizon=measurement_dt*(1:1:num_pred_meas);
                    prediction_sr(i).time=prediction_times;
                    for j=1:num_predictions
                        prediction_sr(i).split_ratios(:,j)=cell2mat(log(j).predicted_splits(i).split_ratios(1))';
                    end
                end
            end
            
            % record the averaege of the predicted flow for each link
            pems_data = nan(num_predictions,numel(good_ml_links));
            pred_data = nan(num_predictions,numel(good_ml_links));
            for i=1:num_predictions
                if size(log(i).predicted_state.flw_vph,1)~=12
                    error('size(log(i).predicted_state.flw_vph,1)~=12') % restrict to one hour predictions
                end
                time = log(i).start_time + (measurement_dt:measurement_dt:prediction_horizon);
                meas_ind = Utils.index_into(time,measured_time);
                if any(meas_ind==0)
                    continue
                end
                pems_data(i,:) = Utils.meanwithnan(measured_flw_vph(:,meas_ind),2)';
                pred_data(i,:) = Utils.meanwithnan(log(i).predicted_state.flw_vph(:,link_ind),1);
            end
            one_hour_flow_error = abs(pems_data-pred_data)./pems_data;
            clear pems_data pred_data meas_ind time
            
            % compare each prediction to data
            prediction_error.flw_vph_Percent = nan(1,num_pred_meas);
            
            prediction_error.flw_vph_RMSD = nan(1,num_pred_meas);
            prediction_error.dty_vpm_RMSD = nan(1,num_pred_meas);
            prediction_error.spd_mph_RMSD = nan(1,num_pred_meas);
            
            prediction_error.flw_vph_MAE = nan(1,num_pred_meas);
            prediction_error.dty_vpm_MAE = nan(1,num_pred_meas);
            prediction_error.spd_mph_MAE = nan(1,num_pred_meas);
            
            prediction_error.flw_vph_Linf = nan(1,num_pred_meas);
            prediction_error.dty_vpm_Linf = nan(1,num_pred_meas);
            prediction_error.spd_mph_Linf = nan(1,num_pred_meas);
            
            for j=1:num_pred_meas
                
                ind = Utils.index_into(prediction_data(j).time,measured_time);
                
                % flw
                p_flw = prediction_data(j).flw_vph(:,ind>0);
                m_flw = measured_flw_vph(:,ind(ind>0));
                
                p_dty = prediction_data(j).dty_vpm(:,ind>0);
                m_dty = measured_dty_vpm(:,ind(ind>0));
                
                % reshape
                n = numel(p_flw);
                p_flw = reshape(p_flw,1,n);
                m_flw = reshape(m_flw,1,n);
                p_dty = reshape(p_dty,1,n);
                m_dty = reshape(m_dty,1,n);
                
                % speed
                p_spd = p_flw./p_dty;
                m_spd = m_flw./m_dty;
                
                % limit predicted speed
                p_spd(p_spd>80) = 80;
                
                % detect suspicious data
                bad_data = m_spd>80 | ...
                    isnan(m_flw) | isinf(m_flw) | ...
                    isnan(m_dty) | isinf(m_dty) | ...
                    isnan(m_spd) | isinf(m_spd) ;
                
                % remove suspicious data
                p_flw = p_flw(~bad_data);
                m_flw = m_flw(~bad_data);
                p_dty = p_dty(~bad_data);
                m_dty = m_dty(~bad_data);
                p_spd = p_spd(~bad_data);
                m_spd = m_spd(~bad_data);
                
                % make column vectors for Utils.meanwithnan
                d_flw = Utils.column_vector(m_flw-p_flw);
                d_dty = Utils.column_vector(m_dty-p_dty);
                d_spd = Utils.column_vector(m_spd-p_spd);
                
                %calculate the link flow percentage error
                if prediction_data(j).horizon<=900
                    prediction_error.flw_vph_Percent(j)=mean(abs(d_flw./(m_flw')*100)<15)*100;
                elseif prediction_data(j).horizon<=1800
                    prediction_error.flw_vph_Percent(j)=mean(abs(d_flw./(m_flw')*100)<20)*100;
                elseif prediction_data(j).horizon<=2700
                    prediction_error.flw_vph_Percent(j)=mean(abs(d_flw./(m_flw')*100)<25)*100;
                else
                    prediction_error.flw_vph_Percent(j)=mean(abs(d_flw./(m_flw')*100)<30)*100;
                end
                
                %RMSD
                prediction_error.flw_vph_RMSD(j) = sqrt(Utils.meanwithnan(d_flw.^2));
                prediction_error.dty_vpm_RMSD(j) = sqrt(Utils.meanwithnan(d_dty.^2));
                prediction_error.spd_mph_RMSD(j) = sqrt(Utils.meanwithnan(d_spd.^2));
                
                %MAE
                prediction_error.flw_vph_MAE(j) = Utils.meanwithnan(abs(d_flw));
                prediction_error.dty_vpm_MAE(j) = Utils.meanwithnan(abs(d_dty));
                prediction_error.spd_mph_MAE(j) = Utils.meanwithnan(abs(d_spd));
                
                %L_inf
                prediction_error.flw_vph_Linf(j) = maxwithnan(abs(d_flw));
                prediction_error.dty_vpm_Linf(j) = maxwithnan(abs(d_dty));
                prediction_error.spd_mph_Linf(j) = maxwithnan(abs(d_spd));
            end
            
        end

    end
    
    methods( Abstract, Access=public )
        [this] = process_result(this);
        [this] = report_result(this,output_file);
    end
    
    methods( Abstract, Access=protected )
        [this] = build_job(this,params);
        [this] = initialize(this);
    end
    
    methods( Static )
        
        function [z] = extract_flow_error(data)
            error_all = horzcat(data.flw_error);
            as_matrix = error_all(1:end-1,:)';
            as_vector = reshape(as_matrix,numel(as_matrix),1);
            as_vector(isnan(as_vector)) = [];
            
            z.as_matrix = as_matrix;
            z.as_vector = as_vector;
        end
        
    end
    
end

