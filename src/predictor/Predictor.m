classdef (Abstract) Predictor < handle

    properties
        
        % config info
        configname
        
        ni                      @NetworkInformation
        day
        start_time
        end_time
        update_dt
        prediction_horizon
        
        measurement_dt
        
        num_predictions             % number of predictions
        num_pred_meas               % number of measurements per horizon
        forced_bad_sensors          % list of sensor ids that should be set 
                                    % to 'bad' on the prediction day
        log
        prediction_error
        
    end
    
    methods(Access=public)
        
        function [this] = Predictor(params)

            % check input
            this.check_inputs(params)
            
            % set values
            this.set_values(params);
            
            % input validation
            this.check_values
            
            % instantiate objects
            this.instantiate_objects(params);
            
        end

        function [this] = reset(this)
            this.log = [];
            this.prediction_error = [];
        end
        
        function [this] = set_day(this,this_day)
            this.day = this_day;
            this.reset;
        end
        
        function [this] = set_force_bad_sensors(this,this_day,x)
            if this.day~=this_day
                error('this.day~=this_day')
            end
            this.forced_bad_sensors = x;
        end
        
        function [prediction_data] = compute_prediction_error(this,error_data_provider,pptfile)
            
            % get good internal links and sensors
            [eval_meas_ids,eval_links] = error_data_provider.get_good_internal_sensors(this.ni,this.day);
            
            % get measured data from data provider
            measured_data = error_data_provider.get_data( ...
                this.day, ...
                eval_meas_ids, ...
                this.start_time , ...
                this.end_time , this.measurement_dt, 'snap');
            measured_time = measured_data(1).time(1:end-1);
            measured_dty_vpm = vertcat(measured_data.dty_vpk)*1.60934;
            measured_flw_vph = vertcat(measured_data.flw_out_vph);
            measured_health = error_data_provider.get_health(this.day,eval_meas_ids);
            clear measured_data eval_meas_ids
            
            % remove bad detectors
            eval_links = eval_links(measured_health);
            num_sensors = numel(eval_links);
            measured_dty_vpm(~measured_health,:) = [];
            measured_flw_vph(~measured_health,:) = [];
            clear measured_health
            
            % organize into matrices per prediction time interval
            prediction_data = repmat(struct('horizon',nan, ...
                'time', [] , ...
                'flw_vph', nan(num_sensors,this.num_predictions) , ...
                'dty_vpm', nan(num_sensors,this.num_predictions) ) , ...
                1,this.num_pred_meas);
            
            prediction_times = [this.log.start_time];
            for j=1:this.num_pred_meas
                prediction_data(j).horizon = j*this.measurement_dt;
                prediction_data(j).time = prediction_times + j*this.measurement_dt;
            end
            
            link_ind = Utils.index_into(eval_links,this.ni.link_ids);
            for i=1:this.num_predictions
                for j=1:this.num_pred_meas
                    % prediction of state at time i, made at time i-j.
                    prediction_data(j).flw_vph(:,i) = this.log(i).predicted_state.flw_vph(j,link_ind)';
                    prediction_data(j).dty_vpm(:,i) = this.log(i).predicted_state.dty_vpm(j+1,link_ind)';
                end
            end
            
            % compare each prediction to data
            this.prediction_error.flw_vph = nan(1,this.num_pred_meas);
            this.prediction_error.dty_vpm = nan(1,this.num_pred_meas);
            this.prediction_error.spd_mph = nan(1,this.num_pred_meas);
            for j=1:this.num_pred_meas
                
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
                
                this.prediction_error.flw_vph(j) = sqrt(Utils.meanwithnan(d_flw.^2));
                this.prediction_error.dty_vpm(j) = sqrt(Utils.meanwithnan(d_dty.^2));
                this.prediction_error.spd_mph(j) = sqrt(Utils.meanwithnan(d_spd.^2));
                
            end
            
            if nargin>=3
                figure('Position',[165 163 866 503])
                pred_time = this.measurement_dt:this.measurement_dt:this.prediction_horizon;
                xlab = 'horizon length [sec]';
                
                ylab = 'error flw [vph]';
                subplotf(3,1,ylab,xlab,pred_time,this.prediction_error.flw_vph,'','b','',2);
                set(gca,'XLim',[0 this.prediction_horizon])
                grid
                
                ylab = 'error dty [vpm]';
                subplotf(3,2,ylab,xlab,pred_time,this.prediction_error.dty_vpm,'','b','',2);
                set(gca,'XLim',[0 this.prediction_horizon])
                set(gca,'YLim',[0 max([1 1.1*max(this.prediction_error.dty_vpm)])])
                grid
                
                ylab = 'error spd [mph]';
                subplotf(3,3,ylab,xlab,pred_time,this.prediction_error.spd_mph,'','b','',2);
                set(gca,'XLim',[0 this.prediction_horizon])
                grid
                
                [ppt,op]=openppt(pptfile,true);
                addslide(op,'Prediction error')
                closeppt(ppt,op)
                
                fprintf('Average flow prediction error: %.1f vph\n',Utils.meanwithnan(this.prediction_error.flw_vph'))
                fprintf('Average density prediction error: %.1f vpm\n',Utils.meanwithnan(this.prediction_error.dty_vpm'))
                fprintf('Average speed prediction error: %.1f mph\n',Utils.meanwithnan(this.prediction_error.spd_mph'))
                
            end
        end
        
        function [error_data_provider] = get_error_data_provider(this,type)
            config = Config.get(this.configname);
            switch type
                case 'pems'
					error_data_provider = Utils.get_pems_dp(this.ni,this.configname);
                case 'simulation'
                    mr = ObjectFactory.beats_model_runner( ...
                                    this.ni , ...
                                    this.sim_dt , ...
                                    this.measurement_dt );
                    error_data_provider = ObjectFactory.sim_data_provider(mr,this.end_time);
            end
        end
        
    end
    
    methods(Access=protected)

        function [this]=set_values(this,params)

            % store inputs
            this.configname = params.configname;
            this.start_time = params.start_time;
            this.end_time = params.end_time;
            this.update_dt = params.update_dt;
            this.prediction_horizon = params.prediction_horizon;
            this.measurement_dt = params.measurement_dt;
            this.num_pred_meas = this.prediction_horizon / this.measurement_dt;
            this.num_predictions = ceil( (this.end_time-this.start_time)/this.update_dt );
            this.log = repmat(struct('start_time',nan, ...
                                     'initial_state',nan, ...
                                     'predicted_demands',nan, ...
                                     'predicted_splits',nan, ...
                                     'predicted_state',nan ),1,this.num_predictions);        
        end        
        
        function [] = check_values(this)
            
            % non-negativity
            if any([ this.start_time ...
                    this.end_time ...
                    this.update_dt ...
                    this.measurement_dt ...
                    this.prediction_horizon] <0 )
                error('violated non-negativity')
            end
            
            % non-zero
            if any([this.update_dt this.prediction_horizon]==0)
                error('violated non-zero')
            end
            
            % range checks
            if this.end_time<=this.start_time
                error('this.end_time<=this.start_time')
            end
            
            if this.prediction_horizon<=this.update_dt
                error('this.prediction_horizon<=this.update_dt')
            end
            
            % mod checks
            if mod(this.end_time-this.start_time,this.update_dt)~=0
                error('mod(this.end_time-this.start_time,this.update_dt)~=0')
            end
            
        end
        
        function [this]=instantiate_objects(this,params)
            
            config = Config.get(params.configname);

            % instantiate network information
            this.ni = ObjectFactory.network_information(config.xml_file);

        end
        
    end
    
    methods(Abstract)
        [this]=run(this);
    end
    
    methods(Access=protected,Static)
        
        function []=check_inputs(params)
            
            if ~isfield(params,'configname')
                error('~isfield(params,''configname'')')
            end
            
            if ~isfield(params,'update_dt')
                error('~isfield(params,''update_dt'')')
            end
            
            if ~isfield(params,'prediction_horizon')
                error('~isfield(params,''prediction_horizon'')')
            end
            
            if ~isfield(params,'measurement_dt')
                error('~isfield(params,''measurement_dt'')')
            end
            
            if ~isfield(params,'start_time')
                error('~isfield(params,''start_time'')')
            end
            
            if ~isfield(params,'end_time')
                error('~isfield(params,''end_time'')')
            end
            
        end
        
    end

end

