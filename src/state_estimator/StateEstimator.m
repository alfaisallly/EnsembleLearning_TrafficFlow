classdef (Abstract) StateEstimator < handle
    
    properties(Access=public)
        model_runner@ModelRunner    % A ModelRunner
        data_provider@DataProvider  % A DataProvider
        day                         % day I am estimating
        link_ids                    % list of links in the network that I am estimating
        xhat_time                   % time stamp for xhat    
        xhat                        % [veh] state on all links in link_ids
        accessible_dp_ids           % list of data provider indices to which I have access
        accessible_link_ind         % link_ids indices corresponding to accessible_dp_ids
        time_init                   % initial time for reset
        xhat_init                   % initial state for reset
    end
    
    methods(Access=public)
        
        function [this] = StateEstimator(model_runner,data_provider,link_ids,use_dp_ids,time_init,xhat_init)
            
            if Utils.isclass(model_runner,'ModelRunner')
                this.model_runner = model_runner.clone('StateEstimator.Constructor');
            end
            this.data_provider = data_provider;
            this.link_ids = Utils.row_vector(link_ids);

            % decide what data can be used
            if nargin<4
                this.accessible_dp_ids = this.data_provider.ids;
            else
                if ischar(use_dp_ids) && strcmp(use_dp_ids,'all')
                    this.accessible_dp_ids = this.data_provider.ids;
                else
                    if ~all(ismember(use_dp_ids, this.data_provider.ids))
                        error('~all(ismember(use_dp_ids, this.data_provider.ids))')
                    end
                    this.accessible_dp_ids = use_dp_ids;
                end
            end
            
            accessible_link_ids = this.data_provider.get_linkids_for_ids(this.accessible_dp_ids);
            this.accessible_link_ind = Utils.index_into(accessible_link_ids,this.link_ids);
            
            if nargin<5
                this.time_init = -inf;
                this.xhat_init = zeros(1,length(this.link_ids));
            else
                if length(xhat_init)~=length(this.link_ids)
                    error('length(xhat_init)~=length(this.link_ids)')
                end
                this.time_init = time_init;
                this.xhat_init = Utils.row_vector(xhat_init);
            end
            
            this.reset();
            
        end
        
        function [this] = set_day(this,this_day)
            
            if ~all(isnan(this.data_provider.days)) && ~ismember(this_day,this.data_provider.days)
                error('~all(isnan(this.data_provider.days)) && ~ismember(this_day,this.data_provider.days)')
            end
            
            this.day = this_day;
            this.reset;
        end

        function [this] = set_force_bad_sensors(this,this_day,bad_vdss)
            this.data_provider.set_force_bad_sensors(this_day,bad_vdss);
            this.set_day(this_day);
        end
        
        function [this] = set_initial_state_veh(this,x_time_init,x_link_ids,x_xhat_init)
            this.time_init = x_time_init;
            this.xhat_init = zeros(1,length(this.link_ids));
            ind = Utils.index_into(x_link_ids,this.link_ids);
            this.xhat_init(ind(ind>0)) = x_xhat_init(ind>0);
        end
            
        function [x] = get_state_for_link_ids(this,lnk_ids)
        % get the current state estimate for a list of link ids
        % (lnk_ids)
        % lnk_ids : list of link ids
        % returns: row vector with link densities in XXX
            x = zeros(1,numel(lnk_ids));
            ind = Utils.index_into(lnk_ids,this.link_ids);
            haveit = ind>0;
            x(haveit) = this.xhat(ind(haveit));
        end

        function [this] = set_accessible_dp_ids(this,x)
        % set the data provider ids to which this state estimator has access
        % (x)
        % x : list of data provider ids
        % returns: nothing
        
            if ~all(ismember(x,this.data_provider.ids))
                error('~all(ismember(x,this.data_provider.ids))')
            end
            this.accessible_dp_ids = x;
            accessible_link_ids = this.data_provider.get_linkids_for_ids(this.accessible_dp_ids);
            this.accessible_link_ind = Utils.index_into(accessible_link_ids,this.link_ids);
        end
        
        function [this] = reset(this)
        % reset the estimator. This resets the mode runner and the state
        % estimate to its original value
        
            this.xhat_time = this.time_init;
            this.xhat = this.xhat_init;
%             if ~isempty(this.model_runner)
%                 this.model_runner.reset();
%             end
        end
       
%         function [internal_dp_ids,internal_link_ids]=get_internal_sensors(this,ni)
%             % get all data provider ids that are internal to the network
%             internal_link_ids = ni.get_internal_link_ids;
%             dp_ids = this.data_provider.get_ids_for_linkids(internal_link_ids);
%             has_dp = cellfun(@(x) ~isempty(x),dp_ids);
%             internal_dp_ids = Utils.cell2array(dp_ids(has_dp));     % internal sensors
%             internal_link_ids = internal_link_ids(has_dp);          % corresponding links                     
%         end
        
    end
    
    methods(Abstract)
        update_estimate(this,new_time)
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)

            dp = params.state_estimator.data_provider;
            
            % get all data provider ids that are internal to the network
            [internal_dp_ids,internal_link_ids] = dp.get_internal_sensors(params.network_info);          
            
            % misc
            do_export = ~isempty(params.ppt_file);
            if do_export
                vis = 'off';
            else
                vis = 'on';
            end
            update_times = params.update_dt:params.update_dt:86400; 
            I = length(internal_dp_ids);
            K = length(update_times);
            
            % allocate for saving estimation results
            X.update_times = update_times;
            X.dp_ids = internal_dp_ids;
            X.dty_estim_veh = nan(I,K,I);        % space X time X removed dp
            
            % for all dp ids
            for i=1:I
                
                % remove one
                remove_dp_id = internal_dp_ids(i);
                
                % reset the state estimator
                params.state_estimator.reset();
                
                % remove remove_vds from the estimator's list of available data sources
                params.state_estimator.set_accessible_dp_ids( setdiff(internal_dp_ids,remove_dp_id,'stable') );

                % run estimations for one day
                fprintf('\testimating without data provider id %d (%d of %d)',remove_dp_id,i,I)
                
                for k=1:K

                    % update estimator
                    params.state_estimator.update_estimate(update_times(k));

                    % record estimated state
                    X.dty_estim_veh(:,k,i) = params.state_estimator.get_state_for_link_ids(internal_link_ids);
                end
                fprintf('\n')
            end
              
            % collect measurements
            X.dty_meas_veh = nan(I,K);
            for k=1:K
                X.dty_meas_veh(:,k) = dp.extract_density_veh( ...
                    dp.get_data('all',internal_dp_ids,update_times(k)) );
            end
            
            % compute errors
            X.error = nan(I,K); % error(removed detector,time)
            for i=1:I
                X.error(i,:) = sum(abs(X.dty_meas_veh - X.dty_estim_veh(:,:,i)),1);
            end
            
            % write to powerpoint -----------------------------------------
            if do_export
                pptfile = fullfile(Folder.reports,params.ppt_file);
                fprintf('\n\tExporting to %s\n',pptfile)
                [ppt,op]=openppt(pptfile,true);
            
                % title slide
                addslideText(op,'State estimation report',...
                    sprintf('Estimator: %s\nConfig: %s\nUpdate dt: %d sec.\n', ...
                    class(params.state_estimator),params.network_info.configfile,params.update_dt) )  
            end
            
            % plot measured state
            figure('Position',[142   236   847   562],'Visible',vis);
            Utils.plot_24hr_contour(X.dty_meas_veh,'dty',internal_dp_ids)
            if do_export
                addslide(op,'Measured density contour','new',[0.5 0.5],'',0.6)
                close
            else
                title('Measured density contour')
            end
            
            % plot errors
            avg_error = mean(X.error,1);
            avg_avg_error = mean(avg_error);
            figure('Position',[249 45 851 420],'Visible',vis);
            
            Utils.plot_24_hour_line(avg_error)
            ylabel('Average one-removed error [veh]')
            tit = sprintf('One-removed errors (avg = %.1f veh)',avg_avg_error);
            if do_export
                addslide(op,tit,'new',[0.5 0.5],'',0.6)
                close
            else
                title(tit)
            end
            
            % plot estimated states
            fprintf('\t')
            for i=1:I
                
                fprintf('.')
                remove_dp_id = internal_dp_ids(i);
       
                figure('Position',[142   236   847   562],'Visible',vis);
                Utils.plot_24hr_contour(X.dty_estim_veh(:,:,i),'dty',internal_dp_ids)
                
                if do_export
                    addslide(op,sprintf('removed dp id=%d',remove_dp_id),'new',[0.5 0.5],'',0.6)
                    close
                else
                    title(sprintf('removed dp id=%d',remove_dp_id))
                end
                
            end
            fprintf('\n')
            
            if do_export
                closeppt(ppt,op)
            end
            
        end
        
    end
    
end
