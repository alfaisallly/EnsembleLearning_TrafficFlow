classdef (Abstract) FlowPredictor < handle

    properties
        data_provider @DataProvider
    end
    
    methods(Access=public)
        
        function [this]=FlowPredictor(data_provider)
            if nargin>0
                this.data_provider = data_provider;
            end
        end
        
        function [b] = check_have_day(this,day)
            b = this.data_provider.check_have_day(day);
        end
        
        function [this] = set_force_bad_sensors(this,day,ids)
            this.data_provider.set_force_bad_sensors(day,ids);
        end
        
        function [flw_vph,time] = predict_flw_for_dpids(this,dpids,day,from,to,dt)
            
            p = this.predict(dpids,day,from,to,dt);
            if isempty(p)
                flw_vph = [];
                time = [];
                return
            else
                flw_vph = nan(length(dpids),length(p(1).flw_out_vph));
                time = p(1).time;
            end
            
            for i=1:length(dpids)
                % use incoming flow if available. 
                % The assumption is that if it is not available it is
                % because this is a sensor, not a link.
                if ~isempty(p(i).flw_in_vph)
                    flw_vph(i,:) = p(i).flw_in_vph;
                else
                    flw_vph(i,:) = p(i).flw_out_vph;
                end
            end
        end

        function [historical_uptonow,historical_later,today_flw,historical_uptonow_cluster,historical_later_cluster] = gather_data_for_prediction(this,day,id,from,dt,to,past_window,lag)
            
            if nargin<6
                to = 86400;
            end
            
            if nargin<7 || isnan(past_window)
                past_window = inf;
            end
            
            if nargin<8
                lag = 0;
            end
            
            if from-past_window-lag<0
                past_window = from-lag;
            end
            
            % get data matrix
            x = this.data_provider.get_data('all',id,0,86400,dt,'interpolate');
            flows = vertcat(x.flw_out_vph);
            clear x
            
            % keep "day" data
            %day_ind = (this.data_provider.days==day )'; %test day is single day
            day_ind = ismember(this.data_provider.days,day )';
            
            if ~any(day_ind)
                error('~any(day_ind)')
            end            

            start_time = from-past_window-lag;
            start_ind = max([1 floor(start_time/dt)]);
            from_lag_ind = max([1 floor((from-lag)/dt)]);
            from_ind = max([1 floor(from/dt)]);
            to_ind = max([1 floor(to/dt)]);
            to_ind = min([to_ind size(flows,2)]);

            % construct historical 
            historical_uptonow = flows(~day_ind , start_ind:from_lag_ind);
            historical_later   = flows(~day_ind , from_ind+1:to_ind);

            % construct today data
            today_flw = flows(day_ind,start_ind:from_lag_ind);
                   
            % get clustered historical data
            historical_uptonow_cluster = this.data_provider.get_representative_data(...
                struct( 'day_of_week',weekday(datenum(day))), ...
                        id, ...
                        from-lag-past_window, ...
                        from-lag , ...
                        dt);
            
            historical_later_cluster = this.data_provider.get_representative_data(...
                struct('day_of_week',weekday(datenum(day))), ...
                    id, ...
                    from - lag , ...
                    to , ...
                    dt);
   
        end
        
           
    end
    
    methods(Abstract)
        [x] = predict(this,ids,days,from,to,dt)
    end
    
    methods(Static)
        
        % (this,ids,update_dt,horizon,sim_data_provider,)
        function [prediction] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,flow_predictor,update_dt,horizon,sim_data_provider}
            
            ppt_file = fullfile(Folder.reports,params.ppt_file);
            ppt_file_error = sprintf('%s_error',ppt_file);
            flow_predictor = params.flow_predictor;
            dp_ids = flow_predictor.data_provider.ids;                              % internal ids to be reported
            link_ids = flow_predictor.data_provider.get_linkids_for_ids(dp_ids);    % corresponding link ids

            dt = flow_predictor.data_provider.measurement_dt;
            update_times = 0:params.update_dt:(86400-params.update_dt);   
            n_horizon = params.horizon/dt;
            
            % allocate prediction array
            prediction = repmat(struct('link_id',nan,'dp_id',nan,'flw',nan(length(update_times),n_horizon)),1,length(dp_ids));
            for i=1:length(dp_ids)
                prediction(i).dp_id = dp_ids(i);
                prediction(i).link_id = link_ids(i);
            end
            
            % collect true flows from simulation data provider
            fprintf('\tCollecting actual flows\n')
            actual = params.sim_data_provider.get_data('all',link_ids,'start','end');
            
            % run predictions
            fprintf('\tRunning %d predictions',length(update_times))
            for k=1:length(update_times)
                fprintf('.')
                start_time = update_times(k);
                x = flow_predictor.predict(dp_ids, ...
                    params.day , ...
                    start_time, ...
                    start_time+params.horizon, ...
                    dt);
                
                lasttime = min([actual(1).time(end)-dt start_time+params.horizon-dt]);
                ind = Utils.index_into(start_time:dt:lasttime,actual(1).time);
                for i=1:length(dp_ids)
                    err = Utils.error(actual(i).flw_out_vph(ind),x(i).flw_out_vph(1:length(ind)));
                    prediction(i).error.mae(k) = err.mae;
                    prediction(i).error.mape(k) = err.mape;
                    prediction(i).error.smape(k) = err.smape;
                    prediction(i).error.rmse(k) = err.rmse;
                    prediction(i).error.linf(k) = err.linf;
                    prediction(i).flw(k,:) = x(i).flw_out_vph;
                end
            end
            fprintf('\n')
                         
            % error table : NOTE: THIS IS AN INCORRECT AVERAGE, CHANGE LATER
            for i=1:length(dp_ids)
                prediction(i).avg_error = Utils.avg_error(prediction(i).error);
            end
            
            if isfield(params,'xls_file')
                X = [prediction.avg_error];
                write(table([prediction.dp_id nan]', ...
                            [[X.mae]';Utils.meanwithnan([X.mae]')], ...
                            [[X.mape]';Utils.meanwithnan([X.mape]')], ...
                            [[X.smape]';Utils.meanwithnan([X.smape]')], ...
                            [[X.rmse]';Utils.meanwithnan([X.rmse]')], ...
                            [[X.linf]';Utils.meanwithnan([X.linf]')], ...
                    'VariableNames',{'id','mae','mape','smape','rmse','linf'}) , ...
                    fullfile(Folder.reports,sprintf('%s.xls',params.xls_file)) );
                clear X
            end

            % collect historical flows from data provider
            fprintf('\tCollecting historical flows\n')
            historical = flow_predictor.data_provider.get_data('all',dp_ids,'start','end');
            
            % loop health from data_provider [vds x days]
            all_is_good = flow_predictor.data_provider.get_health('all',dp_ids);

            % plot predictions for each id
            do_export = ~isempty(params.ppt_file);
            if do_export
                [ppt,op]=openppt(ppt_file,true);
                [ppte,ope]=openppt(ppt_file_error,true);
                
                % title slide
                addslideText(op,'Flow prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d sec.\nPrediction horizon: %d sec.\n', ...
                    class(flow_predictor),class(flow_predictor.data_provider),params.update_dt,params.horizon) )
                
                % title slide
                addslideText(ope,'Flow prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d\nPrediction horizon: %d\n', ...
                    class(flow_predictor),class(flow_predictor.data_provider),params.update_dt,params.horizon) )
            end
            
            for i=1:length(dp_ids)
                
                fprintf('\tSlide %d of %d\n',i,length(dp_ids))
                
                if do_export
                    figure('Position',[279    95   710   571],'Visible','off')
                else
                    figure('Position',[279    95   710   571],'Visible','on')
                end
                
                xlab = 'time [hr]';
                
                % plot good historical data
                vds_is_good = all_is_good(:,i);
                if any(vds_is_good)
                    historical_time = historical(1,1).time(1:end-1)/3600;
                    historical_flw = vertcat(historical(:,i).flw_out_vph);
                    plot(historical_time,historical_flw(vds_is_good,:),'k');
                end
                
                % actual
                hold on
                plot(actual(i).time(1:end-1)/3600,actual(i).flw_out_vph,'r','LineWidth',4)
                
                % prediction
                P = prediction(i);
                for k=1:length(update_times)
                    time = (update_times(k):dt:update_times(k)+params.horizon-dt)/3600;
                    plot(time,P.flw(k,:),'Color',rand(1,3),'LineWidth',2);
                end
                
                xlabel(xlab);
                ylabel('flow [vph]');
                set(gca,'XLim',[0 24])
                grid

                percent_good = 100*sum(vds_is_good)/numel(vds_is_good);
                if do_export
                    addslide(op,sprintf('(good days:%.0f%%) id=%d',percent_good,dp_ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                else
                    title(sprintf('Flow predictions, id=%d (good days:%.0f%%)',dp_ids(i),percent_good))
                end
                
                                
                if do_export
                    figure('Position',[279    95   710   571],'Visible','off')
                else
                    figure('Position',[279    95   710   571],'Visible','on')
                end
                
                subplotf(2,1,'error [vph]',xlab,update_times/3600,[ ...
                                prediction(i).error.linf ; ... 
                                prediction(i).error.rmse ; ...
                                prediction(i).error.mae ],'','','',2);
                set(gca,'XLim',[0 24])
                legend(sprintf('linf %.0f',prediction(i).avg_error.linf), ...
                       sprintf('rmse %.0f',prediction(i).avg_error.rmse), ...
                       sprintf('mae %.0f',prediction(i).avg_error.mae),'Location','Best')
                grid
                
                subplotf(2,2,'error [%]',xlab,update_times/3600,100*[ ...
                                prediction(i).error.mape ; ...
                                prediction(i).error.smape  ],'','','',2);
                set(gca,'XLim',[0 24])
                legend( sprintf('mape %.1f',prediction(i).avg_error.mape), ...
                        sprintf('smape %.1f',prediction(i).avg_error.smape), 'Location','Best')
                grid
                
                if do_export
                    addslide(ope,sprintf('(good days:%.0f%%) id=%d',percent_good,dp_ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                end
                
            end
            
            if do_export
                closeppt(ppt,op)
                closeppt(ppte,ope)
            end
            
        end
    end

end

