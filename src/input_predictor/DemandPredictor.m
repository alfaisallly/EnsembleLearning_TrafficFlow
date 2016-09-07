classdef DemandPredictor < InputPredictor

    properties
        dp_ids      % cell array, data provider internal ids for source links
                    % flow is added for multiple dp ids on a single source link
    end
    
    methods
        
        function [this]=DemandPredictor(network_information,flow_predictor)
            this = this@InputPredictor(flow_predictor);
            this.ids = network_information.get_source_link_ids;
            this.dp_ids = this.flow_predictor.data_provider.get_ids_for_linkids(this.ids);
            has_dp = cellfun(@(x)~isempty(x),this.dp_ids);
            if ~all(has_dp)
                error('~all(has_dp), Try running ''generate_fake_pems_data''')
            end
        end
        
        function [] = set_day(this,day)
            if ~this.check_have_day(day)
                error('~this.check_have_day(day)')
            end
        end
        
        function [x] = predict(this,day,from,to,dt)
            if nargin<5
                dt = nan;
            end
            
            % iterate through all link ids
            x = repmat(struct('link_id',nan,'dp_id',nan,'time',[],'flw_vph',[]),1,length(this.ids));
            for i=1:length(this.ids)
			
				[flw_vph,time] = this.flow_predictor.predict_flw_for_dpids(this.dp_ids{i},day,from,to,dt);
				flw_vph = sum(flw_vph,1);
				
                x(i) = struct( ...
                    'link_id' ,this.ids(i), ...
                    'dp_id'   ,this.dp_ids{i}, ...
                    'time'    ,time, ...
                    'flw_vph' ,flw_vph);
            end
        end
        
        function [x] = predict_given_flows(this,flow_pred)
            
            % iterate through all link ids
            x = repmat(struct('link_id',nan,'dp_id',nan,'time',[],'flw_vph',[]),1,length(this.ids));
            for i=1:length(this.ids)
                ind = ismember([flow_pred.dp_ids],[this.dp_ids{i}]);
                flw_vph = flow_pred(ind).flow;
                flw_vph = sum(flw_vph,1);
                time=flow_pred(ind).time;
                
                x(i) = struct( ...
                    'link_id' ,this.ids(i), ...
                    'dp_id'   ,this.dp_ids{i}, ...
                    'time'    ,time, ...
                    'flw_vph' ,flw_vph);
            end
        end
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,flow_predictor_type,config,output_dt,end_time,update_dt,horizon}
               
            X = [];
            
            config = Config.get(params.config);
            ni = ObjectFactory.network_information(config.xml_file);
            mr = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            sim_dp = ObjectFactory.sim_data_provider(mr,params.end_time);
            day = config.model_day;
            
            % create the flow predictor
            switch params.flow_predictor_type
                case 'simulation'
                    fp = ObjectFactory.simulation_flow_predictor(sim_dp);
                case 'zoh'
                    fp = ObjectFactory.zoh_predictor(sim_dp);
                case 'historical',
					pems_dp = Utils.get_pems_dp(ni,params.config);
                    fp = ObjectFactory.historical_predictor(pems_dp);
                case 'scaled_historical'
					pems_dp = Utils.get_pems_dp(ni,params.config);
                    fp = ObjectFactory.scaled_historical_predictor(sim_dp,pems_dp);
                case 'recursiveARMAX'
                    pems_dp= Utils.get_pems_dp(ni,params.config);
                    fp = ObjectFactory.recursiveARMAX_predictor(pems_dp);
                otherwise
                    error('bad flow_predictor_type')
            end
            
            % create the demand predictor
            dem_pred = DemandPredictor(ni,fp);
            
            % reporting ---------------------------------------------------
            
            ppt_file = fullfile(Folder.reports,params.ppt_file);
            ppt_file_error = sprintf('%s_error',ppt_file);

            dt = dem_pred.flow_predictor.data_provider.measurement_dt;
            update_times = 0:params.update_dt:(86400-params.update_dt);   
            n_horizon = params.horizon/dt;
            dp_ids_set = Utils.flatten_cell_array(dem_pred.dp_ids);
            
            % allocate prediction array
            prediction = repmat(struct('link_id',nan,'dp_ids',nan,'flw',nan(length(update_times),n_horizon)),1,length(dem_pred.ids));
            for i=1:length(dem_pred.ids)
                prediction(i).link_id = dem_pred.ids(i);
                prediction(i).dp_ids = dem_pred.dp_ids{i};
            end
            
            % collect true flows from simulation data provider
            fprintf('\tCollecting actual demands\n')
            actual = sim_dp.get_data('all',dem_pred.ids,'start','end');
            
            % run predictions
            fprintf('\tRunning %d predictions',length(update_times))
            for k=1:length(update_times)

                fprintf('.')
                start_time = update_times(k);
                x = dem_pred.predict( day , ...
                                      start_time, ...
                                      start_time+params.horizon, ...
                                      dt);
                
                lasttime = min([actual(1).time(end)-dt start_time+params.horizon-dt]);
                ind = Utils.index_into(start_time:dt:lasttime,actual(1).time);

                for i=1:length(dem_pred.ids)
                    err = Utils.error(actual(i).flw_in_vph(ind),x(i).flw_vph(1:length(ind)));
                    prediction(i).error.mae(k) = err.mae;
                    prediction(i).error.mape(k) = err.mape;
                    prediction(i).error.smape(k) = err.smape;
                    prediction(i).error.rmse(k) = err.rmse;
                    prediction(i).error.linf(k) = err.linf;
                    prediction(i).flw(k,:) = x(i).flw_vph;
                end

            end
            fprintf('\n')
            
            % error table : NOTE: THIS IS AN INCORRECT AVERAGE, CHANGE LATER
            for i=1:length(dem_pred.ids)
                prediction(i).avg_error = Utils.avg_error(prediction(i).error);
            end
            
            if isfield(params,'xls_file')
                X = [prediction.avg_error];
                write(table([prediction.link_id nan]', ...
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
            historical = dem_pred.flow_predictor.data_provider.get_data('all',dp_ids_set,'start','end');
            
            % loop health from data_provider [vds x days]
            all_is_good = dem_pred.flow_predictor.data_provider.get_health('all',dp_ids_set);

            % plot predictions for each id
            do_export = ~isempty(params.ppt_file);
            if do_export
                [ppt,op]=openppt(ppt_file,true);
                [ppte,ope]=openppt(ppt_file_error,true);
                
                % title slide
                addslideText(op,'Demand prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d sec.\nPrediction horizon: %d sec.\n', ...
                    params.flow_predictor_type,class(dem_pred.flow_predictor.data_provider),params.update_dt,params.horizon) )
                
                % title slide
                addslideText(ope,'Demand prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d\nPrediction horizon: %d\n', ...
                    params.flow_predictor_type,class(dem_pred.flow_predictor.data_provider),params.update_dt,params.horizon) )
            end
            
            for i=1:length(dem_pred.ids)
                
                fprintf('\tSlide %d of %d\n',i,length(dem_pred.ids))
                
                if do_export
                    figure('Position',[279    95   710   571],'Visible','off')
                else
                    figure('Position',[279    95   710   571],'Visible','on')
                end
                
                dp_set_ind = Utils.index_into(dem_pred.dp_ids{i},dp_ids_set);
                
                xlab = 'time [hr]';
                
                % plot good historical data
                vds_is_good = all_is_good(:,dp_set_ind);
                if any(vds_is_good)
                    historical_time = historical(1,1).time(1:end-1)/3600;
                    
                    % WARNING: THIS WILL PROBABLY FAIL FOR MULTIPLE DP IDS
                    % ON A RAMP
                    historical_flw = vertcat(historical(:,dp_set_ind).flw_out_vph);
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
                    addslide(op,sprintf('(good days:%.0f%%) id=%d',percent_good,dem_pred.ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                else
                    title(sprintf('Demand predictions, id=%d (good days:%.0f%%)',dem_pred.ids(i),percent_good))
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
                    addslide(ope,sprintf('(good days:%.0f%%) id=%d',percent_good,dem_pred.ids(i)),'new',[0.5 0.5],'',0.6)
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

