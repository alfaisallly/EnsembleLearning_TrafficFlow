classdef CompleteSensingJobManager < PredictionJobManager
    
    properties
        good_ml_vds_for_day
        num_good_ml_vds_for_day
    end
    
    methods( Access = public )
        
        function [this] = CompleteSensingJobManager(confname,infolder,outfolder)
            this = this@PredictionJobManager(confname,infolder,outfolder);
        end
        
        function [this] = plot_sample_distribution(this)
            figure
            imagesc(this.num_samples_grid)
            xlabel('num decimated ml')
            ylabel('num decimated rp')
            title(sprintf('Total samples = %d',sum(sum(this.num_samples_grid))))
        end
        
        function [result] = process_result(this,results_files)
            
            % Create a data provider for calculating error
            error_data_provider = Utils.get_pems_dp(this.ni,this.configname);
            
            % initialize structures
            num_results = length(results_files);
            good_vdss = repmat(struct('ids',[],'link_ids',[]),num_results,1);
            
            five_minute_prediction = repmat(struct('flw_vph',[]),num_results,1);
            thirty_minute_prediction = repmat(struct('flw_vph',[]),num_results,1);
            sixty_minute_prediction = repmat(struct('flw_vph',[],'flw_error',[]),num_results,1);
            
            estimation_vpm = cell(1,num_results);
            demand_prediction = repmat(struct('demand',[]),num_results,1);
            sr_prediction = repmat(struct('split_ratio',[]),num_results,1);
            bf_prediction = repmat(struct('boundary_flow',[]),num_results,1);
            
            % loop through results
            for i=1:num_results
                
                r = load(fullfile(this.outputfolder,results_files{i}));
                
                % process prediction
                [X,~,Y,vdss,D,SR,BF] = this.process_log_prediction(error_data_provider,r.X.day,r.X.log);
                
                good_vdss(i) = vdss;
                
                ind = [X.horizon]==5*60;
                five_minute_prediction(i).flw_vph = X(ind).flw_vph;
                five_minute_prediction(i).dty_vpm = X(ind).dty_vpm;
                
                ind = [X.horizon]==30*60;
                thirty_minute_prediction(i).flw_vph = X(ind).flw_vph;
                thirty_minute_prediction(i).dty_vpm = X(ind).dty_vpm;
                
                ind = [X.horizon]==60*60;
                sixty_minute_prediction(i).flw_vph = X(ind).flw_vph;
                sixty_minute_prediction(i).dty_vpm = X(ind).dty_vpm;
                sixty_minute_prediction(i).flw_error = Y;
                
                demand_prediction(i).demand=D;
                sr_prediction(i).split_ratio=SR;
                bf_prediction(i).boundary_flow=BF;
                
                % process estimation
                if isfield(r.X.log,'initial_state')
                    estimation_vpm{i} = this.process_log_estimation(r.X.log);
                end
                
            end
            
            % save to result
            result = struct( ...
                'good_vdss', good_vdss , ...
                'five_minute_prediction' , five_minute_prediction , ...
                'thirty_minute_prediction' , thirty_minute_prediction , ...
                'sixty_minute_prediction' , sixty_minute_prediction , ...
                'estimation_vpm' , {estimation_vpm}, ...
                'predicted_demand',demand_prediction, ...
                'predicted_split_ratio',sr_prediction, ...
                'predicted_boundary_flow', bf_prediction );
        end
        
        function [this] = report_result(this,result,output_file)
            
            if nargin<3
                output_file = '';
            end
            
            % resport estimation
            this.report_estimation(result,output_file);
%             
%             % report demand prediction
%             this.report_demand_prediction(result,output_file);
%             
%             % report split ratio prediction
%             this.report_sr_prediction(result,output_file);
            
            % report one hour flow error
            this.report_one_hour_flow_error(result,output_file);
            
            % report predicted contours
            this.report_prediction_contours(result,output_file);
            
        end
    end
    
    methods( Access = protected )
        
        function [this] = initialize(this)
            
            % find removeable internal vdss for each day
            % these are the good mainline vdss on each day
            dp_pems = Utils.get_pems_dp(this.ni,this.configname);
            
            fprintf('Getting good mainline sensors for each day\n');
            for d=1:length(this.days)
                this.good_ml_vds_for_day{d} = dp_pems.get_good_internal_sensors(this.ni,this.days(d));
            end
            this.num_good_ml_vds_for_day = cellfun(@(x) numel(x),this.good_ml_vds_for_day);
            
        end
        
        function [this] = build_job(this,~)

            all_tasks = repmat(struct('day',nan,'remove_ml_vdss',[],'remove_rp_vdss',[]),length(this.days),1);

            for i=1:length(this.days)
                all_tasks(i).day = this.days(i);
            end
            
            % write to file
            this.tasks_file = this.save_job(all_tasks);
        end
        
    end
    
    methods( Access = private )
        
        function [X] = uniform_sampling(this,num_ml_remove,num_rp_remove)
            
            % pick a day with sufficient removable ml vdss
            eligible_days = this.num_good_ml_vds_for_day > num_ml_remove;
            day = datasample(this.days(eligible_days),1);
            today_removable_ml_vds = this.good_ml_vds_for_day{this.days==day};
            
            % select at random the vdss to remove
            if num_ml_remove>0
                decimate_ml = datasample(today_removable_ml_vds,num_ml_remove,'Replace',false);
            else
                decimate_ml = [];
            end
            
            if num_rp_remove>0
                decimate_rp = datasample(this.removable_rp_vds,num_rp_remove,'Replace',false);
            else
                decimate_rp = [];
            end
            
            remove_vdss = [decimate_ml decimate_rp];
            
            % record number of mainline and ramp decimated vdss
            this.record_num_samples(remove_vdss);
            
            % save in struct
%             X = struct('day',day,'remove_vdss',remove_vdss);
            X = struct('day',day,'remove_ml_vdss',decimate_ml,'remove_rp_vdss',decimate_rp);            
            
        end
        
        function [X] = binomial_sampling(this,prob_fail)
            
            % pick a day with sufficient removable ml vdss
            day = datasample(this.days,1);
            today_removable_ml_vds = this.good_ml_vds_for_day{this.days==day};
            
            % all removable vdss
            removable_vdss = [today_removable_ml_vds this.removable_rp_vds];
            
            % select a number of vdss to remove
            num_decimate = binornd(numel(removable_vdss),prob_fail);
            
            % select the vdss to remove, make sure we are not removing all of the mainline sensors.
            not_done = true;
            while not_done
                
                % choose num_decimate from vdss
                remove_vdss = datasample(removable_vdss,num_decimate,'Replace',false);
                
                % finish if you have not selected all removable mainline vdss
                not_done = isempty(setdiff(today_removable_ml_vds,remove_vdss));
                
            end
            
            % record number of mainline and ramp decimated vdss
            this.record_num_samples(remove_vdss);
            
            % save in struct
%             X = struct('day',day,'remove_vdss',remove_vdss);
            decimate_ml = intersect(remove_vdss,today_removable_ml_vds);
            decimate_rp = intersect(remove_vdss,this.removable_rp_vds);
            X = struct('day',day,'remove_ml_vdss',decimate_ml,'remove_rp_vdss',decimate_rp);            

        end
        
        function [this] = record_num_samples(this,remove_vdss)
            num_ml = sum(ismember(remove_vdss,this.all_ml_vds));
            num_rp = sum(ismember(remove_vdss,this.all_rp_vds));
            m = num_ml+1;
            r = num_rp+1;
            this.num_samples_grid(m,r) = this.num_samples_grid(m,r) + 1;
        end
        
        function [this] = report_estimation(this,result,output_file)
            
            % escape if there are no estimations
            if ~isfield(result,'estimation_vpm') || all(cellfun(@(x)isempty(x),result.estimation_vpm))
                return
            end
            
            if all(cellfun(@(x)isempty(x),result.estimation_vpm))
                return
            end
            
            % save to ppt?
            save_ppt = nargin>2 && ~isempty(output_file);
            if save_ppt
                figvis = 'off';
                [ppt,op]=openppt([output_file '_estimation'],true);
            else
                figvis = 'on';
                ppt = [];
                op = [];
            end
            
            % prediction time
            [start_time,update_dt,measurement_dt,end_time,pred_time] = this.get_time_info;
            
            % get pems data provider
            pems_dp = Utils.get_pems_dp(this.ni,this.configname);
            
            for i=1:length(result.estimation_vpm)
                
                fprintf('%d of %d\n',i,length(result.estimation_vpm))
                
                if isempty(result.estimation_vpm{i})
                    continue
                end
                
                % measured data
                good_vdss = result.good_vdss(i).ids;
                meas_data = pems_dp.get_data(this.days(i),good_vdss,start_time,end_time,measurement_dt);
                meas_dty_vpm = vertcat(meas_data.dty_vpk)*1.60934;
                meas_dty_vpm = meas_dty_vpm(:,2:end);
                meas_dty_vpm(isnan(meas_dty_vpm)) = 0;
                n = update_dt/measurement_dt;
                meas_dty_vpm = meas_dty_vpm(:,n:n:end);
                
                Utils.contour_plot( ...
                    figvis , ...
                    pred_time , ...
                    1:numel(good_vdss) , ...
                    meas_dty_vpm , ...
                    sprintf('Measured %s',datestr(this.days(i))) , ...
                    op)
                c = caxis;
                colorbar
                
                good_link_ids = this.ni.vds2linkid(2,Utils.index_into(good_vdss,this.ni.vds2linkid(1,:)));
                ml_ind = Utils.index_into(good_link_ids,this.ni.link_ids);
                
                Utils.contour_plot( ...
                    figvis , ...
                    pred_time , ...
                    1:numel(ml_ind) , ...
                    result.estimation_vpm{i}(:,ml_ind)' , ...
                    sprintf('Estimated %s',datestr(this.days(i)))  , ...
                    op)
                caxis(c);
                colorbar
                
            end
            
            if save_ppt
                closeppt(ppt,op)
            end
            
        end
        
        function [this] = report_one_hour_flow_error(this,result,output_file)
                        
            % save to ppt?
            save_ppt = nargin>2 && ~isempty(output_file);
            if save_ppt
                figvis = 'off';
                [ppt,op]=openppt([output_file '_one_hour_flow_error'],true);
            else
                figvis = 'on';
            end
            
            % parameters
            if save_ppt
                addslideText(op,'Parameters',evalc('disp(this.predictor_params)'))
            else
                fprintf('Parameters:\n')
                disp(this.predictor_params)
            end
            
            % days
            if save_ppt
                z = cell(1,2*length(this.days)-1);
                z(2:2:end) = {'\n'};
                z(1:2:end) = cellfun(@(x)datestr(x),num2cell(this.days),'UniformOutput',false);
                addslideText(op,sprintf('%d days',length(this.days)),sprintf(horzcat(z{:})))
            else
                fprintf('# days = %d\n',length(this.days))
            end
            
%             % VDSs
%             figure('Visible',figvis)
%             z = this.num_good_ml_vds_for_day;
%             Utils.hist(z,'# good mainline VDSs','normal',max(z)-min(z)+1);
%             if save_ppt
%                 addslide(op,'# good mainline VDSs per day',[],[],[],0.5)
%                 close
%             end
            
            % error
            flw_error = PredictionJobManager.extract_flow_error(result.sixty_minute_prediction);
            
            % box plot time-o-day
            figure('Visible',figvis,'Position',[104 267 1228 376])
            boxplot(flw_error.as_matrix,'outliersize',4,'symbol','m.','colors','k','boxstyle','outline','width',0.6)
            set(findobj(gca,'tag','Box'),'LineWidth',2)
            set(findobj(gca,'tag','Median'),'LineWidth',2)
            set(findobj(gca,'tag','Lower Adjacent Value'),'LineWidth',2)
            set(findobj(gca,'tag','Lower Whisker'),'LineWidth',2,'LineStyle','-')
            set(findobj(gca,'tag','Upper Adjacent Value'),'LineWidth',2)
            set(findobj(gca,'tag','Upper Whisker'),'LineWidth',2,'LineStyle','-')
            if save_ppt
                addslide(op,'Error for time of day',[],[],[],0.5)
                close
            end
            
            % histogram of flow errors
            %             Utils.hist(flw_error.as_vector*100,'error [%]','normal');
            Utils.hist(flw_error.as_vector*100,'error [%]','exponential');
            if save_ppt
                addslide(op,'Distribution of the prediction error',[],[],[],0.5)
                close
            end
            
            if save_ppt
                closeppt(ppt,op)
            end
            
        end
        
        function [this] = report_prediction_contours(this,result,output_file)
            
            save_ppt = nargin>2 && ~isempty(output_file);
            
            if save_ppt
                figvis = 'off';
                [ppt,op]=openppt([output_file '_prediction_contours'],true);
            else
                figvis = 'on';
                ppt = [];
                op = [];
            end
            
            % Create a data provider for making data contours
            dp = Utils.get_pems_dp(this.ni,this.configname);
            
            % fwy order
            x = this.ni.scenario.get_ordered_vds;
            ordered_vdss = [x.sensor_vds];
            clear x
            
            [start_time,update_dt,measurement_dt,end_time,pred_time] = this.get_time_info;
            
            for i=1:length(this.days)
                
                fprintf('%d of %d\n',i,length(this.days))                
                
                day = this.days(i);
                vdss = result.good_vdss(i).ids;
                ind = Utils.index_into(ordered_vdss,vdss);
                ind(ind==0) = [];
                ordered_vdss = vdss(ind);
                
                % measured data
                meas_data = dp.get_data(day,ordered_vdss,start_time,end_time,measurement_dt);
                meas_dty_vpm = vertcat(meas_data.dty_vpk)*1.60934;
                meas_dty_vpm = meas_dty_vpm(:,2:end);
                meas_dty_vpm(isnan(meas_dty_vpm)) = 0;
                n = update_dt/measurement_dt;
                meas_dty_vpm = meas_dty_vpm(:,n:n:end);
                
                if save_ppt
                    addslideTitle(op,{datestr(day)})
                end
                
                % plots
                
                Utils.contour_plot( ...
                    figvis , ...
                    pred_time , ...
                    ordered_vdss , ...
                    meas_dty_vpm , ...
                    ['Measured density for ' datestr(day)] , ...
                    op );
                
                Utils.contour_plot( ...
                    figvis , ...
                    pred_time, ...
                    ordered_vdss , ...
                    result.five_minute_prediction(i).dty_vpm(ind,:) , ...
                    ['Five minute density prediction for ' datestr(day)] , ...
                    op );
                
                Utils.contour_plot( ...
                    figvis , ...
                    pred_time, ...
                    ordered_vdss , ...
                    circshift( result.thirty_minute_prediction(i).dty_vpm(ind,:) , 1 , 2 ), ...
                    ['Thirty minute density prediction for ' datestr(day)] , ...
                    op );
                
                Utils.contour_plot( ...
                    figvis , ...
                    pred_time, ...
                    ordered_vdss , ...
                    circshift( result.sixty_minute_prediction(i).dty_vpm(ind,:) , 2 , 2 ), ...
                    ['Sixty minute density prediction for ' datestr(day)] , ...
                    op );
                
            end
            
            if save_ppt
                closeppt(ppt,op)
            end
            
        end
        
        function [this] = report_demand_prediction(this,result,output_file)
            
            % escape if there are no demands
            if ~isfield(result,'predicted_demand') || all(cellfun(@(x)isempty(x),{result.predicted_demand.demand}))
                return
            end
            
            save_ppt = nargin>2 && ~isempty(output_file);
            
            if save_ppt
                figvis = 'off';
            else
                figvis = 'on';
            end

            for i=1:length(this.days)

                if save_ppt
                    [ppt,op]=openppt(sprintf('%d_demand_prediction_%d',output_file,this.days(i)),true);
                else
                    ppt = [];
                    op = [];
                end

                fprintf('%d of %d\n',i,length(this.days))
                
                day = this.days(i);
                
                if save_ppt
                    addslideTitle(op,{datestr(day)})
                end
                
                % plots
                num_ramp=size(result.predicted_demand(i).demand,2);
                for j=1:num_ramp
                    fprintf('\t%d of %d\n',j,num_ramp)
                    Utils.demand_sr_plot( ...
                        figvis , ...
                        result.predicted_demand(i).demand(j).horizon, ...
                        result.predicted_demand(i).demand(j).time , ...
                        result.predicted_demand(i).demand(j).flw_vph , ...
                        ['Predicted demands for ' datestr(day) ' & Link ID: ' num2str(result.predicted_demand(i).demand(j).link_id)] , ...
                        op );
                end
                
                if save_ppt
                    closeppt(ppt,op)
                end
                
            end
            
        end
        
        function [this] = report_sr_prediction(this,result,output_file)
            
            % escape if there are no split ratios
            if ~isfield(result,'predicted_split_ratio') || all(cellfun(@(x)isempty(x),{result.predicted_split_ratio.split_ratio}))
                return
            end
            
            save_ppt = nargin>2 && ~isempty(output_file);
            
            if save_ppt
                figvis = 'off';
            else
                figvis = 'on';
            end
            
            for i=1:length(this.days)
                
                if save_ppt
                    [ppt,op]=openppt(sprintf('%s_split_ratio_prediction_%d',output_file,this.days(i)),true);
                else
                    ppt = [];
                    op = [];
                end
                                
                fprintf('%d of %d\n',i,length(this.days))

                day = this.days(i);
                
                if save_ppt
                    addslideTitle(op,{datestr(day)})
                end
                
                % plots
                num_ramp=size(result.predicted_split_ratio(i).split_ratio,2);
                for j=1:num_ramp
                    Utils.demand_sr_plot( ...
                        figvis , ...
                        result.predicted_split_ratio(i).split_ratio(j).horizon, ...
                        result.predicted_split_ratio(i).split_ratio(j).time , ...
                        result.predicted_split_ratio(i).split_ratio(j).split_ratios , ...
                        ['Predicted split ratios for ' datestr(day) ' & 1st D-Link ID: ' ...
                        num2str(result.predicted_split_ratio(i).split_ratio(j).links_out(1))] , ...
                        op );
                end
                
                if save_ppt
                    closeppt(ppt,op)
                end
            end

        end
        
        function [start_time,update_dt,measurement_dt,end_time,pred_time] = get_time_info(this)
            
            start_time = this.predictor_params.start_time;
            update_dt = this.predictor_params.update_dt;
            measurement_dt = this.predictor_params.measurement_dt;
            end_time = this.predictor_params.end_time;
            pred_time = start_time:update_dt:end_time;
            pred_time = pred_time(2:end)/3600;
            
            %             % aggregation matrix
            %             n = length(pred_time);
            %             m = (end_time-start_time)/measurement_dt;
            %             p = m/n;
            %             aggmatrix = zeros(m,n);
            %             for i=1:n
            %                 aggmatrix(((i-1)*p+1):(i*p),i) = 1/p;
            %             end
            
        end
        
    end
    
end
