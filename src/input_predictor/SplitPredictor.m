classdef SplitPredictor < InputPredictor
    
    properties
        prediction_day  % the prediction day has to be set so that the
        % topology of good sensors can be determined before
        % predictions are made
        dp_info
        all_dps     % list of all data provide ids used for prediction of splits
    end
    
    methods(Access=public)
        
        function [this] = SplitPredictor(ni,flow_predictor)
            
            this = this@InputPredictor(flow_predictor);
            
            node_io = ni.get_node_io_for_multi_output_nodes;
            this.ids = [node_io.id];
            if ~all(arrayfun(@(x)length(x.link_in),node_io)==1)
                error('~all(arrayfun(@(x)length(x.link_in),node_io)==1)')
            end
            
            % the split ratio is calculated as the plus/minus of flows from
            % numerator dps (data providers) divided by the plus/minus of
            % flows from denominator dps.
            this.all_dps = [];
            this.dp_info = repmat(struct('in_struct',[],'out_struct',[]),1,length(this.ids));
            dp_struct = struct('link_id',nan,'has_good_sensors',nan,'my_dps',[],'is_terminal',nan,'plus_dps',[],'minus_dps',[]);
            
            for i=1:length(this.ids)
                
                nIn = length(node_io(i).link_in);
                nOut = length(node_io(i).link_out);
                
                this.dp_info(i).in_struct = repmat(dp_struct,1,nIn);
                this.dp_info(i).out_struct = repmat(dp_struct,1,nOut);
                
                for j=1:nIn
                    link_id = node_io(i).link_in(j);
                    dpids=flow_predictor.data_provider.get_ids_for_linkids(link_id);
                    this.dp_info(i).in_struct(j).my_dps = dpids{1};
                    this.dp_info(i).in_struct(j).link_id = link_id;
                    this.dp_info(i).in_struct(j).is_terminal = ni.get_link_is_terminal(link_id);
                end
                
                for j=1:nOut
                    link_id = node_io(i).link_out(j);
                    dpids=flow_predictor.data_provider.get_ids_for_linkids(link_id);
                    this.dp_info(i).out_struct(j).my_dps = dpids{1};
                    this.dp_info(i).out_struct(j).link_id = link_id;
                    this.dp_info(i).out_struct(j).is_terminal = ni.get_link_is_terminal(link_id);
                end
                
            end
            
        end
        
        function [] = set_day(this,ni,day)
            
            if ~this.check_have_day(day)
                error('~this.check_have_day(day)')
            end
                        
            this.prediction_day = day;
            source_sink_link_ids = setdiff(ni.link_ids,ni.get_internal_link_ids);
            
            % get plus/minus dps corresponding to inputs and outputs of each node
            for i=1:length(this.ids)
                
                for j=1:length(this.dp_info(i).in_struct)
                    link_id = this.dp_info(i).in_struct(j).link_id;
                    is_terminal =  this.dp_info(i).in_struct(j).is_terminal;
                    
                    % check that termainal links have sensors
                    [has_sensors,has_good_sensors] = this.get_has_good_sensors(day,link_id);
                    if is_terminal && ~has_sensors
                        error('is_terminal && ~has_sensors')
                    end
                    
                    % plus minus for non-terminal links
                    [plus_dps,minus_dps] = this.find_linear_combination_dps(ni,day,link_id,source_sink_link_ids);
                    
                    % store
                    this.dp_info(i).in_struct(j).plus_dps = plus_dps;
                    this.dp_info(i).in_struct(j).minus_dps = minus_dps;
                    this.dp_info(i).in_struct(j).has_good_sensors = has_good_sensors;
                    this.all_dps = [this.all_dps plus_dps minus_dps];
                end
                
                for j=1:length(this.dp_info(i).out_struct)
                    link_id = this.dp_info(i).out_struct(j).link_id;
                    is_terminal =  this.dp_info(i).out_struct(j).is_terminal;
                    
                    % check that termainal links have sensors
                    [has_sensors,has_good_sensors] = this.get_has_good_sensors(day,link_id);
                    if is_terminal && ~has_sensors
                        error('is_terminal && ~has_sensors')
                    end
                    
                    % plus minus for non-terminal links
                    [plus_dps,minus_dps] = this.find_linear_combination_dps(ni,day,link_id,source_sink_link_ids);

                    % store
                    this.dp_info(i).out_struct(j).plus_dps = plus_dps;
                    this.dp_info(i).out_struct(j).minus_dps = minus_dps;
                    this.dp_info(i).out_struct(j).has_good_sensors = has_good_sensors;
                    this.all_dps = [this.all_dps plus_dps minus_dps];

                end
                
            end
            
            this.all_dps = unique(this.all_dps);
        end
        
        function [R] = predict(this,day,from,to,dt)
            
            if isempty(this.prediction_day)
                error('set prediction day first')
            end
            
            if isnan(day)
                day=this.prediction_day;
            end
            
            if day~=this.prediction_day
                error('day~=this.prediction_day')
            end
            
            if nargin<5
                dt = nan;
            end
            
            % get predictions for all dp ids in all_dps
            [all_dp_flw,time] = this.flow_predictor.predict_flw_for_dpids(this.all_dps,day,from,to,dt);
            
            R = repmat(struct('node_id',nan,'links_in',[],'links_out',[],'time',[],'split_ratios',{}),1,length(this.ids));
            for i=1:length(this.ids)
                
                nIn = length(this.dp_info(i).in_struct);
                nOut = length(this.dp_info(i).out_struct);
                
                R(i).node_id = this.ids(i);
                R(i).links_in = [this.dp_info(i).in_struct.link_id];
                R(i).links_out = [this.dp_info(i).out_struct.link_id];
                R(i).time = time;
                R(i).split_ratios = cell(nIn,nOut);
                
                % compute split ratio
                sr = this.compute_split_from_dp_info(this.dp_info(i),all_dp_flw);
                
                % distribute
                for j=1:nIn
                    for k=1:nOut
                        R(i).split_ratios{j,k} = sr(k,:);
                    end
                end
                
            end
        end
        
        function [R] = predict_given_flows(this,flow_pred)
           
            % get predictions for all dp ids in all_dps
            time=flow_pred(end).time;
            all_dp_flw=zeros(length(this.all_dps),length(flow_pred(end).flow));
            for i=1:length(this.all_dps)
                ind=ismember([flow_pred.dp_ids],this.all_dps(i));
                all_dp_flw(i,:)=flow_pred(ind).flow;        
            end
            
            R = repmat(struct('node_id',nan,'links_in',[],'links_out',[],'time',[],'split_ratios',{}),1,length(this.ids));
            for i=1:length(this.ids)
                
                nIn = length(this.dp_info(i).in_struct);
                nOut = length(this.dp_info(i).out_struct);
                
                R(i).node_id = this.ids(i);
                R(i).links_in = [this.dp_info(i).in_struct.link_id];
                R(i).links_out = [this.dp_info(i).out_struct.link_id];
                R(i).time = time;
                R(i).split_ratios = cell(nIn,nOut);
                
                % compute split ratio
                sr = this.compute_split_from_dp_info(this.dp_info(i),all_dp_flw);
                
                % distribute
                for j=1:nIn
                    for k=1:nOut
                        R(i).split_ratios{j,k} = sr(k,:);
                    end
                end
                
            end
        end
        
        function [sr] = compute_split_from_dp_info(this,this_dp_info,all_dp_flw)
            
            nOut = length(this_dp_info.out_struct);
            
            % compute in/out flow
            flw_in  = compute_plus_minus_flw(this_dp_info.in_struct);
            flw_out = compute_plus_minus_flw(this_dp_info.out_struct);
            
            % normalize
            tot_in = sum(flw_in,1);
            tot_out = sum(flw_out,1);
            des_tot_flw = mean([tot_in;tot_out],1);
            flw_out = flw_out.*repmat(des_tot_flw./tot_out,nOut,1);
            
            % compute split ratios
            sr = flw_out ./ repmat(des_tot_flw,nOut,1);
            
            function [f] = compute_plus_minus_flw(strct)
                
                f = nan( length(strct),size(all_dp_flw,2) );
                for ii=1:length(strct)
                    
                    % case that the link itself has sensors
                    if strct(ii).has_good_sensors
                        my_ind = Utils.index_into(strct(ii).my_dps,this.all_dps);
                        f(ii,:) = sum(all_dp_flw(my_ind,:),1);
                        continue
                    end
                    
                    % otherwise do the plus/minus calculation
                    plus_ind = Utils.index_into(strct(ii).plus_dps,this.all_dps);
                    minus_ind = Utils.index_into(strct(ii).minus_dps,this.all_dps);
                    net_flw = sum(all_dp_flw(plus_ind,:),1) - sum(all_dp_flw(minus_ind,:),1);
                    
                    if strct(ii).is_terminal
                        f(ii,:) = net_flw;
                    else
                        f(ii,:) = net_flw / 2;
                    end
                end
                
                % eliminate negative values
                f = max(f,0*f);
                
            end
            
        end
        
    end
    
    methods(Access=private)
        
        function [has_sensors,has_good_sensors] = get_has_good_sensors(this,day,start_link)
            dp = this.flow_predictor.data_provider;
            
            % check if this link has a good sensor. If so, stop
            sensor = dp.get_ids_for_linkids(start_link);
            sensor = sensor{1};
            has_sensors = ~isempty(sensor);
            if has_sensors
                has_good_sensors = all(dp.get_health(day,sensor));
            else
                has_good_sensors = false;
            end
            
        end
        
        function [plus_dps,minus_dps] = find_linear_combination_dps(this,ni,day,start_link,source_sink_link_ids)
            
            dp = this.flow_predictor.data_provider;
            
            % otherwise search the network for equivalent set of sensors
            open_leaves = start_link;
            contribution_sign = containers.Map('KeyType','double','ValueType','double');
            SplitPredictor.set_value(contribution_sign,start_link,1);

            plus_dps = [];
            minus_dps = [];
            visited = [];
            while ~isempty(open_leaves)
                                
                visited_neighbors = [];
                unvisited_neighors = [];
                for i=1:length(open_leaves)
                                        
                    this_link = open_leaves(i);
                    
                    % get all links adjacent to this link
                    [up_node,dn_node] = ni.get_begin_end_node_for_link(this_link);
                    [up_enter,up_exit] = ni.get_in_out_links_for_node(up_node);
                    [dn_enter,dn_exit] = ni.get_in_out_links_for_node(dn_node);
                    up_exit  = setdiff(up_exit,this_link);
                    dn_enter = setdiff(dn_enter,this_link);

                    % discard already visited links
                    up_enter = setdiff(up_enter,visited);
                    up_exit  = setdiff(up_exit,visited);
                    dn_enter = setdiff(dn_enter,visited);
                    dn_exit  = setdiff(dn_exit,visited);
                    
                    % get sensors on those links
                    up_enter_sensors = dp.get_ids_for_linkids(up_enter);
                    up_exit_sensors  = dp.get_ids_for_linkids(up_exit);
                    dn_enter_sensors = dp.get_ids_for_linkids(dn_enter);
                    dn_exit_sensors  = dp.get_ids_for_linkids(dn_exit);
                    
                    % get the health of sensors
                    up_enter_health = cellfun(@(x)dp.get_health(day,x),up_enter_sensors,'UniformOutput',false);
                    up_exit_health = cellfun(@(x)dp.get_health(day,x),up_exit_sensors,'UniformOutput',false);
                    dn_enter_health = cellfun(@(x)dp.get_health(day,x),dn_enter_sensors,'UniformOutput',false);
                    dn_exit_health = cellfun(@(x)dp.get_health(day,x),dn_exit_sensors,'UniformOutput',false);
                    
                    % this sensor uses a template
                    up_enter_use_template = cellfun(@(x)dp.uses_template(day,x),up_enter_sensors,'UniformOutput',false);
                    up_exit_use_template = cellfun(@(x)dp.uses_template(day,x),up_exit_sensors,'UniformOutput',false);
                    dn_enter_use_template = cellfun(@(x)dp.uses_template(day,x),dn_enter_sensors,'UniformOutput',false);
                    dn_exit_use_template = cellfun(@(x)dp.uses_template(day,x),dn_exit_sensors,'UniformOutput',false);
                    
                    % good sensor or template
                    has_good_sensor_or_template = @(s,h,t) ~isempty(s) & all(h|t);  % has sensors and all are good
                    up_enter_has_good_sensors_or_template = cellfun(has_good_sensor_or_template, ...
                        up_enter_sensors , ...
                        up_enter_health , ...
                        up_enter_use_template ); 
                    up_exit_has_good_sensors_or_template  = cellfun(has_good_sensor_or_template, ...
                        up_exit_sensors , ...
                        up_exit_health , ...
                        up_exit_use_template);
                    dn_enter_has_good_sensors_or_template = cellfun(has_good_sensor_or_template, ...
                        dn_enter_sensors , ...
                        dn_enter_health , ...
                        dn_enter_use_template);
                    dn_exit_has_good_sensors_or_template  = cellfun(has_good_sensor_or_template, ...
                        dn_exit_sensors , ...
                        dn_exit_health , ...
                        dn_exit_use_template);
                    
                    % check whether any of bad sensors are on sources or sinks.
                    up_enter_is_bndry = ismember(up_enter,source_sink_link_ids);
                    up_exit_is_bndry = ismember(up_exit,source_sink_link_ids);
                    dn_enter_is_bndry = ismember(dn_enter,source_sink_link_ids);
                    dn_exit_is_bndry = ismember(dn_exit,source_sink_link_ids);
                    
                    if any(up_enter_is_bndry & ~up_enter_has_good_sensors_or_template)
                        warning('Boundary link %d has bad sensor and lacks a template.\n', ...
                            up_enter(up_enter_is_bndry & ~up_enter_has_good_sensors_or_template))
                    end
                    
                    if any(up_exit_is_bndry & ~up_exit_has_good_sensors_or_template)
                        warning('Boundary link %d has bad sensor and lacks a template.\n', ...
                            up_exit(up_exit_is_bndry & ~up_exit_has_good_sensors_or_template))
                    end
                    
                    if any(dn_enter_is_bndry & ~dn_enter_has_good_sensors_or_template)
                        warning('Boundary link %d has bad sensor and lacks a template.\n', ...
                            dn_enter(dn_enter_is_bndry & ~dn_enter_has_good_sensors_or_template))
                    end
                    
                    if any(dn_exit_is_bndry & ~dn_exit_has_good_sensors_or_template)
                        warning('Boundary link %d has bad sensor and lacks a template.\n', ...
                            dn_exit(dn_exit_is_bndry & ~dn_exit_has_good_sensors_or_template))
                    end
                    
                    % add good sensors to the plus/minus lists
                    this_plus_dps = Utils.flatten_cell_array( ...
                        [up_enter_sensors(up_enter_has_good_sensors_or_template) ...
                        dn_exit_sensors(dn_exit_has_good_sensors_or_template)]); 
                    
                    this_minus_dps = Utils.flatten_cell_array( ...
                        [up_exit_sensors(up_exit_has_good_sensors_or_template) ...
                        dn_enter_sensors(dn_enter_has_good_sensors_or_template)]);
                    
                    if contribution_sign(this_link)>0   % contributes positively to the flow
                        plus_dps = [plus_dps this_plus_dps];
                        minus_dps = [minus_dps this_minus_dps];
                        SplitPredictor.set_value(contribution_sign,[up_enter dn_exit],1);
                        SplitPredictor.set_value(contribution_sign,[up_exit dn_enter],-1);
                    else                    % contributes negatively to the flow
                        plus_dps = [plus_dps this_minus_dps];
                        minus_dps = [minus_dps this_plus_dps];
                        SplitPredictor.set_value(contribution_sign,[up_enter dn_exit],-1);
                        SplitPredictor.set_value(contribution_sign,[up_exit dn_enter],1);
                    end
                    
                    % add neighbors with good sensors to visited neighbors
                    visited_neighbors = [visited_neighbors ...
                        [up_enter(up_enter_has_good_sensors_or_template) ...
                        up_exit(up_exit_has_good_sensors_or_template) ...
                        dn_enter(dn_enter_has_good_sensors_or_template) ...
                        dn_exit(dn_exit_has_good_sensors_or_template)] ];
                    
                    % add neighbors without detection that are not in visited to unvisited neighbors
                    x = [up_enter(~up_enter_has_good_sensors_or_template) ...
                        up_exit(~up_exit_has_good_sensors_or_template) ...
                        dn_enter(~dn_enter_has_good_sensors_or_template) ...
                        dn_exit(~dn_exit_has_good_sensors_or_template)];
                    unvisited_neighors = [unvisited_neighors setdiff(x,visited)];
                    
                end
                
                % add the visited neighbors previously open leaves to visited
                visited = unique([visited open_leaves visited_neighbors]);
                
                unvisited_neighors = setdiff(unvisited_neighors,visited);
                
                % open leaves are now the unvisited neighbors
                open_leaves = unique(unvisited_neighors);
                
            end
            
            plus_dps = plus_dps(~isnan(plus_dps));
            minus_dps = minus_dps(~isnan(minus_dps));
        end
                
        function [X] = genrerate_report(this,day,ppt_file,ppt_file_error,xls_file,sim_dp,update_dt,horizon,historical)
            
            dt = this.flow_predictor.data_provider.measurement_dt;
            update_times = 0:update_dt:(86400-update_dt);
            
            % allocate stuff
            actual = repmat(struct( 'node_id',nan, ...
                'link_in',nan, ...
                'links_out',[], ...
                'split_ratios',[]), ...
                1,length(this.ids));
            
            prediction = repmat( struct('node_id',nan, ...
                'error', struct() , ...
                'split_ratio', []) , ...
                1,length(this.ids));
            
            % compute true splits from simulation data provider
            for i=1:length(this.ids)
                
                % io information
                prediction(i).node_id = this.ids(i);
                prediction(i).link_in = this.dp_info(i).in_struct.link_id;
                prediction(i).links_out = [this.dp_info(i).out_struct.link_id];
                
                actual(i).node_id = this.ids(i);
                actual(i).link_in = this.dp_info(i).in_struct.link_id;
                actual(i).links_out = [this.dp_info(i).out_struct.link_id];
                
                % compute actual splits from simulation output
                in_link_data = sim_dp.get_data('all',actual(i).link_in,'start','end');
                out_link_data = sim_dp.get_data('all',actual(i).links_out,'start','end');
                
                in_flw = in_link_data.flw_out_vph;
                out_flw = vertcat(out_link_data.flw_in_vph);
                
                actual(i).split_ratios = out_flw./repmat(in_flw,length(actual(i).links_out),1);
            end
            actual_time = in_link_data.time;
            clear in_link_data out_link_data in_flw out_flw
            
            % run predictions
            fprintf('\tRunning %d predictions',length(update_times))
            for k=1:length(update_times)
                
                fprintf('.')
                start_time = update_times(k);
                x = this.predict( day, ...
                    start_time, ...
                    start_time+horizon, ...
                    dt);
                
                lasttime = min([actual_time(end)-dt start_time+horizon-dt]);
                ind = Utils.index_into(start_time:dt:lasttime,actual_time);
                
                for i=1:length(this.ids)
                    
                    a = actual(i).split_ratios(:,ind);
                    b = Utils.stack_cell_array(x(i).split_ratios);
                    b = b(:,1:length(ind));
                    
                    err = Utils.error(a,b);
                    prediction(i).error.mae(k) = err.mae;
                    prediction(i).error.mape(k) = err.mape;
                    prediction(i).error.smape(k) = err.smape;
                    prediction(i).error.rmse(k) = err.rmse;
                    prediction(i).error.linf(k) = err.linf;
                    
                    for j=1:length(x(i).links_out)
                        prediction(i).split_ratio{j}(k,:) = x(i).split_ratios{j};
                    end
                    
                end
                
            end
            fprintf('\n')
            
            % error table : NOTE: THIS IS AN INCORRECT AVERAGE, CHANGE LATER
            for i=1:length(this.ids)
                prediction(i).avg_error = Utils.avg_error(prediction(i).error);
            end
            
            X = [prediction.avg_error];
            write(table([prediction.node_id nan]', ...
                [[X.mae]';Utils.meanwithnan([X.mae]')], ...
                [[X.mape]';Utils.meanwithnan([X.mape]')], ...
                [[X.smape]';Utils.meanwithnan([X.smape]')], ...
                [[X.rmse]';Utils.meanwithnan([X.rmse]')], ...
                [[X.linf]';Utils.meanwithnan([X.linf]')], ...
                'VariableNames',{'id','mae','mape','smape','rmse','linf'}) , ...
                xls_file);
            clear X
            
            % plot predictions for each id
            do_export = ~isempty(ppt_file);
            if do_export
                [ppt,op]=openppt(ppt_file,true);
                [ppte,ope]=openppt(ppt_file_error,true);
                
                % title slide
                addslideText(op,'Split prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d sec.\nPrediction horizon: %d sec.\n', ...
                    class(this.flow_predictor),class(this.flow_predictor.data_provider),update_dt,horizon) )
                
                % title slide
                addslideText(ope,'Split prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d\nPrediction horizon: %d\n', ...
                    class(this.flow_predictor),class(this.flow_predictor.data_provider),update_dt,horizon) )
            end
            
            for i=1:length(this.ids)
                
                fprintf('\tSlide %d of %d\n',i,length(this.ids))
                
                if do_export
                    figure('Position',[279    95   710   571],'Visible','off')
                else
                    figure('Position',[279    95   710   571],'Visible','on')
                end
                
                xlab = 'time [hr]';
                
                numOut = length(prediction(i).links_out);
                
                % plot historical splits
                if ~isempty(historical)
                    for j=1:numOut
                        subplot(numOut,1,j)
                        plot((0:300:86100)/3600,vertcat(historical.split_ratios{i}(:,:,j)),'k');
                    end
                end
                
                % actual
                for j=1:numOut
                    subplot(numOut,1,j)
                    hold on
                    plot(actual_time(1:end-1)/3600,actual(i).split_ratios(j,:),'r','LineWidth',4)
                end
                
                % prediction
                P = prediction(i);
                for k=1:length(update_times)
                    time = (update_times(k):dt:update_times(k)+horizon-dt)/3600;
                    for j=1:numOut
                        subplot(numOut,1,j)
                        plot(time,P.split_ratio{j}(k,:),'Color',rand(1,3),'LineWidth',2);
                    end
                end
                
                for j=1:numOut
                    subplot(numOut,1,j)
                    xlabel(xlab);
                    ylabel('split ratio');
                    set(gca,'XLim',[0 24])
                    grid
                end
                
                if do_export
                    addslide(op,sprintf('id=%d',this.ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                else
                    title(sprintf('Split predictions, id=%d',this.ids(i)))
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
                    addslide(ope,sprintf('id=%d',this.ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                end
                
            end
            
            if do_export
                closeppt(ppt,op)
                closeppt(ppte,ope)
            end
            
        end
        
    end
    
    methods(Static, Access=public)
        
        function [X] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,config,flow_predictor_type,output_dt,update_dt,horizon,end_time}
            
            X = [];
            config = Config.get(params.config);
            
            % run test of model day
            day = config.model_day;
            
            ni = ObjectFactory.network_information(config.xml_file);
            mr = ObjectFactory.beats_model_runner(ni,config.sim_dt,params.output_dt);
            
            % sim_dp provides the 'actual' data
            sim_dp = ObjectFactory.sim_data_provider(mr,params.end_time);
            
            % create flow predictor for all sensors
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
                otherwise
                    error('bad flow predictor type')
            end
            
            % create the split predictor
            split_pred = SplitPredictor(ni,fp);
            split_pred.set_prediction_day(ni,day);
            clear fp
            
            % load historical data
            switch params.flow_predictor_type
                case {'simulation','zoh'}
                    historical = [];
                case {'historical','scaled_historical'}
                    load(fullfile(Folder.data,sprintf('%s_historical_splits',params.config)))
                otherwise
                    error('bad flow predictor type')
            end
            
            split_pred.genrerate_report( ...
                day , ...
                fullfile(Folder.reports,params.ppt_file) , ...
                sprintf('%s_error',fullfile(Folder.reports,params.ppt_file)) , ...
                fullfile(Folder.reports,sprintf('%s.xls',params.xls_file)) , ...
                sim_dp , ...
                params.update_dt , ...
                params.horizon , ...
                historical )
            
        end
        
    end
    
    methods(Static, Access=private)
        
        function []=set_value(map,keys,value)
            for i=1:length(keys)
                key = keys(i);
                if ~map.isKey(key)
                    map(key) = value;
                end
            end
        end
        
    end
end

