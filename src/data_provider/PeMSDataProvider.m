classdef PeMSDataProvider < DataProvider

    properties (Access = public) 
        use_template    % true if the vds does not correspond to an actual PeMS station
        templates
    end
    
    properties(Constant)
        LOOP_HEALTH_THRESHOLD = 90;
        MEASUREMENT_DT = 300;           % 5 minute data
    end
    
    methods( Access = public )
        
        function [this] = PeMSDataProvider(days,internal_ids,use_template,name,link_ids,length_for_ids_km,cluster_method)
            % [this] = PeMSDataProvider(days,internal_ids,use_template,name,link_ids,length_for_ids_km,cluster_method)
            
            this = this@DataProvider(days,internal_ids,link_ids,length_for_ids_km,PeMSDataProvider.MEASUREMENT_DT);
                        
            %  use template for these regardless of health
            this.use_template = ismember(internal_ids,[use_template.id]);
            
            % load traffic data
            pems = PeMS5minData;
            pems.load(Folder.data,this.ids,this.days);

            % load loop health info
            load(fullfile(Folder.data,sprintf('%s_loop_health',name)))

            % check for which I have no pems data
            have_pems_data = cellfun( @(x) ~isempty(x) , pems.data );
            have_pems_data = any(have_pems_data,2)';
            
            % error out if there are vdss for which I have neither pems nor template
            if any(~have_pems_data & ~this.use_template)
                error('any(~have_pems_data & ~this.use_template)')
            end
            
            % error out if there are real pems stations with no health data.
            if any(~ismember(pems.vds(~this.use_template),loop_health.vds))
                error('found pems stations with no health data')
            end

            % load templates data
            load(fullfile(Folder.data,sprintf('%s_template',name)),'template');
            this.templates = template;
            clear template
            
            % extract health data from loop_health
            loop_percent_observed = nan(length(this.ids),length(loop_health.days));
            
            % copy from loop_health for pems vds; 
            ind = Utils.index_into(this.ids,loop_health.vds);
            loop_percent_observed(ind>0,:) = loop_health.percent_observed(ind(ind>0),:);

            % templates are good
            loop_percent_observed(this.use_template,:) = 100;

            % organize days
            loop_percent_observed = loop_percent_observed(:,Utils.index_into(this.days,loop_health.days));
           
            % keep only boolean
            this.vds_is_good = loop_percent_observed' >= PeMSDataProvider.LOOP_HEALTH_THRESHOLD;
            
            % no one has been forced bad
            this.vds_forced_bad = false(size(this.vds_is_good,1),size(this.vds_is_good,2));
            
            clear loop_health loop_percent_observed
            
            % get aggregate pems data
            x = pems.get_data_batch_aggregate(this.ids,this.days,'var',{'flw','dty','spd'},'fill',true);
            
            % distribute into data
            for d=1:length(this.days)
               for i=1:length(this.ids)
                   this.data(d,i).flw_out_vph = x.flw(:,i,d)';
                   this.data(d,i).dty_vpk = ([x.dty(1,i,d) ; x.dty(:,i,d)]')/1.60934;
               end
            end
            clear x
            
            this.time = (0:288)*PeMSDataProvider.MEASUREMENT_DT;
            
            % create a cluster manager
            this.cluster_manager = ClusterManager( ...
                this.ids,this.days,this.time,this.data,cluster_method);
        end
        
        function [x] = get_health(this,getdays,getvds)
            % [x] = get_health(this,getdays,getvds)
            
            if isnan(getdays) & numel(this.days)>1
                error('isnan(getdays) & numel(this.days)>1')
            end
            
            if ischar(getdays) && strcmp(getdays,'all')
                getdays = this.days;
            end
                
            if ischar(getvds) && strcmp(getvds,'all')
                getvds = this.ids;
            end
            
            if isnan(getdays)
                days_ind = 1;
            else
                days_ind = Utils.index_into(getdays,this.days);
            end
            
            vds_ind = Utils.index_into(getvds,this.ids);
            
            x_is_good = this.vds_is_good(days_ind,vds_ind);
            x_forced_bad = this.vds_forced_bad(days_ind,vds_ind);
            
            x = x_is_good & ~x_forced_bad;
            
        end
         
        function [x] = uses_template(this,getdays,getvds)
            
            if numel(getdays)>1
                error('this has not been implemented.')
            end
            
            % use a template if a) the vds appears in use_template, or
            % b) the vds is bad and we have a template
            force_use_template = ismember(getvds,this.ids(this.use_template));            
            is_good = this.get_health(getdays,getvds);
            have_template = ismember(getvds,[this.templates.vds]);
            
            x = force_use_template | (~is_good & have_template);
        end
        
        function [x] = get_data(this,varargin)
            %  days,ids,from,to,dt,method (interpolate or snap)
                
            [getdays,getids,from,to,dt,method] = this.process_get_data_input(varargin);
            x = this.get_data@DataProvider(varargin{:});
            r_time = from:dt:to;
            
            % get health
            health = this.get_health(getdays,getids);
            all_temp_vdss = [this.templates.vds];

            % replace templates into bad days
            for i = 1:length(getids)
                
                force_template = ismember(getids(i),this.ids(this.use_template));
                template_index = Utils.index_into(getids(i),all_temp_vdss);
                
                % no template
                if template_index==0
                    continue
                end
                
                bad_days = ~health(:,i);
                
                if any(bad_days) || force_template
                    
                    % get template values for this vds
                    this_template = this.templates(template_index);
                    flw_vph = [];
                    if numel(r_time)>1
                        switch method
                            case 'interpolate'
                                flw_vph = interp1(this_template.time_sec,this_template.flow_vph,r_time(1:end-1),'previous','extrap');
                            case 'snap'
                                r_time_ind = arrayfun( @(x) find(this_template.time_sec<=x,1,'last') , r_time);
                                n = length(this_template.flow_vph);
                                flw_ind = r_time_ind(1:end-1);
                                flw_ind(flw_ind>n)=n;
                                flw_vph = this_template.flow_vph(flw_ind);
                            otherwise
                                error('unknown method')
                        end
                    end

                    % replace bad days in x with template
                    if force_template
                        for d=1:size(x,1)
                            x(d,i).flw_in_vph = flw_vph;
                            x(d,i).flw_out_vph = flw_vph;
                        end
                    else
                        bad_days_ind = find(bad_days);
                        for d=1:length(bad_days_ind)
                            x(bad_days_ind(d),i).flw_in_vph = flw_vph;
                            x(bad_days_ind(d),i).flw_out_vph = flw_vph;
                        end
                    end
                 
                end
            end

        end
        
        function [] = export_health(this,filename)
            T=array2table(this.vds_is_good , ...
                'VariableNames', arrayfun(@(x) sprintf('vds%d',x),this.ids,'UniformOutput',false) , ...
                'RowNames',arrayfun(@(x) datestr(x,'mmm_dd_yyyy'),this.days,'UniformOutput',false) );
            writetable(T,filename,'WriteRowNames',true)
        end
        
        function [this] = overwrite_force_bad_sensors(this,Z)
            if ~all(size(Z)==size(this.vds_forced_bad))
                error('~all(size(Z)==size(this.vds_forced_bad))')
            end
            this.vds_forced_bad = Z;
        end
            
        function [this] = set_force_bad_sensors(this,this_day,bad_vdss)
            
            this.set_force_bad_sensors@DataProvider(this_day,bad_vdss);
            
            % reset vds_forced_bad
            this.vds_forced_bad = false(size(this.vds_is_good,1),size(this.vds_is_good,2));

            % check day
            day_ind = Utils.index_into(this_day,this.days);
            if day_ind==0
                error('day_ind==0')
            end
            
            % loop through vdss
            for i=1:numel(bad_vdss)
                vds_ind = Utils.index_into(bad_vdss(i),this.ids);
                if any(vds_ind==0)
                    error('any(vds_ind==0)')
                end
                this.vds_forced_bad(day_ind,vds_ind) = true;
            end
            
        end
        
    end
    
end

