classdef (Abstract) DataProvider < handle
    
    properties (Access = public)
        time                % common time vector (in seconds) for all data profiles
        days                % list of days
        ids                 % list of internal data provider ids (link ids or vdss)
        link_ids            % network link ids corresponding to each internal id.
        % This list may contain duplicate values if
        % there are multiple sensors on a link.
        data                % matrix of DataProfile, size #days by #ids
        measurement_dt      % period of measurement in seconds.
        length_for_ids      % [km] link length corresponding to each id (whether vds or link id)
              
        vds_is_good     % vector of health booleans
        vds_forced_bad
        
        cluster_manager     % e.g. historical data profile by day-of-week
    end
    
    methods(Access=protected)
        
        function [this] = DataProvider(days,internal_ids,link_ids,length_for_ids,measurement_dt)
            
            this.days           = Utils.row_vector(days);
            this.ids            = Utils.row_vector(internal_ids);
            this.link_ids       = Utils.row_vector(link_ids);
            this.length_for_ids = Utils.row_vector(length_for_ids);
            
            this.measurement_dt = measurement_dt;
            this.time = [];
            this.data = repmat(DataProfile(),length(days),length(internal_ids));
            for d=1:length(days)
                for i=1:length(internal_ids)
                    this.data(d,i) = DataProfile();
                end
            end
            
        end
        
    end
    
    methods(Access=public)
        
        
        function [this] = set_force_bad_sensors(this,this_day,bad_vdss)
            % do nothing
        end
        
        function [b] = check_have_day(this,day)
            b = ismember(day,this.days);
        end
        
        function [x] = get_days_and_ids(this)
            x.days = this.days;
            x.ids = this.ids;
        end
        
        function [Y] = get_data(this,varargin)
            %  days,ids,from,to,dt,method (interpolate or snap)
            
            [getdays,getids,from,to,dt,method] = this.process_get_data_input(varargin);
            
            % filter the data by days and ids
            if length(this.days)==1 & isnan(getdays)
                days_ind = 1;
            else
                days_ind = Utils.index_into(getdays,this.days);
            end
            ids_ind = Utils.index_into(getids,this.ids);
            X = this.data(days_ind,ids_ind);
            
            % generate output time vector
            r_time = from:dt:to;
            
            switch method
                case 'interpolate'
                    Y = this.get_data_interpolate(X,r_time,getdays,getids);
                case 'snap'
                    Y = this.get_data_snap(X,r_time,getdays,getids);
                otherwise
                    error('unknown method')
            end
            
        end
        
        function [Y] = get_representative_data(this,varargin)
            %  features,ids,from,to,dt,method (interpolate or snap)
            
            features = varargin{1};
            varargin{1} = '';
            [~,getids,from,to,dt,method] = this.process_get_data_input(varargin);
            
            % right now I only know how to process the "day" feature
            if ~isfield(features,'day_of_week')
                error('~isfield(features,''day_of_week'')')
            end
            
            % generate output time vector
            X = repmat(DataProfile(),1,length(getids));
            for i = 1:length(getids)
                c = this.cluster_manager.get_cluster( ...
                    catstruct(features,struct('id',getids(i))));
                X(i) = DataProfile(getids(i),nan,this.cluster_manager.time,[],c.centroid,[]);
            end
            
            % restrict to requested times
            r_time = from:dt:to;
            switch method
                case 'interpolate'
                    Y = this.get_data_interpolate(X,r_time,getdays,getids);
                case 'snap'
                    Y = this.get_data_snap(X,r_time,nan,getids);
                otherwise
                    error('unknown method')
            end
            
        end
        
        function [dty_veh] = extract_density_veh(this,X)
            if size(X,1)>1
                error('size(X,1)>1')
            end
            dty_ids = [X.id];
            use_lengths = this.length_for_ids(Utils.index_into(dty_ids,this.ids));
            d = vertcat(X.dty_vpk);
            dty_veh = d.*repmat(use_lengths',1,size(d,2));
        end
        
        function [x] = get_ids_for_linkids(this,lnkids)
            x = cell(1,length(lnkids));
            for i=1:length(lnkids)
                x{i} = this.ids(lnkids(i)==this.link_ids);
            end
        end
        
        function [x] = get_linkids_for_ids(this,idss)
            x = this.link_ids(Utils.index_into(idss,this.ids));
        end
        
        function [internal_dp_ids,internal_link_ids]=get_internal_sensors(this,ni)
            % get all data provider ids that are internal to the network
            internal_link_ids = ni.get_internal_link_ids;
            dp_ids = this.get_ids_for_linkids(internal_link_ids);
            has_dp = cellfun(@(x) ~isempty(x),dp_ids);
            internal_dp_ids = Utils.cell2array(dp_ids(has_dp));     % internal sensors
            internal_link_ids = internal_link_ids(has_dp);          % corresponding links
        end
        
        function [eval_meas_ids,eval_links]=get_good_internal_sensors(this,ni,day)
            
            [~,internal_link_ids] = this.get_internal_sensors(ni);
            
            meas_ids = this.get_ids_for_linkids(internal_link_ids);
            all_meas_ids = [meas_ids{:}];
            all_meas_health = this.get_health(day,all_meas_ids);
            good_sensors = all_meas_ids(all_meas_health);
            link_has_good_sensor = cellfun( @(x) ~isempty(x) & all(ismember(x,good_sensors)) , meas_ids);
            eval_links = internal_link_ids(link_has_good_sensor);
            eval_meas_ids = meas_ids(link_has_good_sensor);
            
            % check that every link has at most one sensor
            if ~all(cellfun(@(x)length(x)==1,eval_meas_ids))
                error('~all(cellfun(@(x)length(x)==1,eval_meas_ids))')
            end
            
            % Can safely unpack
            eval_meas_ids = [eval_meas_ids{:}];
        end
        
    end
    
    methods( Abstract )
        [x] = get_health(this,getdays,getvds);
        [x] = uses_template(this,getdays,getvds);
        [this] = overwrite_force_bad_sensors(this,Z);
    end
    
    methods( Access = protected )
        
        function [getdays,getids,from,to,dt,method] = process_get_data_input(this,varargin)
            % days, ids, from, to, dt, method
            
            varargin = varargin{1};
            
            % deaults
            if numel(varargin)>0
                getdays = varargin{1};
            else
                getdays = 'all';
            end
            
            if numel(varargin)>1
                getids = varargin{2};
            else
                getids = 'all';
            end
            
            if numel(varargin)>2
                from = varargin{3};
            else
                from = 'start';
            end
            
            if numel(varargin)>3
                to = varargin{4};
            else
                to = from;
            end
            
            if numel(varargin)>4 && ~isnan(varargin{5})
                dt = varargin{5};
            else
                dt = this.measurement_dt;
            end
            
            if numel(varargin)>5
                method = varargin{6};
            else
                method = 'snap';
            end
            
            % interpret string inputs
            if ischar(getdays) && strcmp(getdays,'all')
                getdays = this.days;
            end
            
            if ischar(getids) && strcmp(getids,'all')
                getids = this.ids;
            end
            
            if ischar(from) && strcmp(from,'start')
                from = this.time(1);
            end
            
            if ischar(to) && strcmp(to,'start')
                to = this.time(1);
            end
            
            if ischar(to) && strcmp(to,'end')
                to = this.time(end);
            end
            
            % validate input
            if from<this.time(1)
                error('from<this.time(1)')
            end
            
            if to<from
                error('to<from')
            end
            
            if dt<=0
                error('dt<=0')
            end
            
            if mod(to-from,dt)~=0
                error('mod(to-from,dt)~=0')
            end
            
            % dt should be a whole divisor of measurement_dt
            if mod(this.measurement_dt,dt)~=0
                error('mod(this.measurement_dt,dt)~=0')
            end
            
            % from and to should be multiples of measurement_dt
            if mod(from,this.measurement_dt)~=0 || mod(to,this.measurement_dt)~=0
                error('mod(from,this.measurement_dt)~=0 || mod(to,this.measurement_dt)~=0')
            end
            
            if ~all(ismember(getids,this.ids))
                error('~all(ismember(getids,this.ids))')
            end
            
            if all(isnan(this.days))
                getdays = nan*getdays;
            end
            
%             if ~all(isnan(this.days)) && ~all(ismember(getdays,this.days))
%                 error('~isnan(this.days) && ~all(ismember(getdays,this.days))')
%             end
            
        end
        
        function [Y] = get_data_snap(this,X,r_time,getdays,getids)
            
            r_time_ind = arrayfun( @(x) find(this.time<=x,1,'last') , r_time);
            D = numel(getdays);
            I = numel(getids);
            Y = repmat(DataProfile(),D,I);
            
            n = length(X(1,1).flw_out_vph);
            flw_ind = r_time_ind(1:end-1);
            flw_ind(flw_ind>n)=n;
            for d=1:D
                for i=1:I
                    if numel(r_time)>1
                        r_flw_out_vph = X(d,i).flw_out_vph(flw_ind);
                        if ~isempty(X(d,i).flw_in_vph)
                            r_flw_in_vph = X(d,i).flw_in_vph(flw_ind);
                        else
                            r_flw_in_vph = r_flw_out_vph;
                        end
                    else
                        r_flw_in_vph = [];
                        r_flw_out_vph = [];
                    end
                    %                     if isfield(X(d,i),'dty_vpk')
                    if  ~isempty(X(d,i).dty_vpk)
                        r_dty_vpk = X(d,i).dty_vpk(r_time_ind);
                    else
                        r_dty_vpk = [];
                    end
                    Y(d,i) = DataProfile(getids(i),getdays(d),r_time,r_flw_in_vph,r_flw_out_vph,r_dty_vpk);
                end
            end
        end
        
        function [Y] = get_data_interpolate(this,X,r_time,getdays,getids)
            
            D = numel(getdays);
            I = numel(getids);
            
            Y = repmat(DataProfile(),D,I);
            for d=1:D
                for i=1:I
                    if numel(r_time)>1
                        if ~isempty(X(d,i).flw_in_vph)
                            r_flw_in_vph = interp1(this.time(1:end-1),X(d,i).flw_in_vph,r_time(1:end-1),'previous','extrap');
                        else
                            r_flw_in_vph = [];
                        end
                        r_flw_out_vph = interp1(this.time(1:end-1),X(d,i).flw_out_vph,r_time(1:end-1),'previous','extrap');
                    else
                        r_flw_in_vph = [];
                        r_flw_out_vph = [];
                    end
                    %                     if isfield(X(d,i),'dty_vpk')
                    if ~isempty(X(d,i).dty_vpk)
                        r_dty_vpk = interp1(this.time,X(d,i).dty_vpk,r_time,'linear','extrap');
                    else
                        r_dty_vpk = [];
                    end
                    Y(d,i) = DataProfile(getids(i),getdays(d),r_time,r_flw_in_vph,r_flw_out_vph,r_dty_vpk);
                end
            end
        end
        
    end
    
end

