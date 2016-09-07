classdef ArcadiaDataProvider < DataProvider
    
    properties
    end
    
    properties(Constant)
        BIN_SIZE = 3;       %  3 5-minute intervals
        ONLINE_THRESHOLD = 0.6;     % percentage of ON-LINE needed for use
        MEASUREMENT_DT = 900;           % 15 minute data (aggregated from 5 minutes)
    end
    
    methods( Access = public )
        
        function [this] = ArcadiaDataProvider(id,remove_days)
            
            if nargin<2
                remove_days = [];
            end
               
            if numel(id)>1
                error('works for only 1 detector at a time')
            end
            
            [~,all_ids] = ArcadiaDataProvider.get_metadata();
            if ~all(ismember(id,all_ids))
                error('asdf')
            end
            
            all_days = ArcadiaDataProvider.get_days_for_id(id);
            
            if any(~ismember(remove_days,all_days))
                error('remove day not found in all_days')
            end
            
            % remove days fropm list of all days
            all_days(ismember(all_days,remove_days)) = [];
            
            this = this@DataProvider(all_days,id,[],[],ArcadiaDataProvider.MEASUREMENT_DT);
            
            % load the data
            filename = fullfile(Folder.data,'Arcadia_processed',num2str(this.ids));
            if ~exist([filename '.mat'],'file')
                error('file not found')
            end
            load(filename)
            
            for i=1:length(this.days)
                day = this.days(i);                
                today_data = data([data.day]==day);
                [Bin_X_volume,Bin_DateTime] = ArcadiaDataProvider.VolumeAveragingOneDay(today_data.volume, 0:1/288:(1-1/288), ArcadiaDataProvider.BIN_SIZE);            
                this.data(i).day = today_data.day;
                this.data(i).time = Bin_DateTime;
                this.data(i).flw_out_vph = Bin_X_volume;
            end
            
            % call day filter.
            % 1) remove days with missing records
            % 2) if # positive volumes / day < threshold * size(Bin_X_volume,2)
            
            Bin_X_volume = vertcat(this.data.flw_out_vph);
            
            threshold = 0.6;
            good_days = ArcadiaDataProvider.DayFilter(Bin_X_volume , threshold);
            clear Bin_X_volume  day_ind
            
            this.data(~good_days) = [];
            this.time = 0:ArcadiaDataProvider.MEASUREMENT_DT:86400;
            this.days(~good_days) = [];
            
            this.vds_is_good = true(1,length(this.days));
            
            % create a cluster manager
            cluster_method = 'day_of_week';
            this.cluster_manager = ClusterManager( ...
                this.ids,this.days,this.time,this.data,cluster_method);
        end
        
        function [x] = get_health(this,getdays,~)
            % [x] = get_health(this,getdays,getvds)
            
            if isnan(getdays) & numel(this.days)>1
                error('isnan(getdays) & numel(this.days)>1')
            end
            
            if ischar(getdays) && strcmp(getdays,'all')
                getdays = this.days;
            end
            
            if isnan(getdays)
                days_ind = 1;
            else
                days_ind = Utils.index_into(getdays,this.days);
            end
            
            x = this.vds_is_good(days_ind);
            
        end
        
        function [x] = uses_template(this,getdays,getvds)
            x = false;
        end
        
        function [x] = get_data(this,varargin)
            %  days,ids,from,to,dt,method (interpolate or snap)
            x = this.get_data@DataProvider(varargin{:});
        end
        
        function [this] = overwrite_force_bad_sensors(this,Z)
        end
            
    end
    
    methods( Static, Access = public )
        
        function [days,ids] = get_metadata()
            days = [];
            ids = [];
            load(fullfile(Folder.data,'Arcadia_processed','metadata'))
        end
        
        function [days] = get_days_for_id(id)
            
            filename = fullfile(Folder.data,'Arcadia_processed',num2str(id));
            if ~exist([filename '.mat'],'file')
                error('file not found')
            end
            load(filename)
            days = [data.day];
        end
        
    end
    
    methods( Static, Access = private )
        
        function [Bin_X_volume,Bin_DateTime] = VolumeAveragingOneDay(X_volume, DateTime, BinSize)
            % smoothing the volume data by taking the average in each bin
            % input:
            
            % X_volume: matrix of raw volume data, may contain NaN if the detector
            % fails. Number of rows in X_volume is the number of days we have in the
            % raw data. Number of columns in X_volume is 288 (12*24)
            
            % DateTime: The Date/Time of detection correspond to entries in X_volume, in matlab
            % datenum format.
            
            % BinSize: how many 5 mins interval in each bin, has to be divisable by 12
            
            % output:
            % Bin_X_volume: the volume (#cars/hour) averaged in BinSize * 5 minutes.
            % NaN if a data is not available in the corresponding interval.
            % Bin_DateTime: the corresponding date/time in Bin_X_volume.
            
            if rem(288,BinSize) ~= 0
                error('BinSize must be divisible by 288(24*12)');
            end
            
            Bin_X_volume = NaN(1,288/BinSize);
            Bin_DateTime = NaN(1,288/BinSize);
            
            for i=1:BinSize:288
                Validvolume = 0;
                Sumvolume = 0;
                for k=0:BinSize-1
                    if isnan(X_volume(i+k))~=1
                        Sumvolume = Sumvolume + X_volume(i+k);
                        Validvolume = Validvolume + 1;
                        Bin_DateTime(floor(i/BinSize)+1)=DateTime(i+k);
                    end
                end
                Bin_X_volume(floor(i/BinSize)+1)=Sumvolume/Validvolume;
            end
                        
        end
        
        function [Bin_X_volume,Bin_DateTime] = VolumeAveraging(X_volume, DateTime, BinSize)
            % smoothing the volume data by taking the average in each bin
            % input:
            
            % X_volume: matrix of raw volume data, may contain NaN if the detector
            % fails. Number of rows in X_volume is the number of days we have in the
            % raw data. Number of columns in X_volume is 288 (12*24)
            
            % DateTime: The Date/Time of detection correspond to entries in X_volume, in matlab
            % datenum format.
            
            % BinSize: how many 5 mins interval in each bin, has to be divisable by 12
            
            % output:
            % Bin_X_volume: the volume (#cars/hour) averaged in BinSize * 5 minutes.
            % NaN if a data is not available in the corresponding interval.
            % Bin_DateTime: the corresponding date/time in Bin_X_volume.
            
            if rem(288,BinSize) ~= 0
                error('BinSize must be divisible by 288(24*12)');
            end
            
            Bin_X_volume = NaN(length(X_volume(:,1)),288/BinSize);
            Bin_DateTime = NaN(length(X_volume(:,1)),288/BinSize);
            
            for day = 1:length(X_volume(:,1))
                for i=1:BinSize:288
                    Validvolume = 0;
                    Sumvolume = 0;
                    for k=0:BinSize-1
                        if isnan(X_volume(day,i+k))~=1
                            Sumvolume = Sumvolume + X_volume(day,i+k);
                            Validvolume = Validvolume + 1;
                            Bin_DateTime(day,floor(i/BinSize)+1)=DateTime(day,i+k);
                        end
                    end
                    Bin_X_volume(day,floor(i/BinSize)+1)=Sumvolume/Validvolume;
                end
            end
            clear Sumvolume Validvolume i k day
            
        end
        
        function [good_days] = DayFilter(Bin_X_volume , threshold)
            % filtering bad data according to the following rules:
            % 1. if there are missing data (i.e. NaN) after averaged in day, delete the day
            % 2. if most of the volume records in a day is 0 (nnz < threshold * number
            % of records in a day), delete the day
            
            % input :
            
            % Bin_X_volume: the volume (#cars/hour) after smoothing. Number of rows in Bin_X_volume is the number of days we have in the
            % raw data. Number of columns is the number of intervals after smoothing. NaN if a data is not available in the corresponding interval.
            
            % Bin_DateTime: the corresponding date/time in Bin_X_volume.
            
            % threshold: if the number of nonzero volumes in a row (a day) is less
            % than threshold * (total number of intervals/day ), than delete the record
            % of that day
            
            % output:
            
            % Filtered_Bin_X: matrix of volume after filtering the bad data. each row
            % is a day
            
            % Filtered_DateTime: the corresponding date/time for entries in Filtered_Bin_X
            missingvalues=sum(isnan(Bin_X_volume),2);
            numdays = size(Bin_X_volume,1);
            good_days=true(numdays,1);
            
            % if the fraction of time intervals with positive volume is below
            % threshold, then delete the day
            for day=1:numdays
                good_days(day) = missingvalues(day) == 0  && nnz(Bin_X_volume(day,:)) > threshold * length(Bin_X_volume(1,:));
            end
                        
        end
        
    end
    
end

