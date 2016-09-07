function []=parse_Arcadia_data()
% detector status: 0 - ON_LINE
%                  1 - COMM_ERROR/ON_LINE
%                  2 - COMM_ERROR
%                  3 - FAILURE
%                  4 - COMM_ERROR/FAILURE/ON_LINE
%                  nan - any other unknown status

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_folder = fullfile(root,'data','Arcadia');

d = dir(data_folder);
d([d.isdir])=[];
X = cellfun(@(x)sscanf(x,'DetectorArchive-%d-%d-%d-%d.csv'),{d.name},'UniformOutput',false);
Xh = horzcat(X{:});
clear X

all_ids = Xh(1,:);
all_years = Xh(2,:);
all_months = Xh(3,:);
all_days = Xh(4,:);
unique_ids = unique(Xh(1,:));
clear Xh

% save metada
ids = unique_ids;
days = sort(unique(datenum([all_years;all_months;all_days]')));
save(fullfile(root,'data','Arcadia_processed','metadata'),'ids','days');
clear ids days

day_struct = struct('day',nan,'volume',nan(288,1),'status',nan(288,1));

for i = 1 : length(unique_ids)
    
    id = unique_ids(i);
    display(id)
    
    % find days for this id
    ind = all_ids==id;
    id_days = sort(datenum([all_years(ind);all_months(ind);all_days(ind)]'));

    data = repmat(day_struct,1,0);
    
    for d=1:length(id_days)
        
        day = id_days(d);
        
        filename = fullfile(data_folder,sprintf('DetectorArchive-%d-%s-%0d-%0d.csv', ...
            id, ...
            datestr(day,'yyyy'), ...
            str2double(datestr(day,'mm')), ...
            str2double(datestr(day,'dd'))));
        
        if exist(filename,'file')~=2
            continue
        end
            
        fID = fopen(filename);
        [rawdata]=textscan(fID, '%d %s %s %f %f %f %f %f %f %f %f %f %f',...
            'headerlines', 1,...
            'delimiter',',',...
            'TreatAsEmpty','NA',...
            'EmptyValue', NaN);
        fclose(fID);

        more_data = FormatRawData(rawdata,day_struct); %run script FormatData
        data(end+1:end+length(more_data)) = more_data;
            
    end
    
    % save loaded data as mat file
    save(fullfile(root,'data','Arcadia_processed',num2str(id)),'data');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out_data]=FormatRawData(rawdata,day_struct)

% convert date_time into matlab numeric datenum representation
% e.g: in the raw data "Sun Nov 08 00:00:00 PST 2015"
%                      "Sun Nov 08 00:05:00 PST 2015"

% X_volume is a 288-by-days matrix contains the volume measurements
% don't know how many days with full measurements yet
% expand X_volume while counting days with full measurements
% cell in X_volume with NaN means no record detected in that time

volume = rawdata{1,4};

% transform status
status = rawdata{1,3};
numerical_status = nan(size(status,1),1);
numerical_status(strcmpi(status,'ON_LINE')) = 0;
numerical_status(strcmpi(status,'COMM_ERROR/ON_LINE')) = 1;
numerical_status(strcmpi(status,'COMM_ERROR')) = 2;
numerical_status(strcmpi(status,'FAILURE')) = 3;
numerical_status(strcmpi(status,'COMM_ERROR/FAILURE/ON_LINE')) = 4;
clear status

% transform date/time
date_time = strrep(rawdata{1,2},'PST','');
date_time = strrep(date_time,'PDT','');
date_time = datenum(date_time,'"ddd mmm dd HH:MM:SS  yyyy"');

clear rawdata

days = floor(date_time);
unique_days = unique(days);
num_days = length(unique_days);

out_data = repmat(day_struct,1,num_days);
for d=1:num_days
    
   day = unique_days(d);
   ind = days == day;
   day_volume = volume(ind);
   day_time = date_time(ind);
   day_status = numerical_status(ind);
   
   idx = Utils.index_into(round((day_time-day)*288),0:287);    
   out_data(d).day = day;
   out_data(d).volume(idx) = day_volume;
   out_data(d).status(idx) = day_status;

end
