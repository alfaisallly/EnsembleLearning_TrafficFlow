% Compute and save historical splits for a given configuration file.
% This is done 'off line' because it requires that the topology of sensors
% be recalculted based on the loop health for each day.

clearvars
close all
clc

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(fullfile(root,'src')))

name = '210E';
output_dt = 300;
end_time = 86400;

config = Config.get(name);
ni = ObjectFactory.network_information(config.xml_file);
mr = ObjectFactory.beats_model_runner(ni,config.sim_dt,output_dt);
sim_dp = ObjectFactory.sim_data_provider(mr,end_time);

vdss = ni.get_all_vds;
link_ids = ni.get_link_ids_for_vds(vdss);
link_length = ni.get_lengths_for_linkids_in_km(link_ids);
pems_dp = ObjectFactory.pems_data_provider(config.good_days,vdss,config.vds_use_template,name,link_ids,link_length);
flow_predictor = ObjectFactory.historical_predictor(pems_dp);
                    
% make a split predictor so that I can use dp_info
split_pred = SplitPredictor(ni,flow_predictor);

% go through days 
from = 0;
to = 86400;
dt = 300;
numDays = length(config.good_days);
numIds = length(split_pred.ids);

historical.node_ids = split_pred.ids;
historical.days = config.good_days;
historical.dp_info = cell(numDays,numIds);
historical.split_ratios = cell(1,numIds);

for i=1:numIds
    nOut = length(split_pred.dp_info(i).out_struct);
    historical.split_ratios{i} = nan(numDays,288,nOut);
end

for d=1:numDays
    day = config.good_days(d);
    split_pred.set_prediction_day(ni,day);
    all_dp_flw = split_pred.flow_predictor.predict_flw_for_dpids(split_pred.all_dps,day,from,to,dt);    
    for i=1:numIds
        historical.dp_info{d,i} = split_pred.dp_info(i);
        sr = split_pred.compute_split_from_dp_info(split_pred.dp_info(i),all_dp_flw);
        for j = 1:size(sr,1)
            historical.split_ratios{i}(d,:,j) = sr(j,:);
        end
    end    
end
   
save(fullfile(Folder.data,sprintf('%s_historical_splits',name)),'historical')
