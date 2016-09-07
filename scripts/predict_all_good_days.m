% run 210E prediction for all good days

clearvars
close all
clc

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(fullfile(root,'src')))
config = Config.get('210E');

for d = 1:length(config.good_days)
   
    close all
    day = config.good_days(d);
    fprintf('Generating report for %s\n',datestr(day));
    
    P = ModelFullPredictor( struct(                        ...
        'measurement_dt',       300 ,                      ...
        'model_runner',         'beats',                   ...
        'input_flow_predictor', 'zoh_sim',                 ...
        'state_estimator',      'linear_interpolation',    ...
        'configname' ,          '210E' ,                   ...
        'day' ,                 day ,                      ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ) );

    % create error data provider  
    error_data_provider = P.get_error_data_provider('pems');
    name = sprintf('210E_%d',day);
    
    P.compute_prediction_error(error_data_provider,fullfile(Folder.reports,name))

    save(name,'P')
end