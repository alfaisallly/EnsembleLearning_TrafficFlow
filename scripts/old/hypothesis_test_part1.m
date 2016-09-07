% Part 1 of the hypothesis testing experiment.
% run complete_sensing_experiment for model A, model B, and model_less
% save the result

clear
clc
close all

ObjectFactory.delete_serialized_objects

%% parameters
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
input_folder = fullfile(root,'temp','input');
output_folder = fullfile(root,'temp','output');

model_based_A_params = struct( ...
        'predictor_class',      'ModelFullPredictor' ,     ...
        'measurement_dt',       300 ,                      ...
        'model_runner',         'beats',                   ...
        'input_flow_predictor', 'historical_or_zoh_pems',  ...
        'state_estimator',      'linear_interpolation',    ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ) ;
model_based_A_params.name = 'model_based_A';

model_based_B_params = struct( ...
        'predictor_class',      'ModelFullPredictor' ,     ...
        'measurement_dt',       300 ,                      ...
        'model_runner',         'beats',                   ...
        'input_flow_predictor', 'recursiveARMAX',          ...
        'state_estimator',      'linear_interpolation',    ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ) ;
model_based_B_params.name = 'model_based_B';

model_less_params = struct( ...
        'predictor_class',      'ModelLessPredictor' ,     ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ) ;
model_less_params.name = 'model_less';

all_params = {model_based_B_params,model_based_A_params,model_less_params};

%% create the job manager
mngr = ObjectFactory.complete_sensing_experiment('210E',input_folder,output_folder);

%% execute jobs
results_files = cell(1,numel(all_params));
for i=1:numel(all_params)
    results_files{i} = mngr.execute_job_parallel_build(all_params{i});
end

save hypothesis_testing_experiment

%% process results
for i=1:numel(all_params)
    i
    results = mngr.process_result(results_files{i});
    save(all_params{i}.name,'results')
end
clear results

%% reports
for i=1:numel(all_params)
    load(all_params{i}.name)
    mngr.report_result(results,fullfile(here,all_params{i}.name));
end
clear results



