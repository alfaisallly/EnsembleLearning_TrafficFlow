% Run the modelA/modelB experiment with TrafficPredictionSystemTester

clear
clc
close all

%% user input
experiment = 'complete_sensing';
filters = { @(x) numel(x.remove_ml_vdss)==0  & numel(x.remove_rp_vdss)==0 , ...
            @(x) numel(x.remove_ml_vdss)==15 & numel(x.remove_rp_vdss)==0  , ...
            @(x) numel(x.remove_ml_vdss)==0  & numel(x.remove_rp_vdss)==17 , ...
            @(x) numel(x.remove_ml_vdss)>=27 & numel(x.remove_rp_vdss)==34 };

%% parameters
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
input_folder = fullfile(root,'temp','input');
output_folder = fullfile(root,'temp','output');

ObjectFactory.delete_serialized_objects

model_based_params = struct( ...
        'name',                 'model_based' , ...
        'predictor_class',      'ModelFullPredictor' ,     ...
        'measurement_dt',       300 ,                      ...
        'model_runner',         'beats',                   ...
        'input_flow_predictor', 'historical_or_zoh_pems',  ...
        'state_estimator',      'linear_interpolation',    ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ,                    ...
        'armax_params',         RecursiveArmaxFlowPredictor.get_default_parameters );
    
%% run     
    
fid = fopen('runtimes.txt','w');

% historical or ZOH system
tic
A = TrafficPredictionSystemTester('210E',input_folder,output_folder,model_based_params,experiment);
A.mngr.add_filters(filters);
A.run_experiment('parallel');
A.process_results;
resultA = A.evaluate_hypothesis_test;
lambdaA = A.find_threshold_lambda;
fprintf(fid,'ZOH %f\n',toc);

% ARMAX system
tic
model_based_params.input_flow_predictor = 'recursiveARMAX';
B = TrafficPredictionSystemTester('210E',input_folder,output_folder,model_based_params,experiment);
B.mngr.add_filters(filters);
B.run_experiment('parallel');
B.process_results;
resultB = B.evaluate_hypothesis_test;
lambdaB = B.find_threshold_lambda;
fprintf(fid,'ARMAX %f\n',toc);

fclose(fid);

save(experiment)

%% get info



[rA,cA] = find(~isnan(lambdaA));
lA = nan(1,numel(rA));
mean_ml_A = nan(1,numel(rA));
mean_mb_A = nan(1,numel(rA));
for i=1:numel(rA)
    lA(i) = lambdaA(rA(i),cA(i));
    mean_ml_A(i) = resultA(rA(i),cA(i)).mean_ml;
    mean_mb_A(i) = resultA(rA(i),cA(i)).mean_mb;
end

[rB,cB] = find(~isnan(lambdaB));
mean_ml_B = nan(1,numel(rA));
mean_mb_B = nan(1,numel(rA));
for i=1:numel(rB)
    lB(i) = lambdaB(rB(i),cB(i));
    mean_ml_B(i) = resultB(rB(i),cB(i)).mean_ml;
    mean_mb_B(i) = resultB(rB(i),cB(i)).mean_mb;
end

[lA;lB]
mean_ml_A
mean_mb_A
mean_mb_B

xA = lambdaA;
xA(isnan(xA)) = 0;
xA(xA>0) = 1;
xA(xA<0) = -1;

figure
imagesc(xA,[-1 1])


xB = lambdaB;
xB(isnan(xB)) = 0;
xB(xB>0) = 1;
xB(xB<0) = -1;

figure
imagesc(xB,[-1 1])

