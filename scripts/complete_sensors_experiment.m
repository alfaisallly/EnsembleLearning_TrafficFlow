%  Run complete_sensing_experiment for several variations of mbTPS

clear
clc
close all

ObjectFactory.delete_serialized_objects

%% parameters
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
input_folder = fullfile(root,'temp','input');
output_folder = fullfile(root,'temp','output');

sim_sim_params = struct( ...
        'predictor_class',      'ModelFullPredictor' ,     ...
        'measurement_dt',       300 ,                      ...
        'model_runner',         'beats',                   ...
        'input_flow_predictor', 'simulation',              ...
        'state_estimator',      'simulation',              ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ) ;
sim_sim_params.name = 'sim_sim';
    
model_less_params = struct( ...
        'predictor_class',      'ModelLessPredictor' ,     ...
        'prediction_horizon',   3600 ,                     ...
        'update_dt',            1800 ,                     ...
        'start_time',           0 ,                        ...
        'end_time',             86400 ) ;
model_less_params.name = 'model_less';

lin_sim_params = sim_sim_params;
lin_sim_params.state_estimator = 'linear_interpolation';
lin_sim_params.name = 'lin_sim';

sim_pems_params = sim_sim_params;
sim_pems_params.input_flow_predictor = 'historical_or_zoh_pems';
sim_pems_params.name = 'sim_pems';

lin_pems_params = sim_sim_params;
lin_pems_params.state_estimator = 'linear_interpolation';
lin_pems_params.input_flow_predictor = 'historical_or_zoh_pems';
lin_pems_params.name = 'lin_pems';

lin_armax_pems_params = sim_sim_params;
lin_armax_pems_params.state_estimator = 'linear_interpolation';
lin_armax_pems_params.input_flow_predictor = 'recursiveARMAX';
lin_armax_pems_params.name = 'lin_armax_pems';

all_params = {model_less_params,sim_sim_params,lin_sim_params,sim_pems_params,lin_pems_params,lin_armax_pems_params};
clear model_less_params sim_sim_params lin_sim_params sim_pems_params lin_pems_params lin_armax_pems

%% create the job manager
mngr = ObjectFactory.complete_sensing_experiment('210E',input_folder,output_folder);

%% execute jobs
results_files = cell(1,numel(all_params));
for i=numel(all_params):numel(all_params)%1:numel(all_params)
    results_files{i} = mngr.execute_job_series(all_params{i});
end

% load complete_sensing_experiment

%% process results
for i=numel(all_params):numel(all_params)
    i
    results = mngr.process_result(results_files{i});
    save(all_params{i}.name,'results')
end
% clear results

results.sixty_minute_prediction(11).flw_error=[];
results.sixty_minute_prediction(20).flw_error=[];
%% execute jobs, process results, run reports
for i=numel(all_params):numel(all_params)
%     load(all_params{i}.name)
    mngr.report_result(results,fullfile(here,all_params{i}.name));
end
clear results

% % compare boundary flows
% names = cellfun(@(x)x.name,all_params,'UniformOutput',false);
% Utils.compare_boundary_flows_two_results( mngr.days , ...
%      mngr.process_result(results_files{strcmp(names,'lin_sim')}) , ...
%     'Simulation' , ...
%     mngr.process_result(results_files{strcmp(names,'lin_pems')}) , ...
%     'PeMS' ,fullfile(here,'cmp_sim_pems'));
% % 
% %% comparisons
% for i=1:numel(all_params)
%     load(all_params{i}.name)
%     x=arrayfun(@(x)mean(x.flw_error,2),results.sixty_minute_prediction,'UniformOutput',false);
%     e(i,:) = Utils.meanwithnan(horzcat(x{:}),1);
% end
% figure
% t = 1:35;
% plot(1:35,e,'LineWidth',2)
% legend(regexprep(cellfun(@(x)x.name,all_params,'UniformOutput',false),'_','-'))
% vline(t(mngr.days==datenum(2014,10,13)),'--');
% grid
% 
% myline = repmat(struct('xdata',[],'ydata',[],'mean',nan,'stddev',nan),1,numel(all_params));
% for i=1:numel(all_params)
%     load(all_params{i}.name)
%     flw_error = PredictionJobManager.extract_flow_error(results.sixty_minute_prediction)
%     figure
%     [h,m,s]=Utils.hist(flw_error.as_vector*100,'error [%]','exponential',20,[0 150]);
%     myline(i).xdata=h(2).XData;
%     myline(i).ydata=h(2).YData;
%     myline(i).mean = m;  
%     myline(i).stddev = s;
% %     close
% end
% 
% figure
% for i=1:numel(myline)
%     plot(myline(i).xdata,myline(i).ydata,'LineWidth',2)
%     leg{i} = sprintf('%s (mean %.1f%%)',regexprep(all_params{i}.name,'_',' '),myline(i).mean);
%     hold on
% end
% legend(leg)
% grid
% xlabel('error [%]')
% ylabel('frequency')
% 
