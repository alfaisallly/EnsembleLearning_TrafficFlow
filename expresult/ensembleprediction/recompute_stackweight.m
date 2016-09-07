clear;clc;
% load pre-computed base prediction for each day
load expresult/ensembleprediction/309604_7am8pm.mat
% number of methods
num_basemethods = 5;
% number of previous days used to learn the stack weights
num_days_for_stack = 60;

% store the stack ensembled prediction for each day after the minimum
% number of days required to learn the stack weight
y_stackensemble  = zeros(size(allday_actual) - [num_days_for_stack,0,0]) ;
% store the stack coefficience
stackweight = zeros(num_basemethods , 4 , size(y_stackensemble,1) , length(predictionFrom) ) ;

% compute ensemble weights for each hour, for each day
for t = 1:length(predictionFrom)
    
    [ y_stackensemble(:,:,t) , stackweight(:,:,:,t) ] = ...
        StackEnsembleFlowPrediction.compute_ensemble_only( allday_actual(:,:,t), allday_base(:,:,:,t), num_basemethods , num_days_for_stack);
    
end





