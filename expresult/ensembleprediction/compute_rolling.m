%=====================================================================
% script to run ensemble prediction rollingly for all valid days in a year (roughly half of a year)
%=====================================================================

clear;clc;
station_id =  310215;
arcadia_dp = ArcadiaDataProvider(station_id);

% make one-hour-ahead predictions
predictionFrom = [7,8,9,10,11,12,13,14,15,16,17,18,19] *3600;
predictionTo = [8,9,10,11,12,13,14,15,16,17,18,19,20]*3600;

% need at least 180 days to start the prediction: 
% need previous 120 days to train each base method, need previous 60 days 
% to train the ensemble 
start_date_idx = 181 ;
num_days_for_stack=60;
% five base methods: 1.armax 2.pls 3.svm 4.gaussian process 5.kernel ridge
num_basemethods = 5 ;


allday_actual = zeros(length(arcadia_dp.days)-start_date_idx+1+num_days_for_stack,4,length(predictionFrom));
allday_base = zeros(length(arcadia_dp.days)-start_date_idx+1+num_days_for_stack,4,5,length(predictionFrom));
y_stackensemble = zeros(length(arcadia_dp.days)-start_date_idx+1,4,length(predictionFrom)) ; 

for t = 1:length(predictionFrom)
    
    actual=[];
    base = [];
    for d = start_date_idx:length(arcadia_dp.days)
        
        getdays = arcadia_dp.days(d);
        from = predictionFrom(t);
        to = predictionTo(t);
        
        integrate_predictor = StackEnsembleFlowPrediction(arcadia_dp);
        [y_stackensemble(d-start_date_idx+1,:,t), stackweight , base  , actual ] = integrate_predictor.predict(station_id,getdays,from,to,900 , actual , base  , num_basemethods , num_days_for_stack);
              
    end
    allday_actual(:,:,t) = actual;
    allday_base(:,:,:,t) = base;
    
end





