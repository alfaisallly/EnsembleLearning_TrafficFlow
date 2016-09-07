%=====================================================================
% a running example to make stack ensemble prediction for a single day
%=====================================================================

clear;
station_id =  310215;
arcadia_dp = ArcadiaDataProvider(station_id);

predictionFrom = [7] *3600;
predictionTo = [8]*3600;

dayID = 181 ; % must be greater then 180, the minimum number of past days used for training
getdays = arcadia_dp.days(dayID); % which day to test the prediction 

num_days_for_stack=60;% number of previous days use for training stack ensemble given base predictions
num_basemethods = 5;

for t = 1:length(predictionFrom)
    
    from = predictionFrom(t);
    to = predictionTo(t);
    actual=[];
    base = [];
    integrate_predictor = StackEnsembleFlowPrediction(arcadia_dp);
    [y_today_stackensemble , stackweight , base  , actual ] = integrate_predictor.predict(station_id,getdays,from,to,900 , actual , base  , num_basemethods , num_days_for_stack);
     
    
end


%save('example.mat', 'base', 'actual', 'predictionFrom', 'predictionTo' , 'num_days_for_stack', 'num_basemethods' , 'station_id' )









