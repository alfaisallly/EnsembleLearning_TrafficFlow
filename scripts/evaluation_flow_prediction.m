clear
close all
clc
%% parameters
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
% input_folder = fullfile(root,'temp','input');
% output_folder = fullfile(root,'temp','output');

ObjectFactory.delete_serialized_objects

% Arcadia data provider
station_id = 309426;
arcadia_dp = ArcadiaDataProvider(station_id);
days_ids = arcadia_dp.get_days_and_ids;
getdays = days_ids.days(198);
datestr( getdays)
% arcadia_dp_noday = ArcadiaDataProvider(station_id,getdays);

% create ARMAX predictor
params = struct();
ARMAX_predictor = ObjectFactory.recursiveARMAX_predictor(arcadia_dp,params);

% create PLS predictor
PLS_predictor = PLSFlowPredictor(arcadia_dp);

%create zero order flow predictior
ZOH_predictor = ZOHFlowPredictor(arcadia_dp);

% run prediction
PredictionFrom = linspace(7,23,23-7+1) * 3600 ;
PredictionTo = linspace(8,24,24-8+1) * 3600 ;
%from = 7*3600;
%to = 8*3600;

% store prediction results
pls_result = zeros(length(PredictionFrom)*4,1);
armax_result = zeros(length(PredictionFrom)*4,1);
integrate_result = zeros(length(PredictionFrom)*4,1);
true_result = zeros(length(PredictionFrom)*4,1);
% store error 
pls_RMAE = zeros(length(PredictionFrom),1);
armax_RMAE = zeros(length(PredictionFrom),1);
integrate_RMAE = zeros(length(PredictionFrom),1);
% store integration weight
%weight = zeros(length(PredictionFrom),2);  % common coefficient
predictiondim = 4;
weight = zeros(2,predictiondim,length(PredictionFrom)); % each (to-from)/900 has a pair of coefficients
for t=1:length(PredictionFrom)
    
    from = PredictionFrom(t);
    to = PredictionTo(t);
    
    % PLS
    num_clusters = 1;
    pls_result((t-1)*4+1 : t*4) = PLS_predictor.predict(station_id,getdays,from,to,num_clusters);
    
    %ARMAX
    % need to remove getdays from armax predictor
    armax_flw = ARMAX_predictor.predict(station_id,getdays,from,to,900);
    armax_result((t-1)*4+1 : t*4) = armax_flw.flw_out_vph; 
    
    
    % integrate prediction 
    weight(:,:,t) = IntegratePLS_ARMAX.ensembleweight_multicoefficient(arcadia_dp,station_id,getdays,from,to,1,900) ;
    
    %integrate_result((t-1)*4+1 : t*4) = [pls_result((t-1)*4+1 : t*4),armax_result((t-1)*4+1 : t*4)] * weight(t,:)';
    integrate_result((t-1)*4+1 : t*4) = pls_result((t-1)*4+1 : t*4) .* weight(1,:,t)'  +  armax_result((t-1)*4+1 : t*4) .* weight(2,:,t)' ;

    % compute error
    meas = arcadia_dp.get_data(getdays,station_id,from,to,900,'interpolate');
    true_result((t-1)*4+1 : t*4)  = meas.flw_out_vph; 
    pls_RMAE(t) = norm( pls_result((t-1)*4+1 : t*4) - true_result((t-1)*4+1 : t*4) ,1 )/norm(true_result((t-1)*4+1 : t*4) ,1);
    armax_RMAE(t) = norm( armax_result((t-1)*4+1 : t*4) - true_result((t-1)*4+1 : t*4)  ,1)/norm(true_result((t-1)*4+1 : t*4),1);
    integrate_RMAE(t) = norm( integrate_result((t-1)*4+1 : t*4) - true_result((t-1)*4+1 : t*4)  ,1)/norm(true_result((t-1)*4+1 : t*4),1);
end
%=====================================================================
% time series plot prediction vs actual
figure(1)
P_ts=plot(linspace(1,length(pls_result),length(pls_result)), pls_result ,'*-', ...
     linspace(1,length(pls_result),length(pls_result)), armax_result ,'o-', ... 
     linspace(1,length(integrate_result),length(integrate_result)), integrate_result ,'s-', ... 
     linspace(1,length(pls_result),length(pls_result)), true_result ,'d-'); 
set(P_ts,'LineWidth',2);
LG=legend('PLS','ARMAX','integrate','actual');
%====================================================================
% compute deviation
figure(2)
pls_deviation = pls_result- true_result ;
armax_deviation = armax_result - true_result;
integrate_deviation = integrate_result - true_result;
P_deviation=plot(linspace(1,length(pls_result),length(pls_result)), pls_deviation,'*-', ...
     linspace(1,length(pls_result),length(pls_result)), armax_deviation ,'o-', ... 
     linspace(1,length(integrate_result),length(integrate_result)), integrate_deviation ,'o-', ... 
     linspace(1,length(pls_result),length(pls_result)), zeros(length(pls_result),1) ,'d-'); 
set(P_deviation,'LineWidth',2);
LG=legend('PLS','ARMAX','integrate','actual');
%====================================================================
% compute relative MAE
figure(3)
P_RMAE = plot(linspace(1,length(PredictionFrom),length(PredictionFrom)), pls_RMAE ,'*-',...
              linspace(1,length(PredictionFrom),length(PredictionFrom)), armax_RMAE , 'o-',...
              linspace(1,length(PredictionFrom),length(PredictionFrom)), integrate_RMAE , 'o-',... 
              linspace(1,length(PredictionFrom),length(PredictionFrom)), zeros(length(PredictionFrom),1),'d-' );
set(P_RMAE,'LineWidth',2);
LG=legend('PLS','ARMAX','integrate','actual');
