clear
close all

%==================================================================================
% experiment parameter
PredictionFrom = linspace(7,9,9-7+1) * 3600 ;
PredictionTo = linspace(8,10,10-8+1) * 3600 ;
%PredictionFrom = [7*3600];
%PredictionTo = [8*3600];
num_testdetectors = 10;
num_clusters = 1;
%% parameters
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
% input_folder = fullfile(root,'temp','input');
% output_folder = fullfile(root,'temp','output');

ObjectFactory.delete_serialized_objects

%==================================================================================
%get the list of detector IDs
d = dir([root '/data/Arcadia_processed/']);
alldetectorID = zeros(length(d)-3,1) ;
for i = 1:length(alldetectorID)-1
    alldetectorID(i)=str2num(d(i+2).name(1:6));
end
alldetectorID(end)=9901;
clear d i
%===================================================================================
%testdetectorID  = alldetectorID(randperm(length(alldetectorID),num_testdetectors));
testdetectorID = 309426;


%for i = 1:length(testdetectorID)

%Arcadia data provider
station_id = testdetectorID(1);
arcadia_dp = ArcadiaDataProvider(station_id);
days_ids = arcadia_dp.get_days_and_ids;
if length(days_ids.days) < 181
    error('this detector does not have suffient number of sample days')
end

numdays = length(days_ids.days);
%===================================================================================
pls_result = zeros(4,24,numdays-180);
armax_result = zeros(4,24,numdays-180);
svm_result = zeros(4,24,numdays-180);
integrate_result = zeros(4,24,numdays-180);
true_result = zeros(4,24,numdays-180);
%===================================================================================

for t = 1:length(PredictionFrom)
    
    from = PredictionFrom(t);
    to = PredictionTo(t);
    
    %the first half of the year is for training and validation
    for d = 181:numdays
        
        
        getdays = days_ids.days(d);
        datestr(getdays)
        [training_dp] = construct_trainingdp_includetoday(arcadia_dp,getdays,180);
        
        % obtain ground truth
        meas = training_dp.get_data(getdays,station_id,from,to,900,'interpolate')
        
        
        if norm(meas.flw_out_vph,1) > 0
            
            % PLS
            PLS_predictor = PLSFlowPredictor(training_dp);
            y_pls =  PLS_predictor.predict(station_id,getdays,from,to,num_clusters)
            pls_result(:,PredictionTo(t)/3600,d-180) = y_pls' ;
            
            % ARMAX
            params = struct();
            ARMAX_predictor = ObjectFactory.recursiveARMAX_predictor(training_dp,params);
            y_armax = ARMAX_predictor.predict(station_id,getdays,from,to,900);
            armax_result(:,PredictionTo(t)/3600,d-180) = y_armax.flw_out_vph';
            
            % SVM
            SVM_predictor = SVMFlowPredictor(training_dp);
            y_svm = SVM_predictor.predict(station_id,getdays,from,to,900);
            
            
            
            %             % combine
            %             weight = IntegratePLS_ARMAX.ensembleweight_multicoefficient(arcadia_dp,station_id,getdays,from,to,1,900) ;
            %             y_integrate = y_pls' .* weight(1,:)'  +  y_armax.flw_out_vph' .* weight(2,:)'
            %             integrate_result(:,PredictionTo(t)/3600,d-180)  = y_integrate ;
            
            % true flow
            meas = training_dp.get_data(getdays,station_id,from,to,900,'interpolate');
            true_result(:,PredictionTo(t)/3600,d-180) = meas.flw_out_vph ;
        end
        
        
    end
    
end
%end
%
% cd expresult/exp2/
% save('7to8.mat','PredictionFrom' ,'PredictionTo' ,'error_armax' ,'error_cpls' ,'num_clusters','testdetectorID');
% cd ..