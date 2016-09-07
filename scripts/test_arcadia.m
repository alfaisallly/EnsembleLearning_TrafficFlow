clear
close all

[~,all_ids] = ArcadiaDataProvider.get_metadata();  
       
% -

% test ArcadiaDataProvider
adp = ArcadiaDataProvider(307925);


days_ids = adp.get_days_and_ids;

getdays = days_ids.days(308);
getids = days_ids.ids;

adp.get_health(getdays);

adp.get_health('all');

% x = adp.uses_template(getdays,getids)
from = 10*3600;
to = 11*3600;
x = adp.get_data(getdays,'all',from,to);



pls = PLSFlowPredictor(adp);

num_clusters = 3;
y = pls.predict(nan,getdays,from,to,num_clusters)

 