clear
close all

configname = '210E';
config = Config.get(configname);
ni = ObjectFactory.network_information(config.xml_file);
cluster_method = '';
dp = Utils.get_pems_dp(ni,configname,cluster_method);

        
keys = dp.cluster_manager.get_keys;


cluster = dp.cluster_manager.get_cluster(keys(1))
