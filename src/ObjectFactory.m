classdef ObjectFactory
    
    methods(Static)

        %% model runner
        function [beats] = beats_model_runner(ni,sim_dt,output_dt)      
%             x = ObjectFactory.hash({ni.configfile,sim_dt,output_dt});
%             matfile = fullfile(Folder.objects,sprintf('beats_model_runner_%s',x));
%             if exist([matfile '.mat'],'file')
%                 load(matfile)
%             else
                beats = BeatsModelRunner(ni,sim_dt,output_dt);
%                 save(matfile,'beats')
%             end
        end
        
        %% network information 
        function [ni] = network_information(configfile)
            x = ObjectFactory.hash({configfile});
            matfile = fullfile(Folder.objects,sprintf('network_information_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                ni = NetworkInformation(fullfile(Folder.config,configfile));
                save(matfile,'ni')
            end
        end
        
        function [fwy] = freeway_information(configfile)
            x = ObjectFactory.hash({configfile});
            matfile = fullfile(Folder.objects,sprintf('freeway_information_%s',x));            
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                fwy = FreewayInformation(fullfile(Folder.config,configfile));
                save(matfile,'fwy')
            end
        end
        
        %% data provider
        function [pems_dp] = pems_data_provider(days,internal_ids,use_template,name,link_ids,length_for_ids_km,cluster_method)
            x = ObjectFactory.hash({days,internal_ids,use_template,name,link_ids,length_for_ids_km,cluster_method});
            matfile = fullfile(Folder.objects,sprintf('pems_data_provider_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                pems_dp = PeMSDataProvider(days,internal_ids,use_template,name,link_ids,length_for_ids_km,cluster_method);
                save(matfile,'pems_dp')
            end
        end
        
        function [sim_dp] = sim_data_provider(model_runner,end_time)
            x = ObjectFactory.hash({model_runner.network_information.configfile,end_time});
            matfile = fullfile(Folder.objects,sprintf('sim_data_provider_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                sim_dp = SimDataProvider(model_runner,end_time);
                save(matfile,'sim_dp')
            end
        end
        
        function [sim_vds_dp] = sim_vds_data_provider(model_runner,end_time)
            x = ObjectFactory.hash({model_runner.network_information.configfile,end_time});
            matfile = fullfile(Folder.objects,sprintf('sim_vds_data_provider_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                sim_vds_dp = SimVDSDataProvider(model_runner,end_time);
                save(matfile,'sim_vds_dp')
            end
        end        
                       
        %% flow predictors
        function [simulation_flow_predictor] = simulation_flow_predictor(sim_dp)
            x = ObjectFactory.hash({sim_dp.configfile,sim_dp.end_time});
            matfile = fullfile(Folder.objects,sprintf('simulation_flow_predictor_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                simulation_flow_predictor = SimulationFlowPredictor(sim_dp);
                save(matfile,'simulation_flow_predictor')
            end
        end
        
        function [zoh_predictor] = zoh_predictor_sim(sim_dp)
            x = ObjectFactory.hash({sim_dp.configfile,sim_dp.end_time});
            matfile = fullfile(Folder.objects,sprintf('zoh_sim_predictor_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                zoh_predictor = ZOHFlowPredictor(sim_dp);
                save(matfile,'zoh_predictor')
            end
        end
        
        function [zoh_predictor] = zoh_predictor_pems(pems_dp)
            x = ObjectFactory.hash({pems_dp.days,pems_dp.ids});
            matfile = fullfile(Folder.objects,sprintf('zoh_pems_predictor_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                zoh_predictor = ZOHFlowPredictor(pems_dp);
                save(matfile,'zoh_predictor')
            end
        end
        
        function [historical_predictor] = historical_predictor(pems_dp)
            x = ObjectFactory.hash({pems_dp.days,pems_dp.ids});
            matfile = fullfile(Folder.objects,sprintf('historical_predictor_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                historical_predictor = HistoricalFlowPredictor(pems_dp);
                save(matfile,'historical_predictor')
            end
            
            % set forced bad detectors
            historical_predictor.data_provider.overwrite_force_bad_sensors(pems_dp.vds_forced_bad);
        end
          
        function [historical_or_zoh_predictor] = historical_or_zoh_predictor(pems_dp,sim_dp)
            if nargin>1
                x = ObjectFactory.hash({pems_dp.days,pems_dp.ids,sim_dp.configfile,sim_dp.end_time});
                matfile = fullfile(Folder.objects,sprintf('historical_or_zoh_predictor_sim_%s',x));
            else
                x = ObjectFactory.hash({pems_dp.days,pems_dp.ids});
                matfile = fullfile(Folder.objects,sprintf('historical_or_zoh_predictor_pems_%s',x));
                sim_dp = nan;
            end
            
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                historical_or_zoh_predictor = HistoricalOrZOHFlowPredictor(pems_dp,sim_dp);
                save(matfile,'historical_or_zoh_predictor')
            end
            
            % set forced bad detectors
            historical_or_zoh_predictor.historical_predictor.data_provider.overwrite_force_bad_sensors(pems_dp.vds_forced_bad);
        end
        
        function [recursiveARMAX_predictor] = recursiveARMAX_predictor(pems_dp,params)
            x = ObjectFactory.hash({pems_dp.days,pems_dp.ids,params});
            matfile = fullfile(Folder.objects,sprintf('recursiveARMAX_predictor_pems_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                recursiveARMAX_predictor = RecursiveArmaxFlowPredictor(pems_dp,params);
                save(matfile,'recursiveARMAX_predictor')
            end
            
            % set forced bad detectors
            recursiveARMAX_predictor.data_provider.overwrite_force_bad_sensors(pems_dp.vds_forced_bad);
        end
        
        function [scaled_historical_predictor] = scaled_historical_predictor(sim_dp,pems_dp)
            x = ObjectFactory.hash({sim_dp.configfile,sim_dp.end_time,pems_dp.days,pems_dp.ids});
            matfile = fullfile(Folder.objects,sprintf('scaled_historical_predictor_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                scaled_historical_predictor = ScaledHistoricalFlowPredictor(sim_dp,pems_dp);
                save(matfile,'scaled_historical_predictor')
            end
            
            % set forced bad detectors
            scaled_historical_predictor.data_provider.overwrite_force_bad_sensors(pems_dp.vds_forced_bad);
        end
 
        %% state estimators
        function [simulation_state_estimator] = simulation_state_estimator(mr,sim_dp)
            x = ObjectFactory.hash({mr.network_information,mr.sim_dt,mr.out_dt,sim_dp.configfile,sim_dp.end_time});
            matfile = fullfile(Folder.objects,sprintf('simulation_state_estimator_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                simulation_state_estimator = SimulationStateEstimator(mr,sim_dp);
                save(matfile,'simulation_state_estimator')
            end
        end
        
        function [lininterp_state_estimator] = lininterp_state_estimator(dp,fwy_info)
            x = ObjectFactory.hash({class(dp),dp.days,dp.ids,dp.link_ids,dp.length_for_ids,fwy_info.configfile});
            matfile = fullfile(Folder.objects,sprintf('lininterp_state_estimator_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                lininterp_state_estimator = LinInterpFrewayStateEstimator(dp,fwy_info);
                save(matfile,'lininterp_state_estimator')
            end
            
            % set forced bad detectors
            if ismethod(lininterp_state_estimator.data_provider,'overwrite_force_bad_sensors')
                lininterp_state_estimator.data_provider.overwrite_force_bad_sensors(dp.vds_forced_bad);
            end
        end
        
        %% experimtents
        function [decimation_experiment] = sensor_decimation_experiment(confname,infolder,outfolder)
            x = ObjectFactory.hash({confname,infolder,outfolder});
            matfile = fullfile(Folder.objects,sprintf('decimation_experiment_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                decimation_experiment = SensorDecimationJobManager(confname,infolder,outfolder);
                save(matfile,'decimation_experiment')
            end
        end

        function [complete_sensing_experiment] = complete_sensing_experiment(confname,infolder,outfolder)
            x = ObjectFactory.hash({confname,infolder,outfolder});
            matfile = fullfile(Folder.objects,sprintf('complete_sensing_experiment_%s',x));
            if exist([matfile '.mat'],'file')
                load(matfile)
            else
                complete_sensing_experiment = CompleteSensingJobManager(confname,infolder,outfolder);
                save(matfile,'complete_sensing_experiment')
            end
        end
        
        %% auxiliary
        
        function []=delete_serialized_objects()
        % remove all objects that have been serialized to local .mat files.
            delete(fullfile(Folder.objects,'*.mat'))
        end
        
    end
    
    methods(Static, Access=private)

        function [x] = hash(p)
            x = DataHash(p);
        end
        
        
    end
end