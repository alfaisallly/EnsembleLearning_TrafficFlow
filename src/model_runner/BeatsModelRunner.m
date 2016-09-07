classdef BeatsModelRunner < ModelRunner

    properties
        beats   @BeatsSimulation
    end
    
    properties(Constant)
        vehicle_type_id = 1
    end
    
    % override
    methods (Access = public)

        function [this] = BeatsModelRunner(network_information,sim_dt,output_dt,output_prefix)
            if nargin<5
                param = struct('SIM_DT',sim_dt,'OUTPUT_DT',output_dt);
            else
                param = struct('SIM_DT',sim_dt,'OUTPUT_DT',output_dt,'OUTPUT_PREFIX',output_prefix);
            end
            param.NODE_FLOW_SOLVER =  'tampere';
            this = this@ModelRunner(network_information,sim_dt,output_dt);
            BeatsSimulation.import_beats_classes()
            this.beats = BeatsSimulation;
            this.beats.set_scenario(network_information.scenario);
            this.beats.create_beats_object(param);
            this.beats.intialize_beats_persistent(param);
        end
        
        function [this] = set_state_veh(this,given_linkids,dty_veh)  
            beats_linkids = this.beats.scenario_ptr.get_link_ids;
            if ~all(ismember(given_linkids,beats_linkids)) | ~all(ismember(beats_linkids,given_linkids))
               error('ids dont match')
            end
            ind = Utils.index_into(beats_linkids,given_linkids);
            
            % make sure none are negative
            dty_veh(dty_veh<0) = 0;
            
            this.beats.beats.set.totalDensity(Utils.column_vector(dty_veh(ind)));
        end
        
        function [linkid_dtyveh] = get_state_veh(this)
            link_ids = this.beats.scenario_ptr.get_link_ids();
            dty_in_veh = this.beats.beats.get.totalDensity(NaN);
            linkid_dtyveh = [link_ids' dty_in_veh];
        end
        
        function [this] = set_demand(this,link_id,dt,flw_vph)
            this.beats.beats.set.demand_for_link_si(link_id,dt,flw_vph/3600);
        end
        
        function [start_time,dt,flw_vph] = get_demand(this,link_id)
            dp = this.beats.beats.get.current_demand_for_link(link_id);
            if dp.getDemand.size()~=1
                error('only works for single vehicle type')
            end
            start_time = dp.getStartTime;
            dt = dp.getDt;
            flw_vph = 3600*eval(sprintf('[%s]',char(dp.getDemand.get(0).getContent)));            
        end
        
        function [this] = set_split_ratio(this,node_id,start_time,dt,sr_profiles)
            jaxb = edu.berkeley.path.beats.simulator.JaxbObjectFactory;
            for i=1:length(sr_profiles)
                sr(i) = jaxb.createSplitratio;
                sr(i).setLinkIn(sr_profiles(i).link_in)
                sr(i).setLinkOut(sr_profiles(i).link_out)
                sr(i).setVehicleTypeId(this.vehicle_type_id)
                sr(i).setContent(writecommaformat(sr_profiles(i).profile))
            end            
            this.beats.beats.set.splits_for_node(node_id,start_time,dt,sr);
        end
        
        function [rsrp] = get_split_ratio(this,node_id)
            srp = this.beats.beats.get.splitprofile_for_node(node_id);
            rsrp.node_id = srp.getNodeId;
            rsrp.start_time = srp.getStartTime;
            rsrp.dt = srp.getDt;
            splitratios = srp.getSplitratio();
            rsrp.splitratios = repmat(struct('link_in',nan,'link_out',nan,'profile',[]),1,splitratios.size());
            for i=1:splitratios.size()
                z = splitratios.get(i-1);
                rsrp.splitratios(i).link_in = z.getLinkIn;
                rsrp.splitratios(i).link_out = z.getLinkOut;
                x =  textscan(char(z.getContent),'%f','delimiter',',');
                rsrp.splitratios(i).profile = x{:}';
            end
        end
        
        function [this] = reset(this)
            this.beats.reset_simulation();
        end
        
        function [this] = run_simulation(this,duration)
            this.beats.run_beats_persistent(duration);
        end 
        
        function [that] = clone(this,whoclones)
            fprintf('BeatsModelRunner cloned by %s\n',whoclones)
            that = BeatsModelRunner(this.network_information,this.sim_dt,this.out_dt); 
        end
        
        function [link_ids] = get_link_ids(this)
            link_ids = Utils.row_vector(this.beats.scenario_ptr.get_link_ids);
        end
        
        function [data] = get_data_for_link_id(this,link_id)
            data = this.beats.get_output_for_link_id(link_id);
        end
        
        function [n] = get_num_time(this)
            if this.beats.sim_output_loaded
                n = size(this.beats.density_veh{1},1);
            else 
                n = 0;
            end
        end
        
        function [this] = set_output_dt(this,dt)
            this.out_dt = dt;
            this.beats.out_dt = dt;
            this.beats.beats.set.output_dt(dt);
        end
        
        function [vds_link] = get_vds2linkid(this)
            vds_sensor = this.beats.scenario_ptr.get_sensor_vds2id_map;
            sensor_link = this.beats.scenario_ptr.get_sensor_link_map;
            ind = Utils.index_into(vds_sensor(:,2),sensor_link(:,1));
            vds_link = [vds_sensor(:,1) sensor_link(ind,2)];
        end
          
    end

end

