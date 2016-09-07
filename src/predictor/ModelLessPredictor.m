classdef ModelLessPredictor < Predictor
    
    properties
        r
        good_days
        data_provider   @PeMSDataProvider
    end
    
    methods(Access=public)
        
        function [this] = ModelLessPredictor(params)
            this = this@Predictor(params); 
        end
        
        function [this]=run(this)
                        
            % get mainline vdss and corresponding links
            link_ids = this.ni.link_ids;
            int_link_ids = this.ni.get_internal_link_ids;
            vdss = this.ni.get_vds_for_link_ids(int_link_ids);
            ind = isnan(vdss);
            link_ind = Utils.index_into(int_link_ids(~ind),link_ids);
            vdss(ind) = [];
            clear ind int_link_ids
            
            pred_window = (0:this.measurement_dt:this.prediction_horizon);
            num_vdss = length(vdss);

            % get health information
            health=this.data_provider.get_health(this.day,'all');
            idx=Utils.index_into(vdss,this.data_provider.ids);
            this.r=health(idx);
            clear health idx
            
            % historical data
            X = this.data_provider.get_data(this.good_days,vdss,this.start_time,this.end_time+this.prediction_horizon,nan,'snap');
            
            % average historical data
            hist.time = X(1,1).time;
            hist.dty_vpk = nan(num_vdss,length(hist.time));
            hist.flw_vph = nan(num_vdss,length(hist.time)-1);
            for v = 1:num_vdss
                hist.dty_vpk(v,:) = Utils.meanwithnan(vertcat(X(:,v).dty_vpk),1);
                hist.flw_vph(v,:) = Utils.meanwithnan(vertcat(X(:,v).flw_out_vph),1);
            end
            
            % real time data
            realtime.time = X(1,1).time;
            realtime.dty_vpk = vertcat(X(this.good_days==this.day,:).dty_vpk);
            realtime.flw_vph = vertcat(X(this.good_days==this.day,:).flw_out_vph);
            
            clear X
            
            for p = 1:this.num_predictions
                
                update_time = this.start_time + (p-1)*this.update_dt;
                prediction_times = update_time + pred_window;
                pred_ind = Utils.index_into(prediction_times,hist.time);
                now_ind = pred_ind(1);
                num_times = length(prediction_times);
                ones_row = ones(1,num_times);
                
                % predictions
                pred_dty_vpk = ...
                    (this.r'*ones_row).* repmat(realtime.dty_vpk(:,now_ind),1,num_times) + ...
                    ((1-this.r)'*ones_row) .* hist.dty_vpk(:,pred_ind);
                
                pred_flw_vph = ...
                    (this.r'*ones_row).* repmat(realtime.flw_vph(:,now_ind),1,num_times) + ...
                    ((1-this.r)'*ones_row) .* hist.flw_vph(:,pred_ind);
                pred_flw_vph = pred_flw_vph(:,2:end);
                
                % log
                this.log(p).start_time = update_time;
                this.log(p).predicted_state = struct( ...
                    'flw_vph',nan(num_times-1,numel(link_ids)), ...
                    'dty_vpm',nan(num_times,numel(link_ids)));
                this.log(p).predicted_state.flw_vph(:,link_ind) = pred_flw_vph';
                this.log(p).predicted_state.dty_vpm(:,link_ind) = 1.60934*pred_dty_vpk';
                
            end
            
        end
        
    end
    
    methods(Access=protected)
        
        function [this]=set_values(this,params)
            
            this.set_values@Predictor(params);
            
            % load configuration information
            config = Config.get(params.configname);
            
            this.good_days = config.good_days;
        end
        
    end
    
    methods(Access=protected)
        
        function [this]=instantiate_objects(this,params)
            
            this.instantiate_objects@Predictor(params);
            
            config = Config.get(params.configname);
            
            % create data provider
			this.data_provider = Utils.get_pems_dp(this.ni,params.configname);
            
        end
        
    end
    
end

