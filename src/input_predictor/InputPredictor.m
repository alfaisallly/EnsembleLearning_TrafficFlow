classdef (Abstract) InputPredictor < handle

    properties
        ids     % node ids for SplitPredictor
                % link ids for DemandPredictor
        flow_predictor  @FlowPredictor
    end
    
    methods(Access=public)
       
        function [this]=InputPredictor(flow_predictor)
            this.ids = [];
            this.flow_predictor = flow_predictor;
        end
        
        function [b] = check_have_day(this,day)
            b = this.flow_predictor.check_have_day(day);
        end

        function [this] = set_force_bad_sensors(this,day,ids)
            this.flow_predictor.set_force_bad_sensors(day,ids);
        end
        
    end
    
    methods(Abstract)
        
        % returns structure with predictions for all ids
        [x] = predict(this,day,from,to,dt)
    end
    
end

