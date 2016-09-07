classdef (Abstract) FreewayStateEstimator < StateEstimator
   
    properties
        fwy 
        link_ind
        link_lengths_mile
    end
    
    methods(Access=public)
                
        function [this] = FreewayStateEstimator(model_runner,data_provider,freeway_information)
            % [this] = FreewayStateEstimator(model_runner,data_provider,freeway_information)
            
            this = this@StateEstimator(model_runner,data_provider, ...
                                       freeway_information.get_ordered_ml_links() );                                     
            this.fwy = freeway_information;
            this.link_ind = Utils.index_into(this.link_ids,freeway_information.link_ids);
            this.link_lengths_mile = freeway_information.link_lengths_meters(this.link_ind)/1609.34;
        end
        
    end
    
    methods(Static)
        
        function [X] = run_and_report_error(params)
            X = StateEstimator.run_and_report_error(nan);
        end
        
    end
    
end

