classdef StackEnsembleFlowPrediction < FlowPredictor
    
    properties
    end
    
    methods ( Access = public )
        
        function [this]=StackEnsembleFlowPrediction(dp)
            
            this = this@FlowPredictor(dp);
            
        end
        
        
        function  [y_today_stackensemble , stackweight ,y_uptotoday_base  , y_uptotoday_actual] ...
                = predict(this,id,day,from,to,dt , y_historical_actual, y_historical_base , num_basemethods,num_days_for_stackvalidation)
            % function to compute prediction by base learners in a new day
            % and stack ensemble predictions in a new day
            %====================================================================================
            % output:
            % y_today_stackensemble : stack ensemble prediction for today
            % stackweight: stack weight for today
            % y_uptotoday_base: a 3D tensor, augmenting the 3D tensor y_historical_base by
            % appending base predictions for today to the tensor
            % y_uptotoday_actual: a matrix, augmenting the matrix y_historical_actual by
            % appending today's true flow to the tensor           
            
            %====================================================================================
            % input:
            % y_historical_actual : matrix of historical true flow, number of
            % rows is number of available past days
            %
            % y_historical_base : 3D tensor of prediction by each base method in
            % the past. y_historical_base(:,:,k) are the past predictions by the
            % k-th method.  y_historical_base(i,:,:) are the predictions made in the
            % i-th day
            %====================================================================================
            % number of days in the past to look back for training stacking weights       
            %num_days_for_stackvalidation = 60;  
            % number of days to look back for training each base learner
            num_days_for_trainingbasepredictor = 120 ;
            
            
            if isempty(y_historical_base)==1
                % if y_historical_base and  y_historical_actual are empty matrices, compute the historical base predictions     
                [y_pls, y_armax, y_svm , y_gpr , y_krr , y_actual] = StackEnsembleFlowPrediction.compute_prediction_on_validation(this.data_provider,id,day,...
                    num_days_for_stackvalidation,num_days_for_trainingbasepredictor ,...
                    from,to) ;
                
                y_historical_base = zeros([size(y_pls),5]);
                y_historical_base(:,:,1) = y_pls ;
                y_historical_base(:,:,2) = y_armax ;
                y_historical_base(:,:,3) = y_svm ;
                y_historical_base(:,:,4) = y_gpr ;
                y_historical_base(:,:,5) = y_krr ;
                
                y_historical_actual = y_actual ;
                
            end
            
            % compute weights for stack regression from validation period
            % stack_coefficient is m-by-T.
            % m: number of based methods. T:number of time steps per
            % prediction
            [ stackweight ]  = StackEnsembleFlowPrediction.ensembleweight(  y_historical_actual, y_historical_base , num_basemethods ,num_days_for_stackvalidation ) ;
            
            % get training data for base predictor, TODAY
            [training_dp] = StackEnsembleFlowPrediction.construct_trainingdp_includetoday(this.data_provider,day,num_days_for_trainingbasepredictor);
            
            % PLS, TODAY
            PLS_predictor = PLSFlowPredictor(training_dp);
            y_pls_today =  PLS_predictor.predict(id,day,from,to,1) ;
            % ARMAX, TODAY
            params = struct();
            ARMAX_predictor = ObjectFactory.recursiveARMAX_predictor(training_dp,params) ;
            armax_result = ARMAX_predictor.predict(id,day,from,to,900) ;
            y_armax_today = armax_result.flw_out_vph ;
            % svm, TODAY
            SVM_predictor = SVMFlowPredictor(training_dp);
            y_svm_today = SVM_predictor.predict(id,day,from,to,900) ;
            % gaussian process regression, TODAY
            GPR_predictor = GPRFlowPredictor(training_dp);
            y_gpr_today = GPR_predictor.predict(id,day,from,to,900) ;
            % Kernel ridge regression with gaussian kernel, TODAY
            KRR_predictor = KernelRidgeFlowPredictor(training_dp);
            y_krr_today = KRR_predictor.predict(id,day,from,to,900) ;
            % get actual data, TODAY
            actual =training_dp.get_data(day,id,from,to,900,'interpolate');
            y_actual_today =  actual.flw_out_vph;
                        
            Y_join = [y_pls_today ; y_armax_today ; y_svm_today ; y_gpr_today ; y_krr_today] ;
            y_today_stackensemble = zeros(1,(to-from)/900) ;

            
            for interval = 1: length(y_today_stackensemble)
                y_today_stackensemble(interval) =  Y_join(:,interval)' * stackweight(:,interval) ;            
            end
            
            y_uptotoday_base  = zeros( size(y_historical_base) + [1,0,0] ) ;
            y_uptotoday_base(:,:,1) = [y_historical_base(:,:,1) ; y_pls_today ];
            y_uptotoday_base(:,:,2) = [y_historical_base(:,:,2) ; y_armax_today ];
            y_uptotoday_base(:,:,3) = [y_historical_base(:,:,3) ; y_svm_today ];
            y_uptotoday_base(:,:,4) = [y_historical_base(:,:,4) ; y_gpr_today ];
            y_uptotoday_base(:,:,5) = [y_historical_base(:,:,5) ; y_krr_today ];
            
            y_uptotoday_actual = [y_historical_actual ; y_actual_today] ;
            
            
            
        end
        
    end
    
    methods (Static, Access = public)
        
        function [ Y_stackensemble , stackweight ] = ...
                compute_ensemble_only(  y_historicalactual, y_historicalbase , num_basemethods , num_days_for_stacktraining )
            % function to compute stack ensemble predictions, given the historical base prediction
            %====================================================================================
            % output:
            % Y_stackensemble: matrix of stack ensemble predictions. Each
            % row is the prediction for a day. There are
            % size(y_historicalactual,1)- num_days_for_stacktraining rows
            
            % input:
            % y_historicalactual : matrix of historical true flow, number of
            % rows is number of available past days
            %
            % y_historicalbase : tensor of prediction by each base method in
            % the past. y_historical_base(:,:,k) are the past predictions by the
            % k-th method.  y_historical_base(i,:,:) are the predictions made in the
            % i-th day
            %====================================================================================
            
            if size(y_historicalbase,3) ~= num_basemethods
                error('dimension of y_historicalbase does not match number of base methods \n')
            end
            
            predictionhorizon  = size(y_historicalactual,2) ;
            if size(y_historicalactual,1) - num_days_for_stacktraining < 1
                error('does not satisfy the minimum days needed to compute ensemble')
            end
            
            % need at least num_days_for_stacktraining days for stacking
            % store the stack ensemble prediction for each day after
            % the minimum number of days required to learn the ensemble
            Y_stackensemble = zeros( size(y_historicalactual,1)-num_days_for_stacktraining,predictionhorizon);
            % store the stack coefficient for each method for each day after
            % the minimum number of days required to learn the ensemble 
            stackweight = zeros(num_basemethods , predictionhorizon , size(y_historicalactual,1)-num_days_for_stacktraining ) ;
            
            
            for day_idx = num_days_for_stacktraining + 1: size(y_historicalactual,1)
                
                % extracting the base learner prediction in the past num_days_for_stacktraining days
                Y_actual_for_stack = y_historicalactual(day_idx - num_days_for_stacktraining : day_idx - 1 ,:);
                Y_base_for_stack = y_historicalbase(day_idx - num_days_for_stacktraining : day_idx - 1 ,:,:);
                
                Y_error = zeros(size(Y_base_for_stack)) ;
                for m = 1: num_basemethods
                    Y_error(:,:,m) = Y_actual_for_stack - Y_base_for_stack(:,:,m) ;
                end
                              
                SampleCovariance = zeros(num_basemethods , num_basemethods , predictionhorizon) ;
                for t = 1:predictionhorizon
                    for m = 1: num_basemethods
                        for r= m+1 : num_basemethods
                            SampleCovariance(m,r,t) = Y_error(:,t,m)' * Y_error(:,t,r) ;
                            SampleCovariance(r,m,t) = SampleCovariance(m,r,t) ;
                        end
                        SampleCovariance(m,m,t) = norm(Y_error(:,t,m),2)^2 ;
                    end
                end
                                
                Y_base_today = zeros(num_basemethods , predictionhorizon) ;
                for m = 1:num_basemethods
                    Y_base_today(m,:) = y_historicalbase(day_idx ,:,m);
                end
                
                for t=1:predictionhorizon
                    % compute weight via QP
                    stackweight(:,t,day_idx-num_days_for_stacktraining) = quadprog(SampleCovariance(:,:,t),zeros(num_basemethods,1),[],[],ones(1,num_basemethods),1,zeros(num_basemethods,1),ones(num_basemethods,1)*Inf);                  
                    % linear combination of base predictors
                    Y_stackensemble(day_idx-num_days_for_stacktraining,t) = Y_base_today(:,t)' * stackweight(:,t) ;           
                end
                
                
            end
            
        end
        
    end
    
    methods ( Static , Access = private )
        
        function [training_dp] = construct_trainingdp_includetoday(dp,day,numdays_for_training)
            
            % clone an independent copy of dp; dp is a handle class
            training_dp = feval(class(dp),dp.ids) ;
            if find(training_dp.days == day) <= numdays_for_training
                error('does not satisfy the minimum number of days for training')
            end
            
            today_ind = find(training_dp.days == day) ;
            training_dp.days = training_dp.days(today_ind - numdays_for_training : today_ind);
            training_dp.data = training_dp.data(today_ind - numdays_for_training : today_ind);
            training_dp.vds_is_good =  training_dp.vds_is_good(today_ind - numdays_for_training : today_ind);
            
        end
        
        
        function [validation_days] = construct_validationdays(dp,day,numdays_for_validation)
            day_ind = find(dp.days == day) ;
            if day_ind <= numdays_for_validation
                error('does not satisfy the minimum number of days for valdidation')
            end
            validation_days = dp.days(day_ind-numdays_for_validation:day_ind-1);
            
        end
        
        
        function [y_pls, y_armax, y_svm , y_gpr , y_krr , y_actual] = compute_prediction_on_validation(dp,station_id,testday,...
                numdays_for_validation,numdays_for_training ,...
                from,to)
            
            [validation_days] =  StackEnsembleFlowPrediction.construct_validationdays(dp,testday,numdays_for_validation) ;
            
            y_pls = zeros(length(validation_days) , (to - from) / 900);
            y_armax = zeros(length(validation_days) , (to - from) / 900);
            y_svm = zeros(length(validation_days) , (to - from) / 900);
            y_gpr = zeros(length(validation_days) , (to - from) / 900);
            y_krr = zeros(length(validation_days) , (to - from) / 900);
            y_actual = zeros(length(validation_days) , (to - from) / 900);
            
            for v_day_id = 1:length(validation_days)
                v_day_id
                v_day = validation_days(v_day_id);
                [training_dp] = StackEnsembleFlowPrediction.construct_trainingdp_includetoday(dp,v_day,numdays_for_training);
                
                % gaussian process regression
                GPR_predictor = GPRFlowPredictor(training_dp);
                y_gpr(v_day_id , :) = GPR_predictor.predict(station_id,v_day,from,to,900) ;
                
                % Kernel ridge regression with gaussian kernel
                KRR_predictor = KernelRidgeFlowPredictor(training_dp);
                y_krr(v_day_id , :) = KRR_predictor.predict(station_id,v_day,from,to,900) ;
                
                % svm
                SVM_predictor = SVMFlowPredictor(training_dp);
                y_svm(v_day_id , : ) = SVM_predictor.predict(station_id,v_day,from,to,900) ;
                
                
                % PLS
                PLS_predictor = PLSFlowPredictor(training_dp);
                y_pls(v_day_id , :) =  PLS_predictor.predict(station_id,v_day,from,to,1) ;
                
                % ARMAX
                params = struct();
                ARMAX_predictor = ObjectFactory.recursiveARMAX_predictor(training_dp,params) ;
                armax_result = ARMAX_predictor.predict(station_id,v_day,from,to,900) ;
                y_armax(v_day_id , :) = armax_result.flw_out_vph ;
                
                % get actual data
                actual =training_dp.get_data(v_day,station_id,from,to,900,'interpolate');
                y_actual( v_day_id , :) =  actual.flw_out_vph;
                
                
            end
            
        end
        
        
        function [stackweight ] = ensembleweight(  y_historicalactual, y_historicalbase , num_basemethods , num_days_for_stacktraining  )
            % private function to compute stack coefficients for each base learner
            %====================================================================================
            % input:
            % y_historicalactual : matrix of historical true flow, number of
            % rows is number of available past days
            %
            % y_historicalbase : tensor of prediction by each base method in
            % the past. y_historical_base(:,:,k) are the past predictions by the
            % k-th method.  y_historical_base(i,:,:) are the predictions made in the
            % i-th day
            %====================================================================================
            
            if size(y_historicalbase,3) ~= num_basemethods
                error('dimension of y_historicalbase does not match number of base methods \n')
            end
              
            predictionhorizon  = size(y_historicalactual,2) ;
            
            % extracting true flow in the past num_days_for_stacktraining days
            Y_actual = y_historicalactual(end - num_days_for_stacktraining + 1 : end ,:);
            % extracting prediction by each base learner in the past  num_days_for_stacktraining days
            Y_base = y_historicalbase(end - num_days_for_stacktraining + 1 : end ,:,:);
            Y_error = zeros(size(Y_base)) ;
            for m = 1: num_basemethods
                Y_error(:,:,m) = Y_actual - Y_base(:,:,m) ;
            end
            
            SampleCovariance = zeros(num_basemethods , num_basemethods , predictionhorizon) ;
            
            for t = 1:predictionhorizon
                for m = 1: num_basemethods
                    for r= m+1 : num_basemethods
                        SampleCovariance(m,r,t) = Y_error(:,t,m)' * Y_error(:,t,r) ;
                        SampleCovariance(r,m,t) = SampleCovariance(m,r,t) ;
                    end
                    SampleCovariance(m,m,t) = norm(Y_error(:,t,m),2)^2 ;
                end
            end
            
            stackweight = zeros(num_basemethods , predictionhorizon) ;
            
            for t=1:predictionhorizon
                % compute weight via QP
                stackweight(:,t) = quadprog(SampleCovariance(:,:,t),zeros(num_basemethods,1),[],[],ones(1,num_basemethods),1,zeros(num_basemethods,1),ones(num_basemethods,1)*Inf);
            end
            
        end
        
        
    end
    
    
    
end



