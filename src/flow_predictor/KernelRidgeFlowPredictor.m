classdef KernelRidgeFlowPredictor < FlowPredictor
    
    properties
    end
    
    methods ( Access = public )
        
        function [this]=KernelRidgeFlowPredictor(dp)
            
            this = this@FlowPredictor(dp);
            
        end
        
        function [Y_predicted] = predict(this,id,day,from,to,dt)
            
            % data provider
            training_dp = this.data_provider ;
            
            % construct training data matrix
            % X is D-by-T matrix of training data for up-to-now flows. D is the number of days
            % in training_dp with today removed. 
            X = zeros( nnz(training_dp.days(training_dp.days ~= day)) , from/900 );
            
            % Y is D-by-T' matrix of training data for flows during the prediction horizon. D is the number of days
            % in training_dp with today removed.             
            Y = zeros(nnz(training_dp.days(training_dp.days ~= day)) , ( to - from )/900) ;
            
            idx = find(training_dp.days ~= day) ;
            for i=1:length(idx)
                X(i,:) = training_dp.data(idx(i)).flw_out_vph(1: from/900);
                Y(i,:) =  training_dp.data(idx(i)).flw_out_vph(from/900+1:to/900 );
            end
            % mean shifted training flows up-to-now
            X_centered = X - ones(size(X,1),1)*mean(X,1) ;
            % mean shifted training flows during the prediction horizon
            Y_centered = Y - ones(size(Y,1),1)*mean(Y,1) ;
            
            todayidx = find(training_dp.days == day) ;
            
            lambda = 0.3;
            gaussianwidth = 10;
            Y_predicted = KernelRidgeFlowPredictor.KernelRidge(Y_centered, X_centered , training_dp.data(todayidx).flw_out_vph(1: from/900)-mean(X,1) ,'gaussian', lambda , gaussianwidth ) ;
            % adding back mean flows during the prediction horizon
            Y_predicted = Y_predicted + mean(Y,1); 
            
        end
        
    end
    
    methods ( Static , Access = private )
        
        function [training_dp] = construct_trainingdp_includetoday(dp,day,numdays_for_training)
            
            % clone an independent copy of this; this is a handle class
            training_dp = feval(class(dp),dp.ids) ;
            if find(training_dp.days == day) <= numdays_for_training
                error('does not satisfy the minimum number of days for training')
            end
            training_day_ind = dp.days <= day & dp.days >= day - numdays_for_training;
            training_dp.days = training_dp.days(training_day_ind);
            training_dp.data = training_dp.data(training_day_ind);
            training_dp.vds_is_good =  training_dp.vds_is_good(training_day_ind);
        end
        
        function  y_predict = KernelRidge(y_training, X_training , X_test , kerneloption, lambda , kernelparameter )
            
            if strcmp(kerneloption , 'polynomial')
                % polynomial kernel
                K = (1 + X_training*X_training').^kernelparameter ;
                % kernel between x_test and X_training k(x_test,X_training)
                K_test = (1 + X_test*X_training').^kernelparameter ;
            elseif strcmp(kerneloption , 'gaussian')
                
                %========================================================
                % compute Kernel K(i,j) between all training samples
                X_training_normalize = ( X_training - ones(size(X_training,1),1) * mean(X_training) ) ./ ( ones(size(X_training,1),1) * std(X_training) ); 
                K = zeros(size(X_training,1),size(X_training,1)) ;
                for i=1:size(X_training,1)
                    for j = i+1 : size(X_training,1 )
                        K(i,j) = exp(-norm(X_training_normalize(i,:) - X_training_normalize(j,:) ,2)^2 / (2*kernelparameter^2)  ) ;
                        K(j,i) = K(i,j) ;
                    end
                end
                K = K + diag(ones(size(X_training,1),1)) ;
                %========================================================
                % compute Kernel between the testing data and all
                % training data      
                X_test_normalize = ( X_test - ones(size(X_test,1),1) * mean(X_training) ) ./ ( ones(size(X_test,1),1) * std(X_training) ) ; 
                K_test=zeros(size(X_test,1) , size(X_training,1)) ;
                for i=1:size(X_test,1)
                    for j = 1: size(X_training,1 )
                        K_test(i,j) = exp(-norm(X_test_normalize(i,:) - X_training_normalize(j,:) ,2)^2 / (2*kernelparameter^2)  ) ;
                    end
                end
                %========================================================
            end
            
            % lagrange multiplier
            Kinv = inv((K + lambda*eye(size(K,1))) ) ;
            alpha = Kinv * y_training ;
            % computing prediction
            y_predict = K_test*alpha;
            
        end
        
    end
    
end
