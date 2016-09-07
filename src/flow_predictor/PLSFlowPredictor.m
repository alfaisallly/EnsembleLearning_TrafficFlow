classdef PLSFlowPredictor < FlowPredictor
    
    properties
    end
    
    methods ( Access = public )
        
        function [this]=PLSFlowPredictor(dp)
            
            this = this@FlowPredictor(dp);
            
        end
        
        function [Y_predicted] = predict(this,id,day,from,to,num_clusters,dt)
            
            num_components = 20;
            
            if nargin>6
                error('i dont like it')
            end
            %
            %             % get data matrix
            %             x=this.data_provider.get_data('all','all',0,86400,900,'interpolate');
            %             flows = vertcat(x.flw_out_vph);
            %             clear x
            %
            %             % keep "day" data
            %             day_ind = (this.data_provider.days==day)';
            %             if ~any(day_ind)
            %                 error('~any(day_ind)')
            %             end
            %
            %             Cutoff_Time = floor(from/900);
            %
            %             % construct training
            %             Z_training = flows(~day_ind,1:Cutoff_Time);
            %             Y_training = flows(~day_ind,Cutoff_Time+1:end);
            %
            %             % construct today data (Z_testing)
            %             Z_testing = flows(day_ind,1:Cutoff_Time);
            %
            %             clear day_ind flows
            
            [Z_training,Y_training,Z_testing] = this.gather_data_for_prediction( ...
                day,id,from,900);
                        
            % cluster
            [clusterID,Testing_ClusterID] = ...
                PLSFlowPredictor.kmeans_clustering(Z_training,num_clusters,Z_testing);
            
            % select data from nearest cluster for training
            Z_training = Z_training(clusterID==Testing_ClusterID,:);
            Y_training = Y_training(clusterID==Testing_ClusterID,:);
            
            % PLS Z|Y partition, Z is data before the cutoff time, Y is the data
            % after the cutoff time
            
            % decomposition
            [~,P,C] = PLSFlowPredictor.PLS_Decomposition(Z_training,Y_training,num_components);
            
            % prediction
            Y_predicted = PLSFlowPredictor.PLS_Prediction(P,C,Z_training,Y_training,Z_testing);
            
            % truncation
            to_ind = round((to-from)/900);
            Y_predicted = Y_predicted(1:to_ind);
            Y_predicted = max([Y_predicted;0*Y_predicted]);
            
        end
        
    end
    
    methods ( Static , Access = private )
        
        
        function [clusterID,Testing_ClusterID]=kmeans_clustering(Z_training,K,Z_testing)
            
            if K==1
                clusterID = ones(size(Z_training,1),1);
                Testing_ClusterID = 1;
                return
            end
            
            
            ValidCluster = 0 ;
            trial = 3;
            while ValidCluster == 0
                
                if trial == 0
                    K=K-1;
                    trial = 3;
                    %warningmessage = sprintf('warning: no valid clustering found, automatic decrease K to %d',K);
                    %display(warningmessage);
                end
                

                %[clusterID,clusterCentre]=kmeans(Z_training,K,'EmptyAction','singleton');
                % replacing clustering on distribution for Euclidean space 
                
               Z_VolumeDistribution = Z_training./(sum(Z_training,2) * ones(size(Z_training,2),1)');
               [clusterID,clusterCentre] = kmeans(Z_VolumeDistribution,K,'Distance','cityblock','EmptyAction','singleton');
                
                CheckCluster = zeros(K,1);
                for k = 1:K
                    r = find(clusterID == k);
                    if length(r) > 1/6 * size(Z_training,1)
                        CheckCluster(k)=1;
                    end
                end
                
                if sum(CheckCluster) == K
                    ValidCluster = 1 ;
                end
                
                trial = trial - 1;
            end
            
            % find nearest cluster
            K=size(clusterCentre,1);
            distance = zeros(K,1);
            Z_testing_distribution = Z_testing./(sum(Z_testing,2) * ones(size(Z_testing,2),1)');
            for k=1:K
                distance(k) = norm(Z_testing_distribution-clusterCentre(k,:),1);
            end
            [~,Testing_ClusterID] = min(distance) ;
            
        end
             
  
        function [W,P,C] = PLS_Decomposition(Z,Y,N)
            %==========================================================
            % PLS: doing PLS decomposition between matrix Z and matrix Y. Z and Y
            % must have the same number of rows. In flow prediction, each row is
            % a training day; rows in Z contains the flow volume before the cut off time;
            % Y contains the volume after the cut off time.
            % Decompose Z ~ WP'
            %           Y ~ WC'
            % each columns in W are the projection onto a direction that
            % maximize the covariance matrix Z'Y, subject to mutual orthogonality
            % N is the number of directions
            %==========================================================
            %==========================================================
            % input:
            % Z: matrix of training data before the cutoff time. each row is a day
            % Y: matrix of training data after the cutoff time. each row is a day.
            %==========================================================
            
            % ######### * * * * * * * * * * * * #############
            % !!!!!!!!!!!!! IMPORTANT REMARK !!!!!!!!!!!!!!!!
            % !!!!! before doing the PLS, the mean of each colume will be subtracted.
            % each colume of Z, and each column of Y will be mean-shifted to 0
            % ######### * * * * * * * * * * * * #############
            
            [DZ,T_cutoff]=size(Z);
            [DY,T_remain]=size(Y);
            if DZ ~= DY
                error('error: Z and Y in PLS must have the same number of rows');
            else
                D = DZ;
            end
            
            % shift the mean to zero
            Z_centered = Z - ones(D,1)*mean(Z,1);
            Y_centered = Y - ones(D,1)*mean(Y,1);
            
            % start PLS
            P=zeros(T_cutoff,N);
            C=zeros(T_remain,N);
            W=zeros(D,N);
            for n=1:N
                
                CovarMatrix = Z_centered'*Y_centered;
                [r,d,s]=svds(CovarMatrix,1);
                W(:,n) = Z_centered * r / norm(Z_centered*r,2);  % PLS direction
                P(:,n) = Z_centered' * W(:,n); % projection onto the PLS direction
                C(:,n) = Y_centered' * W(:,n); % projection onto the PLS direction
                Z_centered = Z_centered - W(:,n) * P(:,n)'; % deflation
                Y_centered = Y_centered - W(:,n) * C(:,n)'; % deflation
                
            end
            
            
        end
        
        function Y_predicted = PLS_Prediction(P,C,Z_training,Y_training,Z_testing)
            %==========================================================
            % make predictions using PLS components computed from PLS_decomposition
            % input:
            % P: matrix of PLS predictor components learned
            % C: matrix of PLS predicted components learned
            % Z_training: matrix of training data before the cutoff time. each row is a day
            % Y_training: matrix of training data after the cutoff time. each row is a day.
            % Z_testing: matrix of testing data before the cutoff time. each row is a day
            
            % output:
            % Y_predicted: matrix of prediction, each row is a prediction. Y_predicted
            % has the same number of rows as Z_testing
            
            % ######### * * * * * * * * * * * * #############
            % !!!!!!!!!!!!! IMPORTANT REMARK !!!!!!!!!!!!!!!!
            % the components P, C computed from PLS_decomposition are after subtract
            % the mean of training data (mean shifted to 0).
            % Thus before making predictions, mean of Z_training should be subtracted
            % from Z_testing. After prediction, mean of Y_training should be added back
            % to the prediction
            % ######### * * * * * * * * * * * * #############
            
            [DZ,T_cutoff]=size(Z_training);
            [DY,T_remain]=size(Y_training);
            if DZ ~= DY
                error('error: Z and Y in PLS must have the same number of rows');
            end
            
            % Z_training,Y_training is only used for computing empirical mean
            Z_mean = mean(Z_training,1);
            Y_mean = mean(Y_training,1);
            [D_testing,T_cutoff]=size(Z_testing);
            
            P_inv = pinv(P);
            % compute loading matrix W_predicted for testing data
            W_predictor = P_inv * ( Z_testing - ones(D_testing,1)*Z_mean)';
            % compute prediction
            Y_predicted = C * W_predictor;
            % add mean of training data
            Y_predicted = Y_predicted' + ones(D_testing,1)*Y_mean;
            
        end
        
    end
    
end
