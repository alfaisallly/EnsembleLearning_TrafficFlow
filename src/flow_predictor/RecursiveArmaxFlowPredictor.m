classdef RecursiveArmaxFlowPredictor < FlowPredictor
    properties
        params_armax            % hash map
        
        % for recursive armax
        polyOrder               % order of armax
        initValue               % initial values
        initCovariance          % initial covariance
        estimationMethod        % estimation method
        
        % for armax
        lambda                  % Regularization
        initial_condition       % initial condition
        focus                   % focus
        search_method           % search method
        
        lag                     % smoothing window in seconds
        num_of_pred_Horizon     % to determine the length of representative/real data used
    end
    
    methods ( Access = public )
        
        function [this]=RecursiveArmaxFlowPredictor(pems_dp,params)
        % params is a struct with fields: 
        % polyOrder initValue initCovariance, estimation_method, lambda,
        % num_of_pred_Horizon, lag
            
            this = this@FlowPredictor(pems_dp);
            
            default_params = RecursiveArmaxFlowPredictor.get_default_parameters;
            
            % ARMAX variables, default values
            % A(q) = 1+a_1*q^-1 + ... + a_na*q^-na
            % B(q) = b_0+b_1*q^-1 + ... + b_nb*q^-nb
            % C(q) = 1+c_1*q^-1 + ... + c_nc*q^-nc
            % input-output delay: nk
            % polyOrder=[na nb nc nk]
            if isfield(params,'polyOrder')
                this.polyOrder = params.polyOrder; 
            else
                this.polyOrder = default_params.polyOrder;
            end
            
            if isfield(params,'initValue')
                this.initValue = params.initValue; 
            else
                this.initValue = default_params.initValue;
            end    
            
            % ******************for recursive ARMAX******************** 
            % if you trust your initial stettings, use smaller values, e.g., 0.1.
            % default value is 10000
            if isfield(params,'initCovariance')
                this.initCovariance = params.initCovariance;
            else
                this.initCovariance = default_params.initCovariance;
            end
            
            %estimation method: a structure with fields 'name' and 'value'
            if isfield(params,'estimation_method')
                this.estimationMethod = params.estimation_method;
            else
                this.estimationMethod = default_params.estimationMethod; 
            end
     
            %initialization
            this.params_armax=containers.Map('KeyType','double');
            this.params_armax = containers.Map(pems_dp.ids,...
                repmat([0 this.initValue.a0 this.initValue.b0 this.initValue.c0],1,numel(pems_dp.ids)));
            
            % *******************for ARMAX****************************
            % Regularzation
            if isfield(params,'lambda')
                this.lambda = params.lambda;
            else 
                this.lambda = default_params.lambda;
            end
            
            % Initial condition: auto, zero, estimate, backcast
            this.initial_condition='auto';
            
            % Focus: prediction, simulation, stability
            this.focus='simulation';
            
            % Search method: Gauss-Newton (gn), Adaptive Gauss-Newton (gna)
            % Levenberg-Marquardt (lm), (lsqnonlin), steepest descent(grad)
            % (auto: gn, gna, lm, grad)
            this.search_method='auto';
            
            %**************for prediction horizon and time lag***********
            %look back to the historical data: one prediction horizon
            if isfield(params,'num_of_pred_Horizon')
                this.num_of_pred_Horizon = params.num_of_pred_Horizon;
            else
                this.num_of_pred_Horizon = default_params.num_of_pred_Horizon;  
            end
            
            %smoothing of the data
            if isfield(params,'lag')
                this.lag = params.lag;
            else
                this.lag = default_params.lag; %default: no lag
            end
            
        end
        
        function [y] = predict(this,ids,day,from,to,dt)
            % This function provides predicted flows for given vds ID, day
            % and time period using recursive ARMAX/ ARMAX
            
            if numel(day)~=1
                error('numel(day)~=1')
            end
            
            % get health and template use from pems
            ids_health = this.data_provider.get_health(day,ids);
            uses_template = this.data_provider.uses_template(day,ids);
            
            % recursive ARMAX/ ARMAX
            if(isnan(dt))
                dt=this.data_provider.time(end)-this.data_provider.time(end-1);
            end
            time=(from:dt:to);
            predictionHorizon=time(2:end);
            
            % loop through ids
            y = repmat(DataProfile,1,length(ids));
            for i=1:length(ids)
                [y(i)] = this.recursive_armax(ids(i),day,from,to,dt,ids_health(i),uses_template(i),time, predictionHorizon);
            end
            
        end
        
        function [y]=recursive_armax(this,id,day,from,to,dt,health,template,time,predictionHorizon)
            
            if(template || ~health)
                y = this.data_provider.get_data(day,id,from,to,dt);
                return
            end
            
            %do not use template and healthy, recursive ARMAX/ ARMAX
            if(from<this.num_of_pred_Horizon*(to-from)+this.lag) %less than the default length, return historical profile
                y = this.data_provider.get_representative_data(struct('day_of_week',weekday(datenum(day))),id,from,to,dt);
                return
            end
            
            % get data for prediction
            [ ~,~,today_flw , historical_uptonow_cluster,historical_later_cluster] = ...
              this.gather_data_for_prediction(day,id,from,dt,to,nan,this.lag);
            
            u_match = historical_uptonow_cluster.flw_out_vph;
            u_pred = historical_later_cluster.flw_out_vph;
            y_match = today_flw;
%             if(this.lag>0)
%                 u_match=tsmovavg(u_match,'s',this.lag);
%                 
%                 u_pred=tsmovavg(u_pred,'s',this.lag);
% %                 u_pred=u_pred(this.lag+1:end);
%                 
%                 y_match=tsmovavg(y_match,'s',this.lag);
% %                 y_match=y_match(this.lag+1:end);
%             end
            
            %parameter update and armax prediction
            params_pre = this.get_params_for_ids(id);
            [params_now,w] = this.update_params_for_time_online_simplified(from,u_match,y_match,params_pre);
            this.update_hash_map(id,params_now);
            
            % D-step prediction
            y_pred = this.DStep_prediction_recursive(u_pred, u_match, y_match, w, predictionHorizon, params_now);
            
            % eliminate negative values
            y_pred(y_pred<0) = 0;
            
            inflow=y_pred;
            outflow=inflow;
            y = DataProfile( ...
                id , day , time, ...
                inflow , ...
                outflow , ...
                nan(1,length(time)) );
            
        end
        
        function [ y_pred ] = DStep_prediction(this,u_pred, u, y, predictionHorizon, params_now)
            % Predict y(t+d) by solving the bezout equation
            % This function doesn't need the error input: w
            % But it is slower than the recursive one since it needs to solve the bezout
            % equation and uses transfer functions
            
            % Get A, B, and C 
            A=params_now(2:length(this.initValue.a0)+1);
            B=params_now(2+length(this.initValue.a0):length(this.initValue.a0)+length(this.initValue.b0)+1);
            C=params_now(length(this.initValue.a0)+length(this.initValue.b0)+2:end);
            
            % get transfer functions
            dtSec=predictionHorizon(end)-predictionHorizon(end-1);
            Bzz = tf(B,1,dtSec,'Variable','z^-1');
            Czz = tf(C,1,dtSec,'Variable','z^-1');
            
            % Initialization
            y_pred=nan(1,length(predictionHorizon));
            
            for D = 1: length(predictionHorizon)
                % Bezout equation
                % Solving: C(q^(-1))=F(q^(-1))*A(q^(-1))+q^(-D)*G(q^(-1))
                Delayzz = zeros(1,D);
                Delayzz(end) = 1;
                [F,G] = this.bezout(A,Delayzz,C);

                % Calculate polynomial F and G in predictor
                Fzz = tf(F,1,dtSec,'Variable','z^-1');
                Gzz = tf(G,1,dtSec,'Variable','z^-1');
                
                % Prediction
                % y^hat(t+D|t)=G/C*y(t)+FB/C*u(t)               
                Hy = Gzz/Czz;         % = G/C
                Hu = (Fzz*Bzz)/Czz;   % = FB/C
                
                yGC = lsim(Hy,y); % G/C*y
                uFBC = lsim(Hu,[u,u_pred(1:D)]);    % FB/C*u
                y_pred(D) = yGC(end) + uFBC(end);  % predicted flow
            end
        end
        
        function [ y_pred ] = DStep_prediction_recursive(this,u_pred, u, y, w, predictionHorizon, params_now)
            % Predict y(t+d) using the following recursive function
            % This function is very close to solve the bezout equation
            % The accuracy is impacted by the calculation of w
            
            % y^hat(t+d|t)=-a_1*y^hat(t+d-1|t)-...-a_{d-1}*y^hat(t+1|t)
            %    -a_d*y(t)-...-a_n*y(t+d-n)
            %    +c_d*w(t)+...+c_l*w(t+d-l) + B(q^-1)*u(t+d)
            
            % get A, B, and C
            A=params_now(2:length(this.initValue.a0)+1);
            B=params_now(2+length(this.initValue.a0):length(this.initValue.a0)+length(this.initValue.b0)+1);
            C=params_now(length(this.initValue.a0)+length(this.initValue.b0)+2:end);
            
            %initialization
            y_pred=nan(1,length(predictionHorizon));
            
            %flip y, u, and w
            %time sequence: before: ..., t-2, t-1, t; and after: t, t-1, t-2,....
            y_hat=fliplr(y);
            u_hat=fliplr(u);
            w_hat=fliplr(w);
            
            %D-step recursive prediction
            for D = 1: length(predictionHorizon)                
                u_hat=[u_pred(D) u_hat]; %update u
                
                y_bar=-A(2:end)*y_hat(1:length(A)-1)'; %get the y part
                u_bar =B*u_hat(1:length(B))';  %get the u part
                w_bar=C(2:end)*w_hat(1:length(C)-1)'; %get the w part
                
                y_pred(D)=y_bar+u_bar+w_bar; %predicted y at t+D
                
                y_hat=[y_pred(D) y_hat]; % update y               
                w_hat=[0 w_hat]; %update w
            end
        end
        
        function [params] = get_params_for_ids(this,ids)
            if(isempty(this.params_armax))
                error('isempty(this.params_armax)')
            else
                params= this.params_armax(ids);
            end
        end
        
        function []= update_hash_map(this,ids,params)
            this.params_armax(ids) = params;
        end
        
        function [params_now,w]= update_params_for_time_online_simplified(this,time,u,y,params_pre)
            % This function uses recursive ARMAX to estimate A, B, and C.
            % We use our simplified recursive ARMAX, which is seperated
            % from the Matlab toolbox.
            
            if( params_pre(1)< time-3600*2) %not exist or too old, use default values
                params_pre=[time this.initValue.a0 this.initValue.b0 this.initValue.c0];
            end
            
            % Get initial A, B, and C
            value.a0=params_pre(2:length(this.initValue.a0)+1);
            value.b0=params_pre(2+length(this.initValue.a0):length(this.initValue.a0)+length(this.initValue.b0)+1);
            value.c0=params_pre(length(this.initValue.a0)+length(this.initValue.b0)+2:end);
            
            params=struct(...
                'polyOrder', this.polyOrder,...
                'initValue', value,...
                'estimationMethod', this.estimationMethod,...
                'initCovariance',this.initCovariance);

            % Create a recursiveARMAX object
            obj=recursive_armax(params);
            
            % recursive ARMAX estimation
            % Make u and y of the same length
            if(length(u)>length(y))
                y=[zeros(1,length(u)-length(y)),y];
            else
                u=[zeros(1,length(y)-length(u)),u];
            end
            [A,B,C,y_hat]=obj.estimation(u,y);

            % Get the estimation erriors and return w and params_now
            w=y-y_hat;
            params_now= [time A(end,:) B(end,:) C(end,:)];
        end
        
        function [params_now,w]= update_params_for_time_online(this,time,u,y,params_pre)
            % This function uses recursive ARMAX to estimate A, B, and C
            if( params_pre(1)< time-3600*2) %not exist or too old, use default values
                params_pre=[time this.initValue.a0 this.initValue.b0 this.initValue.c0];
            end
            
            % Make them of the same length
            if(length(u)>length(y))
                y=[zeros(1,length(u)-length(y)),y];
            else
                u=[zeros(1,length(y)-length(u)),u];
            end
            
            % recursiveARMAX
            % Get initial A, B, and C
            a0=params_pre(2:length(this.initValue.a0)+1);
            b0=params_pre(2+length(this.initValue.a0):length(this.initValue.a0)+length(this.initValue.b0)+1);
            c0=params_pre(length(this.initValue.a0)+length(this.initValue.b0)+2:end);
            
            % Create a recursiveARMAX object
            obj=recursiveARMAX(this.polyOrder, a0, b0, c0,'EstimationMethod',this.estimationMethod.name);
            
            % recursive ARMAX options
            if(~isnan(this.initCovariance))
                obj.InitialParameterCovariance=this.initCovariance;
            end
            if(strcmp(this.estimationMethod.name,'ForgettingFactor'))
                obj.ForgettingFactor=this.estimationMethod.value;
            elseif(strcmp(this.estimationMethod.name,'KalmanFilter'))
                obj.ProcessNoiseCovariance=this.estimationMethod.value;
            end
            
            % recursive ARMAX estimation
            EstimatedOutput=nan(size(u));
            for i=1:length(u)
                % A, B, and C will be updated every time
                [A, B, C, EstimatedOutput(i)]=step(obj,y(i),u(i));
            end
            
            % Get the estimation erriors and return w and params_now
            w=y-EstimatedOutput;
            params_now= [time A B C];
        end
        
        function [params_now, w]= update_params_for_time_offline(this,time,u,y,params_pre,dt)
            % This function uses ARMAX to estimate A, B, and C
            % Initial params are not needed in this case
            
            % Make them of the same length
            if(length(u)>length(y))
                y=[zeros(1,length(u)-length(y)),y];
            else
                u=[zeros(1,length(y)-length(u)),u];
            end
            
            % ARMAX
            % Create an iddata object
            data=iddata(y',u',dt);
            
            % ARMAX options
            opt=armaxOptions;
            % Regularization
            if(~isempty(this.lambda))
                opt.Regularization.Lambda=this.lambda;
            end            
            % Initial condition: auto, zero, estimate, backcast
            if(~isempty(this.initial_condition))
                opt.InitialCondition=this.initial_condition;
            end            
            % Focus: prediction, simulation, stability
            if(~isempty(this.focus))
                opt.Focus=this.focus;
            end            
            % Search method: Gauss-Newton (gn), Adaptive Gauss-Newton (gna)
            % Levenberg-Marquardt (lm), (lsqnonlin), steepest descent(grad)
            % (auto: gn, gna, lm, grad)
            if(~isempty(this.search_method))
                opt.SearchMethod=this.search_method;
            end
            
            % Call the ARMAX function
            sys=armax(data,this.polyOrder,opt);
            pred=compare(data,sys);

%             Simulation=armax(data,this.polyOrder,opt);
%             pred_sim=compare(data,Simulation);
%             
%             opt.Focus='prediction';
%             Prediction=armax(data,this.polyOrder,opt);
%             pred_pred=compare(data,Prediction);
%             
%             opt.Focus='stability';
%             Stability=armax(data,this.polyOrder,opt);
%             pred_stab=compare(data,Stability);
%             
%             t=(1:length(u))/10;
%             plot(t, u,'--k');
%             hold on
%             plot(t,y,'k','LineWidth',2)
%             plot(t,pred_pred.y,'b','LineWidth',2)
%             plot(t,pred_stab.y,'g','LineWidth',2)
%             plot(t,pred_sim.y,'r','LineWidth',2)
%             legend('u','y',sprintf('Prediction: MAPE=%.2f %%',mean(abs(pred_pred.y-y')./y'*100)),...
%                 sprintf('Stability: MAPE=%.2f %%',mean(abs(pred_stab.y-y')./y'*100)),sprintf('Simulation: MAPE=%.2f %%',mean(abs(pred_sim.y-y')./y'*100)))
%             xlabel('Time (s)')
%             ylabel('Value')
%             close
            
            % Calculate estimation errors and return w and params_now
            w=(data.y-pred.y)';
            params_now= [time sys.A sys.B sys.C];
        end
        
    end
    
    methods ( Access = private )
        function [R,S] = bezout(this,A,B,C)
            % solution of the Bezout Equation
            % C(q-1) = A(q^(-1)) * R(q^(-1)) + B(q^(-1))* S(q^(-1))
            %
            % A, B and C are polynomials in  q^(-1) = 1/q
            % A and C must be polynomials with 1 as their constant coefficient
            %       (i.e. A(1) = 1 and C(1) = 1)
            % A and B must be co-prime
            %
            % returns polynomials R and S
            
            if (nargin ~=4)
                error('No enough or too many input arguments')
            end
            
            %get compact orders of A, B, and C
            a = this.unpad(A);
            b = this.unpad(B) ;
            c = this.unpad(C);
            
            if c(1) ~= 1 || a(1) ~= 1
                error('either a or c do not have 1 as the constant coefficient')
            end
            n = length(a)-1; % a_1 to a_n
            m = length(b)-1; % b_1 to b_m (Delay term: q^{-m})
            nc = length(c)-1; % c_1 to c_nc
            
            ns = max([n-1,nc-m-1]); %order of S
            nr = m;                 %order of R
            nv = m + ns + 1;        %order of unknowns
            
            ce = [c(2:end) zeros(1,nv - nc )]; %get the same length as the unknowns
            ae = [a(2:end) zeros(1,nv - n )];
            D=zeros(nv,nv);
            
            %fill in the a's
            for jj=1:m
                D(jj:jj+n,jj) = a';
            end
            
            %fill in the b's (unit vectors)
            j=1;
            for jj=m+1:nv
                D(j:j+m,jj) = b';
                j=j+1;
            end
            
            %solve (ce-ae)' = D * v'
            v = (D \ (ce - ae)')';
            
            %extract R and S
            R = [1 v(1:m)];
            S = v(m+1:end);
        end
        
        function vout = unpad(this,v)
            if (size(v,1) > 1)
                v=v';%row vector to column vector
            end
            if (size(v,1) > 1)%not a vector, probably a matrix
                error('not a vector')
            end
            %the following code tries to reduce uncessary zeros (reduce the order) 
            r1 =roots(v); 
            num_zeros = sum(r1==0);
            vout = v(1:end-num_zeros);
        end
    end
    
    methods (Static, Access = public )
        
        function [X] = run_and_report_error(params)
            % fields(params) = {ppt_file,xls_file,configfile,sim_dt,output_dt,end_time,update_dt,horizon}
            
            config = Config.get(params.config);
            ni = ObjectFactory.network_information(config.xml_file);
            pems_dp = Utils.get_pems_dp(ni,params.config);
            
            ppt_file = fullfile(Folder.reports,params.ppt_file);
            ppt_file_error = sprintf('%s_error',ppt_file);
            
            flow_predictor = ObjectFactory.recursiveARMAX_predictor(pems_dp);
            dp_ids = flow_predictor.data_provider.ids;                              % internal ids to be reported
            link_ids = flow_predictor.data_provider.get_linkids_for_ids(dp_ids);    % corresponding link ids
            
            dt = flow_predictor.data_provider.measurement_dt;
            update_times = 0:params.update_dt:(86400-params.update_dt);
            n_horizon = params.horizon/dt;
            
            % allocate prediction array
            prediction = repmat(struct('link_id',nan,'dp_id',nan,'flw',nan(length(update_times),n_horizon)),1,length(dp_ids));
            for i=1:length(dp_ids)
                prediction(i).dp_id = dp_ids(i);
                prediction(i).link_id = link_ids(i);
            end
            
            % collect true flows from pems data provider
            fprintf('\tCollecting actual flows\n')
            actual = pems_dp.get_data(params.day,dp_ids,'start','end');
            actual_rep = pems_dp.get_representative_data(struct('day_of_week',weekday(datenum(params.day))),dp_ids,'start','end');
            
            
            % run predictions
            
            fprintf('\tRunning %d predictions \n',length(update_times))
            flow_predictor.num_of_pred_Horizon=params.num_of_pred_Horizon;
            for k=1:length(update_times)
                fprintf('step=%d \n',k)
                start_time = update_times(k);
                x = flow_predictor.predict(dp_ids, ...
                    params.day , ...
                    start_time, ...
                    start_time+params.horizon, ...
                    dt);
                
                lasttime = min([actual(1).time(end)-dt start_time+params.horizon-dt]);
                ind = Utils.index_into(start_time:dt:lasttime,actual(1).time);
                for i=1:length(dp_ids)
                    err = Utils.error(actual(i).flw_out_vph(ind),x(i).flw_out_vph(1:length(ind)));
                    prediction(i).error.mae(k) = err.mae;
                    prediction(i).error.mape(k) = err.mape;
                    prediction(i).error.smape(k) = err.smape;
                    prediction(i).error.rmse(k) = err.rmse;
                    prediction(i).error.linf(k) = err.linf;
                    prediction(i).flw(k,:) = x(i).flw_out_vph;
                end
            end
            fprintf('\n')
            
            % error table : NOTE: THIS IS AN INCORRECT AVERAGE, CHANGE LATER
            for i=1:length(dp_ids)
                prediction(i).avg_error = Utils.avg_error(prediction(i).error);
            end
            
            if isfield(params,'xls_file')
                X = [prediction.avg_error];
                write(table([prediction.dp_id nan]', ...
                    [[X.mae]';Utils.meanwithnan([X.mae]')], ...
                    [[X.mape]';Utils.meanwithnan([X.mape]')], ...
                    [[X.smape]';Utils.meanwithnan([X.smape]')], ...
                    [[X.rmse]';Utils.meanwithnan([X.rmse]')], ...
                    [[X.linf]';Utils.meanwithnan([X.linf]')], ...
                    'VariableNames',{'id','mae','mape','smape','rmse','linf'}) , ...
                    fullfile(Folder.reports,sprintf('%s.xls',params.xls_file)) );
                clear X
            end
            
            % collect historical flows from data provider
            fprintf('\tCollecting historical flows\n')
            historical = flow_predictor.data_provider.get_data('all',dp_ids,'start','end');
            
            % loop health from data_provider [vds x days]
            all_is_good = flow_predictor.data_provider.get_health('all',dp_ids);
            
            % plot predictions for each id
            do_export = ~isempty(params.ppt_file);
            if do_export
                [ppt,op]=openppt(ppt_file,true);
                [ppte,ope]=openppt(ppt_file_error,true);
                
                % title slide
                addslideText(op,'Flow prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d sec.\nPrediction horizon: %d sec.\n', ...
                    class(flow_predictor),class(flow_predictor.data_provider),params.update_dt,params.horizon) )
                
                % title slide
                addslideText(ope,'Flow prediction report',...
                    sprintf('Flow predictor: %s\nData provider: %s\nPrediction frequency: %d\nPrediction horizon: %d\n', ...
                    class(flow_predictor),class(flow_predictor.data_provider),params.update_dt,params.horizon) )
            end
            
            for i=1:length(dp_ids)
                
                fprintf('\tSlide %d of %d\n',i,length(dp_ids))
                
                if do_export
                    figure('Position',[279    95   710   571],'Visible','off')
                else
                    figure('Position',[279    95   710   571],'Visible','on')
                end
                
                xlab = 'time [hr]';
                
                %                 % plot good historical data
                vds_is_good = all_is_good(:,i);
                %                 if any(vds_is_good)
                %                     historical_time = historical(1,1).time(1:end-1)/3600;
                %                     historical_flw = vertcat(historical(:,i).flw_out_vph);
                %                     plot(historical_time,historical_flw(vds_is_good,:),'k');
                %                 end
                
                % plot representative data
                
                rep_time = actual_rep(1,1).time(1:end-1)/3600;
                rep_flw = actual_rep(1,i).flw_out_vph;
                [h1]=plot(rep_time,rep_flw,'k');
                
                
                % actual
                hold on
                [h2]=plot(actual(i).time(1:end-1)/3600,actual(i).flw_out_vph,'-.b','LineWidth',2);
                
                % prediction
                P = prediction(i);
                for k=1:length(update_times)
                    time = (update_times(k):dt:update_times(k)+params.horizon-dt)/3600;
                    [h3]=plot(time,P.flw(k,:),'r','LineWidth',2);
                end
                
                xlabel(xlab);
                ylabel('flow [vph]');
                set(gca,'XLim',[0 24])
                legend([h1 h2 h3],'Representative','Actual','Predicted','Location','NorthWest')
                grid
                
                percent_good = 100*sum(vds_is_good)/numel(vds_is_good);
                if do_export
                    addslide(op,sprintf('(good days:%.0f%%) id=%d',percent_good,dp_ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                else
                    title(sprintf('Flow predictions, id=%d (good days:%.0f%%)',dp_ids(i),percent_good))
                end
                
                
                if do_export
                    figure('Position',[279    95   710   571],'Visible','off')
                else
                    figure('Position',[279    95   710   571],'Visible','on')
                end
                
                subplotf(2,1,'error [vph]',xlab,update_times/3600,[ ...
                    prediction(i).error.linf ; ...
                    prediction(i).error.rmse ; ...
                    prediction(i).error.mae ],'','','',2);
                set(gca,'XLim',[0 24])
                legend(sprintf('linf %.0f',prediction(i).avg_error.linf), ...
                    sprintf('rmse %.0f',prediction(i).avg_error.rmse), ...
                    sprintf('mae %.0f',prediction(i).avg_error.mae),'Location','Best')
                grid
                
                subplotf(2,2,'error [%]',xlab,update_times/3600,100*[ ...
                    prediction(i).error.mape ; ...
                    prediction(i).error.smape  ],'','','',2);
                set(gca,'XLim',[0 24])
                legend( sprintf('mape %.1f',prediction(i).avg_error.mape), ...
                    sprintf('smape %.1f',prediction(i).avg_error.smape), 'Location','Best')
                grid
                
                if do_export
                    addslide(ope,sprintf('(good days:%.0f%%) id=%d',percent_good,dp_ids(i)),'new',[0.5 0.5],'',0.6)
                    close
                end
                
            end
            
            if do_export
                closeppt(ppt,op)
                closeppt(ppte,ope)
            end
            
        end
        
        function [p] = get_default_parameters()
            p = struct( ...
                'polyOrder',        [1 1 1 0],...
                'initValue',        struct('a0', [1,0.5],'b0', 0.5, 'c0', [1,0.5]),...
                'initCovariance', 0.1,      ...
                'estimationMethod', struct('name', 'KalmanFilter', 'value', 0.01) , ...
                'lambda',1,                ...
                'num_of_pred_Horizon',1,    ...
                'lag',0 );
            
        end
    end
end

