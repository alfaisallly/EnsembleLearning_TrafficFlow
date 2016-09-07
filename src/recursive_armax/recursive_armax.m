classdef recursive_armax
    
    properties
        na                              % Order of A
        nb                              % Order of B
        nc                              % Order of C
        nk                              % Order of the system delay
        nd                              % Currently not needed in ARMAX: set to zero
        nf                              % Currently not needed in ARMAX: set to zero
        
        InitialA                        % Initial coefficients of A
        InitialB                        % Initial coefficients of B
        InitialC                        % Initial coefficients of C
        InitialD                        % Currently not used in ARMAX: set to zero
        InitialF                        % Currently not used in ARMAX: set to zero
        
        u                               % Training data input: array
        y                               % Training data output: array

        EnableAdaptation                % Enable adaptation: logical 0/1
        estimationMethod                % Estimation method:
                                        %   --ForgettingFactor, KalmanFilter, NormalizedGradient, and Gradient
                                        
        ForgettingFactor                % Forgetting method
        ProcessNoiseCovariance          % Kalman filter
        AdaptationGain                  % Normalized gradient/Gradient
        NormalizationBias               % Normalized gradient
        
        InitialParameters               % Initial parameters: column vector of variables in A, B, and C
        Parameters                      % Current parameters: column vector of variables in A, B, and C 
        InitialParameterCovariance      % Initial parameter covariance
        ParameterCovariance             % Current parameter covariance
        
        PastMeasurements                % Past measurements of y, u, and w (err)
        PastJacobianEstimates

        DataType                        % Data type: default "double"       
        AlgorithmParameters             % Algorithm parameters
    end
    
    
    methods (Access = public)
        function this = recursive_armax(params)
            % Recursive ARMAX constructor
            
            % Data type
            if(isfield(params,'DataType'))
                this.DataType=params.DataType;
            else
                % Default
                this.DataType='double';
            end
            
            % Enable Adaptation
            if(isfield(params,'EnableAdaptation'))
                this.EnableAdaptation=params.EnableAdaptation;
            else
                % Default
                this.EnableAdaptation=true();
            end
            
            % Armax orders
            if(isfield(params,'polyOrder'))
                % Should at least contain: na nb nc nk
                this.na=params.polyOrder(1);
                this.nb=params.polyOrder(2);
                this.nc=params.polyOrder(3);
                this.nk=params.polyOrder(4);
                % Have nd ?
                if(length(params.polyOrder)==5)
                    this.nd=params.polyOrder(5);
                else
                    this.nd=0;
                end
                % Have nf
                if(length(params.polyOrder)==6)
                    this.nf=params.polyOrder(6);
                else
                    this.nf=0;
                end 
            else
                % Default
                this.na=1;
                this.nb=1;
                this.nc=1;
                this.nk=0;
                this.nd=0;
                this.nf=0;
            end

            % Armax initial values
            if(isfield(params,'initValue'))
                % Should at least have: A0, B0, and C0
                this.InitialA=params.initValue.a0;
                this.InitialB=params.initValue.b0;
                this.InitialC=params.initValue.c0;
                % Have d0 ?
                if(length(params.polyOrder)==5)
                    this.InitialD=params.initValue.d0;
                else
                    this.InitialD=0;
                end
                % Have f0 ?
                if(length(params.polyOrder)==6)
                    this.InitialF=params.initValue.f0;
                else
                    this.InitialF=0;
                end 
            else
                this.InitialA=[1 0.5];
                this.InitialB=[0.5];
                this.InitialC=[1 0.5];
                this.InitialD=0;
                this.InitialF=0;
            end   
            % Set initial parameters
            this.InitialParameters=this.setInitialParameters;
            
            % Estimation method and the corresponding parameters
            % It is a structure with fields 'name' and 'value'
            dataType=this.DataType;
            if isfield(params,'estimationMethod')
                if(strcmp(params.estimationMethod.name,'ForgettingFactor'))
                    this.estimationMethod = 1;
                    if(~isempty(params.estimationMethod.value))
                        this.ForgettingFactor = cast(params.estimationMethod.value,dataType);
                    else
                        this.ForgettingFactor = 1;
                    end
                elseif(strcmp(params.estimationMethod.name,'KalmanFilter'))
                    this.estimationMethod = 2;
                    if(~isempty(params.estimationMethod.value))
                        this.ProcessNoiseCovariance = this.ScalarToMatrix...
                            (params.estimationMethod.value,length(this.InitialParameters),false());                      
                    else
                        this.ProcessNoiseCovariance = this.ScalarToMatrix...
                            (0.01,length(this.InitialParameters),false());
                    end                   
                elseif(strcmp(params.estimationMethod.name,'NormalizedGradient'))
                    % Need additional input of bias in this case to avoid
                    % jumps with small gradients
                    this.estimationMethod = 3;
                    if(~isempty(params.estimationMethod.value))
                        this.AdaptationGain = cast(params.estimationMethod.value(1),dataType);
                        this.NormalizationBias = cast(params.estimationMethod.value(2),dataType);
                    else
                        this.AdaptationGain = 1;
                        this.NormalizationBias = eps;
                    end   
                elseif(strcmp(params.estimationMethod.name,'Gradient'))
                    this.estimationMethod = 4;
                    if(~isempty(params.estimationMethod.value))
                        this.AdaptationGain = cast(params.estimationMethod.value(1),dataType);
                    else
                        this.AdaptationGain = 1;
                    end  
                else
                    error('Unknown estimation method!')
                end
            else
                % Default: Kalman Filter
                this.estimationMethod = 2;
                this.ProcessNoiseCovariance = 0.1;
            end
            
            % Algorithm parameters
            p = struct(...
                'na', this.na,...
                'nb', this.nb,...
                'nc', this.nc,...
                'nd', this.nd,...
                'nf', this.nf,...
                'nk', this.nk,...
                'estimationMethod',this.estimationMethod);
            % Get the algorithm parameters
            this.AlgorithmParameters = this.getAlgorithmParameters(p);
            
            % Initial parameter covariance
            if isfield(params,'initCovariance')
                initialParameterCovariance=params.initCovariance;
            else
                initialParameterCovariance = 1e4;
            end
            % Create a symmetric matrix for the parameters  
            this.InitialParameterCovariance = this.ScalarToMatrix...
                (initialParameterCovariance,length(this.InitialParameters),true());
                
            % Initialization of past measuments and past Jacobian estimates
            % Past measurements: [1, na+nb+nk+nc]
            this.PastMeasurements = ...
                zeros(1,this.numberOfRequiredPastMeasurements(this.na,this.nb,this.nc,this.nk),dataType);
            % Past Jacobian estimates: [1, max(na,nc)+max(nb+nk,nc+1)+nc]
            this.PastJacobianEstimates = ...
                zeros(1,this.numberOfRequiredPastJacobianEstimates(this.na,this.nb,this.nc,this.nk),dataType);
            
        end
        
        function [A,B,C,y_hat]=estimation(this,u,y)
            % Main function for recursive ARMAX
            %
            %   Inputs
            %       u:   input measurements (array, e.g., 1-by-m ) (training period)
            %       y:   output measurements (array, e.g., 1-by-n) (training period)
            %
            %  Outputs: store all estimated outputs
            %       A:  output estimates of A (na-by-max(m,n))
            %       B:  output estimates of B (nb-by-max(m,n))
            %       C:  output estimates of C (nc-by-max(m,n))
            %   y_hat:  output estimates of y (1-by-max(m,n))
            
            if(isempty(u)||isempty(y))
                error('Not enough inputs for the training data!')
            end
            
            % Make u and y of the same length
            if(length(u)>length(y))
                y=[zeros(1,length(u)-length(y)),y];
            else
                u=[zeros(1,length(y)-length(u)),u];
            end
            
            % Store all estimated A, B, C, and y_hat during the training period
            A=zeros(length(y),this.na+1);
            B=zeros(length(y),this.nb);
            C=zeros(length(y),this.nc+1);
            y_hat=zeros(1,length(y));
            
            % Parameters and parameter covariance
            this.Parameters = this.InitialParameters;
            % ParameterCovariance is only relevant for ForgettingFactor and KalmanFilter algorithms
            if isEstimationMethodForgettingFactor(this) || ...
                    isEstimationMethodKalmanFilter(this)
                this.ParameterCovariance = this.InitialParameterCovariance;
            else
                % Initializing states, even if they are empty, is
                % necessary for codegen
                this.ParameterCovariance = cast([],this.DataType);
            end
            
            % Step along the training data sets
            for i=1:length(y)
                %************ For comparison and testing purposes***********
                %
                % There are three different idrarmaxstep functions in
                % Matlab (C:\Program Files\MATLAB\R2016a\toolbox\ident\idutils)
                % "idrarmaxstep.m", "idrarmaxstep_double_mex.mexw64", and
                % "idrarmaxstep_single_mex.mexw64"
                %
                % What we are using is the function "idrarmaxstep.m"
                % Also, we don't use lower triangular L to compute.
                % Instead, we use P (P=L*L'.) directly.
                %
                %                 if isempty(coder.target)
                %                     if isa(y,'double')
                %                         [A(i,:),B(i,:),C(i,:),e,...
                %                             this.Parameters, this.ParameterCovariance,...
                %                             this.PastMeasurements,this.PastJacobianEstimates] = ...
                %                             idrarmaxstep_double_mex(u(i),y(i),...
                %                             this.EnableAdaptation,...
                %                             this.getAdaptationParameter, ...
                %                             this.getAdaptationParameter2, ...
                %                             this.Parameters, this.ParameterCovariance, ...
                %                             this.PastMeasurements, this.PastJacobianEstimates,...
                %                             this.AlgorithmParameters);
                %                     else
                %                         [A(i,:),B(i,:),C(i,:),e,...
                %                             this.Parameters, this.ParameterCovariance,...
                %                             this.PastMeasurements,this.PastJacobianEstimates] = ...
                %                             idrarmaxstep_single_mex(u(i),y(i),...
                %                             this.EnableAdaptation,...
                %                             this.getAdaptationParameter, ...
                %                             this.getAdaptationParameter2, ...
                %                             this.Parameters, this.ParameterCovariance, ...
                %                             this.PastMeasurements, this.PastJacobianEstimates,...
                %                             this.AlgorithmParameters);
                %                     end
                %                 else
                %                     [A(i,:),B(i,:),C(i,:),e,...
                %                         this.Parameters, this.ParameterCovariance,...
                %                         this.PastMeasurements,this.PastJacobianEstimates] = ...
                %                         idrarmaxstep(u(i),y(i),...
                %                         this.EnableAdaptation,...
                %                         this.getAdaptationParameter, ...
                %                         this.getAdaptationParameter2, ...
                %                         this.Parameters, this.ParameterCovariance, ...
                %                         this.PastMeasurements, this.PastJacobianEstimates,...
                %                         this.AlgorithmParameters);
                %                 end
                %
                %                 if(i==1)
                %                     ParameterCovariance0=ctrlScalarExpansion...
                %                         (this.ParameterCovariance,length(this.InitialParameters),true(),true());
                %                 end
                %                 [A0,B0,C0,e0,...
                %                     Parameters0, ParameterCovariance0,...
                %                     PastMeasurements0,PastJacobianEstimates0] = ...
                %                     idrarmaxstep(u(i),y(i),...
                %                     this.EnableAdaptation,...
                %                     this.getAdaptationParameter, ...
                %                     this.getAdaptationParameter2, ...
                %                     this.Parameters, ParameterCovariance0, ...
                %                     this.PastMeasurements, this.PastJacobianEstimates,...
                %                     this.AlgorithmParameters);
                %************ For comparison and testing purposes***********
                
                [A(i,:),B(i,:),C(i,:),e,...
                    this.Parameters, this.ParameterCovariance,...
                    this.PastMeasurements,this.PastJacobianEstimates] = ...
                    this.idrarmaxstep_simplified(u(i),y(i),...
                    this.EnableAdaptation,...
                    this.getAdaptationParameter, ...
                    this.getAdaptationParameter2, ...
                    this.Parameters, this.ParameterCovariance, ...
                    this.PastMeasurements, this.PastJacobianEstimates,...
                    this.AlgorithmParameters);
                
                % Update the estimates
                y_hat(i) = y(i)-e;
            end
        end
        
    end
    
    methods (Access = private)
        
        function [A,B,C,estimationError,thetaNew,P,phiMemory,psiMemory] = ...
                idrarmaxstep_simplified(this, uMeas,yMeas,isEnabled,adg1,adg2,theta,P,phiMemory,psiMemory,algorithmParams)
            % Recursive estimation of ARMAX models
            % Copied from the Matlab System Identification toolbox
            %
            %   Inputs:
            %     uMeas           - u(t), scalar, representing system inputs
            %     yMeas           - y(t), scalar, representing the system outputs
            %     isEnabled       - Scalar, 1: adaptation enabled 0: hold the coefficients,
            %     adg             - Forgetting factor, noise covariance matrix or the
            %                       adaptation gain, based on the estimation method.
            %                       Signal dimensions depend on the same setting
            %     theta           - n-by-1 vector, contains estimated coefficients
            %                       before the update
            %     P:              - n-by-n matrix, parameter covariance matrix before
            %                       the update
            %     algorithmParams - Structure, contains the necessary parameters from
            %                       the block dialog
            %
            %   Outputs:
            %     A               - The estimated A(q) polynomial
            %     B               - The estimated B(q) polynomial
            %     C               - The estimated C(q) polynomial
            %     estimationError - Scalar, representing the difference between the
            %                       measured output and the estimate
            %     thetaNew        - Internal state: Vectorized version of the A, B and
            %                       C outputs
            %     PNew            - n-ny-n matrix, parameter covariance matrix after
            %                       the update
            %     phiMemoryNew    - Vector containing the IO data necessary for the
            %                       current and future steps of the estimation
            %     psiMemoryNew    - Vector containing the derivative of the estimation
            %                       error with respect to the estimated parameters
            
            % Update the regression data with the latest inputs
            if this.getNumberOfInputs(algorithmParams.nb)>0
                % Store the latest u(t), will be shifted to right recursively
                phiMemory(algorithmParams.b0PosInPhi) = uMeas;  
                
                % Vector containing the derivative of the estimation error with respect to the estimated parameters  
                psiMemory(algorithmParams.b0PosInPsi) = uMeas - ...
                    psiMemory(algorithmParams.b0PosInPsi+1:algorithmParams.b0PosInPsi+algorithmParams.nc) * ...
                    theta(algorithmParams.c1PosInTheta:algorithmParams.c1PosInTheta+algorithmParams.nc-1);
                
                % Testing the difference between psi and phi: results are
                % similar
                %                                 psiMemory(algorithmParams.b0PosInPsi) = uMeas;
            end
            
            % Update Step: Recursive PEM
            psi = psiMemory(algorithmParams.necessaryDataForPsi);
            phi = phiMemory(algorithmParams.necessaryDataForPhi);
            % Estimation error: yMeas - yEstimated
            estimationError = yMeas - phi*theta;
            if isEnabled
                % Adaptation is enabled: update parameters
                [thetaNew,P] = ...
                    this.update_parameters(estimationError,...
                    algorithmParams.estimationMethod,adg1,adg2,theta,psi,P);
            else
                % Adaptation is disabled: hold the previous states
                thetaNew = theta;
            end
            
            % Ensure that we have a stable C polynomial
            C = this.getParametersC(algorithmParams.na, algorithmParams.nb, algorithmParams.nc, thetaNew);
            if isEnabled && ~this.idRecursiveEstimationCheckBistritzCondition(C)
                % New C(q) estimate is unstable. Restore the old C(q)
                thetaNew(algorithmParams.c1PosInTheta:algorithmParams.c1PosInTheta+algorithmParams.nc-1) = ...
                    theta(algorithmParams.c1PosInTheta:algorithmParams.c1PosInTheta+algorithmParams.nc-1);
                C = this.getParametersC( ...
                    algorithmParams.na, algorithmParams.nb, algorithmParams.nc, thetaNew);
            end
            
            epsilon = yMeas - phi*thetaNew;
            
            % Testing the difference between psi and phi: results are
            % similar
            %             psiUpdateA = -yMeas;
            %             psiUpdateC = epsilon;
            
            % Calculate the derivatives w.r.t. a1 and c1
            if algorithmParams.na>0
                psiUpdateA = C*[-yMeas; -psiMemory(1:algorithmParams.nc).'];
            else
                psiUpdateA = -yMeas;
            end
            psiUpdateC = ...
                C*[epsilon; -psiMemory(algorithmParams.c1PosInPsi:algorithmParams.c1PosInPsi+algorithmParams.nc-1).'];
            
            % Data storage for the next time-step
            phiMemory = this.shiftRowVectorToRight(phiMemory);
            psiMemory = this.shiftRowVectorToRight(psiMemory);
            if algorithmParams.na>0
                phiMemory(1) = -yMeas;
                psiMemory(1) = psiUpdateA;
            end
            phiMemory(algorithmParams.c1PosInPhi) = epsilon;
            psiMemory(algorithmParams.c1PosInPsi) = psiUpdateC;
            
            % Block Outputs
            A = this.getParametersA(algorithmParams.na,thetaNew);
            B = this.getParametersB(algorithmParams.na, ...
                algorithmParams.nb,algorithmParams.nk,thetaNew);
            % The C(q) output is calculated right after the update step
        end
        
        function [theta,P] = update_parameters(this,estimationError,...
                estimationMethod,adg1,adg2,theta,phi,P)
            % Update the parameters with a symmetric P
            %
            % Revised from the Matlab System Identification toolbox
            %
            % methods: forgetting factor, kalman filter, normalized gradient/gradient
            %
            % Inputs:
            %   estimationError:  yMeasured[k] - yEstimated[k]
            %   estimationMethod: (uint8) update algorithm. 1: forgetting
            %                     factor, 2: Kalman filter, 3: normalized
            %                     gradient, 4: gradient
            %   adg1:             ForgettingFactor, ProcessNoiseCovariance
            %                     or AdaptationGain, based on estimationMethod
            %   adg2:             NormalizationBias (only used by the
            %                     NormalizedGradient method)
            %   theta:            Current parameter values at time k
            %   phi:              [1 n] regressors/gradient vector
            %   P:                Parameter estimation covariance
            
            % Equations from pg. 367-369, L. Ljung, System
            % Identification: Theory for the User, 1999
            switch estimationMethod
                case 1 % Forgetting-factor: adg=forgetting factor
                    K = P*phi.' / (adg1(1) + phi*P*phi.');
                    P = (P-K*phi*P) / adg1(1);
                                        
                case 2 % Kalman filter: adg= noise covariance matrix
                    K = P*phi.' / (1+phi*P*phi.'); % Eq. (11.29b), R2=1
                    P = P-K*phi*P + adg1; % Eq. (11.29c), R1=adg

                case 3 % Normalized gradient: adg=adaptation gain
                    K = adg1(1)*phi.' / (adg2(1)+phi*phi.');
                    
                case 4 % Gradient: adg=adaptation gain
                    K = adg1(1)*phi.';
                    
                otherwise
                    % This should never happen... But it makes the code more robust
                    % and keeps the MATLAB Coder happy
                    K = zeros(size(theta),'like',theta);
            end
            
            theta = theta+K*estimationError;
        end
       
        function [InitialParameters]=setInitialParameters(this)
            % Initialization of the parameters
            % Copied from the Matlab System Identification toolbox
            
            % Initialization
            InitialParameters=zeros(this.na+sum(this.nb,2,'native')+...
                this.nc+this.nd+sum(this.nf,2,'native'),1);
            % A
            InitialParameters(1:this.na) = this.InitialA(2:1+this.na);
            % B
            InitialParameters(this.na+(1:sum(this.nb,2,'native'))) = this.InitialB;
            % C
            InitialParameters(this.na+sum(this.nb,2,'native')+(1:this.nc)) = this.InitialC(2:1+this.nc);
            % D
            if(this.nd>0)
                InitialParameters(this.na+sum(this.nb,2,'native')+this.nc+(1:this.nd)) = this.InitialD;
            end
            % F
            if(sum(this.nf,2,'native')>0)
                InitialParameters(this.na+sum(this.nb,2,'native')+this.nc+this.nd+(1:sum(this.nf,2,'native'))) = this.InitialF;
            end
        end
        
        function algorithmParams = getAlgorithmParameters(this,p)
            % getAlgorithmParameters Pre-calculate the required parameters
            % for the estimation algortihm. This is used by this system
            % object and the recursive polynomial model estimator simulink
            % block.
            %
            %   algorithmParams = getAlgorithmParameters(p)
            %
            %   Inputs;
            %     p - A structure containing validated parameters
            %
            % Copied from the Matlab System Identification toolbox
            
            % Initial derivatives of the estimation error w.r.t the
            % estimated parameters
            [~,nam,nbm] = ...
                this.numberOfRequiredPastJacobianEstimates(p.na,p.nb,p.nc,p.nk);
            
            % Data memory: phiMemory See localInitializeRarx() for detailed
            % explanation of phi. The elements in Phi that are required in
            % each time step
            % 
            % phiMemory contains more than necessary due to the existence
            % of system delay, nk.
            % necessaryDataForPhi=[y(t-1),...,y(t-na), u(t-nk),...,u(t-nk-nb+1}, w(t-1),...,w(t-nc)]
            %                     |<-... y meas...->|  |<-...  u meas    ...->|  |<-... noise ...->| 
            necessaryDataForPhi = [(1:p.na) ...
                (p.na+p.nk+(1:p.nb)) ...
                (p.na+p.nk+p.nb+(1:p.nc))];
            
            % Derivatives memory: phiMemory The elements in Psi that are
            % required in each time step
            %
            % psiMemory contains more than necessary due to the existence
            necessaryDataForPsi = [(1:p.na) ...
                (nam+p.nk+(1:p.nb)) ...
                (nam+nbm+(1:p.nc))];
            
            % Fcn outputs: Some of these are not needed for arma
            % estimation. However, the mex version of idrarmaxstep() requre
            % that algorithmParams have all these fields defined at all
            % times.
            algorithmParams.estimationMethod = uint8(p.estimationMethod);
            algorithmParams.na = uint16(p.na);
            algorithmParams.nb = uint16(p.nb);
            algorithmParams.nc = uint16(p.nc);
            algorithmParams.nd = uint16(p.nd);
            algorithmParams.nf = uint16(p.nf);
            algorithmParams.nk = uint16(p.nk);
            % Number of parameters
            algorithmParams.nParameters = uint16(p.na+sum(p.nb,2)+p.nc+p.nd+sum(p.nf,2));
            % Necessary data for phi
            algorithmParams.necessaryDataForPhi = uint16(necessaryDataForPhi);
            % Necessary data for psi
            algorithmParams.necessaryDataForPsi = uint16(necessaryDataForPsi);
            % Start position for B(q) parameters
            algorithmParams.b0PosInPhi = uint16(p.na+1);
            algorithmParams.b0PosInPsi = uint16(nam+1);
            % Start position for C(q) parameters
            algorithmParams.c1PosInTheta = uint16(p.na+sum(p.nb,2)+1);
            algorithmParams.c1PosInPhi = uint16(p.na+p.nb+p.nk+1);
            algorithmParams.c1PosInPsi = uint16(nam+nbm+1);
        end
        
        function adg = getAdaptationParameter(this)
            % Get the adaptation parameter for the corresponding estimation
            % method
            %
            % Copied from the Matlab System Identification toolbox
            
            if isEstimationMethodForgettingFactor(this)
                adg = this.ForgettingFactor;
            elseif isEstimationMethodKalmanFilter(this)
                adg = this.ProcessNoiseCovariance;
            elseif isEstimationMethodNormalizedGradient(this) || ...
                    isEstimationMethodGradient(this)
                adg = this.AdaptationGain;
            end
        end
        
        function adg = getAdaptationParameter2(this)
            % As of now only the NormalizedGradient method uses 2
            % parameters, which is the NormalizationBias
            %
            % Copied from the Matlab System Identification toolbox
            
            if isEstimationMethodNormalizedGradient(this)
                adg = this.NormalizationBias;
            else
                adg = cast(0,this.DataType);
            end
        end
        
        function tf = isEstimationMethodForgettingFactor(this)
            %
            % Copied from the Matlab System Identification toolbox
            
            if this.estimationMethod==1
                tf = true();
            else
                tf = false();
            end
        end
        
        function tf = isEstimationMethodKalmanFilter(this)
            %
            % Copied from the Matlab System Identification toolbox
            
            if this.estimationMethod==2
                tf = true();
            else
                tf = false();
            end
        end
        
        function tf = isEstimationMethodNormalizedGradient(this)
            %
            % Copied from the Matlab System Identification toolbox
            
            if this.estimationMethod==3
                tf = true();
            else
                tf = false();
            end
        end
        
        function tf = isEstimationMethodGradient(this)
            %
            % Copied from the Matlab System Identification toolbox
            
            if this.estimationMethod==4
                tf = true();
            else
                tf = false();
            end
        end
    end
    
    methods(Static)
        function y = ScalarToMatrix(u,n,isStrictPosDef)
            %   An n-by-n matrix y is created based on u. If u is a scalar, y has u
            %   on its diagonals. If u is a vector, y has the elements of u on its
            %   diagonals. If u is a matrix y = (u+u.')/2.
            %
            %   Inputs:
            %     u              - (double or single) Scalar, vector or matrix input
            %     n              - size(y,1)
            %     isStrictPosDef - Whether it is strictly positive definite
            %
            %   Output:
            %     y              - Covariance matrix
            %
            % Revised from the Matlab System Identification toolbox

            % If enforcing strict positive definiteness, map non-negative values to
            % eps. Otherwise, to 0.
            if isStrictPosDef
                zeroVariable = cast(eps,'like',u);
            else
                zeroVariable = cast(0,'like',u);
            end
            
            y = zeros(n,n,'like',u);
            if isscalar(u)
                if u(1)<=0
                    u(1) = zeroVariable;
                end
                for kk=uint32(1):uint32(n)
                    y(kk,kk) = u(1);
                end
            elseif isvector(u)
                for kk=uint32(1):uint32(n)
                    if u(kk)>0
                        y(kk,kk) = u(kk);
                    else
                        y(kk,kk) = zeroVariable;
                    end
                end
            elseif ismatrix(u)                
                y = (u + u.')/2;                
            end
         end
        
        function num = numberOfRequiredPastMeasurements(na,nb,nc,nk)
            % Number of required past measurements (y in A, u in B, and w in C)
            % Copied from the Matlab System Identification toolbox
            
            num = na+sum(nb+nk,2)+nc;
        end
        
        function [num,nam,nbm] = numberOfRequiredPastJacobianEstimates(na,nb,nc,nk)
            % Count the number of required past Jacobian estimates
            %
            % Copied from the Matlab System Identification toolbox
            
            if na>0
                nam = max(na,nc);
            else
                nam = 0;
            end
            if nb>0
                nbm = max(nb+nk,nc+1);
            else
                nbm = 0;
            end
            num = nam+nbm+nc;
        end
        
        function val = getNumberOfInputs(nb)
            % Calculate # of inputs in the estimated model based on nb
            % This function considers MI(Multi-Input) systems
            %
            % Copied from the Matlab System Identification toolbox
            
            if isscalar(nb) && nb(1)==0
                val = uint16(0); % time-series
            else
                val = uint16(size(nb,2));
            end
        end
        
        function A = getParametersA(na,allParams)
            % Get A from all params given the size of na
            % A = 1+a_1*q^(-1)+...+a_na*q^(-na)
            % The first term in A should be 1
            %
            % Copied from the Matlab System Identification toolbox
            
            A = zeros(1,na+1,'like',allParams);
            A(1) = 1;
            if na>0
                A(2:na+1) = allParams(1:na);
            end
        end
        
        function B = getParametersB(na,nb,nk,allParams)
            % Get B from all params given the size of na, nb, and nk
            % nk is the system delay
            % B = b_1 + b_2*q^(-1)+...+b_nb*q^(-nb+1)
            % This code can be used for MISO systems
            % But our ARMAX model is a SISO system
            %
            % Copied from the Matlab System Identification toolbox
            
            nInputs = recursive_armax.getNumberOfInputs(nb);
            if nInputs>0
                B = zeros(nInputs,max(nb+nk),'like',allParams);
                % When the delay term is not zero, the first nk term is zero
                cPos = na+1; % Current position (of param b0) in allParams
                for kk=1:nInputs % kk is input index
                    % Coder prefers the double for-loop for MISO: Cannot
                    % determine how to process the assignment
                    % B(kk,nk(kk)+1:nk(kk)+nb(kk)) = allParams(cPos:cPos+nb(kk)-1);
                    % cPos = cPos + nb(kk);
                    for kk2=1:nb(kk)
                        B(kk,nk(kk)+kk2) = allParams(cPos); % Shifted by nk
                        cPos = cPos + 1;
                    end
                end
            else
                % No input, return empty array
                B = cast([],'like',allParams);
            end
        end
        
        function C = getParametersC(na,nb,nc,allParams)
            % Get C from all params given the size of nc
            % C = 1+c_1*q^(-1)+...+c_nc*q^(-na)
            % The first term in C should be 1
            %
            % Copied from the Matlab System Identification toolbox
            
            C = zeros(1,nc+1,'like',allParams);
            C(1) = 1;
            if nc>0
                c1PosInAllParams = na+sum(nb,2)+1;
                C(2:nc+1) = allParams(c1PosInAllParams:c1PosInAllParams+nc-1);
            end
        end
        
        function v = shiftRowVectorToRight(v)
            % Shift the row vector to the right by one step
            % v(2:end) = v(1:end-1);
            %
            % Copied from the Matlab System Identification toolbox
            
            for kk=numel(v):-1:2
                v(kk) = v(kk-1);
            end
        end
        
        function isStable = idRecursiveEstimationCheckBistritzCondition(p)
            % idRecursiveEstimationCheckBistritzCondition Test if all roots of the
            % polynomial with coefficients p are in the unit disk. p must be a
            % row vector of real numbers. p(1) must be positive.
            %
            %   isStable = idRecursiveEstimationCheckBistritzCondition(p)
            %
            %   Algorithm: Yuval Bistritz, Zero Location With Respect To The Unit
            %   Circle Of Discrete-Time Linear System Polynomials, 1984 IEEE
            %
            % Copied from the Matlab System Identification toolbox
            
            %% Check simple, necessary conditions
            if isscalar(p)
                % Static gain
                isStable = true();
                return;
            end
            
            % Necessary conditions with the p(1)>1 assumption (pg. 1134)
            if p(1)<=abs(p(end)) || sum(p,2)<=0
                isStable=false();
                return;
            end
            
            %% Initizalization
            
            % pr: reciprocal polynomial of p
            pr = p(end:-1:1);
            
            % Eq. (12a)
            T_n = p + pr;
            if T_n(1)==0 % sum(p,2)<=0 implies sum(T_n,2)>0 necessary&sufficient condition
                isStable = false();
                return;
            end
            
            % T_nm1=deconv(Tp,[1 -1]) or Tp/(z-1) (Eq. (12b))
            nt = uint16(numel(p));
            T_nm1 = zeros(1,nt-1,'like',p);
            T_nm1(1) = p(1)-pr(1);
            for kk=2:nt-1
                T_nm1(kk) = p(kk)-pr(kk)+T_nm1(kk-1);
            end
            
            % Stability:
            % Necessary condition Eq. (14)
            % Necessary and sufficient condition Eq. (31)
            if T_nm1(1)==0 || sum(T_nm1,2)<=0;
                isStable = false();
                return;
            end
            
            %% Recursion step
            for kk=nt-1:-1:2
                delta = T_n(1)/T_nm1(1); % Eq. (25)
                % delta>0 is a necessary stability condition but skip cheking it
                for mm=1:kk-1
                    T_n(mm) = T_nm1(mm);
                    T_nm1(mm) = delta*(T_nm1(mm)+T_nm1(mm+1))-T_n(mm+1); % Eq. (26)
                end
                T_n(kk) = T_nm1(kk);
                % Stability:
                % Necessary condition Eq. (14)
                % Necessary and sufficient condition Eq. (31)
                sigmaI = T_nm1(1);
                if kk>=2
                    for mm=2:kk-1
                        sigmaI = sigmaI +  T_nm1(mm);
                    end
                end
                if T_nm1(1)==0 || sigmaI<=0; % sign change
                    isStable = false();
                    return;
                end
            end
            isStable = true();
        end
    end
end