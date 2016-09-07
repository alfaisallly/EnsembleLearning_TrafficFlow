classdef TrafficPredictionSystemTester < handle
    
    properties
        
        experiment              % 'complete_sensing' or 'decimation'
        mngr                    % job manager for the experiment
        model_based_params      % parameters of the model based system
        
        model_less_results_files
        model_based_results_files
        
        % experimental results for the model-less and model-based systems.
        % these are cell arrays, each cell contains results for one
        % experiment.
        E_mb
        E_ml
        
    end
    
    methods( Access = public )
        
        function [this] = TrafficPredictionSystemTester(configname,input_folder,output_folder,model_based_params,experiment)
            
            this.model_based_params = model_based_params;
            
            if nargin<5
                this.experiment = 'complete_sensing';
            else
                this.experiment = experiment;
            end
            
            % create manager for a complete sensing experiment
            switch experiment
                case 'complete_sensing'
                    this.mngr = ObjectFactory.complete_sensing_experiment(configname,input_folder,output_folder);
                case 'decimation'
                    this.mngr = ObjectFactory.sensor_decimation_experiment(configname,input_folder,output_folder);
            end
            
        end
        
        function [this] = run_experiment(this,execution)
            
            if nargin<2
                execution = 'series';
            end
            
            model_less_params = struct( ...
                'predictor_class',      'ModelLessPredictor' ,     ...
                'prediction_horizon',   3600 ,                     ...
                'update_dt',            1800 ,                     ...
                'start_time',           0 ,                        ...
                'end_time',             86400 );
            
            % run experiment on model-based predictor
            switch execution
                case 'series'
                    this.model_less_results_files  = this.mngr.execute_job_series( model_less_params );
                    this.model_based_results_files = this.mngr.execute_job_series( this.model_based_params );
                case 'parallel'
                    this.model_less_results_files  = this.mngr.execute_job_parallel_build( model_less_params );
                    this.model_based_results_files = this.mngr.execute_job_parallel_build( this.model_based_params );
            end
            
        end
        
        function [this] = process_results(this)
            
            result_model_less = this.mngr.process_result(this.model_less_results_files);
            result_model_based = this.mngr.process_result(this.model_based_results_files);
            
            switch this.experiment
                case 'complete_sensing'
                    this.E_ml = TrafficPredictionSystemTester.process_complete_sensing({result_model_less.sixty_minute_prediction.flw_error});
                    this.E_mb = TrafficPredictionSystemTester.process_complete_sensing({result_model_based.sixty_minute_prediction.flw_error});
                case 'decimation'
                    this.E_ml = TrafficPredictionSystemTester.process_decimation(result_model_less);
                    this.E_mb = TrafficPredictionSystemTester.process_decimation(result_model_based);
            end
            
        end
        
        function [X] = evaluate_hypothesis_test(this,alpha,lambda)
            
            if nargin<3
                lambda = 0;
            end
            
            if nargin<2
                alpha = 0.05;
            end
            
            % t test
            X = this.my_ttest_all(alpha,lambda);
            
        end
        
        function [lambda_thresh] = find_threshold_lambda(this,alpha)
            
            if nargin<2
                alpha = 0.05;
            end
            
            n1 = size(this.E_mb,1);
            n2 = size(this.E_mb,2);
            
            lambda_thresh = nan(n1,n2);
            
            for i=1:n1
                for j=1:n2
                    
                    e_ml = this.E_ml{i,j};
                    e_mb = this.E_mb{i,j};

                    if isempty(e_ml) || isempty(e_mb)
                        continue
                    end
                    
                    % lower bound lambda = -1
                    lambda_lower = -1;
                    X_lower = TrafficPredictionSystemTester.my_ttest(e_ml,e_mb,alpha,lambda_lower);

                    if X_lower.h~=0
                        error('this should not happen')
                    end

                    % find upper bound lambda
                    lambda_upper = 1;
                    X_upper = TrafficPredictionSystemTester.my_ttest(e_ml,e_mb,alpha,lambda_upper);
                    while X_upper.h~=1
                        lambda_upper = 2*lambda_upper;
                        X_upper = TrafficPredictionSystemTester.my_ttest(e_ml,e_mb,alpha,lambda_upper);
                    end

                    % bisection method
                    d = lambda_upper-lambda_lower;
                    while d > 1e-3
                        lambda_mid = mean([lambda_upper lambda_lower]);
                        X_mid = TrafficPredictionSystemTester.my_ttest(e_ml,e_mb,alpha,lambda_mid);
                        if X_mid.h==0
                            lambda_lower = lambda_mid;
                        else
                            lambda_upper = lambda_mid;
                        end
                        d = lambda_upper-lambda_lower;
                    end
                    
                    lambda_thresh(i,j) = lambda_mid;
                    
                end
            end

        end
        
        function [] = drawhist(d,str,xlim)
            figure
            hist(Utils.column_vector(d),30);
            grid
            h=get(gca,'Children');
            h.FaceColor = 0.5*[1 1 1];
            xlabel('error value')
            set(gca,'YTickLabel',{});
            set(gca,'XLim',xlim)
            textpos(0.7,0.9,0,str,16);
        end
    end

    methods( Access=private )
        
        function [X] = my_ttest_all(this,alpha,lambda)

            n1 = size(this.E_mb,1);
            n2 = size(this.E_mb,2);
            X = repmat( struct( 'mean_mb',nan,...
                                'std_mb',nan,...
                                'mean_ml',nan,...
                                'std_ml',nan,...
                                'mean_ml_lambda',nan,...
                                'h',nan,...
                                's',nan,...
                                'ci',[],...
                                'tstat',nan,...
                                'df',nan,...
                                'sd',nan ) , n1,n2 );

            for i = 1:n1
                for j = 1:n2
                    if ~isempty(this.E_ml{i,j}) && ~isempty(this.E_mb{i,j})
                        X(i,j) = TrafficPredictionSystemTester.my_ttest(this.E_ml{i,j},this.E_mb{i,j},alpha,lambda);
                    end
                end
            end
            
        end

    end
    
    methods( Static, Access = private )
        
        function [X] = my_ttest(e_ml,e_mb,alpha,lambda)
            
            % H0: e_mb=e_ml+mean(e_ml)*lambda vs H1: e_mb<e_ml+mean(e_ml)*lambda
            % h=0: no conclusion
            % h=1: H0 rejected in favor of H1
            [h,s,ci,stats] = ttest(e_mb,e_ml + mean(e_ml)*lambda,'Tail','left','Alpha',alpha);
            
            X.mean_mb = mean(e_mb);
            X.std_mb = std(e_mb);
            X.mean_ml = mean(e_ml);
            X.std_ml = std(e_ml);
            X.mean_ml_lambda = mean(e_ml)*(1+lambda);
            X.h = h;
            X.s = s;
            X.ci = ci;
            X.tstat = stats.tstat;
            X.df = stats.df;
            X.sd = stats.sd;
            
        end
        
        function [x] = process_complete_sensing(d)
            aaa=cellfun(@(x)mean(x,2,'omitnan'),d,'UniformOutput',false);
            x=horzcat(aaa{:});
            x = mean(x,1,'omitnan');
            x = {x};
        end
        
        function [x] = process_decimation(d)
            x = d.error_grid;
        end
        
    end
    
end

