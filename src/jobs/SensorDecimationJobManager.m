classdef SensorDecimationJobManager < PredictionJobManager
    
    properties
        good_ml_vds_for_day
        num_good_ml_vds_for_day
        removable_rp_vds
        num_samples_grid
    end
    
    methods( Access = public )
        
        function [this] = SensorDecimationJobManager(confname,infolder,outfolder)
            this = this@PredictionJobManager(confname,infolder,outfolder);
        end
        
        function [this] = plot_sample_distribution(this)
            figure
            imagesc(this.num_samples_grid)
            xlabel('num decimated ml')
            ylabel('num decimated rp')
            title(sprintf('Total samples = %d',sum(sum(this.num_samples_grid))))
        end
        
        function [result] = process_result(this,results_files)
            
            % Create a data provider for calculating error
            error_data_provider = Utils.get_pems_dp(this.ni,this.configname);
            
            % loop through results
            dim_ml = size(this.num_samples_grid,1);
            dim_rp = size(this.num_samples_grid,2);
            error_grid = cell(dim_ml,dim_rp);
            for i=1:length(results_files)
                r = load(fullfile(this.outputfolder,results_files{i}));
                
                % compute error
                [~,~,one_hour_flow_error] = this.process_log_prediction(error_data_provider,r.X.day,r.X.log);
                
                e = Utils.meanwithnan(Utils.meanwithnan(one_hour_flow_error,1),2);            
                
                % put it into the grid
                num_ml = sum(ismember(r.X.remove_vdss,this.all_ml_vds));
                num_rp = sum(ismember(r.X.remove_vdss,this.all_rp_vds));
                error_grid{num_ml+1,num_rp+1} = [error_grid{num_ml+1,num_rp+1} e];
            end
            
            % statistics on each cell
            error_mean = nan(dim_ml,dim_rp);
            error_var  = nan(dim_ml,dim_rp);
            for i=1:dim_ml
                for j=1:dim_rp
                    error_mean(i,j) = mean(error_grid{i,j});
                    error_var(i,j) = var(error_grid{i,j});
                end
            end
            
            % save to result
            result = struct( ...
                'error_mean' , error_mean , ...
                'error_var' , error_var );
            result.error_grid = error_grid;
            
        end
        
        function [this] = report_result(this,output_file)
            
        end
        
    end
    
    methods( Access = protected )

        function [this] = initialize(this)
            
            % removable ml and ramp vdss
            % remove only vdss that do not already have a template
            removable_vdss = setdiff( this.ni.vds2linkid(1,:) , this.config.vds_use_template_ids );
            this.removable_rp_vds = this.all_rp_vds(ismember(this.all_rp_vds,removable_vdss));
            
            % find removeable internal vdss for each day
            % there are the good mainline vdss on each day		
			dp_pems = Utils.get_pems_dp(this.ni,this.configname);

            fprintf('Getting good mainline sensors for each day\n');
            for d=1:length(this.days)
                this.good_ml_vds_for_day{d} = dp_pems.get_good_internal_sensors(this.ni,this.days(d));
            end
            this.num_good_ml_vds_for_day = cellfun(@(x) numel(x),this.good_ml_vds_for_day);
            
            % initialize the random sampler
            rng('shuffle', 'multFibonacci');
            
        end

        function [this] = build_job(this,params)
            
            if nargin<2
                params = struct();
            end
            
            % minimum number of samples in a cell
            if isfield(params,'min_samp')
                min_samples = params.min_samp;
            else
                min_samples = 10;
            end
            
            % maximum number of samples in a cell
            if isfield(params,'max_samp')
                max_samples = params.max_samp;
            else
                max_samples = 1000;
            end
            
            % total number of samples to take in all cells
            if isfield(params,'tot_samp')
                total_samples = params.tot_samp;
            else
                total_samples = 20000;
            end
            
            % probability of failure for a single sensor
            if isfield(params,'pfail')
                prob_fail = params.pfail;
            else
                prob_fail = 0.1;
            end
            
            if max_samples<min_samples
                error('max_samples<min_samples')
            end
            
            % generate parameters for all runs
            % largest number of ml vds that can be removed; need at least two ml vds for linear interpolation
            max_removable_ml_vds = max(this.num_good_ml_vds_for_day) - 2; 
            num_removable_rp_vds = numel(this.removable_rp_vds);
            
            num_uniform_samples = (1+max_removable_ml_vds)*(1+num_removable_rp_vds)*min_samples;
            num_distribution_samples = total_samples-num_uniform_samples;
            
            if num_distribution_samples<0
                error('num_distribution_samples<0')
            end
            
            % initialize all_tasks and num_samples_grid
%             all_tasks = repmat(struct('day',nan,'remove_vdss',[]),total_samples,1);
            all_tasks = repmat(struct('day',nan,'remove_ml_vdss',[],'remove_rp_vdss',[]),total_samples,1);

            this.num_samples_grid = zeros(max_removable_ml_vds+1,num_removable_rp_vds+1);
            
            % first pass, put min_samples in each cell
            fprintf('First pass: generating %d runs by uniform sampling\n',num_uniform_samples);
            c = 1;
            for num_ml_remove=0:max_removable_ml_vds
                for num_rp_remove=0:num_removable_rp_vds
                    for i=1:min_samples
                        all_tasks(c) = this.uniform_sampling(num_ml_remove,num_rp_remove);
                        c = c+1;
                    end
                end
            end
            
            % second pass, continue sampling according to distribution
            fprintf('Second pass: generating %d runs by binomial sampling\n',num_distribution_samples);
            for i=1:num_distribution_samples
                all_tasks(c) = this.binomial_sampling(prob_fail);
                c = c+1;
            end
            
            % write to file
            this.tasks_file = this.save_job(all_tasks);
            
        end    
    
    end

    methods( Access = private )
        
        function [X] = uniform_sampling(this,num_ml_remove,num_rp_remove)
            
            % pick a day with sufficient removable ml vdss
            eligible_days = this.num_good_ml_vds_for_day >= num_ml_remove+2; %at least two ml vdss left
            day = datasample(this.days(eligible_days),1);
            today_removable_ml_vds = this.good_ml_vds_for_day{this.days==day};
            
            % select at random the vdss to remove
            if num_ml_remove>0
                decimate_ml = datasample(today_removable_ml_vds,num_ml_remove,'Replace',false);
            else
                decimate_ml = [];
            end
            
            if num_rp_remove>0
                decimate_rp = datasample(this.removable_rp_vds,num_rp_remove,'Replace',false);
            else
                decimate_rp = [];
            end
            
            remove_vdss = [decimate_ml decimate_rp];
            
            % record number of mainline and ramp decimated vdss
            this.record_num_samples(remove_vdss);
            
            % save in struct
%             X = struct('day',day,'remove_vdss',remove_vdss);
            X = struct('day',day,'remove_ml_vdss',decimate_ml,'remove_rp_vdss',decimate_rp);
            
        end
        
        function [X] = binomial_sampling(this,prob_fail)
            
            % pick a day with sufficient removable ml vdss
            day = datasample(this.days,1);
            today_removable_ml_vds = this.good_ml_vds_for_day{this.days==day};
            
            % all removable vdss
            removable_vdss = [today_removable_ml_vds this.removable_rp_vds];
            
            % select a number of vdss to remove
            num_decimate = binornd(numel(removable_vdss),prob_fail);
            
            % select the vdss to remove, make sure we are not removing all of the mainline sensors.
            not_done = true;
            while not_done
                
                % choose num_decimate from vdss
                remove_vdss = datasample(removable_vdss,num_decimate,'Replace',false);
                
%                 % finish if you have not selected all removable mainline vdss
%                 not_done = isempty(setdiff(today_removable_ml_vds,remove_vdss));
                
                % finish if you have at least two ml vdss left
                not_done = (sum(ismember(today_removable_ml_vds,remove_vdss))>numel(today_removable_ml_vds)-2);
            end
            
            % record number of mainline and ramp decimated vdss
            this.record_num_samples(remove_vdss);
            
            % save in struct
%             X = struct('day',day,'remove_vdss',remove_vdss);
            decimate_ml = intersect(remove_vdss,today_removable_ml_vds);
            decimate_rp = intersect(remove_vdss,this.removable_rp_vds);
            X = struct('day',day,'remove_ml_vdss',decimate_ml,'remove_rp_vdss',decimate_rp);

        
        end
        
        function [this] = record_num_samples(this,remove_vdss)
            num_ml = sum(ismember(remove_vdss,this.all_ml_vds));
            num_rp = sum(ismember(remove_vdss,this.all_rp_vds));
            m = num_ml+1;
            r = num_rp+1;
            this.num_samples_grid(m,r) = this.num_samples_grid(m,r) + 1;
        end
           
    end
    
end
