classdef Utils
    
    methods(Static)

        %% plotting
        function []=plot_24hr_contour(A,qty,yticks)
            
            if nargin<2
                qty = 'dty';
            end
            
            switch qty
                case 'dty'
                    A = -A;
            end
            
            h=pcolor(A);
            set(h,'EdgeAlpha',0)
            colormap('bone')
            
            % set space axis tick labels
            if nargin>2
                set(gca,'YTick',1:length(yticks))
                set(gca,'YTickLabel',yticks)
            end
            
            % set time axis tick labels
            dt = 86400/size(A,2);
            [timetick,timeticklabel]=Utils.get_time_ticks(dt);
            set(gca,'XTick',timetick)
            set(gca,'XTickLabel',timeticklabel)
            
        end
        
        function []= plot_24_hour_line(A)
            
            dt = 86400/length(A);
            time = dt:dt:86400;            
            fillplot(time,A,{'b'},gcf);
            [timetick,timeticklabel]=Utils.get_time_ticks(dt);
            set(gca,'XTick',time(timetick))
            set(gca,'XTickLabel',timeticklabel)
            grid
        end
        
        function [h,m,s] = hist(data,xstr,fit,numbins,xlim)
            if nargin<3
                fit = '';
            end
            if nargin<4
                numbins = 20;
            end
            [m,s] = normfit(data);
            if ~isempty(fit)
                h=histfit(data,numbins,fit);
                set(h(1),'FaceAlpha',0.5);
                set(h(2),'Color','b');
            else
                hist(data,numbins)
            end
            grid
            xlabel(xstr);
            ylabel('frequency');
            if nargin>4
                set(gca,'XLim',xlim)
            end
            textpos(0.7,0.9,0,sprintf('mean = %.1f%%\nstddev = %.1f%%',m,s),11);
        end
        
        function [] = contour_plot(figvis,time,vdss,data,tit,op)
            
            figure('Visible',figvis)
            h=pcolor(time,1:length(vdss),data);
            set(h,'EdgeAlpha',0)
            timelabelXTick(gca)
            set(gca,'YTick',[])
            ylabel('\rightarrow','FontSize',20)
            title(tit)
            colorbar

            if ~isempty(op)
                addslide(op,tit,[],[],[],0.7)
                close
            end
                
        end
        
        function [] = demand_sr_plot(figvis,horizon,time,data,tit,op)
            
            figure('Visible',figvis)
            h=pcolor(time/3600,horizon,data);
            set(h,'EdgeAlpha',0)
            timelabelXTick(gca)
            xlabel('Time','FontSize',20)
            ylabel('Horizon (sec)','FontSize',20)
            title(tit)
            caxis([0 max([max(data),1])])
            colorbar

            if ~isempty(op)
                addslide(op,tit,[],[],[],0.7)
                close
            end
                
        end
        
        function [] = compare_boundary_flows_two_results(days,r1,name1,r2,name2,output_file)
            
            % escape if there are no boundary flows
            if (~isfield(r1,'predicted_boundary_flow')||~isfield(r2,'predicted_boundary_flow'))
                return
            end
            
            save_ppt = nargin>5 && ~isempty(output_file);
            
            if save_ppt
                figvis = 'off';
                [ppt,op]=openppt([output_file '_boundary_flow'],true);
            else
                figvis = 'on';
                ppt = [];
                op = [];
            end
            
            for i=1:length(days)
                
                day = days(i);
                
                if save_ppt
                    addslideTitle(op,{datestr(day)})
                end
                
                % plots
                num_ramp=size(r1.predicted_boundary_flow(i).boundary_flow,2);
                for j=1:num_ramp
                    Utils.d2( ...
                        figvis , ...
                        r1.predicted_boundary_flow(i).boundary_flow(j).horizon, ...
                        r1.predicted_boundary_flow(i).boundary_flow(j).time , ...
                        r1.predicted_boundary_flow(i).boundary_flow(j).flw_vph , name1 , ...
                        r2.predicted_boundary_flow(i).boundary_flow(j).flw_vph , name2 , ...
                        ['Date: ' datestr(day) ' & Link ID: ' num2str(r1.predicted_boundary_flow(i).boundary_flow(j).link_id) ' & Type: ' char(r1.predicted_boundary_flow(i).boundary_flow(j).link_type) ] , ...
                        op );
                end
            end
            
            if save_ppt
                closeppt(ppt,op)
            end
            
        end
        
        function [] = d2(figvis,horizon,time,d1,name1,d2,name2,tit,op)
            
            figure('Visible',figvis)

            %plot simulation data
            for i=1:length(time)
                tmpTime=(time(i)+horizon)/3600;
                [h1]=plot(tmpTime,d1(:,i),'--r');
                hold on;
            end
            
            %plot PeMS data
            for i=1:length(time)
                tmpTime=(time(i)+horizon)/3600;
                [h2]=plot(tmpTime,d2(:,i),'b');
                hold on;
            end
            
            timelabelXTick(gca)
            xlabel('Time','FontSize',20)
            ylabel('Flow-rate (vph)','FontSize',20)
            legend([h1 h2],name1,name2)
            title(tit)

            if ~isempty(op)
                addslide(op,tit,[],[],[],0.7)
                close
            end
                
        end
        
        %% other
        
        function [y] = row_vector(x)
            
            if isempty(x)
                y = x;
                return
            end
            minsize = min(size(x));
            maxsize = max(size(x));
            y=reshape(x,1,maxsize*minsize);
        end
        
        function [y] = column_vector(x)
            
            if isempty(x)
                y = x;
                return
            end
            minsize = min(size(x));
            maxsize = max(size(x));
            y=reshape(x,maxsize*minsize,1);
        end
                
        function [y] = sum_cell_array(x)
            
            if isempty(x)
                y=nan;
                return
            end
            
            % all have the same length
            if ~all(cellfun(@(z)length(z)==length(x{1}),x))
                error('~all(cellfun(@(z)length(z)==n,x))')
            end
            
            y = x{1};
            for i=2:length(x)
                y = y + x{i};
            end
        end
        
        function [y] = stack_cell_array(x)
            y = [];
            for i=1:length(x)
                y = [y;x{i}];
            end
        end
        
        function [y] = cell2array(x)
            y = nan(size(x,1),size(x,2));
            for i=1:length(x)
                if length(x{i})~=1
                    error('length(x{i})~=1')
                else
                    y(i) = x{i};
                end
            end
        end
        
        function [y] = flatten_cell_array(x)
            y = [];
            for i=1:length(x)
                y = [y Utils.row_vector(x{i})];
            end
        end
        
        function [x] = isclass(obj,str)
            objclass = class(obj);
            
            switch str
                case 'ModelRunner'
                    x = strcmp(objclass,'ModelRunner') || strcmp(objclass,'BeatsModelRunner');
                otherwise
                    x = strcmp(objclass,str);
            end
        end
        
        function [e] = error(a,b)
            e.mae   = mean( mean(abs(a-b)) );
            e.mape  = mean( mean(abs((a-b)./b)) );
            e.smape = mean( mean( 2*abs(a-b)./(abs(a)+abs(b))) );
            e.rmse  = mean( sqrt(mean((a-b).^2)) );
            e.linf  = mean( max(abs(a-b)) );
        end
        
        function [a] = avg_error(e)
            a.mae = mean(e.mae);
            a.mape = mean(e.mape);
            a.smape = mean(e.smape);
            a.rmse = mean(e.rmse);
            a.linf = mean(e.linf);
        end
        
		function [A]=meanwithnan(B,dim)

			if(isempty(B))
				A=[];
				return
			end

			if(nargin==1)
				dim=1;
			end
				
			if(dim==2)
				B=B';
			end

            A = nan(1,size(B,2));
			for i=1:size(B,2)
				if ~all(isnan(B(:,i)))
				A(i) = mean(B(~isnan(B(:,i)),i));
				end
			end

			if(dim==2)
				A=A';
			end
        end
        
        function [A] = varwithnan(B,dim)
            
            if(isempty(B))
                A=[];
                return
            end
            
            if(nargin==1)
                dim=1;
            end
            
            if dim==2
                B=B';
            end
            
            A = nan(1,size(B,2));
            for i=1:size(B,2)
                if ~all(isnan(B(:,i)))
                    A(i) = var(B(~isnan(B(:,i)),i));
                end
            end
            
            if(dim==2)
                A=A';
            end
        end
        
		function [ x ] = index_into( A,B )
		% [~,x]=ismember(A,B);
			[~,x]=ismember(A,B);
		end
		
        %% getters
        function [dp] = get_pems_dp(ni,configname,cluster_method)
            
            if nargin<3
                cluster_method = '';
            end
            
            vdss = ni.get_all_vds;
            link_ids = ni.get_link_ids_for_vds(vdss);
            link_length = ni.get_lengths_for_linkids_in_km(link_ids);
            config = Config.get(configname);
            dp = ObjectFactory.pems_data_provider(config.good_days, ...
                vdss, ...
                config.vds_use_template, ...
                configname, ...
                link_ids, ...
                link_length , ...
				cluster_method );
        end
        
        
    end
    
    methods(Static, Access=private)

        function [timetick,timeticklabel]=get_time_ticks(dt)
            numtick = 13;
            time = 0:dt:86400-dt;
            numtime = length(time);
            timetick = unique(ceil(linspace(1,numtime,min(numtick,numtime))));
            timeticklabel = datestr(time(timetick)/86400,'HH:MM');
        end
        
    end

end

