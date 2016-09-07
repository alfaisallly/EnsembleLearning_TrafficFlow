
stackweight_linearizebytime = zeros(size(stackweight,3)*size(stackweight,4),size(stackweight,2),size(stackweight,1));


for m = 1: num_basemethods
    
   for d = 1:size(stackweight,3)
       
        for h = 1: size(stackweight,4)
            
            for i = 1:4 
                        
                stackweight_linearizebytime((d-1)*size(stackweight,4) + h, i , m ) = stackweight(m,i,d,h) ; 
           
            end 
        end 
        
   end
      
end   
clear i m d
%==================================
% ploting options
% which method: 1.armax 2.pls 3.svm 4.gaussian process 5. kernel ridge
basemethod = 5
% which interval of the hour
intervalofhour = 1
% showing how many days
daystoplot = 7
%==================================
F=plot(stackweight_linearizebytime(1:13*daystoplot,intervalofhour,basemethod),'b-')
set(F,'LineWidth',2)
set(gca,'FontSize',15)
xlim([1,13*daystoplot])
ax = gca;
ax.XTick = [0:1:daystoplot-1]*13 + ones(1,daystoplot);
markerstring = cell(daystoplot,1);
for d = 1:daystoplot
   markerstring{d}=strjoin({'day',num2str(d)});
end    
  
ax.XTickLabel = markerstring;
title('stack coefficients for kernel Ridge regression in 7 days')
