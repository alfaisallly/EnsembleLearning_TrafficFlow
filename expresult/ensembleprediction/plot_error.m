

predictionFrom  = [7,8,9,10,11,12,13,14,15,16,17,18,19] ;
predictionTo = [8,9,10,11,12,13,14,15,16,17,18,19,20] ;
%=====================================================
% compute error
% MAPE
pls_MAPE = zeros(length(predictionFrom) , 1);
armax_MAPE = zeros(length(predictionFrom) , 1);
svm_MAPE = zeros(length(predictionFrom) , 1);
gpr_MAPE = zeros(length(predictionFrom) , 1);
krr_MAPE = zeros(length(predictionFrom) , 1);
ensemble_MAPE = zeros(length(predictionFrom) , 1);
%=====================================================
% compute error
% RMSE
pls_RMSE= zeros(length(predictionFrom) , 1);
armax_RMSE = zeros(length(predictionFrom) , 1);
svm_RMSE = zeros(length(predictionFrom) , 1);
gpr_RMSE = zeros(length(predictionFrom) , 1);
krr_RMSE = zeros(length(predictionFrom) , 1);
ensemble_RMSE = zeros(length(predictionFrom) , 1);
%=====================================================

y_actual = allday_actual(num_days_for_stack+1:end,:,:);
y_base = allday_base(num_days_for_stack+1:end,:,:,:);

for t = 1:length(predictionFrom)
    nz_day=0;
    for d = 1: size(y_actual,1)
        if norm( y_actual(d,:,t) ,1) > 0
            nz_day = nz_day + 1;
            %===================================
            % MAPE
            pls_MAPE(t) = pls_MAPE(t) + norm(y_base(d,:,1,t) - y_actual(d,:,t) ,1)/norm( y_actual(d,:,t) ,1) ;
            armax_MAPE(t) = armax_MAPE(t) + norm(y_base(d,:,2,t) - y_actual(d,:,t) ,1)/norm( y_actual(d,:,t) ,1) ;
            svm_MAPE(t) = svm_MAPE(t) + norm(y_base(d,:,3,t) - y_actual(d,:,t) ,1)/norm( y_actual(d,:,t) ,1) ;
            gpr_MAPE(t) = gpr_MAPE(t) + norm(y_base(d,:,4,t) - y_actual(d,:,t) ,1)/norm( y_actual(d,:,t) ,1) ;
            krr_MAPE(t) = krr_MAPE(t) + norm(y_base(d,:,5,t) - y_actual(d,:,t) ,1)/norm( y_actual(d,:,t) ,1) ;
            
            ensemble_MAPE(t) = ensemble_MAPE(t) + norm(y_stackensemble(d,:,t) - y_actual(d,:,t) ,1)/norm( y_actual(d,:,t) ,1);
            %===================================
            % RMSE
            pls_RMSE(t) = pls_RMSE(t) + norm(y_base(d,:,1,t) - y_actual(d,:,t) ,2) ;
            armax_RMSE(t) = armax_RMSE(t) + norm(y_base(d,:,2,t) - y_actual(d,:,t) ,2) ;
            svm_RMSE(t) = svm_RMSE(t) + norm(y_base(d,:,3,t) - y_actual(d,:,t) ,2) ;
            gpr_RMSE(t) = gpr_RMSE(t) + norm(y_base(d,:,4,t) - y_actual(d,:,t) ,2) ;
            krr_RMSE(t) = krr_RMSE(t) + norm(y_base(d,:,5,t) - y_actual(d,:,t) ,2) ;
            
            ensemble_RMSE(t) = ensemble_RMSE(t) + norm(y_stackensemble(d,:,t) - y_actual(d,:,t) ,2);
        end
    end
    pls_MAPE(t) = pls_MAPE(t) / nz_day ;
    armax_MAPE(t) = armax_MAPE(t) / nz_day ;
    svm_MAPE(t) = svm_MAPE(t) /nz_day ;
    gpr_MAPE(t) = gpr_MAPE(t) / nz_day ;
    krr_MAPE(t) = krr_MAPE(t) / nz_day ;
    
    ensemble_MAPE(t) = ensemble_MAPE(t) / nz_day ;
    
    pls_RMSE(t) = pls_RMSE(t) / nz_day ;
    armax_RMSE(t) = armax_RMSE(t) / nz_day ;
    svm_RMSE(t) = svm_RMSE(t) /nz_day ;
    gpr_RMSE(t) = gpr_RMSE(t) / nz_day ;
    krr_RMSE(t) = krr_RMSE(t) / nz_day ;
    
    ensemble_RMSE(t) = ensemble_RMSE(t) / nz_day ;
    
    
end

pls_RMSE(t) = sqrt(pls_RMSE(t)) ;
armax_RMSE(t) = sqrt(armax_RMSE(t) );
svm_RMSE(t) = sqrt(svm_RMSE(t)) ;
gpr_RMSE(t) = sqrt(gpr_RMSE(t)) ;
krr_RMSE(t) = sqrt(krr_RMSE(t)) ;

ensemble_RMSE(t) = sqrt(ensemble_RMSE(t)) ;

%====================================================================

figure(1)
F = plot( predictionFrom , pls_MAPE , '.-' , ...
    predictionFrom , armax_MAPE , '*-' , ...
    predictionFrom , svm_MAPE , 'o-' , ...
    predictionFrom , gpr_MAPE , 'd-' , ...
    predictionFrom , krr_MAPE , 's-' , ...
    predictionFrom ,ensemble_MAPE, '^-'  ) ;

set(F,'LineWidth',2)
LG = legend('pls','armax','svm','gaussian process','kernel ridge','stacked ensemble')
set(LG,'FontSize',13);

xl = xlabel('Time of the day for prediction');
yl= ylabel('MAPE') ;
set(gca,'FontSize',15)
%====================================================================
figure(2)
F = plot( predictionFrom , pls_RMSE , '.-' , ...
    predictionFrom , armax_RMSE , '*-' , ...
    predictionFrom , svm_RMSE , 'o-' , ...
    predictionFrom , gpr_RMSE , 'd-' , ...
    predictionFrom , krr_RMSE , 's-' , ...
    predictionFrom ,ensemble_RMSE, '^-'  ) ;

set(F,'LineWidth',2)
LG = legend('pls','armax','svm','gaussian process','kernel ridge','stacked ensemble')
set(LG,'FontSize',13);

xl = xlabel('Time of the day for prediction');
yl= ylabel('RMSE') ;
set(gca,'FontSize',15)
% %====================================================================
