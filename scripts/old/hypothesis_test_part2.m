function []=hypothesis_test_part2()
% Part 2 of the hypothesis testing experiment.
% load the result, run ttest, report

% NOTE: REMEMBER TO TRANSPOSE FLW_ERROR

alpha = 0.05; % significance level
load('hypothesis_testing_experiment')
MBA = load('model_based_A');
MBB = load('model_based_B');
ML = load('model_less');

% checks
internal_vdss = mngr.ni.get_vds_for_link_ids(mngr.ni.get_internal_link_ids);
internal_vdss(isnan(internal_vdss)) = [];
num_days = numel(MBA.results.good_vdss);
for i=1:num_days
    if ~all(MBA.results.good_vdss(i).ids==ML.results.good_vdss(i).ids)
        error('~all(MBA.results.good_vdss==ML.results.good_vdss)')
    end
    if ~all(MBB.results.good_vdss(i).ids==ML.results.good_vdss(i).ids)
        error('~all(MBB.results.good_vdss==ML.results.good_vdss)')
    end
    if ~all(ismember(MBA.results.good_vdss(i).ids,internal_vdss))
        error('not all ids in the result are internal')
    end
    if ~all(ismember(MBB.results.good_vdss(i).ids,internal_vdss))
        error('not all ids in the result are internal')
    end
end

% error_all = horzcat(MBA.results.sixty_minute_prediction.flw_error);
% as_matrix = error_all(1:end-1,:)';
% as_vector = reshape(as_matrix,numel(as_matrix),1);
% as_vector(isnan(as_vector)) = [];
% figure, histfit(as_vector,20,'exponential')
% extract errors
e_mb_a = extract_data({MBA.results.sixty_minute_prediction.flw_error} );
e_mb_b = extract_data({MBB.results.sixty_minute_prediction.flw_error} );
e_ml   = extract_data({ML.results.sixty_minute_prediction.flw_error} );

drawhist(e_mb_a,'mbTPS-A',[0,1])
drawhist(e_mb_b,'mbTPS-B',[0,1])
drawhist(e_ml,'model-less',[0,1])
% extract errors
E_mb_a = mean(e_mb_a,1,'omitnan');
E_mb_b = mean(e_mb_b,1,'omitnan');
E_ml   = mean(e_ml,1,'omitnan');
drawhist(E_mb_a,'mbTPS-A',[0,0.25])
drawhist(E_mb_b,'mbTPS-B',[0,0.25])
drawhist(E_ml,'model-less',[0,0.25])
[mean(e_ml) std(e_ml)]
[mean(e_mb_a) std(e_mb_a)]
[mean(e_mb_b) std(e_mb_b)]
% H0: e_mb_a = e_ml vs H1: e_ml > e_mb_a
[p,s,ci,stats] = ttest(E_ml,E_mb_a,'Tail','right','Alpha',alpha);
[p,s]
ci
stats

% H0: e_mb_b = e_ml vs H1: e_ml > e_mb_b
[p,s,ci,stats] = ttest(E_ml,E_mb_b,'Tail','right','Alpha',alpha);
[p,s]
ci
stats

function [x] = extract_data(d)
aaa=cellfun(@(x)mean(x,2,'omitnan'),d,'UniformOutput',false);
x=horzcat(aaa{:});


function [x] = drawhist(d,str,xlim)
figure
hist(Utils.column_vector(d),30);
grid
h=get(gca,'Children');
h.FaceColor = 0.5*[1 1 1];
xlabel('error value')
set(gca,'YTickLabel',{});
set(gca,'XLim',xlim)
textpos(0.7,0.9,0,str,16);
