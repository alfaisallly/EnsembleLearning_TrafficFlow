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