%=====================================================================
%after running the part of the script above, you could save base and actual
%base contains 
clear 
load example.mat

for t = 1:length(predictionFrom)
    
    [ y_stackensemble ] = ...
        StackEnsembleFlowPrediction.compute_ensemble_only( actual, base, num_basemethods , num_days_for_stack);
    
end
