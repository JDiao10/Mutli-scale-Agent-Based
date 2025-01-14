function [P_positive] = Prob_get_test(tau_list,ModelType,Ph,Pl,TimeStep)
    if ModelType ==1        
        P_positive(1:2/TimeStep) = 1- (1-Pl)^(TimeStep);
        P_positive(2/TimeStep+1:length(tau_list)) = 1-(1-Ph)^(TimeStep); 

    elseif ModelType ==2       
        P_positive(1:3/TimeStep) = 1- (1-Pl)^(TimeStep);
        P_positive(3/TimeStep+1:length(tau_list)) = 1-(1-Ph)^(TimeStep);
    end

end