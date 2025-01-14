function [P_testPosit] = Prob_testPosit_VL(VL,minP,range)

    P_testPosit = (VL-0.01)./(max(VL)-0.01)*range+minP;

end