function [T01, T02] = calc_T0(set,temp,tempr,param)
% T0(t) is calculated via T0 = TAU * T_ with
%       [g   1   0   0 ...]
%       [g^2 g   1   0 ...]
% TAU = [g^3 g^2 g   1 ...]
%       [ .   .  .   . ...]
%       [ .   .  .   . ...]    
% g = 1 - 1/tau
% T_ = [T0(1), T(2)/tau, T(3)/tau, ...]
% if n is the size of the steps of our time series, then
% tau = 1 / (1 -g_original^n) with g_original = 1 - 1/tau_original
% 
% [T01 T02] = calc_T0(set,temp,tempr,param)
% 
% INPUT: 
% - set   -> general settings
% - temp  -> temperature draw, dependent on the temperature error model
% - tempr -> original temperature input
% - param -> the parameter set to test [a1, a2, c, tau1, tau2, To(1)]
% 
% OUTPUT: 
% - T01   -> equilibrium temperature 1
% - T02   -> equilibrium temperature 2, only if strcmp(model,'TwoTau')==1

    if set.OptimT0
        T0_rnd = param(6);
    else
        T0_rnd = 0;
    end
    tau1 = param(4);
    tau2 = param(5);
    
    nyrs = size(temp,1);
    
    % Calc T01 with long adjustment time scale, either with burnin in (UseMarT0 
    % -> two different year steps) or without
    if set.UseMarT0
        % number of years over which obs is avereaged (year steps of time series)
        yrs1 = tempr.T0yrs;
        yrs2 = tempr.yrs;
        tau1_1 = tau1/yrs1; 
        tau1_2 = tau1/yrs2;
        G1 = 1-1/tau1_1;
        G2 = 1-1/tau1_2;

        G1_ = (ones(1,nyrs)*G1).^[0:nyrs-1];
        G1_M = toeplitz(G1_,[G1_(1) zeros(1,size(G1_,2)-1)]);
        G2_ = (ones(1,nyrs)*G2).^[0:nyrs-1];
        G2_M = toeplitz(G2_,[G2_(1) zeros(1,size(G2_,2)-1)]);

        G_M1 = [G1_M(:,1:tempr.T0burnin) G2_M(:,tempr.T0burnin+1:end)];
        temp_1 = [temp(1:tempr.T0burnin,:)./tau1_1;temp(tempr.T0burnin+1:end,:)./tau1_2];
        temp_1(1,:) = mean(temp(tempr.T0temp(:,1)<=set.T0period(end),:),1)+ T0_rnd;
    else
        yrs = tempr.yrs;
        tau1 = tau1/yrs;
        G = 1-1/tau1;
        G_ = (ones(1,nyrs)*G).^[0:nyrs-1];
        G_M1 = toeplitz(G_,[G_(1) zeros(1,size(G_,2)-1)]);
        temp_1 = temp/tau1;
        temp_1(1,:) = mean(temp(tempr.T(:,1)<=set.T0period(end),:),1)+ T0_rnd;
    end
    
    % Calc T02 with short time scale
    if strcmp(set.model,'TwoTau')
        if set.UseMarT0
            tau2_1 = tau2/yrs1;
            tau2_2 = tau2/yrs2;
            G1 = 1-1/tau2_1;
            G2 = 1-1/tau2_2;

            G1_ = (ones(1,nyrs)*G1).^[0:nyrs-1];
            G1_M = toeplitz(G1_,[G1_(1) zeros(1,size(G1_,2)-1)]);
            G2_ = (ones(1,nyrs)*G2).^[0:nyrs-1];
            G2_M = toeplitz(G2_,[G2_(1) zeros(1,size(G2_,2)-1)]);

            G_M2 = [G1_M(:,1:tempr.T0burnin) G2_M(:,tempr.T0burnin+1:end)];
            temp_2 = [temp(1:tempr.T0burnin,:)./tau2_1;temp(tempr.T0burnin+1:end,:)./tau2_2];
            temp_2(1,:) = mean(temp(tempr.T0temp(:,1)<=set.T0period(end),:),1);
        else
            tau2 = tau2/yrs;
            G = 1-1/tau2;
            G_ = (ones(1,nyrs)*G).^[0:nyrs-1];
            G_M2 = toeplitz(G_,[G_(1) zeros(1,size(G_,2)-1)]);
            temp_2 = temp/tau2;
            temp_2(1,:) = mean(temp(tempr.T(:,1)<=set.T0period(end),:),1);
        end
    end
    
    % Calc T0 for each of the Tnum temperature realizations
    for i = 1:set.Tnum
       T01(:,i) = G_M1 * temp_1(:,i); 
       if strcmp(set.model,'TwoTau')
            T02(:,i) = G_M2 * temp_2(:,i); 
       end
    end
    if ~strcmp(set.model,'TwoTau')
        T02 = nan(size(T01));
    end

end
