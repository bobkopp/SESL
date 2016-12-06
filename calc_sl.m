function pri2 = calc_sl(set,tempr,temp,T01,T02,param)
    
% Calculates a sea-level time-series, dependent on the model used
% 
% pri2 = calc_sl(set,tempr,temp,T01,T02,param)
% 
% INPUT:
% - set   -> general settings
% - tempr -> original temperature input
% - temp  -> temperature draw, dependent on the temperature error model
% - T01   -> equilibrium temperature 1
% - T02   -> equilibrium temperature 2, only if strcmp(model,'TwoTau')==1
% - param -> the parameter set to test [a1, a2, c, tau1, tau2, To(1)]
% 
% OUTPUT:
% - sl  -> sea level;
% - dsl -> diff(sl);
% - yr  -> year;
% - c   -> c;
% - T   -> temperature;
% - T01 -> equilibrium temperature 1
% - T02 -> equilibrium temperature 2, only if strcmp(model,'TwoTau')==1
     
     
    a1 = param(1);
    a2 = param(2);
    c = param(3);
    tau1 = param(4); 
    tau2 = param(5); 
    
    [year, temp, T01, T02] = resize_T(set,tempr,temp,T01,T02); % make a step-like yearly time series

    % Calc dSL/dt = c + a1*(temp - T01) + a2*(temp - T02);
    if strcmp(set.model,'TwoTau')
        dsea = a1*(temp - T01) + a2*(temp - T02); 
        pri2.T02 = T02;
    elseif strcmp(set.model,'ConstRate')
        dsea = c + a1*(temp - T01);
    elseif strcmp(set.model,'CRdecay')
        % calc. c decaying exponentially on timescale tau2
        g = 1 - 1/tau2;
        G = (ones(size(temp,1),1)*g).^([0:size(temp,1)-1]');
        c = c*G;
        c = repmat(c,1,set.Tnum);
        dsea = c + a1*(temp - T01);
    elseif strcmp(set.model,'CRovTau')
        dsea = c/tau1 + a1*(temp - T01);
    elseif strcmp(set.model,'simpel')
        dsea = a1*(temp - T01);
    end

    sea = cumsum(dsea); 

if real(sea)~=sea
   disp('sea not real')
end
     pri2.sl = sea;
     pri2.dsl = dsea;
     pri2.yr = year;
     pri2.c = c;
     pri2.T = temp;
     pri2.T01 = T01;
    
end % end calc_sl
