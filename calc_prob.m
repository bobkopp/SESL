function prob = calc_prob(S,param,T0500_700,T500_700,c2000,set,dat)
    
% Calculate the Gaussian probability of a certain parameter set.
% 
% prob = calc_prob(S,param,T0500_700,T500_700,c2000,set,dat)
% 
% INPUT: 
% - S         -> output of calc_sl.m
% - param     -> the parameter set to test [a1, a2, c, tau1, tau2, To(1)]
% - T0500_700 -> mean equilibrium temperature between 500-700 CE
% - T500_700  -> mean temperature between 500-700 CE
% - c2000     -> c in the year 2000 CE
% - set       -> general settings
% - dat       -> input data, temperature and sea-level
% 
% OUTPUT: prob is a gaussian probability as defined in Kopp et al. 2016 PNAS

    % Proxy SL data
    P = dat.sea.proxy;

    % calculate optimal Ho between first & last calibration year
    firstyear_calib = set.calibperiod(1);
    lastyear_calib = set.calibperiod(end);
    if set.optHo    
        Ho_offs = find_optim_Ho(set,P,S);
    else
        offset = mean(S.sl(S.yr>=set.baseperiod(1) & S.yr<=set.baseperiod(end),:)); % offset to baseperiod
        offset = repmat(offset,size(S.sl,1),1);
        S.sl = S.sl - offset; 
    end

    % calc residual 'sq_' of simulated and proxy sea level
    S.sl = S.sl(S.yr>=firstyear_calib & S.yr<=lastyear_calib,:);
    S.yr = S.yr(S.yr>=firstyear_calib & S.yr<=lastyear_calib);

    ix1 = ismember(round(P.yr),S.yr);
    P.yr = round(P.yr(ix1));
    P.sl = P.sl(ix1);
    P.sl = repmat(P.sl',1,size(S.sl,2));

    ix2 = ismember(S.yr,P.yr);
    S.sl = S.sl(ix2,:);
    if set.optHo    
        Ho_offs = repmat(Ho_offs,size(S.sl,1),1);
        sq_ = (S.sl - P.sl + Ho_offs);
    else
        sq_ = (S.sl - P.sl);
    end

    % truncate SL covariance matrix to the size of residuals
    if  set.useCov
        C = P.C(ix1,ix1);
    else
        P.slerr = P.slerr(ix1);
        C = diag((P.slerr.^2).*set.fac1);
    end

    % calc likelihood of priors
    if ~set.all_flat_priors

        %param = [a1, a2, c, tau1, tau2, To(1)]
        lik(1) = log(set.lik_a1(param(1)));  
        lik(2) = log(set.lik_a2(param(2)));  
        lik(3) = log(set.lik_c(param(3)));
        if param(4)>0
            if set.TauLogUniform
                lik(4) = log(set.lik_tau1(log(param(4)))); % loguniform
            else
                lik(4) = log(set.lik_tau1(param(4))); % uniform
            end
        else
            lik(4) = log(0);
        end
        if param(5)>0
            if set.TauLogUniform
                lik(5) = log(set.lik_tau2(log(param(5)))); % loguniform
            else
                lik(5) = log(set.lik_tau2(param(5))); %uniform
            end
        else
            lik(5) = log(0);
        end
        if set.OptimT0
            lik(6) = log(set.lik_T01st(param(6)));  
        end

    else
        lik = 0;
    end

    w = nan(1,size(sq_,2));
    for i = 1:size(sq_,2) 
        if ~set.all_flat_priors
            lik(7) = log(set.lik_T02nd(T0500_700(i),T500_700(i))); 
            if strcmp(set.model,'CRdecay')
                lik(8) = log(set.lik_c2000(c2000(i)));
            end
        end
        sq = sq_(:,i);
        w(i) = (-.5 * sq')*(C\sq);
        if set.NormProb
            w(i) = w(i) - .5*log(2*pi*det(C));
        end
        w(i) = w(i) + sum(lik);
        w(i) = exp(w(i));
    end      

    prob = w;

end
