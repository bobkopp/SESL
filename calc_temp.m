function temp = calc_temp(set,tempr)
    
% Represent the temperature uncertainty
% 
% temp = calc_temp(set,tempr)
% 
% INPUT:
% - set   -> general settings
% - tempr -> original temperature input
%
% OUTPUT:
% - temp  -> temperature draws with uncertainty depending on the temperature error model     

    if strcmp(set.T_err,'ar1')
        ar1 = set.tau_ar1;
    elseif strcmp(set.T_err,'ar1ts')
        ar1 = set.tau_ar1;
    else
        ar1 = [];
    end

    if strcmp(set.T_err,'no')
        if set.UseMarT0 %additionally to T_data, Marcott is used to calc. T0
            temp = tempr.T0temp(:,2)';
        else
            temp = tempr.T(:,2)';
        end
    elseif strcmp(set.T_err,'default')
        if set.UseMarT0 %additionally to T_data, Marcott is used to calc. T0
            T = tempr.T0temp;
        else
            T = tempr.T;
        end
        % add noise, sampled from a normal distribution with sig=Terr, to T every ryrs
        ryrs = set.ryrs/tempr.yrs; % Number of years for which to add same noise
        nyrs = length(T);
        Terr = T(:,3);
        r = randn(round(nyrs/ryrs + 1), 1);
        
        mat =  r(1:max(size(Terr(1:ryrs:end))),1) * ones(1,ryrs);
        tempE = reshape(mat',numel(mat),1);
        
        tempE = Terr.* tempE(1:nyrs);
        
        temp = (T(:,2) + tempE)';

    elseif strcmp(set.T_err(1:2),'ar')     
        if set.UseMarT0 %additionally to T_data, Marcott is used to calc. T0
            T = tempr.T0temp;
        else
            T = tempr.T;
        end
        if strcmp(set.T_err,'ar1');
            Cov_ar1 = bsxfun(@times,T(:,3),T(:,3)').* ((ar1*ones(length(T))).^abs(bsxfun(@minus,T(:,1),T(:,1)')));    
        elseif strcmp(set.T_err,'ar1ts'); 
            Cov_ar1 = bsxfun(@times,T(:,3),T(:,3)').* exp(-abs(bsxfun(@minus,T(:,1),T(:,1)')/ar1)); 
        end
        temp = mvnrnd(T(:,2)',Cov_ar1,set.Tnum);   
    end
    
    temp = temp';

end
