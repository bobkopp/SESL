function Ho_offs = find_optim_Ho(set,P,S)
    
% Find optimal offset between simulated and proxy sea level in a least square sense.
%
% Ho_offs = find_optim_Ho(set,P,S)
% 
% INPUT:
% - set -> general settings
% - P   -> sea-level data input
% - S   -> output of calc_sl, i.e. semi-empirical sea-level simulation
% 
% OUTPUT: Ho_offs is the optimal (in a least square sense) offset of semi-empirical
%         and input sea-level data to be subtracted from the latter one.

    firstyear_calib = set.calibperiod(1);
    lastyear_calib = set.calibperiod(end);
    
    % reduce Proxy timeseries to calibration time frame
    ix = P.yr>=firstyear_calib & P.yr<=lastyear_calib; 
    Psea = P.sl(ix);
    Pyear = P.yr(ix);
    Pslerr = P.slerr(ix).^2;
    
    % bring SL_simu to size of proxy (same years)
    ix_ = ismember(S.yr,round(Pyear));
    sea = S.sl(ix_,:);
    year = S.yr(ix_);
    
    % double check size of proxy and andjust to SL_simu if necessary
    ix__ = ismember(round(Pyear),year);
    Psea = Psea(ix__);
    Pslerr = Pslerr(ix__);
    
    % Define error structure
    if exist('C','var') == 1 && set.useCov
        Pcov = P.C(ix,ix); %Pcov = diag(diag(C(ix,ix)));   
        Pcov = Pcov(ix__,ix__);
    else
        Pcov = diag(Pslerr);
    end    
    if ~(size(Pslerr,2)==size(Psea,2)) || ~(size(Pslerr,2)==size(sea,1)) || ~(size(Psea,2)==size(sea,1)) || ~(size(Pslerr,2)==size(Pcov,2)) || ~(size(sea,1)==size(Pcov,2)) || ~(size(Psea,2)==size(Pcov,2)) 
       error('check size') 
    end
    
    % Calc opt. offset
    Psea = repmat(Psea',1,size(sea,2));
    Ho_offs = lscov(ones(size(Pslerr')),Psea-sea,Pcov);

end
