function Calc_SL_from_Param(S,storeNum,Meth,folder)

% Calculates time series of sea level, temperatures and c from the
% parameter set found in SESL.m
%
% Calc_SL_from_Param(S,storeNum,Meth,folder)
%
% INPUT:
% - S        -> Output structure from SESL.m
% - storeNum -> storeNum time series of sea level, temperature, equilibrium 
%               temperature and c are calculated, numbered & saved in 'folder'. 
%               If there is no memory issue storeNum can be set to sample*Tnum  
% - Meth     -> The method used befor by SESL.m
%               'MH' for Metropolis Hastings, 
%               'Slice' for Slice sampling (not working correctly)
% - folder   -> folder to save samples with sizes of 'storeNum' as prescribed below
% 
% OUTPUT: a number of output structures S are saved and numbered in the folder 
%         defined above. Each of these outputs is similar to the SESL.m output 
%         only that the structure S.MH or S.Slice (dependent on the input 
%         Method 'Meth') will be more detailed.
%         Besides Params & alpha, S.MH/Slice includes time series of T01
%         (and T02 if needed), c, temp (temperature), sl (sea level) and
%         the according year counts yrT, yrc, yrsl. The number of samples these
%         time series include is defined by storeNum and is usually smaller
%         than the full set of samples, that is why these outputs later
%         need to be combined by Calc_SESL_Prc.m. This is a memory issue

Tnum = []; % calc. Tnum random temperatures and sea-levels for each param. set and take all. 
           % if isempty(Tnum) take as many random temps. as prescribend in settings

dat = S.data;
set = S.settings;
if strcmp(Meth,'SimAnn')
    smpl = S.(Meth).OptParam;
else
    smpl = S.(Meth).Params;
end

pb = 100; % show a progress bar, progressing every 100 samples
if size(smpl,1)<pb; pb=size(smpl,1);end;
fprintf('%1.0f',zeros(1,size(smpl,1)/pb));fprintf('\n')

ct = 0;
ct2 = 0;
for i = 1:size(smpl,1)
    ct = ct+1;
    
% ------------- Calc temp
    
    if ~isempty(Tnum)
        set.Tnum=Tnum;
    end
    temp = calc_temp(set,dat.temp);
 
% ------------- Calc T0 & SL

    [T01, T02] = calc_T0(set,temp,dat.temp,smpl(i,:));
    sl = calc_sl(set, dat.temp,temp,T01,T02, smpl(i,:));
    Ho = find_optim_Ho(set,dat.sea.proxy,sl);
    Ho = repmat(Ho,size(sl.sl,1),1);
    sl_ = sl.sl+Ho;
    
% ------------- Reduce time series to their true information
        
    T01_((ct-1)*set.Tnum+1:ct*set.Tnum,:) = T01';
    if strcmp(set.model,'TwoTau')
        T02_((ct-1)*set.Tnum+1:ct*set.Tnum,:) = T02';
    end
    temp_((ct-1)*set.Tnum+1:ct*set.Tnum,:) = temp';
    year = sl.yr;
    if strcmp(set.model,'CRdecay')
        if S.settings.UseMarT0
            ixc = ismember(year,S.data.temp.T0temp(:,1));
        else
            ixc = ismember(year,S.data.temp.T(:,1));
        end
        c((ct-1)*set.Tnum+1:ct*set.Tnum,:) = sl.c(ixc,:)';
        yrc = year(ixc);
    else
        c = sl.c;
    end

    if S.settings.UseMarT0
        dyear_ = diff(S.data.temp.T0temp(:,1));
        dyear = [dyear_; dyear_(end)];
        ixsl = ismember(year,S.data.temp.T0temp(:,1)+floor(dyear./2));
    else
        dyear_ = diff(S.data.temp.T(:,1));
        dyear = [dyear_; dyear_(end)];
        ixsl = ismember(year,S.data.temp.T(:,1)+floor(dyear./2));
    end 
    SL((ct-1)*set.Tnum+1:ct*set.Tnum,:) = sl_(ixsl,:)';
    yrsl = year(ixsl);

    % store parts of x (now x==inf) samples
    if (ct*set.Tnum)==storeNum || i==size(smpl,1) 
        ct2 = ct2+1;
        S.(Meth).T01 = T01_;
        if strcmp(set.model,'TwoTau')
            S.(Meth).T02 = T02_;
        end
        S.(Meth).c = c;
        if strcmp(set.model,'CRdecay')
            S.(Meth).yrc = yrc;
        end
        S.(Meth).temp = temp_;
        S.(Meth).sl = SL;
        S.(Meth).yrsl = yrsl;

        if S.settings.UseMarT0
            S.(Meth).yrT = S.data.temp.T0temp(:,1);
        else
            S.(Meth).yrT = S.data.temp.T(:,1);
        end
        save(fullfile(folder, [S.settings.T_data, '_', num2str(ct2)]), 'S'); 
        clear temp_ SL dSL ix year T01_;
        ct = 0;
    end
    
    if mod(i,pb)==0
        fprintf('%1.0f',0);
    end
end

fprintf('\n')

end