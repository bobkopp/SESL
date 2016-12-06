function O = Calc_SESLProjection(S,rcps,dig, NumTmag)

% Calculate sea-level projections under different RCP scenarios.
%
% O = Calc_SESLProjection(S,rcps,dig, NumTmag)
%
% INPUT:
% - S       -> Output of Calc_SESL_Prc.m
% - rcps    -> Cell array filled with strings. RCPs op use: {'RCP3PD', 'RCP45', 'RCP85'}
% - dig     -> number of digits to use, else []
% - NumTmag -> Number of MAGICC RCP realizations to take max&default == 600
%
% OUTPUT: Sea-level time-series for each selected RCP (Psl) and year dates (Pslyr). 
%         The columns in Psl give the percentiles Prc as defined below
    
    usesingle = true; % to save memory 
    
    timeprojected = [2000 2100];
        
    Prc = [5 17 50 83 95]; % Percentiles (68.27=1sigma) 

    folder = 'Data';
    
    if usesingle
        SFnc = @(x)single(x);
    else
        SFnc = @(x)x;
    end
    
    if isempty(NumTmag)
        NumTmag=600;
    end
    Tnum = S.settings.Tnum;
    sample = S.settings.sample;
    
    if ~isempty(dig)
        digits(dig) 
    end
    a = SFnc(S.MH.Params(:,1));
    a = reshape(repmat(a',Tnum,1),Tnum*sample,1);
    a = repmat(a,NumTmag,1);
    if strcmp(S.settings.model,'CRdecay')
        c_2000 = SFnc(S.MH.c_2000);
        c_2000 = repmat(c_2000,NumTmag,1);

        tau_c = SFnc(S.MH.Params(:,5));
        tau_c = reshape(repmat(tau_c',Tnum,1),Tnum*sample,1);
        tau_c = repmat(tau_c,NumTmag,1);
    elseif strcmp(S.settings.model,'simpel')
        c = SFnc(zeros(size(S.MH.Params(:,3))));
        c = reshape(repmat(c',Tnum,1),Tnum*sample,1);
        c = repmat(c,NumTmag,1);
    elseif strcmp(S.settings.model,'CRovTau') 
        c = SFnc(S.MH.Params(:,3)./S.MH.Params(:,4));
        c = reshape(repmat(c',Tnum,1),Tnum*sample,1);
        c = repmat(c,NumTmag,1);
    elseif strcmp(S.settings.model,'ConstRate') 
        c = SFnc(S.MH.Params(:,3));
        c = reshape(repmat(c',Tnum,1),Tnum*sample,1);
        c = repmat(c,NumTmag,1);
    end
    tau = SFnc(S.MH.Params(:,4));
    tau = reshape(repmat(tau',Tnum,1),Tnum*sample,1);
    tau = repmat(tau,NumTmag,1);
    
    T0_2000 = SFnc(S.MH.T01_2000);
    T0_2000 = repmat(T0_2000,NumTmag,1);
    
    if strcmp(S.settings.model,'TwoTau')
        a2 = SFnc(S.MH.Params(:,2));
        a2 = reshape(repmat(a2',Tnum,1),Tnum*sample,1);
        a2 = repmat(a2,NumTmag,1);
        T02_2000 = SFnc(S.MH.T02_2000);
        T02_2000 = repmat(T02_2000,NumTmag,1);
    end
    
    Toff1 = SFnc(S.MH.Tm1970_2000);
    Toff1 = repmat(Toff1,NumTmag,1);
    settings = S.settings;
    clear S
        
    for j = 1:length(rcps) % Count RCPs
        
        sl = SFnc(zeros(size(a)));
        
        if strcmp(settings.model,'CRdecay')
            c = c_2000;
        end
        T0 = T0_2000;
        if strcmp(settings.model,'TwoTau')
            T02 = T0_2000;
        end
        % load MAGICC realizations of RCPS: TIME & DATA 
        if usesingle
            load(fullfile(folder, [cell2mat(rcps(j)), '_single'])); % load rcp 
        else
            load(fullfile(folder, cell2mat(rcps(j)))); % load rcp 
        end

        DATA = DATA(:,1:NumTmag);
        
        Toff2 = repmat(mean(DATA(TIME>=1970 & TIME<=2000,:)),size(DATA,1),1);
        
        DATA = DATA-Toff2;
        
        DATA = DATA(TIME>=timeprojected(1) & TIME<=timeprojected(end),:);
        
        clear TIME T_90prc colheaders textdata Toff2
        
        fprintf('\t %1.0f years of %1s are now calculated: \n',size(DATA,1),cell2mat(rcps(j)));
        
        for i_yr = 1:size(DATA,1); % count years to be projected
            fprintf('%1.0f \t',i_yr);

            Tfc = DATA(i_yr,:);
            Tfc = reshape(repmat(Tfc,Tnum*sample,1),NumTmag*Tnum*sample,1);
            Tfc = Tfc + Toff1;
         
            if i_yr>1 % so it will be relative to 2100
                if ~strcmp(settings.model,'TwoTau')
                    sl = sl + a .* (Tfc-T0) + c;      
                else
                    sl = sl + a .* (Tfc-T0) + a2 .* (Tfc-T02);
                end
            end
    
            T0 = T0 + (1./tau).*(Tfc-T0);
            if strcmp(settings.model,'TwoTau')
                T02 = T02 + (1./tau).*(Tfc-T0);
            end
            if strcmp(settings.model,'CRdecay')
                c = c .* (ones(size(c))-1./tau_c);
            end
            clear Tfc
            
            O.(cell2mat(rcps(j))).Psl(i_yr,:) = prctile(sl,Prc); 
            O.(cell2mat(rcps(j))).Pslyr=timeprojected(1):timeprojected(end);
        end
        fprintf('\n')
    end
end



% sample  Tnum  600
% P1      T1     TR1
% P1      T2     TR1
% P1      T3     TR1
% 
% P2      T1     TR1
% P2      T2     TR1
% P2      T3     TR1
