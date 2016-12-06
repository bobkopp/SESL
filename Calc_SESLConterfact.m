function S = Calc_SESLConterfact(S,CFtemp,column)

% Calculate sea level with different temperature inputs, factual and counterfactual.
%
% S = Calc_SESLConterfact(S,CFtemp,column)
%
% INPUT:    
% - S      -> Output Structure of calc_SESL_Prc.m
% - CFtemp -> string, input temperature to use: 
%             - 'HadCRUT' - HadCRUT4 temperature,
%             - 'CMIP5_mean' - 'historicalNat_CMIP35_tas' mean from Gareth Jones <gareth.s.jones@metoffice.gov.uk>
%             - 'CMIP5' - 'historicalNat_CMIP35_tas' single simu (data column given as input) from Gareth Jones <gareth.s.jones@metoffice.gov.uk>
%             - 'mean' - mean temperature between 500-1800 CE
%             - 'linrate' - linear rate temperature between 500-1800 CE
% - column -> Data column of 'historicalNat_CMIP35_tas' if CFtemp=='CMIP5'
% - Below set the first ['fyr'] and last ['lyr'] year of simulation and the
%   whether to output also the input structure or only the calculation below
%
% OUTPUT: For each in the input selected CFtemp, the calulated sea level
%         (sl), temperature (Tcf), equilibrium temperature (T01) and year 
%         values (time) are given alongside with the general settings. If
%         below output_all == true some of the output of previous scripts 
%         such as data is repeated.

    fyr = 1900; % first forecast year
    lyr = 2015; % last forecastyear
       
    output_all = false; %Output all (true) or only the forecast structure

    sl = []; Tcf = []; T01 = [];
    if strcmp(CFtemp,'HadCRUT')
        Tfc = load(fullfile('Data','T_HC4_gl_2015.mat'));
    elseif strcmp(CFtemp, 'CMIP5_mean')
        Tfc = load(fullfile('Data','historicalNat_CMIP35_tas.mat'));
    elseif strcmp(CFtemp, 'CMIP5_rand')
        Tfc = load(fullfile('Data','historicalNat_CMIP35_tas.mat'));
        Tfc.Tfc = Tfc.T(:,column);
    end
    fprintf('\t "%s" temperature scenario:\n',CFtemp)
    pb = 1000; % show a progress bar, progressing every 1000 sample*Tnum
    if S.settings.sample*S.settings.Tnum<pb; pb=S.settings.sample*S.settings.Tnum;end;
    fprintf('%1.0f',zeros(1,S.settings.sample*S.settings.Tnum/pb));
    fprintf('\n')
    for i_P = 1:S.settings.sample
        for i_T = 1:S.settings.Tnum
            i = (i_P-1)*S.settings.Tnum + i_T;
            if mod(i,pb) == 0
                fprintf('%1.0f',0)
            end
            if strcmp(CFtemp,'mean')
                Tfc_ = S.MH.Tm500_1800(i);
                yrfc_ = fyr:lyr;
                Tfc_ = Tfc_*ones(size(yrfc_));
            elseif strcmp(CFtemp,'linrate')
                po = S.MH.Tpo500_1800(i,:);
                yrfc_ = fyr:lyr;
                Tfc_ = polyval(po,yrfc_);
                Tfc_ = Tfc_ - Tfc_(1) + S.MH.Tm1801_1900(i);
            elseif strcmp(CFtemp,'HadCRUT')
                Tfc_ = Tfc.T(:,2)';
                yrfc_ = Tfc.T(:,1)';
                Tfc_ = Tfc_ - mean(Tfc_(yrfc_>=1850 & yrfc_<=1900)) + S.MH.Tm1850_1900(i);
                Tfc_ = Tfc_(yrfc_>=fyr);
                yrfc_ = yrfc_(yrfc_>=fyr);
            elseif strcmp(CFtemp, 'CMIP5_mean')
                Tfc_ = Tfc.T(:,3)';
                yrfc_ = Tfc.T(:,1);
                Tfc_ = Tfc_ - mean(Tfc_(yrfc_>=1860 & yrfc_<=1900)) + S.MH.Tm1860_1900(i);
                Tfc_ = Tfc_(yrfc_>=fyr);
                yrfc_ = yrfc_(yrfc_>=fyr);
            elseif strcmp(CFtemp, 'CMIP5')         
                Tfc_ = Tfc.Tfc';
                yrfc_ = Tfc.T(:,1);
                Tfc_ = Tfc_ - mean(Tfc_(yrfc_>=1860 & yrfc_<=1900)) + S.MH.Tm1860_1900(i);
                Tfc_ = Tfc_(yrfc_>=fyr);
                yrfc_ = yrfc_(yrfc_>=fyr);
            end

            FC = calc_fc(S,Tfc_,i_P,i);
            sl = [sl;FC.sl];
            Tcf = [Tcf;Tfc_];
            T01 = [T01;FC.T01];

        end
    end
    
    fprintf('\n')
    
    S.proj.sl = sl;
    S.proj.Tcf = Tcf;
    S.proj.T01 = T01;
    S.proj.time = yrfc_;
    
    if lyr<=2000
        % Calculate likelihoods & Differences
        L = calc_lik(S);
        S.proj.Lik = L.Lik;
        S.proj.Lik_yrly = L.Lik_yrly;
        S.proj.Diff = L.Diff;
        S.proj.Diff_yrly = L.Diff_yrly;
        S.proj.Perc = L.Perc;
        S.proj.Perc_yrly = L.Perc_yrly;
        S.proj.slFactual = L.sl;
    end    
    if strcmp(CFtemp, 'CMIP5')
        S.proj.column_CMIP5 = column;
    end
    if ~output_all
        S_ = S.proj;
        S_.settings = S.settings;
        S = S_;
    end
    
end

function fc = calc_fc(S,Tfc,i_P,i)
    
    set = S.settings;
    
    T01 = S.MH.T01_1900(i);  
    if strcmp(S.settings.model,'TwoTau')
        T02 = S.MH.T02_1900(i);
    end
    
    if strcmp(S.settings.model,'CRdecay')
        c = S.MH.c_1900(i);
    else
        c = S.MH.Params(i_P,3);
    end

    a1 = S.MH.Params(i_P,1);
    a2 = S.MH.Params(i_P,2);
    tau1 = S.MH.Params(i_P,4);
    tau2 = S.MH.Params(i_P,5);
        
    nyrs = max(size(Tfc));

    for ii = 1:nyrs-1
        % Compute here T0(sample,year) as the tau-year memory of temp(j)
        T01(ii+1) = T01(ii) + (1/tau1)*(Tfc(ii) - T01(ii));
        if strcmp(set.model,'TwoTau')
            T02(ii+1) = T02(ii) + (1/tau2)*(Tfc(ii) - T02(ii)); 
        elseif strcmp(set.model,'CRdecay')
            c(ii+1) = c(ii) * (1-1/tau2);
        end
    end
    if strcmp(set.model,'TwoTau')
        dsl = a1*(Tfc - T01) + a2*(Tfc - T02); 
        fc.T02 = T02;
    elseif strcmp(set.model,'ConstRate')
        dsl = c + a1*(Tfc - T01);
    elseif strcmp(set.model,'CRdecay')
        dsl = c + a1*(Tfc - T01);
    elseif strcmp(set.model,'CRovTau')
        dsl = c/tau1 + a1*(Tfc - T01);
    elseif strcmp(set.model,'simpel')
        dsl = a1*(Tfc - T01);
    end
    sl = cumsum(dsl);       
    sl = sl-sl(1);

    fc.sl = sl;
    fc.T01 = T01;
end

function L = calc_lik(S)
    
    % which SL to take as factual
    fact_SL = 'calib'; % 'calib': semi-emp calibration 
                       % 'prox': proxy data
                       % 'recalc': semi-empiricaly recalculate with factual temperature
    %---- SL projection under counterfactual T
    sl = S.proj.sl;
    yr = S.proj.time;
               
    %---- SL to calculate likelihood against: either Proxy or semi emp. calibration
    if strcmp(fact_SL,'calib')
        yrPrx = 1900:2000;
        slPrx = S.MH.sl_2000;
    elseif strcmp(fact_SL,'recalc')
        Tfc = S.data.temp.obs';
        yrPrx = S.data.temp.yr;
        ix = yrPrx>=yr(1) & yrPrx<=yr(end);
        fc = calc_fc(S,Tfc(ix),yrPrx(ix),[50,16,84]);
        slPrx = fc.sl;
        yrPrx = yrPrx(ix);
    end       
    % cut proxy to same period as T scenario
    ix = yrPrx>=yr(1) & yrPrx<=yr(end);
    yrPrx = yrPrx(ix);
    slPrx = slPrx(:,ix);
    for i = 1:size(sl,1)%count sample
        Diff_(i,:) = -sl(i,:)+slPrx(i,:);
        Perc_(i,:) = (sl(i,:)./slPrx(i,:))*100;
        for j = 1:size(sl,2)%count years
            L_(i,j) = sum(sl(:,j)>=slPrx(i,j))/size(sl,1);
        end
    end

    % index to give data the spacing of SL & T
    ix = ismember(yrPrx,S.MH.yr_);
    ix(yrPrx==yr(1))=true;
    ix(yrPrx==yr(end))=true;

    Lik(:,1) = yrPrx; % the year at which the likelihoods are calculated
    Lik(:,2) = mean(L_,1);
    L.Lik_yrly = Lik;
    L.Lik = Lik(ix,:);
    
    D = prctile(Diff_,[50 17 83 5 95]);
    L.Diff_yrly(:,1) = yrPrx;
    L.Diff_yrly(:,2:6) = D';
    L.Diff = L.Diff_yrly(ix,:);

    P = prctile(Perc_,[50 17 83 5 95]);
    P(:,1) = 100;
    L.Perc_yrly(:,1) = yrPrx;
    L.Perc_yrly(:,2:6) = P';
    L.Perc = L.Perc_yrly(ix,:);

    L.sl(:,1) = yrPrx(ix);
    L.sl(:,2:6) = prctile(slPrx(:,ix),[50 17 83 5 95])';    
        
end