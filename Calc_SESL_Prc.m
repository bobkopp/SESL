function P = Calc_SESL_Prc(T, num,folder)

% Calculates percentiles of sea-level and temperature (which were previously 
% stored by Calc_SL_from_Param.m) alongside with other values important for 
% projections or counterfactual scenarios
% 
% P = Calc_SESL_Prc(T, num,folder)
% 
% INPUT:
% - T      -> string, name of the temperature used for calibration 
% - num    -> number of files previously stored by Calc_SL_from_Param.m
% - folder -> string, name of the folder where Calc_SL_from_Param.m stored
%             away files
%
% OUTPUT:
% - settings -> settings that were also in previous outputs
% - data     -> original sea-level and temperature data used for calibration
% - MH       -> many items related to temperature T, equilibrium
%               temperature T0, the parameter c and sea level SL . Examples:
%               - T01_2000:      the values of T01 in the year 2000 CE
%               - psl:           percentiles of SL as prescribed by Prc below
%               - T01_1000BC:    the values of T01 in the year 1000 BC
%               - Tm1801_1900:   the mean value of T between 1801-1900 CE
%               - SLpo1800_2000: 1st degree polynomial (polyfit(x,y,1)) of SL 
%                                 between 1800-2000 CE

    calc_Full_T0 = true; % another function at the end of the script to combine all 'num' files for full T0 sample

    Prc = [5 17 50 83 95]; % Percentiles to calculate
    
    S = []; % EM added
    load(fullfile(folder, [T, '_1']));
    P.settings = S.settings;
    P.data = S.data;
    P.MH.Params = S.MH.Params;
    P.MH.alpha = S.MH.alpha;
    P.MH.yrT = S.MH.yrT;
    P.MH.yr = S.MH.yrsl;
    if strcmp(P.settings.model,'CRdecay')
        P.MH.yrc = S.MH.yrc;
    end
    % find proper divisor of the number of sl length
    DivSL = [];
    for x = 1:size(S.MH.sl,2)
        if mod(size(S.MH.sl,2),x)==0
            DivSL = [DivSL,x];
        end
    end
    Nyst = DivSL(2); % choose smallest so it won't run out of memory 
    yst = size(S.MH.sl,2)/Nyst;
    clear S
    
    % Define Output
    sl_ = [];
    psl = [];
    T01_2000 = [];
    T01_1900 = [];
    T01_501 = [];
    T01m501_700 = [];
    T01_1000BC = [];
    T01_2000BC = [];
    if strcmp(P.settings.model,'TwoTau')
        T02_2000 = [];
        T02_1900 = [];
        T02_501 = [];
        T02m501_700 = [];
        T02_1000BC = [];
        T02_2000BC = [];
    end
    if strcmp(P.settings.model,'CRdecay')
        c_2000 = [];
        c_1900 = [];
        c_501 = [];
    end
    sl_20C = [];
    Tpo500_1800 = [];
    Tm1801_1900 = [];
    Tm500_1800 = [];
    Tm501_700 = [];
    Tm2000_1800BC = [];
    Tm1850_1900 = [];
    Tm1860_1900 = [];
    Tm1970_2000 = [];
    SLpo800_1800 = [];
    SLpo1800_2000 = []; 
    SLpo0_1700 = [];
    SLpo0_400 = [];
    SLpo400_800 = [];
    SLpo800_1200 = [];
    SLpo1200_1600 = [];
    SLpo1600_1800 = [];
    SLpo1800_1900 = [];
    SLpo1900_2000 = [];
    
    fprintf('%1.0f',zeros(1,Nyst*num));fprintf('\n');
    for i = 1:Nyst
        for j = 1:num
            fprintf('0')
            load(fullfile(folder, [T, '_', num2str(j)]));
            yrs = S.data.temp.yrs;
            % sl ----------------------------------------------------------
            sl_ = [sl_; S.MH.sl(:,(i-1)*yst+1:i*yst)];
%             
            if i==1
              % T01(2000) -------------------------------------------------
                lyr =2000;
                if sum(S.MH.yrT==lyr)==0
                    T01_2000 = [T01_2000; S.MH.T01(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                else
                    T01_2000 = [T01_2000;S.MH.T01(:,S.MH.yrT == lyr)];
                end
              % T01(1900) -------------------------------------------------
                lyr =1900;
                if sum(S.MH.yrT==lyr)==0
                    T01_1900 = [T01_1900; S.MH.T01(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                else
                    T01_1900 = [T01_1900;S.MH.T01(:,S.MH.yrT == lyr)];
                end
              % T01(501) --------------------------------------------------
                lyr =501;
                if sum(S.MH.yrT==lyr)==0
                    T01_501 = [T01_501; S.MH.T01(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                else
                    T01_501 = [T01_501;S.MH.T01(:,S.MH.yrT == lyr)];
                end
              %T01(1000BC) ------------------------------------------------
                lyr =-1000;
                if sum(S.MH.yrT==lyr)==0
                    T01_1000BC = [T01_1000BC; S.MH.T01(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                else
                    T01_1000BC = [T01_1000BC;S.MH.T01(:,S.MH.yrT == lyr)];
                end
              %T01(2000BC) ------------------------------------------------
                lyr =-2000;
                if sum(S.MH.yrT==lyr)==0
                    T01_2000BC = [T01_2000BC; S.MH.T01(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                else
                    T01_2000BC = [T01_2000BC;S.MH.T01(:,S.MH.yrT == lyr)];
                end
                %T01(501-700) ---------------------------------------------
                ix = S.MH.yrT>=501 & S.MH.yrT<=700;
                T01m501_700 = [T01m501_700; mean(S.MH.T01(:,ix),2)];
                if strcmp(P.settings.model,'TwoTau')
                    % T02(2000) -------------------------------------------------
                    lyr =2000;
                    if sum(S.MH.yrT==lyr)==0
                        T02_2000 = [T02_2000; S.MH.T02(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                    else
                        T02_2000 = [T02_2000;S.MH.T02(:,S.MH.yrT == lyr)];
                    end
                  % T02(1900) -------------------------------------------------
                    lyr =1900;
                    if sum(S.MH.yrT==lyr)==0
                        T02_1900 = [T02_1900; S.MH.T02(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                    else
                        T02_1900 = [T02_1900;S.MH.T02(:,S.MH.yrT == lyr)];
                    end
                  % T02(501) --------------------------------------------------
                    lyr =501;
                    if sum(S.MH.yrT==lyr)==0
                        T02_501 = [T02_501; S.MH.T02(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                    else
                        T02_501 = [T02_501;S.MH.T02(:,S.MH.yrT == lyr)];
                    end
                  %T02(1000BC) ------------------------------------------------
                    lyr =-1000;
                    if sum(S.MH.yrT==lyr)==0
                        T02_1000BC = [T02_1000BC; S.MH.T02(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                    else
                        T02_1000BC = [T02_1000BC;S.MH.T02(:,S.MH.yrT == lyr)];
                    end
                  %T02(2000BC) ------------------------------------------------
                    lyr =-2000;
                    if sum(S.MH.yrT==lyr)==0
                        T02_2000BC = [T02_2000BC; S.MH.T02(:,S.MH.yrT>=lyr-ceil((yrs-1)/2) & S.MH.yrT<=lyr+floor((yrs-1)/2))];
                    else
                        T02_2000BC = [T02_2000BC;S.MH.T02(:,S.MH.yrT == lyr)];
                    end
                    %T02(501-700) ---------------------------------------------
                    ix = S.MH.yrT>=501 & S.MH.yrT<=700;
                    T02m501_700 = [T02m501_700; mean(S.MH.T02(:,ix),2)];
                end
                if strcmp(P.settings.model,'CRdecay')
                  % c(2000) -----------------------------------------------
                    lyr =2000;
                    if sum(S.MH.yrT==lyr)==0
                        c_2000 = [c_2000; S.MH.c(:,S.MH.yrc>=lyr-ceil((yrs-1)/2) & S.MH.yrc<=lyr+floor((yrs-1)/2))];
                    else
                        c_2000 = [c_2000;S.MH.c(:,S.MH.yrc == lyr)];
                    end
                  % c(1900) -----------------------------------------------
                    lyr =1900;
                    if sum(S.MH.yrT==lyr)==0
                        c_1900 = [c_1900; S.MH.c(:,S.MH.yrc>=lyr-ceil((yrs-1)/2) & S.MH.yrc<=lyr+floor((yrs-1)/2))];
                    else
                        c_1900 = [c_1900;S.MH.c(:,S.MH.yrc == lyr)];
                    end
                  % c(501) ------------------------------------------------
                    lyr =501;
                    if sum(S.MH.yrT==lyr)==0
                        c_501 = [c_501; S.MH.c(:,S.MH.yrc>=lyr-ceil((yrs-1)/2) & S.MH.yrc<=lyr+floor((yrs-1)/2))];
                    else
                        c_501 = [c_501;S.MH.c(:,S.MH.yrc == lyr)];
                    end
                end
              % sl(1900-2000) ---------------------------------------------
                sl_20C_ = S.MH.sl(:,S.MH.yrsl>=1900 & S.MH.yrsl<=2000);
                sl_20C =  [sl_20C;sl_20C_ - repmat(sl_20C_(:,1),1,size(sl_20C_,2))];
              % T linear rate 500-1800 ------------------------------------
                ix = S.MH.yrT>=500 & S.MH.yrT<=1800;
                for k = 1: size(S.MH.temp,1)
                    Tpo500_1800_(k,:) = polyfit(S.MH.yrT(ix),S.MH.temp(k,ix)',1);            
                end
                Tpo500_1800 = [Tpo500_1800;Tpo500_1800_]; clear Tpo500_1800_;
              % mean T(1801-1900) -----------------------------------------
                ix = S.MH.yrT>=1801 & S.MH.yrT<=1900;
                Tm1801_1900 = [Tm1801_1900; mean(S.MH.temp(:,ix),2)];
              % mean T(500-1800) ------------------------------------------
                ix = S.MH.yrT>=501 & S.MH.yrT<=1800;
                Tm500_1800 = [Tm500_1800; mean(S.MH.temp(:,ix),2)];
              % mean T(500-700) -------------------------------------------
                ix = S.MH.yrT>=501 & S.MH.yrT<=700;
                Tm501_700 = [Tm501_700; mean(S.MH.temp(:,ix),2)];
              % mean T(-2000--1800) ---------------------------------------
                ix = S.MH.yrT>=-2000 & S.MH.yrT<=-1800;
                Tm2000_1800BC = [Tm2000_1800BC; mean(S.MH.temp(:,ix),2)];
              % mean T(1850-1900) -----------------------------------------
                ix = S.MH.yrT>=1850 & S.MH.yrT<=1900;
                Tm1850_1900 = [Tm1850_1900; mean(S.MH.temp(:,ix),2)];
                %mean T(1860-1900) ----------------------------------------
                ix = S.MH.yrT>=1860 & S.MH.yrT<=1900;
                Tm1860_1900 = [Tm1860_1900; mean(S.MH.temp(:,ix),2)];
              % mean T(1970-2000) -----------------------------------------
                ix = S.MH.yrT>=1970 & S.MH.yrT<=2000;
                Tm1970_2000 = [Tm1970_2000; mean(S.MH.temp(:,ix),2)];       
              % SL linear rate 0-1700 -------------------------------------
                ix = S.MH.yrsl>=0 & S.MH.yrsl<=1700;
                for k = 1: size(S.MH.temp,1)
                    SLpo0_1700_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo0_1700 = [SLpo0_1700;SLpo0_1700_]; clear SLpo0_1700_;
              % SL linear rate 800-1800 -----------------------------------
                ix = S.MH.yrsl>=800 & S.MH.yrsl<=1800;
                for k = 1: size(S.MH.temp,1)
                    SLpo800_1800_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo800_1800 = [SLpo800_1800;SLpo800_1800_]; clear SLpo800_1800_;
              % SL linear rate 1800-2000 ----------------------------------
                ix = S.MH.yrsl>=1800 & S.MH.yrsl<=2000;
                for k = 1: size(S.MH.temp,1)
                    SLpo1800_2000_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo1800_2000 = [SLpo1800_2000;SLpo1800_2000_]; clear SLpo1800_2000_;
              % SL linear rate 0-400 --------------------------------------
                ix = S.MH.yrsl>=0 & S.MH.yrsl<=400;
                for k = 1: size(S.MH.temp,1)
                    SLpo0_400_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo0_400 = [SLpo0_400;SLpo0_400_]; clear SLpo0_400_;
              % SL linear rate 400-800 ------------------------------------
                ix = S.MH.yrsl>=400 & S.MH.yrsl<=800;
                for k = 1: size(S.MH.temp,1)
                    SLpo400_800_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo400_800 = [SLpo400_800;SLpo400_800_]; clear SLpo400_800_;
              % SL linear rate 800-1200 -----------------------------------
                ix = S.MH.yrsl>=800 & S.MH.yrsl<=1200;
                for k = 1: size(S.MH.temp,1)
                    SLpo800_1200_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo800_1200 = [SLpo800_1200;SLpo800_1200_]; clear SLpo800_1200_;
              % SL linear rate 1200-1600 ----------------------------------
                ix = S.MH.yrsl>=1200 & S.MH.yrsl<=1600;
                for k = 1: size(S.MH.temp,1)
                    SLpo1200_1600_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo1200_1600 = [SLpo1200_1600;SLpo1200_1600_]; clear SLpo1200_1600_;
              % SL linear rate 1600-1800 ----------------------------------
                ix = S.MH.yrsl>=1600 & S.MH.yrsl<=1800;
                for k = 1: size(S.MH.temp,1)
                    SLpo1600_1800_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo1600_1800 = [SLpo1600_1800;SLpo1600_1800_]; clear SLpo1600_1800_;
              % SL linear rate 1800-1900 ----------------------------------
                ix = S.MH.yrsl>=1800 & S.MH.yrsl<=1900;
                for k = 1: size(S.MH.temp,1)
                    SLpo1800_1900_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo1800_1900 = [SLpo1800_1900;SLpo1800_1900_]; clear SLpo1800_1900_;
              % SL linear rate 1900-2000 ----------------------------------
                ix = S.MH.yrsl>=1900 & S.MH.yrsl<=2000;
                for k = 1: size(S.MH.temp,1)
                    SLpo1900_2000_(k,:) = polyfit(S.MH.yrsl(ix),S.MH.sl(k,ix),1);            
                end
                SLpo1900_2000 = [SLpo1900_2000;SLpo1900_2000_]; clear SLpo1900_2000_;
            end
            clear S
        end
        psl_ = prctile(sl_,Prc);
        if size(psl_,1) == length(Prc)
            psl = [psl, prctile(sl_,Prc)];
        elseif size(psl_,2) == length(Prc);
            psl = [psl, prctile(sl_,Prc)'];
        end
        sl(:,(i-1)*yst+1:i*yst) = sl_;
        sl_ = [];
    end
    
    P.MH.psl = psl;
    P.MH.sl = sl;
    P.MH.sl_20C = sl_20C;
    P.MH.T01_2000 = T01_2000;
    P.MH.T01_1900 =T01_1900;
    if size(T01_501,2)>1
        P.MH.T01_501 =T01_501(:,2);
    else
        P.MH.T01_501 =T01_501;
    end
    P.MH.T01_1000BC = T01_1000BC;
    P.MH.T01_2000BC = T01_2000BC;
    P.MH.T01m501_700 = T01m501_700;
    if strcmp(P.settings.model,'TwoTau')
        P.MH.T02_2000 = T02_2000;
        P.MH.T02_1900 =T02_1900;
        if size(T02_501,2)>1
            P.MH.T02_501 =T02_501(:,2);
        else
            P.MH.T02_501 =T02_501;
        end
        P.MH.T02_1000BC = T02_1000BC;
        P.MH.T02_2000BC = T02_2000BC;
        P.MH.T02m501_700 = T02m501_700;
    end
    if strcmp(P.settings.model,'CRdecay')
        P.MH.c_2000 = c_2000;
        P.MH.c_1900 = c_1900;
        if size(c_501,2)>1
            P.MH.c_501 = c_501(:,2);
        else
            P.MH.c_501 = c_501;
        end
    end
    P.MH.Tpo500_1800 = Tpo500_1800;
    P.MH.Tm1801_1900 = Tm1801_1900;
    P.MH.Tm500_1800 = Tm500_1800;
    P.MH.Tm501_700 = Tm501_700;
    P.MH.Tm2000_1800BC = Tm2000_1800BC;
    P.MH.Tm1850_1900 = Tm1850_1900;
    P.MH.Tm1860_1900 = Tm1860_1900;
    P.MH.Tm1970_2000 = Tm1970_2000;
    P.MH.SLpo0_1700 = SLpo0_1700;
    P.MH.SLpo800_1800 = SLpo800_1800;
    P.MH.SLpo1800_2000 = SLpo1800_2000;
    P.MH.SLpo0_400 = SLpo0_400 ;
    P.MH.SLpo400_800 = SLpo400_800 ;
    P.MH.SLpo800_1200 = SLpo800_1200 ;
    P.MH.SLpo1200_1600 = SLpo1200_1600 ;
    P.MH.SLpo1600_1800 = SLpo1600_1800 ;
    P.MH.SLpo1800_1900 = SLpo1800_1900 ;
    P.MH.SLpo1900_2000 = SLpo1900_2000 ;
    P.MH.sl = sl;

if calc_Full_T0
    for i = 1:num
        load(fullfile(folder, [T, '_', num2str(i)]))
        if i==1
            Num = size(S.MH.sl,1);
            P.MH.T0 = ones(Num*num, size(S.MH.T01,2));
            if strcmp(P.settings.model,'TwoTau')
                P.MH.T02 = ones(Num*num, size(S.MH.T02,2));
            end
        end
        S.MH.Params = [];
        S.MH.alpha= [];
        S.MH.temp= [];
        S.MH.c= [];
        S.MH.sl= [];
        S.MH.yrsl= [];
        S.MH.yrc= [];
        S.MH.yrT= [];
        P.MH.T0((i-1)*Num+1:i*Num,:) = S.MH.T01;
        if strcmp(P.settings.model,'TwoTau')
            P.MH.T02((i-1)*Num+1:i*Num,:) = S.MH.T02;
        end
    end
end
fprintf('\n')
end