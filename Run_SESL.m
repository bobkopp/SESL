function Run_SESL()

% Start SESL (semi-empirical sea-level) code after defining all settings.
% 
% Run_SESL()
%
% All settings can be defined below. Then the following scripts will be run
% and their outputs saved to 'Output':
% -> SESL.m                (Main SESL script to determine semi-empirical parameter distributions)
% -> Calc_SL_from_Param.m  (Calculates time series of sea level, temperatures and c from the parameter set found in SESL.m)
% -> Calc_SESL_Prc.m       (Calculate percentiles of sea-level, temperature & other important values)
% -> Calc_SESLProjection.m (Calculate sea-level projections under different RCP scenarios)
% -> Calc_SESLConterfact.m (Calculate sea level with different temperature inputs, factual and counterfactual)


start = NaN; % if NaN, starting parameter distribution will be set by simulated annealing,
             % else a parameterset should be given in the following order [a1, a2/b, c, log(tau1), log(tau2), T0(1)]
    
% MAIN SETTINGS
        % Model formulation: 
        % simpel    => dSL/dt = a * (T(t)-T0(t))
        % TwoTau    => dSL/dt = a1 * (T(t)-T0_1(t)) + a2 * (T(t)-T0_2(t))
        % ConstRate => dSL/dt = c + a * (T(t)-T0(t))
        % CRdecay   => dSL/dt = c(t) + a * (T(t)-T0(t))
        % CRovTau    => dSL/dt = c/tau + a * (T(t)-T0(t))
        %   -> dT0(t)/dt = (T-T0)/tau
        %   -> dc(t)/dt = -c(t)/tau
    model = 'CRdecay'; % 'ConstRate' 'CRdecay' 'CRovTau' 'TwoTau' 'simpel'
    % Sea-level data to load has the name 'SL_dat'_'SL_dat2'
    SL_dat = 'GLMW-1ts'; % right now only 'GLMW-1ts'. 'GLMW','GLMW-1amp1ts','GLMW-1ts','GLMW-Gr','GLMW-NC'
    T_data = 'Mann09_11'; % 'Mann08_eiv','Mann08_cps','Mann09','Mann09_11','Marcott13_RegEM-HC3_20','PAGES2k_13'
    
% MH SETTINGS
    JumpDist = 0.005; % Jumping distribution size, default = 0.005
    burning = 1000; % number of monte carlo samples used in the burning-in or spin-up period, default = 1000
    NumSkip = 500; % Thinning: only every NumSkip sample will be selected, default = 500
    sample = 1000; % number of monte carlo samples
    NormProb = false; % Use normalized probabilities -> should make no difference 
    
% TEMPERATURE SETTINGS
    Tnum = 100; % Number of ar1 temperatures to evaluate for likelihood, default = 100
    T_err = 'ar1ts'; % How should the temperature variance be treated ?
        % 'default' -> T + random noise as in Kemp et al. 2011 PNAS
        % 'no'      -> Don't add uncertainty
        % 'ar1'     -> T as AR(1) process with sig as 'default'
        % 'ar1ts'   -> AR(1) Parameter timescale == exp(-abs(t2-t1)/timescale)
    tau_ar1 = 10; % if Terr==ar1ts: timescale of ar1ts uncertainty above (default==10)
                  % if Terr==ar1: AR(1) parameter
    ryrs = 10; % If T_err==default: Number of years for which to add same size noise on temp. in calc_sl
    TerrSc = 1; % scaling of temperature uncertainty 
             
% SEA-LEVEL SETTINGS
    UseCov = true; % if existing, use covariance matrix of SL data to calc likelihood
    CovShrink = 0; % Sea level covariance matrix (C) shrikage parameter f. C gets shrinked to C' while C0 = diag(diag(C)): C' = f*C0 + (1-f)*C
    CovTau = 100; % if ~isnan(CovTau): take the elementwise product of the covariance and a tapering function exp(-delta(t)/CovTau) where CovTau is a time scale (default==100)
    NoNegCov = true; % remove negative Cov Mat entries and replace by 0
    fac1 = 1; % multiplier to SL uncertainty to adjust for autocorrelation ('Kemp et al. PNAS 2011' use 10)
                 
% MODEL SETTINGS
    optHo = true; % for each MC sample, calculate the optimal offset between simulated and calibration sea level via MLS
    OptimT0 = true; % optimize T0(1) (true) or set to mean(T(T0period))
    UseMarT0 = true; % use Marcott temperature to calc. T0 until 'T_data' starts
    T0temp_level = 100; % number of years over which to level Marcott and T_data if UseMarT0==true
    calibperiod = 500:2010; % CE years to include for model calibration. 
    period = -1000:2010; % CE years to include for SL calculation.
    baseperiod = 1400:1800; % % Period for which simulated sl and data are normed. Obsolete if OptimT0==true.
    T0period = -2000:-1800; % period to initialise T0(1)
    T0length = 0; % if==0: T0period is used, else T0(1) = mean(T(first 'T0length' years)) + prior
      
% PRIOR SETTINGS    
    StartDistr =  [0.5, 0.3, 0, log10(50), log10(100) , 0.5]; % starting distribution for simulated annealing: [a1, a2/b, c, log(tau1), log(tau2), T0(1)]
    all_flat_priors = false; % true, false
    AdjustPrior = false; % % adjust the prior in check_settings, according to the model used:true, false
    a1_prior = {'uniform' 0 2};
    a2_prior = {'uniform' 0 2};
    tau1_prior = {'uniform',log(30),log(3000)};
    tau2_prior = {'uniform',log(30),log(20000)};
    c_prior = {'flat', [],[]};
    c2000_prior = {'uniform', -.2,.2};
    T01st_prior = {'uniform',-.6,.6};
    T02nd_prior = {'normal',0,.2};
    TauLogUniform = true; % true => Tau priors are loguniform, else uniform

% SETTINGS FOR CALCULATING TIME SERIES FROM POSTERIOR PARAMETERS
    storeNum = 5000; % storeNum sl, T0, c & temperature curves are calculated from 
                     % the posterior parameter set & saved to the folder 'Simu\' 
                     % default = 5000, maximum = sample*Tnum
                     
% SL PROJECTION SETTINGS
    rcps = {'RCP3PD' 'RCP45' 'RCP85'}; % which RCPs to take - {'RCP3PD' 'RCP45' 'RCP85'}
    dig = []; % number of digits to use in order to save memory, else []
    magicc_realizations = 600; % Number of MAGICC realizations to use for each RCP; max & default == 600
    
% COUNTERFACTUAL/FACTUAL SL SETTINGS
    Tcf = {'mean' 'linrate' 'HadCRUT'}; % which Temperature input to use:'
                     % 'HadCRUT' - HadCRUT4 temperature,
                     % 'CMIP5_mean' - 'historicalNat_CMIP35_tas' mean from Gareth Jones <gareth.s.jones@metoffice.gov.uk>
                     % 'CMIP5' - 'historicalNat_CMIP35_tas' single simu (data column given as input) from Gareth Jones <gareth.s.jones@metoffice.gov.uk>
                     % 'mean' - mean temperature between 500-1800 CE
                     % 'linrate' - linear rate temperature between 500-1800 CE
     col = []; % Data columns of 'historicalNat_CMIP35_tas' to use, if Tcf=='CMIP5'

% RUN MAIN SCRIPT SESL
    fprintf('\n RUN MAIN SCRIPT SESL \n')
    SL = SESL(start,...
        'model', model,...
        'SL_dat', SL_dat,...
        'T_data', T_data,...
        'StartDistr', StartDistr,... 
        'all_flat_priors', all_flat_priors,...
        'AdjustPrior', AdjustPrior,...
        'a1_prior', a1_prior,...
        'a2_prior', a2_prior,...
        'tau1_prior', tau1_prior,...
        'tau2_prior', tau2_prior,...
        'c_prior', c_prior,...
        'c2000_prior', c2000_prior,...
        'T01st_prior', T01st_prior,...
        'T02nd_prior', T02nd_prior,...
        'TauLogUniform', TauLogUniform,...
        'optHo', optHo,...
        'OptimT0', OptimT0,...
        'UseMarT0', UseMarT0,... 
        'T0temp_level', T0temp_level,...
        'calibperiod', calibperiod,...
        'period', period,...
        'baseperiod', baseperiod,...
        'T0period', T0period,...
        'T0length', T0length,...
        'JumpDist', JumpDist,...
        'burning', burning,...
        'NumSkip', NumSkip,...
        'sample', sample,...
        'NormProb', NormProb,...
        'Tnum', Tnum,...
        'T_err', T_err,...
        'tau_ar1', tau_ar1,...
        'ryrs', ryrs,...
        'TerrSc', TerrSc,...
        'UseCov', UseCov,...
        'CovShrink', CovShrink,...
        'CovTau', CovTau,...
        'NoNegCov', NoNegCov,...
        'fac1', fac1); 
    save(fullfile('Out', [T_data '_' SL_dat]), 'SL');
    
% CALCULATE TIME SERIES FROM PARAMETERS
    num =  ceil((SL.settings.sample*SL.settings.Tnum)/storeNum); % number of file to be stored
    fprintf('\n CALCULATE TIME SERIES FROM PARAMETERS \n');
    fprintf('\t %d time series file(s) will be generated\n',num);
    % storeNum sl, T0, c & temperature curves are calculated from the posterior 
    % parameter set & saved to the folder 'Simu\' 
    Calc_SL_from_Param(SL,storeNum,'MH','Simu');                                               
    clear SL

% CALCULATE PERCENTILES OF SL etc.
    fprintf('\n CALCULATE PERCENTILES OF SL etc. \n');
    P = Calc_SESL_Prc(T_data, num,'Simu');
    save(fullfile('Out', [T_data, '_', SL_dat, '_prc']),'P');
    
% CALCULATE SL PROJECTIONS
    fprintf('\n CALCULATE SL PROJECTIONS \n');
    SLpr = Calc_SESLProjection(P,rcps,dig,magicc_realizations);
    save(fullfile('Out', [T_data, '_', SL_dat, '_proj']),'SLpr');
    clear SLpr

% CALCULATE SL COUNTERFACTUALS
    fprintf('\n CALCULATE SL COUNTERFACTUALS \n');
    for i = 1: size(Tcf,2)
        if strcmp(cell2mat(Tcf(i)),'CMIP5')
            for ii = 1:length(col)
                SLcf.([cell2mat(Tcf(i)) '_' num2str(col(ii))]) = Calc_SESLConterfact(P,cell2mat(Tcf(i)),col);
            end
        else
            SLcf.(cell2mat(Tcf(i))) = Calc_SESLConterfact(P,cell2mat(Tcf(i)),col);
        end
    end
    save(fullfile('Out', [T_data, '_', SL_dat, '_cf']),'SLcf');
    clear SLcf
end
