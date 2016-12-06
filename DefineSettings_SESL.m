function set = DefineSettings_SESL(in)

% Defines the general settings.
% 
% set = DefineSettings_SESL(in)
% 
% INPUT: Pairs of one string, defining the setting, and its value.
%        Example: set = DefineSettings_SESL('model','CRdecay','Tnum',100,'useCov',false)
%        
% OUTPUT: Structure with all major settings

P = inputParser;

% Model formulation: 
% simpel    => dSL/dt = a * (T(t)-T0(t))
% TwoTau    => dSL/dt = a1 * (T(t)-T0_1(t)) + a2 * (T(t)-T0_2(t))
% ConstRate => dSL/dt = c + a * (T(t)-T0(t))
% CRdecay   => dSL/dt = c(t) + a * (T(t)-T0(t))
% CRovTau    => dSL/dt = c/tau + a * (T(t)-T0(t))
%   -> dT0(t)/dt = (T-T0)/tau
%   -> dc(t)/dt = -c(t)/tau
defmodel = 'CRdecay';
    valmodel = {'ConstRate' 'CRdecay' 'CRovTau' 'TwoTau' 'simpel'};
    checkmodel = @(x) any(validatestring(x,valmodel));
    addParamValue(P,'model',defmodel,checkmodel)
    
% Prior distribution for Parameters -> important for Likelihood
defall_flat_priors = false;
    addParamValue(P,'all_flat_priors',defall_flat_priors,@islogical);
% adjust the prior in check_settings, according to the model used 
defAdjustPrior = false;
    addParamValue(P,'AdjustPrior',defAdjustPrior,@islogical);
    
defStartDistr = [0.5, 0.3, 0, log10(50), log10(100) , 0]; % Start distribution for simulated annealing
    addParamValue(P,'StartDistr',defStartDistr,@isnumeric);

defa1_prior = {'uniform' 0 2};
    addParamValue(P,'a1_prior',defa1_prior,@iscell);
defa2_prior = {'uniform' 0 2}; % a2(if model=='TwoTau'), b(if useb==true)
    addParamValue(P,'a2_prior',defa2_prior,@iscell);
deftau1_prior = {'uniform',log(30),log(3000)}; % log uniform distribution 
    addParamValue(P,'tau1_prior',deftau1_prior,@iscell);
deftau2_prior = {'uniform',log(1000),log(20000)}; % log uniform distribution 
    addParamValue(P,'tau2_prior',deftau2_prior,@iscell);
defc_prior = {'flat', [],[]}; % flat prior
    addParamValue(P,'c_prior',defc_prior,@iscell);
defc2000_prior = {'uniform', -.2,.2};
    addParamValue(P,'c2000_prior',defc2000_prior,@iscell);
defT01st_prior = {'uniform',-.6,.6}; % T0 at the beginning of T if it is before 500CE 
    addParamValue(P,'T01st_prior',defT01st_prior,@iscell); 
defT02nd_prior = {'normal',0,.2}; % T0(500-700) prior
    addParamValue(P,'T02nd_prior',defT02nd_prior,@iscell); 

defTauLogUniform = true; % true => Tau priors are loguniform, else uniform
    addParamValue(P,'TauLogUniform',defTauLogUniform,@islogical);

% optimize T0(1) (true) or set to mean(T(T0period))
defOptimT0 = true;
    addParamValue(P,'OptimT0',defOptimT0,@islogical);

% Jumping distribution size
defJumpDist = .005;
    addParamValue(P,'JumpDist',defJumpDist,@isnumeric);
    
% Number of ar1 temperatures to evaluate for likelihood
defTnum= 100;
    addParamValue(P,'Tnum',defTnum,@isnumeric);
    
% Sea Level data to use
defSL_dat ='GLMW-1ts'; % which SL proxy should be used? Add dummy data to folder '\Data'
    valSL_dat = {'GLMW','GLMW-1amp1ts','GLMW-1ts','GLMW-Gr','GLMW-NC'};
    checkSL_dat = @(x) any(validatestring(x,valSL_dat));
    addParamValue(P,'SL_dat',defSL_dat,checkSL_dat)
                                 
defCalibperiod = -1000:2010; % years to include for model calibration. KE11[-120 - 2000 AD], Mann Temp. [500 - 2006 AD]
    addParamValue(P,'calibperiod',defCalibperiod,@isnumeric);
defPeriod = -1000:2010; % years of data to include for SL calculation
    addParamValue(P,'period',defPeriod,@isnumeric);
defBaseperiod = 1400:1800; % Period for which simulated sl and data are normed
    addParamValue(P,'baseperiod',defBaseperiod,@isnumeric);
defT0period =-2000:-1800; % period to initialise T0(1)
    addParamValue(P,'T0period',defT0period,@isnumeric);
defT0length = 0; % if==0: T0period is used, else T0(1) = mean(T(first T0length yrs))
    addParamValue(P,'T0length',defT0length,@isnumeric);

defUseMarT0 = true; % use Marcott temperature to calc. T0 until 'T_data' starts
    addParamValue(P,'UseMarT0',defUseMarT0,@islogical);
defT0temp_level = 100; % number of years over which to level Marcott and T_data if UseMarT0==true
    addParamValue(P,'T0temp_level',defT0temp_level,@isnumeric)
    
defOptHo = true; % for each MC sample, calculate the optimal offset between simulated and calibration sea level via MLS
    addParamValue(P,'optHo',defOptHo,@islogical);

defT_data = 'Mann09_11'; % Choose temperature data for calibration
    valT_data = {'Mann08_eiv','Mann08_cps','Mann09','Mann09_11','Marcott13_RegEM-HC3_20','PAGES2k_13'};
    checkT_data = @(x) any(validatestring(x,valT_data));
    addParamValue(P,'T_data',defT_data,checkT_data);
            % 'Mann08_eiv'     -> Mann et al. PNAS 2008 eiv reco (default)
            % 'Mann08_cps'     -> Mann et al. PNAS 2008 cps reco
            % 'Mann09'         -> From Scott Rutherford from Mann et al. 2009 Science                 
            % 'Mann09_11'      -> as above but 11 year averages.
            % 'Marcott13_RegEM'-> Marcott et al. 2013 + HadCRUT3v (20yr mean) starting 9350BC
            % 'PAGES2k13'      -> Pages2K -9-2000AD
            % 'Pages2k'        -> Pages2K -10-1909AD, from 1910-2010AD mean(Pages2k(1850-1909AD))
      
defBurnin = 1000; % number of monte carlo samples used in the burning-in or spin-up period
    addParamValue(P,'burning',defBurnin,@isnumeric);
defNumSkip =100; % only every NumSkip sample will be selected because of autocorrelation
    addParamValue(P,'NumSkip',defNumSkip,@isnumeric);
defSample = 5000; % number of monte carlo samples
    addParamValue(P,'sample',defSample,@isnumeric);
    
defT_err = 'ar1ts'; % How should the temperature variance be treated ?
    valT_err = {'default','no','ar1','ar1ts'};
    checkT_err = @(x) any(validatestring(x,valT_err));
    addParamValue(P,'T_err',defT_err,checkT_err);
        % 'default' -> T + random noise as in KE11
        % 'no'  -> Don't add uncertainty
        % 'ar1'     -> T as AR(1) process with sig as 'default'
        % 'ar1ts'   -> AR(1) Parameter timescale == exp(-abs(t2-t1)/timescale)
deftau_ar1 = 10; % if Terr==ar1ts: timescale of ar1ts uncertainty above
                 % if Terr==ar1: AR(1) parameter
    addParamValue(P,'tau_ar1',deftau_ar1,@isnumeric);
defRyrs = 10; % If T_err==default: Number of years for which to add same noise on temp. in calc_sl
    addParamValue(P,'ryrs',defRyrs,@isnumeric);
        % 10 -> Mann et al. default (Kemp et al. 2011)
        % 20 -> Marcott et al. time series are 20yr means
        % 30 -> PAGES2k time series are 30yr means
defTerrSc = 1; % scaling of temp. uncertainty 
    addParamValue(P,'TerrSc',defTerrSc,@isnumeric);
         
defNormProb = false; % take into account the normfactor 1/sqrt(2*pi*|C|)
    addParamValue(P,'NormProb',defNormProb,@islogical);

defUseCov = true; % if existing, use covariance matrix of SL data to calc likelihood
    addParamValue(P,'useCov',defUseCov,@islogical);
defCovShrink = 0; % Sea level covariance matrix (C) shrikage parameter f. C gets 'shrinked' to C' while C0 = diag(diag(C)): C' = f*C0 + (1-f)*C
    addParamValue(P,'CovShrink',defCovShrink,@isnumeric);
defCovTau = 100; % if ~isnan(CovTau): take the elementwise product of the covariance and a tapering function exp(-delta(t)/CovTau) where CovTau is a time scale
    addParamValue(P,'CovTau',defCovTau,@isnumeric);
defNoNegCov = true; % remove negative Cov Mat entries and replace by 0
    addParamValue(P,'NoNegCov',defNoNegCov,@islogical);

defFac1 = 1; 
    addParamValue(P,'fac1',defFac1,@isnumeric);
        % adjustment factor for autocorrelated NC sea level data as 
        % described in KE11 "Bayesian updating to estimate the model
        % parameters" (default = 10). Gets multiplied to SL uncertainty
    
parse(P,in{:});

set = P.Results;
end

