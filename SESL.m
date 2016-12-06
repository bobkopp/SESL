function out = SESL(P0,varargin) 

% Main semi-empirical sea-level (SESL) script to determine semi-empirical
% parameter distributions
%
% out = SESL(P0,varargin) 
%
% INPUT:
% - P0       -> startin parameter set [a1, a2, c, tau1, tau2, To(1)]. 
%               If P0=NaN it will be asigned via simualted annealing below.
%               Parameters are scaled automatically to have the same order of magnitude
% - varargin -> settings as defined by DefineSettings_SESL
%
% OUTPUT: out.
% - setting     -> all the setting from above plus some internal changes and additions
% - data        -> all data loaded, temperature & sea level 
% - StartParam  -> the parameter set from which the simulated annealing starts
% - lowerBound  -> lower bound for simulated annealing
% - upperBound  -> upper bound for simulated annealing
% - SimAnn      -> the output of the simulated annealing if applied
% - MH          -> Montecarlo-Hastings output if applied
%                  - MH.Params: selected parameters in columns [a1, a2, c, tau1, tau2, To(1)]
%                  - MH.alpha: the MH mixing
% - Slice       -> the output of Slice sampling if applied
% - TimeElapsed -> time needed for calculation in seconds   
    
    calc_lik_only = false; %if true: Calc_SL_SLice_MH2_ only calculates the likelihood 
    % of the last parameter draw from another run with calc_lik_only==false
    
    if ~calc_lik_only
        tic;
        set = DefineSettings_SESL(varargin);
        dat = LoadData_SESL(set);
        set = checksettings(dat,set); % Check that period, calibperiod, baseperiod
        % & T0period are not in conflict with temperature & SL data-set lengths. 
        % Check conflicts of ryrs & data.temp.yrs and data.temp.yrs & tau_lim(1).
        % Make prior pdfs.

        % truncate Mar temp for T0 calculation
        if set.UseMarT0
            ix = dat.temp.T0temp(:,1)>=set.T0period(1);
            dat.temp.T0temp = dat.temp.T0temp(ix,:); 
            dat.temp.T0burnin = sum(dat.temp.T0temp(:,1)<dat.temp.T(1,1)); % size of Mar temp only for T0 burnin in
        end
        printit(set); % output status to shell
        
        out.settings = set;
        out.data = dat;
        if isnan(P0)
            % --------------Bounds & start -----------
            % Starting parameter set 
            [P0, lb, ub, P0offset] = def_Start_Lim(set,dat);
            out.StartParam = P0;
            out.lowerBound = lb;
            out.upperBound = ub;

            %----------------------- simulannealbnd: 1st minimization -----------------
            out.SimAnn = SimulatedAnnealing(P0,P0offset,lb,ub,set,dat);
            out.SimAnn.TimeElapsed = toc;
            P0 = out.SimAnn.OptParam;
            tic;
        else
            fprintf('Params:  a1=%1.3f \n\t\t a2=%1.3f \n\t\t c=%1.3f \n\t\t tau1=%1.0f \n\t\t tau2=%1.0f \n\t\t T_0(1)=%1.2f\n',P0(1),P0(2),P0(3),10^P0(4),10^P0(5),P0(6));
        end
        
        % Check if parameter set has a likelihood > 0
        Lik = zeros(100,1);
        for iLik=1:100
            Lik(iLik) = target_distr([P0(1:3) 10^P0(4) 10^P0(5) P0(6)],set,dat,1);
        end
        if mean(Lik)==0; error('Likelihood of starting parameter set is zero, please find another one.');end
        
        % Start the Metropolis-Hastings or Slice Sampling.
        if set.sample>0
            %----------------------------- MH sampling --------------------------------
                 out.MH = MH_sampler(set,dat,P0);
            %--------------------------- Slice Sampling -------------------------------
        %       out.Slice = Slice_sampler(set,dat,P0); % It never correctly worked if Terr~='no' 
        end
        out.TimeElapsed = toc;
    else
        out = target_distr(P0.MH.Params(end,:),P0.settings,P0.data,0);
    end
end % end main function

function [P0, lb, ub, P0offset] = def_Start_Lim(set,dat)

% Define the starting parameter distribution for simulated annealing, 
% dependent on the semi-empirical model used.
%
% P0 = [a1, a2, c, log10(tau1), log10(tau2), T0(0)]
    P0 = set.StartDistr;
    if strcmp(set.model,'TwoTau')
        %lower and upper bounds for parameters
        lb = [0, 0, -1, log10(dat.temp.yrs), log10(50), -1];
        ub = [10, 10, 1, log10(10000), log10(3000) , 1];
    elseif strcmp(set.model,'ConstRate')
        %lower and upper bounds for parameters
        lb = [0, 0, -1, log10(dat.temp.yrs), log10(dat.temp.yrs), -1];
        ub = [10, 10, 1, log10(10000), log10(300) , 1];
    elseif strcmp(set.model,'CRdecay') % from Bob
        %lower and upper bounds for parameters
        lb = [0, 0, -1, log10(dat.temp.yrs), log10(50), -1];
        ub = [10, 10, 1, log10(5000), log10(20000) , 1];
    elseif strcmp(set.model,'CRovTau')
        %lower and upper bounds for parameters
        lb = [0, 0, -1, log10(dat.temp.yrs), log10(dat.temp.yrs), -1];
        ub = [10, 10, 1, log10(10000), log10(300) , 1];
    elseif strcmp(set.model,'simpel')
        %lower and upper bounds for parameters
        lb = [0, 0, -5, log10(dat.temp.yrs), log10(dat.temp.yrs), -1];
        ub = [10, 10, 5, log10(100000), log10(1000) , 1];
    end
    P0offset = min(lb)*2;
end

function SimAnn = SimulatedAnnealing(P0,P0offset,lb,ub,set,dat)
    
% Start the matlab simulated annealing function simulannealbnd to find an
% maximum likelihood parameter set to start the semi-empirical model from.

    fprintf('\nStart simulated annealing...')
    % -----------Probability function (-Likelihood), bounds & start --------
    set_ = set;

    ProbFunction = @(X) target_distr([X(1:3) 10^X(4) 10^X(5) X(6)],set_,dat,1);
    
    %--------------------- simulated annealing options ------------------------
    options = [];
    options = saoptimset(options,'Display','iter','DisplayInterval',400);
    options = saoptimset(options,'InitialTemperature',500,'TemperatureFcn',@temperaturefast,'MaxFunEval',8000,'TolFun',1e-3);
    
    % Find fval=min(1/Likelihood) and the corresponding parameters x0 with simulated annealing
    [x0,fval,exitFlag,output] = simulannealbnd(@(x) ProbFunction(10.^(x)+P0offset),log10(P0-P0offset),log10(lb-P0offset),log10(ub-P0offset),options);

    x0 = 10.^(x0)+P0offset;
    
    fprintf('The number of iterations was : %d\n', output.iterations);
    fprintf('The number of function evaluations was : %d\n', output.funccount);
    prob = -fval;
    
    fprintf('The maximum likelihood found was : %g\n', prob);
    
    fprintf('Optimal simulannealbnd Parameters: a1    = %g \n',x0(1))
    if strcmp(set.model,'TwoTau')
        fprintf('                                   a2    = %g \n',x0(2));end
    if strcmp(set.model,'ConstRate') || strcmp(set.model,'CRdecay')
        fprintf('                                   c     = %g \n',x0(3));end
    if strcmp(set.model,'CRovTau')
        fprintf('                                   c     = %g \n',x0(3)/10^x0(4));end
    fprintf('                                   tau1  = %g \n',10^x0(4))
    if strcmp(set.model,'TwoTau') || strcmp(set.model,'CRdecay')
        fprintf('                                   tau2  = %g \n',10^x0(5));end
    if set.OptimT0
        fprintf('                                   T0(1) = %g \n',x0(6))
    end
    
    SimAnn.StartParam = P0;
    SimAnn.OptParam = x0; 
    SimAnn.exitFlag = exitFlag;
    SimAnn.output = output;
    SimAnn.MaxLik = prob;
end

function MH = MH_sampler(set,dat,x0)

% Start the slightly adapted Matlab Metropolis-Hastings algorithm mhsample, 
% now called mhsample2 and added as a function below. 
% Here the final parameter set gets calculated.

    fprintf('Start MH Sampling: \n# Burnin: %g \n# Samples: %g \n# Thinning: %g \n',set.burning,set.sample, set.NumSkip)
    
    x0offset = -2; % to avoid taking the log of negative numbers
    % Probability function which outputs Likelihood
    ProbFunction = @(x) target_distr([x(1:3) 10^x(4) 10^x(5) x(6)],set,dat,0);

    pr = set.JumpDist;    
    proprnd = @(x) x + randn*pr;

    samplerset = struct('pdf', @(xx) ProbFunction(10.^(xx)+x0offset), ...
        'proprnd',proprnd, ...
        'burnin',set.burning, ...
        'thin',set.NumSkip, ...
        'symmetric',1);

    [smpl, accept] = mhsample2(log10(x0-x0offset),set.sample,samplerset);
   
    smpl = 10.^smpl+x0offset;
    smpl(:,4) = 10.^smpl(:,4);
    smpl(:,5) = 10.^smpl(:,5);
    MH.Params = smpl;
    MH.alpha = accept;  
end

function Slice = Slice_sampler(set,dat,x0)

% A theoretical alternative to the MH-sampler is the slice samplig method
% which did by now not work in this context.

    fprintf('Start Slice Sampling: \n# Burnin: %g \n# Samples: %g \n# Thinning: %g \n',set.burning,set.sample, set.NumSkip)
    % Probability function which outputs Likelihood
    ProbFunction = @(x) target_distr(x,set,dat,0);
    w = .005;% width -> default 10
     
    [rnd, neval] = slicesample2(x0,set.sample,'pdf',ProbFunction,'width',w,'burnin',set.burning,'thin',set.NumSkip);
    
    Slice.Params = rnd;
    Slice.neval = neval;
end

function prob = target_distr(param,set,dat,negLik)

% Here the target distribution for the above chosen sampling method gets
% evaluated by calculating a semi-empirical sea-level (calc_sl.m) and
% comparing it to the input sea-level data (calc_prob).

    % Add uncertainty to the temperature input
    temp = calc_temp(set,dat.temp);   
    % Calculate the equilibrium temperature
    [T01, T02] = calc_T0(set,temp,dat.temp,param);
    % Calculate the semi-empirical sea-level and a few parameters
    pri = calc_sl(set, dat.temp,temp,T01,T02,param);
    if strcmp(set.model,'CRdecay') && sum(pri.yr==2000) == 1
        c2000 = pri.c(pri.yr==2000,:);
    else
        c2000 = pri.c; % c is fixed
    end
    T0500_700 =mean(pri.T01(pri.yr>=500 & pri.yr<=700,:));
    T500_700 =mean(pri.T(pri.yr>=500 & pri.yr<=700,:));
    % Compare sl data input to semi-emp. sl -> Gaussian likelihood
    prob = calc_prob(pri,param,T0500_700,T500_700,c2000,set,dat);
    if negLik % simulated annealing is looking for a minimum
        prob = -mean(prob); 
    else
        prob = mean(prob);
    end

end % end calc_prop

function printit(settings)

% Print a few of the settings to the command window

    fprintf('\n The MODEL used is: ')
    if strcmp(settings.model,'simpel')
        fprintf('dSL/dt = a * (T(t)-T0(t)) \n')
    elseif strcmp(settings.model,'ConstRate')
        fprintf('dSL/dt = c + a * (T(t)-T0(t)) \n')
    elseif strcmp(settings.model,'CRdecay')
        fprintf('dSL/dt = c(t) + a * (T(t)-T0(t)) \n')
        fprintf('                    dc(t)/dt = c(t)/tau_ \n')
    elseif strcmp(settings.model,'CRovTau')
        fprintf('dSL/dt = c/tau + a * (T(t)-T0(t)) \n')
    elseif strcmp(settings.model,'TwoTau')
        fprintf('dSL/dt = a1 * (T(t)-T0_1(t)) + a2 * (T(t)-T0_2(t)) \n')
    end
    fprintf('                    dT0(t)/dt=(T-T0)/tau \n')
    fprintf(' -> Calib. period: %4d - %4d \n',settings.calibperiod(1),settings.calibperiod(end))
    fprintf(' -> Calc. period : %4d - %4d \n',settings.period(1),settings.period(end))
    fprintf(' -> Base  period : %4d - %4d \n',settings.baseperiod(1),settings.baseperiod(end))
    fprintf(' -> To(1) period : %4d - %4d \n',settings.T0period(1),settings.T0period(end))

    if ~settings.OptimT0
        fprintf('    T0(1) is not optimized T0(1) = mean(T(T0period))\n')
    end

    fprintf('-> "%s" temperature + "%s" error \n',settings.T_data,settings.T_err);
    if settings.TerrSc>1
        fprintf('    -> T unc. scaled by %1.0f !!!\n',settings.TerrSc);
    end
    if strcmp(settings.T_err,'default')
        fprintf('   added every %2d years \n',settings.ryrs);
    elseif strcmp(settings.T_err,'ar1ts')
        fprintf('   with ar1 timescale = %2d years \n',settings.tau_ar1);
    end
    
    fprintf('-> "%s" SL proxy \n',settings.SL_dat);
    
    if ~settings.optHo
        fprintf('-> H_0 is NOT optimal \n')
    end
    
    if settings.useCov
        fprintf('-> COV MAT for MC update\n');
        if settings.CovShrink ~=0
            fprintf('   Cov Mat scaling f: %4.0f (C_ = f*C0 + (1-f)*C) \n',settings.CovShrink);
        end
        if ~isnan(settings.CovTau)
            fprintf('   Cov Mat scaling tau: %4.0f (C_ = C.*exp(-delta(t)/CovTau)) \n',settings.CovTau);
        end
        if settings.NoNegCov
            fprintf('   Negative terms will be REMOVED from the Cov Mat \n')
        end
    end
end

function set = checksettings(dat,set)

% check that period, calibperiod, baseperiod & T0period are not in
% conflict with temperature & SL data. Make Parameter prior PDFs.
    
    % CHECK LENGTH OF PERIOD
    if set.UseMarT0
        Tpath = 'T0temp';
    else
        Tpath = 'T';
    end
    if set.period(1)<round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2))
        period(1) = round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2));
    else
        period(1) = set.period(1);
    end
    
    if set.period(end)>round(dat.temp.T(end,1)+ceil((dat.temp.yrs-1)/2))
        period = period(1):round(dat.temp.T(end,1)+ceil((dat.temp.yrs-1)/2));
    else
        period = period(1):set.period(end);
    end
    if ~(isequal(period,set.period))
       fprintf('\n!!!WARNING!!! "period" changed from %4.0f-%4.0f CE to %4.0f-%4.0f CE',set.period(1),set.period(end),period(1),period(end))
       set.period = period;
    end

    % CHECK LENGTH OF T0-PERIOD
    if set.T0length==0
        if set.T0period(1)<round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2))
            T0period(1) = round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2));
        else
            T0period(1) = set.T0period(1);
        end
        if set.T0period(end)>round(dat.temp.(Tpath)(end,1)+ceil((dat.temp.yrs-1)/2))
            T0period = T0period(1):round(dat.temp.(Tpath)(end,1)+ceil((dat.temp.yrs-1)/2));
        elseif set.T0period(end)<round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2))
            T0period = T0period(1):T0period(1)+length(set.T0period);
        else
            T0period = T0period(1):set.T0period(end);
        end
        if ~(isequal(T0period,set.T0period))
           fprintf('\n!!!WARNING!!! "T0period" changed from %4.0f-%4.0f CE to %4.0f-%4.0f CE',set.T0period(1),set.T0period(end),T0period(1),T0period(end))
           set.T0period = T0period;
        end
    else
        T0period = (round(dat.temp.(Tpath)(1,1)):round(dat.temp.(Tpath)(1,1))+set.T0length)-floor((dat.temp.yrs-1)/2);
        fprintf('\n!!!WARNING!!! "T0period" changed according to T0 length from %4.0f-%4.0f CE to %4.0f-%4.0f CE',set.T0period(1),set.T0period(end),T0period(1),T0period(end))
        set.T0period = T0period;
    end
    
    % CHECK LENGTH OF CALIBPERIOD
    if set.calibperiod(1)<max([round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2)),round(dat.sea.proxy.yr(1))]);
        calibperiod(1) = max([round(dat.temp.(Tpath)(1,1)-floor((dat.temp.yrs-1)/2)),round(dat.sea.proxy.yr(1))]);
    else
        calibperiod(1) = set.calibperiod(1);
    end

    if set.calibperiod(end)>min([round(dat.temp.(Tpath)(end,1)+ceil((dat.temp.yrs-1)/2)),round(dat.sea.proxy.yr(end))]);
        calibperiod = calibperiod(1):min([round(dat.temp.(Tpath)(end,1)+ceil((dat.temp.yrs-1)/2)),round(dat.sea.proxy.yr(end))]);
    else
        calibperiod = calibperiod(1):set.calibperiod(end);
    end
    if ~(isequal(calibperiod,set.calibperiod))
       fprintf('\n!!!WARNING!!! "calibperiod" changed from %4.0f-%4.0f CE to %4.0f-%4.0f CE',set.calibperiod(1),set.calibperiod(end),calibperiod(1),calibperiod(end))
       set.calibperiod = calibperiod;
    end
    
    %CHECK LENGTH OF BASEPERIOD
    if set.baseperiod(1)<round(dat.sea.proxy.yr(1))
        baseperiod(1) = round(dat.sea.proxy.yr(1));
    else
        baseperiod(1) = set.baseperiod(1);
    end
    if set.baseperiod(end)>round(dat.sea.proxy.yr(end))
        baseperiod = baseperiod(1):round(dat.sea.proxy.yr(end));
    else
        baseperiod = baseperiod(1):set.baseperiod(end);
    end
    if ~(isequal(baseperiod,set.baseperiod))
       fprintf('\n!!!WARNING!!! "baseperiod" changed from %4.0f-%4.0f CE to %4.0f-%4.0f CE',set.baseperiod(1),set.baseperiod(end),baseperiod(1),baseperiod(end))
       set.baseperiod = baseperiod;
    end    
    
    %CHECK Tnum & T_err
    if (strcmp(set.T_err,'no') || strcmp(set.T_err,'default')) && set.Tnum>1
        fprintf('\n!!!WARNING!!! "Tnum" changed from %4.0f to %4.0f',set.Tnum,1)
        set.Tnum = 1;
    end    
    
    % ryrs should not be smaller than temperature time steps
    if ~(isequal(set.ryrs/dat.temp.yrs,round(set.ryrs/dat.temp.yrs)))
        if set.ryrs<dat.temp.yrs;
            ryrs = dat.temp.yrs;
        else
            if mod(set.ryrs,dat.temp.yrs)>=set.ryrs/2
                ryrs = set.ryrs+mod(set.ryrs,dat.temp.yrs);
            else
                ryrs = set.ryrs-mod(set.ryrs,dat.temp.yrs);
            end
        end
        fprintf('\n!!!WARNING!!! "ryrs" changed from %4.0f to %4.0f',set.ryrs,ryrs)
        set.ryrs = ryrs;
    end
    
    % Depending on the temp. uncertainty setting, tau_ar1 needs to be >= or <= 1. 
    % One time it is the AR(1) timescale the other time it represents the AR(1) parameter 
    if strcmp(set.T_err,'ar1ts') && set.tau_ar1<1
        error('For the temperature uncertainty setting ar1ts, tau_ar1 needs to be >=1. Please adjust!')
    elseif strcmp(set.T_err,'ar1') && set.tau_ar1>1
        error('For the temperature uncertainty setting ar1, tau_ar1 needs to be <=1. Please adjust!')
    end

    % Adjust param prior according to model
    if set.AdjustPrior
        fprintf('\n!!!WARNING!!! The PRIORS get ADJUSTED to the model used.')
        if strcmp(set.model,'TwoTau')
            set.a1_prior = {'uniform' 0 10};
            set.a2_prior = {'uniform' 0 10};
            if set.TauLogUniform
                set.tau1_prior = {'uniform',log(1000),log(10000)}; % log uniform distribution 
                set.tau2_prior = {'uniform',log(dat.temp.yrs),log(3000)}; % log uniform distribution 
            else
                set.tau1_prior = {'uniform',301,10000}; %  uniform distribution 
                set.tau2_prior = {'uniform',dat.temp.yrs,300}; %  uniform distribution 
            end
            set.c_prior = {'flat', [],[]};
        elseif strcmp(set.model,'ConstRate')
            set.a1_prior = {'uniform' 0 10};
            set.a2_prior = {'flat' [] []};
            if set.TauLogUniform
                set.tau1_prior = {'uniform',log(dat.temp.yrs),log(10000)}; % log uniform distribution 
                set.tau2_prior = {'flat',[],[]}; % log uniform distributio
            else
                set.tau1_prior = {'uniform',dat.temp.yrs,10000}; %  uniform distribution
                set.tau2_prior = {'flat',[],[]};  
            end
            set.c_prior = {'uniform',-1,1};      
        elseif strcmp(set.model,'CRdecay')
            set.a1_prior = {'uniform' 0 2};
            set.a2_prior = {'flat' [] []};
            if set.TauLogUniform
                set.tau1_prior = {'uniform',log(30),log(3000)}; % log uniform distribution 
                set.tau2_prior = {'uniform',log(1000),log(20000)}; % log uniform distribution 
            else
                set.tau1_prior = {'uniform',30,3000}; %  uniform distribution 
                set.tau2_prior = {'uniform',1000,10000}; %  uniform distribution 
            end
        elseif strcmp(set.model,'CRovTau')
            set.a1_prior = {'uniform' 0 10};
            set.a2_prior = {'flat' [] []};
            if set.TauLogUniform
                set.tau1_prior = {'uniform',log(dat.temp.yrs),log(10000)}; % log uniform distribution 
                set.tau2_prior = {'flat',[],[]}; 
            else
                set.tau1_prior = {'uniform',dat.temp.yrs,10000}; %  uniform distribution 
                set.tau2_prior = {'flat',[],[]}; 
            end
            set.c_prior = {'uniform', -1, 1}; 
        elseif strcmp(set.model,'simpel')
            set.a1_prior = {'uniform' 0 10};
            set.a2_prior = {'flat' [] []};
            if set.TauLogUniform
                set.tau1_prior = {'uniform',log(dat.temp.yrs),log(10000)}; % log uniform distribution 
                set.tau2_prior = {'flat',[],[]}; % log uniform distribution 
            else
                set.tau1_prior = {'uniform',dat.temp.yrs,10000}; % uniform distribution 
                set.tau2_prior = {'flat',[],[]}; 
            end
            set.c_prior = {'flat',[],[]}; 
        end           
            set.T0_prior = {'normal',0,.2}; 
    end
    
    % Check that adjustment time scale is not smaller than timesteps 
    if set.TauLogUniform
        if cell2mat(set.tau1_prior(2))<log(dat.temp.yrs)
            set.tau1_prior{2} = log(dat.temp.yrs);
        end
        if cell2mat(set.tau2_prior(2))<log(dat.temp.yrs)
            set.tau2_prior{2} = log(dat.temp.yrs);
        end
    else
        if cell2mat(set.tau1_prior(2))<(dat.temp.yrs)
            set.tau1_prior{2} = (dat.temp.yrs);
        end
        if cell2mat(set.tau2_prior(2))<(dat.temp.yrs)
            set.tau2_prior{2} = (dat.temp.yrs);
        end
    end
    
    % define Prior likelihood functions (not log likelihoods!)
    if strcmp(set.a1_prior(1),'flat')
        set.lik_a1 = @(x) 1;
    else
        set.lik_a1 = @(x) pdf(cell2mat(set.a1_prior(1)),x,cell2mat(set.a1_prior(2)),cell2mat(set.a1_prior(3)));
    end
    if strcmp(set.a2_prior(1),'flat')
        set.lik_a2 = @(x) 1;
    else
        set.lik_a2 = @(x) pdf(cell2mat(set.a2_prior(1)),x,cell2mat(set.a2_prior(2)),cell2mat(set.a2_prior(3)));
    end
    if strcmp(set.c_prior(1),'flat')
        set.lik_c = @(x) 1;
    else
        set.lik_c = @(x) pdf(cell2mat(set.c_prior(1)),x,cell2mat(set.c_prior(2)),cell2mat(set.c_prior(3)));
    end
    if strcmp(set.c2000_prior(1),'flat')
        set.lik_c2000 = @(x) 1;
    else
        set.lik_c2000 = @(x) pdf(cell2mat(set.c2000_prior(1)),x,cell2mat(set.c2000_prior(2)),cell2mat(set.c2000_prior(3)));
    end
    if strcmp(set.tau1_prior(1),'flat')
        set.lik_tau1 = @(x) 1;
    else
        set.lik_tau1 = @(x) pdf(cell2mat(set.tau1_prior(1)),x,cell2mat(set.tau1_prior(2)),cell2mat(set.tau1_prior(3)));
    end
    if strcmp(set.tau2_prior(1),'flat')
        set.lik_tau2 = @(x) 1;
    else
        set.lik_tau2 = @(x) pdf(cell2mat(set.tau2_prior(1)),x,cell2mat(set.tau2_prior(2)),cell2mat(set.tau2_prior(3)));
    end
    if strcmp(set.T01st_prior(1),'flat')
        set.lik_T01st = @(x) 1;
    else
        set.lik_T01st = @(x) pdf(cell2mat(set.T01st_prior(1)),x,cell2mat(set.T01st_prior(2)),cell2mat(set.T01st_prior(3)));
    end
    if strcmp(set.T02nd_prior(1),'flat')
        set.lik_T02nd = @(x,y) 1;
    else
        set.lik_T02nd = @(x, y) pdf(cell2mat(set.T02nd_prior(1)),x,y,cell2mat(set.T02nd_prior(3)));
    end

    
    % if parameters are not needed for certain models, remove prior.
    if (strcmp(set.model,'ConstRate') || strcmp(set.model,'CRdecay') || strcmp(set.model,'CRovTau') || strcmp(set.model,'simpel'))
        set.lik_a2 = @(x) 1; end
    if strcmp(set.model,'TwoTau') || strcmp(set.model,'simpel')
        set.lik_c = @(x) 1; end
    if strcmp(set.model,'ConstRate') || strcmp(set.model,'CRovTau') || strcmp(set.model,'simpel')
        set.lik_tau2 = @(x) 1; end
         
    fprintf('\n')
end

function [smpl,accept] = mhsample2(start,nsamples,optionstruct)

pdf = optionstruct.pdf;
proprnd = optionstruct.proprnd;
burnin = optionstruct.burnin;
thin = optionstruct.thin;
sym = optionstruct.symmetric;

if ~isfield(optionstruct, 'nchain')
    nchain = 1;
end

% log density is preferred for numerical stability
logpdf = @(x) mylog(pdf(x));
logproppdf = @(x,y) mylog(proppdf(x,y));

outclass = superiorfloat(start); % single or double

% Put the replicates dimension second.
distnDims = size(start,2);
smpl = zeros([nsamples,nchain,distnDims],outclass);

x0 = start;  %x0  is the place holder for the current value
accept =zeros(nchain,1,outclass);
% Metropolis-Hasting Algorithm.
U = log(rand(nchain,nsamples*thin+burnin));
for i = 1-burnin:nsamples*thin
    y = zeros(size(x0));
    for ii = 1:size(x0,2)
        y(ii) = proprnd(x0(ii));
    end
    rho =nan;
    while isnan(rho)
        if ~sym
            q1 = logproppdf(x0,y);
            q2 = logproppdf(y,x0);
            % this is a generic formula.
            rho = (q1+logpdf(y))-(q2+logpdf(x0));
        else
            %save the evaluation time for symmetric proposal dist'n
            rho = logpdf(y)-logpdf(x0);
        end
    end
    % Accept or reject the proposal.
    Ui = U(:,i+burnin);
    acc = Ui<= min(rho,0);
    x0(acc,:) = y(acc,:); % preserves x's shape.
    accept = accept+(acc);
    if mod(i,100)==0 && i<0
        fprintf('%1.0f  ',i)
        if mod(i,1000)==0
            fprintf('\n')
        end
    end
    if i>0 && mod(i,thin)==0; % burnin and thin
        if mod(i,100)==0
            fprintf('%1.0f  ',i)
            if mod(i,1000)==0
                fprintf('\n')
            end
        end
        smpl(i/thin,:,:) = x0;
    end
end
fprintf('\n')
% Accept rate can be used to optimize the choice of scale parameters in
% random walk MH sampler. See for example Roberts, Gelman and Gilks (1997).
accept = accept/(nsamples*thin+burnin);

% Move the replicates dimension to the end to make samples easier to
% manipulate.
smpl = permute(smpl,[1 3 2]);

end
function  y = mylog(x)
% my log function is to define to avoid the warnings.
y = -Inf(size(x));
y(x>0) = log(x(x>0));
end
