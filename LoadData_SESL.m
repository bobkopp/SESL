function d = LoadData_SESL(set)

% Load input temperature and sea-level data.
% 
% d = LoadData_SESL(set)
% 
% INPUT: general settings from DefineSettings_SESL.m
% 
% OUTPUT: data structure


    folder = 'Data';
    
    %%%%%%%%%%%%%%%%%%%%%%% SEA LEVEL in cm %%%%%%%%%%%%%%%%%%%%%%%

    SL_dat = set.SL_dat;
    
    load(fullfile(folder,[SL_dat, '.mat'])); 

    d.sea.proxy.yr = sl(:,1)';
    d.sea.proxy.sl = sl(:,2)'/10;
    d.sea.proxy.slerr = sl(:,3)'/10;

    C = C/100;
    C = C + eye(length(C))*eps;

    if set.useCov
        % Scale Cov Mat and remove negative elements
        if ~isnan(set.CovTau) && exist('C','var')==1
            Csc = exp(-abs(bsxfun(@minus,d.sea.proxy.yr,d.sea.proxy.yr')/set.CovTau));
            d.sea.proxy.C = C.*Csc;
        else
            d.sea.proxy.C = set.CovShrink*diag(diag(C))+ (1-set.CovShrink)*C;
        end
        if set.NoNegCov && exist('C','var')==1
            ix = d.sea.proxy.C<0;
            d.sea.proxy.C(ix) = 0;
        end
    end
    
    % offset proxy SL data
    offset = mean(d.sea.proxy.sl(d.sea.proxy.yr>=set.baseperiod(1) & d.sea.proxy.yr<=set.baseperiod(end))); % offset to baseperiod
    d.sea.proxy.sl = d.sea.proxy.sl - offset; 
    
    
    %%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%
    
    temperature = load(fullfile(folder, set.T_data)); % Load Temp 
    d.temp.T = temperature.T;

    dyr = diff(d.temp.T(:,1));
    if dyr(1)==mean(dyr)
        d.temp.yrs = dyr(1);
    else
        error('Temperature not recorded in constant steps')
    end
    d.temp.T(:,3) = d.temp.T(:,3)*set.TerrSc;
    
    if set.UseMarT0 %Use Marcott RegEM temperature to spin up T0
        temperature = load(fullfile(folder, 'Marcott13_RegEM-HC3_20'));

        d.temp.T0yrs = 20; % year steps
        T = temperature.T;
        
        % bring Marcott temp. (for T0 calculation) to one level whith the
        % other temperature's first 'T0temp_level' years of data
        ix = T(:,1)>=d.temp.T(1,1) & T(:,1)<=d.temp.T(1,1)+set.T0temp_level;
        offs = mean(d.temp.T(d.temp.T(:,1)<=d.temp.T(1,1)+set.T0temp_level,2));
        T(:,2) = T(:,2) - mean(T(ix,2)) + offs;
        
        % merge Mar and other temperature
        ix = T(:,1)<d.temp.T(1,1)-floor((d.temp.yrs-1)/2);
        T = T(ix,:);     
        d.temp.T0temp = [T; d.temp.T(:,1:3)];
     end

    
end % end loaddata