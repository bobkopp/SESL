function [year, temp, T01, T02] = resize_T(set,tempr,temp, T01, T02)

% If the input temperature is not given in yearly values, here it will, 
% alongside with the equilibrium temperatures, be resized to a step-like 
% yearly time series
%
% [year temp T01 T02] = resize_T(set,tempr,temp, T01, T02)
% 
% INPUT:
% - set   -> general settings
% - tempr -> original temperature input
% - temp  -> temperature draw, dependent on the temperature error model
% - T01   -> equilibrium temperature 1
% - T02   -> equilibrium temperature 2, only if strcmp(model,'TwoTau')==1
% 
% OUTPUT:
% - year  -> year;
% - temp -> step-like yearly temperature;
% - T01  -> step-like yearly equilibrium temperature 1
% - T02  -> step-like yearly equilibrium temperature 2, only if strcmp(model,'TwoTau')==1

    if set.UseMarT0
        Tpath = 'T0temp';
    else
        Tpath = 'T';
    end

    % Turn temp & T0 to steplike yearly data and truncate to desired length 
    fyr = max([tempr.T(1,1) min([set.period(1) set.calibperiod(1)])]);
    lyr = max([set.period(end) set.calibperiod(end)]);
    
    year = repmat(tempr.T(tempr.T(:,1)>=fyr & tempr.T(:,1)<=lyr,1),tempr.yrs,1);
    [year, ix] = sort(year);
    year1 = [year(1)-floor((tempr.yrs-1)/2) : year(end)+ceil((tempr.yrs-1)/2)];%tempr.yr;
    
    temp1 = temp(tempr.(Tpath)(:,1)>=max([tempr.T(1,1) fyr]),:);
    temp1 = repmat(temp1,tempr.yrs,1);
    temp1 = temp1(ix,:);
    
    T01_1 = T01(tempr.(Tpath)(:,1)>=max([tempr.T(1,1) fyr]),:);
    T01_1 = repmat(T01_1,tempr.yrs,1);
    T01_1 = T01_1(ix,:);
    if strcmp(set.model,'TwoTau')
        T02_1 = T02(tempr.(Tpath)(:,1)>=max([tempr.T(1,1) fyr]),:);
        T02_1 = repmat(T02_1,tempr.yrs,1);
        T02_1 = T02_1(ix,:);
    end
    
    % Turn temp & T0 to steplike yearly data and truncate to desired length 
    fyr_ = max([tempr.(Tpath)(1,1),min([set.period(1) set.calibperiod(1)])]);
    lyr_ = fyr-1;
    
    if fyr~=fyr_ && set.UseMarT0 % if Marcott temperature is used to spin up T0
        year = repmat(tempr.T0temp(tempr.T0temp(:,1)>=fyr_ & tempr.T0temp(:,1)<=lyr_,1),tempr.T0yrs,1);
        [year, ix] = sort(year);
        year2 = [year(1)-floor((tempr.T0yrs-1)/2) : year(end)+ceil((tempr.T0yrs-1)/2)];%tempr.yr;

        temp2 = temp(tempr.(Tpath)(:,1)>=fyr_& tempr.(Tpath)(:,1)<=lyr_,:);
        temp2 = repmat(temp2,tempr.T0yrs,1);
        temp2 = temp2(ix,:);

        T01_2 = T01(tempr.(Tpath)(:,1)>=fyr_& tempr.(Tpath)(:,1)<=lyr_,:);
        T01_2 = repmat(T01_2,tempr.T0yrs,1);
        T01_2 = T01_2(ix,:);
        if strcmp(set.model,'TwoTau')
            T02_2 = T02(tempr.(Tpath)(:,1)>=fyr_& tempr.(Tpath)(:,1)<=lyr_,:);
            T02_2 = repmat(T02_2,tempr.T0yrs,1);
            T02_2 = T02_2(ix,:);
        end
       
        temp = [temp2(year2<year1(1),:); temp1];

        T01 = [T01_2(year2<year1(1),:); T01_1];
        if strcmp(set.model,'TwoTau')
            T02 = [T02_2(year2<year1(1),:); T02_1];
        end
        year = [year2(year2<year1(1)), year1];
    else
        temp = temp1;
        year = year1;
        T01 = T01_1;
        if strcmp(set.model,'TwoTau')
            T02 = T02_1;
        end
        
    end
    
end
