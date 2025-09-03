%% PROCESS AND FILTER
clear all
cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/0_results')
load('grid4_results_apr29_v1.mat')
outname = 'grid4_processed_monthly_jun26.mat';

cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/radiative_kernel_data');
load('era5_ker_grid4_v1.mat');


for p = 1:2
    if p == 1; output = lakeoutput; end
    if p == 2; output = landoutput; end

%STEP 0: convert to NaN
for i = 1:length(output)
    t = output(i).snowa_ts;
    output(i).snowa_ts(t == -999) = NaN;
    t = output(i).watera_ts;
    output(i).watera_ts(t == -999) = NaN;
    t = output(i).snowa_std;
    output(i).snowa_std(t == -999) = NaN;
    t = output(i).watera_std;
    output(i).watera_std(t == -999) = NaN;
end


%STEP 0: ADD IN PERCENTAGES
    for i = 1:length(output)
        clear r
        r(:,1) = output(i).cloud_ts;
        r(:,2) = output(i).snow_ts;
        r(:,3) = output(i).water_ts;
        output(i).all_ts = r;
    end

    for i = 1:length(output)
        output(i).sum = sum(output(i).all_ts,2);
        r = output(i).all_ts;
        output(i).cloudp = 100*r(:,1)./output(i).sum;
        output(i).snowp = 100*r(:,2)./sum(r(:,2:3),2);
        output(i).waterp = 100*r(:,3)./sum(r(:,2:3),2);
    end


%STEP 1: REMOVE OBVIOUS OUTLIERS
for i = 1:length(output)
    ts = output(i).snowa_ts;
    II = find(ts < 0 | ts > 1);
    if isempty(II) == 0
    output(i).snowa_ts(II) = NaN;
    output(i).snowa_std(II) = NaN;
    end
    ts = output(i).watera_ts;
    II = find(ts < 0 | ts > 1);
    if isempty(II) == 0
    output(i).watera_ts(II) = NaN;
    output(i).watera_std(II) = NaN;
    end
end

%STEP 2: interpolate to daily after filtering
for i = 1:length(output)
    if max(max(output(i).all_ts)) >= 4 && nansum(output(i).snowp > 50) > 14 %needs at least 4 pixels and 2 weeks of snow cover
    dates = output(i).days;
    cloud = output(i).cloudp;
    maxpixels = max(sum(output(i).all_ts,2));
    %snow fraction
    ts = output(i).snowp;
    snowfrac = filter_ts(ts,cloud,dates,75);

    %water fraction
    ts = output(i).waterp;
    waterfrac = filter_ts(ts,cloud,dates,75);

    %snow albedo
    ts = output(i).snowa_ts;
    vc = output(i).snow_ts;
    snowalb = filter_ts_albedo(ts,cloud,vc,maxpixels,dates,75,10);

    %water albedo
    ts = output(i).watera_ts;
    vc = output(i).water_ts;
    wateralb = filter_ts_albedo(ts,cloud,vc,maxpixels,dates,75,10);
    
    fdata(i).snowfrac = snowfrac;
    fdata(i).waterfrac = waterfrac;
    fdata(i).snowalb = snowalb;
    fdata(i).wateralb = wateralb;
    else
        fdata(i).snowfrac = NaN;
        fdata(i).waterfrac = NaN;
        fdata(i).snowalb = NaN;
        fdata(i).wateralb = NaN;
    end
end

dates = datetime(2001,1,1):datetime(2022,12,31);
%STEP 3: Convert to monthly time series
for i =1:length(fdata)
    if length(fdata(i).snowfrac) > 2
    tsin = fdata(i).snowfrac;
    [btsmean,btsmedian] = get_monthly_values(tsin,dates);
    fdata(i).snowfracbw_mean = btsmean;
    fdata(i).snowfracbw_median = btsmedian;
    
    tsin = fdata(i).waterfrac;
    [btsmean,btsmedian] = get_monthly_values(tsin,dates);
    fdata(i).waterfracbw_mean = btsmean;
    fdata(i).waterfracbw_median = btsmedian;

    tsin = fdata(i).snowalb;
    [btsmean,btsmedian] = get_monthly_values(tsin,dates);
    fdata(i).snowalbbw_mean = btsmean;
    fdata(i).snowalbbw_median = btsmedian;

    tsin = fdata(i).wateralb;
    [btsmean,btsmedian] = get_monthly_values(tsin,dates);
    fdata(i).wateralbbw_mean = btsmean;
    fdata(i).wateralbbw_median = btsmedian;
    end
end

%STEP 4: Calculate average snowy days and snow-off date for lake and land
for i = 1:length(fdata)
if length(fdata(i).snowfrac) > 2
tsin = fdata(i).snowfrac;
y = 2002:2022;
for n = 1:length(y)
    years = year(dates);
    ind = find(years == y(n));
    ts = tsin(ind);
    %number of days with snow > 25%
    snow25(n,1) = sum(ts > 25);
    %average first day with snow below 25 (and stays below 25 for the next
    %20 days)
    snowf(n,1) = get_first_snowfree_day(ts,25);
end
fdata(i).snow = snow25;
fdata(i).snowf = snowf;
fdata(i).meansnow = nanmean(snow25);
fdata(i).meansnowf = nanmean(snowf);
else
    fdata(i).meansnow = 0;
    fdata(i).meansnowf = 0;
end


if p == 1; lakefdata = fdata; end
if p == 2; landfdata = fdata; end
end
end

   
%% Convert to radiative forcing

%Start for lakes
fdata = lakefdata;
for i = 1:length(fdata)
        %VARIABLES: snow fraction, snow albedo, water albedo
        snowfrac = 0.01*fdata(i).snowfracbw_mean;
        snowalb = fdata(i).snowalbbw_mean;
        snowfreealb = fdata(i).wateralbbw_mean;
        
        if isempty(snowfrac) == 0
        snowfrac(snowfrac < 0.001) = 0;
        %loop through years
        for j = 1:22
            %get water albedo       
            s = snowfrac(j,:);
            w = snowfreealb(j,:);
            ind = find(s < 0.1); %first month where snow < 10%
            if isempty(ind) == 0
                 sfa(j,1) = nanmean(w(ind(1))); 
            else
                 sfa(j,1) = NaN;
            end
        end

        sfa = nanmean(sfa);
        %correct for if no month < 10%, choose july august mean
        if isnan(sfa)
            w = snowfreealb(:,[7,8]);
            sfa = nanmean(nanmean(w));
        end


        %get era5 kernel values
        ker = 100*era5_ker(i,:);
        %calculate CRE (monthly)
        for j = 1:22
            for n = 1:12
                cre(j,n) = snowfrac(j,n).*(snowalb(j,n)-sfa).*ker(n);
            end
        end
        out(i).lakescre = cre;
        out(i).lakessnowalb = snowalb;
        out(i).lakessnowfrac = snowfrac;
        out(i).lakessnowfreealb = sfa;
        out(i).lakessnow = fdata(i).snow;
        out(i).lakessnowf = fdata(i).snowf;
        else
        out(i).lakescre = 0;
        out(i).lakessnowalb = 0;
        out(i).lakessnowfrac = 0;
        out(i).lakessnowfreealb = 0;
        out(i).lakessnow = 0;
        out(i).lakessnowf = 0;
        end
end


% Land 
fdata = landfdata;
for i = 1:length(fdata)
   %VARIABLES: snow fraction, snow albedo, water albedo
        snowfrac = 0.01*fdata(i).snowfracbw_mean;
        snowalb = fdata(i).snowalbbw_mean;
        snowfreealb = fdata(i).wateralbbw_mean;
        
        if isempty(snowfrac) == 0
        snowfrac(snowfrac < 0.001) = 0;

        %loop through years
        for j = 1:22
            %get water albedo       
            s = snowfrac(j,:);
            w = snowfreealb(j,:);
            ind = find(s < 0.1); %first month where snow < 10%
            if isempty(ind) == 0
                 sfa(j,1) = nanmean(w(ind(1))); 
            else
                 sfa(j,1) = NaN;
            end
        end

        sfa = nanmean(sfa);

        %correct for if no month < 10%, choose july and august mean
        if isnan(sfa)
            w = snowfreealb(:,[7,8]);
            sfa = nanmean(nanmean(w));
        end
            

        %get era5 kernel values
        ker = 100*era5_ker(i,:);
        %calculate CRE (monthly)
        for j = 1:22
            for n = 1:12
                cre(j,n) = snowfrac(j,n).*(snowalb(j,n)-sfa).*ker(n);
            end
        end
        out(i).landcre = cre;
        out(i).landsnowalb = snowalb;
        out(i).landsnowfrac = snowfrac;
        out(i).landsnowfreealb = sfa;
        out(i).landsnow = fdata(i).snow;
        out(i).landsnowf = fdata(i).snowf;
        else
        out(i).landcre = 0;
        out(i).landsnowalb = 0;
        out(i).landsnowfrac = 0;
        out(i).landsnowfreealb = 0;
        out(i).landsnow = 0;
        out(i).landsnowf = 0;
        end    

end


cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/0_results');
save(outname,'out','lakefdata','landfdata');

