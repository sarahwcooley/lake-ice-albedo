% COMPLETE ANALYSIS 

%% Step 1: read in files, calculate variables, process
clear all
for n = 1:5
    clear out lakefdata landfdata o
    cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/0_results')
    if n == 1; load('grid1_processed_monthly_jun26.mat'); end
    if n == 2; load('grid2_processed_monthly_jun26.mat'); end
    if n == 3; load('grid3_processed_monthly_jun26.mat'); end
    if n == 4; load('grid4_processed_monthly_jun26.mat'); end
    if n == 5; load('gridland_processed_monthly_jun26.mat'); end

    cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/shapefiles/0_grids_dec17');
    if n == 1; roi = shaperead('grids_lakes_1_dec17.shp'); end
    if n == 2; roi = shaperead('grids_lakes_2_dec17.shp'); end
    if n == 3; roi = shaperead('grids_lakes_3_dec17.shp'); end
    if n == 4; roi = shaperead('grids_lakes_4_dec17.shp'); end
    if n == 5; roi = shaperead('grids_land_all_dec17_v3.shp'); end

    for i = 1:length(out)
        
        %1. Spatial Variables
        o(i).lakefrac = roi(i).lake_fract/(1-(0.01*roi(i).glacier_fr)); %correcting lake frac for glacier frac
        o(i).glacierfrac = roi(i).glacier_fr;
        o(i).lat = roi(i).lat;
        o(i).area = (1-0.01*roi(i).glacier_fr).*roi(i).area; %correcting area for glacier area

        %2. Snow Off Variables
        landsnowf = out(i).landsnowf;
        lakesnowf = out(i).lakessnowf;
        o(i).landsnowoff = calculate_statistics(landsnowf);
        o(i).lakesnowoff = calculate_statistics(lakesnowf);
        if nansum(nansum(lakesnowf)) ~= 0 && nansum(nansum(landsnowf)) ~= 0
            o(i).snowoffdif = calculate_statistics(lakesnowf-landsnowf);
        else
            o(i).snowoffdif = nan(1,5);
        end
        o(i).monthlylandsnow = calculate_statistics_monthly(out(i).landsnowfrac);
        o(i).monthlylakesnow = calculate_statistics_monthly(out(i).lakessnowfrac);
        
        %3. Snow Albedo Variables
        %land
        if nansum(nansum(out(i).landsnowfrac)) ~= 0
            sf  = out(i).landsnowfrac;
            snowalb = out(i).landsnowalb;
            snowfreealb = out(i).landsnowfreealb;
            snowalb(sf < 0.25) = NaN; %only care about snowalb values where sf > 0.25
            snowalb = nanmean(snowalb(:,2:5),2);
            o(i).landsnowalb = calculate_statistics(snowalb);
            o(i).landsnowfreealb = snowfreealb;
            landsnow = snowalb - snowfreealb;
        else
            o(i).landsnowalb = nan(1,5);
            o(i).landsnowfreealb = NaN;
        end

        %lakes
        if nansum(nansum(out(i).lakessnowfrac)) ~= 0
            sf  = out(i).lakessnowfrac;
            snowalb = out(i).lakessnowalb;
            snowfreealb = out(i).lakessnowfreealb;
            snowalb(sf < 0.25) = NaN; %only care about snowalb values where sf > 0.25
            snowalb = nanmean(snowalb(:,2:5),2);
            o(i).lakesnowalb = calculate_statistics(snowalb);
            o(i).lakesnowfreealb = snowfreealb;
            lakesnow = snowalb - snowfreealb;
            o(i).albdif = calculate_statistics(lakesnow-landsnow);
        else
            o(i).lakesnowalb = nan(1,5);
            o(i).lakesnowfreealb = NaN;
            o(i).albdif = nan(1,5);
        end
      
        %4. CRE Variables
        %land 
        landcre = out(i).landcre;
        o(i).landmonthlycre = calculate_statistics_monthly(landcre);
        cre = nansum(landcre,2)./12;
        o(i).landtotalcre = calculate_statistics(cre);

        %lakes
        lakecre = out(i).lakescre;
        o(i).lakemonthlycre = calculate_statistics_monthly(lakecre);
        cre = nansum(lakecre,2)./12;
        o(i).laketotalcre = calculate_statistics(cre);

        %total CRE
        if nansum(nansum(lakecre)) ~= 0
            lf = 0.01*o(i).lakefrac;
            totalcre = lf*lakecre + (1-lf)*landcre;
        else
            totalcre = landcre;
        end
        o(i).totalmonthlycre = calculate_statistics_monthly(totalcre);
        cre = sum(totalcre,2)./12;
        o(i).totalcre = calculate_statistics(cre);

        %5. Lake Contribution
        %lake contribution
        if nansum(nansum(lakecre)) ~= 0; 
            lf = 0.01*o(i).lakefrac;
            lc = 100*((lf*lakecre)./totalcre);
            o(i).monthlylakecont = calculate_statistics_monthly(lc);
            var = 100*((lf*nansum(lakecre,2))./nansum(totalcre,2));
            o(i).lakecont = calculate_statistics(var);
        else
            if nansum(nansum(landcre)) ~= 0;
                o(i).monthlylakecont = zeros(12,5);
                o(i).lakecont = zeros(1,5);
                lc = zeros(22,12);
            else
                o(i).monthlylakecont = nan(12,5);
                o(i).lakecont = nan(1,5);
                lc = nan(22,12);
            end
        end
          
        %lake fraction increase
        lfi = 100*lc./o(i).lakefrac;
        o(i).monthlylakefracinc = calculate_statistics_monthly(lc);
        var = 100*(100*((lf*nansum(lakecre,2))./nansum(totalcre,2)))./o(i).lakefrac;
        o(i).lakefracinc = calculate_statistics(var);

        %lake increase
        if nansum(nansum(lakecre)) ~= 0
            li = 100*(lakecre - landcre)./landcre;
            o(i).monthlylakeinc = calculate_statistics_monthly(li);
            var = 100*(nansum(lakecre,2) - nansum(landcre,2))./nansum(landcre,2);
            o(i).lakeinc = calculate_statistics(var);
        else
            if nansum(nansum(landcre)) ~= 0;
                o(i).monthlylakeinc = zeros(12,5);
                o(i).lakeinc = zeros(1,5);
            else
                o(i).monthlylakeinc = nan(12,5);
                o(i).lakeinc = nan(1,5);
            end
        end
              
    end
    if n < 5; roi = rmfield(roi,'idall'); end
    if n == 1; outvar = o; outout = out; roio = roi; else; outvar = cat(2,outvar,o); outout = cat(2,outout,out); roio = cat(1,roio,roi); end
end

% Add Biome data
cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/0_results');
load('ecodata_may9.mat')
o = outvar;
for i = 1:length(o)
    o(i).id = i;
    outout(i).id = i;
    o(i).biome = ecodata(i).biome_1;
    outout(i).biome = ecodata(i).biome_1;
end

save('complete_results_jul3.mat','out','o');
op = o;

%% STEP 2: GLOBAL and BIOME STATISTICS
%globally, calculate statistics in key variables
%write function for calculating land weighting globally and by biome

%KEY VARIABLES: 
% Snow-off day lakes and land, snow-off day difference
% Albedo difference lakes and land
% Lake CRE, Land CRE, Total CRE (monthly and annual)
% Lake Contribution, Lake Cont Inc (two ways) monthly and annual)

%For each variable, calculate (1) mean, (2) median, (3) std, (4) IQR
%Statistics calculated 2 ways - (1) one where we take the mean of all stats
%variables, other where we calculate stats variables on the mean

biome = [op.biome]';
bb = zeros(size(biome));
%reordering biomes from largest to smallest
bb(biome == 6) = 1;
bb(biome == 4) = 2;
bb(biome == 11) = 3;
bb(biome == 8) = 4;
bb(biome == 13) = 5;
bb(biome == 10) = 6;
bb(biome == 5) = 7;

for n = 1:8
    if n == 1
       o = op; 
    else
       o = op(bb == n-1);
    end
    area = [o.area]';
    clear var var2

% % 1. Snow-Off Day %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lf = 0.01*[o.lakefrac]';
%land snow off timing
for i = 1:length(o); var(i,:) = o(i).landsnowoff; end
stats1.landsnowoff = get_global_weighted_lc_stats(var,(1-lf).*area);
stats2.landsnowoff = get_global_weighted_stats(var(:,1),(1-lf).*area);

%lake snow off timing
for i = 1:length(o); var(i,:) = o(i).lakesnowoff; end
stats1.lakesnowoff = get_global_weighted_lc_stats(var,lf.*area);
stats2.lakesnowoff = get_global_weighted_stats(var(:,1),lf.*area);

%snow off dif
for i = 1:length(o); var(i,:) = o(i).snowoffdif; end
stats1.snowoffdif = get_global_weighted_lc_stats(var,area);
stats2.snowoffdif = get_global_weighted_stats(var(:,1),area);

% % 2. Albedo Difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%land
for i = 1:length(o); sfa = o(i).landsnowfreealb; var(i,:) = o(i).landsnowalb - [sfa, sfa, 0, sfa, sfa]; end
stats1.landalbdif = get_global_weighted_lc_stats(var,(1-lf).*area);
stats2.landalbdif = get_global_weighted_stats(var(:,1),(1-lf).*area);

%lakes
for i = 1:length(o); sfa = o(i).lakesnowfreealb; var(i,:) = o(i).lakesnowalb - [sfa, sfa, 0, sfa, sfa]; end
stats1.lakealbdif = get_global_weighted_lc_stats(var,lf.*area);
stats2.lakealbdif = get_global_weighted_stats(var(:,1),lf.*area);

%difference
for i = 1:length(o); var(i,:) = o(i).albdif; end
stats1.albdif = get_global_weighted_lc_stats(var,area);
stats2.albdif = get_global_weighted_stats(var(:,1),area);


%3. MONTHLY CRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%land
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).landmonthlycre(p,:); end
    stats1.landmonthlycre(p,:)  = get_global_weighted_lc_stats(var,(1-lf).*area);
    stats2.landmonthlycre(p,:) = get_global_weighted_stats(var(:,1),(1-lf).*area);
end

%lakes
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).lakemonthlycre(p,:); end
    stats1.lakemonthlycre(p,:) = get_global_weighted_lc_stats(var,lf.*area);
    stats2.lakemonthlycre(p,:) = get_global_weighted_stats(var(:,1),lf.*area);
end

%total
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).totalmonthlycre(p,:); end
    stats1.totalmonthlycre(p,:) = get_global_weighted_lc_stats(var,area);
    stats2.totalmonthlycre(p,:) = get_global_weighted_stats(var(:,1),area);
end

%4. ANNUAL CRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%land
for i = 1:length(o); var(i,:) = o(i).landtotalcre; end
stats1.landtotalcre = get_global_weighted_lc_stats(var,(1-lf).*area);
stats2.landtotalcre = get_global_weighted_stats(var(:,1),(1-lf).*area);

%lake
for i = 1:length(o); var(i,:) = o(i).laketotalcre; end
stats1.laketotalcre = get_global_weighted_lc_stats(var,lf.*area);
stats2.laketotalcre = get_global_weighted_stats(var(:,1),lf.*area);

%total
for i = 1:length(o); var(i,:) = o(i).totalcre; end
stats1.totalcre = get_global_weighted_lc_stats(var,area);
stats2.totalcre = get_global_weighted_stats(var(:,1),area);

% % 5. Lake Contribution Monthly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lake contribution
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylakecont(p,:); var2(i,:) = o(i).totalmonthlycre(p,:); end
    stats1.monthlylakecont(p,:) = get_global_weighted_lc_stats(var,area.*abs(var2));
    stats2.monthlylakecont(p,:) = get_global_weighted_stats(var(:,1),area.*abs(var2(:,1)));
end

%lake frac inc
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylakefracinc(p,:); var2(i,:) = o(i).totalmonthlycre(p,:); end
    stats1.monthlylakefracinc(p,:) = get_global_weighted_lc_stats(var,area);
    stats2.monthlylakefracinc(p,:) = get_global_weighted_stats(var(:,1),area);
end

%lake inc
for p = 1:12
    for i = 1:length(o)
        var(i,:) = o(i).monthlylakeinc(p,:);
        var2(i,:) = o(i).totalmonthlycre(p,:);
        II = isoutlier(var(i,:),'median',10);
        var(II) = NaN;
    end
    stats1.monthlylakeinc(p,:) = get_global_weighted_lc_stats(var,area);
    stats2.monthlylakeinc(p,:) = get_global_weighted_stats(var(:,1),area);
end

% % 6. Lake Contribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%total lake cont
for i = 1:length(o); var(i,:) = o(i).lakecont; var2(i,:) = o(i).totalcre; end
stats1.lakecont = get_global_weighted_lc_stats(var,area.*abs(var2));
stats2.lakecont = get_global_weighted_stats(var(:,1),area.*abs(var2(:,1)));

%total lake frac inc
for i = 1:length(o); var(i,:) = o(i).lakefracinc; var2(i,:) = o(i).totalcre; end
stats1.lakefracinc = get_global_weighted_lc_stats(var,area);
stats2.lakefracinc = get_global_weighted_stats(var(:,1),area);

%total lake inc
for i = 1:length(o); var(i,:) = o(i).lakeinc; var2(i,:) = o(i).totalcre; end
II = var(:,1) ~= 0;
stats1.lakeinc = get_global_weighted_lc_stats(var(II,:),area(II));
stats2.lakeinc = get_global_weighted_stats(var(II,1),area(II));

% % 7. Lake Contribution other calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lake contribution monthly
clear landcre lakecre lakefrac
for p = 1:12
    for i = 1:length(o); landcre(i,:) = o(i).landmonthlycre(p,:); lakecre(i,:) = o(i).lakemonthlycre(p,:); lakefrac(i,:) = o(i).lakefrac; end
    lcstats.monthlylakecont(p,:) = get_global_weighted_lakecont_stats(landcre,lakecre,lakefrac,area);
    [lcstats.monthlylakefracinc(p,:),lcstats.monthlylakeinc(p,:)] = get_global_weighted_lakeinc_stats(landcre,lakecre,lakefrac,area);
end

%lake contribution annual
clear landcre lakecre
for i = 1:length(o); landcre(i,:) = o(i).landtotalcre; lakecre(i,:) = o(i).laketotalcre; lakefrac(i,:) = o(i).lakefrac; end
lcstats.lakecont = get_global_weighted_lakecont_stats(landcre,lakecre,lakefrac,area);
[lcstats.lakefracinc,lcstats.lakeinc] = get_global_weighted_lakeinc_stats(landcre,lakecre,lakefrac,area);

% % Convert into tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n == 1 %global tables
    [t1,t2,t3,t4] = get_global_stats_tables(stats1,stats2,lcstats);
else %biome tables
    [b1(n-1,:),b12(n-1,:),b13(n-1,:),b2(n-1,:)] = get_biome_stats_tables(stats1,stats2,lcstats,lakefrac,area);
end

% % store values for figures
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylandsnow(p,:); end
    stats1.landmonthlysnow(p,:)  = get_global_weighted_lc_stats(var,area);
    stats2.landmonthlysnow(p,:) = get_global_weighted_stats(var(:,1),area);
end

for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylakesnow(p,:); end
    stats1.lakemonthlysnow(p,:)  = get_global_weighted_lc_stats(var,area);
    stats2.lakemonthlysnow(p,:) = get_global_weighted_stats(var(:,1),area);
end

s1(n) = stats1;
s2(n) = stats2;
lcs(n) = lcstats;


end

%calculate proportion in each biome
areaper = 100*b2(:,2)./sum(b2(:,2));
landarea = b2(:,2).*(1-0.01.*b2(:,1));
landper = 100*b1(:,7).*landarea./sum(b1(:,7).*landarea);
lakearea = b2(:,2).*0.01.*b2(:,1);
lakeareaper = 100*lakearea./sum(lakearea);
lakeper = 100*b1(:,8).*lakearea./sum(b1(:,8).*lakearea);
totalper = 100*b1(:,9).*b2(:,2)./(sum(b1(:,9).*b2(:,2)));
pertable = cat(2,areaper,totalper,lakeareaper,lakeper);

%% Radiative Forcing (Sensitivity) Analysis
cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/ancillary_data/era5_temp');
load('temp_data_jun20.mat')

%Run sensitivity analyses (1) globally, (2) by biome, (3) individually by
%grid cell

%1. Global Lake and Global Land Mean CRE
for w = 1:8
    if w == 1
       o = op; 
       td = tempdata;
       ot = outout;
    else
       o = op(bb == w-1);
       td = tempdata(bb == w-1);
       ot = outout(bb == w-1);
    end
    area = [o.area]';
    lf = 0.01*[o.lakefrac]';

    lakecre = nan(22,length(o));
    landcre = nan(22,length(o));
    temp = nan(22,length(o));
    
    %get temperature data
    for n = 1:12
        for i = 1:length(td); temp(:,i) = td(i).monthly(:,n); end
        for j = 1:22; ts = temp(j,:); II = ~isnan(ts); 
            tempmeanlake(j,1) = mean(ts(II),Weight=lf(II).*area(II)); 
            tempmeanland(j,1) = mean(ts(II),Weight=(1-lf(II)).*area(II)); 
        end

    %get CRE data
    %clear lakecre landcre
    for i = 1:length(ot)
        if length(ot(i).lakescre) > 1; lakecre(:,i) = ot(i).lakescre(:,n); end
        if length(ot(i).landcre) > 1; landcre(:,i) = ot(i).landcre(:,n); end
    end
    for j = 1:22
        ts = lakecre(j,:)';
        II = ~isnan(ts);
        lake(j,1) = mean(ts(II),Weight=lf(II).*area(II));
        ts = landcre(j,:)';
        II = ~isnan(ts);
        land(j,1) = mean(ts(II),Weight=(1-lf(II)).*area(II));
    end

    %regress lakes and land
    x = tempmeanlake;
    y = lake;
    II = ~isnan(y);
    [r,p] = corrcoef(x(II),y(II));
    P1=polyfit(x(II),y(II),1); 
    mdl = fitlm(x(II),y(II),'linear');
    tempresults(n).lake = [r(2),p(2),P1(1)];
    tempresults(n).lakemdl = mdl.Coefficients;

    x = tempmeanland;
    y = land;
    II = ~isnan(y);
    [r,p] = corrcoef(x(II),y(II));
    P2=polyfit(x(II),y(II),1);
    mdl = fitlm(x(II),y(II),'linear');
    tempresults(n).land = [r(2),p(2),P2(1)];
    tempresults(n).landmdl = mdl.Coefficients;
    end
    
    if w == 1
        [global_sens_table] = get_sensitivity_table(tempresults);
    else
        [c] = get_sensitivity_table(tempresults);
        biomeresults(w-1).results = c;
    end
end

%get biome results table
for w = 1:7
    cc = 1;
    for n = 3:6
        for p = 3:6; biome_sens_table(w,cc) = biomeresults(w).results(p,n); cc = cc+1; end
    end
end

%3. Calculate individually by grid cell

for i = 1:length(op)
    temp = tempdata(i).monthly;
    lakecre = outout(i).lakescre;
    landcre = outout(i).landcre;
    if length(lakecre) > 1
        op(i).lakesens = calculate_sensitivity(temp,lakecre);
        [m,ind] = nanmax(op(i).lakesens);
        op(i).maxlakesens = [m,ind];
    else
        op(i).lakesens = nan(12,1);
        op(i).maxlakesens = nan(1,2);
    end
    if length(landcre) > 1
        op(i).landsens = calculate_sensitivity(temp,landcre);
        [m,ind] = nanmax(op(i).landsens);
        op(i).maxlandsens = [m,ind];
    else
        op(i).landsens = nan(12,1);
        op(i).maxlandsens = nan(1,2);
    end
end

%% Test for significance of difference over 22 years

%variables = lakecre, landcre, snowofflakes, snowoffland, albedolakes,
%albedoland
for w = 1:8
    if w == 1
       o = op; 
       td = tempdata;
       ot = outout;
    else
       o = op(bb == w-1);
       td = tempdata(bb == w-1);
       ot = outout(bb == w-1);
    end
    area = [o.area]';
    lf = 0.01*[o.lakefrac]';

    lakecre = nan(22,length(o));
    landcre = nan(22,length(o));
    lakessf = nan(21,length(o));
    landsf = nan(21,length(o));
    lakealb = nan(22,length(o));
    landalb = nan(22,length(o));

  %get data

    for i = 1:length(ot)
        if length(ot(i).lakescre) > 1; lakecre(:,i) = nanmean(ot(i).lakescre,2); end
        if length(ot(i).landcre) > 1; landcre(:,i) = nanmean(ot(i).landcre,2); end
        if length(ot(i).lakessnowf) > 1; lakessf(:,i) = ot(i).lakessnowf; end
        if length(ot(i).landsnowf) > 1; landsf(:,i) = ot(i).landsnowf; end

        if length(ot(i).lakessnowalb) > 1
            for n = 1:22
            sf  = ot(i).lakessnowfrac(n,:);
            snowalb = ot(i).lakessnowalb(n,:);
            snowfreealb = ot(i).lakessnowfreealb;
            snowalb(sf < 0.25) = NaN; %only care about snowalb values where sf > 0.25
            snowa(n,1) = nanmean(snowalb(:,2:5));
            end
            lakealb(:,i) = snowa;
        end
        if length(ot(i).landsnowalb) > 1
            for n = 1:22
            sf  = ot(i).landsnowfrac(n,:);
            snowalb = ot(i).landsnowalb(n,:);
            snowfreealb = ot(i).landsnowfreealb;
            snowalb(sf < 0.25) = NaN; %only care about snowalb values where sf > 0.25
            snowa(n,1) = nanmean(snowalb(:,2:5));
            end
            landalb(:,i) = snowa;
        end
    end


    %testing lakes vs land distribution
    for jj = 1:3
        if jj == 1; varlake = lakecre; varland = landcre; end
        if jj == 2; varlake = lakessf; varland = landsf; end
        if jj == 3; varlake = lakealb; varland = landalb; end

    clear lake land
    if jj == 2; jt = 21; else; jt = 22; end
    for j = 1:jt
        ts = varlake(j,:)';
        II = ~isnan(ts);
        lake(j,1) = mean(ts(II),Weight=lf(II).*area(II));
        ts = varland(j,:)';
        II = ~isnan(ts);
        land(j,1) = mean(ts(II),Weight=(1-lf(II)).*area(II));
    end
    [h,p,ci] = ttest2(lake,land);
    if jj ==1; tres(w).cre = [h,p]; tres(w).creci = ci; end
    if jj == 2; tres(w).snow = [h,p]; tres(w).snowci = ci; end
    if jj == 3; tres(w).alb = [h,p]; tres(w).albci = ci; end
    end
end


%% Write everything to shapefile

%variables: lake fraction, landsnowff, lakesnowoff, snowoffdif, albdif
%land total cre, lake total cre, lake cont, lakefracinc, lakeinc
%max lake sens, max land sens, lake sens dif w/m2, lake sens dif %

for i = 1:length(op)
    os(i).Geometry = 'Polygon';
    os(i).id = i;
    os(i).X = roio(i).X;
    os(i).Y = roio(i).Y;
    os(i).lakefrac = op(i).lakefrac;
    os(i).landsnowoff = op(i).landsnowoff(1);
    os(i).lakesnowoff = op(i).lakesnowoff(1);
    os(i).snowoffdif = op(i).snowoffdif(1);
    os(i).albdif = op(i).albdif(1);
    os(i).landcre = op(i).landtotalcre(1);
    os(i).lakecre = op(i).laketotalcre(1);
    os(i).lakecont = op(i).lakecont(1);
    os(i).lakefracinc = op(i).lakefracinc(1);
    os(i).lakeinc = op(i).lakeinc(1);
    os(i).landsens = op(i).maxlandsens(1);
    os(i).lakesens = op(i).maxlakesens(1);
    os(i).sensdif = op(i).maxlakesens(1) - op(i).maxlandsens(1);
    os(i).sensdifp = 100*(op(i).maxlakesens(1) - op(i).maxlandsens(1))./op(i).maxlandsens(1);
end
cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/0_results/shapefiles')
shapewrite(os,'complete_results_jul7_lakeweighted.shp')


%% FIGURES 

%1. Panel plot showing (1) monthly snow fraction, (2) albedo contrast
%(3) total CRE, (4) monthly CRE, (5) monthly lake contribution (maybe add line showing lake
%fraction?) - GLOBAL AND FOR EACH BIOME (combine into one plot?)

%2. Sensitivity figures - monthly sensitivity globally and by biome? 


%Figure 1: Snow Fraction, albedo contrast, monthly CRE
map = {'#6CBA7D','#3C93C2'};
cmap = validatecolor(map,'multiple');
titles = {'Northern Hemisphere','Boreal Forest','Temperate Broadleaf Forest','Tundra','Temperate Grassland','Desert','Montane Grassland','Temperate Conifer Forest'};
for i = 1:8

    %Monthly CRE
    figure(1)
    subplot(4,2,i)
    [mts1,x,inBetween] = variables_for_fillplot_v2(s1(i).landmonthlycre);
    hold off
    h = fill(x,inBetween,cmap(1,:));
    h.FaceAlpha = 0.25;
    h.EdgeAlpha = 0;
    hold on
    [mts2,x,inBetween] = variables_for_fillplot_v2(s1(i).lakemonthlycre);
    h = fill(x,inBetween,cmap(2,:));
    h.FaceAlpha = 0.25;
    h.EdgeAlpha = 0;
    h1 = plot(1:12,mts1,'LineWidth',2,'Color',cmap(1,:));
    h2 = plot(1:12,mts2,'LineWidth',2,'Color',cmap(2,:));
    ylabel('CrRE (W/m^2)')
    xticks([1:12]);
    xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'})
    xlim([0.5 12.5])
    set(gca,'FontSize',14)
    if i == 1; legend([h1,h2],{'Land','Lakes'}); end
        title(titles{i});
    xlim([1 12])

end
    %

% FIGURE 2: CRE (lake/land), lake frac vs. lake cont, snow-off timing,
% albedo contrast
 clear var1 var1l var1h var2

for i = 1:8
    var1(i,1) = s1(i).landtotalcre(1);
    var1h(i,1) = s1(i).landtotalcre(1)+s1(i).landtotalcre(3);
    var1l(i,1) = s1(i).landtotalcre(1)-s1(i).landtotalcre(3);
    var1(i,2) = s1(i).laketotalcre(1);
    var1h(i,2) = s1(i).laketotalcre(1)+s1(i).laketotalcre(3);
    var1l(i,2) = s1(i).laketotalcre(1)-s1(i).laketotalcre(3);
end

    figure(2)
    t = tiledlayout(2,2,'TileSpacing','tight');
    nexttile
    hold off
    colororder(cmap);
    bar(1:8,var1,1)
    hold on
    x = [0.85 1.85 2.85 3.85 4.85 5.85 6.85 7.85];
    errorbar(x,var1(:,1),var1(:,1)-var1l(:,1),var1h(:,1)-var1(:,1),'LineStyle','none','LineWidth',1,'Color','k');
    x = [1.15, 2.15, 3.15, 4.15, 5.15, 6.15, 7.15, 8.15];
    errorbar(x,var1(:,2),var1(:,2)-var1l(:,2),var1h(:,2)-var1(:,2),'LineStyle','none','LineWidth',1,'Color','k');
    plot([1.5 1.5],[-26 0],'k-')
    xticklabels(titles)
    ylabel('Mean Annual CrRE_t (W/m^2)')
    legend({'Land','Lakes'})
    set(gca,'FontSize',14)
    tcre1 = var1;
    tcre1h = var1h;
    tcre1l = var1l;
    ylim([-26 0])


 clear var1 var1l var1h var2
 for i = 1:8
     var1(i,1) = lcs(i).lakecont(1);
     var1h(i,1) = lcs(i).lakecont(1)+s1(i).lakecont(3);
     var1l(i,1) = lcs(i).lakecont(1)-s1(i).lakecont(3);
     if i == 1; vv(i,1) = 3.70; else
        vv(i,1) = b2(i-1,1); end
 end
    nexttile
    hold off
    cm = [0.45 0.45 0.45;0.75 0.75 0.75];
    b = bar(1:8,[vv,var1],1,'FaceColor','flat');
    for k = 1:2
        b(k).CData = cm(k,:);
    end
    hold on
    x = [1.15, 2.15, 3.15, 4.15, 5.15, 6.15, 7.15, 8.15];
    errorbar(x,var1(:,1),var1(:,1)-var1l(:,1),var1h(:,1)-var1(:,1),'LineStyle','none','LineWidth',1,'Color','k');
    plot([1.5 1.5],[0 11],'k-')
    xticklabels(titles)
    ylabel('Lake %')
    legend({'Lake Fraction','Lake CrRE_t Contribution'})
    set(gca,'FontSize',14)
    ylim([0 11])
    lf1 = var1;
    lf1h = var1h;
    lf1l = var1l;

for i = 1:8
     var1(i,1) = s1(i).landsnowoff(1);
     var1h(i,1) = s1(i).landsnowoff(1)+s1(i).landsnowoff(3);
     var1l(i,1) = s1(i).landsnowoff(1)-s1(i).landsnowoff(3);
     var1(i,2) = s1(i).lakesnowoff(1);
     var1h(i,2) = s1(i).lakesnowoff(1)+s1(i).lakesnowoff(3);
     var1l(i,2) = s1(i).lakesnowoff(1)-s1(i).lakesnowoff(3);
 end 
   
    figure(2)
    nexttile
    hold off
    colororder(cmap);
    bar(1:8,var1,1)
    hold on
    x = [0.85 1.85 2.85 3.85 4.85 5.85 6.85 7.85];
    errorbar(x,var1(:,1),var1(:,1)-var1l(:,1),var1h(:,1)-var1(:,1),'LineStyle','none','LineWidth',1,'Color','k');
    x = [1.15, 2.15, 3.15, 4.15, 5.15, 6.15, 7.15, 8.15];
    errorbar(x,var1(:,2),var1(:,2)-var1l(:,2),var1h(:,2)-var1(:,2),'LineStyle','none','LineWidth',1,'Color','k');
    plot([1.5 1.5],[1 182],'k-')
    xticklabels(titles)
    ylabel('Snow-Off Timing')
    legend({'Land','Lakes'})
    set(gca,'FontSize',14)
    ylim([1 182])
    yticks([1 32 60 91 121 152 182])
    yticklabels({'Jan 1','Feb 1','Mar 1','Apr 1','May 1','Jun 1','Jul 1'})
    so1 = var1;
    so1h = var1h;
    so1l = var1l;

for i = 1:8
    var1(i,1) = s1(i).landalbdif(1);
    var1h(i,1) = s1(i).landalbdif(1)+s1(i).landalbdif(3);
    var1l(i,1) = s1(i).landalbdif(1)-s1(i).landalbdif(3);
    var1(i,2) = s1(i).lakealbdif(1);
    var1h(i,2) = s1(i).lakealbdif(1)+s1(i).lakealbdif(3);
    var1l(i,2) = s1(i).lakealbdif(1)-s1(i).lakealbdif(3);
end
   
    figure(2)
    nexttile
    hold off
    colororder(cmap);
    bar(1:8,var1,1)
    hold on
    x = [0.85 1.85 2.85 3.85 4.85 5.85 6.85 7.85];
    errorbar(x,var1(:,1),var1(:,1)-var1l(:,1),var1h(:,1)-var1(:,1),'LineStyle','none','LineWidth',1,'Color','k');
    x = [1.15, 2.15, 3.15, 4.15, 5.15, 6.15, 7.15, 8.15];
    errorbar(x,var1(:,2),var1(:,2)-var1l(:,2),var1h(:,2)-var1(:,2),'LineStyle','none','LineWidth',1,'Color','k');
    plot([1.5 1.5],[0 0.62],'k-')
    xticklabels(titles)
    ylabel('Seasonal Albedo Contrast')
    legend({'Land','Lakes'})
    set(gca,'FontSize',14)
    ylim([0 0.62])
    sa1 = var1;
    sa1h = var1h;
    sa1l = var1l;


% % FIGURE 3: Sensitivity By Month %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%March, April, May, June
monthtitles = {'February','March','April','May','June'};
clear var varh varl
figure(3)
t = tiledlayout(2,4,'TileSpacing','compact');
for i = 1:8
    for j = 1:5
        if i == 1
            tt = global_sens_table; 
        else
            tt = biomeresults(i-1).results;
        end
        var(j,1) = tt(3,j+1);
        var(j,2) = tt(4,j+1);
        varh(j,1) = tt(7,j+1);
        varh(j,2) = tt(8,j+1);
    end
    varl = varh;
    nexttile
   
    hold off
    colororder(cmap);
    bar(1:5,var,1)
    hold on
    x = [0.85 1.85 2.85 3.85 4.85]; % 5.85 6.85 7.85];
    errorbar(x,var(:,1),varl(:,1),varh(:,1),'LineStyle','none','LineWidth',1,'Color','k');
    x = [1.15, 2.15, 3.15, 4.15, 5.15]; % 6.15, 7.15, 8.15];
    errorbar(x,var(:,2),varl(:,2),varh(:,2),'LineStyle','none','LineWidth',1,'Color','k');
    %title(titles{i})
    if i == 1 || i == 5; 
    ylabel('CrRF_t (W/m^2/°C)')
    end
    if i ==1; legend({'Land','Lakes'}); end
    set(gca,'FontSize',14)
    if i > 4
    xticklabels(monthtitles)
    else
    xticklabels({});
    end

    ylim([0 7])

end


%Figure 4: Alternate plotting of all variables %%%%%%%%%%%%%%%%%%%%

figure(4)
map = {'#6CBA7D','#3C93C2'};
cmap = validatecolor(map,'multiple');
titles = {'Northern Hemisphere','Boreal Forest','Temperate Broadleaf Forest','Tundra','Temperate Grassland','Desert','Montane Grassland','Temperate Conifer Forest'};
figure(4)
t = tiledlayout(3,3,'TileSpacing','compact');   
for i = 1:8
    inds = [3 4 5 6 7 8 9];
    %Monthly CRE (plots 
    if i == 1; nexttile([1 2]); %subplot(9,6,[1:3,7:9]); 
    else; nexttile %subplot(9,6,[(13+(i-2)*6):(13+(i-2)*6)+2])
    end
    [mts1,x,inBetween] = variables_for_fillplot_v2(s1(i).landmonthlycre);
    hold off
    h = fill(x,inBetween,cmap(1,:));
    h.FaceAlpha = 0.25;
    h.EdgeAlpha = 0;
    hold on
    [mts2,x,inBetween] = variables_for_fillplot_v2(s1(i).lakemonthlycre);
    h = fill(x,inBetween,cmap(2,:));
    h.FaceAlpha = 0.25;
    h.EdgeAlpha = 0;
    h1 = plot(1:12,mts1,'LineWidth',2,'Color',cmap(1,:));
    h2 = plot(1:12,mts2,'LineWidth',2,'Color',cmap(2,:));
    if i == 1 || i == 3 || i == 6
    ylabel('CrRE_t (W/m^2)')
    end
    xticks([1:12]);
    xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'})
    xlim([0.5 12.5])
    ylim([-100 0])
    set(gca,'FontSize',14)
    if i == 1; legend([h1,h2],{'Land','Lakes'}); ylim([-100 0]); end
        %title(titles{i});
        xlim([1 12])
 


end

% Save figures
cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/figures/aug5')
%gcf = figure(1);
%saveas(gcf,'figure1_lw_v2.png');
%gcf = figure(2);
%saveas(gcf,'figure2_lw_v3.png');
%gcf = figure(3);
%saveas(gcf,'figure3_lw_v3_nt.png');
%gcf = figure(4);
%saveas(gcf,'figure1alt_lw_v2_3_nt.png');

%% MAP FIGURES
ss = os;
for i = 1:length(os);
    ss(i).Lat = os(i).Y;
    ss(i).Lon = os(i).X;
end

ss = rmfield(ss,'X');
ss = rmfield(ss,'Y');

%Snow Off Day Difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.snowoffdif]';
thresh = [nanmin(var) -25 -15 -10 -5 0 5 10 15 25 nanmax(var)];
cblabel = 'Snow-Off Timing Δ (days)';
name = 'map1';
cmap = brewermap(length(thresh)-1,'PRGn');
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%Albedo Difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.albdif]';
thresh = [nanmin(var) -0.1 -0.05 -0.025 0 0.025 0.05 0.1 nanmax(var)];
cblabel = 'Seasonal Albedo Contrast Δ';
name = 'map2';
cmap = brewermap(length(thresh)-1,'BrBG');
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%Land CRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.landcre]';
thresh = [nanmin(var) -25 -20 -15 -10 -6 -3 -1 nanmax(var)];
cblabel = 'Land CrRE_t (W/m^2)';
name = 'map3';
cmap = brewermap(length(thresh)-1,'OrRd'); cmap = flipud(cmap);
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%Lake CRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.lakecre]';
thresh = [nanmin(var) -25 -20 -15 -10 -6 -3 -1 nanmax(var)];
cblabel = 'Lake CrRE_t (W/m^2)';
name = 'map4';
cmap = brewermap(length(thresh)-1,'OrRd'); cmap = flipud(cmap);
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%Lake Cont %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.lakecont]';
thresh = [nanmin(var) 0.1 1 5 10 15 20 25 nanmax(var)];
cblabel = 'Lake Contribution to CrRE_t (%)';
name = 'map5';
cmap = brewermap(length(thresh)-1,'Blues');
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%Lake Inc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.lakeinc]';
var1 = [os.lakefracinc]';
var(isnan(var1)) = NaN; % correcting for NaN for tiles with no lakes
thresh = [nanmin(var) 0 10 25 50 75 100 nanmax(var)];
cblabel = 'Lake CrRE_t Increase (%)';
name = 'map6';
cmap = brewermap(length(thresh)-1,'Reds');
cmap(1,:) = [0.9 0.9 0.4];
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%Lake Fraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var = [os.lakefrac]';
thresh = [nanmin(var) 0.1 1 5 10 15 20 25 nanmax(var)];
cblabel = 'Lake Fraction (%)';
name = 'map7';
cmap = brewermap(length(thresh)-1,'Blues');
[h] = make_map_plot(ss,var,thresh,cblabel,name,cmap);

%% Results for ROIs
cd('/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/ancillary_data');
load('rois_jul9.mat');

%globally, calculate statistics in key variables
%write function for calculating land weighting globally and by biome

%KEY VARIABLES: 
% Snow-off day lakes and land, snow-off day difference
% Albedo difference lakes and land
% Lake CRE, Land CRE, Total CRE (monthly and annual)
% Lake Contribution, Lake Cont Inc (two ways) monthly and annual)

%For each variable, calculate (1) mean, (2) median, (3) std, (4) IQR
%Statistics calculated 2 ways - (1) one where we take the mean of all stats
%variables, other where we calculate stats variables on the mean
roi = [rois.roi]';

for n = 1:4
  
    o = op(roi == n);
    area = [o.area]';
    clear var var2

% % 1. Snow-Off Day %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lf = 0.01*[o.lakefrac]';
%land snow off timing
for i = 1:length(o); var(i,:) = o(i).landsnowoff; end
rstats1.landsnowoff = get_global_weighted_lc_stats(var,(1-lf).*area);
rstats2.landsnowoff = get_global_weighted_stats(var(:,1),(1-lf).*area);

%lake snow off timing
for i = 1:length(o); var(i,:) = o(i).lakesnowoff; end
rstats1.lakesnowoff = get_global_weighted_lc_stats(var,lf.*area);
rstats2.lakesnowoff = get_global_weighted_stats(var(:,1),lf.*area);

%snow off dif
for i = 1:length(o); var(i,:) = o(i).snowoffdif; end
rstats1.snowoffdif = get_global_weighted_lc_stats(var,area);
rstats2.snowoffdif = get_global_weighted_stats(var(:,1),area);

% % 2. Albedo Difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%land
for i = 1:length(o); sfa = o(i).landsnowfreealb; var(i,:) = o(i).landsnowalb - [sfa, sfa, 0, sfa, sfa]; end
rstats1.landalbdif = get_global_weighted_lc_stats(var,(1-lf).*area);
rstats2.landalbdif = get_global_weighted_stats(var(:,1),(1-lf).*area);

%lakes
for i = 1:length(o); sfa = o(i).lakesnowfreealb; var(i,:) = o(i).lakesnowalb - [sfa, sfa, 0, sfa, sfa]; end
rstats1.lakealbdif = get_global_weighted_lc_stats(var,lf.*area);
rstats2.lakealbdif = get_global_weighted_stats(var(:,1),lf.*area);

%difference
for i = 1:length(o); var(i,:) = o(i).albdif; end
rstats1.albdif = get_global_weighted_lc_stats(var,area);
rstats2.albdif = get_global_weighted_stats(var(:,1),area);


%3. MONTHLY CRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%land
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).landmonthlycre(p,:); end
    rstats1.landmonthlycre(p,:)  = get_global_weighted_lc_stats(var,(1-lf).*area);
    rstats2.landmonthlycre(p,:) = get_global_weighted_stats(var(:,1),(1-lf).*area);
end

%lakes
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).lakemonthlycre(p,:); end
    rstats1.lakemonthlycre(p,:) = get_global_weighted_lc_stats(var,lf.*area);
    rstats2.lakemonthlycre(p,:) = get_global_weighted_stats(var(:,1),lf.*area);
end

%total
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).totalmonthlycre(p,:); end
    rstats1.totalmonthlycre(p,:) = get_global_weighted_lc_stats(var,area);
    rstats2.totalmonthlycre(p,:) = get_global_weighted_stats(var(:,1),area);
end

%4. ANNUAL CRE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%land
for i = 1:length(o); var(i,:) = o(i).landtotalcre; end
rstats1.landtotalcre = get_global_weighted_lc_stats(var,(1-lf).*area);
rstats2.landtotalcre = get_global_weighted_stats(var(:,1),(1-lf).*area);

%lake
for i = 1:length(o); var(i,:) = o(i).laketotalcre; end
rstats1.laketotalcre = get_global_weighted_lc_stats(var,lf.*area);
rstats2.laketotalcre = get_global_weighted_stats(var(:,1),lf.*area);

%total
for i = 1:length(o); var(i,:) = o(i).totalcre; end
rstats1.totalcre = get_global_weighted_lc_stats(var,area);
rstats2.totalcre = get_global_weighted_stats(var(:,1),area);

% % 5. Lake Contribution Monthly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lake contribution
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylakecont(p,:); var2(i,:) = o(i).totalmonthlycre(p,:); end
    rstats1.monthlylakecont(p,:) = get_global_weighted_lc_stats(var,area.*abs(var2));
    rstats2.monthlylakecont(p,:) = get_global_weighted_stats(var(:,1),area.*abs(var2(:,1)));
end

%lake frac inc
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylakefracinc(p,:); var2(i,:) = o(i).totalmonthlycre(p,:); end
    rstats1.monthlylakefracinc(p,:) = get_global_weighted_lc_stats(var,area);
    rstats2.monthlylakefracinc(p,:) = get_global_weighted_stats(var(:,1),area);
end

%lake inc
for p = 1:12
    for i = 1:length(o)
        var(i,:) = o(i).monthlylakeinc(p,:);
        var2(i,:) = o(i).totalmonthlycre(p,:);
        II = isoutlier(var(i,:),'median',10);
        var(II) = NaN;
    end
    rstats1.monthlylakeinc(p,:) = get_global_weighted_lc_stats(var,area);
    rstats2.monthlylakeinc(p,:) = get_global_weighted_stats(var(:,1),area);
end

% % 6. Lake Contribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%total lake cont
for i = 1:length(o); var(i,:) = o(i).lakecont; var2(i,:) = o(i).totalcre; end
rstats1.lakecont = get_global_weighted_lc_stats(var,area.*abs(var2));
rstats2.lakecont = get_global_weighted_stats(var(:,1),area.*abs(var2(:,1)));

%total lake frac inc
for i = 1:length(o); var(i,:) = o(i).lakefracinc; var2(i,:) = o(i).totalcre; end
rstats1.lakefracinc = get_global_weighted_lc_stats(var,area);
rstats2.lakefracinc = get_global_weighted_stats(var(:,1),area);

%total lake inc
for i = 1:length(o); var(i,:) = o(i).lakeinc; var2(i,:) = o(i).totalcre; end
II = var(:,1) ~= 0;
rstats1.lakeinc = get_global_weighted_lc_stats(var(II,:),area(II));
rstats2.lakeinc = get_global_weighted_stats(var(II,1),area(II));

% % 7. Lake Contribution other calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lake contribution monthly
clear landcre lakecre lakefrac
for p = 1:12
    for i = 1:length(o); landcre(i,:) = o(i).landmonthlycre(p,:); lakecre(i,:) = o(i).lakemonthlycre(p,:); lakefrac(i,:) = o(i).lakefrac; end
    rlcstats.monthlylakecont(p,:) = get_global_weighted_lakecont_stats(landcre,lakecre,lakefrac,area);
    [rlcstats.monthlylakefracinc(p,:),rlcstats.monthlylakeinc(p,:)] = get_global_weighted_lakeinc_stats(landcre,lakecre,lakefrac,area);
end

%lake contribution annual
clear landcre lakecre
for i = 1:length(o); landcre(i,:) = o(i).landtotalcre; lakecre(i,:) = o(i).laketotalcre; lakefrac(i,:) = o(i).lakefrac; end
rlcstats.lakecont = get_global_weighted_lakecont_stats(landcre,lakecre,lakefrac,area);
[rlcstats.lakefracinc,rlcstats.lakeinc] = get_global_weighted_lakeinc_stats(landcre,lakecre,lakefrac,area);

% % Convert into tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r1(n,:),r12(n,:),r13(n,:),r2(n,:)] = get_biome_stats_tables(rstats1,rstats2,rlcstats,lakefrac,area);


% % store values for figures
for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylandsnow(p,:); end
    rstats1.landmonthlysnow(p,:)  = get_global_weighted_lc_stats(var,area);
    rstats2.landmonthlysnow(p,:) = get_global_weighted_stats(var(:,1),area);
end

for p = 1:12
    for i = 1:length(o); var(i,:) = o(i).monthlylakesnow(p,:); end
    rstats1.lakemonthlysnow(p,:)  = get_global_weighted_lc_stats(var,area);
    rstats2.lakemonthlysnow(p,:) = get_global_weighted_stats(var(:,1),area);
end

rs1(n) = rstats1;
rs2(n) = rstats2;
rlcs(n) = rlcstats;


end
