%READ IN CSV DATA

folder = '/Volumes/SWC_SSD/0_Lake_Ice_TimeSeries/0_final_grid1';
ids = 1:2050;

%process data
[lakeoutput,landoutput] = read_in_csv_time_series_v1(folder,ids);

outfolder = '/Users/sc961/Library/CloudStorage/GoogleDrive-cooleysarahw@gmail.com/My Drive/Research/Lake_Ice_Albedo_project/0_results/';
cd(outfolder)
save('grid1_results_apr29_v1.mat','landoutput');

%check for missing data (sees if any grid cells are missing years)

[missing_lakes] = check_for_missing_data(lakeoutput);
[missing_land] = check_for_missing_data(landoutput);

ind = find([missing_lakes.sum]' > 56);
missing_lakes = missing_lakes(ind);

ind = find([missing_land.sum]' > 56);
missing_land = missing_land(ind);

for i = 1:length(missing_land)
    m = missing_land(i).myears;
    ind = find(m(:,2) > 360);
    missing_land(i).fullyears = m(ind,1);
end

count = 1;
for i = 1:length(missing_land);
    f = missing_land(i).fullyears;
    for j = 1:length(f)
        fyears(count,1) = missing_land(i).id;
        fyears(count,2) = f(j);
        count = count + 1;
    end
end

unyears = unique(fyears(:,2));
for i = 1:length(unyears)
    ind = find(fyears(:,2) == unyears(i));
    od(i).year = unyears(i);
    od(i).ids = fyears(ind,1);
end