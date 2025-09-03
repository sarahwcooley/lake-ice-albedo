function [lakeoutput,landoutput] = read_in_csv_time_series_v1(folder,ids)
%order of CSVs: Cloud, Snow, Water, Snow Albedo, Water Albedo, Snow Albedo
%sigma, Water Albedo sigma, Date
cd(folder)
%read in lakes
for n = 1:length(ids)

    %LAKES %%%%%
    fname = ['grid1*_lakes_' num2str(ids(n)) '.csv'];
    files = dir(fname);
    for j = 1:length(files)
        input = readcell(files(j).name,'Delimiter',',');
        c1 = input(:,4);
        c2 = input(:,5);
        ind = find([c1{:}] == 3 & [c2{:}] == 4);
        ind = ind';
        input(ind,:) = [];
        clear o d
        for i = 1:length(input
            for t = 1:7
                o(i,t) = input{i,t}; 
            end
            d(i,1) = input{i,8};
        end
        if j == 1
            lakes = o;
            datelakes = d;
        else
            lakes = cat(1,lakes,o);
            datelakes = cat(1,datelakes,d);
        end
    end    

    %LAND %%%%%
    fname = ['grid1*_land_' num2str(ids(n)) '.csv'];
    files = dir(fname);
    if isempty(files) == 0
    for j = 1:length(files)
        input = readcell(files(j).name,'Delimiter',',');
        c1 = input(:,4);
        c2 = input(:,5);
        ind = find([c1{:}] == 3 & [c2{:}] == 4);
        input(ind,:) = [];
        clear o d
        for i = 1:length(input)
            for t = 1:7
                o(i,t) = input{i,t}; 
            end
            d(i,1) = input{i,8};
        end
        if j == 1
            land = o;
            dateland = d;
        else
            land = cat(1,land,o);
            dateland = cat(1,dateland,d);
        end
    end    

%correct for missing dates
    tdays = datetime(2001,1,1):datetime(2022,12,31);
    for i = 1:length(tdays)
        td = tdays(i);
        ind = find(dateland == td);
        if length(ind) > 1
            ind = ind(1);
        end
        if isempty(ind) == 0
            outland(i,:) = land(ind,:);
        else
            outland(i,:) = NaN;
        end
    end


    for i = 1:length(tdays)
        td = tdays(i);
        ind = find(datelakes == td);
        if length(ind) > 1
            ind = ind(1);
        end
        if isempty(ind) == 0
            outlakes(i,:) = lakes(ind,:);
        else
            outlakes(i,:) = NaN;
        end
    end

    for pp = 4:7
        ind = find(outland(:,pp) ~= -999);
        outland(ind,pp) = 0.001*outland(ind,pp);
        ind = find(outlakes(:,pp) ~= -999);
        outlakes(ind,pp) = 0.001*outlakes(ind,pp);
    end


    landoutput(n).id = ids(n);
    landoutput(n).cloud_ts = outland(:,1);
    landoutput(n).snow_ts = outland(:,2);
    landoutput(n).water_ts = outland(:,3);
    landoutput(n).snowa_ts = outland(:,4);
    landoutput(n).watera_ts = outland(:,5);
    landoutput(n).snowa_std = outland(:,6);
    landoutput(n).watera_std = outland(:,7);
    landoutput(n).days = tdays';

    lakeoutput(n).id = ids(n);
    lakeoutput(n).cloud_ts = outlakes(:,1);
    lakeoutput(n).snow_ts = outlakes(:,2);
    lakeoutput(n).water_ts = outlakes(:,3);
    lakeoutput(n).snowa_ts = outlakes(:,4);
    lakeoutput(n).watera_ts = outlakes(:,5);
    lakeoutput(n).snowa_std = outlakes(:,6);
    lakeoutput(n).watera_std = outlakes(:,7);
    lakeoutput(n).days = tdays';
    disp(['Finished ' num2str(n)])
    else
    landoutput(n).id = ids(n);
    landoutput(n).cloud_ts = NaN;
    landoutput(n).snow_ts = NaN;
    landoutput(n).water_ts = NaN;
    landoutput(n).snowa_ts = NaN;
    landoutput(n).watera_ts = NaN;
    landoutput(n).snowa_std = NaN;
    landoutput(n).watera_std = NaN;
    landoutput(n).days = NaN;

    lakeoutput(n).id = ids(n);
    lakeoutput(n).cloud_ts = NaN;
    lakeoutput(n).snow_ts = NaN;
    lakeoutput(n).water_ts = NaN;
    lakeoutput(n).snowa_ts = NaN;
    lakeoutput(n).watera_ts = NaN;
    lakeoutput(n).snowa_std = NaN;
    lakeoutput(n).watera_std = NaN;
    lakeoutput(n).days = tdays';
    disp(['Missing file ' num2str(n)])

    end
end