function [btsmean,btsmedian] = get_monthly_values(tsin,dates)
years = year(dates);
unyears = unique(years);
for i = 1:length(unyears)
    ind = find(years == unyears(i));
    days = dates(ind);
    ts = tsin(ind);

    months = month(days);
    unmonths = unique(months);
    count = 1;
    for j = 1:12
        ind2 = find(months == unmonths(j));
        tsm = ts(ind2);
        btsmean(i,j) = nanmean(tsm);
        btsmedian(i,j) = nanmedian(tsm);
    end

end
