function [tsout] = filter_ts_albedo(tsin,cloud,valcounts,maxpixels,dates,cloudthresh,pixelthresh)

%filter clouds
II = cloud > cloudthresh & 100*valcounts./maxpixels < pixelthresh;
tsin(II) = NaN;
II = isnan(tsin);
tsy = tsin;
tsy(II) = [];
dates(II) = [];
qdays = datetime(2001,1,1):datetime(2022,12,31);
daily = interp1(dates,tsy,qdays,'linear');
daily = medfilt1(daily,5,'truncate');
tsout = daily';
end


