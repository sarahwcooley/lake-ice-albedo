function varout = calculate_statistics_monthly(var)  

if nansum(nansum(var)) ~= 0
     for p = 1:12
         v = var(:,p);
         v = rmoutliers(v,'median');
         varout(p,1) = nanmean(v);
         varout(p,2) = nanmedian(v);
         varout(p,3) = nanstd(v);
         varout(p,4) = prctile(v,25);
         varout(p,5) = prctile(v,75);
      end
else
    varout = nan(12,5);
end
