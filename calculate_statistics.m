function out = calculate_statistics(var)
    if nansum(nansum(var)) ~= 0
        var = rmoutliers(var,'median');
        out = [nanmean(var),nanmedian(var),nanstd(var),prctile(var,25),prctile(var,75)];
    else
        out = nan(1,5);
    end
end