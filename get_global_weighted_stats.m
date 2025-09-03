function [output] = get_global_weighted_stats(var,area)
    II = ~isnan(var);
    var = var(II);
    area = area(II);
    output(1,1) = mean(var,Weights=area);
    output(1,2) = median(var,Weights=area);
    output(1,3)  = std(var,area);
    output(1,4:5)  = wprctile(var,[25 75],area,5);
end
    
