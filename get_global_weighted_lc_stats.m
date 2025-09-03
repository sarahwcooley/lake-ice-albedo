function [var] = get_global_weighted_lc_stats(varin,area)
    for i = 1:5
        lc = varin(:,i);
        II = ~isnan(lc) & ~isinf(lc);
        lc = lc(II);
        a = area(II);
        var(1,i) = mean(lc,Weights=a);
    end
end


    
