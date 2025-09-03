function [out] = get_global_weighted_lakecont_stats(landcre,lakecre,lakefrac,area)
    for p = 1:5
        for i = 1:length(area)
            lf = 0.01*lakefrac(i);
            landtotal(i,1) = (1-lf)*area(i)*landcre(i,p);
            laketotal(i,1) = lf*area(i)*lakecre(i,p);
        end
    out(p,1) = 100*nansum(laketotal)./(nansum(laketotal)+nansum(landtotal));
    end
end
