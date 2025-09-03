function [out1,out2] = get_global_weighted_lakeinc_stats(landcre,lakecre,lakefrac,area)
    for p = 1:2
        for i = 1:length(area)
            lf = 0.01*lakefrac(i);
            landtotal(i,1) = (1-lf)*area(i)*landcre(i,p);
            laketotal(i,1) = lf*area(i)*lakecre(i,p);
            lakearea(i,1) = lf*area(i);
        end
        lc = (nansum(laketotal)./(nansum(laketotal)+nansum(landtotal)));
        lf = nansum(lakearea)./nansum(area);
        out1(p,1) = 100*lc./lf;
        II = lakefrac >= 0.1;
        lakec = lakecre(:,p);
        landc = landcre(:,p);
        out2(p,1) = 100*(nansum(lakec(II).*area(II))-nansum(landc(II).*area(II)))./nansum(landc(II).*area(II));
    end

end
    
 