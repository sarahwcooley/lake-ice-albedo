function [b1,b12,b13,b2] = get_biome_stats_tables(stats1,stats2,lcstats,lakefrac,area)
n = 1;
%table b1
b1(n,1) = stats1.landsnowoff(1);
b1(n,2) = stats1.lakesnowoff(1);
b1(n,3) = stats1.snowoffdif(1);
b1(n,4) = stats1.landalbdif(1);
b1(n,5) = stats1.lakealbdif(1);
b1(n,6) = stats1.albdif(1);
b1(n,7) = stats1.landtotalcre(1);
b1(n,8) = stats1.laketotalcre(1);
b1(n,9) = stats1.totalcre(1);

%table b12
b12(n,1) = stats1.landsnowoff(3);
b12(n,2) = stats1.lakesnowoff(3);
b12(n,3) = stats1.snowoffdif(3);
b12(n,4) = stats1.landalbdif(3);
b12(n,5) = stats1.lakealbdif(3);
b12(n,6) = stats1.albdif(3);
b12(n,7) = stats1.landtotalcre(3);
b12(n,8) = stats1.laketotalcre(3);
b12(n,9) = stats1.totalcre(3);

%table b13
b13(n,1) = stats2.landsnowoff(3);
b13(n,2) = stats2.lakesnowoff(3);
b13(n,3) = stats2.snowoffdif(3);
b13(n,4) = stats2.landalbdif(3);
b13(n,5) = stats2.lakealbdif(3);
b13(n,6) = stats2.albdif(3);
b13(n,7) = stats2.landtotalcre(3);
b13(n,8) = stats2.laketotalcre(3);
b13(n,9) = stats2.totalcre(3);

%table b2
b2(n,1) = 100*nansum(0.01*lakefrac.*area)./nansum(area);
b2(n,2) = 10^-6*nansum(area);
b2(n,3) = lcstats.lakecont(1);
b2(n,4) = stats1.lakecont(3);
b2(n,5) = lcstats.lakefracinc(1);
b2(n,6) = lcstats.lakeinc(1);

end