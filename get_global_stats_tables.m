function [t1,t2,t3,t4] = get_global_stats_tables(stats1,stats2,lcstats)
    %Stats 1 table
    clear t
    s = stats1;
    t(1,:) = s.landsnowoff;
    t(2,:) = s.lakesnowoff;
    t(3,:) = s.snowoffdif;
    t(4,:) = s.landalbdif;
    t(5,:) = s.lakealbdif;
    t(6,:) = s.albdif;
    t(7,:) = s.landtotalcre;
    t(8,:) = s.laketotalcre;
    t(9,:) = s.totalcre;
    t1 = t;

    %Stats 2 table
    s = stats2;
    t(1,:) = s.landsnowoff;
    t(2,:) = s.lakesnowoff;
    t(3,:) = s.snowoffdif;
    t(4,:) = s.landalbdif;
    t(5,:) = s.lakealbdif;
    t(6,:) = s.albdif;
    t(7,:) = s.landtotalcre;
    t(8,:) = s.laketotalcre;
    t(9,:) = s.totalcre;
    t2 = t;

    %Table for Lake Cont
    clear t
    t(1,:) = stats1.lakecont;
    t(2,:) = lcstats.lakecont;
    t(3,:) = stats1.lakefracinc;
    t(4,1:2) = lcstats.lakefracinc;
    t(5,:) = stats1.lakeinc;
    t(6,:) = stats2.lakeinc;
    t(7,1:2) = lcstats.lakeinc;
    t3 = t;

    %Table for Monthly Data
    clear t
    t(:,1) = stats1.landmonthlycre(:,1);
    t(:,2) = stats1.landmonthlycre(:,3);
    t(:,3) = stats2.landmonthlycre(:,3);
    t(:,4) = stats1.lakemonthlycre(:,1);
    t(:,5) = stats1.lakemonthlycre(:,3);
    t(:,6) = stats2.lakemonthlycre(:,3);
    t(:,7) = stats1.totalmonthlycre(:,1);
    t(:,8) = stats1.totalmonthlycre(:,3);
    t(:,9) = stats2.totalmonthlycre(:,3);
    t(:,10) = lcstats.monthlylakecont(:,1);
    t(:,11) = lcstats.monthlylakefracinc(:,1);
    t(:,12) = lcstats.monthlylakeinc(:,1);
    t4 = t';
end