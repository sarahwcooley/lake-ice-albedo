function [c] = get_sensitivity_table(tempresults)
    c = nan(8,12);
    for i = 1:length(tempresults)
        c(1,i) = tempresults(i).land(2);
        c(2,i) = tempresults(i).lake(2);
        if c(1,i) < 0.05; c(3,i) = tempresults(i).land(3); end 
        if c(2,i) < 0.05; c(4,i) = tempresults(i).lake(3); end
        if c(1,i) && c(2,i) < 0.05; c(5,i) = tempresults(i).lake(3)-tempresults(i).land(3); end
        if c(1,i) && c(2,i) < 0.05; c(6,i) = 100*(tempresults(i).lake(3)-tempresults(i).land(3))./tempresults(i).land(3); end
        if c(1,i) < 0.05; c(7,i) = table2array(tempresults(i).landmdl(2,2)); end
        if c(2,i) < 0.05; c(8,i) = table2array(tempresults(i).lakemdl(2,2)); end
    end
end