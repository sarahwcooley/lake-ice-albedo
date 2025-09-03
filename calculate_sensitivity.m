function out = calculate_sensitivity(temp,cre)

for i = 1:12
    x = temp(:,i);
    y = cre(:,i);
    JJ = isoutlier(y,'median');
    II = ~isnan(y) & ~JJ;
    [r,p] = corrcoef(x(II),y(II));
    P=polyfit(x(II),y(II),1); 
    if p(2) < 0.05
        out(i,1) = P(1);
    else
        out(i,1) = NaN;
    end
end


