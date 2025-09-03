function [mts,x,inBetween] = variables_for_fillplot_v2(var)
    mts = var(:,1)';
    mthigh = var(:,1)'+ var(:,3)';
    mtlow = var(:,1)' - var(:,3)';
    x1 = (1:12);
    x1(isnan(mts)) = [];
    x = [x1,fliplr(x1)];
    inBetween = [mthigh(~isnan(mts)),fliplr(mtlow(~isnan(mts)))];
end
  