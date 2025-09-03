function [h] = make_map_plot(ss,var,thresh,cblabel,name,cmap)

figure(5)
hold off
axesm('MapProjection','robinson','MapLatLimit',[25 90],'MapLonLimit',[-180 180]);
axis off; %framem on; framem('FLineWidth',1); 
tightmap; set(gca, 'Position', [0 0 1 1]);
load coastlines
geoshow(coastlat,coastlon,'DisplayType','line','Color','k')
hold on
geoshow(ss(isnan(var)),'DisplayType', 'polygon', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor',[0 0 0],'LineWidth',0.0005);
[ids, coloridx] = get_color_idx(var,thresh);
s = ss(ids);
for k = 1:length(thresh)-1
    s1 = s(coloridx == k);
    faceColor = cmap(k, :);
    geoshow(s1, 'DisplayType', 'polygon', 'FaceColor', faceColor, 'EdgeColor',[0 0 0],'LineWidth',0.0005);
end
gcf = figure(5);
fig = gcf;
clim([1 length(thresh)]);
colormap(cmap);
cb = colorbar;
cb.Ticks = [1:length(thresh)];
for p = 1:length(thresh)
    if p == 1; cb.TickLabels{p,1} = ' '; else; cb.TickLabels{p,1} = num2str(thresh(p));end
    if p == length(thresh); cb.TickLabels{p,1} = ' '; end
end
set(cb, 'Position', [0.91 0.15 0.015 0.57]);  
cb.Label.String = cblabel;
set(gca,'FontSize',10)
set(fig, 'Units', 'inches', 'Position', [1 1 11 2.5]);
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [11 2.5]);
set(fig, 'PaperPosition', [0 0 11 2.5]); 
print(fig, name, '-dpng', '-r500')
h = 1;
end