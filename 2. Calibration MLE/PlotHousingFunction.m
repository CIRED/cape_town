%% Print and export figure for construction function

% Create figure
fig = figure('units','centimeters','position',[0,0,30,15]);

colorSimul = 0.9*[0 1 0];

% Subplot 1
subplot(1,2,1)
hold on 
scatter1 = scatter(data.sp2011Distance(selectedDensity), dataDensity(selectedDensity)./10000, 100, '.','MarkerEdgeColor', 0.6*[1 1 1]);
% Fit of Data
plot1 = plot(1:35, f1(1:35), 'LineWidth',2.3, 'Color', 0.3*[1 1 1]);
p1 = predint(f1,1:35,0.95,'functional','on');
plot(1:35, p1, 'LineWidth',0.5, 'Color', 0.3*[1 1 1])
% Fit of CD
plot2 = plot(1:35, f2(1:35), 'LineWidth',3.4, 'Color', colorSimul);
p2 = predint(f2,1:35,0.95,'functional','on');
plot(1:35, p2, ':', 'LineWidth',1.5, 'Color', 0.8.*colorSimul)
xlim([0 35])
ylim([0 0.6])
legend([scatter1 plot1 plot2],'Data','Fit of data', 'Fit of simulated value')
ylabel('Floor area ratio')
xlabel('Distance from the center') 
ax = gca;
ax.FontSize = 14;
title('Cobb-Douglas function')

% Subplot 2
subplot(1,2,2)
hold on 
scatter1 = scatter(data.sp2011Distance(selectedDensity), dataDensity(selectedDensity)./10000, 100, '.','MarkerEdgeColor', 0.6*[1 1 1]);
% Fit of Data
plot1 = plot(1:35, f1(1:35), 'LineWidth',2.3, 'Color', 0.3*[1 1 1]);
p1 = predint(f1,1:35,0.95,'functional','on');
plot(1:35, p1, 'LineWidth',0.5, 'Color', 0.3*[1 1 1])
% Fit of CES
plot2 = plot(1:35, f3(1:35), 'LineWidth',3.4, 'Color', colorSimul);
p3 = predint(f3,1:35,0.95,'functional','on');
plot(1:35, p3, ':', 'LineWidth',1.5, 'Color', 0.8.*colorSimul)
xlim([0 35])
ylim([0 0.6])
legend([scatter1 plot1 plot2],'Data','Fit of data', 'Fit of simulated value')
xlabel('Distance from the center') 
ax = gca;
ax.FontSize = 14;
title('CES function')

nameFile = strcat('.', slash, nameFolder, slash, 'housingProductionFunction_CD_CES');
format_d = '-dpng';
format2 = '.png';
print('-f1', format_d, [nameFile,format2]);
saveas(1,[nameFile,'.fig']);


close all
    

