function [] = make_spider_plot(inputs, variable_name, scenario_name, index_file_courbes)
% Displays a diamond plot from inputs and legends

global slash

scenario_name(1) = strcat('Reference - ', ' ', scenario_name(1));

nb_scenario = size(inputs,1);
nb_variable = size(inputs,2);

% First scenario is the reference scenario (normalized to 1)
r_data = zeros(size(inputs));
r_data(1,:) = ones(1,nb_variable);
r_data(2:size(r_data,1),:) = inputs(2:size(inputs,1), :) ./ repmat(inputs(1,:), size(inputs,1) - 1, 1);
max_range = max(max(r_data));

angle = 2 .* pi .* (0:nb_variable - 1) ./ nb_variable;

X = (r_data .* cos(angle))';
Y = (r_data .* sin(angle))';

graduation = [0.5; 1.5; 2];
nameGraduation = {'-50%', '+50%','+100%'};
xGraduation = (graduation .* cos(angle))';
yGraduation = (graduation .* sin(angle))';
xGraduation = [xGraduation; graduation' .* cos(angle(1))];
yGraduation = [yGraduation; graduation' .* sin(angle(1))];


%colours = colormap(parula);
%colours = colours(round((1:(nb_scenario-1))./nb_scenario.*size(colours,1)),:);
colours = summer(nb_scenario - 1);

% Contours
X_out = [X; X(1,:)];
Y_out = [Y; Y(1,:)];

% Plot
h.ax = axes('xlim', [-max_range max_range], 'ylim', [-max_range max_range], 'dataaspectratio', [1 1 1], 'visible', 'off');
hold on
for i = 2:nb_scenario
    h.poly = patch(X(:,i), Y(:,i), colours(i-1,:));
    set(h.poly, 'EdgeColor', 'none', 'FaceAlpha', 0.08, 'HandleVisibility','off')
end
plot([zeros(1,nb_variable); max_range.*cos(angle)], [zeros(1,nb_variable); max_range.*sin(angle)], ':k', 'HandleVisibility','off');
plot(xGraduation, yGraduation, 'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility','off');
r = plot(X_out(:,1), Y_out(:,1), 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.5);
k = plot(X_out(:,2:end), Y_out(:,2:end), 'LineWidth', 1.5);
set(k, {'color'}, num2cell(colours,2));
rk = [r; k];
lgd = legend(rk, scenario_name);
lgd.Location = 'westoutside';
lgd.FontSize = 14;
legend('boxoff');
% Add name of the variables
h.axlbl = text(1.2*max_range*cos(angle'), 1.2*max_range*sin(angle'), variable_name, 'horiz', 'center', 'vert', 'top', 'FontSize', 14);
h.axlbl2 = text(graduation, 0.*graduation + 0.15, nameGraduation, 'horiz', 'center', 'vert', 'top', 'FontSize', 10);
set(gcf,'units','centimeters','position',[0,0,30,15])
hold off 

sauvegarde(strcat('.', slash, index_file_courbes, slash, 'spider_plot'))

end

function sauvegarde(nom)



    format_d='-dpng ';
    format2='.png';
            
    print('-f1', format_d, [nom,format2]);
    saveas(1,[nom,'.fig']);

end