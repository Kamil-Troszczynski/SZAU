%%%% Charakterystyka statyczna modelu obiektu

clear; clc; close all;

%% Parametry i punkt pracy
[model, F1_0, FD_0, h1_0, h2_0] = zad1_obiekt.parametry();

% Zakresy wejść
F1v = linspace(0, 2*F1_0, 200);      % sterowanie F1
FDv = linspace(0, 2*FD_0, 200);      % zakłócenie FD

% Folder na wykresy
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

%% Stan ustalony:

h1_stat = @(F1,FD) ((F1 + FD)/model.a1).^2;
h2_stat = @(F1,FD) ((F1 + FD)/model.a2).^2;

%% Wykres 1: h1, h2 w funkcji F1 dla kilku wartości FD
FD_list = FD_0 * [0 0.5 1 1.5 2];
f = figure('Name','Charakterystyki statyczne h(F1)', 'Color', 'w', 'Visible', 'off');
for i = 1:2
    subplot(1,2,i); hold on; grid on;
    for FD = FD_list
        if i == 1
            plot(F1v, h1_stat(F1v, FD), 'LineWidth', 1.5);
        else
            plot(F1v, h2_stat(F1v, FD), 'LineWidth', 1.5);
        end
    end
    if i == 1
        plot(F1_0, h1_0, 'ko', 'MarkerFaceColor', 'k');
        ylabel('h_1 [cm]');
        title('h_1(F_1)');
    else
        plot(F1_0, h2_0, 'ko', 'MarkerFaceColor', 'k');
        ylabel('h_2 [cm]'); 
        title('h_2(F_1)');
    end
    xlabel('F_1 [cm^3/s]');
    legend([compose('F_D = %.1f', FD_list), {'punkt pracy'}], 'Location', 'northwest');
end
exportgraphics(f, fullfile(folder, 'zad1_char_stat_h_F1.png'), 'Resolution', 300); close(f);

%% Wykres charakterystyki statycznej h2(F1in, Fd)
[F1g, FDg] = meshgrid(F1v, FDv);
f = figure('Name','Charakterystyka statyczna h2(F1,FD)', 'Color', 'w', 'Visible', 'off');
surf(F1g, FDg, h2_stat(F1g, FDg), 'EdgeColor', 'none');
hold on;
plot3(F1_0, FD_0, h2_0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('F_1 [cm^3/s]'); 
ylabel('F_D [cm^3/s]'); 
zlabel('h_2 [cm]');
title('h_2(F_1, F_D)'); 
colorbar; 
grid on; 
view(-35, 30);
exportgraphics(f, fullfile(folder, 'zad1_char_stat_h2_F1_FD.png'), 'Resolution', 300); close(f);