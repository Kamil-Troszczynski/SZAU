%% Parametry
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;

% Punkt pracy
F1_0 = 55;
FD_0 = 10;

% Zakresy wejść
F1v = linspace(0, 2*F1_0, 200);      % sterowanie F1
FDv = linspace(0, 2*FD_0, 200);      % zakłócenie FD

%% Stan ustalony:

h1_stat = @(F1,FD) ((F1 + FD)/p.a1).^2;
h2_stat = @(F1,FD) ((F1 + FD)/p.a2).^2;

h1_0 = h1_stat(F1_0, FD_0);
h2_0 = h2_stat(F1_0, FD_0);

%% Wykres 1: h1, h2 w funkcji F1 dla kilku wartości FD
FD_list = FD_0 * [0 0.5 1 1.5 2];
figure('Name','Charakterystyki statyczne h(F1)');
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
        ylabel('h_1 [cm]'); title('h_1(F_1)');
    else
        plot(F1_0, h2_0, 'ko', 'MarkerFaceColor', 'k');
        ylabel('h_2 [cm]'); title('h_2(F_1)');
    end
    xlabel('F_1 [cm^3/s]');
    legend([compose('F_D = %.1f', FD_list), {'punkt pracy'}], 'Location', 'northwest');
end

%% Wykres 2: powierzchnia h2(F1, FD)
[F1g, FDg] = meshgrid(F1v, FDv);
figure('Name','Charakterystyka statyczna h2(F1,FD)');
surf(F1g, FDg, h2_stat(F1g, FDg), 'EdgeColor', 'none');
hold on;
plot3(F1_0, FD_0, h2_0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
xlabel('F_1 [cm^3/s]'); ylabel('F_D [cm^3/s]'); zlabel('h_2 [cm]');
title('h_2(F_1, F_D)'); colorbar; grid on; view(-35, 30);