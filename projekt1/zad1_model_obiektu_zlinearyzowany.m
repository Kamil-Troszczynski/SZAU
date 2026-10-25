%%%% Linearyzacja modelu nieliniowego obiektu
clear; clc; close all;

%% Parametry
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

% Folder na wykresy
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

%% Linearyzacja w punkcie pracy
lin.dFh1_dh1 = -p.a1/(4*p.C1*h1_0^(3/2));
lin.dFh1_dF1 =  1/(2*p.C1*h1_0);
lin.dFh1_dFD =  1/(2*p.C1*h1_0);
lin.dFh2_dh1 =  p.a1/(2*p.A2*sqrt(h1_0));
lin.dFh2_dh2 = -p.a2/(2*p.A2*sqrt(h2_0));

%% Dyskretyzacja
Ts = 0.5;                 % krok [s]
Tk = 3000;                % czas symulacji [s]
N  = round(Tk/Ts);
d  = round(p.tau/Ts);     % opóźnienie w próbkach
t    = (0:N)*Ts;

%% Eksperyment 1: te same zmiany wejść co w scenariuszu DMC (+/-50%)
F1in = F1_0*ones(1,N+1);
F1in(t >= 300)  = 1.5*F1_0;
F1in(t >= 1200) = 0.5*F1_0;
FD = FD_0*ones(1,N+1);
FD(t >= 700)  = 0.5*FD_0;
FD(t >= 2000) = 1.5*FD_0;

[h2, h2n] = sym_modele(p, lin, F1in, FD, Ts, d, F1_0, FD_0, h1_0, h2_0);

f = figure('Color', 'w', 'Visible', 'off');
subplot(2,1,1); plot(t, h2n, t, h2, '--'); grid on; ylabel('h_2 [cm]');
legend('Model nieliniowy',' Model liniowy','Location','best');
subplot(2,1,2); stairs(t, F1in); hold on; stairs(t, FD, '--'); hold off;
grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
legend('F_{1in}', 'F_D', 'Location', 'best');
exportgraphics(f, fullfile(folder, 'zad1_liniowy_vs_nieliniowy.png'), 'Resolution', 300); close(f);

%% Eksperyment 2: skoki F1in w dwóch miejscach - wzrost (+25%, +50%) i spadek (-25%, -50%)
kierunki = [1 -1];
tytuly   = {'Wzrost F_{1in}', 'Spadek F_{1in}'};
FD = FD_0*ones(1,N+1);

f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 700]);
for ik = 1:numel(kierunki)
    F1in = F1_0*ones(1,N+1);
    F1in(t >= 100)  = (1 + kierunki(ik)*0.25)*F1_0;
    F1in(t >= 1500) = (1 + kierunki(ik)*0.5)*F1_0;

    [h2, h2n] = sym_modele(p, lin, F1in, FD, Ts, d, F1_0, FD_0, h1_0, h2_0);

    subplot(2,2,ik);
    plot(t, h2n, t, h2, '--'); grid on; ylabel('h_2 [cm]'); title(tytuly{ik});
    legend('Model nieliniowy', 'Model liniowy', 'Location', 'best');
    subplot(2,2,ik+2);
    stairs(t, F1in); grid on; ylabel('F_{1in} [cm^3/s]'); xlabel('t [s]');
    ylim([0.4 1.6]*F1_0);
end
exportgraphics(f, fullfile(folder, 'zad1_liniowy_vs_nieliniowy_skoki.png'), 'Resolution', 300); close(f);

%% Eksperyment 3: pojedyncze skoki wejść o różnej wielkości i kierunku z punktu pracy (tylko tabela)
skoki_F1 = [-0.5 -0.25 -0.1 0.1 0.25 0.5];   % względem F1_0
skoki_FD = [-0.5 0.5];                       % względem FD_0
t_skok   = 100;

wejscia = [repmat({'F_1'}, 1, numel(skoki_F1)), repmat({'F_D'}, 1, numel(skoki_FD))];
skoki   = [skoki_F1, skoki_FD];
E   = zeros(size(skoki));          % suma kwadratów różnicy h2 liniowy - nieliniowy
dh2 = zeros(2, numel(skoki));      % zmiana h2 w stanie ustalonym: [nieliniowy; liniowy]

for is = 1:numel(skoki)
    F1in = F1_0*ones(1,N+1);
    FD   = FD_0*ones(1,N+1);
    if strcmp(wejscia{is}, 'F_1')
        F1in(t >= t_skok) = (1 + skoki(is))*F1_0;
    else
        FD(t >= t_skok) = (1 + skoki(is))*FD_0;
    end

    [h2, h2n] = sym_modele(p, lin, F1in, FD, Ts, d, F1_0, FD_0, h1_0, h2_0);
    E(is) = sum((h2 - h2n).^2);
    dh2(:, is) = [h2n(end) - h2_0; h2(end) - h2_0];
end

%% Tabela: jakość przybliżenia liniowego w zależności od wielkości skoku
fprintf('\nPorównanie modelu liniowego z nieliniowym (wyjście h2)\n');
fprintf('%8s %8s %14s %14s %12s\n', 'wejście', 'skok', 'dh2 nieliniowy', 'dh2 liniowy', 'E');
for is = 1:numel(skoki)
    fprintf('%8s %+7g%% %14.3f %14.3f %12.2f\n', wejscia{is}, 100*skoki(is), dh2(1,is), dh2(2,is), E(is));
end

%% Symulacja modelu liniowego i nieliniowego dla tych samych wejść
function [h2, h2n] = sym_modele(p, lin, F1in, FD, Ts, d, F1_0, FD_0, h1_0, h2_0)
    N = numel(F1in) - 1;
    h1  = h1_0*ones(1,N+1);  h2  = h2_0*ones(1,N+1);    % model liniowy
    h1n = h1_0*ones(1,N+1);  h2n = h2_0*ones(1,N+1);    % model nieliniowy
    for k = 1:N
        if k > d
            F1 = F1in(k-d);
        else
            F1 = F1_0;
        end

        % Model zlinearyzowany
        h1(k+1) = h1(k) + Ts*( lin.dFh1_dh1*(h1(k)-h1_0) + lin.dFh1_dF1*(F1-F1_0) + lin.dFh1_dFD*(FD(k)-FD_0) );
        h2(k+1) = h2(k) + Ts*( lin.dFh2_dh1*(h1(k)-h1_0) + lin.dFh2_dh2*(h2(k)-h2_0) );

        % Model nieliniowy
        F2 = p.a1*sqrt(h1n(k));
        F3 = p.a2*sqrt(h2n(k));
        h1n(k+1) = h1n(k) + Ts*(F1 + FD(k) - F2)/(2*p.C1*h1n(k));
        h2n(k+1) = h2n(k) + Ts*(F2 - F3)/p.A2;

        h1(k+1)  = max(h1(k+1),  1e-6);   h2(k+1)  = max(h2(k+1),  0);
        h1n(k+1) = max(h1n(k+1), 1e-6);   h2n(k+1) = max(h2n(k+1), 0);
    end
end
