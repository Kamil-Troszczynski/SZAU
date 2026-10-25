%%%% Regulacja DMC zaprojektowana na modelu liniowym sterująca obiektem nieliniowym
clear; clc; close all;

%% Parametry obiektu
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

% Folder na wykresy
folder = fullfile(fileparts(mfilename('fullpath')), 'zad1_plots');

%% Parametry regulatora DMC
Tp_reg   = 10;
D        = 60;
N        = 30;
Nu       = 5;
lambda   = 0.7;
limit_u  = [0 110];
limit_du = [-10 10];

%% Parametry do badania wpływu nastaw
lambda_lista = [0.1 0.7 2 5];
N_lista      = [10 20 30 60];
Nu_lista     = [1 2 5 10];
D_lista      = [30 45 60 90];

%% Linearyzacja w punkcie pracy
dFh1_dh1 = -p.a1/(4*p.C1*h1_0^(3/2));
dFh1_dF1 =  1/(2*p.C1*h1_0);
dFh2_dh1 =  p.a1/(2*p.A2*sqrt(h1_0));
dFh2_dh2 = -p.a2/(2*p.A2*sqrt(h2_0));

%% Dyskretyzacja
Ts   = 0.5;
Tk   = 3000;
N_sim = round(Tk/Ts);
d    = round(p.tau/Ts);
t    = (0:N_sim)*Ts;

%% Odpowiedź skokowa modelu liniowego (liczona do najdłuższego badanego D)
co = round(Tp_reg/Ts);
D_max = max([D D_lista]);
s  = zeros(D_max,1);
dh1 = 0; dh2 = 0;
for k = 1:D_max*co
    du  = double(k > d);
    dh1_new = dh1 + Ts*( dFh1_dh1*dh1 + dFh1_dF1*du );
    dh2 = dh2 + Ts*( dFh2_dh1*dh1 + dFh2_dh2*dh2 );
    dh1 = dh1_new;
    if mod(k, co) == 0
        s(k/co) = dh2;
    end
end

% Wzmocnienie statyczne modelu liniowego - sprawdzenie, czy D obejmuje stan ustalony
K_stat = (dFh1_dF1/(-dFh1_dh1)) * (dFh2_dh1/(-dFh2_dh2));
fprintf('Odpowiedź skokowa po D = %d krokach: %.2f%% wartości ustalonej\n', D, 100*s(D)/K_stat);

f = figure('Color', 'w', 'Visible', 'off');
stairs((1:D_max)*Tp_reg, s); 
hold on;
yline(K_stat, '--'); 
xline(D*Tp_reg, ':'); 
hold off;
grid on; 
ylabel('s'); 
xlabel('t [s]');
legend('Odpowiedź skokowa - model liniowy', 'wartość ustalona', sprintf('D = %d', D), 'Location', 'best');
exportgraphics(f, fullfile(folder, 'zad1_dmc_odp_skokowa.png'), 'Resolution', 300); close(f);

%% Scenariusz testowy
yzad = h2_0*ones(1,N_sim+1);
yzad(t >= 300)  = 1.5*h2_0;
yzad(t >= 1200) = 0.5*h2_0;
FDc = FD_0*ones(1,N_sim+1);
FDc(t >= 700) = 0.5*FD_0;
FDc(t >= 2000) = 1.5*FD_0;

sc.p = p; sc.yzad = yzad; sc.FD = FDc; sc.Ts = Ts; sc.co = co; sc.d = d;
sc.F1_0 = F1_0; sc.h1_0 = h1_0; sc.h2_0 = h2_0;
sc.limit_u = limit_u; sc.limit_du = limit_du;

%% Regulacja DMC obiektu nieliniowego
[ke, ku] = dmc_macierze_ke_ku(s(1:D), N, Nu, D, lambda);
[h2c, u, E] = sym_dmc(ke, ku, D, sc);
fprintf('Wskaźnik jakości regulacji DMC: %.2f\n', E);

%% Wykresy
ip = 1:co:N_sim+1;

f = figure('Color', 'w', 'Visible', 'off');
subplot(2,1,1);
stairs(t(ip), h2c(ip)); hold on;
stairs(t(ip), yzad(ip), '--'); hold off;
grid on; ylabel('h_2 [cm]');
legend('wyjście h_2','wartość zadana','Location','best');

subplot(2,1,2);
stairs(t(ip), u(ip)); 
hold on;
stairs(t(ip), FDc(ip), '--'); 
hold off;
grid on; 
ylabel('F [cm^3/s]'); 
xlabel('t [s]');
legend('sterowanie F_{1in}','zakłócenie F_D','Location','best');
exportgraphics(f, fullfile(folder, 'zad1_dmc_regulacja.png'), 'Resolution', 300); close(f);

%% Wpływ nastaw na jakość regulacji - zmiana jednego parametru
nastawy = [lambda N Nu D];
badania = {'\lambda', 'lambda', lambda_lista;
           'N',       'N',      N_lista;
           'N_u',     'Nu',     Nu_lista;
           'D',       'D',      D_lista};

for ib = 1:size(badania, 1)
    wart = badania{ib,3};
    E_b  = zeros(size(wart));
    f = figure('Color', 'w', 'Visible', 'off');
    for iw = 1:numel(wart)
        n = nastawy; n(ib) = wart(iw);       % n = [lambda N Nu D]
        [ke_b, ku_b] = dmc_macierze_ke_ku(s(1:n(4)), n(2), n(3), n(4), n(1));
        [h2b, ub, E_b(iw)] = sym_dmc(ke_b, ku_b, n(4), sc);
        subplot(2,1,1); 
        stairs(t(ip), h2b(ip)); 
        hold on;
        
        subplot(2,1,2); 
        stairs(t(ip), ub(ip));  
        hold on;
    end
    etyk = arrayfun(@(w, e) sprintf('%s = %g (E = %.1f)', badania{ib,1}, w, e), wart, E_b, 'UniformOutput', false);
    subplot(2,1,1); 
    stairs(t(ip), yzad(ip), 'k--'); 
    hold off;
    grid on; 
    ylabel('h_2 [cm]');
    legend([etyk, {'wartość zadana'}], 'Location', 'best');

    subplot(2,1,2); 
    hold off;
    grid on; 
    ylabel('F_{1in} [cm^3/s]'); 
    xlabel('t [s]');
    exportgraphics(f, fullfile(folder, sprintf('zad1_dmc_strojenie_%s.png', badania{ib,2})), 'Resolution', 300); close(f);

    fprintf('\nWpływ %s na wskaźnik jakości E:\n', badania{ib,2});
    fprintf('%10g  E = %.2f\n', [wart; E_b]);
end 

%% Macierze regulatora DMC w wersji analitycznej
function [ke, ku] = dmc_macierze_ke_ku(s, N, Nu, D, lambda)
    M = zeros(N, Nu);
    for i = 1:N
        for j = 1:min(i, Nu)
            M(i,j) = s(min(i-j+1, D));
        end
    end
    Mp = zeros(N, D-1);
    for i = 1:N
        for j = 1:D-1
            Mp(i,j) = s(min(i+j, D)) - s(j);
        end
    end
    K  = (M'*M + lambda*eye(Nu)) \ M';
    ke = sum(K(1,:));
    ku = K(1,:) * Mp;
end

% Wyznaczenie wartości sterującej w jednym kroku oraz zaktualizowanie
% wektora przechowującego zmiany wartości sterujących
function [u, dUp] = dmc_krok(yzad, y, u_prev, dUp, ke, ku, limit_u, limit_du)
    du = ke*(yzad - y) - ku*dUp;
    du = min(max(du, limit_du(1)), limit_du(2));
    u = min(max(u_prev + du, limit_u(1)), limit_u(2));
    dUp = [u - u_prev; dUp(1:end-1)];
end

% Symulacja układu regulacji DMC z obiektem nieliniowym, E - suma kwadratów
% błędów w chwilach próbkowania regulatora
function [h2, u, E] = sym_dmc(ke, ku, D, sc)
    p = sc.p;
    n = numel(sc.yzad) - 1;
    h1 = sc.h1_0*ones(1,n+1);  h2 = sc.h2_0*ones(1,n+1);
    u = sc.F1_0*ones(1,n+1);
    u_akt = sc.F1_0;
    dUp = zeros(D-1,1);

    for k = 1:n
        if mod(k-1, sc.co) == 0
            [u_akt, dUp] = dmc_krok(sc.yzad(k), h2(k), u_akt, dUp, ke, ku, sc.limit_u, sc.limit_du);
        end
        u(k) = u_akt;
        if k > sc.d
            F1 = u(k-sc.d);
        else
            F1 = sc.F1_0;
        end

        % Obiekt nieliniowy
        F2 = p.a1*sqrt(h1(k));
        F3 = p.a2*sqrt(h2(k));
        h1(k+1) = h1(k) + sc.Ts*(F1 + sc.FD(k) - F2)/(2*p.C1*h1(k));
        h2(k+1) = h2(k) + sc.Ts*(F2 - F3)/p.A2;

        h1(k+1) = max(h1(k+1), 1e-6);
        h2(k+1) = max(h2(k+1), 0);
    end
    u(end) = u(end-1);

    ip = 1:sc.co:n+1;
    E = sum((sc.yzad(ip) - h2(ip)).^2);
end
