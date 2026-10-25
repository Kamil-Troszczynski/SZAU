%%%% Rozmyty regulator DMC w wersji analitycznej z ograniczeniami (podejście PDC),
%%%% porównanie z konwencjonalnym DMC z zadania 1 oraz strojenie funkcji przynależności
%%%% regulatora. Wymaga pliku zad2_model_rozmyty.mat tworzonego przez zad2_modele_rozmyte.m
clear; clc; close all;

%% Parametry obiektu
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

% Folder na wykresy i plik z funkcjami przynależności modelu
katalog     = fileparts(mfilename('fullpath'));
folder      = fullfile(katalog, 'zad2_plots');
plik_modelu = fullfile(katalog, 'zad2_model_rozmyty.mat');
if ~isfile(plik_modelu)
    error('Brak pliku %s - najpierw uruchom zad2_modele_rozmyte.m', plik_modelu);
end
dane = load(plik_modelu);          % mf_wyb, mf_popr, h2_min, h2_max

%% Parametry regulatorów
Tp_reg   = 10;
limit_u  = [0 110];
limit_du = [-10 10];

% Konwencjonalny DMC - nastawy jak w zadaniu 1
konw.D = 60; konw.N = 30; konw.Nu = 5; konw.lambda = 0.7;

% Rozmyty DMC - te same N, Nu, lambda; dłuższe D, bo przy wysokim h2 obiekt jest wolniejszy
rozm.D = 100; rozm.N = 30; rozm.Nu = 5; rozm.lambda = 0.7;

% PDC: początkowe funkcje przynależności regulatora takie same jak w modelu obiektu
mf_start = dane.mf_popr;

% Warianty kształtu funkcji przynależności sprawdzane przy strojeniu regulatora
g_lista    = [1 1.25 1.5 2 2.5];
szer_lista = [0.3 0.5 0.7 1];

%% Dyskretyzacja
Ts    = 0.5;
Tk    = 3000;
N_sim = round(Tk/Ts);
d     = round(p.tau/Ts);
t     = (0:N_sim)*Ts;
co    = round(Tp_reg/Ts);
ip    = 1:co:N_sim+1;

%% Scenariusz testowy (jak w zadaniu 1)
yzad = h2_0*ones(1,N_sim+1);
yzad(t >= 300)  = 1.5*h2_0;
yzad(t >= 1200) = 0.5*h2_0;
FDc = FD_0*ones(1,N_sim+1);
FDc(t >= 700)  = 0.5*FD_0;
FDc(t >= 2000) = 1.5*FD_0;

sc.p = p; sc.yzad = yzad; sc.FD = FDc; sc.Ts = Ts; sc.co = co; sc.d = d;
sc.F1_0 = F1_0; sc.h1_0 = h1_0; sc.h2_0 = h2_0;
sc.limit_u = limit_u; sc.limit_du = limit_du;

%% Regulator konwencjonalny: jeden model lokalny w punkcie pracy
mf_lin   = struct('L', 1, 'typ', 'trojkatne', 'g', 1, 'szer', 0, 'c', h2_0);
reg_konw = projekt_dmc(p, mf_lin, konw, Ts, co, d);
[h2_k, u_k, E_k] = sym_dmc_rozmyty(reg_konw, sc);

%% Regulator rozmyty (PDC): lokalne regulatory DMC zaprojektowane w centrach zbiorów
reg_start = projekt_dmc(p, mf_start, rozm, Ts, co, d);
[h2_r, u_r, E_r] = sym_dmc_rozmyty(reg_start, sc);

fprintf('DMC konwencjonalny:                  E = %.2f\n', E_k);
fprintf('DMC rozmyty (%s): E = %.2f\n', zad2_opis_mf(mf_start), E_r);

% Sprawdzenie, czy horyzont D obejmuje stan ustalony każdego modelu lokalnego
fprintf('\nOdpowiedzi skokowe modeli lokalnych po D = %d krokach:\n', rozm.D);
for i = 1:numel(mf_start.c)
    m = zad2_lokalny(p, mf_start.c(i));
    K = (m.b/(-m.d11))*(m.d21/(-m.d22));
    fprintf('h2 = %5.2f cm: K = %.4f, s(D) = %.2f%% K\n', m.h2, K, 100*reg_start.s(end,i)/K);
end

% Odpowiedzi skokowe modeli lokalnych
f = figure('Color', 'w', 'Visible', 'off');
stairs((1:rozm.D)*Tp_reg, reg_start.s); grid on; ylabel('s'); xlabel('t [s]');
legend(arrayfun(@(c) sprintf('model lokalny h_2 = %.2f cm', c), mf_start.c, 'UniformOutput', false), ...
       'Location', 'best');
exportgraphics(f, fullfile(folder, 'zad2_dmc_odp_skokowe_lokalne.png'), 'Resolution', 300); close(f);

rysuj_regulacje(t, ip, yzad, FDc, {h2_k, h2_r}, {u_k, u_r}, ...
    {sprintf('DMC konwencjonalny (E = %.1f)', E_k), sprintf('DMC rozmyty (E = %.1f)', E_r)}, ...
    fullfile(folder, 'zad2_dmc_konw_vs_rozmyty.png'));

%% Strojenie funkcji przynależności regulatora (ta sama liczba zbiorów)
warianty = zad2_warianty_mf(numel(mf_start.c), g_lista, szer_lista, dane.h2_min, dane.h2_max);
nW       = numel(warianty);
E_war    = zeros(1, nW);
for iw = 1:nW
    reg = projekt_dmc(p, warianty{iw}, rozm, Ts, co, d);
    [~, ~, E_war(iw)] = sym_dmc_rozmyty(reg, sc);
end
[E_best, iw_best] = min(E_war);
if E_best < E_r
    mf_stroj = warianty{iw_best};
else
    mf_stroj = mf_start;           % żaden wariant nie poprawił jakości regulacji
end
reg_stroj = projekt_dmc(p, mf_stroj, rozm, Ts, co, d);
[h2_s, u_s, E_s] = sym_dmc_rozmyty(reg_stroj, sc);

fprintf('\nStrojenie funkcji przynależności regulatora rozmytego\n');
for iw = 1:nW
    fprintf('%-40s E = %10.2f\n', zad2_opis_mf(warianty{iw}), E_war(iw));
end
fprintf('Przed strojeniem: %-40s E = %.2f\n', zad2_opis_mf(mf_start), E_r);
fprintf('Po strojeniu:     %-40s E = %.2f\n', zad2_opis_mf(mf_stroj), E_s);
fprintf('DMC konwencjonalny:                                        E = %.2f\n', E_k);

% Funkcje przynależności regulatora przed i po strojeniu
hh = linspace(dane.h2_min - 1, dane.h2_max + 1, 500);
f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 400]);
subplot(1,2,1); plot(hh, wartosci_mf(hh, mf_start)); grid on;
xlabel('h_2 [cm]'); ylabel('\mu'); title(['Przed: ' zad2_opis_mf(mf_start)]);
subplot(1,2,2); plot(hh, wartosci_mf(hh, mf_stroj)); grid on;
xlabel('h_2 [cm]'); ylabel('\mu'); title(['Po: ' zad2_opis_mf(mf_stroj)]);
exportgraphics(f, fullfile(folder, 'zad2_dmc_funkcje_przynaleznosci.png'), 'Resolution', 300); close(f);

rysuj_regulacje(t, ip, yzad, FDc, {h2_k, h2_r, h2_s}, {u_k, u_r, u_s}, ...
    {sprintf('DMC konwencjonalny (E = %.1f)', E_k), ...
     sprintf('DMC rozmyty przed strojeniem (E = %.1f)', E_r), ...
     sprintf('DMC rozmyty po strojeniu (E = %.1f)', E_s)}, ...
    fullfile(folder, 'zad2_dmc_przed_po_strojeniu.png'));


%% Funkcje
% Projekt rozmytego DMC: dla każdego zbioru model lokalny, jego odpowiedź skokowa
% i współczynniki analitycznego regulatora DMC (ke, ku)
function reg = projekt_dmc(p, mf, nast, Ts, co, d)
    L = numel(mf.c);
    reg.mf = mf;
    reg.D  = nast.D;
    reg.KE = zeros(L, 1);
    reg.KU = zeros(L, nast.D-1);
    reg.s  = zeros(nast.D, L);
    for i = 1:L
        m = zad2_lokalny(p, mf.c(i));
        s = odp_skokowa(m, nast.D, Ts, co, d);
        [reg.KE(i), reg.KU(i,:)] = dmc_ke_ku(s, nast.N, nast.Nu, nast.D, nast.lambda);
        reg.s(:,i) = s;
    end
end

% Odpowiedź skokowa modelu lokalnego próbkowana z okresem regulatora (jak w zadaniu 1)
function s = odp_skokowa(m, D, Ts, co, d)
    s  = zeros(D,1);
    dh1 = 0; dh2 = 0;
    for k = 1:D*co
        du  = double(k > d);
        dh1_new = dh1 + Ts*( m.d11*dh1 + m.b*du );
        dh2 = dh2 + Ts*( m.d21*dh1 + m.d22*dh2 );
        dh1 = dh1_new;
        if mod(k, co) == 0
            s(k/co) = dh2;
        end
    end
end

% Współczynniki regulatora DMC w wersji analitycznej
function [ke, ku] = dmc_ke_ku(s, N, Nu, D, lambda)
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

% Symulacja układu regulacji rozmytym DMC z obiektem nieliniowym. Przyrost sterowania
% to suma przyrostów regulatorów lokalnych ważona funkcjami przynależności zmierzonego h2.
% E - suma kwadratów błędów w chwilach próbkowania regulatora
function [h2, u, E] = sym_dmc_rozmyty(reg, sc)
    p = sc.p;
    n = numel(sc.yzad) - 1;
    h1 = sc.h1_0*ones(1,n+1);  h2 = sc.h2_0*ones(1,n+1);
    u = sc.F1_0*ones(1,n+1);
    u_akt = sc.F1_0;
    dUp = zeros(reg.D-1,1);

    for k = 1:n
        if mod(k-1, sc.co) == 0
            w  = zad2_przynaleznosc(h2(k), reg.mf);
            du = w*(reg.KE*(sc.yzad(k) - h2(k)) - reg.KU*dUp);
            % Ograniczenia przyrostu i wartości sterowania
            du = min(max(du, sc.limit_du(1)), sc.limit_du(2));
            u_nowe = min(max(u_akt + du, sc.limit_u(1)), sc.limit_u(2));
            dUp = [u_nowe - u_akt; dUp(1:end-1)];
            u_akt = u_nowe;
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

% Wartości funkcji przynależności dla wektora h2 (wiersz = punkt, kolumna = zbiór)
function W = wartosci_mf(hh, mf)
    W = zeros(numel(hh), numel(mf.c));
    for k = 1:numel(hh)
        W(k,:) = zad2_przynaleznosc(hh(k), mf);
    end
end

% Wykres przebiegów regulacji: wyjście z wartością zadaną i sterowanie z zakłóceniem
function rysuj_regulacje(t, ip, yzad, FD, h2, u, nazwy, plik)
    f = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 700]);
    subplot(2,1,1); hold on; box on;
    for i = 1:numel(h2)
        stairs(t(ip), h2{i}(ip));
    end
    stairs(t(ip), yzad(ip), 'k--'); hold off;
    grid on; ylabel('h_2 [cm]');
    legend([nazwy, {'wartość zadana'}], 'Location', 'best');

    subplot(2,1,2); hold on; box on;
    for i = 1:numel(u)
        stairs(t(ip), u{i}(ip));
    end
    stairs(t(ip), FD(ip), 'k--'); hold off;
    grid on; ylabel('F [cm^3/s]'); xlabel('t [s]');
    legend([nazwy, {'zakłócenie F_D'}], 'Location', 'best');
    exportgraphics(f, plik, 'Resolution', 300); close(f);
end
