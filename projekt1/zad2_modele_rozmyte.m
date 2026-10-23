%%%% Modele rozmyte Takagi-Sugeno 2, 3, 4, 5 modeli lokalnych oraz porównanie z modelem liniowym i nieliniowym
clear; clc; close all;

%% Parametry obiektu
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

%% Parametry modeli rozmytych
L_lista = [2 3 4 5];        % ilość modeli
typ_mf  = 'trojkatne';      % ['trojkatne', 'gauss']
szer    = 0.5;              % szerokość gaussa

% Zmienna harmonogramująca: wyjście h2, podzielone równomiernie
% w zakresie wynikającym ze zmian F1 o +/-50% od punktu pracy
h2_min = ((0.5*F1_0 + FD_0)/p.a2)^2;
h2_max = ((1.5*F1_0 + FD_0)/p.a2)^2;

%% Dyskretyzacja
Ts   = 0.5;
Tk   = 3000;
N = round(Tk/Ts);
d    = round(p.tau/Ts);
t    = (0:N)*Ts;

%% Model liniowy = model TS z jednym modelem lokalnym w punkcie pracy
lin0 = lokalny(p, h2_0);

%% Funkcje przynależności
figure;
hh = linspace(h2_min - 1, h2_max + 1, 500);
for il = 1:numel(L_lista)
    c = linspace(h2_min, h2_max, L_lista(il));
    W = zeros(numel(hh), numel(c));
    for k = 1:numel(hh), W(k,:) = przynaleznosc(hh(k), c, typ_mf, szer); end
    subplot(2,2,il); plot(hh, W); grid on;
    xlabel('h_2 [cm]'); ylabel('\mu');
    title(sprintf('%d modele lokalne', L_lista(il)));
end

%% Eksperymenty: skoki F1 o różnej wielkości i kierunku
skoki = [-0.5 -0.25 0.25 0.5]*F1_0;
E = zeros(numel(skoki), numel(L_lista) + 1);   % błąd: [liniowy, TS2..TS5]
FD = FD_0*ones(1, N+1);

figure;
for is = 1:numel(skoki)
    F1in = F1_0*ones(1, N+1);
    F1in(t >= 100) = F1_0 + skoki(is);

    h2n = sym_nielin(p, F1in, FD, Ts, d, F1_0, h1_0, h2_0);
    h2l = sym_ts(lin0, h2_0, typ_mf, szer, F1in, FD, Ts, d, F1_0, h1_0, h2_0);
    E(is,1) = sum((h2l - h2n).^2);

    subplot(2,2,is);
    plot(t, h2n, 'k', t, h2l, '--'); hold on;
    for il = 1:numel(L_lista)
        c   = linspace(h2_min, h2_max, L_lista(il));
        lok = arrayfun(@(ci) lokalny(p, ci), c);
        h2f = sym_ts(lok, c, typ_mf, szer, F1in, FD, Ts, d, F1_0, h1_0, h2_0);
        E(is, il+1) = sum((h2f - h2n).^2);
        plot(t, h2f);
    end
    hold off; grid on; xlabel('t [s]'); ylabel('h_2 [cm]');
    title(sprintf('Skok F_1 o %+.1f cm^3/s', skoki(is)));
    if is == 1                         % jedna legenda wystarczy dla całej figury
        legend([{'nieliniowy', 'liniowy'}, ...
                arrayfun(@(L) sprintf('TS, %d modele', L), L_lista, 'UniformOutput', false)], ...
                'Location', 'best');
    end
end

%% Tabela błędów (suma kwadratów różnicy względem modelu nieliniowego)
fprintf('\nBłąd modeli względem modelu nieliniowego (suma kwadratów)\n');
fprintf('%10s %10s', 'skok F1', 'liniowy');
kol = arrayfun(@(L) sprintf('TS%d', L), L_lista, 'UniformOutput', false);
fprintf('%10s', kol{:});
fprintf('\n');
for is = 1:numel(skoki)
    fprintf('%10.1f', skoki(is)); fprintf('%10.2f', E(is,:)); fprintf('\n');
end

%% Funkcje
% Model lokalny: linearyzacja w punkcie równowagi o poziomie h2i
function m = lokalny(p, h2i)
    m.h2  = h2i;
    m.F   = p.a2*sqrt(h2i);             % całkowity dopływ F1 + FD w równowadze
    m.h1  = (m.F/p.a1)^2;
    m.d11 = -p.a1/(4*p.C1*m.h1^(3/2));
    m.b   =  1/(2*p.C1*m.h1);
    m.d21 =  p.a1/(2*p.A2*sqrt(m.h1));
    m.d22 = -p.a2/(2*p.A2*sqrt(h2i));
end

% Znormalizowane funkcje przynależności (centra c, zmienna y = h2)
function w = przynaleznosc(y, c, typ, szer)
    L = numel(c);
    if L == 1, w = 1; return; end
    y = min(max(y, c(1)), c(end));      % poza zakresem działa model skrajny
    switch typ
        case 'trojkatne'
            w = zeros(1, L);
            for i = 1:L
                if i > 1 && y >= c(i-1) && y <= c(i)
                    w(i) = (y - c(i-1))/(c(i) - c(i-1));
                end
                if i < L && y >= c(i) && y <= c(i+1)
                    w(i) = (c(i+1) - y)/(c(i+1) - c(i));
                end
            end
        case 'gauss'
            sig = szer*mean(diff(c));
            w = exp(-(y - c).^2/(2*sig^2));
    end
    w = w/sum(w);
end

% Model rozmyty TS: ważona suma pochodnych modeli lokalnych
function h2 = sym_ts(lok, c, typ, szer, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
    N = numel(F1in) - 1;
    h1 = h1_0; h2 = h2_0*ones(1, N+1);
    for k = 1:N
        if k > d, F1 = F1in(k-d); else, F1 = F1_0; end
        w = przynaleznosc(h2(k), c, typ, szer);
        f1 = 0; f2 = 0;
        for i = 1:numel(lok)
            m = lok(i);
            f1 = f1 + w(i)*( m.d11*(h1 - m.h1) + m.b*(F1 + FD(k) - m.F) );
            f2 = f2 + w(i)*( m.d21*(h1 - m.h1) + m.d22*(h2(k) - m.h2) );
        end
        h1 = h1 + Ts*f1;
        h2(k+1) = h2(k) + Ts*f2;
    end
end

% Model nieliniowy
function h2 = sym_nielin(p, F1in, FD, Ts, d, F1_0, h1_0, h2_0)
    N = numel(F1in) - 1;
    h1 = h1_0; h2 = h2_0*ones(1, N+1);
    for k = 1:N
        if k > d, F1 = F1in(k-d); else, F1 = F1_0; end
        F2 = p.a1*sqrt(h1); F3 = p.a2*sqrt(h2(k));
        h1 = max(h1 + Ts*(F1 + FD(k) - F2)/(2*p.C1*h1), 1e-6);
        h2(k+1) = max(h2(k) + Ts*(F2 - F3)/p.A2, 0);
    end
end