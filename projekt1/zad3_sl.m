% Rozmyty regulator DMC typu SL (Successive Linearization) w wersji numerycznej
% z ograniczeniami sterowania.
%
% W każdej chwili próbkowania regulatora:
%   1. model rozmyty jest linearyzowany w bieżącym punkcie pracy: M(k) = sum_i w_i*M_i, Mp(k) = sum_i w_i*Mp_i,
%   2. wyznaczana jest odpowiedź swobodna Y0 = y(k) + Mp(k)*dUp,
%   3. rozwiązywane jest zadanie QP:
%         min ||Yzad - Y0 - M(k)*dU||^2 + lambda*||dU||^2
%         przy du_min <= dU <= du_max  oraz  u_min <= u(k-1) + J*dU <= u_max,
%   4. do obiektu trafia pierwszy element dU.
%
% Pole regulator.typ:
%   'analityczny' - rozmyty DMC analityczny (PDC) z zadania 2, ograniczenia przez przycięcie,
%   'sl'          - SL numeryczny, ograniczenia w zadaniu QP.

classdef zad3_sl
    methods (Static)
        function regulator = projekt(model, mf, nastawy, Ts, co, d, typ)
            regulator = zad2_rozmyte.projekt_dmc(model, mf, nastawy, Ts, co, d);
            regulator.typ = typ;
            regulator.N = nastawy.N;
            regulator.Nu = nastawy.Nu;
            regulator.lambda = nastawy.lambda;

            % Macierze M_i, Mp_i modeli lokalnych
            L = numel(mf.c);
            regulator.M  = zeros(nastawy.N, nastawy.Nu, L);
            regulator.Mp = zeros(nastawy.N, nastawy.D-1, L);
            for i = 1:L
                [regulator.M(:,:,i), regulator.Mp(:,:,i)] = zad1_dmc.macierze(regulator.s(:,i), nastawy.N, nastawy.Nu, nastawy.D);
            end

            if strcmp(typ, 'sl')
                regulator.opcje = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
            end
        end

        % Model rozmyty zlinearyzowany w punkcie pracy y
        function [M, Mp] = linearyzacja(regulator, y)
            w = regulator.wagi(y);
            M  = zeros(regulator.N, regulator.Nu);
            Mp = zeros(regulator.N, regulator.D-1);
            for i = 1:numel(w)
                M  = M  + w(i)*regulator.M(:,:,i);
                Mp = Mp + w(i)*regulator.Mp(:,:,i);
            end
        end

        % Plan przyrostów sterowania dU (Nu x 1); do obiektu trafia dU(1)
        function dU = plan(regulator, y, yzad, dUp, u_prev, limit_u, limit_du)
            if strcmp(regulator.typ, 'analityczny')
                dU = regulator.wagi(y)*(regulator.KE*(yzad - y) - regulator.KU*dUp);
                return;
            end

            [M, Mp] = zad3_sl.linearyzacja(regulator, y);
            Y0 = y + Mp*dUp;
            Yzad = yzad*ones(regulator.N, 1);
            H = 2*(M'*M + regulator.lambda*eye(regulator.Nu));
            f = -2*M'*(Yzad - Y0);

            dU = zad3_sl.qp(regulator, H, f, u_prev, limit_u, limit_du);
        end

        function dU = qp(regulator, H, f, u_prev, limit_u, limit_du)
            Nu = regulator.Nu;
            J = tril(ones(Nu));
            jedynki = ones(Nu, 1);
            A = [J; -J];                                                  % ograniczenia wartości u
            b = [(limit_u(2) - u_prev)*jedynki; (u_prev - limit_u(1))*jedynki];
            dU_min = limit_du(1)*jedynki;                                 % ograniczenia przyrostu
            dU_max = limit_du(2)*jedynki;
            [dU, ~, flaga] = quadprog(H, f, A, b, [], [], dU_min, dU_max, zeros(Nu, 1), regulator.opcje);
            if flaga <= 0
                error('quadprog nie znalazł rozwiązania (flaga %d)', flaga);
            end
        end

        % Symulacja układu regulacji z obiektem nieliniowym (jak zad1_dmc.symulacja).
        % czas - średni czas wyznaczenia sterowania [ms]
        function [h2, u, E, czas] = symulacja(regulator, scenariusz)
            n = numel(scenariusz.yzad) - 1;
            h1 = scenariusz.h1_0*ones(1, n+1);  h2 = scenariusz.h2_0*ones(1, n+1);
            u = scenariusz.F1_0*ones(1, n+1);
            u_akt = scenariusz.F1_0;
            dUp = zeros(regulator.D-1, 1);
            czasy = [];

            for k = 1:n
                if mod(k-1, scenariusz.co) == 0
                    tic;
                    dU = zad3_sl.plan(regulator, h2(k), scenariusz.yzad(k), dUp, u_akt, scenariusz.limit_u, scenariusz.limit_du);
                    czasy(end+1) = toc; %#ok<AGROW>
                    [u_akt, dUp] = zad1_dmc.ogranicz(dU(1), u_akt, dUp, scenariusz.limit_u, scenariusz.limit_du);
                end
                u(k) = u_akt;
                if k > scenariusz.d
                    F1 = u(k-scenariusz.d);
                else
                    F1 = scenariusz.F1_0;
                end
                [h1(k+1), h2(k+1)] = zad1_obiekt.krok(scenariusz.model, h1(k), h2(k), F1, scenariusz.FD(k), scenariusz.Ts);
            end
            u(end) = u(end-1);

            ip = 1:scenariusz.co:n+1;
            E = sum((scenariusz.yzad(ip) - h2(ip)).^2);
            czas = 1e3*mean(czasy);
        end

        % Wskaźniki jakości liczone w chwilach próbkowania regulatora:
        % TV - całkowita zmienność sterowania sum|du|, n_ogr - liczba kroków z sterowaniem na ograniczeniu
        function [TV, n_ogr] = wskazniki(u, scenariusz)
            u = u(1:scenariusz.co:end);
            du = diff(u);
            TV = sum(abs(du));
            tol = 1e-6;
            na_ogr_u  = abs(u(2:end) - scenariusz.limit_u(1)) < tol | abs(u(2:end) - scenariusz.limit_u(2)) < tol;
            na_ogr_du = abs(du - scenariusz.limit_du(1)) < tol | abs(du - scenariusz.limit_du(2)) < tol;
            n_ogr = sum(na_ogr_u | na_ogr_du);
        end

        % Największe przeregulowanie po skokach wartości zadanej, liczone w kierunku skoku
        function p = przeregulowanie(h2, yzad)
            skoki = [find(diff(yzad) ~= 0) + 1, numel(yzad) + 1];
            p = 0;
            for i = 1:numel(skoki) - 1
                kierunek = sign(yzad(skoki(i)) - yzad(skoki(i) - 1));
                odcinek = skoki(i):skoki(i+1) - 1;
                p = max(p, max(kierunek*(h2(odcinek) - yzad(odcinek))));
            end
        end
    end
end
