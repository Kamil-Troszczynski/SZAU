% Rozmyty regulator DMC typu NPL (Nonlinear Prediction and Linearization) w wersji
% numerycznej z ograniczeniami sterowania.
%
% W każdej chwili próbkowania regulatora:
%   1. model rozmyty jest linearyzowany w bieżącym punkcie pracy: M(k) = sum_i w_i*M_i
%      (odpowiedź wymuszona - tak samo jak w SL),
%   2. odpowiedź swobodna Y0 jest liczona iteracyjnie rozmytym modelem TS obiektu
%      (zad2_rozmyte.sym_ts) przy zamrożonym sterowaniu u(k|k) = ... = u(k-1);
%      wagi modeli lokalnych są przeliczane w każdym kroku predykcji z przewidywanego h2,
%      uwzględniane są sterowania już wysłane do obiektu, a jeszcze "w drodze" (opóźnienie tau),
%   3. Y0 jest korygowana zakłóceniem typu DMC: d(k) = h2(k) - h2_model(k), gdzie model TS
%      jest symulowany równolegle z obiektem (h1 nie jest mierzone - bierze się je z modelu),
%   4. rozwiązywane jest zadanie QP jak w SL (zad3_sl.qp), do obiektu trafia pierwszy element dU.
%
% Zakłócenie FD jest traktowane jako niemierzalne - w modelu wewnętrznym FD = FD_0,
% a jego wpływ kompensuje d(k).
%
% Pole regulator.typ:
%   'analityczny' - DMC analityczny (konwencjonalny lub rozmyty PDC), przycięcie,
%   'sl'          - rozmyty DMC-SL numeryczny (QP),
%   'npl'         - rozmyty DMC-NPL numeryczny (QP).

classdef zad4_npl
    methods (Static)
        function regulator = projekt(model, mf, nastawy, Ts, co, d, typ, F1_0, FD_0)
            if strcmp(typ, 'analityczny')
                regulator = zad3_sl.projekt(model, mf, nastawy, Ts, co, d, typ);
                return;
            end
            regulator = zad3_sl.projekt(model, mf, nastawy, Ts, co, d, 'sl');   % M_i, Mp_i, opcje quadprog
            regulator.typ = typ;
            if ~strcmp(typ, 'npl'), return; end

            % Parametry modeli lokalnych modelu TS w postaci wektorów (szybki krok modelu)
            modele_lok = zad2_rozmyte.modele_lokalne(model, mf.c);
            regulator.P = struct('h1', [modele_lok.h1], 'h2', [modele_lok.h2], 'F', [modele_lok.F], ...
                                 'd11', [modele_lok.d11], 'b', [modele_lok.b], 'd21', [modele_lok.d21], 'd22', [modele_lok.d22]);
            regulator.Ts = Ts;
            regulator.co = co;
            regulator.d = d;
            regulator.FD_model = FD_0;

            % Stan ustalony modelu TS dla F1_0, FD_0 - punkt startowy modelu wewnętrznego
            % (model TS nie musi mieć równowagi dokładnie w (h1_0, h2_0))
            x = [((F1_0 + FD_0)/model.a1)^2, ((F1_0 + FD_0)/model.a2)^2];
            for k = 1:round(5000/Ts)
                [x(1), x(2)] = zad4_npl.krok_ts(regulator, x(1), x(2), F1_0, FD_0);
            end
            regulator.x0_model = x;
        end

        % Jeden krok Eulera rozmytego modelu TS (to samo co w zad2_rozmyte.sym_ts)
        function [h1, h2] = krok_ts(regulator, h1, h2, F1, FD)
            P = regulator.P;
            w = zad2_rozmyte.przynaleznosc(h2, regulator.mf);
            f1 = w*(P.d11.*(h1 - P.h1) + P.b.*(F1 + FD - P.F))';
            f2 = w*(P.d21.*(h1 - P.h1) + P.d22.*(h2 - P.h2))';
            h1 = max(h1 + regulator.Ts*f1, 1e-6);
            h2 = max(h2 + regulator.Ts*f2, 0);
        end

        % Odpowiedź swobodna (bez korekty d(k)) z modelu TS startującego ze stanu xm w chwili k.
        % u - dotychczasowe sterowania (u(1:k-1) znane), u_prev - sterowanie zamrożone na horyzoncie
        function Y0 = odp_swobodna(regulator, xm, k, u, u_prev, F1_0)
            co = regulator.co;  d = regulator.d;
            h1 = xm(1);  h2 = xm(2);
            Y0 = zeros(regulator.N, 1);
            for j = k:k + regulator.N*co - 1
                jd = j - d;
                if jd < 1
                    F1 = F1_0;            % sterowanie sprzed początku symulacji
                elseif jd <= k - 1
                    F1 = u(jd);           % sterowanie już wysłane, jeszcze "w drodze"
                else
                    F1 = u_prev;          % przyszłe sterowanie zamrożone (dU = 0)
                end
                [h1, h2] = zad4_npl.krok_ts(regulator, h1, h2, F1, regulator.FD_model);
                if mod(j - k + 1, co) == 0
                    Y0((j - k + 1)/co) = h2;
                end
            end
        end

        % Plan przyrostów sterowania NPL
        function dU = plan_npl(regulator, y, yzad, xm, k, u, u_prev, F1_0, limit_u, limit_du)
            dk = y - xm(2);                                          % zakłócenie typu DMC
            Y0 = zad4_npl.odp_swobodna(regulator, xm, k, u, u_prev, F1_0) + dk;
            M  = zad3_sl.linearyzacja(regulator, y);                 % odpowiedź wymuszona jak w SL
            Yzad = yzad*ones(regulator.N, 1);
            H = 2*(M'*M + regulator.lambda*eye(regulator.Nu));
            f = -2*M'*(Yzad - Y0);
            dU = zad3_sl.qp(regulator, H, f, u_prev, limit_u, limit_du);
        end

        % Symulacja układu regulacji z obiektem nieliniowym dla dowolnego typu regulatora.
        % Dla NPL równolegle z obiektem symulowany jest model TS (to samo F1, FD = FD_model).
        % czas - średni czas wyznaczenia sterowania [ms]
        function [h2, u, E, czas] = symulacja(regulator, scenariusz)
            n = numel(scenariusz.yzad) - 1;
            h1 = scenariusz.h1_0*ones(1, n+1);  h2 = scenariusz.h2_0*ones(1, n+1);
            u = scenariusz.F1_0*ones(1, n+1);
            u_akt = scenariusz.F1_0;
            dUp = zeros(regulator.D-1, 1);
            npl = strcmp(regulator.typ, 'npl');
            if npl, xm = regulator.x0_model; end
            czasy = [];

            for k = 1:n
                if mod(k-1, scenariusz.co) == 0
                    tic;
                    if npl
                        dU = zad4_npl.plan_npl(regulator, h2(k), scenariusz.yzad(k), xm, k, u, u_akt, ...
                                               scenariusz.F1_0, scenariusz.limit_u, scenariusz.limit_du);
                    else
                        dU = zad3_sl.plan(regulator, h2(k), scenariusz.yzad(k), dUp, u_akt, scenariusz.limit_u, scenariusz.limit_du);
                    end
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
                if npl
                    [xm(1), xm(2)] = zad4_npl.krok_ts(regulator, xm(1), xm(2), F1, regulator.FD_model);
                end
            end
            u(end) = u(end-1);

            ip = 1:scenariusz.co:n+1;
            E = sum((scenariusz.yzad(ip) - h2(ip)).^2);
            czas = 1e3*mean(czasy);
        end
    end
end
