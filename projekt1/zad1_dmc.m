% Regulator DMC w wersji analitycznej z ograniczeniami (konwencjonalny) oraz rozmyty po odpowiedniej konfiguracji

classdef zad1_dmc
    methods (Static)
        % Odpowiedź skokowa modelu liniowego
        function s = odp_skokowa(model_zlinearyzowany, D, Ts, co, d)
            s  = zeros(D, 1);
            dh1 = 0; dh2 = 0;
            for k = 1:D*co
                du  = double(k > d);
                dh1_new = dh1 + Ts*( model_zlinearyzowany.d11*dh1 + model_zlinearyzowany.b*du );
                dh2 = dh2 + Ts*( model_zlinearyzowany.d21*dh1 + model_zlinearyzowany.d22*dh2 );
                dh1 = dh1_new;
                if mod(k, co) == 0
                    s(k/co) = dh2;
                end
            end
        end

        % Współczynniki ke, ku regulatora DMC w wersji analitycznej
        function [ke, ku] = wspolczynniki(s, N, Nu, D, lambda)
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

        % Projektowanie regulatora dla modeli lokalnych modele_lok. Dla jednego modelu wagi = 1;
        % dla wielu modeli pole regulator.wagi trzeba ustawić na funkcję przynależności.
        function regulator = projekt(modele_lok, nastawy, Ts, co, d)
            L = numel(modele_lok);
            regulator.D = nastawy.D;
            regulator.KE = zeros(L, 1);
            regulator.KU = zeros(L, nastawy.D-1);
            regulator.s = zeros(nastawy.D, L);
            regulator.wagi = @(y) 1;
            for i = 1:L
                regulator.s(:,i) = zad1_dmc.odp_skokowa(modele_lok(i), nastawy.D, Ts, co, d);
                [regulator.KE(i), regulator.KU(i,:)] = zad1_dmc.wspolczynniki(regulator.s(:,i), nastawy.N, nastawy.Nu, nastawy.D, nastawy.lambda);
            end
        end

        function [u, dUp] = ogranicz(du, u_prev, dUp, limit_u, limit_du)
            du  = min(max(du, limit_du(1)), limit_du(2));
            u   = min(max(u_prev + du, limit_u(1)), limit_u(2));
            dUp = [u - u_prev; dUp(1:end-1)];
        end

        function scenariusz = utworz_scenariusz(model, yzad, FD, Ts, co, d, F1_0, h1_0, h2_0, limit_u, limit_du)
            scenariusz = struct('model', model, 'yzad', yzad, 'FD', FD, 'Ts', Ts, 'co', co, 'd', d, 'F1_0', F1_0, 'h1_0', h1_0, 'h2_0', h2_0, 'limit_u', limit_u, 'limit_du', limit_du);
        end

        % Symulacja układu regulacji z obiektem nieliniowym.
        % E - suma kwadratów błędów w chwilach próbkowania regulatora
        function [h2, u, E] = symulacja(regulator, scenariusz)
            n = numel(scenariusz.yzad) - 1;
            h1 = scenariusz.h1_0*ones(1, n+1);  h2 = scenariusz.h2_0*ones(1, n+1);
            u = scenariusz.F1_0*ones(1, n+1);
            u_akt = scenariusz.F1_0;
            dUp = zeros(regulator.D-1, 1);

            for k = 1:n
                if mod(k-1, scenariusz.co) == 0
                    w  = regulator.wagi(h2(k));
                    du = w*(regulator.KE*(scenariusz.yzad(k) - h2(k)) - regulator.KU*dUp);
                    [u_akt, dUp] = zad1_dmc.ogranicz(du, u_akt, dUp, scenariusz.limit_u, scenariusz.limit_du);
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
        end
    end
end
