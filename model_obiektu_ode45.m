%% Parametry
p.A2 = 430; p.C1 = 0.85; p.a1 = 18; p.a2 = 24; p.tau = 70;
F1_0 = 55; FD_0 = 10;

%% Punkt pracy
h1_0 = ((F1_0 + FD_0)/p.a1)^2;
h2_0 = ((F1_0 + FD_0)/p.a2)^2;

%% Wejścia: skok F1in w t=100 s, skok zakłócenia w t=1500 s
F1in = @(t) F1_0 + 5*(t >= 100);
FD   = @(t) FD_0 + 5*(t >= 1500);
F1   = @(t) F1in(t - p.tau);

%% Rozwiązanie
opts = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep',1);
[t, h] = ode45(@(t,h) zbiorniki(t,h,p,F1,FD), [0 3000], [h1_0; h2_0], opts);

%% Wykresy
subplot(2,1,1); plot(t, h(:,1)); grid on; ylabel('h_1 [cm]');
subplot(2,1,2); plot(t, h(:,2)); grid on; ylabel('h_2 [cm]'); xlabel('t [s]');

%% Model
function dh = zbiorniki(t, h, p, F1, FD)
h1 = max(h(1), 1e-6);
h2 = max(h(2), 0);
F2 = p.a1*sqrt(h1);
F3 = p.a2*sqrt(h2);
dh = [ (F1(t) + FD(t) - F2) / (2*p.C1*h1) ;
    (F2 - F3) / p.A2 ];
end