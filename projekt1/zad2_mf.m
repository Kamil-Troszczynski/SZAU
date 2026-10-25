% Funkcje przynależności dla zmiennej harmonogramującej h2
%   L    - liczba zbiorów rozmytych (modeli lokalnych)
%   typ  - 'trojkatne' lub 'gauss'
%   g    - rozkład centrów: g = 1 podział równomierny, g > 1 centra gęściej przy małych h2
%   szer - szerokość funkcji Gaussa względem odstępu między sąsiednimi centrami
function mf = zad2_mf(L, typ, g, szer, h2_min, h2_max)
    x = linspace(0, 1, L);
    mf.L    = L;
    mf.typ  = typ;
    mf.g    = g;
    mf.szer = szer;
    mf.c    = h2_min + (h2_max - h2_min)*x.^g;
end
