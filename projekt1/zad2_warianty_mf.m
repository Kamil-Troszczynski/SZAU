% Warianty kształtu funkcji przynależności do strojenia: dla każdego rozkładu
% centrów g funkcje trójkątne oraz funkcje Gaussa o szerokościach z szer_lista
function war = zad2_warianty_mf(L, g_lista, szer_lista, h2_min, h2_max)
    war = cell(1, numel(g_lista)*(1 + numel(szer_lista)));
    iw = 0;
    for g = g_lista
        iw = iw + 1;
        war{iw} = zad2_mf(L, 'trojkatne', g, 0, h2_min, h2_max);
        for sz = szer_lista
            iw = iw + 1;
            war{iw} = zad2_mf(L, 'gauss', g, sz, h2_min, h2_max);
        end
    end
end
