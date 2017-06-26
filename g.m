function G = g(V)
    E = 1.17391304;
    Ga = -50/66*1e-3;
    Gb = -9/22*1e-3;
    Emax = 8.1818;
    Gc = 4.591e-3;
    if V<=-Emax
        G = Gc*V + Emax*(Gc-Gb) + E*(Gb-Ga);
    elseif V<=E
        G = Gb*V + (Gb-Ga)*E;
    elseif V<E
        G = Ga*V;
    elseif V<Emax
        G = Gb*V + (Ga-Gb)*E;
    else
        G = Gc*V + Emax*(Gb-Gc) + E*(Ga-Gb);
    end
end