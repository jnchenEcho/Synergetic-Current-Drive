function Res = FM(g,Te)
mc2 = 0.510998e6;           % eV
w = (g - 1)*mc2;
Res = exp( - w/Te);
end
