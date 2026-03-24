
function Res = CouLog(g, ne, c_ue)
omegap = 5.641*10*sqrt(ne);     % s^(-1)
hbar = 6.5821195e-16;       % eV*s
mc2 = 0.510998e6;       % eV

u = c_ue*sqrt(g^2 - 1);
p = u/c_ue;
Lambda = 2^(1/4) * mc2 * ( 1 + p^2)^(3/4) / (hbar*omegap);
CouLog = log(Lambda);
Res = CouLog;
end
