function Res = DiffChi0(y, g, n_para, Zeff , fc, IAR, h, c, N, ue)
% h = B/Bmax;
up = c*(g - y)/n_para;
uv2 = uVer2(y, g, n_para, c);
uv = sqrt(uv2);
u = sqrt(up^2 + uv^2);
Bmax_B = 1/h;

if u == 0
    Res = 0;
    return
end

L = Bmax_B * uv2 / u^2;

Res = sign(up) * ( g*ue^2/u  * DFu(u, Zeff, fc, c, N, ue) * HL(L, IAR, N) ...
    + 2 * Bmax_B * up*ue^2/u^3 * (g*up/u - n_para*u/c)*Fu(u, Zeff, fc, c, N, ue)*DHL(L, IAR) );

end



