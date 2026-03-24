
function res = DiffDeltaChi(y, gamma, n_para, Zeff, fc, fn, IAR, h, h1,  c, Dlh0, dlh, up1, up2, Ck, k, N,ue)
% h = B/Bmax;
up = c*(gamma - y)/n_para;
uv2 = uVer2(y, gamma, n_para, c);
uv = sqrt(uv2);
u = sqrt(up^2+uv^2);

Bmax_B = 1/h;

if u == 0
    res = 0;
    return
end

L = Bmax_B*uv2/u^2;

 res = sign(up) * ...
    (  ( gamma*ue^2/u * EL(L, IAR, k, N) *   DSuR_y(u, up, Zeff, fc, fn, h1, c, Dlh0, dlh, up1, up2, Ck, k,N, ue) ) ...            % зЂвт
     + ( 2 * Bmax_B * ue^2*up /u^3 * ( (gamma * up /u - n_para * u/c) * SuR_y(u, up, Zeff,fc, fn, h1, c, Dlh0, dlh, up1, up2, Ck,k,N,ue) * DEL(L, IAR, k) ) ) ...
     );       
  

end
