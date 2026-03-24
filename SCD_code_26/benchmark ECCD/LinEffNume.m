function Res = LinEffNume(y, ga, gb, Te, n_para, eLL, Zeff, fc, IAR, h, c, N, ue)
fun = @(g) uVer2(y, g, n_para, c)^eLL * FM(g, Te) * DiffChi0(y, g, n_para, Zeff , fc, IAR, h, c, N, ue);
IntRes = GaussInt(fun, ga, gb, N);
Res = IntRes;
end
