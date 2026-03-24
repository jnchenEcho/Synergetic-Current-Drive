function Res = EffDeno(y, ga, gb, Te, n_para, eLL, c, N)
fun = @(g) uVer2(y, g, n_para, c)^eLL * FM(g, Te) ;
IntRes = GaussInt(fun, ga, gb, N);
Res = IntRes;
end
