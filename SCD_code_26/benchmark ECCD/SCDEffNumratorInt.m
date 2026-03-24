function f = SCDEffNumratorInt(y, a, b, Te, n_para, eLL, Zeff,  fc, fn, IAR, h, h1, c, Dlh0, dlh, up1, up2, Ck, k, N, ue) 
fun = @(gamma) uVer2(y, gamma, n_para, c)^eLL * FM(gamma, Te) ...
    * ( DiffDeltaChi(y, gamma, n_para, Zeff, fc, fn, IAR, h, h1,  c, Dlh0, dlh, up1, up2, Ck, k, N,ue) + DiffChi0(y, gamma, n_para, Zeff , fc, IAR, h, c, N,ue) ) ;

 IntRes = GaussInt(fun, a, b, N);
 f = IntRes;
end

