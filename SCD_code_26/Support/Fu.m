function Res = Fu(u, Zeff, fc, c, N, ue)
rho = (Zeff + 1)/fc;
fun = @(x) x^(rho + 3) * (1 + (u*x/c)^2)^(-3/2) * ( (1+sqrt(1+(u/c)^2)) / (1+sqrt(1+(u*x/c)^2)) )^rho;
IntRes = GaussInt(fun, 0, 1, N);
Res = (u/ue)^4/fc * IntRes;
end
