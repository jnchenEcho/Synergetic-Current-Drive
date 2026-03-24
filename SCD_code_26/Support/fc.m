function Res = fc(IAR, N)
e = IAR;
h2 = (1 - e)^2/( 1 - e^2)^(1/2);
fun = @(L) HL(L, IAR, N);
IntRes = GaussInt(fun, 0, 1, N);
Res = 3/2*h2*IntRes;
end
