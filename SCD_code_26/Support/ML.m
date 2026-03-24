function Res = ML(L, IAR, k, N)
if L > 1
    Res = 0;
    return
end

e = IAR;
h1 = 1 - e;
h2 = (1 - e)^2/( 1 - e^2)^(1/2);
h3 = 1/pi * ( (1+e)*asin(sqrt(2*e/(1+e))) + sqrt(2*e*(1 + e)) );
c2 =  - 1/4 * (h2 - h1^2) / ((1-h1)^2);
c4 = h3^2/(1 - h1) - 1;
s = L*(1-h1)/(1-L*h1);

if s == 1
    Res = 0;
    return
end

fun = @(z)  ( k*(1-h1) + 3*h1*z )*( (1 - h1 + z*h1)^(5/2) * (1 + z^2*( c2 * (1-z^2) + c4 * z^2))^(1/2) )^(-1);
IntRes = GaussInt(fun, s, 1, N);
Res = 0.5*sqrt(1-h1)/( k - (k-3)*L*h1 )*IntRes;
end