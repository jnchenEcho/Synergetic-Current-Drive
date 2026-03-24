function Res = HL(L, IAR, N)
if L > 1
    Res = 0;
    return
end

e = IAR;
h1 = 1 - e;
h2 = (1 - e)^2/( 1 - e^2)^(1/2);
h3 = 1/pi * ( (1+e)*asin(sqrt(2*e/(1+e))) + sqrt(2*e*(1 + e)) );
% h3 = 1/pi * ( (1+e)*(sin(sqrt(2*e/(1+e))))^(-1) + sqrt(2*e*(1 + e)) );
c2 =  - 1/4 * (h2 - h1^2) / ((1-h1)^2);
c4 = h3^2/(1 - h1) - 1;
s = L*(1-h1)/(1-L*h1);

if s == 1
    Res = 0;
    return
end

fun = @(z) ( (1 - h1 + z*h1)^(3/2) * (1 + z^2*( c2 * (1-z^2) + c4 * z^2))^(1/2) )^(-1);
IntRes = GaussInt(fun, s, 1, N);
ap = sqrt(1 - h1)/h1 * (1/sqrt(1 - (1-s)*h1) - 1);
C1 = 0.0;
Res = 1/2*sqrt(1 - h1)*IntRes + C1*ap;
end

