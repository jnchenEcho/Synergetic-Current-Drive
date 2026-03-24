function Res = DEL(L, IAR, k)
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
Res = -0.5 * (1 - L*h1)^(-1/2) * (k - (k-3)*L*h1) * (1 + s^2*(c2*(1-s^2) + c4*s^2))^(-1/2);        % 这里做了改动，在第二项中加入了 lambda
end
