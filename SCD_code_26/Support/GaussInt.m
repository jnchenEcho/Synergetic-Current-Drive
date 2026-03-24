function Res = GaussInt(f, a, b, N)
[x, w]=lgwt(N,a,b);
sum = 0;
for i = 1 : N
    res = w(i) * f(x(i));
    sum = sum + res;
end
Res = sum;
end
