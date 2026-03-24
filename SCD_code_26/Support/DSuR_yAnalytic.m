function f = DSuR_yAnalytic(u, Zeff,  IAR,theta_p, fc, c_ue, Dlh0, up1, up2, N)
C0 = (Zeff+1)/fc;
fun = @(y) y^(C0) * DPartu(y, u, Zeff, fc, c_ue,  N)* Int_L( y, u, IAR, theta_p, up1, up2, N);
IntRes = GaussInt(fun, 0, 1, N);
h2 = (1-IAR)^2/(1-IAR^2)^(1/2);  
Bmax_B = (1+IAR*cos(theta_p))/(1-IAR);
f = - 3*h2*Bmax_B*Dlh0*IntRes;
end

function Res = DPartu(y, u, Zeff, fc, c_ue,  N)
C0  = (Zeff+1)/fc;
Res = Partu1(y, u, C0, c_ue) * DPartu2(y, u, Zeff, fc, c_ue, N) + DPartu1(y, u, C0, c_ue) * Partu2(y, u, Zeff, fc, c_ue, N);
end

function Res = Partu1(y, u, C0, c_ue)
Res = u * ( (Gamma(u, c_ue)+1)/(Gamma_p(y, u, c_ue)+1) )^(C0) * Gamma_p(y, u, c_ue)^(-2);
end

function Res = DPartu1(y, u, C0, c_ue)
c = c_ue;
p = C0;
Res = (-2).*c.^(-2).*u.^2.*y.^2.*(1+c.^(-2).*u.^2.*y.^2).^(-2).*((1+(1+ ...
  c.^(-2).*u.^2).^(1/2)).*(1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).^(-1)) ...
  .^p+(1+c.^(-2).*u.^2.*y.^2).^(-1).*((1+(1+c.^(-2).*u.^2).^(1/2)).* ...
  (1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).^(-1)).^p+p.*u.*(1+c.^(-2).* ...
  u.^2.*y.^2).^(-1).*((1+(1+c.^(-2).*u.^2).^(1/2)).*(1+(1+c.^(-2).* ...
  u.^2.*y.^2).^(1/2)).^(-1)).^((-1)+p).*((-1).*c.^(-2).*u.*(1+(1+ ...
  c.^(-2).*u.^2).^(1/2)).*y.^2.*(1+c.^(-2).*u.^2.*y.^2).^(-1/2).*(1+ ...
  (1+c.^(-2).*u.^2.*y.^2).^(1/2)).^(-2)+c.^(-2).*u.*(1+c.^(-2).* ...
  u.^2).^(-1/2).*(1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).^(-1));
end

function Res = Partu2(y, u, Zeff, fc, c_ue, N)
Res = Fu_y(y, u, Zeff, fc, c_ue, N);
end

function Res = DPartu2(y, u, Zeff, fc, c_ue, N)
C0 = (Zeff+1)/fc;
c = c_ue;
p = C0;
fun = @(x) (-3).*c.^(-2).*fc.^(-1).*u.^5.*x.^(5+p).*y.^6.*(1+c.^(-2).*u.^2.* ...
  x.^2.*y.^2).^(-5/2).*((1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).*(1+(1+ ...
  c.^(-2).*u.^2.*x.^2.*y.^2).^(1/2)).^(-1)).^p+4.*fc.^(-1).*u.^3.* ...
  x.^(3+p).*y.^4.*(1+c.^(-2).*u.^2.*x.^2.*y.^2).^(-3/2).*((1+(1+c.^( ...
  -2).*u.^2.*y.^2).^(1/2)).*(1+(1+c.^(-2).*u.^2.*x.^2.*y.^2).^(1/2)) ...
  .^(-1)).^p+fc.^(-1).*p.*u.^4.*x.^(3+p).*y.^4.*(1+c.^(-2).*u.^2.* ...
  x.^2.*y.^2).^(-3/2).*((1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).*(1+(1+ ...
  c.^(-2).*u.^2.*x.^2.*y.^2).^(1/2)).^(-1)).^((-1)+p).*((-1).*c.^( ...
  -2).*u.*x.^2.*y.^2.*(1+c.^(-2).*u.^2.*x.^2.*y.^2).^(-1/2).*(1+(1+ ...
  c.^(-2).*u.^2.*y.^2).^(1/2)).*(1+(1+c.^(-2).*u.^2.*x.^2.*y.^2).^( ...
  1/2)).^(-2)+c.^(-2).*u.*y.^2.*(1+c.^(-2).*u.^2.*y.^2).^(-1/2).*(1+ ...
  (1+c.^(-2).*u.^2.*x.^2.*y.^2).^(1/2)).^(-1));
IntRes = GaussInt(fun, 0, 1, N);
Res = IntRes;
end

function f = Fu_y(y, u, Zeff, fc, c_ue, N)
C0 = (Zeff+1)/fc;
fun = @(x) x^(C0+3)/(1+(u*x*y/c_ue)^2)^(3/2) * ( (1+sqrt(1+(u*y/c_ue)^2)) / (1+sqrt(1+(u*x*y/c_ue)^2)) )^(C0);
IntRes = GaussInt(fun, 0, 1, N);
f = 1/fc * u^4 * y^4 * IntRes;
end

function f = Gamma_p(y, u, c_ue)
f = sqrt(1 + (y*u/c_ue)^2);
end


function f = Int_L(y, u, IAR, theta_p, up1, up2, N)
fun = @(L) Int_L_Fun(L, y, u, IAR, theta_p, up1, up2);
IntRes = GaussInt(fun, 0, 1, N);
f = IntRes;
end

function f = Int_L_Fun(L, y, u, IAR, theta_p, up1, up2)
Bmax_B = (1+IAR*cos(theta_p))/(1-IAR);
h = 1/Bmax_B;
up = u*sqrt(1-h*L);

DHLnum = DHL(L, IAR);
DDHLnum = DDHL(L, IAR);

m1 = DdLH(up, up1, up2)  * DHLnum * up;
     
m2 = dLH(up, up1, up2) * (DHLnum + DDHLnum*(2*(L - Bmax_B))/y^2 );
     
 f = m1 + m2;
end