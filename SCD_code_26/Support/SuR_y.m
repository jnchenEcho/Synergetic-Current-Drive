%%
function f = SuR_y(u, u_para, Zeff,fc, fn, h1,  c, Dlh0, dlh, up1, up2, Ck, k, N,ue)
C3 = (k-1) * Dlh0 * dlh * fc/fn * Ck;
rho = (Zeff+1)/fn;
if u_para >= 0 && u_para < up1
%     A1 = 4.6327;
%     x = (u/c_ue)^2;
%     f = u^(-rho)*exp(rho*x*(1/4 - 3/32*x^2) + A1) ;
    f = 0;
    return
elseif u_para >= up1 && u_para <= up2
    fun = @(y) y^(rho + k)/( Gamma_p(y, u, c) + 1)^rho /(Gamma_p(y, u, c))^2 ;
    IntRes = GaussInt(fun, 0, 1, N);
    f = C3 * ( Gamma(u, c) + 1 )^(rho)* (u/ue)^(k+1) * IntRes;
    return
elseif u_para > up2
    f = 0;
    return
end
end

%%
function f = Gamma_p(y, u, c)
f = sqrt(1 + (y*u/c)^2);
end
%%
% function f = Fu_y(y, u, Zeff, fc, c_ue, N)
% C0 = (Zeff+1)/fc;
% fun = @(x) x^(C0+3)/(1+(u*x*y/c_ue)^2)^(3/2) * ( (1+sqrt(1+(u*y/c_ue)^2)) / (1+sqrt(1+(u*x*y/c_ue)^2)) )^(C0);
% IntRes = GaussInt(fun, 0, 1, N);
% f = 1/fc * u^4 * y^4 * IntRes;
% end

%%
% function f = SuR_y(u, u_para, Zeff,fc, fn, h1,  c_ue, Dlh0, dlh, up1, up2, N)
% C3 = 3 * Dlh0 * dlh * h1* fc/fn/(Zeff + 1 + 4*fc);
% rho = (Zeff+1)/fn;
% if u_para >= 0 && u_para < up1
%     A1 = 4.6327;
%     x = (u/c_ue)^2;
%     f = u^(-rho)*exp(rho*x*(1/4 - 3/32*x^2) + A1) ;
%     f = 0;
%     return
% elseif u_para >= up1 && u_para <= up2
%     n = 3.12;
%     f = u^(3.12);
%     return
% elseif u_para > up2
%     f = 0;
%     return
% end
% end



