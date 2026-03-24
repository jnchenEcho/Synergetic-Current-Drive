%%
function Res = DSuR_y(u, u_para, Zeff, fc, fn, h1, c, Dlh0, dlh, up1, up2, Ck, k, N,ue)


if u_para >= 0 && u_para < up1
    Res = 0;
    return
elseif u_para >= up1 && u_para <= up2
    rho = (Zeff+1)/fn;    % 公式中的rho
    C3 = (k-1) * Dlh0 * dlh * fc/fn * Ck;
    fun = @(y) (-1).*c.^(-2).*rho.*u.*(1+(1+c.^(-2).*u.^2).^(1/2)).^rho.*(u.* ...
  ue.^(-1)).^(1+k).*y.^(2+k+rho).*(1+c.^(-2).*u.^2.*y.^2).^(-3/2).*( ...
  1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).^((-1)+(-1).*rho)+(-2).*c.^(-2) ...
  .*u.*(1+(1+c.^(-2).*u.^2).^(1/2)).^rho.*(u.*ue.^(-1)).^(1+k).*y.^( ...
  2+k+rho).*(1+c.^(-2).*u.^2.*y.^2).^(-2).*(1+(1+c.^(-2).*u.^2.* ...
  y.^2).^(1/2)).^((-1).*rho)+c.^(-2).*rho.*u.*(1+c.^(-2).*u.^2).^( ...
  -1/2).*(1+(1+c.^(-2).*u.^2).^(1/2)).^((-1)+rho).*(u.*ue.^(-1)).^( ...
  1+k).*y.^(k+rho).*(1+c.^(-2).*u.^2.*y.^2).^(-1).*(1+(1+c.^(-2).* ...
  u.^2.*y.^2).^(1/2)).^((-1).*rho)+(1+k).*(1+(1+c.^(-2).*u.^2).^( ...
  1/2)).^rho.*(u.*ue.^(-1)).^k.*ue.^(-1).*y.^(k+rho).*(1+c.^(-2).* ...
  u.^2.*y.^2).^(-1).*(1+(1+c.^(-2).*u.^2.*y.^2).^(1/2)).^((-1).*rho);
    IntRes = GaussInt(fun, 0, 1, N);
    Res = C3*IntRes;
    return
elseif u_para > up2
    Res = 0;
    return

end

%%
% function Res = DSuR_y(u, u_para, Zeff, fc, fn, h1, c_ue, Dlh0, dlh, up1, up2, N)
% 
% 
% if u_para >= 0 && u_para < up1
%     Res = 0;
%     return
% elseif u_para >= up1 && u_para <= up2
%     C = (Zeff+1)/fn;    % 公式中的rho
%     C3 = 3 * Dlh0 * dlh * h1* fc/fn/(Zeff + 1 + 4*fc);
%     c = c_ue;
%     ue = 1; % 相当于归一化到 ue上
%     n = 3.12;
%     Res = n*u^(n-1);
%     return
% elseif u_para > up2
%     Res = 0;
%     return
% 
% end




