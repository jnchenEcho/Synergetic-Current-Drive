function f = dLH(up, up1, up2, dlh)
% 目前这个还是有错误，因为没有导出dLH的正确形式。     20230729

% Cofficience
dup1 = 0.5;
dup2 = 1;

upl = up1 - dup1^2/(2*up1);
upr = up2 - dup2^2/(2*up2);


Al = exp(dup1^2/(4*up1^2));
Ar = up1/up2*exp(dup2^2/(4*up2^2));

% Al = 0;
% Ar = 0;

ue = 1;

if  up < 0
    f = 0;
    return
elseif up >= 0 && up < up1
    dLH1 = Al.*exp(1).^((-1).*dup1.^(-2).*(ue.^(-1).*up+(-1).*upl).^2);
    f = dLH1;
    return
elseif up >= up1 && up <= up2
    f = dlh;
    return
elseif up >= up2
    dLH3 = Ar.*exp(1).^((-1).*dup2.^(-2).*(ue.^(-1).*up+(-1).*upr).^2);
    f = dLH3;
    return
end
end