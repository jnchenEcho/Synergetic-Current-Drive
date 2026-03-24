function f = DdLH(up, up1, up2, dlh)
% 目前这个还是有错误，因为没有导出dLH的正确形式。     20230729
% Cofficience

dup1 = 0.5;
dup2 = 1;    % 这里会影响导数,2倍处为零。例如，dup2 = 1，up2 = 5,则导数在 up = 7 处为 0。

upl = up1 - dup1^2/(2*up1);
upr = up2 - dup2^2/(2*up2);


Al = exp(dup1^2/(4*up1^2));
Ar = up1/up2*exp(dup2^2/(4*up2^2));

% Al = 0;
% Ar = 0;


if up < 0 
    f = 0;
    return
elseif up >= 0 && up < up1
    DdLH1 = Al*exp(-(up - upl)^2/dup1^2)*(-2/dup1^2)*(up - upl);
    f = DdLH1;
    return
elseif  up >= up1 && up <= up2 
    f = 0; 
    return
elseif up > up2 
    DdLH3 =Ar*exp(-(up - upr)^2/dup2^2) * (-2/dup2^2) *(up - upr);
    f = DdLH3;
    return
end
end