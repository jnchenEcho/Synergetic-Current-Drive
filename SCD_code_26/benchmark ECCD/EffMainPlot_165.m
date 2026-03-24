% ECCD & SCD Result
%% 初始化
clc
clear
dbstop if error
format long
TimeStart = tic;
%% 参数设置
%20250402 N = 200, yn = 200, 10 min
%20250715 N = 300, yn = 250, 30 min

N = 300;     % 高斯积分点数
yn = 150;       % 计算点数

thetaP = 11*pi/12;      % 极向角
thetaP_deg = rad2deg(thetaP);
fprintf('极向角 = %3.1f \n', thetaP_deg);

Te = 2000;      % eV
eLL = 2;        % 谐波数
Zeff = 1.6;     
IAR = 0.2;      % inverse aspect ratio
ne = 2e19;      % m^(-3)
ve = 5.931*1e5*sqrt(Te);        % 最可几速度，m/s
c = 3e8;        % 光速
ue = 5.931*1e5*sqrt(Te);  
CouLog = 29.27 + 3/2*log(Te) - 1/2*log(ne);
mc2 = 0.510998e6;
tau = Te/mc2;
%% 低杂波设置
% Dlh = 0.2 20250402
% Dlh = 0.5, 20250715
Dlh0 = 0.2;     
up1 = 2*ue;
up2 = 4*ue;
dlh = 1;
%% Fu数值拟合设置
% KOrder = 4;     % Lin, Fu 非相对论
KOrder = 3.12;      % Fu 相对论数值拟合

%% 磁场设置
Bt0 = 1;
Bp0 = 1;
B0 = sqrt(Bt0^2 + Bp0^2);
Bt = Bt0/(1 + IAR*cos(thetaP));
Bp = Bp0/(1 + IAR*cos(thetaP));
B = B0/(1 + IAR*cos(thetaP));
Bmax = B0/(1 - IAR);
h = B/Bmax;
h1 = 1 - IAR;

%% 通行粒子
fcnum = fc(IAR, N);
fcpart = fcnum*100;
fprintf('通行粒子所占比例: %4.2f%% \n',fcpart);
fnnum = fn(IAR, KOrder, N);

%% Fu数值拟合设置
% Ck = 1/(Zeff+1+4*fcnum);  % Lin Fu 非相对论
Ck = 1;     % 相对论数值拟合

%% 捕获区
[~,TrappedAngle] = TrappedRegion(thetaP, IAR); 
k = tan(pi/2-deg2rad(TrappedAngle));
%% 数组设置
%  n_para_Mtx = [0.3 0.5 0.9];
 n_para_Mtx = [0.3 0.5 0.7];
% n_para_Mtx = [0.2 0.4 0.6];
% n_para_Mtx = [0.4];
[~, n_para_length] = size(n_para_Mtx);
% fprintf('共有 %d 条数据线,电子回旋波平行折射率分别为 %2.2f, %2.2f, %2.2f \n', n_para_length, n_para_Mtx(1), n_para_Mtx(2), n_para_Mtx(3));
LinEffNumeArray = zeros();
LinEffStarArray = zeros();      % ECCD
SCDEffNumeArray = zeros();
SCDEffStarArray = zeros();
EffDenoArray = zeros();
ya = zeros(n_para_length);
%%
for j = 1 : n_para_length
fprintf('正在计算第 %d / %d \n', j, n_para_length);
n_para = n_para_Mtx(j);
ymin = sqrt(1 - n_para^2);      % b^2 - 4ac > 0
yanum = 0.90;
ya(j) = yanum;      % 左区间
if ymin > ya(j)
    ymin = ceil(ymin*100)/100;
    ya(j) = ymin;
end
yb = 1;         % 右区间

dy = (yb - ya(j))/yn;
yy(j,:) = ya(j) : dy : yb;          % 共有 (yn+1) 个数值
i = 1;
for y = ya(j) : dy : yb
    % 确定积分上下限，即共振曲线左右两端
    syms g
    fung = g^2 - 1 - ((g - y)/n_para)^2;
    gRoot = solve(fung == 0);
    ga = double(gRoot(1));
    gb = double(gRoot(2));
    
    gammamin = ga;
    
    % 加热截断
    HeatingIndex = 50;      % 类似于Toray中的 etcutoff 
    gcut = gammamin + HeatingIndex * tau;
    
    % 捕获区与共振线的交界
    syms up 
    funh = (k^2+1)*up^2 - c^2*( (n_para*up/c + y)^2 - 1);
    hroot = solve(funh == 0);
    gi = hroot.*n_para./c + y;
    gc = double(gi(1));
    gd = double(gi(2));

    
    % 分母代表吸收功率，在求解gamma区域内全部积分
    % 分子代表驱动电流，在可加热区间内积分
    if isreal(gc) == 1 && isreal(gd) == 1    % 判定共振线是否穿过捕获区
        if gcut >= gd 
            gammamax = min([gcut gb]);
            LinEffNumeArray(1, i) = LinEffNume(y, gammamin, gc, Te, n_para, eLL, Zeff, fcnum, IAR, h, c, N,ue) + LinEffNume(y, gd, gammamax, Te, n_para, eLL, Zeff, fcnum, IAR, h, c, N,ue);
            SCDEffNumeArray(1, i) = SCDEffNumratorInt(y, gammamin, gc, Te, n_para, eLL, Zeff, fcnum, fnnum, IAR, h, h1, c, Dlh0, dlh, up1, up2, Ck, KOrder, N,ue) +  SCDEffNumratorInt(y, gd, gammamax, Te, n_para, eLL, Zeff, fcnum, fnnum, IAR, h, h1, c, Dlh0, dlh, up1, up2, Ck, KOrder, N,ue);      
            EffDenoArray(1, i) = EffDeno(y, gammamin, gammamax, Te, n_para, eLL, c, N);
        elseif gcut < gd && gcut >= gc
            gammamax = gcut;
            LinEffNumeArray(1, i) = LinEffNume(y, gammamin, gc, Te, n_para, eLL, Zeff, fcnum, IAR, h, c, N,ue) ;
            SCDEffNumeArray(1, i) = SCDEffNumratorInt(y, gammamin, gc, Te, n_para, eLL, Zeff, fcnum, fnnum, IAR, h, h1, c, Dlh0, dlh, up1, up2, Ck,KOrder, N,ue);
            EffDenoArray(1, i) = EffDeno(y, gammamin, gammamax, Te, n_para, eLL, c, N);
        elseif gcut < gc && gcut >= ga
            gammamax = gcut;
            LinEffNumeArray(1, i) = LinEffNume(y, gammamin, gammamax, Te, n_para, eLL, Zeff, fcnum, IAR, h, c, N,ue) ;
            SCDEffNumeArray(1, i) = SCDEffNumratorInt(y, gammamin, gammamax, Te, n_para, eLL, Zeff, fcnum, fnnum, IAR, h, h1, c, Dlh0, dlh, up1, up2, Ck,KOrder, N,ue);
            EffDenoArray(1, i) = EffDeno(y, gammamin, gammamax, Te, n_para, eLL, c, N);
        end
    else
        gammamax = min([gcut gb]);
        LinEffNumeArray(1, i) =  LinEffNume(y, gammamin, gammamax, Te, n_para, eLL, Zeff, fcnum, IAR, h, c, N, ue);
        SCDEffNumeArray(1, i) = SCDEffNumratorInt(y, gammamin, gammamax, Te, n_para, eLL, Zeff, fcnum, fnnum, IAR, h, h1, c, Dlh0, dlh, up1, up2, Ck, KOrder, N,ue);
        EffDenoArray(1, i) = EffDeno(y, gammamin, gammamax, Te, n_para, eLL, c, N);
    end

    LinEffStarArray(j, i) = - 4/CouLog * h1 * LinEffNumeArray(1, i) / EffDenoArray(1, i);
    SCDEffStarArray(j, i) = - 4/CouLog * h1 * SCDEffNumeArray(1, i) / EffDenoArray(1, i);
    i = i + 1;
end

end
fprintf('计算完成 \n')
%% 几何因子
GemotryFactor = -1;
LinEffArray = GemotryFactor .* LinEffStarArray;
SCDEffArray = GemotryFactor .* SCDEffStarArray;
%% 效率画图
figure

ymin = yanum;
ymax = yb;
dy = 0.02;
Emin = -0.1;
dE = 0.1;

ColorMatrix = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];         % 支持三种颜色
for k = 1 : n_para_length
    hold on
    plot(yy(k,:), LinEffArray(k,:), ...
        'Color', [ColorMatrix(k,:)],'LineWidth', 1.8, 'LineStyle', '-.','DisplayName', ['ECCD, $n_{\parallel}=$', num2str(n_para_Mtx(k))]);
    Cur(k) = plot(yy(k,:), SCDEffArray(k,:), ...
        'Color', [ColorMatrix(k,:)],'LineWidth', 1.8, 'LineStyle', '-', 'DisplayName', ['SCD, $n_{\parallel}=$', num2str(n_para_Mtx(k))]);
end
lgd = legend('Interpreter','latex', 'Location', 'best','FontSize',10,'FontWeight','bold');
plot([ymin ymax], [0 0], 'Color', 'black','HandleVisibility','off')

box on
xlabel('2 \omega_{c} / \omega');
ylabel('CD Efficiency');

if thetaP == pi/12
    Emax = 0.3;
    text(ymin+0.092, Emin + 0.03 ,'(a)','FontSize',14,'FontName','Times New Roman','FontWeight', 'bold');
elseif thetaP == pi/2
    Emax = 0.3;
    text(ymin+0.092, Emin + 0.03 ,'(b)','FontSize',14,'FontName','Times New Roman','FontWeight', 'bold');
elseif thetaP == 11*pi/12
    Emax = 0.3;
    text(ymin+0.092, Emin + 0.03,'(c)','FontSize',14,'FontName','Times New Roman','FontWeight', 'bold');
end


axis([ymin ymax Emin Emax]);
set(gca,'FontName','Times New Roman','FontSize',14, 'FontWeight', 'bold', ...
    'XTick', ymin : dy : ymax,'XMinorTick','off','YTick', Emin : dE : Emax,'YMinorTick','off', ...
    'TickDir','in','TickLength',[0.02 0.035]);
%%
saveas(gcf, 'CDefficiency_benchmarkECCD_thetaP165_DLH0.2_region2-4ue.fig');
exportgraphics(gcf, 'CDefficiency_benchmarkECCD_thetaP165_DLH0.2_region2-4ue.eps');

%%
TimeEnd = toc(TimeStart);
fprintf('运行时间为 %2.1f 分钟 \n', TimeEnd/60);

















