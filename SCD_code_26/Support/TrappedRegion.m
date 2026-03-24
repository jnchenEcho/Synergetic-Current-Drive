
function [thetaC, thetaT] = TrappedRegion(theta_p, IAR)

B0 = 1;
e = IAR;

B = B0/(1+e*cos(theta_p));
Bm = B0/(1-e);
yita = Bm/B;

theta_c = asin(sqrt(1/yita));            % radian, for circulating particle
thetaC = rad2deg(theta_c);             % the degree of angle, for circulating particle

theta_t = pi/2 - theta_c;       % radian, for trapped particle
thetaT = rad2deg(theta_t);
end



