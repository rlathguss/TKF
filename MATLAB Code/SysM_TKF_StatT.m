% "Capture Region of Tactical Missile Equipped with Semi-Active Laser Seeker Using Tobit Kalman Filter"
% TKF VALIDATED
% System Model is composed with Look Angle & Relative Range

function [E_range, E_sigma, Cov_P, Ez_total, pl, ph, pus, bias_lambda, RR, HH, Ez_us, KK, K_EKF, x_EKF, xpp, x_nobias, HH_nobias, RR_nobias, Ez_total_nobias, K_nobias] = SysM_TKF_StatT(z, dt, missile, accel_m, count, sat_upper, sat_lower) % Add 'count' if necessary.

persistent A H B Q R sigma_v sigma_w %sat_upper sat_lower 
persistent x P
persistent firstRun

if isempty(firstRun)
    firstRun = 1;
    %sat_upper = deg2rad(+7.5);  %7.5            % MODIFY IF NECESSARY.
    %sat_lower = -sat_upper;

    % m: # of Measurement components , n: # of State components
    % m = 1, n = 2
    % Output matrix mxn
    H = [ 1 0 ];
  
    x = [ deg2rad(15) 10000 ]';  % Initial guess 15 10000    % MODIFY WHEN INITIAL M&T GEOMETRY CHANGES.
    % 5, 7000
    %A = Ajacobian(x, missile, dt);

    % System(Process) noise nxn
    %Q = [0.1^2 0; 0 0.1^2];
    sigma_w = 0.1;
    %Q = sigma_w^2 * [dt^3/3 dt^2/2; dt^2/2 dt];
    %Q = [dt^3/3 dt^2/2; dt^2/2 dt];
    %Q = Qupdate(A, missile, dt, sigma_w);

    % Measurement noise mxm
    sigma_v = 0.0052;
    R = sigma_v^2;          % 0.0052[rad] = 0.298[deg]
    %R = 0.0052^2;

    % Error Covariance nxn
    P = [0.0052 0;
            0   50];

    B = [1/missile.V 0]';

end

%xp = A*x + B*u ;               % Current state prediction(k)c
A = Ajacobian(x, missile, dt);                  % Linearize

xp = fxp(x, dt, B, accel_m, missile, sigma_w);  % Prediction
xpp = xp(1);

Q = Qupdate(A, missile, dt, sigma_w);           % Process Noise Covariance Matrix

Pp = A*P*A' + Q;                % Current error covariance prediction(k) = A*(Previous error covariance)*A' + System noise


% ------------------ Ez, HH, RR CALCULATION -------------------
[Ez_total, HH, RR, pl, ph, pus, bias_lambda, Ez_us, HH_nobias, RR_nobias, Ez_total_nobias] = zExpectation(H, xp, sat_upper, sat_lower, sigma_v, count);
% -------------------------------------------------------------

K = Pp * HH' / (HH*Pp*HH' + RR);
KK = K;
% Original) K = Pp*H'/(H*Pp*H' + R);        % Current Kalman gain(k) = func(Current P prediction(k), H, R)
K_EKF = Pp*H'/(H*Pp*H'+R);

x = xp + K*(z - Ez_total);
%x = xp + K*(z - HH * xp);
% Original) x = xp + K*(z - H*xp);          % Current state estimation(k) = Current state prediction(k) + Current Kalman gain(k)*(Current Measurement(k) -H*Current state prediction(k))
x_EKF = xp + K_EKF*(z - H*xp);
x_EKF = x_EKF(1);                           % Look Angle 

K_nobias = Pp * HH_nobias' / (HH_nobias*Pp*HH_nobias' + RR_nobias);
x_nobias = xp + K_nobias*(z - Ez_total_nobias);
% Nobias Estimation: x = xp + K_EKF*(z - Ez_total_nobias);


% -------------------------------------------------------------
P = Pp - K*HH*Pp;                % Current error covariance calculation(k) = Current error covariance prediction(k) - K*H*Pp
Cov_P = P;
%==============================================================
E_sigma = x(1);             % Estimated Look Angle
E_range = x(2);             % Estimated Range
%==============================================================
function A = Ajacobian(x, missile, dt)

A = zeros(2, 2);
A(1,1) = missile.V*cos(x(1))/x(2);
A(1,2) = -missile.V*sin(x(1))/x(2)^2;
A(2,1) = missile.V*sin(x(1));
A(2,2) = 0;
A = eye(2) + A*dt;

%==============================================================
function Q = Qupdate(A, missile, dt, sigma_w)
Q = zeros(2,2);
Q(1,1) = ((A(1,1))^2*dt^2)/3 + dt*A(1,1) + 1;
Q(1,2) = (A(1,1)*A(2,1)*dt^2)/3 + dt*A(2,1)/2;
Q(2,1) = Q(1,2);
Q(2,2) = dt^2 * (A(2,1))^2 /3;
Q = (Q * sigma_w^2 * dt) / missile.V^2;
% %============================================================
function xp = fxp(x, dt, B, accel_m, missile, sigma_w)

xdot = zeros(2,1);
xdot(1) = missile.V*sin(x(1))/x(2);
xdot(2) = - missile.V*cos(x(1));
u = accel_m + randn*sigma_w;        % LATAX HAS ADDIVE UNCERTAINTY.

xp = x + xdot*dt + B*u*dt;
xp(1) = xp(1);%+randn*deg2rad(sigma_w);
%================================
function [Ez_total, HH, RR, pl, ph, pus, bias_lambda, Ez_us, HH_nobias, RR_nobias, Ez_total_nobias] = zExpectation(H, xp, sat_upper, sat_lower, sigma_v, count) % Add 'count' if necessary.

%count_nd = fitdist(count, 'normal');   % count normal distribution
           % OUTPUTs: count_nd.mu & count_nd.sigma 
           % IF SPECIFIC mu & sigma ARE GIVEN, THERE IS ONLY ONE NORMAL DISTRIBUTION CORRESPONDING TO IT.
           
          % count_nd.sigma = sigma_v;

alpha = (sat_lower - H*xp) / sigma_v;
beta  = (sat_upper - H*xp) / sigma_v;

%pdf_l = pdf(count_nd, alpha);                 % PDF value @ alpha   : phi(alpha)
%pdf_h = pdf(count_nd, beta);                  % PDF value @ beta    : phi(beta)
pdf_l = normpdf(alpha);
pdf_h = normpdf(beta);

%pl = cdf(count_nd, alpha);                    % CDF value [-inf, alpha] : Probability for lower Meas.Sat.
%ph = 1 - cdf(count_nd, beta);                 % CDF value [beta, +inf]  : Probability for upper Meas.Sat.
pl = cdf('normal', alpha, 0, 1);
ph = 1 - cdf('normal', beta, 0, 1);

if ph==1
    ph = 0.99999;
end

if pl==1
    pl = 0.99999;
end

%pus = 1 - pl - ph;                        % Probability for Meas.Unsat.
%pus = mean(count,"omitnan");
%pus = cdf('normal', beta, 0, 1) - cdf('normal', alpha, 0, 1);
pus = 1 - ph - pl;

bias_lambda = (pdf_h - pdf_l) / pus;             % BIAS TERM to avoid selection bias(thus biased estimate), Inverse Mill's Ratio

Ez_total = (H*xp - sigma_v*bias_lambda)*pus +...
                               sat_lower*pl +...
                               sat_upper*ph;        % Ez = Ez[UNSAT]*pus + Ez[SAT_l]*pl + Ez[SAT_h]*ph       

% Ez_total = H*xp*pus - sigma_v*(pdf_h-pdf_l) + pl*sat_lower + ph*sat_upper;

% Ez_total WITHOUT BIAS
Ez_total_nobias = H*xp*pus + pl*sat_lower + ph*sat_upper;


HH = pus*H;                                   % Modified H
RR = sigma_v^2*(1-bias_lambda^2)+...          % Modified R
     sigma_v*H*xp*bias_lambda+...
     sigma_v*(sat_lower*pdf_l - sat_upper*pdf_h)/pus;

HH_nobias = pus*H;
RR_nobias = sigma_v^2 + 0 + sigma_v*(sat_lower*pdf_l - sat_upper*pdf_h)/pus;

Ez_us = H*xp*pus - sigma_v*(pdf_h - pdf_l);   % E[z|l<z<h]*pus
%==================================











