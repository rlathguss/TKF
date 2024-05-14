% "Capture Region of Tactical Missile Equipped with Semi-Active Laser Seeker Using Tobit Kalman Filter"
% "Look Angle Constrained Impact Angle Control Guidance Law for Homing Missiles with Bearings-Only Measurements"
% TKF Validated.
% Guidance Law : IACG with Bearings-Only Measurement - Sliding Mode Control
% ONLY VALID FOR STATIONARY TARGET CASE
% Measure Look angle ==> Estimate Look angle rate & Range rate ==> Calculate LOS angle rate.

% Fixed amf, varying k1 parametric study case
%function [Max_true_sigma, SatReg_E, SatReg_T, t_seg_E, t_seg_T, SatT_E, SatT_T, amf, Cov_P] = IACGBOM_TKF(k1)

% Fixed k1, varying amf parametric study case
%function [Max_true_sigma, SatReg_E, SatReg_T, SatT_E, SatT_T, amf, Cov_P] = IACGBOM_TKF(amf)

% Both k1 & amf parametric study case
%function [Max_true_sigma, SatReg_E, SatReg_T, SatT_E, SatT_T, amf, Cov_P] = IACGBOM_TKF(amf, k1)

clear all; close all% COMMENT WHEN PERFORMING PARAMETRIC STUDY

% Navigation Constant
N = 3;

% Missile init params
missile.y   =  0;%7000*sind(27);                       % y pos, m 
missile.x   =  0;%-7000*cosd(27);                      % x pos, m
% am0 = deg2rad(30);
missile.HDG = deg2rad(30); %am0; 
sigma_0 = missile.HDG;                              % Initial Look Angle
%sigma_0 = deg2rad(5);                              % Initial Look Angle

missile.V   =  250;%200;                                 % velocity magnitude, m/s 
missile.yV  =  missile.V * sin(missile.HDG);        % y velocity, m/s 
missile.xV  =  missile.V * cos(missile.HDG);        % x velocity, m/s 

accel_m = 0;                                        % Initial Missile Acceleration, m/s^2
%sigma_max = deg2rad(+45);                           % rad, Seeker's Look Angle Limit 

% Target init params
target.y   =  0;                                    % y pos, m
target.x   =  10000;%0;                                    % x pos, m
target.HDG = deg2rad(0);                            % heading, rad
target.yV  =  0;                                    % y velocity, m/s
target.xV  =  0;                                    % x velocity, m/s
target.V   = sqrt(target.yV^2 + target.xV^2);       % velocity magnitude, m/s

% =========================================================================
% sigma_d = deg2rad(7.5); % rad, Desired Look Angle / Commanded Look Angle
% =========================================================================

sat_upper = deg2rad(20);
sat_lower = -sat_upper;

accel_max = 100; % m/s, acceleration limit

amf = deg2rad(-160); %159 rad, Desired Impact Angle
%amf = deg2rad(amf);  % rad, COMMENT OFF WHEN PARAMETRIC STUDY

% ==================== USER CHOSEN INITIAL PARAMETERS =====================
Rf = 0.01;           % m
phi1 = 0.15;        % Unitless
k2 = deg2rad(10);   % rad

% ================ k1 plays the role of DESIRED LOOK ANGLE ================
% Originally, sigma_max(FOV limit) - 0.01 [deg]
k1 = deg2rad(9) %is the value when Saturated region actually starts appearing.
%k1 = deg2rad(50);     % rad
%k1 = sigma_max - deg2rad(0.01);            % Conservative Case: Pure IACGBOM
%k1 = deg2rad(k1);     % COMMENT OFF WHEN PARAMETRIC STUDY
% 8.6 is the value where Max(true_sigma) is closest to saturation.
% =========================================================================

% Sim params
Simtime  = 100;                        % max sim duration, seconds
dt = 0.01;                             % time-step size, seconds 
Niter = Simtime/dt;                    % max iteration num

% Pre-allocate logger
logger.t   = nan(1, Niter);    % elapsed time
logger.my  = nan(1, Niter);    % pursuer north pos
logger.mx  = nan(1, Niter);    % pursuer east pos
logger.myV = nan(1, Niter);    % pursuer north vel
logger.mxV = nan(1, Niter);    % pursuer east vel
logger.m.HDG = nan(1, Niter);    % pursuer heading angle
logger.M_HDG_dot = nan(1, Niter); % Pursuer heading angle derivative

logger.ty  = nan(1, Niter);    % target north y pos
logger.tx  = nan(1, Niter);    % target east x pos
logger.tyV = nan(1, Niter);    % target north y vel
logger.txV = nan(1, Niter);    % target east x vel

logger.z_sigma = nan(1, Niter); % Meausrement : Look angle
logger.v = nan(1,Niter);        % Conditioned Measurement Noise: Bernoulli distribution

logger.E_range = nan(1, Niter);
logger.true_range = nan(1, Niter);
logger.est_error = nan(1,Niter);
logger.lambda = nan(1, Niter); % pursuer/target bearing
logger.lambda_dot = nan(1, Niter);
logger.true_lambda = nan(1, Niter);
logger.true_lambda_dot = nan(1,Niter);
logger.sigma = nan(1, Niter);  % Look angle 
logger.true_sigma = nan(1, Niter);
logger.E_sigma = nan(1, Niter);
logger.S = nan(1,Niter);
logger.e1 = nan(1,Niter);
logger.e2 = nan(1,Niter);

logger.TKF_lookerror = nan(1, Niter);   % TKF Look Angle Error
logger.EKF_lookerror = nan(1, Niter);    % EKF Look Angle Error

logger.M_Accel = nan(1, Niter);
logger.M_Accel_x = nan(1, Niter);
logger.M_Accel_y = nan(1, Niter);

logger.normKalman = nan(1,Niter); % norm(Kalman)

logger.count = nan(1,Niter)';
logger.Ez_total = nan(1,Niter)';
logger.pl = nan(1,Niter)';
logger.ph = nan(1,Niter)';
logger.pus = nan(1,Niter)';
logger.bias = nan(1,Niter)';
logger.RR = nan(1,Niter)';
logger.HH = nan(1,Niter)';
logger.R = nan(1,Niter);
logger.H = nan(1,Niter);

logger.Ez_l = nan(1,Niter)';
logger.Ez_h = nan(1,Niter)';
logger.Ez_us = nan(1,Niter)';
logger.Ez_total_nobias = nan(1,Niter)';

logger.KK = nan(1,Niter);
logger.K_EKF = nan(1,Niter)';
logger.x_EKF = nan(1,Niter)';
logger.xpp = nan(1,Niter)';

logger.E_sigma_nobias = nan(1,Niter)';
logger.HH_nobias = nan(1, Niter)';
logger.RR_nobias = nan(1, Niter)';
logger.K_nobias = nan(1, Niter)';

logger.tkf_rhsterm = nan(1,Niter);
logger.ekf_rhsterm = nan(1,Niter);

logger.Cov_P = nan(1,Niter)';

%logger.nonG_v = nan(1,Niter);

%--------------------------------------------------------------------------
% Init sim
%--------------------------------------------------------------------------
RTP_last = [(target.x - missile.x);...          % delta y % RTP: Relative Target Position w.r.t Missile; 2x1 Vector
    (target.y - missile.y)];	                % delta x

true_sigma = sigma_0;                           % Initial true_sigma is validated as sigma_0

% Target pos
target.y = target.y + target.yV*dt;             % dt시간 뒤의 target.y pos 
target.x = target.x + target.xV*dt;             % dt시간 뒤의 target.x pos 

% Missile pos
missile.y = missile.y + missile.yV*dt;          % dt시간 뒤의 missile.y pos 
missile.x = missile.x + missile.xV*dt;          % dt시간 뒤의 missile.x pos 

%--------------------------------------------------------------------------
% Run sim
%--------------------------------------------------------------------------
for k = 1:Niter

% SELECT ONE FILTER AMONG THREE BELOW : TKF ||  EKF w/ SAT as input ||  EKF w/o SAT. as input
% ================================================== TOBIT KALMAN FILTER ====================================================

[z, v, count] = MeasM_TKF_StatT(sigma_0, true_sigma, sat_upper, sat_lower);       % MEASUREMENT MODEL

logger.count(k) = count;                                    % Save Bernoulli Distribution
count = logger.count;

[E_range, E_sigma, Cov_P, Ez_total, pl, ph, pus, bias_lambda, RR, HH, Ez_us, KK, K_EKF, x_EKF, xpp, x_nobias, HH_nobias, RR_nobias, Ez_total_nobias, K_nobias] = SysM_TKF_StatT(z, dt, missile, accel_m, count, sat_upper, sat_lower); % SYSTEM MODEL & GUIDANCE FILTER

    logger.Ez_total(k) = Ez_total;
    logger.pl(k) = pl;
    logger.ph(k) = ph;
    logger.pus(k) = pus;
    bias = bias_lambda * 0.0052;                % bias = bias_lambda * "Measurement Noise Deviation(sigma_v)" 0.0052
    logger.bias(k) = bias';
    logger.RR(k) = norm(RR);
    logger.HH(k) = norm(HH);
    logger.R(k) = 0.0052^2;                     % MODIFY: sigma_v
    logger.H(k) = 1;                            % MODIFY: norm(H)
    logger.Ez_us(k) = Ez_us;
    logger.Ez_l(k) = pl*sat_lower;
    logger.Ez_h(k) = ph*sat_upper;

    logger.Ez_total_nobias(k) = Ez_total_nobias;     % Modify sigma_v if necessary. 0.0052
    logger.KK(k) = KK(1);
    logger.K_EKF(k) = K_EKF(1);
    logger.x_EKF(k) = x_EKF;

    logger.E_sigma_nobias(k) = x_nobias(1);
    logger.RR_nobias(k) = norm(RR_nobias);
    logger.HH_nobias(k) = norm(HH_nobias);
    logger.K_nobias(k) = K_nobias(1);

    logger.tkf_rhsterm(k) = KK(1)*(z - Ez_total);
    logger.ekf_rhsterm(k) = K_EKF(1)*(z - xpp);

    logger.Cov_P(k) = Cov_P(1,1);

    if pus==0
        disp 'UNSAT Probability is 0. All the parameters become NaN from now on.'
        break;
    end


% ====================== EXTENDED KALMAN FILTER with SATURATED MEASUREMENTS CONSIDERED WHEN ESTIMATING ======================

% [z] = MeasM_TKF_EKF_zSAT_StatT(sigma_0, true_sigma);
% 
% [E_range, E_sigma, Cov_P, K_EKF, xpp] = SysM_TKF_EKF_zSAT_StatT(z, dt, missile, accel_m);
% 
% logger.K_EKF(k) = K_EKF;
% logger.ekf_rhsterm(k) = K_EKF*(z - xpp);

% ======================================== TOBIT KALMAN FILTER WITHOUT BIAS TERM ============================================

% [z, v, count] = MeasM_TKF_StatT(sigma_0, true_sigma);       % MEASUREMENT MODEL: It does not depend on the bias term.
% 
% logger.count(k) = count;                                    % Save Bernoulli Distribution
% count = logger.count;
% 
% [E_range, E_sigma, Cov_P, Ez_total, pl, ph, pus, bias_lambda, RR, HH, Ez_us, xpp, x_nobias, HH_nobias, RR_nobias, Ez_total_nobias, K_nobias]...
%             = SysM_TKF_StatT_withoutBIAS(z, dt, missile, accel_m, count);    % SYSTEM MODEL & GUIDANCE FILTER
% 
%     logger.Ez_total(k) = Ez_total;
%     logger.pl(k) = pl;
%     logger.ph(k) = ph;
%     logger.pus(k) = pus;
%     bias = bias_lambda * 0.0052;                % bias = bias_lambda * "Measurement Noise Deviation(sigma_v)" 0.0052
%     logger.bias(k) = bias';
%     logger.RR(k) = norm(RR);
%     logger.HH(k) = norm(HH);
%     logger.R(k) = 0.0052^2;                     % MODIFY: sigma_v
%     logger.H(k) = 1;                            % MODIFY: norm(H)
%     logger.Ez_us(k) = Ez_us;
%     logger.Ez_l(k) = pl*sat_lower;
%     logger.Ez_h(k) = ph*sat_upper;
% 
%     logger.Ez_total_nobias(k) = Ez_total_nobias;     % Modify sigma_v if necessary. 0.0052
%     %logger.KK(k) = KK(1);
%     %logger.K_EKF(k) = K_EKF(1);
%     %logger.x_EKF(k) = x_EKF;
% 
%     logger.E_sigma_nobias(k) = x_nobias(1);
%     logger.RR_nobias(k) = norm(RR_nobias);
%     logger.HH_nobias(k) = norm(HH_nobias);
%     logger.K_nobias(k) = K_nobias(1);
% 
%     %logger.tkf_rhsterm(k) = KK(1)*(z - Ez_total);
%     %logger.ekf_rhsterm(k) = K_EKF(1)*(z - xpp);
% 
%     if pus==0
%         disp 'UNSAT Probability is 0. All the parameters become NaN from now on.'
%         break;
%     end

% =================== EXTENDED KALMAN FILTER with SATURATED MEASUREMENTS "NOT" CONSIDERED WHEN ESTIMATING ===================
% ====================== If measured values are Saturated, Estimation is done only with PREDICTED STATE =====================

% [z] = MeasM_TKF_EKF_zSAT_NOT_StatT(sigma_0, true_sigma);
% 
% [E_range, E_sigma] = SysM_TKF_EKF_zSAT_NOT_StatT(z, dt, missile, accel_m);

% ===========================================================================================================================

    % True RTP(Range = LOSdistance)
    RTP_true = [(target.x - missile.x); (target.y - missile.y)];      
    true_range = norm(RTP_true);
    
    % Error btw true_range & E_range
    est_error = true_range - E_range;
    
    % Missile velocity, m/s
    Vm = sqrt(missile.yV^2 + missile.xV^2);
    
    % Target velocity, m/s
    %Vt = sqrt(target.yV^2 + target.xV^2);
    
    % Line-of-sight (LOS) angle, rad w.r.t Ref.
    true_lambda = atan2(RTP_true(2), RTP_true(1));
    lambda = missile.HDG - E_sigma;
    %lambda_nobias = missile.HDG - x_nobias(1);

    % LOS angle time derivatives (d/dt lambda), rad/s
    lambda_dot = - missile.V * sin(E_sigma) / E_range;
    true_lambda_dot = - missile.V * sin(missile.HDG - true_lambda) / true_range;
    %lambda_dot_nobias = - missile.V * sin(x_nobias(1)) / E_range;

    % True Look Angle, rad
    true_sigma = missile.HDG - true_lambda; % Look angle, rad
    
    % Missile Heading Angle time derivative (d/dt missile.HDG), rad/s
    % if k==1
    %     missile.HDG_dot = 0;
    % else
    %     missile.HDG_dot = (missile.HDG - missile.HDG_prev)/dt;
    % end

    % =======================================================================================
    e1 = lambda - amf;              % rad
    %e1 = true_lambda - amf;
    e2 = E_sigma;                   % rad
    %e2 = true_sigma;
    S = e2 - k1*(e1/sqrt(e1^2+phi1^2));     % rad
    f2 = (1+k1*phi1^2/((e1^2+phi1^2)^1.5))*abs(sin(E_sigma));
    %f2 = (1+k1*phi1^2/((e1^2+phi1^2)^1.5))*abs(sin(true_sigma));

    % IACCG, accel_m: Missile Acceleration
    accel_m = -(Vm*f2/Rf + k2)*Vm*tanh(0.01*S);            % Change tanh constant if necessary. Original: 10 but Oscillation occurs.
    
    % Acceleration Saturation Condition
    if ( abs(accel_m) > accel_max && accel_m < 0 )% m/s^2
        accel_m = -accel_max;
    elseif ( abs(accel_m) > accel_max && accel_m > 0 )
            accel_m = accel_max;

    else
            accel_m;
    end
    % =======================================================================================

    % Terminate sim at interception
    % if (abs(true_range) <= 10)||(lambda < amf) % interception 기준
    %     disp('True Range is less than 10m or LOS angle achieved desired impact angle')
    %     break;
    % end
    
    if abs(true_range) <= 3
    disp('True Range is less than 3m')
    %rad2deg(amf - logger.lambda(k-1))
        break;
    end

    if (norm(RTP_true)-norm(RTP_last))/dt > 0
    disp('Range is increasing. Look into the LOS angle. Engagement may have failed.')
        break;
    end

    % Update missile pos for time-step(dt) and apply current slip
    missile.y = missile.y + missile.yV*dt; %+ current.yV*dt;  
    missile.x = missile.x + missile.xV*dt; %+ current.xV*dt;  
    
    % Compute the y/x acceleration commands
    % In pure pro-nav accel commands are applied normal to missile's velocity vector
    M_accel_y = +accel_m * cos(missile.HDG);
    M_accel_x = -accel_m * sin(missile.HDG);
    
    % Update missile y/x velocities
    missile.yV = missile.yV + M_accel_y*dt;
    missile.xV = missile.xV + M_accel_x*dt;
    
    % Update missile heading
    missile.HDG_prev = missile.HDG;
    missile.HDG = atan2(missile.yV, missile.xV);    % renewed
    
    RTP_last = RTP_true;

    % Update target pos for time step
    % if
    %     target.y = target.y + target.yV*dt; % + current.yV*dt;
    %     target.x = target.x + target.xV*dt; % + current.xV*dt;
    % end

    % RTP_last = RTP;     % renewed
    
    %-------------------------------------
    % Log time-step data
    logger.t(k)   = k*dt;
    logger.my(k)  = missile.y;
    logger.mx(k)  = missile.x;
    logger.myV(k) = missile.yV;
    logger.mxV(k) = missile.xV;
    logger.m.HDG(k) = missile.HDG;
    %logger.m.HDG_dot(k) = missile.HDG_dot;

    logger.ty(k)  = target.y;
    logger.tx(k)  = target.x;
    logger.tyV(k) = target.yV;
    logger.txV(k) = target.xV;

    logger.E_range(k) = E_range;
    logger.true_range(k) = true_range;

    logger.est_error(k) = est_error;

    logger.lambda(k) = lambda;
    logger.lambda_dot(k) = lambda_dot;
    logger.true_lambda(k) = true_lambda;
    logger.true_lambda_dot(k) = true_lambda_dot;
    logger.E_sigma(k) = E_sigma;
    logger.true_sigma(k) = true_sigma;
    logger.z_sigma(k) = z;
    % logger.v(k) = v;
    logger.xpp(k) = xpp;
    logger.TKF_lookerror(k) = true_sigma - E_sigma;
    % logger.EKF_lookerror(k) = true_sigma - x_EKF;

    logger.M_Accel(k) = accel_m;
    logger.M_Accel_x(k) = M_accel_x;
    logger.M_Accel_y(k) = M_accel_y;
    logger.S(k) = S;
    logger.e1(k) = e1;
    logger.e2(k) = e2;

    % logger.normKalman(k) = norm(K);

%Cov_P
% 
end

% =================== Saturated Region Area Calculation ===================
index_E_upper = find(logger.E_sigma >= sat_upper);            % Upper Saturated index
index_E_lower = find(logger.E_sigma <= sat_lower);            % Lower Saturated index
t_seg_E_upper = logger.t(index_E_upper);
t_seg_E_lower = logger.t(index_E_lower);
E_sigma_seg_E_upper = logger.E_sigma(index_E_upper);                % Saturated Estimated look angles
E_sigma_seg_E_lower = logger.E_sigma(index_E_lower);


index_T_upper = find(logger.true_sigma >= sat_upper);  
index_T_lower = find(logger.true_sigma <= sat_lower);
t_seg_T_upper = logger.t(index_T_upper);
t_seg_T_lower = logger.t(index_T_lower);
true_sigma_seg_T_upper = logger.true_sigma(index_T_upper);          % Saturated True look angles
true_sigma_seg_T_lower = logger.true_sigma(index_T_lower);

SatReg_E_upper = trapz(t_seg_E_upper, rad2deg(E_sigma_seg_E_upper)) - ((max(t_seg_E_upper)-min(t_seg_E_upper))*rad2deg(sat_upper));   % Estimation Saturated Region Area
SatReg_E_lower = -trapz(t_seg_E_lower, rad2deg(E_sigma_seg_E_lower)) + ((max(t_seg_E_lower)-min(t_seg_E_lower))*rad2deg(sat_lower));
if isempty(SatReg_E_upper)       
    SatReg_E_upper = 0;
end
if isempty(SatReg_E_lower)
    SatReg_E_lower = 0;
end
SatReg_E = SatReg_E_upper + SatReg_E_lower;             % SatReg_E automatically becomes [](NaN) if UNSAT whole time.

SatReg_T_upper = trapz(t_seg_T_upper, rad2deg(true_sigma_seg_T_upper)) - ((max(t_seg_T_upper)-min(t_seg_T_upper))*rad2deg(sat_upper));% True Saturated Region Area
SatReg_T_lower = -trapz(t_seg_T_lower, rad2deg(true_sigma_seg_T_lower)) + ((max(t_seg_T_lower)-min(t_seg_T_lower))*rad2deg(sat_lower));
if isempty(SatReg_T_upper)
    SatReg_T_upper = 0;
end
if isempty(SatReg_T_lower)
    SatReg_T_lower = 0;
end
SatReg_T = SatReg_T_upper + SatReg_T_lower;             % SatReg_T automatically becomes [](NaN) if UNSAT whole time.

SatT_E = dt * (length(t_seg_E_upper) + length(t_seg_E_lower));                          % Estimation Saturated duration
SatT_T = dt * (length(t_seg_T_upper) + length(t_seg_T_lower));                          % True Saturated duration

Max_true_sigma = rad2deg(max(logger.true_sigma));       % deg, Achievable max true look angle 


%--------------------------------------------------------------------------
% Visualize results
%--------------------------------------------------------------------------
close all;

% Impact index
[M,I] = min(logger.E_range);

% Range
%-------------------------------------
figure;
plot(logger.t, logger.E_range); hold on;
plot(logger.t, logger.true_range,'linewidth',2); hold on;
scatter(k*dt, M, 'filled')
title(['IACCG(DPP&PPN), N = ' num2str(N) ])
legend('Estimated Range', 'True Range',...
    ['Intercept: r = ' num2str(E_range) ' m, t = ' num2str(k*dt) ' s'],...
    'Location','nw')
ylabel('Range (m)')
xlabel('Elapsed time (sec)')
set(gca, 'YDir', 'reverse')

% Error between Estimation and True Range
figure;
scatter(logger.t, logger.est_error,'filled'); hold on;
title(['IACCG(DPP&PPN), N = ' num2str(N) ])
legend('Range Error')
ylabel('Range Error btw Estimated and True Range [m])')
xlabel('Elapsed time (sec)')
% 



% Estimated and True Lambda (LOS angle)
% figure;
% plot(logger.t, rad2deg(logger.lambda), 'linewidth', 2); hold on;
% plot(logger.t, rad2deg(logger.true_lambda),'linewidth',2); hold on;
% title('Estimated and True LOS Angle')
% legend('Estimated LOS angle', 'True LOS angle')
% xlabel('Elapsed time(sec)')
% ylabel('LOS Angle(deg)')



% 
% % Lambda dot True and Estimated
% %--------------------------------------
% figure;
% plot(logger.t, rad2deg(logger.lambda_dot)); hold on;
% plot(logger.t, rad2deg(logger.true_lambda_dot),'linewidth',2); hold on;
% title(['IACCG(DPP&PPN), N = ' num2str(N) ])
% legend('Estimated Lambda dot', 'True Lambda dot')
% ylabel('Lambda dot(deg/s)')
% xlabel('Elapsed time (sec)')
% 



% Estimated, True and Measured Sigma(Look angle)
%-------------------------------------
figure;
plot(logger.t, rad2deg(logger.E_sigma),'linewidth', 2); hold on;
plot(logger.t, rad2deg(logger.true_sigma),'linewidth',2); hold on;
scatter(logger.t, rad2deg(logger.z_sigma), 5); hold on;
%plot(logger.t, rad2deg(logger.xpp)); hold on;
title('Estimated and True Look Angle')
legend('Estimated Look Angle', 'True Look Angle', 'Measured Look Angle')%, 'Predicted Look Angle')
xlabel('Elapsed time(sec)')
ylabel('deg')





% Missile Heading angle derivative
% -------------------------------------
figure;
plot(logger.t, rad2deg(logger.M_HDG_dot),'linewidth',1.5); hold on;
title('Missile Heading angle derivative')
legend('Missile HDG Derivative')
xlabel('Elapsed time(sec)')
ylabel('Missile HDG angle derivative (deg/s)')

% Heading
%-------------------------------------
figure;
plot(logger.t, rad2deg(logger.m.HDG)); hold on;
plot(logger.t, rad2deg(logger.lambda)); hold on;
title('Modified DPP + PPN')
legend('Missile heading', 'LOS angle', 'Location', 'nw' )
ylabel('Heading angle (deg)')
xlabel('Elapsed time (sec)')





% Position
%-------------------------------------
figure;
scatter(logger.mx, logger.my, 'filled'); hold on;
scatter(logger.tx, logger.ty, 'filled');
title('M&T Trajectory')
legend('Pursuer', 'Target', 'Location', 'southoutside',...
    'Orientation','horizontal')
xlabel('+x (m)')
ylabel('+y (m)')





% Acceleration
%--------------------------------------
% figure;
% plot(logger.t, logger.M_Accel, 'linewidth', 1.8); hold on;
% title('IACCG(DPP+PPN)')
% legend('Missile Acceleration')
% ylabel('Missile Acceleration (m/s^2)')
% xlabel('Elapsed time (sec)')
% 
% % Kalman Gain
% % figure;
% % %plot(logger.t, logger.Kalman1); hold on;
% % %plot(logger.t, logger.Kalman2); hold on;
% % plot(logger.t, logger.normKalman); hold on;
% % title(['Pure Proportional Navigation, N = ' num2str(N) ])
% % legend('norm(Kalman)')
% % xlabel('Elapsed time(sec)')
% 
% % Expectation Ez
% figure;
% plot(logger.t, rad2deg(logger.Ez_total), 'linewidth', 1.8); hold on;
% plot(logger.t, rad2deg(logger.Ez_l), 'linewidth', 1.8); hold on;
% plot(logger.t, rad2deg(logger.Ez_h), 'linewidth', 1.8); hold on;
% plot(logger.t, rad2deg(logger.Ez_us), 'linewidth', 1.8); hold on;
% title('IACCG TKF Expectation(Ez) Components');
% legend('Ez_{total}', 'Ez_l', 'Ez_h', 'Ez_{us}')
% xlabel('Elapsed time (sec)')
% ylabel('deg')
% % 
% Probabilities pl, ph, pus
% figure;
% plot(logger.t, logger.pl, 'linewidth', 1.8); hold on;
% plot(logger.t, logger.ph, 'linewidth', 1.8); hold on;
% plot(logger.t, logger.pus, 'linewidth', 1.8); hold on;
% title('Probabilities')
% legend('p_l', 'p_h', 'p_{us}')
% ylabel('Probability')
% xlabel('Elapsed time(sec)')

% Bias term
% figure;
% plot(logger.t, rad2deg(logger.bias), 'linewidth', 1.2);
% title('Bias Term')
% legend('Bias')
% ylabel('Bias Value')
% xlabel('Elapsed time(sec)')
% 
% % Ez with Bias VS. No bias
% figure;
% plot(logger.t, rad2deg(logger.Ez_total), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.Ez_total_nobias), 'linewidth', 1.2); hold on;
% title('Ez_{total} with Bias VS. No bias')
% legend('Ez_{total} with Bias', 'Ez_{total} without Bias')
% ylabel('deg')
% xlabel('Elapsed time(sec)')
% % 
% % H VS. HH VS. HH w/o bias
% figure; 
% plot(logger.t, logger.H, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.HH, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.HH_nobias, 'linewidth', 1.2); hold on;
% title('H_{EKF} VS. H_{TKF} VS. H_{TKF} w/o bias')
% legend('H_{EKF}', 'H_{TKF}', 'H_{TKF} w/o bias')
% ylabel('Norm of the matrices')
% xlabel('Elapsed time(sec)')
% 
% % R VS. RR VS. RR w/o bias
% figure;
% plot(logger.t, logger.R, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.RR, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.RR_nobias, 'linewidth', 1.2); hold on;
% title('R_{EKF} VS. R_{TKF} VS. R_{TKF} w/o bias')
% legend('R_{EKF}', 'R_{TKF}', 'R_{TKF} w/o bias')
% ylabel('Norm of the matrices')
% xlabel('Elapsed time(sec)')
% 
% % TKF Kalman gain VS. EKF Kalman gain
% figure;
% plot(logger.t, logger.KK, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.K_EKF, 'linewidth', 1.2); hold on;
% title('TKF K_k VS. EKF K_k')
% legend('TKF K_k', 'EKF K_k')
% ylabel('Norm of the Kalman gain')
% xlabel('Elaped time(sec)')

% TKF Kalman Gain VS. TKF_nobias Kalman Gain
% figure;
% plot(logger.t, logger.KK, 'linewidth', 1.8); hold on;
% plot(logger.t, logger.K_nobias, 'linewidth', 1.8); hold on;
% title('K_{TKF} VS. K_{TKF} w/o bias')
% legend('K_{TKF}', 'K_{TKF} w/o bias')
% ylabel('Norm of the Kalman gains')
% xlabel('Elapsed time(sec')

% 
% % TKF E_sigma VS. EKF E_sigma Comparison
% figure;
% plot(logger.t, rad2deg(logger.x_EKF), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.E_sigma), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.true_sigma), 'linewidth', 1.1); hold on;
% title('Look Angle | TKF Estimate VS. EKF Estimate')
% legend('EKF Estimate', 'TKF Estimate', 'True Look Angle')
% ylabel('deg')
% xlabel('Elapsed time(sec)')
% 
% Look Angle Estimation Error
% figure;
% plot(logger.t, rad2deg(logger.TKF_lookerror)); hold on;
% plot(logger.t, rad2deg(logger.EKF_lookerror)); hold on;
% title('Look Angle Error')
% legend('TKF Look Angle Error', 'EKF Look Angle Error')
% ylabel('deg')
% xlabel('Elapsed time(sec)')




% Look Angle Estimation Error & 3-sigma(P Mtx) bound
figure; grid on;
plot(logger.t, rad2deg(logger.TKF_lookerror), 'linewidth', 1.5); hold on;
plot(logger.t, 3*rad2deg(sqrt(logger.Cov_P)), 'linewidth', 1.5, 'color', 'r'); hold on;
plot(logger.t, -3*rad2deg(sqrt(logger.Cov_P)), 'linewidth', 1.5, 'color', 'r'); hold on;
title('Look Angle Error Trajectory')
legend('TKF Estimation Error', '3-sigma bound')
ylabel('deg')
xlabel('Time (s)')





% 
% % Hxp VS. Ez_total VS. Ez_total_nobias (in Estimation equation)
% figure;
% plot(logger.t, rad2deg(logger.xpp)); hold on;
% plot(logger.t, rad2deg(logger.Ez_total)); hold on;
% plot(logger.t, rad2deg(logger.Ez_total_nobias)); hold on;
% title('H*xp VS. Ez_{total} VS. Ez_{total} w/o bias')
% legend('H*xp','Ez_{total}','Ez_{total} w/o bias')
% ylabel('deg')
% xlabel('Elapsed time(sec')
% 
% % TKF and EKF Estimate equation RHS term comparison
% figure;
% plot(logger.t, rad2deg(logger.tkf_rhsterm)); hold on;
% plot(logger.t, rad2deg(logger.ekf_rhsterm)); hold on;
% title('TKF&EKF Equation RHS term Comparison')
% legend('TKF RHS', 'EKF RHS')
% ylabel('deg')
% xlabel('Elapsed time(sec)')

% figure;
% plot(logger.t, rad2deg(logger.E_sigma), 'linewidth', 2); hold on;
% plot(logger.t, rad2deg(logger.E_sigma_nobias),'linewidth', 2); hold on;
% plot(logger.t, rad2deg(logger.true_sigma), 'linewidth', 2); hold on;
% title('Look Angle Trajectory')
% legend('Estimated Look Angle', 'Estimated Look Angle w/o bias', 'True Look Angle')
% ylabel('deg')
% xlabel('Elapsed time (sec)')
% 
% % Sliding Surface
% %-------------------------------------
% figure;
% plot(logger.t, logger.S,'linewidth', 1.8,'Marker','v','MarkerIndices',1:500:length(logger.t)); hold on;
% title('Sliding Surface')
% legend('Sliding Surface Variable')
% ylabel('unit')
% xlabel('Elapsed time(sec)')
% 
% % Error Variables
% %-------------------------------------
% figure;
% plot(logger.t, rad2deg(logger.e1),'linewidth', 1.8,'Marker','v','MarkerIndices',1:500:length(logger.t)); hold on;
% plot(logger.t, rad2deg(logger.e2),'linewidth', 1.8,'Marker','v','MarkerIndices',1:500:length(logger.t)); hold on;
% title('Error Variables')
% legend('e1', 'e2')
% ylabel('deg')
% xlabel('Elapsed time (sec)')

% =====================================================================================