% "Capture Region of Tactical Missile Equipped with Semi-Active Laser Seeker Using Tobit Kalman Filter"
% "Composite Guidance Scheme for Impact Angle Control Against a Nonmaneuvering Moving Target"
% TKF Validated.
% Guidance Law : IACCG(DPP+PPN)
% ONLY VALID FOR STATIONARY TARGET CASE
% Measure Look angle ==> Estimate Look angle rate & Range rate ==> Calculate LOS angle rate.

clear all; close all

% Navigation Constant
N = 3;
index = 0;

% Missile init params
missile.y   =  7000*sind(27);                       % y pos, m 
missile.x   =  -7000*cosd(27);                      % x pos, m  
missile.HDG = deg2rad(-22);
%sigma_0 = missile.HDG;                             % Initial Look Angle
sigma_0 = deg2rad(5);                               % Initial Look Angle

missile.V   =  200;                                 % velocity magnitude, m/s 
missile.yV  =  missile.V * sin(missile.HDG);        % y velocity, m/s 
missile.xV  =  missile.V * cos(missile.HDG);        % x velocity, m/s 

accel_m = 0;                                        % Initial Missile Acceleration, m/s^2
accel_max = 100;                                    % m/s, acceleration limit
sigma_max = deg2rad(+60);                           % rad, Seeker's Look Angle Limit 

% Target init params
target.y   =  0;                                    % y pos, m
target.x   =  0;                                    % x pos, m
target.HDG = deg2rad(0);                            % heading, rad
target.yV  =  0;                                    % y velocity, m/s
target.xV  =  0;                                    % x velocity, m/s
target.V   = sqrt(target.yV^2 + target.xV^2);       % velocity magnitude, m/s

% Other Initial Conditions and Parameters
eta = target.V / missile.V; 

% sigma_d Selection in order to bring out lambda_switch and also compare whether initial setup is good or not.
sigma_d_lower =  asin(eta); % rad, Lower bound for sigma_d
sigma_d_upper =  min(acos(2*eta), sigma_max); % rad, Upper bound for sigma_d

sigma_d = deg2rad(30); % rad, Desired Look Angle / Commanded Look Angle

sat_upper = deg2rad(7.5);
sat_lower = -sat_upper;

if not ((abs(sigma_d) > sigma_d_lower)&&(abs(sigma_d) <= sigma_d_upper))
    disp 'Check sigma_d initial condition. Must be within a range.'
    return
end

if sigma_0 > abs(sigma_d)
    disp 'sigma_0 must be equal or less than sigma_d'
    return
end

zeta1 = eta/sin(sigma_d);
zeta2 = eta/cos(sigma_d);


amf = deg2rad(-45); % rad, Desired Impact Angle
% labmda_switch and amf are functions of each other, thus if you decide amf, then lambda_switch can be calculated itself.
lambda_switch = (( (atan2(sin(amf)-eta*sin(target.HDG), cos(amf)-eta*cos(target.HDG)))) - ((amf) - (sigma_d)) / N ) * (N/(N-1));
% rad, lambda_switch is specific LOS angle where the guidance gain changes.

% K Gain Selection
K = 230; % Select K value as you wish, IF YOU CAN.
R0 = norm([(target.x - missile.x), (target.y - missile.y)]); % Initialized R0 in order to calculate K_upper  

K_upper = (accel_max * R0 + missile.V^2 * sin(sigma_0)) / (R0*(abs(sigma_d) - sigma_0));

if K > K_upper
    disp 'K must be less than K_upper'
    return
end

% ===================================================IN REPAIR=====================================================
% accel_max comparison thus determining whether abs(amf) =< abs(amf_max)
% if zeta1 < 1
%     Denom_int = (2/sqrt(1-zeta1^2))*(atan2((tan(lambda_switch/2)+zeta1),sqrt(1-zeta1^2))-atan2(zeta1, sqrt(1-zeta1^2)));
% elseif zeta1 > 1
%     Denom_int = (1/sqrt(zeta1^2-1))*log((tan(lambda_switch)*(zeta1+sqrt(zeta1^2-1))+1) / (tan(lambda_switch)*(zeta1-sqrt(zeta1^2-1))+1));
% else % zeta1 = 1
%     Denom_int = (1-tan(pi/4 + lambda_switch/2));
% end
% % Denom_int is a integral part of the Denominator
% Denom_1 = R0/(1+zeta1*sin(lambda_switch)); % in order to simplify the accel_max calculation code.
% Numerator = N * missile.V^2 * (sin(sigma_d)+eta*sin(lambda_switch));
% 
% gamma = (Numerator) / (Denom_1 * exp(Denom_int / tan(sigma_d))); % Gamma is calculated accel_max via designated amf(& lambda_switch)
% 
% if gamma > accel_max    % Finding amf_max via code is not validated due to limitation of my own.
%     disp 'Required accel_max is over than designated accel_max. Reduce amf.'
%     return
% end
% ==================================================================================================================
% =========Computing amf_max from EQ.32 using Newton-Raphson Method(or other method) is not implemented!!===========

% Sim params
S  = 100;                        % max sim duration, seconds
dt = 0.001;                       % time-step size, seconds 
Niter = S/dt;                    % max iteration num

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

true_sigma = sigma_0;

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

[z, v, count] = MeasM_TKF_StatT(sigma_0, true_sigma, sat_upper,sat_lower);       % MEASUREMENT MODEL

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
    if k==1
        missile.HDG_dot = 0;
    else
        missile.HDG_dot = (missile.HDG - missile.HDG_prev)/dt;
    end

    K = 300;
    % =======================================================================================
    % IACCG, accel_m: Missile Acceleration
    accel_m = Vm * lambda_dot + K * (sigma_d - E_sigma);
    
    if (abs(lambda) >= abs(lambda_switch))||(index)
        accel_m = N * Vm * lambda_dot;
        index = 1;
    end
   
    % Acceleration Saturation Condition
    if ( abs(accel_m) > 100 && accel_m < 0 )% m/s^2
        accel_m = -100;
    elseif ( abs(accel_m) > 100 && accel_m > 0 )
            accel_m = 100;

    else
            accel_m;
    end
    % =======================================================================================

    % Terminate sim at interception
    if abs(E_range) <= 30 % interception 기준
        disp('Estimated Range is less than 30m')
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
    logger.m.HDG_dot(k) = missile.HDG_dot;

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
    %logger.xpp(k) = xpp;
    logger.TKF_lookerror(k) = true_sigma - E_sigma;
    %logger.EKF_lookerror(k) = true_sigma - x_EKF;

    logger.M_Accel(k) = accel_m;
    logger.M_Accel_x(k) = M_accel_x;
    logger.M_Accel_y(k) = M_accel_y;
    % logger.normKalman(k) = norm(K);

%Cov_P
% 
end



%--------------------------------------------------------------------------
% Visualize results
%--------------------------------------------------------------------------
close all;

% Impact index
[M,I] = min(logger.E_range);

% Range
%-------------------------------------
figure();
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

% % Error between Estimation and True Range
% figure;
% scatter(logger.t, logger.est_error,'filled'); hold on;
% title(['IACCG(DPP&PPN), N = ' num2str(N) ])
% legend('Range Error')
% ylabel('Range Error btw Estimated and True Range [m])')
% xlabel('Elapsed time (sec)')

% % Estimated and True Lambda (LOS angle)
% figure;
% plot(logger.t, rad2deg(logger.lambda)); hold on;
% plot(logger.t, rad2deg(logger.true_lambda),'linewidth',2); hold on;
% title('Estimated and True LOS Angle')
% legend('Estimated LOS angle', 'True LOS angle')
% xlabel('Elapsed time(sec)')
% ylabel('LOS Angle(deg)')

% Lambda dot True and Estimated
%--------------------------------------
figure;
plot(logger.t, rad2deg(logger.lambda_dot)); hold on;
plot(logger.t, rad2deg(logger.true_lambda_dot),'linewidth',2); hold on;
title(['IACCG(DPP&PPN), N = ' num2str(N) ])
legend('Estimated Lambda dot', 'True Lambda dot')
ylabel('Lambda dot(deg/s)')
xlabel('Elapsed time (sec)')

% Estimated, True and Measured Sigma(Look angle)
%-------------------------------------
figure;
plot(logger.t, rad2deg(logger.E_sigma)); hold on;
plot(logger.t, rad2deg(logger.true_sigma),'linewidth',2); hold on;
scatter(logger.t, rad2deg(logger.z_sigma)); hold on;
plot(logger.t, rad2deg(logger.xpp)); hold on;
title('Estimated and True Look Angle')
legend('Estimated Look Angle', 'True Look Angle', 'Measured Look Angle', 'Predicted Look Angle')
xlabel('Elapsed time(sec)')
ylabel('Look Angle (deg)')

% % Missile Heading angle derivative
% % -------------------------------------
% figure;
% plot(logger.t, rad2deg(logger.M_HDG_dot),'linewidth',1.5); hold on;
% title('Missile Heading angle derivative')
% legend('Missile HDG Derivative')
% xlabel('Elapsed time(sec)')
% ylabel('Missile HDG angle derivative (deg/s)')

% % Heading
% %-------------------------------------
% figure;
% plot(logger.t, rad2deg(logger.m.HDG)); hold on;
% plot(logger.t, rad2deg(logger.lambda)); hold on;
% title('Modified DPP + PPN')
% legend('Missile heading', 'LOS angle', 'Location', 'nw' )
% ylabel('Heading angle (deg)')
% xlabel('Elapsed time (sec)')

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

% % Acceleration
% %--------------------------------------
% figure;
% plot(logger.t, logger.M_Accel); hold on;
% title('IACCG(DPP+PPN)')
% legend('Missile Acceleration')
% ylabel('Missile Acceleration (m/s^2)')
% xlabel('Elapsed time (sec)')

% % Kalman Gain
% figure;
% plot(logger.t, logger.Kalman1); hold on;
% plot(logger.t, logger.Kalman2); hold on;
% plot(logger.t, logger.normKalman); hold on;
% title(['Pure Proportional Navigation, N = ' num2str(N) ])
% legend('norm(Kalman)')
% xlabel('Elapsed time(sec)')
% 
% % Expectation Ez
% figure;
% plot(logger.t, rad2deg(logger.Ez_total), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.Ez_l), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.Ez_h), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.Ez_us), 'linewidth', 1.2); hold on;
% title('IACCG TKF Expectation(Ez) Components');
% legend('Ez_{total}', 'Ez_l', 'Ez_h', 'Ez_{us}')
% xlabel('Elapsed time (sec)')
% ylabel('deg')

% Probabilities pl, ph, pus
figure;
plot(logger.t, logger.pl, 'linewidth', 1.2); hold on;
plot(logger.t, logger.ph, 'linewidth', 1.2); hold on;
plot(logger.t, logger.pus, 'linewidth', 1.2); hold on;
title('Probabilities')
legend('p_l', 'p_h', 'p_{us}')
ylabel('Probability')
xlabel('Elapsed time(sec)')

% Bias term
figure;
plot(logger.t, rad2deg(logger.bias), 'linewidth', 1.2);
title('Bias Term')
legend('Bias')
ylabel('Bias Value')
xlabel('Elapsed time(sec)')

% % Ez with Bias VS. No bias
% figure;
% plot(logger.t, rad2deg(logger.Ez_total), 'linewidth', 1.2); hold on;
% plot(logger.t, rad2deg(logger.Ez_total_nobias), 'linewidth', 1.2); hold on;
% title('Ez_{total} with Bias VS. No bias')
% legend('Ez_{total} with Bias', 'Ez_{total} without Bias')
% ylabel('deg')
% xlabel('Elapsed time(sec)')
% 
% % H VS. HH VS. HH w/o bias
% figure; 
% plot(logger.t, logger.H, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.HH, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.HH_nobias, 'linewidth', 1.2); hold on;
% title('H_{EKF} VS. H_{TKF} VS. H_{TKF} w/o bias')
% legend('H_{EKF}', 'H_{TKF}', 'H_{TKF} w/o bias')
% ylabel('Norm of the matrices')
% xlabel('Elapsed time(sec)')

% % R VS. RR VS. RR w/o bias
% figure;
% plot(logger.t, logger.R, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.RR, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.RR_nobias, 'linewidth', 1.2); hold on;
% title('R_{EKF} VS. R_{TKF} VS. R_{TKF} w/o bias')
% legend('R_{EKF}', 'R_{TKF}', 'R_{TKF} w/o bias')
% ylabel('Norm of the matrices')
% xlabel('Elapsed time(sec)')

% % TKF Kalman gain VS. EKF Kalman gain
% figure;
% plot(logger.t, logger.KK, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.K_EKF, 'linewidth', 1.2); hold on;
% title('TKF K_k VS. EKF K_k')
% legend('TKF K_k', 'EKF K_k')
% ylabel('Norm of the Kalman gain')
% xlabel('Elaped time(sec)')

% % TKF Kalman Gain VS. TKF_nobias Kalman Gain
% figure;
% plot(logger.t, logger.KK, 'linewidth', 1.2); hold on;
% plot(logger.t, logger.K_nobias, 'linewidth', 1.2); hold on;
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

% % Look Angle Estimation Error
% figure;
% plot(logger.t, rad2deg(logger.TKF_lookerror)); hold on;
% plot(logger.t, rad2deg(logger.EKF_lookerror)); hold on;
% title('Look Angle Error')
% legend('TKF Look Angle Error', 'EKF Look Angle Error')
% ylabel('deg')
% xlabel('Elapsed time(sec)')
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
% 
% figure;
% plot(logger.t, rad2deg(logger.E_sigma), 'linewidth', 2); hold on;
% plot(logger.t, rad2deg(logger.E_sigma_nobias),'linewidth', 2); hold on;
% plot(logger.t, rad2deg(logger.true_sigma), 'linewidth', 2); hold on;
% title('Look Angle Trajectory')
% legend('Estimated Look Angle', 'Estimated Look Angle w/o bias', 'True Look Angle')
% ylabel('deg')
% xlabel('Elapsed time (sec)')