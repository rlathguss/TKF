% "Look Angle Constrained Impact Angle Control Guidance Law for Homing Missiles with Bearings-Only Measurements"
% IACGBOM 

clear; close all

% Pro-Nav gain, (typically set between 3-5)
% N = 3;

% Missile initial parameters
missile.y   =  0; % North pos, m
missile.x   =  0; % East pos, m 
am0 = deg2rad(15); % Alpha_M_0: Initial M HDG: Fixed constant.
missile.HDG =  am0; % Initial M HDG: updates in the loop
missile.V   =  250; % velocity magnitude, m/s 
missile.yV  =  missile.V * sin(missile.HDG); % y velocity, m/s 
missile.xV  =  missile.V * cos(missile.HDG); % x velocity, m/s 
% sigma_0 = deg2rad(5);

% Target init params
%target.anchored = false; % fix target to init pos when true
target.y   =  0; % y pos, m
target.x   =  10000;   % x pos, m
target.yV  =  0; % y velocity, m/s
target.xV  =  0; % x velocity, m/s
target.V   = sqrt(target.yV^2 + target.xV^2); % velocity magnitude, m/s

% Desired engagement parameters
amf = deg2rad(-32);       % Alpha_M_Final: Desired impact angle [-pi, 0], rad, cw dir
sig_max = deg2rad(+45); % Maximum look angle, rad, symmetric to M.body axis
% sig_max must be more or equal than am0

% ====================== INITIAL PARAMETERS =============================
Rf = 0.01; % m      % Rf=0.5
%phi1 = deg2rad(0.15);
k2 = deg2rad(10);
phi1 = (0.15);
%k2 = (10);

k1 = deg2rad(9);
%k1 = deg2rad(rad2deg(sig_max) - 0.01);        % Originally, sigma_max(FOV limit) - 0.01 [deg]
%k1 = sig_max - deg2rad(0.01);
% Sim params
Simtime  = 100;                        % max sim duration, seconds
dt = 0.01;                       % time-step size, seconds 
Niter = Simtime/dt;                    % max iteration num

% Pre-allocate logger
logger.t   = nan(1, Niter);    % elapsed time
logger.my  = nan(1, Niter);    % pursuer north pos
logger.mx  = nan(1, Niter);    % pursuer east pos
logger.myV = nan(1, Niter);    % pursuer north vel
logger.mxV = nan(1, Niter);    % pursuer east vel
logger.m.HDG = nan(1, Niter);    % pursuer heading angle

logger.ty  = nan(1, Niter);    % target north y pos
logger.tx  = nan(1, Niter);    % target east x pos
logger.tyV = nan(1, Niter);    % target north y vel
logger.txV = nan(1, Niter);    % target east x vel

logger.R   = nan(1, Niter);    % pursuer/target range
logger.Lambda = nan(1, Niter); % pursuer/target bearing
logger.Lambdadot = nan(1, Niter); % Lambda dot
logger.sigma = nan(1, Niter);  % Look angle 
logger.S = nan(1,Niter);
logger.e1 = nan(1,Niter);
logger.e2 = nan(1,Niter);

logger.M_Accel = nan(1, Niter);
logger.M_Accel_x = nan(1, Niter);
logger.M_Accel_y = nan(1, Niter);

%--------------------------------------------------------------------------
% Init sim
%--------------------------------------------------------------------------
RTP_last = [(target.x - missile.x);... % delta x % RTP: Relative Target Position w.r.t Missile; 2x1 Vector
    (target.y - missile.y)];	       % delta y

% Target pos
target.y = target.y + target.yV*dt; % dt시간 뒤의 target.y pos
target.x = target.x + target.xV*dt; % dt시간 뒤의 target.x pos 

% Missile pos
missile.y = missile.y + missile.yV*dt; % dt시간 뒤의 missile.y pos 
missile.x = missile.x + missile.xV*dt; % dt시간 뒤의 missile.x pos

%--------------------------------------------------------------------------
% Run sim
%--------------------------------------------------------------------------
for k = 1:Niter 
    
    % Relative position in the inertial frame, m ; Relative Target Position w.r.t Missile
    RTP = [(target.x - missile.x);... % delta x
        (target.y - missile.y)];      % delta y
    
    % Range(LOSdistance) to target 
    LOSdistance = norm(RTP);
    
    % Relative velocity in the inertial frame, m/s
    % VTP = [(target.yV - missile.yV);... % delta y velocity
    %        (target.xV - missile.xV)];   % delta x velocity
    VTP = (RTP - RTP_last) ./ dt;  % ./ 는 A의 각 요소를 B의 각 해당하는 요소로 나누는 cmd
    
    % Closing velocity, m/s
    Vc = abs(-RTP'*VTP / LOSdistance); 
    
    % Missile velocity, m/s
    Vm = sqrt(missile.yV^2 + missile.xV^2);
    
    % Target velocity, m/s
    %Vt = sqrt(target.yV^2 + target.xV^2);
    
    % Line-of-sight (LOS) angle, rad w.r.t Ref.
    lambda = atan2(RTP(2), RTP(1)); % lambda: Ghose책 기준 Theta
    
    % LOS angle time derivative (d/dt lambda), rad/s
    lambda_dot = (RTP(1)*VTP(2) - RTP(2)*VTP(1)) / LOSdistance^2;
    
    % Look Angle, rad
    sigma = missile.HDG - lambda; % Look angle, rad
    
    % Lead angle, rad % Lead angle
    %L = asin(Vt * sin(alpha_t + lambda) / Vm);
    
    % True Proportional Navigation, rad/s2
    % nc = N * Vc * lambda_dot;
    
    % =======================================================================================
    e1 = lambda - amf;      %rad
    e2 = sigma;             %rad
    S = e2 - k1*(e1/sqrt(e1^2+phi1^2));
    f2 = ((1+k1*phi1^2/((e1^2+phi1^2)^1.5))*abs(sin(sigma)));

    % IACCG, accel_m: Missile Acceleration
    accel_m = -(Vm*f2/Rf + k2)*Vm*tanh(0.01*S); % a=10
    
    % Acceleration Saturation Condition
    if ( abs(accel_m) > 100 && accel_m < 0 )% m/s^2
        accel_m = -100;
    elseif ( abs(accel_m) > 100 && accel_m > 0 )
            accel_m = 100;

    else
            accel_m;
    end
    % =======================================================================================
    
    % Terminate sim at intercept
    if (abs(LOSdistance)<=3) && (abs(rad2deg(amf - logger.Lambda(k-1)))<2.5) % interception 기준
        disparg = ['Range is less than 3m. Final range is ', num2str(norm(RTP)) '. LOS angle error is ', num2str(rad2deg(amf - logger.Lambda(k-1)))];
        disp(disparg)
        break;
    end

    if (norm(RTP)-norm(RTP_last))/dt > 0
        disparg = ['Range is increasing(FAILED). ' ...
        'Final range & LOS angle error are ', num2str(LOSdistance),' & ' num2str(rad2deg(amf-logger.Lambda(k-1))), ' Desired impact angle is ', num2str(rad2deg(amf))];
        disp(disparg)
        break;
    end
    
    % Update missile pos for time-step(dt) and apply current slip
    missile.y = missile.y + missile.yV*dt; %+ current.yV*dt;  % renewed
    missile.x = missile.x + missile.xV*dt; %+ current.xV*dt;  % renewed
    
    % Compute the y/x acceleration commands
    % In pure pro-nav accel commands are applied normal to missile's velocity vector
    M_accel_y = +accel_m * cos(missile.HDG);
    M_accel_x = -accel_m * sin(missile.HDG);
    
    % Update missile y/x velocities
    missile.yV = missile.yV + M_accel_y*dt;
    missile.xV = missile.xV + M_accel_x*dt;
    
    % Update missile heading
    missile.HDG = atan2(missile.yV, missile.xV);    % renewed
    
    RTP_last = RTP;     % renewed
    
    %-------------------------------------
    % Log time-step data
    logger.t(k)   = k*dt;
    logger.my(k)  = missile.y;
    logger.mx(k)  = missile.x;
    logger.myV(k) = missile.yV;
    logger.mxV(k) = missile.xV;
    logger.m.HDG(k) = missile.HDG;
    logger.ty(k)  = target.y;
    logger.tx(k)  = target.x;
    logger.tyV(k) = target.yV;
    logger.txV(k) = target.xV;
    logger.R(k)   = LOSdistance;
    logger.Lambda(k) = lambda;
    logger.Lambdadot(k) = lambda_dot;
    logger.sigma(k) = sigma;
    logger.M_Accel(k) = accel_m;
    logger.M_Accel_x(k) = M_accel_x;
    logger.M_Accel_y(k) = M_accel_y;
    logger.S(k) = S;
    logger.e1(k) = e1;
    logger.e2(k) = e2;
    
end

%--------------------------------------------------------------------------
% Visualize results
%--------------------------------------------------------------------------
close all;

% Impact index
[M,I] = min(logger.R);

% Range
%-------------------------------------
% figure; grid on;
% plot(logger.t,logger.R); hold on;
% scatter(k*dt, M, 'filled')
% title('Relative Range')
% legend('Range',...
%     ['Intercept: r = ' num2str(LOSdistance) ' m, t = ' num2str(k*dt) ' s'],...
%     'Location','nw')
% ylabel('Range (m)')
% xlabel('Elapsed time (sec)')
% set(gca, 'YDir', 'reverse')

% Missile Heading , LOS angle, Look angle
%-------------------------------------
figure; grid on;
plot(logger.t, rad2deg(logger.m.HDG),'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
plot(logger.t, rad2deg(logger.Lambda),'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
plot(logger.t, rad2deg(logger.sigma),'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
title('Missile Heading Angle & LOS Angle & Look Angle')
legend('Missile heading', 'LOS angle', 'Look Angle' )
ylabel('(deg)')
xlabel('Elapsed time (sec)')

% Position
%-------------------------------------
figure; grid on;
%scatter(logger.mx, logger.my, 'filled'); hold on;
plot(logger.mx, logger.my, 'linewidth',1.8, 'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
scatter(logger.tx, logger.ty, 'filled');
%set(gca, 'DataAspectRatio',[1 1 1]);
title('M&T Trajectory')
legend('Pursuer', 'Target', 'Location', 'southoutside',...
    'Orientation','horizontal')
xlabel('+x (m)')
ylabel('+y (m)')

% Acceleration
% --------------------------------------
% figure; grid on;
% plot(logger.t, logger.M_Accel,'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
% title('Missile Acceleration')
% legend('Missile Acceleration')
% ylabel('Missile Acceleration (m/s^2)')
% xlabel('Elapsed time (sec)')

% Look Angle
%-------------------------------------
% figure; grid on;
% plot(logger.t, rad2deg(logger.sigma),'linewidth', 1.8); hold on;
% title('Look Angle')
% legend('Look Angle')
% ylabel('Look Angle sigma (deg)')
% xlabel('Elapsed time (sec)')

% Lambda dot
%-------------------------------------
% figure; grid on;
% plot(logger.t, logger.Lambdadot); hold on;
% title('LOS angle rate')
% legend('Lambda dot')
% ylabel('Lambda dot')
% xlabel('Elapsed time(sec)')

%Sliding Surface
%-------------------------------------
% figure;
% plot(logger.t, logger.S,'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
% title('Sliding Surface')
% legend('Sliding Surface Variable')
% ylabel('unit')
% xlabel('Elapsed time(sec)')

%Error Variables
%-------------------------------------
figure;
plot(logger.t, rad2deg(logger.e1),'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
plot(logger.t, rad2deg(logger.e2),'linewidth', 1.8,'Marker','o','MarkerIndices',1:3000:length(logger.t)); hold on;
title('Error Variables')
legend('e1', 'e2')
ylabel('deg')
xlabel('Elapsed time (sec)')
