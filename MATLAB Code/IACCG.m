% Composite Guidance Scheme for IAC Against a Nonmnaneuvering Moving Target
% B.G.Park, H.H.Kwon, and Y.H.Kim
% Journal of Guidance, Control, and Dynamics, Vol.39, No.5, May 2016
clear; close all

% Pro-Nav gain, unitless (typically set between 3-5) / N >= 3
N = 3;

% Missile init params
missile.y   =  0;%7000*sind(27); % north pos, m 
missile.x   =  0;%-7000*cosd(27); % east pos, m  
%missile.HDG =  deg2rad(+30); % heading, rad  
am0 = deg2rad(30);
missile.HDG = am0;%deg2rad(-22);
sigma_0 = missile.HDG; % Initial Look Angle
%sigma_0 = deg2rad(5);
missile.V   =  250; % velocity magnitude, m/s 
missile.yV  =  missile.V * sin(missile.HDG); % y velocity, m/s 
missile.xV  =  missile.V * cos(missile.HDG); % x velocity, m/s 

sigma_max = deg2rad(+45); % rad, Seeker's Look Angle Limit --> Symmetry is not reflected, multiply 2 if wanted.

% Target init params
target.anchored = false; % fix target to init pos when true
target.y   =  0; % y pos, m
target.x   =  5000;%0; %5000;   % x pos, m
target.HDG = deg2rad(0); % heading, rad
target.yV  =  0; % y velocity, m/s
target.xV  =  0; %50; % x velocity, m/s
target.V   = sqrt(target.yV^2 + target.xV^2); % velocity magnitude, m/s

% Other Initial Conditions and Parameters
eta = target.V / missile.V; 
% ==========> Must be positioned in for loop for reality: NO. Even realistic one did use average speed Mach 0.7

% sigma_d Selection in order to bring out lambda_switch and also compare whether initial setup is good or not.
sigma_d_lower =  asin(eta); % rad, Lower bound for sigma_d
sigma_d_upper =  min(acos(2*eta), sigma_max); % rad, Upper bound for sigma_d

sigma_d = deg2rad(45); % rad, Desired Look Angle

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

accel_max = 100; % m/s, acceleration limit ==========> OFF for reality and revise as real-time one in the for loop.
% Guess it doesn't have to be in the for loop bcuz the realistic one also determined the accel_max as a constant.

amf = deg2rad(-120); % rad, Desired Impact Angle
% labmda_switch and amf are functions of each other, thus if you decide amf, then lambda_switch can be calculated itself.
lambda_switch = (( (atan2(sin(amf)-eta*sin(target.HDG), cos(amf)-eta*cos(target.HDG)))) - ((amf) - (sigma_d)) / N ) * (N/(N-1));
% rad, lambda_switch is specific LOS angle where the guidance gain changes, thus logic switch occurs. 

% K Gain Selection
K = 300; % Select K value as you wish, IF YOU CAN.
R0 = norm([(target.x - missile.x), (target.y - missile.y)]); % Initialized R0 in order to calculate K_upper  

K_upper = (accel_max * R0 + missile.V^2 * sin(sigma_0)) / (R0*(abs(sigma_d) - sigma_0));

if K > K_upper
    disp 'K must be less than K_upper'
    return
end

% accel_max comparison thus determining whether abs(amf) =< abs(amf_max)
if zeta1 < 1
    Denom_int = (2/sqrt(1-zeta1^2))*(atan2((tan(lambda_switch/2)+zeta1),sqrt(1-zeta1^2))-atan2(zeta1, sqrt(1-zeta1^2)));
elseif zeta1 > 1
    Denom_int = (1/sqrt(zeta1^2-1))*log((tan(lambda_switch)*(zeta1+sqrt(zeta1^2-1))+1) / (tan(lambda_switch)*(zeta1-sqrt(zeta1^2-1))+1));
else % zeta1 = 1
    Denom_int = (1-tan(pi/4 + lambda_switch/2));
end
% Denom_int is a integral part of the Denominator
Denom_1 = R0/(1+zeta1*sin(lambda_switch)); % in order to simplify the accel_max calculation code.
Numerator = N * missile.V^2 * (sin(sigma_d)+eta*sin(lambda_switch));

gamma = (Numerator) / (Denom_1 * exp(Denom_int / tan(sigma_d))); % Gamma is calculated accel_max via designated amf(& lambda_switch)

if gamma > accel_max    % Finding amf_max via code is not validated due to limitation of my own.
    disp 'Required accel_max is over than designated accel_max. Reduce amf.'
    return
end
% %=========Computing amf_max from EQ.32 using Newton-Raphson Method(or other method) is not implemented!!============


% Current prams
% Applied as a disturbance to missile's and target's velocity
%current.yV =  0.0; % y velocity, m/s
%current.xV =  0.0; % x velocity, m/s

% Sim params
S  = 100;    % max sim duration, seconds
dt = 0.0001;     % time-step size, seconds 
Niter = S/dt; % max iteration num

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
logger.lambda_dot = nan(1, Niter);
logger.sigma = nan(1, Niter);  % Look angle 

logger.lambda_dot_test = nan(1,Niter);

logger.M_Accel = nan(1, Niter);
logger.M_Accel_x = nan(1, Niter);
logger.M_Accel_y = nan(1, Niter);

%--------------------------------------------------------------------------
% Init sim
%--------------------------------------------------------------------------
RTP_last = [(target.x - missile.x);... % delta y % RTP: Relative Target Position w.r.t Missile; 2x1 Vector
    (target.y - missile.y)];	       % delta x

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
    VTP = (RTP - RTP_last) ./ dt;  
    
    % Closing velocity, m/s
    Vc = abs(-RTP'*VTP / LOSdistance); 
    
    % Missile velocity, m/s
    Vm = sqrt(missile.yV^2 + missile.xV^2);
    
    % Target velocity, m/s
    %Vt = sqrt(target.yV^2 + target.xV^2);
    
    % Line-of-sight (LOS) angle, rad w.r.t Ref.
    if k==1
        lambda = atan2(RTP(2), RTP(1));
    else
    lambda_prev = lambda;
    lambda = atan2(RTP(2), RTP(1)); % lambda: Ghose책 기준 Theta
    end

    
    % LOS angle time derivative (d/dt lambda), rad/s
    if k==1
    lambda_dot = (RTP(1)*VTP(2) - RTP(2)*VTP(1)) / LOSdistance^2;
    else
    lambda_dot_test = (lambda - lambda_prev)/dt;
    lambda_dot = (RTP(1)*VTP(2) - RTP(2)*VTP(1)) / LOSdistance^2;
    end

    % Look Angle, rad
    sigma = missile.HDG - lambda; % Look angle, rad
    
    % Target heading, rad
    %beta = atan2(target.yV, target.xV);
    
    % Lead angle, rad % Lead angle
    %L = asin(Vt * sin(beta + lambda) / Vm);
    
    % True Proportional Navigation, rad/s2
    % nc = N * Vc * lambda_dot;
    
    %-------------------------------------
    % Pure Proportional Navigation, accel_m: Missile Acceleration
    accel_m = Vm * lambda_dot + K * (sigma_d - sigma); % Initial Phase accel_m, Modified DPP
    
    if abs(lambda) >= abs(lambda_switch)
        accel_m = N * Vm * lambda_dot; % Final Phase accel_m, PPN
    end
   
    % Acceleration Saturation Condition
    if ( abs(accel_m) > 100 && accel_m < 0 )% m/s^2
        accel_m = -100;

    elseif ( abs(accel_m) > 100 && accel_m > 0 )
            accel_m = 100;

    else
            accel_m;
    end
    %-------------------------------------
    
    % Terminate sim at intercept
    if abs(LOSdistance) <= 0.3 % interception 기준
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
    missile.HDG = atan2(missile.yV, missile.xV);    % renewed
    
    % Update target pos for time step
    if(~target.anchored)
        target.y = target.y + target.yV*dt; % + current.yV*dt;
        target.x = target.x + target.xV*dt; % + current.xV*dt;
    end
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
    logger.lambda_dot(k) = lambda_dot;
    logger.sigma(k) = rad2deg(sigma);
    logger.M_Accel(k) = accel_m;
    logger.M_Accel_x(k) = M_accel_x;
    logger.M_Accel_y(k) = M_accel_y;
    
    if k==1
    
    else
    logger.lambda_dot_test(k) = lambda_dot_test;
    end

end

%--------------------------------------------------------------------------
% Visualize results
%--------------------------------------------------------------------------
close all;

% Impact index
[M,I] = min(logger.R);

% Range
%-------------------------------------
figure;
plot(logger.t,logger.R); hold on;
scatter(k*dt, M, 'filled')
title('Modified DPP + PPN @ lambda switch')
legend('Range',...
    ['Intercept: r = ' num2str(LOSdistance) ' m, t = ' num2str(k*dt) ' s'],...
    'Location','nw')
ylabel('Range (m)')
xlabel('Elapsed time (sec)')
set(gca, 'YDir', 'reverse')

% Heading
%-------------------------------------
rad2deg = @(x) x*180/pi;  % helper function

figure;
plot(logger.t, rad2deg(logger.m.HDG)); hold on;
plot(logger.t, rad2deg(logger.Lambda)); hold on;
title('Modified DPP + PPN @ lambda switch')
legend('Missile heading', 'LOS angle', 'Location', 'nw' )
ylabel('Heading angle (deg)')
xlabel('Elapsed time (sec)')

% Position
%-------------------------------------
figure;
scatter(logger.mx, logger.my, 'filled'); hold on;
scatter(logger.tx, logger.ty, 'filled');
set(gca, 'DataAspectRatio',[1 1 1]);
title('Modified DPP + PPN @ lambda switch')
legend('Pursuer', 'Target', 'Location', 'southoutside',...
    'Orientation','horizontal')
xlabel('+x (m)')
ylabel('+y (m)')

% Acceleration
%--------------------------------------
figure;
plot(logger.t, logger.M_Accel); hold on;
title('Modified DPP + PPN @ lambda switch')
legend('Missile Acceleration')
ylabel('Missile Acceleration (m/s^2)')
xlabel('Elapsed time (sec)')

% Look Angle
%-------------------------------------
figure; grid on;
plot(logger.t, logger.sigma); hold on;
title('Modified DPP + PPN @ lambda switch')
legend('Look Angle')
ylabel('Look Angle sigma (deg)')
xlabel('Elapsed time (sec)')

% Lambda dot (Strictly increasing or decreasing)
figure; grid on;
plot(logger.t, logger.lambda_dot); hold on;
%plot(logger.t, logger.lambda_dot_test); 
title('Modified DPP + PPN @ lambda switch')
legend('Lambda dot')
ylabel('Lambda dot (rad/s)')
xlabel('Elapsed time (sec)')
% Position Animtaion
%-------------------------------------

%figure(4);
%ax = gca;
%axis(ax,'manual');
%axis(ax, [min([logger.mx; logger.tx]) max([logger.mx; logger.tx]) min([logger.my; logger.ty]) max([logger.my; logger.ty])]);

%line1=animatedline(ax, 'Color', 'b', 'Linewidth', 2);
%line2=animatedline(ax, 'Color', 'r', 'Linewidth', 2);

%xlabel('X location');
%ylabel('Y location');
%title('M&T Trajectory Animation');
%legend('Dataset 1', 'Dataset 2');
%grid on;

%for ii = 1:length(logger.mx)
%    addpoints(line1, logger.mx(ii), logger.tx(ii));
%    addpoints(line2, logger.my(ii), logger.ty(ii));
%    drawnow;
%    pause(2);
%end

