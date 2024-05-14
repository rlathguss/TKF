% "Capture Region of Tactical Missile Equipped with Semi-Active Laser Seeker Using Tobit Kalman Filter"
% TKF VALIDATED
% Measure ONLY Look Angle

function [z, v, count] = MeasM_TKF_StatT(sigma_0, true_sigma, sat_upper, sat_lower)

persistent sigma %sat_upper sat_lower
persistent firstRun

if isempty(firstRun)               
    sigma = sigma_0;
    % sat_upper = deg2rad(+7.5); % 7.5
    % sat_lower = -sat_upper;
end

v = 0.0052*randn;       % rad, Covariance(sigma_v) = 0.0052 , 3*Covariance(sigma_v) = 99% Region

if firstRun==1
    z_sigma = true_sigma + v;           % rad, Look angle measurement = True + Noise: Idealized Seeker
else % firstRun is empty(k=1: Initial run)
    z_sigma = sigma + v;
    firstRun = 1;
end

count = 1;      % count follows the Bernoulli distribution: UNSAT.Meas.: count=1 or SAT.Meas.: count=0

if ( abs(z_sigma) >= sat_upper && z_sigma > 0)       % Upper saturated
    % v = 0;
    count = 0;
    z_sigma = sat_upper;
    
elseif ( abs(z_sigma) >= sat_upper && z_sigma < 0)   % Lower saturated
    % v = 0;
    count = 0;
    z_sigma = sat_lower;
    
end

z = z_sigma;