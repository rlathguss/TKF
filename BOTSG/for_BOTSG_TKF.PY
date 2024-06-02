import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from SysM import sysm_tkf
from MeasM import measm_tkf


def run_simulation(sat_upper1, amf1, sigma_desired_deg1):

    sat_upper = np.deg2rad(sat_upper1)             
    sat_lower = -sat_upper

    am0 = np.deg2rad(+30);                      # Alpha_M_0: Initial M HDG: Fixed constant.
    sigma_0 = am0                               # rad, Seeker's Initial Look Angle 
    sigma_desired = sat_upper + np.deg2rad(sigma_desired_deg1)  # rad, Seeker's Look Angle Limit


    amf = np.deg2rad(amf1)       # rad, Desired Impact Angle
    accel_m = 0                 # Initial Missile Acceleration, m/s^2
    accel_max = 1000            # m/s^2, acceleration limit

    # Missile init params
    missile = {'y': 0,
            'x': 0,
            'HDG': am0,
            'V': 300}                                   # velocity magnitude, m/s
    missile['yV'] = missile['V'] * np.sin(missile['HDG'])  # y velocity, m/s
    missile['xV'] = missile['V'] * np.cos(missile['HDG'])  # x velocity, m/s


    # Target init params
    target = {'y': 0,
            'x': 5000,
            'HDG': 0,              # heading, rad
            'yV': 0,               # y velocity, m/s
            'xV': 0}               # x velocity, m/s
    target['V'] = np.sqrt(target['yV']**2 + target['xV']**2)

    # Verifying Guidance Gains N=(k1- k2*lambda): N is defined in the PPN loop
    beta = 4*(sigma_desired - am0)/amf;                                  # beta is less or equal than 0
    k1 = (beta - np.sqrt(beta**2 + 4 * beta * am0 / amf))/2 + 1;         # Guidance gain k1
    k2 = (amf*(beta-np.sqrt(beta**2+4*beta*am0/amf))+2*am0)/amf**2;      # Guidance gain k2


    # ================ Desired Impact Angle Achievement Survey ================
    if sigma_desired < am0:
        print('aaaa')
        return 0
    if amf < -2 * am0:
        if sigma_desired < -amf ** 2 / (4 * (amf + am0)):
            print('bbbb')
            return 0
    else:
        if sigma_desired < am0:
            print('cccc')
            return 0

    # =========================================================================

    S = 50                  # 최대 시뮬레이션 지속 시간, 초
    dt = 0.01               # 시간 단계 크기, 초
    Niter = int(S / dt)     # 최대 반복 횟수
    logger = {
        'true_range': [],    
    }

    #--------------------------------------------------------------------------
    # Init sim
    #--------------------------------------------------------------------------
    true_sigma = sigma_0

    # Target pos update
    target['y'] = target['y'] + target['yV'] * dt
    target['x'] = target['x'] + target['xV'] * dt

    # Missile pos update
    missile['y'] = missile['y'] + missile['yV'] * dt
    missile['x'] = missile['x'] + missile['xV'] * dt

    #--------------------------------------------------------------------------
    # Run sim
    #--------------------------------------------------------------------------
    for k in range(Niter):
        RTP_true = np.array([target['x'] - missile['x'], target['y'] - missile['y']])           # True RTP (Range = LOS distance)
        true_range = np.linalg.norm(RTP_true)
        Vm = np.sqrt(missile['yV']**2 + missile['xV']**2)                                       # Missile velocity, m/s
        true_lambda = np.arctan2(RTP_true[1], RTP_true[0])                                      # Line-of-sight (LOS) angle, rad w.r.t Ref.
        true_lambda_dot = - missile['V'] * np.sin(missile['HDG'] - true_lambda) / true_range    # LOS angle time derivatives (d/dt lambda), rad/s
        true_sigma = missile['HDG'] - true_lambda                                               # True Look Angle, rad


        # ================================================== TOBIT KALMAN FILTER ====================================================
        z, v = measm_tkf(sigma_0, true_sigma, sat_upper, sat_lower,k)  # MEASUREMENT MODEL

        # SYSTEM MODEL & GUIDANCE FILTER
        E_range, E_sigma, Cov_P, Ez_total, pdf_l, pdf_h, pl, ph, pus, bias_lambda, Ez_us, RR, HH, K, xpp = sysm_tkf(z, dt, missile, accel_m, sat_upper, sat_lower,k)

        E_lambda = missile['HDG'] - E_sigma                         # Line-of-sight (LOS) angle, rad w.r.t Ref.
        E_lambda_dot = - missile['V'] * np.sin(E_sigma) / E_range   # LOS angle time derivatives (d/dt lambda), rad/s


        # ================================================== Missile Acceleration ====================================================   
        N = k1 - k2 * E_lambda
        accel_m = N * Vm * E_lambda_dot

        # Acceleration Saturation Condition
        if abs(accel_m) > accel_max: 
            if accel_m < 0:
                accel_m = -accel_max
            elif accel_m > 0:
                accel_m = accel_max
            else:
                pass

        # ================================================== Terminate sim condition ====================================================   
        interception_range = 2.5
        # Terminate sim at interception
        if RTP_true[1]> -2 and abs(RTP_true[0]) <= interception_range:
            impact_angle_error = abs(np.rad2deg(amf - missile['HDG']))
            if impact_angle_error <= 2.5:
                print('zzzzzzzzzzzzz')
                return 1
            else:
                print('dddd')
                return 0
        if len(logger['true_range'])>=2:
            if logger['true_range'][-1] > logger['true_range'][-2] and (abs(true_range) <= 4 * interception_range):
                print('eeee')
                return 0


        # Calculate accelerations
        M_accel_y = +accel_m * np.cos(missile['HDG'])
        M_accel_x = -accel_m * np.sin(missile['HDG'])

        # Update velocities
        missile['yV'] += M_accel_y * dt 
        missile['xV'] += M_accel_x * dt

        # Update positions
        if true_range < 100:
            noise = 0.1
        else:
            noise = 0.3
        
        missile['y'] += missile['yV'] * dt + noise * np.random.randn()
        missile['x'] += missile['xV'] * dt + noise * np.random.randn()


        # 로깅
        logger['true_range'].append(true_range)
    return 0


# 격자점 생성
sat_range = np.arange(30, 90 + 10,10)
amf_range = np.arange(-180, 0, 10)

print(sat_range)
print(amf_range)

success_count = np.zeros((len(sat_range), len(amf_range)))
success_probabilities = []
total_runs_per_grid = 1

# 시그마 desired 변화에 따른 시뮬레이션
for sigma_desired_deg in [3, 6]:  # 5도와 10도에 대해 각각 시뮬레이션 실행
    for i, sat in enumerate(sat_range):
        for j, amf in enumerate(amf_range):
            for _ in range(total_runs_per_grid):  # 각 격자점에서 total_runs_per_grid만큼의 시뮬레이션 실행
                success_count[i, j] += run_simulation(sat, amf, sigma_desired_deg)


    # 성공 확률 계산
    success_probability = success_count / total_runs_per_grid    # total_runs는 시뮬레이션 횟수
    
    # 시그마 desired에 대한 성공 확률 저장
    success_probabilities.append(success_probability)





print(success_probabilities)

# 한 번에 플로팅
plt.figure(figsize=(12, 6))
for idx, sigma_desired_deg in enumerate([5, 10]):
    plt.subplot(1, 2, idx+1)
    plt.imshow(success_probabilities[idx], extent=[amf_range[0], amf_range[-1], 
                                                   sat_range[-1], sat_range[0]],
               cmap='viridis', origin='lower', aspect='auto')
    plt.colorbar(label='Success Probability')
    plt.xlabel('Alpha_M_F (degrees)')
    plt.ylabel('Saturation (degrees)')
    plt.title(f'Success Probability Heatmap (Sigma Desired: {sigma_desired_deg} degrees)')
    plt.grid(True)

plt.tight_layout()
plt.show()