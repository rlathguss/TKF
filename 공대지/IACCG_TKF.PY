import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from SysM import sysm_tkf
from MeasM import measm_tkf


# Navigation Constant
N = 3
index = 0

# Missile init params
missile = {'y': 7000*np.sin(np.deg2rad(27)),  # y pos, m
           'x': -7000*np.cos(np.deg2rad(27)),  # x pos, m
           'HDG': np.deg2rad(-22),
           'HDG_prev' : 0,
           'HDG_dot' : 0,
           'V': 200}  # velocity magnitude, m/s
missile['yV'] = missile['V'] * np.sin(missile['HDG'])  # y velocity, m/s
missile['xV'] = missile['V'] * np.cos(missile['HDG'])  # x velocity, m/s


# Target init params
target = {'y': 0,  # y pos, m
          'x': 0,  # x pos, m
          'HDG': np.deg2rad(0),  # heading, rad
          'yV': 0,  # y velocity, m/s
          'xV': 0}  # x velocity, m/s
target['V'] = np.sqrt(target['yV']**2 + target['xV']**2)  # velocity magnitude, m/s


# Initial Look Angle
sigma_0 = np.deg2rad(5)     # rad
sigma_max = np.deg2rad(60)  # rad, Seeker's Look Angle Limit
accel_m = 0                 # Initial Missile Acceleration, m/s^2
accel_max = 100             # m/s^2, acceleration limit

sat_upper = np.deg2rad(7.5)
sat_lower = -sat_upper

# Other Initial Conditions and Parameters
eta = target['V'] / missile['V']

# sigma_d Selection in order to bring out lambda_switch and also compare whether initial setup is good or not.
sigma_d_lower = np.arcsin(eta)  # rad, Lower bound for sigma_d
sigma_d_upper = min(np.arccos(2*eta), sigma_max)  # rad, Upper bound for sigma_d
sigma_d = np.deg2rad(10)  # rad, Desired Look Angle / Commanded Look Angle



if not ((abs(sigma_d) > sigma_d_lower) and (abs(sigma_d) <= sigma_d_upper)):
    print('Check sigma_d initial condition. Must be within a range.')
    exit()

if sigma_0 > abs(sigma_d):
    print('sigma_0 must be equal or less than sigma_d')
    exit()

zeta1 = eta / np.sin(sigma_d)
zeta2 = eta / np.cos(sigma_d)


amf = np.deg2rad(-50)  # rad, Desired Impact Angle
lambda_switch = ((np.arctan2(np.sin(amf) - eta * np.sin(target['HDG']), np.cos(amf) - eta * np.cos(target['HDG']))) - ((amf) - (sigma_d)) / N) * (N / (N - 1))


# K Gain Selection
K = 230  # K 값을 원하는 대로 선택하세요.
# R0 초기화를 위해 target과 missile의 위치 차이의 노름을 계산합니다.
R0 = np.linalg.norm([target['x'] - missile['x'], target['y'] - missile['y']])

# 상한 K 값을 계산합니다.
K_upper = (accel_max * R0 + missile['V']**2 * np.sin(sigma_0)) / (R0 * (abs(sigma_d) - sigma_0))

if K > K_upper:
    print('K must be less than K_upper')
else:
    # 계속 진행할 코드를 여기에 작성하세요.
    pass  # 'pass'를 사용하여 빈 블록을 표시할 수 있습니다.

# Sim params
S = 100  # 최대 시뮬레이션 지속 시간, 초
dt = 0.01  # 시간 단계 크기, 초
Niter = int(S / dt)  # 최대 반복 횟수

logger = {
    't': [],
    'my': [],
    'mx': [],
    'myV': [],
    'mxV': [],
    'm_HDG': [],
    'm_HDG_dot': [],
    'ty': [],
    'tx': [],
    'tyV': [],
    'txV': [],
    'z_sigma': [],
    'v': [],
    'E_range': [],
    'true_range': [],
    'est_error': [],
    'lambda': [],
    'lambda_dot': [],
    'true_lambda': [],
    'true_lambda_dot': [],
    'sigma': [],
    'true_sigma': [],
    'E_sigma': [],
    'TKF_lookerror': [],
    'EKF_lookerror': [],
    'M_Accel': [],
    'M_Accel_x': [],
    'M_Accel_y': [],
    'normKalman': [],
    'count': [],
    'Ez_total': [],
    'pdf_h' :[],
    'pdf_l' :[],
    'pl': [],
    'ph': [],
    'pus': [],
    'bias': [],
    'RR': [],
    'HH': [],
    'R': [],
    'H': [],
    'Ez_l': [],
    'Ez_h': [],
    'Ez_us': [],
    'Ez_total_nobias': [],
    'KK': [],
    'K_EKF': [],
    'x_EKF': [],
    'xpp': [],
    'E_sigma_nobias': [],
    'HH_nobias': [],
    'RR_nobias': [],
    'K_nobias': [],
    'tkf_rhsterm': [],
    'ekf_rhsterm': [],
    'Cov_P': [],
}


#--------------------------------------------------------------------------
# Init sim
#--------------------------------------------------------------------------
# delta y, delta x
RTP_last = np.array([target['x'] - missile['x'], target['y'] - missile['y']])

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
    # SELECT ONE FILTER AMONG THREE BELOW : TKF || EKF w/ SAT as input || EKF w/o SAT. as input
    # ================================================== TOBIT KALMAN FILTER ====================================================
    z, v, count = measm_tkf(sigma_0, true_sigma, sat_upper, sat_lower,k)  # MEASUREMENT MODEL
    logger['count'].append(count)
    
    # SYSTEM MODEL & GUIDANCE FILTER
    E_range, E_sigma, Cov_P, Ez_total, pdf_l, pdf_h, pl, ph, pus, bias_lambda, RR, HH, Ez_us, KK, K_EKF, x_EKF, xpp, x_nobias, HH_nobias, RR_nobias, Ez_total_nobias, K_nobias = sysm_tkf(z, dt, missile, accel_m, count, sat_upper, sat_lower,k)
    bias = bias_lambda * 0.0052  # bias = bias_lambda * "Measurement Noise Deviation(sigma_v)" 0.0052
    logger['bias'].append(bias)
    logger['RR'].append(np.linalg.norm(RR))
    logger['HH'].append(np.linalg.norm(HH))
    logger['R'].append(0.0052**2)  # MODIFY: sigma_v
    logger['H'].append(1)  # MODIFY: norm(H)
    logger['pus'].append(pus)
    logger['ph'].append(ph)
    logger['pl'].append(pl)
    logger['pdf_h'].append(pdf_h)
    logger['pdf_l'].append(pdf_l)

    logger['Ez_us'].append(Ez_us)
    logger['Ez_l'].append(pl * sat_lower)
    logger['Ez_h'].append(ph * sat_upper)
    logger['Ez_total'].append(Ez_total)
    logger['Ez_total_nobias'].append(Ez_total_nobias)
    logger['KK'].append(KK[0])
    logger['K_EKF'].append(K_EKF[0])
    logger['x_EKF'].append(x_EKF)
    logger['E_sigma_nobias'].append(x_nobias[0])
    logger['RR_nobias'].append(np.linalg.norm(RR_nobias))
    logger['HH_nobias'].append(np.linalg.norm(HH_nobias))
    logger['K_nobias'].append(K_nobias[0])
    logger['tkf_rhsterm'].append(KK[0] * (z - Ez_total))
    logger['ekf_rhsterm'].append(K_EKF[0] * (z - xpp))
    logger['Cov_P'].append(Cov_P[0, 0])
    
    if pus == 0:
        print('UNSAT Probability is 0. All the parameters become NaN from now on.')
        break

    # True RTP (Range = LOS distance)
    RTP_true = np.array([target['x'] - missile['x'], target['y'] - missile['y']])
    true_range = np.linalg.norm(RTP_true)
    est_error = true_range - E_range                        # Error btw true_range & E_range
    Vm = np.sqrt(missile['yV']**2 + missile['xV']**2)       # Missile velocity, m/s
    
    # Line-of-sight (LOS) angle, rad w.r.t Ref.
    true_lambda = np.arctan2(RTP_true[1], RTP_true[0])
    lambda_ = missile['HDG'] - E_sigma
    
    # LOS angle time derivatives (d/dt lambda), rad/s
    lambda_dot = - missile['V'] * np.sin(E_sigma) / E_range
    true_lambda_dot = - missile['V'] * np.sin(missile['HDG'] - true_lambda) / true_range
    
    # True Look Angle, rad
    true_sigma = missile['HDG'] - true_lambda
    
    # Missile Heading Angle time derivative (d/dt missile.HDG), rad/s
    if k == 0:
        missile_HDG_dot = 0
    else:
        missile_HDG_dot = (missile['HDG'] - missile['HDG_prev']) / dt
    K = 300

    # Missile Acceleration 계산
    accel_m = Vm * lambda_dot + K * (sigma_d - E_sigma)

    if abs(lambda_) >= abs(lambda_switch) or index:
        accel_m = N * Vm * lambda_dot
        index = 1

    # Acceleration Saturation Condition
    if abs(accel_m) > 100:
        if accel_m < 0:
            accel_m = -100
        elif accel_m > 0:
            accel_m = 100
        else:
            pass

    # Terminate sim at interception
    if abs(E_range) <= 30:  # interception 기준
        print('Estimated Range is less than 30m')
        break

    # Missile 위치 업데이트
    missile['y'] += missile['yV'] * dt + 1 * np.random.randn()
    missile['x'] += missile['xV'] * dt + 1 * np.random.randn()

    M_accel_y = +accel_m * np.cos(missile['HDG'])
    M_accel_x = -accel_m * np.sin(missile['HDG'])

    # Missile y/x 속도 업데이트
    missile['yV'] += M_accel_y * dt
    missile['xV'] += M_accel_x * dt

    # Missile 헤딩 업데이트
    missile['HDG_prev'] = missile['HDG']
    missile['HDG'] = np.arctan2(missile['yV'], missile['xV'])

    # Target 위치 업데이트는 주석 처리되어 있으므로 여기서는 생략합니다.

    # 로깅
    logger['t'].append(k*dt)
    logger['my'].append(missile['y'])
    logger['mx'].append(missile['x'])
    logger['myV'].append(missile['yV'])
    logger['mxV'].append(missile['xV'])
    logger['m_HDG'].append(missile['HDG'])
    logger['m_HDG_dot'].append(missile['HDG_dot'])

    logger['ty'].append(target['y'])
    logger['tx'].append(target['x'])
    logger['tyV'].append(target['yV'])
    logger['txV'].append(target['xV'])

    logger['E_range'].append(E_range)
    logger['true_range'].append(true_range)

    logger['est_error'].append(est_error)

    logger['lambda'].append(lambda_)
    logger['lambda_dot'].append(lambda_dot)
    logger['true_lambda'].append(true_lambda)
    logger['true_lambda_dot'].append(true_lambda_dot)
    logger['E_sigma'].append(E_sigma)
    logger['true_sigma'].append(true_sigma)
    logger['TKF_lookerror'].append(true_sigma - E_sigma)
    logger['z_sigma'].append(z)
    logger['xpp'].append(xpp)

    logger['M_Accel'].append(accel_m)
    logger['M_Accel_x'].append(M_accel_x)
    logger['M_Accel_y'].append(M_accel_y)





#--------------------------------------------------------------------------
# Sim result
#--------------------------------------------------------------------------
# Impact index
M, I = np.min(logger['E_range']), np.argmin(logger['E_range'])
print(f' M : {M} ')
print(f' I : {I} ')


fig, axs = plt.subplots(2, 3, figsize=(15, 10))
# Position
# axs[0, 0].scatter(logger['mx'], logger['my'], label='Missile', s=0.7)
# axs[0, 0].scatter(logger['tx'], logger['ty'], label='Target', s=2)
# axs[0, 0].set_title('M&T Trajectory')
# axs[0, 0].legend()
# axs[0, 0].set_xlabel('+x (m)')
# axs[0, 0].set_ylabel('+y (m)')

mx_array = np.array(logger['mx'])
my_array = np.array(logger['my'])
tx_array = np.array(logger['tx'])
ty_array = np.array(logger['ty'])
e_range = np.array(logger['E_range'])

# 희석 인덱스 생성
sparse_index = np.arange(0, len(mx_array), 100)
if sparse_index[-1] != len(mx_array) - 1:
    sparse_index = np.append(sparse_index, len(mx_array) - 1)

# 데이터 희석
sparse_mx = mx_array[sparse_index]
sparse_my = my_array[sparse_index]
sparse_tx = tx_array[sparse_index]
sparse_ty = ty_array[sparse_index]
sparse_e_range = e_range[sparse_index]

# 충돌 시점 찾기 (희석된 데이터 사용)
collision_index = np.where(sparse_e_range <= 3)[0]
if collision_index.size > 0:
    last_frame = collision_index[0] 
else:
    last_frame = len(sparse_e_range)


# 애니메이션 업데이트 함수
def update(frame):
    axs[0, 0].clear()
    axs[0, 0].plot(sparse_mx[:frame], sparse_my[:frame], 'r', label='Missile')
    # axs[0, 0].plot(sparse_tx[:frame], sparse_ty[:frame], 'b', label='Target')
    axs[0, 0].plot(sparse_mx[frame-1], sparse_my[frame-1], 'ro')
    axs[0, 0].plot(sparse_tx[frame-1], sparse_ty[frame-1], 'bo', label='Target')

    axs[0, 0].legend()
    axs[0, 0].set_title('M&T Trajectory')
    axs[0, 0].set_xlabel('+x (m)')
    axs[0, 0].set_ylabel('+y (m)')

ani = FuncAnimation(fig, update, frames=range(last_frame + 1), interval=50, repeat = False)


# Estimated, True and Measured Sigma(Look angle)
axs[0, 1].plot(logger['t'], np.rad2deg(logger['E_sigma']), linewidth=0.8, color='steelblue', label='Estimated Look Angle')
axs[0, 1].plot(logger['t'], np.rad2deg(logger['true_sigma']), linewidth=1, color='gold', label='True Look Angle')
axs[0, 1].scatter(logger['t'], np.rad2deg(logger['z_sigma']), s=0.8, c='orange', label='Measured Look Angle')
axs[0, 1].plot(logger['t'], np.rad2deg(logger['xpp']), label='Predicted Look Angle')
axs[0, 1].set_title('Estimated and True Look Angle')
axs[0, 1].legend()
axs[0, 1].set_xlabel('Elapsed time (sec)')
axs[0, 1].set_ylabel('Look Angle (deg)')

# Look Angle Estimation Error & 3-sigma(P Mtx) bound
axs[0, 2].plot(logger['t'], np.rad2deg(logger['TKF_lookerror']), linewidth=1.5, label='TKF Estimation Error')
axs[0, 2].plot(logger['t'], 3*np.rad2deg(np.sqrt(logger['Cov_P'][1:])), linewidth=2, color='r', label='3-sigma bound')
axs[0, 2].plot(logger['t'], -3*np.rad2deg(np.sqrt(logger['Cov_P'][1:])), linewidth=2, color='r')
axs[0, 2].set_title('Look Angle Error Trajectory')
axs[0, 2].legend()
axs[0, 2].set_xlabel('Time (s)')
axs[0, 2].set_ylabel('deg')

# Range Error between Estimation and True Range
ax1 = axs[1, 0]
t_index_value = logger['t'][I]
ax1.plot(logger['t'], logger['E_range'], label='Estimated Range')
ax1.plot(logger['t'], logger['true_range'], linewidth=2, label='True Range')
ax1.scatter(logger['t'][I], M, color='red', label=f'Intercept: r = {M} m, t = {t_index_value} s', s=20)
ax1.set_ylabel('Range (m)')
ax1.set_xlabel('Elapsed time (sec)')
ax1.set_title(f'IACCG(DPP&PPN), N = {N}')
# 두 번째 y축 설정
ax2 = ax1.twinx()  # 두 번째 y축을 위한 설정
ax2.scatter(logger['t'], logger['est_error'], color='blue', label='Range Error', s=0.8)
ax2.set_ylabel('Range Error btw Estimated and True Range [m])')

# 레전드 핸들과 레이블 수집
handles, labels = [], []
for ax in [ax1, ax2]:
    for handle, label in zip(*ax.get_legend_handles_labels()):
        handles.append(handle)
        labels.append(label)
plt.legend(handles, labels)

# Pus // pl // ph
axs[1, 1].plot(logger['t'],  logger['pl'][1:], label='pl')
axs[1, 1].plot(logger['t'],  logger['ph'][1:],  label='ph')
axs[1, 1].plot(logger['t'],   logger['pus'][1:], label='pus')
axs[1, 1].set_title('Probabilities')
axs[1, 1].legend()
axs[1, 1].set_xlabel('Elapsed time(sec)')
axs[1, 1].set_ylabel('Probability')

# pdf_l // pdf_h
axs[1, 2].plot(logger['t'],  logger['pdf_l'][1:], label='pdf_l')
axs[1, 2].plot(logger['t'],  logger['pdf_h'][1:], label='pdf_h')
axs[1, 2].set_title('PDF_H // PDF_L')
axs[1, 2].legend()
axs[1, 2].set_xlabel('Elapsed time(sec)')
axs[1, 2].set_ylabel('PDF')
plt.tight_layout()
plt.show()

# Acceleration
plt.figure()
plt.plot(logger['t'], logger['M_Accel'], label='M_Accel')
plt.title(f'IACCG(DPP+PPN)')
plt.legend('Missile Acceleration')
plt.ylabel('Missile Acceleration (m/s^2)')
plt.xlabel('Elapsed time (sec)')
plt.tight_layout()
plt.show()




# # Bias
# converted_bias = np.rad2deg(logger['bias'][1:])
# axs[1, 0].plot(logger['t'],  converted_bias.flatten(), label='bias')
# axs[1, 0].set_title('Bias Term')
# axs[1, 0].legend()
# axs[1, 0].set_xlabel('Elapsed time(sec)')
# axs[1, 0].set_ylabel('Bias Value')


# # Range
# plt.figure()
# t_index_value = logger['t'][I]
# plt.plot(logger['t'], logger['E_range'], label='Estimated Range')
# plt.plot(logger['t'], logger['true_range'], linewidth=2, label='True Range')
# plt.scatter(logger['t'][I], M, color='red', label=f'Intercept: r = {M} m, t = {t_index_value} s')
# plt.title(f'IACCG(DPP&PPN), N = {N}')
# plt.legend(loc='upper left')
# plt.ylabel('Range (m)')
# plt.xlabel('Elapsed time (sec)')
# plt.gca().invert_yaxis()


# # 거리 오차
# # Error between Estimation and True Range
# plt.figure()
# plt.scatter(logger['t'], logger['est_error'], color='blue', label='Range Error')
# plt.title(f'IACCG(DPP&PPN), N = {N}')
# plt.legend()
# plt.ylabel('Range Error btw Estimated and True Range [m])')
# plt.xlabel('Elapsed time (sec)')


# # LOS ANGLE
# # Estimated and True Lambda (LOS angle)
# plt.figure()
# plt.plot(logger['t'], np.rad2deg(logger['lambda']), label='Estimated LOS angle')
# plt.plot(logger['t'], np.rad2deg(logger['true_lambda']), linewidth=2, label='True LOS angle')
# plt.title('Estimated and True LOS Angle')
# plt.legend()
# plt.xlabel('Elapsed time (sec)')
# plt.ylabel('LOS Angle (deg)')

# # LOS ANGLE DOT
# # Lambda dot True and Estimated
# plt.figure()
# plt.plot(logger['t'], np.rad2deg(logger['lambda_dot']), label='Estimated Lambda dot')
# plt.plot(logger['t'], np.rad2deg(logger['true_lambda_dot']), linewidth=2, label='True Lambda dot')
# plt.title(f'IACCG(DPP&PPN), N = {N}')
# plt.legend()
# plt.ylabel('Lambda dot (deg/s)')
# plt.xlabel('Elapsed time (sec)')