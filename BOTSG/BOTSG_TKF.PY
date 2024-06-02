import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from SysM import sysm_tkf
from MeasM import measm_tkf


message = None

sat_upper = np.deg2rad(40)                  # lim
sat_lower = -sat_upper

am0 = np.deg2rad(+30)                       # Alpha_M_0: Initial M HDG: Fixed constant.
sigma_0 = am0                               # rad, Seeker's Initial Look Angle 
sigma_desired = np.deg2rad(+50)             # sat_upper + np.deg2rad(5)  # rad, Seeker's Look Angle Limit


amf = np.deg2rad(-80)       # rad, Desired Impact Angle

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


sigma_v = 0.0052;   # rad, ~= 0.298deg, Measurement Noise Standard Deviation
sigma_w = 0.1;  

# Verifying Guidance Gains N=(k1- k2*lambda): N is defined in the PPN loop
beta = 4*(sigma_desired - am0)/amf;                                  # beta is less or equal than 0
k1 = (beta - np.sqrt(beta**2 + 4 * beta * am0 / amf))/2 + 1;    # Guidance gain k1
k2 = (amf*(beta-np.sqrt(beta**2+4*beta*am0/amf))+2*am0)/amf**2;  # Guidance gain k2


# ================ Desired Impact Angle Achievement Survey ================
if sigma_desired < am0:
    print('Check initial am0 & sig_max conditions: Initial M.HDG is above the FOV limit. k1 & k2 become complex.')
    raise SystemExit

# if amf < -2 * am0:
#     if sigma_desired < -amf ** 2 / (4 * (amf + am0)):
#         print('Check initial conditions. Out of feasible solution area. k1 & k2 become complex.')
#         raise SystemExit
# # else:
#     if sigma_desired < am0:
#         print('Initial M.HDG > FOV limit. Out of feasible solution area. k1 & k2 become complex.')
#         raise SystemExit
# =========================================================================

S = 100                 # 최대 시뮬레이션 지속 시간, 초
dt = 0.01               # 시간 단계 크기, 초
Niter = int(S / dt)     # 최대 반복 횟수

logger = {
    't': [],
    'my': [],
    'mx': [],
    'myV': [],
    'mxV': [],
    'm_HDG': [],
    'z_sigma': [],  

    'ty': [],
    'tx': [],
    'tyV': [],
    'txV': [],

    'true_lambda': [],
    'true_lambda_dot': [],
    'true_sigma': [],
    'true_range': [],    
    
    'E_lambda': [],
    'E_lambda_dot': [],
    'E_sigma': [],
    'E_range': [], 

    'TKF_range_error': [],  # 실제 거리 - 예측 거리
    'TKF_look_error': [],

    'M_Accel': [],
    'M_Accel_x': [],
    'M_Accel_y': [],
    
    'xpp': [],

    'bias': [],
    'Ez_total': [],
    'Ez_l': [],
    'Ez_h': [],
    'Ez_us': [],

    'pdf_h' :[],
    'pdf_l' :[],
    'pl': [],
    'ph': [],
    'pus': [],

    'Cov_P1': [],
    'Cov_P2': [],
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

    # Update heading
    missile['HDG'] = np.arctan2(missile['yV'], missile['xV']) # +  np.deg2rad(sigma_w * np.random.randn())
        
    RTP_true = np.array([target['x'] - missile['x'], target['y'] - missile['y']])           # True RTP (Range = LOS distance)
    true_range = np.linalg.norm(RTP_true) #+ 0.1 * np.random.randn()
    Vm = np.sqrt(missile['yV']**2 + missile['xV']**2)                                       # Missile velocity, m/s
    true_lambda = np.arctan2(RTP_true[1], RTP_true[0])                                      # Line-of-sight (LOS) angle, rad w.r.t Ref.
    true_lambda_dot = - missile['V'] * np.sin(missile['HDG'] - true_lambda) / true_range    # LOS angle time derivatives (d/dt lambda), rad/s
    true_sigma = missile['HDG'] - true_lambda                                               # True Look Angle, rad


    # ================================================== TOBIT KALMAN FILTER ====================================================
    z, v = measm_tkf(sigma_0, true_sigma, sat_upper, sat_lower,k)  # MEASUREMENT MODEL

    # SYSTEM MODEL & GUIDANCE FILTER
    E_range, E_sigma, Cov_P, Ez_total, pdf_l, pdf_h, pl, ph, pus, bias_lambda, Ez_us, RR, HH, K, xpp = sysm_tkf(z, dt, missile, accel_m, sat_upper, sat_lower,k)
    bias = bias_lambda * 0.0052  # bias = bias_lambda * "Measurement Noise Deviation(sigma_v)" 0.0052
    logger['bias'].append(bias)

    logger['Ez_total'].append(Ez_total)
    logger['Ez_us'].append(Ez_us)
    logger['Ez_l'].append(pl * sat_lower)
    logger['Ez_h'].append(ph * sat_upper)

    logger['pus'].append(pus)
    logger['ph'].append(ph)
    logger['pl'].append(pl)
    logger['pdf_h'].append(pdf_h)
    logger['pdf_l'].append(pdf_l)

    logger['Cov_P1'].append(Cov_P[0, 0])
    logger['Cov_P2'].append(Cov_P[1, 1])
    

    TKF_range_error = true_range - E_range                      # Error btw true_range & E_range
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
    if RTP_true[1]> -2 and abs(RTP_true[0])<=interception_range:
        impact_angle_error = abs(np.rad2deg(amf - missile['HDG']))
        if impact_angle_error < 2.5:
            message = (f"Target intercepted successfully with desired impact angle.\n"
                f"Range is less than {interception_range:.2f}m. Final range: {true_range:.2f}m.\n"
                f"Missile heading error: {impact_angle_error:.2f} degrees.")
        else:
            message = (f"Target intercepted, but without the desired impact angle.\n"
                f"Range is less than {interception_range:.2f}m. Final range: {true_range:.2f}m.\n"
                f"Missile heading error: {impact_angle_error:.2f} degrees.")
        print(message)
        break

    if len(logger['true_range'])>=2:
        if logger['true_range'][-1] > logger['true_range'][-2] and (abs(true_range) <= 4 * interception_range):
            message = ("Missile missed the target. Range is increasing, indicating divergence.")
            print(message)
            break



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
    logger['t'].append(k*dt)
    logger['my'].append(missile['y'])
    logger['mx'].append(missile['x'])
    logger['myV'].append(missile['yV'])
    logger['mxV'].append(missile['xV'])
    logger['m_HDG'].append(missile['HDG'])

    logger['ty'].append(target['y'])
    logger['tx'].append(target['x'])
    logger['tyV'].append(target['yV'])
    logger['txV'].append(target['xV'])
    logger['z_sigma'].append(z)

    logger['true_range'].append(true_range)
    logger['true_lambda'].append(true_lambda)
    logger['true_lambda_dot'].append(true_lambda_dot)
    logger['true_sigma'].append(true_sigma)

    logger['E_lambda'].append(E_lambda)
    logger['E_lambda_dot'].append(E_lambda_dot)
    logger['E_sigma'].append(E_sigma)
    logger['E_range'].append(E_range)

    logger['TKF_range_error'].append(TKF_range_error)
    logger['TKF_look_error'].append(true_sigma - E_sigma)

    logger['M_Accel'].append(accel_m)
    logger['M_Accel_x'].append(M_accel_x)
    logger['M_Accel_y'].append(M_accel_y)

    logger['xpp'].append(xpp)





#--------------------------------------------------------------------------
# Sim result
#--------------------------------------------------------------------------
# Impact index
M, I = np.min(logger['E_range']), np.argmin(logger['E_range'])
print(f' M : {M} ')
print(f' I : {I} ')


fig, axs = plt.subplots(2, 3, figsize=(14, 10))
fig.suptitle(message)
# Position
axs[0, 0].plot(logger['mx'], logger['my'], label='Missile')
axs[0, 0].scatter(logger['tx'], logger['ty'], label='Target', c ='r',s=2)
axs[0, 0].set_title('Missile & Target trajectory')
axs[0, 0].legend()
axs[0, 0].set_xlabel('+x (m)')
axs[0, 0].set_ylabel('+y (m)')
# axs[0, 0].grid(True)

# mx_array = np.array(logger['mx'])
# my_array = np.array(logger['my'])
# tx_array = np.array(logger['tx'])
# ty_array = np.array(logger['ty'])
# e_range = np.array(logger['E_range'])

# # 희석 인덱스 생성
# sparse_index = np.arange(0, len(mx_array), 100)
# if sparse_index[-1] != len(mx_array) - 1:
#     sparse_index = np.append(sparse_index, len(mx_array) - 1)

# # 데이터 희석
# sparse_mx = mx_array[sparse_index]
# sparse_my = my_array[sparse_index]
# sparse_tx = tx_array[sparse_index]
# sparse_ty = ty_array[sparse_index]
# sparse_e_range = e_range[sparse_index]

# # 충돌 시점 찾기 (희석된 데이터 사용)
# collision_index = np.where(sparse_e_range <= 3)[0]
# if collision_index.size > 0:
#     last_frame = collision_index[0] 
# else:
#     last_frame = len(sparse_e_range)


# # 애니메이션 업데이트 함수
# def update(frame):
#     axs[0, 0].clear()
#     axs[0, 0].plot(sparse_mx[:frame], sparse_my[:frame], 'r', label='Missile')
#     #axs[0, 0].plot(sparse_tx[:frame], sparse_ty[:frame], 'b', label='Target')
#     axs[0, 0].plot(sparse_mx[frame-1], sparse_my[frame-1], 'ro')
#     axs[0, 0].plot(sparse_tx[frame-1], sparse_ty[frame-1], 'bo', label='Target')

#     axs[0, 0].legend()
#     axs[0, 0].set_title('Missile & Target Trajectory')
#     axs[0, 0].set_xlabel('+x (m)')
#     axs[0, 0].set_ylabel('+y (m)')
#     axs[0, 0].grid(True)

# ani = FuncAnimation(fig, update, frames=range(last_frame + 1), interval=50, repeat = False)


# Estimated, True and Measured Sigma(Look angle)
axs[0, 1].plot(logger['t'], np.rad2deg(logger['E_sigma']), linewidth=0.8, color='steelblue', label='Estimated Look Angle')
axs[0, 1].plot(logger['t'], np.rad2deg(logger['true_sigma']), linewidth=1, color='gold', label='True Look Angle')
axs[0, 1].scatter(logger['t'], np.rad2deg(logger['z_sigma']), s=0.8, c='orange', label='Measured Look Angle')
axs[0, 1].plot(logger['t'], np.rad2deg(logger['xpp']), label='Predicted Look Angle')
axs[0, 1].set_title('Estimated and True Look Angle')
axs[0, 1].legend()
axs[0, 1].set_xlabel('Elapsed time (sec)')
axs[0, 1].set_ylabel('Look Angle (deg)')
# axs[0, 1].grid(True)

# Look Angle Estimation Error & 3-sigma(P Mtx) bound
axs[0, 2].plot(logger['t'], np.rad2deg(logger['TKF_look_error']), linewidth=1.5, label='TKF Estimation Error')
axs[0, 2].plot(logger['t'], 3*np.rad2deg(np.sqrt(logger['Cov_P1'][1:])), linewidth=2, color='r', label='3-sigma bound')
axs[0, 2].plot(logger['t'], -3*np.rad2deg(np.sqrt(logger['Cov_P1'][1:])), linewidth=2, color='r')
axs[0, 2].set_title('Look Angle Error')
axs[0, 2].legend()
axs[0, 2].set_xlabel('Time (s)')
axs[0, 2].set_ylabel('deg')
axs[0, 2].grid(True)

# Range Error between Estimation and True Range
ax1 = axs[1, 0]
t_index_value = logger['t'][I]
ax1.plot(logger['t'], logger['E_range'], label='Estimated Range')
ax1.plot(logger['t'], logger['true_range'], label='True Range')
ax1.scatter(logger['t'][I], M, color='red', label=f'Intercept: r = {M:.2F} m, t = {t_index_value:.2F} s', s=20)
ax1.set_ylabel('Range (m)')
ax1.set_xlabel('Elapsed time (sec)')
ax1.set_title(f'Estimated and True Range Error')
# 두 번째 y축 설정
ax2 = ax1.twinx()  # 두 번째 y축을 위한 설정
ax2.scatter(logger['t'], logger['TKF_range_error'], color='blue', label='Range Error', s=0.8)
# ax2.plot(logger['t'], 3*np.sqrt(logger['Cov_P2'][1:]), linewidth=2, color='r', label='3-sigma bound')
# ax2.plot(logger['t'], -3*np.sqrt(logger['Cov_P2'][1:]), linewidth=2, color='r')
ax2.set_ylabel('Range Error (m)')

# 레전드 핸들과 레이블 수집
handles, labels = [], []
for ax in [ax1, ax2]:
    for handle, label in zip(*ax.get_legend_handles_labels()):
        handles.append(handle)
        labels.append(label)
plt.legend(handles, labels)



# print(logger['m_HDG'])
axs[1, 1].plot(logger['t'],  np.rad2deg(logger['m_HDG']), label='Missile Heading')
axs[1, 1].plot(logger['t'],  np.rad2deg(logger['true_lambda']),  label='LOS Angle')
axs[1, 1].plot(logger['t'],  np.rad2deg(logger['true_sigma']), label='Look Angle')
axs[1, 1].legend()
axs[1, 1].set_title('Missile HDG, LOS Angle, Look Angel')
axs[1, 1].set_xlabel('Elapsed time(sec)')
axs[1, 1].set_ylabel('Degree')
axs[1, 1].grid(True)


# Acceleration
axs[1, 2].plot(logger['t'],  logger['M_Accel'])
axs[1, 2].set_title(f'Missile Acceleration')
axs[1, 2].set_xlabel('Elapsed time(sec)')
axs[1, 2].set_ylabel('Missile Acceleration (m/s^2)')
axs[1, 2].grid(True)
plt.tight_layout()
plt.show()





# # Pus // pl // ph
# axs[1, 1].plot(logger['t'],  logger['pl'][1:], label='pl')
# axs[1, 1].plot(logger['t'],  logger['ph'][1:],  label='ph')
# axs[1, 1].plot(logger['t'],   logger['pus'][1:], label='pus')
# axs[1, 1].set_title('Probabilities')
# axs[1, 1].legend()
# axs[1, 1].set_xlabel('Elapsed time(sec)')
# axs[1, 1].set_ylabel('Probability')


# # Bias
# converted_bias = np.rad2deg(logger['bias'][1:])
# axs[1, 2].plot(logger['t'],  converted_bias.flatten(), label='bias')
# axs[1, 2].set_title('Bias Term')
# axs[1, 2].legend()
# axs[1, 2].set_xlabel('Elapsed time(sec)')
# axs[1, 2].set_ylabel('Bias Value')


# # pdf_l // pdf_h
# axs[1, 2].plot(logger['t'],  logger['pdf_l'][1:], label='pdf_l')
# axs[1, 2].plot(logger['t'],  logger['pdf_h'][1:], label='pdf_h')
# axs[1, 2].set_title('PDF_H // PDF_L')
# axs[1, 2].legend()
# axs[1, 2].set_xlabel('Elapsed time(sec)')
# axs[1, 2].set_ylabel('PDF')
# plt.tight_layout()
# plt.show()

# # Acceleration
# plt.figure()
# plt.plot(logger['t'], logger['M_Accel'], label='M_Accel')
# plt.title(f'IACCG(DPP+PPN)')
# plt.legend('Missile Acceleration')
# plt.ylabel('Missile Acceleration (m/s^2)')
# plt.xlabel('Elapsed time (sec)')
# plt.tight_layout()
# plt.show()

# # Bias
# plt.figure()
# converted_bias = np.rad2deg(logger['bias'][1:])
# plt.plot(logger['t'],  converted_bias.flatten(), label='bias')
# plt.title('Bias Term')
# plt.legend()
# plt.xlabel('Elapsed time(sec)')
# plt.ylabel('Bias Value')
# plt.tight_layout()
# plt.show()



plt.figure()
fig, ax1 = plt.subplots()

p1, = ax1.plot(logger['t'], logger['E_range'], label='Estimated Range')
p2, = ax1.plot(logger['t'], logger['true_range'], linewidth=2, label='True Range')
ax1.set_ylabel('Range (m)')
ax1.set_xlabel('Elapsed time (sec)')

# 두 번째 y축 (오른쪽)
ax2 = ax1.twinx()
p3, = ax2.plot(logger['t'], np.rad2deg(logger['E_sigma']), linewidth=0.8, color='steelblue', label='Estimated Look Angle')
p4, = ax2.plot(logger['t'], np.rad2deg(logger['true_sigma']), linewidth=1, color='gold', label='True Look Angle')
p5 = ax2.scatter(logger['t'], np.rad2deg(logger['z_sigma']), s=0.8, c='orange', label='Measured Look Angle')
p6, = ax2.plot(logger['t'], np.rad2deg(logger['xpp']), label='Predicted Look Angle')
ax2.set_ylabel('Look Angle (degrees)')

# 모든 레전드를 하나로 합침
lines = [p1, p2, p3, p4, p5, p6]
labels = [l.get_label() for l in lines[:4]] + ['Measured Look Angle', 'Predicted Look Angle']
ax1.legend(lines, labels)

plt.title('Estimated and True Look Angle // Range')
plt.tight_layout()
plt.show()




# axs[0, 1].plot(logger['t'], np.rad2deg(logger['E_sigma']), linewidth=0.8, color='steelblue', label='Estimated Look Angle')
# axs[0, 1].plot(logger['t'], np.rad2deg(logger['true_sigma']), linewidth=1, color='gold', label='True Look Angle')
# axs[0, 1].scatter(logger['t'], np.rad2deg(logger['z_sigma']), s=0.8, c='orange', label='Measured Look Angle')
# axs[0, 1].plot(logger['t'], np.rad2deg(logger['xpp']), label='Predicted Look Angle')
# axs[0, 1].set_title('Estimated and True Look Angle')
# axs[0, 1].legend()
# axs[0, 1].set_xlabel('Elapsed time (sec)')
# axs[0, 1].set_ylabel('Look Angle (deg)')



# # 거리 오차
# # Error between Estimation and True Range
# plt.figure()
# plt.scatter(logger['t'], logger['TKF_range_error'], color='blue', label='Range Error')
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