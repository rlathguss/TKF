import numpy as np
from scipy.stats import norm

x = None
P = None

def sysm_tkf(z, dt, missile, accel_m, count, sat_upper, sat_lower, firstRun):
    global x, P
    H = np.array([[1, 0]]) # 측정 행렬
    sigma_w = 0.1  # System(Process) noise
    sigma_v = 0.0052  # Measurement noise
    R = sigma_v ** 2  #측정 노이즈 공분산
    # 제어 입력 행렬
    B = np.array([[1/ missile['V'], 0]]).T

    if firstRun == 0:
        # 공대지 초기 추정값    
        x = np.array([[np.deg2rad(5)], [7000]])    
        
        # 오차 공분산 초기 설정
        P = np.array([[0.0052, 0],[0, 50]])

    # 상태 예측
    A = Ajacobian(x, missile, dt)  # Linearize
    xp = fxp(x, dt, B, accel_m, missile, sigma_w) # 상태예측
    xpp = xp[0] # 예측된 각도 
    
    Q = Qupdate(A, missile, dt, sigma_w)  # Process Noise Covariance Matrix 업데이트
    Pp = A @ P @ A.T + Q  # 오차 공분산 예측
    
    # Ez, HH, RR 계산
    Ez_total, HH, RR, pdf_l, pdf_h, pl, ph, pus, bias_lambda, Ez_us, HH_nobias, RR_nobias, Ez_total_nobias = zExpectation(H, xp, sat_upper, sat_lower, sigma_v, count)
    
    # 칼만 이득 계산 및 상태 업데이트
    K = Pp @ HH.T @ np.linalg.inv(HH @ Pp @ HH.T + RR)
    K_EKF = Pp @ H.T @ np.linalg.inv(H @ Pp @ H.T + R)
    
    x = xp + K @ (z - Ez_total)
    x_EKF = xp + K_EKF @ (z - H @ xp)
    x_EKF = x_EKF[0]  # Look Angle
    
    K_nobias = Pp @ HH_nobias.T @ np.linalg.inv(HH_nobias @ Pp @ HH_nobias.T + RR_nobias)
    x_nobias = xp + K_nobias @ (z - Ez_total_nobias)
    
    P = Pp - K @ HH @ Pp  # 오차 공분산 계산
    
    E_sigma = x[0,0]  # 추정된 Look Angle
    E_range = x[1,0]  # 추정된 범위

    return E_range, E_sigma, P, Ez_total, pdf_l, pdf_h, pl, ph, pus, bias_lambda, RR, HH, Ez_us, K, K_EKF, x_EKF, xpp, x_nobias, HH_nobias, RR_nobias, Ez_total_nobias, K_nobias







class SysM_TKF_StatT:
    def __init__(self):
        self.firstRun = True

    def run(self, z, dt, missile, accel_m, count, sat_upper, sat_lower):
        if self.firstRun:
            self.firstRun = False
            
            self.H = np.array([[1, 0]]) # 측정 행렬
            #self.x = np.array([[np.deg2rad(15)], [10000]])  # 지대지 초기 추정값
            self.x = np.array([[np.deg2rad(5)], [7000]])  # 공대지 초기 추정값

            self.sigma_w = 0.1  # System(Process) noise
            self.sigma_v = 0.0052  # Measurement noise
            self.R = self.sigma_v ** 2  #측정 노이즈 공분산
            
            # 오차 공분산 초기 설정
            self.P = np.array([[0.0052, 0],
                               [0, 50]])
            # 제어 입력 행렬
            self.B = np.array([[1/ missile['V'], 0]]).T

        # 상태 예측
        A = Ajacobian(self.x, missile, dt)  # Linearize
        xp = fxp(self.x, dt, self.B, accel_m, missile, self.sigma_w) # 상태예측
        xpp = xp[0] # 예측된 각도 
        
        Q = Qupdate(A, missile, dt, self.sigma_w)  # Process Noise Covariance Matrix 업데이트
        Pp = A @ self.P @ A.T + Q  # 오차 공분산 예측
        
        # Ez, HH, RR 계산
        # 측정값의 기대치 계산 및 측정 업데이트
        Ez_total, HH, RR, pdf_l, pdf_h, pl, ph, pus, bias_lambda, Ez_us, HH_nobias, RR_nobias, Ez_total_nobias = zExpectation(self.H, xp, sat_upper, sat_lower, self.sigma_v, count)
        
        # 칼만 이득 계산 및 상태 업데이트
        K = Pp @ HH.T @ np.linalg.inv(HH @ Pp @ HH.T + RR)
        K_EKF = Pp @ self.H.T @ np.linalg.inv(self.H @ Pp @ self.H.T + self.R)
        
        self.x = xp + K @ (z - Ez_total)
        x_EKF = xp + K_EKF @ (z - self.H @ xp)
        x_EKF = x_EKF[0]  # Look Angle
        
        K_nobias = Pp @ HH_nobias.T @ np.linalg.inv(HH_nobias @ Pp @ HH_nobias.T + RR_nobias)
        x_nobias = xp + K_nobias @ (z - Ez_total_nobias)
        
        self.P = Pp - K @ HH @ Pp  # 오류 공분산 계산
        
        E_sigma = self.x[0,0]  # 추정된 Look Angle
        E_range = self.x[1,0]  # 추정된 범위
        # print(E_range)
        return E_range, E_sigma, self.P, Ez_total, pdf_l, pdf_h, pl, ph, pus, bias_lambda, RR, HH, Ez_us, K, K_EKF, x_EKF, xpp, x_nobias, HH_nobias, RR_nobias, Ez_total_nobias, K_nobias



def Ajacobian(x, missile, dt):
    A = np.zeros((2, 2))
    A[0, 0] =  missile['V'] * np.cos(x[0]) / x[1]
    A[0, 1] = - missile['V'] * np.sin(x[0]) / x[1]**2
    A[1, 0] =  missile['V'] * np.sin(x[0])
    A[1, 1] = 0
    A = np.eye(2) + A * dt
    return A


def Qupdate(A, missile, dt, sigma_w):
    Q = np.zeros((2,2))
    Q[0, 0] = ((A[0, 0])**2 * dt**2) / 3 + dt * A[0, 0] + 1
    Q[0, 1] = (A[0, 0] * A[1, 0] * dt**2) / 3 + dt * A[1, 0] / 2
    Q[1, 0] = Q[0, 1]
    Q[1, 1] = dt**2 * (A[1, 0])**2 / 3
    Q = (Q * sigma_w**2 * dt) /  missile['V']**2
    return Q

def fxp(x, dt, B, accel_m, missile, sigma_w):
    xdot = np.zeros((2,1))
    xdot[0] = missile['V'] * np.sin(x[0]) / x[1]
    xdot[1] = -missile['V'] * np.cos(x[0])
    
    u = accel_m + np.random.randn() * sigma_w  # 평균 0, 표준편차 1 정규분포를 따르는 노이즈 생성
    xp = x + xdot * dt + B * u * dt

    
    return xp



def zExpectation(H, xp, sat_upper, sat_lower, sigma_v, count):

    # 상한과 하한까지의 정규화된 거리 계산
    alpha = (sat_lower - np.sum(H @ xp)) / sigma_v
    beta = (sat_upper - np.sum(H @ xp)) / sigma_v
    # alpha와 beta에서의 확률 밀도 함수(PDF) 값 계산
    pdf_l = norm.pdf(alpha)
    pdf_h = norm.pdf(beta)

    # alpha와 beta에서의 누적 분포 함수(CDF) 값 계산
    pl = norm.cdf(alpha)  # 하한 이하일 확률
    ph = 1 - norm.cdf(beta)  # 상한 이상일 확률

    if ph == 1:
        ph = 0.99999

    if pl ==1:
        pl = 0.99999

    # 측정값이 포화 상태가 아닐 확률
    pus = 1 - ph - pl

    # 포화 상태가 아닌 측정값에 대한 편향 보정 계수 계산
    bias_lambda = (pdf_h - pdf_l) / pus
    
    # 편향 보정을 포함한 측정값의 예상 값 계산
    Ez_total = (H @ xp - sigma_v * bias_lambda) * pus + sat_lower * pl + sat_upper * ph

    # 편향 보정 없이 측정값의 예상 값 계산
    Ez_total_nobias = H @ xp * pus + pl * sat_lower + ph * sat_upper
    
    # 측정 행렬과 잡음 공분산 행렬 업데이트
    HH = pus * H
    RR = sigma_v**2 * (1 - bias_lambda**2) + sigma_v * H @ xp * bias_lambda + sigma_v * (sat_lower * pdf_l - sat_upper * pdf_h) / pus
    HH_nobias = pus * H
    RR_nobias = sigma_v**2 + 0 + sigma_v * (sat_lower * pdf_l - sat_upper * pdf_h) / pus
    
    # 포화 상태가 아닌 측정값의 예상 값 계산
    Ez_us = H @ xp * pus - sigma_v * (pdf_h - pdf_l)
    
    # 계산된 값들 반환
    return Ez_total, HH, RR, pdf_l, pdf_h, pl, ph, pus, bias_lambda, Ez_us, HH_nobias, RR_nobias, Ez_total_nobias
