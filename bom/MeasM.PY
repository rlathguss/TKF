import numpy as np


def measm_tkf(sigma_0, true_sigma, sat_upper, sat_lower, firstRun):
    sigma_v = 0.0052  # rad, 공분산(sigma_v) = 0.0052, 3*공분산(sigma_v) = 99% 영역
    v = sigma_v * np.random.randn()  # 평균 0, 표준편차 1 정규분포를 따르는 노이즈 생성

    if firstRun == 0:
        z_sigma = sigma_0 + v
    else:
        z_sigma = true_sigma + v

    count = 1  # count는 측정값이 포화 상태인지 여부를 나타냄

    # 측정값이 상한값을 초과하는 경우
    if abs(z_sigma) >= sat_upper and z_sigma > 0:
        count = 0  # 포화 상태로 설정
        z_sigma = sat_upper  # 측정값을 상한값으로 조정
    # 측정값이 하한값 이하인 경우
    elif abs(z_sigma) >= sat_upper and z_sigma < 0:
        count = 0  # 포화 상태로 설정
        z_sigma = sat_lower  # 측정값을 하한값으로 조정

    return z_sigma, v, count


class MeasM_TKF_StatT:
    def __init__(self, sigma_0, true_sigma, sat_upper, sat_lower):
        self.sigma = sigma_0  # 초기 시그마 값 설정
        self.true_sigma = true_sigma
        self.sat_upper = sat_upper  # 상한값 설정
        self.sat_lower = sat_lower  # 하한값 설정
        self.firstRun = True  # 첫 실행 여부 플래그

    def run(self):
        sigma_v = 0.0052  # rad, 공분산(sigma_v) = 0.0052, 3*공분산(sigma_v) = 99% 영역
        v = sigma_v * np.random.randn()  # 평균 0, 표준편차 1 정규분포를 따르는 노이즈 생성

        # 첫 실행인 경우 초기 시그마 값에 노이즈를 더해 측정값 설정
        # 그렇지 않은 경우 실제 시그마 값에 노이즈를 더해 측정값 설정
        if self.firstRun:
            z_sigma = self.sigma + v
            self.firstRun = False
        else:
            z_sigma = self.true_sigma + v

        count = 1  # count는 측정값이 포화 상태인지 여부를 나타냄

        # 측정값이 상한값을 초과하는 경우
        if abs(z_sigma) >= self.sat_upper and z_sigma > 0:
            count = 0  # 포화 상태로 설정
            z_sigma = self.sat_upper  # 측정값을 상한값으로 조정
        # 측정값이 하한값 이하인 경우
        elif abs(z_sigma) >= self.sat_upper and z_sigma < 0:
            count = 0  # 포화 상태로 설정
            z_sigma = self.sat_lower  # 측정값을 하한값으로 조정

        return z_sigma, v, count

