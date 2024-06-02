import numpy as np


def measm_tkf(sigma_0, true_sigma, sat_upper, sat_lower, firstRun):
    sigma_v = 0.0052  # rad, 공분산(sigma_v) = 0.0052, 3*공분산(sigma_v) = 99% 영역
    v = sigma_v * np.random.randn()  # 평균 0, 표준편차 1 정규분포를 따르는 노이즈 생성

    if firstRun == 0:
        z_sigma = sigma_0 + v
    else:
        z_sigma = true_sigma + v

    # 측정값이 상한값을 초과하는 경우
    if abs(z_sigma) >= sat_upper and z_sigma > 0:

        z_sigma = sat_upper  # 측정값을 상한값으로 조정
    # 측정값이 하한값 이하인 경우
    elif abs(z_sigma) >= sat_upper and z_sigma < 0:

        z_sigma = sat_lower  # 측정값을 하한값으로 조정

    return z_sigma, v