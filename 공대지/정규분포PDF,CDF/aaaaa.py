import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad



# 설정값
true_value = np.deg2rad(7.8)  # 실제 값 (라디안으로 변환)
measurement_limit = np.deg2rad(7.5)  # 관측치의 최대 한계 (라디안으로 변환)
measurement_std = 0.0052  # 관측치의 표준편차 (라디안)

# 가우시안 분포의 x값 범위 설정
x = np.linspace(true_value - 5 * measurement_std, measurement_limit + 5 * measurement_std, 1000)
x_deg = np.rad2deg(x)  # 라디안에서 도로 변환

# 정상적인 가우시안 분포
normal_pdf = norm.pdf(x, loc=true_value, scale=measurement_std)

# 중도 절단 분포
saturated_pdf = np.copy(normal_pdf)
saturated_pdf[x > measurement_limit] = 0  # 한계값 이상은 확률 0으로 설정

# 그래프 그리기
plt.figure(figsize=(10, 6))
plt.plot(x_deg, normal_pdf,color='blue', label='Normal Distribution')
# plt.fill_between(x_deg, 0, normal_pdf, alpha=0.2)
plt.plot(x_deg, saturated_pdf, color='r', label='Saturated Distribution', linestyle='--')
# plt.fill_between(x_deg, 0, saturated_pdf, alpha=0.2)
plt.axvline(x=np.rad2deg(measurement_limit), color='orange', label='Measurement Limit')
plt.title('Gaussian Distribution with Censoring at Measurement Limit')
plt.xlabel('Degrees')
plt.ylabel('Probability Density')
plt.legend()
plt.show()



# 설정값
true_value = np.deg2rad(7.3)  # 실제 값 (라디안으로 변환)
measurement_limit = np.deg2rad(7.5)  # 관측치의 최대 한계 (라디안으로 변환)
measurement_std = 0.0052  # 관측치의 표준편차 (라디안)

# 가우시안 분포의 x값 범위 설정
x = np.linspace(true_value - 5 * measurement_std, measurement_limit + 5 * measurement_std, 1000)
dx = x[1] - x[0]  # x값 간격

# 정상적인 가우시안 분포
normal_pdf = norm.pdf(x, loc=true_value, scale=measurement_std)

# 포화된 분포 생성
saturated_pdf = norm.pdf(x, loc=true_value, scale=measurement_std)
saturated_pdf[x > measurement_limit] = 0  # 한계값 이상은 확률 0으로 설정

# 정상적인 분포의 평균과 표준편차 계산
normal_mean = np.sum(x * normal_pdf) * dx
normal_var = np.sum((x - normal_mean)**2 * normal_pdf) * dx
normal_std = np.sqrt(normal_var)

# 포화된 분포의 평균과 표준편차 계산
saturated_mean = np.sum(x * saturated_pdf) * dx
saturated_var = np.sum((x - saturated_mean)**2 * saturated_pdf) * dx
saturated_std = np.sqrt(saturated_var)

# 결과 출력
print("Normal Distribution Mean (rad):", normal_mean)
print("Normal Distribution Std (rad):", normal_std)
print("Saturated Distribution Mean (rad):", saturated_mean)
print("Saturated Distribution Std (rad):", saturated_std)

# 도 단위로 변환
normal_mean_deg = np.rad2deg(normal_mean)
normal_std_deg = np.rad2deg(normal_std)
saturated_mean_deg = np.rad2deg(saturated_mean)
saturated_std_deg = np.rad2deg(saturated_std)

print("\nNormal Distribution Mean (deg):", normal_mean_deg)
print("Normal Distribution Std (deg):", normal_std_deg)
print("Saturated Distribution Mean (deg):", saturated_mean_deg)
print("Saturated Distribution Std (deg):", saturated_std_deg)





import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# 설정값
true_value = np.deg2rad(7.3)  # 실제 값 (라디안으로 변환)
measurement_limit = np.deg2rad(7.5)  # 관측치의 최대 한계 (라디안으로 변환)
measurement_std = 0.0052  # 관측치의 표준편차 (라디안)

# 가우시안 분포의 x값 범위 설정
x = np.linspace(true_value - 5 * measurement_std, measurement_limit + 5 * measurement_std, 1000)
x_deg = np.rad2deg(x)  # 라디안에서 도로 변환

# 포화된 분포 생성
saturated_pdf = norm.pdf(x, loc=true_value, scale=measurement_std)
saturated_pdf[x > measurement_limit] = 0  # 한계값 이상은 확률 0으로 설정

# 그래프 그리기
plt.figure(figsize=(10, 6))

# 포화된 분포 그리기
plt.plot(x_deg, saturated_pdf, label='Saturated Distribution', linestyle='--')

plt.title('Saturated Gaussian Distribution')
plt.xlabel('Degrees')
plt.ylabel('Probability Density')
plt.legend()
plt.show()








# ############################## PDF, CDF ################################
# # 데이터 생성
# x = np.linspace(-4, 4, 1000)
# pdf = norm.pdf(x, 0, 1)  # 표준정규분포의 확률밀도함수
# cdf = norm.cdf(x, 0, 1)  # 표준정규분포의 누적확률밀도함수


# plt.figure(figsize=(12, 6))
# # 확률밀도함수 (PDF)
# plt.subplot(1, 2, 1)
# plt.plot(x, pdf, label='PDF', color='blue')
# plt.title('Probability Density Function (PDF)')
# plt.xlabel('x')
# plt.ylabel('Probability Density')

# # 누적확률밀도함수 (CDF)
# plt.subplot(1, 2, 2)
# plt.plot(x, cdf, label='CDF', color='red')
# plt.title('Cumulative Distribution Function (CDF)')
# plt.xlabel('x')
# plt.ylabel('Cumulative Probability')
# plt.tight_layout()
# plt.show()



# ############################## sat_lower ################################
# # 정규분포의 평균과 표준편차 설정
# mu = 0
# sigma = 1
# x = np.linspace(-5, 5, 1000)

# # 정규분포 확률밀도함수 계산
# pdf = norm.pdf(x, mu, sigma)

# a, b, c = -np.Infinity, np.deg2rad(-7.5 + 7.8)/0.0052, np.deg2rad(7.5 + 7.8)/0.0052
# print(f'probability pus : { norm.cdf(c)- norm.cdf(b)}')
# print(f'probability pl : { norm.cdf(b)}')
# print(f'probability ph : { 1 - norm.cdf(c)}')
# # 특정 범위에서 확률밀도함수의 적분 계산
# integral, _ = quad(norm.pdf, a, b, args=(mu, sigma))

# # 그래프 그리기
# plt.figure(figsize=(10, 6))
# plt.plot(x, pdf, label='PDF')
# plt.fill_between(x, 0, pdf, where=(x >= a) & (x <= b), color='gray', alpha=0.5, label=f'Integral from -∞ to {b:.2f} = {integral:.4f}')
# plt.title('Probability Density Function')
# plt.xlabel('x')
# plt.ylabel('Density')
# plt.legend()
# plt.grid()
# plt.show()

# ############################## sat_upper ################################
# # 정규분포의 평균과 표준편차 설정
# mu = 0
# sigma = 1
# x = np.linspace(-5, 5, 1000)

# # 정규분포 확률밀도함수 계산
# pdf = norm.pdf(x, mu, sigma)

# # 특정 범위 설정
# a = np.deg2rad(7.5-10)/0.0052
# b = np.Infinity

# # 특정 범위에서 확률밀도함수의 적분 계산
# integral, _ = quad(norm.pdf, a, np.inf, args=(mu, sigma))

# # 그래프 그리기
# plt.figure(figsize=(10, 6))
# plt.plot(x, pdf, label='PDF')

# # 특정 범위에서 색칠
# x_fill = np.linspace(a, 5, 1000)
# y_fill = norm.pdf(x_fill, mu, sigma)
# plt.fill_between(x_fill, y_fill, alpha=0.5, color='gray', label=f'Integral from {a:.2f} to ∞ = {integral:.4f}')

# # 반대쪽 범위에서 색칠
# x_fill_left = np.linspace(-5, a, 1000)
# y_fill_left = norm.pdf(x_fill_left, mu, sigma)
# plt.fill_between(x_fill_left, y_fill_left, alpha=0.3, color='blue', label=f'Integral from -∞ to {a:.2f} = {(1-integral):.4f}')

# # 제목, 라벨, 범례 추가
# plt.title('Probability Density Function')
# plt.xlabel('x')
# plt.ylabel('Density')
# plt.legend()
# plt.grid()
# plt.show()
# print(f'Integral from {a:.2f} to ∞: {integral:.4f}')