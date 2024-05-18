import numpy as np

# 임의의 F_k 행렬과 타우(τ) 정의
F_k = np.array([[1, 2], [3, 4]])  # 예시 F_k
tau = 0.1  # 예시 타우 값

# 단위 행렬 I 정의
I = np.eye(F_k.shape[0])

# (I + F_k*τ) 계산
expr1 = I + F_k*tau

# (I + F_k*τ)의 전치 계산
expr1_T = expr1.T

# 두 표현식 곱하기
result = np.dot(expr1, expr1_T)

print(result)