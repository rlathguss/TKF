import numpy as np
gt=np.array([1/10,0])
ww=np.array([[10,10],[10,10]])
g=gt.T
print(gt.shape)
print(ww.shape)
print(g.shape)

print(ww@gt)
print(g@ww@gt)


# print(np.array([[1/10],[0]])@np.array([[10,10],[10,10]])@np.array([1/10,0]))