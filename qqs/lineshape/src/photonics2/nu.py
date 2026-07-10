import numpy as np

def direct_mul(a,b):
    la=len(a)
    lb=len(b)
    return np.matmul(a.reshape(la,1),b.reshape(1,lb))


a=np.array([1,2,3])
b=np.array([1,10,100])#.reshape(3,1)
print(b)
print(np.kron(a,b))