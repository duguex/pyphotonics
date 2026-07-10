import numpy as np
a=np.array([1,2,3]).reshape(3,1)
b=np.array([4,5,6]).reshape(1,3)
print(np.matmul(a,b))