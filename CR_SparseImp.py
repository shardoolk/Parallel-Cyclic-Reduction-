


import scipy
import scipy.sparse
import scipy.linalg
import matplotlib.pyplot as plt
import numpy as np
import time

N = int(input("Enter the value of N:  "))

p = []

for k in range(25):
    p.append(pow(2,k)-1)

p = np.array(p)

if N in p:
    A = scipy.sparse.diags([-1, 2, -1], [-1, 0, 1], shape=(N,N)).toarray()
    b = np.ones((N,1)).reshape(N,1)
    x = np.zeros((N,1))
else:
    N1 = p[p>=N][0]
    A = np.zeros((N1,N1))
    b = np.zeros((N1,1))
    A[:N,:N] = scipy.sparse.diags([-1, 2, -1], [-1, 0, 1], shape=(N,N)).toarray()
    b = np.ones((N1,1)).reshape(N1,1)
    x = np.zeros((N1,1))
    for i in range(N,N1):
        A[i][i] = 1    
        
subD = np.diag(A, k = -1).copy()
supD = np.diag(A, k = 1).copy()
mainD = np.diag(A, k = 0).copy()   
newArr = np.zeros(len(mainD))
subD = np.insert(subD,0,0)
supD = np.append(supD,0)

xexact = scipy.linalg.solve(A[:N,:N],b[:N])
N = len(A)

start = time.time()


for i in range(0,int(np.log2(N+1))):
    for j in range(pow(2,i+1)-1,N,pow(2,i+1)):
        offset = pow(2,i)
        index1 = j - offset
        index2 = j + offset
        
        alpha = subD[j]/mainD[index1]
        gamma = supD[j]/mainD[index2]
        
        subD[j] = -subD[index1]*(alpha)
        mainD[j] = mainD[j] - supD[index1]*alpha - subD[index2]*gamma
        supD[j] = -supD[index2]*(gamma)
        b[j] = b[j] - b[index1] * alpha - b[index2] * gamma

index = int((N-1)/2)
x[index] = b[index]/mainD[index]

for i in range(int(np.log2(N+1)),-1,-1):
    for j in range(pow(2,i+1)-1,N,pow(2,i+1)):
        offset = pow(2,i)
        index1 = j - offset
        index2 = j + offset
        
        if (j != index1):
            if (index1 - offset < 0):
                #print("index1 - offset = ",index1 - offset)
                x[index1] = (b[index1]- supD[index1]*x[index1+offset])/mainD[index1]
                
            else:
                x[index1] = (b[index1] - subD[index1]*x[index1-offset] - supD[index1]*x[index1+offset])/mainD[index1]
                
        if(j != index2):
            if(index2 + offset >= N ):
                x[index2] = (b[index2] - subD[index2]*x[index2-offset])/mainD[index2]
                
            else:
               x[index2] = (b[index2] - subD[index2]*x[index2-offset] - supD[index2]*x[index2+offset])/mainD[index2]
               
        




             
end = time.time()

error = xexact - x[:int(len(xexact))]

mse = (np.square(error)).mean(axis=0)
print(mse[0])
#plt.plot(xexact,'*')
#plt.plot(x[:int(len(xexact))])
#plt.grid()
#plt.legend(['Exact', 'Cyclic Reduction Solution'])    

#print("time = ", (end-start))
        