from scipy.fftpack import dct, idct
from scipy.sparse import rand, eye
import numpy as np
from numpy.random import randint, random
import pylab as pl

N = 10
ERROR_MIN = 1e-2
PERC_MAX = 100
PERC_STEP = 2

def f(A,B):
    #relative error in frobenius norm
    s = 0.0
    for i in range(N):
        for j in range(N):
            s += (A[i][j] - B[i][j])**2
    return np.sqrt(s)/N
    
def imshowM(A):
    pl.imshow(A, interpolation="nearest")
    pl.colorbar()
    pl.show()

def constructSparseMatrix(i):
    A = np.zeros([N,N], float)
    numNonZeros = i*N*N/100
    for i in range(numNonZeros):
        i = randint(N)
        j = randint(N)
        while True:
            if (A[i][j] == 0):
                A[i][j] = 2.0*(random()-.5)
                break
            else:
                i = randint(N)
                j = randint(N)
    return A
         
    
def main():
    count = 0
    
    percSampled = []
    perc = []
    # for i = percentage of non-zero elements in matrix A
    percRange = np.arange(1,PERC_MAX,PERC_STEP)
    percSampled = np.zeros(len(percRange))
    pl.ion()
    for p in percRange:
        print "working on " + str(p) + " out of " + str(PERC_MAX)
        
        Sampled = np.zeros([N,N], int) #record of which indeces have been sampled
        A = constructSparseMatrix(p)   #true matrix in sparse basis
        B = dct(A, norm='ortho')       #true matrix in dense basis
        Ap = np.zeros([N,N], float)    #sampled matrix in sparse basis
        Bp = np.zeros([N,N], float)    #sampled matrix in dense basis
        error = 1.0                    #starting error    
        numSteps = 0                   #count of sampling iterations
        #imshowM(A)
        #imshowM(B)
        ax3 = pl.subplot(2, 2, 3)
        ax3.clear()
        pl.xlabel("percentage non-zero elements")
        pl.ylabel("num iter for err. " + str(ERROR_MIN))

        ax3.plot(percRange,percSampled)
        #txt3 = ax3.text(1,1, 'iteration: ' + str(numSteps))
        pl.draw()
        
        while (error > ERROR_MIN):
            i = randint(N)                 
            j = randint(N)
            while True:
                #make sure this spot has not been sampled (can remove as condition later to test in practice)
                if (Sampled[i][j] == 1):
                    i = randint(N)
                    j = randint(N)
                else:
                    Bp[i][j] = B[i][j]            #sample the true dense matrix
                    Ap = idct(Bp, norm='ortho')    #get the sparse matrix guess from the sampled dense matrix
                    
                    if numSteps%100 == 0:
                        
                        ax2 = pl.subplot(2, 2, 2)
                        ax2.imshow(Ap, interpolation="nearest")
                        txt2 = ax2.text(1,1, 'iteration: ' + str(numSteps))
                        pl.draw()
                        ax1 = pl.subplot(2, 2, 1)
                        txt1 = ax1.text(1,1, 'perc. non-0: ' + str(p))
                        ax1.imshow(A, interpolation="nearest")
                        pl.draw()
                        txt2.remove()
                        txt1.remove()
                        
                    error = f(Ap,A)               #error of sparse matrix calculated
                    numSteps += 1                 #update the number of counts
                    break
        percSampled[count] = numSteps
        perc.append(p)
        count += 1
                
                
                
                



if __name__ == "__main__":
    main()



