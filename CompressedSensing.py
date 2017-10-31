from scipy.fftpack import dct, idct
from scipy.sparse import rand, eye
import numpy as np
from numpy.random import randint, random
import pylab as pl


N = 20
ERROR_MIN = 1e-1
PERC_MAX = 100
PERC_STEP = 3
NUM_ITERS = 5
NUM_FRAMES = 10

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
    
def calcAPriori(percRange):
    a = []
    for p in percRange:
        a.append(p*N*N/100.0)
    return a 
   
def calcNaive(percRange):
    a = []
    for p in percRange:
        a.append(N*N)
    return a

def newAverage(average, percSampled, n):
    #calc running average given new data
    if (n != 0):
        for i in range(len(average)):
            s = average[i] * (n - 1)
            s += percSampled[i]
        
            average[i] = s/n
    return average
        
    

def main():
    percRange = np.arange(1,PERC_MAX,PERC_STEP)
    aPriori = calcAPriori(percRange)
    naive = calcNaive(percRange)
    percSampled = np.zeros(len(percRange))    
    average = np.zeros(len(percRange), float)
    for n in range(NUM_ITERS):
        average = newAverage(average, percSampled, n)
        count = 0
    
        percSampled = []
        perc = []
        # for i = percentage of non-zero elements in matrix A
        aPriori = calcAPriori(percRange)
        naive = calcNaive(percRange)
        percSampled = np.zeros(len(percRange))
        pl.ion()
        
        ax4 = pl.subplot(2, 2, 4)
        ax4.clear()
        pl.xlabel("percentage non-zero elements")
        pl.ylabel("average iters ")
        pl.ylim([0, N**2 + 5])
        ax4.plot(percRange, average)
        ax4.plot(percRange, aPriori, "g-")
        ax4.plot(percRange, naive, "r.")
        ax4.legend(['compressed','naive','perfect info'], loc=4, prop={'size':8})
        
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
            pl.ylim([0, N**2 + 5])
            ax3.plot(percRange,percSampled)
            ax3.plot(percRange, aPriori, "g-")
            ax3.plot(percRange, naive, "r.")
            ax3.legend(['compressed','naive','perfect info'], loc=4, prop={'size':8})
        
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
                        Sampled[i][j] = 1
                        Ap = idct(Bp, norm='ortho')    #get the sparse matrix guess from the sampled dense matrix
                    
                        if numSteps%(N*N/NUM_FRAMES) == 0:
                        
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
                        break
                numSteps += 1                 #update the number of counts
            percSampled[count] = numSteps
            perc.append(p)
            count += 1
                
                
                



if __name__ == "__main__":
    main()



