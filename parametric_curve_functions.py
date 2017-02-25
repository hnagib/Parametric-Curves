import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import *
import seaborn as sns

def bezier_poly_interact(ctrl_points, pt, t=0.5, P1x=-1, P1y=2, P2x=3, P2y=2):
    
    if type(ctrl_points) == int:
        P_i = np.array([[-2,-1,3,4,1,0,3], [0,2,2,0,1,0.5,-0.5]])
        
        # Control points
        P_i[0][1] = P1x
        P_i[1][1] = P1y
        P_i[0][2] = P2x
        P_i[1][2] = P2y
        
    else:
        P_i = ctrl_points
    
    # Parameter space
    h = 0.01
    u = np.arange(0,1,h)
    t = int(t/h)
    

    
    # Select N control points from P_i
    P_i = P_i[:,0:pt]
    
    # Initialize emp
    P_u = np.zeros((2,u.size))
    Bi = []
    
    n = P_i.shape[1]-1; 
    for i in range(n+1):
        Bi.append(comb(n,i)*(u**i)*((1-u)**(n-i))) 
        P_u += Bi[i]*P_i[:,i].reshape(2,1)
           
    plt.subplot(1,2,1)
    ax = plt.gca()
    plt.plot(P_i[0,:],P_i[1,:],'k--')
    for i in range(n+1): 
        plt.plot(P_i[0,i],P_i[1,i],'o')
    plt.plot(P_u[0,:],P_u[1,:],'k'), plt.plot(P_u[0,t],P_u[1,t],'ko')
    ax.set_xlim(-3, 5)
    ax.set_ylim(-1, 3)
    plt.title('Bezier Curves')
    
    plt.subplot(1,2,2)
    for i in range(n+1): 
        plt.plot(u,Bi[i]), plt.plot(u[t],Bi[i][t],'ko')
    plt.title('Blending functions')
    
def non_preiodic_knot(n,k):
    u_i = np.zeros(n+k+1)
    for i in range(n+k+1):
        if 0<=i<k:
            u_i[i] = 0
        elif k<=i<=n:
            u_i[i] = i-k+1
        elif n<i<=(n+k):
            u_i[i] = n-k+2
    return u_i

def blending_func(u,u_i,k):
    N =[]
    for k in range(k+1):

        N.append([])

        if k == 0:
            for i in range(len(u_i)-k):
                N[k].append(u_i[i])

        elif k == 1:
            for i in range(len(u_i)-k):
                N[k].append(((u >= u_i[i]) & (u < u_i[i+1]) ).astype(int))

        else:
            for i in range(len(u_i)-k):

                if (u_i[i+k-1]-u_i[i]) == 0:
                    N_ik_1 = np.zeros(u.size)
                else:
                    N_ik_1 = ((u-u_i[i])*N[k-1][i])/(u_i[i+k-1]-u_i[i])

                if (u_i[i+k]-u_i[i+1]) == 0:
                    N_ik_2 = np.zeros(u.size)
                else:
                    N_ik_2 = ((u_i[i+k]-u)*N[k-1][i+1])/(u_i[i+k]-u_i[i+1])

                N_ik = N_ik_1 + N_ik_2

                N[k].append(N_ik)
                
    return N


def B_spline(N,u,u_i,P_i,k):  
    P_u = np.zeros((2,u.size))
    for i in range(len(u_i)-k):
        P_u += N[k][i]*P_i[:,i].reshape(2,1)
    
    return P_u

def plot_bsline_blendfunc(t, N, P_i, P_u, u, u_i,k):   
    plt.subplot(3,2,1)
    ax = plt.gca()
    plt.plot(P_u[0],P_u[1],'k')
    plt.plot(P_u[0][t],P_u[1][t],'ko')
    
    plt.plot(P_i[0,:],P_i[1,:],'k--')
    for i in range(P_i.shape[1]): 
        plt.plot(P_i[0,i],P_i[1,i],'o')  
    ax.set_xlim(-3, 5)
    ax.set_ylim(-1, 3)
    plt.title('B-Spline Curve')
 
    for k in np.arange(2,k+1):
        plt.subplot(3,2,k)
        plt.axvline(u[t],color='k', linestyle='dotted')
        for i in range(len(u_i)-k):
            plt.plot(u,N[k][i])
            plt.plot(u[t],N[k][i][t],'ko')
            plt.title('Blending functions of degree %d' %(k-1))
            
def B_spline_interact(ctrl_points, pt=8, k=4, t=0.5, P1x=-1,P1y=2):

    if type(ctrl_points) == int:
        P_i = np.array([[-2,-1,3,4,1,0,2,0], [0,2,2,0,1,0.5,-0.5,-0.5]])

        # Control points
        P_i[0][1] = P1x
        P_i[1][1] = P1y
     
    else:
        P_i = ctrl_points

    t = float(t); 
    if k == '':
        k = 2
    else:
        k = int(k)
        
    P_i = P_i[:,0:pt]
    n = P_i.shape[1]-1; k = k; h = 0.01
    
    if k > (n+1):
        print("Can't draw degreee = %d curve through %d control points"%(k,n+1))
    elif k == 0:
        print("k must be greater than 0")
    else:
        u_i = non_preiodic_knot(n,k)        
        u   = np.arange(u_i[k-1],u_i[n+1],h)
        t = int((t/h)*u_i[-1])
        N   = blending_func(u,u_i,k)
        P_u = B_spline(N,u,u_i,P_i,k)

        plot_bsline_blendfunc(t,N, P_i, P_u, u, u_i,k);