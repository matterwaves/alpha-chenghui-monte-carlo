import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

def GenerateEllipse():
    XY = np.array([(np.cos(i), np.cos(i+0.5)) for i in range(20)])
    return XY

def DirectEllipseFit(XY):
    centroid = np.mean(XY, axis=0)
    Z = np.array([XY[:,0]**2, XY[:,0]*XY[:,1], XY[:,1]**2, XY[:,0], XY[:,1]]).T
    D1 = np.array([(XY[:,0]-centroid[0])**2, (XY[:,0]-centroid[0])*(XY[:,1]-centroid[1]),\
      (XY[:,1]-centroid[1])**2]).T
    D2 = np.array([XY[:,0]-centroid[0], XY[:,1]-centroid[1], np.ones(np.size(XY,axis=0))]).T
    S1 = np.dot(D1.T,D1)
    S2 = np.dot(D1.T,D2)
    S3 = np.dot(D2.T,D2)
    T = -np.dot(LA.inv(S3),S2.T)
    M = S1 + np.dot(S2,T)
    M = np.array([M[2,:]/2, -M[1,:], M[0,:]/2])
    val, vec = LA.eig(M)
    veclist = [vec[:,i] for i in range(3) if 4*vec[0,i]*vec[2,i]-vec[1,i]**2 > 0]
    A1 = veclist[0]
    A = np.concatenate((A1,np.dot(T,A1)))
    A3 = A[3]-2*A[0]*centroid[0]-A[1]*centroid[1]
    A4 = A[4]-2*A[2]*centroid[1]-A[1]*centroid[0]
    A5 = A[5]+A[0]*centroid[0]**2+A[2]*centroid[1]**2+\
    A[1]*centroid[0]*centroid[1]-A[3]*centroid[0]-A[4]*centroid[1]
    A[3] = A3
    A[4] = A4
    A[5] = A5
    A = -A/A[5]
    Contrast = np.zeros(2)
    Center = np.zeros(2)

    if A[0] < 0:
        Phase = np.arccos(0.5*A[1]/np.sqrt(A[0]*A[2]))
    else:
        Phase = np.arccos(-0.5*A[1]/np.sqrt(A[0]*A[2]))
    Ratio = np.sqrt(A[0]/A[2])
    Center[0] = (A[3]/Ratio/np.cos(Phase)+A[4])/A[1]/np.tan(Phase)**2
    Center[1] = (A[4]*Ratio/np.cos(Phase)+A[3])/A[1]/np.tan(Phase)**2
    Weight = (A[5]-A[0]*Center[0]**2-A[1]*Center[0]*Center[1]-A[2]*Center[1]**2)/(-np.sin(Phase)**2);
    Contrast[0] = np.sqrt(abs(Weight/A[0]));
    Contrast[1] = np.sqrt(abs(Weight/A[2]));

    #PhaseErr = abs(cot(Phase))*sqrt(Std(2,2)/A(2)^2+Std(1,1)/4/A(1)^2+Std(3,3)/A(3)^2-...
    #Std(1,2)/A(2)/A(1)-Std(2,3)/A(2)/A(3)+Std(1,3)/A(1)/A(3)/2);
    #sigma = sqrt(sum((Z*A(1:5)-1).^2)/(size(XY,1)-5));
    #MM = Z'*Z;
    #Std = sigma*MM.^(-0.5);
    
    return A, Phase, Contrast, Center

def DisplayEllipse(XY, Phase, Contrast, Center):
    Theta = np.arange(1000)/1000*2*np.pi;
    FitX = Contrast[0]*np.cos(Theta)+Center[0];
    FitY = Contrast[1]*np.cos(Theta+Phase)+Center[1];
    plt.plot(XY[:,0], XY[:,1], 'o', FitX, FitY);
