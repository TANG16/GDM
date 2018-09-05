# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 12:12:28 2018

@author: Erdem Varol
"""
import numpy as np
from scipy import linalg

def gdr_primal(X,Y,covar,lam):
    dim = X.shape
    C = np.concatenate((np.ones((dim[0],1)),covar),axis=1)
    CCC = C @ linalg.inv(C.T @ C)
    CCCC = CCC @ C.T
    J = linalg.inv((X.T @ X) - X.T @ CCCC @ X + lam*(Y.T @ Y) * np.eye(dim[1]) - lam * (Y.T @ CCCC @ Y) * np.eye(dim[1])) @ (2*lam * (X.T @ Y - X.T @ CCCC @ Y))
    A0 = X.T @ CCC - J @ Y.T @ CCC
    W0 = CCC.T @ Y - CCC.T @ X @ J
    Q = linalg.inv((X.T @ X) - X.T @ CCCC @ X + lam*(Y.T @ Y) * np.eye(dim[1]) - lam * (Y.T @ CCCC @ Y) * np.eye(dim[1])) @ (2*lam * (X.T - X.T @ CCCC))
    return {'J':J,'A0':A0,'W0':W0,'Q':Q}
    
def gdr_dual(X,Y,covar,lam):
    dim = X.shape
    dimk = covar.shape
    C = np.concatenate((np.ones((dim[0],1)),covar),axis=1)
    CCC = C @ linalg.inv(C.T @ C)
    CCCC = CCC @ C.T
    K = X @ X.T
    M2 = np.concatenate((np.concatenate(((-K / (lam*(Y.T @ Y) - lam*(Y.T @ CCCC @ Y)) - np.eye(dim[0])),C),axis=1),np.concatenate((C.T, np.zeros((dimk[1]+1,dimk[1]+1))),axis=1)),axis=0)
    b2 = np.concatenate((np.eye(dim[0]) + (lam*K @ (CCCC - np.eye(dim[0]))) / (lam*Y.T@Y - lam*Y.T @ CCCC @ Y), np.zeros((dimk[1]+1,dim[0]))),axis=0)
    iM = linalg.inv(M2)
    
    v = iM @ b2 @ Y
    L = v[0:dim[0]]
    J = (lam*X.T@Y - lam*X.T @ CCCC @ Y - X.T @ L) @ linalg.inv(lam*(Y.T@Y) - lam*Y.T @ CCCC @ Y)
    A0 = X.T @ CCC - J @ Y.T @ CCC
    W0 = CCC.T @ Y - CCC.T @ X @ J
    Q = (lam*X.T - lam*X.T @ CCCC - X.T @ iM[0:dim[0],0:dim[0]] @ b2[0:dim[0],0:dim[0]]) / (lam*(Y.T @ Y) - lam*(Y.T @ CCCC @ Y))
    return {'J':J,'A0':A0,'W0':W0,'Q':Q}