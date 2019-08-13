import numpy as np
from sklearn.metrics import *

def slope(pred,true):
    return ( np.sum(true*pred) / np.sum(pred*pred) )

def Rsquared0(pred,true):
    true_mean = np.mean(true)
    yr0 = pred * slope(pred,true)
    first_term = (true - yr0)**2 #(true - yr0)
    second_term = (true-true_mean)**2 #(true-true_mean)
    return (1 - ( np.sum(first_term) / np.sum(second_term) ) )

#squared correlation coefficient
def Rsquared(pred,true):
    true_mean = np.mean(true)
    pred_mean = np.mean(pred)
    first_term = np.sum( (true-true_mean) * (pred - pred_mean) )
    second_term = np.sqrt( np.sum( (true-true_mean)**2) * np.sum( (pred-pred_mean)**2) )
    div = first_term / second_term
    return div*2

#crossvalidated coefficient of determination
def Qsquared(pred,true):
    true_mean = np.mean(true)
    first_term = (np.abs(pred-true))**2                                                                                            
    second_term = (np.abs(true-true_mean))**2
    ret = 1 - (np.sum(first_term) / np.sum(second_term))
    return ret

def RMSE(pred,true):
    rmse = np.sqrt(mean_squared_error(true,pred))
    return rmse

class ValidationMetrics:
      def __init__(self,pred,true):
          self.RMSE = RMSE(pred,true)
          self.Qsquared = Qsquared(pred,true)
          self.Rsquared = Rsquared(pred,true)
          self.Rsquared0 = Rsquared0(pred,true)
          self.slope = slope(pred,true)

