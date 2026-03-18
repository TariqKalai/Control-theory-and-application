import math
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error, median_absolute_error, explained_variance_score

def Broida_complicated(t1, t2):

    T = 5.5*(t2 - t1)
    theta = 2.8*t1 - 1.8*t2

    return (T, theta)


def Broida_simple(Tg, Tu):

    T = Tg
    theta = Tu
    return (T, theta)




def Van_der_Grinten(a, Tu, Tg):
    e = math.exp(1)

    T1 = Tg*((3*a*e - 1)/(1 + a*e))
    T2 = Tg*((1 - a*e)/(1 + a*e))
    theta = Tu - ((T1*T2)/(T1 + 3*T2))

    return (T1, T2, theta)


def Strejc(Tu, Tg):

    ratio = Tu/Tg

    match (ratio) :
        case n if 0 < n <= 0.1 :
            a_n = 0
            order = 1
            b_n = 1

        case n if 0.1 < n <= 0.22 :
            a_n = 0.1
            order =2
            b_n = 2.72

        case n if 0.22 < n <= 0.32 :
            a_n = 0.22
            order =3
            b_n = 3.69
        case _ :
            raise ValueError("error Tu/Tg ratio is too high")
    
    T = Tg / b_n
    Tuth = a_n * Tg
    theta = Tu -Tuth

    return (T,theta, order)




def Error_quant(y_true, y_pred):
    return {
        "R2": r2_score(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": mean_squared_error(y_true, y_pred),
        "MedianAE": median_absolute_error(y_true, y_pred)
    }
