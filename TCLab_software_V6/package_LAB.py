import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from package_DBR import Process, Bode

  

def LL_RT(MV,Kp,Tlead, Tlag, Ts , PV , PVInit=0,method='EBD'):
    
    """
    The function "FO_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "FO_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    
    K = Ts/Tlag
    if len(PV) == 0:
        PV.append(PVInit)
    else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
        if method == 'EBD':
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+ (Tlead/Ts) )*MV[-1] )- ((Tlead/Ts)*MV[-2])) )


        elif method == 'EFD':
            PV.append((1-K)*PV[-1] + K*Kp* ((Tlead/Ts)*MV[-1] + (1- (Tlead/Ts))*MV[-2]))


        # elif method == 'TRAP':
        #     PV.append((1/(2*Tlag+Ts))*((2*Tlag-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))      


        else:
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+ (Tlead/Ts) )*MV[-1] )- ((Tlead/Ts)*MV[-2])) )
   



def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E , ManFF = False, PVinit=0, method="EBD"):

    # initialisation de E 

    if len(PV) == 0 :
        E.append(SP[-1] - PVinit)
    else :
        E.append(SP[-1] - PV[-1])

    # Partie MV process 

    MVP.append(Kc*E[-1])

    # Partie MV intergal
    if len(MVI) == 0:
        MVI.append(Kc*(Ts/Ti)*E[-1])
    else:
        MVI.append(MVI[-1] + Kc*(Ts/Ti)*E[-1])

    

    # Partie MV derivée
    TFD = Td * alpha
    if len(MVD) != 0:
        if len(E) == 1 :
            MVD.append( (((TFD)/(TFD+Ts)) * MVD [-1] ) + ((Kc*Td)/(TFD+Ts))*(E[-1])) 
        else :
            MVD.append( (((TFD)/(TFD+Ts)) * MVD [-1] ) + ((Kc*Td)/(TFD+Ts))*(E[-1] - E[-2]))
    else :
        if len(E) == 1 :
            MVD.append(  ((Kc*Td)/(TFD+Ts))*(E[-1])) 
        else :
            MVD.append(  ((Kc*Td)/(TFD+Ts))*(E[-1] - E[-2]))

## pas normal j comprend pas ce qu il se passe MVFF ne change pas grand chose, mais MVMan n est pas respecté? 
    if Man[-1] == True :
        if ManFF :
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else :
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]

      


    MV.append( MVD[-1] +MVI[-1] + MVP[-1] + MVFF[-1])

    if MV[-1] > MVMax :

        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF [-1]
        MV[-1] =  MVD[-1] +MVI[-1] + MVP[-1] + MVFF[-1]


    if MV[-1] < MVMin :

        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF [-1]
        MV[-1] =  MVD[-1] +MVI[-1] + MVP[-1] + MVFF[-1]
    


def IMC_Tuning(Kp, gamma, theta, T1, T2=0, order=2 ):

    if T1 > 10*T2:
        Tc = T1 *gamma
    elif T2 > 10*T1:
        Tc = T2*gamma
    else :
        Tc = (T1+T2)*gamma
    
    if order == 1 :
        Kc  = (T1 + (theta/2)) / ((Tc+ (theta/2)) * Kp )
        Ti = T1 + theta/2
        Td = (T1*theta)/ (2*T1 +theta)
        
        return (Kc, Ti, Td)
    
    else :

        Kc = (T1 + T2)/((Tc + theta)*Kp)
        Ti = T1 + T2
        Td = (T1*T2)/(T1 + T2)

        return (Kc, Ti, Td)
    

def Margin_error(Process : Process):
    omega = np.logspace(-4, 1, 10000)
    bode_result = Bode(Process, omega, Show=False)

    gain = np.abs(bode_result)
    phase = np.degrees(np.unwrap(np.angle(bode_result)))

    # Gain margin
    if np.min(phase) > -180:
        Gain_margin = np.inf
    else:
        index_GainMargin = np.argmin(np.abs(phase + 180))
        
        # Formule : GM = 1/Gain au croisement, converti en dB
        Gain_margin = 20 * np.log10(1 / gain[index_GainMargin])
    
    # Phase margin 
    if np.max(gain) < 1:
        Phase_margin = np.inf
    else:
        index_PhaseMargin = np.argmin(np.abs(gain - 1))
        # Formule : Distance entre la phase et -180°
        Phase_margin = phase[index_PhaseMargin] + 180

    return Gain_margin, Phase_margin