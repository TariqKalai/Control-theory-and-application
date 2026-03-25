import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from package_DBR import Process, Bode

  

def LL_RT(MV,Kp,Tlead, Tlag, Ts , PV , PVInit=0,method='EBD'):
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "LL_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that takes into consideration lead and lag time constants.
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

    """
    The function "PID_RT" needs to be included in a "for or while loop".

    :SP: Set Point vector
    :PV: Process Value/output vector
    :Man: Manual controller mode vector (True or False)
    :MVMan: Manual Value for MV vector
    :MVFF: Feedforward vector
    :Kc: Controller gain
    :Ti: Integral time constant [s]
    :Td: Derivative time constant [s]
    :alpha: TFD = Td * alpha, with TFD being the derivative filter time constant [s]
    :Ts: Sampling period [s]
    :MVMin: Minimum Value of MV (used for saturation and anti wind-up)
    :MVMax: Maximum Value of MV (used for saturation and anti wind-up)
    :MV: Manipulated Value/input vector
    :MVP: Proportional part of MV vector
    :MVI: Integral part of MV vector
    :MVD: Derivative part of MV vector
    :E: Control Error vector
    :ManFF: FF in Manual Mode (optional: default boolean value is False)
    :PVinit: Initial value for PV (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD".
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    The saturation of "MV" is implemented with anti wind-up. 
    """


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
        if ManFF[-1] :
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

    """
    Computes the optimised PID controller settings for FOPDT and SOPDT processes.

    :Kp: Process gain
    :gamma: Loop response time as a ratio of T1    
    :theta: Process delay.
    :T1: First order lag time constant.
    :T2: (SOPDT only) Second order lag time constant.
    :order:
        1 : First Order Plus Dead Time for PID control (IMC tuning case H)
        2 : Second Order Plus Dead Time for PID control (IMC tuning case I)
        
    return : PID controller parameters Kc, Ti and Td in a tuple (Kc, Ti, Td)
    """

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
    

# def Margin_error(Process : Process):

#     """
#     Computes the gain and phase margin and analyses the robustness of the PID controller.

#     :Process: Process as defined by the class "Process".
#         Use the following command to define the default process which is simply a unit gain process:
#             P = Process({})
#     """

#     omega = np.logspace(-4, 1, 10000)
#     bode_result = Bode(Process, omega, Show=False)

#     gain = np.abs(bode_result)
#     phase = np.degrees(np.unwrap(np.angle(bode_result)))

#     # Gain margin
#     if np.min(phase) > -180:
#         Gain_margin = np.inf
#     else:
#         index_GainMargin = np.argmin(np.abs(phase + 180))
        
#         # Formule : GM = 1/Gain au croisement, converti en dB
#         Gain_margin = 20 * np.log10(1 / gain[index_GainMargin])
    
#     # Phase margin 
#     if np.max(gain) < 1:
#         Phase_margin = np.inf
#     else:
#         index_PhaseMargin = np.argmin(np.abs(gain - 1))
#         # Formule : Distance entre la phase et -180°
#         Phase_margin = phase[index_PhaseMargin] + 180

#     return Gain_margin, Phase_margin

def Margin_error(Process: Process):

    """
    Computes and displays on a graph the gain and phase margin and analyses the robustness of the PID controller.

    :Process: Process as defined by the class "Process".
        Use the following command to define the default process which is simply a unit gain process:
            P = Process({})

    :returns: gain margin and phase margin in a tuple (gm_val, pm_val)
    """
    omega = np.logspace(-4, 1, 10000)
    bode_result = Bode(Process, omega, Show=False)

    gain_linear = np.abs(bode_result)
    gain_db = 20 * np.log10(gain_linear)
    phase = np.degrees(np.unwrap(np.angle(bode_result)))

    idx_180 = np.where(np.diff(np.sign(phase + 180)))[0]
    
    idx_0dB = np.where(np.diff(np.sign(gain_linear - 1)))[0]

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(18, 12), sharex=True)

    ax_mag.semilogx(omega, gain_db, color='blue')
    ax_mag.axhline(0, color='black', linestyle='--', lw=1)
    
    if len(idx_180) > 0:
        w_180 = omega[idx_180[0]]
        gm_val = -gain_db[idx_180[0]] 
        ax_mag.annotate('', xy=(w_180, 0), xytext=(w_180, gain_db[idx_180[0]]),
                        arrowprops=dict(arrowstyle='<->', color='green', lw=2))
        ax_mag.text(w_180, gain_db[idx_180[0]]/2, f' gm: {gm_val:.1f} dB', color='green', fontweight='bold')
    
    ax_phase.semilogx(omega, phase, color='orange')
    ax_phase.axhline(-180, color='black', linestyle='--', lw=1)
    
    if len(idx_0dB) > 0:
        w_0dB = omega[idx_0dB[0]]
        pm_val = phase[idx_0dB[0]] + 180
        ax_phase.annotate('', xy=(w_0dB, -180), xytext=(w_0dB, phase[idx_0dB[0]]),
                          arrowprops=dict(arrowstyle='<->', color='red', lw=2))
        ax_phase.text(w_0dB, (phase[idx_0dB[0]] - 180)/2, f' pm: {pm_val:.1f}°', color='red', fontweight='bold')
    else:
        pm_val = np.inf
        # ax_phase.text(0.1, -10, "Marge de Phase Infinie (Gain < 0dB)", color='red', transform=ax_phase.transAxes)

    ax_phase.set_ylim([-270, 0]) 
    ax_mag.grid(True, which="both")
    ax_phase.grid(True, which="both")
    plt.show()

    return gm_val, pm_val