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
    if len(PV) == 0 :
        PV.append(PVInit)
    else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
        if method == 'EBD':
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+ (Tlead/Ts) )*MV[-1] )- ((Tlead/Ts)*MV[-2])) )


        elif method == 'EFD':
            PV.append((1-K)*PV[-1] + K*Kp* ((Tlead/Ts)*MV[-1] + (1- (Tlead/Ts))*MV[-2]))


        elif method == 'TRAP':
            PV.append(PV[-1]*((2-K)/(2+K)) + ((Kp*K)/(2+K)) * (((2*Tlead/Ts)+1)*MV[-1] + (1-(2*Tlead/Ts) )*MV[-2]))

        else: #default : EBD
            PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+ (Tlead/Ts) )*MV[-1] )- ((Tlead/Ts)*MV[-2])) )
   

def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E , ManFF = False, PVinit=0, method="EBD-EBD"):

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
    - EBD-EBD: EBD for integral action and EBD for derivative action 
    - EBD-TRAP: EBD for integral action and TRAP for derivative action 
    - TRAP-EBD: TRAP for integral action and EBD for derivative action 
    - TRAP-TRAP: TRAP for integral action and TRAP for derivative action

    The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD".
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    The saturation of "MV" is implemented with anti wind-up. 
    """

    method_I, method_D= method.split('-')
    # initialisation de E 

    if len(PV) == 0 :
        E.append(SP[-1] - PVinit)
    else :
        E.append(SP[-1] - PV[-1])

    # Partie MV process 

    MVP.append(Kc*E[-1])

    # Partie MV intergal
    if method_I == 'TRAP' : 
        if len(MVI) == 0 :
            MVI.append(0.5*Kc*Ts/Ti * (E[-1]))
        else : 
            MVI.append(MVI[-1] + 0.5*Kc*Ts/Ti * (E[-1] + E[-2]))
    else :
        if len(MVI) == 0 :
            MVI.append((Kc*(Ts/Ti))*E[-1])
        else : 
            MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])


    # Partie MV derivée
    TFD = Td * alpha
    if len(MVD) != 0:
        if method_D =="EBD" : 
            if len(E) == 1 :
            #elif car au debut on peut pas aller chercher des erreurs qui n'existe pas
                MVD.append( (((TFD)/(TFD+Ts)) * MVD [-1] ) + ((Kc*Td)/(TFD+Ts))*(E[-1])) 
            else :
                MVD.append( (((TFD)/(TFD+Ts)) * MVD [-1] ) + ((Kc*Td)/(TFD+Ts))*(E[-1] - E[-2]))
        else :
            if len(E) == 1 : 
                MVD.append(((TFD-Ts/2)/(TFD+Ts/2))*MVD[-1] + ((Kc*Td)/(TFD + Ts/2))*(E[-1]))
            else :
                MVD.append(((TFD-Ts/2)/(TFD+Ts/2))*MVD[-1] + ((Kc*Td)/(TFD + Ts/2))*(E[-1]- E[-2]))

    else :
        if method_D =="EBD" : 
            if len(E) == 1 :
                MVD.append( ((Kc*Td)/(TFD+Ts))*(E[-1])) 
            else :
                MVD.append( ((Kc*Td)/(TFD+Ts))*(E[-1] - E[-2]))


        else :
            if len(E) == 1 :
                MVD.append(((Kc*Td)/(TFD + Ts/2))*(E[-1]))
            else : 
                MVD.append(((Kc*Td)/(TFD + Ts/2))*(E[-1]- E[-2]))

    # anti Wind up

    # mode manuel
    if Man[-1] == True :
        if ManFF :
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else :
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]

      


    MV.append( MVD[-1] +MVI[-1] + MVP[-1] + MVFF[-1])

    #Saturation
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

    if T1 > 5*T2:
        Tc = T1 *gamma
    elif T2 > 5*T1:
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
    
class PID:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kp'] if 'Kp' in parameters else 0
        self.parameters['Ti'] = parameters['theta'] if 'theta' in parameters else 1
        self.parameters['Td'] = parameters['Tlead1'] if 'Tlead1' in parameters else 0
        self.parameters['alpha'] = parameters['Tlead2'] if 'Tlead2' in parameters else 0


def Margin_error(P: Process, C: PID, omega, title =""):
    """
    The function "Margin_error" computes and plots the Bode diagram of the
    open-loop transfer function L(s) = C(s) * P(s) and returns the gain
    and phase margins.

    :P: Process object containing the process transfer function parameters imported from package DBR
    :C: PID object containing the controller parameters:
        - Kc   : controller gain
        - Ti   : integral time constant [s]
        - Td   : derivative time constant [s]
        - alpha: derivative filter coefficient
    :omega: frequency vector [rad/s] over which the Bode diagram is computed
    :title: (optional: default value is "")
        title of the Bode diagram figure

    The function "Margin_error" returns:
        :gm_val: gain margin [dB]
            Amount of gain increase required to make the open-loop gain
            equal to 0 dB at the -180° phase crossing frequency.
            Returns np.inf if no -180° phase crossing is found.
        :pm_val: phase margin [°]
            Amount of additional phase lag required to bring the phase
            to -180° at the 0 dB gain crossover frequency.
            Returns np.inf if no 0 dB gain crossing is found.

    The function also displays a two-panel Bode plot showing:
        - Magnitude [dB] vs frequency for P(s), C(s) and L(s)
        - Phase [°] vs frequency for P(s), C(s) and L(s)
        - Gain margin annotated on the magnitude plot (green arrow)
        - Phase margin annotated on the phase plot (red arrow)
    """
    s = 1j * omega

    # 1. Calcul du contrôleur C(s)
    Kc = C.parameters['Kc']
    Ti = C.parameters['Ti']
    Td = C.parameters['Td']
    alpha = C.parameters['alpha']
    
    Cs_complex = Kc * (1 + 1 / (Ti * s) + (Td * s) / (alpha * Td * s + 1))
    gain_C_db = 20 * np.log10(np.abs(Cs_complex))
    phase_C = np.degrees(np.angle(Cs_complex))

    # 2. Calcul du processus P(s)
    bode_result = Bode(P, omega, Show=False)
    gain_P_db = 20 * np.log10(np.abs(bode_result))
    
    # Correction : On évite l'unwrap infini pour le plot, ou on limite l'axe Y après
    phase_P = np.degrees(np.unwrap(np.angle(bode_result)))

    # 3. Boucle Ouverte L(s)
    gain_L_db = gain_P_db + gain_C_db
    phase_L = phase_P + phase_C
    gain_L_linear = 10**(gain_L_db / 20)

    # --- Recherche des indices critiques ---
    # On cherche le croisement à -180
    idx_180 = np.where(np.diff(np.sign(phase_L + 180)))[0]
    idx_0dB = np.where(np.diff(np.sign(gain_L_linear - 1)))[0]

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)

    fig.suptitle(title, y=0.98)
   

    # Magnitude Plot
    ax_mag.semilogx(omega, gain_P_db, 'g--', alpha=0.4, label='P(s)')
    ax_mag.semilogx(omega, gain_C_db, 'r--', alpha=0.4, label='C(s)')
    ax_mag.semilogx(omega, gain_L_db, 'b', lw=2, label='L(s) = CP')
    ax_mag.axhline(0, color='black', label ="1dB", lw=1.5)
    
    # Phase Plot
    ax_phase.semilogx(omega, phase_P, 'g--', alpha=0.4, label='P(s)')
    ax_phase.semilogx(omega, phase_C, 'r--', alpha=0.4, label='C(s)')
    ax_phase.semilogx(omega, phase_L, 'orange', lw=2, label = 'L(s) = CP')
    ax_phase.axhline(-180, color='black',label = "-180° " ,lw=1.5, linestyle='--')

    # --- Marges ---
    gm_val = np.inf
    if len(idx_180) > 0:
        idx = idx_180[0] # Premier croisement
        w_180 = omega[idx]
        gm_val = -gain_L_db[idx]
        ax_mag.annotate('', xy=(w_180, 0), xytext=(w_180, gain_L_db[idx]),
                        arrowprops=dict(arrowstyle='<->', color='green', lw=2))
        ax_mag.text(w_180, gain_L_db[idx]/2, f' GM: {gm_val:.1f} dB', color='green', weight='bold')

    pm_val = np.inf
    if len(idx_0dB) > 0:
        idx = idx_0dB[-1] # Dernier croisement (souvent le plus critique)
        w_0dB = omega[idx]
        pm_val = phase_L[idx] + 180
        ax_phase.annotate('', xy=(w_0dB, -180), xytext=(w_0dB, phase_L[idx]),
                          arrowprops=dict(arrowstyle='<->', color='red', lw=2))
        ax_phase.text(w_0dB, (phase_L[idx]-180)/2, f' PM: {pm_val:.1f}°', color='red', weight='bold')

    # --- RÉGLAGE CRUCIAL DU ZOOM ---
    ax_phase.set_ylim([-300, 45]) # On centre sur la zone -180 / 0
    
    ax_mag.grid(True, which="both", alpha=0.3)
    ax_phase.grid(True, which="both", alpha=0.3)
    ax_mag.legend()
    ax_phase.legend()
    plt.tight_layout()
    plt.show()

    return gm_val, pm_val