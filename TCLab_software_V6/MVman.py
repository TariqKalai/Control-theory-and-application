import os
import pandas as pd
from datetime import datetime

def add_man_mvman_to_file(input_path, scenario, MV0, TSim):
    """
    Lit un fichier CSV existant, ajoute les colonnes Man et MVMan
    en fonction du scénario, et sauvegarde un nouveau fichier nommé
    scenario + nom_original.
    
    :input_path: chemin vers le fichier CSV original
    :scenario: "CLP+FF", "CLP+noFF", "OLP+noFF", "OLP+FF"
    :MV0: valeur initiale de MV
    :TSim: temps de simulation total
    """

    # Lecture du fichier
    df = pd.read_csv(input_path)
    t = df['t'].values

    # Définition des paths selon le scénario
    ManPaths = {
        "CLP+FF":   {0: True,  150: False, TSim: False},
        "CLP+noFF": {0: True,  150: False, TSim: False},
        "OLP+noFF": {0: True,  TSim: True},
        "OLP+FF":   {0: True,  TSim: True},
    }

    MVManPaths = {
        "CLP+FF":   {0: MV0 + 15, TSim: MV0 + 15},
        "CLP+noFF": {0: MV0 + 15, TSim: MV0 + 15},
        "OLP+noFF": {0: MV0,      TSim: MV0},
        "OLP+FF":   {0: MV0,      TSim: MV0},
    }

    if scenario not in ManPaths:
        raise ValueError(f"Scénario inconnu : '{scenario}'. "
                         f"Choisir parmi {list(ManPaths.keys())}")

    # Construction des colonnes Man et MVMan
    man_col   = []
    mvman_col = []

    for time in t:
        man_val   = None
        mvman_val = None
        for key in ManPaths[scenario]:
            if time >= key:
                man_val = ManPaths[scenario][key]
        for key in MVManPaths[scenario]:
            if time >= key:
                mvman_val = MVManPaths[scenario][key]
        man_col.append(man_val)
        mvman_col.append(mvman_val)

    df['Man']   = man_col
    df['MVMan'] = mvman_col

    # Réordonnement des colonnes selon le header demandé
    cols = ['t', 'MV', ' MVP', ' MVI', ' MVD', ' SP', 'PV', 'DV', 'Man', 'MVMan']
    # Garde uniquement les colonnes disponibles dans le bon ordre
    cols = [c for c in cols if c in df.columns]
    print(cols)
    df = df[cols]

    # Construction du nom de sortie : scenario + nom_original
    input_dir      = os.path.dirname(input_path)
    input_filename = os.path.basename(input_path)
    output_filename = f"{scenario}_{input_filename}"
    output_path     = os.path.join("Code_data", output_filename)

    df.to_csv(output_path, index=False)
    print(f"Fichier sauvegardé : {output_path}")

    return output_path


add_man_mvman_to_file(
    input_path = "Data/TCLAB_PID_FF_Test_on_2026-04-01-09h10.txt",
    scenario   = "CLP+noFF",
    MV0        = 50,
    TSim       = 1700
)
# → CLP+FF_TCLAB_PID_FF_Test_on_2026-04-01-09h10.csv