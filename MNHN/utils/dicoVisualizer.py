import numpy as np
import pandas as pd

def dico_visualizer(path_dico):
    dico = np.load(path_dico, allow_pickle='TRUE').item()    
    df_dico = np.transpose(pd.DataFrame.from_dict(dico))  
    print(path_dico)
    print(df_dico)
    return df_dico



#dico_visualizer("/Users/pauline/Desktop/Projet_MNHN/PF00198.26.pid_v2.npy")
#dico_visualizer("/Users/pauline/Desktop/Projet_MNHN/PF00198.26.pid.npy")







# dico_visualizer("/Users/pauline/Desktop/MNHN_EvolProt/PF00244.23.pid.npy")

# dico_visualizer("/Users/pauline/Desktop/MNHN_EvolProt/PF09847.12.pid.npy")

# dico_visualizer("/Users/pauline/Desktop/MNHN_EvolProt/PF10417.12.pid.npy")

# dico_visualizer("/Users/pauline/Desktop/MNHN_EvolProt/PF10922.11.pid.npy")

# dico_visualizer("/Users/pauline/Desktop/MNHN_EvolProt/PF12574.11.pid.npy")