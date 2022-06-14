import numpy as np
import pandas as pd



########### 
# blosum part
###########
dico_blosum = {"A":{"A": 0.35, "C": 0.65},
        "C":{"A": 0.15, "C":0.85}}
blosum_path = "/Users/pauline/Desktop/cube_manu/blosum_test"
np.save(blosum_path, dico_blosum)

def dico_visualizer(path_dico):
    dico = np.load(path_dico, allow_pickle='TRUE').item()    
    df_dico = np.transpose(pd.DataFrame.from_dict(dico))  
    print(path_dico)
    print(df_dico)
    return df_dico

dico_visualizer("/Users/pauline/Desktop/cube_manu/blosum_test.npy")


########### 
# cubes part
###########

cube_1 = {"A":{"A": {"A": 0.2, "C": 0.3},  "C": {"A": 0.8, "C": 0.7}},   # REVOIR COMMENT MES DICO ONT ÉTÉ CONSTUITS, possible erreur?
          "C":{"A": {"A": 0.1, "C": 0.4},  "C": {"A": 0.9, "C": 0.6}}}
cube_1_path = "/Users/pauline/Desktop/cube_manu/proba_cond_(-1,k)"
np.save(cube_1_path , cube_1)
print(cube_1["A"]["C"]["C"])
dico_visualizer("/Users/pauline/Desktop/cube_manu/proba_cond_(-1,k).npy")

cube_2 = {"A":{"A": {"A": 0.1, "C": 0.2},  "C": {"A": 0.9, "C": 0.8}},
          "C":{"A": {"A": 0.5, "C": 0.1},  "C": {"A": 0.5, "C": 0.9}}}
cube_2_path = "/Users/pauline/Desktop/cube_manu/proba_cond_(-2,k)"
np.save(cube_2_path , cube_2)
dico_visualizer("/Users/pauline/Desktop/cube_manu/proba_cond_(-2,k).npy")


########### 
# seed part
###########

path_Seed = "/Users/pauline/Desktop/cube_manu/Seed"

########### 
# faux pid part
###########
pid_dico = {"prot1":{"prot1": 100, "prot2": 100},
            "prot2":{"prot1": 100, "prot2": 100}}
pid_path = "/Users/pauline/Desktop/cube_manu/pid"
np.save(pid_path, pid_dico)
dico_visualizer("/Users/pauline/Desktop/cube_manu/pid.npy")
