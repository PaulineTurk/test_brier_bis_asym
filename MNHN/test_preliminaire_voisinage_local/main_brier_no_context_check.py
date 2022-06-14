import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]  # 0: meme niveau, 1: 1 niveau d'écart etc.
sys.path.append(str(package_root_directory_MNHN))


import MNHN.brier.brierScore as brierScore
import MNHN.brier.predictor as predictor

path_folder_fasta = "/Users/pauline/Desktop/Test_preliminaire/data_test_split/Pfam_split/seed_test_10" # les environs 1000 seeds tests
path_folder_pid = "/Users/pauline/Desktop/Overfitting_test/PID_couple"  # tous les anciens pid
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]



path_cond_proba = "/Users/pauline/Desktop/Overfitting_test/test_1/BlosumRes/BlosumRes_50_A/Blosum_proba_cond.npy" # l'ancienne variante blosum proba
unit_Brier = predictor.predictor_blosum(path_cond_proba, list_residu)
print("\npredictor Blosum")
brier_score = brierScore.brier_score_matrix(path_folder_fasta, path_folder_pid, unit_Brier, list_residu, pid_inf = 62)
print(brier_score)
