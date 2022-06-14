import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [2]
root_path = file.parents[3]
sys.path.append(str(package_root_directory_MNHN))
import MNHN.localNeighbour.localNeighbourfonction as localNeighbourfonction

path_folder_pid = "DATA_MNHN/PID" # chemin des pid
path_folder_fasta = "DATA_MNHN/Pfam_Upper"  # chemin des seeds d'entrainement
path_NeighborRes = "DATA_MNHN/data_Result/Pfam_split/CUBE" # folder a créer avant de lancer le calcul des cubes
delay_num = -1
kp_SeqChoice = "k"
list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
pid_inf = 0
triplet_count, name_folder_fasta = localNeighbourfonction.multi_triplet_count_for_cube(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf)
cond_proba, path_proba_cond = localNeighbourfonction.cube(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu)