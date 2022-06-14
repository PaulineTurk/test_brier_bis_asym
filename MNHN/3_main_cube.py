import sys  
from pathlib import Path  
file = Path(__file__). resolve()  
package_root_directory_MNHN = file.parents [1]
root_path = file.parents[2]
sys.path.append(str(package_root_directory_MNHN))
import MNHN.localNeighbour.localNeighbourfonction as localNeighbourfonction
import MNHN.utils.folder as folder

#DATA_MNHN = "DATA_MNHN_9id_1deca"
#DATA_MNHN = "DATA_MNHN_10id"
#DATA_MNHN = "DATA_MNHN_1id_9deca"
#DATA_MNHN = "DATA_MNHN_1id_9deca_random"
#DATA_MNHN = "DATA_MNHN_5id_5deca"
#DATA_MNHN = "DATA_MNHN_10deca"
#DATA_MNHN = "DATA_MNHN_1id_9deca_mini"
#DATA_MNHN = "Mini_Pfam"
DATA_MNHN = "DATA_MNHN_1id_9deca"

# DATA_MNHN = "DATA_MNHN_1id_2deca_short"
# list_standard_aa = ["A", "R", "N", "D"]
# list_residu = ["A", "R", "N", "D"]

path_folder_pid = f"{DATA_MNHN}/PID" # chemin des pid
path_folder_fasta = f"{DATA_MNHN}/Pfam_Upper"  # chemin des seeds d'entrainement
path_NeighborRes = f"{DATA_MNHN}/data_Result/CUBE"
folder.creat_folder(path_NeighborRes) # folder a créer avant de lancer le calcul des cubes

list_residu = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

#list_residu = ["A", "Y", "V"]

pid_inf = 0


for delay_num in range(1, 5):
    for kp_SeqChoice in ["k", "p"]:
        print(f"\n{delay_num}, {kp_SeqChoice }\n")
        triplet_count, name_folder_fasta = localNeighbourfonction.multi_triplet_count_for_cube(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf)
        cond_proba, path_proba_cond = localNeighbourfonction.cube(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu)

for delay_num in [-1, -2, -3, -4]:
    for kp_SeqChoice in ["k", "p"]:
        print(f"\n{delay_num}, {kp_SeqChoice }\n")
        triplet_count, name_folder_fasta = localNeighbourfonction.multi_triplet_count_for_cube(path_folder_fasta, path_folder_pid, delay_num, kp_SeqChoice, list_residu, pid_inf)
        cond_proba, path_proba_cond = localNeighbourfonction.cube(triplet_count, name_folder_fasta, path_NeighborRes, delay_num, kp_SeqChoice, list_residu)

