import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[1]
root_path = file.parents[2]
sys.path.append(str(package_root_directory_MNHN))


import MNHN.utils.folder as folder
import MNHN.treatment.description as description
import MNHN.treatment.stockholm as stockholm
import MNHN.treatment.capitalizer as capitalizer
import MNHN.treatment.split as split
import MNHN.treatment.pid as pid
import MNHN.treatment.redundancy as redundancy


#DATA_MNHN = "DATA_MNHN_9id_1deca"
#DATA_MNHN = "DATA_MNHN_10id"
#DATA_MNHN = "DATA_MNHN_1id_9deca_asym_v2"
#DATA_MNHN = "DATA_MNHN_1id_9deca_random"
#DATA_MNHN = "DATA_MNHN_5id_5deca"
#DATA_MNHN = "DATA_MNHN_10deca"
#DATA_MNHN = "DATA_MNHN_1id_9deca_mini"
#DATA_MNHN = "Mini_Pfam"
DATA_MNHN = "DATA_MNHN_1id_9deca"

# DATA_MNHN = "DATA_MNHN_1id_2deca_short"
# list_standard_aa = ["A", "R", "N", "D"]
# list_residu = ["A", "R", "N", "D"]

# Séparation multi-Stockholm en un fichier Fasta par seed
list_standard_aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

#list_standard_aa = ["A", "Y", "V"]

# path_file = "SEED_majorite_identite"  # seed from Pfam in multi_Stockholm

# path_file = f"{root_path}/DATA_MNHN/Pfam-A.seed"
# path_folder_save = "DATA_MNHN/Pfam_Stockholm" # chemin à choisir
# stockholm.stockholm_separator(path_file, path_folder_save)
# path_folder_stockholm = path_folder_save

path_folder_fasta = f"{DATA_MNHN}/SEED"  # chemin à choisir
# stockholm.multi_stockholm_to_fasta(path_folder_stockholm, path_folder_fasta)
path_data = path_folder_fasta
# print(f"\nVisualisation of {path_data}")
# list_residu = list_standard_aa
# residu_count, total_residu, character_count, total_character = description.data_count(path_data, list_residu)
# description.bar_plot_data_count(path_data, residu_count, total_residu, "Standard amino-acid")
# description.bar_plot_data_count(path_data, character_count, total_character , "Character")


# Capitalisation des caractères de Pfam
path_data_corrected = f"{DATA_MNHN}/Pfam_Upper"  # chemin à choisir
print(path_data_corrected)
capitalizer.multi_capitalization(path_data, path_data_corrected)
path_data = path_data_corrected
print(f"\nVisualisation of {path_data}")
list_residu = list_standard_aa
residu_count, total_residu, character_count, total_character = description.data_count(path_data, list_residu)
description.bar_plot_data_count(path_data, residu_count, total_residu, "Standard amino-acid")
description.bar_plot_data_count(path_data, character_count, total_character , "Character")


# Calcul des PID
print("\nCalcul de PID")
path_folder_fasta = path_data_corrected # upper correction only
path_folder_pid = f"{DATA_MNHN}/PID"    # chemin à choisir
list_inclusion = list_standard_aa
pid.save_pid(path_folder_fasta, path_folder_pid, list_inclusion)

# Clustering
# print("\nClustering")
# path_folder_fasta = path_data_corrected # still upper correction only
# path_folder_fasta_nonRedondant = f"{root_path}/DATA_MNHN/Pfam_nonRedondant"  # chemin à choisir
# list_residu = list_standard_aa
# redundancy.multi_non_redundancy_correction(path_folder_fasta, path_folder_fasta_nonRedondant, list_residu, pid_sup = 99)
# residu_count, total_residu, character_count, total_character = description.data_count(path_folder_fasta_nonRedondant, list_residu)
# description.bar_plot_data_count(path_folder_fasta_nonRedondant, residu_count, total_residu, "Standard amino-acid")
# description.bar_plot_data_count(path_folder_fasta_nonRedondant, character_count, total_character , "Character")

# Split des seeds en seed train et seed test (50%/50% de Pfam) 
# path_folder_data = path_data_corrected     # data to split train/test
# result_folder = f"{DATA_MNHN}/data_Result"   # chemin à choisir
# folder.creat_folder(result_folder)
# path_folder_data_split = f"{DATA_MNHN}/data_Result/Pfam_split" 
# folder.creat_folder(path_folder_data_split)
# percentage_A = 50
# name_data_A, name_data_B = "Pfam_train", "Pfam_test"
# split.data_split(path_folder_data, path_folder_data_split, percentage_A, name_data_A, name_data_B)