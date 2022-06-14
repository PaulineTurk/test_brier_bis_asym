import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[2]
sys.path.append(str(package_root_directory_MNHN))


import MNHN.utils.folder as folder
import MNHN.treatment.split as split



# Split des seeds en seed train et seed test (95%/5% de Pfam) 
path_folder_fasta_nonRedondant = "/Users/pauline/Desktop/Test_preliminaire/Pfam_B"  # tous les fasta
path_folder_data = path_folder_fasta_nonRedondant      # redondant corrected
result_folder = "/Users/pauline/Desktop/Test_preliminaire/data_test_split"   # chemin à choisir
folder.creat_folder(result_folder)
path_folder_data_split = f"{result_folder}/Pfam_split" 
percentage_A = 90
name_data_A, name_data_B = "seed_test_90", "seed_test_10"
split.data_split(path_folder_data, path_folder_data_split, percentage_A, name_data_A, name_data_B)