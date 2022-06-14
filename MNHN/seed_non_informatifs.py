import sys  
from pathlib import Path  
file = Path(__file__).resolve()  
package_root_directory_MNHN = file.parents[1]
sys.path.append(str(package_root_directory_MNHN))


import MNHN.utils.fastaReader as fastaReader
from MNHN.utils.timer import Timer


def nbreSeed(path_folder_fasta):   # étape de pré-traitement à ajouter pour éviter certains calculs inutils à des étapes suivantes
    
    t = Timer()
    t.start()
    
    files_in_path_folder_fasta = Path(path_folder_fasta).iterdir()

    print("\nPath folder fasta tested:", path_folder_fasta)

    #list_non_informative_seed = []
    count = 0
    print("\nNon informative seeds:")

    for path_file_fasta in files_in_path_folder_fasta:
        len_seed = len(fastaReader.read_multi_fasta(path_file_fasta))
        if len_seed <= 1:
                print(f"{path_file_fasta}, len = {len_seed}")
                count += 1
                #list_non_informative_seed.append(path_file_fasta)
                

    print("\nTotal non informative seeds:", count)
    t.stop("Identification of non-informative seeds")

            








if __name__ == '__main__': 
    path_folder_fasta = "/Users/pauline/Desktop/Overfitting_test/Pfam_fasta_99"
    nbreSeed(path_folder_fasta)