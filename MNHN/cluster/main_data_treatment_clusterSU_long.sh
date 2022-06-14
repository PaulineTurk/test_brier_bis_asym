#!/bin/bash

#PBS -q beta
#PBS -l select=2:ncpus=24:mpiprocs=24
#PBS -l walltime=60:00:00
#PBS -N pid
module load python/3.9

cd /scratchbeta/turkp

path_file=(/scratchbeta/turkp/Pfam_Sample)                               # chemin fichier Stockholm téléchargé depuis Pfam
path_folder_save=(/scratchbeta/turkp/Pfam_Stockholm)                     # chemin a fixer ou mettre les fichier stockholms (1/seed)
path_folder_fasta=(/scratchbeta/turkp/Pfam_Fasta)                        # chemin ou mettre les fichiers fasta (1/seed)
path_data_corrected=(/scratchbeta/turkp/Pfam_Upper)                      # chemin ou mettre les fichiers fasta avec tous les caractères en majuscule
path_folder_pid=(/scratchbeta/turkp/PID)                                 # chemin ou mettre les numpy de pid de chaque seed
path_folder_fasta_nonRedondant=(/scratchbeta/turkp/Pfam_nonRedondant)    # chemin ou mettre les seed après clustering et retrait de redondance
result_folder=(/scratchbeta/turkp/data_Result)                           # chemin ou mettre Pfam splittée à 50/50 par rapport aux seeds entre train/test



python ./Stage_MNHN_EvolProt/MNHN/cluster/main_data_treatment_cluster.py $path_file $path_folder_save $path_folder_fasta $path_data_corrected $path_folder_pid $path_folder_fasta_nonRedondant $result_folder