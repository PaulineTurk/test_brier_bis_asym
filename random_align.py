



def write_fasta(list_standard_aa, len_align, general_filename, number_file):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
    """

    import textwrap

    for i in range(1,number_file+1):
        dictionary = {">seq_1": "".join(random.choices(list_standard_aa, k=len_align)),
                      ">seq_2": "".join(random.choices(list_standard_aa, k=len_align))}


        with open(f"{general_filename}.{i}", "w") as outfile:
            for key, value in dictionary.items():
                outfile.write(key + "\n")
                outfile.write("\n".join(textwrap.wrap(value, 60)))
                outfile.write("\n")

        print("Success! File written")

import random


list_standard_aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
general_filename = "/home/pauline/Bureau/test_brier_bis_asym/DATA_MNHN_1id_9deca_random/SEED/seed_random"
len_align = 20_000
number_file = 1


write_fasta(list_standard_aa, len_align,  general_filename, number_file)