import json
import numpy as np

# Nested dictionary
path_2d = "/home/pauline/Bureau/test_brier_bis_asym/DATA_MNHN_1id_9deca/data_Result/Blosum_proba_cond.npy"
table_2d = np.load(path_2d, allow_pickle='TRUE').item()
print(json.dumps(table_2d , indent=4))
