from pymongo import MongoClient
from collections import Counter
import pandas as pd

# -------------------------
# MongoDB connection
# -------------------------
client = MongoClient("mongodb://localhost:27017/")

GAMDB = client["GAM_w15_ncc"]
GANMDB = client["GANM_w25_ncc"]
LPMDB = client["LPM_w15_ncc"]
LPNMDB = client["LPNM_w20_ncc"]
OPMDB = client["OPM_w30_ncc"]
OPNMDB = client["OPNM_w20_ncc"]

ConservedDB = client["Conserved_cluster_inter"]
conserved_col = ConservedDB["inter_specie_ncc"]

databases = [GAMDB, GANMDB, LPMDB, LPNMDB, OPMDB, OPNMDB]

AA = list("ACDEFGHIKLMNPQRSTVWY")

results = []

# -------------------------
# Iterate through conserved clusters
# -------------------------
for conserved_doc in conserved_col.find():

    cluster = conserved_doc["cluster"]

    sequences = []

    for db in databases:
        for specie in db.list_collection_names():

            doc = db[specie].find_one({"cluster": cluster})

            if doc and "sequence" in doc:
                sequences.append(doc["sequence"])

    if not sequences:
        continue

    # concatenate all sequences belonging to this cluster
    merged_seq = "".join(sequences)

    aa_count = Counter(merged_seq)
    total = len(merged_seq)

    row = {
        "cluster": cluster,
        "n_sequences": len(sequences),
        "total_length": total
    }

    for aa in AA:
        row[aa] = round(100 * aa_count.get(aa, 0) / total, 3)

    results.append(row)

# -------------------------
# Save results
# -------------------------
df = pd.DataFrame(results)

df.to_csv("amino_acid_composition.csv", index=False)

print(df.head())
print(f"\nProcessed {len(df)} clusters")
print("Results saved to amino_acid_composition.csv")