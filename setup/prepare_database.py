#
# %%
import os
path = "/home/mngs/genbank_viruses_dehydrated/ncbi_dataset/data"
combined_fasta_path = os.path.join(path, "combined.fasta")
seq2genome_map = {}
# find and cat
for genome_acc in os.listdir(path):
    if "GC" in genome_acc:
        for fasta_file_name in os.listdir(os.path.join(path, genome_acc)):
            
            with open(os.path.join(path, genome_acc, fasta_file_name)) as f:
                    contents = f.read()
            with open(combined_fasta_path, "a") as combined_fasta:
                combined_fasta.write(contents)
            for line in contents.split("\n"):
                if line.startswith(">"):
                    seq_id = line[1:].split(" ")[0]
                    seq2genome_map[seq_id] = genome_acc
# %%
#write seq2genome_map as json
import json
with open(os.path.join(path, "seq2genome_map.json"), "w") as f:
    json.dump(seq2genome_map, f)
# %%
