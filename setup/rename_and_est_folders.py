read_dir = "/media/mngs/48TBRAID5HDD/20261008_230815/Fastq"
samples = [
           "RSV01",
           "RSV02",
           "RSV03",
           "RSV04",
           "RSV05",
           "RSV06",
           "RSV07",
           "RSV08",
           "RSV09",
           "RSV10",
           "RSV11",
           "RSV12"
]

new_dir = "/media/mngs/48TBRAID5HDD/20261008_230815_Renamed"

#make directories for the new folder structure
import os
import shutil
os.makedirs(new_dir, exist_ok=True)
for sample in samples:
    os.makedirs(os.path.join(new_dir, sample), exist_ok=True)


#parse and identify read files and move them
for file in os.listdir(read_dir):
    if "fastq.gz" in file:
        sample_name = file.split("_")[0]
        if sample_name in samples:
            read_1 = "R1" in file
            read_2 = "R2" in file
            if read_1:
                new_file_path = os.path.join(new_dir, sample_name, f"{sample_name}_1.fastq.gz")
                old_file_path = os.path.join(read_dir, file)
            if read_2:
                new_file_path = os.path.join(new_dir, sample_name, f"{sample_name}_2.fastq.gz")
                old_file_path = os.path.join(read_dir, file)
            shutil.copyfile(old_file_path, new_file_path)