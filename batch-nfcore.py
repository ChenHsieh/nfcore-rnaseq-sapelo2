import os
import shutil
import pandas as pd
import numpy as np
import subprocess

# TODO - create directories for user, if they don't exist
sub_directories = ["species", "clean", "count",
                   "fastq", "map", "miscelaneous", "reference"]

table = pd.read_csv("./input/sample_table.csv")
directory_ids = list(set(table["folder_id"].to_numpy()))
myID = table["myID"][0]
start_sh_list = []


def make_directories():
    """
    Creates directory structure to store genome data

    @param species_name: the species to sequence
    """

    print("### creating directories ###")
    # TODO - table["species"] can be used to write readme files for human to read in each directory
    species = table["species"]
    root_dir = "./data/"

    if not (os.path.isdir(root_dir)):
        os.mkdir(root_dir)
    for id in directory_ids:
        if not (os.path.isdir(root_dir + id)):
            print(f"creating directory " + id)
            os.mkdir(root_dir + id)
        for dir in sub_directories:
            if not (os.path.isdir(root_dir + id + '/' + dir)):
                print(f"creating {root_dir + id + '/' + dir}")
                os.mkdir(root_dir + id + "/" + dir + "/")


def transfer_data():
    """
    As of now, extracts data from Phytozome and
    moves it to specific species directory

    @param species_name: the species to sequence
    """
    print("### extracting data ###")
    # copy genome reference, fasta and gff3 respectively
    for i in range(len(table)):
        if not os.path.isfile("./data/" + table.at[i, "folder_id"] + "/reference/" + table.at[i, "folder_id"] + "_genome.fa"):
            print(
                f'copying {table.at[i, "genome_fa"]} to {"./data/" + table.at[i,"folder_id"] + "/reference/" + table.at[i, "folder_id"] + "_genome.fa"}')
            shutil.copyfile(table.at[i, "genome_fa"], "./data/" + table.at[i,
                                                                           "folder_id"] + "/reference/" + table.at[i, "folder_id"] + "_genome.fa")
        # copy gff3, must change the extension to gff
        if not os.path.isfile("./data/" + table.at[i, "folder_id"] + "/reference/" + table.at[i, "folder_id"] + ".gff3"):
            print(
                f'copying {table.at[i, "gene_gff3"]} to {"./data/" + table.at[i,"folder_id"] + "/reference/" + table.at[i, "folder_id"] + ".gff"}')
            shutil.copyfile(table.at[i, "gene_gff3"], "./data/" + table.at[i,
                                                                           "folder_id"] + "/reference/" + table.at[i, "folder_id"] + ".gff")
        # copy fastq files
        print(
            f'copying {table.at[i, "fastq_dir"] + table.at[i, "sample_name"] + ".R1.fastq"}')
        shutil.copyfile(table.at[i, "fastq_dir"] + table.at[i, "sample_name"] + ".R1.fastq",
                        "./data/" + table.at[i, "folder_id"] + "/fastq/" + table.at[i, "sample_name"] + ".R1.fastq")
        shutil.copyfile(table.at[i, "fastq_dir"] + table.at[i, "sample_name"] + ".R2.fastq",
                        "./data/" + table.at[i, "folder_id"] + "/fastq/" + table.at[i, "sample_name"] + ".R2.fastq")




def generate_scripts():
    print("### generating scripts ###")
    for element in directory_ids:
        """
        makes a script to prepare for array mapping
        """
        element = element.strip()
        start_sh_list.append(f"cd ./data/{element}")
        start_sh_list.append(f"sbatch run_nfcore_{element}.sh")
        start_sh_list.append("cd ../..")
        with open(f"./data/{element}/run_nfcore_{element}.sh", "wt") as script:
            script.write("#!/bin/sh\n")
            script.write(f"#SBATCH -J {element}_nfcore\n")
            script.write("#SBATCH --partition highmem_p\n")
            script.write("#SBATCH --nodes=1\n")
            script.write("#SBATCH --ntasks-per-node=12\n")
            script.write("#SBATCH --time=8:00:00\n")
            script.write("#SBATCH --mem=130gb\n")
            script.write(f"#SBATCH --mail-user={myID}@uga.edu\n")
            script.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
            script.write("date\n\n")
            script.write("cd $SLURM_SUBMIT_DIR\n")

            script.write("ml Nextflow\n")
            # TODO correct the parameters
            script.write(f"""nextflow run nf-core/rnaseq \
--input sample_sheet.csv \
--fasta ./reference/{element}_genome.fa \
--gff ./reference/{element}.gff \
-profile singularity --max_memory 128GB""")
#  --fc_group_features 'gene_id' was removed temporary
        with open(f"./data/{element}/sample_sheet.csv", 'w') as sample_sheet:
            sample_sheet.write("sample,fastq_1,fastq_2,strandedness\n")
            sample_list = table[table["folder_id"] == element].sample_name.to_list()
            count = 1
            for sample in sample_list:
                sample_sheet.write(
                    f"{element}_rep{count},./fastq/{element}.R1.fastq,./fastq/{element}.R2.fastq,unstranded\n")
                count += 1


def submit_scripts():
    """
    Submit all scripts to the cluster
    """
    print("### Submitting scripts to the cluster ###")
    with open(f"submit_all_scripts.sh", "w") as submit_script:
        submit_script.write("\n".join(start_sh_list))
    bashCommand = f"bash ./submit_all_scripts.sh"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print(f"subprocess output: {output}")
    print(f"subprocess error: {error}")


make_directories()
transfer_data()
generate_scripts()
submit_scripts()
