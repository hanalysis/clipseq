import glob
import pandas as pd

log_files = glob.glob("*tetranscripts.log")

for file in log_files:
    with open(file, "rt") as logfile:
        contents = logfile.read()

ind = 1
# don't include the control lib
final_ind = contents.find("In library test.rna.paired_end")

ann_reads_dict = {}
multi_reads_dict = {}
unann_reads_dict = {}

while contents.find("In library ", ind) < final_ind:
    sample_1 = contents.find("In library ", ind)
    sample_2 = contents.find(".bam", sample_1)
    sample_name = contents[(sample_1 + 11):(sample_2 - 18)]

    ann_reads = contents.find("Total annotated reads = ", ind)
    multi_reads = contents.find("Total non-uniquely mapped reads = ", ind)
    unann_reads = contents.find("Total unannotated reads = ", ind)

    ann_reads_dict[sample_name] = ann_reads + len("Total annotated reads = ")
    multi_reads_dict[sample_name] = multi_reads + len("Total non-uniquely mapped reads = ")
    unann_reads_dict[sample_name] = unann_reads + len("Total unannotated reads = ")

    # move onto next sample
    ind = contents.find("In library ", sample_2)

# Annotated reads
ann_reads_val_dict = {}

for sample, pos in ann_reads_dict.items():
    i = 0
    while contents[pos+i].isdigit():
        i = i + 1
    if not contents[pos+i].isdigit():
        ann_reads_val_dict[sample] = contents[pos:pos+i]

# Multi reads
multi_reads_val_dict = {}

for sample, pos in multi_reads_dict.items():
    i = 0
    while contents[pos+i].isdigit():
        i = i + 1
    if not contents[pos+i].isdigit():
        multi_reads_val_dict[sample] = contents[pos:pos+i]

#Unann reads
unann_reads_val_dict = {}

for sample, pos in unann_reads_dict.items():
    i = 0
    while contents[pos+i].isdigit():
        i = i + 1
    if not contents[pos+i].isdigit():
        unann_reads_val_dict[sample] = contents[pos:pos+i]

## Put into df

df = pd.DataFrame({
    "Total annotated reads": ann_reads_val_dict,
    "Total multi-mapping reads": multi_reads_val_dict,
    "Total unannotated reads": unann_reads_val_dict
})


df.to_csv("tetranscripts_qc.tsv", sep = "\t")