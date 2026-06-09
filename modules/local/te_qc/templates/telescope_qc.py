import glob
import pandas as pd

log_files = glob.glob("*telescope.log")

phrase_list = [" fragments mapped to reference",
               " had one unique alignment",
               " had multiple alignments",
               " fragments overlapped annotation",
               " map to one locus",
               " map to multiple loci"]

assign_list = ["Fragments mapped to reference",
               "Unique alignments",
               "Multiple alignments",
               "Fragments mapped to TE annotation",
               "One locus",
               "Multiple loci"]

extract_dict = dict(zip(phrase_list, assign_list))

# Empty rows for populating later
ref_row = []
ann_row = []

# Find position of each number need to extract
for file in log_files:
    with open(file, "rt") as logfile:
        contents = logfile.read()

    # make new dict for positions
    positions_dict = {}

    for phrase, label in extract_dict.items():
        positions_dict[label] = contents.find(phrase)


    # Extract number for each
    value_dict = {}

    for label, pos in positions_dict.items():
        i = 1
        while contents[pos -i].isdigit():
            i = i + 1
        if not contents[pos -i].isdigit():
            value_dict[label] = contents[pos-i:pos]

    print(value_dict)

    # Make one .tsv file for each of the two metrics

    all_df = pd.DataFrame({"label": list(value_dict.keys()), "value": list(value_dict.values())})

    ref_df = all_df[0:3]
    ann_df = all_df[3:6]

    # Export to .tsv

    sample_name = file.replace(".log", "")
    ref_row.append({"Sample": sample_name, **{label: value_dict[label] for label in assign_list[:3]}})

    ann_row.append({"Sample": sample_name, **{label: value_dict[label] for label in assign_list[3:]}})

ref_df = pd.DataFrame(ref_row).set_index("Sample")
ann_df = pd.DataFrame(ann_row).set_index("Sample")

ref_df.to_csv("ref_genome_telescope_qc.tsv", sep="\t")
ann_df.to_csv("te_ann_telescope_qc.tsv", sep="\t")
