import os
from collections import Counter

H = "seq known t5 t3 mut add correct mapped amb"


def read_bam(fn_name):
    counts = Counter()
    data = {}
    pairs = {}
    with open(fn_name) as in_handle:
        for line in in_handle:
            cols = line.strip().split("\t")
            attr = cols[20].split(";")
            mir_ann = attr[2].split("=")[1]
            is_correct = False
            if cols[14] == "miRNA":
                ann = get_name(cols[3])
                if ann[0] == mir_ann:
                    is_correct = True
                ann = " ".join(ann)
                pair = (cols[3], mir_ann)
                if pair not in pairs:
                    counts[cols[3]] += 1
                pairs[pair] = 0
                data[(cols[3], is_correct)] = "%s %s %s %s" % (cols[3], ann, str(is_correct), mir_ann)
    return data, counts


def print_output(data, counts, out_file):
    with open(out_file, 'w') as out_handle:
        out_handle.write(H + "\n")
        for k in counts:
            if counts[k] > 1 and (k, True) in data:
                out_handle.write("%s %s\n" % (data[(k, True)],counts[k]))
            elif counts[k] > 1:
                out_handle.write("%s %s\n" % (data[(k, False)],counts[k]))
            elif counts[k] == 1 and (k, True) in data:
                out_handle.write("%s %s\n" % (data[(k, True)],counts[k]))
            elif counts[k] == 1:
                out_handle.write("%s %s\n" % (data[(k, False)],counts[k]))


def get_name(name):
    f = name.strip().split("_")
    is_add, is_mut, is_5trim, is_3trim = False, False, False, False
    if not f[4].endswith("null"):
        is_mut = True
    if not f[5].endswith("null"):
        is_add = True
    if not f[3].endswith("0"):
        is_3trim = True
    if not f[3].startswith("0"):
        is_5trim = True
    return [f[1], str(is_5trim), str(is_3trim), str(is_mut), str(is_add)]


def main():
    data, counts = read_bam("sim.20.hsa.bam.anno")
    print_output(data, counts,  "sim.20.hsa.dat")
    data, counts = read_bam("sim.20.hsa.primary.bam.anno")
    print_output(data, counts,  "sim.20.hsa.primary.dat")


if __name__ == "__main__":
    main()
