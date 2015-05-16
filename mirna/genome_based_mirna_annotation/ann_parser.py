from collections import Counter
import pysam

H = "seq known t5 t3 mut add correct mapped amb tool"


def read_sim_fa(fn):
    sim = {}
    with open(fn) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                name = line.strip()[1:]
                sim[name] = get_name(name)
    return sim

def read_bam(fn_name):
    mapped = {}
    fn = pysam.AlignmentFile(fn_name, 'rb')
    for record in fn:
        mapped[record.query_name] = get_name(record.query_name)
    return mapped

def read_ann(fn_name):
    counts = Counter()
    data = {}
    pairs = {}
    with open(fn_name) as in_handle:
        for line in in_handle:
            cols = line.strip().split()
            # attr = cols[20].split(";")
            mir_ann = cols[15]
            is_correct = False
            if cols[16] == "miRNA":
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

def print_output(data, counts, sim, mapped, out_file, prefix):
    with open(out_file, 'w') as out_handle:
        # out_handle.write(H + "\n")
        for k in counts:
            if counts[k] > 1 and (k, True) in data:
                out_handle.write("%s %s %s\n" % (data[(k, True)], counts[k], prefix))
            elif counts[k] > 1:
                out_handle.write("%s %s %s\n" % (data[(k, False)], counts[k], prefix))
            elif counts[k] == 1 and (k, True) in data:
                out_handle.write("%s %s %s\n" % (data[(k, True)], counts[k], prefix))
            elif counts[k] == 1:
                out_handle.write("%s %s %s\n" % (data[(k, False)], counts[k], prefix))

        for mirna in mapped:
            if mirna not in counts:
                out_handle.write("%s %s False Non-mirna 0 %s\n" % (mirna, " ".join(sim[mirna]), prefix))
                counts[mirna] = 0

        for mirna in sim:
            if mirna not in counts:
                out_handle.write("%s %s False False 0 %s\n" % (mirna, " ".join(sim[mirna]), prefix))


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
