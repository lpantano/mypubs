import logging
from argparse import ArgumentParser

logger = logging.getLogger()
name_data = {}

class realSeq:
    def __init__(self, seq):
        self.seq = seq
        self.miraligner = ""
        self.bench = ""


def _read_miraligner(in_file):
    """read miraligner results"""
    global name_data
    seen = {}
    with open(in_file) as in_handle:
        in_handle.readline()
        for line in in_handle:
            cols = line.strip().split("\t")
            logger.debug(line)
            ok_mir, ok_add, ok_mut, ok_t5, ok_t3 = _check_pos_mirlaigner_strict(cols[1], cols[2], cols[5], cols[6], cols[7], cols[8])
            print "%s %s %s %s %s %s %s True Yes miraligner" % (cols[1], _is_iso(cols[1]), ok_mir, ok_mut, ok_add, ok_t5, ok_t3)
            logger.debug("miraligner %s %s %s %s %s %s True" % (cols[1], ok_mir, ok_mut, ok_add, ok_t5, ok_t3))
            seen[cols[0]] = 0
    for name in name_data:
        if name not in seen:
            print "%s %s False False False False False False None miraligner" % (name, _is_iso(_get_name(name)))


def _check_pos_miraligner(name, chr, mut, add, t5, t3):
    """compare name of the seq with position map
    just compatible with miraligner"""
    ok_mir, ok_add, ok_mut, ok_t5, ok_t3 = False, False, False, False, False
    cref = name.split("_")[1].lower()
    cquery = chr.lower()
    logger.debug("%s %s" % (cref, cquery))
    t5_ref, t3_ref = name.split("_")[3].split(":")
    mut_ref = name.split("_")[4].split(":")[1]
    add_ref = name.split("_")[5].split(":")[1]
    if (t3 != "0" and t3_ref != "0") or (t3 == t3_ref):
        ok_t3 = True
    if (t5 != "0" and t5_ref != "0") or (t5 == t5_ref):
        ok_t5 = True
    if mut != "0" and mut_ref != "null":
        ok_mut = True
    if mut == "0" and mut_ref == "null":
        ok_mut = True
    if (add != "null" and add_ref != "null") or (add == "0" and add_ref == "null"):
        ok_add = True
    if cref == cquery:
        ok_mir = True
    return ok_mir, ok_add, ok_mut, ok_t5, ok_t3


def _check_pos_mirlaigner_strict(name, chr, mut, add, t5, t3):
    """compare name of the seq with position map
    just compatible with miraligner"""
    ok_mir, ok_add, ok_mut, ok_t5, ok_t3 = False, False, False, False, False
    cref = name.split("_")[1].lower()
    cquery = chr.lower()
    logger.debug("%s %s" % (cref, cquery))
    t5_ref, t3_ref = name.split("_")[3].split(":")
    mut_ref = name.split("_")[4].split(":")[1]
    add_ref = name.split("_")[5].split(":")[1]
    if (len(t3) == abs(int(t3_ref))) or (t3 == t3_ref):
        ok_t3 = True
    if (len(t5) == abs(int(t5_ref))) or (t5 == t5_ref):
        ok_t5 = True
    if mut != "0" and mut_ref != "null":
        ok_mut = True
    if mut == "0" and mut_ref == "null":
        ok_mut = True
    if len(add) == len(add_ref) or (add == "0" and add_ref == "null"):
        ok_add = True
    if cref == cquery:
        ok_mir = True
    return ok_mir, ok_add, ok_mut, ok_t5, ok_t3


def _read_srnabench(in_file):
    """read srnabench results"""
    global name_data
    seen = {}
    with open(in_file) as in_handle:
        in_handle.readline()
        for line in in_handle:
            cols = line.strip().split("\t")
            logger.debug(line)
            ok_mir, ok_add, ok_mut, ok_t5, ok_t3, info = _check_pos_srnabench(_get_name(cols[0]), cols[1], cols[3], cols[4])
            print "%s %s %s %s %s %s %s %s Yes srnabench" % (_get_name(cols[0]), _is_iso(_get_name(cols[0])), ok_mir, ok_mut, ok_add, ok_t5, ok_t3, info)
            logger.debug("srnabench %s %s %s %s %s %s %s" % (cols[1], ok_mir, ok_mut, ok_add, ok_t5, ok_t3, info))
            seen[cols[0]] = 0
    for name in name_data:
        if name not in seen:
            print "%s %s False False False False False False None srnabench" % (name, _is_iso(_get_name(name)))


def _get_name(seqs):
    """get the name of the sequence"""
    global name_data
    return name_data[seqs]


def _is_iso(name):
    """Identify if it is an isomir"""
    logger.debug("_is_iso: %s" % name)
    t5_ref, t3_ref = name.split("_")[3].split(":")
    mut_ref = name.split("_")[4].split(":")[1]
    add_ref = name.split("_")[5].split(":")[1]
    if t5_ref != "0" or t3_ref != "0" or mut_ref != "null" or add_ref != "null":
        return True
    else:
        return False


def  _check_pos_srnabench(name, chr, isoclass, nucvar):
    """compare name of the seq with position map
    just compatible with srnabench"""
    ok_mir, ok_add, ok_mut, ok_t5, ok_t3, info = False, False, False, False, False, True
    t5_ref, t3_ref = name.split("_")[3].split(":")
    logger.debug("_check_pos_srnabench: name %s" % name)
    logger.debug("_check_pos_srnabench: isoclass %s" % isoclass)
    logger.debug("_check_pos_srnabench: nucvar %s" % nucvar)
    cref = name.split("_")[1].lower()
    cquery = chr.lower()
    mut_ref = name.split("_")[4].split(":")[1]
    add_ref = name.split("_")[5].split(":")[1]
    isoclass = isoclass
    if "|" in isoclass:
        var_type = isoclass.split("|")[1]
        logger.debug("_check_pos_srnabench: var_type %s" % var_type)
        if var_type.startswith("lv3p") and t3_ref != "0":
            ok_t3 = True
        if var_type.startswith("lv5p") and t5_ref != "0":
            ok_t5 = True
        if var_type.startswith("nta") and add_ref != "null":
            ok_add = True
    if not ok_t3 and t3_ref == "0":
        ok_t3 = True
    if not ok_t5 and t5_ref == "0":
        ok_t5 = True
    if isoclass == "mv" and sum([t3_ref != "0", t5_ref != "0", add_ref != "null"]) > 1:
        ok_t3, ok_t5, ok_add = True, True, True
        info = False
    if nucvar != "-" and mut_ref != "null":
        ok_mut = True
    if nucvar == "-" and mut_ref == "null":
        ok_mut = True
    if cquery == cref:
        ok_mir = True
    logger.debug("_check_pos_srnabench: info %s ok_t3 %s ok_t5 %s" % (info, ok_t5, ok_t3))
    return ok_mir, ok_add, ok_mut, ok_t5, ok_t3, info


def _create_seqs_name_dict(in_file):
    """read fasta file and populate dict"""
    global name_data
    name = ""
    with open(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith(">"):
                name_data.update({line.strip(): name})
            else:
                name = line.strip().replace(">", "")


if __name__ == "__main__":
        parser = ArgumentParser(description="check annotation")
        parser.add_argument("--log", help="debug mode", action='store_true')
        args = parser.parse_args()
        if not args.log:
            numeric_level = getattr(logging, "INFO", None)
        else:
            numeric_level = getattr(logging, "DEBUG", None)
        logging.basicConfig(level=numeric_level)
        print "name is_iso is_mir is_mut is_add is_t5 is_t3 info find tool"
        _create_seqs_name_dict("sim.21.hsa.fa")
        _read_srnabench("srnabench/srnabench/miRBase_isoAnnotation.txt")
        _read_miraligner("miraligner/sim.21.hsa.mirna")
