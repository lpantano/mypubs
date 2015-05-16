"""
Run multiple tools to measure accuracy of isomiRs genome based annotation
"""

from argparse import ArgumentParser
import os
import contextlib

from os.path import exists as is_there
from os.path import abspath as full
from bcbio.provenance import do
from bcbio.utils import safe_makedir
import ann_parser as res

@contextlib.contextmanager
def ch_directory(dir):
    cur_dir = os.getcwd()
    os.chdir(dir)
    yield dir
    os.chdir(cur_dir)


def _stats(ann, fasta, prefix):
    output = prefix + ".tsv"
    if not is_there(output):
        sim_data = res.read_sim_fa(fasta)
        data, counts = res.read_bam(ann)
        res.print_output(data, counts, sim_data, output, prefix)
    return output


def _annotate(input, mirbase):
    output = "mirbase.bed"
    cmd = ("bedtools intersect -bed -wo -s -f 0.80 -abam"
           " {input} -b {mirbase} >| {output}")
    if not is_there(output):
        do.run(cmd.format(**locals()), "")
    return full(output)


def _star(input, index, mirbase):
    safe_makedir("star")
    with ch_directory("star"):
        output = "star_map.bam"
        cmd = ("STAR --genomeDir {index} --readFilesIN {input}"
               " --outFilterMultimapNmax 50"
               " --outSAMattributes NH HI NM"
               " --alignIntronMax 1")
        cmd_bam = "samtools view -Sbh Aligned.out.sam >| {output}"

        if not is_there("Aligned.out.sam"):
            do.run(cmd.format(**locals()), "")
        if not is_there(output):
            do.run(cmd_bam.format(**locals()), "")

        mirbase_output = _annotate(output, mirbase)
    return mirbase_output


def _bowtie2(input, index, mirbase):
    safe_makedir("bowtie2")
    with ch_directory("bowtie2"):
        output = "bowtie2_map.bam"
        cmd = ("bowtie2 -f -k 50 -L 18 -x {index}"
               " -U {input}"
               " >| hits.sam")
        cmd_bam = "samtools view -Sbh hits.sam >| {output}"

        if not is_there("hits.sam"):
            do.run(cmd.format(**locals()), "")
        if not is_there(output):
            do.run(cmd_bam.format(**locals()), "")

        mirbase_output = _annotate(output, mirbase)
    return mirbase_output


def _hisat(input, index, mirbase):
    safe_makedir("hisat")
    with ch_directory("hisat"):
        output = "hisat_map.bam"
        cmd = ("hisat -f -k 50 -L 18 -x {index}"
               " -U {input}"
               " >| hits.sam")
        cmd_bam = "samtools view -Sbh hits.sam >| {output}"

        if not is_there("hits.sam"):
            do.run(cmd.format(**locals()), "")
        if not is_there(output):
            do.run(cmd_bam.format(**locals()), "")

        mirbase_output = _annotate(output, mirbase)
    return mirbase_output


if __name__ == "__main__":
    parser = ArgumentParser(description="Run different tools in simulated fasta file")
    parser.add_argument("--fasta", required=True, help="short reads")
    parser.add_argument("--mirbase", required=True, help="bed file with mirbase annotation")
    parser.add_argument("--star", help="star index")
    parser.add_argument("--bowtie2", help="bowtie2 index")
    parser.add_argument("--hisat", help="hisat index")
    args = parser.parse_args()

    outputs = {}
    if args.star:
        print "doing STAR"
        outputs.update({"star": _star(full(args.fasta), full(args.star), full(args.mirbase))})
    if args.bowtie2:
        print "doing bowtie2"
        outputs.update({"bowtie2": _bowtie2(full(args.fasta), full(args.bowtie2), full(args.mirbase))})
    if args.hisat:
        print "doing hisat"
        outputs.update({"hisat": _hisat(full(args.fasta), full(args.hisat), full(args.mirbase))})

    os.remove("summary.tsv") if os.path.exists("summary.tsv") else None
    with open("summary.tsv", 'w') as out_handle:
        out_handle.write(res.H + "\n")
    for tool, stat in outputs.items():
        stat_file = _stats(stat, args.fasta, tool)
        do.run("cat %s >> summary.tsv" % stat_file, "merging %s" % tool)
