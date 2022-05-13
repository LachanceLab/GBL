#!/usr/bin/python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import argparse
import pandas as pd


def parse_genome(fasta_file):
    records = SeqIO.parse(fasta_file, 'fasta')
    return records


def write_cds_records(records, output_file):
    SeqIO.write(records, output_file, "fasta")


def translate_cds(records):
    translated_records = []
    for rec in records:
        if len(rec.seq) % 3 != 0:
            print("Length is not a multiple of 3. Skip {}".format(rec.id))
            continue
        rec_aa = rec.translate()
        rec_aa.id = rec.id
        rec_aa.description = rec.description.replace('mRNA', 'pep')
        translated_records.append(rec_aa)
    return translated_records


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Nucleotide CDS in fasta format')
    parser.add_argument('-o', '--output_file', help='Output file; amino acid sequences of CDS')
    args = parser.parse_args()
    cds_records = parse_genome(args.input_file)
    translated_records = translate_cds(cds_records)
    write_cds_records(translated_records, args.output_file)


if __name__ == '__main__':
    main(sys.argv[1:])
