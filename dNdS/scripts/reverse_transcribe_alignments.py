#!/usr/bin/env python3
import sys
import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO


def read_alignment(input_file):
    """
    Read alignment file in phylip format
    @params input_file: str, input alignment file
    @return: list, MSA
    """
    alignment = list(AlignIO.read(input_file, 'phylip'))
    return alignment


def parse_nucleotides(input_file):
    """
    Parse fasta file
        @params input_file: str, input alignment file
    @return: dict, nucleotide sequences
    """
    nt_seqs = SeqIO.parse(input_file, 'fasta')
    nt_seqs_dict = {}
    for seq in nt_seqs:
        # phylip format only shows first 10 symbols --> only need the first 10 symbols for matching ids
        nt_seqs_dict[seq.id[:10]] = seq
    return nt_seqs_dict


def reverse_translate(alignment, nt_seqs, species, output):
    """
    Reverse translate multiple sequence alignments of proteins
    @params alignment: list, MSA of protein sequences
    @params nt_seqs: dict, dictionary of corresponding nucleotide sequences
    @params species: list, of species to use in order of alignment
    @params output: str, output filename
    """
    aligned_nt_seqs = []
    # iterate over protein alignments
    for i, prot in enumerate(alignment):
        # extract corresponding nucleotide sequence
        nt = nt_seqs[prot.id]
        a = 0
        n = 0
        aligned_nt = ''
        while a < len(prot):
            # gap
            if prot.seq[a] == '-':
                aligned_nt += '---'
            # get current codon
            else:
                # assert it translates to the right aa
                assert nt.seq[n * 3: n * 3 + 3].translate() == prot[a], 'Wrong translation'
                aligned_nt += nt.seq[n * 3: n * 3 + 3]
                n += 1
            a += 1
        # create SeqRecords
        aligned_nt = SeqRecord(aligned_nt, id=nt.id, name=nt.name, description=nt.description)
        # asign ID based on species
        aligned_nt.id = species[i]
        # add to list
        aligned_nt_seqs.append(aligned_nt)
    # write nucleotide MSA in pseudo phylip format that is accepted by codeml
    with open(output, 'w') as out:
        out.write(f"{len(alignment)} {len(aligned_nt)}")
        
        for seq in aligned_nt_seqs:
            out.write('\n')
            out.write(f'{seq.id}  {seq.seq}')
    out.close()


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alignments', nargs='+', help='MSA files in phylip format', required=True)
    parser.add_argument('-c', '--cds', nargs='+', help='Corresponding files with nucleotide CDS in fasta format',
                        required=True)
    parser.add_argument('-s', '--species', help='List of species in order of alignments alignments', nargs='+',
                        required=True)
    parser.add_argument('-o', '--output', nargs='+', help='Outfiles in pseudo phylip format that is accept by codeml',
                        required=True)
    args = parser.parse_args()
    for alignments, cds, output in zip(args.alignments, args.cds, args.output):
        alignment = read_alignment(alignments)
        nt_seqs = parse_nucleotides(cds)
        reverse_translate(alignment, nt_seqs, args.species, output)


if __name__ == '__main__':
   main(sys.argv[1:])
