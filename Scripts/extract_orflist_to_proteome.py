from Bio import SeqIO
import sys

def revcomplement(seq):
    seq = seq[::-1]
    complement = str()
    for AA in seq:
        if AA == 'A': complement += 'T'
        elif AA == 'T': complement += 'A'
        elif AA == 'G': complement += 'C'
        elif AA == 'C': complement += 'G'
    return complement


codon = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
        "TAT":"Y", "TAC":"Y", "TAA":"", "TAG":"",
        "TGT":"C", "TGC":"C", "TGA":"", "TGG":"W",
        "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

with open(sys.argv[1], 'r') as seq:
    handle = SeqIO.read(seq, "fasta")
    sequence = (handle.seq._data).upper()
    revseq = revcomplement(sequence)

outfile = sys.argv[1]+".pfa"
with open(outfile, 'w') as out:
    with open(sys.argv[2], 'r') as fasta:
        fasta = fasta.read().splitlines()
        orfname = fasta[0].split()[0]
        for line in fasta[1:]:
            orfno, start, stop, length, frame = line.split()
            start = int(start)
            stop = int(stop)
            if frame[0] == '+':
                name = orfname+"_"+orfno
                frame = int(frame[1])
                if frame == 1:
                    theseq = sequence[0:]
                    writeseq = theseq[start:stop]
                    protseq = str()
                    # testcount = 0
                    for nt in range(0,len(writeseq),3):
                        protseq += codon[writeseq[nt:nt+3]]
                        # testcount += 1
                        # print (codon[writeseq[nt:nt+3]])
                        # if testcount == 10: sys.exit()
                    out.write("%s\n%s\n" % (name, protseq))
                elif frame == 2:
                    theseq = sequence[1:]
                    writeseq = theseq[start:stop]
                    protseq = str()
                    for nt in range(0,len(writeseq),3):
                        protseq += codon[writeseq[nt:nt+3]]
                    out.write("%s\n%s\n" % (name, protseq))
                elif frame == 3:
                    theseq = sequence[1:]
                    writeseq = theseq[start:stop]
                    protseq = str()
                    for nt in range(0,len(writeseq),3):
                        protseq += codon[writeseq[nt:nt+3]]
                    out.write("%s\n%s\n" % (name, protseq))

            elif frame[0] == '-':
                name = orfname+"_"+orfno+"_rev"
                frame = int(frame[1])
                if frame == 1:
                    theseq = revseq[0:]
                    writeseq = theseq[start:stop]
                    protseq = str()
                    for nt in range(0,len(writeseq),3):
                        protseq += codon[writeseq[nt:nt+3]]
                    out.write("%s\n%s\n" % (name, protseq))
                elif frame == 2:
                    theseq = revseq[1:]
                    writeseq = theseq[start:stop]
                    protseq = str()
                    for nt in range(0,len(writeseq),3):
                        protseq += codon[writeseq[nt:nt+3]]
                    out.write("%s\n%s\n" % (name, protseq))
                elif frame == 3:
                    theseq = revseq[1:]
                    writeseq = theseq[start:stop]
                    protseq = str()
                    for nt in range(0,len(writeseq),3):
                        protseq += codon[writeseq[nt:nt+3]]
                    out.write("%s\n%s\n" % (name, protseq))
print ("Output: %s" % outfile)
