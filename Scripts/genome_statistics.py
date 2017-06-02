import glob
from Bio.Blast import NCBIXML
from Bio import SeqIO
import numpy, math, time

def statistics(sequence):
    sequence = sequence.upper()
    gc_count = 0
    a, g, t, c = 0, 0, 0, 0
    cc, gg, aa, tt = 0, 0, 0, 0
    cg, ca, ct = 0, 0, 0
    gc, ga, gt = 0, 0, 0
    ac, ag, at = 0, 0, 0
    tc, tg, ta = 0, 0, 0
    length = 0
    for nt in range(len(sequence)):
        if sequence[nt] != 'N':
            length += 1
            if sequence[nt] == 'C':
                c += 1
                gc_count += 1
            elif sequence[nt] == 'G':
                g += 1
                gc_count += 1
            elif sequence[nt] == 'A':
                a += 1
            elif sequence[nt] == 'T':
                t += 1

            if nt < len(sequence)-1:
                if sequence[nt:nt+2] == 'CC': cc += 1
                elif sequence[nt:nt+2] == 'GG': gg += 1
                elif sequence[nt:nt+2] == 'AA': aa += 1
                elif sequence[nt:nt+2] == 'TT': tt += 1

                elif sequence[nt:nt+2] == 'CG': cg += 1
                elif sequence[nt:nt+2] == 'CA': ca += 1
                elif sequence[nt:nt+2] == 'CT': ct += 1
                elif sequence[nt:nt+2] == 'GC': gc += 1

                elif sequence[nt:nt+2] == 'GA': ga += 1
                elif sequence[nt:nt+2] == 'GT': gt += 1
                elif sequence[nt:nt+2] == 'AC': ac += 1
                elif sequence[nt:nt+2] == 'AG': ag += 1

                elif sequence[nt:nt+2] == 'AT': at += 1
                elif sequence[nt:nt+2] == 'TC': tc += 1
                elif sequence[nt:nt+2] == 'TG': tg += 1
                elif sequence[nt:nt+2] == 'TA': ta += 1
    gc_count = [gc_count/length]
    nt_stats = [a/length, g/length, c/length, t/length]
    dnt_stats = [cc, gg, aa, tt, cg, ca, ct, gc,
                    ga, gt, ac, ag, at, tc, tg, ta]
    dnt_stats = [float(a/(length-1)) for a in dnt_stats]
    the_stats = [gc_count, nt_stats, dnt_stats]
    return the_stats

dnt_names = ['CC', 'GG', 'AA', 'TT',
            'CG', 'CA', 'CT', 'GC',
            'GA', 'GT', 'AC', 'AG',
            'AT', 'TC', 'TG', 'TA']
agtc = ["A","G","T","C"]

with open(sys.argv[1], 'r') as fa:
    print ("Parsing {}...".format(sys.arg[1]))
    starttime = time.time()
    fa = SeqIO.read(fa, "fasta")
    sequence = str()
    for i in fa.seq._data:
        if i in agtc:
            sequence += i
    stats = statistics(sequence)
    with open(sys.argv[1]+".stats", 'w') as fs:
        fs.write("NUCLEOTIDE STATISTICS FOR: {}\n".format(fastafile))
        fs.write("GC Content: {0: >7}%\n".format(str(stats[0][0])))
        fs.write("Nucleotide Frequency:\n")
        for i in range(len(agtc)):
            fs.write("{}: {}\n".format(agtc[i], stats[1][i]))
        fs.write("Dinucleotide Frequency:\n")
        for c,d in zip(dnt_names, stats[2]):
            fs.write("{}: {}\n".format(c,d))
    print ("Finished processing in {} s.".format(round(time.time()-starttime),3))
