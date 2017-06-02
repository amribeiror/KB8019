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

def gc_content(sequence):

    gc_content = 0
    tot_nucleotides = 0

    seq = sequence.upper()

    for nucleotide in seq:

        if nucleotide != 'N':
            tot_nucleotides += 1

            if nucleotide == 'G':
            	gc_content += 1

            if nucleotide == 'C':
            	gc_content += 1

    fraction = float(gc_content/tot_nucleotides)*100
    # fraction = float(gc_content / tot_nucleotides)

    #return ('GC content: %.2f, total nucleotides: %.2f, fraction: %.2f %%' % (gc_content, tot_nucleotides, fraction))
    return (fraction)

def count_nucleotide(sequence):
    sequence = sequence.upper()
    A, G, T, C = 0, 0, 0, 0
    for nt in sequence:
        if nt == "A": A+= 1
        elif nt == "G": G+= 1
        elif nt == "T": T+= 1
        elif nt == "C": C+= 1
    nt_count = [A,G,T,C]
    nt_count = [q/len(sequence) for q in nt_count]
    return (nt_count)

def dinucleotides(sequence):
    print ("\tCounting dinucleotides...")
    cc = 0
    gg = 0
    aa = 0
    tt = 0
    cg = 0
    ca = 0
    ct = 0
    gc = 0
    ga = 0
    gt = 0
    ac = 0
    ag = 0
    at = 0
    tc = 0
    tg = 0
    ta = 0

    seq = sequence.upper()

    for n in range(len(seq)):
        print ("{}/{}".format(n, len(seq)))
        if (n) < len(seq)-1:
            if seq[n:n+1] == 'CC': cc += 1
            elif seq[n:n+1] == 'GG': gg += 1
            elif seq[n:n+1] == 'AA': aa += 1
            elif seq[n:n+1] == 'TT': tt += 1

            elif seq[n:n+1] == 'CG': cg += 1
            elif seq[n:n+1] == 'CA': ca += 1
            elif seq[n:n+1] == 'CT': ct += 1
            elif seq[n:n+1] == 'GC': gc += 1

            elif seq[n:n+1] == 'GA': ga += 1
            elif seq[n:n+1] == 'GT': gt += 1
            elif seq[n:n+1] == 'AC': ac += 1
            elif seq[n:n+1] == 'AG': ag += 1

            elif seq[n:n+1] == 'AT': at += 1
            elif seq[n:n+1] == 'TC': tc += 1
            elif seq[n:n+1] == 'TG': tg += 1
            elif seq[n:n+1] == 'TA': ta += 1

    dnt_stats = [cc, gg, aa, tt, cg, ca, ct, gc,
                    ga, gt, ac, ag, at, tc, tg, ta]
    # tot = (cc+gg+aa+tt+cg+ca+ct+gc+ga+gt+ac+ag+at+tc+tg+ta)
    tot = sum(dnt_stats)
    dnt_stats = [a/tot for a in dnt_stats]
    # gc = (cc+gg+gc+cg)
    # dnt_stats = [cc, gg, aa, tt, cg, ca, ct, gc,
    #                 ga, gt, ac, ag, at, tc, tg, ta]
    # print (cc, gg, aa, tt, cg, ca, ct, gc, ga, gt, ac, ag, at, tc, tg, ta)
    # fraction = gc / tot
    # return fraction
    return dnt_stats



dnt_names = ['CC', 'GG', 'AA', 'TT',
            'CG', 'CA', 'CT', 'GC',
            'GA', 'GT', 'AC', 'AG',
            'AT', 'TC', 'TG', 'TA']
agtc = ["A","G","T","C"]
# for fastafile in glob.glob("../raw_sequences/*.fa"):
#     records = list(SeqIO.parse(fastafile, "fasta"))
#     #print (type(records[0].seq))
#     seqId = records[0].id
#     seq = records[0].seq
#     #print (dinucleotides(seq[::-1]))
#     #print (dinucleotides(seq))
#     gc = gc_content(seq)
#     gc_pairs = dinucleotides(seq)
#
#     with open('gc_content.txt', 'a') as f:
#         f.write('%s\nGC content: %.2f %%\n' % (seqId, gc))

for fastafile in glob.glob("../raw_sequences/??.fa"):
    print ("Parsing {}...".format(fastafile))
    starttime = time.time()
    with open(fastafile, 'r') as fa:
        fa = SeqIO.read(fa, "fasta")
        sequence = str()
        for i in fa.seq._data:
            if i in agtc:
                sequence += i
        stats = statistics(sequence)
        fastafile = fastafile.replace(".fa", "")
        with open(fastafile+".stats", 'w') as fs:
            fs.write("NUCLEOTIDE STATISTICS FOR: {}\n".format(fastafile))
            fs.write("GC Content: {0: >7}%\n".format(str(stats[0][0])))
            fs.write("Nucleotide Frequency:\n")
            for i in range(len(agtc)):
                fs.write("{}: {}\n".format(agtc[i], stats[1][i]))
            fs.write("Dinucleotide Frequency:\n")
            for c,d in zip(dnt_names, stats[2]):
                fs.write("{}: {}\n".format(c,d))
    print ("Finished processing in {} s.".format(round(time.time()-starttime),3))
