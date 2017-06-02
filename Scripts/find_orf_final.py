from Bio import SeqIO
import sys, re, time
from operator import itemgetter
import itertools


################################################################################
##### PARAMETERS ###############################################################
startcodon = ['ATG', 'GTG', 'TTG'] #GLIMMER also uses gtg,ttg
stopcodon = ['TGA', 'TAA', 'TAG']
purine = ["A", "G"]
consensusSeq10 = list("TATAAT")
consensusSeq35 = list("TTGACA")
Shine_Dalgarno = list("AGGAGG")
Kozak = list("GCCRCC")

################################################################################
##### SETTINGS #################################################################
transcription = 0
translation = 1
################################################################################

def revcomplement(seq):
    seq = seq[::-1]
    complement = str()
    for AA in seq:
        if AA == 'A': complement += 'T'
        elif AA == 'T': complement += 'A'
        elif AA == 'G': complement += 'C'
        elif AA == 'C': complement += 'G'
    return complement

def forward_prokaryote(sequence):
    theORFs = []
    ### FORWARD STRAND #########################################################
    for rf in [0,1,2]:
        print ("Checking reading frame +%s..." %(rf+1))
        orf_start = []
        orf_end = []
        frame = []
        startmode = True
        for nt in range(rf, len(sequence), 3):
            if nt+3 < len(sequence):
                triplets = sequence[nt:nt+3]
                frame.append("+%s" %str(rf+1))
                if triplets in startcodon and startmode == True:
                    orf_start.append(nt)
                    startmode = False
                elif triplets in stopcodon and startmode == False:
                    orf_end.append(nt+3)
                    startmode = True
        startmode = False
        ### Removes if there are more start codons than there are stop ones.
        if len(orf_start) > len(orf_end):
            diff = len(orf_end)-len(orf_start)
            del (orf_start[diff:])
        orf = list(zip(orf_start, orf_end))
        ### Collects the nucleotide positions for spacers - between ORFs
        inbetween = []
        inbetween.insert(0, (0, orf_start[0]))
        for a,b in zip(orf_start[1:], orf_end):
            inbetween.append((b,a))
        # inbetween.insert(len(inbetween), (orf_end[-1], len(seq)))
        ########################################################################
        transcript = []
        for z, sep in enumerate(inbetween):
            tt = []
            ### Only transcription (-35 and Pribnow boxes)
            checkboxes = sequence[sep[0]:(sep[1]-16)]
            if len(checkboxes) >= 45:
                box35_final_score = 0
                temp_n = 0
                for b in range(len(checkboxes)):
                    box35_subscore = 0
                    box35_seq = checkboxes[b:b+6]
                    if len(box35_seq) == 6:
                        for nt in range(len(box35_seq)):
                            if box35_seq[nt] == consensusSeq35[nt]:
                                box35_subscore += 1
                        if (box35_subscore > box35_final_score and
                            box35_subscore >= 3):
                            box35_final_score = box35_subscore
                            temp_n = b
                pribbnow_region = checkboxes[temp_n+6:temp_n+6+25]
                if len(pribbnow_region) >= 6:
                    pribnow_final_score = 0
                    for c in range(len(pribbnow_region)):
                        pribnow_subscore = 0
                        pribnow_check_box = pribbnow_region[c:c+6]
                        for i in range(len(pribnow_check_box)):
                            if pribnow_check_box[i] == consensusSeq10[i]:
                                pribnow_subscore += 1
                        if (pribnow_subscore > pribnow_final_score and
                            pribnow_subscore >= 4):
                            pribnow_final_score = pribnow_subscore
                if box35_final_score >= 3 and pribnow_final_score >= 4:
                    if (transcription == True and
                        translation == False):
                        transcript.append(z)
                    elif translation == True and translation == True:
                        tt.append(z)

            ### Only translation (Shine_Dalgarno)
            sd_final_score = 0
            checkseq = sequence[(-16+sep[1]):sep[1]]
            for nt in range(0,len(checkseq)):
                sd_seq = checkseq[nt:nt+6]
                sd_subscore = 0
                if len(sd_seq) == 6:
                    for i in range(len(sd_seq)):
                        if sd_seq[i] == Shine_Dalgarno[i]:
                            sd_subscore += 1
                if sd_subscore > sd_final_score and sd_subscore >= 4:
                    sd_final_score = sd_subscore
            if sd_final_score >= 4:
                if translation == True and transcription == False:
                    transcript.append(z)
                elif translation == True and transcription == True:
                    tt.append(z)

            if len(tt) == 2:
                transcript.append(z)

        theorf = [orf[q] for q in transcript
                    if orf[q][1]-orf[q][0] >= 90]
        orf_start = [theorf[a][0] for a in range(len(theorf))]
        orf_end = [theorf[a][1] for a in range(len(theorf))]
        orflen = [theorf[i][1]-theorf[i][0] for i in range(len(theorf))]
        ziporf = set(list(zip(orf_start, orf_end, orflen, frame)))
        theORFs.extend(ziporf)
    return theORFs

def reverse_prokaryote(sequence):
    theORFs = []
    ### REVERSE STRAND #########################################################
    for rf in [0,1,2]:
        print ("Checking reading frame -%s..." %(rf+1))
        orf_start = []
        orf_end = []
        frame = []
        startmode = True
        for nt in range(rf, len(sequence), 3):
            if nt+3 < len(sequence):
                triplets = sequence[nt:nt+3]
                frame.append("+%s" %str(rf+1))
                if triplets in startcodon and startmode == True:
                    orf_start.append(nt)
                    startmode = False
                elif triplets in stopcodon and startmode == False:
                    orf_end.append(nt+3)
                    startmode = True
        startmode = False
        ### Removes if there are more start codons than there are stop ones.
        if len(orf_start) > len(orf_end):
            diff = len(orf_end)-len(orf_start)
            del (orf_start[diff:])
        orf = list(zip(orf_start, orf_end))
        ### Collects the nucleotide positions for spacers - between ORFs
        inbetween = []
        inbetween.insert(0, (0, orf_start[0]))
        for a,b in zip(orf_start[1:], orf_end):
            inbetween.append((b,a))
        # inbetween.insert(len(inbetween), (orf_end[-1], len(seq)))

        ########################################################################
        transcript = []
        for z, sep in enumerate(inbetween):
            tt = []
            ### Only transcription (-35 and Pribnow boxes)
            checkboxes = sequence[sep[0]:(sep[1]-16)]
            if len(checkboxes) >= 45:
                box35_final_score = 0
                temp_n = 0
                for b in range(len(checkboxes)):
                    box35_subscore = 0
                    box35_seq = checkboxes[b:b+6]
                    if len(box35_seq) == 6:
                        for nt in range(len(box35_seq)):
                            if box35_seq[nt] == consensusSeq35[nt]:
                                box35_subscore += 1
                        if (box35_subscore > box35_final_score and
                            box35_subscore >= 3):
                            box35_final_score = box35_subscore
                            temp_n = b
                pribbnow_region = checkboxes[temp_n+6:temp_n+6+25]
                if len(pribbnow_region) >= 6:
                    pribnow_final_score = 0
                    for c in range(len(pribbnow_region)):
                        pribnow_subscore = 0
                        pribnow_check_box = pribbnow_region[c:c+6]
                        for i in range(len(pribnow_check_box)):
                            if pribnow_check_box[i] == consensusSeq10[i]:
                                pribnow_subscore += 1
                        if (pribnow_subscore > pribnow_final_score and
                            pribnow_subscore >= 4):
                            pribnow_final_score = pribnow_subscore
                if box35_final_score >= 3 and pribnow_final_score >= 4:
                    if (transcription == True and
                        translation == False):
                        transcript.append(z)
                    elif translation == True and translation == True:
                        tt.append(z)

            ### Only translation (Shine_Dalgarno)
            sd_final_score = 0
            checkseq = sequence[(-16+sep[1]):sep[1]]
            for nt in range(0,len(checkseq)):
                sd_seq = checkseq[nt:nt+6]
                sd_subscore = 0
                if len(sd_seq) == 6:
                    for i in range(len(sd_seq)):
                        if sd_seq[i] == Shine_Dalgarno[i]:
                            sd_subscore += 1
                if sd_subscore > sd_final_score and sd_subscore >= 4:
                    sd_final_score = sd_subscore
            if sd_final_score >= 4:
                if translation == True and transcription == False:
                    transcript.append(z)
                elif translation == True and transcription == True:
                    tt.append(z)

            if len(tt) == 2:
                transcript.append(z)

        theorf = [orf[q] for q in transcript
                    if orf[q][1]-orf[q][0] >= 90]
        orf_start = [theorf[a][0] for a in range(len(theorf))]
        orf_end = [theorf[a][1] for a in range(len(theorf))]
        orflen = [theorf[i][1]-theorf[i][0] for i in range(len(theorf))]
        ziporf = set(list(zip(orf_start, orf_end, orflen, frame)))
        theORFs.extend(ziporf)
    return theORFs

def forward_eukaryote(sequence):
    theORFs = []
    ### FORWARD STRAND #########################################################
    for rf in [0,1,2]:
        print ("Checking reading frame +%s..." %(rf+1))
        orf_start = []
        orf_end = []
        frame = []
        startmode = True
        for nt in range(rf, len(sequence), 3):
            if nt+3 < len(sequence):
                triplets = sequence[nt:nt+3]
                frame.append("+%s" %str(rf+1))
                if triplets == "ATG" and startmode == True:
                    orf_start.append(nt)
                    startmode = False
                elif triplets in stopcodon and startmode == False:
                    orf_end.append(nt+3)
                    startmode = True
        startmode = False
        ### Removes if there are more start codons than there are stop ones.
        if len(orf_start) > len(orf_end):
            diff = len(orf_end)-len(orf_start)
            del (orf_start[diff:])
        orf = list(zip(orf_start, orf_end))
        ### Collects the nucleotide positions for spacers - between ORFs
        inbetween = []
        inbetween.insert(0, (0, orf_start[0]))
        for a,b in zip(orf_start[1:], orf_end):
            inbetween.append((b, a))
        # inbetween.insert(len(inbetween), (orf_end[-1], len(seq)))
        ########################################################################
        transcript = []
        for z, sep in enumerate(inbetween):
            ### Only translation (Kozak)
            kozak_final_score = 0
            checkseq = sequence[(-3+sep[1]):sep[1]+4]
            # print (len(pattern))
            if len(checkseq) == 7:
                if checkseq[0] == "A":
                # if checkseq[-1] == "T":
                    if translation == 1 and transcription == 0:
                        transcript.append(z)

        theorf = [orf[q] for q in transcript
                        if q < len(orf)
                        and orf[q][1]-orf[q][0] >= 90]
        orf_start = [theorf[a][0] for a in range(len(theorf))]
        orf_end = [theorf[a][1] for a in range(len(theorf))]
        orflen = [theorf[i][1]-theorf[i][0] for i in range(len(theorf))]
        ziporf = set(list(zip(orf_start, orf_end, orflen, frame)))
        theORFs.extend(ziporf)
    return theORFs

def reverse_eukaryote(sequence):
    theORFs = []
    ### FORWARD STRAND #########################################################
    for rf in [0,1,2]:
        print ("Checking reading frame -%s..." %(rf+1))
        orf_start = []
        orf_end = []
        frame = []
        # pattern = r"([AG]{2}[AGTC]ATG{90,}[ATGC][T(GA|AA|AG)])"
        startmode = True
        for nt in range(rf, len(sequence), 3):
            if nt+3 < len(sequence):
                triplets = sequence[nt:nt+3]
                frame.append("-%s" %str(rf+1))
                if triplets == "ATG" and startmode == True:
                    orf_start.append(nt)
                    startmode = False
                elif triplets in stopcodon and startmode == False:
                    orf_end.append(nt+3)
                    startmode = True
        startmode = False
        ### Removes if there are more start codons than there are stop ones.
        if len(orf_start) > len(orf_end):
            diff = len(orf_end)-len(orf_start)
            del (orf_start[diff:])
        orf = list(zip(orf_start, orf_end))
        ### Collects the nucleotide positions for spacers - between ORFs
        inbetween = []
        inbetween.insert(0, (0, orf_start[0]))
        for a,b in zip(orf_start[1:], orf_end):
            inbetween.append((b, a))
        # inbetween.insert(len(inbetween), (orf_end[-1], len(seq)))
        ########################################################################
        transcript = []
        for z, sep in enumerate(inbetween):
            ### Only translation (Kozak)
            kozak_final_score = 0
            checkseq = sequence[(-3+sep[1]):sep[1]+4]
            # print (len(pattern))
            if len(checkseq) == 7:
                if checkseq[0] == "A":
                # if checkseq[-1] == "T":
                    if translation == 1 and transcription == 0:
                        transcript.append(z)

        theorf = [orf[q] for q in transcript
                        if q < len(orf)
                        and orf[q][1]-orf[q][0] >= 90]
        orf_start = [theorf[a][0] for a in range(len(theorf))]
        orf_end = [theorf[a][1] for a in range(len(theorf))]
        orflen = [theorf[i][1]-theorf[i][0] for i in range(len(theorf))]
        ziporf = set(list(zip(orf_start, orf_end, orflen, frame)))
        theORFs.extend(ziporf)
    return theORFs

with open(sys.argv[1], 'r') as seq:
    starttime = time.time()
    seq = seq.read().splitlines()
    sequence = seq[1].replace("N","")
    revseq = revcomplement(sequence)
    print ("Length of sequence %s: %r bp." %(sys.argv[1], len(sequence)))
    ORF = []

    while True:
        kingdom = int(input("Is the genome prokaryotic[1] or eukaryotic[2]:\t"))
        if kingdom == 1 or kingdom == 2:
            break

    if kingdom == 1:
        ORF.extend(forward_prokaryote(sequence))
        ORF.extend(reverse_prokaryote(revseq))
    if kingdom == 2:
        ORF.extend(forward_eukaryote(sequence))
        ORF.extend(reverse_eukaryote(revseq))
    print (len(ORF))

    if translation == True and transcription == False:
        extraname = ".translate"
    elif transcription == True and translation == False:
        extraname = ".trancription"
    elif transcription == True and translation == True:
        extraname = ".transcript_translate"

    with open(sys.argv[1]+extraname+".orf_coords", 'w') as oc:
        oc.write(">{} ORFcount:{}\n".format(sys.argv[1], str(len(ORF))))
        for z, ea in enumerate(sorted(ORF, key=itemgetter(0)), start=1):
            oc.write("orf{:05d}  ".format(z))
            oc.write("{: >7} {: >7} {: >7} {: >7}\n".format(*ea))

    print ("Finished processing %s in %s seconds.\n" %
        (sys.argv[1], round(time.time()-starttime,3)))
