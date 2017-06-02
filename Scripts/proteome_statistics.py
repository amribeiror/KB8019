import sys, time
from Bio import SeqIO
from itertools import product
from collections import OrderedDict as odict

proteome = sys.argv[1]

aa = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S",
            "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
daa = ['GG', 'GA', 'GL', 'GM', 'GF', 'GW', 'GK', 'GQ', 'GE', 'GS',
        'GP', 'GV', 'GI', 'GC', 'GY', 'GH', 'GR', 'GN', 'GD', 'GT',
        'AG', 'AA', 'AL', 'AM', 'AF', 'AW', 'AK', 'AQ', 'AE', 'AS',
        'AP', 'AV', 'AI', 'AC', 'AY', 'AH', 'AR', 'AN', 'AD', 'AT',
        'LG', 'LA', 'LL', 'LM', 'LF', 'LW', 'LK', 'LQ', 'LE', 'LS',
        'LP', 'LV', 'LI', 'LC', 'LY', 'LH', 'LR', 'LN', 'LD', 'LT',
        'MG', 'MA', 'ML', 'MM', 'MF', 'MW', 'MK', 'MQ', 'ME', 'MS',
        'MP', 'MV', 'MI', 'MC', 'MY', 'MH', 'MR', 'MN', 'MD', 'MT',
        'FG', 'FA', 'FL', 'FM', 'FF', 'FW', 'FK', 'FQ', 'FE', 'FS',
        'FP', 'FV', 'FI', 'FC', 'FY', 'FH', 'FR', 'FN', 'FD', 'FT',
        'WG', 'WA', 'WL', 'WM', 'WF', 'WW', 'WK', 'WQ', 'WE', 'WS',
        'WP', 'WV', 'WI', 'WC', 'WY', 'WH', 'WR', 'WN', 'WD', 'WT',
        'KG', 'KA', 'KL', 'KM', 'KF', 'KW', 'KK', 'KQ', 'KE', 'KS',
        'KP', 'KV', 'KI', 'KC', 'KY', 'KH', 'KR', 'KN', 'KD', 'KT',
        'QG', 'QA', 'QL', 'QM', 'QF', 'QW', 'QK', 'QQ', 'QE', 'QS',
        'QP', 'QV', 'QI', 'QC', 'QY', 'QH', 'QR', 'QN', 'QD', 'QT',
        'EG', 'EA', 'EL', 'EM', 'EF', 'EW', 'EK', 'EQ', 'EE', 'ES',
        'EP', 'EV', 'EI', 'EC', 'EY', 'EH', 'ER', 'EN', 'ED', 'ET',
        'SG', 'SA', 'SL', 'SM', 'SF', 'SW', 'SK', 'SQ', 'SE', 'SS',
        'SP', 'SV', 'SI', 'SC', 'SY', 'SH', 'SR', 'SN', 'SD', 'ST',
        'PG', 'PA', 'PL', 'PM', 'PF', 'PW', 'PK', 'PQ', 'PE', 'PS',
        'PP', 'PV', 'PI', 'PC', 'PY', 'PH', 'PR', 'PN', 'PD', 'PT',
        'VG', 'VA', 'VL', 'VM', 'VF', 'VW', 'VK', 'VQ', 'VE', 'VS',
        'VP', 'VV', 'VI', 'VC', 'VY', 'VH', 'VR', 'VN', 'VD', 'VT',
        'IG', 'IA', 'IL', 'IM', 'IF', 'IW', 'IK', 'IQ', 'IE', 'IS',
        'IP', 'IV', 'II', 'IC', 'IY', 'IH', 'IR', 'IN', 'ID', 'IT',
        'CG', 'CA', 'CL', 'CM', 'CF', 'CW', 'CK', 'CQ', 'CE', 'CS',
        'CP', 'CV', 'CI', 'CC', 'CY', 'CH', 'CR', 'CN', 'CD', 'CT',
        'YG', 'YA', 'YL', 'YM', 'YF', 'YW', 'YK', 'YQ', 'YE', 'YS',
        'YP', 'YV', 'YI', 'YC', 'YY', 'YH', 'YR', 'YN', 'YD', 'YT',
        'HG', 'HA', 'HL', 'HM', 'HF', 'HW', 'HK', 'HQ', 'HE', 'HS',
        'HP', 'HV', 'HI', 'HC', 'HY', 'HH', 'HR', 'HN', 'HD', 'HT',
        'RG', 'RA', 'RL', 'RM', 'RF', 'RW', 'RK', 'RQ', 'RE', 'RS',
        'RP', 'RV', 'RI', 'RC', 'RY', 'RH', 'RR', 'RN', 'RD', 'RT',
        'NG', 'NA', 'NL', 'NM', 'NF', 'NW', 'NK', 'NQ', 'NE', 'NS',
        'NP', 'NV', 'NI', 'NC', 'NY', 'NH', 'NR', 'NN', 'ND', 'NT',
        'DG', 'DA', 'DL', 'DM', 'DF', 'DW', 'DK', 'DQ', 'DE', 'DS',
        'DP', 'DV', 'DI', 'DC', 'DY', 'DH', 'DR', 'DN', 'DD', 'DT',
        'TG', 'TA', 'TL', 'TM', 'TF', 'TW', 'TK', 'TQ', 'TE', 'TS',
        'TP', 'TV', 'TI', 'TC', 'TY', 'TH', 'TR', 'TN', 'TD', 'TT']

def proteome_stat(metagene):
    G, A, L, M, F = 0,0,0,0,0
    W, K, Q, E, S = 0,0,0,0,0
    P, V, I, C, Y = 0,0,0,0,0
    H, R, N, D, T = 0,0,0,0,0
    length = 0
    for i in metagene:
        length += 1
        if i == "G": G += 1
        elif i == "A": A+= 1
        elif i == "L": L += 1
        elif i == "M": M += 1
        elif i == "F": F += 1
        elif i == "W": W += 1
        elif i == "K": K += 1
        elif i == "Q": Q += 1
        elif i == "E": E += 1
        elif i == "S": S += 1
        elif i == "P": P += 1
        elif i == "V": V += 1
        elif i == "I": I += 1
        elif i == "C": C += 1
        elif i == "Y": Y += 1
        elif i == "H": H += 1
        elif i == "R": R += 1
        elif i == "N": N += 1
        elif i == "D": D += 1
        elif i == "T": T += 1

    daalist1 = [(metagene[d:d+2]) for d in range(0,len(metagene),2)
                if len(metagene[d:d+2]) == 2]
    daacount1 = odict([[x, daalist1.count(x)] for x in daa])

    daalist2 = [(metagene[d:d+2]) for d in range(1,len(metagene),2)
                if len(metagene[d:d+2]) == 2]
    daacount2 = odict([[x, daalist2.count(x)] for x in daa])

    daacount = []
    for i in daa:
        addition = (daacount1[i]+daacount2[i])/(length-1)
        daacount.append(addition)

    aastat = [G, A, L, M, F, W, K, Q, E, S,
                P, V, I, C, Y, H, R, N, D, T]
    aastat = [float(i/length) for i in aastat]
    thestats = [aastat, daacount]
    return thestats

with open(proteome, 'r') as pt:
    print ("Parsing file {}...".format(proteome))
    metagene = list()
    records = list(SeqIO.parse(pt, 'fasta'))
    for i in range(len(records)):
        metagene.extend(records[i].seq._data)
    metagene = ''.join(metagene)
    metagene = metagene.upper()
    stats = proteome_stat(metagene)
    with open(proteome+".stats", 'w') as ps:
        ps.write("AMINO ACID STATISTICS FOR: {}\n".format(proteome))
        ps.write("Amino Acid Frequency:\n")
        for i in range(len(aa)):
            ps.write("{}: {}\n".format(aa[i], stats[0][i]))
        ps.write("Checksum: {}\n".format(sum(stats[1])))
        ps.write("Diamino Acid Frequency:\n")
        for i in range(len(daa)):
            ps.write("{}: {}\n".format(daa[i], stats[1][i]))
        ps.write("Checksum: {}".format(sum(stats[0])))

    print ("Finished parsing file {}.".format(proteome))
