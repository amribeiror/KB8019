import glob, re, math, sys, time
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity as cs_angle

seq_no = []
aa_count = []
daa_count = []

def aa_distance(proteome1, proteome2):
    summation = [(a-b)**2 for a,b in zip(proteome1, proteome2)]
    distance = round(math.sqrt(sum(summation)),4)
    return distance

def daa_distance(proteome1, proteome2):
    summation = [(a-b)**2 for a,b in zip(proteome1, proteome2)]
    distance = round(math.sqrt(sum(summation)),4)
    return distance

def matrix_files(thematrix, filename):
    location = "distance_matrices/"
    with open(location+filename+'.txt', 'w') as st:
        for numbers in seq_no:
            st.write("{0: >8} ".format(numbers))
        st.write("\n")
        for eachline in thematrix:
            for eachel in eachline:
                st.write("{0: >8} ".format(eachel))
            st.write("\n")
    return True

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


for filename in glob.glob("??.fa.noN.pfa.stats"):
    seq_no.append(re.match(".*([0-9][0-9]+)", filename).group(1))
    with open(filename, 'r') as f:
        f = f.read().splitlines()
        aa_range = []
        for ea in f[2:22]:
            aa_range.append(float((ea.split())[1]))
        aa_count.append(aa_range)
        daa_range = []
        for ea in f[24:424]:
            daa_range.append(float((ea.split())[1]))
        daa_count.append(daa_range)

aa_matrix = []
for i in aa_count:
    for j in aa_count:
        aa_matrix.append(aa_distance(i, j))
aa_matrix = np.array(aa_matrix).reshape((c,c))
matrix_files(aa_matrix, "aa_matrix")

daa_matrix = []
for i in daa_count:
    for j in daa_count:
        daa_matrix.append(daa_distance(i, j))
daa_matrix = np.array(daa_matrix).reshape((c,c))
matrix_files(daa_matrix, "daa_matrix")

aa_angle_matrix = []
for i in aa_count:
    i = np.array([i])
    for j in aa_count:
        j = np.array([j])
        aa_angle_matrix.append(cs_angle(i, j))
aa_angle_matrix = np.array(aa_angle_matrix).reshape((c,c))
aa_angle_matrix = np.around(aa_angle_matrix, decimals=4)
matrix_files(aa_angle_matrix, "aa_angle_matrix")

daa_angle_matrix = []
for i in daa_count:
    i = np.array([i])
    for j in daa_count:
        j = np.array([j])
        daa_angle_matrix.append(cs_angle(i, j))
daa_angle_matrix = np.array(daa_angle_matrix).reshape((c,c))
daa_angle_matrix = np.around(daa_angle_matrix, decimals=4)
matrix_files(daa_angle_matrix, "daa_angle_matrix")

print ("Finished.")
