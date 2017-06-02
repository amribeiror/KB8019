import glob, re, math, sys, time
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity as cs_angle

print ("Executing script {}...".format(sys.argv[0]))
starttime = time.time()
seq_no = [] # eachrow
gc_content = []
nt_count = []
dnt_count = []

def gc_distance(gc_genome1, gc_genome2):
	distance = math.sqrt((gc_genome1-gc_genome2)**2)
	return distance

def nucleotide_distance(NT_genome1, NT_genome2):
	"""A,G,T,C --> 0,1,2,3"""
	distance = round(math.sqrt((NT_genome1[0]-NT_genome2[0])**2 +
	                    (NT_genome1[1]-NT_genome2[1])**2 +
	                    (NT_genome1[2]-NT_genome2[2])**2 +
	                    (NT_genome1[3]-NT_genome2[3])**2),4)
	return distance

def dinucleotide_distance(DNT_genome1, DNT_genome2):
	summation = [(a-b)**2 for a,b in zip(DNT_genome1, DNT_genome2)]
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

dnt_names = ['CC', 'GG', 'AA', 'TT',
	        'CG', 'CA', 'CT', 'GC',
	        'GA', 'GT', 'AC', 'AG',
	        'AT', 'TC', 'TG', 'TA']
agtc = ["A","G","T","C"]


c = 0
for filename in glob.glob("??.stats"):
	if not "36" in filename:
		c += 1
		seq_no.append((re.match(".*([0-9][0-9]+)", filename)).group(1))
		with open(filename, 'r') as f:
			f = f.read().splitlines()
			gc_content.append(float((f[1].replace("%", "").split())[2]))
			nt_range = []
			for ea in f[3:7]:
			    nt_range.append(float((ea.split())[1]))
			nt_count.append(nt_range)
			dnt_range = []
			for ea in f[8:24]:
			    dnt_range.append(float((ea.split())[1]))
			dnt_count.append(dnt_range)

gc_matrix = []
for i in gc_content:
	for j in gc_content:
	    gc_matrix.append(round(math.sqrt((i-j)**2),4))
gc_matrix = np.array(gc_matrix).reshape((c,c))
matrix_files(gc_matrix, "gc_matrix")

nt_matrix = []
for i in nt_count:
	for j in nt_count:
	    nt_matrix.append(nucleotide_distance(i, j))
nt_matrix = np.array(nt_matrix).reshape((c,c))
matrix_files(nt_matrix, "nt_matrix")

dnt_matrix = []
for i in dnt_count:
	for j in dnt_count:
	    dnt_matrix.append(dinucleotide_distance(i, j))
dnt_matrix = np.array(dnt_matrix).reshape((c,c))
matrix_files(dnt_matrix, "dnt_matrix")

nt_angle_matrix = []
for i in nt_count:
	i = np.array([i])
	for j in nt_count:
	    j = np.array([j])
	    nt_angle_matrix.append(cs_angle(i, j))
nt_angle_matrix = np.array(nt_angle_matrix).reshape((c,c))
nt_angle_matrix = np.around(nt_angle_matrix, decimals=4)
matrix_files(nt_angle_matrix, "nt_angle_matrix")

dnt_angle_matrix = []
for i in dnt_count:
	i = np.array([i])
	for j in dnt_count:
	    j = np.array([j])
	    dnt_angle_matrix.append(cs_angle(i, j))
dnt_angle_matrix = np.array(dnt_angle_matrix).reshape((c,c))
dnt_angle_matrix = np.around(dnt_angle_matrix, decimals=4)
matrix_files(dnt_angle_matrix, "dnt_angle_matrix")


print ("Finished execution of script {} in {} second(s).".format(
	    sys.argv[0], round(time.time()-starttime,2)))
