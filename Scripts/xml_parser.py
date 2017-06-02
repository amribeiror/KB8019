import sys, re
from Bio.Blast import NCBIXML

handle = open(sys.argv[1])
print (sys.argv[1])
records = NCBIXML.parse(handle)

counthit, count_nohit = 0, 0
id_eval_more0 = []
id_eval_less0 = []
less0, more0 = 0, 0

for eachrecord in records:
    if len(eachrecord.alignments) > 0:
        counthit += 1
        querydef = eachrecord.query
        queryid = eachrecord.query_id
        besthit = float(eachrecord.alignments[0].hsps[0].expect)
        if besthit <= 0.001:
            id_eval_less0.append((queryid,besthit))
            less0 += 1
        elif besthit > 0.001:
            id_eval_more0.append((queryid,besthit))
            more0 += 1
    else:
        count_nohit += 1


print ("Total number of queries: %i" %(count_nohit+counthit))
print ("Hit found: %i\t No hits: %i" %(counthit, count_nohit))
print ("<= 0.001: %i hits\t >0.001: %i hits\n" %(less0, more0))
