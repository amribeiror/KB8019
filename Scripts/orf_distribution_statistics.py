import glob
import matplotlib.pyplot as plt
from statistics import median, mean

colours = ['r','c','b','m','k']
number = []
orfbin = []
orf_stat = []
bins = [x*1000 for x in range(1,16)]
for filename in glob.glob("??.fa.translate.orf_coords"):
    no = filename.replace(".fa.translate.orf_coords", "")
    number.append(no)
    orf_sub_bin = []
    with open(filename, 'r') as fr:
        fr = fr.read().splitlines()
        for line in fr:
            if not ">" in line:
                line = line.split()
                orf_sub_bin.append(int(line[-2]))
    orfbin.append(orf_sub_bin)
    orf_sub_stat = [max(orf_sub_bin), min(orf_sub_bin),
                    median(orf_sub_bin), round(mean(orf_sub_bin),0)]
    orf_stat.append(orf_sub_stat)

plt.figure()
for i in range(len(number)):
    plt.hist(orfbin[i], histtype='step', color=colours[i], label=number[i], bins=bins)
plt.xlabel("Size (bp)")
plt.ylabel("Count")
plt.title("Distribution of Predicted ORF Size")
plt.legend(loc="upper right")
# plt.show()
plt.savefig("ORF_Stat_fig.png")
plt.savefig("ORF_Stat_fig.svg")

print ("Finished.")
