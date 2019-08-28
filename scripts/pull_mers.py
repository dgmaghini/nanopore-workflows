import sys

inputfile = sys.argv[1]
outputname = sys.argv[2]
mercount = int(sys.argv[3])

# with open("mer_counts_dumps_101.fa", "r") as f:
# 	with open("fivemers.fa", "w") as output:
with open(inputfile, "r") as f:
	with open(outputname, "w") as output:
		line = f.readline().strip()
		counter = 0
		while line != "":
			if line[0] == ">" and int(line[1:]) >= mercount:
				output.write('>contig_' + line[1:] + "_" + str(counter) + "\n")
				counter += 1
				line = f.readline().strip()
				output.write(line + "\n")
				line = f.readline().strip()
			else:
				line = f.readline().strip()
