with open("treponema_succ_bin114.5repeats101.fa", "r") as f:
	with open("treponema_kmer_intensities.tsv", "w") as output:
		for line in f: 
			if line[0] == ">":
				output.write(line.strip() + "\t" + "1" "\n")
			
