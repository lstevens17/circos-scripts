#!/usr/bin/env python3
from Bio import SeqIO
import argparse
import sys

def print_karyotype_file(fasta, species):
	with open(fasta, 'r') as fasta:
		with open(species + ".kt", 'w') as outfile:
			print("[+] Writing " + species + ".kt")
			chr2num = {}
			seqnum = 1
			for record in SeqIO.parse(fasta, "fasta"):
				if species == "sp1":
					outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr", "-", record.id, seqnum, "0", len(record.seq), seqnum))
					chr2num[record.id] = seqnum
					seqnum += 1
				else:
					outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("chr", "-", record.id, record.id, "0", len(record.seq), record.id))
	return chr2num


def parse_table(table_file):
	with open(table_file, 'r') as table:
		table_dict = {}
		for line in table:
			if not line.startswith("#"):
				cols = line.rstrip("\n").split()
				buscoID, status = cols[0], cols[1]
				if status == 'Complete':
					chr, start, stop = cols[2], int(cols[3]), int(cols[4])
					table_dict[buscoID] = [chr, start, stop]
	return table_dict

def print_links_file(t1_dict, t2_dict, chr2num):
	with open("links.txt", 'w') as outfile:
		print("[+] Writing links.txt")
		for buscoID, t1_list in t1_dict.items():
			t1_chr, t1_start, t1_stop = t1_list[0], t1_list[1], t1_list[2]
			num = str(chr2num[t1_chr.split(":")[0]])
			try:
				t2_list = t2_dict[buscoID]
				t2_chr, t2_start, t2_stop = t2_list[0], t2_list[1], t2_list[2]
				outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (t1_chr, t1_start, t1_stop, t2_chr, t2_start, t2_stop, "color=chr"+num))
			except KeyError:
				pass

if __name__ == "__main__":
	SCRIPT = "proccess_busco.py"
	# argument set up
	parser = argparse.ArgumentParser()
	parser.add_argument("-f1", "--fasta1", type=str, help = "genome FASTA for sp1", required=True)
	parser.add_argument("-f2", "--fasta2", type=str, help = "genome FASTA for sp2", required=True)
	parser.add_argument("-t1", "--table1", type=str, help = "full_table.tsv from BUSCO for sp1", required=True)
	parser.add_argument("-t2", "--table2", type=str, help = "full_table.tsv from BUSCO for sp1", required=True)
	args = parser.parse_args()
	f1_file = args.fasta1
	f2_file = args.fasta2
	t1_file = args.table1
	t2_file = args.table2
	chr2num = print_karyotype_file(f1_file, "sp1")
	print_karyotype_file(f2_file, "sp2")
	t1_dict = parse_table(t1_file)
	t2_dict = parse_table(t2_file)
	print_links_file(t1_dict, t2_dict, chr2num)
