import os
import sys
import argparse
import itertools 

parser = argparse.ArgumentParser(description='Process some inputs')

parser.add_argument('-r', '--ref', type=str, help='Path to reference genome (with prefix)')
parser.add_argument('-f', '--fastq_dir', type=str, help='Path to fastq_files')
parser.add_argument('-o', '--out_dir', type=str, help='output directory')
parser.add_argument('-s', '--samples_file', type=str, help='samples.txt file (having sample names and metadata)')
parser.add_argument('-c', '--controls_file', type=str, help='controls.txt file (having control names and metadata)')
parser.add_argument('-v', '--snp_file', type=str, help='snp file required for SNPsplit')
args = parser.parse_args()

ref = args.ref
fastq_dir = args.fastq_dir
out_dir = args.out_dir
samples_file = args.samples_file
controls_file = args.controls_file
snp_file = args.snp_file

sample_names = open(samples_file,'r')
control_names = open(controls_file, 'r')

sample_control = []
controls = []
trimming = []
encodings = []
libraries = []
peak_types = []
effective_genomes = []
allele_specifics = []
alleles_1 = []
alleles_2 = []

total_samples = 0

for line in sample_names:
	total_samples += 1
	fields = line.rstrip().split("\t")
	sample_control.append(fields[0].rstrip())
	controls.append(fields[1].rstrip())
	trimming.append(fields[2].rstrip())
	encodings.append(fields[3].rstrip())
	libraries.append(fields[4].rstrip())
	effective_genomes.append(fields[5].rstrip())
	allele_specifics.append(fields[6].rstrip())
	alleles_1.append(fields[7].rstrip())
	alleles_2.append(fields[8].rstrip())
	peak_types.append(fields[9].rstrip())


for line in control_names:
	fields = line.rstrip().split("\t")
	sample_control.append(fields[0].rstrip())
	trimming.append(fields[1].rstrip())
	encodings.append(fields[2].rstrip())
	libraries.append(fields[3].rstrip())
	effective_genomes.append(fields[4].rstrip())
	allele_specifics.append(fields[5].rstrip())
	alleles_1.append(fields[6].rstrip())
	alleles_2.append(fields[7].rstrip())
	peak_types.append(fields[8].rstrip())


sample_names.close()
control_names.close()

##########################################################################################################################################################################################################
################################################################################## bowtie alignment pbs script creation ##################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a STAR-align script\n#PBS -N STAR-align\n#PBS -l nodes=1:ppn=28\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o STAR-align_pbs.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	qual_encoding = encodings[i]
	library = libraries[i]
	trimmed = trimming[i]

	pbs_file_name = out_dir + folder_name + "/STAR_align.pbs"

	os.system("mkdir -p " + out_dir + folder_name)

	cmd_1 = "module load STAR\n\n"
	cmd_2 = "cd " + out_dir + folder_name + "\n\n"

	if (library == "paired"):
		if (trimmed == "no"):
			cmd_3 = "( time STAR --runMode alignReads --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --runThreadN 24 --genomeDir " + ref + " --readFilesIn " + fastq_dir + "raw/" + folder_name + "_1.fastq.gz " + fastq_dir + "raw/" + folder_name + "_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix " + folder_name + " --outSAMtype BAM SortedByCoordinate ) > STAR_align_out.log 2>&1\n\n"
		elif (trimmed == "yes"):
			cmd_3 = "( time STAR --runMode alignReads --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --runThreadN 24 --genomeDir " + ref + " --readFilesIn " + fastq_dir + "trimmed/" + folder_name + "_1_val_1.fq.gz " + fastq_dir + "trimmed/" + folder_name + "_2_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix " + folder_name + " --outSAMtype BAM SortedByCoordinate ) > STAR_align_out.log 2>&1\n\n"	
	elif (library == "single"):
		if (trimmed == "no"):
			cmd_3 = "( time STAR --runMode alignReads --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --runThreadN 24 --genomeDir " + ref + " --readFilesIn " + fastq_dir + "raw/" + folder_name + "_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix " + folder_name + " --outSAMtype BAM SortedByCoordinate ) > STAR_align_out.log 2>&1\n\n"
		elif (trimmed == "yes"):
			cmd_3 = "( time STAR --runMode alignReads --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --runThreadN 24 --genomeDir " + ref + " --readFilesIn " + fastq_dir + "trimmed/" + folder_name + "_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix " + folder_name + " --outSAMtype BAM SortedByCoordinate  ) > STAR_align_out.log 2>&1\n\n"

	cmd_4 = "qsub samtools_markdup.pbs\n"
	pbs_script = open(pbs_file_name,'w+')
	pbs_script.write('%s\n'%(pbs_header + cmd_1 + cmd_2 + cmd_3 + cmd_4))
	pbs_script.close()


##########################################################################################################################################################################################################
################################################################################## samtools markdup pbs script creation ##################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a samtools duplicate marking pipeline script\n#PBS -N mark_dup\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o samtools_markdup_pbs.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	allele_specific = allele_specifics[i]
	library = libraries[i]
	allele_1 = alleles_1[i]
	allele_2 = alleles_2[i]
	pbs_file_name = out_dir + folder_name + "/samtools_markdup.pbs"
	fixmate = ""

	cmd_1 = "module load samtools\n\n"
	cmd_2 = "cd " + out_dir + folder_name + "\n\n"
	cmd_7 = "( time samtools view -q 20 -@ 13 -O BAM -o " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam " + folder_name + "_NS" + fixmate + "_CS_dedup.bam ) > samtools_view_filter_out.log 2>&1\n\n"

	cmd_3 = "( time samtools sort -@ 13 -n -O BAM -o " + folder_name + "_NS.bam " + folder_name + "Aligned.sortedByCoord.out.bam ) > samtools_name_sort_out.log 2>&1\n"
	cmd_4 = ""
	if (library == "paired"):
		fixmate = "_fixmate"
		cmd_4 = "( time samtools fixmate -@ 13 -m -O BAM " + folder_name + "_NS.bam " + folder_name + "_NS" + fixmate + ".bam ) > samtools_fixmate_out.log 2>&1\n"
		cmd_7 = "( time samtools view -q 20 -f 3 -@ 13 -O BAM -o " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam " + folder_name + "_NS" + fixmate + "_CS_dedup.bam ) > samtools_view_filter_out.log 2>&1\n\n"
	cmd_5 = "( time samtools sort -@ 13 -O BAM -o " + folder_name + "_NS" + fixmate + "_CS.bam " + folder_name + "_NS" + fixmate + ".bam ) > samtools_coord_sort_out.log 2>&1\n"
	cmd_6 = "( time samtools markdup -@ 13 -r -O BAM " + folder_name + "_NS" + fixmate + "_CS.bam " + folder_name + "_NS" + fixmate + "_CS_dedup.bam ) > samtools_markdup_out.log 2>&1\n"
	cmd_8 = ("rm " + folder_name + "_NS.bam " + folder_name + "_NS_fixmate.bam " + folder_name + "_NS" + fixmate + "_CS.bam " + folder_name + "_NS_fixmate_CS_dedup.bam\n\n"  )
	cmd_9 = ("samtools index " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam\n")
	cmd_10 = ("samtools idxstats " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam > " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered_idxstats.txt\n\n")

	if (allele_specific == "yes"):
		cmd_11 = "qsub SNPsplit.pbs\n"
		cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4 + cmd_5 + cmd_6 + cmd_7 + cmd_8 + cmd_9 + cmd_10 + cmd_11
	else:
		cmd = cmd_1 + cmd_2 + cmd_3 + cmd_5 + cmd_6 + cmd_7 + cmd_8 + cmd_9 + cmd_10
		cmd = cmd + "cd " + out_dir + "\ncd ..\n\n"
		cmd = cmd + "python scripts/bam_normalize.py " + out_dir + " " + folder_name + " " + allele_specific + " " + allele_1 + " " + allele_2 + " " + ref + " " + library + "\n\n"
		cmd = cmd + "cd " + out_dir + folder_name + "\n\n"
		cmd = cmd + "qsub bam_normalize.pbs\n"
		cmd = cmd + "qsub htseq-count.pbs\n"

	pbs_script = open(pbs_file_name,'w+')
	pbs_script.write('%s\n'%(pbs_header + cmd))
	pbs_script.close()


##########################################################################################################################################################################################################
################################################################################## SNP_split pbs script creation #########################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a SNPsplit script\n#PBS -N SNPsplit\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o SNPsplit_pbs.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	library = libraries[i]
	allele_specific = allele_specifics[i]
	allele_1 = alleles_1[i]
	allele_2 = alleles_2[i]

	if (allele_specific == "yes"):
		pbs_file_name = out_dir + folder_name + "/SNPsplit.pbs"
		cmd_1 = "module load SNPsplit\nmodule load samtools\n\n"
		cmd_2 = "cd " + out_dir + folder_name + "\n\n"

		if (library == "paired"):
			cmd_3 = "( time SNPsplit --paired --no_sort --snp_file "  + snp_file + " " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam ) > SNPsplit_out.log 2>&1\n\n"
		elif (library == "single"):
			cmd_3 = "( time SNPsplit --bam --no_sort --snp_file "  + snp_file + " " + folder_name + "_NS_fixmate_CS_dedup_filtered.bam ) > SNPsplit_out.log 2>&1\n\n"

		cmd_4 = "rm " + folder_name + "_NS_fixmate_CS_dedup_filtered.allele_flagged.bam\n\n"
		cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4
		cmd = cmd + "samtools index " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bam\n"
		cmd = cmd + "samtools index " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bam\n\n"

		cmd = cmd + "cd " + out_dir + "\ncd ..\n\n"
		cmd = cmd + "python scripts/bam_normalize.py " + out_dir + " " + folder_name + " " + allele_specific + " " + allele_1 + " " + allele_2 + " " + ref + " " + library + "\n\n"
		cmd = cmd + "cd " + out_dir + folder_name + "\n\n"
		cmd = cmd + "qsub bam_normalize.pbs\n"
		cmd = cmd + "qsub htseq-count.pbs\n"

		pbs_script = open(pbs_file_name,'w+')
		pbs_script.write('%s\n'%(pbs_header + cmd))
		pbs_script.close()


##########################################################################################################################################################################################################
################################################################################## htseq-count pbs script creation #######################################################################################
##########################################################################################################################################################################################################

pbs_header = "# This is a htseq-count script\n#PBS -N htseq\n#PBS -l nodes=1:ppn=13\n#PBS -l walltime=10:00:00\n#PBS -q submit\n#PBS -j oe\n#PBS -o htseq-count.out\n#PBS -m abe\n#PBS -M semwal.a@wehi.edu.au\n\n"

for i in range(len(sample_control)):

	folder_name = sample_control[i]
	library = libraries[i]
	allele_specific = allele_specifics[i]
	allele_1 = alleles_1[i]
	allele_2 = alleles_2[i]

	fixmate = ""

	if (library == "paired"):
		fixmate = "_fixmate"

	if (allele_specific == "yes"):
		pbs_file_name = out_dir + folder_name + "/htseq-count.pbs"
		cmd_1 = "cd " + out_dir + folder_name + "\n\n"
		cmd_2 = "( time htseq-count -f bam -r pos --additional-attr=gene_name " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.bam /wehisan/general/academic/seq_data/Ayush/genomes/GRCm38/gencode.vM23.chr_patch_hapl_scaff.annotation_2.gtf >> " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome1.gene_counts.txt) > htseq-count_out.log 2>&1\n\n"
		cmd_3 = "( time htseq-count -f bam -r pos --additional-attr=gene_name " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.bam /wehisan/general/academic/seq_data/Ayush/genomes/GRCm38/gencode.vM23.chr_patch_hapl_scaff.annotation_2.gtf >> " + allele_1 + "_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.genome2.gene_counts.txt) > htseq-count_out.log 2>&1\n\n"

	cmd_4 = "( time htseq-count -f bam -r pos --additional-attr=gene_name " + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.bam /wehisan/general/academic/seq_data/Ayush/genomes/GRCm38/gencode.vM23.chr_patch_hapl_scaff.annotation_2.gtf >> comp_" + folder_name + "_NS" + fixmate + "_CS_dedup_filtered.gene_counts.txt) > htseq-count_out.log 2>&1\n\n"

	cmd = cmd_1 + cmd_2 + cmd_3 + cmd_4

	pbs_script = open(pbs_file_name,'w+')
	pbs_script.write('%s\n'%(pbs_header + cmd))
	pbs_script.close()


