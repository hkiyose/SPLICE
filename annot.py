import datetime
import sys
import os 
import tools
from os.path import join, basename, exists
import luigi
from luigi.util import requires
import subprocess
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="query sequence data (FASTQ)", required=True)
parser.add_argument("-g", "--genome", help="reference genome sequence (FASTA)", required=True)
parser.add_argument("-r", "--ref", help="reference transcriptome file (output file of ref_exonnum.py)", required=True)
parser.add_argument("-p", "--prefix", help="prefix for the output files", required=True)
parser.add_argument("-o", "--output", help="output directory", default="./out/annot")
parser.add_argument("--tmp", help="temporary directory", default="./tmp/annot_tmp")
parser.add_argument("-w", "--workers", help="number of threads", default="1")
parser.add_argument("--minimap2", help="execution file of minimap2", default="minimap2")

parser.add_argument("--bq_filt", help="read quality cutoff. Minimum average base quality score (15)", default=15)
parser.add_argument("--min_sc_len", help="minimum length of the softclip region to be remapped (60)", default=60)
parser.add_argument("--mq_filt", help="mapping quality cutoff (0)", default=0)
parser.add_argument("--min_fusion_dist", help="minimum distance of each transcript in the fusion transcript (200000)", default=200000)
parser.add_argument("--max_fusion_bp_merge", help="maximam distance to merge fusion gene breakpoints (5)", default=5)
parser.add_argument("--min_fusion_read", help="minimum number of support reads for fusion transcripts (1)", default=1)
parser.add_argument("--min_fusion_freq", help="minimum frequency of the fusion transcript in the total amount of the gene (percentage) (0.1)", default=0.1)
parser.add_argument("--mq_filt_novel_exon", help="mapping quality cutoff of novel exon (1)", default=1)

args = parser.parse_args()

if not os.path.exists(args.output):
	os.makedirs(args.output)

if not os.path.exists(args.tmp):
	os.makedirs(args.tmp)

# Base quality filtering (1/20)
class BqFilt(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.tmp, basename(args.input)+".BqFilt"))
		
	def run(self):
		with tools.SetIO(self.output().path):	
			tools.bq_filter(args.input, args.bq_filt)

# Mapping to reference genome (2/20)		
@requires(BqFilt)
class Minimap2_1(luigi.Task):
	task_namespace = 'tasks'

	def output(self):
		return luigi.LocalTarget(self.input().path + ".sam")

	def run(self):
		subprocess.run([args.minimap2 + ' -ax splice --cs ' + args.genome + ' ' + self.input().path + ' -o ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Convert the SAM(genome) format (3/20)
@requires(Minimap2_1)
class ConvSam(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".ConvSam")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.sam_convert(self.input().path)

@requires(ConvSam, BqFilt)
class ConvSam2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + "2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.sam_convert2(self.input()[0].path, self.input()[1].path, 0.3)

# Re-mapping the softclip region (4/20)
@requires(ConvSam2, BqFilt)
class GetSoftclip(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".softclip.fa")
		
	def run(self):
		with tools.SetIO(self.output().path):
			tools.get_softclip(self.input()[0].path, self.input()[1].path, args.min_sc_len)

@requires(GetSoftclip)
class Minimap2_2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".sam")
	
	def run(self):
		subprocess.run([args.minimap2 + ' -ax splice --cs ' + args.genome + ' ' + self.input().path + ' -o ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Integrating mapped softclip regions into converted SAM (5/20)
@requires(Minimap2_2)
class ConvSam3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".ConvSam3")
		
	def run(self):
		with tools.SetIO(self.output().path):
			tools.sam_convert(self.input().path)

@requires(ConvSam2, ConvSam3, BqFilt, Minimap2_2)
class ConvSam4(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[3].path + ".ConvSam4")
	
	def run(self):
		subprocess.run(["cat " + self.input()[1].path + " >> " + self.input()[1].path + ".MergeFile"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		subprocess.run(["cat " + self.input()[0].path + " >> " + self.input()[1].path + ".MergeFile"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		subprocess.run(["sort -k1 " + self.input()[1].path + ".MergeFile > " + self.input()[1].path + ".MergeFile2"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.sam_convert3(self.input()[1].path + ".MergeFile2", self.input()[2].path, 0.3)
		subprocess.run(["rm " + self.input()[1].path + ".MergeFile*"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		
# Gene annotation (6/20)
@requires(ConvSam4)
class AnnotExonnum(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".annot")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.annot_exonnum(args.ref, self.input().path, args.mq_filt, "1", "1")

# Splicing variant annotation (7/20)
@requires(AnnotExonnum)
class AnnotTranscript(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + "2")
	
	def run(self):	
		with tools.SetIO(self.output().path):
			tools.annot_transcript(args.ref, self.input().path)

#  Mapping to reference transcriptome (8/20)
@requires(BqFilt)
class Minimap2_3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".cDNA.sam")
	
	def run(self):
		subprocess.run([args.minimap2 + ' -a --cs ' + args.ref + ".fa" + ' ' + self.input().path + ' -o ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Add the number of alignment matches (Converted SAM(Genome)/SAM(Transcriptome)) (9/20)	
@requires(Minimap2_3)
class GetMatchRate(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".match_num")
		
	def run(self):
		with tools.SetIO(self.output().path):
			tools.match_rate(self.input().path)

@requires(AnnotTranscript)
class GetMatchRate_2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".sj")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.calc_error(self.input().path)

@requires(GetMatchRate, GetMatchRate_2)
class GetMatchRate_3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.match_rate2(self.input()[0].path, self.input()[1].path)

# Table (10/20)
@requires(GetMatchRate_3)
class ChangeFormat(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".table")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.format_change(self.input().path)

# Merge at position inside the TSS side exon and inside the TES side exon (11/20)
@requires(ChangeFormat, GetMatchRate_3)
class MergeTssTes(luigi.Task):
	task_namespace = 'tasks'

	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".conv")
	
	def run(self):
		subprocess.run(["sort -k3 " + self.input()[0].path + " > " + self.input()[0].path + ".sort"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.merge_transcript(self.input()[1].path, self.input()[0].path + ".sort")

@requires(MergeTssTes, GetMatchRate_3)
class MergeTssTes_2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".gap")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_transcript2(args.ref, self.input()[0].path, self.input()[1].path)

@requires(MergeTssTes_2)
class MergeTssTes_3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".mrg")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_transcript3(self.input().path)

# Detect novel exon combination candidates (12/20)
@requires(MergeTssTes_3)
class DetNovelExonComb(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".mid_comb")
	
	def run(self):
		subprocess.run(["sort -nr " + self.input().path + " > " + self.input().path + ".sort"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.det_novel(args.ref, self.input().path + ".sort", "mid_comb")

# Detect novel exon length candidates (13/20)
@requires(MergeTssTes_3, DetNovelExonComb)
class DetNovelExonLen(luigi.Task):
	task_namespace = 'tasks'	
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".mid_sl")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_novel(args.ref, self.input()[0].path + ".sort", "mid_sl")

#  Detect novel TES candidates (14/20)
@requires(MergeTssTes_2)
class DetNovelTes(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(self.input().path + ".mrg2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_novel_tes(self.input().path)

@requires(DetNovelTes)
class DetNovelTes_2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".known")
	
	def run(self):
		subprocess.run(["sort -nr " + self.input().path + " > " + self.input().path + ".sort"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.det_novel_tes2(args.ref, self.input().path + ".sort")

@requires(DetNovelTes_2, GetMatchRate_3)
class DetNovelTes_3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + "2")
	
	def run(self):
		subprocess.run(["sort -k3 " + self.input()[0].path + " > " + self.input()[0].path + ".sort"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.det_novel_tes3(self.input()[1].path, self.input()[0].path + ".sort")

@requires(DetNovelTes_3, DetNovelTes_2)
class DetNovelTes_4(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_novel_tes4(self.input()[0].path)

@requires(DetNovelTes_4, DetNovelTes_2)
class DetNovelTes_5(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "4")	
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_novel_tes5(args.ref, self.input()[0].path)

@requires(DetNovelTes_5, GetMatchRate_3, BqFilt)
class DetNovelTes_6(luigi.Task):	
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".pA")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_novel_tes6(args.ref, self.input()[1].path, self.input()[2].path, "15", "0.8", self.input()[0].path)

@requires(DetNovelTes_6, GetMatchRate_3)
class DetNovelTes_7(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".TES")
	
	def run(self):
		subprocess.run(["sort -nr -k2 " + self.input()[0].path + " > " + self.input()[0].path + ".sort"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.det_novel_tes7(args.ref, self.input()[0].path + ".sort", self.input()[1].path)

# Convert the SAM(transcriptome) format (15/20)
@requires(Minimap2_3)
class ConvSam5(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".conv") 
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.sam_convert(self.input().path)

@requires(ConvSam5, BqFilt)
class ConvSam6(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + "2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			 tools.sam_convert2(self.input()[0].path, self.input()[1].path, "0.3")

# Detect fusion gene candidates (16/20)
@requires(AnnotTranscript, ConvSam6, BqFilt)
class DetFusion(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[0].path + ".fusion")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_fusion(args.ref, self.input()[0].path, self.input()[1].path, self.input()[2].path, "15", "0.8")

@requires(DetFusion)
class DetFusion_2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		 return luigi.LocalTarget(self.input().path + "2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_fusion2(self.input().path, args.min_fusion_dist)

@requires(DetFusion_2, DetFusion)
class DetFusion_3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_fusion3(self.input()[0].path, args.max_fusion_bp_merge, args.min_fusion_read, args.min_fusion_freq)

@requires(DetFusion_3)
class DetFusion_4(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.output, args.prefix + ".fusion"))
	
	def run(self):
		subprocess.run(["sort -nr " + self.input().path + " > " + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Add information on the number of alignment matches (Converted SAM(Transcriptome)) (17/20)
@requires(ConvSam6)
class ConvSam7(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".map_rate")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.sam_convert4(self.input().path)

# Detect novel position candidates (18/20)
@requires(ConvSam7, AnnotTranscript, BqFilt)
class DetNovelExon(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.output, args.prefix + ".novel_exon"))
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.det_novel2(self.input()[1].path, args.mq_filt_novel_exon, self.input()[0].path, "0.8", "60", self.input()[2].path, "15", "0.8")

# Expression table (19/20)		
@requires(DetNovelExonComb, DetNovelExonLen, DetNovelTes_7, DetNovelTes_5, GetMatchRate_3, DetFusion_4, DetNovelExon)
class GetExpTable(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(join(args.output, args.prefix + ".annot"))
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.get_table(self.input()[0].path, self.input()[1].path, self.input()[2].path, self.input()[3].path, self.input()[4].path)

# Analysis report (20/20)
@requires(BqFilt, ConvSam2, GetExpTable)
class GetReport(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(join(args.output, args.prefix +".stat"))

	def run(self):
		with tools.SetIO(self.output().path):
			tools.get_report(args.input, self.input()[0].path, self.input()[1].path)		

if __name__ == '__main__':
	luigi.run(['tasks.GetReport', '--workers', args.workers, '--local-scheduler'])

	
