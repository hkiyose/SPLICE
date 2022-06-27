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

parser.add_argument("-i", "--input", help="input directory (output directory of annot.py)", required=True)
parser.add_argument("-r", "--ref", help="reference transcriptome file (output file of ref_exonnum.py)", required=True)
parser.add_argument("-o", "--output", help="output directory ('./out/exp')", default="./out/exp")
parser.add_argument("--tmp", help="temporary directory ('./tmp/exp_tmp')", default="./tmp/exp_tmp")
parser.add_argument("-w", "--workers", help="number of threads", default="1")

parser.add_argument("--min_read_num", help="minimum number of support reads (3)", default=3)
parser.add_argument("--min_read_freq", help="minimum frequency of the transcript in the total amount of the gene (percentage) (1)", default=1)
parser.add_argument("--range_sj_eva", help="maximum change of novel exon length for evaluate the error rate (20)", default=20)
parser.add_argument("--err_rate_filt", help="mapping error rate cutoff at splicing junctino sites (percentage) (20)", default=20)
parser.add_argument("--min_novel_len_gap", help="minimum change of novel exon length (5)", default=5)
parser.add_argument("--min_novel_exon_len", help="minimum length of novel exon (60)", default=60)


args = parser.parse_args()

if not os.path.exists(args.output):
	os.makedirs(args.output)

if not os.path.exists(args.tmp):
	os.makedirs(args.tmp)

#  Merging tables (1/14)
class MergeTable(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.tmp, "exp.table"))
		
	def run(self):
		with tools.SetIO(self.output().path):	
			tools.merge_table(args.input)

# Filtering out low-frequency Novel splicing variants (2/14)	
@requires(MergeTable)
class FreqFilt(luigi.Task):
	task_namespace = 'tasks'

	def output(self):
		return luigi.LocalTarget(self.input().path + "2")

	def run(self):
		with tools.SetIO(self.output().path):
			tools.freq_filter(self.input().path, args.min_read_num, args.min_read_freq)

# Change the format, evaluate the splicing junction error rate and remove unreliable novel splicing variants (3/14)
@requires(FreqFilt, MergeTable)
class ChangeFormat2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.format_change2("./tools/gene_type/genecode.txt", "./tools/gene_type/refseq.txt",  args.ref, self.input()[0].path, "30", args.err_rate_filt, args.min_novel_len_gap, args.range_sj_eva, args.range_sj_eva)

# Merge each candidate by category (4/14)
# Known candidates
@requires(ChangeFormat2)
class MergeCategoryKnown(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".known")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_known(self.input().path)

# Novel exon length candidates
@requires(ChangeFormat2)
class MergeCategoryLen(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".len")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_len_comb(args.ref, self.input().path, "LENGTH")
	
@requires(MergeCategoryLen)
class MergeCategoryLen2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + "2")
	
	def run(self):
		subprocess.run(['sort -k 4,4 -k 7nr,7 ' + self.input().path + ' > ' + self.input().path + '.sort'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.merge_len_comb2(self.input().path + '.sort')

@requires(MergeCategoryLen2, MergeCategoryLen)
class MergeCategoryLen3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_len_comb3(self.input()[0].path)
			
# Novel exon combination candidates
@requires(ChangeFormat2)
class MergeCategoryComb(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".comb")
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_len_comb(args.ref, self.input().path, "COMBINATION")

@requires(MergeCategoryComb)
class MergeCategoryComb2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + "2")

	def run(self):
		subprocess.run(['sort -k 4,4 -k 7nr,7 ' + self.input().path + ' > ' + self.input().path + '.sort'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.merge_len_comb2(self.input().path + '.sort')

@requires(MergeCategoryComb2, MergeCategoryComb)
class MergeCategoryComb3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_len_comb3(self.input()[0].path)

# Merging novel unannotated exon tables (5/14)
class MergeNovelExon(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(join(args.tmp, "novel_exon.table"))

	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_novel_exon(args.input)

@requires(MergeNovelExon)
class MergeNovelExon2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + "2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_novel_exon2(self.input().path)

@requires(MergeNovelExon2, MergeNovelExon)
class MergeNovelExon3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_novel_exon3(self.input()[0].path)

@requires(MergeNovelExon3, MergeNovelExon)
class MergeNovelExon4(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "4")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_novel_exon4(self.input()[0].path, self.input()[1].path)

@requires(MergeNovelExon4)
class MergeNovelExon5(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".conv")	
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_novel_exon5(self.input().path)

# Filtering out low-frequency Novel unannot exon candidates (6/14)
@requires(MergeTable, MergeNovelExon5, MergeNovelExon)
class freqFilt2(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[2].path + "5")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.freq_filter2(self.input()[0].path, self.input()[1].path)

# Extracting Novel unannot exon candidates detected together with the gene region (7/14)
@requires(freqFilt2, MergeNovelExon)
class GetNovelExon(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "6") 
		
	def run(self):
		with tools.SetIO(self.output().path):
			tools.novel_exon_filter("./tools/gene_type/genecode.txt", "./tools/gene_type/refseq.txt", self.input()[0].path)

# Removal of Fusion gene candidates containing Novel exon (8/14)
@requires(GetNovelExon, MergeNovelExon)
class RmvFusion(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "7")
	
	def run(self):
		subprocess.run(['sort -k4 ' + self.input()[0].path + ' > ' + self.input()[0].path + '.sort'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.novel_exon_filter2(args.ref, self.input()[0].path + '.sort')

# Merging novel unannot exon candidates (9/14)
@requires(RmvFusion, MergeNovelExon)
class MergeNovelExon6(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "8")
	
	def run(self):
		subprocess.run(['sort -k 4,4 -k 7nr,7 ' + self.input()[0].path + ' > ' + self.input()[0].path + '.sort'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		with tools.SetIO(self.output().path):
			tools.merge_novel_exon6(self.input()[0].path + '.sort')

@requires(MergeNovelExon6, MergeNovelExon)
class MergeNovelExon7(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "9")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_len_comb3(self.input()[0].path)

# Conversion of novel exon location information to the most frequent location, filtering out unreliable candidates (10/14)
@requires(MergeNovelExon)
class MergeNovelExon8(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input().path + ".sort")
		
	def run(self):
		with tools.SetIO(self.output().path):
			tools.novel_exon_sort(self.input().path)

@requires(MergeNovelExon7, MergeNovelExon)
class MergeNovelExon9(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "10")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.novel_exon_convert(self.input()[0].path)

@requires(MergeNovelExon7, MergeNovelExon8, MergeNovelExon9, MergeNovelExon)
class MergeNovelExon10(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(self.input()[3].path + "11")

	def run(self):
		with tools.SetIO(self.output().path):
			tools.novel_exon_convert2(self.input()[0].path, self.input()[1].path, args.ref, self.input()[2].path)
	
@requires(MergeNovelExon10, MergeNovelExon7, MergeNovelExon)
class MergeNovelExon11(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[2].path + "12")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.novel_exon_convert3(args.ref, self.input()[0].path, args.min_novel_exon_len, args.min_novel_exon_len, self.input()[1].path)

# Integrate each table (11/14)
@requires(MergeCategoryKnown, MergeNovelExon11, MergeCategoryComb3, MergeCategoryLen3, ChangeFormat2)
class IntegrateTable(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[4].path + ".merge")
	
	def run(self):
		subprocess.run(['cp ' + self.input()[0].path + ' ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		subprocess.run(['cat ' + self.input()[1].path + ' >> ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		subprocess.run(['cat ' + self.input()[2].path + ' >> ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		subprocess.run(['cat ' + self.input()[3].path + ' >> ' + self.output().path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Add information such as splicing variant frequency and total read count. Exclude low-frequency candidates (12/14)
@requires(MergeTable, IntegrateTable)
class changeFormat3(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "2")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.format_change3(self.input()[0].path, self.input()[1].path)

@requires(changeFormat3, IntegrateTable)
class changeFormat4(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "3")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.format_change4(self.input()[0].path, args.min_read_freq, args.min_read_num)

# Integrate 3'UTR information (13/14)
@requires(changeFormat4, IntegrateTable)
class Merge3utr(luigi.Task):
	task_namespace = 'tasks'
		
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "4")
	
	def run(self):
		with tools.SetIO(self.output().path):
			tools.merge_3utr(self.input()[0].path)

# Modify the coding type and generate the final file (14/14)
@requires(Merge3utr, IntegrateTable)
class changeFormat5(luigi.Task):
	task_namespace = 'tasks'	
	
	def output(self):
		return luigi.LocalTarget(self.input()[1].path + "5")

	def run(self):
		with tools.SetIO(self.output().path):
			tools.format_change5("./tools/gene_type/genecode.txt", "./tools/gene_type/refseq.txt", args.ref, self.input()[0].path)

@requires(changeFormat5)
class changeFormat6(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.output, "expression.tsv"))
		
	def run(self):
		with tools.SetIO(self.output().path):
			tools.format_change6(self.input().path) 

if __name__ == '__main__':
	luigi.run(['tasks.changeFormat6', '--workers', '1', '--local-scheduler'])

	
