import datetime
import sys
import os 
import tools
from os.path import join, basename
import luigi
from luigi.util import requires
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="reference transcriptome (UCSC table browser output file)", required=True)
parser.add_argument("-g", "--genome", help="reference genome sequence (FASTA)", required=True)
parser.add_argument("-o", "--output", help="output directory", default="./out/ref_exonnum")
parser.add_argument("-p", "--prefix", help="prefix for the output files", default="RefGenes")
parser.add_argument("-w", "--workers", help="number of threads", default="1")
parser.add_argument("--tmp", help="temporary directory", default="./tmp/ref_exonnum_tmp")

args = parser.parse_args()

class AddExonNum(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.tmp, basename(args.input)+".AddExonNum"))
		
	def run(self):
		with tools.SetIO(self.output().path): # 標準出力を一時的にfile出力に変更し,その後console出力に戻すための処理
			tools.add_exonnum(args.input)
		
@requires(AddExonNum)
class RemoveConjoined(luigi.Task):
	task_namespace = 'tasks'

	def output(self):
		return luigi.LocalTarget(self.input().path + ".RemoveConjoined") #self.input().path: 依存先のタスク(AddExonNum)のoutput file

	def run(self):
		with tools.SetIO(self.output().path): #self.output().path: 現在のタスク(RemoveConjoined)のoutput file
			tools.rm_conjoined(self.input().path)

@requires(RemoveConjoined)
class ConvertRef(luigi.Task):
	task_namespace = 'tasks'

	def output(self):
		return luigi.LocalTarget(join(args.output, args.prefix))
	
	def run(self):
		with tools.SetIO(self.output().path):	
			tools.ref_convert(self.input().path)

class ConvertFa(luigi.Task):
	task_namespace = 'tasks'

	def output(self):
		return luigi.LocalTarget(join(args.tmp, basename(args.genome)+".ConvertFa"))
	
	def run(self):
		with tools.SetIO(self.output().path):	
			tools.convert_fa(args.genome)

@requires(ConvertFa, ConvertRef)#複数依存タスクが存在する場合、カンマ区切りでおｋ
class MakeRefSeq(luigi.Task):
	task_namespace = 'tasks'
	
	def output(self):
		return luigi.LocalTarget(join(args.output, args.prefix + ".fa"))

	def run(self):
		with tools.SetIO(self.output().path):	
			tools.make_ref_seq_cp(self.input()[0].path, self.input()[1].path) #依存タスクが複数存在する場合、self.input()[int].pathで入力ファイルを選択する


if __name__ == '__main__':
    # luigiの実行
#	luigi.configuration.LuigiConfigParser.add_config_path('./luigi.cfg')#config fileの場所を指定する
	luigi.run(['tasks.MakeRefSeq', '--workers', args.workers, '--local-scheduler'])

	
