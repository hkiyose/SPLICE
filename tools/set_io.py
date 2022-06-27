import sys
import os 

#標準出力を一時的にfile出力に変更し、その後console出力に戻す
class SetIO():
	def __init__(self, filename: str):
		self.filename = filename
	def __enter__(self):
		sys.stdout = _STDLogger(out_file=self.filename)
	def __exit__(self, *args):
		sys.stdout = sys.__stdout__

class _STDLogger():
	def __init__(self, out_file='out.log'):
		self.log = open(out_file, "a+")
	def write(self, message):
		self.log.write(message)
	def flush(self):
		pass

