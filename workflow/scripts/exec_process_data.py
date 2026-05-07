import os
import sys
from subprocess import run, PIPE
from hummuspy.utils import *

wrkdir = sys.argv[1]
config_file = sys.argv[2]
log_path = os.path.join(wrkdir, "execution_log.txt")

ready_dir = os.path.join( wrkdir, "ready")
if( not os.path.exists(ready_dir) ):
	os.makedirs(ready_dir)

@report_execution_information("Data Processing step", log_path)
def execute_process_data():
	rscript = os.path.join( os.path.dirname( os.path.abspath(__file__) ), "process_data.R")
	out = run( ["Rscript", rscript, "-d", wrkdir, "-f", config_file], stdout=PIPE, stderr=PIPE )
	if( out.returncode != 0 ):
		raise Exception(out.stderr.decode())
	else:
		with open(log_path, 'a') as f:
			f.write("\n"+out.stdout.decode())

		ready_file = os.path.join( ready_dir, "analysis_step1_done.txt")
		f = open( ready_file, "w" )
		f.close()

execute_process_data()