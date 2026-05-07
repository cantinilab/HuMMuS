import os
import sys
from subprocess import run, PIPE
from hummuspy.utils import *

wrkdir = sys.argv[1]
config_file = sys.argv[2]
log_path = os.path.join(wrkdir, "execution_log.txt")

ready_dir = os.path.join( wrkdir, "ready")

@report_execution_information("Layers and networks generation step", log_path)
def execute_generate_networks_layers():
	rscript = os.path.join( os.path.dirname( os.path.abspath(__file__) ), "generate_networks_layers.R")
	out = run( ["Rscript", rscript, "-d", wrkdir, "-f", config_file], stdout=PIPE, stderr=PIPE )
	if( out.returncode != 0 ):
		raise Exception(out.stderr.decode())
	else:
		ready_file = os.path.join( ready_dir, "analysis_step2_done.txt")
		f = open( ready_file, "w" )
		f.close()

execute_generate_networks_layers()