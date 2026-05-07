# importing libraries
import os
import psutil
import datetime
import time

from hummuspy.utils import *
import logging
logger = logging.getLogger("profiler")

# inner psutil function
def process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss
    
# decorator function
def profile_memory(func):
    def wrapper(*args, **kwargs):
        mem_before = process_memory()
        result = func(*args, **kwargs)
        mem_after = process_memory()
        
        delta = mem_after - mem_before
        msg = "Function {}: Memory Trace - Before: {} | After: {} | Increment: {}".format( func.__name__, mem_before, mem_after, delta)
        
        logger.info( msg )
        return result
    return wrapper

def report_execution_information(stage, log_path):
    def decorator(func):
        name = func.__name__
        def wrapper(*args, **kwargs):
            before = datetime.datetime.now()
            mem_before = process_memory()
            result = func(*args, **kwargs)
            mem_after = process_memory()
            after = datetime.datetime.now()  
            
            datefmt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            delta_time = after-before
            formatted = str(datetime.timedelta(seconds = delta_time.seconds))
            
            delta = mem_after - mem_before
            
            msg = "[HuMMuS - Time trace - %s] Stage: %s >> Function: %s - Execution time: %s\n" %(datefmt, stage, name, formatted)
            msg += "[HuMMuS - Memory trace - %s] Stage: %s >> Function: %s - Before: %s | After: %s | Increment: %s\n" %(datefmt, stage, name, mem_before, mem_after, delta)
            with open( log_path,"a") as f:
                f.write( msg )
                
            return result
        wrapper.__name__ = func.__name__
        return wrapper
    return decorator

# --------- Test area
@report_execution_information("Initialization", "test.txt")
@profile_memory
def test_func():
    x = [1] * (10 ** 7)
    y = [2] * (4 * 10 ** 8)
    del x
    time.sleep(2)
    return 'done'
    
if( __name__ == "__main__" ):
    test_func()
