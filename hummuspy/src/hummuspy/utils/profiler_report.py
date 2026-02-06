# importing libraries
import os
import psutil

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
        logger.info("Function {}: Memory Trace - Before: {}, After: {}, Increment: {}".format( func.__name__, mem_before, mem_after, mem_after - mem_before) )
        return result
    return wrapper

@profile_memory
def func():
    x = [1] * (10 ** 7)
    y = [2] * (4 * 10 ** 8)
    del x
    return 'done'
    
print(__name__)
if( __name__ == "__main__" ):

    func()
