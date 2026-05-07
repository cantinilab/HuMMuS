try:
    import logging
    import absl.logging
    logging.root.removeHandler(absl.logging._absl_handler)
    absl.logging._warn_preinit_stderr = False
except Exception as e:
    pass
from hummuspy.utils.logging.hummus_logging import *
from autologging import logged
from hummuspy.utils.profiler_report import *
