"""Find default logging config file."""
import os as _os
import logging.config as _logging

_log_conf = _os.path.join(_os.path.dirname( _os.path.abspath(__file__)), 'logging_conf.ini')
if not _os.path.isfile(_log_conf):
    print("ERROR initializing log: %s not found!" % _log_conf)
else:
    _logging.fileConfig(_log_conf)
