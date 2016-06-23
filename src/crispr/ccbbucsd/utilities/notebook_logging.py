# standard libraries
import logging
import sys

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def set_stdout_info_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    str_handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(str_handler)