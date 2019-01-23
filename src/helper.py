#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/22 3:16 PM
import os
import logging
from logging.handlers import RotatingFileHandler


def set_logging(log_name):
    """
    return a two logger object, one for stream and another for the log file

    :param log_name:
    :return:
    """
    formatter = logging.Formatter(
        fmt="[%(asctime)s] - [%(levelname)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # streamhandler
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    sh.setLevel(logging.INFO)

    # filehandler

    fh = RotatingFileHandler(
        os.path.join("{}_error.log".format(log_name)),
        mode='a+',
        encoding="utf-8",
        maxBytes=20 * 1024 * 1024,
        backupCount=2,
        delay=0
    )

    fh.setFormatter(formatter)
    fh.setLevel(logging.ERROR)

    log = logging.getLogger(log_name)
    log.setLevel(logging.DEBUG)
    log.addHandler(sh)
    log.addHandler(fh)

    return log
