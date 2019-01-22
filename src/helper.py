#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/22 3:16 PM
import logging



def set_logging(log_name):
    """
    a DEBUG logger
    :param log_name:
    :return:
    """
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        fmt="[%(asctime)s] - [%(levelname)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    handler.setFormatter(formatter)
    handler.setLevel(logging.DEBUG)

    log = logging.getLogger(log_name)
    log.setLevel(logging.DEBUG)
    log.addHandler(handler)

    return log
