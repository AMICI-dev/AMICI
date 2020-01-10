from __future__ import absolute_import
import logging
import platform
import socket
import amici
import os
import warnings
import pysb

SECONDS_IN_HOUR = 3600
LOG_LEVEL_ENV_VAR = 'AMICI_LOG'
BASE_LOGGER_NAME = 'amici'
NAMED_LOG_LEVELS = {'NOTSET': logging.NOTSET,
                    'DEBUG': logging.DEBUG,
                    'INFO': logging.INFO,
                    'WARNING': logging.WARNING,
                    'ERROR': logging.ERROR,
                    'CRITICAL': logging.CRITICAL}


def _setup_logger(level=logging.WARNING, console_output=True,
                  file_output=False, capture_warnings=True):
    """Set up a new logging.Logger for AMICI logging
    Arguments:

        level: Logging level, typically using a constant like logging.INFO or
            logging.DEBUG

        console_output: Set up a default console log handler if True (default)

        file_output: Supply a filename to copy all log output to that file, or
            set to False to disable (default)

        capture_warnings: Capture warnings from Python's warnings module
            if True (default)
    Returns:
        A logging.Logger object for PySB logging. Note that other PySB modules
        should use a logger specific to their namespace instead by calling
        :func:`get_logger`.

    Raises:
    """
    log = logging.getLogger(BASE_LOGGER_NAME)

    # Logging level can be overridden with environment variable
    if LOG_LEVEL_ENV_VAR in os.environ:
        try:
            level = int(os.environ[LOG_LEVEL_ENV_VAR])
        except ValueError:
            # Try parsing as a name
            level_name = os.environ[LOG_LEVEL_ENV_VAR]
            if level_name in NAMED_LOG_LEVELS.keys():
                level = NAMED_LOG_LEVELS[level_name]
            else:
                raise ValueError(f'Environment variable {LOG_LEVEL_ENV_VAR} '
                                 f'contains an invalid value "{level_name}".'
                                 f' If set, its value must be one of '
                                 f'{", ".join(NAMED_LOG_LEVELS.keys())}'
                                 f' (case-sensitive) or an integer log level.')

    log.setLevel(level)

    # Remove default logging handler
    log.handlers = []

    log_fmt = logging.Formatter('%(asctime)s.%(msecs).3d - %(name)s - '
                                '%(levelname)s - %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S')

    if console_output:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(log_fmt)
        log.addHandler(stream_handler)

    if file_output:
        file_handler = logging.FileHandler(file_output)
        file_handler.setFormatter(log_fmt)
        log.addHandler(file_handler)

    log.info('Logging started on AMICI version %s', amici.__version__)

    log.debug('OS Platform: %s', platform.platform())
    log.debug('Python version: %s', platform.python_version())
    log.debug('Hostname: %s', socket.getfqdn())

    logging.captureWarnings(capture_warnings)

    return log


def get_logger(logger_name=BASE_LOGGER_NAME, log_level=None, **kwargs):
    """ Returns (if extant) or creates an AMICI logger
    If the AMICI base logger has already been set up, this method will
    return it or any of its descendant loggers without overriding the
    settings - i.e. any values supplied as kwargs will be ignored.
    Parameters:

        logger_name : Get a logger for a specific namespace, typically __name__
            for code outside of classes or self.__module__ inside a class

        log_level : Override the default or preset log level for the requested
            logger. None or False uses the default or preset value. True
            evaluates to logging.DEBUG. Any integer is used directly.

        **kwargs : kwargs
            Keyword arguments to supply to :func:`setup_logger`. Only used when
            the AMICI logger hasn't been set up yet (i.e. there have been no
            calls to this function or :func:`get_logger` directly).

    Returns:
        A logging.Logger object with the requested name

    Raises:

    """
    if BASE_LOGGER_NAME not in logging.Logger.manager.loggerDict.keys():
        _setup_logger(**kwargs)
    elif kwargs:
        warnings.warn('AMICI logger already exists, ignoring keyword '
                      'arguments to setup_logger')

    logger = logging.getLogger(logger_name)

    if log_level is not None and log_level is not False:
        if isinstance(log_level, bool):
            log_level = logging.DEBUG
        elif not isinstance(log_level, int):
            raise ValueError('log_level must be a boolean, integer or None')

        if logger.getEffectiveLevel() != log_level:
            logger.debug('Changing log_level from %d to %d' % (
                logger.getEffectiveLevel(), log_level))
            logger.setLevel(log_level)

    return logger
