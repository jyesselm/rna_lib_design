import logging
import colorlog
import sys

APP_LOGGER_NAME = "RLD"

log_colors = {
    "DEBUG": "cyan",
    "WARNING": "yellow",
    "ERROR": "red",
    "CRITICAL": "bold_red",
}


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=False, file_name=None):
    """
    Set up the logger for the app
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    log_format = "%(levelname)-8s %(name)-12s %(message)s"
    formatter = colorlog.ColoredFormatter("%(log_color)s"+log_format, log_colors=log_colors)

    # pylint: disable=C0103
    sh = colorlog.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(logging.Formatter(log_format))
        logger.addHandler(fh)

    return logger


def get_logger(module_name):
    """
    Get the logger for the module
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
