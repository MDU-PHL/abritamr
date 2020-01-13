import logging
# create logger with 'spam_application'
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# # create file handler which logs even debug messages
fh = logging.FileHandler('abritamr.log')
fh.setLevel(logging.INFO)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(ch)
logger.addHandler(fh)
