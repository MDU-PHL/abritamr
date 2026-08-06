import logging
import sys

logging.basicConfig(format = '[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', level=logging.INFO) 
handler = logging.StreamHandler(sys.stderr)

log = logging.getLogger(__name__)
log.addHandler(handler)