import logging
import sys
from pathlib import Path

from sqs_wrapper import SQSWrapper

logger = logging.getLogger('sqs_writer')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')
formatter.default_msec_format = '%s.%03d'
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

sqs = SQSWrapper('default', 'vep')

if len(sys.argv) != 2:
    print('Input file required')
    sys.exit(-1)

inputfile = Path(sys.argv[1])

if not inputfile.exists():
    print('Input file does not exist')
    sys.exit(-1)

logger.info('Attempting to queue items')

with open(inputfile, 'r') as f:
    items = f.readlines()
    for item in items:
        work_item = item.strip()
        sqs.queue_work_item(work_item)

# logger.info('Queueing item')
# work_items = ['LIHC;TCGA-2V-A95S-01A']
# for work_item in work_items:
#     sqs.queue_work_item(work_item)

