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

sqs = SQSWrapper('default', 'igv')

if len(sys.argv) != 2:
    print('Input file required')
    sys.exit(-1)

inputfile = Path(sys.argv[1])

# if not inputfile.exists():
#     print('Input file does not exist')
#     sys.exit(-1)

logger.info('Attempting to queue items')

# with open(inputfile, 'r') as f:
#     items = f.readlines()
#     for item in items:
#         work_item = item.strip()
#         sqs.queue_work_item(work_item)

logger.info('Queueing item')
work_items = ['SKCM|E|TCGA-BF-A5EQ-01A|chr10|101002946|101003691|chr10_101002946_101003691_D_chr10:101002688-101002689', 
'SKCM|E|TCGA-BF-AAP2-01A|chr10|123154485|123154942|chr10_123154485_123154942_D_chr10:123154933-123154934',
'SKCM|E|TCGA-EB-A44O-01A|chr10|13600348|13610039|chr10_13600348_13610039_NDA_chr10:13610100-13610101',
'SKCM|E|TCGA-EB-A41B-01A|chr10|17236379|17237312|chr10_17236379_17237312_D_chr10:17234785-17234786']
for work_item in work_items:
    sqs.queue_work_item(work_item)

