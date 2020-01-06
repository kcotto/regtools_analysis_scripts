import datetime
import json
import logging

from sqs_wrapper import SQSWrapper

logger = logging.getLogger('sqs_status')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')
formatter.default_msec_format = '%s.%03d'
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

sqs = SQSWrapper('default')
logger.info('Awaiting status messages')
for status_item in sqs.status_items():
    work_item = status_item['work_item']
    message = status_item['message']
    is_fatal = status_item['is_fatal']
    is_complete = status_item['is_complete']
    hostname = status_item['hostname']
    logger.info(f'Got status item "{status_item}"')
    if status_item['is_fatal']:
        with open('failures_v5.tsv', 'w') as f:
            f.write(f'{datetime.datetime.now()}\t{work_item}\t{message}\n')