import logging

from sqs_wrapper import SQSWrapper

logger = logging.getLogger('sqs_reader')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')
formatter.default_msec_format = '%s.%03d'
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

sqs = SQSWrapper('default')

for work_item in sqs.work_items():
    logger.info(f'Got work item "{work_item}"')
    sqs.report_status(work_item, f'Reporting status for item {work_item}', is_fatal=False, is_complete=True)
