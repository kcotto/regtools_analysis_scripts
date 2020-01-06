"""Simple wrapper for AWS SQS"""
import logging
import socket
import json
import boto3

"""Simple wrapper for AWS SQS"""
import logging
import socket
import json
import boto3
class SQSWrapper():
    """
    Simple wrapper class for AWS SQS
    Creates or attaches to two queues named f'{base_queue_name}_work.fifo' and f'{base_queue_name}_status
    Note: the 'base_queue_name' parameter MUST be provided in this version otherwise an exception will
          be raised
    """
    def __init__(self, profile_name: str, base_queue_name: str = None, log_level: int = logging.INFO):
        # Get the service resource
        self.__logger = logging.getLogger('SQSWrapper')
        self.__logger.setLevel(log_level)
        if base_queue_name is None:
            raise Exception('SQSWrapper: base_queue_name must be passed to the constructor')
        self.__base_queue_name = base_queue_name
        self.__work_queue_name = f'{self.__base_queue_name}_work.fifo'
        self.__status_queue_name = f'{self.__base_queue_name}_status'
        formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)s:%(message)s')
        formatter.default_msec_format = '%s.%03d'
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        self.__logger.addHandler(stream_handler)
        self.__logger.debug(f'Initializing SQS wrapper from profile "{profile_name}"')
        self.__session = boto3.Session(profile_name=profile_name)
        self.__sqs = self.__session.resource('sqs')
        self.__hostname = socket.gethostname()
        try:
            self.__work_item_fifo = self.__sqs.get_queue_by_name(QueueName=self.__work_queue_name)
            self.__logger.debug(f'Found existing queue "{self.__work_queue_name}"')
        except self.__sqs.meta.client.exceptions.QueueDoesNotExist:
            self.__work_item_fifo = self.__sqs.create_queue(QueueName=self.__work_queue_name, Attributes={'DelaySeconds': '5', 'FifoQueue': 'true', 'VisibilityTimeout': '43200'})
            self.__logger.debug(f'Created new FIFO queue "{self.__work_queue_name}"')
        try:
            self.__status_queue = self.__sqs.get_queue_by_name(QueueName=self.__status_queue_name)
            self.__logger.debug(f'Found existing queue "{self.__status_queue_name}"')
        except self.__sqs.meta.client.exceptions.QueueDoesNotExist:
            self.__status_queue = self.__sqs.create_queue(QueueName=self.__status_queue_name, Attributes={'VisibilityTimeout': '10'})
            self.__logger.debug(f'Created new queue "{self.__status_queue_name}"')
    def queue_work_item(self, work_item_name: str) ->str:
        """Put a work item on the queue"""
        response = self.__work_item_fifo.send_message(MessageBody=work_item_name, MessageDeduplicationId=work_item_name, MessageGroupId='WorkItems')
        self.__logger.info(f'Queued work item "{work_item_name}" with MessageId: {response["MessageId"]}')
        return response['MessageId']
    def work_items(self):
        """Get the next work item from the queue"""
        while True:
            for message in self.__work_item_fifo.receive_messages(WaitTimeSeconds=20):
                if message is None:
                    self.__logger.debug('No work items currently available')
                    continue
                self.__logger.debug(f'Got work item "{message.body}"')
                self.__logger.debug(f'Deleting work item "{message.body}"')
                body = message.body
                # message.delete()
                yield body
    def report_status(self, work_item: str, message: str, is_fatal: bool = False, is_complete: bool = False):
        """Report a status message to the status queue"""
        body = {
            'work_item': work_item,
            'message': message,
            'is_fatal': is_fatal,
            'is_complete': is_complete,
            'hostname': self.__hostname
        }
        self.__logger.debug(f'Reporting status: "{message}", is_fatal: {is_fatal}')
        response = self.__status_queue.send_message(MessageBody=json.dumps(body))
        return response['MessageId']
    def status_items(self):
        """Get the next status item from the queue"""
        while True:
            for message in self.__status_queue.receive_messages(WaitTimeSeconds=20):
                if message is None:
                    self.__logger.debug('No status items currently available')
                    continue
                self.__logger.debug(f'Got status item "{message.body}"')
                yield json.loads(message.body)
                self.__logger.debug(f'Deleting status item "{message.body}"')
                message.delete()