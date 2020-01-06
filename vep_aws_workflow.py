"""Process a bunch of stuff"""
import subprocess
from pathlib import Path
import os
import sys
import csv
import gzip
import time
import shutil
import tarfile
import logging
import datetime
import traceback
import requests
import fnmatch

from sqs_wrapper import SQSWrapper


class RegtoolsWorkflowException(Exception):
    """Exception raised when the workflow detects an error"""

    def __init__(self, sample_id: str = None, message: str = None):
        """Initialize a RegtoolsWorkflowException"""
        super().__init__()
        self.message = message
        self.sample_id = sample_id


class RegtoolsWorkflow:
    """Container class for the regtools workflow"""

    def __init__(self, sample_id: str = None, filesystem_path: str = None, logger=None,
                 s3_token_download_url: str = None, s3_archive_upload_url: str = None):
        if sample_id is None or filesystem_path is None or logger is None:
            raise RegtoolsWorkflowException(
                'RegtoolsWorkflow.__init__(): sample_id, filesystem_path, and logger must ALL be specified')
        logger.info(f'Initializing RegtoolsWorkflow. sample_id: "{sample_id}", filesystem_path: "{filesystem_path}"')
        self.sample_id = sample_id
        self.filesystem_path = Path(filesystem_path)
        self.logger = logger
        self.cwd = os.getcwd()
        self.s3_token_download_url = s3_token_download_url
        self.s3_archive_upload_url = s3_archive_upload_url

        # uid and gid are to (hopefully) fix the ownership of the files created by Docker
        self.uid = os.getuid()
        self.gid = os.getgid()

        self.cohort = self.sample_id.split(';')[0]
        self.sample = self.sample_id.split(';')[1]

        # Make sure the working directory exists
        if not self.filesystem_path.is_dir():
            self.filesystem_path.mkdir()

        # Make sure the working directory is empty
        # self.cleanup()

    def cleanup(self):
        """Clean up the working directory"""
        if self.filesystem_path.is_dir():
            shutil.rmtree(str(self.filesystem_path), ignore_errors=True)

    def pull_docker_image(self, image_name: str):
        """Pull a docker image with retries on failure"""
        self.logger.info(f'Pulling Docker image {image_name}')
        for retry_number in range(5):
            try:
                subprocess.run(f'docker pull {image_name}', shell=True, check=True)
                return
            except subprocess.CalledProcessError as e:
                self.logger.exception(
                    f'An error occurred while initializing the workflow environment: Command {e.cmd} returned non-zero exit status {e.returncode}')
                time.sleep(2 * (retry_number + 1))
        raise RegtoolsWorkflowException('initialization',
                                        f'Retry count exceeded while pulling Docker image {image_name}')

    def run_docker_image_as_current_user(self, image_name_and_all_args: str, stdout=None):
        docker_args = image_name_and_all_args.split()
        self.logger.info(f'Running docker container command "{docker_args[0]} {docker_args[1]}"')
        subprocess.run(f'docker run --user {self.uid}:{self.gid} -v /vep_data:/opt/vep/.vep:Z {image_name_and_all_args}',
                       shell=True, check=True, stdout=stdout)

    def download_from_s3_and_run_vep(self):
        self.logger.info('Downloading sample from s3')
        subprocess.run(
            f'aws s3 cp {self.s3_archive_upload_url}/{self.cohort}/{self.sample}.tar.gz {self.sample}.tar.gz',
            shell=True,
            check=True)
        subprocess.run(f'tar xzf {self.sample}.tar.gz {self.sample}/{self.sample}_master.vcf', shell=True, check=True)
        os.remove(f'{self.sample}.tar.gz')
        subprocess.run(f'mv {self.sample}/*.vcf vep_data/', shell=True, check=True)
        self.run_docker_image_as_current_user(
            f'ensemblorg/ensembl-vep ./vep --input_file=/opt/vep/.vep/{self.sample_id}.vcf --output_file=/opt/vep/.vep/{self.sample_id}_master.vep.vcf --vcf --everything --cache --dir_cache=/opt/vep/.vep/vep/  --force_overwrite --per_gene --format=vcf --assembly=GRCh38 --offline  --fasta=/opt/vep/.vep/GRCh38.d1.vd1.fa')
        subprocess.run(
            f'aws s3 cp vep_data/{self.sample_id}_master.vep.vcf {self.s3_archive_upload_url}/VEP_vcfs',
            shell=True,
            check=True)
        shutil.rmtree(f'{self.sample}/')
        files_to_remove_from_vep = glob.glob('vep_data/*.vcf')
        for file in files_to_remove_from_vep:
            os.remove(file)


logger = logging.getLogger('regtools_reader')
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter('%(asctime)s:%(levelname)7s:%(name)s:%(message)s')
log_formatter.default_msec_format = '%s.%03d'
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_formatter)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

sqs = SQSWrapper('default', 'vep')

try:
    logger.info('Waiting for work item')
    for sample_id in sqs.work_items():
        filesystem_path = 'regtools_wd_' + sample_id
        try:
            workflow = RegtoolsWorkflow(sample_id=sample_id,
                                        filesystem_path=filesystem_path,
                                        logger=logger,
                                        s3_token_download_url='s3://regtools-cwl-sharedfiles/gdc-user-token.txt',
                                        s3_archive_upload_url='s3://regtools-results-unstranded')

            file_handler = logging.FileHandler(f'{filesystem_path}/{sample_id}.log')
            file_handler.setFormatter(log_formatter)
            file_handler.setLevel(logging.DEBUG)
            logger.addHandler(file_handler)
            logger.info(f'Starting sample "{sample_id}"')
            sqs.report_status(sample_id, f'Starting sample {sample_id}')

            workflow.download_from_s3_and_run_vep()
            sqs.report_status(sample_id, f'Downloaded sample {sample_id}')


            start_time = datetime.datetime.now()
            workflow.download_from_s3_and_run_vep()
            sqs.report_status(sample_id, f'Downloaded sample {sample_id}')
            end_time = datetime.datetime.now()
            elapsed_time = end_time - start_time
            logger.info(f'Elapsed time: {elapsed_time.total_seconds()} seconds')


            logger.info(f'Completed sample "{sample_id}"')
            sqs.report_status(sample_id, f'Completed sample {sample_id}', is_complete=True)
        except RegtoolsWorkflowException as e:
            logger.exception(f'An error occurred while processing sample "{sample_id}"')
            sqs.report_status(e.sample_id, e.message, is_fatal=True)
        except subprocess.CalledProcessError as e:
            logger.exception(f'An error occurred while processing sample "{sample_id}"')
            sqs.report_status(sample_id, f'Command {e.cmd} returned non-zero exit status {e.returncode}', is_fatal=True)
        except:
            logger.exception(f'An unexpected exception occurred while processing sample "{sample_id}"')
        finally:
            file_handler.close()
            logger.removeHandler(file_handler)
            logger.info('Cleaning up scratch directory')
            # workflow.cleanup()
            logger.info(f'Closed log file "{sample_id}.log"')
            logger.info('Waiting for work item')
except KeyboardInterrupt:
    logger.info('Exiting at user request')
    sqs.report_status('exiting', 'Exiting at user request')
except:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    tb = traceback.format_exception(exc_type, exc_value, exc_traceback)
    logger.exception('Unhandled exception')
    sqs.report_status('unhandled_exception', ''.join(tb), is_fatal=True)
