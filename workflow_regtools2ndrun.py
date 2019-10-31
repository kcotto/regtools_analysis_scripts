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

        # sample_id_base = self.sample_id.rsplit('-', 1)[0]
        # self.vcf_sample_id = sample_id_base + '-01A'
        self.cohort = self.sample_id.split(';')[0]
        self.sample = self.sample_id.split(';')[1]

        # uid and gid are to (hopefully) fix the ownership of the files created by Docker
        self.uid = os.getuid()
        self.gid = os.getgid()

        # define the intermediate state variables we will be using later, to keep pylint happy...
        self.vcf_gz_files = None
        self.vcf_files = None
        self.bam_manifest_path = None
        self.vcf_manifest_path = None
        self.gdc_token = Path('gdc-user-token.txt')

        # Make sure the working directory exists
        if not self.filesystem_path.is_dir():
            self.filesystem_path.mkdir()

        # Make sure the working directory is empty
        # self.cleanup()

    def cleanup(self):
        """Clean up the working directory"""
        if self.filesystem_path.is_dir():
            shutil.rmtree(str(self.filesystem_path), ignore_errors=True)
        shutil.rmtree(self.sample, ignore_errors=True)
        os.remove(f'{self.sample}.tar.gz')

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
        subprocess.run(f'docker run --user {self.uid}:{self.gid} -v {self.cwd}:/manifests {image_name_and_all_args}',
                       shell=True, check=True, stdout=stdout)


    def get_gdc_token(self):
        self.logger.info('Getting GDC token file')
        if self.s3_token_download_url:
            subprocess.run(f'aws s3 cp {self.s3_token_download_url} ./{self.gdc_token}', shell=True, check=True)
        else:
            self.logger.warning('S3 Token Download URL is not set - new token not downloaded.')

        if not self.gdc_token.exists():
            raise RegtoolsWorkflowException(self.sample_id,
                                            'GDC Token file does not exist.  Unable to download manifest files')

        os.chmod(self.gdc_token, 0o600)
        with open(self.gdc_token, 'r') as myfile:
            raw_token = myfile.read() + '========'



    def run_regtools(self, identify_args: str, file_suffix: str):
        self.logger.info(
            f'About to run regtools griffithlab/regtools regtools cis-splice-effects associate {identify_args} /manifests/all_variants_sorted.vcf /manifests/{self.sample}/{self.sample}_extracted.bed /manifests/GRCh38.d1.vd1.fa /manifests/gencode.v29.annotation.gtf -o /manifests/{self.sample}/output/cse_identify_filtered_compare_{file_suffix}.tsv -v /manifests/{self.sample}/output/cse_identify_filtered_compare_{file_suffix}.vcf -j /manifests/{self.sample}/output/cse_identify_filtered_compare_{file_suffix}.bed')
        self.run_docker_image_as_current_user(
            f'griffithlab/regtools regtools cis-splice-effects associate {identify_args} /manifests/{self.sample}/all_variants_sorted.vcf /manifests/{self.sample}/{self.sample}_extracted.bed /manifests/GRCh38.d1.vd1.fa /manifests/gencode.v29.annotation.gtf -o /manifests/{self.sample}/output/cse_identify_filtered_compare_{file_suffix}.tsv -v /manifests/{self.sample}/output/cse_identify_filtered_compare_{file_suffix}.vcf -j /manifests/{self.sample}/output/cse_identify_filtered_compare_{file_suffix}.bed')

    def download_from_s3(self):
        self.logger.info('Downloading sample from s3')
        subprocess.run(f'aws s3 cp {self.s3_archive_upload_url}/{self.cohort}/{self.sample}.tar.gz {self.sample}.tar.gz', shell=True,
                       check=True)
        subprocess.run(f'tar xzf {self.sample}.tar.gz', shell=True, check=True)
        os.rename(filesystem_path, sample_id.split(';')[1])
        subprocess.run(
            f'aws s3 cp {self.s3_archive_upload_url}/{self.cohort}/all_variants_sorted.vcf {self.sample}/',
            shell=True,
            check=True)


    def identify_cis_splice_effects(self):
        self.logger.info('Running regtools cis-splice-effects identify')
        os.mkdir(f'{self.sample}/output')
        # default regtools run #
        self.run_regtools('', 'default')
        # i50e5 regtools run #
        self.run_regtools('-i 50 -e 5', 'i50e5')
        # E regtools run #
        self.run_regtools('-E', 'E')
        # I regtools run #
        self.run_regtools('-I', 'I')

    def check_files_and_save(self):
        self.logger.info('Archiving result and saving to S3')
        os.remove(f'{self.sample}/all_variants_sorted.vcf')
        subprocess.run(f'tar czf {self.sample}.tar.gz {self.sample}', shell=True, check=True)

        if self.s3_archive_upload_url:
            subprocess.run(f'aws s3 cp {self.sample}.tar.gz {self.s3_archive_upload_url}/{self.cohort}/{self.sample}.tar.gz', shell=True,
                           check=True)
        else:
            self.logger.error('S3 Archive Upload URL is not set - archive file not uploaded')


logger = logging.getLogger('regtools_reader')
logger.setLevel(logging.DEBUG)
log_formatter = logging.Formatter('%(asctime)s:%(levelname)7s:%(name)s:%(message)s')
log_formatter.default_msec_format = '%s.%03d'
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(log_formatter)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

sqs = SQSWrapper('default')

try:
    logger.info('Waiting for work item')
    # for sample_id in sqs.work_items():
    while True:
        sample_id = 'BLCA;TCGA-2F-A9KQ-01A'
        filesystem_path = 'regtools_wd_' + sample_id.split(';')[1]
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

            # grab the updated user token file from S3
            workflow.get_gdc_token()

            sqs.report_status(sample_id, f'Running regtools on sample {sample_id}')
            start_time = datetime.datetime.now()
            workflow.download_from_s3()
            workflow.identify_cis_splice_effects()
            end_time = datetime.datetime.now()
            elapsed_time = end_time - start_time
            logger.info(f'Elapsed time: {elapsed_time.total_seconds()} seconds')

            workflow.check_files_and_save()

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
            workflow.cleanup()
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
