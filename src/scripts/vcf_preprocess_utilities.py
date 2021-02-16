#! /usr/bin/env python
__author__ = 'BinglanLi'

import os
import shutil
import urllib.request as request
from contextlib import closing
import  vcf_preprocess_exceptions as Exceptions


def obtain_vcf_file_prefix(path):
    vcf_file_full_name = os.path.split(path)[1].split('.')
    if (vcf_file_full_name[-2] == "vcf") and (vcf_file_full_name[-1] == "gz"):
        vcf_file_prefix = '.'.join(vcf_file_full_name[:len(vcf_file_full_name) - 2])
        return vcf_file_prefix
    else:
        raise Exceptions.InappropriateVCFSuffix(path)

def download_from_url(url, download_to_dir, save_to_file = None):
    local_filename = os.path.join(download_to_dir, url.split('/')[-1]) if not save_to_file else save_to_file
    with closing(request.urlopen(url)) as r:
        with open(local_filename, 'wb') as f:
            print('Downloading %s' %local_filename)
            shutil.copyfileobj(r, f)
    return local_filename
