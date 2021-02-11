__author__ = 'BinglanLi'

import os
import  vcf_preprocess_exceptions as Exceptions

def obtain_folder_path(path):
    return os.path.split(path)[0]

def obtain_vcf_file_prefix(path):
    vcf_file_full_name = os.path.split(path)[1].split('.')
    if (vcf_file_full_name[-2] == "vcf") and (vcf_file_full_name[-1] == "gz"):
        vcf_file_prefix = '.'.join(vcf_file_full_name[:len(vcf_file_full_name) - 2])
        return vcf_file_prefix
    else:
        raise Exceptions.InappropriateVCFSuffix(path)

def setOutputPath(folder):
    return os.getcwd() if not folder else folder

def AppendCompressedVCFSuffix(name):
    return name + ".vcf.gz"  


