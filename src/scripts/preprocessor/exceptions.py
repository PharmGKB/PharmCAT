__author__ = 'BinglanLi'

from pathlib import Path

from typing import Union


# partial credits to Heroico (GitHub User Name) and https://github.com/hakyimlab/MetaXcan


class ReportableException(Exception):
    """Simple exception with message(s)"""
    def __init__(self, msg: str):
        self.msg: str = msg

    def __str__(self):
        return self.msg


class InappropriateVCFSuffix(ReportableException):
    """Inappropriate Input VCF suffix (not ending with .vcf.gz)"""
    def __init__(self, msg: Union[Path, str]):
        super(InappropriateVCFSuffix, self).__init__("Inappropriate VCF suffix (not ending with '.vcf.gz'): %s"
                                                     % str(msg))


class InvalidURL(ReportableException):
    """Inappropriate URL. No downloadable content found."""
    def __init__(self, url: str):
        super(InvalidURL, self).__init__("Invalid URL: %s" % url)
