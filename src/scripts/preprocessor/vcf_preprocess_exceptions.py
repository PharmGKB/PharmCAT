__author__ = 'BinglanLi'
# partial credits to Heroico (GitHub User Name) and https://github.com/hakyimlab/MetaXcan

class ReportableException(Exception):
    """Simple exeception with message(s)"""
    def __init__(self, msg):
        self.msg = msg

class InappropriateVCFSuffix(ReportableException):
    """Inappropriate Input VCF suffix (not ending with .vcf.gz)"""
    def __init__(self, msg):
        super(InappropriateVCFSuffix, self).__init__("Inappropriate VCF suffix (not ending with '.vcf.gz'): %s" % (msg))

class InvalidURL(ReportableException):
    """Inappropriate URL. No downloadable content found."""
    def __init__(self, msg):
        super(InvalidURL, self).__init__("Invalid downloading URL: %s" % (msg))
