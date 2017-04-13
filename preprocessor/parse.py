import sys
import json
from pprint import pprint


def run(json_file):
    with open(json_file) as data_file:
        data = json.load(data_file)
    # pprint(data)
    summary = data['metadata']['inputFilename']
    for result in data['results']:
        matches = ""
        for diplotyope in result['diplotypes']:
            matches = "%s %s:%s" %(matches, diplotyope['name'], diplotyope['score'])
        summary = "%s %s:%s\t" % (summary, result['gene'], matches)
        print "PHARMCAT:\t", data['metadata']['inputFilename'], "\t", result['gene'], "\t", matches
    print "PHARMCAT SUMMARY: ", summary
    pass


if __name__ == '__main__':
    input_json = sys.argv[1]
    # print 'Running on:', input_json
    run(input_json)
