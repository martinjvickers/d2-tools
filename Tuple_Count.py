#!/usr/bin/env python
import os,sys
import string
import optparse
import SeqSign
import copy

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-l", "--SampleFileList", action = "store", type = "string", dest = "samples_files_list",
                  help = "list of reads files of Samples")
parser.add_option("-k", "--kofKTuple", action = "store", type = "string", dest = "k_of_KTuple",
                  help = "the value k of KTuple")
parser.add_option("-t", "--tcp", action = "store", type = "string", dest = "tc_path", default='./',
                  help = "the folder to save the TupleCount files")

(options, args) = parser.parse_args()
if (options.samples_files_list is None or
    options.k_of_KTuple is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(0)
FileList=options.samples_files_list
k=int(options.k_of_KTuple)
TC_Path=options.tc_path
if TC_Path[-1] is not '/':
        TC_Path+='/'


bk=open(FileList)
iter=0
title_line='k-tuple'+'\t'+'number_of_occurrences'
for fileName in bk.readlines():
	iter=iter+1
        fileName=fileName[0:-1]
        temp1=fileName.split('/')
        temp2=temp1[-1].rsplit('.',1)
        sampleName=temp2[0]
        cv=SeqSign.countVector(fileName, k)
	outFileName=TC_Path+sampleName+'_k'+str(k)+'_tupleCount.txt'
	out_ktuple=open(outFileName, 'w')
	if iter==1:
                KTuple_list=cv.KTuple_count.keys()
	print >>out_ktuple, title_line
	#print >>out_ktuple, 'total'+'\t'+str(cv.total)
	for ktuple in KTuple_list:
		write_lines=ktuple+'\t'+str(cv.KTuple_count[ktuple])
		print >>out_ktuple, write_lines
	out_ktuple.close()
bk.close()
