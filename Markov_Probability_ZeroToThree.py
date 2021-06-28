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
parser.add_option("-m", "--mpp", action = "store", type = "string", dest = "mp_path", default='./',
                  help = "the folder to save the Markov Provility(order is from 0~3) files")

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
MP_Path=options.mp_path
if MP_Path[-1] is not '/':
        MP_Path+='/'
bk=open(FileList)
iter=0
title_line='k-tuple'+'\t'+'number_of_occurrences'+'\t'+'M0_possibility'+'\t'+'M1_possibility'+'\t'+'M2_possibility'+'\t'+'M3_possibility'
for fileName in bk.readlines():
	iter=iter+1
        fileName=fileName[0:-1]
        temp1=fileName.split('/')
        temp2=temp1[-1].rsplit('.',1)
        sampleName=temp2[0]
	KTuple_count={}
	r_count_dict={}
	r_total_dict={}
	for i in range(1,5):
		tempFileName=TC_Path+sampleName+'_k'+str(i)+'_tupleCount.txt'
		temp_file=open(tempFileName, 'r')
		total=0
		r_count_dict[i]={}
		line=temp_file.readline()
		for lines in temp_file.readlines():
                	fields=lines[0:-1].split('\t')
                	ktuple=fields[0]
                	count=int(fields[1])
                	r_count_dict[i][ktuple]=count
                	total+=count
		r_total_dict[i]=total
	#print r_count_dict.keys()
	#print r_total_dict.keys()
	tupleFileName=TC_Path+sampleName+'_k'+str(k)+'_tupleCount.txt'
	tupleFile=open(tupleFileName, 'r')
	line=tupleFile.readline()
	for lines in tupleFile.readlines():
		fields=lines[0:-1].split('\t')
		ktuple=fields[0]
		count=int(fields[1])
		KTuple_count[ktuple]=count
	mp0=SeqSign.MarkovPossibility(r_total_dict[1], r_count_dict[1], r_count_dict[1], KTuple_count, 0)
	mp1=SeqSign.MarkovPossibility(r_total_dict[1], r_count_dict[1], r_count_dict[2], KTuple_count, 1)
	mp2=SeqSign.MarkovPossibility(r_total_dict[2], r_count_dict[2], r_count_dict[3], KTuple_count, 2)
	mp3=SeqSign.MarkovPossibility(r_total_dict[3], r_count_dict[3], r_count_dict[4], KTuple_count, 3)
	outFileName=MP_Path+sampleName+'_k'+str(k)+'_wordcount_pw.txt'
	out_ktuple=open(outFileName, 'w')
	if iter==1:
                KTuple_list=KTuple_count.keys()
	print >>out_ktuple, title_line
	for ktuple in KTuple_list:
		n=KTuple_count[ktuple]
		if n==0:
			write_lines="%s\t%E" %(ktuple,1e-10)
		else:
			write_lines=ktuple+'\t'+str(KTuple_count[ktuple])
		write_lines+="\t%s" % (mp0.KTuple_possibility[ktuple])
		write_lines+="\t%s" % (mp1.KTuple_possibility[ktuple])
		write_lines+="\t%s" % (mp2.KTuple_possibility[ktuple])
		write_lines+="\t%s" % (mp3.KTuple_possibility[ktuple])
		print >>out_ktuple, write_lines
	out_ktuple.close()
bk.close()
