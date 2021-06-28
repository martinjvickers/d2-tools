#!/usr/bin/env python
import os,sys
import string
import optparse
import SeqSign
import math
import copy

prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-l", "--SampleFileList", action = "store", type = "string", dest = "samples_files_list",
                  help = "list of reads files of Samples")
parser.add_option("-k", "--kofKTuple", action = "store", type = "string", dest = "k_of_KTuple",
                  help = "the value k of KTuple")
parser.add_option("-r", "--order", action = "store", type = "string", dest = "Markov_model_order",default=0,
                  help = "the order of markov model(0,1,2,3,k-2)")
parser.add_option("-d", "--distance", action = "store", type = "string", dest = "function_of_distance",
                  help = "the function of distance calculating")                                
parser.add_option("-m", "--mpp", action = "store", type = "string", dest = "mp_path", default='./',
                  help = "the folder to save the Markov Provility(order is from 0~3) files") 
parser.add_option("-o", "--output", action = "store", type = "string", dest = "output_prefix",
                  help = "prefix of output files")

(options, args) = parser.parse_args()
if (options.samples_files_list is None or
    options.k_of_KTuple is None or
    options.Markov_model_order is None or
    options.output_prefix is None):
    print prog_base + ": error: missing required command-line argument."
    parser.print_help()
    sys.exit(0)
fd_list=['Eu','Ma','Ch','d2', 'd2S', 'd2Star','Hao','S2']    
if options.function_of_distance is None:
	fd='d2S'
else:
	fd=options.function_of_distance
	if fd not in fd_list:
		print "Error: distance calculating function should be Eu, Ma, Ch, d2, d2S, d2Star"
		sys.exit(1)

FileList=options.samples_files_list
k=int(options.k_of_KTuple)
r=int(options.Markov_model_order)
outPrefix=options.output_prefix
MP_Path=options.mp_path
if MP_Path[-1] is not '/':
        MP_Path+='/'

if r>3 and r<0:
        print "Error: for d2Meta and this version, the markov model order just can be 0,1,2,3."
        sys.exit(1)
if r >= k:
        print "Error: the k-tuple size must be more than markov model order."
        sys.exit(1)

bk=open(FileList)
sample_file_list=[]
sample_name_list=[]
sample_count={}
sample_possibilty={}
sample_total={}
for fileName in bk.readlines():
	fileName=fileName[0:-1]
	#ample_file_list.append(fileName)
	temp1=fileName.split('/')
	temp2=temp1[-1].rsplit('.', 1)
	sampleName=temp2[0]
	sample_name_list.append(sampleName)	
	kfileName=MP_Path+sampleName+'_k'+str(k)+'_wordcount_pw.txt'
	#rfileName=sampleName+'.k'+str(r)+'_ktuple_matrix.txt'
	#r1fileName=sampleName+'.k'+str(r+1)+'_ktuple_matrix.txt'
	try:
		kfile=open(kfileName)
	except:
		print 'File: '+kfileName+' does not exist!'
		sys.exit(1)
	KTuple_count={}
	KTuple_possibility={}
        total=0
        lines=kfile.readline()
        for lines in kfile.readlines():
           	fields=lines[0:-1].split('\t')
               	ktuple=fields[0]
                count=eval(fields[1])
		probability=eval(fields[2+r])
                KTuple_count[ktuple]=count
		KTuple_possibility[ktuple]=probability
                total+=count
	kfile.close()
	#####################################################################
	#mp=SeqSign.MarkovPossibility(fileName,KTuple_count,r)
	sample_count[sampleName]=copy.deepcopy(KTuple_count)
	del KTuple_count
	sample_possibilty[sampleName]=copy.deepcopy(KTuple_possibility)
	del KTuple_possibility
	sample_total[sampleName]=total
bk.close()
OUT_matrix=open(outPrefix+'.dissimilarity_matrix.txt', 'w')
title_line=''
for sampleName in sample_name_list:
	title_line+='\t'
	title_line+=sampleName
print >>OUT_matrix, title_line

nsample=len(sample_name_list)
for i in range(nsample):
	sampleName=sample_name_list[i]
	write_line=sampleName
	sampleX=sample_name_list[i]
	cX=sample_count[sampleX]
	pX=sample_possibilty[sampleX]
	nX=sample_total[sampleX]
	for j in range(0,i):
		sampleY=sample_name_list[j]
		cY=sample_count[sampleY]
		pY=sample_possibilty[sampleY]
		nY=sample_total[sampleY]
		if fd=='d2':
        		dist=SeqSign.d2(cX, cY)
		elif fd=='Eu':
			dist=SeqSign.Eu(cX, cY, nX, nY)
		elif fd=='Ma':
                        dist=SeqSign.Ma(cX, cY, nX, nY)
		elif fd=='Ch':
                        dist=SeqSign.Ch(cX, cY, nX, nY)
		elif fd=='d2S':
        		dist=SeqSign.d2S(cX, cY, pX, pY, nX, nY)
		elif fd=='d2Star':
        		dist=SeqSign.d2Star(cX, cY, pX, pY, nX, nY)
		elif fd=='Hao':
        		dist=SeqSign.Hao(cX, cY, pX, pY, nX, nY)
		elif fd=='S2':
                        dist=SeqSign.S2(cX, cY, pX, pY, nX, nY)
		write_line+='\t'+str(round(dist, 5))
	for j in range(i,nsample):
                write_line+='\t'+'0'
	print >>OUT_matrix, write_line
del sample_count
del sample_possibilty
del sample_total
OUT_matrix.close()
'''
cv1=SeqSign.countVector(readsFileName1, k)
mp1=SeqSign.MarkovPossibility(readsFileName1, cX, r)
cv2=SeqSign.countVector(readsFileName2, k)
mp2=SeqSign.MarkovPossibility(readsFileName2, cY, r)

if fd=='d2':
	dist=SeqSign.d2(cX, cY)
elif fd=='dS2':	
	dist=SeqSign.dS2(cX, cY, pX, pY, nX, nY)
elif fd=='dStar2':
	dist=SeqSign.dStar2(cX, cY, pX, pY, nX, nY)
elif fd=='Hao':
	dist=SeqSign.Hao(cX, cY, pX, pY, nX, nY)
print dist
'''
##############################################DONE##################################################################################	
