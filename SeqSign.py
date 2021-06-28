'''
Created on 2012-09-28

@author: liulin
'''

import os,sys
import string
import optparse
import math
##############################################dissimilarity functions#############################################
########################################################dAI#######################################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
pX is k-tuple possibility under Markov model of sample X
pY is k-tuple possibility under Markov model of sample Y
nX is total number of k-tuples of sample X
'''
def S2(cX, cY, pX, pY, nX, nY):
	S2kr=0.0;n=0
	tempX=0.0;tempY=0.0
	for ktuple in cX:
		n+=1
		if ktuple in cY:
			cXi=cX[ktuple];cYi=cY[ktuple]
			pXi=pX[ktuple];pYi=pY[ktuple]
			fXi=cXi/float(nX);fYi=cYi/float(nY)
			cXi_sigma=cXi*pXi;cYi_sigma=cYi*pYi
			if cXi_sigma==cYi_sigma:
				temp1=0;temp2=0
			else:
				temp=float(cXi_sigma+cYi_sigma)
				temp3=2*cXi_sigma/temp;temp1=cXi_sigma*math.log(temp3)
				temp4=2*cYi_sigma/temp;temp2=cYi_sigma*math.log(temp4)
			tempX+=temp1;tempY+=temp2
		else:
			print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!"
			sys.exit(1)
	if tempX==0 and tempY==0:
		S2kr=0.0
	else:
		S2kr=(tempX+tempY)/n+2*math.log(2.0)
	return S2kr
##################################################################################################################
#################################################Euclidean Distance###############################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
'''
def Eu(cX, cY, nX, nY):
        d=0.0
        D=0.0
	temp=0.0
        for ktuple in cX:
                if ktuple in cY:
                        fXi=cX[ktuple]/float(nX)
                        fYi=cY[ktuple]/float(nY)
			temp=(fXi-fYi)**(2)
                        D+=temp
                else:
                        print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!" % (ktuple)
                        sys.exit(1)
        d=D**(0.5)########sqrt
        return d
##################################################################################################################
#################################################Manhattan Distance###############################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
'''
def Ma(cX, cY, nX, nY):
        d=0.0
        D=0.0
        temp=0.0
        for ktuple in cX:
                if ktuple in cY:
                        fXi=cX[ktuple]/float(nX)
                      	fYi=cY[ktuple]/float(nY)
                        temp=abs(fXi-fYi)
                        D+=temp
                else:
                        print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!" % (ktuple)
                        sys.exit(1)
        d=D########sqrt
        return d
##################################################################################################################
#################################################Chebyshev Distance###############################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
'''
def Ch(cX, cY, nX, nY):
        d=0.0
        D=0.0
        temp=0.0
        for ktuple in cX:
                if ktuple in cY:
                        fXi=cX[ktuple]/float(nX)
                        fYi=cY[ktuple]/float(nY)
                        temp=abs(fXi-fYi)
                        if temp>D:
				D=temp
                else:
                        print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!" % (ktuple)
                        sys.exit(1)
        d=D
        return d
##################################################################################################################

#########################################################d2#######################################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
'''
def d2(cX, cY):
	d_2=0.0
	D_2=0.0
	tempX=0.0
	tempY=0.0
	for ktuple in cX:
		if ktuple in cY:
			cXi=cX[ktuple]
			cYi=cY[ktuple]
			D_2+=cXi*cYi
			tempX+=cXi*cXi
			tempY+=cYi*cYi
		else:
			print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!" % (ktuple)
			sys.exit(1)
	tempX=tempX**(0.5)########sqrt
	tempY=tempY**(0.5)
	temp=D_2/(tempX*tempY)
	d_2=0.5*(1-temp)
	return d_2
##################################################################################################################
########################################################dS2#######################################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
pX is k-tuple possibility under Markov model of sample X
pY is k-tuple possibility under Markov model of sample Y
nX is total number of k-tuples of sample X
'''
def d2S(cX, cY, pX, pY, nX, nY):
	d_2S=0.0
	D_2S=0.0
	tempX=0.0
	tempY=0.0
	n=0
	for ktuple in cX:
		n=n+1
		if ktuple in cY:
			cXi=cX[ktuple]
			cYi=cY[ktuple]
			pXi=pX[ktuple]
			pYi=pY[ktuple]
			cXi_bar=cXi-nX*pXi
			cYi_bar=cYi-nY*pYi
			temp1=cXi_bar**2
			temp2=cYi_bar**2
			temp3=(temp1+temp2)**0.5
			if temp3==0:
				temp3=1.0	
			D_2S+=cXi_bar*cYi_bar/temp3
			tempX+=cXi_bar*cXi_bar/temp3
			tempY+=cYi_bar*cYi_bar/temp3
		else:
			print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!"
			sys.exit(1)
	tempX=tempX**0.5
	tempY=tempY**0.5
	temp=D_2S/(tempX*tempY)
	d_2S=0.5*(1-temp)
	return d_2S
##################################################################################################################
######################################################dStar2######################################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
pX is k-tuple possibility under Markov model of sample X
pY is k-tuple possibility under Markov model of sample Y
nX is total number of k-tuples of sample X
'''
def d2Star(cX, cY, pX, pY, nX, nY):
	d_2Star=0.0
	D_2Star=0.0
	tempX=0.0
	tempY=0.0
	for ktuple in cX:
		if ktuple in cY:
			cXi=cX[ktuple]
			cYi=cY[ktuple]
			pXi=pX[ktuple]
			pYi=pY[ktuple]
			temp1=nX*pXi
			temp2=nY*pYi
			cXi_bar=cXi-temp1
			cYi_bar=cYi-temp2
			temp3=(temp1*temp2)**0.5
			#if temp3==0
			if temp1==0:
				temp1=1.0
                       		temp3=1.0
			if temp2==0:
				temp2=1.0
				temp3=1.0
			D_2Star+=cXi_bar*cYi_bar/temp3
			tempX+=cXi_bar*cXi_bar/temp1
			tempY+=cYi_bar*cYi_bar/temp2
		else:
			print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!"
			sys.exit(1)
	tempX=tempX**0.5
	tempY=tempY**0.5
	temp=D_2Star/(tempX*tempY)
	d_2Star=0.5*(1-temp)
	return d_2Star
##################################################################################################################
########################################################Hao#######################################################
'''
cX is k-tuple count vector of sample X
cY is k-tuple count vector of sample Y
pX is k-tuple possibility under Markov model of sample X
pY is k-tuple possibility under Markov model of sample Y
nX is total number of k-tuples of sample X
'''
def Hao(cX, cY, pX, pY, nX, nY):
	d_Hao=0.0
	tempXY=0.0
	tempX=0.0
	tempY=0.0
	for ktuple in cX:
		if ktuple in cY:
			cXi=cX[ktuple]
			cYi=cY[ktuple]
			pXi=pX[ktuple]
			pYi=pY[ktuple]
			fXi=float(cXi)/nX
			fYi=float(cYi)/nY
			if pXi==0:
				temp1=-1
			else:
				temp1=(fXi/pXi)-1.0
			if pYi==0:
				temp2=-1
			else:
				temp2=(fYi/pYi)-1.0
			tempXY+=temp1*temp2
			tempX+=temp1*temp1
			tempY+=temp2*temp2
		else:
			print >>sys.stderr, "Error: the k-tuple %s in sampleX, but not in samleY!"
			sys.exit(1)
	tempX=tempX**0.5
	tempY=tempY**0.5
	temp=tempXY/(tempX*tempY)
	d_Hao=(1-temp)/2
	return d_Hao
##################################################################################################################
#####################################class countVector: caculate count vector#####################################
class countVector:
	####################initial function####################
	def __init__(self, readsFile, k, down=0):
		self.fileName = str(readsFile)
		self.kbp=int(k)
		self.KTuple_count = {}
		self.total=0
		self.ACGT_pair={}
                self.ACGT_pair['A']='T'
                self.ACGT_pair['C']='G'
                self.ACGT_pair['G']='C'
                self.ACGT_pair['T']='A'
		self.ACGT_pair['N']='N'
		self.seq_range=int(down)#sequence range
		#self.Get_Count_Vector()
		if k==1:
			self.Get_Count_Vector_one()
		else:
			self.Get_Count_Vector()	
			
	#########define the function to get all k-tuples########
	def Get_KTuple_sets(self, position, S_ktuple, diction):
		if position < 0:
			key_temp = ''.join(S_ktuple)
			diction[key_temp] = 0
			return
		else:
			S_ktuple[position] = 'A'
			self.Get_KTuple_sets(position-1,S_ktuple,diction)
			S_ktuple[position] = 'C'
			self.Get_KTuple_sets(position-1,S_ktuple,diction)
			S_ktuple[position] = 'G'
			self.Get_KTuple_sets(position-1,S_ktuple,diction)
			S_ktuple[position] = 'T'
			self.Get_KTuple_sets(position-1,S_ktuple,diction)
			return
	def Get_reads_complements(self, read_seq, n):
		read_complement=''
		read_reverse=read_seq[::-1]
		for w in read_reverse:
			try:
				read_complement+=self.ACGT_pair[w]
			except:
				print n
				print read_seq
		return read_complement
			 			
	##############read files to caculate the k-tuple count vector for k > 1#########
	def Get_Count_Vector(self):
		ktemp=self.kbp
		sr=self.seq_range
		tempSep=self.fileName.split('.')
		fileType=tempSep[-1]
		if fileType=='fastq' or fileType=='fq':
			startChar='@'
			mod_n=4
			tar_n=2
		elif fileType=='fna' or fileType=='fa' or fileType=='fasta':
			startChar='>'
			mod_n=2
                        tar_n=0
		else:
			print "File format is wrong! Please input FA or FQ files."
			sys.exit(1)
		bk=open(self.fileName)
		KTuple_intialStr=['N']*ktemp
		self.Get_KTuple_sets(ktemp-1, KTuple_intialStr, self.KTuple_count)
		line_count=0
		n=0
		for seq in bk.xreadlines():
			line_count+=1
			n+=1
			if seq[0]==startChar:
				line_count=1
			if (line_count%mod_n) != tar_n:
				continue
			seq=seq.upper()
			seq_forward=seq[0:-1-sr]
			read_seq_in=seq[sr:-1]	
			seq_complement=self.Get_reads_complements(read_seq_in, n)
			length=len(seq_forward)
			##########calculate the k-tuple count vector################
			for st in range(length-ktemp+1):
				ktuple=seq_forward[st:st+ktemp]
				if ktuple in self.KTuple_count:
					self.KTuple_count[ktuple]+=1
					self.total+=1
				ktuple=seq_complement[st:st+ktemp]
				if ktuple in self.KTuple_count:
					self.KTuple_count[ktuple]+=1
					self.total+=1
		bk.close()
		
	##############read files to caculate the k-tuple count vector for k = 1#########		
	def Get_Count_Vector_one(self):
		if self.kbp != 1:
			print "k is not equal to 1!"
			return
		tempSep=self.fileName.split('.')
                fileType=tempSep[-1]
                if fileType=='fastq' or fileType=='fq':
			startChar='@'
                        mod_n=4
                        tar_n=2
                elif fileType=='fna' or fileType=='fa' or fileType=='fasta':
			startChar='>'
                        mod_n=2
                        tar_n=0
                else:
                        print "File format is wrong! Please input FA or FQ files."
                        sys.exit(1)
		bk = open(self.fileName)
		line_count=0
		self.KTuple_count['A']=0
		self.KTuple_count['C']=0
		self.KTuple_count['G']=0
		self.KTuple_count['T']=0
		for seq in bk.xreadlines():
			line_count+=1
			if seq[0]==startChar:
                                line_count=1
			if (line_count%mod_n) != tar_n:
				continue
			seq=seq.upper()
			seq_forward=seq[0:-1]
			A_count=seq_forward.count('A')
			C_count=seq_forward.count('C')
			G_count=seq_forward.count('G')
			T_count=seq_forward.count('T')
			temp=A_count+C_count+G_count+T_count
			self.KTuple_count['A']+=A_count
                        self.KTuple_count['C']+=C_count
                        self.KTuple_count['G']+=G_count
                        self.KTuple_count['T']+=T_count
			self.total+=temp
			if self.seq_range==1:
				to_del=seq_forward[0]
				if to_del in self.KTuple_count:
					self.KTuple_count[to_del]-=1
					self.total-=1
				to_del=self.ACGT_pair[seq_forward[-1]]
				if to_del in self.KTuple_count:
					self.KTuple_count[to_del]-=1
					self.total-=1
		bk.close()
		
	##############function: write the vector to a file#########
	def write_ktuple_to_file(self, outPrefix):
		OUT=open(outPrefix+".k-tuple_count_vector.txt", 'w')
		for key in self.KTuple_count:
			print >>OUT, str(key)+'\t'+str(self.KTuple_count[key])
		OUT.close()
#####################################################class end####################################################

###################class MarkovPossibility: calculate k-tuple possibility under markov modle######################
class MarkovPossibility:
	def __init__(self, total, r_count, r1_count, ktuple_count, r_order):
                self.total = int(total)
                self.order = int(r_order)
		if r_order!=0:
			self.r_count = r_count
		self.r1_count = r1_count
                self.KTuple_count = ktuple_count
                self.KTuple_possibility = {}
                self.Get_Markov_Possibility()
	def Get_Markov_Possibility(self):
		r=self.order
		if r==0:
			f_base={}
			f_base['A']=self.r1_count['A']/float(self.total)
                        f_base['C']=self.r1_count['C']/float(self.total)
                        f_base['G']=self.r1_count['G']/float(self.total)
                        f_base['T']=self.r1_count['T']/float(self.total)
			for ktuple in self.KTuple_count:
                                p_temp=1
                                for w in ktuple:
                                        p_temp*=f_base[w]
				self.KTuple_possibility[ktuple]=format(p_temp,'e')
                                #self.KTuple_possibility[ktuple]=round(p_temp, 5)
                        del f_base
		else:
                        for ktuple in self.KTuple_count:
				k=len(ktuple)
				if k<=r:
					self.KTuple_possibility[ktuple]='NA'
					continue
                                w_start=ktuple[0:r]
                                p_temp=self.r_count[w_start]/float(self.total)
                                for j in range(k-r):
                                        w1=ktuple[j:j+r+1]
					w2=ktuple[j:j+r]
					w2A=''.join([w2,'A'])
					w2C=''.join([w2,'C'])
					w2G=''.join([w2,'G'])
					w2T=''.join([w2,'T'])
                                        nk1=self.r1_count[w1]
                                        nk2=self.r1_count[w2A]+self.r1_count[w2C]+self.r1_count[w2G]+self.r1_count[w2T]
                                        if nk2==0:
                                                p_temp=0
						self.KTuple_possibility[ktuple]=format(p_temp,'e')
                                                continue
                                        else:
                                                p_temp*=nk1/float(nk2)
                               			self.KTuple_possibility[ktuple]=format(p_temp,'e')
					#self.KTuple_possibility[ktuple]=round(p_temp, 5)

	def Get_Markov_Possibility_d2Meta(self):
		inputFile=self.fileName
		#####M0,r=0#####################
		obj=countVector(inputFile,1)
                obj.KTuple_count['A']=obj.KTuple_count['A']/float(obj.total)
                obj.KTuple_count['C']=obj.KTuple_count['C']/float(obj.total)
                obj.KTuple_count['G']=obj.KTuple_count['G']/float(obj.total)
                obj.KTuple_count['T']=obj.KTuple_count['T']/float(obj.total)
		###M1,M2,M3,r=1,2,3######################
		obj1=countVector(inputFile, 1, 1)
                obj2=countVector(inputFile, 2, 1)
		obj3=countVector(inputFile, 3, 1)
		obj4=countVector(inputFile, 4, 1)
		for ktuple in self.KTuple_count:
			k=len(ktuple)
			#######################M0#######################
                        p_temp=1
                        for w in ktuple:
                                p_temp*=obj.KTuple_count[w]
			self.KTuple_possibility[ktuple]=format(p_temp,'e')
                        #self.M0_possibility[ktuple]=round(p_temp,5)
			#######################M1#######################
			if k<=1:
				self.M1_possibility[ktuple]='NA'
			else:
				w_start=ktuple[0:1]
				p_temp=obj1.KTuple_count[w_start]/float(obj1.total)
				for j in range(k-1):
					w1=ktuple[j:j+1]
					w2=ktuple[j:j+2]
					nk1=obj2.KTuple_count[w2]
					nk2=obj1.KTuple_count[w1]
					if nk2==0:
						p_temp=0
						continue
					else:
						p_temp*=nk1/float(nk2)
				self.KTuple_possibility[ktuple]=format(p_temp,'e')
				#self.M1_possibility[ktuple]=round(p_temp,5)
			#######################M2#######################
			if k<=2:
				self.M2_possibility[ktuple]='NA'
			else:
				w_start=ktuple[0:2]
	                        p_temp=obj2.KTuple_count[w_start]/float(obj2.total)
        	                for j in range(k-2):
                	                w1=ktuple[j:j+2]
                        	        w2=ktuple[j:j+3]
                                	nk1=obj3.KTuple_count[w2]
                               		nk2=obj2.KTuple_count[w1]
                             	   	if nk2==0:
                                        	p_temp=0
                                        	continue
                                	else:
                                        	p_temp*=nk1/float(nk2)
				self.KTuple_possibility[ktuple]=format(p_temp,'e')
                        	#self.M2_possibility[ktuple]=round(p_temp,5)
			#######################M3#######################
			if k<=3:
                                self.M3_possibility[ktuple]='NA'
			else:
                        	w_start=ktuple[0:3]
                        	p_temp=obj3.KTuple_count[w_start]/float(obj3.total)
                        	for j in range(k-3):
                                	w1=ktuple[j:j+3]
                                	w2=ktuple[j:j+4]
                                	nk1=obj4.KTuple_count[w2]
                                	nk2=obj3.KTuple_count[w1]
                                	if nk2==0:
                                        	p_temp=0
                                        	continue
                                	else:
                                        	p_temp*=nk1/float(nk2)
				self.KTuple_possibility[ktuple]=format(p_temp,'e')
                        	#self.M3_possibility[ktuple]=round(p_temp,5)
		############################################################################	
                del obj
		del obj1
		del obj2
		del obj3
		del obj4
#######################main function##################
