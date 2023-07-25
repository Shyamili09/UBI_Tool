#!/user/bin/python
#from Bio.SeqUtils.ProtParam import ProteinAnalysis
#>>> X = ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT")
#>>> print(X.count_amino_acids()['A']) 
import pprint
import re
import sys
import csv
import numpy
import operator
import os
import fnmatch
import subprocess
global ssplist
###################################   AA composition  ##############################################
def weightage(lys,kmerAAfreq):
	neweachAA=[]
	inpwei = {}
	with open('weightage_for_AA13mer.txt', 'rb') as csv_file:
		next(csv_file)
    		for row in csv.reader(csv_file, delimiter='\t'):
        		inpwei[row[0]] = row[1:]
	for AAfli in kmerAAfreq[1:]:
		for eachAA in AAfli:  #list with AA and freq
			for key in inpwei:
				if eachAA[0] == key:
					maxx=inpwei[key][0]
					wei=float(0.000000)
					limits=numpy.linspace(0,float(maxx)+3.5,5)
					if int(inpwei[key][2]) == 1:			#enriched AA
						if 0.07 <= float(eachAA[1]) < limits[1]:	#freq very high within range
							wei = 1*float(inpwei[key][3])	#weightage multiplied by 4 (high) and by weightage factor based on p value
						elif  limits[1]<= float(eachAA[1]) < limits[2]:	#weightage multiplied by 2 (low) and by weightage 
							wei = 2*float(inpwei[key][3])
						elif  limits[2]<= float(eachAA[1]) < limits[3]:	#weightage multiplied by 2 (low) and by weightage 
							wei = 3*float(inpwei[key][3])
						elif  limits[3]<= float(eachAA[1]) <= limits[4]:	#weightage multiplied by 2 (low) and by weightage 
							wei = 4*float(inpwei[key][3])
					elif int(inpwei[key][2]) == -1: 			#depleted AA	
						if 0 <= float(eachAA[1]) < limits[1]:	#freq very high within range - low threshold
							wei = 4*float(inpwei[key][3])	#weightage multiplied by 2 (low) and by weightage 
						elif limits[1] <= float(eachAA[1]) < limits[2]:	#weightage multiplied by 4 (high) and by weightage 
							wei = 3*float(inpwei[key][3])
						elif limits[2] <= float(eachAA[1]) < limits[3]:	#weightage multiplied by 4 (high) and by weightage 
							wei = 2*float(inpwei[key][3])
						elif limits[3] <= float(eachAA[1]) < limits[4]:	#weightage multiplied by 4 (high) and by weightage 
							wei = 1*float(inpwei[key][3])
					
			eachAA.append(wei)
			neweachAA.append(eachAA)
	return(kmerAAfreq)
					
				
		
	
def weightageouter(lys,kmerAAfreq):
	neweachAA=[]
	inpwei = {}
	with open('weightage_for_AAouter13mer.txt', 'rb') as csv_file:
		next(csv_file)
    		for row in csv.reader(csv_file, delimiter='\t'):
        		inpwei[row[0]] = row[1:]
	for eachAA in kmerAAfreq:  #list with AA and freq
		for key in inpwei:
			if eachAA[0] == key:
				maxx=inpwei[key][0]
				wei=float(0.000000)
				limits=numpy.linspace(0,float(maxx),5)
				if int(inpwei[key][2]) == 1:			#enriched AA
					if 0.07 <= float(eachAA[1]) < limits[1]:	#freq very high within range
						wei = 1*float(inpwei[key][3])	#weightage multiplied by 4 (high) and by weightage factor based on p value
					elif  limits[1]<= float(eachAA[1]) < limits[2]:	#weightage multiplied by 2 (low) and by weightage 
						wei = 2*float(inpwei[key][3])
					elif  limits[2]<= float(eachAA[1]) < limits[3]:	#weightage multiplied by 2 (low) and by weightage 
						wei = 3*float(inpwei[key][3])
					elif  limits[3]<= float(eachAA[1]) <= limits[4]:	#weightage multiplied by 2 (low) and by weightage 
						wei = 4*float(inpwei[key][3])
				elif int(inpwei[key][2]) == -1: 			#depleted AA	
					if 0 <= float(eachAA[1]) < limits[1]:	#freq very high within range - low threshold
						wei = 4*float(inpwei[key][3])	#weightage multiplied by 2 (low) and by weightage 
					elif limits[1] <= float(eachAA[1]) < limits[2]:	#weightage multiplied by 4 (high) and by weightage 
						wei = 3*float(inpwei[key][3])
					elif limits[2] <= float(eachAA[1]) < limits[3]:	#weightage multiplied by 4 (high) and by weightage 
						wei = 2*float(inpwei[key][3])
					elif limits[3] <= float(eachAA[1]) < limits[4]:	#weightage multiplied by 4 (high) and by weightage 
						wei = 1*float(inpwei[key][3])
				
		eachAA.append(wei)
	return(kmerAAfreq)
def AAcomp(lys,mer):	#passing each lysine and  its 13mer, AA count is done
	aa=['A','C','E','D','G','F','I','H','K','M','L','N','Q','P','S','R','T','W','V','Y']
	AAcount=[]
	for a in aa:
		val=mer.count(a)
		val=round(float(val)/13.00000,3)
		temp=[a,val]		#each AA and its count is appended to  list AAcount
		AAcount.append(temp)
	return (AAcount)

#############################Secondary structure prediction#######################
###################run_psipred####################################
def run_psipred(proteinid):
	results_paths = ["/home/shyamili/shyamili/fasta_files/"]
	for path in results_paths:
	    os.chdir(path)
	    files = os.listdir(path)
	    for file in files:
		if fnmatch.fnmatch(file, "*horiz"):
		    uniprotid = file.rstrip(".horiz")
		    if uniprotid == proteinid:
		    	horizfile=file
			sspf=open(horizfile,'r')
			global ssplist
			ssplist=sspf.readlines()
			return(ssplist)
			
def extract(ssplist):
	pred=[]
	seq=''
	for line in ssplist:
		if 'Pred:' in line:
			line=line.replace('Pred:','')
			line=line.strip(' ')
			pred.append(line.rstrip('\n'))
	sspff=''.join(pred)
	return (sspff)
def countssp(sspmot):
	gre=''
	val=0.00
	sspweigh=0.00
	types=['H','E','T']
	for i in types:
		Hc=sspmot.count('H')
		Ec=sspmot.count('E')
		Tc=sspmot.count('C')
	if Hc > 0 and Ec == 0 and Tc == 0:
		val="100"
		sspweigh=-0.05
	elif Hc == 0 and Ec > 0 and Tc == 0:
		val="010"
		sspweigh=-0.05
	elif Hc == 0 and Ec == 0 and Tc > 0:
		val="001"
	elif Hc > 0 and Ec > 0 and Tc == 0:
		val="110"
		sspweigh=0.05
		if Hc > Ec:
			gre=12
		if Ec > Hc:
			gre=21
	elif Hc > 0 and Ec == 0 and Tc > 0:
		val="101"
		if Hc > Tc:
			gre=13
		if Tc > Hc:
			gre=31
	elif Hc == 0 and Ec > 0 and Tc > 0:
		val="011"	
		if Ec > Tc:
			gre=23
		if Tc > Ec:
			gre=32
	elif Hc > 0 and Ec > 0 and Tc > 0:
		val="111"
		l={1:Hc,2:Ec,3:Tc}
		sspweigh=0.5
		from collections import OrderedDict
		l_sor = OrderedDict(sorted(l.items(), key=lambda x: x[1],reverse=True))
		for s in l_sor:	
			gre+=str(s)
	else:
		val="000"
	fssplist=[(float(Hc)),(float(Ec)),(float(Tc)),val,gre,sspweigh]
	return fssplist
def logopredict(keyy,valuee):
	q1={'K':3,'D':2,'A':1}
	q2={'G':1,'Y':1,'Q':1}
	q3={'R':-3,'H':-1,'C':-2,'S':-1}
	q4={'C':-2,'H':-1,'K':-3}
	q5={'R':3}
	q6={'R':2}
	q78={'S':-3}
	posq3k={'K':[1,2,3,4,5]}
	negq3k={'K':[7,8,9,10]}
	posq3r={'R':[4,5]}
	posq3s={'S':[0,2,3]}
	posor={'R':[2,3,4]}
	negos={'S':[0,1,2,9,10,11]}
	negop={'P':[1]}
	posq4={'C':[0,2],'H':[1,1],'K':[3,3]}
	negq7={'L':[3,-1]}
	negq3={'L':[1,-0.5]}
	negq4={'L':[4,-1]}
	negq8={'L':[0,-1],'L':[2,-1]}
	enrich={'D':2,'G':3,'L':2,'Q':3,'A':2}
	deple={'K':-3,'C':-3,'H':-3,'P':-3,'M':-2}
	en=de=0
	q1v=q2v=q3v=q4v=q5v=q6v=q7v=q8v=q9v=q10v=q11v=q12v=q13v=q14v=q15v=q16v=0
	for i in range(len(valuee)):
		for k, v in enrich.iteritems():
			en+=valuee[0].count(k)*int(v)
		for k,v in deple.iteritems():
			de+=valuee[0].count(k)*int(v)
		for k,v in q1.iteritems():
			q1v+=valuee[2].count(k)*int(v)
		for k,v in q2.iteritems():
			q2v+=valuee[3].count(k)*v
		for k,v in q3.iteritems():
			q3v+=valuee[2].count(k)*v
		for k,v in q4.iteritems():
			q4v+=valuee[3].count(k)*v	
		for k,v in q5.iteritems():
			q5v+=valuee[7].count(k)*v
		for k,v in q6.iteritems():
			q6v+=valuee[8].count(k)*v
		for k,v in q78.iteritems():
			q7v+=valuee[7].count(k)*v
		for k,v in q78.iteritems():
			q8v+=valuee[8].count(k)*v
		q9v=valuee[4].count('L')*2
		q10v=valuee[5].count('F')*1
		q11v=valuee[6].count('P')*-1
		for k,v in posq3k.iteritems():
			for i in v:
				if valuee[0][i] == k:
					q12v+=1
		for k,v in negq3k.iteritems():
			for i in v:
				if valuee[0][i] == k:
					q12v-=4
		for k,v in posq3r.iteritems():
			for i in v:
				if valuee[0][i] == k:
					q12v-=3
		for k,v in posq3s.iteritems():
			for i in v:
				if valuee[0][i] == k:
					q12v-=3
		for k,v in posor.iteritems():
			for i in v:
				if valuee[1][i] == k:
					q12v+=2
		for k,v in negos.iteritems():
			for i in v:
				if valuee[1][i] == k:
					q12v-=3
		for k,v in negop.iteritems():
			for i in v:
				if valuee[1][i] == k:
					q12v-=1
		for k,v in negq7.iteritems():
			if valuee[2][v[0]] == k:
				q13v+=v[1] 
		for k,v in negq3.iteritems():
			if valuee[2][v[0]] == k:
				q14v+=v[1]
		for k,v in negq4.iteritems():
			if valuee[3][v[0]] == k:
				q15v+=v[1]
		for k,v in negq8.iteritems():
			if valuee[8][v[0]] == k:
				q16v+=v[1]   
		count=q1v+q2v+q3v+q4v+q5v+q6v+q7v+q8v+q9v+q10v+q11v+q12v
		score=en+de+q12v
		return(score)
############################main progrma########################
		
def runkmer(uni,seqq):			#protein seq in string format
	print "MAIN PROGRAM PROTEIN ID",uni
	kmer={}
	kmer2={}
	for lys in range(len(seqq)):
		if seqq[lys] == "K":
			ins=(int(lys)-12)
			off=(int(lys)+13)
			inset=(int(lys)-1)-5
			offset=(int(lys)-1)+8
			upstream=''
			downstream=''
			if inset < 0:
				temp=inset
				inset = 0
				add='-'*abs(int(temp))
				upstream=add
			if offset > len(seqq):
				temp=offset-len(seqq)
				offset=len(seqq)
				downstream='-'*abs(temp)	
			ubikmer=(''.join((upstream+seqq[inset:offset]+downstream)))
			kmer[lys]=ubikmer
			if inset< 0 or ins < 0:
				temp=ins
				ins = 0
				add='-'*abs(int(temp))
				upstream=add
			if offset > len(seqq) or off > len(seqq):
				temp=off-len(seqq)
				off=len(seqq)
				downstream='-'*abs(temp)	
			ubifullmer=(''.join((upstream+seqq[ins:off]+downstream)))
			ubioutermer=ubifullmer[:6]+ubifullmer[19:25]
			kmer[lys]=[ubikmer,ubioutermer]
			kmer2[lys]=[ubikmer,ubioutermer,ubikmer[:6],ubikmer[7:],ubikmer[2:11],ubikmer[7:9],ubikmer[5],ubioutermer[:6],ubioutermer[6:]]

	for key,value in kmer2.iteritems():	
		count=logopredict(key,value)
		value.append(count)
############################################## computing AA composition #####################################
	for key in (kmer):		#key is lysine value is 13mer
		AAlist=AAcomp(key,kmer[key][0])
		kmer[key]=[kmer[key],AAlist]	#appended AAcounts to the respective keys (lysines kmers in dict)

	for key in kmer:
		AAweigh=weightage(key,kmer[key])
		kmer[key]=AAweigh

	for key,value in kmer.iteritems():
		a=0	
		freqlist=kmer[key][1]
		for l in freqlist:
			a+=(l[2])
			sumofallfreq=a
		value.append(sumofallfreq)
################################# outer k-mer composiotin ###############################
	for key,value in kmer.iteritems():		#key is lysine value is 13mer
		AAlist=AAcomp(key,kmer[key][0][1])
		AAouterweigh=weightageouter(key,AAlist)
		value.append(AAouterweigh)

	for key,value in kmer.iteritems():
		a=0	
		freqlist=kmer[key][3]
		for l in freqlist:
			a+=(l[2])
			sumofallfreq=a
		value.append(sumofallfreq)
	for key ,value in kmer2.iteritems():
		for k,v in  kmer.iteritems():
			if k==key:
				logoscore=value[9]
		value.append(logoscore)
	
	for key,val in kmer.iteritems():
		#Pscore=val[2]+val[3][5]+val[4][2][5]+val[4][3]+val[4][6]+  #aascore,disordered,sspscore,sspscore2,outer13merscore
		#print "values",val
		#print val[4]
		#print val[2],val[3][5],val[4][0][2][5],val[4][0][3],
		#print val[4][3],val[6]
		#Pscore=val[2]+val[3][5]+val[4][0][2][5]+val[4][0][3]+val[4][3]+val[6]
		#Pscore=	val[4][3]+val[2]+val[6]+val[3][5]+val[4][0][2][5]+val[4][0][3]	
		for k,v in kmer2.iteritems():
			#print "before",k,key,uni
			if key == k:
				Pscore=v[9]+val[4]
		val.append(Pscore)
	print "writing output file"
	print uni+'.txt'
	sorted_kmer2 = sorted(kmer2.items(), key=lambda x: (float(x[1][10])),reverse=True)  #sorted based on weightage
	
	sorted_kmer = sorted(kmer.items(), key=lambda x: (float(x[1][5])),reverse=True)  #sorted based on weightage
	#print "pscoresorted",uni,sorted_kmer	

	f8=open("/home/shyamili/shyamili/Tool/ubiresults.txt",'a')
	top8=sorted_kmer[:8]
	for key,val in top8:
		for k,v in sorted_kmer2:
			if key == k:
				if v[9] >= 0:
					f8.write(uni+"\t"+str(val[0])+"\t"+str(int(key)+1)+"\t"+str(val[5])+"\n")
	#print "sortes_kmer",sorted_kmer
		
		
#####################################opening input file#####################################



fastadict={}
with open(sys.argv[1]) as f:
    lines=f.read()
    lines=lines.split('>')
    # print "spmnsd",lines
    lines=['>'+x for x in lines[1:]]
    #print lines
    for x in lines:
        file_name=x.split('\n')[0][1:]  #use this variable to create the new file
	file_name=file_name.split("|")[1]
        fil=open("/home/shyamili/shyamili/Tool/fasta_files/"+file_name+'.fasta','w')
        fil.write(x)
        fil.close()
	temp=[]
	temp=x.split("\n")[1:]
	#print temp
	seq=''.join(temp)
	seq=seq.strip("\n")
    	fastadict[file_name]=seq
#print fastadict


      
'''cmd = 'python blast_dir_files.py'
code=subprocess.call(cmd, shell=True)
print "exit status: "+str(code)
import subprocess
cmd = 'python run_psipred40.py'
code=subprocess.call(cmd, shell=True)
print "exit status: "+str(code)
print "dictionary fasta file",fastadict'''
for i in fastadict:
	runkmer(i,fastadict[i])
		
	


