blosum=open(r'C:\course work\689\A2\BLOSUM62.txt')
blsm=[]

for line in blosum:
    blsm.append(line.split())

import numpy as np
def score(a,b):
    z=float(blsm[blsm[0].index(a)+1][blsm[0].index(b)+1])
    return(z);
    
seq1='GLSDGEWQLVLNVWGKVEADVAGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGNTVLTALGGILKKKGHHEAELTPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMSKALELFRNDMAAKYKELG'
seq3='GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELG'
seq2='GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDRFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELG'
pen=-8
m=len(seq1)
n=len(seq2)
mat=np.zeros((m+1,n+1))
for i in range(0,m+1):
    mat[i][0]=pen*i
for j in range(0,n+1):
    mat[0][j]=pen*j
for i in range(1,m+1):
    for j in range(1,n+1):
        subs=mat[i-1][j-1]+score(seq1[i-1],seq2[j-1])
        dele=mat[i-1][j]+pen
        insrt=mat[i][j-1]+pen
        mat[i][j]=max(subs,insrt,dele)
aln1,aln2='',''
i,j=m,n
while i>0 and j>0:
    current=mat[i][j]
    diaganol=mat[i-1][j-1]
    up=mat[i][j-1]
    left=mat[i-1][j]
    if current==diaganol+score(seq1[i-1],seq2[j-1]):
        aln1+=seq1[i-1]
        aln2+=seq2[j-1]
        i,j=i-1,j-1
		
    elif current==left+pen:
        aln1+=seq1[i-1]
        aln2+='-'
        i-=1
    elif current==up+pen:
        aln1+='-'
        aln2+=seq2[j-1]
        j-=1
while i>0:
    aln1+=seq1[i-1]
    aln2+='-'
    i-=1
while j>0:
    aln1+='-'
    aln2+=seq2[j-1]
    j-=1		
aln1=aln1[::-1]
aln2=aln2[::-1]
i,j=0,0
seqscore=0
identity=0
symbol=''
for i in range(0,len(aln1)):
	if aln1[i]==aln2[i]:
		identity+=1
		seqscore+=score(aln1[i],aln2[i])
		symbol+= '|'
	elif aln1[i]=='-' or aln2[i]=='-':
		seqscore+=pen
		symbol+=' '
	elif aln1[i]!=aln2[i] and aln1[i]!='-' and aln2[i]!='-':
		seqscore+=score(aln1[i],aln2[i])
		symbol+=' '
				
identity= float(identity)/len(aln1) *100
	
print('Sequence similarity:', seqscore)
print ('Identity Score', identity,'%')
print (aln1)
print (symbol)
print (aln2)

