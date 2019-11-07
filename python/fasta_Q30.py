import os,re,sys

infile=file("test.fq.txt")
outfile=file("out.txt","w")

def qual_stat(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(q) - 33
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30

arr=[]
for line in infile:
	line=line.strip()
	arr.append(line)
j=0
agct=[]
a,g,c,t,n,dna=0,0,0,0,0,0
total_count = 0
q20_count = 0
q30_count = 0

for i in range(len(arr)):
	if arr[i].startswith("@") and re.search("^A|^G|^C|^T",arr[i+1]):
		j+=1
		#dna+=arr[i+1].count("A")+arr[i+1].count("G")+arr[i+1].count("C")+arr[i+1].count("T")
		dna+=len(arr[i+1])
		a+=arr[i+1].count("A")
		g+=arr[i+1].count("G")
		c+=arr[i+1].count("C")
		t+=arr[i+1].count("T")
		n+=arr[i+1].count("N")
		#total_count += len(arr[i+3])
		q20, q30 = qual_stat(arr[i+3])
		#q20_count += q20
		q30_count += q30
ap=float(a)*100/float(dna)
gp=float(g)*100/float(dna)
cp=float(c)*100/float(dna)
tp=float(t)*100/float(dna)
np=float(n)*100/float(dna)
q30=float(q30_count)*100/float(dna)
print "Total Reads:",j
print "Total Bases:",dna
print "Q30 Bases:",q30_count
print "q30 percents:",("%.2f" % q30)+"%"
print "A:",a,("%.2f" % ap)+"%"
print "G:",g,("%.2f" % gp)+"%"
print "C:",c,("%.2f" % cp)+"%"
print "T:",t,("%.2f" % tp)+"%"
print "N:",n,("%.2f" % np)+"%"

#print("total bases:", total_count)
#print("q20 bases:", q20_count)

#print("q20 percents:", 100 * float(q20_count)/float(total_count))


