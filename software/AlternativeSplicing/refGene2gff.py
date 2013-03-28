#!/usr/bin/python
#####################################################################################
#Author: Qi Liu
#Email: liuqi_1984@yahoo.com.cn
#Description: covert the RefGene file of UCSC into GFF format file, and sort the gene
#             by the chromosomes and location on each chromosome.
#####################################################################################
import sys,os
import re
import string

def getGff(file):
   """convert the format the file to gff format"""
   try:
      FI = open(file,"r") 
   except IOError, e:
      print >> sys.stderr, "File read error! ",e.message
      sys.exit(2)
   gffResult = {}
   for line in FI:
      line = line.strip()
      tem = line.split("\t")
      mRnaStart = str(string.atoi(tem[4])+1)
      mRnaEnd = str(string.atoi(tem[5]))
      firstCdsStart = str(string.atoi(tem[6])+1)
      lastUTRStart = str(string.atoi(tem[7])+1)
      exonStarts = tem[9].split(",")
      exonStarts.pop()
      exonEnds = tem[10].split(",")
      exonEnds.pop()
      tss = tem[15].split(",")
      resultOneLine = ""
      resultOneLine += tem[2]+"\t"+"refGene"+"\t"+"mRNA"+"\t"+mRnaStart+"\t"+mRnaEnd+"\t"+"."+"\t"+tem[3]+"\t"+"."+"\t"+"ID="+tem[1]+"; name="+tem[12]+";"+"\n"
      if string.atoi(tem[8]) == 1:
          resultOneLine += tem[2]+"\t"+"refGene"+"\t"+"exon"+"\t"+str(string.atoi(exonStarts[0])+1)+"\t"+exonEnds[0]+"\t"+"."+"\t"+tem[3]+"\t"+"."+"\t"+"Parent="+tem[1]+";"
      else:
        for i in range(0,string.atoi(tem[8])):
         exonStart = str(string.atoi(exonStarts[i])+1)
         exonEnd = str(string.atoi(exonEnds[i]))
         if i == 0:
                 resultOneLine += tem[2]+"\t"+"refGene"+"\t"+"exon"+"\t"+exonStart+"\t"+exonEnd+"\t"+"."+"\t"+tem[3]+"\t"+"."+"\t"+"Parent="+tem[1]+";"+"\n"
         elif i == string.atoi(tem[8])-1:
                 resultOneLine += tem[2]+"\t"+"refGene"+"\t"+"exon"+"\t"+exonStart+"\t"+exonEnd+"\t"+"."+"\t"+tem[3]+"\t"+"."+"\t"+"Parent="+tem[1]+";"
         else:
                 resultOneLine += tem[2]+"\t"+"refGene"+"\t"+"exon"+"\t"+exonStart+"\t"+exonEnd+"\t"+"."+"\t"+tem[3]+"\t"+"."+"\t"+"Parent="+tem[1]+";"+"\n"
      if not gffResult.has_key(tem[2]):
           gffResult[tem[2]] = {string.atoi(mRnaStart):[resultOneLine]}
      else:
           if not gffResult[tem[2]].has_key(string.atoi(mRnaStart)):
               gffResult[tem[2]][string.atoi(mRnaStart)] = [resultOneLine]
           elif resultOneLine not in gffResult[tem[2]][string.atoi(mRnaStart)]:
               gffResult[tem[2]][string.atoi(mRnaStart)].append(resultOneLine)
   FI.close()
   printGff(gffResult)

def getEmbedderNumber(str):
      pat = re.compile(r'chr(\d+)')
      pat2 = re.compile(r'chr(\w)')
      a = pat.search(str)
      if a:
         number = a.group(1)
         return string.atoi(number)
      else:
         b = pat2.search(str)
         if b:
            if b.group(1) == 'X':
              return 23
            elif b.group(1) == 'Y':
              return 24
            elif b.group(1) == 'M':
              return 25
 
def sortByNumber(list):
      """sort the gene by the chromosomes and the location on each chromosome"""
      an = [(getEmbedderNumber(e),e) for e in list]
      an.sort()
      return [ s for f, s in an ]

def printGff(gffResult): 
     for chr in sortByNumber(gffResult.keys()):
        poss = gffResult[chr].keys() 
        poss.sort()
        for pos in poss:
            for ele in gffResult[chr][pos]:
                 print ele
         


if __name__ == "__main__":
  if len(sys.argv) < 2:
     print "Usage: python %s <UCSC_RefGene_File>" %sys.argv[0]
     sys.exit(2)
  getGff(sys.argv[1])
