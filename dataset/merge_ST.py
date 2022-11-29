#!/usr/bin/python
import sys
oldfile = open(sys.argv[1],'r')
newfile = open(sys.argv[2],'r')

newst = {}
newfile.readline()
for line in newfile:
    c = line.split("\t")
    newst[c[0]]=c[1]

newdata = {}
for line in oldfile:
    c = line.split(" ")
    d = c[0].split("_")
    idc = d[0]
    if not idc in newst:
#        print("coucou")
        newdata[idc] = c[1]
    else:
#        print("not coucou")
        newdata[idc] = newst[idc]

for i in newdata:
    print(i+"\t"+newdata[i][:-1])


