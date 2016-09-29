#!/usr/bin/env python

import sys

if __name__=='__main__':
    if(len(sys.argv)!=4):
        raise Exception('Invalid number of arguments!');
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    extension = sys.argv[3]
    
    f = open(inFile)
    x = f.readlines()
    f.close()
    
    aList = []
    for i in x:
        aList.append(i.rstrip('\n')+extension)
    
    # write out the new file
    f = open(outFile, 'w')
    for i in aList:
        f.write('%s\n' % i)
    f.close()