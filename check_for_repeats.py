# Author: Roshan M Regy
# EmailID: roshanm.regy@tamu.edu
import re
import pandas as pd
import sys
from colorama import Fore
from colorama import Style

# Function to list substrings of 2 or more 
# characters 
def FindSubStr(value):
    substr=[]
    for i in range(2,int(len(value)/2)):
        for start in range(len(value)-i):
            sub = value[start:start+i]
            if sub not in substr:
                substr.append(sub)
    return substr

# Function to match substrings to 
# original sequence
def PatternMatch(substr,value):
    output = []
    for sub in substr:
        check = re.findall(sub,value)
        if len(check)>1:
            splitvalue = value.split(sub)
            # For pretty printing on console
            # Only keep writestr lines if 
            # printing is not needed
            printstr = ""
            writestr= ""
            if value[:len(sub)]==sub:
                printstr+=Fore.RED
                printstr+=sub
                writestr+='|'+sub+'|'
                printstr+=Style.RESET_ALL
            for  s,part in enumerate(splitvalue):
                if s<len(splitvalue)-1:
                    printstr+=part
                    writestr+=part
                    printstr+=Fore.RED
                    printstr+=sub
                    writestr+='|'+sub+'|'
                    printstr+=Style.RESET_ALL
                else:
                    printstr+=part
                    writestr+=part
            if value[-len(sub):]==sub:
                printstr+=Fore.RED
                printstr+=sub
                writestr+='|'+sub+'|'
                printstr+=Style.RESET_ALL
            print (Fore.RED+sub+Style.RESET_ALL,f'%s'%printstr)
            output.append([sub,len(check),writestr])
    return output

seq = open(sys.argv[1]).read().strip()
#seq="MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS"
substr = FindSubStr(seq)
output = PatternMatch(substr,seq)
df=pd.DataFrame(output,columns=['repeat','# of repeats','prettyprint'])
df = df.sort_values('# of repeats',ascending=False)
df.to_csv('%s_repeatcheck.csv'%(sys.argv[1][:-4]),index=False)
