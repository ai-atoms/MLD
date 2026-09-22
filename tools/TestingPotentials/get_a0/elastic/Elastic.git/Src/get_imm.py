import re
import sys

def get_imm_from_din(dinNDM):
 f=open(dinNDM,'r')
 for line in f:
   if 'imm' in line:
    word=line.replace("="," ").split()
 f.close()
 return word[1];


imm=get_imm_from_din(sys.argv[1])
print(imm)



