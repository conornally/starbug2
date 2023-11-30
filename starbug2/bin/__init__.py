import os
from starbug2.utils import printf,perror

EXIT_SUCCESS=0
EXIT_FAIL   =1
EXIT_EARLY  =2
EXIT_MIXED  =3

def usage(docstring,verbose=0):
    if verbose: perror(docstring)
    else: perror("%s\n"%docstring.split('\n')[1])
    return 1

def parsecmd(args):
    cmd=os.path.basename(args[0])
    return cmd,args[1:]

