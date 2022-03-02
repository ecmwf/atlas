#!/usr/bin/env python3

#
# This script checks all example-grids and updates the uid to a newly calculated uid
# The only argument to this script is the path to the "atlas-grids" executable
#
# !!! WARNING !!! 
#   Execution overwrites existing ".yml" files in the same directory as this script!
#

def run_program_get_error(file):
    import os
    from subprocess import check_output, CalledProcessError

    command = [program, file, '--check' ]
    
    from subprocess import Popen, PIPE
    
    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    output, error = p.communicate()
    if p.returncode != 0:
       return str(error.strip(),'ascii')
    return None


def process_file( file, newfile ):

    import string
    print( "File: " , file )

    f = open(file,"r+")
    text = f.read()
    f.close()

    newtext = text

    lines = run_program_get_error( file )
    if lines:
        for line in lines.splitlines():
            import re
            matched = re.match( r"Check failed: grid uid (.*) expected to be (.*)$", line )
            if matched :
                print( line )
                replace = matched.group(1)
                uid = re.findall( r"^  uid : (.*)$", text, re.MULTILINE)[0]
                print("Search/Replace ",uid, " to ", replace )
                newtext = text.replace( uid, replace )
    else:
        print( "Nothing to replace, copy-only" )

    print( "New file: ", newfile)
    f = open(newfile,"w")
    f.write(newtext)
    f.close()


import sys
program = sys.argv[1]

print( "Updating uids of grids with program [",program,"]")

import os
directory = os.path.dirname(os.path.realpath(__file__))

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith('.yml'):
            file = root+"/"+file
            newfile = file
            process_file(file,newfile)
