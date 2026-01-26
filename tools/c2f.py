#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# (C) Copyright 1996-2017 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


from argparse import ArgumentParser
import re
import time

INDENT = "    "

supported = [
    'void',
    'char',
    'short',
    'int',
    'long',
    'size_t',
    'float',
    'double',
    '_Bool',
]

c2ftype = {
    'char':'character(c_char)',
    'short':'integer(c_short)',
    'int':'integer(c_int)',
    'long':'integer(c_long)',
    'size_t':'integer(c_size_t)',
    'float':'real(c_float)',
    'double':'real(c_double)',
    '_Bool':'logical(c_bool)',
}

ftype2use = {
    'character(c_char)':'c_char',
    'integer(c_short)':'c_short',
    'integer(c_int)':'c_int',
    'integer(c_long)':'c_long',
    'integer(c_size_t)':'c_size_t',
    'real(c_float)':'c_float',
    'real(c_double)':'c_double',
    'logical(c_bool)':'c_bool'
}

function_signature = re.compile(r'''^
(                                            #1 Type
  (const\s+)?                                #2 Leading const
  ([^\*\&\s\[]+)                             #3
  \s*
  (\&|\*)?                                   #4
  \s*
  (const)?                                   #5
  \s*
  (\&|\*)?                                   #6
  \s*
)
(.+?)\s*                                   #7 Function name
\((.*?)\)$                                 #8 Argument list
''',re.VERBOSE)

class ParsingFailed(Exception):
    def __init__(self,msg):
        Exception.__init__(self, msg)

class Argument:

    arg_signature = re.compile(r'''^
    (                                            #1 Type
      (const\s+)?                                #2 Leading const
      ([^\*\&\s\[]+)                             #3
      \s*
      (\&|\*)?                                   #4
      \s*
      (const)?                                   #5
      \s*
      (\&|\*)?                                   #6
      \s+?
      (\&|\*)?                                   #7
    )
    (.+?)      #8
    (\[\])?$   #9
    ''',re.VERBOSE)

    def __init__(self,arg,debug=False):
        self.c = arg
        self.type=""
        self.value=""
        self.ref=False
        
        match = Argument.arg_signature.match(arg)
        if not match:
            raise ParsingFailed("Could not parse argument "+arg)

        self.type  = str(match.group(1)).strip()
        self.base_type  = str(match.group(3)).strip()
        self.name = str(match.group(8)).strip()
        if( match.group(9) and self.base_type in supported ):
            self.array = True
        else:
            self.array = False
        
        if( "&" in self.type ):
            self.ref = True
        else:
            self.ref = False
        
        if( "*" in self.type ):
            self.ptr = True
        else:
            self.ptr = False
        
        if debug:
            print ("name = ",self.name)
            print ("type = ",self.type)
            print ("base_type = ",self.base_type)
            print ("array = ",self.array)
            print ("ref = ",self.ref)
            print ("ptr = ",self.ptr)
 
        self.use = set()
        self.fortran = ""
        if( self.base_type == "char" and self.ptr and not self.ref  ) :
            self.fortran += "character(c_char), dimension(*)"
            self.use.add("c_char")
        else:
            if(self.ptr):
                self.fortran += "type(c_ptr)"
                self.use.add("c_ptr")
            else:
                try:
                    ftype = c2ftype[self.base_type]
                    self.fortran += ftype
                    self.use.add(ftype2use[ftype])
                except KeyError:
                    raise ParsingFailed("Could not parse argument "+arg)

            if( self.array ):
                self.fortran += ", dimension(*)"

            elif( not self.ref ):
                self.fortran += ", value"
                
        self.fortran += " :: "+self.name
        
        if( debug ):
            print( self.fortran )

class ReturnType:

    def __init__(self,line):

        match = function_signature.match(line)
        if( not match ):
            raise ParsingFailed("Could not parse return type for line "+line)
        self.type = str(match.group(1)).strip()
        self.base_type = str(match.group(3)).strip()
        self.ptr = match.group(4) or match.group(6)
        self.use = set()
        self.ftype = ""

        if( self.type != "void" ):
            
            if( self.ptr ):
                self.ftype = "type(c_ptr)"
                self.use.add("c_ptr")
            else:
                if( not self.type in supported ):
                    raise ParsingFailed("Cannot return type "+self.type+" for statement "+str(line))
                try:
                    self.ftype = c2ftype[self.base_type]
                    self.use.add(ftype2use[self.ftype])
                except KeyError:
                    raise ParsingFailed("Could not parse return type for statement "+line)

    def __bool__(self): # python 3
        return self.type != "void"
        
class Function:
    
    def __init__(self,line):
        self.line = line
        self.return_type = ReturnType(line)
        self.name = self.extract_function_name(line)
        self.arguments = [ Argument(function_arg) for function_arg in self.extract_function_args(line) ]

    def extract_function_name(self,line):
        match = function_signature.match(line)
        if match:
            fn = str(match.group(7))
            return fn

    def extract_function_args(self,line):
        match = function_signature.match(line)
        if match:
            args = str(match.group(8)).strip()
            if args:
                return [arg.strip() for arg in args.split(',')]
        return []

    def fix_fortran_linewidth(self,code,linewidth):
        newcode = ""
        for line in code.split('\n'):
            comment = line.strip()[0] == "!"
            if comment :
                width = linewidth
                end = ''
                begin = '!   '
            else :
                width = linewidth-1
                begin = "  "+"&"
                end = "&"
            while len(line) > width:
                newcode += line[:width] + end + "\n"
                line = begin + line[width:]
            newcode += line + "\n"
        return newcode

    def fortran_interface(self,linewidth=80):
        sep = '!'+'-'*(linewidth-1)
        intf = sep+"\n"
        intf += "! "+self.line+"\n"
        intf += sep+"\n"
        if( self.return_type ):
            intf += "function "
        else:
            intf += "subroutine "
        intf += self.name+"("

        if self.arguments:
            intf += " "+", ".join([arg.name for arg in self.arguments])+" "
        intf += ")"+" bind(C,name=\""+self.name+"\")\n"
        
        self.use = set()
        self.use.update(self.return_type.use)
        for arg in self.arguments:
            self.use.update(arg.use)
        if( self.use ):
            intf += INDENT+"use iso_c_binding, only: "+", ".join(self.use)+"\n"

        
        if( self.return_type ):
            intf += INDENT+self.return_type.ftype + " :: " + self.name + "\n"
        for arg in self.arguments:
            intf += INDENT+arg.fortran + "\n"
        
        if( self.return_type ):
            intf += "end function\n"
        else:
            intf += "end subroutine\n"

        intf += sep
        
        intf = self.fix_fortran_linewidth(intf,linewidth)
        return intf

class Code:
    
    def __init__(self):
        self.content=''
    
    def readline(self,line):
        self.content+=line
    
    def statements(self):
        self.statements = [ statement.strip() for statement in self.content.split(';')[:-1] ]
        return self.statements

# -----------------------------------------------------------------------------

parser = ArgumentParser()
parser.add_argument("file", type=str, help="file to parse")
parser.add_argument("-o", "--output", type=str, help="output")
parser.add_argument("-m", "--module", type=str, help="module")
parser.add_argument("-t", "--types", type=str, help='e.g. {"idx_t":"integer(c_int)","gidx_t":"integer(c_long)"}' )
parsed = parser.parse_args()
input  = parsed.file
output = parsed.output
module = parsed.module

if( parsed.types ) :
    import json
    extra_types = json.loads( str(parsed.types) )
    for k in extra_types:
        supported.append(k)
        c2ftype[k] = extra_types[k]

if not module:
    module = output
    module = str(re.sub(r'(.+)(\..+)',r'\1',module))
    module.replace('-','_')

code = Code()

with open(input,'r') as file:
    content = file.read()
    externC = re.compile(r'^extern\s+"C"(\s*{)?')
    
    regex_externC = [re.compile(r'^extern\s+"C"(\s*{)?'),re.compile('^{')]
    in_externC = [False,False]

    for line in content.splitlines():
        line = line.strip()
        
        # Go in extern "C" region
        if not in_externC[0] :
            match = regex_externC[0].match(line)
            if match:
                in_externC[0]=True
                if match.group(1):
                    in_externC[1] = True
                    continue
        if in_externC[0] and not in_externC[1] :
            match = regex_externC[1].match(line)
            if match:
                in_externC[1] = True
                continue

        # Go out of extern "C" region
        if re.search("^}",line):
            in_externC[0] = False
            in_externC[1] = False

        # Skip empty lines and comment lines
        if re.search("^//|^$",line):
            continue
        
        # Line to parse
        if in_externC[1]:
            code.readline(line)

intf = ""
intf += "! Do not modify!\n"
intf += "! Automatically generated file that contains Fortran bindings\n"
intf += "! of C funtions in file "+input+" contained within extern \"C\" scope.\n"
intf += "! Time of creation: "+str(time.asctime(time.localtime()))+"\n\n"
intf += "module "+module+"\n"
intf += "implicit none\n"
intf += "interface\n\n"
for statement in code.statements():
    try:
        intf += Function(statement).fortran_interface()+"\n"
    except ParsingFailed as e:
        print("\n\n"+"-"*80+"\n"+"Automatic generation of Fortran bindings failed for file\n\n"+
              "    "+input+"\n\n"+"ParsingFailed for statement:\n\n"+
              "    "+statement+"\n\nError: "+str(e)+"\n"+"-"*80+"\n")
        raise
intf += "end interface\n"
intf += "end module\n"

with open(output,'w') as out:
    out.write(intf)
