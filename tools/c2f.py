# (C) Copyright 1996-2014 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

#! /usr/bin/env python
from argparse import ArgumentParser
import re
import os
parser = ArgumentParser()
parser.add_argument("files", type=str, nargs='*', help="files to parse")
parser.add_argument("-o", "--output", type=str, help="output")
parsed = parser.parse_args()
files = parsed.files
output = parsed.output
used_types = []

# A tabulation:
TAB = "  "

# -------------------------------------------------------------------------
# These dictionaries give the Fortran type and its KIND for each C type:
# -------------------------------------------------------------------------
# One word types:
TYPES_DICT = { 
    "int":("integer(c_int)","c_int"),
    "long":("integer(c_long)","c_long"),
    "short":("integer(c_short)","c_short"),
    "_Bool":("logical(c_bool)","c_bool"),
    "boolean":("logical(c_bool)","c_bool"),
    "double": ("real(c_double)","c_double"),
    "float":("real(c_float)","c_float"),
    "size_t":  ("integer(c_size_t)","c_size_t"),
     }

# Two words types:
TYPES2_DICT = {
    "long double": ("real(c_long_double)","c_long_double"),
    "unsigned long":("integer(c_long)","c_long"),
    "unsigned short":("integer(c_short)","c_short"),
    "unsigned int":("integer(c_int)","c_int"),
    }

class ParsingError(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

#---------------------------------------------------------------------------
# Regular expressions used to identify the different parts of a C c_line:
#---------------------------------------------------------------------------
REGEX_RETURNED_TYPE = re.compile( "^ *([_0-9a-zA-Z ]+ *\**&?)" )
REGEX_FUNCTION_NAME = re.compile( "([0-9a-zA-Z_]+) *\(" )
REGEX_ARGUMENTS = re.compile( "\(([&0-9a-zA-Z_\s\,\*\[\]]*)\).*;$" )
REGEX_ARGS = re.compile( " *([&0-9a-zA-Z_\s\*\[\]]+),?" )
REGEX_VAR_TYPE = re.compile( " *([_0-9a-zA-Z]+)[ |\*]" )
REGEX_TYPE = re.compile( "^ *((const )?\w+)[ \*]?" )
REGEX_VAR_NAME = re.compile( "[ |\*&]([_0-9a-zA-Z]+)(?:\[\])?$" )

def multiline(line, maxlength): 
    """Split a long line in a multiline, following Fortran syntax."""
    result = ""
    while len(line) > maxlength-1:
        result += line[0:maxlength-1] + "&\n"
        line = "&"+ line[maxlength-1:]
    result += line
    return result

def iso_c_binding(declaration, returned):
    """ Returns the Fortran type corresponding to a C type in the 
        ISO_C_BINDING module and the KIND type """
    global REGEX_TYPE
    global TYPES_DICT

    try:
        c_type = REGEX_TYPE.search(declaration).group(1)
    except AttributeError:
        return "?", "?"    # error

    # Is it an array ?
    if declaration.find("[") != -1:
        array = ", dimension(*)"
    else:
        array = ""

    # Is it a pointer ?
    if declaration.find("*") != -1:
        # Is it a string (char or gchar array) ?
        # TODO: what about "unsigned char"   "guchar" gunichar ?
        if ((c_type.find("char") != -1) or (c_type.find("char*") != -1)) and (not returned):
            if declaration.find("**") != -1:
                return "type(c_ptr), dimension(*)", "c_ptr"
            else:
                return "character(kind=c_char), dimension(*)", "c_char"
        else:
            return "type(c_ptr)", "c_ptr"

    # Other cases:
    if len(declaration.split()) >= 3:   # Two words type
        for each in TYPES2_DICT:
            if set(each.split()).issubset(set(declaration.split())):
                return TYPES2_DICT[each][0] + array, TYPES2_DICT[each][1]
    else:  # It is a one word type
        for each in TYPES_DICT:
            if each in c_type.split():
                return TYPES_DICT[each][0] + array, TYPES_DICT[each][1]

    # It is finally an unknown type:
    return "?", "?"

def C_F_binding(c_line):

    type_returned = REGEX_RETURNED_TYPE.search(c_line)
    try:
        function_type = type_returned.group(1)
    except AttributeError:
        raise ParsingError, "Returned type not found "+ c_line

    # Will it be a Fortran function or a subroutine ?
    if (function_type.find("void") != -1) and (function_type.find("*") == -1):
        f_procedure = "subroutine "
        f_the_end   = "end subroutine"
        isfunction  = False
        f_use       = ""
    else:
        f_procedure = "function "
        f_the_end   = "end function"
        isfunction  = True
        returned_type, iso_c = iso_c_binding(type_returned.group(1), True)
        f_use = iso_c
        if returned_type.find("?") != -1:
            error_flag = True
            raise ParsingError, "Unknown data type:    " + type_returned.group(1) + "  " + c_line 

    # f_name is the name of the function in fortran:
    function_name = REGEX_FUNCTION_NAME.search(c_line)
    try:
        f_name = function_name.group(1)
    except AttributeError:
        raise ParsingError, "Function name not found "+c_line
                           
    arguments = REGEX_ARGUMENTS.search(c_line)
    try:
        args = REGEX_ARGS.findall(arguments.group(1))
    except AttributeError:
        raise ParsingError, "Arguments not found " + c_line

    # Each argument of the function is analyzed:
    declarations = ""
    args_list = ""
    for arg in args:
        if arg != "void":
            try:
                var_type = REGEX_VAR_TYPE.search(arg).group(1)
            except AttributeError:
                raise ParsingError, "Variable type not found " + c_line

            f_type, iso_c = iso_c_binding(arg, False)
            if iso_c not in used_types:
                used_types.append(iso_c)
                        
            if f_type.find("c_") != -1:
                if f_use == "":
                    f_use = iso_c
                else:
                    # each iso_c must appear only once:
                    REGEX_ISO_C = re.compile( "("+iso_c+")"+"([^\w]|$)" )
                    if REGEX_ISO_C.search(f_use) == None:
                        f_use += ", " + iso_c 
            elif f_type.find("?") != -1:
                raise ParsingError, "Unknown data type:    " + arg + "  " + c_line

            
            try:
                var_name = REGEX_VAR_NAME.search(arg).group(1)
            except AttributeError:
                raise ParsingError,"Variable name not found "+c_line+"  "+arg
    
            # Array with unknown dimension are passed by adress,
            # the others by value:
            if (f_type.find("(*)") != -1):
                passvar = ""
            else:
                if (arg.find("&") !=  -1):
                    passvar = ""
                else:
                    passvar = ", value"
                
            declarations += 1*TAB + f_type + passvar + " :: " + var_name + "\n"
            if args_list == "":
                args_list = var_name
            else:
                args_list += ", " + var_name
    interface = 0*TAB + "! " + c_line + "\n"
    first_line = 0*TAB + f_procedure + f_name + "(" + args_list + ') bind(c,name="{f_name}")'.format(f_name=f_name)
    
    interface += multiline(first_line, 79) + " \n"

    interface += 1*TAB + "use iso_c_binding, only: " + f_use + "\n"
    if isfunction:
        interface += 1*TAB + returned_type + " :: " + f_name + "\n"
    interface += declarations
    interface += 0*TAB + f_the_end + "\n\n" 

    return interface

def parse_file(c_file):
    file_name, file_ext = os.path.splitext(c_file)
    with open(c_file,'r') as open_file:
        text = open_file.read()
        repl = r""
        text = re.sub('(?ms).*?extern "C"\s*{', '', text)
        text = re.sub('(?ms)}.*$', '', text)
        text = re.sub('(?ms)^\s+', '', text)
        lines = [ line.strip() for line in text.strip().splitlines() ]
    f_file_path = output#file_name+'-auto.f90'
    module_name = file_name.split('/')[-1]+"_c_binding"
    with open(f_file_path,'w') as f_file:
        import time
        f_file_header =  "! Automatically generated by c2f.py on {time}\n".format(time=time.asctime(time.localtime()))
        f_file_header += "! Based on file {filename}\n".format(filename=c_file)
        f_file_header += "! Please do not modify.\n"

        f_file.write(f_file_header+"\nmodule {module_name}\nimplicit none\ninterface\n\n".format(module_name=module_name))
        
        for line in lines:
          try:
            f_file.write( C_F_binding(line) + '\n\n')
          except ParsingError, e:
            print("\n----------------------------------------\n"+
                  "Parsing failed for file\n    "+c_file+" : \n"+str(e)+
                  "\n----------------------------------------\n"+text)
        f_file.write("end interface\nend module {module_name}\n".format(module_name=module_name))



for file in files:
    parse_file(file)