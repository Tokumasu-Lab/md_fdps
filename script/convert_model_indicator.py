#!/usr/bin/python
#--------------------------------------------------------------------------------------------------
# This script use python
#   "md_fdps_enum_model.hpp" file generator.
#   source: ***.mol2 files in ./model
#
# developed: Kiyoto Kawai, Tokumasu Lab in Tohoku Univ.
#    21 July 2017
#--------------------------------------------------------------------------------------------------
# coding: UTF-8

import subprocess
import os

target_locate = "./model/"
target_ext  = ".mol2"
result_file = "./src/md_fdps_enum_model.hpp"

def load_atom_name(model, atom_list):
    model_file = model + target_ext
    f = open(model_file)
    line_file = f.readlines()
    mode = "header"
    for line in line_file:

        #--- skip comment line
        if len(line) <= 2:
            continue
        if line[0] == "!":
            continue
        if line[:2] == "//":
            continue

        #--- convert string to list
        line_list = line.split()

        #--- decode header
        if line_list[0] == "@<TRIPOS>MOLECULE":
            if mode != "header":
                print "ERROR: each ***.mol2 file must have 1 molecule."
                sys.exit()
            mode = "header"
            continue
        if line_list[0] == "@<TRIPOS>ATOM":
            mode = "atom"
            continue
        if line_list[0] == "@<TRIPOS>BOND":
            mode = "bond"
            continue

        #--- decode atom data
        if mode == "atom":
            atom = line_list[1]
            if not (atom in atom_list):
                atom_list.append(atom)

    f.close


#--- enum class definition
def make_enum_class(typename, value_list, result):
    result.append("enum class " + typename + " : uint_fast32_t {")
    for value in value_list:
        result.append(" "*4 + value + ",")
    result.append("};")
    result.append("")


#--- "whatis(enum value)" function (return std::string)
def make_whatis(typename, value_list, result):
    result.append(" "*4 + "inline std::string whatis(const " + typename + " &e){")

    result.append(" "*8 + "switch (e) {")

    for value in value_list:
        result.append(" "*12 + "case " + typename + "::" + value + ":")
        result.append(" "*16 + 'return "' + value + '";')
        result.append(" "*12 + "break;")
        result.append("")

    result.append(" "*12 + "default:")
    result.append(" "*16 + 'throw std::out_of_range("undefined enum value in ' + typename + '");')

    result.append(" "*8 + "}")

    result.append(" "*4 + "}")
    result.append("")


#--- "which_***(std::string)" function (return enum class value)
def make_which(typename, value_list, result):
    result.append(" "*4 + "inline " + typename + " which_" + typename + "(const std::string &str){")

    result.append(" "*8 + 'if(str == "' + value_list[0] + '"){')
    result.append(" "*12 + 'return ' + typename + '::' + value_list[0] + ';')
    result.append(" "*8 + "}")

    for value in value_list[1:]:
        result.append(" "*8 + 'else if(str == "' + value + '"){')
        result.append(" "*12 + 'return ' + typename + '::' + value + ';')
        result.append(" "*8 + "}")

    result.append(" "*8 + "else {")
    result.append(" "*12 + 'std::cerr << "  ' + typename + ': input = " << str << std::endl;')
    result.append(" "*12 + 'throw std::out_of_range("undefined enum value in ' + typename + '");')
    result.append(" "*8 + "}")

    result.append(" "*4 + "}")
    result.append("")


#--- "size_***()" function (return the number of defined value)
def make_size(typename, value_list, result):
    result.append(" "*4 + "inline size_t size_" + typename + "(){")

    result.append(" "*8 + "return " + str(len(value_list)) + ";")

    result.append(" "*4 + "}")
    result.append("")


#--- "whatis(std::tuple<...>)" function (return std::string)
def make_whatis_tuple(type_list, result):
    #--- string converter: int and double
    result.append(" "*4 + "template<typename T>")
    result.append(" "*4 + "std::string whatis_string_Impl(T const &v){")
    result.append(" "*8 + "return std::to_string(v);")
    result.append(" "*4 + "}")
    result.append("")

    #--- string converter: ENUM
    for t in type_list:
        result.append(" "*4 + "inline std::string whatis_string_Impl(" + t + " const &e){")
        result.append(" "*8 + "return whatis(e);")
        result.append(" "*4 + "}")
        result.append("")

    #--- recursive part
    result.append(" "*4 + "template<typename Tuple, size_t Index = std::tuple_size<Tuple>::value-1>")
    result.append(" "*4 + "struct whatis_Impl{")
    result.append(" "*8 + "static void apply(std::string &str, Tuple const &tuple){")
    result.append(" "*12 + "whatis_Impl<Tuple, Index-1>::apply(str, tuple);")
    result.append(" "*12 + 'str += ", " + whatis_string_Impl(std::get<Index>(tuple));')
    result.append(" "*8 + "}")
    result.append(" "*4 + "};")
    result.append("")

    #--- terminator part
    result.append(" "*4 + "template<typename Tuple>")
    result.append(" "*4 + "struct whatis_Impl<Tuple, 0>{")
    result.append(" "*8 + "static void apply(std::string &str, Tuple const &tuple){")
    result.append(" "*12 + "str = whatis_string_Impl(std::get<0>(tuple));")
    result.append(" "*8 + "}")
    result.append(" "*4 + "};")
    result.append("")

    #--- interface
    result.append(" "*4 + "template<typename Tuple>")
    result.append(" "*4 + "inline std::string whatis(Tuple const &tt){")
    result.append(" "*8 + "std::string str{""};")
    result.append(" "*8 + "whatis_Impl<Tuple>::apply(str, tt);")
    result.append(" "*8 + 'return "(" + str + ")";')
    result.append(" "*4 + "}")
    result.append("")

#--- function of "std::cout << enum class::(value)"
def make_output_ostream(type_list, result):
    for t in type_list:
        result.append(        "inline std::ostream& operator << (std::ostream& s, const " + t + " &e){")
        result.append(" "*4 + "s << ENUM::whatis(e);")
        result.append(" "*4 + "return s;")
        result.append(        "}")
        result.append("")


if __name__ == "__main__":

    #--- get result of ls
    pwd = os.getcwd()
    os.chdir(target_locate)
    cmd = "ls"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    os.chdir(pwd)
    ls_result = p.communicate()[0]
    ls_files = ls_result.split()
    #print ls_files

    #--- get target file name
    model_list = []
    for f in ls_files:
        if f[(len(f)-len(target_ext)):] == target_ext:
            model_list.append(f[:(len(f)-len(target_ext))])

    #--- import all atom names
    atom_list = []
    print "  model list : " + str(model_list)
    for model in model_list:
        load_atom_name(target_locate + model, atom_list)

    print "  atom list  : " + str(atom_list)


    #--- output C++ header file (definition of enum class)
    result = []

    #------ file header
    result.append("//" + "-"*88)
    result.append("//  This file is enum class of atom types and molecular types.")
    result.append("//    generated by ./convert_model_indicator.py")
    result.append("//" + "-"*88)
    result.append("#pragma once")
    result.append("")
    result.append("#include <cstdint>")
    result.append("#include <string>")
    result.append("#include <tuple>")
    result.append("#include <stdexcept>")
    result.append("")
    result.append("")

    #------ enum class of model name (using file name)
    make_enum_class("MolName", model_list, result)
    #------ enum class of atom name
    make_enum_class("AtomName", atom_list, result)

    #--- adding interface for enum class value
    result.append("namespace ENUM {")
    result.append("")
    #------ for MolType
    result.append(" "*4 + '//--- basic interface for "MolType"')
    make_whatis("MolName", model_list, result)
    make_which( "MolName", model_list, result)
    make_size(  "MolName", model_list, result)
    #------ for AtomType
    result.append(" "*4 + '//--- basic interface for "AtomType"')
    make_whatis("AtomName", atom_list, result)
    make_which( "AtomName", atom_list, result)
    make_size(  "AtomName", atom_list, result)

    type_list = []
    type_list.append("MolName")
    type_list.append("AtomName")
    result.append(" "*4 + '//--- convert std::tuple<...> to std::string')
    make_whatis_tuple(type_list, result)
    result.append("}")
    result.append("")

    #--- add to global namespace
    result.append('//--- output function as "std::cout << (enum class::value)"')
    make_output_ostream(type_list, result)
    result.append("")


    #--- output file
    print "  output file: " + result_file
    f = open(result_file, "w")
    for line in result:
        f.write(line + "\n")
    f.close
