#!/usr/bin/python
#--------------------------------------------------------------------------------------------------
#  This script uses python.
#  "pdb" and "crd" file converter for VMD from the output file of "md_fdps_main.cpp".
#
#    2016/3/3  Kiyoto Kawai, Tokumasu Lab., Tohoku Univ.
#--------------------------------------------------------------------------------------------------
# coding: UTF-8

import sys
import glob
import os.path

from optparse import OptionParser

#--- output pdb file setting
out_pdb_name = "vmd_movie.pdb"

#--- output crd file setting
out_crd_name = "vmd_movie.crd"

#--- input file setting
#------ file name of AMBER VMD files
locate_pdb  = "./pdb/"
in_pdb_head = "vmd"
in_pdb_ext  = ".pdb"

#------ file name of pos files
locate_pos  = "./posdata/"
in_pos_head = "pos"
in_pos_ext  = ".dat"

#------ format of pos files
index_id    = 0
index_pos_x = 4
index_pos_y = 5
index_pos_z = 6


#------ functions ---------------------------------------------------------------------------------
def extract_number_part_(file, header, ext):
    #--- extract number part
    tmp = file.rstrip(ext)
    index  = tmp.rindex(header)
    index += len(header)
    tmp = tmp[index:]

    return tmp


def search_numbered_file_(files, header, ext, tgt_num):
    #--- search target number file
    for file_tmp in files:
        #--- extract number part
        num_tmp = extract_number_part_(file_tmp, header, ext)

        #--- check argument tgt_num
        if len(num_tmp) < len(str(tgt_num)):
            print " ERROR: invalid target number digit."
            print "   file digit:", len(num_tmp), " | target digit:", len(str(tgt_num))
            print "   target digit must be lower than file digit."
            sys.exit()

        #--- adjust digits
        tgt_tmp = "0"*( len(num_tmp) - len(str(tgt_num)) ) + str(tgt_num)

        #--- judge the file_tmp is target or not
        if num_tmp == tgt_tmp:
            return file_tmp

    #--- target was not found
    print "ERROR: target pdb file was not found."
    sys.exit()


def write_sorted_pdb_(target_pdb):
    #--- input target file
    pdb_dic = {}
    for line in open(target_pdb, "r"):
        list = line.split()
        id = int(list[1])
        pdb_dic[id] = line

    #--- get n_atom
    n_atom = len(pdb_dic)

    #--- sort the contents
    tmp = ""
    for i in range(len(pdb_dic)):
        tmp += pdb_dic[i]

    #--- output
    f = open(out_pdb_name, "w")
    f.write(tmp)
    f.close

    return n_atom


def estimate_start_(files, header, ext):
    #--- make number part list
    num_list = []
    for file in files:
        tmp = extract_number_part_(file, header, ext)
        num_list.append(int(tmp))

    return min(num_list)


def estimate_end_(files, header, ext):
    #--- make number part list
    num_list = []
    for file in files:
        tmp = extract_number_part_(file, header, ext)
        num_list.append(int(tmp))

    return max(num_list)


def estimate_cycle_(files, header, ext, i_from):
    #--- make number part list
    num_list = []
    for file in files:
        tmp = extract_number_part_(file, header, ext)
        num_list.append(int(tmp))

    num1 = sys.maxint
    for num in num_list:
        if num1 > num:
            if num > i_from:
                num1 = num

    cycle = num1 - i_from

    return cycle


def load_pos_data_(target_pos):
    #--- input target file
    pos_dic = {}
    for line in open(target_pos, "r"):
        list = line.split()

        #--- skip header
        if len(list) < index_pos_z:
            continue
        if not list[0].isdigit():
            continue

        #--- clear tmp
        pos_tmp = []

        #--- get value
        id   = int(list[index_id])
        pos_tmp.append(float(list[index_pos_x]))
        pos_tmp.append(float(list[index_pos_y]))
        pos_tmp.append(float(list[index_pos_z]))

        #--- add to dictionary
        pos_dic[id] = pos_tmp

    return pos_dic



#------ main sequence -----------------------------------------------------------------------------
if __name__ == "__main__":
    #--- get execute options
    parser = OptionParser()
    parser = OptionParser(usage="%prog [-f (start number)] [-t (end number)] [-c (cycle)]")
    parser.add_option("-f","--from",action="store",type="int",dest="i_from",
                    help="the start file of 'pos(number).dat' with 'vmd(number).dat'")
    parser.add_option("-t","--till",action="store",type="int",dest="i_till",
                    help="the end file of 'pos(number).dat'")
    parser.add_option("-c","--cycle",action="store",type="int",dest="cycle",
                    help="increment of file number, 'pos(number).dat'")

    parser.set_defaults(i_from =  0,
                        i_till = -1,
                        cycle  = -1)

    (options, args) = parser.parse_args()

    i_from = options.i_from
    i_till = options.i_till
    cycle  = options.cycle


    #--- find target pdb file
    pdb_files = glob.glob(locate_pdb+"*"+in_pdb_ext)
    if len(pdb_files) <= 0:
        print "ERROR: there are no '*"+in_pdb_ext+"' files in '"+locate_pdb+"' directory."
        sys.exit()

    #--- estimate i_from
    if i_from <= 0:
        i_from = estimate_start_(pdb_files, in_pdb_head, in_pdb_ext)
        print "start file was estimated = "+str(i_from)

    #--- set target pdb file
    target_pdb = search_numbered_file_(pdb_files, in_pdb_head, in_pdb_ext, i_from)

    #--- make sorted pdb file
    n_atom = write_sorted_pdb_(target_pdb)


    #--- get list of pos file
    pos_files = glob.glob(locate_pos+"*"+in_pos_ext)

    #--- error check
    if len(pos_files) <= 0:
        print "ERROR: there are no '*"+in_pos_ext+"' files in '"+locate_pos+"' directory."
        sys.exit()

    if len(pos_files) < 2:
        print "ERROR: number of '*"+ext+"' files in '"+locate_pos+"' is insufficient."
        print "   more than 2 files are required."
        print pos_files


    #--- estimate i_till
    if i_till <= 0:
        i_till = estimate_end_(pos_files, in_pos_head, in_pos_ext)
        print "end file was estimated = "+str(i_till)

    #--- estimate cycle
    if cycle <= 0:
        cycle = estimate_cycle_(pos_files, in_pos_head, in_pos_ext, i_from)
        print "cycle was estimated = "+str(cycle)

    #--- load first file
    target_pos = search_numbered_file_(pos_files, in_pos_head, in_pos_ext, i_from)
    pos_ref = load_pos_data_(target_pos)


    #--- convert pos file to AMBER trajectory
    #------ header (atom number)
    f = open(out_crd_name, "w")
    f.write(str(n_atom))
    f.write("\n")

    x_max = 0.0
    y_max = 0.0
    z_max = 0.0
    i_step   = i_from
    tot_step = int( (i_till - i_from)/cycle )
    for i in range(tot_step):
        #--- increment
        i_step += cycle

        #--- load target file
        target_pos = search_numbered_file_(pos_files, in_pos_head, in_pos_ext, i_step)
        print "  loading the file '"+target_pos+"'..."
        pos_now = load_pos_data_(target_pos)

        #--- convert 1 dimention value sequence
        pos_list = []
        for j in range(len(pos_now)):
            tmp_pos = pos_now[j]
            for tmp in tmp_pos:
                pos_list.append(tmp)

            x_max = max(tmp_pos[0], x_max)
            y_max = max(tmp_pos[1], y_max)
            z_max = max(tmp_pos[2], z_max)

        #--- convert AMBER trajectory data
        value_count = 0
        for pos in pos_list:
            f.write(str(" %9.4f" % pos))
            value_count += 1
            if value_count >= 6:
                f.write("\n")
                value_count = 0
            else:
                f.write(" ")

    #--- footer
    f.write( "  " + str(x_max) + "  " + str(y_max) + "  " + str(z_max) + "  90.0  90.0  90.0\n" )

    f.close
    print " ...'"+out_crd_name+"' file was generated."
