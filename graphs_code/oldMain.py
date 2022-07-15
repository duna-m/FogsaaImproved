# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import subprocess
import os
import sys
import time

from needlemanwunsch import *
from metrics import *



def run(fasta1, fasta2,match, mismatch, insertion, deletion, similarity=0.3):
    print("match:" + str(match)+ "mismatch: "+ str(mismatch)+ "ins: "+str(insertion)+ "dele: "+str(deletion)+ "similarity: "+ str(similarity))
    args = ["./fogsaa.exe", fasta1, fasta2, "1", "0", str(match), str(mismatch), str(insertion), str(deletion), str(similarity)]
    #args = ["./fogsaa.exe", fasta1, fasta2, "1", "0", "0", "-1", "-2", "-3", str(similarity)]
    #args = ["./fogsaa.exe", fasta1, fasta2, "1", "0", "0", "-1", "-2"]
    ret = subprocess.call(args, stdin=None, stdout=None, stderr=None, shell=False)


def peek_line(file):
    pos = file.tell()
    line = file.readline()
    file.seek(pos)
    return line




def file_to_fasta(strands_and_refs,metrics_file):
    line1 = strands_and_refs.readline()
    line2 = strands_and_refs.readline()


    strands_and_refs.readline()
    if line1 and line2:
        fasta1 = open("fasta1.txt","w")
        fasta2 = open("fasta2.txt","w")
        fasta1.write(">seq1\n")
        fasta1.write(line1)
        fasta2.write(">seq2\n")
        fasta2.write(line2)
        if peek_line(strands_and_refs):
            return True
        else:
            return 2
    return False

def fogsaa_loop_runner(args):
    metrics_file = open("aligned_strands_fogsaa.txt", "w")
    strands_and_refs = open(args[0], "r")
    res = True
    num = 0
    avg_time = 0
    while res:
        res = file_to_fasta(strands_and_refs, metrics_file)

        starting_time = time.time()
        run("fasta1.txt", "fasta2.txt", args[1], args[2], args[3], args[4], args[5])
        num += 1
        avg_time += time.time() - starting_time

        res_file1 = open("new_alseq1.txt", "r")
        res_file2 = open("new_alseq2.txt", "r")
        line1_aligned = res_file1.readline()
        line2_aligned = res_file2.readline()
        if res == 2:
            metrics_file.writelines([line1_aligned, "\n", line2_aligned])
            res = False
        else:
            metrics_file.writelines([line1_aligned, "\n", line2_aligned, "\n", "\n"])
    metrics_file.close()

    avg_time /= num
    return avg_time



def old_main(args):

    args_wunsch = [args[0]]

    wunsch_avg_time = wunsch_loop_runner(args_wunsch)

    total = args[1]+args[2]+args[3]
    mis = args[1]/total*3-1
    ins = args[2]/total*3-1
    dele = args[3]/total*3-1

    print(mis,ins,dele)

    args_fogsaa = ["orig_strands.txt", 0, mis, ins, dele, args[4]]
    fogsaa_avg_time = fogsaa_loop_runner(args_fogsaa)

    args_metrics_wunsch = ["orig_strands.txt", "aligned_strands_wunsch.txt", "36","needleman-metrics"]
    print("needleman wunsch metrics are:")
    wunsch_scs, wunsch_lcs, wunsch_lev = metrics_runner(args_metrics_wunsch)

    args_metrics_fogsaa = ["orig_strands.txt", "aligned_strands_fogsaa.txt", "36", "fogsaa_metrics"]
    print("fogsaa metrics are:")
    fogsaa_scs, fogsaa_lcs, fogsaa_lev = metrics_runner(args_metrics_fogsaa)

    return wunsch_avg_time, wunsch_scs, wunsch_lcs, wunsch_lev, fogsaa_avg_time, fogsaa_scs, fogsaa_lcs, fogsaa_lev






    #args = ["./fogsaa.exe", "fasta1.txt", "fasta2.txt", "1", "0", "1", "-1", "-2"]
    #ret = subprocess.call(args, stdin=None, stdout=None, stderr=None, shell=False)



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
