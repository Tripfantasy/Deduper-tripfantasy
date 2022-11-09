#!/usr/bin/env python3.10
import argparse
import re
from itertools import islice

def get_args():
    parser = argparse.ArgumentParser(description = "Program to Dedupe SAM Files")
    parser.add_argument("-f","--filename", help ="Input sorted SAM File name", type = str, required = True)
    parser.add_argument("-u","--umi",help = "File containing known UMIs",type = str, required = True)
    parser.add_argument("-o","--output", help = "Output de-duped SAM file.", type = str, required = False)
    return parser.parse_args()
args = get_args()


#Function chunk: get_umis() , strandflag(), clipcig()  
def get_umis():
    "Function to read in a file containing known UMIS and add them to a list"
    umis = []
    for line in umi_file:
            line = line.strip()
            umis.append(line)
            umi_len = str(len(umis))
    print("UMI list made with: " + umi_len + " UMIs!")
    return(umis)

def clipcig(CIGAR: str, POSITION:int, strand: str) -> int:
    "Interprets cigar string (=CIGAR) and provides modification to position"
    if strand == "+":                                   # + strand and - strand are clipped differently 
        S = re.findall(r'(\d+)S', CIGAR)                # Regular expression to find digits preceeding a specific character in a string. **Stores as list of STR**
        if S == []:                                     # If the regex returns empty, move on. 
            pass
        else:                                           # If there is soft clipping, position adjusts by subtracting start_clip
            start_clip = int(S[0])
            POSITION -= start_clip
    else:
        S = re.findall(r'(\d+)S', CIGAR)                # - strand adds the end-clip. 
        if S == []:
            pass
        else:
            start_clip = int(S[1])
            POSITION += start_clip
        M = re.findall(r'(\d+)M', CIGAR)                # - strand adds the M sum. 
        if M == []:
            pass
        else:
            M = list(map(int, M))
            POSITION += sum(M)
        D = re.findall(r'(\d+)D', CIGAR)                # - strand adds D sum
        if D == []:
            pass
        else:
            D = list(map(int, D))
        POSITION += sum(D)
        N = re.findall(r'(\d+)N', CIGAR)               # - strand adds N sum 
        if N == []:
            pass
        else:
            N = list(map(int, N))
        POSITION += sum(N)
    return POSITION

def cigarmutate(CIGAR: str) -> list:
    "Turn cigar string into a list parsed by letter"
    cigar_info = "[0-9]+[MIDNSHPX]+"
    ciglist = re.findall(cigar_info,CIGAR)
    return ciglist

def clipcig2(ciglist:list, POSITION:int, strand:str)-> int:
    "Interprets cigar str using cigarmutate as base input. (Alt method)"
    m = 0
    i = 0
    d = 0
    sp = 0
    ls = 0
    rs = 0
    adj_pos = 0
    for item in ciglist:
        letter = item[-1]
        digit = int(item[0:-1])
        if letter == "M":
            m += digit
        elif letter == "D":
            d += digit
        elif letter == "N":
            sp += digit
        elif letter == "S":
            if letter in ciglist[0]:
                ls += digit
            elif letter in ciglist[-1]:
                rs += digit
    if strand =="+": 
        adj_pos = (int(POSITION) - ls)
    elif strand == "-":
        adj_pos = (int(POSITION) + m + d + sp + rs)
    return adj_pos


def field_grep(line: str)->set:
    "Pull relevant de-duplification info from SAM records"
    if line.strip().startswith("@"):
        # args.output.write(line)
        pass
    else:
        fields = line.split()
        QNAME = re.findall("([a-zA-Z]{8})$",fields[0])
        UMI = QNAME[0]
        FLAG = str(fields[1])
        strand = "+"
        if((int(FLAG) & 16) == 16):
            strand = "-"
        RNAME = str(fields[2])
        CIGAR = str(fields[5])
        POSITION = int(fields[3])
        #POSITION = clipcig(CIGAR, POSITION, strand)
        return [UMI, strand, RNAME, POSITION, CIGAR]


#----------------------Time consuming method [ Attempt 1] Resulted in ~1mil+ excess good reads comparing to consensus----------------------#

# with open(args.filename, "r") as sam, open(args.umi,"r") as umi_file, open(args.output,"w") as outfile: 
#     umis = get_umis()
#     record_checker = []
#     umi_error = 0
#     duplicates = 0
#     good_records = 0
#     header_lines = 0
#     for line in sam:
#         if line[0] == "@":
#             outfile.write(line)
#             header_lines += 1
#         else:
#             record = field_grep(line) 
#             ciglist = cigarmutate(record[4])
#             if record[0] in umis:
#                 adj_pos = clipcig2(ciglist, record[3],record[1])
#                 record[3] = adj_pos
#                 record_ID = ''.join(map(str, record))
#                 if record_ID not in record_checker:
#                     record_checker.append(record_ID)
#                     if len(record_checker) > 25000:
#                         record_checker.pop(0)
#                     good_records += 1
#                     outfile.write(line)
#                 else: 
#                     duplicates += 1
#             else: 
#                 umi_error += 1

#----------------------Dictionary implementation and modified clipcig2 function [Attempt 2] approx 2 min and ~10k+ excess good reads----------------------
with open(args.filename, "r") as sam, open(args.umi,"r") as umi_file, open(args.output,"w") as outfile: 
    umis = get_umis()
    good_reads = {}
    umi_error = 0
    duplicates = 0
    good_records = 0
    header_lines = 0
    for line in sam:                      
        if line.startswith("@"):                                            #Write header lines to output file    
            header_lines += 1
            outfile.write(line)
        else:
            record = field_grep(line)                                       #Pull relevant sam fields record = [UMI, strand, RNAME, POSITION, CIGAR]
            if record[0] not in umis:                                       #Check for bad/unknown umis
                umi_error += 1
                continue                                                    #Process further if umi is good
            ciglist = cigarmutate(record[4])                                #generate cigar string list
            adj_pos = clipcig2(ciglist, record[3],record[1])                #Adjust the position based on ciglist2 function and replace with adjusted value
            record[3] = adj_pos
            record_ID = [record[0], record[1], record[2], record[3]]        #Make key of relevant info for dictionary check
            record_ID = ''.join(map(str, record_ID))                        #Simplify key as a string (Removed ~50 more duplicates)

            if record_ID not in good_reads:                                 #Check for duplicates based on dictionary key
                good_records += 1
                good_reads[record_ID] = 1
                outfile.write(line)
            else:
                duplicates += 1
                                                                            #Print the counts! 
print("PCR Duplicates:" + str(duplicates))
print("UMI Error: " + str(umi_error))
print("Good Records: " + str(good_records))
print("Header Lines: " + str(header_lines))
                    

