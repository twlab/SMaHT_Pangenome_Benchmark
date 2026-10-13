#!/usr/bin/env python3

################################################
#
#   Add read groups to a BAM file using pysam
#    and Samtools to multi-thread compression
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################

################################################
#   Libraries
################################################
import sys, argparse, os
import subprocess
import copy as cp
import pysam as ps

################################################
#   Functions
################################################
def header_as_str(header):
    str_header = ''
    for key in header:
        if type(header[key]) is list:
            for d in header[key]:
                s = '@' + str(key) + '\t'
                for k, v in d.items():
                    s += str(k) + ':' + str(v) + '\t'
                str_header += s.rstrip() + '\n'
        elif type(header[key]) is dict:
            s = '@' + str(key) + '\t'
            for k, v in header[key].items():
                s += str(k) + ':' + str(v) + '\t'
            str_header += s.rstrip() + '\n'
    return str_header

def get_read_group_illumina(query_name):
    qn_as_list = query_name.split(':')
    l = len(qn_as_list) # number of fields in query_name
    if l == 7: # NEW style read name
               #  <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        instrument_run_flowcell = '_'.join(qn_as_list[:3])
        lane = qn_as_list[3]
        return f'{instrument_run_flowcell}.{lane}'
    elif l == 5: # OLD style read name
                 #  H0164ALXX140820:2:1101:10003:23460
                 #  H0164____________ <flowcell ID>
                 #  _____ALXX140820__ <barcode or index in a multiplexed run>
                 #  _______________:2 <lane>
        return '.'.join(qn_as_list[:2])
    elif l == 1: # modified read name, just use a placeholder
        qn_as_list_modified = query_name.split('.')
        modified_read_group = qn_as_list_modified[0]
        return f'{modified_read_group}'
    else:
        sys.exit('\nFORMAT ERROR: read name format {0} not recognized\n'
                    .format(query_name))

def check_EOF(filename):
    EOF_hex = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    size = os.path.getsize(filename)
    fb = open(filename, "rb")
    fb.seek(size - 28)
    EOF = fb.read(28)
    fb.close()
    if EOF != EOF_hex:
        sys.stderr.write('EOF is missing, adding EOF\n')
        fb = open(filename, "ab")
        fb.write(EOF_hex)
        fb.close()
    else:
        sys.stderr.write('EOF is present, skipping\n')

def main(args):

    # Variables
    directory = args['directory'] if args['directory'] else '.'
    platform = args['platform'] if args['platform'] else 'ILLUMINA'
    library = args['library'] if args['library'] else 'LIBRARY'
    samplename = args['samplename']
    threads = int(args['threads']) if args['threads'] else 1

    # Open input BAM file using Samtools and send to pipe
    pipe_in = subprocess.Popen(['samtools', 'view', '--no-PG', '-h', '-@ {0}'.format(threads), args['inputfile']], stdout=subprocess.PIPE)

    # Data structures
    QNAMEs = set()

    # Read header and all possible QNAMEs
    samfile = ps.AlignmentFile(pipe_in.stdout, 'r')

    header = cp.deepcopy(dict(samfile.header))
    header.setdefault('RG', [])

    for read in samfile:
        QNAME = get_read_group_illumina(read.query_name)
        QNAMEs.add(QNAME)

    # Update header with read groups
    {header['RG'].append({'ID': f'{samplename}.{QNAME}', 'SM': samplename, 'PL': platform, 'PU': QNAME, 'LB': f'{samplename}.{library}'}) for QNAME in QNAMEs}

    # Open output file
    filename = directory + '/' + args['inputfile'].split('/')[-1].split('.')[0] + '_rg' + '.bam'
    bamfile = open(filename, 'w')

    # Buffer to stream output
    pipe_out = subprocess.Popen(['samtools', 'view', '--no-PG', '-b', '-S', '-h', '-@ {0}'.format(threads), '-'], stdin=subprocess.PIPE, stdout=bamfile)

    # Write header
    pipe_out.stdin.write(header_as_str(header).encode())

    # Add read groups to alignments
    pipe_in = subprocess.Popen(['samtools', 'view', '--no-PG', '-h', '-@ {0}'.format(threads), args['inputfile']], stdout=subprocess.PIPE)
    samfile = ps.AlignmentFile(pipe_in.stdout, 'r')
    for read in samfile:
        QNAME = get_read_group_illumina(read.query_name)
        ID = f'{samplename}.{QNAME}'
        # not using pysam to add read tag because somehow it is doing
        # something extremely slow and locking all the multi-threading
        read_str = read.tostring() + '\t' + 'RG:Z:{0}'.format(ID) + '\n'
        pipe_out.stdin.write(read_str.encode())

    pipe_out.communicate()
    bamfile.close()

    # Check if output file has EOF, if not add EOF
    check_EOF(filename)

    # Create index for output file
    if args['index']:
        ps.index(filename)

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add read groups to a BAM file for a single sample')

    parser.add_argument('-i','--inputfile', help='Input BAM file', required=True)
    parser.add_argument('-s','--samplename', help='Name of the sample', required=True)
    parser.add_argument('-l','--library', help='Identifier for the sequencing library preparation [LIBRARY]', required=False)
    parser.add_argument('-p','--platform', help='Name of the sequencing platform (ILLUMINA) [ILLUMINA]', required=False)
    parser.add_argument('-d','--directory', help='Directory to use to write results [.]', required=False)
    parser.add_argument('-t','--threads', help='Number of threads to use for compression/decompression [1]', required=False)
    parser.add_argument('-x','--index', action='store_true', help='Create index for the output file, the input file must be sorted by coordinates')

    args = vars(parser.parse_args())

    main(args)

#end if
