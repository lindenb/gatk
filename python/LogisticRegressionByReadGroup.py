#!/usr/bin/env python

import os,sys

def exit(msg,errorlevel):
    "Exit the program with the specified message and error code."
    print msg
    sys.exit(errorlevel)

def open_source_file(source_filename):
    "Open the source file with the given name.  Make sure it's readable."
    if not os.access(source_filename,os.R_OK):
        exit("Unable to read covariate counts file '" + sys.argv[1] + "'",1)
    if not source_filename.endswith('.csv'):
        exit("Source file is in incorrect format.  Must be csv.")
    return open(source_filename,'r')

def read_header(source_file):
    "Read the header from the given source file.  Do basic validation."
    header = source_file.readline().split(',')
    if header[0] != 'rg' or header[1] != 'dn':
        exit("Input file is in invalid format.  First two columns should be <read group>,<dinucleotide>",1)
    return header

def create_intermediate_file(source_file,header,read_group,dinuc):
    "Create an intermediate file for a particular read group / dinuc from a given source file"
    base = source_file.name[:source_file.name.rfind('.csv')]
    intermediate_filename = "%s.%s.%s.csv" % (base,read_group,dinuc)
    intermediate_file = open(intermediate_filename,"w")
    intermediate_file.write(','.join(header[2:]))
    return intermediate_file

def open_target_file(target_filename):
    "Open a target file and write out the header."
    target_file = open(target_filename,'w')
    target_file.write("rg\tdinuc\t")
    for p in range(5):
        for q in range(5):
            target_file.write("%d,%d\t" % (p,q))
    target_file.write("\n")
    return target_file

def process_file(source_file,read_group,dinuc,target_file,R_exe,logistic_regression_script):    
    "Process the contents of an intermediate file.  An intermediate file is the specific data arranged per read group, per dinuc."
    base = source_file.name[:source_file.name.rfind('.csv')] + '.' + read_group
    regression_command = ' '.join((R_exe,logistic_regression_script,base,base,dinuc))
    result = os.system(regression_command)
    if result != 0:
        exit('Unable to run linear regression; command was %s, error code was %d' % (regression_command, result),1)
    parameters_filename = '.'.join((base,dinuc,'parameters'))
    if not os.access(parameters_filename,os.R_OK):
        exit("Unable to read output of R from file " + parameters_filename)
    parameters_file = open(parameters_filename,'r')
    parameters = ' '.join([line.rstrip('\n') for line in parameters_file]).split(' ')
    target_file.write('\t'.join([read_group,dinuc]+parameters)+'\n')
    # optionally remove the input and output files.  For now, keep them around for debugging purposes.
    # os.remove('.'.join((base,dinuc,'csv')))
    # os.remove('.'.join((base,dinuc,'parameters')))
    
def compute_logistic_regression(source_filename,target_filename,R_exe,logistic_regression_script):
    "Do the heavy lifting.  Group covariant data into individual output files, compute the logistic regression, and format the results."
    source_file = open_source_file(source_filename)
    target_file = open_target_file(target_filename)

    header = read_header(source_file)

    intermediate_file = None
    read_group = None
    dinuc = None

    for data_line in source_file:
        data_line = data_line.strip()
        if len(data_line) == 0:
            continue
        data = data_line.split(',')
        if read_group != data[0] or dinuc != data[1]:
            if intermediate_file:
                intermediate_file.close()
                process_file(source_file,read_group,dinuc,target_file,R_exe,logistic_regression_script)
            read_group,dinuc = data[0:2]
            intermediate_file = create_intermediate_file(source_file,header,read_group,dinuc)
        intermediate_file.write(','.join(data[2:])+'\n')

    if intermediate_file:
        intermediate_file.close()
        process_file(source_file,read_group,dinuc,target_file,R_exe,logistic_regression_script)

    source_file.close()
    target_file.close()

if __name__ == "__main__":
    # set up sensible defaults for the locations of R and the logistic regression script
    R_exe="/broad/tools/apps/R-2.6.0/bin/Rscript"
    logistic_regression_script="resources/logistic_regression.R"

    if len(sys.argv) < 3:
        exit("Usage: logistic_regression <covariate count input csv> <parameters csv>",1)

    compute_logistic_regression(sys.argv[1],sys.argv[2],R_exe,logistic_regression_script)
