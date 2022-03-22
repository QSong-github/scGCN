import argparse
#import threading
import time
import os
import configparser
import logging
import coloredlogs
import sys
import subprocess
import multiprocessing
import pandas as pd

class NewConfigParser(configparser.ConfigParser):
    def optionxform(self, optionstr):
        return optionstr

def err_exit():
    sys.exit('\033[1;31;47m!!The program exited abnormally, please check the log file !!\033[0m')

def main():
    pwd= os.getcwd()
    parser = argparse.ArgumentParser(description='scGCN pipline')
    parser.add_argument('-i', type=str, required=True, metavar='input_list',help='The count matix and cell type.')
    parser.add_argument('-o', type=str, metavar='outputdir', default=os.sep.join([pwd, 'output']),
                        help='The output dir.default:' + os.sep.join([pwd, 'output']))
    parser.add_argument('--log', type=str, metavar='logfile',help='log file name,or log will print on stdout')
    parser.add_argument('--debug', action='store_true', default=False, help='The log file will contain the output of each software itself, which is convenient for finding errors (-log is required)')

    args = parser.parse_args()

    if not args.log:
        coloredlogs.install(
                fmt='%(asctime)s: %(levelname)s\t%(message)s',
                evel='info'.upper(), stream=sys.stdout
            )
    else:
        if args.debug:
            l='debug'
        else:
            l='info'
        print('Logs will be written to {}'.format(args.log))
        if os.path.exists(args.log):
            os.remove(args.log)
        logging.basicConfig(
                    filename=args.log,
                    filemode='a',
                    format='%(asctime)s: %(levelname)s\t%(message)s',
                    datefmt='%H:%M:%S',
                    level=l.upper()
                    )
    home_dir=os.path.dirname(os.path.abspath(__file__))

    logging.info("Start checking the input file..")
    try:
        input_list=pd.read_csv(args.i, sep='\t', header=0, index_col='Species')
    except FileNotFoundError:
        logging.error('Checking input file :[ERR] -- The input file :%s  does not exist' %(args.input_list))
        err_exit()

    if input_list.columns.tolist()!=["specie","count","type","dbclass"]:
        logging.error('The input file :%s format error.' % (args.input_list))
        err_exit()
    for f in input_list.count:
        if not os.path.exists(f):
            logging.error('The input count file :%s does not exist.'%f)
            err_exit()
    for f in input_list.type:
        if not os.path.exists(f):
            logging.error('The cell type count file :%s does not exist.'%f)
            err_exit()

    os.mkdir(args.o)
    currentPath = os.path.realpath(__file__)
    for ref in input_list[input_list.dbclass=="ref"].itertuples():
        for query in input_list[input_list.dbclass=="query"].itertuples():
            os.mkdir(os.sep.join([args.o,ref.Species,query.Species]))
            logging.info('Start processing %s to %s  type predictions.' % (ref.Species, query.Species))
            logging.info('Data pre-processing...')
            os.chdir(os.sep.join([args.o,ref.Species,query.Species]))
            stdout=subprocess.run(['Rscript',os.sep.join([currentPath,'data_preprocess.R']),ref.count,ref.type,query.count,query.type],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logging.debug('{} stdout:\n'.format('Data pre-processing') +stdout.stdout.decode('utf-8'))
            logging.debug('{} stderr:\n'.format('Data pre-processing') +stdout.stderr.decode('utf-8'))
            logging.info("Data pre-processing finish.")

            logging.info("Cell type predictions...")
            stdout = subprocess.run(['python',os.sep.join([currentPath,'train.py'])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logging.debug('{} stdout:\n'.format('Data pre-processing') +stdout.stdout.decode('utf-8'))
            logging.debug('{} stderr:\n'.format('Data pre-processing') +stdout.stderr.decode('utf-8'))
            logging.info("Cell type predictions finish.")

    logging.info("ALL finish.")

if __name__ == '__main__':
    main()