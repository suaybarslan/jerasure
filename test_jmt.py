#!/usr/bin/env python
# Copyright (C) 2014 Quantum Corporation <contact@quantum.com>
#
# Author: Suayb S. Arslan <suayb.arslan@quantum.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library Public License for more details.
#
# This function assumes that you have succesfully installed the libraries
# GF-Complete and Jerasure version @ https://github.com/suaybarslan/jerasure
# and have the following executables ready:
# 
# 1. encoder (Jerasure original)
# 2. decoder (Jerasure original)
# 3. encoderMT2 (Multi-threaded Jerasure)
# 4. decoderMT2 (Multi-threaded Jerasure)
#
# Additionally, since we use numpy and matplotlib packages,
# you may need to install them before succesfully running this script or
# else you can take out the code snippets below by choosing the appropriate 
# inputs. 
#

import os
import sys
import random
from subprocess import Popen, PIPE
try:
	import numpy as np
	print "INFO: numpy is succesfully imported."
except ImportError:
	print "INFO: numpy is not installed, skipping..."
try:
	import matplotlib.pyplot as plt
	print "INFO: pyplot is succesfully imported."
#except ImportError:
#	print "INFO: pyplot is not installed, skipping..."
except Exception, err:
        sys.stderr.write('ERROR: %s\n' % str(err))

import math
import thread
import filecmp

P_BLOCK_SIZE = 2000  		# in bytes. this is the increment of the packetsize.
B_BLOCK_SIZE = 10000		# in bytes. this is the increment of the buffersize.
R_INC        = 100  		# This is the range increment for the buffer size. 
							# the greater it is, the shorter the simulation time. 
PLOT		 = False        # choose whether you want results to be plotted. 
                  
							# Parameters for testing:
NUMBER_OF_RUN_TIMES = 10 	# number of run times for averaging.
K                   = 16  	# Coding parameter K
M                   = 8 	# Coding parameter M
sumKM               = K + M # Sum of K and M which is the codeword length N

GF_BITS             = 8 	# number of bits per symbols.
numOfErrors         = M 	# number of errors allowed.

def run_encode(input_file, K, M, GF_BITS, PSIZE, BSIZE): 	# func for running encoding.
    CMD = ['encoder', input_file, str(K), str(M), 'reed_sol_van', str(GF_BITS), str(PSIZE), str(BSIZE)]
    process = Popen(CMD, stdout = PIPE, stderr = PIPE)
    stdout, stderr = process.communicate()
    return (stdout, stderr)
    
    
def run_encode_MT(input_file, K, M, GF_BITS, PSIZE, BSIZE): # func for running multi-threaded encoding.
    CMD = ['encoderMT2', input_file, str(K), str(M), str(GF_BITS), str(PSIZE), str(BSIZE)]
    process = Popen(CMD, stdout = PIPE, stderr = PIPE)
    stdout, stderr = process.communicate()
    return (stdout, stderr)
    
def run_decode(input_file):									# func for running decoding.
    CMD = ['decoder', input_file]
    process = Popen(CMD, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    return (stdout, stderr)
    
def run_decode_MT(input_file):								# func for running multi-threaded decoding.
    CMD = ['decoderMT2', input_file]
    process = Popen(CMD, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    return (stdout, stderr)
    
def getErasedFile(path):									# choose file to erase
    files = os.listdir(path)
    index = random.randrange(0, len(files))
    return files[index]

def define_range(start, end, step): 						# define the range of values based on R_INC. 
    while start <= end:
        yield start
        start += step 
    
    
def main(arg):
	#input_file = 'InputBIN.bin'  # input function you may change the input function.
	input_file = arg[1]
	# check if the file exists:
	if not os.path.exists(input_file):
		print "File does not exist. Check the file name."
		sys.exit(1)
	print "file is found.......................................[OK]."
    
	x_range = []
	y_range_enc = []
	y_range_enc_mt = []
	y_range_dec = []
	y_range_dec_mt = []
    
	range_values = define_range(30, 1500, R_INC) # generator object
	length_values = len(list(range_values))
	cnt = 0
	for s in define_range(5, 1500, R_INC):
		max_avg_enc_speed = 0; max_avg_dec_speed = 0; 
		max_avg_enc_speed_mt = 0; max_avg_dec_speed_mt = 0; 
		avg_enc_best = 0; avg_dec_best = 0;
		avg_enc_mt_best = 0; avg_dec_mt_best = 0;
        
		for nn in range(1,11):
			PBS = nn*P_BLOCK_SIZE;
			avg_enc = 0; avg_dec = 0;
			avg_enc_mt = 0; avg_dec_mt = 0;
			for j in range(NUMBER_OF_RUN_TIMES):
				try:
					# ENCODE/DECODE w/ original Jerasure encoder/decoder functions.
					stdout, stderr = run_encode(input_file, K, M, GF_BITS, PBS, B_BLOCK_SIZE*s)
                    
					# print float(stdout.splitlines()[0].split()[2])
					# avg_enc += float(stdout.splitlines()[0][-15:].lstrip())
					avg_enc += float(stdout.splitlines()[0].split()[2])
					# now erase the appropriate set of data files:
					for kk in range(1,numOfErrors+1):
						while True:
							ErasedFile = getErasedFile('Coding/')
							if not ErasedFile.startswith(os.path.splitext(input_file)[0]+'_m'):
								if os.path.exists('Coding/'+ErasedFile):
									break
						os.remove('Coding/'+ErasedFile)
                        
					stdout, stderr = run_decode(input_file)
					
					# print float(stdout.splitlines()[0].split()[2])
					# avg_dec += float(stdout.splitlines()[0][-15:])
					# print float(stdout.splitlines()[0].split()[2])
					avg_dec += float(stdout.splitlines()[0].split()[2])
					# compare the files for correct decoding...
                    
					decodedFile = 'Coding/'+os.path.splitext(input_file)[0]+'_decoded.txt';
					if os.path.exists(decodedFile):
						if not filecmp.cmp(decodedFile,input_file):
							print "Incorrect decoding..............[ORJ2.0]";
							sys.exit(1);
						os.remove(decodedFile) 

					# ENCODE/DECODE w/ MT version of Jerasure encoder/decoder functions.
					stdout, stderr = run_encode_MT(input_file, K, M, GF_BITS, PBS, B_BLOCK_SIZE*s)
                                
					# avg_enc_mt += float(stdout.splitlines()[0][-15:])
					avg_enc_mt += float(stdout.splitlines()[0].split()[2])
					# now erase the appropriate set of data files:
					for kk in range(1,numOfErrors+1):
						while True:
							ErasedFile = getErasedFile('Coding/')
							if not ErasedFile.startswith(os.path.splitext(input_file)[0]+'_m'):
								if os.path.exists('Coding/'+ErasedFile):
									break
						os.remove('Coding/'+ErasedFile)
                       
					stdout, stderr = run_decode_MT(input_file)
					# avg_dec_mt += float(stdout.splitlines()[0][-15:])
					avg_dec_mt += float(stdout.splitlines()[0].split()[2])
                    
					decodedFile = 'Coding/'+os.path.splitext(input_file)[0]+'_decoded.txt';
					if os.path.exists(decodedFile):
						if not filecmp.cmp(decodedFile,input_file):
							print "Incorrect decoding..............[MTJ2.0]";
							sys.exit(1);
						os.remove(decodedFile) 
                     
				except:
					print "Error in encoding/decoding module. Check file size."
					sys.exit(1)
                 
			# find minimums:       
			if avg_enc_best < avg_enc:
				avg_enc_best = avg_enc;
			if avg_enc_mt_best < avg_enc_mt:
				avg_enc_mt_best = avg_enc_mt;
			if avg_dec_best < avg_dec:
				avg_dec_best = avg_dec;
			if avg_dec_mt_best < avg_dec_mt:
				avg_dec_mt_best = avg_dec_mt;
       
		# find averages:    
		max_avg_enc_speed = str(avg_enc_best/NUMBER_OF_RUN_TIMES)     
		max_avg_enc_speed_mt = str(avg_enc_mt_best/NUMBER_OF_RUN_TIMES)     
		max_avg_dec_speed = str(avg_dec_best/NUMBER_OF_RUN_TIMES)     
		max_avg_dec_speed_mt = str(avg_dec_mt_best/NUMBER_OF_RUN_TIMES)   
       
		# accumulate averages:     
		x_range.append(s)
		y_range_enc.append(max_avg_enc_speed)
		y_range_enc_mt.append(max_avg_enc_speed_mt)
		y_range_dec.append(max_avg_dec_speed)
		y_range_dec_mt.append(max_avg_dec_speed_mt)
        
        # print test results:
		print "-- BSIZE:" + str(s*B_BLOCK_SIZE) + " (BLOCKSIZE X " + str(s) + ")"    
		print "Avg. Encode Rate (ORJ2.0): " +str(max_avg_enc_speed)+ " MB/sec"
		print "Avg. Encode Rate (MTJ2.0):" +str(max_avg_enc_speed_mt)+ " MB/sec"
		print "Avg. Decode Rate (ORJ2.0): " +str(max_avg_dec_speed)+ " MB/sec"
		print "Avg. Decode Rate (MTJ2.0):" +str(max_avg_dec_speed_mt)+ " MB/sec"
        
    
	if PLOT:   
		print "The simulation is terminated succesfully, ploting the results..."
		print(x_range)
		print(map(float,y_range_enc))
		print(map(float,y_range_enc_mt))
		# Plot encoding results here:
		plt.figure()
		plt.plot(x_range, y_range_enc, label="Jerasure 2.0-Encoding ("+str(K)+","+str(M)+")", lw=2)
		plt.plot(x_range, y_range_enc_mt, label="MT_J2.0-Encoding ("+str(K)+","+str(M)+")", lw=2)
		plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True)
		plt.ylabel("Speed in MB/sec")
		plt.xlabel("Buffersize: BLOCK_SIZE X (x value)")
		plt.grid(True)
		plt.ioff()
		plt.show()
		#plt.show(block=True)
		
		print(x_range)
		print(map(float,y_range_dec))
		print(map(float,y_range_dec_mt))
		# Plot decoding results here:
		plt.figure()
		plt.plot(x_range, y_range_dec, label="Jerasure 2.0-Decoding ("+str(K)+","+str(M)+")", lw=2)
		plt.plot(x_range, y_range_dec_mt, label="MT_J2.0-Decoding ("+str(K)+","+str(M)+")", lw=2)
		plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, shadow=True)
		plt.ylabel("Speed in MB/sec")
		plt.xlabel("Buffersize: BLOCK_SIZE X (x value)")
		plt.grid(True)
		plt.ioff()
		plt.show()
		#plt.show(block=True)
    
if __name__ == '__main__':

	# check the existance of the input:
	if (len(sys.argv) is not 2):
		print "USAGE: python test_jmt.py <filename>"
		sys.exit(1);
	else:
		main(sys.argv)
    
    
    
    
    
    
    
    
    
