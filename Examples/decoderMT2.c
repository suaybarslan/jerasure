/* *
 * Copyright (c) 2014, James S. Plank and Kevin Greenan
 * All rights reserved.
 *
 * Jerasure - A C/C++ Library for a Variety of Reed-Solomon and RAID-6 Erasure
 * Coding Techniques
 *
 * Revision 2.0: Galois Field backend now links to GF-Complete
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  - Neither the name of the University of Tennessee nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* Jerasure's authors:

   Revision 2.x - 2014: James S. Plank and Kevin M. Greenan.
   Revision 1.2 - 2008: James S. Plank, Scott Simmerman and Catherine D. Schuman.
   Revision 1.0 - 2007: James S. Plank.
 */

/* 
This program takes as input an inputfile, k, m, a coding
technique, w, and packetsize.  It is the companion program
of encoder.c, which creates k+m files.  This program assumes 
that up to m erasures have occurred in the k+m files.  It
reads in the k+m files or marks the file as erased. It then
recreates the original file and creates a new file with the
suffix "decoded" with the decoded contents of the file.

This program does not error check command line arguments because 
it is assumed that encoder.c has been called previously with the
same arguments, and encoder.c does error check.
*/

/*********************************************************************************
	Beta Version (Revision 1.0): Suayb S. Arslan and Hoa Le 
	Jerasure 2.x0 Multi-Thread Decoder Implementation using shared-memory method.
*/
/*
Preface: The original library has different set of erasure coding techniques. 
Based on our testing and simulations, we figured the coding technqiue "reed_sol_van"
performed the best among others and it turns out that its construction aligns more 
friendly with the MT-implementation. Major changes with respect to original encoder.c
are detailed below with each code line item.

There are multiple ways of implementing multi-threading for Jerasure. We took 
a shared-memory (no-mutex) approach. In this approach, threads share the same buffer 
space. Since threads are only allowed to read from but not write to the buffer, we 
have not used mutexes and hence secured good performance. 

Multi-threaded decoder implementation (version 1.0) deals with data and parity blocks 
one after another because parity computation depends on the data blocks being  
available. Exluding this sequential operation, erased data block as well as parity 
block computations are fully multi-threaded (operated upon in parallel). Below is a 
case for instance how data block recovery is made fully multi-threaded. 

In original decoder.c, the last lost data block is recovered using simple XOR for 
performance reasons. However this creates dependency upon the presence of all data
blocks. So we do not invoke XOR for full-threaded implementation.  

SSA @ QTM, 10/18/16.
	                     
Reference Doc: See PPT file Jerasure2,0_setup_math_perf at Quantum/Box repositories. 
**********************************************************************************/

#include <pthread.h> //important to include.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <signal.h>
#include <unistd.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "galois.h"
#include "cauchy.h"
#include "liberation.h"
#include "timing.h"

/* #define MULTIPROCESS */

/* N is not used */
#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

enum Coding_Technique {Reed_Sol_Van, No_Coding};
/* We decided to implement reed_sol_van only due to performance reasons */
/* Global variables for signal handler */
/* No method variable is needed */
int readins, n;
/* No signal handler is used */

/* This is for multi-threading implementation (original decoder does not have it) */
/* Better to put into a header file. We decided to make minimum changes though */
typedef struct{
	int k;
	int w;
	int *matrix_row;
	int *src_ids;
	int dest_id;
	char** data;
	char** coding;
	int blocksize;
}jinput_args;

/* The following function is a modified version of jerasure_matrix_dotprod_mt 
to utilize multi-thread capabilities of modern CPUs. Note that it is a bit 
different from encoderMT2's version.*/
void *jerasure_matrix_dotprod_mt(void *x_void_ptr){
	jinput_args *x_ptr = (jinput_args *)x_void_ptr;    
	int k = x_ptr->k;
	int w = x_ptr->w;
	int *matrix_row = x_ptr->matrix_row;
	int *src_ids = x_ptr->src_ids;
	int dest_id = x_ptr->dest_id;  
	char** data_ptrs = x_ptr->data;
	char** coding_ptrs = x_ptr->coding;
	int size = x_ptr->blocksize;
	
	int init = 0;
	char *dptr, *sptr;
	int i;
    
	dptr = (dest_id < k) ? data_ptrs[dest_id] : coding_ptrs[dest_id-k];
	/* First copy or xor any data that does not need to be multiplied by a factor */
	for (i = 0; i < k; i++) {
		if (matrix_row[i] == 1) {
			if (src_ids == NULL) {
				sptr = data_ptrs[i];
			} else if (src_ids[i] < k) {
			sptr = data_ptrs[src_ids[i]];
			} else {
			sptr = coding_ptrs[src_ids[i]-k];
			}
			if (init == 0) {
				memcpy(dptr, sptr, size);
				init = 1;
			} else {
				galois_region_xor(sptr, dptr, size);
			}
		}
	}

	/* Now do the data that needs to be multiplied by a factor */
	for (i = 0; i < k; i++) {
		if (matrix_row[i] != 0 && matrix_row[i] != 1) {
			if (src_ids == NULL) {
				sptr = data_ptrs[i];
			} else if (src_ids[i] < k) {
				sptr = data_ptrs[src_ids[i]];
			} else {
				sptr = coding_ptrs[src_ids[i]-k];
			}
			switch (w) {
				case 8:  galois_w08_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
				case 16: galois_w16_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
				case 32: galois_w32_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
			}
			init = 1;
		}
	}
	return NULL;
}

int main (int argc, char **argv) {
	FILE *fp;				// File pointer

	/* Jerasure arguments */
	char **data;
	char **coding;
	int *erasures;
	int *erased;
	int *matrix; 
	/* No bitmatrix is defined for reed_sol_van */
	
	/* Parameters */
	int k, m, w, packetsize, buffersize;
	int tech;
	char *c_tech;
	
	int i, j;				// loop control variable, s
	int blocksize = 0;			// size of individual files
	unsigned int origsize;	// size of file before padding - we use unsigned for file sizes > 1GB.
	unsigned int total;		// used to write data, not padding to file - we use unsigned for file sizes > 1GB.
	struct stat status;		// used to find size of individual files
	int numerased;			// number of erased files
		
	/* Used to recreate file names */
	char *temp;
	char *cs1, *cs2, *extension;
	char *fname;
	int md;
	char *curdir;
	
	/* Used to time decoding */
	struct timing t1, t2, t3, t4;
	struct timeval tstart; /* extra timing parameters for measuring performance */
	double tsec;
	double totalsec;
	
	/* Internal variables for Jerasure decoding. Original file do not have them. 
	We need them to be able to implement multi-threaded Jerasure.
	edd: number of erasures in data section
	edp: number of erasures in parity section */
	int edd, edp, ind, *decoding_matrix, *dm_ids, *tmpids;
	int lastdrive;
	
	matrix = NULL;
	dm_ids = NULL; /* Initialization of internal arrays */
	decoding_matrix = NULL; /* Initialization of internal arrays */
	
	totalsec = 0.0;
	#if defined(MULTIPROCESS)
		double sstart; /* extra timing parameters for measuring performance */
	#endif
	
	/* Start timing */
	timing_set(&t1);
	
	/* Error checking parameters */
	if (argc != 2) {
		fprintf(stderr, "usage: inputfile\n");
		exit(0);
	}
	/* Get current working directory for construction of file names */
	curdir = (char*)malloc(sizeof(char)*1000);	
	assert(curdir == getcwd(curdir, 1000));
	
	/* Begin recreation of file names */
	cs1 = (char*)malloc(sizeof(char)*strlen(argv[1]));
	cs2 = strrchr(argv[1], '/');
	if (cs2 != NULL) {
		cs2++;
		strcpy(cs1, cs2);
	}
	else {
		strcpy(cs1, argv[1]);
	}
	cs2 = strchr(cs1, '.');
	if (cs2 != NULL) {
                extension = strdup(cs2);
		*cs2 = '\0';
	} else {
           extension = strdup("");
        }	
	fname = (char *)malloc(sizeof(char*)*(100+strlen(argv[1])+20));

	/* Read in parameters from metadata file */
	sprintf(fname, "%s/Coding/%s_meta.txt", curdir, cs1);
	
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		fprintf(stderr, "Error: no metadata file %s\n", fname);
		exit(1);
	}
	temp = (char *)malloc(sizeof(char)*(strlen(argv[1])+20));
	if (fscanf(fp, "%s", temp) != 1) {
		fprintf(stderr, "Metadata file - bad format\n");
		exit(0);
	}
	
	if (fscanf(fp, "%u", &origsize) != 1) { /* Need u% for unsigned int */
		fprintf(stderr, "Original size is not valid\n");
		exit(0);
	}
	if (fscanf(fp, "%d %d %d %d %d", &k, &m, &w, &packetsize, &buffersize) != 5) {
		fprintf(stderr, "Parameters are not correct\n");
		exit(0);
	}
	c_tech = (char *)malloc(sizeof(char)*(strlen(argv[1])+20));
	if (fscanf(fp, "%s", c_tech) != 1) {
		fprintf(stderr, "Metadata file - bad format\n");
		exit(0);
	}
	if (fscanf(fp, "%d", &tech) != 1) {
		fprintf(stderr, "Metadata file - bad format\n");
		exit(0);
	}
	/* we did not have to use method since we only implement reed_sol_van */
	if (fscanf(fp, "%d", &readins) != 1) {
		fprintf(stderr, "Metadata file - bad format\n");
		exit(0);
	}
	fclose(fp);	
	
	/* Allocate memory */
	erased = (int *)malloc(sizeof(int)*(k+m));
	for (i = 0; i < k+m; i++)
		erased[i] = 0;
	erasures = (int *)malloc(sizeof(int)*(k+m));

	data = (char **)malloc(sizeof(char *)*k);
	coding = (char **)malloc(sizeof(char *)*m);
	if (buffersize != origsize) {
		for (i = 0; i < k; i++) {
			data[i] = (char *)malloc(sizeof(char)*(buffersize/k));
		}
		for (i = 0; i < m; i++) {
			coding[i] = (char *)malloc(sizeof(char)*(buffersize/k));
		}
		blocksize = buffersize/k;
	}

	sprintf(temp, "%d", k);
	md = strlen(temp);
	timing_set(&t3);

	/* Create coding matrix or bitmatrix */ /* We use reed_sol_van only */
	matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
	timing_set(&t4);
	totalsec += timing_delta(&t3, &t4);
	
	
	/* In the original decoder, the erased file indexes are calculated within the while
	loop. This is examined to be inefficient. So we take it out by placing internal 
	variables inside the encoderMT2.c source code (edd, edp, etc.). */
	numerased = 0;
	lastdrive = k;
	edd = 0; edp = 0;
	/* Determine the missing files and erasures. */
	for (i = 1; i <= k+m; i++) {
		if (i <= k){
			sprintf(fname, "%s/Coding/%s_k%0*d%s", curdir, cs1, md, i, extension);
		}else{
			sprintf(fname, "%s/Coding/%s_m%0*d%s", curdir, cs1, md, i - k, extension);
		}    
		fp = fopen(fname, "rb");
		if (fp == NULL) {
			erased[i-1] = 1;
			erasures[numerased] = i-1;
			numerased++;
			if(i <= k){
				edd++;
				lastdrive = i-1;
			}else{
				edp++;
			}
		}else{
			fclose(fp);
		}
	}
	/* Check whether the decoder is able to decode. */	
	/* In the original encoder, this check is done after encoding attempt. In this version,
	we perform the checking before any attempt is made. */
	if (erased[k]) lastdrive = k;
	if(numerased > m){
		printf("Unsuccesful!\n");
		exit(0);
	}	
	/* Check on w */
	if (w != 8 && w != 16 && w != 32){
		fprintf(stderr, "Decoding cannot be terminated successfully!\n");
		exit(0);
	}
    
	/* Generate the decoding matrix if need be. */
	timing_set(&t3); /* Calculation of inverse matrix is part of pure timing estimation */
	if (edd > 1 || (edd > 0 && erased[k])) {
		dm_ids = talloc(int, k);
		decoding_matrix = talloc(int, k*k);
		if (jerasure_make_decoding_matrix(k, m, w, matrix, erased, decoding_matrix, dm_ids) < 0) {
			free(dm_ids);
			free(decoding_matrix);
			fprintf(stderr, "Decoding cannot be terminated successfully!\n");
			exit(0);
		}	
	}
	timing_set(&t4);
	totalsec += timing_delta(&t3, &t4);

	/* Helping Array structs for MT-implementation */
	jinput_args args4d[edd];
	jinput_args args4p[edp];
	pthread_t jmed_mt_thread[edd];
	pthread_t jmep_mt_thread[edp];

	/* Begin decoding process */
	total = 0;
	n = 1;	
	while (n <= readins) {
		/* Open files, check for erasures, read in data/coding */	
		/* Original encoder uses two different for loops. It turns out that we can combine those
		for loops into a single one to reduce the number of instuctions */
		for (i = 1; i <= k+m; i++) {
			if (i <= k){
				sprintf(fname, "%s/Coding/%s_k%0*d%s", curdir, cs1, md, i, extension);
			}else{
				sprintf(fname, "%s/Coding/%s_m%0*d%s", curdir, cs1, md, i - k, extension);
			}
			fp = fopen(fname, "rb");

			if (fp != NULL){
				if (buffersize == origsize) {
					stat(fname, &status);
					blocksize = status.st_size;
					if (i <= k){
						data[i-1] = (char *)malloc(sizeof(char)*blocksize);
						assert(blocksize == fread(data[i-1], sizeof(char), blocksize, fp));
					}else{
						coding[i-k-1] = (char *)malloc(sizeof(char)*blocksize);
						assert(blocksize == fread(coding[i-k-1], sizeof(char), blocksize, fp));
					}			
				}
				else {
					fseek(fp, blocksize*(n-1), SEEK_SET); 
					if (i <= k){
						assert(buffersize/k == fread(data[i-1], sizeof(char), buffersize/k, fp));
					}else{
						assert(blocksize == fread(coding[i-k-1], sizeof(char), blocksize, fp));
					}
				}
				fclose(fp);
			}
		}
    	/* Finish allocating data/coding if needed */
		if (n == 1) {
			for (i = 0; i < numerased; i++) {
				if (erasures[i] < k) {
					data[erasures[i]] = (char *)malloc(sizeof(char)*blocksize);
				}
				else {
					coding[erasures[i]-k] = (char *)malloc(sizeof(char)*blocksize);
				}
			}
		}
	    
		/* DEcode and REpair starts here. */ 
		if (n==1){ /* This is for my own timing */
			gettimeofday(&tstart, NULL);
			#if defined(MULTIPROCESS)
				sstart = (double) tstart.tv_sec + ((double) tstart.tv_usec) / 1000000;
			#endif		
		}
		timing_set(&t3);
		if (tech == Reed_Sol_Van){		
			if(edp == 0 && edd == 0)
				goto edp0edd0; /* If no error is detected, there is no need to run the rest of the code */
			ind = 0;
			if (edd > 1 || (edd > 0 && erased[k])) {
				/* recover raw data: - DECODE process: */
				for (i = 0; i < k; i++) {
					if (erased[i]) {
						args4d[ind].k = k; args4d[ind].w = w; args4d[ind].matrix_row = decoding_matrix+(i*k);
						args4d[ind].src_ids = dm_ids; args4d[ind].dest_id = i; args4d[ind].data = data;
						args4d[ind].coding = coding; args4d[ind].blocksize = blocksize;
						if(pthread_create(&jmed_mt_thread[ind], NULL, jerasure_matrix_dotprod_mt, &args4d[ind])){
							fprintf(stderr, "Error creating data thread/s %d\n", i);
							return 1;
						}
					ind++;
					}
				}  
				for (i = 0; i < ind; i++){
					if(pthread_join(jmed_mt_thread[i], NULL)) {
						fprintf(stderr, "Error joining data thread/s %d\n", i);
						return 2;
					}
				}
			}
			if(edd == 1 && erased[k] == 0){
				tmpids = talloc(int, k);
				for (i = 0; i < k; i++) {
					tmpids[i] = (i < lastdrive) ? i : i+1;
				}
				jerasure_matrix_dotprod(k, w, matrix, tmpids, lastdrive, data, coding, blocksize);
				free(tmpids);
	        }
	        /* recover parity data: - REPAIR process: */
			if(edp == 1){
				for (i = 0; i < m; i++) {
					if (erased[k+i])
						jerasure_matrix_dotprod(k, w, matrix+(i*k), NULL, i+k, data, coding, blocksize);
				}
			}else if(edp>1){
				ind = 0;
				for (i = 0; i < m; i++) {
					if (erased[k+i]) {
						args4p[ind].k = k; args4p[ind].w = w; args4p[ind].matrix_row = matrix+(i*k);
						args4p[ind].src_ids = NULL; args4p[ind].dest_id = i+k; args4p[ind].data = data;
						args4p[ind].coding = coding; args4p[ind].blocksize = blocksize;
						if(pthread_create(&jmep_mt_thread[ind], NULL, jerasure_matrix_dotprod_mt, &args4p[ind])){
							fprintf(stderr, "Error creating coding thread/s %d\n", i);
							return 1;
						}
						ind++;       
					}
				}
				for (i = 0; i < edp; i++){
					if(pthread_join(jmep_mt_thread[i], NULL)) {
						fprintf(stderr, "Error joining coding thread/s %d\n", i);
						return 2;
					}
				}
			}         
		}else{
			fprintf(stderr, "Not a valid coding technique.\n");
			exit(0);
		}
		/* DEcode and Repair stops here. */
		edp0edd0:
		timing_set(&t4);	
	  
		/* Create decoded file */
		sprintf(fname, "%s/Coding/%s_decoded%s", curdir, cs1, extension);
		if (n == 1) {
			fp = fopen(fname, "wb");
			if(fp == NULL){ /* file handler checks were missing in encoder.c */
				printf("Error handling opening the file %s", fname);
			}
		}
		else {
			fp = fopen(fname, "ab");
			if(fp == NULL){ /* file handler checks were missing in encoder.c */
				printf("Error handling opening the file %s", fname);
			}			
		}
		for (i = 0; i < k; i++) {
			if (total+blocksize <= origsize) {
				fwrite(data[i], sizeof(char), blocksize, fp);
				total+= blocksize;
			}
			else {
				for (j = 0; j < blocksize; j++) {
					if (total < origsize) {
						fprintf(fp, "%c", data[i][j]);
						total++;
					}
					else {
						break;
					}
					
				}
			}
		}
		n++;
		fclose(fp);	
		totalsec += timing_delta(&t3, &t4);
	} 
	
	/* Free allocated memory */
	free(cs1);
	free(extension);
	free(fname);
	free(data);
	free(coding);
	free(erasures);
	free(erased);
	if (dm_ids != NULL) free(dm_ids);
	if (decoding_matrix != NULL) free(decoding_matrix);
	
	/* Stop timing and print time */
	timing_set(&t2);
	tsec = timing_delta(&t1, &t2);
	#if defined(MULTIPROCESS)
		/* Here is an extra timing computation/display for multi-process environments */
		printf("%0.6f %0.6f\n", sstart, sstart + totalsec);
	#endif
	printf("Decoding (MB/sec): %0.10f\n", (((double) origsize)/1024.0/1024.0)/totalsec);
	printf("De_Total (MB/sec): %0.10f\n\n", (((double) origsize)/1024.0/1024.0)/tsec);
	
	return 0;
}
	
	
