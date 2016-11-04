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
technique, w, and packetsize.  It creates k+m files from 
the original file so that k of these files are parts of 
the original file and m of the files are encoded based on 
the given coding technique. The format of the created files 
is the file name with "_k#" or "_m#" and then the extension.  
(For example, inputfile test.txt would yield file "test_k1.txt".)
*/

/*********************************************************************************
	Beta Version (Revision 1.0): Suayb S. Arslan and Hoa Le 
	Jerasure 2.x Multi-Thread Encoder Implementation using shared-memory method.
*/
/*
Preface: The original library has different set of erasure coding techniques. 
Based on our testing and simulations, we figured the coding technqiue "reed_sol_van"
performed the best among others and it turns out that it construction aligns more 
friendly with the MT-implementation. Major changes with respect to original encoder.c
are detailed below with each code line item.

There are multiple ways of implementing multi-threading for Jerasure. We took 
a shared-memory (no-mutex) approach. In this approach, threads share the same buffer 
space. Since threads are only allowed to read from but not write to the buffer, we 
have not used mutexes and hence secured good performance. 
	                     
Reference Doc: See PPT file Jerasure2,0_setup_math_perf at Quantum/Box repositories. 
**********************************************************************************/

#include <assert.h>
#include <pthread.h> //important to include.
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>
#include <gf_rand.h>
#include <unistd.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "cauchy.h"
#include "liberation.h"
#include "galois.h"

/* #define MULTIPROCESS */
enum Coding_Technique {Reed_Sol_Van, Reed_Sol_R6_Op, Cauchy_Orig, Cauchy_Good, Liberation, Blaum_Roth, Liber8tion, RDP, EVENODD, No_Coding};
/* We decided to implement reed_sol_van only due to performance reasons */
/* Global variables */
int readins, n;
/* No signal handler is used */

/* No need for function prototypes:*/
/* 1. is_prime is NOT used by reed_sol_van, */
/* 2. No need for ctrl handler. */
double get_time_usec() /* This is for debugging purposes - time optimization (not needed for production) */
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return 1000000.0 * tv.tv_sec + tv.tv_usec;
}

/* This is for multi-threading implementation (original encoder does not have it) */
/* Better to put into a header file. We decided to make minimum changes though */
typedef struct{ 
    int k;
    int w;
    int *matrix_row;
    int dest_id;
    char** data;
    char** coding;
    int blocksize;
}jinput_args;

/* The following function is a modified version of jerasure_matrix_dotprod_mt 
to utilize multi-thread capabilities of modern CPUs.*/
void *jerasure_matrix_dotprod_mt(void *x_void_ptr){ 
    jinput_args *x_ptr = (jinput_args *)x_void_ptr;
    
    int k = x_ptr->k;
    int w = x_ptr->w;
    int dest_id = x_ptr->dest_id;
    int *matrix_row = x_ptr->matrix_row;
    char** data_ptrs = x_ptr->data;
    char** coding_ptrs = x_ptr->coding;
    int blocksize = x_ptr->blocksize;
    
    int init = 0;
    char *dptr, *sptr;
    int i;
    int j = dest_id-k;
    
    dptr = coding_ptrs[j];
    for (i = 0; i < k; i++) {
        sptr = data_ptrs[i];
        if (matrix_row[i+j*k] == 1) {
            if (init == 0) {
                memcpy(dptr, sptr, blocksize);
                init = 1;
            }else{
                galois_region_xor(sptr, dptr, blocksize);
            } 
        }
        if (matrix_row[i+j*k] > 1) {
             switch (w) {
                case 8:  galois_w08_region_multiply(sptr, matrix_row[i+j*k], blocksize, dptr, init); break;
                case 16: galois_w16_region_multiply(sptr, matrix_row[i+j*k], blocksize, dptr, init); break;
                case 32: galois_w32_region_multiply(sptr, matrix_row[i+j*k], blocksize, dptr, init); break;
            }
            init = 1;   
        }
    }
    
    return NULL;
}


int jfread(void *ptr, int size, int nmembers, FILE *stream)
{
  if (stream != NULL) return fread(ptr, size, nmembers, stream);

  MOA_Fill_Random_Region(ptr, size);
  return size;
}


int main (int argc, char **argv) {
	FILE *fp, *fp2;					// file pointers
	char *block;					// padding file
	unsigned int size, newsize;		// size of file and temp size - we made it unsigned to cover large files.
	struct stat status;				// finding file size

	
	enum Coding_Technique tech;		// coding technique (parameter)
	int k, m, w, packetsize;		// parameters
	int buffersize;					// paramter
	int i;							// loop control variables
	int blocksize;					// size of k+m files
	unsigned int total;				/* we made it unsigned to cover large files. */
	int extra;

	/* Jerasure Arguments */
	char **data;				
	char **coding;
	int *matrix;
	/* We do not need bitmatrix for reed_sol_van */
	/* We do not need scheduling for reed_sol_van */
	
	/* Creation of file name variables */
	char temp[5];  
	char *s1, *s2, *extension;
	char *fname;
	int md;
	char *curdir;

	/* Timing variables */
	double t1, t2, t3, t4;
	struct timeval tstart, tfinish; /* extra timing parameters for measuring performance */
	double tsec;
	double totalsec;
	#if defined(MULTIPROCESS)
		double sstart; /* measuring time in multi-process environment */
	#endif
	
	/* Find buffersize */
	int up, down;
	
	/* Start timing */
	t1 = (double)get_time_usec();
	
	//Initialization
	totalsec = 0.0;
	matrix = NULL;
	/* We do not need bitmatrix for reed_sol_van */
	/* We do not need scheduling for reed_sol_van */
	
	/* Error check Arguments*/
	if (argc != 7) {
		fprintf(stderr,  "usage: inputfile k m w packetsize buffersize\n");
		/* We have shortened the comments due to using only reed_sol_van coding technique */
		exit(0);
	}
	/* Conversion of parameters and error checking */
	if (sscanf(argv[2], "%d", &k) == 0 || k <= 0) {
		fprintf(stderr,  "Invalid value for k\n");
		exit(0);
	}
	if (sscanf(argv[3], "%d", &m) == 0 || m < 0) {
		fprintf(stderr,  "Invalid value for m\n");
		exit(0);
	}
	if (sscanf(argv[4],"%d", &w) == 0 || w <= 0) { /* we have one less parameter in the argument */
		fprintf(stderr,  "Invalid value for w.\n");
		exit(0);
	}
	if (argc == 5) { /* we have one less parameter in the argument */
		packetsize = 0;	
	}
	else {
		if (sscanf(argv[5], "%d", &packetsize) == 0 || packetsize < 0) {
			fprintf(stderr,  "Invalid value for packetsize.\n");
			exit(0);
		}
	}
	if (argc != 7) { /* we have one less parameter in the argument */
	    buffersize = 0;
	}
	else{ 
		if (sscanf(argv[6], "%d", &buffersize) == 0 || buffersize < 0) {
			fprintf(stderr, "Invalid value for buffersize\n");
			exit(0);
		}
	}
	/* Determine proper buffersize by finding the closest valid buffersize to the input value  */
	if (buffersize != 0) {
		if (packetsize != 0 && buffersize%(sizeof(long)*w*k*packetsize) != 0) { 
			up = buffersize;
			down = buffersize;
			while (up%(sizeof(long)*w*k*packetsize) != 0 && (down%(sizeof(long)*w*k*packetsize) != 0)) {
				up++;
				if (down == 0) {
					down--;
				}
			}
			if (up%(sizeof(long)*w*k*packetsize) == 0) {
				buffersize = up;
			}
			else {
				if (down != 0) {
					buffersize = down;
				}
			}
		}
		else if (packetsize == 0 && buffersize%(sizeof(long)*w*k) != 0) {
			up = buffersize;
			down = buffersize;
			while (up%(sizeof(long)*w*k) != 0 && down%(sizeof(long)*w*k) != 0) {
				up++;
				down--;
			}
			if (up%(sizeof(long)*w*k) == 0) {
				buffersize = up;
			}
			else {
				buffersize = down;
			}
		}
	}
	
	/* Setting of coding technique and error checking */
	// We only use Vandermonde matrices as they provide sufficiently good performance.
	
	
	
	
	
	tech = Reed_Sol_Van; 
	if (w != 8 && w != 16 && w != 32) {
		fprintf(stderr,  "w must be one of {8, 16, 32}\n");
		exit(0);
	}
	/* No need for method variable when we use single code choice:reed_sol_van */
	/* Get current working directory for construction of file names */
	curdir = (char*)malloc(sizeof(char)*1000);	
	assert(curdir == getcwd(curdir, 1000));

	if (argv[1][0] != '-') {

		/* Open file and error check */
		fp = fopen(argv[1], "rb");
		if (fp == NULL) {
			fprintf(stderr,  "Unable to open file.\n");
			exit(0);
		}
	
		/* Create Coding directory */
		i = mkdir("Coding", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Coding directory.\n");
			exit(0);
		}
	
		/* Determine original size of file */
		stat(argv[1], &status);	
		size = status.st_size;
        } else {
        	if (sscanf(argv[1]+1, "%d", &size) != 1 || size <= 0) {
                	fprintf(stderr, "Files starting with '-' should be sizes for randomly created input\n");
			exit(1);
		}
        	fp = NULL;
		MOA_Seed(time(0));
	}

	newsize = size;
	
	/* Find new size by determining next closest multiple */
	if (packetsize != 0) {
		if (size%(k*w*packetsize*sizeof(long)) != 0) {
			while (newsize%(k*w*packetsize*sizeof(long)) != 0) 
				newsize++;
		}
	}
	else {
		if (size%(k*w*sizeof(long)) != 0) {
			while (newsize%(k*w*sizeof(long)) != 0) 
				newsize++;
		}
	}
	
	if (buffersize != 0) {
		while (newsize%buffersize != 0){ 
			newsize++;
		}
	}
	
	
	/* Determine size of k+m files */
	blocksize = newsize/k;

	/* Allow for buffersize and determine number of read-ins */
	if (size > buffersize && buffersize != 0) {
		/* We figured the conditional statement here has no effect */
		readins = newsize/buffersize;
		block = (char *)malloc(sizeof(char)*buffersize);
		blocksize = buffersize/k;
	}
	else {
		readins = 1;
		buffersize = size;
		block = (char *)malloc(sizeof(char)*newsize);
	}
	
	/* Break inputfile name into the filename and extension */	
	s1 = (char*)malloc(sizeof(char)*(strlen(argv[1])+20));
	s2 = strrchr(argv[1], '/');
	if (s2 != NULL) {
		s2++;
		strcpy(s1, s2);
	}
	else {
		strcpy(s1, argv[1]);
	}
	s2 = strchr(s1, '.');
	if (s2 != NULL) {
		extension = strdup(s2);
		*s2 = '\0';
	} else {
		extension = strdup("");
	}
	
	/* Allocate for full file name */
	fname = (char*)malloc(sizeof(char)*(strlen(argv[1])+strlen(curdir)+20));
	sprintf(temp, "%d", k);
	md = strlen(temp);
	
	/* Allocate data and coding */
	data = (char **)malloc(sizeof(char*)*k);
	coding = (char **)malloc(sizeof(char*)*m);
	for (i = 0; i < m; i++) {
		coding[i] = (char *)malloc(sizeof(char)*blocksize);
                if (coding[i] == NULL) { perror("malloc"); exit(1); }
	}
	
	
	
	/* Create coding matrix or bitmatrix and schedule */
	t3 = (double)get_time_usec();
	matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
	t4 = (double)get_time_usec();	
	if ((t4 - t3) < 0) /* accurate calculation of timing ( with sanity check ) */
		totalsec += (double)(unsigned)(~0);
	else
		totalsec += t4 - t3;
    
	/* Read in data until finished */ 
	n = 1;
	total = 0;
	/* input values for MT-implementation */
	jinput_args args[m];
	pthread_t jme_mt_thread[m];
	
	while (n <= readins) {
		/* Check if padding is needed, if so, add appropriate 
		   number of zeros */
		if (total < size && total+buffersize <= size) {
			total += jfread(block, sizeof(char), buffersize, fp);
		}
		else if (total < size && total+buffersize > size) {
			extra = jfread(block, sizeof(char), buffersize, fp);
			for (i = extra; i < buffersize; i++) {
				block[i] = '0';
			}
		}
		else if (total == size) {
			for (i = 0; i < buffersize; i++) {
				block[i] = '0';
			}
		}
	
		/* Set pointers to point to file data */
		for (i = 0; i < k; i++) {
			data[i] = block+(i*blocksize);
		}
		
		t3 = (double)get_time_usec();	
		/* Since we use a subfunction of jerasure_matrix_encode function in multi-threaded mode, we
		needed to change this part of the code using a for loop. */
		if (n==1){
		    gettimeofday(&tstart, NULL);
		    #if defined(MULTIPROCESS)
			sstart = (double) tstart.tv_sec + ((double) tstart.tv_usec) / 1000000;
		    #endif
		}
		for (i = 0; i < m; i++){
			args[i].k = k; args[i].w = w;args[i].matrix_row = matrix;args[i].dest_id = k+i;
			args[i].data = data;args[i].coding = coding;args[i].blocksize = blocksize;
			if(pthread_create(&jme_mt_thread[i], NULL, jerasure_matrix_dotprod_mt, &args[i])){
				fprintf(stderr, "Error creating thread %d\n", i);
				return 1;
				}
		}	
		for (i = 0; i < m; i++){
			if(pthread_join(jme_mt_thread[i], NULL)) {
				fprintf(stderr, "Error joining thread %d\n", i);
				return 2;
			}
		}
		gettimeofday(&tfinish, NULL);
		t4 = (double)get_time_usec();	
		if ((t4 - t3) < 0)
			totalsec += (double)(unsigned)(~0);
		else
			totalsec += t4 - t3;
		
		/* Write data and encoded data to k+m files */
		/* This change does not affect the MT performance. We make a single for loop for reducing the 
		number of instructions. Particularly for big files and small buffersize, since this loop is 
		inside the while loop this can lead to a change in the performance.*/
		for	(i = 1; i <= k+m; i++) {
		    if (fp == NULL) {
		    	bzero(data[i-1], blocksize);
		    } else {
				if (i <= k){
					sprintf(fname, "%s/Coding/%s_k%0*d%s", curdir, s1, md, i, extension);
				}else{
					sprintf(fname, "%s/Coding/%s_m%0*d%s", curdir, s1, md, i - k, extension);
				}
				if (n == 1) {
					fp2 = fopen(fname, "wb");
				}
				else {
					fp2 = fopen(fname, "ab");
				}

				if (i <= k){
					fwrite(data[i-1], sizeof(char), blocksize, fp2);
				}else{
					fwrite(coding[i - k - 1], sizeof(char), blocksize, fp2);
				}
				fclose(fp2);
			}
		}
		n++;		
	}
	/* Create metadata file */
	if(fp != NULL){
		sprintf(fname, "%s/Coding/%s_meta.txt", curdir, s1);
		fp2 = fopen(fname, "wb");
		fprintf(fp2, "%s\n", argv[1]); 
		fprintf(fp2, "%u\n", size); /* u% for unsigned numbers */
		fprintf(fp2, "%d %d %d %d %d\n", k, m, w, packetsize, buffersize);
		fprintf(fp2, "reed_sol_van\n");
		fprintf(fp2, "%d\n", tech);
		fprintf(fp2, "%d\n", readins);
		fclose(fp2);
	}
    
	/* Free allocated memory */
	free(s1);
	free(fname);
	free(block);
	free(curdir);
	
	/* Calculate rate in MB/sec and print */
	t2 = (double)get_time_usec();
	tsec = t2 - t1;	
 	#if defined(MULTIPROCESS)	
		/* Here is an extra timing computation/display for multi-process environments */
		printf("%0.6f %0.6f\n", sstart, sstart + (double)totalsec/1000000);
		//printf("%lu %lu", 1000000 * tstart.tv_sec + tstart.tv_usec, 1000000 * tfinish.tv_sec + tfinish.tv_usec);
 	#endif	
	printf("Encoding (MB/sec): %0.10f\n", 1000000.0*(((double) size)/1024.0/1024.0)/totalsec);
	printf("En_Total (MB/sec): %0.10f\n", 1000000.0*(((double) size)/1024.0/1024.0)/tsec);
    
	return 0;
}
