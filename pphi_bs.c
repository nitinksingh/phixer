#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>

/* Pragma for User defined values */
#define NROW 10
#define TSAMPLE_COUNT 10
#define NBOOTSTRAP 1
/* Compute other pragmas based on above */
#define NBINS (int) floor(0.5*sqrt(TSAMPLE_COUNT))
//#define NBINS 3
#define SAMPLE_COUNT (int) floor(TSAMPLE_COUNT*0.95)
#define NUM_THREADS 1

int compare_floats(const void *a, const void *b);
char *gnu_basename(char *path);

/* Function for computing phi-mixing coefficient
 * @input:
 * 	CSV file where rows are genes (variables) and columns
 * 	are samples (independent observation)
 * @output:
 * 	Pruned graph in MATLAB sparse format
 */
int
main(int argc, char *argv[]) {
	// Input and output file pointers, and vars needed for file reading
	FILE *infp, *outfp;

	size_t read, len = 0;
	char *line = NULL, *tok;

	// Matrices holding phimix computation vals
	float phimix[NROW][NROW], sorted_input[NROW][SAMPLE_COUNT], input_all[NROW][TSAMPLE_COUNT], percentile_edges[NROW][NBINS-1];
	float input[NROW][SAMPLE_COUNT];
	float phimix_pruned_avg[NROW][NROW] = {{0.0}};
	float phimix_pruned[NROW][NROW];

	float x, delta;
	int rank;

	int i, bs, row = 0, col = 0, p, thread_id, nthreads;

	if (argc != 2) {
		printf("Help: %s <input_file>\n", argv[0]);
		return 1;
	}

	/* Open and read input file */
	infp = fopen(argv[1], "r");

	if (infp == 0) {
		printf("Could not open %s input file to read\n", argv[1]);
		return 1;
	}

	// Read each line, tokenize and create input arrays
	row = 0;
	while ((read = getline(&line, &len, infp)) != -1) {
		col = 0;
		tok = strtok(line, ",");
		while (tok != NULL) {
			input_all[row][col] = atof(tok);
			tok = strtok(NULL, ",");
			col = col + 1;
		}
		row = row + 1;
		//printf("row:%d\n", row);
	}
	free(line);
	/* Handle error */
	if (fclose(infp) || ferror(infp)) {
		fprintf(stderr, "Error in file closing or reading!\n");
		exit(1);
	}
	printf("Input file loaded successfully\n");

	// Initialize random seed
	srand(time(0));
	/******************************************
	// Test code: print input matrix
	for (row=0; row<NROW; row++) {
		for (col=0; col<TSAMPLE_COUNT; col++)
			printf("%f\t", input_all[row][col]);
		printf("\n");
	}
	return 0;
	******************************************/
	for (bs=0; bs<NBOOTSTRAP; bs++) {
		/* Generate Randomly sampled inputs from the original input */
		for (i=0; i<SAMPLE_COUNT; i++) {
			// Get a random number and copy that sample
			col = rand() / (RAND_MAX / SAMPLE_COUNT + 1);
			for (row=0; row<NROW; row++) {
				input[row][col] = input_all[row][col];
				sorted_input[row][col] = input_all[row][col];
			}
		}

		// Sort the input data for percentile bining
		for (row=0; row<NROW; row++) {
			qsort(sorted_input[row], SAMPLE_COUNT, sizeof(float), compare_floats);
			// Percentile implementation. Note that its not a generic implementation
			for (p=1; p<NBINS; p++) {
				x = p*(SAMPLE_COUNT+1)/NBINS;
				// Since I know that x shall be always positive, typecast would
				// return same as floor but save floor function call
				rank = (int) x;
				delta = x - rank;

				// Now save the edges of percentile values
				if (rank+1 >= SAMPLE_COUNT) {
					percentile_edges[row][p-1] = sorted_input[row][SAMPLE_COUNT-1];
				} else {
					percentile_edges[row][p-1] = sorted_input[row][rank] + sorted_input[row][rank+1] * delta;
				}
			}
		}


		/******************************************
		// Test code: Print the percentile edges for each data row
		printf("percentile edeges\n");
		for (row=0; row < NROW; row++){
			printf("row %d:",row);
			for(col=0; col < NBINS-1; col++){
				printf("%f ",percentile_edges[row][col]);
			}

			printf("\n");
		}
		printf("\n");
		*****************************************/
		printf("Percentile Edges Calculated\n");

		/******************************************
		// Test code: Print sorted input
		for (row=0; row<NROW; row++) {
			for (col=0; col<SAMPLE_COUNT; col++)
				printf("%f\t", sorted_input[row][col]);
			printf("\n");
		}

		******************************************/

		/******************************************
		// Test code: Print percentile edged sorting
		for (row=0; row<9; row++) {
			for (col=0; col<NBINS-1; col++)
				printf("%f\t", percentile_edges[row][col]);
			printf("\n");
		}
		******************************************/


		/* BINNING code: Bin to compute marginal and joint PDFs */

		// Set number of theads to be used
		omp_set_num_threads(NUM_THREADS);

		// Begin Parallel block
		#pragma omp parallel
		{
			int j=0, k=0, l=0, xbin, ybin, b, rind, cind;
			float x_cord, y_cord;
			float phi[NBINS][NBINS], alpha[NBINS], beta[NBINS];

			nthreads = omp_get_num_threads();
			thread_id = omp_get_thread_num();
			printf("thread id: %d of %d\n", thread_id, nthreads);

			// Parallel for loops
			#pragma omp for schedule(static)
			for (i=0; i<NROW; i++) {
				for (j=0; j<NROW; j++) {
					// Initialize phi, alpha and beta
					for (k=0; k<NBINS; k++) {
						for (l=0; l<NBINS; l++) {
							phi[k][l] = 0;
						}
						alpha[k]=0;
						beta[k]=0;
					}

					// Assign each sample to appropriate bin as defined by percentile_edges
					// Thats how we'll get joint dist (phi) and marginals: alpha and beta
					for (k=0; k<SAMPLE_COUNT; k++) {
						// x-cord, y-cord
						x_cord = input[i][k];
						y_cord = input[j][k];

						// Get bin number for x-cordinate of sample.
						// Default is 0 ie first bin then iterate from right and if x_cord exceeds
						// from percentile edge then assign it to that bin
						xbin = 0;
						for (b=NBINS-2; b >=0; b--) {
							if (x_cord >=  percentile_edges[i][b]) {
								xbin = b+1;
								break;
							}
						}
						// Similarily get bin for y-cordinate of sample
						ybin = 0;
						for (b=NBINS-2; b >=0; b--) {
							if (y_cord >=  percentile_edges[j][b]) {
								ybin = b+1;
								break;
							}
						}

						// Increase the sample freq count in respective bin
						alpha[xbin] += 1.0/SAMPLE_COUNT;
						beta[ybin] += 1.0/SAMPLE_COUNT;
						phi[xbin][ybin] += 1.0/SAMPLE_COUNT;

						// Test code to see if samples are going into right bins
						//if (i==0 && j == 8)
						//	printf("x:%f y:%f %d %d %d \n", x_cord, y_cord, k, xbin, ybin);

					}


					/***********************************************************
					// Test code: Print alpha beta phi
					float asum = 0, bsum = 0, phisum = 0;
					if (i==0 && j==0) {
						printf("alpha\n");
						for (k=0; k<NBINS; k++) {
							printf("%f ", alpha[k]);
							asum += alpha[k];
						}
						printf("\nbeta\n");
						for (l=0; l<NBINS; l++) {
							printf("%f ", beta[l]);
							bsum += beta[l];
						}
						printf("\nphi\n");

						for (rind=0; rind < NBINS; rind++){
							for(cind=0; cind < NBINS; cind++){
								printf("%f ",phi[rind][cind]);
								phisum += phi[rind][cind];

							}
							printf("\n");
						}
						printf("i: %d, j: %d, asum: %f, bsum; %f, ,phisum; %f\n", i,j, asum, bsum, phisum);
					}
					************************************************************/


					// Now phi, alpha and beta are ready. Compute phi-mix
					float abs_colsum, temp_phi = 0.0, phi_coeff = 0.0;

					if (i==j) {
						phi_coeff = 1.0;
						phimix[i][j] =  phi_coeff;
						phimix_pruned[i][j] = phi_coeff;
						continue;
					}


					for (cind=0; cind < NBINS; cind++) {
						abs_colsum = 0.0;
						// Compute the sum of each column of phi - alpha*beta.
						for (rind=0; rind < NBINS; rind++){
							abs_colsum += fabs(phi[rind][cind]-alpha[rind]*beta[cind]);
						}

					 	temp_phi = 0.5 * (abs_colsum/beta[cind]);
						// Take max over col wise phimix
						if (temp_phi > phi_coeff){
							phi_coeff = temp_phi;
						}
					}

					// Assign computed phi-coefficeint to corresponding value
					// Also, initialize phimix pruning matrix as well
					phimix[i][j] =  phi_coeff;
					phimix_pruned[i][j] = phi_coeff;

				}
				//printf("row: %d finished\n", i+1);
			}

		// Parallel block ends here
		}
		printf("Finished phimixing computation.\n");
		printf("Beginning pruning block....\n");
		// Begin Parallel block
		#pragma omp parallel
		{
			int j=0, k=0;
			float ij=0, kj=0, ik=0, min=0;
			nthreads = omp_get_num_threads();
			thread_id = omp_get_thread_num();
			printf("thread id: %d of %d\n", thread_id, nthreads);

			// Parallel for loops
			#pragma omp for schedule(static)
			for (i=0; i<=NROW-1; i++) {
				for (j=0; j<=NROW-1; j++) {
					for (k=0; k<=NROW-1; k++){
						ij = phimix[i][j];
						ik = phimix[i][k];
						kj = phimix[k][j];

						min = (ik < kj) ? ik : kj;

						if (ij < min) {
							phimix_pruned[i][j] = 0;
							break;
						}

					}
				}

				//printf("row: %d finished\n", i);
			}

		// Parallel block ends here
		}
		/***********************************************
		//Test code: print epsilon values
		for ( row=0; row < NROW; row++){
			for(col=0; col < NROW; col++){
				printf("%f ",epsilon[row][col]);
			}
			printf("\n");
		}

		printf("\n");
		***********************************************/
		for (row=0; row<=NROW-1; row++) {
			for (col=0; col<=NROW-1; col++) {
				phimix_pruned_avg[row][col] += phimix_pruned[row][col]/NBOOTSTRAP;
			}
		}

		printf("Phimixing pruning completed for bs: %d\n", bs);
	}

 	// Commented phimix file writing
	/***********************************************
	// Open output files
	char phimix_file[100];
	strcpy(phimix_file, "phimix_");
	strcat(phimix_file, gnu_basename(argv[1]));
	outfp = fopen(phimix_file, "w");

	if (outfp == 0) {
		printf("Could not open %s output file to write\n", phimix_file);
		return 1;
	}

	printf("Writing phimix values to file...\n");
	// Dump the  phi matrix into MATLAB sparse format output file
	for (row=0; row<=NROW-1; row++) {
		for (col=0; col<=NROW-1; col++) {
			fprintf(outfp, "%d\t%d\t%f\n", row+1, col+1, phimix[row][col]);
		}
	}

	fclose(outfp);
	***********************************************/


	/* Dump the  pruned phi matrix into MATLAB sparse format output file */
	char pruned_file[100], binsize[30];
	strcpy(pruned_file, "pruned_nw_bs");
	sprintf(binsize, "%d_bins%d_", NBOOTSTRAP, NBINS);
	strcat(pruned_file, binsize);
	strcat(pruned_file, gnu_basename(argv[1]));
	outfp = fopen(pruned_file, "w");
	if (outfp == 0) {
		printf("Could not open %s output file to write\n", pruned_file);
		return 1;
	}

	/* Currently a_ij is a(i|j) i.e j-->i; transposinig so that a_ij is i-->j */
	printf("Writing PRUNED phimixing matrix to file...\n");
	for (row=0; row<=NROW-1; row++) {
		for (col=0; col<=NROW-1; col++) {
			if (phimix_pruned_avg[row][col] > 0) {
				fprintf(outfp, "%d,%d,%f\n", col+1, row+1, phimix_pruned_avg[row][col]);
			}
		}
	}
	fclose(outfp);

	printf("Finished writing output file\\s...\n");

	return 0;
}

/* gnu basename function implementation -- copied */
char *
gnu_basename(char *path)
{
    char *base = strrchr(path, '/');
    return base ? base+1 : path;
}

/*
 * Compare function needed for qsort copied from qsort documentation
 */
int
compare_floats(const void *a, const void *b)
     {
       const float *da = (const float *) a;
       const float *db = (const float *) b;

       return (*da > *db) - (*da < *db);
     }

/*
 * Percentile implementation reference
 * http://commons.apache.org/math/apidocs/org/apache/commons/math/stat/descriptive/rank/Percentile.html
*/
