####################################
# About Phixer
####################################
Phixer is an algorithm to generate Gene Interaction Network
from gene-expression data. The resulting GIN graph is directed
and may contain cycles.

The core of the Phixer is implemented in C which depends
on OpenMP library for carrying out computations in parallel.

A high level description of the algorithm is given below.

1. Read input gene-expression CSV files. Rows are genes and columns
are number of samples. The input file should be preprocessed, and should
not contain NA or missing values. In other words, the input CSV file 
should be a m-by-n data matrix (m genes, n samples). There should
not be any header row or gene-id column.

2. For each resampled bootstrap runs, compute phi-mixing coefficient in parallel.

3. Prune the resulting GIN graph in parallel.

4. Write the pruned GIN in Matlab sparse format for downstream processing.

5. Perform thresholding in Matlab.


####################################
# Instructions for running Phixer 
####################################
The following instructions are meant for Linux operating system
with GCC and OpenMP installed.

# Step 1 
Edit pphi_bs.c Set NROW, TSAMPLE_COUNT, BOOTSTRAPS, NUM_THREADS

# Step 2: Compile from bash shell
```bash
gcc -Wall pphi_bs.c -fopenmp -o phixer.out
```
# Step 3: Execute
Increase memory and stack size. Depending on machine hardwars/OS settings,
this step may be needed for practical sized 
problems i.e number of genes greater than 15,000. 
The actual values will vary from machine to machine.
```bash
ulimit -s 1300000000

export GOMP_STACKSIZE=2000000
```
Compute phi-mixing coefficients, do prunning and write prunned graph.
Note that output file shall be created in the directory where executable
phixer.out resides.

```bash
time ./phixer.out input_data_csv_file 
```
# Step 4: Threshold (MATLAB)
run phixer_threshold.m that takes the output of step 3 as an input.

####################################
# Demo
####################################
# Step 1
The current pphi_bs.c has been edited for the demo. 

# Step 2
```bash
gcc -Wall pphi_bs.c -fopenmp -o phixer_demo.out
```
# Step 3
```bash
./phixer_demo.out test_100_585.txt
```
# Step 4 -- FROM Matlab
```matlab
nw = phixer_threshold(step3_output_file);
```
