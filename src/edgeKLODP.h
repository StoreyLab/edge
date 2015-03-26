#include <stdio.h>/* --------- printf */
#include <stdlib.h>/* -------- malloc(), free(), etc. */
#include <math.h>/* ---------- log(), exp(), etc. */
#include <R.h>/* ------------- R functions */
#include <Rinternals.h>

#define NR_END 1
#define FREE_ARG void*

/***********************************************************************
  EDGE-specific functions
***********************************************************************/
void odpScoreCluster(double *, double *, double *, int *, int *, int *, int *, int *, double *);
void kldistance(double *, double *, double *, double *, int *, int *, int *, double *);

/***********************************************************************
  utility functions 
***********************************************************************/
  void sortQK(int, int, int, double *);
void swapQK(int, int, double *);
double *vector(int, int);
void free_vector(double *, int, int);
int *ivector(int, int);
void free_ivector(int *, int, int);
double **matrix(int, int, int, int);
void free_matrix(double **, int, int, int, int);
int **imatrix(int, int, int, int);
void free_imatrix(int **, int, int, int, int);
