#include "edgeKLODP.h"

/********************************************************************************************
functions for KLODP:
  odpScoreCluster: compute sum of normal densities to be used as numerator or denominator in score
  with c.member.
********************************************************************************************/

void odpScoreCluster(double *sumDat, double *mu, double *sigma, int *m, int *n, int *p, int *null, int *cluster, double *scr) {
  int i, j, g;
  double *first, *middle;

  /* if alternative component, set up a couple of vectors */
  /* allocate memory */
   first = vector(0, *m - 1);

  /* initialize to zero */
	for(i = 0; i < *m; i++)
      first[i] = 0.0;

  if(*null == 0) {
    /* allocate memory */
    middle = vector(0, *p - 1);

    /* initialize to zero */
	for(i = 0; i < *p; i++) {
      middle[i] = 0.0;
    }
  }

  for(i = 0; i < *m; i++) {
	for(j=0; j< *n ; j++){
		  first[i] += sumDat[j + i * *n]*sumDat[j + i * *n];
	}
  }

  for(i = 0; i < *m; i++) {
    scr[i] = 0.0;

    for(g = 0; g < *p; g++) {			/* g scans genes */
      /* alternative component */
      if(*null == 0) {
         /* middle[j] += 2 * sumDat[i + (l + 1) * *m] * mu[g + l * *m];*/
		    for(j=0; j< *n ; j++){
			    middle[g] += 2 * sumDat[j + i * *n]*sumDat[j + g * *n + *n * *m];
		    }
  		  /*last[g] += nGrp[l] * mu[g + l * *m] * mu[g + l * *m];*/
		    scr[i] += pow(1 / sigma[g], *n) * exp(-0.5 / sigma[g] / sigma[g] * (first[i] - middle[g] + mu[g])) * cluster[g];
      } else /* null component */
        scr[i] += pow(1 / sigma[g], *n) * exp(-0.5 / sigma[g] / sigma[g] * first[i]) * cluster[g];
    }
	    /* reset vectors to zero, if necessary */
    if(*null == 0) {
      for(g = 0; g < *p; g++) {
        middle[g] = 0.0;
      }
    }

  }

  /* free memory, if necessary */
    free_vector(first, 0, *m - 1);

  if(*null == 0) {
    free_vector(middle, 0, *p - 1);
  }
}

void kldistance(double *centerFit, double *centerVar, double *fit, double *var, int *m, int *nc, int *n, double *kldd) {
  int i, j, l;
  double sum;

  for(i = 0; i < *m; i++) {
    for(j = 0; j < *nc; j++) {			/* l scans clusters */
	kldd[j + i* *nc] = 0.0;
	sum = 0.0;
	for(l=0; l< *n ; l++){
		sum += pow((centerFit[l + j* *n]-fit[l + i* *n]),2);
	}
		kldd[j + i* *nc] =  (sum * (1 / centerVar[j] + 1 / var[i]))/2 + *n * (centerVar[j] / var[i] + var[i] / centerVar[j])/2 - *n;
    }
  }
}



/* quicksort routine */
void sortQK(int low, int high, int n, double *w) {
  if(low < high) {
    int lo = low, hi = high + 1;
    double elem = w[low];
    for (;;) {
      while ((lo < n) && (w[++lo] < elem));
      while ((hi >= 0) && (w[--hi] > elem));
      if (lo < hi) swapQK(lo, hi, w);
      else break;
    }

    swapQK(low, hi, w);
    sortQK(low, hi - 1, n, w);
    sortQK(hi + 1, high, n, w);
  }
}

/* swap function for use with sortQK() */
void swapQK(int i, int j, double *w) {
  double tmp = w[i];

  w[i] = w[j];
  w[j] = tmp;
}

/* allocate a int vector with subscript range v[nl...nh] */
int *ivector(int nl, int nh) {
  int *v;

  v = (int *) malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
  if(!v) Rprintf("\n allocation failure in ivector()\n");
  return v - nl + NR_END;
}

/* free a int vector allocated with ivector() */
void free_ivector(int *v, int nl, int nh) {
  free((FREE_ARG) (v + nl - NR_END));
}

/* allocate a int matrix with subscript ranges m[nrl...nrh][ncl...nch] */
int **imatrix(int nrl, int nrh, int ncl, int nch) {
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc((size_t)((nrow + NR_END) * sizeof(int*)));
  if(!m) Rprintf("%s", "allocation fialure\n");

  m += NR_END;
  m -= nrl;

  /* set pointer to rows */
  m[nrl] = (int *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
  if(!m[nrl]) Rprintf("%s", "allocation fialure\n");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
  return m;
}

/* free int matrix allocated with imatrix() */
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch) {
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

/* allocate a double matrix with subscript ranges m[nrl...nrh][ncl...nch] */
double **matrix(int nrl, int nrh, int ncl, int nch) {
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t)((nrow + NR_END) * sizeof(double*)));
  if(!m) Rprintf("%s", "allocation fialure\n");

  m += NR_END;
  m -= nrl;

  /* set pointer to rows */
  m[nrl] = (double *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if(!m[nrl]) Rprintf("%s", "allocation fialure\n");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
  return m;
}

/* free double matrix allocated with matrix() */
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch) {
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

/* allocate a double vector with subscript range v[nl...nh] */
double *vector(int nl, int nh) {
  double *v;

  v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
  if(!v) Rprintf("\n allocation failure in vector()\n");
  return v - nl + NR_END;
}

/* free double vector allocated with vector() */
void free_vector(double *v, int nl, int nh) {
  free((FREE_ARG) (v + nl - NR_END));
}
