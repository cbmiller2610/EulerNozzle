#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <signal.h>

static volatile sig_atomic_t keep_running=1;

int main (void)
{
  return 0;
}

static void sig_handler (int _)
{
  (void)_;
  keep_running = 0;
}

int heat_solver (double *U, int xnodes, double dx, int tnodes, double dt, double Nu, double r, double C)
{
  signal(SIGINT, sig_handler);
  while (keep_running) 
  {
    int i, j, status;
    gsl_matrix_view U_view = gsl_matrix_view_array(U, tnodes, xnodes);
    gsl_vector_view U_old_row, U_new_row;
    double U_old_back, U_old_mid, U_old_fwd;
    gsl_vector *diag = gsl_vector_alloc(xnodes);
    gsl_vector *e = gsl_vector_alloc(xnodes-1);
    gsl_vector *f = gsl_vector_alloc(xnodes-1);
    gsl_vector *b = gsl_vector_alloc(xnodes);
    
    /*Initialize tridiagonal coefficient vectors for Ax=b
    Diagonal vector*/
    for (i=0; i<xnodes; i++) {
      gsl_vector_set(diag, i, 2+2*r);
    }
    gsl_vector_set(diag, xnodes-1, 2+2*r+2*r*dx*Nu);

    /*above diagonal vector*/
    gsl_vector_set(e, 0, -2*r);
    for (i=1; i<xnodes-1; i++) {
      gsl_vector_set(e, i, -r);
    }
    
    /*below diagonal vector*/
    for (i=0; i<xnodes-2; i++) {
      gsl_vector_set(f, i, -r);
    }
    gsl_vector_set(e, xnodes-2, -2*r);

    /*March in time updating the b vector using values from the current timestep
     *  and calculating x (U vector) with GSL tridiagonal solver*/
    for (i=0; i<tnodes-1; i++) {
      U_old_row = gsl_matrix_row(&U_view.matrix, i);
      U_new_row = gsl_matrix_row(&U_view.matrix, i+1);
      /*Set b_0*/
      U_old_mid = gsl_vector_get(&U_old_row.vector, 0);
      U_old_fwd = gsl_vector_get(&U_old_row.vector, 1);
      /*Zero Neumann BC*/
      /*gsl_vector_set(b, 0, (2-2*r)*U_old_mid + 2*r*U_old_fwd);*/
      /*Constant Neumann BC*/
      gsl_vector_set(b, 0, (2-2*r)*U_old_mid + 2*r*U_old_fwd-4*r*dx*C);
      /*Set internal nodes*/
      for (j=1; j<xnodes-1; j++) {
	U_old_back = gsl_vector_get(&U_old_row.vector, j-1); 
	U_old_mid = gsl_vector_get(&U_old_row.vector, j);
	U_old_fwd = gsl_vector_get(&U_old_row.vector, j+1);
	gsl_vector_set(b, j, (2-2*r)*U_old_mid+r*U_old_fwd+r*U_old_back);
      }
      /*Set b_end*/
      U_old_back = gsl_vector_get(&U_old_row.vector, xnodes-2);
      U_old_mid = gsl_vector_get(&U_old_row.vector, xnodes-1);
      gsl_vector_set(b, xnodes-1, (2-2*r-2*r*dx*Nu)*U_old_mid+2*r*U_old_back);
      /*Tridiagonal Solve, storing results to new U vector*/
      status = gsl_linalg_solve_tridiag(diag, e, f, b, &U_new_row.vector);
      if (status != 0) {
	printf("Something went wrong!!!\n");
      }
    }
    gsl_vector_free(e);
    gsl_vector_free(f);
    gsl_vector_free(b);
    gsl_vector_free(diag);
    keep_running=0;
  }

  return EXIT_SUCCESS;
}
