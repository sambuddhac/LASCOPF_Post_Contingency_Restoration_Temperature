//solver.c
/* Produced by CVXGEN, 2016-05-03 19:35:14 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */
#include "gensolverFirst.h"
#include <array>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;
double GensolverFirst::eval_gap(void) {
  int i;
  double gap;
  gap = 0;
  for (i = 0; i < 26; i++)
    gap += work.z[i]*work.s[i];
  return gap;
}
void GensolverFirst::set_defaults(void) {
  settings.resid_tol = 1e-6;
  settings.eps = 1e-4;
  settings.max_iters = 25;
  settings.refine_steps = 1;
  settings.s_init = 1;
  settings.z_init = 1;
  settings.debug = 0;
  settings.verbose = 0;
  settings.verbose_refinement = 0;
  settings.better_start = 1;
  settings.kkt_reg = 1e-7;
}
void GensolverFirst::setup_pointers(void) {
  work.y = work.x + 23;
  work.s = work.x + 23;
  work.z = work.x + 49;
  vars.Pg = work.x + 0;
  vars.PgNext = work.x + 1;
  vars.Thetag = work.x + 12;
}
void GensolverFirst::setup_indexing(void) {
  setup_pointers();
}
void GensolverFirst::set_start(void) {
  int i;
  for (i = 0; i < 23; i++)
    work.x[i] = 0;
  for (i = 0; i < 0; i++)
    work.y[i] = 0;
  for (i = 0; i < 26; i++)
    work.s[i] = (work.h[i] > 0) ? work.h[i] : settings.s_init;
  for (i = 0; i < 26; i++)
    work.z[i] = settings.z_init;
}
double GensolverFirst::eval_objv(void) {
  int i;
  double objv;
  /* Borrow space in work.rhs. */
  multbyP(work.rhs, work.x);
  objv = 0;
  for (i = 0; i < 23; i++)
    objv += work.x[i]*work.rhs[i];
  objv *= 0.5;
  for (i = 0; i < 23; i++)
    objv += work.q[i]*work.x[i];
  objv += params.c0[0]+work.quad_608032620544[0]+work.quad_900708769792[0]+work.frac_121674190848*(work.quad_324888432640[0]+work.quad_331507691520[0]);
  return objv;
}
void GensolverFirst::fillrhs_aff(void) {
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 23;
  r3 = work.rhs + 49;
  r4 = work.rhs + 75;
  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, work.y);
  multbymGT(work.buffer, work.z);
  for (i = 0; i < 23; i++)
    r1[i] += work.buffer[i];
  multbyP(work.buffer, work.x);
  for (i = 0; i < 23; i++)
    r1[i] -= work.buffer[i] + work.q[i];
  /* r2 = -z. */
  for (i = 0; i < 26; i++)
    r2[i] = -work.z[i];
  /* r3 = -Gx - s + h. */
  multbymG(r3, work.x);
  for (i = 0; i < 26; i++)
    r3[i] += -work.s[i] + work.h[i];
  /* r4 = -Ax + b. */
  multbymA(r4, work.x);
  for (i = 0; i < 0; i++)
    r4[i] += work.b[i];
}
void GensolverFirst::fillrhs_cc(void) {
  int i;
  double *r2;
  double *ds_aff, *dz_aff;
  double mu;
  double alpha;
  double sigma;
  double smu;
  double minval;
  r2 = work.rhs + 23;
  ds_aff = work.lhs_aff + 23;
  dz_aff = work.lhs_aff + 49;
  mu = 0;
  for (i = 0; i < 26; i++)
    mu += work.s[i]*work.z[i];
  /* Don't finish calculating mu quite yet. */
  /* Find min(min(ds./s), min(dz./z)). */
  minval = 0;
  for (i = 0; i < 26; i++)
    if (ds_aff[i] < minval*work.s[i])
      minval = ds_aff[i]/work.s[i];
  for (i = 0; i < 26; i++)
    if (dz_aff[i] < minval*work.z[i])
      minval = dz_aff[i]/work.z[i];
  /* Find alpha. */
  if (-1 < minval)
      alpha = 1;
  else
      alpha = -1/minval;
  sigma = 0;
  for (i = 0; i < 26; i++)
    sigma += (work.s[i] + alpha*ds_aff[i])*
      (work.z[i] + alpha*dz_aff[i]);
  sigma /= mu;
  sigma = sigma*sigma*sigma;
  /* Finish calculating mu now. */
  mu *= 0.038461538461538464;
  smu = sigma*mu;
  /* Fill-in the rhs. */
  for (i = 0; i < 23; i++)
    work.rhs[i] = 0;
  for (i = 49; i < 75; i++)
    work.rhs[i] = 0;
  for (i = 0; i < 26; i++)
    r2[i] = work.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}
void GensolverFirst::refine(double *target, double *var) {
  int i, j;
  double *residual = work.buffer;
  double norm2;
  double *new_var = work.buffer2;
  for (j = 0; j < settings.refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 75; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
#ifndef ZERO_LIBRARY_MODE
    if (settings.verbose_refinement) {
      if (j == 0)
        printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
      else
        printf("After refinement we get squared norm %.6g.\n", norm2);
    }
#endif
    /* Solve to find new_var = KKT \ (target - A*var). */
    ldl_solve(residual, new_var);
    /* Update var += new_var, or var += KKT \ (target - A*var). */
    for (i = 0; i < 75; i++) {
      var[i] -= new_var[i];
    }
  }
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 75; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }
    if (j == 0)
      printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
    else
      printf("After refinement we get squared norm %.6g.\n", norm2);
  }
#endif
}
double GensolverFirst::calc_ineq_resid_squared(void) {
  /* Calculates the norm ||-Gx - s + h||. */
  double norm2_squared;
  int i;
  /* Find -Gx. */
  multbymG(work.buffer, work.x);
  /* Add -s + h. */
  for (i = 0; i < 26; i++)
    work.buffer[i] += -work.s[i] + work.h[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 26; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
double GensolverFirst::calc_eq_resid_squared(void) {
  /* Calculates the norm ||-Ax + b||. */
  double norm2_squared;
  int i;
  /* Find -Ax. */
  multbymA(work.buffer, work.x);
  /* Add +b. */
  for (i = 0; i < 0; i++)
    work.buffer[i] += work.b[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 0; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
void GensolverFirst::better_start(void) {
  /* Calculates a better starting point, using a similar approach to CVXOPT. */
  /* Not yet speed optimized. */
  int i;
  double *x, *s, *z, *y;
  double alpha;
  work.block_33[0] = -1;
  /* Make sure sinvz is 1 to make hijacked KKT system ok. */
  for (i = 0; i < 26; i++)
    work.s_inv_z[i] = 1;
  fill_KKT();
  ldl_factor();
  fillrhs_start();
  /* Borrow work.lhs_aff for the solution. */
  ldl_solve(work.rhs, work.lhs_aff);
  /* Don't do any refinement for now. Precision doesn't matter too much. */
  x = work.lhs_aff;
  s = work.lhs_aff + 23;
  z = work.lhs_aff + 49;
  y = work.lhs_aff + 75;
  /* Just set x and y as is. */
  for (i = 0; i < 23; i++)
    work.x[i] = x[i];
  for (i = 0; i < 0; i++)
    work.y[i] = y[i];
  /* Now complete the initialization. Start with s. */
  /* Must have alpha > max(z). */
  alpha = -1e99;
  for (i = 0; i < 26; i++)
    if (alpha < z[i])
      alpha = z[i];
  if (alpha < 0) {
    for (i = 0; i < 26; i++)
      work.s[i] = -z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 26; i++)
      work.s[i] = -z[i] + alpha;
  }
  /* Now initialize z. */
  /* Now must have alpha > max(-z). */
  alpha = -1e99;
  for (i = 0; i < 26; i++)
    if (alpha < -z[i])
      alpha = -z[i];
  if (alpha < 0) {
    for (i = 0; i < 26; i++)
      work.z[i] = z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 26; i++)
      work.z[i] = z[i] + alpha;
  }
}
void GensolverFirst::fillrhs_start(void) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 23;
  r3 = work.rhs + 49;
  r4 = work.rhs + 75;
  for (i = 0; i < 23; i++)
    r1[i] = -work.q[i];
  for (i = 0; i < 26; i++)
    r2[i] = 0;
  for (i = 0; i < 26; i++)
    r3[i] = work.h[i];
  for (i = 0; i < 0; i++)
    r4[i] = work.b[i];
}
long GensolverFirst::solve(void) {
  int i;
  int iter;
  double *dx, *ds, *dy, *dz;
  double minval;
  double alpha;
  work.converged = 0;
  setup_pointers();
  pre_ops();
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose)
    printf("iter     objv        gap       |Ax-b|    |Gx+s-h|    step\n");
#endif
  fillq();
  fillh();
  fillb();
  if (settings.better_start)
    better_start();
  else
    set_start();
  for (iter = 0; iter < settings.max_iters; iter++) {
    for (i = 0; i < 26; i++) {
      work.s_inv[i] = 1.0 / work.s[i];
      work.s_inv_z[i] = work.s_inv[i]*work.z[i];
    }
    work.block_33[0] = 0;
    fill_KKT();
    ldl_factor();
    /* Affine scaling directions. */
    fillrhs_aff();
    ldl_solve(work.rhs, work.lhs_aff);
    refine(work.rhs, work.lhs_aff);
    /* Centering plus corrector directions. */
    fillrhs_cc();
    ldl_solve(work.rhs, work.lhs_cc);
    refine(work.rhs, work.lhs_cc);
    /* Add the two together and store in aff. */
    for (i = 0; i < 75; i++)
      work.lhs_aff[i] += work.lhs_cc[i];
    /* Rename aff to reflect its new meaning. */
    dx = work.lhs_aff;
    ds = work.lhs_aff + 23;
    dz = work.lhs_aff + 49;
    dy = work.lhs_aff + 75;
    /* Find min(min(ds./s), min(dz./z)). */
    minval = 0;
    for (i = 0; i < 26; i++)
      if (ds[i] < minval*work.s[i])
        minval = ds[i]/work.s[i];
    for (i = 0; i < 26; i++)
      if (dz[i] < minval*work.z[i])
        minval = dz[i]/work.z[i];
    /* Find alpha. */
    if (-0.99 < minval)
      alpha = 1;
    else
      alpha = -0.99/minval;
    /* Update the primal and dual variables. */
    for (i = 0; i < 23; i++)
      work.x[i] += alpha*dx[i];
    for (i = 0; i < 26; i++)
      work.s[i] += alpha*ds[i];
    for (i = 0; i < 26; i++)
      work.z[i] += alpha*dz[i];
    for (i = 0; i < 0; i++)
      work.y[i] += alpha*dy[i];
    work.gap = eval_gap();
    work.eq_resid_squared = calc_eq_resid_squared();
    work.ineq_resid_squared = calc_ineq_resid_squared();
#ifndef ZERO_LIBRARY_MODE
    if (settings.verbose) {
      work.optval = eval_objv();
      printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter+1, work.optval, work.gap, sqrt(work.eq_resid_squared),
          sqrt(work.ineq_resid_squared), alpha);
    }
#endif
    /* Test termination conditions. Requires optimality, and satisfied */
    /* constraints. */
    if (   (work.gap < settings.eps)
        && (work.eq_resid_squared <= settings.resid_tol*settings.resid_tol)
        && (work.ineq_resid_squared <= settings.resid_tol*settings.resid_tol)
       ) {
      work.converged = 1;
      work.optval = eval_objv();
      return iter+1;
    }
  }
  return iter;
}

//util.c
/* Produced by CVXGEN, 2016-05-03 19:35:15 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: util.c. */
/* Description: Common utility file for all cvxgen code. */

void GensolverFirst::tic(void) {
  tic_timestart = clock();
}
float GensolverFirst::toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.2f.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
float GensolverFirst::tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
void GensolverFirst::printmatrix(char *name, double *A, int m, int n, int sparse) {
  int i, j;
  printf("%s = [...\n", name);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      if ((sparse == 1) && (A[i+j*m] == 0))
        printf("         0");
      else
        printf("  % 9.4f", A[i+j*m]);
    printf(",\n");
  }
  printf("];\n");
}
double GensolverFirst::unif(double lower, double upper) {
  return lower + ((upper - lower)*rand())/RAND_MAX;
}
/* Next function is from numerical recipes in C. */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float GensolverFirst::ran1(long*idum, int reset) {
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (reset) {
    iy = 0;
  }
  if (*idum<=0||!iy) {
    if (-(*idum)<1)*idum=1;
    else *idum=-(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum<0)*idum+=IM;
      if (j<NTAB)iv[j]=*idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum<0)*idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy)> RNMX) return RNMX;
  else return temp;
}
/* Next function is from numerical recipes in C. */
float GensolverFirst::randn_internal(long *idum, int reset) {
  static int iset=0;
  static float gset;
  float fac, rsq, v1, v2;
  if (reset) {
    iset = 0;
  }
  if (iset==0) {
    do {
      v1 = 2.0*ran1(idum, reset)-1.0;
      v2 = 2.0*ran1(idum, reset)-1.0;
      rsq = v1*v1+v2*v2;
    } while(rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
double GensolverFirst::randn(void) {
  return randn_internal(&global_seed, 0);
}
void GensolverFirst::reset_rand(void) {
  srand(15);
  global_seed = 1;
  randn_internal(&global_seed, 1);
}

//matrix_support.c
/* Produced by CVXGEN, 2016-05-03 19:35:14 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */

void GensolverFirst::multbymA(double *lhs, double *rhs) {
}
void GensolverFirst::multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
  lhs[10] = 0;
  lhs[11] = 0;
  lhs[12] = 0;
  lhs[13] = 0;
  lhs[14] = 0;
  lhs[15] = 0;
  lhs[16] = 0;
  lhs[17] = 0;
  lhs[18] = 0;
  lhs[19] = 0;
  lhs[20] = 0;
  lhs[21] = 0;
  lhs[22] = 0;
}
void GensolverFirst::multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[0]*(1);
  lhs[2] = -rhs[0]*(params.ones[0])-rhs[1]*(-1);
  lhs[3] = -rhs[0]*(params.ones[1])-rhs[2]*(-1);
  lhs[4] = -rhs[0]*(params.ones[2])-rhs[3]*(-1);
  lhs[5] = -rhs[0]*(params.ones[3])-rhs[4]*(-1);
  lhs[6] = -rhs[0]*(params.ones[4])-rhs[5]*(-1);
  lhs[7] = -rhs[0]*(params.ones[5])-rhs[6]*(-1);
  lhs[8] = -rhs[0]*(params.ones[6])-rhs[7]*(-1);
  lhs[9] = -rhs[0]*(params.ones[7])-rhs[8]*(-1);
  lhs[10] = -rhs[0]*(params.ones[8])-rhs[9]*(-1);
  lhs[11] = -rhs[0]*(params.ones[9])-rhs[10]*(-1);
  lhs[12] = -rhs[0]*(params.ones[10])-rhs[11]*(-1);
  lhs[13] = -rhs[0]*(-params.ones[0])-rhs[1]*(1);
  lhs[14] = -rhs[0]*(-params.ones[1])-rhs[2]*(1);
  lhs[15] = -rhs[0]*(-params.ones[2])-rhs[3]*(1);
  lhs[16] = -rhs[0]*(-params.ones[3])-rhs[4]*(1);
  lhs[17] = -rhs[0]*(-params.ones[4])-rhs[5]*(1);
  lhs[18] = -rhs[0]*(-params.ones[5])-rhs[6]*(1);
  lhs[19] = -rhs[0]*(-params.ones[6])-rhs[7]*(1);
  lhs[20] = -rhs[0]*(-params.ones[7])-rhs[8]*(1);
  lhs[21] = -rhs[0]*(-params.ones[8])-rhs[9]*(1);
  lhs[22] = -rhs[0]*(-params.ones[9])-rhs[10]*(1);
  lhs[23] = -rhs[0]*(-params.ones[10])-rhs[11]*(1);
  lhs[24] = -rhs[0]*(-1);
  lhs[25] = -rhs[0]*(1);
}
void GensolverFirst::multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(params.ones[0])-rhs[3]*(params.ones[1])-rhs[4]*(params.ones[2])-rhs[5]*(params.ones[3])-rhs[6]*(params.ones[4])-rhs[7]*(params.ones[5])-rhs[8]*(params.ones[6])-rhs[9]*(params.ones[7])-rhs[10]*(params.ones[8])-rhs[11]*(params.ones[9])-rhs[12]*(params.ones[10])-rhs[13]*(-params.ones[0])-rhs[14]*(-params.ones[1])-rhs[15]*(-params.ones[2])-rhs[16]*(-params.ones[3])-rhs[17]*(-params.ones[4])-rhs[18]*(-params.ones[5])-rhs[19]*(-params.ones[6])-rhs[20]*(-params.ones[7])-rhs[21]*(-params.ones[8])-rhs[22]*(-params.ones[9])-rhs[23]*(-params.ones[10])-rhs[24]*(-1)-rhs[25]*(1);
  lhs[1] = -rhs[2]*(-1)-rhs[13]*(1);
  lhs[2] = -rhs[3]*(-1)-rhs[14]*(1);
  lhs[3] = -rhs[4]*(-1)-rhs[15]*(1);
  lhs[4] = -rhs[5]*(-1)-rhs[16]*(1);
  lhs[5] = -rhs[6]*(-1)-rhs[17]*(1);
  lhs[6] = -rhs[7]*(-1)-rhs[18]*(1);
  lhs[7] = -rhs[8]*(-1)-rhs[19]*(1);
  lhs[8] = -rhs[9]*(-1)-rhs[20]*(1);
  lhs[9] = -rhs[10]*(-1)-rhs[21]*(1);
  lhs[10] = -rhs[11]*(-1)-rhs[22]*(1);
  lhs[11] = -rhs[12]*(-1)-rhs[23]*(1);
  lhs[12] = 0;
  lhs[13] = 0;
  lhs[14] = 0;
  lhs[15] = 0;
  lhs[16] = 0;
  lhs[17] = 0;
  lhs[18] = 0;
  lhs[19] = 0;
  lhs[20] = 0;
  lhs[21] = 0;
  lhs[22] = 0;
}
void GensolverFirst::multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848*work.quad_877793640448[0]));
  lhs[1] = rhs[1]*(2*work.frac_479733190656);
  lhs[2] = rhs[2]*(2*work.frac_479733190656);
  lhs[3] = rhs[3]*(2*work.frac_479733190656);
  lhs[4] = rhs[4]*(2*work.frac_479733190656);
  lhs[5] = rhs[5]*(2*work.frac_479733190656);
  lhs[6] = rhs[6]*(2*work.frac_479733190656);
  lhs[7] = rhs[7]*(2*work.frac_479733190656);
  lhs[8] = rhs[8]*(2*work.frac_479733190656);
  lhs[9] = rhs[9]*(2*work.frac_479733190656);
  lhs[10] = rhs[10]*(2*work.frac_479733190656);
  lhs[11] = rhs[11]*(2*work.frac_479733190656);
  lhs[12] = rhs[12]*(2*work.frac_121674190848);
  lhs[13] = rhs[13]*(2*work.frac_121674190848);
  lhs[14] = rhs[14]*(2*work.frac_121674190848);
  lhs[15] = rhs[15]*(2*work.frac_121674190848);
  lhs[16] = rhs[16]*(2*work.frac_121674190848);
  lhs[17] = rhs[17]*(2*work.frac_121674190848);
  lhs[18] = rhs[18]*(2*work.frac_121674190848);
  lhs[19] = rhs[19]*(2*work.frac_121674190848);
  lhs[20] = rhs[20]*(2*work.frac_121674190848);
  lhs[21] = rhs[21]*(2*work.frac_121674190848);
  lhs[22] = rhs[22]*(2*work.frac_121674190848);
}
void GensolverFirst::fillq(void) {
  work.q[0] = params.c1[0]-2*params.PgNu[0]*work.frac_479733190656+params.gamma[0]*(params.B[0]+params.B[1]+params.B[2]+params.B[3]+params.B[4]+params.B[5]+params.B[6]+params.B[7]+params.B[8]+params.B[9]+params.B[10])+params.lambda_1[0]+params.lambda_1[1]+params.lambda_1[2]+params.lambda_1[3]+params.lambda_1[4]+params.lambda_1[5]+params.lambda_1[6]+params.lambda_1[7]+params.lambda_1[8]+params.lambda_1[9]+params.lambda_1[10]+2*work.frac_121674190848*((-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*params.ones[0]+(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])*params.ones[1]+(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])*params.ones[2]+(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])*params.ones[3]+(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])*params.ones[4]+(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])*params.ones[5]+(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])*params.ones[6]+(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])*params.ones[7]+(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])*params.ones[8]+(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])*params.ones[9]+(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10])*params.ones[10]);
  work.q[1] = -2*work.frac_479733190656*params.PgNextNu[0]+params.gamma[0]*params.D[0]+params.lambda_2[0];
  work.q[2] = -2*work.frac_479733190656*params.PgNextNu[1]+params.gamma[0]*params.D[1]+params.lambda_2[1];
  work.q[3] = -2*work.frac_479733190656*params.PgNextNu[2]+params.gamma[0]*params.D[2]+params.lambda_2[2];
  work.q[4] = -2*work.frac_479733190656*params.PgNextNu[3]+params.gamma[0]*params.D[3]+params.lambda_2[3];
  work.q[5] = -2*work.frac_479733190656*params.PgNextNu[4]+params.gamma[0]*params.D[4]+params.lambda_2[4];
  work.q[6] = -2*work.frac_479733190656*params.PgNextNu[5]+params.gamma[0]*params.D[5]+params.lambda_2[5];
  work.q[7] = -2*work.frac_479733190656*params.PgNextNu[6]+params.gamma[0]*params.D[6]+params.lambda_2[6];
  work.q[8] = -2*work.frac_479733190656*params.PgNextNu[7]+params.gamma[0]*params.D[7]+params.lambda_2[7];
  work.q[9] = -2*work.frac_479733190656*params.PgNextNu[8]+params.gamma[0]*params.D[8]+params.lambda_2[8];
  work.q[10] = -2*work.frac_479733190656*params.PgNextNu[9]+params.gamma[0]*params.D[9]+params.lambda_2[9];
  work.q[11] = -2*work.frac_479733190656*params.PgNextNu[10]+params.gamma[0]*params.D[10]+params.lambda_2[10];
  work.q[12] = 2*work.frac_121674190848*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0]);
  work.q[13] = 2*work.frac_121674190848*(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1]);
  work.q[14] = 2*work.frac_121674190848*(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2]);
  work.q[15] = 2*work.frac_121674190848*(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3]);
  work.q[16] = 2*work.frac_121674190848*(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4]);
  work.q[17] = 2*work.frac_121674190848*(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5]);
  work.q[18] = 2*work.frac_121674190848*(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6]);
  work.q[19] = 2*work.frac_121674190848*(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7]);
  work.q[20] = 2*work.frac_121674190848*(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8]);
  work.q[21] = 2*work.frac_121674190848*(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9]);
  work.q[22] = 2*work.frac_121674190848*(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10]);
}
void GensolverFirst::fillh(void) {
  work.h[0] = -params.PgMin[0];
  work.h[1] = params.PgMax[0];
  work.h[2] = -params.RgMin[0]*params.ones[0];
  work.h[3] = -params.RgMin[0]*params.ones[1];
  work.h[4] = -params.RgMin[0]*params.ones[2];
  work.h[5] = -params.RgMin[0]*params.ones[3];
  work.h[6] = -params.RgMin[0]*params.ones[4];
  work.h[7] = -params.RgMin[0]*params.ones[5];
  work.h[8] = -params.RgMin[0]*params.ones[6];
  work.h[9] = -params.RgMin[0]*params.ones[7];
  work.h[10] = -params.RgMin[0]*params.ones[8];
  work.h[11] = -params.RgMin[0]*params.ones[9];
  work.h[12] = -params.RgMin[0]*params.ones[10];
  work.h[13] = params.RgMax[0]*params.ones[0];
  work.h[14] = params.RgMax[0]*params.ones[1];
  work.h[15] = params.RgMax[0]*params.ones[2];
  work.h[16] = params.RgMax[0]*params.ones[3];
  work.h[17] = params.RgMax[0]*params.ones[4];
  work.h[18] = params.RgMax[0]*params.ones[5];
  work.h[19] = params.RgMax[0]*params.ones[6];
  work.h[20] = params.RgMax[0]*params.ones[7];
  work.h[21] = params.RgMax[0]*params.ones[8];
  work.h[22] = params.RgMax[0]*params.ones[9];
  work.h[23] = params.RgMax[0]*params.ones[10];
  work.h[24] = -(params.RgMin[0]+params.PgPrev[0]);
  work.h[25] = -(-params.PgPrev[0]-params.RgMax[0]);
}
void GensolverFirst::fillb(void) {
}
void GensolverFirst::pre_ops(void) {
  work.frac_479733190656 = params.beta[0];
  work.frac_479733190656 /= 2;
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_877793640448[0] = params.ones[0]*params.ones[0]+params.ones[1]*params.ones[1]+params.ones[2]*params.ones[2]+params.ones[3]*params.ones[3]+params.ones[4]*params.ones[4]+params.ones[5]*params.ones[5]+params.ones[6]*params.ones[6]+params.ones[7]*params.ones[7]+params.ones[8]*params.ones[8]+params.ones[9]*params.ones[9]+params.ones[10]*params.ones[10];
  work.quad_608032620544[0] = params.PgNu[0]*work.frac_479733190656*params.PgNu[0];
  work.quad_900708769792[0] = work.frac_479733190656*(params.PgNextNu[0]*params.PgNextNu[0]+params.PgNextNu[1]*params.PgNextNu[1]+params.PgNextNu[2]*params.PgNextNu[2]+params.PgNextNu[3]*params.PgNextNu[3]+params.PgNextNu[4]*params.PgNextNu[4]+params.PgNextNu[5]*params.PgNextNu[5]+params.PgNextNu[6]*params.PgNextNu[6]+params.PgNextNu[7]*params.PgNextNu[7]+params.PgNextNu[8]*params.PgNextNu[8]+params.PgNextNu[9]*params.PgNextNu[9]+params.PgNextNu[10]*params.PgNextNu[10]);
  work.quad_324888432640[0] = ((-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])+(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])*(-params.Pg_N_init[1]+params.Pg_N_avg[1]+params.ug_N[1])+(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])*(-params.Pg_N_init[2]+params.Pg_N_avg[2]+params.ug_N[2])+(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])*(-params.Pg_N_init[3]+params.Pg_N_avg[3]+params.ug_N[3])+(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])*(-params.Pg_N_init[4]+params.Pg_N_avg[4]+params.ug_N[4])+(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])*(-params.Pg_N_init[5]+params.Pg_N_avg[5]+params.ug_N[5])+(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])*(-params.Pg_N_init[6]+params.Pg_N_avg[6]+params.ug_N[6])+(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])*(-params.Pg_N_init[7]+params.Pg_N_avg[7]+params.ug_N[7])+(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])*(-params.Pg_N_init[8]+params.Pg_N_avg[8]+params.ug_N[8])+(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])*(-params.Pg_N_init[9]+params.Pg_N_avg[9]+params.ug_N[9])+(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10])*(-params.Pg_N_init[10]+params.Pg_N_avg[10]+params.ug_N[10]));
  work.quad_331507691520[0] = ((-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])+(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1])*(-params.Vg_N_avg[1]-params.Thetag_N_avg[1]+params.vg_N[1])+(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2])*(-params.Vg_N_avg[2]-params.Thetag_N_avg[2]+params.vg_N[2])+(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3])*(-params.Vg_N_avg[3]-params.Thetag_N_avg[3]+params.vg_N[3])+(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4])*(-params.Vg_N_avg[4]-params.Thetag_N_avg[4]+params.vg_N[4])+(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5])*(-params.Vg_N_avg[5]-params.Thetag_N_avg[5]+params.vg_N[5])+(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6])*(-params.Vg_N_avg[6]-params.Thetag_N_avg[6]+params.vg_N[6])+(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7])*(-params.Vg_N_avg[7]-params.Thetag_N_avg[7]+params.vg_N[7])+(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8])*(-params.Vg_N_avg[8]-params.Thetag_N_avg[8]+params.vg_N[8])+(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9])*(-params.Vg_N_avg[9]-params.Thetag_N_avg[9]+params.vg_N[9])+(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10])*(-params.Vg_N_avg[10]-params.Thetag_N_avg[10]+params.vg_N[10]));
}

//ldl.c
/* Produced by CVXGEN, 2016-05-03 19:35:14 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: ldl.c. */
/* Description: Basic test harness for solver.c. */

/* Be sure to place ldl_solve first, so storage schemes are defined by it. */
void GensolverFirst::ldl_solve(double *target, double *var) {
  int i;
  /* Find var = (L*diag(work.d)*L') \ target, then unpermute. */
  /* Answer goes into var. */
  /* Forward substitution. */
  /* Include permutation as we retrieve from target. Use v so we can unpermute */
  /* later. */
  work.v[0] = target[12];
  work.v[1] = target[13];
  work.v[2] = target[14];
  work.v[3] = target[15];
  work.v[4] = target[16];
  work.v[5] = target[17];
  work.v[6] = target[18];
  work.v[7] = target[19];
  work.v[8] = target[20];
  work.v[9] = target[21];
  work.v[10] = target[22];
  work.v[11] = target[23];
  work.v[12] = target[24];
  work.v[13] = target[25];
  work.v[14] = target[26];
  work.v[15] = target[27];
  work.v[16] = target[28];
  work.v[17] = target[29];
  work.v[18] = target[30];
  work.v[19] = target[31];
  work.v[20] = target[32];
  work.v[21] = target[33];
  work.v[22] = target[34];
  work.v[23] = target[35];
  work.v[24] = target[36];
  work.v[25] = target[37];
  work.v[26] = target[38];
  work.v[27] = target[39];
  work.v[28] = target[40];
  work.v[29] = target[41];
  work.v[30] = target[42];
  work.v[31] = target[43];
  work.v[32] = target[44];
  work.v[33] = target[45];
  work.v[34] = target[46];
  work.v[35] = target[47];
  work.v[36] = target[48];
  work.v[37] = target[49]-work.L[0]*work.v[11];
  work.v[38] = target[50]-work.L[1]*work.v[12];
  work.v[39] = target[73]-work.L[2]*work.v[35];
  work.v[40] = target[74]-work.L[3]*work.v[36];
  work.v[41] = target[1];
  work.v[42] = target[51]-work.L[4]*work.v[13]-work.L[5]*work.v[41];
  work.v[43] = target[62]-work.L[6]*work.v[24]-work.L[7]*work.v[41]-work.L[8]*work.v[42];
  work.v[44] = target[2];
  work.v[45] = target[52]-work.L[9]*work.v[14]-work.L[10]*work.v[44];
  work.v[46] = target[63]-work.L[11]*work.v[25]-work.L[12]*work.v[44]-work.L[13]*work.v[45];
  work.v[47] = target[3];
  work.v[48] = target[53]-work.L[14]*work.v[15]-work.L[15]*work.v[47];
  work.v[49] = target[64]-work.L[16]*work.v[26]-work.L[17]*work.v[47]-work.L[18]*work.v[48];
  work.v[50] = target[4];
  work.v[51] = target[54]-work.L[19]*work.v[16]-work.L[20]*work.v[50];
  work.v[52] = target[65]-work.L[21]*work.v[27]-work.L[22]*work.v[50]-work.L[23]*work.v[51];
  work.v[53] = target[5];
  work.v[54] = target[55]-work.L[24]*work.v[17]-work.L[25]*work.v[53];
  work.v[55] = target[66]-work.L[26]*work.v[28]-work.L[27]*work.v[53]-work.L[28]*work.v[54];
  work.v[56] = target[6];
  work.v[57] = target[56]-work.L[29]*work.v[18]-work.L[30]*work.v[56];
  work.v[58] = target[67]-work.L[31]*work.v[29]-work.L[32]*work.v[56]-work.L[33]*work.v[57];
  work.v[59] = target[7];
  work.v[60] = target[57]-work.L[34]*work.v[19]-work.L[35]*work.v[59];
  work.v[61] = target[68]-work.L[36]*work.v[30]-work.L[37]*work.v[59]-work.L[38]*work.v[60];
  work.v[62] = target[8];
  work.v[63] = target[58]-work.L[39]*work.v[20]-work.L[40]*work.v[62];
  work.v[64] = target[69]-work.L[41]*work.v[31]-work.L[42]*work.v[62]-work.L[43]*work.v[63];
  work.v[65] = target[9];
  work.v[66] = target[59]-work.L[44]*work.v[21]-work.L[45]*work.v[65];
  work.v[67] = target[70]-work.L[46]*work.v[32]-work.L[47]*work.v[65]-work.L[48]*work.v[66];
  work.v[68] = target[10];
  work.v[69] = target[60]-work.L[49]*work.v[22]-work.L[50]*work.v[68];
  work.v[70] = target[71]-work.L[51]*work.v[33]-work.L[52]*work.v[68]-work.L[53]*work.v[69];
  work.v[71] = target[0]-work.L[54]*work.v[37]-work.L[55]*work.v[38]-work.L[56]*work.v[39]-work.L[57]*work.v[40]-work.L[58]*work.v[42]-work.L[59]*work.v[43]-work.L[60]*work.v[45]-work.L[61]*work.v[46]-work.L[62]*work.v[48]-work.L[63]*work.v[49]-work.L[64]*work.v[51]-work.L[65]*work.v[52]-work.L[66]*work.v[54]-work.L[67]*work.v[55]-work.L[68]*work.v[57]-work.L[69]*work.v[58]-work.L[70]*work.v[60]-work.L[71]*work.v[61]-work.L[72]*work.v[63]-work.L[73]*work.v[64]-work.L[74]*work.v[66]-work.L[75]*work.v[67]-work.L[76]*work.v[69]-work.L[77]*work.v[70];
  work.v[72] = target[11];
  work.v[73] = target[61]-work.L[78]*work.v[23]-work.L[79]*work.v[71]-work.L[80]*work.v[72];
  work.v[74] = target[72]-work.L[81]*work.v[34]-work.L[82]*work.v[71]-work.L[83]*work.v[72]-work.L[84]*work.v[73];
  /* Diagonal scaling. Assume correctness of work.d_inv. */
  for (i = 0; i < 75; i++)
    work.v[i] *= work.d_inv[i];
  /* Back substitution */
  work.v[73] -= work.L[84]*work.v[74];
  work.v[72] -= work.L[80]*work.v[73]+work.L[83]*work.v[74];
  work.v[71] -= work.L[79]*work.v[73]+work.L[82]*work.v[74];
  work.v[70] -= work.L[77]*work.v[71];
  work.v[69] -= work.L[53]*work.v[70]+work.L[76]*work.v[71];
  work.v[68] -= work.L[50]*work.v[69]+work.L[52]*work.v[70];
  work.v[67] -= work.L[75]*work.v[71];
  work.v[66] -= work.L[48]*work.v[67]+work.L[74]*work.v[71];
  work.v[65] -= work.L[45]*work.v[66]+work.L[47]*work.v[67];
  work.v[64] -= work.L[73]*work.v[71];
  work.v[63] -= work.L[43]*work.v[64]+work.L[72]*work.v[71];
  work.v[62] -= work.L[40]*work.v[63]+work.L[42]*work.v[64];
  work.v[61] -= work.L[71]*work.v[71];
  work.v[60] -= work.L[38]*work.v[61]+work.L[70]*work.v[71];
  work.v[59] -= work.L[35]*work.v[60]+work.L[37]*work.v[61];
  work.v[58] -= work.L[69]*work.v[71];
  work.v[57] -= work.L[33]*work.v[58]+work.L[68]*work.v[71];
  work.v[56] -= work.L[30]*work.v[57]+work.L[32]*work.v[58];
  work.v[55] -= work.L[67]*work.v[71];
  work.v[54] -= work.L[28]*work.v[55]+work.L[66]*work.v[71];
  work.v[53] -= work.L[25]*work.v[54]+work.L[27]*work.v[55];
  work.v[52] -= work.L[65]*work.v[71];
  work.v[51] -= work.L[23]*work.v[52]+work.L[64]*work.v[71];
  work.v[50] -= work.L[20]*work.v[51]+work.L[22]*work.v[52];
  work.v[49] -= work.L[63]*work.v[71];
  work.v[48] -= work.L[18]*work.v[49]+work.L[62]*work.v[71];
  work.v[47] -= work.L[15]*work.v[48]+work.L[17]*work.v[49];
  work.v[46] -= work.L[61]*work.v[71];
  work.v[45] -= work.L[13]*work.v[46]+work.L[60]*work.v[71];
  work.v[44] -= work.L[10]*work.v[45]+work.L[12]*work.v[46];
  work.v[43] -= work.L[59]*work.v[71];
  work.v[42] -= work.L[8]*work.v[43]+work.L[58]*work.v[71];
  work.v[41] -= work.L[5]*work.v[42]+work.L[7]*work.v[43];
  work.v[40] -= work.L[57]*work.v[71];
  work.v[39] -= work.L[56]*work.v[71];
  work.v[38] -= work.L[55]*work.v[71];
  work.v[37] -= work.L[54]*work.v[71];
  work.v[36] -= work.L[3]*work.v[40];
  work.v[35] -= work.L[2]*work.v[39];
  work.v[34] -= work.L[81]*work.v[74];
  work.v[33] -= work.L[51]*work.v[70];
  work.v[32] -= work.L[46]*work.v[67];
  work.v[31] -= work.L[41]*work.v[64];
  work.v[30] -= work.L[36]*work.v[61];
  work.v[29] -= work.L[31]*work.v[58];
  work.v[28] -= work.L[26]*work.v[55];
  work.v[27] -= work.L[21]*work.v[52];
  work.v[26] -= work.L[16]*work.v[49];
  work.v[25] -= work.L[11]*work.v[46];
  work.v[24] -= work.L[6]*work.v[43];
  work.v[23] -= work.L[78]*work.v[73];
  work.v[22] -= work.L[49]*work.v[69];
  work.v[21] -= work.L[44]*work.v[66];
  work.v[20] -= work.L[39]*work.v[63];
  work.v[19] -= work.L[34]*work.v[60];
  work.v[18] -= work.L[29]*work.v[57];
  work.v[17] -= work.L[24]*work.v[54];
  work.v[16] -= work.L[19]*work.v[51];
  work.v[15] -= work.L[14]*work.v[48];
  work.v[14] -= work.L[9]*work.v[45];
  work.v[13] -= work.L[4]*work.v[42];
  work.v[12] -= work.L[1]*work.v[38];
  work.v[11] -= work.L[0]*work.v[37];
  /* Unpermute the result, from v to var. */
  var[0] = work.v[71];
  var[1] = work.v[41];
  var[2] = work.v[44];
  var[3] = work.v[47];
  var[4] = work.v[50];
  var[5] = work.v[53];
  var[6] = work.v[56];
  var[7] = work.v[59];
  var[8] = work.v[62];
  var[9] = work.v[65];
  var[10] = work.v[68];
  var[11] = work.v[72];
  var[12] = work.v[0];
  var[13] = work.v[1];
  var[14] = work.v[2];
  var[15] = work.v[3];
  var[16] = work.v[4];
  var[17] = work.v[5];
  var[18] = work.v[6];
  var[19] = work.v[7];
  var[20] = work.v[8];
  var[21] = work.v[9];
  var[22] = work.v[10];
  var[23] = work.v[11];
  var[24] = work.v[12];
  var[25] = work.v[13];
  var[26] = work.v[14];
  var[27] = work.v[15];
  var[28] = work.v[16];
  var[29] = work.v[17];
  var[30] = work.v[18];
  var[31] = work.v[19];
  var[32] = work.v[20];
  var[33] = work.v[21];
  var[34] = work.v[22];
  var[35] = work.v[23];
  var[36] = work.v[24];
  var[37] = work.v[25];
  var[38] = work.v[26];
  var[39] = work.v[27];
  var[40] = work.v[28];
  var[41] = work.v[29];
  var[42] = work.v[30];
  var[43] = work.v[31];
  var[44] = work.v[32];
  var[45] = work.v[33];
  var[46] = work.v[34];
  var[47] = work.v[35];
  var[48] = work.v[36];
  var[49] = work.v[37];
  var[50] = work.v[38];
  var[51] = work.v[42];
  var[52] = work.v[45];
  var[53] = work.v[48];
  var[54] = work.v[51];
  var[55] = work.v[54];
  var[56] = work.v[57];
  var[57] = work.v[60];
  var[58] = work.v[63];
  var[59] = work.v[66];
  var[60] = work.v[69];
  var[61] = work.v[73];
  var[62] = work.v[43];
  var[63] = work.v[46];
  var[64] = work.v[49];
  var[65] = work.v[52];
  var[66] = work.v[55];
  var[67] = work.v[58];
  var[68] = work.v[61];
  var[69] = work.v[64];
  var[70] = work.v[67];
  var[71] = work.v[70];
  var[72] = work.v[74];
  var[73] = work.v[39];
  var[74] = work.v[40];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared norm for solution is %.8g.\n", check_residual(target, var));
  }
#endif
}
void GensolverFirst::ldl_factor(void) {
  work.d[0] = work.KKT[0];
  if (work.d[0] < 0)
    work.d[0] = settings.kkt_reg;
  else
    work.d[0] += settings.kkt_reg;
  work.d_inv[0] = 1/work.d[0];
  work.v[1] = work.KKT[1];
  work.d[1] = work.v[1];
  if (work.d[1] < 0)
    work.d[1] = settings.kkt_reg;
  else
    work.d[1] += settings.kkt_reg;
  work.d_inv[1] = 1/work.d[1];
  work.v[2] = work.KKT[2];
  work.d[2] = work.v[2];
  if (work.d[2] < 0)
    work.d[2] = settings.kkt_reg;
  else
    work.d[2] += settings.kkt_reg;
  work.d_inv[2] = 1/work.d[2];
  work.v[3] = work.KKT[3];
  work.d[3] = work.v[3];
  if (work.d[3] < 0)
    work.d[3] = settings.kkt_reg;
  else
    work.d[3] += settings.kkt_reg;
  work.d_inv[3] = 1/work.d[3];
  work.v[4] = work.KKT[4];
  work.d[4] = work.v[4];
  if (work.d[4] < 0)
    work.d[4] = settings.kkt_reg;
  else
    work.d[4] += settings.kkt_reg;
  work.d_inv[4] = 1/work.d[4];
  work.v[5] = work.KKT[5];
  work.d[5] = work.v[5];
  if (work.d[5] < 0)
    work.d[5] = settings.kkt_reg;
  else
    work.d[5] += settings.kkt_reg;
  work.d_inv[5] = 1/work.d[5];
  work.v[6] = work.KKT[6];
  work.d[6] = work.v[6];
  if (work.d[6] < 0)
    work.d[6] = settings.kkt_reg;
  else
    work.d[6] += settings.kkt_reg;
  work.d_inv[6] = 1/work.d[6];
  work.v[7] = work.KKT[7];
  work.d[7] = work.v[7];
  if (work.d[7] < 0)
    work.d[7] = settings.kkt_reg;
  else
    work.d[7] += settings.kkt_reg;
  work.d_inv[7] = 1/work.d[7];
  work.v[8] = work.KKT[8];
  work.d[8] = work.v[8];
  if (work.d[8] < 0)
    work.d[8] = settings.kkt_reg;
  else
    work.d[8] += settings.kkt_reg;
  work.d_inv[8] = 1/work.d[8];
  work.v[9] = work.KKT[9];
  work.d[9] = work.v[9];
  if (work.d[9] < 0)
    work.d[9] = settings.kkt_reg;
  else
    work.d[9] += settings.kkt_reg;
  work.d_inv[9] = 1/work.d[9];
  work.v[10] = work.KKT[10];
  work.d[10] = work.v[10];
  if (work.d[10] < 0)
    work.d[10] = settings.kkt_reg;
  else
    work.d[10] += settings.kkt_reg;
  work.d_inv[10] = 1/work.d[10];
  work.v[11] = work.KKT[11];
  work.d[11] = work.v[11];
  if (work.d[11] < 0)
    work.d[11] = settings.kkt_reg;
  else
    work.d[11] += settings.kkt_reg;
  work.d_inv[11] = 1/work.d[11];
  work.L[0] = (work.KKT[12])*work.d_inv[11];
  work.v[12] = work.KKT[13];
  work.d[12] = work.v[12];
  if (work.d[12] < 0)
    work.d[12] = settings.kkt_reg;
  else
    work.d[12] += settings.kkt_reg;
  work.d_inv[12] = 1/work.d[12];
  work.L[1] = (work.KKT[14])*work.d_inv[12];
  work.v[13] = work.KKT[15];
  work.d[13] = work.v[13];
  if (work.d[13] < 0)
    work.d[13] = settings.kkt_reg;
  else
    work.d[13] += settings.kkt_reg;
  work.d_inv[13] = 1/work.d[13];
  work.L[4] = (work.KKT[16])*work.d_inv[13];
  work.v[14] = work.KKT[17];
  work.d[14] = work.v[14];
  if (work.d[14] < 0)
    work.d[14] = settings.kkt_reg;
  else
    work.d[14] += settings.kkt_reg;
  work.d_inv[14] = 1/work.d[14];
  work.L[9] = (work.KKT[18])*work.d_inv[14];
  work.v[15] = work.KKT[19];
  work.d[15] = work.v[15];
  if (work.d[15] < 0)
    work.d[15] = settings.kkt_reg;
  else
    work.d[15] += settings.kkt_reg;
  work.d_inv[15] = 1/work.d[15];
  work.L[14] = (work.KKT[20])*work.d_inv[15];
  work.v[16] = work.KKT[21];
  work.d[16] = work.v[16];
  if (work.d[16] < 0)
    work.d[16] = settings.kkt_reg;
  else
    work.d[16] += settings.kkt_reg;
  work.d_inv[16] = 1/work.d[16];
  work.L[19] = (work.KKT[22])*work.d_inv[16];
  work.v[17] = work.KKT[23];
  work.d[17] = work.v[17];
  if (work.d[17] < 0)
    work.d[17] = settings.kkt_reg;
  else
    work.d[17] += settings.kkt_reg;
  work.d_inv[17] = 1/work.d[17];
  work.L[24] = (work.KKT[24])*work.d_inv[17];
  work.v[18] = work.KKT[25];
  work.d[18] = work.v[18];
  if (work.d[18] < 0)
    work.d[18] = settings.kkt_reg;
  else
    work.d[18] += settings.kkt_reg;
  work.d_inv[18] = 1/work.d[18];
  work.L[29] = (work.KKT[26])*work.d_inv[18];
  work.v[19] = work.KKT[27];
  work.d[19] = work.v[19];
  if (work.d[19] < 0)
    work.d[19] = settings.kkt_reg;
  else
    work.d[19] += settings.kkt_reg;
  work.d_inv[19] = 1/work.d[19];
  work.L[34] = (work.KKT[28])*work.d_inv[19];
  work.v[20] = work.KKT[29];
  work.d[20] = work.v[20];
  if (work.d[20] < 0)
    work.d[20] = settings.kkt_reg;
  else
    work.d[20] += settings.kkt_reg;
  work.d_inv[20] = 1/work.d[20];
  work.L[39] = (work.KKT[30])*work.d_inv[20];
  work.v[21] = work.KKT[31];
  work.d[21] = work.v[21];
  if (work.d[21] < 0)
    work.d[21] = settings.kkt_reg;
  else
    work.d[21] += settings.kkt_reg;
  work.d_inv[21] = 1/work.d[21];
  work.L[44] = (work.KKT[32])*work.d_inv[21];
  work.v[22] = work.KKT[33];
  work.d[22] = work.v[22];
  if (work.d[22] < 0)
    work.d[22] = settings.kkt_reg;
  else
    work.d[22] += settings.kkt_reg;
  work.d_inv[22] = 1/work.d[22];
  work.L[49] = (work.KKT[34])*work.d_inv[22];
  work.v[23] = work.KKT[35];
  work.d[23] = work.v[23];
  if (work.d[23] < 0)
    work.d[23] = settings.kkt_reg;
  else
    work.d[23] += settings.kkt_reg;
  work.d_inv[23] = 1/work.d[23];
  work.L[78] = (work.KKT[36])*work.d_inv[23];
  work.v[24] = work.KKT[37];
  work.d[24] = work.v[24];
  if (work.d[24] < 0)
    work.d[24] = settings.kkt_reg;
  else
    work.d[24] += settings.kkt_reg;
  work.d_inv[24] = 1/work.d[24];
  work.L[6] = (work.KKT[38])*work.d_inv[24];
  work.v[25] = work.KKT[39];
  work.d[25] = work.v[25];
  if (work.d[25] < 0)
    work.d[25] = settings.kkt_reg;
  else
    work.d[25] += settings.kkt_reg;
  work.d_inv[25] = 1/work.d[25];
  work.L[11] = (work.KKT[40])*work.d_inv[25];
  work.v[26] = work.KKT[41];
  work.d[26] = work.v[26];
  if (work.d[26] < 0)
    work.d[26] = settings.kkt_reg;
  else
    work.d[26] += settings.kkt_reg;
  work.d_inv[26] = 1/work.d[26];
  work.L[16] = (work.KKT[42])*work.d_inv[26];
  work.v[27] = work.KKT[43];
  work.d[27] = work.v[27];
  if (work.d[27] < 0)
    work.d[27] = settings.kkt_reg;
  else
    work.d[27] += settings.kkt_reg;
  work.d_inv[27] = 1/work.d[27];
  work.L[21] = (work.KKT[44])*work.d_inv[27];
  work.v[28] = work.KKT[45];
  work.d[28] = work.v[28];
  if (work.d[28] < 0)
    work.d[28] = settings.kkt_reg;
  else
    work.d[28] += settings.kkt_reg;
  work.d_inv[28] = 1/work.d[28];
  work.L[26] = (work.KKT[46])*work.d_inv[28];
  work.v[29] = work.KKT[47];
  work.d[29] = work.v[29];
  if (work.d[29] < 0)
    work.d[29] = settings.kkt_reg;
  else
    work.d[29] += settings.kkt_reg;
  work.d_inv[29] = 1/work.d[29];
  work.L[31] = (work.KKT[48])*work.d_inv[29];
  work.v[30] = work.KKT[49];
  work.d[30] = work.v[30];
  if (work.d[30] < 0)
    work.d[30] = settings.kkt_reg;
  else
    work.d[30] += settings.kkt_reg;
  work.d_inv[30] = 1/work.d[30];
  work.L[36] = (work.KKT[50])*work.d_inv[30];
  work.v[31] = work.KKT[51];
  work.d[31] = work.v[31];
  if (work.d[31] < 0)
    work.d[31] = settings.kkt_reg;
  else
    work.d[31] += settings.kkt_reg;
  work.d_inv[31] = 1/work.d[31];
  work.L[41] = (work.KKT[52])*work.d_inv[31];
  work.v[32] = work.KKT[53];
  work.d[32] = work.v[32];
  if (work.d[32] < 0)
    work.d[32] = settings.kkt_reg;
  else
    work.d[32] += settings.kkt_reg;
  work.d_inv[32] = 1/work.d[32];
  work.L[46] = (work.KKT[54])*work.d_inv[32];
  work.v[33] = work.KKT[55];
  work.d[33] = work.v[33];
  if (work.d[33] < 0)
    work.d[33] = settings.kkt_reg;
  else
    work.d[33] += settings.kkt_reg;
  work.d_inv[33] = 1/work.d[33];
  work.L[51] = (work.KKT[56])*work.d_inv[33];
  work.v[34] = work.KKT[57];
  work.d[34] = work.v[34];
  if (work.d[34] < 0)
    work.d[34] = settings.kkt_reg;
  else
    work.d[34] += settings.kkt_reg;
  work.d_inv[34] = 1/work.d[34];
  work.L[81] = (work.KKT[58])*work.d_inv[34];
  work.v[35] = work.KKT[59];
  work.d[35] = work.v[35];
  if (work.d[35] < 0)
    work.d[35] = settings.kkt_reg;
  else
    work.d[35] += settings.kkt_reg;
  work.d_inv[35] = 1/work.d[35];
  work.L[2] = (work.KKT[60])*work.d_inv[35];
  work.v[36] = work.KKT[61];
  work.d[36] = work.v[36];
  if (work.d[36] < 0)
    work.d[36] = settings.kkt_reg;
  else
    work.d[36] += settings.kkt_reg;
  work.d_inv[36] = 1/work.d[36];
  work.L[3] = (work.KKT[62])*work.d_inv[36];
  work.v[11] = work.L[0]*work.d[11];
  work.v[37] = work.KKT[63]-work.L[0]*work.v[11];
  work.d[37] = work.v[37];
  if (work.d[37] > 0)
    work.d[37] = -settings.kkt_reg;
  else
    work.d[37] -= settings.kkt_reg;
  work.d_inv[37] = 1/work.d[37];
  work.L[54] = (work.KKT[64])*work.d_inv[37];
  work.v[12] = work.L[1]*work.d[12];
  work.v[38] = work.KKT[65]-work.L[1]*work.v[12];
  work.d[38] = work.v[38];
  if (work.d[38] > 0)
    work.d[38] = -settings.kkt_reg;
  else
    work.d[38] -= settings.kkt_reg;
  work.d_inv[38] = 1/work.d[38];
  work.L[55] = (work.KKT[66])*work.d_inv[38];
  work.v[35] = work.L[2]*work.d[35];
  work.v[39] = work.KKT[67]-work.L[2]*work.v[35];
  work.d[39] = work.v[39];
  if (work.d[39] > 0)
    work.d[39] = -settings.kkt_reg;
  else
    work.d[39] -= settings.kkt_reg;
  work.d_inv[39] = 1/work.d[39];
  work.L[56] = (work.KKT[68])*work.d_inv[39];
  work.v[36] = work.L[3]*work.d[36];
  work.v[40] = work.KKT[69]-work.L[3]*work.v[36];
  work.d[40] = work.v[40];
  if (work.d[40] > 0)
    work.d[40] = -settings.kkt_reg;
  else
    work.d[40] -= settings.kkt_reg;
  work.d_inv[40] = 1/work.d[40];
  work.L[57] = (work.KKT[70])*work.d_inv[40];
  work.v[41] = work.KKT[71];
  work.d[41] = work.v[41];
  if (work.d[41] < 0)
    work.d[41] = settings.kkt_reg;
  else
    work.d[41] += settings.kkt_reg;
  work.d_inv[41] = 1/work.d[41];
  work.L[5] = (work.KKT[72])*work.d_inv[41];
  work.L[7] = (work.KKT[73])*work.d_inv[41];
  work.v[13] = work.L[4]*work.d[13];
  work.v[41] = work.L[5]*work.d[41];
  work.v[42] = work.KKT[74]-work.L[4]*work.v[13]-work.L[5]*work.v[41];
  work.d[42] = work.v[42];
  if (work.d[42] > 0)
    work.d[42] = -settings.kkt_reg;
  else
    work.d[42] -= settings.kkt_reg;
  work.d_inv[42] = 1/work.d[42];
  work.L[8] = (-work.L[7]*work.v[41])*work.d_inv[42];
  work.L[58] = (work.KKT[75])*work.d_inv[42];
  work.v[24] = work.L[6]*work.d[24];
  work.v[41] = work.L[7]*work.d[41];
  work.v[42] = work.L[8]*work.d[42];
  work.v[43] = work.KKT[76]-work.L[6]*work.v[24]-work.L[7]*work.v[41]-work.L[8]*work.v[42];
  work.d[43] = work.v[43];
  if (work.d[43] > 0)
    work.d[43] = -settings.kkt_reg;
  else
    work.d[43] -= settings.kkt_reg;
  work.d_inv[43] = 1/work.d[43];
  work.L[59] = (work.KKT[77]-work.L[58]*work.v[42])*work.d_inv[43];
  work.v[44] = work.KKT[78];
  work.d[44] = work.v[44];
  if (work.d[44] < 0)
    work.d[44] = settings.kkt_reg;
  else
    work.d[44] += settings.kkt_reg;
  work.d_inv[44] = 1/work.d[44];
  work.L[10] = (work.KKT[79])*work.d_inv[44];
  work.L[12] = (work.KKT[80])*work.d_inv[44];
  work.v[14] = work.L[9]*work.d[14];
  work.v[44] = work.L[10]*work.d[44];
  work.v[45] = work.KKT[81]-work.L[9]*work.v[14]-work.L[10]*work.v[44];
  work.d[45] = work.v[45];
  if (work.d[45] > 0)
    work.d[45] = -settings.kkt_reg;
  else
    work.d[45] -= settings.kkt_reg;
  work.d_inv[45] = 1/work.d[45];
  work.L[13] = (-work.L[12]*work.v[44])*work.d_inv[45];
  work.L[60] = (work.KKT[82])*work.d_inv[45];
  work.v[25] = work.L[11]*work.d[25];
  work.v[44] = work.L[12]*work.d[44];
  work.v[45] = work.L[13]*work.d[45];
  work.v[46] = work.KKT[83]-work.L[11]*work.v[25]-work.L[12]*work.v[44]-work.L[13]*work.v[45];
  work.d[46] = work.v[46];
  if (work.d[46] > 0)
    work.d[46] = -settings.kkt_reg;
  else
    work.d[46] -= settings.kkt_reg;
  work.d_inv[46] = 1/work.d[46];
  work.L[61] = (work.KKT[84]-work.L[60]*work.v[45])*work.d_inv[46];
  work.v[47] = work.KKT[85];
  work.d[47] = work.v[47];
  if (work.d[47] < 0)
    work.d[47] = settings.kkt_reg;
  else
    work.d[47] += settings.kkt_reg;
  work.d_inv[47] = 1/work.d[47];
  work.L[15] = (work.KKT[86])*work.d_inv[47];
  work.L[17] = (work.KKT[87])*work.d_inv[47];
  work.v[15] = work.L[14]*work.d[15];
  work.v[47] = work.L[15]*work.d[47];
  work.v[48] = work.KKT[88]-work.L[14]*work.v[15]-work.L[15]*work.v[47];
  work.d[48] = work.v[48];
  if (work.d[48] > 0)
    work.d[48] = -settings.kkt_reg;
  else
    work.d[48] -= settings.kkt_reg;
  work.d_inv[48] = 1/work.d[48];
  work.L[18] = (-work.L[17]*work.v[47])*work.d_inv[48];
  work.L[62] = (work.KKT[89])*work.d_inv[48];
  work.v[26] = work.L[16]*work.d[26];
  work.v[47] = work.L[17]*work.d[47];
  work.v[48] = work.L[18]*work.d[48];
  work.v[49] = work.KKT[90]-work.L[16]*work.v[26]-work.L[17]*work.v[47]-work.L[18]*work.v[48];
  work.d[49] = work.v[49];
  if (work.d[49] > 0)
    work.d[49] = -settings.kkt_reg;
  else
    work.d[49] -= settings.kkt_reg;
  work.d_inv[49] = 1/work.d[49];
  work.L[63] = (work.KKT[91]-work.L[62]*work.v[48])*work.d_inv[49];
  work.v[50] = work.KKT[92];
  work.d[50] = work.v[50];
  if (work.d[50] < 0)
    work.d[50] = settings.kkt_reg;
  else
    work.d[50] += settings.kkt_reg;
  work.d_inv[50] = 1/work.d[50];
  work.L[20] = (work.KKT[93])*work.d_inv[50];
  work.L[22] = (work.KKT[94])*work.d_inv[50];
  work.v[16] = work.L[19]*work.d[16];
  work.v[50] = work.L[20]*work.d[50];
  work.v[51] = work.KKT[95]-work.L[19]*work.v[16]-work.L[20]*work.v[50];
  work.d[51] = work.v[51];
  if (work.d[51] > 0)
    work.d[51] = -settings.kkt_reg;
  else
    work.d[51] -= settings.kkt_reg;
  work.d_inv[51] = 1/work.d[51];
  work.L[23] = (-work.L[22]*work.v[50])*work.d_inv[51];
  work.L[64] = (work.KKT[96])*work.d_inv[51];
  work.v[27] = work.L[21]*work.d[27];
  work.v[50] = work.L[22]*work.d[50];
  work.v[51] = work.L[23]*work.d[51];
  work.v[52] = work.KKT[97]-work.L[21]*work.v[27]-work.L[22]*work.v[50]-work.L[23]*work.v[51];
  work.d[52] = work.v[52];
  if (work.d[52] > 0)
    work.d[52] = -settings.kkt_reg;
  else
    work.d[52] -= settings.kkt_reg;
  work.d_inv[52] = 1/work.d[52];
  work.L[65] = (work.KKT[98]-work.L[64]*work.v[51])*work.d_inv[52];
  work.v[53] = work.KKT[99];
  work.d[53] = work.v[53];
  if (work.d[53] < 0)
    work.d[53] = settings.kkt_reg;
  else
    work.d[53] += settings.kkt_reg;
  work.d_inv[53] = 1/work.d[53];
  work.L[25] = (work.KKT[100])*work.d_inv[53];
  work.L[27] = (work.KKT[101])*work.d_inv[53];
  work.v[17] = work.L[24]*work.d[17];
  work.v[53] = work.L[25]*work.d[53];
  work.v[54] = work.KKT[102]-work.L[24]*work.v[17]-work.L[25]*work.v[53];
  work.d[54] = work.v[54];
  if (work.d[54] > 0)
    work.d[54] = -settings.kkt_reg;
  else
    work.d[54] -= settings.kkt_reg;
  work.d_inv[54] = 1/work.d[54];
  work.L[28] = (-work.L[27]*work.v[53])*work.d_inv[54];
  work.L[66] = (work.KKT[103])*work.d_inv[54];
  work.v[28] = work.L[26]*work.d[28];
  work.v[53] = work.L[27]*work.d[53];
  work.v[54] = work.L[28]*work.d[54];
  work.v[55] = work.KKT[104]-work.L[26]*work.v[28]-work.L[27]*work.v[53]-work.L[28]*work.v[54];
  work.d[55] = work.v[55];
  if (work.d[55] > 0)
    work.d[55] = -settings.kkt_reg;
  else
    work.d[55] -= settings.kkt_reg;
  work.d_inv[55] = 1/work.d[55];
  work.L[67] = (work.KKT[105]-work.L[66]*work.v[54])*work.d_inv[55];
  work.v[56] = work.KKT[106];
  work.d[56] = work.v[56];
  if (work.d[56] < 0)
    work.d[56] = settings.kkt_reg;
  else
    work.d[56] += settings.kkt_reg;
  work.d_inv[56] = 1/work.d[56];
  work.L[30] = (work.KKT[107])*work.d_inv[56];
  work.L[32] = (work.KKT[108])*work.d_inv[56];
  work.v[18] = work.L[29]*work.d[18];
  work.v[56] = work.L[30]*work.d[56];
  work.v[57] = work.KKT[109]-work.L[29]*work.v[18]-work.L[30]*work.v[56];
  work.d[57] = work.v[57];
  if (work.d[57] > 0)
    work.d[57] = -settings.kkt_reg;
  else
    work.d[57] -= settings.kkt_reg;
  work.d_inv[57] = 1/work.d[57];
  work.L[33] = (-work.L[32]*work.v[56])*work.d_inv[57];
  work.L[68] = (work.KKT[110])*work.d_inv[57];
  work.v[29] = work.L[31]*work.d[29];
  work.v[56] = work.L[32]*work.d[56];
  work.v[57] = work.L[33]*work.d[57];
  work.v[58] = work.KKT[111]-work.L[31]*work.v[29]-work.L[32]*work.v[56]-work.L[33]*work.v[57];
  work.d[58] = work.v[58];
  if (work.d[58] > 0)
    work.d[58] = -settings.kkt_reg;
  else
    work.d[58] -= settings.kkt_reg;
  work.d_inv[58] = 1/work.d[58];
  work.L[69] = (work.KKT[112]-work.L[68]*work.v[57])*work.d_inv[58];
  work.v[59] = work.KKT[113];
  work.d[59] = work.v[59];
  if (work.d[59] < 0)
    work.d[59] = settings.kkt_reg;
  else
    work.d[59] += settings.kkt_reg;
  work.d_inv[59] = 1/work.d[59];
  work.L[35] = (work.KKT[114])*work.d_inv[59];
  work.L[37] = (work.KKT[115])*work.d_inv[59];
  work.v[19] = work.L[34]*work.d[19];
  work.v[59] = work.L[35]*work.d[59];
  work.v[60] = work.KKT[116]-work.L[34]*work.v[19]-work.L[35]*work.v[59];
  work.d[60] = work.v[60];
  if (work.d[60] > 0)
    work.d[60] = -settings.kkt_reg;
  else
    work.d[60] -= settings.kkt_reg;
  work.d_inv[60] = 1/work.d[60];
  work.L[38] = (-work.L[37]*work.v[59])*work.d_inv[60];
  work.L[70] = (work.KKT[117])*work.d_inv[60];
  work.v[30] = work.L[36]*work.d[30];
  work.v[59] = work.L[37]*work.d[59];
  work.v[60] = work.L[38]*work.d[60];
  work.v[61] = work.KKT[118]-work.L[36]*work.v[30]-work.L[37]*work.v[59]-work.L[38]*work.v[60];
  work.d[61] = work.v[61];
  if (work.d[61] > 0)
    work.d[61] = -settings.kkt_reg;
  else
    work.d[61] -= settings.kkt_reg;
  work.d_inv[61] = 1/work.d[61];
  work.L[71] = (work.KKT[119]-work.L[70]*work.v[60])*work.d_inv[61];
  work.v[62] = work.KKT[120];
  work.d[62] = work.v[62];
  if (work.d[62] < 0)
    work.d[62] = settings.kkt_reg;
  else
    work.d[62] += settings.kkt_reg;
  work.d_inv[62] = 1/work.d[62];
  work.L[40] = (work.KKT[121])*work.d_inv[62];
  work.L[42] = (work.KKT[122])*work.d_inv[62];
  work.v[20] = work.L[39]*work.d[20];
  work.v[62] = work.L[40]*work.d[62];
  work.v[63] = work.KKT[123]-work.L[39]*work.v[20]-work.L[40]*work.v[62];
  work.d[63] = work.v[63];
  if (work.d[63] > 0)
    work.d[63] = -settings.kkt_reg;
  else
    work.d[63] -= settings.kkt_reg;
  work.d_inv[63] = 1/work.d[63];
  work.L[43] = (-work.L[42]*work.v[62])*work.d_inv[63];
  work.L[72] = (work.KKT[124])*work.d_inv[63];
  work.v[31] = work.L[41]*work.d[31];
  work.v[62] = work.L[42]*work.d[62];
  work.v[63] = work.L[43]*work.d[63];
  work.v[64] = work.KKT[125]-work.L[41]*work.v[31]-work.L[42]*work.v[62]-work.L[43]*work.v[63];
  work.d[64] = work.v[64];
  if (work.d[64] > 0)
    work.d[64] = -settings.kkt_reg;
  else
    work.d[64] -= settings.kkt_reg;
  work.d_inv[64] = 1/work.d[64];
  work.L[73] = (work.KKT[126]-work.L[72]*work.v[63])*work.d_inv[64];
  work.v[65] = work.KKT[127];
  work.d[65] = work.v[65];
  if (work.d[65] < 0)
    work.d[65] = settings.kkt_reg;
  else
    work.d[65] += settings.kkt_reg;
  work.d_inv[65] = 1/work.d[65];
  work.L[45] = (work.KKT[128])*work.d_inv[65];
  work.L[47] = (work.KKT[129])*work.d_inv[65];
  work.v[21] = work.L[44]*work.d[21];
  work.v[65] = work.L[45]*work.d[65];
  work.v[66] = work.KKT[130]-work.L[44]*work.v[21]-work.L[45]*work.v[65];
  work.d[66] = work.v[66];
  if (work.d[66] > 0)
    work.d[66] = -settings.kkt_reg;
  else
    work.d[66] -= settings.kkt_reg;
  work.d_inv[66] = 1/work.d[66];
  work.L[48] = (-work.L[47]*work.v[65])*work.d_inv[66];
  work.L[74] = (work.KKT[131])*work.d_inv[66];
  work.v[32] = work.L[46]*work.d[32];
  work.v[65] = work.L[47]*work.d[65];
  work.v[66] = work.L[48]*work.d[66];
  work.v[67] = work.KKT[132]-work.L[46]*work.v[32]-work.L[47]*work.v[65]-work.L[48]*work.v[66];
  work.d[67] = work.v[67];
  if (work.d[67] > 0)
    work.d[67] = -settings.kkt_reg;
  else
    work.d[67] -= settings.kkt_reg;
  work.d_inv[67] = 1/work.d[67];
  work.L[75] = (work.KKT[133]-work.L[74]*work.v[66])*work.d_inv[67];
  work.v[68] = work.KKT[134];
  work.d[68] = work.v[68];
  if (work.d[68] < 0)
    work.d[68] = settings.kkt_reg;
  else
    work.d[68] += settings.kkt_reg;
  work.d_inv[68] = 1/work.d[68];
  work.L[50] = (work.KKT[135])*work.d_inv[68];
  work.L[52] = (work.KKT[136])*work.d_inv[68];
  work.v[22] = work.L[49]*work.d[22];
  work.v[68] = work.L[50]*work.d[68];
  work.v[69] = work.KKT[137]-work.L[49]*work.v[22]-work.L[50]*work.v[68];
  work.d[69] = work.v[69];
  if (work.d[69] > 0)
    work.d[69] = -settings.kkt_reg;
  else
    work.d[69] -= settings.kkt_reg;
  work.d_inv[69] = 1/work.d[69];
  work.L[53] = (-work.L[52]*work.v[68])*work.d_inv[69];
  work.L[76] = (work.KKT[138])*work.d_inv[69];
  work.v[33] = work.L[51]*work.d[33];
  work.v[68] = work.L[52]*work.d[68];
  work.v[69] = work.L[53]*work.d[69];
  work.v[70] = work.KKT[139]-work.L[51]*work.v[33]-work.L[52]*work.v[68]-work.L[53]*work.v[69];
  work.d[70] = work.v[70];
  if (work.d[70] > 0)
    work.d[70] = -settings.kkt_reg;
  else
    work.d[70] -= settings.kkt_reg;
  work.d_inv[70] = 1/work.d[70];
  work.L[77] = (work.KKT[140]-work.L[76]*work.v[69])*work.d_inv[70];
  work.v[37] = work.L[54]*work.d[37];
  work.v[38] = work.L[55]*work.d[38];
  work.v[39] = work.L[56]*work.d[39];
  work.v[40] = work.L[57]*work.d[40];
  work.v[42] = work.L[58]*work.d[42];
  work.v[43] = work.L[59]*work.d[43];
  work.v[45] = work.L[60]*work.d[45];
  work.v[46] = work.L[61]*work.d[46];
  work.v[48] = work.L[62]*work.d[48];
  work.v[49] = work.L[63]*work.d[49];
  work.v[51] = work.L[64]*work.d[51];
  work.v[52] = work.L[65]*work.d[52];
  work.v[54] = work.L[66]*work.d[54];
  work.v[55] = work.L[67]*work.d[55];
  work.v[57] = work.L[68]*work.d[57];
  work.v[58] = work.L[69]*work.d[58];
  work.v[60] = work.L[70]*work.d[60];
  work.v[61] = work.L[71]*work.d[61];
  work.v[63] = work.L[72]*work.d[63];
  work.v[64] = work.L[73]*work.d[64];
  work.v[66] = work.L[74]*work.d[66];
  work.v[67] = work.L[75]*work.d[67];
  work.v[69] = work.L[76]*work.d[69];
  work.v[70] = work.L[77]*work.d[70];
  work.v[71] = work.KKT[141]-work.L[54]*work.v[37]-work.L[55]*work.v[38]-work.L[56]*work.v[39]-work.L[57]*work.v[40]-work.L[58]*work.v[42]-work.L[59]*work.v[43]-work.L[60]*work.v[45]-work.L[61]*work.v[46]-work.L[62]*work.v[48]-work.L[63]*work.v[49]-work.L[64]*work.v[51]-work.L[65]*work.v[52]-work.L[66]*work.v[54]-work.L[67]*work.v[55]-work.L[68]*work.v[57]-work.L[69]*work.v[58]-work.L[70]*work.v[60]-work.L[71]*work.v[61]-work.L[72]*work.v[63]-work.L[73]*work.v[64]-work.L[74]*work.v[66]-work.L[75]*work.v[67]-work.L[76]*work.v[69]-work.L[77]*work.v[70];
  work.d[71] = work.v[71];
  if (work.d[71] < 0)
    work.d[71] = settings.kkt_reg;
  else
    work.d[71] += settings.kkt_reg;
  work.d_inv[71] = 1/work.d[71];
  work.L[79] = (work.KKT[142])*work.d_inv[71];
  work.L[82] = (work.KKT[143])*work.d_inv[71];
  work.v[72] = work.KKT[144];
  work.d[72] = work.v[72];
  if (work.d[72] < 0)
    work.d[72] = settings.kkt_reg;
  else
    work.d[72] += settings.kkt_reg;
  work.d_inv[72] = 1/work.d[72];
  work.L[80] = (work.KKT[145])*work.d_inv[72];
  work.L[83] = (work.KKT[146])*work.d_inv[72];
  work.v[23] = work.L[78]*work.d[23];
  work.v[71] = work.L[79]*work.d[71];
  work.v[72] = work.L[80]*work.d[72];
  work.v[73] = work.KKT[147]-work.L[78]*work.v[23]-work.L[79]*work.v[71]-work.L[80]*work.v[72];
  work.d[73] = work.v[73];
  if (work.d[73] > 0)
    work.d[73] = -settings.kkt_reg;
  else
    work.d[73] -= settings.kkt_reg;
  work.d_inv[73] = 1/work.d[73];
  work.L[84] = (-work.L[82]*work.v[71]-work.L[83]*work.v[72])*work.d_inv[73];
  work.v[34] = work.L[81]*work.d[34];
  work.v[71] = work.L[82]*work.d[71];
  work.v[72] = work.L[83]*work.d[72];
  work.v[73] = work.L[84]*work.d[73];
  work.v[74] = work.KKT[148]-work.L[81]*work.v[34]-work.L[82]*work.v[71]-work.L[83]*work.v[72]-work.L[84]*work.v[73];
  work.d[74] = work.v[74];
  if (work.d[74] > 0)
    work.d[74] = -settings.kkt_reg;
  else
    work.d[74] -= settings.kkt_reg;
  work.d_inv[74] = 1/work.d[74];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared Frobenius for factorization is %.8g.\n", check_factorization());
  }
#endif
}
double GensolverFirst::check_factorization(void) {
  /* Returns the squared Frobenius norm of A - L*D*L'. */
  double temp, residual;
  /* Only check the lower triangle. */
  residual = 0;
  temp = work.KKT[141]-1*work.d[71]*1-work.L[54]*work.d[37]*work.L[54]-work.L[55]*work.d[38]*work.L[55]-work.L[58]*work.d[42]*work.L[58]-work.L[60]*work.d[45]*work.L[60]-work.L[62]*work.d[48]*work.L[62]-work.L[64]*work.d[51]*work.L[64]-work.L[66]*work.d[54]*work.L[66]-work.L[68]*work.d[57]*work.L[68]-work.L[70]*work.d[60]*work.L[70]-work.L[72]*work.d[63]*work.L[72]-work.L[74]*work.d[66]*work.L[74]-work.L[76]*work.d[69]*work.L[76]-work.L[59]*work.d[43]*work.L[59]-work.L[61]*work.d[46]*work.L[61]-work.L[63]*work.d[49]*work.L[63]-work.L[65]*work.d[52]*work.L[65]-work.L[67]*work.d[55]*work.L[67]-work.L[69]*work.d[58]*work.L[69]-work.L[71]*work.d[61]*work.L[71]-work.L[73]*work.d[64]*work.L[73]-work.L[75]*work.d[67]*work.L[75]-work.L[77]*work.d[70]*work.L[77]-work.L[56]*work.d[39]*work.L[56]-work.L[57]*work.d[40]*work.L[57];
  residual += temp*temp;
  temp = work.KKT[71]-1*work.d[41]*1;
  residual += temp*temp;
  temp = work.KKT[78]-1*work.d[44]*1;
  residual += temp*temp;
  temp = work.KKT[85]-1*work.d[47]*1;
  residual += temp*temp;
  temp = work.KKT[92]-1*work.d[50]*1;
  residual += temp*temp;
  temp = work.KKT[99]-1*work.d[53]*1;
  residual += temp*temp;
  temp = work.KKT[106]-1*work.d[56]*1;
  residual += temp*temp;
  temp = work.KKT[113]-1*work.d[59]*1;
  residual += temp*temp;
  temp = work.KKT[120]-1*work.d[62]*1;
  residual += temp*temp;
  temp = work.KKT[127]-1*work.d[65]*1;
  residual += temp*temp;
  temp = work.KKT[134]-1*work.d[68]*1;
  residual += temp*temp;
  temp = work.KKT[144]-1*work.d[72]*1;
  residual += temp*temp;
  temp = work.KKT[0]-1*work.d[0]*1;
  residual += temp*temp;
  temp = work.KKT[1]-1*work.d[1]*1;
  residual += temp*temp;
  temp = work.KKT[2]-1*work.d[2]*1;
  residual += temp*temp;
  temp = work.KKT[3]-1*work.d[3]*1;
  residual += temp*temp;
  temp = work.KKT[4]-1*work.d[4]*1;
  residual += temp*temp;
  temp = work.KKT[5]-1*work.d[5]*1;
  residual += temp*temp;
  temp = work.KKT[6]-1*work.d[6]*1;
  residual += temp*temp;
  temp = work.KKT[7]-1*work.d[7]*1;
  residual += temp*temp;
  temp = work.KKT[8]-1*work.d[8]*1;
  residual += temp*temp;
  temp = work.KKT[9]-1*work.d[9]*1;
  residual += temp*temp;
  temp = work.KKT[10]-1*work.d[10]*1;
  residual += temp*temp;
  temp = work.KKT[11]-1*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[13]-1*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[15]-1*work.d[13]*1;
  residual += temp*temp;
  temp = work.KKT[17]-1*work.d[14]*1;
  residual += temp*temp;
  temp = work.KKT[19]-1*work.d[15]*1;
  residual += temp*temp;
  temp = work.KKT[21]-1*work.d[16]*1;
  residual += temp*temp;
  temp = work.KKT[23]-1*work.d[17]*1;
  residual += temp*temp;
  temp = work.KKT[25]-1*work.d[18]*1;
  residual += temp*temp;
  temp = work.KKT[27]-1*work.d[19]*1;
  residual += temp*temp;
  temp = work.KKT[29]-1*work.d[20]*1;
  residual += temp*temp;
  temp = work.KKT[31]-1*work.d[21]*1;
  residual += temp*temp;
  temp = work.KKT[33]-1*work.d[22]*1;
  residual += temp*temp;
  temp = work.KKT[35]-1*work.d[23]*1;
  residual += temp*temp;
  temp = work.KKT[37]-1*work.d[24]*1;
  residual += temp*temp;
  temp = work.KKT[39]-1*work.d[25]*1;
  residual += temp*temp;
  temp = work.KKT[41]-1*work.d[26]*1;
  residual += temp*temp;
  temp = work.KKT[43]-1*work.d[27]*1;
  residual += temp*temp;
  temp = work.KKT[45]-1*work.d[28]*1;
  residual += temp*temp;
  temp = work.KKT[47]-1*work.d[29]*1;
  residual += temp*temp;
  temp = work.KKT[49]-1*work.d[30]*1;
  residual += temp*temp;
  temp = work.KKT[51]-1*work.d[31]*1;
  residual += temp*temp;
  temp = work.KKT[53]-1*work.d[32]*1;
  residual += temp*temp;
  temp = work.KKT[55]-1*work.d[33]*1;
  residual += temp*temp;
  temp = work.KKT[57]-1*work.d[34]*1;
  residual += temp*temp;
  temp = work.KKT[59]-1*work.d[35]*1;
  residual += temp*temp;
  temp = work.KKT[61]-1*work.d[36]*1;
  residual += temp*temp;
  temp = work.KKT[12]-work.L[0]*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[14]-work.L[1]*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[16]-work.L[4]*work.d[13]*1;
  residual += temp*temp;
  temp = work.KKT[18]-work.L[9]*work.d[14]*1;
  residual += temp*temp;
  temp = work.KKT[20]-work.L[14]*work.d[15]*1;
  residual += temp*temp;
  temp = work.KKT[22]-work.L[19]*work.d[16]*1;
  residual += temp*temp;
  temp = work.KKT[24]-work.L[24]*work.d[17]*1;
  residual += temp*temp;
  temp = work.KKT[26]-work.L[29]*work.d[18]*1;
  residual += temp*temp;
  temp = work.KKT[28]-work.L[34]*work.d[19]*1;
  residual += temp*temp;
  temp = work.KKT[30]-work.L[39]*work.d[20]*1;
  residual += temp*temp;
  temp = work.KKT[32]-work.L[44]*work.d[21]*1;
  residual += temp*temp;
  temp = work.KKT[34]-work.L[49]*work.d[22]*1;
  residual += temp*temp;
  temp = work.KKT[36]-work.L[78]*work.d[23]*1;
  residual += temp*temp;
  temp = work.KKT[38]-work.L[6]*work.d[24]*1;
  residual += temp*temp;
  temp = work.KKT[40]-work.L[11]*work.d[25]*1;
  residual += temp*temp;
  temp = work.KKT[42]-work.L[16]*work.d[26]*1;
  residual += temp*temp;
  temp = work.KKT[44]-work.L[21]*work.d[27]*1;
  residual += temp*temp;
  temp = work.KKT[46]-work.L[26]*work.d[28]*1;
  residual += temp*temp;
  temp = work.KKT[48]-work.L[31]*work.d[29]*1;
  residual += temp*temp;
  temp = work.KKT[50]-work.L[36]*work.d[30]*1;
  residual += temp*temp;
  temp = work.KKT[52]-work.L[41]*work.d[31]*1;
  residual += temp*temp;
  temp = work.KKT[54]-work.L[46]*work.d[32]*1;
  residual += temp*temp;
  temp = work.KKT[56]-work.L[51]*work.d[33]*1;
  residual += temp*temp;
  temp = work.KKT[58]-work.L[81]*work.d[34]*1;
  residual += temp*temp;
  temp = work.KKT[60]-work.L[2]*work.d[35]*1;
  residual += temp*temp;
  temp = work.KKT[62]-work.L[3]*work.d[36]*1;
  residual += temp*temp;
  temp = work.KKT[63]-work.L[0]*work.d[11]*work.L[0]-1*work.d[37]*1;
  residual += temp*temp;
  temp = work.KKT[65]-work.L[1]*work.d[12]*work.L[1]-1*work.d[38]*1;
  residual += temp*temp;
  temp = work.KKT[74]-work.L[4]*work.d[13]*work.L[4]-1*work.d[42]*1-work.L[5]*work.d[41]*work.L[5];
  residual += temp*temp;
  temp = work.KKT[81]-work.L[9]*work.d[14]*work.L[9]-1*work.d[45]*1-work.L[10]*work.d[44]*work.L[10];
  residual += temp*temp;
  temp = work.KKT[88]-work.L[14]*work.d[15]*work.L[14]-1*work.d[48]*1-work.L[15]*work.d[47]*work.L[15];
  residual += temp*temp;
  temp = work.KKT[95]-work.L[19]*work.d[16]*work.L[19]-1*work.d[51]*1-work.L[20]*work.d[50]*work.L[20];
  residual += temp*temp;
  temp = work.KKT[102]-work.L[24]*work.d[17]*work.L[24]-1*work.d[54]*1-work.L[25]*work.d[53]*work.L[25];
  residual += temp*temp;
  temp = work.KKT[109]-work.L[29]*work.d[18]*work.L[29]-1*work.d[57]*1-work.L[30]*work.d[56]*work.L[30];
  residual += temp*temp;
  temp = work.KKT[116]-work.L[34]*work.d[19]*work.L[34]-1*work.d[60]*1-work.L[35]*work.d[59]*work.L[35];
  residual += temp*temp;
  temp = work.KKT[123]-work.L[39]*work.d[20]*work.L[39]-1*work.d[63]*1-work.L[40]*work.d[62]*work.L[40];
  residual += temp*temp;
  temp = work.KKT[130]-work.L[44]*work.d[21]*work.L[44]-1*work.d[66]*1-work.L[45]*work.d[65]*work.L[45];
  residual += temp*temp;
  temp = work.KKT[137]-work.L[49]*work.d[22]*work.L[49]-1*work.d[69]*1-work.L[50]*work.d[68]*work.L[50];
  residual += temp*temp;
  temp = work.KKT[147]-work.L[78]*work.d[23]*work.L[78]-1*work.d[73]*1-work.L[79]*work.d[71]*work.L[79]-work.L[80]*work.d[72]*work.L[80];
  residual += temp*temp;
  temp = work.KKT[76]-work.L[6]*work.d[24]*work.L[6]-1*work.d[43]*1-work.L[7]*work.d[41]*work.L[7]-work.L[8]*work.d[42]*work.L[8];
  residual += temp*temp;
  temp = work.KKT[83]-work.L[11]*work.d[25]*work.L[11]-1*work.d[46]*1-work.L[12]*work.d[44]*work.L[12]-work.L[13]*work.d[45]*work.L[13];
  residual += temp*temp;
  temp = work.KKT[90]-work.L[16]*work.d[26]*work.L[16]-1*work.d[49]*1-work.L[17]*work.d[47]*work.L[17]-work.L[18]*work.d[48]*work.L[18];
  residual += temp*temp;
  temp = work.KKT[97]-work.L[21]*work.d[27]*work.L[21]-1*work.d[52]*1-work.L[22]*work.d[50]*work.L[22]-work.L[23]*work.d[51]*work.L[23];
  residual += temp*temp;
  temp = work.KKT[104]-work.L[26]*work.d[28]*work.L[26]-1*work.d[55]*1-work.L[27]*work.d[53]*work.L[27]-work.L[28]*work.d[54]*work.L[28];
  residual += temp*temp;
  temp = work.KKT[111]-work.L[31]*work.d[29]*work.L[31]-1*work.d[58]*1-work.L[32]*work.d[56]*work.L[32]-work.L[33]*work.d[57]*work.L[33];
  residual += temp*temp;
  temp = work.KKT[118]-work.L[36]*work.d[30]*work.L[36]-1*work.d[61]*1-work.L[37]*work.d[59]*work.L[37]-work.L[38]*work.d[60]*work.L[38];
  residual += temp*temp;
  temp = work.KKT[125]-work.L[41]*work.d[31]*work.L[41]-1*work.d[64]*1-work.L[42]*work.d[62]*work.L[42]-work.L[43]*work.d[63]*work.L[43];
  residual += temp*temp;
  temp = work.KKT[132]-work.L[46]*work.d[32]*work.L[46]-1*work.d[67]*1-work.L[47]*work.d[65]*work.L[47]-work.L[48]*work.d[66]*work.L[48];
  residual += temp*temp;
  temp = work.KKT[139]-work.L[51]*work.d[33]*work.L[51]-1*work.d[70]*1-work.L[52]*work.d[68]*work.L[52]-work.L[53]*work.d[69]*work.L[53];
  residual += temp*temp;
  temp = work.KKT[148]-work.L[81]*work.d[34]*work.L[81]-1*work.d[74]*1-work.L[82]*work.d[71]*work.L[82]-work.L[83]*work.d[72]*work.L[83]-work.L[84]*work.d[73]*work.L[84];
  residual += temp*temp;
  temp = work.KKT[67]-work.L[2]*work.d[35]*work.L[2]-1*work.d[39]*1;
  residual += temp*temp;
  temp = work.KKT[69]-work.L[3]*work.d[36]*work.L[3]-1*work.d[40]*1;
  residual += temp*temp;
  temp = work.KKT[64]-1*work.d[37]*work.L[54];
  residual += temp*temp;
  temp = work.KKT[66]-1*work.d[38]*work.L[55];
  residual += temp*temp;
  temp = work.KKT[75]-1*work.d[42]*work.L[58];
  residual += temp*temp;
  temp = work.KKT[82]-1*work.d[45]*work.L[60];
  residual += temp*temp;
  temp = work.KKT[89]-1*work.d[48]*work.L[62];
  residual += temp*temp;
  temp = work.KKT[96]-1*work.d[51]*work.L[64];
  residual += temp*temp;
  temp = work.KKT[103]-1*work.d[54]*work.L[66];
  residual += temp*temp;
  temp = work.KKT[110]-1*work.d[57]*work.L[68];
  residual += temp*temp;
  temp = work.KKT[117]-1*work.d[60]*work.L[70];
  residual += temp*temp;
  temp = work.KKT[124]-1*work.d[63]*work.L[72];
  residual += temp*temp;
  temp = work.KKT[131]-1*work.d[66]*work.L[74];
  residual += temp*temp;
  temp = work.KKT[138]-1*work.d[69]*work.L[76];
  residual += temp*temp;
  temp = work.KKT[142]-work.L[79]*work.d[71]*1;
  residual += temp*temp;
  temp = work.KKT[72]-work.L[5]*work.d[41]*1;
  residual += temp*temp;
  temp = work.KKT[79]-work.L[10]*work.d[44]*1;
  residual += temp*temp;
  temp = work.KKT[86]-work.L[15]*work.d[47]*1;
  residual += temp*temp;
  temp = work.KKT[93]-work.L[20]*work.d[50]*1;
  residual += temp*temp;
  temp = work.KKT[100]-work.L[25]*work.d[53]*1;
  residual += temp*temp;
  temp = work.KKT[107]-work.L[30]*work.d[56]*1;
  residual += temp*temp;
  temp = work.KKT[114]-work.L[35]*work.d[59]*1;
  residual += temp*temp;
  temp = work.KKT[121]-work.L[40]*work.d[62]*1;
  residual += temp*temp;
  temp = work.KKT[128]-work.L[45]*work.d[65]*1;
  residual += temp*temp;
  temp = work.KKT[135]-work.L[50]*work.d[68]*1;
  residual += temp*temp;
  temp = work.KKT[145]-work.L[80]*work.d[72]*1;
  residual += temp*temp;
  temp = work.KKT[77]-1*work.d[43]*work.L[59]-work.L[8]*work.d[42]*work.L[58];
  residual += temp*temp;
  temp = work.KKT[84]-1*work.d[46]*work.L[61]-work.L[13]*work.d[45]*work.L[60];
  residual += temp*temp;
  temp = work.KKT[91]-1*work.d[49]*work.L[63]-work.L[18]*work.d[48]*work.L[62];
  residual += temp*temp;
  temp = work.KKT[98]-1*work.d[52]*work.L[65]-work.L[23]*work.d[51]*work.L[64];
  residual += temp*temp;
  temp = work.KKT[105]-1*work.d[55]*work.L[67]-work.L[28]*work.d[54]*work.L[66];
  residual += temp*temp;
  temp = work.KKT[112]-1*work.d[58]*work.L[69]-work.L[33]*work.d[57]*work.L[68];
  residual += temp*temp;
  temp = work.KKT[119]-1*work.d[61]*work.L[71]-work.L[38]*work.d[60]*work.L[70];
  residual += temp*temp;
  temp = work.KKT[126]-1*work.d[64]*work.L[73]-work.L[43]*work.d[63]*work.L[72];
  residual += temp*temp;
  temp = work.KKT[133]-1*work.d[67]*work.L[75]-work.L[48]*work.d[66]*work.L[74];
  residual += temp*temp;
  temp = work.KKT[140]-1*work.d[70]*work.L[77]-work.L[53]*work.d[69]*work.L[76];
  residual += temp*temp;
  temp = work.KKT[143]-work.L[82]*work.d[71]*1;
  residual += temp*temp;
  temp = work.KKT[73]-work.L[7]*work.d[41]*1;
  residual += temp*temp;
  temp = work.KKT[80]-work.L[12]*work.d[44]*1;
  residual += temp*temp;
  temp = work.KKT[87]-work.L[17]*work.d[47]*1;
  residual += temp*temp;
  temp = work.KKT[94]-work.L[22]*work.d[50]*1;
  residual += temp*temp;
  temp = work.KKT[101]-work.L[27]*work.d[53]*1;
  residual += temp*temp;
  temp = work.KKT[108]-work.L[32]*work.d[56]*1;
  residual += temp*temp;
  temp = work.KKT[115]-work.L[37]*work.d[59]*1;
  residual += temp*temp;
  temp = work.KKT[122]-work.L[42]*work.d[62]*1;
  residual += temp*temp;
  temp = work.KKT[129]-work.L[47]*work.d[65]*1;
  residual += temp*temp;
  temp = work.KKT[136]-work.L[52]*work.d[68]*1;
  residual += temp*temp;
  temp = work.KKT[146]-work.L[83]*work.d[72]*1;
  residual += temp*temp;
  temp = work.KKT[68]-1*work.d[39]*work.L[56];
  residual += temp*temp;
  temp = work.KKT[70]-1*work.d[40]*work.L[57];
  residual += temp*temp;
  return residual;
}
void GensolverFirst::matrix_multiply(double *result, double *source) {
  /* Finds result = A*source. */
  result[0] = work.KKT[141]*source[0]+work.KKT[64]*source[49]+work.KKT[66]*source[50]+work.KKT[75]*source[51]+work.KKT[82]*source[52]+work.KKT[89]*source[53]+work.KKT[96]*source[54]+work.KKT[103]*source[55]+work.KKT[110]*source[56]+work.KKT[117]*source[57]+work.KKT[124]*source[58]+work.KKT[131]*source[59]+work.KKT[138]*source[60]+work.KKT[142]*source[61]+work.KKT[77]*source[62]+work.KKT[84]*source[63]+work.KKT[91]*source[64]+work.KKT[98]*source[65]+work.KKT[105]*source[66]+work.KKT[112]*source[67]+work.KKT[119]*source[68]+work.KKT[126]*source[69]+work.KKT[133]*source[70]+work.KKT[140]*source[71]+work.KKT[143]*source[72]+work.KKT[68]*source[73]+work.KKT[70]*source[74];
  result[1] = work.KKT[71]*source[1]+work.KKT[72]*source[51]+work.KKT[73]*source[62];
  result[2] = work.KKT[78]*source[2]+work.KKT[79]*source[52]+work.KKT[80]*source[63];
  result[3] = work.KKT[85]*source[3]+work.KKT[86]*source[53]+work.KKT[87]*source[64];
  result[4] = work.KKT[92]*source[4]+work.KKT[93]*source[54]+work.KKT[94]*source[65];
  result[5] = work.KKT[99]*source[5]+work.KKT[100]*source[55]+work.KKT[101]*source[66];
  result[6] = work.KKT[106]*source[6]+work.KKT[107]*source[56]+work.KKT[108]*source[67];
  result[7] = work.KKT[113]*source[7]+work.KKT[114]*source[57]+work.KKT[115]*source[68];
  result[8] = work.KKT[120]*source[8]+work.KKT[121]*source[58]+work.KKT[122]*source[69];
  result[9] = work.KKT[127]*source[9]+work.KKT[128]*source[59]+work.KKT[129]*source[70];
  result[10] = work.KKT[134]*source[10]+work.KKT[135]*source[60]+work.KKT[136]*source[71];
  result[11] = work.KKT[144]*source[11]+work.KKT[145]*source[61]+work.KKT[146]*source[72];
  result[12] = work.KKT[0]*source[12];
  result[13] = work.KKT[1]*source[13];
  result[14] = work.KKT[2]*source[14];
  result[15] = work.KKT[3]*source[15];
  result[16] = work.KKT[4]*source[16];
  result[17] = work.KKT[5]*source[17];
  result[18] = work.KKT[6]*source[18];
  result[19] = work.KKT[7]*source[19];
  result[20] = work.KKT[8]*source[20];
  result[21] = work.KKT[9]*source[21];
  result[22] = work.KKT[10]*source[22];
  result[23] = work.KKT[11]*source[23]+work.KKT[12]*source[49];
  result[24] = work.KKT[13]*source[24]+work.KKT[14]*source[50];
  result[25] = work.KKT[15]*source[25]+work.KKT[16]*source[51];
  result[26] = work.KKT[17]*source[26]+work.KKT[18]*source[52];
  result[27] = work.KKT[19]*source[27]+work.KKT[20]*source[53];
  result[28] = work.KKT[21]*source[28]+work.KKT[22]*source[54];
  result[29] = work.KKT[23]*source[29]+work.KKT[24]*source[55];
  result[30] = work.KKT[25]*source[30]+work.KKT[26]*source[56];
  result[31] = work.KKT[27]*source[31]+work.KKT[28]*source[57];
  result[32] = work.KKT[29]*source[32]+work.KKT[30]*source[58];
  result[33] = work.KKT[31]*source[33]+work.KKT[32]*source[59];
  result[34] = work.KKT[33]*source[34]+work.KKT[34]*source[60];
  result[35] = work.KKT[35]*source[35]+work.KKT[36]*source[61];
  result[36] = work.KKT[37]*source[36]+work.KKT[38]*source[62];
  result[37] = work.KKT[39]*source[37]+work.KKT[40]*source[63];
  result[38] = work.KKT[41]*source[38]+work.KKT[42]*source[64];
  result[39] = work.KKT[43]*source[39]+work.KKT[44]*source[65];
  result[40] = work.KKT[45]*source[40]+work.KKT[46]*source[66];
  result[41] = work.KKT[47]*source[41]+work.KKT[48]*source[67];
  result[42] = work.KKT[49]*source[42]+work.KKT[50]*source[68];
  result[43] = work.KKT[51]*source[43]+work.KKT[52]*source[69];
  result[44] = work.KKT[53]*source[44]+work.KKT[54]*source[70];
  result[45] = work.KKT[55]*source[45]+work.KKT[56]*source[71];
  result[46] = work.KKT[57]*source[46]+work.KKT[58]*source[72];
  result[47] = work.KKT[59]*source[47]+work.KKT[60]*source[73];
  result[48] = work.KKT[61]*source[48]+work.KKT[62]*source[74];
  result[49] = work.KKT[12]*source[23]+work.KKT[63]*source[49]+work.KKT[64]*source[0];
  result[50] = work.KKT[14]*source[24]+work.KKT[65]*source[50]+work.KKT[66]*source[0];
  result[51] = work.KKT[16]*source[25]+work.KKT[74]*source[51]+work.KKT[75]*source[0]+work.KKT[72]*source[1];
  result[52] = work.KKT[18]*source[26]+work.KKT[81]*source[52]+work.KKT[82]*source[0]+work.KKT[79]*source[2];
  result[53] = work.KKT[20]*source[27]+work.KKT[88]*source[53]+work.KKT[89]*source[0]+work.KKT[86]*source[3];
  result[54] = work.KKT[22]*source[28]+work.KKT[95]*source[54]+work.KKT[96]*source[0]+work.KKT[93]*source[4];
  result[55] = work.KKT[24]*source[29]+work.KKT[102]*source[55]+work.KKT[103]*source[0]+work.KKT[100]*source[5];
  result[56] = work.KKT[26]*source[30]+work.KKT[109]*source[56]+work.KKT[110]*source[0]+work.KKT[107]*source[6];
  result[57] = work.KKT[28]*source[31]+work.KKT[116]*source[57]+work.KKT[117]*source[0]+work.KKT[114]*source[7];
  result[58] = work.KKT[30]*source[32]+work.KKT[123]*source[58]+work.KKT[124]*source[0]+work.KKT[121]*source[8];
  result[59] = work.KKT[32]*source[33]+work.KKT[130]*source[59]+work.KKT[131]*source[0]+work.KKT[128]*source[9];
  result[60] = work.KKT[34]*source[34]+work.KKT[137]*source[60]+work.KKT[138]*source[0]+work.KKT[135]*source[10];
  result[61] = work.KKT[36]*source[35]+work.KKT[147]*source[61]+work.KKT[142]*source[0]+work.KKT[145]*source[11];
  result[62] = work.KKT[38]*source[36]+work.KKT[76]*source[62]+work.KKT[77]*source[0]+work.KKT[73]*source[1];
  result[63] = work.KKT[40]*source[37]+work.KKT[83]*source[63]+work.KKT[84]*source[0]+work.KKT[80]*source[2];
  result[64] = work.KKT[42]*source[38]+work.KKT[90]*source[64]+work.KKT[91]*source[0]+work.KKT[87]*source[3];
  result[65] = work.KKT[44]*source[39]+work.KKT[97]*source[65]+work.KKT[98]*source[0]+work.KKT[94]*source[4];
  result[66] = work.KKT[46]*source[40]+work.KKT[104]*source[66]+work.KKT[105]*source[0]+work.KKT[101]*source[5];
  result[67] = work.KKT[48]*source[41]+work.KKT[111]*source[67]+work.KKT[112]*source[0]+work.KKT[108]*source[6];
  result[68] = work.KKT[50]*source[42]+work.KKT[118]*source[68]+work.KKT[119]*source[0]+work.KKT[115]*source[7];
  result[69] = work.KKT[52]*source[43]+work.KKT[125]*source[69]+work.KKT[126]*source[0]+work.KKT[122]*source[8];
  result[70] = work.KKT[54]*source[44]+work.KKT[132]*source[70]+work.KKT[133]*source[0]+work.KKT[129]*source[9];
  result[71] = work.KKT[56]*source[45]+work.KKT[139]*source[71]+work.KKT[140]*source[0]+work.KKT[136]*source[10];
  result[72] = work.KKT[58]*source[46]+work.KKT[148]*source[72]+work.KKT[143]*source[0]+work.KKT[146]*source[11];
  result[73] = work.KKT[60]*source[47]+work.KKT[67]*source[73]+work.KKT[68]*source[0];
  result[74] = work.KKT[62]*source[48]+work.KKT[69]*source[74]+work.KKT[70]*source[0];
}
double GensolverFirst::check_residual(double *target, double *multiplicand) {
  /* Returns the squared 2-norm of lhs - A*rhs. */
  /* Reuses v to find the residual. */
  int i;
  double residual;
  residual = 0;
  matrix_multiply(work.v, multiplicand);
  for (i = 0; i < 23; i++) {
    residual += (target[i] - work.v[i])*(target[i] - work.v[i]);
  }
  return residual;
}
void GensolverFirst::fill_KKT(void) {
  work.KKT[141] = 2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848*work.quad_877793640448[0]);
  work.KKT[71] = 2*work.frac_479733190656;
  work.KKT[78] = 2*work.frac_479733190656;
  work.KKT[85] = 2*work.frac_479733190656;
  work.KKT[92] = 2*work.frac_479733190656;
  work.KKT[99] = 2*work.frac_479733190656;
  work.KKT[106] = 2*work.frac_479733190656;
  work.KKT[113] = 2*work.frac_479733190656;
  work.KKT[120] = 2*work.frac_479733190656;
  work.KKT[127] = 2*work.frac_479733190656;
  work.KKT[134] = 2*work.frac_479733190656;
  work.KKT[144] = 2*work.frac_479733190656;
  work.KKT[0] = 2*work.frac_121674190848;
  work.KKT[1] = 2*work.frac_121674190848;
  work.KKT[2] = 2*work.frac_121674190848;
  work.KKT[3] = 2*work.frac_121674190848;
  work.KKT[4] = 2*work.frac_121674190848;
  work.KKT[5] = 2*work.frac_121674190848;
  work.KKT[6] = 2*work.frac_121674190848;
  work.KKT[7] = 2*work.frac_121674190848;
  work.KKT[8] = 2*work.frac_121674190848;
  work.KKT[9] = 2*work.frac_121674190848;
  work.KKT[10] = 2*work.frac_121674190848;
  work.KKT[11] = work.s_inv_z[0];
  work.KKT[13] = work.s_inv_z[1];
  work.KKT[15] = work.s_inv_z[2];
  work.KKT[17] = work.s_inv_z[3];
  work.KKT[19] = work.s_inv_z[4];
  work.KKT[21] = work.s_inv_z[5];
  work.KKT[23] = work.s_inv_z[6];
  work.KKT[25] = work.s_inv_z[7];
  work.KKT[27] = work.s_inv_z[8];
  work.KKT[29] = work.s_inv_z[9];
  work.KKT[31] = work.s_inv_z[10];
  work.KKT[33] = work.s_inv_z[11];
  work.KKT[35] = work.s_inv_z[12];
  work.KKT[37] = work.s_inv_z[13];
  work.KKT[39] = work.s_inv_z[14];
  work.KKT[41] = work.s_inv_z[15];
  work.KKT[43] = work.s_inv_z[16];
  work.KKT[45] = work.s_inv_z[17];
  work.KKT[47] = work.s_inv_z[18];
  work.KKT[49] = work.s_inv_z[19];
  work.KKT[51] = work.s_inv_z[20];
  work.KKT[53] = work.s_inv_z[21];
  work.KKT[55] = work.s_inv_z[22];
  work.KKT[57] = work.s_inv_z[23];
  work.KKT[59] = work.s_inv_z[24];
  work.KKT[61] = work.s_inv_z[25];
  work.KKT[12] = 1;
  work.KKT[14] = 1;
  work.KKT[16] = 1;
  work.KKT[18] = 1;
  work.KKT[20] = 1;
  work.KKT[22] = 1;
  work.KKT[24] = 1;
  work.KKT[26] = 1;
  work.KKT[28] = 1;
  work.KKT[30] = 1;
  work.KKT[32] = 1;
  work.KKT[34] = 1;
  work.KKT[36] = 1;
  work.KKT[38] = 1;
  work.KKT[40] = 1;
  work.KKT[42] = 1;
  work.KKT[44] = 1;
  work.KKT[46] = 1;
  work.KKT[48] = 1;
  work.KKT[50] = 1;
  work.KKT[52] = 1;
  work.KKT[54] = 1;
  work.KKT[56] = 1;
  work.KKT[58] = 1;
  work.KKT[60] = 1;
  work.KKT[62] = 1;
  work.KKT[63] = work.block_33[0];
  work.KKT[65] = work.block_33[0];
  work.KKT[74] = work.block_33[0];
  work.KKT[81] = work.block_33[0];
  work.KKT[88] = work.block_33[0];
  work.KKT[95] = work.block_33[0];
  work.KKT[102] = work.block_33[0];
  work.KKT[109] = work.block_33[0];
  work.KKT[116] = work.block_33[0];
  work.KKT[123] = work.block_33[0];
  work.KKT[130] = work.block_33[0];
  work.KKT[137] = work.block_33[0];
  work.KKT[147] = work.block_33[0];
  work.KKT[76] = work.block_33[0];
  work.KKT[83] = work.block_33[0];
  work.KKT[90] = work.block_33[0];
  work.KKT[97] = work.block_33[0];
  work.KKT[104] = work.block_33[0];
  work.KKT[111] = work.block_33[0];
  work.KKT[118] = work.block_33[0];
  work.KKT[125] = work.block_33[0];
  work.KKT[132] = work.block_33[0];
  work.KKT[139] = work.block_33[0];
  work.KKT[148] = work.block_33[0];
  work.KKT[67] = work.block_33[0];
  work.KKT[69] = work.block_33[0];
  work.KKT[64] = -1;
  work.KKT[66] = 1;
  work.KKT[75] = params.ones[0];
  work.KKT[82] = params.ones[1];
  work.KKT[89] = params.ones[2];
  work.KKT[96] = params.ones[3];
  work.KKT[103] = params.ones[4];
  work.KKT[110] = params.ones[5];
  work.KKT[117] = params.ones[6];
  work.KKT[124] = params.ones[7];
  work.KKT[131] = params.ones[8];
  work.KKT[138] = params.ones[9];
  work.KKT[142] = params.ones[10];
  work.KKT[72] = -1;
  work.KKT[79] = -1;
  work.KKT[86] = -1;
  work.KKT[93] = -1;
  work.KKT[100] = -1;
  work.KKT[107] = -1;
  work.KKT[114] = -1;
  work.KKT[121] = -1;
  work.KKT[128] = -1;
  work.KKT[135] = -1;
  work.KKT[145] = -1;
  work.KKT[77] = -params.ones[0];
  work.KKT[84] = -params.ones[1];
  work.KKT[91] = -params.ones[2];
  work.KKT[98] = -params.ones[3];
  work.KKT[105] = -params.ones[4];
  work.KKT[112] = -params.ones[5];
  work.KKT[119] = -params.ones[6];
  work.KKT[126] = -params.ones[7];
  work.KKT[133] = -params.ones[8];
  work.KKT[140] = -params.ones[9];
  work.KKT[143] = -params.ones[10];
  work.KKT[73] = 1;
  work.KKT[80] = 1;
  work.KKT[87] = 1;
  work.KKT[94] = 1;
  work.KKT[101] = 1;
  work.KKT[108] = 1;
  work.KKT[115] = 1;
  work.KKT[122] = 1;
  work.KKT[129] = 1;
  work.KKT[136] = 1;
  work.KKT[146] = 1;
  work.KKT[68] = -1;
  work.KKT[70] = 1;
}

//testsolver.c
/* Produced by CVXGEN, 2016-05-03 19:35:15 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

//Vars vars;
//Params params;
//Workspace work;
//Settings settings;
GensolverFirst::GensolverFirst( int numberOfCases, double cost2, double cost1, double cost0, double PgenMax, double PgenMin, double RgenMax, double RgenMin, double Beta, double Gamma, double PgenPrevious ) {
  params.caseNumber = numberOfCases;
  params.contNum[0] = 0.0;
  params.c2[0] = cost2;
  params.c1[0] = cost1;
  params.c0[0] = cost0;
  params.beta[0] = Beta;
  params.PgNu[0] = 0.0;
  params.gamma[0] = Gamma;
  params.rho[0] = 1.0;
  params.PgMin[0] = PgenMin;
  params.PgMax[0] = PgenMax;
  params.RgMin[0] = RgenMin;
  params.RgMax[0] = RgenMax;
  params.PgPrev[0] = PgenPrevious;
for ( int i = 0; i <= 10; ++i ) {
  params.ones[i] = 0.0;
  params.B[i] = 0.0;
  params.D[i] = 0.0;
  params.lambda_1[i] = 0.0;
  params.lambda_2[i] = 0.0;
  params.PgNextNu[i] = 0.0;
  params.Pg_N_init[i] = 0.0;
  params.Pg_N_avg[i] = 0.0;
  params.ug_N[i] = 0.0;
  params.Vg_N_avg[i] = 0.0;
  params.Thetag_N_avg[i] = 0.0;
  params.vg_N[i] = 0.0;
}
  global_seed = 1;
  Piterate = NULL;
  Thiterate = NULL;
  PgNextiterate = NULL;
  /*for ( int i = 0; i <= 3; ++i ) {
   Thiterate[i] = 0.0;
  }*/
  set_defaults();
  setup_indexing();
}
GensolverFirst::~GensolverFirst() {
}
#define NUMTESTS 0
void GensolverFirst::mainsolve(int APPItCount, double gsRho, double Pgprev[], double Pgavrg[], double Pprice[], double Apriceavg[], double Aavg[], double Aprice[], double PgenAPP, double PgenNextAPP[], double BAPP[], double DAPP[], double LambAPP1[], double LambAPP2[]) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data(APPItCount, gsRho, Pgprev, Pgavrg, Pprice, Apriceavg, Aavg, Aprice, PgenAPP, PgenNextAPP, BAPP, DAPP, LambAPP1, LambAPP2);
  /* Solve problem instance for the record. */
  settings.verbose = 0;
  num_iters = solve();
#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;
  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();
  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;
  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif
}

void GensolverFirst::load_default_data(int APPItCount, double gsRho, double Pgenprev[], double Pgenavg[], double Powerprice[], double Angpriceavg[], double Angavg[], double Angprice[], double PgenAPP, double PgenNextAPP[], double BAPP[], double DAPP[], double LambAPP1[], double LambAPP2[]) {
 if (APPItCount==1) {
  params.rho[0] = gsRho;
  params.PgNu[0] = 0.0;
 }
 else {
  params.rho[0] = gsRho;
  params.PgNu[0] = PgenAPP;
 }
for ( int i = 0; i <= params.caseNumber; ++i ) {
  params.ones[i] = 1.0;
  params.B[i] = BAPP[i];
  params.D[i] = DAPP[i];
  params.PgNextNu[i] = PgenNextAPP[i];
  params.Pg_N_init[i] = Pgenprev[i];
  params.Pg_N_avg[i] = Pgenavg[i];
  params.lambda_1[i] = LambAPP1[i];
  params.lambda_2[i] = LambAPP2[i];
  params.ug_N[i] = Powerprice[i];
  params.Vg_N_avg[i] = Angpriceavg[i];
  params.Thetag_N_avg[i] = Angavg[i];
  params.vg_N[i] = Angprice[i];
}
  set_defaults();
  setup_indexing();
}
double GensolverFirst::getPSol(void) {
 Piterate = (vars.Pg);
 return *Piterate;
}

double GensolverFirst::getPPrevSol(void) {
 return params.PgPrev[0];
}

double *GensolverFirst::getPNextSol(void) {
 PgNextiterate = (vars.PgNext);
 return PgNextiterate;
}

double GensolverFirst::getObj(void) {
 return (params.c2[0])*(*Piterate)*(*Piterate)+(params.c1[0])*(*Piterate)+(params.c0[0]);
}

double *GensolverFirst::getThetaPtr(void){
Thiterate = (vars.Thetag);
return Thiterate;
}

double GensolverFirst::getRMax() // start of getPMax function
{
	return params.RgMax[0];
} // end of getMax

double GensolverFirst::getRMin() // start of getPMin function
{
	return params.RgMin[0];
} // end of getPMin

double GensolverFirst::getPMax() // start of getPMax function
{
	return params.PgMax[0];
} // end of getMax

double GensolverFirst::getPMin() // start of getPMin function
{
	return params.PgMin[0];
} // end of getPMin

double GensolverFirst::getQuadCoeff() // start of getQuadCoeff function
{
	return params.c2[0]; 
} // end of getQuadCoeff

double GensolverFirst::getLinCoeff() // start of getLinCoeff function
{
	return params.c1[0]; 
} // end of getLinCoeff

double GensolverFirst::getConstCoeff() // start of getConstCoeff function
{
	return params.c0[0]; 
} // end of getConstCoeff
