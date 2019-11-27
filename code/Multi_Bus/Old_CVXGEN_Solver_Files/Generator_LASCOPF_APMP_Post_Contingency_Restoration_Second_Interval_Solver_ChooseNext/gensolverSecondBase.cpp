//solver.pp
/* Produced by CVXGEN, 2019-01-20 18:36:23 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */
#include "solver.h"
double eval_gap(void) {
  int i;
  double gap;
  gap = 0;
  for (i = 0; i < 6; i++)
    gap += work.z[i]*work.s[i];
  return gap;
}
void set_defaults(void) {
  settings.resid_tol = 1e-6;
  settings.eps = 1e-4;
  settings.max_iters = 25;
  settings.refine_steps = 1;
  settings.s_init = 1;
  settings.z_init = 1;
  settings.debug = 0;
  settings.verbose = 1;
  settings.verbose_refinement = 0;
  settings.better_start = 1;
  settings.kkt_reg = 1e-7;
}
void setup_pointers(void) {
  work.y = work.x + 3;
  work.s = work.x + 3;
  work.z = work.x + 9;
  vars.Pg = work.x + 0;
  vars.PgPrev = work.x + 1;
  vars.Thetag = work.x + 2;
}
void setup_indexing(void) {
  setup_pointers();
}
void set_start(void) {
  int i;
  for (i = 0; i < 3; i++)
    work.x[i] = 0;
  for (i = 0; i < 0; i++)
    work.y[i] = 0;
  for (i = 0; i < 6; i++)
    work.s[i] = (work.h[i] > 0) ? work.h[i] : settings.s_init;
  for (i = 0; i < 6; i++)
    work.z[i] = settings.z_init;
}
double eval_objv(void) {
  int i;
  double objv;
  /* Borrow space in work.rhs. */
  multbyP(work.rhs, work.x);
  objv = 0;
  for (i = 0; i < 3; i++)
    objv += work.x[i]*work.rhs[i];
  objv *= 0.5;
  for (i = 0; i < 3; i++)
    objv += work.q[i]*work.x[i];
  objv += params.c0[0]+work.quad_399913824256[0]+work.quad_608032620544[0]+work.quad_155015135232[0]+work.quad_497818112000[0];
  return objv;
}
void fillrhs_aff(void) {
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 3;
  r3 = work.rhs + 9;
  r4 = work.rhs + 15;
  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, work.y);
  multbymGT(work.buffer, work.z);
  for (i = 0; i < 3; i++)
    r1[i] += work.buffer[i];
  multbyP(work.buffer, work.x);
  for (i = 0; i < 3; i++)
    r1[i] -= work.buffer[i] + work.q[i];
  /* r2 = -z. */
  for (i = 0; i < 6; i++)
    r2[i] = -work.z[i];
  /* r3 = -Gx - s + h. */
  multbymG(r3, work.x);
  for (i = 0; i < 6; i++)
    r3[i] += -work.s[i] + work.h[i];
  /* r4 = -Ax + b. */
  multbymA(r4, work.x);
  for (i = 0; i < 0; i++)
    r4[i] += work.b[i];
}
void fillrhs_cc(void) {
  int i;
  double *r2;
  double *ds_aff, *dz_aff;
  double mu;
  double alpha;
  double sigma;
  double smu;
  double minval;
  r2 = work.rhs + 3;
  ds_aff = work.lhs_aff + 3;
  dz_aff = work.lhs_aff + 9;
  mu = 0;
  for (i = 0; i < 6; i++)
    mu += work.s[i]*work.z[i];
  /* Don't finish calculating mu quite yet. */
  /* Find min(min(ds./s), min(dz./z)). */
  minval = 0;
  for (i = 0; i < 6; i++)
    if (ds_aff[i] < minval*work.s[i])
      minval = ds_aff[i]/work.s[i];
  for (i = 0; i < 6; i++)
    if (dz_aff[i] < minval*work.z[i])
      minval = dz_aff[i]/work.z[i];
  /* Find alpha. */
  if (-1 < minval)
      alpha = 1;
  else
      alpha = -1/minval;
  sigma = 0;
  for (i = 0; i < 6; i++)
    sigma += (work.s[i] + alpha*ds_aff[i])*
      (work.z[i] + alpha*dz_aff[i]);
  sigma /= mu;
  sigma = sigma*sigma*sigma;
  /* Finish calculating mu now. */
  mu *= 0.16666666666666666;
  smu = sigma*mu;
  /* Fill-in the rhs. */
  for (i = 0; i < 3; i++)
    work.rhs[i] = 0;
  for (i = 9; i < 15; i++)
    work.rhs[i] = 0;
  for (i = 0; i < 6; i++)
    r2[i] = work.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}
void refine(double *target, double *var) {
  int i, j;
  double *residual = work.buffer;
  double norm2;
  double *new_var = work.buffer2;
  for (j = 0; j < settings.refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 15; i++) {
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
    for (i = 0; i < 15; i++) {
      var[i] -= new_var[i];
    }
  }
#ifndef ZERO_LIBRARY_MODE
  if (settings.verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 15; i++) {
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
double calc_ineq_resid_squared(void) {
  /* Calculates the norm ||-Gx - s + h||. */
  double norm2_squared;
  int i;
  /* Find -Gx. */
  multbymG(work.buffer, work.x);
  /* Add -s + h. */
  for (i = 0; i < 6; i++)
    work.buffer[i] += -work.s[i] + work.h[i];
  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 6; i++)
    norm2_squared += work.buffer[i]*work.buffer[i];
  return norm2_squared;
}
double calc_eq_resid_squared(void) {
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
void better_start(void) {
  /* Calculates a better starting point, using a similar approach to CVXOPT. */
  /* Not yet speed optimized. */
  int i;
  double *x, *s, *z, *y;
  double alpha;
  work.block_33[0] = -1;
  /* Make sure sinvz is 1 to make hijacked KKT system ok. */
  for (i = 0; i < 6; i++)
    work.s_inv_z[i] = 1;
  fill_KKT();
  ldl_factor();
  fillrhs_start();
  /* Borrow work.lhs_aff for the solution. */
  ldl_solve(work.rhs, work.lhs_aff);
  /* Don't do any refinement for now. Precision doesn't matter too much. */
  x = work.lhs_aff;
  s = work.lhs_aff + 3;
  z = work.lhs_aff + 9;
  y = work.lhs_aff + 15;
  /* Just set x and y as is. */
  for (i = 0; i < 3; i++)
    work.x[i] = x[i];
  for (i = 0; i < 0; i++)
    work.y[i] = y[i];
  /* Now complete the initialization. Start with s. */
  /* Must have alpha > max(z). */
  alpha = -1e99;
  for (i = 0; i < 6; i++)
    if (alpha < z[i])
      alpha = z[i];
  if (alpha < 0) {
    for (i = 0; i < 6; i++)
      work.s[i] = -z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 6; i++)
      work.s[i] = -z[i] + alpha;
  }
  /* Now initialize z. */
  /* Now must have alpha > max(-z). */
  alpha = -1e99;
  for (i = 0; i < 6; i++)
    if (alpha < -z[i])
      alpha = -z[i];
  if (alpha < 0) {
    for (i = 0; i < 6; i++)
      work.z[i] = z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 6; i++)
      work.z[i] = z[i] + alpha;
  }
}
void fillrhs_start(void) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;
  r1 = work.rhs;
  r2 = work.rhs + 3;
  r3 = work.rhs + 9;
  r4 = work.rhs + 15;
  for (i = 0; i < 3; i++)
    r1[i] = -work.q[i];
  for (i = 0; i < 6; i++)
    r2[i] = 0;
  for (i = 0; i < 6; i++)
    r3[i] = work.h[i];
  for (i = 0; i < 0; i++)
    r4[i] = work.b[i];
}
long solve(void) {
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
    for (i = 0; i < 6; i++) {
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
    for (i = 0; i < 15; i++)
      work.lhs_aff[i] += work.lhs_cc[i];
    /* Rename aff to reflect its new meaning. */
    dx = work.lhs_aff;
    ds = work.lhs_aff + 3;
    dz = work.lhs_aff + 9;
    dy = work.lhs_aff + 15;
    /* Find min(min(ds./s), min(dz./z)). */
    minval = 0;
    for (i = 0; i < 6; i++)
      if (ds[i] < minval*work.s[i])
        minval = ds[i]/work.s[i];
    for (i = 0; i < 6; i++)
      if (dz[i] < minval*work.z[i])
        minval = dz[i]/work.z[i];
    /* Find alpha. */
    if (-0.99 < minval)
      alpha = 1;
    else
      alpha = -0.99/minval;
    /* Update the primal and dual variables. */
    for (i = 0; i < 3; i++)
      work.x[i] += alpha*dx[i];
    for (i = 0; i < 6; i++)
      work.s[i] += alpha*ds[i];
    for (i = 0; i < 6; i++)
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

// ldl.cpp
/* Produced by CVXGEN, 2019-01-20 18:36:21 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: ldl.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
/* Be sure to place ldl_solve first, so storage schemes are defined by it. */
void ldl_solve(double *target, double *var) {
  int i;
  /* Find var = (L*diag(work.d)*L') \ target, then unpermute. */
  /* Answer goes into var. */
  /* Forward substitution. */
  /* Include permutation as we retrieve from target. Use v so we can unpermute */
  /* later. */
  work.v[0] = target[2];
  work.v[1] = target[3];
  work.v[2] = target[4];
  work.v[3] = target[5];
  work.v[4] = target[6];
  work.v[5] = target[7];
  work.v[6] = target[8];
  work.v[7] = target[9]-work.L[0]*work.v[1];
  work.v[8] = target[10]-work.L[1]*work.v[2];
  work.v[9] = target[11]-work.L[2]*work.v[3];
  work.v[10] = target[12]-work.L[3]*work.v[4];
  work.v[11] = target[0]-work.L[4]*work.v[7]-work.L[5]*work.v[8]-work.L[6]*work.v[9]-work.L[7]*work.v[10];
  work.v[12] = target[1];
  work.v[13] = target[13]-work.L[8]*work.v[5]-work.L[9]*work.v[11]-work.L[10]*work.v[12];
  work.v[14] = target[14]-work.L[11]*work.v[6]-work.L[12]*work.v[11]-work.L[13]*work.v[12]-work.L[14]*work.v[13];
  /* Diagonal scaling. Assume correctness of work.d_inv. */
  for (i = 0; i < 15; i++)
    work.v[i] *= work.d_inv[i];
  /* Back substitution */
  work.v[13] -= work.L[14]*work.v[14];
  work.v[12] -= work.L[10]*work.v[13]+work.L[13]*work.v[14];
  work.v[11] -= work.L[9]*work.v[13]+work.L[12]*work.v[14];
  work.v[10] -= work.L[7]*work.v[11];
  work.v[9] -= work.L[6]*work.v[11];
  work.v[8] -= work.L[5]*work.v[11];
  work.v[7] -= work.L[4]*work.v[11];
  work.v[6] -= work.L[11]*work.v[14];
  work.v[5] -= work.L[8]*work.v[13];
  work.v[4] -= work.L[3]*work.v[10];
  work.v[3] -= work.L[2]*work.v[9];
  work.v[2] -= work.L[1]*work.v[8];
  work.v[1] -= work.L[0]*work.v[7];
  /* Unpermute the result, from v to var. */
  var[0] = work.v[11];
  var[1] = work.v[12];
  var[2] = work.v[0];
  var[3] = work.v[1];
  var[4] = work.v[2];
  var[5] = work.v[3];
  var[6] = work.v[4];
  var[7] = work.v[5];
  var[8] = work.v[6];
  var[9] = work.v[7];
  var[10] = work.v[8];
  var[11] = work.v[9];
  var[12] = work.v[10];
  var[13] = work.v[13];
  var[14] = work.v[14];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared norm for solution is %.8g.\n", check_residual(target, var));
  }
#endif
}
void ldl_factor(void) {
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
  work.L[0] = (work.KKT[2])*work.d_inv[1];
  work.v[2] = work.KKT[3];
  work.d[2] = work.v[2];
  if (work.d[2] < 0)
    work.d[2] = settings.kkt_reg;
  else
    work.d[2] += settings.kkt_reg;
  work.d_inv[2] = 1/work.d[2];
  work.L[1] = (work.KKT[4])*work.d_inv[2];
  work.v[3] = work.KKT[5];
  work.d[3] = work.v[3];
  if (work.d[3] < 0)
    work.d[3] = settings.kkt_reg;
  else
    work.d[3] += settings.kkt_reg;
  work.d_inv[3] = 1/work.d[3];
  work.L[2] = (work.KKT[6])*work.d_inv[3];
  work.v[4] = work.KKT[7];
  work.d[4] = work.v[4];
  if (work.d[4] < 0)
    work.d[4] = settings.kkt_reg;
  else
    work.d[4] += settings.kkt_reg;
  work.d_inv[4] = 1/work.d[4];
  work.L[3] = (work.KKT[8])*work.d_inv[4];
  work.v[5] = work.KKT[9];
  work.d[5] = work.v[5];
  if (work.d[5] < 0)
    work.d[5] = settings.kkt_reg;
  else
    work.d[5] += settings.kkt_reg;
  work.d_inv[5] = 1/work.d[5];
  work.L[8] = (work.KKT[10])*work.d_inv[5];
  work.v[6] = work.KKT[11];
  work.d[6] = work.v[6];
  if (work.d[6] < 0)
    work.d[6] = settings.kkt_reg;
  else
    work.d[6] += settings.kkt_reg;
  work.d_inv[6] = 1/work.d[6];
  work.L[11] = (work.KKT[12])*work.d_inv[6];
  work.v[1] = work.L[0]*work.d[1];
  work.v[7] = work.KKT[13]-work.L[0]*work.v[1];
  work.d[7] = work.v[7];
  if (work.d[7] > 0)
    work.d[7] = -settings.kkt_reg;
  else
    work.d[7] -= settings.kkt_reg;
  work.d_inv[7] = 1/work.d[7];
  work.L[4] = (work.KKT[14])*work.d_inv[7];
  work.v[2] = work.L[1]*work.d[2];
  work.v[8] = work.KKT[15]-work.L[1]*work.v[2];
  work.d[8] = work.v[8];
  if (work.d[8] > 0)
    work.d[8] = -settings.kkt_reg;
  else
    work.d[8] -= settings.kkt_reg;
  work.d_inv[8] = 1/work.d[8];
  work.L[5] = (work.KKT[16])*work.d_inv[8];
  work.v[3] = work.L[2]*work.d[3];
  work.v[9] = work.KKT[17]-work.L[2]*work.v[3];
  work.d[9] = work.v[9];
  if (work.d[9] > 0)
    work.d[9] = -settings.kkt_reg;
  else
    work.d[9] -= settings.kkt_reg;
  work.d_inv[9] = 1/work.d[9];
  work.L[6] = (work.KKT[18])*work.d_inv[9];
  work.v[4] = work.L[3]*work.d[4];
  work.v[10] = work.KKT[19]-work.L[3]*work.v[4];
  work.d[10] = work.v[10];
  if (work.d[10] > 0)
    work.d[10] = -settings.kkt_reg;
  else
    work.d[10] -= settings.kkt_reg;
  work.d_inv[10] = 1/work.d[10];
  work.L[7] = (work.KKT[20])*work.d_inv[10];
  work.v[7] = work.L[4]*work.d[7];
  work.v[8] = work.L[5]*work.d[8];
  work.v[9] = work.L[6]*work.d[9];
  work.v[10] = work.L[7]*work.d[10];
  work.v[11] = work.KKT[21]-work.L[4]*work.v[7]-work.L[5]*work.v[8]-work.L[6]*work.v[9]-work.L[7]*work.v[10];
  work.d[11] = work.v[11];
  if (work.d[11] < 0)
    work.d[11] = settings.kkt_reg;
  else
    work.d[11] += settings.kkt_reg;
  work.d_inv[11] = 1/work.d[11];
  work.L[9] = (work.KKT[22])*work.d_inv[11];
  work.L[12] = (work.KKT[23])*work.d_inv[11];
  work.v[12] = work.KKT[24];
  work.d[12] = work.v[12];
  if (work.d[12] < 0)
    work.d[12] = settings.kkt_reg;
  else
    work.d[12] += settings.kkt_reg;
  work.d_inv[12] = 1/work.d[12];
  work.L[10] = (work.KKT[25])*work.d_inv[12];
  work.L[13] = (work.KKT[26])*work.d_inv[12];
  work.v[5] = work.L[8]*work.d[5];
  work.v[11] = work.L[9]*work.d[11];
  work.v[12] = work.L[10]*work.d[12];
  work.v[13] = work.KKT[27]-work.L[8]*work.v[5]-work.L[9]*work.v[11]-work.L[10]*work.v[12];
  work.d[13] = work.v[13];
  if (work.d[13] > 0)
    work.d[13] = -settings.kkt_reg;
  else
    work.d[13] -= settings.kkt_reg;
  work.d_inv[13] = 1/work.d[13];
  work.L[14] = (-work.L[12]*work.v[11]-work.L[13]*work.v[12])*work.d_inv[13];
  work.v[6] = work.L[11]*work.d[6];
  work.v[11] = work.L[12]*work.d[11];
  work.v[12] = work.L[13]*work.d[12];
  work.v[13] = work.L[14]*work.d[13];
  work.v[14] = work.KKT[28]-work.L[11]*work.v[6]-work.L[12]*work.v[11]-work.L[13]*work.v[12]-work.L[14]*work.v[13];
  work.d[14] = work.v[14];
  if (work.d[14] > 0)
    work.d[14] = -settings.kkt_reg;
  else
    work.d[14] -= settings.kkt_reg;
  work.d_inv[14] = 1/work.d[14];
#ifndef ZERO_LIBRARY_MODE
  if (settings.debug) {
    printf("Squared Frobenius for factorization is %.8g.\n", check_factorization());
  }
#endif
}
double check_factorization(void) {
  /* Returns the squared Frobenius norm of A - L*D*L'. */
  double temp, residual;
  /* Only check the lower triangle. */
  residual = 0;
  temp = work.KKT[21]-1*work.d[11]*1-work.L[4]*work.d[7]*work.L[4]-work.L[5]*work.d[8]*work.L[5]-work.L[6]*work.d[9]*work.L[6]-work.L[7]*work.d[10]*work.L[7];
  residual += temp*temp;
  temp = work.KKT[24]-1*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[0]-1*work.d[0]*1;
  residual += temp*temp;
  temp = work.KKT[1]-1*work.d[1]*1;
  residual += temp*temp;
  temp = work.KKT[3]-1*work.d[2]*1;
  residual += temp*temp;
  temp = work.KKT[5]-1*work.d[3]*1;
  residual += temp*temp;
  temp = work.KKT[7]-1*work.d[4]*1;
  residual += temp*temp;
  temp = work.KKT[9]-1*work.d[5]*1;
  residual += temp*temp;
  temp = work.KKT[11]-1*work.d[6]*1;
  residual += temp*temp;
  temp = work.KKT[2]-work.L[0]*work.d[1]*1;
  residual += temp*temp;
  temp = work.KKT[4]-work.L[1]*work.d[2]*1;
  residual += temp*temp;
  temp = work.KKT[6]-work.L[2]*work.d[3]*1;
  residual += temp*temp;
  temp = work.KKT[8]-work.L[3]*work.d[4]*1;
  residual += temp*temp;
  temp = work.KKT[10]-work.L[8]*work.d[5]*1;
  residual += temp*temp;
  temp = work.KKT[12]-work.L[11]*work.d[6]*1;
  residual += temp*temp;
  temp = work.KKT[13]-work.L[0]*work.d[1]*work.L[0]-1*work.d[7]*1;
  residual += temp*temp;
  temp = work.KKT[15]-work.L[1]*work.d[2]*work.L[1]-1*work.d[8]*1;
  residual += temp*temp;
  temp = work.KKT[17]-work.L[2]*work.d[3]*work.L[2]-1*work.d[9]*1;
  residual += temp*temp;
  temp = work.KKT[19]-work.L[3]*work.d[4]*work.L[3]-1*work.d[10]*1;
  residual += temp*temp;
  temp = work.KKT[27]-work.L[8]*work.d[5]*work.L[8]-1*work.d[13]*1-work.L[9]*work.d[11]*work.L[9]-work.L[10]*work.d[12]*work.L[10];
  residual += temp*temp;
  temp = work.KKT[28]-work.L[11]*work.d[6]*work.L[11]-1*work.d[14]*1-work.L[12]*work.d[11]*work.L[12]-work.L[13]*work.d[12]*work.L[13]-work.L[14]*work.d[13]*work.L[14];
  residual += temp*temp;
  temp = work.KKT[14]-1*work.d[7]*work.L[4];
  residual += temp*temp;
  temp = work.KKT[16]-1*work.d[8]*work.L[5];
  residual += temp*temp;
  temp = work.KKT[18]-1*work.d[9]*work.L[6];
  residual += temp*temp;
  temp = work.KKT[20]-1*work.d[10]*work.L[7];
  residual += temp*temp;
  temp = work.KKT[22]-work.L[9]*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[25]-work.L[10]*work.d[12]*1;
  residual += temp*temp;
  temp = work.KKT[23]-work.L[12]*work.d[11]*1;
  residual += temp*temp;
  temp = work.KKT[26]-work.L[13]*work.d[12]*1;
  residual += temp*temp;
  return residual;
}
void matrix_multiply(double *result, double *source) {
  /* Finds result = A*source. */
  result[0] = work.KKT[21]*source[0]+work.KKT[14]*source[9]+work.KKT[16]*source[10]+work.KKT[18]*source[11]+work.KKT[20]*source[12]+work.KKT[22]*source[13]+work.KKT[23]*source[14];
  result[1] = work.KKT[24]*source[1]+work.KKT[25]*source[13]+work.KKT[26]*source[14];
  result[2] = work.KKT[0]*source[2];
  result[3] = work.KKT[1]*source[3]+work.KKT[2]*source[9];
  result[4] = work.KKT[3]*source[4]+work.KKT[4]*source[10];
  result[5] = work.KKT[5]*source[5]+work.KKT[6]*source[11];
  result[6] = work.KKT[7]*source[6]+work.KKT[8]*source[12];
  result[7] = work.KKT[9]*source[7]+work.KKT[10]*source[13];
  result[8] = work.KKT[11]*source[8]+work.KKT[12]*source[14];
  result[9] = work.KKT[2]*source[3]+work.KKT[13]*source[9]+work.KKT[14]*source[0];
  result[10] = work.KKT[4]*source[4]+work.KKT[15]*source[10]+work.KKT[16]*source[0];
  result[11] = work.KKT[6]*source[5]+work.KKT[17]*source[11]+work.KKT[18]*source[0];
  result[12] = work.KKT[8]*source[6]+work.KKT[19]*source[12]+work.KKT[20]*source[0];
  result[13] = work.KKT[10]*source[7]+work.KKT[27]*source[13]+work.KKT[22]*source[0]+work.KKT[25]*source[1];
  result[14] = work.KKT[12]*source[8]+work.KKT[28]*source[14]+work.KKT[23]*source[0]+work.KKT[26]*source[1];
}
double check_residual(double *target, double *multiplicand) {
  /* Returns the squared 2-norm of lhs - A*rhs. */
  /* Reuses v to find the residual. */
  int i;
  double residual;
  residual = 0;
  matrix_multiply(work.v, multiplicand);
  for (i = 0; i < 3; i++) {
    residual += (target[i] - work.v[i])*(target[i] - work.v[i]);
  }
  return residual;
}
void fill_KKT(void) {
  work.KKT[21] = 2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848);
  work.KKT[24] = 2*work.frac_479733190656;
  work.KKT[0] = 2*work.frac_121674190848;
  work.KKT[1] = work.s_inv_z[0];
  work.KKT[3] = work.s_inv_z[1];
  work.KKT[5] = work.s_inv_z[2];
  work.KKT[7] = work.s_inv_z[3];
  work.KKT[9] = work.s_inv_z[4];
  work.KKT[11] = work.s_inv_z[5];
  work.KKT[2] = 1;
  work.KKT[4] = 1;
  work.KKT[6] = 1;
  work.KKT[8] = 1;
  work.KKT[10] = 1;
  work.KKT[12] = 1;
  work.KKT[13] = work.block_33[0];
  work.KKT[15] = work.block_33[0];
  work.KKT[17] = work.block_33[0];
  work.KKT[19] = work.block_33[0];
  work.KKT[27] = work.block_33[0];
  work.KKT[28] = work.block_33[0];
  work.KKT[14] = -1;
  work.KKT[16] = 1;
  work.KKT[18] = params.selectZero[0];
  work.KKT[20] = -params.selectZero[0];
  work.KKT[22] = -1;
  work.KKT[25] = 1;
  work.KKT[23] = 1;
  work.KKT[26] = -1;
}

// matrix_support.cpp
/* Produced by CVXGEN, 2019-01-20 18:36:21 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
}
void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1);
  lhs[1] = -rhs[0]*(1);
  lhs[2] = -rhs[0]*(params.selectZero[0]);
  lhs[3] = -rhs[0]*(-params.selectZero[0]);
  lhs[4] = -rhs[0]*(-1)-rhs[1]*(1);
  lhs[5] = -rhs[0]*(1)-rhs[1]*(-1);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(-1)-rhs[1]*(1)-rhs[2]*(params.selectZero[0])-rhs[3]*(-params.selectZero[0])-rhs[4]*(-1)-rhs[5]*(1);
  lhs[1] = -rhs[4]*(1)-rhs[5]*(-1);
  lhs[2] = 0;
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*(params.c2[0]+work.frac_479733190656+work.frac_121674190848));
  lhs[1] = rhs[1]*(2*work.frac_479733190656);
  lhs[2] = rhs[2]*(2*work.frac_121674190848);
}
void fillq(void) {
  work.q[0] = params.c1[0]-2*params.PgNu[0]*work.frac_479733190656+params.gammaSC[0]*(params.BSC[0]+params.BSC[1]+params.BSC[2]+params.BSC[3]+params.BSC[4]+params.BSC[5]+params.BSC[6]+params.BSC[7]+params.BSC[8]+params.BSC[9]+params.BSC[10]+params.BSC[11]+params.BSC[12]+params.BSC[13]+params.BSC[14]+params.BSC[15]+params.BSC[16]+params.BSC[17]+params.BSC[18]+params.BSC[19]+params.BSC[20]+params.BSC[21]+params.BSC[22]+params.BSC[23]+params.BSC[24]+params.BSC[25]+params.BSC[26]+params.BSC[27]+params.BSC[28]+params.BSC[29]+params.BSC[30]+params.BSC[31]+params.BSC[32]+params.BSC[33]+params.BSC[34]+params.BSC[35]+params.BSC[36]+params.BSC[37]+params.BSC[38]+params.BSC[39]+params.BSC[40]+params.BSC[41]+params.BSC[42]+params.BSC[43]+params.BSC[44]+params.BSC[45]+params.BSC[46]+params.BSC[47]+params.BSC[48]+params.BSC[49]+params.BSC[50]+params.BSC[51]+params.BSC[52]+params.BSC[53]+params.BSC[54]+params.BSC[55]+params.BSC[56]+params.BSC[57]+params.BSC[58]+params.BSC[59]+params.BSC[60]+params.BSC[61]+params.BSC[62]+params.BSC[63]+params.BSC[64]+params.BSC[65]+params.BSC[66]+params.BSC[67]+params.BSC[68]+params.BSC[69]+params.BSC[70]+params.BSC[71]+params.BSC[72]+params.BSC[73]+params.BSC[74]+params.BSC[75]+params.BSC[76]+params.BSC[77]+params.BSC[78]+params.BSC[79]+params.BSC[80]+params.BSC[81]+params.BSC[82]+params.BSC[83]+params.BSC[84]+params.BSC[85]+params.BSC[86]+params.BSC[87]+params.BSC[88]+params.BSC[89]+params.BSC[90]+params.BSC[91]+params.BSC[92]+params.BSC[93]+params.BSC[94]+params.BSC[95]+params.BSC[96]+params.BSC[97]+params.BSC[98]+params.BSC[99]+params.BSC[100]+params.BSC[101]+params.BSC[102]+params.BSC[103]+params.BSC[104]+params.BSC[105]+params.BSC[106]+params.BSC[107]+params.BSC[108]+params.BSC[109]+params.BSC[110]+params.BSC[111]+params.BSC[112]+params.BSC[113]+params.BSC[114]+params.BSC[115]+params.BSC[116]+params.BSC[117]+params.BSC[118]+params.BSC[119]+params.BSC[120]+params.BSC[121]+params.BSC[122]+params.BSC[123]+params.BSC[124]+params.BSC[125]+params.BSC[126]+params.BSC[127]+params.BSC[128]+params.BSC[129]+params.BSC[130]+params.BSC[131]+params.BSC[132]+params.BSC[133]+params.BSC[134]+params.BSC[135]+params.BSC[136]+params.BSC[137]+params.BSC[138]+params.BSC[139]+params.BSC[140]+params.BSC[141]+params.BSC[142]+params.BSC[143]+params.BSC[144]+params.BSC[145]+params.BSC[146]+params.BSC[147]+params.BSC[148]+params.BSC[149]+params.BSC[150]+params.BSC[151]+params.BSC[152]+params.BSC[153]+params.BSC[154]+params.BSC[155]+params.BSC[156]+params.BSC[157]+params.BSC[158]+params.BSC[159]+params.BSC[160]+params.BSC[161]+params.BSC[162]+params.BSC[163]+params.BSC[164]+params.BSC[165]+params.BSC[166]+params.BSC[167]+params.BSC[168]+params.BSC[169]+params.BSC[170]+params.BSC[171]+params.BSC[172]+params.BSC[173]+params.BSC[174]+params.BSC[175]+params.BSC[176]+params.BSC[177]+params.BSC[178]+params.BSC[179]+params.BSC[180]+params.BSC[181]+params.BSC[182]+params.BSC[183]+params.BSC[184]+params.BSC[185]+params.BSC[186]+params.BSC[187]+params.BSC[188]+params.BSC[189]+params.BSC[190]+params.BSC[191]+params.BSC[192]+params.BSC[193]+params.BSC[194]+params.BSC[195]+params.BSC[196]+params.BSC[197]+params.BSC[198]+params.BSC[199]+params.BSC[200]+params.BSC[201]+params.BSC[202]+params.BSC[203]+params.BSC[204]+params.BSC[205]+params.BSC[206]+params.BSC[207]+params.BSC[208]+params.BSC[209]+params.BSC[210]+params.BSC[211]+params.BSC[212]+params.BSC[213]+params.BSC[214]+params.BSC[215]+params.BSC[216]+params.BSC[217]+params.BSC[218]+params.BSC[219]+params.BSC[220]+params.BSC[221]+params.BSC[222]+params.BSC[223]+params.BSC[224]+params.BSC[225]+params.BSC[226]+params.BSC[227]+params.BSC[228]+params.BSC[229]+params.BSC[230]+params.BSC[231]+params.BSC[232]+params.BSC[233]+params.BSC[234]+params.BSC[235]+params.BSC[236]+params.BSC[237]+params.BSC[238]+params.BSC[239]+params.BSC[240]+params.BSC[241]+params.BSC[242]+params.BSC[243]+params.BSC[244]+params.BSC[245]+params.BSC[246]+params.BSC[247]+params.BSC[248]+params.BSC[249]+params.BSC[250]+params.BSC[251]+params.BSC[252]+params.BSC[253]+params.BSC[254]+params.BSC[255]+params.BSC[256]+params.BSC[257]+params.BSC[258]+params.BSC[259]+params.BSC[260]+params.BSC[261]+params.BSC[262]+params.BSC[263]+params.BSC[264]+params.BSC[265]+params.BSC[266]+params.BSC[267]+params.BSC[268]+params.BSC[269]+params.BSC[270]+params.BSC[271]+params.BSC[272]+params.BSC[273]+params.BSC[274]+params.BSC[275]+params.BSC[276]+params.BSC[277]+params.BSC[278]+params.BSC[279]+params.BSC[280]+params.BSC[281]+params.BSC[282]+params.BSC[283]+params.BSC[284]+params.BSC[285]+params.BSC[286]+params.BSC[287]+params.BSC[288]+params.BSC[289]+params.BSC[290]+params.BSC[291]+params.BSC[292]+params.BSC[293]+params.BSC[294]+params.BSC[295]+params.BSC[296]+params.BSC[297]+params.BSC[298]+params.BSC[299]+params.BSC[300]+params.BSC[301]+params.BSC[302]+params.BSC[303]+params.BSC[304]+params.BSC[305]+params.BSC[306]+params.BSC[307]+params.BSC[308]+params.BSC[309]+params.BSC[310]+params.BSC[311]+params.BSC[312]+params.BSC[313]+params.BSC[314]+params.BSC[315]+params.BSC[316]+params.BSC[317]+params.BSC[318]+params.BSC[319]+params.BSC[320]+params.BSC[321]+params.BSC[322]+params.BSC[323]+params.BSC[324]+params.BSC[325]+params.BSC[326]+params.BSC[327]+params.BSC[328]+params.BSC[329]+params.BSC[330]+params.BSC[331]+params.BSC[332]+params.BSC[333]+params.BSC[334]+params.BSC[335]+params.BSC[336]+params.BSC[337]+params.BSC[338]+params.BSC[339]+params.BSC[340]+params.BSC[341]+params.BSC[342]+params.BSC[343]+params.BSC[344]+params.BSC[345]+params.BSC[346]+params.BSC[347]+params.BSC[348]+params.BSC[349]+params.BSC[350]+params.BSC[351]+params.BSC[352]+params.BSC[353]+params.BSC[354]+params.BSC[355]+params.BSC[356]+params.BSC[357]+params.BSC[358]+params.BSC[359]+params.BSC[360]+params.BSC[361]+params.BSC[362]+params.BSC[363]+params.BSC[364]+params.BSC[365]+params.BSC[366]+params.BSC[367]+params.BSC[368]+params.BSC[369]+params.BSC[370]+params.BSC[371]+params.BSC[372]+params.BSC[373]+params.BSC[374]+params.BSC[375]+params.BSC[376]+params.BSC[377]+params.BSC[378]+params.BSC[379]+params.BSC[380]+params.BSC[381]+params.BSC[382]+params.BSC[383]+params.BSC[384]+params.BSC[385]+params.BSC[386]+params.BSC[387]+params.BSC[388]+params.BSC[389]+params.BSC[390]+params.BSC[391]+params.BSC[392]+params.BSC[393]+params.BSC[394]+params.BSC[395]+params.BSC[396]+params.BSC[397]+params.BSC[398]+params.BSC[399]+params.BSC[400]+params.BSC[401]+params.BSC[402]+params.BSC[403]+params.BSC[404]+params.BSC[405]+params.BSC[406]+params.BSC[407]+params.BSC[408]+params.BSC[409]+params.BSC[410]+params.BSC[411]+params.BSC[412]+params.BSC[413]+params.BSC[414]+params.BSC[415]+params.BSC[416]+params.BSC[417]+params.BSC[418]+params.BSC[419]+params.BSC[420]+params.BSC[421]+params.BSC[422]+params.BSC[423]+params.BSC[424]+params.BSC[425]+params.BSC[426]+params.BSC[427]+params.BSC[428]+params.BSC[429]+params.BSC[430]+params.BSC[431]+params.BSC[432]+params.BSC[433]+params.BSC[434]+params.BSC[435]+params.BSC[436]+params.BSC[437]+params.BSC[438]+params.BSC[439]+params.BSC[440]+params.BSC[441]+params.BSC[442]+params.BSC[443]+params.BSC[444]+params.BSC[445]+params.BSC[446]+params.BSC[447]+params.BSC[448]+params.BSC[449]+params.BSC[450]+params.BSC[451]+params.BSC[452]+params.BSC[453]+params.BSC[454]+params.BSC[455]+params.BSC[456]+params.BSC[457]+params.BSC[458]+params.BSC[459]+params.BSC[460]+params.BSC[461]+params.BSC[462]+params.BSC[463]+params.BSC[464]+params.BSC[465]+params.BSC[466]+params.BSC[467]+params.BSC[468]+params.BSC[469]+params.BSC[470]+params.BSC[471]+params.BSC[472]+params.BSC[473]+params.BSC[474]+params.BSC[475]+params.BSC[476]+params.BSC[477]+params.BSC[478]+params.BSC[479]+params.BSC[480]+params.BSC[481]+params.BSC[482]+params.BSC[483]+params.BSC[484]+params.BSC[485]+params.BSC[486]+params.BSC[487]+params.BSC[488]+params.BSC[489]+params.BSC[490]+params.BSC[491]+params.BSC[492]+params.BSC[493]+params.BSC[494]+params.BSC[495]+params.BSC[496]+params.BSC[497]+params.BSC[498]+params.BSC[499])+params.lambda_1SC[0]+params.lambda_1SC[1]+params.lambda_1SC[2]+params.lambda_1SC[3]+params.lambda_1SC[4]+params.lambda_1SC[5]+params.lambda_1SC[6]+params.lambda_1SC[7]+params.lambda_1SC[8]+params.lambda_1SC[9]+params.lambda_1SC[10]+params.lambda_1SC[11]+params.lambda_1SC[12]+params.lambda_1SC[13]+params.lambda_1SC[14]+params.lambda_1SC[15]+params.lambda_1SC[16]+params.lambda_1SC[17]+params.lambda_1SC[18]+params.lambda_1SC[19]+params.lambda_1SC[20]+params.lambda_1SC[21]+params.lambda_1SC[22]+params.lambda_1SC[23]+params.lambda_1SC[24]+params.lambda_1SC[25]+params.lambda_1SC[26]+params.lambda_1SC[27]+params.lambda_1SC[28]+params.lambda_1SC[29]+params.lambda_1SC[30]+params.lambda_1SC[31]+params.lambda_1SC[32]+params.lambda_1SC[33]+params.lambda_1SC[34]+params.lambda_1SC[35]+params.lambda_1SC[36]+params.lambda_1SC[37]+params.lambda_1SC[38]+params.lambda_1SC[39]+params.lambda_1SC[40]+params.lambda_1SC[41]+params.lambda_1SC[42]+params.lambda_1SC[43]+params.lambda_1SC[44]+params.lambda_1SC[45]+params.lambda_1SC[46]+params.lambda_1SC[47]+params.lambda_1SC[48]+params.lambda_1SC[49]+params.lambda_1SC[50]+params.lambda_1SC[51]+params.lambda_1SC[52]+params.lambda_1SC[53]+params.lambda_1SC[54]+params.lambda_1SC[55]+params.lambda_1SC[56]+params.lambda_1SC[57]+params.lambda_1SC[58]+params.lambda_1SC[59]+params.lambda_1SC[60]+params.lambda_1SC[61]+params.lambda_1SC[62]+params.lambda_1SC[63]+params.lambda_1SC[64]+params.lambda_1SC[65]+params.lambda_1SC[66]+params.lambda_1SC[67]+params.lambda_1SC[68]+params.lambda_1SC[69]+params.lambda_1SC[70]+params.lambda_1SC[71]+params.lambda_1SC[72]+params.lambda_1SC[73]+params.lambda_1SC[74]+params.lambda_1SC[75]+params.lambda_1SC[76]+params.lambda_1SC[77]+params.lambda_1SC[78]+params.lambda_1SC[79]+params.lambda_1SC[80]+params.lambda_1SC[81]+params.lambda_1SC[82]+params.lambda_1SC[83]+params.lambda_1SC[84]+params.lambda_1SC[85]+params.lambda_1SC[86]+params.lambda_1SC[87]+params.lambda_1SC[88]+params.lambda_1SC[89]+params.lambda_1SC[90]+params.lambda_1SC[91]+params.lambda_1SC[92]+params.lambda_1SC[93]+params.lambda_1SC[94]+params.lambda_1SC[95]+params.lambda_1SC[96]+params.lambda_1SC[97]+params.lambda_1SC[98]+params.lambda_1SC[99]+params.lambda_1SC[100]+params.lambda_1SC[101]+params.lambda_1SC[102]+params.lambda_1SC[103]+params.lambda_1SC[104]+params.lambda_1SC[105]+params.lambda_1SC[106]+params.lambda_1SC[107]+params.lambda_1SC[108]+params.lambda_1SC[109]+params.lambda_1SC[110]+params.lambda_1SC[111]+params.lambda_1SC[112]+params.lambda_1SC[113]+params.lambda_1SC[114]+params.lambda_1SC[115]+params.lambda_1SC[116]+params.lambda_1SC[117]+params.lambda_1SC[118]+params.lambda_1SC[119]+params.lambda_1SC[120]+params.lambda_1SC[121]+params.lambda_1SC[122]+params.lambda_1SC[123]+params.lambda_1SC[124]+params.lambda_1SC[125]+params.lambda_1SC[126]+params.lambda_1SC[127]+params.lambda_1SC[128]+params.lambda_1SC[129]+params.lambda_1SC[130]+params.lambda_1SC[131]+params.lambda_1SC[132]+params.lambda_1SC[133]+params.lambda_1SC[134]+params.lambda_1SC[135]+params.lambda_1SC[136]+params.lambda_1SC[137]+params.lambda_1SC[138]+params.lambda_1SC[139]+params.lambda_1SC[140]+params.lambda_1SC[141]+params.lambda_1SC[142]+params.lambda_1SC[143]+params.lambda_1SC[144]+params.lambda_1SC[145]+params.lambda_1SC[146]+params.lambda_1SC[147]+params.lambda_1SC[148]+params.lambda_1SC[149]+params.lambda_1SC[150]+params.lambda_1SC[151]+params.lambda_1SC[152]+params.lambda_1SC[153]+params.lambda_1SC[154]+params.lambda_1SC[155]+params.lambda_1SC[156]+params.lambda_1SC[157]+params.lambda_1SC[158]+params.lambda_1SC[159]+params.lambda_1SC[160]+params.lambda_1SC[161]+params.lambda_1SC[162]+params.lambda_1SC[163]+params.lambda_1SC[164]+params.lambda_1SC[165]+params.lambda_1SC[166]+params.lambda_1SC[167]+params.lambda_1SC[168]+params.lambda_1SC[169]+params.lambda_1SC[170]+params.lambda_1SC[171]+params.lambda_1SC[172]+params.lambda_1SC[173]+params.lambda_1SC[174]+params.lambda_1SC[175]+params.lambda_1SC[176]+params.lambda_1SC[177]+params.lambda_1SC[178]+params.lambda_1SC[179]+params.lambda_1SC[180]+params.lambda_1SC[181]+params.lambda_1SC[182]+params.lambda_1SC[183]+params.lambda_1SC[184]+params.lambda_1SC[185]+params.lambda_1SC[186]+params.lambda_1SC[187]+params.lambda_1SC[188]+params.lambda_1SC[189]+params.lambda_1SC[190]+params.lambda_1SC[191]+params.lambda_1SC[192]+params.lambda_1SC[193]+params.lambda_1SC[194]+params.lambda_1SC[195]+params.lambda_1SC[196]+params.lambda_1SC[197]+params.lambda_1SC[198]+params.lambda_1SC[199]+params.lambda_1SC[200]+params.lambda_1SC[201]+params.lambda_1SC[202]+params.lambda_1SC[203]+params.lambda_1SC[204]+params.lambda_1SC[205]+params.lambda_1SC[206]+params.lambda_1SC[207]+params.lambda_1SC[208]+params.lambda_1SC[209]+params.lambda_1SC[210]+params.lambda_1SC[211]+params.lambda_1SC[212]+params.lambda_1SC[213]+params.lambda_1SC[214]+params.lambda_1SC[215]+params.lambda_1SC[216]+params.lambda_1SC[217]+params.lambda_1SC[218]+params.lambda_1SC[219]+params.lambda_1SC[220]+params.lambda_1SC[221]+params.lambda_1SC[222]+params.lambda_1SC[223]+params.lambda_1SC[224]+params.lambda_1SC[225]+params.lambda_1SC[226]+params.lambda_1SC[227]+params.lambda_1SC[228]+params.lambda_1SC[229]+params.lambda_1SC[230]+params.lambda_1SC[231]+params.lambda_1SC[232]+params.lambda_1SC[233]+params.lambda_1SC[234]+params.lambda_1SC[235]+params.lambda_1SC[236]+params.lambda_1SC[237]+params.lambda_1SC[238]+params.lambda_1SC[239]+params.lambda_1SC[240]+params.lambda_1SC[241]+params.lambda_1SC[242]+params.lambda_1SC[243]+params.lambda_1SC[244]+params.lambda_1SC[245]+params.lambda_1SC[246]+params.lambda_1SC[247]+params.lambda_1SC[248]+params.lambda_1SC[249]+params.lambda_1SC[250]+params.lambda_1SC[251]+params.lambda_1SC[252]+params.lambda_1SC[253]+params.lambda_1SC[254]+params.lambda_1SC[255]+params.lambda_1SC[256]+params.lambda_1SC[257]+params.lambda_1SC[258]+params.lambda_1SC[259]+params.lambda_1SC[260]+params.lambda_1SC[261]+params.lambda_1SC[262]+params.lambda_1SC[263]+params.lambda_1SC[264]+params.lambda_1SC[265]+params.lambda_1SC[266]+params.lambda_1SC[267]+params.lambda_1SC[268]+params.lambda_1SC[269]+params.lambda_1SC[270]+params.lambda_1SC[271]+params.lambda_1SC[272]+params.lambda_1SC[273]+params.lambda_1SC[274]+params.lambda_1SC[275]+params.lambda_1SC[276]+params.lambda_1SC[277]+params.lambda_1SC[278]+params.lambda_1SC[279]+params.lambda_1SC[280]+params.lambda_1SC[281]+params.lambda_1SC[282]+params.lambda_1SC[283]+params.lambda_1SC[284]+params.lambda_1SC[285]+params.lambda_1SC[286]+params.lambda_1SC[287]+params.lambda_1SC[288]+params.lambda_1SC[289]+params.lambda_1SC[290]+params.lambda_1SC[291]+params.lambda_1SC[292]+params.lambda_1SC[293]+params.lambda_1SC[294]+params.lambda_1SC[295]+params.lambda_1SC[296]+params.lambda_1SC[297]+params.lambda_1SC[298]+params.lambda_1SC[299]+params.lambda_1SC[300]+params.lambda_1SC[301]+params.lambda_1SC[302]+params.lambda_1SC[303]+params.lambda_1SC[304]+params.lambda_1SC[305]+params.lambda_1SC[306]+params.lambda_1SC[307]+params.lambda_1SC[308]+params.lambda_1SC[309]+params.lambda_1SC[310]+params.lambda_1SC[311]+params.lambda_1SC[312]+params.lambda_1SC[313]+params.lambda_1SC[314]+params.lambda_1SC[315]+params.lambda_1SC[316]+params.lambda_1SC[317]+params.lambda_1SC[318]+params.lambda_1SC[319]+params.lambda_1SC[320]+params.lambda_1SC[321]+params.lambda_1SC[322]+params.lambda_1SC[323]+params.lambda_1SC[324]+params.lambda_1SC[325]+params.lambda_1SC[326]+params.lambda_1SC[327]+params.lambda_1SC[328]+params.lambda_1SC[329]+params.lambda_1SC[330]+params.lambda_1SC[331]+params.lambda_1SC[332]+params.lambda_1SC[333]+params.lambda_1SC[334]+params.lambda_1SC[335]+params.lambda_1SC[336]+params.lambda_1SC[337]+params.lambda_1SC[338]+params.lambda_1SC[339]+params.lambda_1SC[340]+params.lambda_1SC[341]+params.lambda_1SC[342]+params.lambda_1SC[343]+params.lambda_1SC[344]+params.lambda_1SC[345]+params.lambda_1SC[346]+params.lambda_1SC[347]+params.lambda_1SC[348]+params.lambda_1SC[349]+params.lambda_1SC[350]+params.lambda_1SC[351]+params.lambda_1SC[352]+params.lambda_1SC[353]+params.lambda_1SC[354]+params.lambda_1SC[355]+params.lambda_1SC[356]+params.lambda_1SC[357]+params.lambda_1SC[358]+params.lambda_1SC[359]+params.lambda_1SC[360]+params.lambda_1SC[361]+params.lambda_1SC[362]+params.lambda_1SC[363]+params.lambda_1SC[364]+params.lambda_1SC[365]+params.lambda_1SC[366]+params.lambda_1SC[367]+params.lambda_1SC[368]+params.lambda_1SC[369]+params.lambda_1SC[370]+params.lambda_1SC[371]+params.lambda_1SC[372]+params.lambda_1SC[373]+params.lambda_1SC[374]+params.lambda_1SC[375]+params.lambda_1SC[376]+params.lambda_1SC[377]+params.lambda_1SC[378]+params.lambda_1SC[379]+params.lambda_1SC[380]+params.lambda_1SC[381]+params.lambda_1SC[382]+params.lambda_1SC[383]+params.lambda_1SC[384]+params.lambda_1SC[385]+params.lambda_1SC[386]+params.lambda_1SC[387]+params.lambda_1SC[388]+params.lambda_1SC[389]+params.lambda_1SC[390]+params.lambda_1SC[391]+params.lambda_1SC[392]+params.lambda_1SC[393]+params.lambda_1SC[394]+params.lambda_1SC[395]+params.lambda_1SC[396]+params.lambda_1SC[397]+params.lambda_1SC[398]+params.lambda_1SC[399]+params.lambda_1SC[400]+params.lambda_1SC[401]+params.lambda_1SC[402]+params.lambda_1SC[403]+params.lambda_1SC[404]+params.lambda_1SC[405]+params.lambda_1SC[406]+params.lambda_1SC[407]+params.lambda_1SC[408]+params.lambda_1SC[409]+params.lambda_1SC[410]+params.lambda_1SC[411]+params.lambda_1SC[412]+params.lambda_1SC[413]+params.lambda_1SC[414]+params.lambda_1SC[415]+params.lambda_1SC[416]+params.lambda_1SC[417]+params.lambda_1SC[418]+params.lambda_1SC[419]+params.lambda_1SC[420]+params.lambda_1SC[421]+params.lambda_1SC[422]+params.lambda_1SC[423]+params.lambda_1SC[424]+params.lambda_1SC[425]+params.lambda_1SC[426]+params.lambda_1SC[427]+params.lambda_1SC[428]+params.lambda_1SC[429]+params.lambda_1SC[430]+params.lambda_1SC[431]+params.lambda_1SC[432]+params.lambda_1SC[433]+params.lambda_1SC[434]+params.lambda_1SC[435]+params.lambda_1SC[436]+params.lambda_1SC[437]+params.lambda_1SC[438]+params.lambda_1SC[439]+params.lambda_1SC[440]+params.lambda_1SC[441]+params.lambda_1SC[442]+params.lambda_1SC[443]+params.lambda_1SC[444]+params.lambda_1SC[445]+params.lambda_1SC[446]+params.lambda_1SC[447]+params.lambda_1SC[448]+params.lambda_1SC[449]+params.lambda_1SC[450]+params.lambda_1SC[451]+params.lambda_1SC[452]+params.lambda_1SC[453]+params.lambda_1SC[454]+params.lambda_1SC[455]+params.lambda_1SC[456]+params.lambda_1SC[457]+params.lambda_1SC[458]+params.lambda_1SC[459]+params.lambda_1SC[460]+params.lambda_1SC[461]+params.lambda_1SC[462]+params.lambda_1SC[463]+params.lambda_1SC[464]+params.lambda_1SC[465]+params.lambda_1SC[466]+params.lambda_1SC[467]+params.lambda_1SC[468]+params.lambda_1SC[469]+params.lambda_1SC[470]+params.lambda_1SC[471]+params.lambda_1SC[472]+params.lambda_1SC[473]+params.lambda_1SC[474]+params.lambda_1SC[475]+params.lambda_1SC[476]+params.lambda_1SC[477]+params.lambda_1SC[478]+params.lambda_1SC[479]+params.lambda_1SC[480]+params.lambda_1SC[481]+params.lambda_1SC[482]+params.lambda_1SC[483]+params.lambda_1SC[484]+params.lambda_1SC[485]+params.lambda_1SC[486]+params.lambda_1SC[487]+params.lambda_1SC[488]+params.lambda_1SC[489]+params.lambda_1SC[490]+params.lambda_1SC[491]+params.lambda_1SC[492]+params.lambda_1SC[493]+params.lambda_1SC[494]+params.lambda_1SC[495]+params.lambda_1SC[496]+params.lambda_1SC[497]+params.lambda_1SC[498]+params.lambda_1SC[499]+params.gamma[0]*params.B[0]-params.lambda_4[0]+2*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*work.frac_121674190848;
  work.q[1] = -2*params.PgPrevNu[0]*work.frac_479733190656+params.gamma[0]*params.A[0]-params.lambda_3[0];
  work.q[2] = 2*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*work.frac_121674190848;
}
void fillh(void) {
  work.h[0] = -params.PgMin[0];
  work.h[1] = params.PgMax[0];
  work.h[2] = -(params.RgMin[0]-params.selectZero[0]*params.PgNext[0]);
  work.h[3] = -(params.selectZero[0]*params.PgNext[0]-params.RgMax[0]);
  work.h[4] = -params.RgMin[0];
  work.h[5] = params.RgMax[0];
}
void fillb(void) {
}
void pre_ops(void) {
  work.frac_479733190656 = params.beta[0];
  work.frac_479733190656 /= 2;
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_399913824256[0] = params.PgPrevNu[0]*work.frac_479733190656*params.PgPrevNu[0];
  work.quad_608032620544[0] = params.PgNu[0]*work.frac_479733190656*params.PgNu[0];
  work.quad_155015135232[0] = (-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*work.frac_121674190848*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0]);
  work.quad_497818112000[0] = (-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*work.frac_121674190848*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0]);
}

// util.cpp
/* Produced by CVXGEN, 2019-01-20 18:36:23 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: util.c. */
/* Description: Common utility file for all cvxgen code. */
#include "solver.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
long global_seed = 1;
static clock_t tic_timestart;
void tic(void) {
  tic_timestart = clock();
}
float toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.2f.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
float tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}
void printmatrix(char *name, double *A, int m, int n, int sparse) {
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
double unif(double lower, double upper) {
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
float ran1(long*idum, int reset) {
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
float randn_internal(long *idum, int reset) {
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
double randn(void) {
  return randn_internal(&global_seed, 0);
}
void reset_rand(void) {
  srand(15);
  global_seed = 1;
  randn_internal(&global_seed, 1);
}

// testsolver.cpp
/* Produced by CVXGEN, 2019-01-20 18:36:23 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"
Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 0
int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif
  set_defaults();
  setup_indexing();
  load_default_data();
  /* Solve problem instance for the record. */
  settings.verbose = 1;
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
  return 0;
}
void load_default_data(void) {
  params.c2[0] = 1.101595805149151;
  params.c1[0] = 1.4162956452362097;
  params.c0[0] = 0.5818094778258887;
  params.beta[0] = 1.021655210395326;
  params.PgPrevNu[0] = 1.5717878173906188;
  params.PgNu[0] = 1.5851723557337523;
  params.gammaSC[0] = 0.2511706209276725;
  params.BSC[0] = -1.171028487447253;
  params.BSC[1] = -1.7941311867966805;
  params.BSC[2] = -0.23676062539745413;
  params.BSC[3] = -1.8804951564857322;
  params.BSC[4] = -0.17266710242115568;
  params.BSC[5] = 0.596576190459043;
  params.BSC[6] = -0.8860508694080989;
  params.BSC[7] = 0.7050196079205251;
  params.BSC[8] = 0.3634512696654033;
  params.BSC[9] = -1.9040724704913385;
  params.BSC[10] = 0.23541635196352795;
  params.BSC[11] = -0.9629902123701384;
  params.BSC[12] = -0.3395952119597214;
  params.BSC[13] = -0.865899672914725;
  params.BSC[14] = 0.7725516732519853;
  params.BSC[15] = -0.23818512931704205;
  params.BSC[16] = -1.372529046100147;
  params.BSC[17] = 0.17859607212737894;
  params.BSC[18] = 1.1212590580454682;
  params.BSC[19] = -0.774545870495281;
  params.BSC[20] = -1.1121684642712744;
  params.BSC[21] = -0.44811496977740495;
  params.BSC[22] = 1.7455345994417217;
  params.BSC[23] = 1.9039816898917352;
  params.BSC[24] = 0.6895347036512547;
  params.BSC[25] = 1.6113364341535923;
  params.BSC[26] = 1.383003485172717;
  params.BSC[27] = -0.48802383468444344;
  params.BSC[28] = -1.631131964513103;
  params.BSC[29] = 0.6136436100941447;
  params.BSC[30] = 0.2313630495538037;
  params.BSC[31] = -0.5537409477496875;
  params.BSC[32] = -1.0997819806406723;
  params.BSC[33] = -0.3739203344950055;
  params.BSC[34] = -0.12423900520332376;
  params.BSC[35] = -0.923057686995755;
  params.BSC[36] = -0.8328289030982696;
  params.BSC[37] = -0.16925440270808823;
  params.BSC[38] = 1.442135651787706;
  params.BSC[39] = 0.34501161787128565;
  params.BSC[40] = -0.8660485502711608;
  params.BSC[41] = -0.8880899735055947;
  params.BSC[42] = -0.1815116979122129;
  params.BSC[43] = -1.17835862158005;
  params.BSC[44] = -1.1944851558277074;
  params.BSC[45] = 0.05614023926976763;
  params.BSC[46] = -1.6510825248767813;
  params.BSC[47] = -0.06565787059365391;
  params.BSC[48] = -0.5512951504486665;
  params.BSC[49] = 0.8307464872626844;
  params.BSC[50] = 0.9869848924080182;
  params.BSC[51] = 0.7643716874230573;
  params.BSC[52] = 0.7567216550196565;
  params.BSC[53] = -0.5055995034042868;
  params.BSC[54] = 0.6725392189410702;
  params.BSC[55] = -0.6406053441727284;
  params.BSC[56] = 0.29117547947550015;
  params.BSC[57] = -0.6967713677405021;
  params.BSC[58] = -0.21941980294587182;
  params.BSC[59] = -1.753884276680243;
  params.BSC[60] = -1.0292983112626475;
  params.BSC[61] = 1.8864104246942706;
  params.BSC[62] = -1.077663182579704;
  params.BSC[63] = 0.7659100437893209;
  params.BSC[64] = 0.6019074328549583;
  params.BSC[65] = 0.8957565577499285;
  params.BSC[66] = -0.09964555746227477;
  params.BSC[67] = 0.38665509840745127;
  params.BSC[68] = -1.7321223042686946;
  params.BSC[69] = -1.7097514487110663;
  params.BSC[70] = -1.2040958948116867;
  params.BSC[71] = -1.3925560119658358;
  params.BSC[72] = -1.5995826216742213;
  params.BSC[73] = -1.4828245415645833;
  params.BSC[74] = 0.21311092723061398;
  params.BSC[75] = -1.248740700304487;
  params.BSC[76] = 1.808404972124833;
  params.BSC[77] = 0.7264471152297065;
  params.BSC[78] = 0.16407869343908477;
  params.BSC[79] = 0.8287224032315907;
  params.BSC[80] = -0.9444533161899464;
  params.BSC[81] = 1.7069027370149112;
  params.BSC[82] = 1.3567722311998827;
  params.BSC[83] = 0.9052779937121489;
  params.BSC[84] = -0.07904017565835986;
  params.BSC[85] = 1.3684127435065871;
  params.BSC[86] = 0.979009293697437;
  params.BSC[87] = 0.6413036255984501;
  params.BSC[88] = 1.6559010680237511;
  params.BSC[89] = 0.5346622551502991;
  params.BSC[90] = -0.5362376605895625;
  params.BSC[91] = 0.2113782926017822;
  params.BSC[92] = -1.2144776931994525;
  params.BSC[93] = -1.2317108144255875;
  params.BSC[94] = 0.9026784957312834;
  params.BSC[95] = 1.1397468137245244;
  params.BSC[96] = 1.8883934547350631;
  params.BSC[97] = 1.4038856681660068;
  params.BSC[98] = 0.17437730638329096;
  params.BSC[99] = -1.6408365219077408;
  params.BSC[100] = -0.04450702153554875;
  params.BSC[101] = 1.7117453902485025;
  params.BSC[102] = 1.1504727980139053;
  params.BSC[103] = -0.05962309578364744;
  params.BSC[104] = -0.1788825540764547;
  params.BSC[105] = -1.1280569263625857;
  params.BSC[106] = -1.2911464767927057;
  params.BSC[107] = -1.7055053231225696;
  params.BSC[108] = 1.56957275034837;
  params.BSC[109] = 0.5607064675962357;
  params.BSC[110] = -1.4266707301147146;
  params.BSC[111] = -0.3434923211351708;
  params.BSC[112] = -1.8035643024085055;
  params.BSC[113] = -1.1625066019105454;
  params.BSC[114] = 0.9228324965161532;
  params.BSC[115] = 0.6044910817663975;
  params.BSC[116] = -0.0840868104920891;
  params.BSC[117] = -0.900877978017443;
  params.BSC[118] = 0.608892500264739;
  params.BSC[119] = 1.8257980452695217;
  params.BSC[120] = -0.25791777529922877;
  params.BSC[121] = -1.7194699796493191;
  params.BSC[122] = -1.7690740487081298;
  params.BSC[123] = -1.6685159248097703;
  params.BSC[124] = 1.8388287490128845;
  params.BSC[125] = 0.16304334474597537;
  params.BSC[126] = 1.3498497306788897;
  params.BSC[127] = -1.3198658230514613;
  params.BSC[128] = -0.9586197090843394;
  params.BSC[129] = 0.7679100474913709;
  params.BSC[130] = 1.5822813125679343;
  params.BSC[131] = -0.6372460621593619;
  params.BSC[132] = -1.741307208038867;
  params.BSC[133] = 1.456478677642575;
  params.BSC[134] = -0.8365102166820959;
  params.BSC[135] = 0.9643296255982503;
  params.BSC[136] = -1.367865381194024;
  params.BSC[137] = 0.7798537405635035;
  params.BSC[138] = 1.3656784761245926;
  params.BSC[139] = 0.9086083149868371;
  params.BSC[140] = -0.5635699005460344;
  params.BSC[141] = 0.9067590059607915;
  params.BSC[142] = -1.4421315032701587;
  params.BSC[143] = -0.7447235390671119;
  params.BSC[144] = -0.32166897326822186;
  params.BSC[145] = 1.5088481557772684;
  params.BSC[146] = -1.385039165715428;
  params.BSC[147] = 1.5204991609972622;
  params.BSC[148] = 1.1958572768832156;
  params.BSC[149] = 1.8864971883119228;
  params.BSC[150] = -0.5291880667861584;
  params.BSC[151] = -1.1802409243688836;
  params.BSC[152] = -1.037718718661604;
  params.BSC[153] = 1.3114512056856835;
  params.BSC[154] = 1.8609125943756615;
  params.BSC[155] = 0.7952399935216938;
  params.BSC[156] = -0.07001183290468038;
  params.BSC[157] = -0.8518009412754686;
  params.BSC[158] = 1.3347515373726386;
  params.BSC[159] = 1.4887180335977037;
  params.BSC[160] = -1.6314736327976336;
  params.BSC[161] = -1.1362021159208933;
  params.BSC[162] = 1.327044361831466;
  params.BSC[163] = 1.3932155883179842;
  params.BSC[164] = -0.7413880049440107;
  params.BSC[165] = -0.8828216126125747;
  params.BSC[166] = -0.27673991192616;
  params.BSC[167] = 0.15778600105866714;
  params.BSC[168] = -1.6177327399735457;
  params.BSC[169] = 1.3476485548544606;
  params.BSC[170] = 0.13893948140528378;
  params.BSC[171] = 1.0998712601636944;
  params.BSC[172] = -1.0766549376946926;
  params.BSC[173] = 1.8611734044254629;
  params.BSC[174] = 1.0041092292735172;
  params.BSC[175] = -0.6276245424321543;
  params.BSC[176] = 1.794110587839819;
  params.BSC[177] = 0.8020471158650913;
  params.BSC[178] = 1.362244341944948;
  params.BSC[179] = -1.8180107765765245;
  params.BSC[180] = -1.7774338357932473;
  params.BSC[181] = 0.9709490941985153;
  params.BSC[182] = -0.7812542682064318;
  params.BSC[183] = 0.0671374633729811;
  params.BSC[184] = -1.374950305314906;
  params.BSC[185] = 1.9118096386279388;
  params.BSC[186] = 0.011004190697677885;
  params.BSC[187] = 1.3160043138989015;
  params.BSC[188] = -1.7038488148800144;
  params.BSC[189] = -0.08433819112864738;
  params.BSC[190] = -1.7508820783768964;
  params.BSC[191] = 1.536965724350949;
  params.BSC[192] = -0.21675928514816478;
  params.BSC[193] = -1.725800326952653;
  params.BSC[194] = -1.6940148707361717;
  params.BSC[195] = 0.15517063201268;
  params.BSC[196] = -1.697734381979077;
  params.BSC[197] = -1.264910727950229;
  params.BSC[198] = -0.2545716633339441;
  params.BSC[199] = -0.008868675926170244;
  params.BSC[200] = 0.3332476609670296;
  params.BSC[201] = 0.48205072561962936;
  params.BSC[202] = -0.5087540014293261;
  params.BSC[203] = 0.4749463319223195;
  params.BSC[204] = -1.371021366459455;
  params.BSC[205] = -0.8979660982652256;
  params.BSC[206] = 1.194873082385242;
  params.BSC[207] = -1.3876427970939353;
  params.BSC[208] = -1.106708108457053;
  params.BSC[209] = -1.0280872812241797;
  params.BSC[210] = -0.08197078070773234;
  params.BSC[211] = -1.9970179118324083;
  params.BSC[212] = -1.878754557910134;
  params.BSC[213] = -0.15380739340877803;
  params.BSC[214] = -1.349917260533923;
  params.BSC[215] = 0.7180072150931407;
  params.BSC[216] = 1.1808183487065538;
  params.BSC[217] = 0.31265343495084075;
  params.BSC[218] = 0.7790599086928229;
  params.BSC[219] = -0.4361679370644853;
  params.BSC[220] = -1.8148151880282066;
  params.BSC[221] = -0.24231386948140266;
  params.BSC[222] = -0.5120787511622411;
  params.BSC[223] = 0.3880129688013203;
  params.BSC[224] = -1.4631273212038676;
  params.BSC[225] = -1.0891484131126563;
  params.BSC[226] = 1.2591296661091191;
  params.BSC[227] = -0.9426978934391474;
  params.BSC[228] = -0.358719180371347;
  params.BSC[229] = 1.7438887059831263;
  params.BSC[230] = -0.8977901479165817;
  params.BSC[231] = -1.4188401645857445;
  params.BSC[232] = 0.8080805173258092;
  params.BSC[233] = 0.2682662017650985;
  params.BSC[234] = 0.44637534218638786;
  params.BSC[235] = -1.8318765960257055;
  params.BSC[236] = -0.3309324209710929;
  params.BSC[237] = -1.9829342633313622;
  params.BSC[238] = -1.013858124556442;
  params.BSC[239] = 0.8242247343360254;
  params.BSC[240] = -1.753837136317201;
  params.BSC[241] = -0.8212260055868805;
  params.BSC[242] = 1.9524510112487126;
  params.BSC[243] = 1.884888920907902;
  params.BSC[244] = -0.0726144452811801;
  params.BSC[245] = 0.9427735461129836;
  params.BSC[246] = 0.5306230967445558;
  params.BSC[247] = -0.1372277142250531;
  params.BSC[248] = 1.4282657305652786;
  params.BSC[249] = -1.309926991335284;
  params.BSC[250] = 1.3137276889764422;
  params.BSC[251] = -1.8317219061667278;
  params.BSC[252] = 1.4678147672511939;
  params.BSC[253] = 0.703986349872991;
  params.BSC[254] = -0.2163435603565258;
  params.BSC[255] = 0.6862809905371079;
  params.BSC[256] = -0.15852598444303245;
  params.BSC[257] = 1.1200128895143409;
  params.BSC[258] = -1.5462236645435308;
  params.BSC[259] = 0.0326297153944215;
  params.BSC[260] = 1.4859581597754916;
  params.BSC[261] = 1.71011710324809;
  params.BSC[262] = -1.1186546738067493;
  params.BSC[263] = -0.9922787897815244;
  params.BSC[264] = 1.6160498864359547;
  params.BSC[265] = -0.6179306451394861;
  params.BSC[266] = -1.7725097038051376;
  params.BSC[267] = 0.8595466884481313;
  params.BSC[268] = -0.3423245633865686;
  params.BSC[269] = 0.9412967499805762;
  params.BSC[270] = -0.09163346622652258;
  params.BSC[271] = 0.002262217745727657;
  params.BSC[272] = -0.3297523583656421;
  params.BSC[273] = -0.8380604158593941;
  params.BSC[274] = 1.6028434695494038;
  params.BSC[275] = 0.675150311940429;
  params.BSC[276] = 1.1553293733718686;
  params.BSC[277] = 1.5829581243724693;
  params.BSC[278] = -0.9992442304425597;
  params.BSC[279] = 1.6792824558896897;
  params.BSC[280] = 1.4504203490342324;
  params.BSC[281] = 0.02434104849994556;
  params.BSC[282] = 0.27160869657612263;
  params.BSC[283] = -1.5402710478528858;
  params.BSC[284] = 1.0484633622310744;
  params.BSC[285] = -1.3070999712627054;
  params.BSC[286] = 0.13534416402363814;
  params.BSC[287] = -1.4942507790851232;
  params.BSC[288] = -1.708331625671371;
  params.BSC[289] = 0.436109775042258;
  params.BSC[290] = -0.03518748153727991;
  params.BSC[291] = 0.6992397389570906;
  params.BSC[292] = 1.1634167322171374;
  params.BSC[293] = 1.9307499705822648;
  params.BSC[294] = -1.6636772756932747;
  params.BSC[295] = 0.5248484497343218;
  params.BSC[296] = 0.30789958152579144;
  params.BSC[297] = 0.602568707166812;
  params.BSC[298] = 0.17271781925751872;
  params.BSC[299] = 0.2294695501208066;
  params.BSC[300] = 1.4742185345619543;
  params.BSC[301] = -0.1919535345136989;
  params.BSC[302] = 0.13990231452144553;
  params.BSC[303] = 0.7638548150610602;
  params.BSC[304] = -1.6420200344195646;
  params.BSC[305] = -0.27229872445076087;
  params.BSC[306] = -1.5914631171820468;
  params.BSC[307] = -1.4487604283558668;
  params.BSC[308] = -1.991497766136364;
  params.BSC[309] = -1.1611742553535152;
  params.BSC[310] = -1.133450950247063;
  params.BSC[311] = 0.06497792493777155;
  params.BSC[312] = 0.28083295396097263;
  params.BSC[313] = 1.2958447220129887;
  params.BSC[314] = -0.05315524470737154;
  params.BSC[315] = 1.5658183956871667;
  params.BSC[316] = -0.41975684089933685;
  params.BSC[317] = 0.97844578833777;
  params.BSC[318] = 0.2110290496695293;
  params.BSC[319] = 0.4953003430893044;
  params.BSC[320] = -0.9184320124667495;
  params.BSC[321] = 1.750380031759156;
  params.BSC[322] = 1.0786188614315915;
  params.BSC[323] = -1.4176198837203735;
  params.BSC[324] = 0.149737479778294;
  params.BSC[325] = 1.9831452222223418;
  params.BSC[326] = -1.8037746699794734;
  params.BSC[327] = -0.7887206483295461;
  params.BSC[328] = 0.9632534854086652;
  params.BSC[329] = -1.8425542093895406;
  params.BSC[330] = 0.986684363969033;
  params.BSC[331] = 0.2936851199350441;
  params.BSC[332] = 0.9268227022482662;
  params.BSC[333] = 0.20333038350653299;
  params.BSC[334] = 1.7576139132046351;
  params.BSC[335] = -0.614393188398918;
  params.BSC[336] = 0.297877839744912;
  params.BSC[337] = -1.796880083990895;
  params.BSC[338] = 0.21373133661742738;
  params.BSC[339] = -0.32242822540825156;
  params.BSC[340] = 1.9326471511608059;
  params.BSC[341] = 1.7824292753481785;
  params.BSC[342] = -1.4468823405675986;
  params.BSC[343] = -1.8335374338761512;
  params.BSC[344] = -1.5172997317243713;
  params.BSC[345] = -1.229012129120719;
  params.BSC[346] = 0.9046719772422094;
  params.BSC[347] = 0.17591181415489432;
  params.BSC[348] = 0.13970133814112584;
  params.BSC[349] = -0.14185208214985234;
  params.BSC[350] = -1.9732231264739348;
  params.BSC[351] = -0.4301123458221334;
  params.BSC[352] = 1.9957537650387742;
  params.BSC[353] = 1.2811648216477893;
  params.BSC[354] = 0.2914428437588219;
  params.BSC[355] = -1.214148157218884;
  params.BSC[356] = 1.6818776980374155;
  params.BSC[357] = -0.30341101038214635;
  params.BSC[358] = 0.47730909231793106;
  params.BSC[359] = -1.187569373035299;
  params.BSC[360] = -0.6877370247915531;
  params.BSC[361] = -0.6201861482616171;
  params.BSC[362] = -0.4209925183921568;
  params.BSC[363] = -1.9110724537712471;
  params.BSC[364] = 0.6413882087807936;
  params.BSC[365] = -1.3200399280087032;
  params.BSC[366] = 0.41320105301312626;
  params.BSC[367] = 0.4783213861392275;
  params.BSC[368] = 0.7916189857293743;
  params.BSC[369] = -0.8322752558146558;
  params.BSC[370] = -0.8318720537426154;
  params.BSC[371] = 1.0221179076113445;
  params.BSC[372] = -0.4471032189262627;
  params.BSC[373] = -1.3901469561676985;
  params.BSC[374] = 1.6210596051208572;
  params.BSC[375] = -1.9476687601912737;
  params.BSC[376] = 1.5459376306231292;
  params.BSC[377] = -0.830972896191656;
  params.BSC[378] = -0.47269983955176276;
  params.BSC[379] = 1.913620609584223;
  params.BSC[380] = -0.25329703423935124;
  params.BSC[381] = 0.8635279149674653;
  params.BSC[382] = -0.35046893227111564;
  params.BSC[383] = 1.6541432486772365;
  params.BSC[384] = 0.8779619968413503;
  params.BSC[385] = -0.07723284625844862;
  params.BSC[386] = -1.6631134040635196;
  params.BSC[387] = -0.54546452868516;
  params.BSC[388] = -0.03757319061095998;
  params.BSC[389] = -0.864543266194465;
  params.BSC[390] = 0.13856203767859343;
  params.BSC[391] = -1.1613957272733684;
  params.BSC[392] = -0.022681697832835024;
  params.BSC[393] = 0.11202078062843634;
  params.BSC[394] = 0.6934385624164641;
  params.BSC[395] = 0.9814633803279791;
  params.BSC[396] = 0.9198949681022897;
  params.BSC[397] = -0.3035363988458051;
  params.BSC[398] = -0.1761906755724203;
  params.BSC[399] = 1.4940284058791686;
  params.BSC[400] = -0.5488483097174393;
  params.BSC[401] = 0.9521313238305416;
  params.BSC[402] = 1.9762689267600413;
  params.BSC[403] = 1.6992335341478482;
  params.BSC[404] = 0.1969474711697119;
  params.BSC[405] = -0.7795544525014559;
  params.BSC[406] = 0.4892505434034007;
  params.BSC[407] = 0.7372066729248594;
  params.BSC[408] = 0.10784901966517557;
  params.BSC[409] = -0.6340934767066218;
  params.BSC[410] = -0.17829371464242083;
  params.BSC[411] = -1.6728370279392784;
  params.BSC[412] = -0.8348711800042916;
  params.BSC[413] = -1.4204129800590897;
  params.BSC[414] = 0.6659229232859376;
  params.BSC[415] = 1.8369365661533168;
  params.BSC[416] = -1.371061267737546;
  params.BSC[417] = -1.8868237125934915;
  params.BSC[418] = 0.9654286768651104;
  params.BSC[419] = -0.5833420409292005;
  params.BSC[420] = 0.02386510653728502;
  params.BSC[421] = -1.7558076992858345;
  params.BSC[422] = -1.2889402130475411;
  params.BSC[423] = 0.7820251677632606;
  params.BSC[424] = 0.4208424784688227;
  params.BSC[425] = 1.4136448896755982;
  params.BSC[426] = 1.8516928541530757;
  params.BSC[427] = -0.5615396035790421;
  params.BSC[428] = 0.4809940266433177;
  params.BSC[429] = -0.20929035114697303;
  params.BSC[430] = 0.022387850798402553;
  params.BSC[431] = -0.43399296564115764;
  params.BSC[432] = 1.9095769077945013;
  params.BSC[433] = 0.4945512698336847;
  params.BSC[434] = -1.4324582900293557;
  params.BSC[435] = 0.790913765746676;
  params.BSC[436] = 1.8630250293383734;
  params.BSC[437] = 1.5793975466121069;
  params.BSC[438] = 0.2320163334712646;
  params.BSC[439] = -1.9411408650055968;
  params.BSC[440] = 1.2221853270725478;
  params.BSC[441] = 1.7274453600045607;
  params.BSC[442] = 0.9357159281665783;
  params.BSC[443] = -0.2841874908331623;
  params.BSC[444] = -0.4766355664552626;
  params.BSC[445] = 0.9784190546201912;
  params.BSC[446] = -1.5685956114005477;
  params.BSC[447] = 1.1387833891036;
  params.BSC[448] = -0.004779126480003892;
  params.BSC[449] = -1.7195239474925414;
  params.BSC[450] = 1.2921808565147272;
  params.BSC[451] = -0.43317009071966606;
  params.BSC[452] = -1.572940257279357;
  params.BSC[453] = -1.3048062231674988;
  params.BSC[454] = 1.4377304631579175;
  params.BSC[455] = -1.3090328020145874;
  params.BSC[456] = 1.1370018620707785;
  params.BSC[457] = 1.2164644012668289;
  params.BSC[458] = -1.6539274174499985;
  params.BSC[459] = -0.25845368809725544;
  params.BSC[460] = 1.1486358936399745;
  params.BSC[461] = -0.03975647517318137;
  params.BSC[462] = 1.4640632749164326;
  params.BSC[463] = -0.48111499989733186;
  params.BSC[464] = 0.5132576752843594;
  params.BSC[465] = -1.1459189400462249;
  params.BSC[466] = 1.3690255364554855;
  params.BSC[467] = 1.3574291456003253;
  params.BSC[468] = 0.26333733823037253;
  params.BSC[469] = -0.7076462135286032;
  params.BSC[470] = -0.6097272363453645;
  params.BSC[471] = 0.37873096815108465;
  params.BSC[472] = -1.4863636934585411;
  params.BSC[473] = 0.04189135833804869;
  params.BSC[474] = -0.8182949160834703;
  params.BSC[475] = -0.6336865828985854;
  params.BSC[476] = -0.7126437991119396;
  params.BSC[477] = 1.3381487344587226;
  params.BSC[478] = -1.2979975504895949;
  params.BSC[479] = -1.0542097271412714;
  params.BSC[480] = -1.3421003125955435;
  params.BSC[481] = -1.9395969070507038;
  params.BSC[482] = -0.29758108058547306;
  params.BSC[483] = 1.3757899684264032;
  params.BSC[484] = 1.6109970296148042;
  params.BSC[485] = -0.050537352418498216;
  params.BSC[486] = -0.3144945653528741;
  params.BSC[487] = 1.4726689240031474;
  params.BSC[488] = 0.11397910876468265;
  params.BSC[489] = 0.19466869962815858;
  params.BSC[490] = 0.5972476722406035;
  params.BSC[491] = -1.6815490772221828;
  params.BSC[492] = 1.3540223072599735;
  params.BSC[493] = -1.577027832358222;
  params.BSC[494] = 0.12928618615237353;
  params.BSC[495] = 1.704038169667271;
  params.BSC[496] = 0.19482725189070793;
  params.BSC[497] = -0.6311686254597215;
  params.BSC[498] = 0.9065234706582928;
  params.BSC[499] = 1.604058201281767;
  params.lambda_1SC[0] = 0.4649414640474294;
  params.lambda_1SC[1] = -1.7764554290993346;
  params.lambda_1SC[2] = 1.5152343936830337;
  params.lambda_1SC[3] = -1.9280901945449935;
  params.lambda_1SC[4] = 0.7129569482366098;
  params.lambda_1SC[5] = 1.6001840923928201;
  params.lambda_1SC[6] = -1.3702177446733126;
  params.lambda_1SC[7] = 0.11266051920028186;
  params.lambda_1SC[8] = 0.8202183589903962;
  params.lambda_1SC[9] = -1.297953481011172;
  params.lambda_1SC[10] = -1.0192096617939002;
  params.lambda_1SC[11] = -1.7337200441949867;
  params.lambda_1SC[12] = -1.3639899659742465;
  params.lambda_1SC[13] = -1.5273517222086332;
  params.lambda_1SC[14] = -0.8374302703303731;
  params.lambda_1SC[15] = 1.00229367551592;
  params.lambda_1SC[16] = 0.7747378843920099;
  params.lambda_1SC[17] = 1.0504096866871468;
  params.lambda_1SC[18] = 0.638655773812761;
  params.lambda_1SC[19] = 1.176936790033046;
  params.lambda_1SC[20] = -1.4041747524796162;
  params.lambda_1SC[21] = 0.21725437512222667;
  params.lambda_1SC[22] = -1.9141609882936188;
  params.lambda_1SC[23] = -0.03334441105363828;
  params.lambda_1SC[24] = 1.3736673884387467;
  params.lambda_1SC[25] = -0.11085150689269163;
  params.lambda_1SC[26] = -0.8176560931958075;
  params.lambda_1SC[27] = -0.9013799953302866;
  params.lambda_1SC[28] = -0.42583422050124753;
  params.lambda_1SC[29] = 1.6552920005330618;
  params.lambda_1SC[30] = 1.8971842560697287;
  params.lambda_1SC[31] = 0.9935321777966784;
  params.lambda_1SC[32] = 1.9500402929402196;
  params.lambda_1SC[33] = 1.0489535977170181;
  params.lambda_1SC[34] = -0.8630392743714372;
  params.lambda_1SC[35] = -0.25967183338596733;
  params.lambda_1SC[36] = 0.8925966402843359;
  params.lambda_1SC[37] = 0.8373600738876834;
  params.lambda_1SC[38] = 0.7125001994938436;
  params.lambda_1SC[39] = -0.048447588572545275;
  params.lambda_1SC[40] = -1.4274714856193604;
  params.lambda_1SC[41] = 1.8385542904833923;
  params.lambda_1SC[42] = -1.1195070325474288;
  params.lambda_1SC[43] = 1.9175373793884956;
  params.lambda_1SC[44] = -1.49030500627704;
  params.lambda_1SC[45] = 1.9213425364706396;
  params.lambda_1SC[46] = -0.49553546476315047;
  params.lambda_1SC[47] = 1.2437464435895134;
  params.lambda_1SC[48] = -1.970831509470568;
  params.lambda_1SC[49] = -0.219996830259797;
  params.lambda_1SC[50] = -1.0042329091607591;
  params.lambda_1SC[51] = 0.7781008085794774;
  params.lambda_1SC[52] = 0.65210699599452;
  params.lambda_1SC[53] = -0.152326999732443;
  params.lambda_1SC[54] = 0.8265434509993406;
  params.lambda_1SC[55] = 1.9130464561754126;
  params.lambda_1SC[56] = -1.6270096836882288;
  params.lambda_1SC[57] = 0.2507042290048189;
  params.lambda_1SC[58] = 0.7038441998600256;
  params.lambda_1SC[59] = 0.5328743207925606;
  params.lambda_1SC[60] = -0.9509907719589208;
  params.lambda_1SC[61] = 1.499815178589135;
  params.lambda_1SC[62] = -1.0178753663037017;
  params.lambda_1SC[63] = 1.3798461831617561;
  params.lambda_1SC[64] = -0.11708553759234386;
  params.lambda_1SC[65] = -1.4276299186218124;
  params.lambda_1SC[66] = 1.296518419303864;
  params.lambda_1SC[67] = -1.6872707956138546;
  params.lambda_1SC[68] = 1.1799585157870145;
  params.lambda_1SC[69] = 0.4000488706320535;
  params.lambda_1SC[70] = 1.506638004200894;
  params.lambda_1SC[71] = 1.2128180682740366;
  params.lambda_1SC[72] = -0.39211699471717854;
  params.lambda_1SC[73] = -1.4592313874139302;
  params.lambda_1SC[74] = -0.9352340128154211;
  params.lambda_1SC[75] = -1.994709862977336;
  params.lambda_1SC[76] = 0.6136129920637026;
  params.lambda_1SC[77] = -1.6579503948780245;
  params.lambda_1SC[78] = -1.2828456921062488;
  params.lambda_1SC[79] = -1.0200938896697522;
  params.lambda_1SC[80] = -0.3755900704115436;
  params.lambda_1SC[81] = 0.747199791836243;
  params.lambda_1SC[82] = -0.22212974213441683;
  params.lambda_1SC[83] = 0.015082263441096089;
  params.lambda_1SC[84] = -1.6271688108937168;
  params.lambda_1SC[85] = -0.6472903955867526;
  params.lambda_1SC[86] = -1.1733258209627806;
  params.lambda_1SC[87] = 0.9565501943340924;
  params.lambda_1SC[88] = -1.929389541307601;
  params.lambda_1SC[89] = 0.4671837668673531;
  params.lambda_1SC[90] = 0.7915477026785647;
  params.lambda_1SC[91] = 0.018572068486599758;
  params.lambda_1SC[92] = -1.8220899973808726;
  params.lambda_1SC[93] = -0.995629851336445;
  params.lambda_1SC[94] = -1.0486975119711213;
  params.lambda_1SC[95] = -0.9289312699596386;
  params.lambda_1SC[96] = -0.9472402942019333;
  params.lambda_1SC[97] = 1.8908619466142156;
  params.lambda_1SC[98] = 1.164645007668001;
  params.lambda_1SC[99] = 1.5636429264767182;
  params.lambda_1SC[100] = 0.8540115800503387;
  params.lambda_1SC[101] = -0.6133530465568309;
  params.lambda_1SC[102] = 1.7674136894457204;
  params.lambda_1SC[103] = -0.06217940181271242;
  params.lambda_1SC[104] = -1.2582602406204213;
  params.lambda_1SC[105] = 0.9179968784775836;
  params.lambda_1SC[106] = -0.9627796203753647;
  params.lambda_1SC[107] = 1.2911416493727805;
  params.lambda_1SC[108] = 0.9619156621267284;
  params.lambda_1SC[109] = -0.8391987363014124;
  params.lambda_1SC[110] = -0.16142857857315818;
  params.lambda_1SC[111] = 0.8603892868304936;
  params.lambda_1SC[112] = 0.672061858055037;
  params.lambda_1SC[113] = 0.10631385676272265;
  params.lambda_1SC[114] = -1.1434283104802896;
  params.lambda_1SC[115] = -0.7024280087663541;
  params.lambda_1SC[116] = 0.7791723379458899;
  params.lambda_1SC[117] = -0.17206766925671246;
  params.lambda_1SC[118] = 0.8714406054415362;
  params.lambda_1SC[119] = 0.7364640800101268;
  params.lambda_1SC[120] = -0.577393318625969;
  params.lambda_1SC[121] = -1.603607371381821;
  params.lambda_1SC[122] = 0.7231454736647596;
  params.lambda_1SC[123] = -0.5776666119800344;
  params.lambda_1SC[124] = 0.25985922282642804;
  params.lambda_1SC[125] = -1.500019293846674;
  params.lambda_1SC[126] = -1.41591503759888;
  params.lambda_1SC[127] = -0.30464385789747794;
  params.lambda_1SC[128] = 0.677515340905404;
  params.lambda_1SC[129] = -1.5301412809058377;
  params.lambda_1SC[130] = 1.097788736551506;
  params.lambda_1SC[131] = 1.4054563154505293;
  params.lambda_1SC[132] = 0.6904915185274869;
  params.lambda_1SC[133] = 0.9984361169236493;
  params.lambda_1SC[134] = -1.0460788838474921;
  params.lambda_1SC[135] = -1.5989319614177124;
  params.lambda_1SC[136] = -0.6834813660758638;
  params.lambda_1SC[137] = -1.4978328637140224;
  params.lambda_1SC[138] = -0.3340404173113156;
  params.lambda_1SC[139] = 1.044497402438696;
  params.lambda_1SC[140] = -0.875611719278079;
  params.lambda_1SC[141] = 1.4233779191761733;
  params.lambda_1SC[142] = -0.1880612910960302;
  params.lambda_1SC[143] = -1.3523791242997114;
  params.lambda_1SC[144] = 0.5691200673315562;
  params.lambda_1SC[145] = -0.24590364206081272;
  params.lambda_1SC[146] = -0.6790819241936314;
  params.lambda_1SC[147] = 0.06554105230580287;
  params.lambda_1SC[148] = 1.9642897275976492;
  params.lambda_1SC[149] = 1.0075323744403706;
  params.lambda_1SC[150] = -0.8257682212557649;
  params.lambda_1SC[151] = 1.7097592474973915;
  params.lambda_1SC[152] = -1.2633370473270409;
  params.lambda_1SC[153] = -0.3674317265957998;
  params.lambda_1SC[154] = -0.5096221670767425;
  params.lambda_1SC[155] = 1.9427867788797188;
  params.lambda_1SC[156] = -1.9819265272693376;
  params.lambda_1SC[157] = 0.59706237040941;
  params.lambda_1SC[158] = 0.03464508712113412;
  params.lambda_1SC[159] = -1.3762622976910213;
  params.lambda_1SC[160] = 0.5667125534930704;
  params.lambda_1SC[161] = 1.2314327073557654;
  params.lambda_1SC[162] = 1.067926609423298;
  params.lambda_1SC[163] = -1.6894983403969057;
  params.lambda_1SC[164] = 1.1803977296204375;
  params.lambda_1SC[165] = -0.6163708084181705;
  params.lambda_1SC[166] = 1.410387054314267;
  params.lambda_1SC[167] = 0.2483811753322378;
  params.lambda_1SC[168] = -1.6802703802890249;
  params.lambda_1SC[169] = 0.979002476278382;
  params.lambda_1SC[170] = 0.25437753015395215;
  params.lambda_1SC[171] = -0.19429046032083797;
  params.lambda_1SC[172] = 0.5538337109514373;
  params.lambda_1SC[173] = -0.22516774177478682;
  params.lambda_1SC[174] = -1.8522987901467385;
  params.lambda_1SC[175] = 1.6508179296419305;
  params.lambda_1SC[176] = -0.8487785288135496;
  params.lambda_1SC[177] = -0.8826237616160246;
  params.lambda_1SC[178] = -1.8011033482653875;
  params.lambda_1SC[179] = 0.8695812040274169;
  params.lambda_1SC[180] = 1.6378165319423914;
  params.lambda_1SC[181] = -0.1481829774264356;
  params.lambda_1SC[182] = 1.2564969039671574;
  params.lambda_1SC[183] = -0.684739163169866;
  params.lambda_1SC[184] = -0.4610590690471823;
  params.lambda_1SC[185] = -0.8331048225539304;
  params.lambda_1SC[186] = -0.8202459090544569;
  params.lambda_1SC[187] = 0.24515655655414248;
  params.lambda_1SC[188] = -0.6962337393386173;
  params.lambda_1SC[189] = -1.7368155738752113;
  params.lambda_1SC[190] = 0.7610728396291946;
  params.lambda_1SC[191] = -0.765449239026013;
  params.lambda_1SC[192] = -0.43085918845486093;
  params.lambda_1SC[193] = -1.2548272353144259;
  params.lambda_1SC[194] = -1.8025410663253707;
  params.lambda_1SC[195] = 0.6171321070899451;
  params.lambda_1SC[196] = -0.3071517457288193;
  params.lambda_1SC[197] = 1.7295616144888504;
  params.lambda_1SC[198] = -1.924302635086176;
  params.lambda_1SC[199] = 1.299583715888006;
  params.lambda_1SC[200] = -1.9258432649268848;
  params.lambda_1SC[201] = -1.308131173573865;
  params.lambda_1SC[202] = -1.944655670458073;
  params.lambda_1SC[203] = 0.03753627429395401;
  params.lambda_1SC[204] = -0.5851814302714002;
  params.lambda_1SC[205] = -0.7266056329999255;
  params.lambda_1SC[206] = 1.8975656076388434;
  params.lambda_1SC[207] = 0.3121902737027811;
  params.lambda_1SC[208] = 1.7053730032601906;
  params.lambda_1SC[209] = 1.555838365783658;
  params.lambda_1SC[210] = -0.641731508455837;
  params.lambda_1SC[211] = 0.6141829877795306;
  params.lambda_1SC[212] = -1.2197849909613159;
  params.lambda_1SC[213] = -0.5380071919857308;
  params.lambda_1SC[214] = 0.739143852117047;
  params.lambda_1SC[215] = 0.2752842971358995;
  params.lambda_1SC[216] = -1.6601930270596967;
  params.lambda_1SC[217] = -1.6874428445363736;
  params.lambda_1SC[218] = 1.0840077884635324;
  params.lambda_1SC[219] = 0.8042461357116335;
  params.lambda_1SC[220] = 0.5943031843200117;
  params.lambda_1SC[221] = -1.1167184421998662;
  params.lambda_1SC[222] = 1.1352733713960155;
  params.lambda_1SC[223] = 0.3077638771814004;
  params.lambda_1SC[224] = -1.4115482222137845;
  params.lambda_1SC[225] = -1.3353078694260594;
  params.lambda_1SC[226] = 0.48320878332770656;
  params.lambda_1SC[227] = 1.9130931062707663;
  params.lambda_1SC[228] = 0.34428617274287987;
  params.lambda_1SC[229] = -0.819331379598593;
  params.lambda_1SC[230] = 0.8644114302172246;
  params.lambda_1SC[231] = -1.1871445028885255;
  params.lambda_1SC[232] = 0.5654759806166565;
  params.lambda_1SC[233] = -0.34438329868703654;
  params.lambda_1SC[234] = -1.8767908666912492;
  params.lambda_1SC[235] = 0.06559655821858801;
  params.lambda_1SC[236] = 0.030042529289834974;
  params.lambda_1SC[237] = 1.6476947083696136;
  params.lambda_1SC[238] = 0.843474742904093;
  params.lambda_1SC[239] = -0.7846732095313289;
  params.lambda_1SC[240] = 0.5806097744148899;
  params.lambda_1SC[241] = -0.9999926236517398;
  params.lambda_1SC[242] = 0.616600769909577;
  params.lambda_1SC[243] = 0.40576562725410703;
  params.lambda_1SC[244] = 0.4579434796746926;
  params.lambda_1SC[245] = -1.3554959003596077;
  params.lambda_1SC[246] = 0.842737370393269;
  params.lambda_1SC[247] = 0.15089405206073003;
  params.lambda_1SC[248] = 0.9692197132519018;
  params.lambda_1SC[249] = -1.8213178505691485;
  params.lambda_1SC[250] = 0.9320547684014473;
  params.lambda_1SC[251] = 1.8010468183262986;
  params.lambda_1SC[252] = 0.34260384914401065;
  params.lambda_1SC[253] = -1.439404381078433;
  params.lambda_1SC[254] = 1.6950932713659075;
  params.lambda_1SC[255] = 0.22279836918864637;
  params.lambda_1SC[256] = -0.8903349770781586;
  params.lambda_1SC[257] = 0.48856177553753355;
  params.lambda_1SC[258] = 1.4355721493163385;
  params.lambda_1SC[259] = -1.566804863345856;
  params.lambda_1SC[260] = 0.23865167716993874;
  params.lambda_1SC[261] = 0.5309468550013152;
  params.lambda_1SC[262] = 1.7044705300707474;
  params.lambda_1SC[263] = -1.847248469391193;
  params.lambda_1SC[264] = 0.02972819554506545;
  params.lambda_1SC[265] = 0.9304968547048014;
  params.lambda_1SC[266] = 0.6387704531503036;
  params.lambda_1SC[267] = 1.4456832184180497;
  params.lambda_1SC[268] = 0.625745180667252;
  params.lambda_1SC[269] = 0.5265177505358345;
  params.lambda_1SC[270] = 0.43961206692696697;
  params.lambda_1SC[271] = -1.3412493780273906;
  params.lambda_1SC[272] = 1.418721847116943;
  params.lambda_1SC[273] = 0.06342988689319062;
  params.lambda_1SC[274] = -1.9145956108108106;
  params.lambda_1SC[275] = 1.2818619663415407;
  params.lambda_1SC[276] = 0.754165816794254;
  params.lambda_1SC[277] = 0.6296241802331286;
  params.lambda_1SC[278] = -1.792949927289066;
  params.lambda_1SC[279] = 0.8217878144889776;
  params.lambda_1SC[280] = -0.6479435408117138;
  params.lambda_1SC[281] = 1.2579034508883895;
  params.lambda_1SC[282] = -1.3757785464314858;
  params.lambda_1SC[283] = -0.6582651266611776;
  params.lambda_1SC[284] = -1.5875771458755419;
  params.lambda_1SC[285] = -1.9265755661365582;
  params.lambda_1SC[286] = -0.7057348420149179;
  params.lambda_1SC[287] = 0.7975180092519181;
  params.lambda_1SC[288] = -0.32088268527672037;
  params.lambda_1SC[289] = -0.8082275050375718;
  params.lambda_1SC[290] = 1.7450103338687706;
  params.lambda_1SC[291] = 0.9068020986809042;
  params.lambda_1SC[292] = 1.8865025396757575;
  params.lambda_1SC[293] = -0.7882162112934408;
  params.lambda_1SC[294] = -1.2239019904562953;
  params.lambda_1SC[295] = 1.873426511441803;
  params.lambda_1SC[296] = -1.682198810826204;
  params.lambda_1SC[297] = 1.1419698359283421;
  params.lambda_1SC[298] = 1.7891316972075328;
  params.lambda_1SC[299] = 0.8722740536456279;
  params.lambda_1SC[300] = 0.2977725076084803;
  params.lambda_1SC[301] = -0.5241861000703971;
  params.lambda_1SC[302] = -0.8178996238468783;
  params.lambda_1SC[303] = -1.7316157946122606;
  params.lambda_1SC[304] = -0.4845032252049113;
  params.lambda_1SC[305] = 0.011693701254325006;
  params.lambda_1SC[306] = -1.0047763318758451;
  params.lambda_1SC[307] = 1.8784654853736424;
  params.lambda_1SC[308] = 1.9941864300472405;
  params.lambda_1SC[309] = -0.8420756262594038;
  params.lambda_1SC[310] = 0.40100822473216713;
  params.lambda_1SC[311] = 1.5311441352200075;
  params.lambda_1SC[312] = 0.33483870227787005;
  params.lambda_1SC[313] = -0.5925395671108;
  params.lambda_1SC[314] = 1.38348005534032;
  params.lambda_1SC[315] = 0.11490978184210254;
  params.lambda_1SC[316] = 0.47708847918959174;
  params.lambda_1SC[317] = 1.027763410227907;
  params.lambda_1SC[318] = -0.8468811136558307;
  params.lambda_1SC[319] = 0.7447704042795364;
  params.lambda_1SC[320] = -0.9338755969722978;
  params.lambda_1SC[321] = 0.6098347185312418;
  params.lambda_1SC[322] = -0.9735649686562615;
  params.lambda_1SC[323] = 1.1849625913858515;
  params.lambda_1SC[324] = -0.6107093125803296;
  params.lambda_1SC[325] = 1.8417016330277591;
  params.lambda_1SC[326] = -1.3064755864737192;
  params.lambda_1SC[327] = 0.8573085847685804;
  params.lambda_1SC[328] = 1.419896827798702;
  params.lambda_1SC[329] = 1.4850286501940344;
  params.lambda_1SC[330] = 0.08010077466299359;
  params.lambda_1SC[331] = 0.19744445712508796;
  params.lambda_1SC[332] = -0.8231563858898774;
  params.lambda_1SC[333] = 1.4265140784690264;
  params.lambda_1SC[334] = 1.049128669818165;
  params.lambda_1SC[335] = 0.05345768194173983;
  params.lambda_1SC[336] = 0.011929446531647514;
  params.lambda_1SC[337] = 1.6778383073204326;
  params.lambda_1SC[338] = 0.2435581699233258;
  params.lambda_1SC[339] = -1.6246124734146075;
  params.lambda_1SC[340] = -1.5673027875725514;
  params.lambda_1SC[341] = 1.0664660584889862;
  params.lambda_1SC[342] = 0.10332493781246344;
  params.lambda_1SC[343] = 0.22544507411930503;
  params.lambda_1SC[344] = -0.26831321995271074;
  params.lambda_1SC[345] = -1.9092064697978284;
  params.lambda_1SC[346] = 1.212948652923116;
  params.lambda_1SC[347] = -0.4738331727256955;
  params.lambda_1SC[348] = -0.14322528944343116;
  params.lambda_1SC[349] = -0.23716036261764994;
  params.lambda_1SC[350] = -1.4651836630834838;
  params.lambda_1SC[351] = 1.8528675440253215;
  params.lambda_1SC[352] = 1.8376546118537322;
  params.lambda_1SC[353] = -1.960900547345425;
  params.lambda_1SC[354] = -1.7817327704522716;
  params.lambda_1SC[355] = -0.4552931154965534;
  params.lambda_1SC[356] = 1.9271279075012147;
  params.lambda_1SC[357] = -0.07280939033984346;
  params.lambda_1SC[358] = 1.690650811162882;
  params.lambda_1SC[359] = -0.6332117170608655;
  params.lambda_1SC[360] = -0.018669079412589884;
  params.lambda_1SC[361] = 1.315234926662395;
  params.lambda_1SC[362] = 0.6344738466833548;
  params.lambda_1SC[363] = -0.7737875816724431;
  params.lambda_1SC[364] = -0.17341331528078197;
  params.lambda_1SC[365] = 0.49644539915075914;
  params.lambda_1SC[366] = -0.3964610484283231;
  params.lambda_1SC[367] = -1.2070063596450797;
  params.lambda_1SC[368] = -1.1009836657191028;
  params.lambda_1SC[369] = 0.9623725771782534;
  params.lambda_1SC[370] = -1.1582253709307806;
  params.lambda_1SC[371] = -1.2580585493348209;
  params.lambda_1SC[372] = -0.8509947608983248;
  params.lambda_1SC[373] = 1.9561937247097654;
  params.lambda_1SC[374] = 1.9657515960355676;
  params.lambda_1SC[375] = -0.08071916045506011;
  params.lambda_1SC[376] = -1.1279544365648775;
  params.lambda_1SC[377] = -0.22637881027549733;
  params.lambda_1SC[378] = 1.531730731854938;
  params.lambda_1SC[379] = 0.31862671186881597;
  params.lambda_1SC[380] = 1.602374423783099;
  params.lambda_1SC[381] = 1.5213334094065325;
  params.lambda_1SC[382] = 1.2566660352704577;
  params.lambda_1SC[383] = 1.5357418564739316;
  params.lambda_1SC[384] = 0.05154245184013817;
  params.lambda_1SC[385] = -1.324753811510325;
  params.lambda_1SC[386] = 1.1179260556284074;
  params.lambda_1SC[387] = 1.9694585876958852;
  params.lambda_1SC[388] = 0.6491674951260187;
  params.lambda_1SC[389] = -1.0114223834972633;
  params.lambda_1SC[390] = 1.7124968047182398;
  params.lambda_1SC[391] = 1.6948655432277633;
  params.lambda_1SC[392] = 1.696354213900042;
  params.lambda_1SC[393] = 0.6519303632129039;
  params.lambda_1SC[394] = 1.7532595015771681;
  params.lambda_1SC[395] = -1.0689277585787904;
  params.lambda_1SC[396] = 0.2999711052481464;
  params.lambda_1SC[397] = 1.9736064630176338;
  params.lambda_1SC[398] = 1.979501700160962;
  params.lambda_1SC[399] = 1.2039783558004364;
  params.lambda_1SC[400] = 0.3479105891743006;
  params.lambda_1SC[401] = 0.5106971584398714;
  params.lambda_1SC[402] = -0.4307115965384396;
  params.lambda_1SC[403] = -1.1091317585098595;
  params.lambda_1SC[404] = 1.5860176306348546;
  params.lambda_1SC[405] = 0.34032204461048954;
  params.lambda_1SC[406] = -0.5646458063571318;
  params.lambda_1SC[407] = -0.6907891858235171;
  params.lambda_1SC[408] = -1.3101847145724244;
  params.lambda_1SC[409] = -1.2598104838651492;
  params.lambda_1SC[410] = 0.6705825076020981;
  params.lambda_1SC[411] = -0.8546192381265021;
  params.lambda_1SC[412] = -0.7325468788548601;
  params.lambda_1SC[413] = 1.3803554562999998;
  params.lambda_1SC[414] = -1.6256843580507114;
  params.lambda_1SC[415] = 1.7965441279272834;
  params.lambda_1SC[416] = -0.3399339055070154;
  params.lambda_1SC[417] = 1.3097965458493785;
  params.lambda_1SC[418] = 1.2956610507014727;
  params.lambda_1SC[419] = 0.7300739114640376;
  params.lambda_1SC[420] = -1.7888378883354057;
  params.lambda_1SC[421] = 0.6188282208921394;
  params.lambda_1SC[422] = -0.2137568032106767;
  params.lambda_1SC[423] = -0.7582359932208433;
  params.lambda_1SC[424] = 0.07552442078430044;
  params.lambda_1SC[425] = 0.6582807647803639;
  params.lambda_1SC[426] = -0.4395825380378424;
  params.lambda_1SC[427] = -0.36462729392519533;
  params.lambda_1SC[428] = -1.138369869414038;
  params.lambda_1SC[429] = -1.3340903092039142;
  params.lambda_1SC[430] = 0.49102820219859833;
  params.lambda_1SC[431] = 1.898898721456833;
  params.lambda_1SC[432] = -1.6076018815980047;
  params.lambda_1SC[433] = 0.6277675334260895;
  params.lambda_1SC[434] = 0.3987769566194541;
  params.lambda_1SC[435] = 0.9905862265950569;
  params.lambda_1SC[436] = -0.7234378898375802;
  params.lambda_1SC[437] = -0.9535061418892306;
  params.lambda_1SC[438] = -0.8515323336274423;
  params.lambda_1SC[439] = 1.351262314290167;
  params.lambda_1SC[440] = 1.395297466152738;
  params.lambda_1SC[441] = 1.2510630987219682;
  params.lambda_1SC[442] = -0.8855777049876719;
  params.lambda_1SC[443] = -0.4886371844463899;
  params.lambda_1SC[444] = -0.8389918878672464;
  params.lambda_1SC[445] = 0.7454772770396576;
  params.lambda_1SC[446] = -1.0595538407433325;
  params.lambda_1SC[447] = 0.9721423047498363;
  params.lambda_1SC[448] = 1.1717797996750612;
  params.lambda_1SC[449] = 0.13341997483211987;
  params.lambda_1SC[450] = -1.3019613911703187;
  params.lambda_1SC[451] = -0.9498985431049976;
  params.lambda_1SC[452] = 0.26884959101944306;
  params.lambda_1SC[453] = 0.9405861503538735;
  params.lambda_1SC[454] = -1.196653066182705;
  params.lambda_1SC[455] = 0.5928416070331406;
  params.lambda_1SC[456] = 1.6343988763781745;
  params.lambda_1SC[457] = 0.9289064189489133;
  params.lambda_1SC[458] = 1.553667698452747;
  params.lambda_1SC[459] = 1.835102931426607;
  params.lambda_1SC[460] = -0.902070249980941;
  params.lambda_1SC[461] = -1.1678394203688107;
  params.lambda_1SC[462] = 1.6024713591344644;
  params.lambda_1SC[463] = -0.44192252766618223;
  params.lambda_1SC[464] = -1.5683492873332128;
  params.lambda_1SC[465] = 0.4939387908625674;
  params.lambda_1SC[466] = -0.08586482539630502;
  params.lambda_1SC[467] = -1.2318581361203762;
  params.lambda_1SC[468] = 1.664718902921094;
  params.lambda_1SC[469] = -0.08074373943701518;
  params.lambda_1SC[470] = 1.9647077429046287;
  params.lambda_1SC[471] = -0.8747876532172518;
  params.lambda_1SC[472] = -1.761476745360381;
  params.lambda_1SC[473] = 1.6305706546063696;
  params.lambda_1SC[474] = -0.8150291080939795;
  params.lambda_1SC[475] = 1.695214585600739;
  params.lambda_1SC[476] = -0.06154179073234811;
  params.lambda_1SC[477] = 0.773371010154682;
  params.lambda_1SC[478] = -1.5350540911294388;
  params.lambda_1SC[479] = -0.9415755162983497;
  params.lambda_1SC[480] = 0.7155316230529167;
  params.lambda_1SC[481] = 1.8091069849498709;
  params.lambda_1SC[482] = 1.280009351521683;
  params.lambda_1SC[483] = 1.618102617058292;
  params.lambda_1SC[484] = 1.098765036733019;
  params.lambda_1SC[485] = 1.8880892752962968;
  params.lambda_1SC[486] = -1.0477239549841721;
  params.lambda_1SC[487] = -1.5104935754779705;
  params.lambda_1SC[488] = -1.0379189379087528;
  params.lambda_1SC[489] = 0.32032752790065455;
  params.lambda_1SC[490] = -0.5172073222873852;
  params.lambda_1SC[491] = -1.8593384370605137;
  params.lambda_1SC[492] = -1.501224960877439;
  params.lambda_1SC[493] = -0.8752575006728813;
  params.lambda_1SC[494] = -0.6583101718910673;
  params.lambda_1SC[495] = -0.7297827107849395;
  params.lambda_1SC[496] = 0.26204251997878325;
  params.lambda_1SC[497] = 0.621736081396008;
  params.lambda_1SC[498] = 1.3792584896597941;
  params.lambda_1SC[499] = -0.48040661895150727;
  params.gamma[0] = 1.8759900709408046;
  params.A[0] = -0.4757318797181278;
  params.B[0] = -0.30661897767384616;
  params.lambda_3[0] = -0.34295343044091897;
  params.lambda_4[0] = 1.0586466044456668;
  params.rho[0] = 1.4755818267568077;
  params.Pg_N_init[0] = -0.31061806296636574;
  params.Pg_N_avg[0] = -0.6294403219130329;
  params.ug_N[0] = -1.6701576009286305;
  params.Vg_N_avg[0] = -0.9847849033565876;
  params.Thetag_N_avg[0] = -1.6978969797336316;
  params.vg_N[0] = -0.6696178462723301;
  params.PgMin[0] = 0.3787827614029411;
  params.PgMax[0] = 1.291613023537978;
  params.RgMin[0] = -1.4925354476392583;
  params.selectZero[0] = -0.04594985808100738;
  params.PgNext[0] = 0.9650762635535861;
  params.RgMax[0] = 0.6891187182955274;
}
