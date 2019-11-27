/* Produced by CVXGEN, 2019-02-10 21:04:07 -0500.  */
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
  lhs[0] = rhs[0]*(2*(params.c2[0]+work.frac_479733190656+work.frac_957443481600+work.frac_121674190848));
  lhs[1] = rhs[1]*(2*work.frac_479733190656);
  lhs[2] = rhs[2]*(2*work.frac_121674190848);
}
void fillq(void) {
  work.q[0] = params.c1[0]-2*params.PgNu[0]*work.frac_479733190656-2*params.PgAPPSC[0]*work.frac_957443481600+params.BSC[0]*params.gammaSC[0]-params.lambda_2SC[0]+params.gamma[0]*params.B[0]-params.lambda_4[0]+2*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*work.frac_121674190848;
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
  work.frac_957443481600 = params.betaSC[0];
  work.frac_957443481600 /= 2;
  work.frac_121674190848 = params.rho[0];
  work.frac_121674190848 /= 2;
  work.quad_399913824256[0] = params.PgPrevNu[0]*work.frac_479733190656*params.PgPrevNu[0];
  work.quad_608032620544[0] = params.PgNu[0]*work.frac_479733190656*params.PgNu[0];
  work.quad_583541456896[0] = params.PgAPPSC[0]*work.frac_957443481600*params.PgAPPSC[0];
  work.quad_155015135232[0] = (-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0])*work.frac_121674190848*(-params.Pg_N_init[0]+params.Pg_N_avg[0]+params.ug_N[0]);
  work.quad_497818112000[0] = (-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0])*work.frac_121674190848*(-params.Vg_N_avg[0]-params.Thetag_N_avg[0]+params.vg_N[0]);
}
