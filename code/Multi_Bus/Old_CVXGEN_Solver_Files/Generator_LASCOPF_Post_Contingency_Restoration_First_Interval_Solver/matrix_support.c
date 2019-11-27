/* Produced by CVXGEN, 2016-05-03 19:35:14 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
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
void multbymG(double *lhs, double *rhs) {
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
void multbymGT(double *lhs, double *rhs) {
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
void multbyP(double *lhs, double *rhs) {
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
void fillq(void) {
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
void fillh(void) {
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
void fillb(void) {
}
void pre_ops(void) {
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
