/* Produced by CVXGEN, 2019-02-10 21:04:09 -0500.  */
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
  params.betaSC[0] = 0.2511706209276725;
  params.PgAPPSC[0] = -1.171028487447253;
  params.gammaSC[0] = 0.10293440660165976;
  params.BSC[0] = -0.23676062539745413;
  params.lambda_2SC[0] = -1.8804951564857322;
  params.gamma[0] = 0.9136664487894222;
  params.A[0] = 0.596576190459043;
  params.B[0] = -0.8860508694080989;
  params.lambda_3[0] = 0.7050196079205251;
  params.lambda_4[0] = 0.3634512696654033;
  params.rho[0] = 0.04796376475433073;
  params.Pg_N_init[0] = 0.23541635196352795;
  params.Pg_N_avg[0] = -0.9629902123701384;
  params.ug_N[0] = -0.3395952119597214;
  params.Vg_N_avg[0] = -0.865899672914725;
  params.Thetag_N_avg[0] = 0.7725516732519853;
  params.vg_N[0] = -0.23818512931704205;
  params.PgMin[0] = 0.31373547694992654;
  params.PgMax[0] = 1.0892980360636895;
  params.RgMin[0] = -1.560629529022734;
  params.selectZero[0] = -0.774545870495281;
  params.PgNext[0] = 0.4439157678643628;
  params.RgMax[0] = 0.7759425151112975;
}
