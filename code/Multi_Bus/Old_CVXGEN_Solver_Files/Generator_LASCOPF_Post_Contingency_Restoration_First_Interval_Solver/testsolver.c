/* Produced by CVXGEN, 2016-05-03 19:35:15 -0400.  */
/* CVXGEN is Copyright (C) 2006-2012 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2012 Jacob Mattingley. */
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
  params.PgNu[0] = 1.5717878173906188;
  params.PgNextNu[0] = 1.5851723557337523;
  params.PgNextNu[1] = -1.497658758144655;
  params.PgNextNu[2] = -1.171028487447253;
  params.PgNextNu[3] = -1.7941311867966805;
  params.PgNextNu[4] = -0.23676062539745413;
  params.PgNextNu[5] = -1.8804951564857322;
  params.PgNextNu[6] = -0.17266710242115568;
  params.PgNextNu[7] = 0.596576190459043;
  params.PgNextNu[8] = -0.8860508694080989;
  params.PgNextNu[9] = 0.7050196079205251;
  params.PgNextNu[10] = 0.3634512696654033;
  params.gamma[0] = 0.04796376475433073;
  params.B[0] = 0.23541635196352795;
  params.B[1] = -0.9629902123701384;
  params.B[2] = -0.3395952119597214;
  params.B[3] = -0.865899672914725;
  params.B[4] = 0.7725516732519853;
  params.B[5] = -0.23818512931704205;
  params.B[6] = -1.372529046100147;
  params.B[7] = 0.17859607212737894;
  params.B[8] = 1.1212590580454682;
  params.B[9] = -0.774545870495281;
  params.B[10] = -1.1121684642712744;
  params.D[0] = -0.44811496977740495;
  params.D[1] = 1.7455345994417217;
  params.D[2] = 1.9039816898917352;
  params.D[3] = 0.6895347036512547;
  params.D[4] = 1.6113364341535923;
  params.D[5] = 1.383003485172717;
  params.D[6] = -0.48802383468444344;
  params.D[7] = -1.631131964513103;
  params.D[8] = 0.6136436100941447;
  params.D[9] = 0.2313630495538037;
  params.D[10] = -0.5537409477496875;
  params.lambda_1[0] = -1.0997819806406723;
  params.lambda_1[1] = -0.3739203344950055;
  params.lambda_1[2] = -0.12423900520332376;
  params.lambda_1[3] = -0.923057686995755;
  params.lambda_1[4] = -0.8328289030982696;
  params.lambda_1[5] = -0.16925440270808823;
  params.lambda_1[6] = 1.442135651787706;
  params.lambda_1[7] = 0.34501161787128565;
  params.lambda_1[8] = -0.8660485502711608;
  params.lambda_1[9] = -0.8880899735055947;
  params.lambda_1[10] = -0.1815116979122129;
  params.lambda_2[0] = -1.17835862158005;
  params.lambda_2[1] = -1.1944851558277074;
  params.lambda_2[2] = 0.05614023926976763;
  params.lambda_2[3] = -1.6510825248767813;
  params.lambda_2[4] = -0.06565787059365391;
  params.lambda_2[5] = -0.5512951504486665;
  params.lambda_2[6] = 0.8307464872626844;
  params.lambda_2[7] = 0.9869848924080182;
  params.lambda_2[8] = 0.7643716874230573;
  params.lambda_2[9] = 0.7567216550196565;
  params.lambda_2[10] = -0.5055995034042868;
  params.rho[0] = 1.336269609470535;
  params.ones[0] = -0.6406053441727284;
  params.ones[1] = 0.29117547947550015;
  params.ones[2] = -0.6967713677405021;
  params.ones[3] = -0.21941980294587182;
  params.ones[4] = -1.753884276680243;
  params.ones[5] = -1.0292983112626475;
  params.ones[6] = 1.8864104246942706;
  params.ones[7] = -1.077663182579704;
  params.ones[8] = 0.7659100437893209;
  params.ones[9] = 0.6019074328549583;
  params.ones[10] = 0.8957565577499285;
  params.Pg_N_init[0] = -0.09964555746227477;
  params.Pg_N_init[1] = 0.38665509840745127;
  params.Pg_N_init[2] = -1.7321223042686946;
  params.Pg_N_init[3] = -1.7097514487110663;
  params.Pg_N_init[4] = -1.2040958948116867;
  params.Pg_N_init[5] = -1.3925560119658358;
  params.Pg_N_init[6] = -1.5995826216742213;
  params.Pg_N_init[7] = -1.4828245415645833;
  params.Pg_N_init[8] = 0.21311092723061398;
  params.Pg_N_init[9] = -1.248740700304487;
  params.Pg_N_init[10] = 1.808404972124833;
  params.Pg_N_avg[0] = 0.7264471152297065;
  params.Pg_N_avg[1] = 0.16407869343908477;
  params.Pg_N_avg[2] = 0.8287224032315907;
  params.Pg_N_avg[3] = -0.9444533161899464;
  params.Pg_N_avg[4] = 1.7069027370149112;
  params.Pg_N_avg[5] = 1.3567722311998827;
  params.Pg_N_avg[6] = 0.9052779937121489;
  params.Pg_N_avg[7] = -0.07904017565835986;
  params.Pg_N_avg[8] = 1.3684127435065871;
  params.Pg_N_avg[9] = 0.979009293697437;
  params.Pg_N_avg[10] = 0.6413036255984501;
  params.ug_N[0] = 1.6559010680237511;
  params.ug_N[1] = 0.5346622551502991;
  params.ug_N[2] = -0.5362376605895625;
  params.ug_N[3] = 0.2113782926017822;
  params.ug_N[4] = -1.2144776931994525;
  params.ug_N[5] = -1.2317108144255875;
  params.ug_N[6] = 0.9026784957312834;
  params.ug_N[7] = 1.1397468137245244;
  params.ug_N[8] = 1.8883934547350631;
  params.ug_N[9] = 1.4038856681660068;
  params.ug_N[10] = 0.17437730638329096;
  params.Vg_N_avg[0] = -1.6408365219077408;
  params.Vg_N_avg[1] = -0.04450702153554875;
  params.Vg_N_avg[2] = 1.7117453902485025;
  params.Vg_N_avg[3] = 1.1504727980139053;
  params.Vg_N_avg[4] = -0.05962309578364744;
  params.Vg_N_avg[5] = -0.1788825540764547;
  params.Vg_N_avg[6] = -1.1280569263625857;
  params.Vg_N_avg[7] = -1.2911464767927057;
  params.Vg_N_avg[8] = -1.7055053231225696;
  params.Vg_N_avg[9] = 1.56957275034837;
  params.Vg_N_avg[10] = 0.5607064675962357;
  params.Thetag_N_avg[0] = -1.4266707301147146;
  params.Thetag_N_avg[1] = -0.3434923211351708;
  params.Thetag_N_avg[2] = -1.8035643024085055;
  params.Thetag_N_avg[3] = -1.1625066019105454;
  params.Thetag_N_avg[4] = 0.9228324965161532;
  params.Thetag_N_avg[5] = 0.6044910817663975;
  params.Thetag_N_avg[6] = -0.0840868104920891;
  params.Thetag_N_avg[7] = -0.900877978017443;
  params.Thetag_N_avg[8] = 0.608892500264739;
  params.Thetag_N_avg[9] = 1.8257980452695217;
  params.Thetag_N_avg[10] = -0.25791777529922877;
  params.vg_N[0] = -1.7194699796493191;
  params.vg_N[1] = -1.7690740487081298;
  params.vg_N[2] = -1.6685159248097703;
  params.vg_N[3] = 1.8388287490128845;
  params.vg_N[4] = 0.16304334474597537;
  params.vg_N[5] = 1.3498497306788897;
  params.vg_N[6] = -1.3198658230514613;
  params.vg_N[7] = -0.9586197090843394;
  params.vg_N[8] = 0.7679100474913709;
  params.vg_N[9] = 1.5822813125679343;
  params.vg_N[10] = -0.6372460621593619;
  params.PgMin[0] = 0.12934639598056652;
  params.PgMax[0] = 1.7282393388212876;
  params.RgMin[0] = -0.581744891658952;
  params.RgMax[0] = 1.4821648127991252;
  params.PgPrev[0] = 0.316067309402988;
}
