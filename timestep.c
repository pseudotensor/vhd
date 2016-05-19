#include "global.h"
#include "defs.h"

void timestep(void)
{
  FTYPE dtother;
  int i, j, k, l;
  FTYPE idt2[NUMDTCHECKS + 1];
  int ks[NUMDTCHECKS + 1], js[NUMDTCHECKS + 1], is[NUMDTCHECKS + 1];
  FTYPE dt2inv_max[NUMDTCHECKS + 1];
  int didfail, didfail_full;
  FTYPE finaln;
  static int firsttime = 1;
  static FTYPE ttimestep = 0;
  static FTYPE dtlast;
  // for slow idtcreate change
  FTYPE bxa, bya, dv;
  FTYPE rho, u;
  FTYPE odx1, odx2, ods, odl;
  FTYPE valphen, valphen2, cs2;
  int reall, viscl, nonvl;
  FTYPE ftemp;
  FTYPE vel1, vel2;
  int gosub, gosup;
  FTYPE dt2invl[3];
  static FTYPE dtotherlowest;
  static int laststep;



  if (visc_real == 1) {
    nu_compute();
  }


  for (l = 2; l <= NUMDTCHECKS; l++) {
    dt2inv_max[l] = 0.;
    ks[l] = js[l] = is[l] = 0;
  }
  dtlast = dt;

  didfail = 0;
  didfail_full = 0;
  if (firsttime == 1) {
    ttimestep = t - 1.E-12;
    laststep = 0;
  }

  LOOP {
#if(TS0CHECK)
    if (s[2][k][j][i] < 0) {	// actually detects nan too
      sprintf(tempc1, "%3f", s[2][k][j][i]);
      if (tempc1[0] == 'n') {
	fprintf(fail_file,
		"nan internal energy density error: k: %3d j: %3d i: %3d u: %15.10g \n",
		k, j, i, s[2][k][j][i]);
      } else if (s[2][k][j][i] < 0)
	fprintf(fail_file,
		"negative internal energy density error: k: %3d j: %3d i: %3d en: %15.10g \n",
		k, j, i, s[2][k][j][i]);
      didfail = 1;
    }
    if (s[1][k][j][i] < 0) {
      sprintf(tempc1, "%3f", s[1][k][j][i]);
      if (tempc1[0] == 'n') {
	fprintf(fail_file,
		"nan mass density error: k: %3d j: %3d i: %3d rho: %15.10g \n",
		k, j, i, s[1][k][j][i]);
      } else if (s[1][k][j][i] < 0)
	fprintf(fail_file,
		"negative mass density error: k: %3d j: %3d i: %3d rho: %15.10g \n",
		k, j, i, s[1][k][j][i]);
      didfail = 1;
    }
#endif



    // not inlining this function for some reason, so slow in loop, 
#include "timestep1.h"

    for (l = 2; l <= NUMDTCHECKS; l++) {
      ftemp = idt2[l];
      if (ftemp > dt2inv_max[l]) {
	dt2inv_max[l] = ftemp;
	ks[l] = k;
	js[l] = j;
	is[l] = i;
      }
#if(CHECKDTLOW==1)
      if (ftemp > SQIDTLOWEST) {
	timecheck(-l, idt2, k, j, i, 0);
	didfail = 1;
	fflush(fail_file);
      }
#endif
    }



  }				// end loop over domain


#if(CHECKDTLOW==1)
  // check if any cpu has failure
  if (numprocs > 1) {
  } else {
    didfail_full = didfail;
  }
  if (didfail_full) {
    if (myid <= 0) {
      fprintf(log_file, "timestep failure\n");
    }
    if (DOGENDIAG) {
      diag(2);
    }
    myexit(5);
  }
#endif



  // find lowest constrainer on dt
  reall = 2;
  for (l = 3; l <= NUMDTCHECKS; l++) {
    if (dt2inv_max[l] > dt2inv_max[reall]) {
      reall = l;
    }
  }

#if(DODTDIAG)
  // do check up on dominates of timestep for each type
  if (t > ttimestep) {		// per cpu pure dt data
    for (l = 2; l <= NUMDTCHECKS; l++) {
      timecheck(l, idt2, ks[l], js[l], is[l], reall);
    }
    fflush(logdt_file);
  }
#endif

#if(DOTSTEPDIAG)
  if (t > ttimestep) {
    // output timescales (create own DTtimescale later)
    timescale();
  }
#endif

#if((DODTDIAG)||(DOTSTEPDIAG))
  ttimestep = t + DTtimestep;
#endif


  // find lowest constrainer on dt due to visc
  if (dt2inv_max[8] > dt2inv_max[9]) {
    viscl = 8;
  } else
    viscl = 9;

  // find lowest constrainer on dt of non-viscosity type (next
  // highest
  // dt^2)
  nonvl = 2;
  for (l = 3; l <= NUMDTCHECKS; l++) {
    if (l == 8)
      l = 10;			// skip viscosity
    if (dt2inv_max[l] > dt2inv_max[nonvl]) {
      nonvl = l;
    }
  }
  // communicate the lowest dt values to all cpus
  if (numprocs > 1) {
  } else {
    dt2invl[0] = dt2inv_max[reall];
    dt2invl[1] = dt2inv_max[viscl];
    dt2invl[2] = dt2inv_max[nonvl];
  }
  
  dt = 1.0 / sqrt(dt2invl[0]);	// normal case of no subcycling

  if (analoutput == 6) {
    // for checking visc code
    dt = pow(invcour2 * alpha_real / (dx[2][1][0] * dx[2][1][0]), -1.0);
    ftemp =
	pow(invcour2 * alpha_real /
	    (x[2][1][0] * x[2][1][0] * dx[2][2][N2 / 2] *
	     dx[2][2][N2 / 2]), -1.0);
    if (ftemp < dt)
      ftemp = dt;
  }
  // because first time step is bad if e.g. viscosity on and no v at
  // first
  if (firsttime == 1) {		// tweak for given problem so starts
    // out
    // good
    if (dt > 1.E-5)
      dt = 1.E-5;
    firsttime = 0;
  }
  // don't increase timestep by too much on subcycle or normal cycle.
  // don't check if supercycle since need to force non-visc back to
  // visc 
  // time.
  if (subcyclen >= 0) {
    if (dt > 1.3 * dtlast)
      dt = 1.3 * dtlast;
  }

  /* don't step beyond end of run */
  if (t + dt >= tf) {
    // last timestep
    laststep = 1;
    if (subcyclen == 1) {
      dt = tf - t;
      reallaststep = 1;
    }
    if (subcyclen == -1) {
      fprintf(fail_file,
	      "shouldn't be here at end of run at super cycle\n");
      myexit(1);
    }
    if (subcyclen >= 2) {
      // just end subcycle and let next timestep() figure final
      // dt,
      // never subcycling again
      dt = (t - tscyclefrom);
      subcyclen = -1;
      tscycleto = t;
      t = tscyclefrom;
      nthsubcycle = 0;
      reallaststep = 0;
    }
    // make sure don't get pathological case of dt=0 on last step
    if (dt < SSMALL) {
      reallaststep = 1;
      laststep = 1;
      dt = SSMALL;
    }
  }
}



void idtcreate(FTYPE * idt2, int k, int j, int i)
{
  FTYPE bxa, bya, dv;
  FTYPE rho, u;
  FTYPE odx1, odx2, ods, odl;
  FTYPE valphen, valphen2, cs2;
  FTYPE ftemp;
  FTYPE vel1, vel2;

#include "timestep1.h"
}


#define NUMPATHS 5

void timescale(void)
{
  int i, j, k, l, m;
  FTYPE bxa, bya, dv;
  FTYPE rho, u;
  FTYPE odx1, odx2, ods, odl;
  FTYPE valphen, valphen2, cs2;
  FTYPE ftemp;
  FTYPE vel1, vel2;
  FTYPE idt2[NUMDTCHECKS + 1];
  static int firsttime = 1;
  static FTYPE paths[NUMPATHS];
  static FTYPE limits[2];	// inner and outer x1 limits
  FTYPE timescales[NUMDTCHECKS + 1][NUMPATHS];	// 5 time scale
  // paths
  static char filename[MAXFILENAME];
  char temps[50];
  static FILE *timescale_file;
  static int startm, lastm;
  int gotfirst, gotlast;

  if (firsttime == 1) {

    // limits on x1-range
    limits[0] = L[1][1];
    limits[1] = R0;		// torus center or injection center

    // locations on x2
    paths[0] = M_PI * 0.5 - M_PI * 0.5 * 2.0 / 3.0;
    paths[1] = M_PI * 0.5 - M_PI * 0.5 / 3.0;
    paths[2] = M_PI * 0.5;
    paths[3] = M_PI * 0.5 + M_PI * 0.5 / 3.0;
    paths[4] = M_PI * 0.5 + M_PI * 0.5 * 2.0 / 3.0;

    // determine what cpu gets what trajectories
    gotfirst = -1;
    gotlast = -1;
    startm = 0;
    lastm = 4;
    k = 0;
    for (m = 0; m <= NUMPATHS - 1; m++) {
      for (j = 0; j < N2; j++) {

	// 0.6 might get 2 trajectories, but so close that's ok
	// anyways since overwritten by 2nd one
	if (fabs(x[2][2][j] - paths[m]) < 0.6 * dx[1][2][j]) {	// if 
	  // 
	  // within 
	  // zone

	  if (gotfirst < 0)
	    gotfirst = m;	// gets which m is first m for
	  // this cpu
	  gotlast = m;		// gets last gotten

	}			// end over this trajectory if found

      }				// end over x2-dir seeking trajectories
    }				// end over all trajectories

    // now determine which cpus should do what trajectories
    if (numprocs > 1) {
      lastm = gotlast;
      if (myid == 0) {
	startm = gotfirst;
      }
    } else {
      startm = gotfirst;
      lastm = gotlast;
    }
    fprintf(log_file, "startm: %d lastm: %d\n", startm, lastm);
    fflush(log_file);
    sprintf(filename, "%s0_timescales%s", DATADIR, DATEXT);


    // write file header
    if (myid <= 0) {
      if ((timescale_file = fopen(filename, WRITETYPE)) == NULL) {	// naively 
									// 
	// 
	// append 
	// if 
	// appendold=1
	fprintf(fail_file, "timescales: Cannot open: %s\n", filename);
	myexit(1);
      }
      if (appendold == 0) {
	fprintf(timescale_file, "#%10s\n%10d %10d\n", "TSVER",
		TSVER, TSTYPE);
	fprintf(timescale_file, "#%15s", "time");
	for (m = startm; m <= lastm; m++) {
	  for (l = 2; l <= NUMDTCHECKS; l++) {
	    sprintf(temps, "p%1d-c%2d", m, l);
	    fprintf(timescale_file, " %15s", temps);
	  }
	}
      }
      fclose(timescale_file);
    }
  }				// endif firsttime==1


  // now find different timescales along different trajectories
  // reset timescales
  for (l = 2; l <= NUMDTCHECKS; l++) {
    for (m = 0; m <= NUMPATHS - 1; m++) {
      timescales[l][m] = 0;
    }
  }
  k = 0;
  for (m = startm; m <= lastm; m++) {
    for (j = 0; j < N2; j++) {

      if (fabs(x[2][2][j] - paths[m]) < 0.6 * dx[1][2][j]) {	// if 
	// 
	// within 
	// zone

	for (i = 0; i < N1; i++) {

	  if ((x[2][1][i] > limits[0]) && (x[2][1][i] < limits[1])) {	// if 
									// 
	    // 
	    // within 
	    // bounds 
	    // of 
	    // trajectory

#include "timestep2.h"

	    for (l = 2; l <= NUMDTCHECKS; l++) {
	      timescales[l][m] += 1.0 / sqrt(SSMALL + idt2[l]);	// add 
	      // 
	      // up 
	      // dt 
	      // for 
	      // this 
	      // zone
	    }

	  }			// end if within bounds of trajectory

	}			// end loop over this trajectory

      }				// end over this trajectory if found

    }				// end over x2-dir seeking trajectories

  }				// end over all trajectories

  // now write timescales to file
  for (i = 0; i < numprocs; i++) {
    if (myid == i) {
      if ((timescale_file = fopen(filename, "at")) == NULL) {
	fprintf(fail_file, "timescales: Cannot open: %s\n", filename);
	myexit(1);
      }
      if (i == 0) {
	fprintf(timescale_file, " %15.10g", t);
      }
      for (m = startm; m <= lastm; m++) {
	for (l = 2; l <= NUMDTCHECKS; l++) {
	  fprintf(timescale_file, " %15.10g", timescales[l][m]);
	}
      }
      if (i == numprocs - 1) {
	fprintf(timescale_file, "\n");
      }
      fclose(timescale_file);
    }
  }

  firsttime = 0;

}

// failmode==0 fail
// failmode=# which to check
void timecheck(int failmode, FTYPE * idt2, int k, int j, int i,
	       int reall)
{
  FTYPE ftemp;
  int wfail, l;
  FTYPE dtinv;
  FTYPE valphen, bxa, bya;
  FILE *out;
  static int firsttime = 1;
  FTYPE vel1, vel2;
  FTYPE cs2;

  if (failmode > 0) {
    idtcreate(idt2, k, j, i);
    dtinv = sqrt(idt2[failmode]);
  } else {
    // no need to create idt's if failmode<0 since just created it
    // in
    // loop
    dtinv = sqrt(idt2[-failmode]);
    failmode = 0;		// reactivate fail mode
  }


  if (wgam)
    cs2 = gam * (gam - 1.) * s[2][k][j][i] / s[1][k][j][i];
  else
    cs2 = cs * cs;
  cs = sqrt(cs2);		// need cs now
  vel1 = (0.5 * (v[1][1][k][j][i] + v[1][1][k][j][i + 1]) - vg[1]);
  vel2 = (0.5 * (v[1][2][k][j][i] + v[1][2][k][j + 1][i]) - vg[2]);
  if (vel1 >= 0) {
    vel1 = vel1 + cs;
  } else {
    vel1 = vel1 - cs;
  }
  if (vel2 >= 0) {
    vel2 = vel2 + cs;
  } else {
    vel2 = vel2 - cs;
  }

  /* alfven velocity */
  bxa = 0.5 * (v[2][1][k][j][i] + v[2][1][k][j][i + 1]);
  bya = 0.5 * (v[2][2][k][j][i] + v[2][2][k][j + 1][i]);
  valphen =
      sqrt((bxa * bxa + bya * bya +
	    v[2][3][k][j][i] * v[2][3][k][j][i]) / s[1][k][j][i]);

  // see who failed or see what dominates for checkup
  wfail = 2;
  ftemp = fabs(idt2[wfail]);
  for (l = 3; l <= NUMDTCHECKS; l++) {
    if (fabs(idt2[l]) > ftemp) {
      ftemp = fabs(idt2[l]);
      wfail = l;
    }
  }
  ftemp = 1.0 / (sqrt(ftemp + SSMALL));	// real dt of failure 

  if (wfail > NUMDTCHECKS) {
    fprintf(fail_file,
	    "unexpected failure in timestep.c when checking dt\n");
    myexit(1);
  }


  if ((failmode == 0) || (DODTDIAG == 0)) {
    out = fail_file;
    fprintf(fail_file,
	    "dt has dropped below set lowest threshold.  dt=%15.10g\n",
	    1.0 / dtinv);
  } else
    out = logdt_file;


  if (failmode == 0) {
    fprintf(out, "Time check at t=%15.10g #%d(%d) has dt=%15.10g\n", t,
	    wfail, reall, ftemp);
#if( (COORD==3)&&(SAMPLED>0) )
    fprintf(out, "Velocities at k=%d j=%d i=%d x=%15.10g z=%15.10g\n",
	    k, j, i, x[1][1][i] * sin(x[1][2][j]),
	    x[1][1][i] * cos(x[1][2][j]));
#else
    fprintf(out,
	    "Velocities at  k=%d j=%d i=%d x1=%15.10g x2=%15.10g\n", k,
	    j, i, x[1][1][i], x[1][2][j]);
#endif

    fprintf(out,
	    "#  Vel Checked   %15s : %15s\n"
	    "-------------------------------------\n"
	    "1: sound speed : %15.10g : %15.10g\n"
	    "2: x1-vel+-cs  : %15.10g : %15.10g\n"
	    "3: x2-vel+-cs  : %15.10g : %15.10g\n"
	    "4: Alven vel   : %15.10g : %15.10g\n"
	    "5: Linear visc : %15.10g : %15.10g\n"
	    "6: visc(x1)    : %15.10g : %15.10g\n"
	    "7: visc(x2)    : %15.10g : %15.10g\n"
	    "8: rvisc(x1):  : %15.10g : %15.10g\n"
	    "9: rvisc(x2):  : %15.10g : %15.10g\n"
	    "10: resist:    : %15.10g : %15.10g\n\n", "dt", "v||dv",
	    SSMALL, cs, 1. / (fabs(sqrt(idt2[2])) + SSMALL), vel1,
	    1. / fabs(sqrt(idt2[3]) + SSMALL), vel2,
	    1. / fabs(sqrt(idt2[4]) + SSMALL), valphen,
	    1. / fabs(sqrt(idt2[5]) + SSMALL), invcour2 * nu_l * cs,
	    1. / fabs(sqrt(idt2[6]) + SSMALL),
	    v[1][1][k][j][i + 1] - v[1][1][k][j][i],
	    1. / fabs(sqrt(idt2[7]) + SSMALL),
	    v[1][2][k][j + 1][i] - v[1][2][k][j][i],
	    1. / fabs(sqrt(idt2[8]) + SSMALL),
	    v[1][1][k][j][i + 1] - v[1][1][k][j][i],
	    1. / fabs(sqrt(idt2[9]) + SSMALL),
	    v[1][2][k][j + 1][i] - v[1][2][k][j][i],
	    1. / fabs(sqrt(idt2[10]) + SSMALL), invcour2 * resist * cs);

  } else {
    if ((firsttime == 1) && (appendold == 0)) {
      fprintf(out, "#%10s\n%10d %10d\n", "LOGDTVER", LOGDTVER,
	      LOGDTTYPE);
      fprintf(out,
	      "#%15s %15s %5s %5s %5s %5s %5s %5s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
	      "t", "dt", "r", "w", "l", "k", "j", "i", "x1", "x2",
	      "cs_dt", "cs_v", "x1v_dt", "x1v_(v+-cs)", "x2v_dt",
	      "x2v_(v+-cs)", "Bv_dt", "Bv_v", "lv_dt", "lv_dv",
	      "vx1_dt", "vx1_dv", "vx2_dt", "vx2_dv", "rvx1_dt",
	      "rvx1_nu", "rvx2_dt", "rvx2_nu", "resist_dt", "resist_v");
    }

    fprintf(out, " %15.10g %15.10g %5d %5d %5d", t,
	    1. / (fabs(sqrt(idt2[failmode])) + SSMALL), reall, wfail,
	    failmode);
#if( (COORD==3)&&(SAMPLED>0) )
    fprintf(out, " %5d %5d %5d %15.10g %15.10g", k, j, i,
	    x[1][1][i] * sin(x[1][2][j]), x[1][1][i] * cos(x[1][2][j]));
#else
    fprintf(out, " %5d %5d %5d %15.10g %15.10g", k, j, i, x[1][1][i],
	    x[1][2][j]);
#endif
    fprintf(out,
	    " %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
	    SSMALL, cs, 1. / (fabs(sqrt(idt2[2])) + SSMALL), vel1,
	    1. / fabs(sqrt(idt2[3]) + SSMALL), vel2,
	    1. / fabs(sqrt(idt2[4]) + SSMALL), valphen,
	    1. / fabs(sqrt(idt2[5]) + SSMALL), invcour2 * nu_l * cs,
	    1. / fabs(sqrt(idt2[6]) + SSMALL),
	    v[1][1][k][j][i + 1] - v[1][1][k][j][i],
	    1. / fabs(sqrt(idt2[7]) + SSMALL),
	    v[1][2][k][j + 1][i] - v[1][2][k][j][i],
	    1. / fabs(sqrt(idt2[8]) + SSMALL), nu_real[k][j][i],
	    1. / fabs(sqrt(idt2[9]) + SSMALL), nu_real[k][j][i],
	    1. / fabs(sqrt(idt2[10]) + SSMALL), invcour2 * resist * cs);

  }


  firsttime = 0;
}
