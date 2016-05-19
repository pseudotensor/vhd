#include "global.h"
#include "defs.h"

#include <signal.h>
#include <time.h>
int main(int argc, char *argv[], char *envp[]
    )
{
  int dostepout=0;
  int j, i;
  SFTYPE tlasttime;
  SFTYPE tnext, tnextavg;
  static long subnstep = 0, supnstep = 0, nornstep = 0;
  int error = 0;
  FILE *gogo_file;
  FILE *fileit=NULL;
  char stemp[50];
  int flooroutc;
  SFTYPE tfloor;
  char goch;
  int goend, goend_full=0;
  int avgdone = 0;
  time_t timestart, timestop;
  time_t gtimestart, gtimestop;
  SFTYPE walltime = 0, walltimelocal = 0;
  static int logcount = 0;
  static int avgcount = 0;
  static int floorcount = 0;
  int numzones=0;
  FILE *cpuout;
  SFTYPE comptstart;
  int imagescaleflag = 0;

  // ///////////////// INITIALIZATION

  // uncomment if get SIGFPE 8 error at run with floats and no result
  // problem
  // signal(8,SIG_IGN);

  myid = 0;			// defines single process run
  numprocs = 1;



  if (!(cpuout = fopen("numcpus.txt", "wt"))) {
    printf("Can't open numcpus.txt\n");
    myexit(1);
  } else {
    fprintf(cpuout, "%d\n", numprocs);
    fclose(cpuout);
  }


  // start interps fresh no matter what.
  globalinterpmod = 1;
  /* Initialize ALL */
  init_genfiles(0);		// setup general log files
  if (myid <= 0) {
    fprintf(logfull_file, "#Start Initialization.\n");
    fflush(logfull_file);
  }
  fprintf(stderr, "proc: %02d Start Initialization.\n", myid);
  if ((error = init(argc, argv, envp)) > 0) {
    fprintf(fail_file, "\nError: initialize_par: %d\n", error);
    myexit(1);
  }

  floorcount = 0;
  logcount = 0;
  nstep = 0;
  flooroutc = floor_start;
  avgdone = 0;
  avgcount = 0;
  logcount = 0;
  tnext = tstart;
  tnext = tstart + (FTYPE) (logcount) * DTl;
  tnextavg = tavgi;
  tlasttime = 0;
  comptstart = t;

  fprintf(stderr, "proc: %02d   End Initialization.\n", myid);
  if (myid <= 0) {
    fprintf(logfull_file, "#End Initialization.\n");
    fflush(logfull_file);
  }
  // write version header info for perf and step logs
  if ((myid <= 0) && (appendold == 0)) {
    fprintf(logstep_file, "#%10s\n%10d %10d\n", "STEPVER", STEPVER,
	    STEPTYPE);
    fprintf(logperf_file, "#%10s\n%10d %10d\n", "PERFVER", PERFVER,
	    PERFTYPE);
  }
  /* do initial diagnostics */
#if(DOLOSSDIAG)
  diag(-1);			// must come before normal diag
#endif
#if(DOGENDIAG)
  diag(0);
#endif
  logcount++;
  tnext = tstart + (FTYPE) (logcount) * DTl;

#if(DOFLOORDIAG>=1)
  tfloor = t - 1.E-12;
  for (i = 0; i < NUMFLOOROUT; i++) {
    for (j = 1; j <= NUMFLOORVAR; j++) {
      floorcnt[i][j] = 0;
    }
  }
#endif



  // //////////////// BEGIN COMPUTATION

  fprintf(stderr, "Starting Computation on proc: %02d . . .\n", myid);
  if (myid <= 0) {
    fprintf(logfull_file, "#Starting Computation\n");
    fflush(logfull_file);
  }
  time(&timestart);
  time(&gtimestart);


  goend = 0;
  subcyclen = 1;		// start fresh
  reallaststep = 0;
  if (myid <= 0) {
    fprintf(logstep_file, "#");
  }
  // while(t < tf) {
  while (reallaststep == 0) {

    /* find timestep */
    timestep();

    /* step variables forward in time */
    stepvar();

    tdep_compute();		// find new time dep stuff
    if (DOSPDIAG) {
      sp_compute();		// check on sonic point
    }
    analsolve(0);		// find new analytic values


    nstep++;			// any type of step
    if (subcyclen == -1)
      supnstep++;
    if (subcyclen == 1)
      nornstep++;

    /* Perform required diagnostics that need to be done per timestep */
    // only do diags if just supercycled or normal cycle
    if (subcyclen <= 1) {
#if(DOLOSSDIAG)
      diag(-1);
#endif

      // perform general diagnostics 
#if(DOGENDIAG)
      diag(1);
#endif




      // floor diagnostics
#if(DOFLOORDIAG>=1)
      if (t > tfloor) {
	flooroutc++;
	floorcount++;
	tfloor = tstart + (FTYPE) (floorcount) * DTfloor;

	fprintf(logfl_file, "floorcnt: t=%15.10g\n", t);
	for (i = 0; i < NUMFLOOROUT; i++) {
	  for (j = 1; j <= NUMFLOORVAR; j++) {
	    fprintf(logfl_file, "%12d ", floorcnt[i][j]);
	  }
	  fprintf(logfl_file, "\n");
	}
	for (j = 1; j <= NUMFLOORVAR; j++) {
	  fprintf(logfl_file, "lowest %d %15.10g\n",
		  wherelowest[j], floorlowest[j]);
	}


      }
#if(FLUSHFAILDT)
      fflush(logfl_file);
#endif
#endif

    }				// end if normal or super cycle

    // check if time to output step/time/dt info
    // setup so can plot in sm
    if ((!(nstep % (NDTCCHECK / 100)))) {
      if (myid <= 0) {
	fprintf(logstep_file, ".");
	fflush(logstep_file);
      }
    }
    if ((!(nstep % NDTCCHECK)) || (t >= tf - 1.0E-7)
	|| (t <= tstart + 1.0E-7)) {
      if (myid <= 0) {
	fprintf(logstep_file, "\n");
      }
      for (i = 1; i <= 1; i++) {	// ==0 not done since really
	// not
	// needed
	if (i == 0) {
	  fileit = log_file;
	  dostepout = 1;
	} else if (i == 1) {
	  fileit = logstep_file;
	  if (myid <= 0)
	    dostepout = 1;
	  else
	    dostepout = 0;
	}
	if (dostepout) {
	  fprintf(fileit, "#step,norm,sup,sub,t,dt,upto,i,N\n"
		  "%10ld %10ld %10ld %10ld %15.10g %15.10g ",
		  nstep, nornstep, supnstep, subnstep, t, dt);
	  if (subcyclen > 1)
	    fprintf(fileit, "%15.10g %5d %5d\n", tscycleto,
		    nthsubcycle, subcyclen);
	  else if (subcyclen == 1)
	    fprintf(fileit, "%15.10g 0 0\n", t + dt);
	  else if (subcyclen == -1)
	    fprintf(fileit, "%15.10g 0 0\n", tscycleto);
	  fflush(fileit);
	}
      }
      if (myid <= 0) {
	fprintf(logstep_file, "#");
      }
      // check image scale too
      if (myid <= 0) {
	if (imagescaleflag == 0) {
	  if (t >= timagescale) {
	    fprintf(logfull_file,
		    "#t>timagescale hit: t: %21.15g timagescale: %21.15g\n",
		    t, timagescale);
	  }
	  imagescaleflag = 1;
	}
      }
    }
#if(FLUSHFAILDT)
    fflush(fail_file);
    // fflush(log_file);
#endif

    // check if user wants to stop or not(go.go)
    if (!(nstep % NGOCHECK)) {

      sprintf(stemp, "%sgo.go%s", DATADIR, myidtxt);

      if ((gogo_file = fopen(stemp, "rt")) == NULL) {
	fprintf(fail_file, "Could not open go file: %s\n", stemp);
	myexit(1);
      }
      goch = fgetc(gogo_file);
      // fprintf(stderr,"got here: %c\n",goch);
      fflush(stdout);
      if ((goch == 'n') || (goch == 'N')) {
	goend = 1;
	fprintf(stderr, "proc: %02d: go.go called\n", myid);
	if (myid <= 0) {
	  fprintf(logfull_file, "#go.go called\n");
	  fflush(logfull_file);
	}
	fflush(stdout);
      }
      fclose(gogo_file);
      if (numprocs > 1) {
      } else {
	goend_full = goend;
      }
      if (goend_full != goend) {
	if (goend == 1)
	  fprintf(stderr, "proc: %02d: go.go call rejected\n", myid);
	else if (goend == 0)
	  fprintf(stderr, "proc: %02d: go.go call not active\n", myid);
	goend = 0;		// restore original state
      } else {
	if (goend_full == 1) {
	  reallaststep = 1;	// if really all see goend=1,
	  // then 
	  // is really last step, all have
	  // goend=1
	}
      }
    }
    // speed check
    // setup so can plot in sm
    if (!(nstep % NZCCHECK)) {
      time(&gtimestop);
      // running average
      // running average zonecycle rate
      walltime = (SFTYPE) difftime(gtimestop, timestart);
      if (walltime < 1)
	walltime = 1;
      // local zonecycle rate
      walltimelocal = (SFTYPE) difftime(gtimestop, gtimestart);
      if (walltimelocal < 1)
	walltimelocal = 1;
      for (i = 1; i <= 1; i++) {	// don't really want perf for
	// i==0
	if (i == 0) {
	  fileit = log_file;
	  dostepout = 1;
	  numzones = N1 * N2 * N3;
	  strcpy(stemp, "");
	} else if (i == 1) {
	  fileit = logperf_file;
	  if (myid <= 0)
	    dostepout = 1;
	  else
	    dostepout = 0;
	  numzones = totalzones;
	  strcpy(stemp, "all");
	}
	if (dostepout) {
	  fprintf(fileit,
		  "#t, ete, n, wt, zc, tu/hr,  lete, ln, lwt, lzc, ltu/hr\n");
	  fprintf(fileit,
		  "%15.10g %15.10g %10ld %15.10g %10d %10d  %15.10g %5d %15.10g %10d %10d\n",
		  t,
		  ((tf - t + 1.0E-6) / (t - comptstart +
					1.0E-6) * walltime *
		   2.777777E-4)
		  , nstep, walltime * 2.777777E-4,
		  (int) ((FTYPE) (numzones) * (FTYPE) (nstep) /
			 walltime)
		  , (int) ((t - comptstart) / (walltime * 2.777777E-4))

		  ,
		  ((tf - t + 1.0E-6) / (t - tlasttime +
					1.0E-6) * walltimelocal *
		   2.777777E-4)
		  , NZCCHECK, walltimelocal * 2.777777E-4,
		  (int) ((FTYPE) (numzones) *
			 (FTYPE) (NZCCHECK) / walltimelocal)
		  ,
		  (int) ((t -
			  tlasttime) / (walltimelocal * 2.777777E-4))
	      );
	  fflush(fileit);
	}
      }
      time(&gtimestart);
      tlasttime = t;
    }				// end if output speed

  }				// end over all time


  // /////////////// END COMPUTATION




  time(&timestop);

#if(DOGENDIAG)
  /* do final diagnostics */
  diag(2);
#endif
#if(DOAVGDIAG)
  if (avgdone == 0) {		// finish average if didn't actually
    // get
    // past final average time
    diagavg(2);
    avgdone = 1;
  }
#endif

  if (goend == 1) {
    fprintf(log_file, "proc: %02d Go end called(go.go)\n", myid);
    if (myid <= 0) {
      fprintf(logfull_file, "#Go end called(go.go)\n");
      fflush(logfull_file);
    }
  }
  fprintf(log_file, "nstep: %ld\n", nstep);
  fprintf(log_file, "subnstep: %ld\n", subnstep);


  // running average zonecycle rate
  walltime = (FTYPE) difftime(timestop, timestart);
  if (walltime < 1)
    walltime = 1;

  if (myid <= 0) {
    fprintf(logfull_file,
	    "#allproc: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g\n",
	    nstep, walltime * 2.777777E-4,
	    (int) ((FTYPE) (totalzones) * (FTYPE) nstep / walltime),
	    (t - comptstart));
    fprintf(logperf_file,
	    "#done: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g tu/hour: %10.5g\n",
	    nstep, walltime * 2.777777E-4,
	    (int) ((FTYPE) (totalzones) * (FTYPE) nstep / walltime),
	    (t - comptstart),
	    (t - comptstart) / (walltime * 2.777777E-4));
    fprintf(logstep_file,
	    "#done: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g\n",
	    nstep, walltime * 2.777777E-4,
	    (int) ((FTYPE) (totalzones) * (FTYPE) nstep / walltime),
	    (t - comptstart));
  }
  // end cut pasted from in loop

#if(DOFLOORDIAG>=1)
  fprintf(logfl_file, "floorcnt: t=%15.10g\n", t);
  for (i = 0; i < NUMFLOOROUT; i++) {
    for (j = 1; j <= NUMFLOORVAR; j++) {
      fprintf(logfl_file, "%12d ", floorcnt[i][j]);
    }
    fprintf(logfl_file, "\n");
  }
  for (j = 1; j <= NUMFLOORVAR; j++) {
    fprintf(logfl_file, "lowest %d %15.10g\n", wherelowest[j],
	    floorlowest[j]);
  }
#endif


  return (myexit(0));
}
