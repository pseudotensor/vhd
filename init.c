#include "global.h"
#include "defs.h"

#if(SENSITIVE==1)
#define INPUT1 "%lf %d %d"
#define INPUT1OLD "%lf %d"
#define HEADER3_S2 "%lf %lf %lf %lf %lf %lf\n"
#define HEADER3_S3 "%lf %lf %lf %lf %lf %lf %lf\n"
#define HEADER3_S4 "%lf %lf %lf %lf %lf\n"
#define HEADER3_S5 "%lf %lf\n"
#define HEADER3_S6 "%lf %lf %lf %lf %d %lf\n"
#define HEADER3_S7 "%lf %lf %lf %lf %lf %lf %lf %lf\n"
#define HEADER3_S8 "%lf %lf\n"
#define HEADER3_S9 "%lf %lf %lf\n"
#else
#define INPUT1 "%f %d %d"
#define INPUT1OLD "%f %d"
#define HEADER3_S2 "%f %f %f %f %f %f\n"
#define HEADER3_S3 "%f %f %f %f %f %f %f\n"
#define HEADER3_S4 "%f %f %f %f %f\n"
#define HEADER3_S5 "%f %f\n"
#define HEADER3_S6 "%f %f %f %f %d %f\n"
#define HEADER3_S7 "%f %f %f %f %f %f %f %f\n"
#define HEADER3_S8 "%f %f\n"
#define HEADER3_S9 "%f %f %f\n"
#endif

// 1/2/3-d valid
int init_params(int seed, FTYPE beta, FTYPE nj)
{
  SFTYPE tvisc;
  SFTYPE numvtime;
  /* BEGIN: Assign global parameters */

  runtype = 0;
  directinput = 0;		// only applies if runtype>0

  // 
  // runtype:
  // change numbers at botton of function too!
  // 0: nonreentrant
  // 1: data only reentrant
  // 2: par only reentrant
  // 3: data and par reentrant
  // check for runtime parameter file, and if exists, use it.
  // Think binary, 0=scratch, 1=do from file
  // if 00=0 or no number(argc==1) -> par->scratch data->scratch
  // for rest, argc==2
  // if 01=1 -> par->scratch data->file
  // if 10=2 -> par->file data->scratch
  // if 11=3 -> par->file data->file
  // runtype is the decimal equiv of that binary code
  // 
  // 
  // directinput:
  // 0: user uses .in files for all inputs(good for distinct start
  // from
  // other data set as initial data set) (currently, data:
  // 0_pdump.in.xx 
  // when DUMPSM==1, 0_dump.in.xx when DUMPSM==0)
  // 1: use given numbers for dumps in init.c and normal .par ext for
  // inputs, directly starting where left off when doing runtype>0
  // 2: Total reentrance based upon assumption that want last data
  // outputted to be used (most useful for startup after crash/go.go
  // stop, or continuing after run done but want to continue(extend
  // time)
  // 3: like #1, but specify time instead of dump # (most useful for
  // random startup)

  timereenter = 2.0E4;
  // time to reenter calculation
  // time only used when directinput==3

  if (runtype == 0) {
    deleteolddat = 1;		// 0: keep old data 1: delete them
    deleteoldimg = 1;		// 0: keep old images 1: delete them
    deleteoldpar = 1;		// 0: keep old par files 1: delete them
  } else {
    deleteolddat = 0;		// assume want to keep when reentrant
    deleteoldimg = 0;
    deleteoldpar = 0;
  }
  if (runtype == 0) {
    directinput = 0;		// NO choice
  }
  // assume that directinput>0 means normally want to continue
  // ener/loss 
  // data sets too, and that directinput==0 means starting fresh from
  // seperate data set as initial data, so no appending wanted
  if (runtype == 0) {
    appendold = 0;		// 0: no append time series data file
    // 1:
    // do append
  } else {
    if (directinput > 0) {
      appendold = 1;		// assume want to append if doing
      // reentrance of any kind
    } else {
      appendold = 0;		// assume no append wanted when
      // inputting
      // files using directinput==0
    }
  }



  // -1=test solution 1
  // 0=no analytic solution output
  // 5=tori1
  // 8=inject: initial setup of just atmosphere
  analoutput = 8;



  press = 1;			// 0->no pressure at all 1->pressure at
  // all (includes curvature)
  trans = 1;			// 0->no trans 1->do trans step at all

  transx1 = 1;			// 0->no trans 1->do trans step in
  // x1-dir
  transrhox1 = 1;		// 0->do not transport mass density in
  // x1-dir 1->do
  transiex1 = 1;		// 0->do not transport internal energy
  // density in x1-dir 1->do
  transv1x1 = 1;		// 0->do not transport velocity1 in
  // x1-dir 
  // 1->do
  transv2x1 = 1;		// 0->do not transport velocity2 in
  // x1-dir 
  // 1->do
  transv3x1 = 1;		// 0->do not transport velocity3 in
  // x1-dir 
  // 1->do


  transx2 = 1;			// 0->no trans 1->do trans step in
  // x2-dir
  transrhox2 = 1;		// 0->do not transport mass density in
  // x2-dir 1->do
  transiex2 = 1;		// 0->do not transport internal energy
  // density in x2-dir 1->do
  transv1x2 = 1;		// 0->do not transport velocity1 in
  // x2-dir 
  // 1->do
  transv2x2 = 1;		// 0->do not transport velocity2 in
  // x2-dir 
  // 1->do
  transv3x2 = 1;		// 0->do not transport velocity3 in
  // x2-dir 
  // 1->do

  ie = 1;			// 0-> no ie step 1->do ie step
  // (controls
  // art viscosity steps too)

  visc_art = 1;			// 0->no artificial visc step 1->do
  // visc
  // step

  visc_real = 1;		// 0-> do not do real visc 1-> do
  vreal = 4;
  // 1: gammie1(visce13 on)
  // 2: igumenshchev1(all visce terms on)
  // 3: stone et al. 1(visce13,visce23 on)
  // 4: mg (igu with sin)
  vischeat = 1;


  mdotin = 1;			// 0->no mass influx on grid, 1: yes
  // (see
  // injection())
  mag = 0;			// not used

  res = 0;			// not used


  advint = 1;
  // 0=DONOR
  // 1=VANLEER

  kever = 0;
  // 0: use native fluxes when possible for ke (seems to be best)
  // 1: use mass flux as basis for flux of ke

  // zone to start/stop on, in which to integrate over for integrated
  // values, and flux terms
  // used to avoid problems with boundary conditions
  // NOTE: Outer boundary is loop ending value where <INTOX, not
  // <=INTOX 
  // !
  // restrict boundary to avoid error due to outflow boundary
  // condition, 
  // factor of 20 for injection case in error rate decrease
  // global, need to correct per cpu
  intix1 = 2;
  intox1 = N1;			// should keep close to outer edge
  // since
  // otherwise too different for different
  // res, and injection won't add up
  intix2 = 0;
  intox2 = N2;
  intix3 = 0;
  intox3 = N3;

  nonunigridx1 = 5;		// 0=unigrid 1+=nonunigrid
  nonunigridx2 = 0;		// 0=unigrid 1+=nonunigrid
  nonunigridx3 = 0;		// 0=unigrid 1+=nonunigrid
  // x1: 1 : Cos in r(puts highest res in middle)
  // 2 : Log in r (puts highest res in inner r-edge)
  // 3 : Power law per decade
  // 4 : Power law per decade (dx/x=constant)
  // 5 : Power law per decade equal, spacing in log(r-rgp)
  // x2: 1 : Cos in theta(puts highest res in middle of theta=pi/2
  // 2 : Log split at middle
  // 3 : Power law per decade




  // stuff only needed if not reading in parameter data
  if ((runtype == 1) || (runtype == 0)) {

    simplebc = 1;		// whether to use simple boundary
    // conditions or not
    // those simple preferences
    bcix1 = 4;
    bcox1 = 4;
    bcix2 = 1;
    bcox2 = 1;
    bcix3 = 1;
    bcox3 = 1;

    // //////////// no need to modify/////
    // /
    // both inner and outer must be 5 above!
    periodicx1 = ((bcix1 == 5) || (bcox1 == 5) ? 1 : 0);	// 0:
								// not
    // periodic 
    // on
    // x1-dir
    // 1: is
    // if 1, implies global reflection on that boundary
    reflectix1 = ((bcix1 == 1) || (bcix1 == 2) ? 1 : 0);	// if
								// r=0
    // is
    // inner
    // edge,
    // so
    // reflecting 
    // in spc
    reflectox1 = ((bcox1 == 1) || (bcox1 == 2) ? 1 : 0);	// outer
    // radial
    // reflection
    // skips used to avoid computing inner edge values for that
    // direction.  Needed to ensure compute all relevant stuff and
    // only relevant stuff.  Generally 0 for periodic and 1
    // otherwise.
    skipix1 = ((periodicx1 == 1) ? 0 : 1);	// what zone to start
    // at
    // on inner edge
    // allows not to compute inner r-edge zones.  Set to 0 if not
    // reflecting at r=0 OR 1 if reflecting at r=0 or if don't want 
    // to 
    // calculate inner r-edge zone because it's a boundary
    // zone(e.g.
    // outflow)

    periodicx2 = ((bcix2 == 5) || (bcox2 == 5) ? 1 : 0);
    // if 1, implies global reflection on that boundary
    reflectix2 = ((bcix2 == 1) || (bcix2 == 2) ? 1 : 0);	// as
    // above
    // but
    // with x2 
    // dir
    reflectox2 = ((bcox2 == 1) || (bcox2 == 2) ? 1 : 0);	// as
    // above
    // but
    // with x2 
    // dir
    skipix2 = ((periodicx2 == 1) ? 0 : 1);	// what zone to start
    // at
    // on inner edge of
    // x2-grid
    // skip inner x2 edge when is boundary zone.


    periodicx3 = ((bcix3 == 5) || (bcox3 == 5) ? 1 : 0);
    reflectix3 = ((bcix3 == 1) || (bcix3 == 2) ? 1 : 0);	// for
    // completeness
    reflectox3 = ((bcox3 == 1) || (bcox3 == 2) ? 1 : 0);	// for
    // completeness
    skipix3 = ((periodicx3 == 1) ? 0 : 1);	// for completeness
    // /
    // //////////// no need to modify/////


    // remember to set REFLECTIX1=1 if you set below to 0 for spc.
    // and set SKIPIX1 to 1
    rg = 2.0;			// rg=2GM/c^2
    rgp = rg;			// PW-potential gravitational
    // radius(set
    // by units==2 if chosen) (rg==0 if
    // Newtonian pot)
    // rgp = 0.0; // newtonian pot
    x1in = (1.4 * rg);
    x1out = ((20.0 + 1.4) * rg);
    x2in = (0.0);
    x2out = (M_PI);
    x3in = (0.0);
    x3out = (2.0 * M_PI);

    gam = (5.0 / 3.0);		// adiabatic index
    alpha_real0 = 0.01;		// alpha(see step.c:tdep_compute())
    n_real = 2.0;		// power of sin(theta) for Length scale
    // for Gammie visc
    coolfact = 0;		// not used
    alpha = 5.0 / 2.0;		// equipartition parameter(for gam=1
    // case)
    cs = (1.0);			// sound speed, only used for gam=1.0

    tstart = 0.0;		// start unless reading in start time
    // tf = 10.0;
    // tf = 4000.0;
    // tf = .245;
    // tf = 0.240869213846978;
    // tf = .279;
    // tf = 50.0*(2.*M_PI);
    // tf = 3.0+4.4*(2.*M_PI) ; 
    // //tf = 3.0+2.3*(2.*M_PI) ; // end time
    // tf = 100.0*(2.*M_PI) ; // end time
    tvisc = pow(x1out, 1.5) / alpha_real0;	// viscous time scales
    numvtime = 4.0;
    tf = numvtime * tvisc;	// N viscous time scales
    // tf = 1E5;
    // tf=50.0;

    timagescale = tvisc;

    // #define DTDUMP (1.0)
    // #define DTCONS (0.1)
    // DTl = (1./1000.)*(2.*M_PI) ; /* logging period(lower limit
    // on
    // all except DTtimestep, this is also the ener/loss file
    // output
    // period */

    DTl = (5.0 / 1.0);

    // strictly for purposes of reentrant ability, should have
    // dump/floor same DT
    DTfloor = DTd = DTl * 300.0;	/* dumping period(data and
					   analytic) */
    // DTfloor=DTd = DTDUMP;
    // strictly for reentract ability, DTi should be integral
    // number
    // of DTds.
    DTpd = 10.0 * DTd;		// always reentrant data
    // DTpd = DTDUMP;
    DTi = DTl * 30.0;		/* image period */
    // DTi=DTDUMP;

    // below not restricted by DTl calls, just each self
    DTener = DTl * 1.0;		// negative means do every time step
    // (must 
    // be multiple of DTl)
    // f2's FFT shows could do 20.0 here(was 4.0)
    DTloss = DTl * 1.0;		// negative means do every time step
    // (must 
    // be multiple of DTl
    // DTloss=DTCONS;
    DTtimestep = DTl * 400.0;	// how often to output dominate
    // term in timestep to failfile
    DTsp = DTl * 300.0;		// how often to output sonic point info
    // DTtimestep=DTCONS;

    // //////
    // generally below don't change


    // x1
    L[1][1] = x1in;
    // L[1][1]= -1.0;
    L[2][1] = (x1out - x1in);	// width of x1-direction
    // L[2][1]= 2.0;

    // x2
    // remember to set reflectix2=1 if you set below to 0 for spc.
    L[1][2] = x2in;		// x2 inner edge position
    // L[1][2]= -1.0 ;
    L[2][2] = (x2out - x2in);	// width of x2-direction
    // L[2][2]= 2.0;

    // x3
    L[1][3] = x3in;		// x3 inner edge position
    L[2][3] = (x3out - x3in);	// width of x3-direction


    // nu_l==0.05(min) to .1 to .2(max) works best for oscillations
    // near shocks
    // the bigger you make this, the lower the oscillations, but
    // the
    // longer time it takes and the worse errors will become around
    // shocks due to non-oscilatting parts.  Can let some
    // oscillations 
    // appear for sake of those non-oscillating errors.
    // nu_vnr==3.0 works best for shocks
    nu_vnr = 3.0;		// Von Neuman Rictmeyer artificial
    // viscosity
    nu_l = 0.00;		// linear viscosity(not needed unless
    // want 
    // pretty shocks)(see global.h)
    nu_ten = 3.0;		// tensor viscosity
    GRAVC = 1.0;		// Gravitational Constant
    MASSBH = 1.0;		// Mass of central object
    GM = 1.0;			// as units, should be set in
    // analsol.c,
    // not here
    dt = 1.e10;			// start time step, used as reference,
    // keep large!
    cour = 0.5;			// courant condition
    cour2 = 0.25;		// diffuse courant condition

    resist = 0.0;		// not used
    nu_sh = 0.0;		// not used
    vg[1] = 0.0;		// velocity of grid in x1-dir, not used
    vg[2] = 0.0;		// velocity of grid in x2-dir, not used
    vg[3] = 0.0;		// velocity of grid in x3-dir, not used
  }
  /* END: global params */
  return (0);
}

// 1/2/3-d valid
int init_reentrance2(SFTYPE time)
{
  if (fabs(time - tstart) < 1.0E-8 * time + 1.0E-10) {
    ireenter = 0;
  } else if (directinput == 2)
    ireenter = 1;
  else if (directinput == 3)
    ireenter = 0;

  // must have time to do this, so call after time set
  pdump_start = (int) ((time - tstart) / DTpd) + ireenter;
  dump_start = (int) ((time - tstart) / DTd) + ireenter;
  npdump_start = (int) ((time - tstart) / DTd) + ireenter;
  adump_start = (int) ((time - tstart) / DTd) + ireenter;
  floor_start = (int) ((time - tstart) / DTfloor) + ireenter;
  image_start = (int) ((time - tstart) / DTi) + ireenter;

  if (ADUMPFLAG == -1) {
    adump_start = 0;		// force since only 1
  }
  if (PDUMPFLAG == 0)
    pdump_start = 0;
  if (DUMPFLAG == 0)
    dump_start = 0;
  if (NPDUMPFLAG == 0)
    npdump_start = 0;
  if (FLOORDUMPFLAG == 0)
    floor_start = 0;
  if (ADUMPFLAG == 0)
    adump_start = 0;
  if (IMAGEFLAG == 0)
    image_start = 0;

  if (myid <= 0) {
    fprintf(logfull_file, "starts: %d %d %d %d %d %d\n", pdump_start,
	    dump_start, npdump_start, adump_start, floor_start,
	    image_start);
    fflush(logfull_file);
  }
  return (0);
}

// 1/2/3-d valid
int init_reentrance(void)
{
  FILE *in;
  int numpdumps, numdumps, numadumps, numfloordumps, numimages,
      numnpdumps;
  char temps[MAXFILENAME];


  if (directinput == 3) {
    init_reentrance2(timereenter);
    if (myid <= 0) {
      fprintf(logfull_file, "Level 3 reentrance\n");
    }
  } else if (directinput == 2) {

    // determine number of primitive dumps
    sprintf(temps, "%s0_numpdumps%s", DATADIR, DATEXT);
    if ((in = fopen(temps, "r")) == NULL) {
      fprintf(fail_file, "error opening dump output file %s\n", temps);
      myexit(1);
    }
    while (fgetc(in) != '\n');
    fscanf(in, "%d", &numpdumps);
    fclose(in);
    if (myid <= 0) {
      fprintf(logfull_file, "number of pdumps: %d\n", numpdumps);
    }
    // determine number of normal dump files
    sprintf(temps, "%s0_numdumps%s", DATADIR, DATEXT);
    if ((in = fopen(temps, "r")) == NULL) {
      fprintf(fail_file, "error opening dump output file %s\n", temps);
      myexit(1);
    }
    while (fgetc(in) != '\n');
    fscanf(in, "%d", &numdumps);
    fclose(in);
    if (myid <= 0) {
      fprintf(logfull_file, "number of dumps: %d\n", numdumps);
    }

    sprintf(temps, "%s0_numnpdumps%s", DATADIR, DATEXT);
    if ((in = fopen(temps, "r")) == NULL) {
      fprintf(fail_file, "error opening dump output file %s\n", temps);
      myexit(1);
    }
    while (fgetc(in) != '\n');
    fscanf(in, "%d", &numnpdumps);
    fclose(in);

    if (myid <= 0) {
      fprintf(logfull_file, "number of np dumps: %d\n", numnpdumps);
    }

    sprintf(temps, "%s0_numadumps%s", DATADIR, DATEXT);
    if ((in = fopen(temps, "r")) == NULL) {
      fprintf(fail_file, "error opening dump output file %s\n", temps);
      myexit(1);
    }
    while (fgetc(in) != '\n');
    fscanf(in, "%d", &numadumps);
    fclose(in);

    if (myid <= 0) {
      fprintf(logfull_file, "number of adumps: %d\n", numadumps);
    }

    sprintf(temps, "%s0_numfloordumps%s", DATADIR, DATEXT);
    if ((in = fopen(temps, "r")) == NULL) {
      fprintf(fail_file, "error opening dump output file %s\n", temps);
      myexit(1);
    }
    while (fgetc(in) != '\n');
    fscanf(in, "%d", &numfloordumps);
    fclose(in);

    if (myid <= 0) {
      fprintf(logfull_file, "number of floor dumps: %d\n",
	      numfloordumps);
    }
    // get # of images
    sprintf(temps, "%s%s0_numimages%s", DATADIR, "i/", DATEXT);
    if ((in = fopen(temps, "r")) == NULL) {
      fprintf(fail_file, "error opening dump output file %s\n", temps);
      myexit(1);
    }
    while (fgetc(in) != '\n');	// skip comment line
    fscanf(in, "%d", &numimages);
    fclose(in);

    if (myid <= 0) {
      fprintf(logfull_file, "number of image dumps: %d\n", numimages);
    }
    // first guess
    pdump_start = numpdumps - 1;
    dump_start = numdumps - 1;
    npdump_start = numnpdumps - 1;
    adump_start = numadumps - 1;
    floor_start = numfloordumps - 1;
    image_start = numimages - 1;

    // use count as reference(this is the file one reads)
    if (DUMPSM == 1) {
      pdump_start = pdump_start;
    } else {
      dump_start = dump_start;
    }
  } else {
    // only exactly valid for good DT if match t on dump and
    // images,
    // and only can do that if dumps occur at image dump interval.
    // So dump DT needs to be integral of image DT.

    pdump_start = 0;
    dump_start = npdump_start = 0;
    adump_start = 0;
    floor_start = 0;
    image_start = 0;
  }


  /* END: reentrance params */
  fflush(log_file);
  if (myid <= 0) {
    fflush(logfull_file);
  }
  return (0);
}



int init_otherparams(void)	// stuff not stored to gparm file
{
  if ((runtype == 3) || (runtype == 1)) {
    // t set by input of data
  } else {
    t = tstart;			// start time
  }
  // compute things that always need but formula doesn't change
  invcour = 1.0 / cour;		// inverse of courant number
  invcour2 = 1.0 / cour2;	// inverse of diffuse courant number
  DTavg = 0;
  wgam53 = (int) (fabs(gam - 5. / 3.) < ERR);	// switch to see if
  // gam==5/3
  wgam1 = (int) (fabs(gam - 1.) < ERR);	// switch to see if gam==1
  wgam = (int) (fabs(gam - 1.) > ERR);	// switch to see if need
  // to compute pressure
  wpw = (int) (fabs(rgp - 0.) > ERR);	// switch to see if newtonian
  // or
  // PW pot

  // here since don't want to dump and still want reentrance
  IOBound[0] = 0.1;		// first transition from out to inflow
  IOBound[1] = 0.9;		// second transition from inflow to
  // outflow
  IOBound[2] = M_PI - IOBound[1];	// third transition from inflow 
					// to 
  // outflow
  IOBound[3] = M_PI - IOBound[0];	// fourth transition from
					// inflow
  // to outflow

  return (0);
}


int init_mainbc(int px1, int six1, int rix1, int rox1, int px2,
		int six2, int rix2, int rox2, int px3, int six3,
		int rix3, int rox3)
{
  // intix2 stuff not really right if account region goes inside
  // inner
  // cpus, but this should never be done anyways!


  periodicx1 = px1;
  skipix1 = six1;
  reflectix1 = rix1;
  reflectox1 = rox1;
  if (skipix1 == 1) {
    if (intix1 <= 1) {
      skipintix1 = 1;
    } else
      skipintix1 = intix1;
  } else
    skipintix1 = intix1;

  // x2
  if (numprocs == 1) {
    periodicx2 = px2;
    skipix2 = six2;
    reflectix2 = rix2;
    reflectox2 = rox2;
    if (skipix2 == 1) {
      if (intix2 <= 1) {
	skipintix2 = 1;
      } else
	skipintix2 = intix2;
    } else
      skipintix2 = intix2;
  }
  if (numprocs > 1) {
  }
  // x3
  periodicx3 = px3;
  skipix3 = six3;
  reflectix3 = rix3;
  reflectox3 = rox3;
  if (skipix3 == 1) {
    if (intix3 <= 1) {
      skipintix3 = 1;
    } else
      skipintix3 = intix3;
  } else
    skipintix3 = intix3;

  return (0);
}



// 1/2/3-d valid
int init_bc(int simple, int ix1, int ox1, int ix2, int ox2, int ix3,
	    int ox3)
{
  int i, j, k, l, m, n;
  int N[3 + 1];
  int numbc[3 + 1];
  int skipi[3 + 1];
  int itemp;

  N[1] = N1;
  N[2] = N2;
  N[3] = N3;
  numbc[1] = N1BND;
  numbc[2] = N2BND;
  numbc[3] = N3BND;
  skipi[1] = skipix1;
  skipi[2] = skipix2;
  skipi[3] = skipix3;

  for (l = 1; l <= NUMSCA - 1; l++) {
    LOOP {
      /* Stick here the assignment of boundary conditions for active
         zones */
      bcs[l][1][k][j][i] = 0;	// other entries do not matter
    }
  }

  for (l = 1; l <= NUMVEC - 1; l++) {
    LOOP {
      /* Stick here the assignment of boundary conditions for active
         zones */
      bcv[l][1][k][j][i] = 0;	// other entries do not matter
    }
  }

  /* Now assign bc for boundary zones--good for normal boundary zones
     in any situation - just change type from here with conditions */

  // go over all except potential
  for (l = 1; l <= NUMSCA - NOBOUNDPOT; l++) {
    for (m = 1; m <= 3; m++) {
      for (j = -numbc[3 - (4 - m) % 3];
	   j < N[3 - (4 - m) % 3] + numbc[3 - (4 - m) % 3]; j++) {
	for (i = -numbc[m % 3 + 1];
	     i < N[m % 3 + 1] + numbc[m % 3 + 1]; i++) {
	  /* m%3+1 gives next 1->2,2->3,3->1 3-(4-m)%3 gives previous
	     1->3,2->1,3->2 */
	  for (n = 0; n < numbc[m]; n++) {
	    if (m == 1) {	/* Assign over x=const boundaries */

	      if (simple == 1)
		itemp = ix1;
	      else
		itemp = 4;
	      // inner boundary
	      bcs[l][1][j][i][-numbc[m] + n] = itemp;
	      bcs[l][2][j][i][-numbc[m] + n] = m;
	      bcs[l][3][j][i][-numbc[m] + n] = 1;

	      // outer boundary
	      if (IOBOUNDARY == 0) {
		if (simple == 1)
		  itemp = ox1;
		else
		  itemp = 4;
	      } else {
		if (x[2][2][j] < IOBound[0]) {
		  itemp = 3;
		} else if ((x[2][2][j] >= IOBound[0])
			   && (x[2][2][j] < IOBound[1])) {
		  itemp = 4;
		} else if ((x[2][2][j] >= IOBound[1])
			   && (x[2][2][j] < IOBound[2])) {
		  itemp = 3;
		} else if ((x[2][2][j] >= IOBound[2])
			   && (x[2][2][j] < IOBound[3])) {
		  itemp = 4;
		} else if (x[2][2][j] >= IOBound[3]) {
		  itemp = 3;
		}
	      }

	      bcs[l][1][j][i][N[m] + n] = itemp;
	      bcs[l][2][j][i][N[m] + n] = m;
	      bcs[l][3][j][i][N[m] + n] = -1;

	    }
	    if (m == 2) {	/* Assign over y=const boundaries */

	      // setup for global grid, only need to
	      // change
	      // this

	      // inner boundary
	      if (simple == 1)
		itemp = ix2;
	      else
		itemp = 1;
	      bcs[l][1][i][-numbc[m] + n][j] = itemp;
	      bcs[l][2][i][-numbc[m] + n][j] = m;
	      bcs[l][3][i][-numbc[m] + n][j] = 1;

	      // outer boundary
	      if (simple == 1)
		itemp = ox2;
	      else
		itemp = 1;
	      bcs[l][1][i][N[m] + n][j] = itemp;
	      bcs[l][2][i][N[m] + n][j] = m;
	      bcs[l][3][i][N[m] + n][j] = -1;

	      // no need to change remaining local
	      // assignments

	      if (numprocs > 1) {
	      }
	    }
	    if (m == 3) {	/* Assign over z=const boundaries */

	      // inner x3
	      if (simple == 1)
		itemp = ix3;
	      else
		itemp = 1;
	      bcs[l][1][-numbc[m] + n][j][i] = itemp;
	      bcs[l][2][-numbc[m] + n][j][i] = m;
	      bcs[l][3][-numbc[m] + n][j][i] = 1;
	      // outer x3
	      if (simple == 1)
		itemp = ox3;
	      else
		itemp = 1;
	      bcs[l][1][N[m] + n][j][i] = itemp;
	      bcs[l][2][N[m] + n][j][i] = m;
	      bcs[l][3][N[m] + n][j][i] = -1;
	    }
	  }
	}
      }
    }
    // make corner zones special
    bcs[l][1][0][-1][-1] = -1;
    bcs[l][1][0][-1][N1] = -1;
    bcs[l][1][0][N2][-1] = -1;
    bcs[l][1][0][N2][N1] = -1;

  }

  if (NOBOUNDPOT == 1) {
    // potential is static everywhere for all time
    l = 3;
    LOOPF {
      /* Stick here the assignment of boundary conditions for active
         zones */
      bcs[l][1][k][j][i] = 0;	// other entries do not matter
    }
  }
  // now do vectors
  for (l = 1; l <= NUMVEC; l++) {
    for (m = 1; m <= 3; m++) {
      for (j = -numbc[3 - (4 - m) % 3];
	   j < N[3 - (4 - m) % 3] + numbc[3 - (4 - m) % 3]; j++) {
	for (i = -numbc[m % 3 + 1];
	     i < N[m % 3 + 1] + numbc[m % 3 + 1]; i++) {
	  /* m%3+1 gives next 1->2,2->3,3->1 3-(4-m)%3 gives previous
	     1->3,2->1,3->2 */
	  for (n = 0; n < numbc[m]; n++) {
	    if (m == 1) {	/* Assign over x=const boundaries */

	      // inner boundary
	      if (simple == 1)
		itemp = ix1;
	      else
		itemp = 4;
	      bcv[l][1][j][i][-numbc[m] + n] = itemp;
	      bcv[l][2][j][i][-numbc[m] + n] = m;
	      bcv[l][3][j][i][-numbc[m] + n] = 1;

	      // outer boundary
	      if (IOBOUNDARY == 0) {
		if (simple == 1)
		  itemp = ox1;
		else
		  itemp = 4;
	      } else {
		if (x[2][2][j] < IOBound[0]) {
		  itemp = 3;
		} else if ((x[2][2][j] >= IOBound[0])
			   && (x[2][2][j] < IOBound[1])) {
		  itemp = 4;
		} else if ((x[2][2][j] >= IOBound[1])
			   && (x[2][2][j] < IOBound[2])) {
		  itemp = 3;
		} else if ((x[2][2][j] >= IOBound[2])
			   && (x[2][2][j] < IOBound[3])) {
		  itemp = 4;
		} else if (x[2][2][j] >= IOBound[3]) {
		  itemp = 3;
		}
	      }
	      bcv[l][1][j][i][N[m] + n] = itemp;
	      bcv[l][2][j][i][N[m] + n] = m;
	      bcv[l][3][j][i][N[m] + n] = -1;
	    }
	    if (m == 2) {	/* Assign over y=const boundaries */

	      // setup for global grid, only need to
	      // change
	      // this

	      // inner boundary
	      if (simple == 1)
		itemp = ix2;
	      else
		itemp = 1;
	      bcv[l][1][i][-numbc[m] + n][j] = itemp;
	      bcv[l][2][i][-numbc[m] + n][j] = m;
	      bcv[l][3][i][-numbc[m] + n][j] = 1;
	      // outer boundary
	      if (simple == 1)
		itemp = ox2;
	      else
		itemp = 1;
	      bcv[l][1][i][N[m] + n][j] = itemp;
	      bcv[l][2][i][N[m] + n][j] = m;
	      bcv[l][3][i][N[m] + n][j] = -1;


	      // no need to change remaining local
	      // assignments

	      if (numprocs > 1) {
	      }

	    }
	    if (m == 3) {	/* Assign over z=const boundaries */

	      // inner x3
	      if (simple == 1)
		itemp = ix3;
	      else
		itemp = 1;
	      bcv[l][1][-numbc[m] + n][j][i] = itemp;
	      bcv[l][2][-numbc[m] + n][j][i] = m;
	      bcv[l][3][-numbc[m] + n][j][i] = 1;

	      // outer x3
	      if (simple == 1)
		itemp = ox3;
	      else
		itemp = 1;
	      bcv[l][1][N[m] + n][j][i] = itemp;
	      bcv[l][2][N[m] + n][j][i] = m;
	      bcv[l][3][N[m] + n][j][i] = -1;
	    }
	  }
	}
      }
    }
    // make corner zones special
    bcv[l][1][0][-1][-1] = -1;
    bcv[l][1][0][-1][N1] = -1;
    bcv[l][1][0][N2][-1] = -1;
    bcv[l][1][0][N2][N1] = -1;
  }



  return (0);
}



int init(int argc, char *argv[], char *envp[])
{
  int k, j, i;
  int seed;
  FTYPE beta, nj;
  int error = 0;
  char temps[100];
  int itemp;
  char extension[MAXFILENAME];


  // left over parameters
  seed = 13;
  ranc(seed);			/* fire up random number generator */
  beta = 0.1;
  nj = 1.0;


  error += init_mem();
  error += init_pointers();
  // initialize general stuff
  init_general();
  error += init_params(seed, beta, nj);	// must come before most
  // everything else

  // Setup directory structure in case do not exist
  if (myid <= 0) {
    sprintf(temps, "mkdir %s%s", DATADIR, "i/");
    system(temps);
    if (deleteolddat) {
      sprintf(temps, "rm %s*", DATADIR);
      strcat(temps, DATEXT);
      strcat(temps, "*");
      system(temps);
    }
    if (deleteoldpar) {
      sprintf(temps, "rm %s*", DATADIR);
      strcat(temps, PAREXT);
      strcat(temps, "*");
      system(temps);
    }
  }
  // get or set parameters needed for data(get/set)
  if ((runtype == 3) || (runtype == 2)) {
    error += init_runpar(directinput);	// function for par files
    error += init_otherparams();	// other parameter init
    if (directinput > 0) {
      error += init_reentrance2(t);
    }
  }
  if ((runtype == 1) || (runtype == 0)) {
    error += init_otherparams();
    error +=
	init_mainbc(periodicx1, skipix1, reflectix1, reflectox1,
		    periodicx2, skipix2, reflectix2, reflectox2,
		    periodicx3, skipix3, reflectix3, reflectox3);
    if (directinput > 0) {
      error += init_reentrance2(t);
    }
    error += init_dx(N1, N2, N3, N1BND, N2BND, N3BND, startx, 0, 0);
    error += init_x(N1, N2, N3, N1BND, N2BND, N3BND, startx, 0, 0);
    error += init_diffs();
    init_compsave();
    error +=
	init_bc(simplebc, bcix1, bcox1, bcix2, bcox2, bcix3, bcox3);
  }
  // get data
  if ((runtype == 3) || (runtype == 1)) {
    if (directinput > 0) {
      error += init_reentrance();
    }
    if (DUMPSM == 0) {		// only valid to read in dump files if
      // not 
      // in SM interpolated format(in general)
      if (directinput > 0)
	itemp = dump_start;
      else
	itemp = -1;
      error += init_rundat(itemp, 0);	// function for dump dat
      // files
    } else {
      if (directinput > 0)
	itemp = pdump_start;
      else
	itemp = -1;
      error += init_rundat(itemp, 100);	// function for dump dat
      // files
    }
    // no need to input adump right now
    // if(directinput>0) itemp=dump_start;
    // error+=init_rundat(itemp,1);//function for adump dat files
    // no input of floor needed anymore
    // if(directinput>0) itemp=floor_start; else itemp=-1;
    // error+=init_rundat(itemp,2);//function for floor dat files
  }
  // set data
  if ((runtype == 0) || (runtype == 2)) {
    error += init_data();
    // only init below 2 if not reentrant on data
    init_inflows();
    init_radiations();
  } else {
    // need to solve analytic solution for many cases: need bc
    // data,
    // etc...just as with nonloading case.
    analsolve(0);
    tdep_compute();		// compute time dep stuff
    // place any overrides here.  e.g. file loads potential, but
    // potential is not bounded since static and needs proper
    // static
    // dependence on boundary zones, so force s[3] to be sanal[3],
    // ignoring file data, for wherever you assign the override.
    // Can
    // use init_data as template.

    LOOPF {
      s[3][k][j][i] = sanal[3][k][j][i];
    }
  }

  if (DOPARDIAG) {
    error += init_outgparm(-1);	// -1 means out all within define
    // params
  }


  /* enforce boundary conditions on all scalars and vectors */
  bound(NULL, NULL, -1, -1);

  // init final things
  accountstoreset();		// must be before any accounting
  init_floor();
  init_visc();
  init_loss();


  if (appendold == 1) {
    sprintf(WRITETYPE, "at+");
  } else {
    sprintf(WRITETYPE, "wt");
  }
  strcpy(extension, OUTEXT);


  // setup computational log files
  if (DODTDIAG) {
    sprintf(temps, "%s0_logdt%s%s", DATADIR, extension, myidtxt);

    if ((logdt_file = fopen(temps, WRITETYPE)) == NULL) {	// naive
      // append
      // if
      // appendold==1
      fprintf(stderr, "dtdiag: Cannot open: %s\n", temps);
      exit(1);
    }
  }
  if (DOLOGSTEP) {
    if (myid <= 0) {
      sprintf(temps, "%s0_logstep%s", DATADIR, extension);

      if ((logstep_file = fopen(temps, WRITETYPE)) == NULL) {	// naive 
								// 
	// 
	// append 
	// if 
	// appendold==1
	fprintf(stderr, "logstep: Cannot open: %s\n", temps);
	exit(1);
      }
    }
  }
  if (DOLOGPERF) {
    if (myid <= 0) {
      sprintf(temps, "%s0_logperf%s", DATADIR, extension);

      if ((logperf_file = fopen(temps, WRITETYPE)) == NULL) {	// naive 
								// 
	// 
	// append 
	// if 
	// appendold==1
	fprintf(stderr, "logperf: Cannot open: %s\n", temps);
	exit(1);
      }
    }
  }

  if (DOFLOORDIAG > 0) {

    sprintf(temps, "%s0_logfl%s%s", DATADIR, extension, myidtxt);

    if ((logfl_file = fopen(temps, WRITETYPE)) == NULL) {	// naive
      // append
      // if
      // appendold==1
      fprintf(stderr, "floordiag: Cannot open: %s\n", temps);
      exit(1);
    }
  }


  return (error);
}


// expensive for 1-D problems to do each time step
void init_loss(void)
{
  int i, j, k, l, m;

  // initialize loss data
  for (i = 1; i <= NUMSCA; i++) {
    for (j = 1; j <= 3; j++) {
      for (k = 0; k <= 1; k++) {
	for (l = 0; l < NBIG; l++) {
	  losss[i][j][k][l] = 0.0;
	}
      }
    }
  }
  for (i = 1; i <= NUMVEC; i++) {
    for (m = 0; m <= 3; m++) {	// 0 for KE
      for (j = 1; j <= 3; j++) {
	for (k = 0; k <= 1; k++) {
	  for (l = 0; l < NBIG; l++) {
	    lossv[i][m][j][k][l] = 0.0;
	  }
	}
      }
    }
  }
  for (i = 1; i <= 1; i++) {
    for (j = 1; j <= 3; j++) {
      for (k = 0; k <= 1; k++) {
	for (l = 0; l < NBIG; l++) {
	  lossvisc[i][j][k][l] = 0.0;
	}
      }
    }
  }

}

void init_general(void)		// as is setup for a tori case.  Needed 
				// 
																// to 
																// 
				// avoid general floating point issues
				// in
				// image routines if not assigned
{
  int i, j;


  // scalars

  j = 0;			// normal output
  i = 0;			// view large

  mms[i][j][1][0] = (9. / 10.) * 1.e-06;
  mms[i][j][1][1] = 1.02012908;

  mms[i][j][2][0] = (9. / 10.) * 1.33333333e-11;
  mms[i][j][2][1] = 0.7428376443;

  mms[i][j][3][0] = -3.265618839;
  mms[i][j][3][1] = -0.1881626381;

  j = 1;			// second type of comp
  i = 0;			// view large

  mms[i][j][1][0] = (9. / 10.) * 1.e-06;
  mms[i][j][1][1] = 1.02012908;

  mms[i][j][2][0] = (9. / 10.) * 1.33333333e-11;
  mms[i][j][2][1] = 0.7428376443;

  mms[i][j][3][0] = -3.265618839;
  mms[i][j][3][1] = -0.1881626381;

  j = 0;			// normal comp
  i = 1;			// view zoom
  // for alpha=.01 with both terms
  mms[i][j][1][0] = .0001;
  mms[i][j][1][1] = .2;

  mms[i][j][2][0] = .0001;
  mms[i][j][2][1] = 0.7428376443;

  mms[i][j][3][0] = -3.265618839;
  mms[i][j][3][1] = -0.1881626381;

  j = 1;			// 2nd comp
  i = 1;			// view zoom
  // for alpha=.01 with both terms
  mms[i][j][1][0] = .0001;
  mms[i][j][1][1] = .1;

  mms[i][j][2][0] = .0001;
  mms[i][j][2][1] = 0.7428376443;

  mms[i][j][3][0] = -3.265618839;
  mms[i][j][3][1] = -0.1881626381;


  // vectors

  i = 0;			// normal view
  j = 0;			// 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0] = 0.0;
  mmv[i][j][1][0][1] = 10.0;
  // vx1
  mmv[i][j][1][1][0] = -2.0;
  mmv[i][j][1][1][1] = .1;
  // vx2
  mmv[i][j][1][2][0] = -1.5;
  mmv[i][j][1][2][1] = 1.5;
  // vx3
  mmv[i][j][1][3][0] = -.1;
  mmv[i][j][1][3][1] = 10.0;

  // B0
  mmv[i][j][2][0][0] = 0;
  mmv[i][j][2][0][1] = 0;
  // Bx1
  mmv[i][j][2][1][0] = 0;
  mmv[i][j][2][1][1] = 0;
  // Bx2
  mmv[i][j][2][2][0] = 0;
  mmv[i][j][2][2][1] = 0;
  // Bx3
  mmv[i][j][2][3][0] = 0;
  mmv[i][j][2][3][1] = 0;

  i = 0;
  j = 1;
  // rho times v or B
  // must get from sm plots

  // rho*v0
  mmv[i][j][1][0][0] = 0.0;
  mmv[i][j][1][0][1] = 1.0;
  // rho*vx1
  mmv[i][j][1][1][0] = -0.01;
  mmv[i][j][1][1][1] = 0.01;
  // rho*vx2
  mmv[i][j][1][2][0] = -0.01;
  mmv[i][j][1][2][1] = .01;
  // rho*vx3
  mmv[i][j][1][3][0] = -.1;
  mmv[i][j][1][3][1] = 10.0;

  // rho*B0
  mmv[i][j][2][0][0] = 0;
  mmv[i][j][2][0][1] = 0;
  // rho*Bx1
  mmv[i][j][2][1][0] = 0;
  mmv[i][j][2][1][1] = 0;
  // rho*Bx2
  mmv[i][j][2][2][0] = 0;
  mmv[i][j][2][2][1] = 0;
  // rho*Bx3
  mmv[i][j][2][3][0] = 0;
  mmv[i][j][2][3][1] = 0;




  i = 1;			// zoom view (vectors)
  j = 0;			// 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0] = 0.0;
  mmv[i][j][1][0][1] = 10.0;
  // vx1
  mmv[i][j][1][1][0] = -2.0;
  mmv[i][j][1][1][1] = .1;
  // vx2
  mmv[i][j][1][2][0] = -1.5;
  mmv[i][j][1][2][1] = 1.5;
  // vx3
  mmv[i][j][1][3][0] = -.1;
  mmv[i][j][1][3][1] = 10.0;

  // B0
  mmv[i][j][2][0][0] = 0;
  mmv[i][j][2][0][1] = 0;
  // Bx1
  mmv[i][j][2][1][0] = 0;
  mmv[i][j][2][1][1] = 0;
  // Bx2
  mmv[i][j][2][2][0] = 0;
  mmv[i][j][2][2][1] = 0;
  // Bx3
  mmv[i][j][2][3][0] = 0;
  mmv[i][j][2][3][1] = 0;

  i = 1;			// zoom view (vectors)
  j = 1;
  // rho times v or B
  // must get from sm plots

  // rho*v0
  mmv[i][j][1][0][0] = 0.0;
  mmv[i][j][1][0][1] = 1.0;
  // rho*vx1
  mmv[i][j][1][1][0] = -0.01;
  mmv[i][j][1][1][1] = 0.01;
  // rho*vx2
  mmv[i][j][1][2][0] = -0.01;
  mmv[i][j][1][2][1] = .01;
  // rho*vx3
  mmv[i][j][1][3][0] = -.1;
  mmv[i][j][1][3][1] = 1.0;

  // rho*B0
  mmv[i][j][2][0][0] = 0;
  mmv[i][j][2][0][1] = 0;
  // rho*Bx1
  mmv[i][j][2][1][0] = 0;
  mmv[i][j][2][1][1] = 0;
  // rho*Bx2
  mmv[i][j][2][2][0] = 0;
  mmv[i][j][2][2][1] = 0;
  // rho*Bx3
  mmv[i][j][2][3][0] = 0;
  mmv[i][j][2][3][1] = 0;




  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for (i = 0; i < ITYPES; i++) {	// both views
    for (j = 0; j < CTYPES; j++) {	// both comps
      outerdefs[i][j][1] = mms[i][j][1][0];	// rho
      outerdefs[i][j][2] = mms[i][j][2][0];	// en
      outerdefs[i][j][3] = mms[i][j][3][0];	// pot

      outerdefv[i][j][1][0] = mmv[i][j][1][0][0];	// magnitude of 
							// 
      // v
      outerdefv[i][j][1][1] = mmv[i][j][1][1][0];	// v1
      outerdefv[i][j][1][2] = mmv[i][j][1][2][0];	// v2
      outerdefv[i][j][1][3] = mmv[i][j][1][3][0];	// v3

      outerdefv[i][j][2][0] = mmv[i][j][2][0][0];
      outerdefv[i][j][2][1] = mmv[i][j][2][1][0];
      outerdefv[i][j][2][2] = mmv[i][j][2][2][0];
      outerdefv[i][j][2][3] = mmv[i][j][2][3][0];
    }
  }
  // for dumps(never use 2nd comp type, and always same for any view)
  douterdefs[1] = mms[0][0][1][0];	// rho
  douterdefs[2] = mms[0][0][2][0];	// en
  douterdefs[3] = mms[0][0][3][0];	// pot

  douterdefv[1][0] = 0.0;	// magnitude of v
  douterdefv[1][1] = 0.0;	// v1
  douterdefv[1][2] = 0.0;	// v2
  douterdefv[1][3] = 0.0;	// v3

  douterdefv[2][0] = 0.0;
  douterdefv[2][1] = 0.0;
  douterdefv[2][2] = 0.0;
  douterdefv[2][3] = 0.0;

}



void init_visc(void)
{
  int i, j, k;

  // init visc param(part the never changes)(must be consistent with
  // step.c)
  LOOPF {
    if (vreal == 1) {
      // alpha-parameterization
      // $\nu = \alpha r sin(\theta)^(n_real)
      nu_fact[k][j][i] = x[2][1][i] * pow(g[2][4][j], n_real);
    }

    if (vreal == 2) {
      // Igumenshchev alpha param: $\nu=\alpha*cs^2/Omega_{k}$
      // here by Keplerian we mean spherical rotation, even if
      // fluid 
      // is meant to rotate on cylinders.
      // nu_fact[k][j][i]=pow(x[2][1][i],1.5)/sqrt(GM);
      // below was being used for all IGU runs, when IGU really
      // uses 
      // above!
      // nu_fact[k][j][i]=pow(x[2][1][i]*g[2][4][j],1.5)/sqrt(GM);
      // below for PW pot, use IGU like below
      nu_fact[k][j][i] =
	  pow(x[2][1][i], 1.5) * (1.0 - rgp / x[2][1][i]) / sqrt(GM);
    }

    if (vreal == 3) {
      // Stone et al.
      // Run A through J with different alpha (J diff nu_real)
      nu_fact[k][j][i] = R0 * R0 * Omega0 / rho0;

      // Run K with different alpha
      // nu_fact[k][j][i]=alpha_real*sqrt(x[2][1][i]*g[2][4][j])
      // ;
    }
    if (vreal == 4) {
      // below was being used for all IGU runs, when IGU really
      // uses 
      // above!
      nu_fact[k][j][i] = pow(x[2][1][i] * g[2][4][j], 1.5) / sqrt(GM);
      // below for PW pot, use IGU like below
      // nu_fact[k][j][i]=pow(x[2][1][i],1.5)*(1.0-rgp/x[2][1][i])/sqrt(GM);
    }
#if(analoutput==6)
    // for test of the visc real code
    nu_fact[k][j][i] = 1.0;
#endif

  }
}


void init_floor(void)
{
  int i, j, k, l, m;

  // init floor variables
  LOOP {
    for (l = 1; l <= NUMSCA; l++) {
      floorvars[l][k][j][i] = 0.0;
    }
    for (l = 1; l <= NUMVEC; l++) {
      floorvar0[l][k][j][i] = 0.0;
      for (m = 1; m <= 3; m++) {
	floorvarv[l][m][k][j][i] = 0.0;
      }
    }
  }
  // initialize floor checks
  for (i = 0; i <= NUMFLOORVAR; i++) {
    floors[i] = 0.0;
    floorlowest[i] = 1.0;
    wherelowest[i] = -1;
  }
}

void init_inflows(void)
{
  int i;

  // init inflows variables
  for (i = 0; i <= NUMLOSSVAR; i++) {
    inflows[i] = 0.0;
  }
}
void init_radiations(void)
{
  int i;

  // init inflows variables
  for (i = 0; i <= NUMLOSSVAR; i++) {
    radiations[i] = 0.0;
  }
}


// 1/2/3-d valid
int init_mem(void)
{

  /* END: Allocate memory for parameter structure */
  return (0);
}

void init_genfiles(int gopp)
{
  char temps[MAXFILENAME];
  char extension[MAXFILENAME];


  if (gopp == 1) {
    strcpy(extension, PPEXT);
  } else if (gopp == 0) {
    strcpy(extension, OUTEXT);
  }
  // always have fail and general log open

  sprintf(temps, "%s0_fail%s%s", DATADIR, extension, myidtxt);


  if ((fail_file = fopen(temps, "wt")) == NULL) {
    fprintf(stderr, "fail: Cannot open: %s\n", temps);
    exit(1);
  }

  sprintf(temps, "%s0_log%s%s", DATADIR, extension, myidtxt);

  if ((log_file = fopen(temps, "wt")) == NULL) {
    fprintf(stderr, "log: Cannot open: %s\n", temps);
    exit(1);
  }
  if (myid <= 0) {
    sprintf(temps, "%s0_logfull%s", DATADIR, extension);

    if ((logfull_file = fopen(temps, "wt")) == NULL) {
      fprintf(stderr, "logfull: Cannot open: %s\n", temps);
      exit(1);
    }
  }


}

int init_runpar(int gopp)
{
  FILE *in[3];
  char fname[200];
  int dumi[50];
  char ch;
  int i, j, k, l, m;
  char extension[MAXFILENAME];

  // assume: grid1, grid2, and gparam for parameter data

  // remember to make same changes to this and output in
  // diag.c/init.c!!

  if (gopp >= 1)
    gopp = 1;

  if (gopp == 1) {
    strcpy(extension, PAREXT);
  } else if (gopp == 0) {
    strcpy(extension, INEXT);
  }


  for (i = 0; i < numprocs; i++) {
    if (i == myid) {
      sprintf(fname, "%s0_gparam%s", DATADIR, extension);

      if ((in[0] = fopen(fname, "rt")) == NULL) {
	fprintf(fail_file, "gparam: Cannot open: %s\n", fname);
	myexit(1);
      }


      l = 0;
      ch = fgetc(in[l]);
      if (ch == '#') {
	while (fgetc(in[l]) != '\n');
      } else {
	fprintf(fail_file,
		"1: gparm: global parm file doesn't match expected format\n");
	myexit(1);
      }

      fscanf(in[l], HEADER3_S0, &dumi[0], &dumi[1]);
      if ((dumi[0] != PVER) || (dumi[1] != PTYPE)) {
	fprintf(fail_file,
		"Expected pver/ptype: %d %d got %d %d\n", PVER,
		PTYPE, dumi[0], dumi[1]);
	myexit(6);
      }
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S1, &dumi[0], &dumi[1], &dumi[2]);
      if ((dumi[0] != N1) || (dumi[1] != totalsize[2])
	  || (dumi[2] != N3)) {
	fprintf(fail_file,
		"2: gparm: Gridsize in gparm file needs to match global.h file\n");
	fprintf(fail_file, "got: N1: %d N2: %d N3: %d\n", dumi[0],
		dumi[1], dumi[2]);
	fprintf(fail_file, "expected: N1: %d N2: %d N3: %d\n", N1,
		totalsize[2], N3);
	myexit(1);
      }
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S2, &L[1][1], &L[1][2], &L[1][3],
	     &L[2][1], &L[2][2], &L[2][3]);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S3, &rg, &rgp, &cs, &coolfact, &gam,
	     &alpha_real0, &n_real);
      // fprintf(stdout,"%15.10g %15.10g %15.10g %15.10g %15.10g
      // %15.10g
      // %15.10g\n",rg,rgp,cs,coolfact,gam,alpha_real0,n_real);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S4, &nu_vnr, &nu_l, &nu_ten, &cour, &cour2);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S5, &GRAVC, &MASSBH);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S6, &tstart, &tf, &tf, &tf,
	     &tf, &timagescale);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S7, &DTl, &DTd, &DTi, &DTloss, &DTfloor,
	     &DTtimestep, &DTpd, &DTener);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S8, &resist, &nu_sh);
      while ((ch = fgetc(in[l])) != '\n');
      fscanf(in[l], HEADER3_S9, &vg[1], &vg[2], &vg[3]);

      fclose(in[0]);
    }
  }

  sprintf(fname, "%s0_grid1%s%s", DATADIR, extension, myidtxt);


  if ((in[1] = fopen(fname, "rt")) == NULL) {
    fprintf(fail_file, "grid1: Cannot open: %s\n", fname);
    myexit(1);
  }
  sprintf(fname, "%s0_grid2%s%s", DATADIR, extension, myidtxt);

  if ((in[2] = fopen(fname, "rt")) == NULL) {
    fprintf(fail_file, "grid2: Cannot open: %s\n", fname);
    myexit(1);
  }


  for (l = 1; l <= NUMGRID; l++) {
    ch = fgetc(in[l]);
    if (ch == '#') {
      while (fgetc(in[l]) != '\n');
    } else {
      fprintf(fail_file,
	      "1: grid data file doesn't match expected format: %d\n",
	      l);
      myexit(1);
    }
    fscanf(in[l], "%d %d\n", &dumi[0], &dumi[1]);
    if ((dumi[0] != GRIDVER) || (dumi[1] != GRIDTYPE)) {
      fprintf(fail_file,
	      "expected gridver/gridtype: %d %d got %d %d\n",
	      GRIDVER, GRIDTYPE, dumi[0], dumi[1]);
      myexit(6);
    }
    ch = fgetc(in[l]);
    if (ch == '#') {
      while (fgetc(in[l]) != '\n');
    } else {
      fprintf(fail_file,
	      "1.1: grid data file doesn't match expected format: %d\n",
	      l);
      myexit(1);
    }

    LOOPF {
      fscanf(in[l], HEADER4_S,
	     &dumi[0], &dumi[1], &dumi[2], &dumi[3],
	     &dx[l][1][i], &dx[l][2][j], &dx[l][3][k],
	     &x[l][1][i], &x[l][2][j], &x[l][3][k],
	     &g[l][1][i], &dg[l][1][i], &g[l][2][i], &dg[l][2][i],
	     &g[l][3][i], &dg[l][3][i], &g[l][4][j], &dg[l][4][j],
	     &cotan[l][j], &dvl[l][1][i], &dvl[l][2][j], &dvl[l][3][k]);
      for (m = 1; m <= NUMSCA; m++) {
	fscanf(in[l], "%d ", &dumi[4 + (m - 1) * 3]);
	fscanf(in[l], "%d ", &dumi[5 + (m - 1) * 3]);
	fscanf(in[l], "%d ", &dumi[6 + (m - 1) * 3]);
      }
      for (m = 1; m <= NUMVEC; m++) {
	fscanf(in[l], "%d ", &dumi[4 + NUMSCA * 3 + (m - 1) * 3]);
	fscanf(in[l], "%d ", &dumi[5 + NUMSCA * 3 + (m - 1) * 3]);
	fscanf(in[l], "%d", &dumi[6 + NUMSCA * 3 + (m - 1) * 3]);
	if (m == NUMVEC)
	  fscanf(in[l], "\n");
	else
	  fscanf(in[l], " ");
      }
      // must do it this way to read in properly and get memory
      // writting correct too
      for (m = 1; m <= NUMSCA; m++) {
	bcs[m][1][k][j][i] = (short) dumi[4 + (m - 1) * 3];
	bcs[m][2][k][j][i] = (short) dumi[5 + (m - 1) * 3];
	bcs[m][3][k][j][i] = (short) dumi[6 + (m - 1) * 3];
      }
      for (m = 1; m <= NUMVEC; m++) {
	bcv[m][1][k][j][i] = (short) dumi[4 + NUMSCA * 3 + (m - 1) * 3];
	bcv[m][2][k][j][i] = (short) dumi[5 + NUMSCA * 3 + (m - 1) * 3];
	bcv[m][3][k][j][i] = (short) dumi[6 + NUMSCA * 3 + (m - 1) * 3];
      }

      if ((dumi[0] != l) || (dumi[1] != k) || (dumi[2] != j + startj)
	  || (dumi[3] != i)) {
	fprintf(fail_file,
		"2: grid data file doesn't match expected format: %d %d %d %d\n",
		l, k, j + startj, i);
	fprintf(fail_file, "%d %d %d %d\n", dumi[0], dumi[1],
		dumi[2], dumi[3]);
	myexit(1);
      }
    }
  }

  for (l = 1; l <= NUMGRID; l++) {
    fclose(in[l]);
  }

  return (0);
}



int init_rundat(int dump_cnt, int which)
{
  FILE *data;
  char fname[200];
  int dumi[50];
  SFTYPE dumfs[10];
  char ch;
  int i, j, k, l, m;
  FTYPE(*sin)[N3M][N2M][N1M];
  FTYPE(*vin)[3][N3M][N2M][N1M];
  char dfheader[MAXFILENAME];
  static FTYPE(*sigma)[3][N3M][N2M][N1M];	// -2*rho*nu*e_{ij}=
  // sigma_{ij}
  int gopp;

  // for now, assume filename is: dump for data

  // remember to make same changes to this and output in
  // diag.c/init.c!!
  sin=NULL;
  vin=NULL;

  // if post process call
  if (which == 9) {
    gopp = 1;
    which = 3;
  } else
    gopp = 0;

  // can only read in non-interp'ed data
  if (which == 0) {
    sin = s;
    vin = v;
    strcpy(dfheader, "dump");
  } else if (which == 50) {	// for post proc feature(doexpand)
    sin = s;
    vin = v;
    strcpy(dfheader, "toedump");
    which = 0;			// otherwise treat as dump
  } else if (which == 100) {
    sin = s;
    vin = v;
    strcpy(dfheader, "pdump");
    which = 0;			// otherwise treat as dump
  } else if (which == 1) {
    sin = sanal;
    vin = vanal;
    strcpy(dfheader, "adump");
  } else if (which == 2) {
    sin = floorvars;
    vin = floorvarv;
    // floorvar0 used natively
    strcpy(dfheader, "floor");
  } else if (which == 3) {
    sigma = workt1;		// assume workt1 holds sigma data in
    // all
    // cases
    // use workt1 as sigma(assume dump will use same)
    strcpy(dfheader, "npdump");
  } else {
    fprintf(fail_file, "not setup for which=%d in init_rundat\n",
	    which);
    myexit(1);
  }

  if (dump_cnt == -1) {
    sprintf(fname, "%s0_%s%s%s", DATADIR, dfheader, INEXT, myidtxt);
  } else if (dump_cnt >= 0) {
    if (which < 10) {
      sprintf(fname, "%s%s%04d%s%s", DATADIR, dfheader, dump_cnt,
	      DATEXT, myidtxt);
    } else if (which == 10) {	// no number for average file
      sprintf(fname, "%s%s%s%s", DATADIR, dfheader, DATEXT, myidtxt);
    }
  }
  fprintf(log_file, "proc:%02d: Start Reading in %s\n", myid, fname);

  if ((data = fopen(fname, "rt")) == NULL) {
    fprintf(fail_file, "data: Cannot open: %s\n", fname);
    myexit(1);
  }
  // check version info
  ch = fgetc(data);
  if (ch == '#') {
    while (fgetc(data) != '\n');
  } else {
    fprintf(fail_file,
	    "1:read: data file doesn't match expected format\n");
    myexit(1);
  }
  // version header
  fscanf(data, "%d %d\n", &dumi[0], &dumi[1]);

  if ((which == 0) || (which == 1)) {
    if ((dumi[0] != DVER) || (dumi[1] != DTYPE)) {
      fprintf(fail_file,
	      "read: Expected DVER/DTYPE %d %d got %d %d\n", DVER,
	      DTYPE, dumi[0], dumi[1]);
      myexit(6);
    }
  } else if (which == 2) {
    if ((dumi[0] != FLVER) || (dumi[1] != FLTYPE)) {
      fprintf(fail_file,
	      "read: Expected FLVER/FLTYPE %d %d got %d %d\n", FLVER,
	      FLTYPE, dumi[0], dumi[1]);
      myexit(6);
    }
  } else if (which == 3) {
    if ((dumi[0] != NPVER) || (dumi[1] != NPTYPE)) {
      fprintf(fail_file,
	      "read: Expected NPVER/NPTYPE %d %d got %d %d\n", NPVER,
	      NPTYPE, dumi[0], dumi[1]);
      myexit(6);
    }
  }
  // continue with data file and normal header
  ch = fgetc(data);
  if (ch == '#') {
    while (fgetc(data) != '\n');
  } else {
    fprintf(fail_file,
	    "1.1:read: data file doesn't match expected format\n");
    myexit(1);
  }
  if (OLDSCHOOL == 1) {
    if ((which == 0) || (which == 10))
      fscanf(data, INPUT1OLD, &t, &dumi[0]);
    else
      fscanf(data, INPUT1OLD, &dumfs[0], &dumi[0]);
  } else if (OLDSCHOOL == 0) {
    if ((which == 0) || (which == 10))
      fscanf(data, INPUT1, &t, &dumi[0], &dumi[1]);
    else
      fscanf(data, INPUT1, &dumfs[0], &dumi[0], &dumi[1]);
  }

  fprintf(log_file, "restarting at t=%15.10g\n", t);
  if (myid <= 0) {
    fprintf(logfull_file, "restarting at t=%15.10g\n", t);
  }

  if (dumi[0] != 0) {
    fprintf(log_file,
	    "1.5: Warning: cannot take in interpolated data: expected: 0 got: %d\n",
	    dumi[0]);
    if (myid <= 0) {
      fprintf(logfull_file,
	      "1.5: Warning: cannot take in interpolated data: expected: 0 got: %d\n",
	      dumi[0]);
    }
    // myexit(1);
  }
  if (OLDSCHOOL == 0) {
    if (dumi[1] != 0) {
      fprintf(log_file,
	      "1.6: Warning, data is not gridded correctly: expected: 0 got: %d\n",
	      dumi[1]);
      if (myid <= 0) {
	fprintf(logfull_file,
		"1.6: Warning, data is not gridded correctly: expected: 0 got: %d\n",
		dumi[1]);
      }
    }
  }
  while (fgetc(data) != '\n');	// skip to next line

  ch = fgetc(data);
  if (ch == '#') {
    while (fgetc(data) != '\n');	// skip comment line
  } else {
    fprintf(fail_file, "2: data file doesn't match expected format\n");
    myexit(1);
  }

  fscanf(data, "%d %d %d", &dumi[0], &dumi[1], &dumi[2]);
  if ((dumi[0] != N1) || (dumi[1] != N2) || (dumi[2] != N3)) {
    fprintf(fail_file, "2.4: wrong data size\n");
    myexit(1);
  }
  while (fgetc(data) != '\n');	// skip to next line



  ch = fgetc(data);
  if (ch == '#') {
    while (fgetc(data) != '\n');	// skip comment line
  } else {
    fprintf(fail_file,
	    "2.5: data file doesn't match expected format\n");
    myexit(1);
  }

#if(FULLINPUT)
  LOOPF
#else
  LOOP
#endif
  {
    if (which < 3) {
      for (l = 1; l <= NUMSCA; l++) {
	fscanf(data, INPUT2, &sin[l][k][j][i]);
      }
      // only input first vector(v) for now
      for (l = 1; l <= NUMVEC - 1; l++) {
	if (which == 2)
	  fscanf(data, INPUT2, &floorvar0[l][k][j][i]);
	for (m = 1; m <= 3; m++) {
	  fscanf(data, INPUT3I, &vin[l][m][k][j][i]);
	  if ((l == NUMVEC - 1) && (m == 3))
	    fscanf(data, "\n");
	  else
	    fscanf(data, " ");
	}
      }
    }				// end if which<3
    else if (which == 3) {	// sigma read(only used during post
      // process)
      for (l = 1; l <= 3; l++) {
	for (m = 1; m <= 3; m++) {
	  fscanf(data, INPUT2, &sigma[l][m][k][j][i]);
	}
      }
    }
  }				// end over loop of data

  fclose(data);

  // check inputted data for correctness
  LOOP {
    if (s[1][k][j][i] < 1E-20) {
      fprintf(fail_file,
	      "input data has mass density <1E-20: %d %d %d %15.10g\n",
	      k, j, i, s[1][k][j][i]);
      myexit(1);
    }
    if (s[2][k][j][i] < 1E-20) {
      fprintf(fail_file,
	      "input data has ie density <1E-20: %d %d %d %15.10g\n",
	      k, j, i, s[2][k][j][i]);
      myexit(1);
    }
  }

  fprintf(log_file, "Done Reading in %s file\n", dfheader);
  fflush(log_file);
  return (0);
}

// mapping from image space to function space
#define CMAPSCA1(x) pow(10.0,x)
#define CMAPSCA2(x) pow(10.0,x)
#define CMAPSCA3(x) (x)
#define CMAPSCAGEN(x) (x)

#define CMAPVEC1(x) (x)
#define CMAPVEC2(x) (x)

// mapping from function to image space in value(color)
#define FMAPSCA1(x) log10(x)
#define FMAPSCA2(x) log10(x)
#define FMAPSCA3(x) (x)
#define FMAPSCAGEN(x) (x)

#define FMAPVEC1(x) (x)
#define FMAPVEC2(x) (x)

#define LINLOGFIX 0		// whether to use both linear and log10
				// input data to recreate old data

// used to input image data(much extracted from image()
// assume one doesn't care about data scaling, just used to interpolate
// during post process run.  Only need to invert color mapping for ppm
int init_runimage(int im_cnt, int wsca, int wvec, int call_code,
		  int outtype)
{
  int i, j, l;
  char ifheader[MAXFILENAME];
  SFTYPE liq, a=0, b=0, c=0;
  unsigned char tempuch;
  int iii, dualgo;
  int ll;
  int q;
  SFTYPE ftempfix=0;

  FILE *im_file;
  static char ifnam[MAXFILENAME], temp[MAXFILENAME];

  int im_cnts[ITYPES][NUMSCA + 1];
  int im_cntv[ITYPES][NUMVEC + 1];

  char temps[MAXFILENAME];

  // for now assume all im_cnt same
  for (i = 1; i <= NUMSCA; i++) {
    im_cnts[outtype][i] = im_cnt;
  }
  for (i = 1; i <= NUMVEC; i++) {
    im_cntv[outtype][i] = im_cnt;
  }


  // ///////// SCALARS
  if (wsca != 0) {
    for (l = 1; l <= NUMSCA; l++) {
      /* if not to do all, pick */
      if (wsca != -1) {
	ll = wsca;
      } else
	ll = l;

      // assume samplei originally 0, no interp from interp'ed
      // data

      dualgo = LINLOGFIX;	// need to read in linear for linear
      // range 
      // and log for lower detail range.
      for (iii = 0; iii <= dualgo; iii++) {
	// setup file input
	sprintf(temps, "%s%s", DATADIR, "i/");
	strcpy(ifheader, "imx");
	sprintf(temps, "%simx%01d-%01d-s%01d/", temps, outtype,
		iii, ll);
	sprintf(ifnam, "%s%s%01d-%01d-s%01d-%04d%s%s", temps,
		ifheader, outtype, iii, ll, im_cnts[outtype][ll],
		DATEXT, myidtxt);

	if (IMAGEFORMATINPUT == 0) {
	  strcat(ifnam, ".r8");
	  if (GZIPIMAGEINPUT == 0) {
	    im_file = fopen(ifnam, "rb");
	  } else if (GZIPIMAGEINPUT > 0) {
	    strcat(ifnam, ".gz");

	    sprintf(temp, "gzip -d < %s", ifnam);
	    im_file = popen(temp, "r");
	    // fprintf(stdout,"%s %s\n",ifnam,temp);
	  }
	  if (im_file == NULL) {
	    fprintf(fail_file, "error opening image file: %s\n", ifnam);
	    myexit(2);
	  }
	  // skip 4 comment lines
#if(OLDSCHOOL2==0)
	  while (fgetc(im_file) != '\n');	// skip to next
	  // line
	  while (fgetc(im_file) != '\n');	// skip to next
	  // line
	  while (fgetc(im_file) != '\n');	// skip to next
	  // line
	  while (fgetc(im_file) != '\n');	// skip to next
	  // line
#endif
	}
	if (IMAGEFORMATINPUT == 1) {
	  fprintf(fail_file,
		  "Can't input this format of image, since not simple to reverse lookup function value to interpolate\n");
	  myexit(1);
	}
	// now set function map using min/max
	// same code as in image()
	if (ll == 1) {
	  if (iii == 0) {
	    b = FMAPSCA1(mms[outtype][iii][ll][0]);
	    a = 255. / (FMAPSCA1(mms[outtype][iii][ll][1]) - b);
	  }
	  if (iii == 1) {
	    b = FMAPSCAGEN(mms[outtype][iii][ll][0]);
	    a = 255. / (FMAPSCAGEN(mms[outtype][iii][ll][1]) - b);
	  }
	} else if (ll == 2) {
	  if (iii == 0) {
	    b = FMAPSCA2(mms[outtype][iii][ll][0]);
	    a = 255. / (FMAPSCA2(mms[outtype][iii][ll][1]) - b);
	  }
	  if (iii == 1) {
	    b = FMAPSCAGEN(mms[outtype][iii][ll][0]);
	    a = 255. / (FMAPSCAGEN(mms[outtype][iii][ll][1]) - b);
	  }
	} else if (ll == 3) {
	  b = FMAPSCA3(mms[outtype][iii][ll][0]);
	  a = 255. / (FMAPSCA3(mms[outtype][iii][ll][1]) - b);
	}
	c = -a * b;

	LOOPINI {
	  /* read value */
	  fread(&tempuch, sizeof(unsigned char), 1, im_file);
	  // tempuch=(unsigned char)fgetc(im_file);
	  liq = (SFTYPE) ((int) (tempuch)) / a + b;	// inverse 
	  // 
	  // now invert from color to function space
	  if (ll == 1) {
	    if (iii == 0) {
	      ftempfix = CMAPSCA1(liq);
	    }
	    if (iii == 1) {
	      ftempfix = CMAPSCAGEN(liq);
	    }
	  } else if (ll == 2) {
	    if (iii == 0) {
	      ftempfix = CMAPSCA2(liq);
	    }
	    if (iii == 1) {
	      ftempfix = CMAPSCAGEN(liq);
	    }
	  } else if (ll == 3) {
	    ftempfix = CMAPSCA3(liq);
	  }
	  liq = ftempfix;

	  if (LINLOGFIX) {
	    // only modify in a way which uses the linear
	    // and
	    // log10 data best
	    if ((iii == 0)
		&& (liq < mms[outtype][iii][ll][1] / (256.0))) {
	      s[ll][0][j][i] = liq;	// assume interp
	      // only uses k=0
	    }
	    if ((iii == 1)
		&& (liq >= mms[outtype][iii][ll][1] / (256.0))) {
	      s[ll][0][j][i] = liq;	// assume interp
	      // only uses k=0
	    }
	  } else {
	    if (iii == 0) {
	      s[ll][0][j][i] = liq;
	    }
	  }
	}

	if (GZIPIMAGEINPUT == 0) {
	  fclose(im_file);
	} else if (GZIPIMAGEINPUT > 0) {
	  pclose(im_file);
	}
      }				// over dualgo
      /* cut short loop if only to do one */
      if (wsca != -1)
	l = NUMSCA;
    }				// over scalars
  }				// if any scalars




  // ////////// VECTORS
  if (wvec != 0) {
    for (l = 1; l <= NUMVEC - 1; l++) {
      /* if not to do all, pick */
      if (wvec != -1) {
	ll = wvec;
      } else
	ll = l;


      for (q = 1; q <= 3; q++) {	// over components (no need to
	// do
	// q=0)


	dualgo = 0;		// force for input, no need for other
	// data 
	// since can just deduce this from linear
	// rho/v in pp
	for (iii = 0; iii <= dualgo; iii++) {
	  // setup file output
	  sprintf(temps, "%s%s", DATADIR, "i/");

	  strcpy(ifheader, "imx");
	  sprintf(temps, "%simx%01d-%01d-v%01d-%01d/", temps,
		  outtype, iii, ll, q);
	  sprintf(ifnam, "%s%s%01d-%01d-v%01d-%01d-%04d%s%s",
		  temps, ifheader, outtype, iii, ll, q,
		  im_cnts[outtype][ll], DATEXT, myidtxt);
	  if (IMAGEFORMATINPUT == 0) {
	    strcat(ifnam, ".r8");
	    if (GZIPIMAGEINPUT == 0) {
	      im_file = fopen(ifnam, "rb");
	    } else if (GZIPIMAGEINPUT > 0) {
	      strcat(ifnam, ".gz");
	      sprintf(temp, "gzip -d < %s", ifnam);
	      strcpy(ifnam, temp);	// for below
	      // fprintf
	      im_file = popen(ifnam, "r");
	    }
	    if (im_file == NULL) {
	      fprintf(fail_file,
		      "error opening image file: %s\n", ifnam);
	      myexit(2);
	    }
	    // skip 4 comment lines
#if(OLDSCHOOL2==0)
	    while (fgetc(im_file) != '\n');	// skip to next
	    // line
	    while (fgetc(im_file) != '\n');	// skip to next
	    // line
	    while (fgetc(im_file) != '\n');	// skip to next
	    // line
	    while (fgetc(im_file) != '\n');	// skip to next
	    // line
#endif
	  }
	  if (IMAGEFORMATINPUT == 1) {
	    fprintf(fail_file,
		    "Can't input this format of image, since not simple to reverse lookup function value to interpolate\n");
	    myexit(1);
	  }
	  // setup map
	  // same code as in image()
	  if (ll == 1) {
	    if (iii == 0) {
	      // warning: silent fail if 0,1 are both 0
	      // for
	      // q=0 ?
	      b = FMAPVEC1(mmv[outtype][iii][ll][q][0]);
	      a = 255. / (FMAPVEC1(mmv[outtype][iii][ll][q][1]) - b);
	    }
	    if (iii == 1) {	// rho*v // only good if v
	      // range
	      // is - to + values
	      b = FMAPVEC1(mmv[outtype][iii][ll][q][0]);
	      a = 255. / (FMAPVEC1(mmv[outtype][iii][ll][q][1]) - b);
	    }
	  }
	  if (ll == 2) {
	    if (iii == 0) {	// B
	      b = FMAPVEC2(mmv[outtype][iii][ll][q][0]);
	      a = 255. / (FMAPVEC2(mmv[outtype][iii][ll][q][1]) - b);
	    }
	    if (iii == 1) {	// B*v // only good if v range
	      // is
	      // - to + values
	      b = FMAPVEC2(mmv[outtype][iii][ll][q][0]);
	      a = 255. / (FMAPVEC2(mmv[outtype][iii][ll][q][1]) - b);
	    }
	  }
	  c = -a * b;

	  // input image without map
	  LOOPINI {

	    fread(&tempuch, sizeof(unsigned char), 1, im_file);
	    liq = (SFTYPE) ((int) (tempuch)) / a + b;	// inverse 
	    // 
	    // value

	    // invert
	    if (ll == 1) {
	      ftempfix = CMAPVEC1(liq);
	    } else if (ll == 2) {
	      ftempfix = CMAPVEC2(liq);
	    }
	    liq = ftempfix;

	    v[ll][q][0][j][i] = liq;	// assume interp
	    // only uses k=0 
	  }


	  if (GZIPIMAGEINPUT == 0) {
	    fclose(im_file);
	  } else if (GZIPIMAGEINPUT > 0) {
	    pclose(im_file);
	  }
	}			// dualgo
      }				// over components of the vector
      /* cut short loop if only to do one */
      if (wvec != -1)
	l = NUMVEC - 1;
    }				// over vectors
  }				// if any vectors
  return (0);
}



// which: -1: all 0: params 1: grid no interp 2: grid interp 3: interp
// params
int init_outgparm(int which)
{
  int i, j, k, l, m;
  FILE *out;
  char filename[MAXFILENAME];
  char temps[MAXFILENAME];

  // output global parameters
  if (myid <= 0) {
    if ((which == -1) || (which == 0)) {
      sprintf(filename, "%s0_gparam%s", DATADIR, PAREXT);

      if ((out = fopen(filename, "wt")) == NULL) {
	fprintf(fail_file, "Can't open %s for writting\n", filename);
	myexit(1);
      } else {
	fprintf(out, HEADER3_P, "PVER", "TYPE", PVER, PTYPE, "N1",
		"N2", "N3", N1, totalsize[2], N3, "x1-start",
		"x2-start", "x3-start", "x1-length", "x2-length",
		"x3-length", L[1][1], L[1][2], L[1][3], L[2][1],
		L[2][2], L[2][3]
		, "rg", "rgp", "cs", "coolfact", "gam",
		"alpha_real0", "n_real", rg, rgp, cs, coolfact,
		gam, alpha_real0, n_real, "nu_vnr", "nu_l",
		"nu_ten", "cour", "cour2", nu_vnr, nu_l, nu_ten,
		cour, cour2, "GRAVC", "MASSBH", GRAVC, MASSBH,
		"tstart", "tf", "tavgi", "tavgf", "numavg",
		"timagescale", tstart, tf, tf, tf, tf,
		timagescale, "DTl", "DTd", "DTi", "DTloss",
		"DTfloor", "DTtimestep", "DTpd", "DTener", DTl,
		DTd, DTi, DTloss, DTfloor, DTtimestep, DTpd,
		DTener, "resist", "nu_sh", resist, nu_sh, "vg-x1",
		"vg-x2", "vg-x3", vg[1], vg[2], vg[3]);
      }

      fclose(out);
    }
  }
  // always need general grid data
  // output dx x diffs to file--active and ghost zones
  if ((which == -1) || (which == 1)) {
    for (l = 1; l <= NUMGRID; l++) {
      sprintf(temps, "%s0_grid", DATADIR);

      sprintf(filename, "%s%01d%s%s", temps, l, PAREXT, myidtxt);


      if ((out = fopen(filename, "wt")) == NULL) {
	fprintf(fail_file, "Can't open %s for writting\n", filename);
	myexit(1);
      } else {
	fprintf(out, "#%10s\n%10d %10d\n", "GRIDVER", GRIDVER,
		GRIDTYPE);
	fprintf(out,
		"#%4s %4s %4s %4s " "%21s %21s %21s "
		"%21s %21s %21s "
		"%21s %21s %21s %21s %21s %21s %21s %21s %21s "
		"%21s %21s %21s ", "grid", "k", "j", "i", "dx1",
		"dx2", "dx3", "x1", "x2", "x3", "g1", "dg1", "g2",
		"dg2", "g3", "dg3", "g4", "dg4", "cot", "dvl1",
		"dvl2", "dvl3");
	for (m = 1; m <= NUMSCA; m++) {
	  fprintf(out, "%6s%02d %6s%02d %6s%02d ", "bsty", m,
		  "bsdm", m, "bsdr", m);
	}
	for (m = 1; m <= NUMVEC; m++) {
	  fprintf(out, "%6s%02d %6s%02d %6s%02d", "bvty", m,
		  "bvdm", m, "bvdr", m);
	  if (m == NUMVEC)
	    fprintf(out, "\n");
	  else
	    fprintf(out, " ");
	}

	LOOPF {
	  fprintf(out, HEADER4_P,
		  l, k, startj + j, i,
		  dx[l][1][i], dx[l][2][j], dx[l][3][k],
		  x[l][1][i], x[l][2][j], x[l][3][k],
		  g[l][1][i], dg[l][1][i], g[l][2][i],
		  dg[l][2][i], g[l][3][i], dg[l][3][i],
		  g[l][4][j], dg[l][4][j], cotan[l][j],
		  dvl[l][1][i], dvl[l][2][j], dvl[l][3][k]
	      );
	  for (m = 1; m <= NUMSCA; m++) {
	    fprintf(out, "%6d ", bcs[m][1][k][j][i]);
	    fprintf(out, "%6d ", bcs[m][2][k][j][i]);
	    fprintf(out, "%6d ", bcs[m][3][k][j][i]);
	  }
	  for (m = 1; m <= NUMVEC; m++) {
	    fprintf(out, "%6d ", bcv[m][1][k][j][i]);
	    fprintf(out, "%6d ", bcv[m][2][k][j][i]);
	    fprintf(out, "%6d", bcv[m][3][k][j][i]);
	    if (m == NUMVEC)
	      fprintf(out, "\n");
	    else
	      fprintf(out, " ");
	  }
	}
      }
      fclose(out);
    }

    if (FULLOUTPUT == 0) {	// only need if not doing fulloutput
      // output dx x diffs to file--active zones only (only diff
      // is
      // LOOPF->LOOP from above)
      for (l = 1; l <= NUMGRID; l++) {
	sprintf(temps, "%s0_gridact", DATADIR);

	sprintf(filename, "%s%01d%s%s", temps, l, PAREXT, myidtxt);

	if ((out = fopen(filename, "wt")) == NULL) {
	  fprintf(fail_file, "Can't open %s for writting\n", filename);
	  myexit(1);
	} else {
	  fprintf(out, "#%10s\n%10d %10d\n", "GRIDVER", GRIDVER,
		  GRIDTYPE);
	  fprintf(out,
		  "#%4s %4s %4s %4s " "%21s %21s %21s "
		  "%21s %21s %21s "
		  "%21s %21s %21s %21s %21s %21s %21s %21s %21s "
		  "%21s %21s %21s ", "grid", "k", "j", "i",
		  "dx1", "dx2", "dx3", "x1", "x2", "x3", "g1",
		  "dg1", "g2", "dg2", "g3", "dg3", "g4", "dg4",
		  "cot", "dvl1", "dvl2", "dvl3");
	  for (m = 1; m <= NUMSCA; m++) {
	    fprintf(out, "%6s%02d %6s%02d %6s%02d ", "bsty", m,
		    "bsdm", m, "bsdr", m);
	  }
	  for (m = 1; m <= NUMVEC; m++) {
	    fprintf(out, "%6s%02d %6s%02d %6s%02d", "bvty", m,
		    "bvdm", m, "bvdr", m);
	    if (m == NUMVEC)
	      fprintf(out, "\n");
	    else
	      fprintf(out, " ");
	  }

	  LOOP {
	    fprintf(out, HEADER4_P,
		    l, k, startj + j, i,
		    dx[l][1][i], dx[l][2][j], dx[l][3][k],
		    x[l][1][i], x[l][2][j], x[l][3][k],
		    g[l][1][i], dg[l][1][i], g[l][2][i],
		    dg[l][2][i], g[l][3][i], dg[l][3][i],
		    g[l][4][j], dg[l][4][j], cotan[l][j],
		    dvl[l][1][i], dvl[l][2][j], dvl[l][3][k]
		);
	    for (m = 1; m <= NUMSCA; m++) {
	      fprintf(out, "%6d ", bcs[m][1][k][j][i]);
	      fprintf(out, "%6d ", bcs[m][2][k][j][i]);
	      fprintf(out, "%6d ", bcs[m][3][k][j][i]);
	    }
	    for (m = 1; m <= NUMVEC; m++) {
	      fprintf(out, "%6d ", bcv[m][1][k][j][i]);
	      fprintf(out, "%6d ", bcv[m][2][k][j][i]);
	      fprintf(out, "%6d", bcv[m][3][k][j][i]);
	      if (m == NUMVEC)
		fprintf(out, "\n");
	      else
		fprintf(out, " ");
	    }
	  }
	}
	fclose(out);
      }
    }
  }

  return (0);
}



// 1/2/3-d valid
int init_pointers(void)
{

  dx = (FTYPE(*)[3][NBIGM]) (&(dxa[-1][-1][NBIGBND]));
  x = (FTYPE(*)[3][NBIGM]) (&(xa[-1][-1][NBIGBND]));

  idx = (FTYPE(*)[3][INTNBIG]) (&(idxa[-1][-1][0]));
  ix = (FTYPE(*)[3][INTNBIG]) (&(ixa[-1][-1][0]));

  s = (FTYPE(*)[N3M][N2M][N1M]) (&(sa[-1][N3BND][N2BND][N1BND]));
  v = (FTYPE(*)[3][N3M][N2M][N1M]) (&(va[-1][-1][N3BND][N2BND][N1BND]));

  accountstore =
      (short (*)[N2M][N1M]) (&(accountstorea[N3BND][N2BND][N1BND]));

  g = (FTYPE(*)[NUMMETRIC][NBIGM]) (&(ga[-1][-1][NBIGBND]));
  dg = (FTYPE(*)[NUMMETRIC][NBIGM]) (&(dga[-1][-1][NBIGBND]));
  dvl = (FTYPE(*)[3][NBIGM]) (&(dvla[-1][-1][NBIGBND]));
  cotan = (FTYPE(*)[NBIGM]) (&(cotana[-1][NBIGBND]));


  DS = (FTYPE(*)[3][N3M][N2M][N1M]) (&
				     (DSa[-1][-1][N3BND][N2BND]
				      [N1BND]));
  OARCL = (FTYPE(*)[3][N3M][N2M][N1M]) (&(OARCLa[-1][-1][N3BND][N2BND]
					  [N1BND]));
  OVOL = (FTYPE(*)[N3M][N2M][N1M]) (&(OVOLa[-1][N3BND][N2BND][N1BND]));


  bcs = (short (*)[3][N3M][N2M][N1M]) (&(bcsa[-1][-1][N3BND][N2BND]
					 [N1BND]));
  bcv = (short (*)[3][N3M][N2M][N1M]) (&(bcva[-1][-1][N3BND][N2BND]
					 [N1BND]));

  worksbc = (FTYPE(*)[2][N2BND * N1M]) (&(worksbca[-1][-1][0]));
  workvbc = (FTYPE(*)[2][3 * N2BND * N1M]) (&(workvbca[-1][-1][0]));


  sanal =
      (FTYPE(*)[N3M][N2M][N1M]) (&(sanala[-1][N3BND][N2BND][N1BND]));
  vanal = (FTYPE(*)[3][N3M][N2M][N1M]) (&(vanala[-1][-1][N3BND][N2BND]
					  [N1BND]));

  floorvars =
      (FTYPE(*)[N3M][N2M][N1M]) (&
				 (floorvarsa[-1][N3BND][N2BND][N1BND]));
  floorvarv =
      (FTYPE(*)[3][N3M][N2M][N1M]) (&(floorvarva[-1][-1][N3BND][N2BND]
				      [N1BND]));
  floorvar0 =
      (FTYPE(*)[N3M][N2M][N1M]) (&
				 (floorvar0a[-1][N3BND][N2BND][N1BND]));

  losss = (SFTYPE(*)[3][2][NBIG]) (&(losssa[-1][-1][0][0]));
  lossv = (SFTYPE(*)[3 + 1][3][2][NBIG]) (&(lossva[-1][0][-1][0][0]));
  lossvisc = (SFTYPE(*)[3][2][NBIG]) (&(lossvisca[-1][-1][0][0]));

  nu_fact = (FTYPE(*)[N2M][N1M]) (&(nu_facta[N3BND][N2BND][N1BND]));
  nu_real = (FTYPE(*)[N2M][N1M]) (&(nu_reala[N3BND][N2BND][N1BND]));


  rhoinject = (FTYPE(*)[N2M][N1M]) (&(rhoinjecta[N3BND][N2BND][N1BND]));
  eninject = (FTYPE(*)[N2M][N1M]) (&(eninjecta[N3BND][N2BND][N1BND]));

  work1 = (FTYPE(*)[N2M][N1M]) (&(work1a[N3BND][N2BND][N1BND]));
  work2 = (FTYPE(*)[N2M][N1M]) (&(work2a[N3BND][N2BND][N1BND]));
  work3 = (FTYPE(*)[N2M][N1M]) (&(work3a[N3BND][N2BND][N1BND]));
  work4 = (FTYPE(*)[N2M][N1M]) (&(work4a[N3BND][N2BND][N1BND]));
  work5 = (FTYPE(*)[N2M][N1M]) (&(work5a[N3BND][N2BND][N1BND]));
  work6 = (FTYPE(*)[N2M][N1M]) (&(work6a[N3BND][N2BND][N1BND]));
  work7 = (FTYPE(*)[N2M][N1M]) (&(work7a[N3BND][N2BND][N1BND]));
  work8 = (FTYPE(*)[N2M][N1M]) (&(work8a[N3BND][N2BND][N1BND]));
  work9 = (FTYPE(*)[N2M][N1M]) (&(work9a[N3BND][N2BND][N1BND]));
  work10 = (FTYPE(*)[N2M][N1M]) (&(work10a[N3BND][N2BND][N1BND]));

  works1 = (SFTYPE(*)[N2M][N1M]) (&(works1a[N3BND][N2BND][N1BND]));

  workv1 =
      (FTYPE(*)[N3M][N2M][N1M]) (&(workv1a[-1][N3BND][N2BND][N1BND]));
  workv2 =
      (FTYPE(*)[N3M][N2M][N1M]) (&(workv2a[-1][N3BND][N2BND][N1BND]));
  workv3 =
      (FTYPE(*)[N3M][N2M][N1M]) (&(workv3a[-1][N3BND][N2BND][N1BND]));
  workv4 =
      (FTYPE(*)[N3M][N2M][N1M]) (&(workv4a[-1][N3BND][N2BND][N1BND]));
  workv5 =
      (FTYPE(*)[N3M][N2M][N1M]) (&(workv5a[-1][N3BND][N2BND][N1BND]));

  workt1 = (FTYPE(*)[3][N3M][N2M][N1M]) (&(workt1a[-1][-1][N3BND][N2BND]
					   [N1BND]));
  workt2 = (FTYPE(*)[3][N3M][N2M][N1M]) (&(workt2a[-1][-1][N3BND][N2BND]
					   [N1BND]));
  workt3 = (FTYPE(*)[3][N3M][N2M][N1M]) (&(workt3a[-1][-1][N3BND][N2BND]
					   [N1BND]));

  workiq = (FTYPE(*)[INTN2][INTN1]) (&(workiqa[-1][0][0]));
  work0iq = (FTYPE(*)[INTN2][INTN1]) (&(work0iqa[-1][0][0]));
  // below big waste of memory
  workviq = (FTYPE(*)[3][INTN2][INTN1]) (&(workviqa[-1][-1][0][0]));


  // ptraddr(-1);

  // END: Setup arrays giving NBND boundary zones--makes 0,0 first
  // active zone 
  return (0);
}


// Notes on dx(dvl):
// -----
// Without any interpolation, must setup dx to be consistent at
// boundary
// zones with spatial structure.
// So I use BC interpolation to avoid worrying about this :) Kinda
// ingore
// below
// ------
// Periodic-> spatial condition.  dx must be periodic->
// dx(ghost)=dx(what
// copied), dx as scalar
// Reflect-> spatial condition. dx(-1)=dx(0), etc., dx as scalar
// AOS-> as with Reflect, dx as scalar, do not invert dx->-dx for
// z-axis
// 
// Rest do not give spatial condition, all spatial info indeterminate.
// Outflow: Soft Edge: ~Neuman, Like bigger domain existed or true
// outflow 
// to edge.
// Inflow: Hard Edge: Dirichlet: or not really hard edge, like bigger
// domain but with influx of matter.

// Outflow using copy needs Reflect-like spatial condition(i.e. local
// copy 
// of dx to boundary)
// Outflow using interpolation using position does not care what dx is
// on
// bc.  Probably reasonable to make sure dx does not put middle of bc
// zone 
// at x=0 for non-cart coords.

// Inflow is just Dirichlet, does not care what dx is, as with outflow
// interpolated.  Make sure to use correct position in
// condition(center/edge of zone, etc.).

// Below I have chosen out/in-flow to be treated same as with reflect
// for
// dx/dvl unless I decide otherwise
// Might be interesting to set dx large for outflow, as if domain
// existed, 
// but only as average value.  Probably need v to be linearly
// interpolated 
// within zone to be really meaninglful.

// Necessary dxs:
// besides active grid, must at least have dx on:
// for advection terms(for fluxes(sweepy only)) and dq for sweepx:
// advection(sweepy):
// dx1a(i) <> i=0..N1-1, j=-1..N2
// dx1b(i) <> i=0..N1-1, j=0..N2
// dq(sweepx):
// dx1b(i) <> i=0..N1, j=0..N2-1
// dx1a(i) <> i=-1..N1-1, j=0..N2-1
// dq(sweepy):
// dx2b(i) <> i=0..N1-1, j=0..N2
// dx2a(i) <> i=0..N1-1, j=-1..N2-1


// 1/2/3-d valid
// which: 0 = real grid 1=interp grid
// outtype: for interp grid only: what kind of grid:
// 0=native type
// 1=linear cart
int init_dx(int n1, int n2, int n3, int n1b, int n2b, int n3b,
	    SFTYPE * sx, int which, int outtype)
{
  SFTYPE frac;
  int i, j, k, m;
  int N[3 + 1], NB[3 + 1];
  int NC;
  SFTYPE dellogx;
  SFTYPE cc, cw, dx1=0, dx2=0, Ne1, Ne2, dxe, dx0=0;
  SFTYPE bcoef, acoef, ccoef, dcoef, truea;
  SFTYPE acoefi;
  SFTYPE acoeff;
  SFTYPE ftemp=0;
  SFTYPE xasofar;
  SFTYPE ftemp1, ftemp2;
  int itemp;
  int realsize2=0;


  N[1] = n1;
  N[2] = n2;
  N[3] = n3;
  NB[1] = n1b;
  NB[2] = n2b;
  NB[3] = n3b;
  /* BEGIN: Apply global parameters to cell structure */


  /* Assign dxs--uniform grid */
  /* Note this data would be dynamic if dynamic grid */

  // a-grid


  // X1
  // add up sizes for total size of grid
  if (which == 0)
    totalsize[1] = N[1];
  else if (which == 1)
    itotalsize[1] = N[1];

  if ((which == 0) || (outtype == 0)) {

    if (which == 1) {
      iL[1][1] = L[1][1];
      iL[2][1] = L[2][1];
    }

    sx[1] = L[1][1];
    if (nonunigridx1 == 0) {
      for (i = -NB[1]; i < N[1] + NB[1]; i++) {
	ftemp = L[2][1] / N[1];
	if (which == 0)
	  dx[1][1][i] = ftemp;
	else if (which == 1)
	  idx[1][1][i] = ftemp;
      }
    } else if (nonunigridx1 == 1) {
      // use mathematica page to get these values or pick from
      // domain determined below:
      NC = N[1] / 3;		// specify how many zones for cosine
      // region
      if (NC != N[1]) {
	cw = 1.0;
	cc = 1.0;
	dx1 = (cc - cw * .5) - L[1][1];
	dx2 = L[2][1] + L[1][1] - (cs + cw * .5);
	dxe = (dx1 + dx2) / ((SFTYPE) (N[1] - NC));
	Ne1 = rint(dx1 / dxe);
	Ne2 = rint(dx2 / dxe);
	if (cw > L[2][1]) {
	  fprintf(fail_file,
		  "Can't have coswidth greater than real width\n");
	  myexit(1);
	}
	if ((cc - cw / 2 < L[1][1])
	    || (cc + cw / 2 > L[1][1] + L[2][1])) {
	  fprintf(fail_file,
		  "Can't have edges of cos outside domain!\n");
	  myexit(1);
	}
	// partition rest of space evenly
	for (i = -NB[1]; i < (int) (Ne1); i++) {
	  ftemp = dxe;
	  if (which == 0)
	    dx[1][1][i] = ftemp;
	  else if (which == 1)
	    idx[1][1][i] = ftemp;
	}
	for (i = N[1] - 1 + NB[1]; i >= N[1] - (int) (Ne2); i--) {
	  ftemp = dxe;
	  if (which == 0)
	    dx[1][1][i] = ftemp;
	  else if (which == 1)
	    idx[1][1][i] = ftemp;
	}

	// no acoef here(use l21-dx2-dx1 instead of cw for fix
	// on
	// aliasing accuracy
	bcoef =
	    ((SFTYPE) (NC) * dxe -
	     (L[2][1] - dx1 - dx2)) / ((SFTYPE) (NC - 1.));
	// fprintf(stderr,"%15.10g %15.10g %15.10g %15.10g
	// %15.10g 
	// %15.10g\n",dx1,dx2,dxe,Ne1,Ne2,bcoef);
	if (bcoef < 0) {
	  fprintf(fail_file,
		  "bcoef<0 does not produce intended results!\n");
	  myexit(1);
	}
	if (bcoef > dxe * 0.5) {
	  fprintf(fail_file,
		  "bcoef>dxe/2 does not produce intended results!\n");
	  myexit(1);
	}


	for (i = Ne1; i < N[1] - Ne2; i++) {
	  ftemp =
	      (dxe - bcoef) +
	      bcoef * cos(2. * M_PI * (SFTYPE) (i - Ne1) /
			  ((SFTYPE) (NC - 1)));
	  if (which == 0)
	    dx[1][1][i] = ftemp;
	  else if (which == 1)
	    idx[1][1][i] = ftemp;
	}
      } else {
	acoeff = (L[2][1] / ((SFTYPE) N[1]));
	acoefi = (L[2][1] / ((SFTYPE) (N[1] + 1)));
	// must have a>b for nonnegative dx1
	acoef = (acoeff - acoefi) * .1 + acoefi;	// add
	// some%
	// off
	// inner
	// acoef
	// to get
	// smallest 
	// dx in
	// center
	bcoef = L[2][1] - (SFTYPE) (N[1]) * acoef;

	for (i = -NB[1]; i < N[1] + NB[1]; i++) {
	  ftemp =
	      acoef +
	      bcoef * cos(2. * M_PI * i / ((SFTYPE) (N[1] - 1)));
	  if (which == 0)
	    dx[1][1][i] = ftemp;
	  else if (which == 1)
	    idx[1][1][i] = ftemp;
	}
      }
    } else if (nonunigridx1 == 2) {
      ccoef = 0.0;		// leave at 0.0
      // bcoef starts at >3.000000 and goes to
      // infinity(100.0~uniform) 
      bcoef = 3.1;		// free
      acoef =
	  (L[2][1] - ccoef * N[1]) / (alnfact(bcoef + N[1] - 2) -
				      alnfact(bcoef - 2));

      for (i = -NB[1]; i < N[1] + NB[1]; i++) {
	ftemp = ccoef + acoef * log(i + bcoef);
	if (which == 0)
	  dx[1][1][i] = ftemp;
	else if (which == 1)
	  idx[1][1][i] = ftemp;
      }
    } else if ((nonunigridx1 == 3) || (nonunigridx1 == 4)) {
      acoef = 10.0;		// Factor by which each larger decade
      // dx
      // is to last decade.
      ftemp1 = L[1][1] + L[2][1];
      ftemp2 = L[1][1];
      bcoef = (SFTYPE) (N[1]) / (log10(ftemp1) - log10(ftemp2));	// Number 
									// 
      // 
      // of 
      // zones 
      // per 
      // decade

      ccoef =
	  (pow(acoef, (SFTYPE) (N[1]) / (bcoef - 1.)) -
	   1.) / (pow(acoef, 1. / (bcoef - 1.)) - 1.);

      i = 0;
      ftemp = L[2][1] / ccoef;
      if (which == 0)
	dx[1][1][i] = dx0 = ftemp;
      else if (which == 1)
	idx[1][1][i] = dx0 = ftemp;

      if (nonunigridx1 == 4) {
	sx[1] = dx0 * (bcoef - 1.0) / log(acoef);
      }
      if (myid <= 0) {
	fprintf(logfull_file,
		"Factor: %15.10g Nr: %15.10g ccoef: %15.10g dx0: %15.10g\n",
		acoef, bcoef, ccoef, dx0);
      }
      // just copy down to boundary zones
      // for(i=-1;i>=-NB[1];i--){
      // if(which==0) dx[1][1][i]=dx0;
      // else if(which==1) idx[1][1][i]=dx0;
      // }

      for (i = -NB[1]; i < N[1] + NB[1]; i++) {	// compute as
	// factor of dx0
	// so don't
	// accumulate
	// error
	ftemp = pow(acoef, (SFTYPE) (i) / (bcoef - 1.)) * dx0;
	if (which == 0)
	  dx[1][1][i] = ftemp;
	else if (which == 1)
	  idx[1][1][i] = ftemp;
      }
    } else if (nonunigridx1 == 5) {
      acoef = 10.0;		// Factor by which each larger decade
      // dx
      // is to last decade.
      bcoef = 10.0;		// definition of decade in this base in
      // log(useless)
      ftemp1 = log(L[1][1] - rgp) / log(acoef);	// xin
      ftemp2 = log(L[1][1] + L[2][1] - rgp) / log(acoef);	// xout

      ccoef = -log((L[1][1] - rgp) / (L[1][1] + L[2][1] - rgp)) / log(bcoef);	// Number 
										// 
      // 
      // of 
      // decades

      dellogx = (ftemp2 - ftemp1) / ((SFTYPE) (N[1]));	// dx in
      // log

      if (myid <= 0) {
	fprintf(logfull_file,
		"xin: %15.10g xout: %15.10g F: %15.10g D: %15.10g Nd: %15.10g dellogx: %15.10g\n",
		ftemp1, ftemp2, acoef, bcoef, ccoef, dellogx);
	fflush(logfull_file);
      }

      for (i = -NB[1]; i < N[1] + NB[1]; i++) {
	ftemp =
	    pow(acoef,
		ftemp1 + i * dellogx) * (pow(acoef, dellogx) - 1.0);
	if (which == 0)
	  dx[1][1][i] = ftemp;
	else if (which == 1)
	  idx[1][1][i] = ftemp;
      }
    } else {
      fprintf(fail_file, "No known X1 grid type specified\n");
      myexit(1);
    }


    // X2
    if (which == 1) {
      iL[1][2] = L[1][2];
      iL[2][2] = L[2][2];
    }

    if (numprocs > 1) {
    } else {
      i = 0;
      itemp = N[2];
      if (which == 0)
	sizes[2][i] = itemp;
      else if (which == 1)
	isizes[2][i] = itemp;
    }

    // now find place on grid(only needed for real grid)
    if (which == 0) {
      startj = 0;
      for (i = 0; i < myid; i++) {
	startj += sizes[2][i];
      }
      endj = 0;
      for (i = 0; i <= myid; i++) {
	endj += sizes[2][i];
      }
      endj -= 1;		// actual j-location from sizes
    } else if (which == 1) {
      startj = 0;
      for (i = 0; i < myid; i++) {
	startj += isizes[2][i];
      }
      endj = 0;
      for (i = 0; i <= myid; i++) {
	endj += isizes[2][i];
      }
      endj -= 1;		// actual j-location from sizes
    }
    // add up sizes for total size of grid
    if (which == 0)
      totalsize[2] = 0;
    else if (which == 1)
      itotalsize[2] = 0;
    for (i = 0; i < numprocs; i++) {
      if (which == 0)
	totalsize[2] += sizes[2][i];
      else if (which == 1)
	itotalsize[2] += isizes[2][i];
    }
    if (which == 0)
      realsize2 = totalsize[2];
    else if (which == 1)
      realsize2 = itotalsize[2];
    else {
      fprintf(fail_file, "realsize2: no which: %d here\n", which);
      myexit(1);
    }


    if (nonunigridx2 == 0) {
      xasofar = L[1][2];
      for (j = -NB[2]; j < realsize2 + NB[2]; j++) {

	ftemp = L[2][2] / ((SFTYPE) (realsize2));


	if (j == 0)
	  sx[2] = xasofar;
	if (which == 0)
	  dx[1][2][j] = ftemp;
	else if (which == 1)
	  idx[1][2][j] = ftemp;

	if ((j >= 0) && (j < realsize2)) {
	  xasofar += ftemp;
	}
      }
    } else if (nonunigridx2 == 1) {
      acoeff = (L[2][2] / ((SFTYPE) realsize2));
      acoefi = (L[2][2] / ((SFTYPE) (realsize2 + 1)));
      // must have a>b for nonnegative dx2
      acoef = (acoeff - acoefi) * .9 + acoefi;	// add some% off
      // inner acoef to
      // get smallest dx 
      // in center
      bcoef = L[2][2] - (SFTYPE) (realsize2) * acoef;

      xasofar = L[1][2];
      for (j = -NB[2]; j < realsize2 + NB[2]; j++) {
	ftemp =
	    acoef + bcoef * cos(2. * M_PI * j / ((SFTYPE) (N[2] - 1)));

	if (j == 0)
	  sx[2] = xasofar;
	if (which == 0)
	  dx[1][2][j] = ftemp;
	else if (which == 1)
	  idx[1][2][j] = ftemp;

	if ((j >= 0) && (j < realsize2)) {
	  xasofar += ftemp;
	}
      }
    } else if (nonunigridx2 == 2) {
      ccoef = 0.0;		// leave at 0.0
      // bcoef starts at >2.000000 and goes to
      // infinity(100.0~uniform) 
      bcoef = 2.1;		// free
      acoef =
	  (0.52 * L[2][2] -
	   0.5 * ccoef * realsize2) / (alnfact(bcoef +
					       0.5 * realsize2 -
					       1.) - alnfact(bcoef -
							     1.));

      // printf("acoef: %21.10g\n",acoef);
      xasofar = L[1][2];
      for (j = -NB[2]; j < realsize2 + NB[2]; j++) {
	if (j >= realsize2 / 2) {
	  ftemp =
	      ccoef + acoef * log((SFTYPE) j -
				  (SFTYPE) N[2] *
				  (SFTYPE) numprocs * 0.5 - 1.0 +
				  bcoef);
	} else if (j <= realsize2 / 2 - 1) {
	  ftemp =
	      ccoef +
	      acoef * log((SFTYPE) N[2] * (SFTYPE) numprocs *
			  0.5 - 1.0 - (SFTYPE) j - 1.0 + bcoef);
	}

	if (j == 0)
	  sx[2] = xasofar;
	if (which == 0)
	  dx[1][2][j] = ftemp;
	else if (which == 1)
	  idx[1][2][j] = ftemp;

	if ((j >= 0) && (j < realsize2)) {
	  xasofar += ftemp;
	}
      }
    } else if (nonunigridx2 == 3) {
      // N[2] needs to be factor of 2 right now

      // F
      acoef = 4.0;		// Factor by which each larger decade
      // dx
      // is to last decade. 
      dcoef = (SFTYPE) (realsize2) * 0.5;	// number of zones per
      // decade for half the
      // region (usually N/2 or
      // something less)
      truea = pow(acoef, 1.0 / (dcoef - 1.));	// true factor of
      // increase per
      // entire range
      // N
      bcoef = (SFTYPE) (realsize2);	// number of zones in Pi

      ccoef =
	  L[2][2] * 0.5 * (truea - 1.) / (pow(truea, bcoef * 0.5) -
					  1.0);
      xasofar = L[1][2];
      for (j = -NB[2]; j < realsize2 + NB[2]; j++) {

	// real dxa(j) grid value as if whole grid here 
	if (j <= realsize2 / 2 - 1) {	// if in -j half
	  // ftemp=ccoef*pow(acoef,(bcoef*0.5-1.0-(SFTYPE)j)/(bcoef*0.5-1.0));
	  ftemp = ccoef * pow(truea, bcoef * 0.5 - 1.0 - (SFTYPE) j);
	} else if (j >= realsize2 / 2) {
	  // ftemp=ccoef*pow(acoef,(j-bcoef*0.5)/(bcoef*0.5));
	  ftemp = ccoef * pow(truea, (SFTYPE) (j) - bcoef * 0.5);
	}
	// can use this below part in general for any grid!
	// Just
	// make sure above calcs are based upon full grid, not
	// each cpu

	if (j == 0)
	  sx[2] = xasofar;
	if (which == 0)
	  dx[1][2][j] = ftemp;
	else if (which == 1)
	  idx[1][2][j] = ftemp;

	if ((j >= 0) && (j < realsize2)) {
	  xasofar += ftemp;
	}
      }
    } else {
      fprintf(fail_file, "No known X2 grid type specified\n");
      myexit(1);
    }


    // X3
    if (which == 1) {
      iL[1][3] = L[1][3];
      iL[2][3] = L[2][3];
    }
    // add up sizes for total size of grid
    if (which == 0)
      totalsize[3] = N[3];
    else if (which == 1)
      itotalsize[3] = N[3];
    sx[3] = L[1][3];
    if (nonunigridx3 == 0) {
      for (k = -NB[3]; k < N[3] + NB[3]; k++) {
	ftemp = L[2][3] / N[3];
	if (which == 0)
	  dx[1][3][k] = ftemp;
	else if (which == 1)
	  idx[1][3][k] = ftemp;
      }
    } else {
      fprintf(fail_file, "No known X3 grid type specified\n");
      myexit(1);
    }

  }				// end if normal grid or interp wants
  // native grid type
  else {
    if (outtype == 1) {
      frac = 0.03;
      if (COORD == 3) {		// z vs x, y=0 plane
	startix[1] = 0.0;
	startix[2] = -frac * (L[1][1] + L[2][1]);
	startix[3] = L[1][3];
	iL[1][1] = startix[1];
	// iL[2][1]=2.0*frac*(L[1][1]+L[2][1]); // choose for
	// square output for square image size
	iL[2][1] = frac * (L[1][1] + L[2][1]);	// choose for
	// rect/square
	// aspect ratio
	// with rect image 
	// size
	iL[1][2] = startix[2];
	iL[2][2] = 2.0 * frac * (L[1][1] + L[2][1]);
	iL[1][3] = startix[3];
	iL[2][3] = 2.0 * M_PI;
	dx1 = iL[2][1] / N[1];
	dx2 = iL[2][2] / N[2];
      } else {
	startix[1] = L[1][1];
	startix[2] = L[1][2];
	startix[3] = L[1][3];
	iL[1][1] = startix[1];
	iL[2][1] = L[2][1];
	iL[1][2] = startix[2];
	iL[2][2] = L[2][2];
	iL[1][3] = startix[3];
	iL[2][3] = L[2][3];
	dx1 = iL[2][1] / N[1];
	dx2 = iL[2][2] / N[2];
      }
    } else {
      fprintf(fail_file, "no outtype: %d here\n", outtype);
      myexit(1);
    }				// if outtype==1
    for (i = -NB[1]; i < N[1] + NB[1]; i++) {
      idx[1][1][i] = dx1;
    }
    for (j = -NB[2]; j < N[2] + NB[2]; j++) {
      idx[1][2][j] = dx2;
    }
    k = 0;
    idx[1][3][k] = L[1][3] + L[2][3];
  }				// if other outtypes for interp grid

  // now assign b-grid, which merely depends on a-grid except for
  // dxb(-1) which is never used, but set anyways

  if (which == 0) {
    for (m = 1; m <= 3; m++) {
      for (i = -NB[m] + 1; i < N[m] + NB[m]; i++) {
	dx[2][m][i] = 0.5 * (dx[1][m][i] + dx[1][m][i - 1]);
      }
      if (m != 3) {
	dx[2][m][-NB[m]] = dx[2][m][-NB[m] + 1];	// for
	// definiteness, 
	// really
	// not
	// defined 
	// from
	// a-grid
      } else {
	dx[2][m][0] = dx[1][m][0];	// since x3-dir has no
	// boundary zones and has
	// only 1 grid zone for
	// now
      }
    }
  } else if (which == 1) {
    for (m = 1; m <= 3; m++) {
      for (i = -NB[m] + 1; i < N[m] + NB[m]; i++) {
	idx[2][m][i] = 0.5 * (idx[1][m][i] + idx[1][m][i - 1]);
      }
      if (m != 3) {
	idx[2][m][-NB[m]] = idx[2][m][-NB[m] + 1];	// for
	// definiteness, 
	// really
	// not
	// defined 
	// from
	// a-grid
      } else {
	idx[2][m][0] = idx[1][m][0];	// since x3-dir has no
	// boundary zones and has
	// only 1 grid zone for
	// now
      }
    }
  }

  if (which == 0)
    totalzones = totalsize[1] * totalsize[2] * totalsize[3];
  else if (which == 1)
    itotalzones = itotalsize[1] * itotalsize[2] * itotalsize[3];


  return (0);
}


#define RGMOVE (.001)		// how far to move beyond rgp if start
																// off
				// below rgp
#define STARTF (1.0/4.0)	// what factor or less 2nd condition
																// below
				// should be

/* Once dxs fixed, below applies to any grid */
/* BEGIN: Setup x1,x2,x3 positions from dx data */
/* I break this loop from the previous because x depends on diffs in dx 
 */
/* That is to say I do not have to worry about order problems */
/* Note this data would be dynamic if dynamic grid */
// 1/2/3-d valid
int init_x(int n1, int n2, int n3, int n1b, int n2b, int n3b,
	   SFTYPE * sx, int which, int outtype)
{
  int i, j, k, l, m;
  SFTYPE startbc[3 + 1] = { 0 };
  int N[3 + 1], NB[3 + 1];
  SFTYPE ftemp=0;


  N[1] = n1;
  N[2] = n2;
  N[3] = n3;
  NB[1] = n1b;
  NB[2] = n2b;
  NB[3] = n3b;

  // determine first BC's starting a-grid position
  startbc[1] = sx[1];
  if ((which == 0) || (outtype == 0)) {

    for (i = -1; i >= -NB[1]; i--) {
      if (which == 0)
	startbc[1] -= dx[1][1][i];
      else if (which == 1)
	startbc[1] -= idx[1][1][i];
    }
    if ((COORD == 3)
	&& ((reflectix1 == 0)
	    || ((reflectix1 == 1) && (L[1][1] > 1.0E-6)))) {
      if (startbc[1] <= rgp) {
	if (myid <= 0) {
	  fprintf(logfull_file,
		  "r coordinate assigned at boundary value to less than rg!\n");
	}
	if (nonunigridx1 == 4) {
	  fprintf(fail_file,
		  "This grid requires sx[1] to be set and not changed\n");
	  myexit(1);
	} else {
	  if (myid <= 0) {
	    fprintf(logfull_file,
		    "... moving from: %21.10g to ", startbc[1]);
	  }
	  // so force to be good
	  if (which == 0)
	    startbc[1] = dx[2][1][-NB[1]] / (STARTF * .9999) + rgp;
	  else if (which == 1)
	    startbc[1] = idx[2][1][-NB[1]] / (STARTF * .9999) + rgp;
	  // sx[1]+=(rgp-startbc[1]+RGMOVE);
	  // startbc[1]=sx[1];
	  // for(i=-1;i>=-NB[1];i--){
	  // if(which==0) startbc[1]-=dx[1][1][i];
	  // else if(which==1) startbc[1]-=idx[1][1][i];
	  // }
	  // now get comp grid start value
	  sx[1] = startbc[1];
	  for (i = -1; i >= -NB[1]; i--) {
	    if (which == 0)
	      sx[1] += dx[1][1][i];
	    else if (which == 1)
	      sx[1] += idx[1][1][i];
	  }
	  if (startbc[1] <= rgp) {
	    fprintf(fail_file,
		    "Failed to correct starting position(1)\n");
	    myexit(1);
	  } else {
	    if (myid <= 0) {
	      fprintf(logfull_file, " %21.10g\n", startbc[1]);
	    }
	  }
	  if (myid <= 0) {
	    fprintf(logfull_file,
		    "Can only have negative (r-rgp) with reflecting BC and rgp=0\n");
	  }
	}
      }
      // now check more conditions (zone i=0 used since bzones
      // may
      // be artificially low and)
      if (which == 0)
	ftemp = dx[2][1][0] / (startbc[1] - rgp);
      else if (which == 1)
	ftemp = idx[2][1][0] / (startbc[1] - rgp);
      else {
	fprintf(fail_file, "no which: %d here\n", which);
	myexit(1);
      }
      if (myid <= 0) {
	if (ftemp > STARTF) {
	  fprintf(logfull_file,
		  "Bad R=dr/(r-rgp) inner boundary ratio: %15.10g\n",
		  ftemp);
	  fprintf(logfull_file,
		  "This bad R has known to produce instabilities with gam=5/3\n");
	} else {
	  fprintf(logfull_file,
		  "Good R=dr[0]/(r[0]-rgp) inner boundary ratio: %15.10g\n",
		  ftemp);
	}
      }
    }
    // only if no CPU chop at different fixed radii
    L[1][1] = sx[1];
  }



  startbc[2] = sx[2];
  if ((which == 0) || (outtype == 0)) {
    for (j = -1; j >= -NB[2]; j--) {
      if (which == 0)
	startbc[2] -= dx[1][2][j];
      else if (which == 1)
	startbc[2] -= idx[1][2][j];
    }
  }
  startbc[3] = sx[3];

  if ((which == 0) || (outtype == 0)) {
    for (k = -1; k >= -NB[3]; k--) {
      if (which == 0)
	startbc[3] -= dx[1][3][k];
      else if (which == 1)
	startbc[3] -= idx[1][3][k];
    }
    for (l = 1; l <= 3; l++) {
      fprintf(log_file, "sx[%d]: %21.10g startbc[%d]: %21.10g\n", l,
	      sx[l], l, startbc[l]);
    }
  }

  if (which == 0) {
    // allow r<0 so volume factor is correct?
    for (m = 1; m <= 3; m++) {
      for (i = -NB[m]; i < N[m] + NB[m]; i++) {
	if (i == -NB[m]) {
	  x[1][m][i] = startbc[m];	/* a-grid */
	  x[2][m][i] = x[1][m][i] + dx[1][m][i] * 0.5;	/* b-grid */
	}
	/* Applying rule for a-grid: x(i)=x(i-1)+dx(i-1) */
	/* Applying rule for b-grid: x(i)=x(i-1)+dx(i) */
	else {
	  x[1][m][i] = x[1][m][i - 1] + dx[1][m][i - 1];
	  x[2][m][i] = x[2][m][i - 1] + dx[2][m][i];
	}
      }
    }				/* END: Setup x1,x2,x3 positions from
				   dx data */
  } else if (which == 1) {
    // allow r<0 so volume factor is correct?
    for (m = 1; m <= 3; m++) {
      for (i = -NB[m]; i < N[m] + NB[m]; i++) {
	if (i == -NB[m]) {
	  ix[1][m][i] = startbc[m];	/* a-grid */
	  ix[2][m][i] = ix[1][m][i] + idx[1][m][i] * 0.5;	/* b-grid 
								 */
	}
	/* Applying rule for a-grid: x(i)=x(i-1)+dx(i-1) */
	/* Applying rule for b-grid: x(i)=x(i-1)+dx(i) */
	else {
	  ix[1][m][i] = ix[1][m][i - 1] + idx[1][m][i - 1];
	  ix[2][m][i] = ix[2][m][i - 1] + idx[2][m][i];
	}
      }
    }				/* END: Setup x1,x2,x3 positions from
				   dx data */

  }

  return (0);
}



/* Currently assume g1=1 in rest of code, h=high, l=low */
#if(COORD==1)
/* x1=x, x2=y, x3=z */
#define g1(x1) 1
#define dg1(x1) 0		// always zero currently, not used
				// currently
#define g2(x1) 1
#define dg2(x1) 0		// par(g2i,dx1)
#define g31(x1) 1
#define dg31(x1) 0		// par(g31,dx1)
#define g32(x2) 1
#define dg32(x2) 0		// par(g32,dx2)

#define dvl1(x1h,x1l) (x1h-x1l)
#define dvl2(x2h,x2l) (x2h-x2l)
#define dvl3(x2h,x2l) (x2h-x2l)

#endif
#if(COORD==2)
/* x1=z, x2=r, x3=phi */
#define g1(x1) 1
#define dg1(x1) 0
#define g2(x1) 1
#define dg2(x1) 0
#define g31(x1) 1
#define dg31(x1) 0
#define g32(x2) fabs(x2)
#define dg32(x2) 1

#define dvl1(x1h,x1l) (x1h-x1l)
#define dvl2(x2h,x2l) (0.5*(x2h*x2h-x2l*x2l))
#define dvl3(x2h,x2l) (x2h-x2l)

#endif
#if(COORD==3)
/* x1=r, x2=theta, x3=phi */
#define g1(x1) 1
#define dg1(x1) 0
#define g2(x1) fabs(x1)		// fabs meant to make so metric isn't
																// <0?
#define dg2(x1) 1
#define g31(x1) fabs(x1)
#define dg31(x1) 1
#define g32(x2) fabs(sin(x2))
#define cotangent(x2) (1./tan(x2))
#define dg32(x2) (cos(x2))	// not sure, should it be periodic due
																// to
				// g32?

#define dvl1(x1h,x1l) (THIRD*(x1h*x1h*x1h-x1l*x1l*x1l))
#define dvl2(x2h,x2l) (-(cos(x2h)-cos(x2l)))
#define dvl3(x2h,x2l) (x2h-x2l)

#endif


// 1/2/3-d valid

// Note: As with dxs, only necessary to have g and dg and dvl defined
// on:
// For advection terms(for fluxes):
// for 1-dir(sweepx):
// g2a: i=0..N1 , j=0..N2-1
// g2b: i=-1..N1-1, j=0..N2-1
// g31a: i=0..N1 , j=0..N2-1
// g31b: i=-1..N1-1, j=0..N2-1
// g32a: not needed
// g32b: i=0..N1, j=0..N2-1
// dvl2a: i=-1..N1, j=0..N2-1
// dvl2b: i=0..N1, j=0..N2-1
// for 2-dir(sweepy):
// g2a: not needed
// g2b: i=0..N1-1, j=-1..N2-1
// g31a: i=0..N1-1 , j=0..N2
// g31b: i=0..N1-1, j=0..N2
// g32a: i=0..N1-1, j=0..N2
// g32b: i=0..N1-1, j=-1..N2-1
// rest only need to be on active grid


int init_diffs(void)
{

  int i, l, m;
  SFTYPE xm, xml=0, xmh=0;

  int N[3 + 1], NB[3 + 1];

  N[1] = N1;
  N[2] = N2;
  N[3] = N3;
  NB[1] = N1BND;
  NB[2] = N2BND;
  NB[3] = N3BND;

  /* BEGIN: assign values for matric scales and volume factors */
  /* I break this loop from the previous because metric may need diffs
     in x */
  /* That is to say I do not have to worry about order problems */
  /* Note this data would be dynamic if dynamic grid */

  // Not general!! --must modify below for different coords than
  // given
  // above
  for (l = 1; l <= NUMGRID; l++) {
    for (m = 1; m <= 2; m++) {	// no need for 3
      for (i = -NB[m]; i < N[m] + NB[m]; i++) {
	/* Setup metric */
	xm = x[l][m][i];

	if (m == 1) {
	  g[l][1][i] = g1(xm);
	  dg[l][1][i] = dg1(xm);
	  g[l][2][i] = g2(xm);
	  dg[l][2][i] = dg2(xm);
	  g[l][3][i] = g31(xm);
	  dg[l][3][i] = dg31(xm);
	} else if (m == 2) {
	  g[l][4][i] = g32(xm);
#if(COORD==3)
	  if (((i == 0) && (l == 1)) || ((i == N2) && (l == 1)))
	    cotan[l][i] = 0;	// just for printability,
	  // should never use this
	  // value!!
	  else
	    cotan[l][i] = cotangent(xm);
#else
	  cotan[l][i] = 1.0;
#endif

	  dg[l][4][i] = dg32(xm);
	}
      }
    }
  }
  /* END: assign values for matric scales and volume factors */

  /* Assign volume differences */
  for (l = 1; l <= NUMGRID; l++) {
    for (m = 1; m <= 3; m++) {
      for (i = -NB[m]; i < N[m] + NB[m]; i++) {
	if (l == 1) {
	  if (i == N[m] + NB[m] - 1) {
	    xml = x[l][m][i];
	    xmh = xml + dx[l][m][i];
	  } else {
	    xml = x[l][m][i];
	    xmh = x[l][m][i + 1];
	  }
	} else if (l == 2) {
	  if (i == -NB[m]) {
	    xmh = x[l][m][i];
	    xml = xmh - dx[l][m][i];
	  } else {
	    xml = x[l][m][i - 1];
	    xmh = x[l][m][i];
	  }
	}

	if (m == 1)
	  dvl[l][m][i] = fabs(dvl1(xmh, xml));	// fabs
	// for
	// boundary 
	// zones.
	// Shouldn't 
	// have
	// negative 
	// volumes
	else if (m == 2)
	  dvl[l][m][i] = fabs(dvl2(xmh, xml));
	else if (m == 3)
	  dvl[l][m][i] = fabs(dvl3(xmh, xml));
      }

    }
  }

  // correct for volume diff across theta=0,Pi boundary on b-grid
  // when
  // reflecting...r is already fine
  // below needed for x1-dir vy-adv at theta=0,theta=Pi
  // but not really used for v, just a place holder for advecting vy
  // in
  // x-dir.  Should really take care of there instead.

  // volume at theta edges are nonzero!
  if ((COORD == 3) && (reflectix2 == 1) && (L[1][2] < 1.0E-6)) {
    if (myid <= 0) {
      dvl[2][2][0] =
	  fabs(dvl2(x[2][2][0], 0.0)) + fabs(dvl2(0, x[2][2][-1]));
    }
  }
  if ((COORD == 3) && (reflectox2 == 1)
      && (L[1][2] + L[2][2] > M_PI - 1.0E-6)) {
    if (myid == numprocs - 1) {
      dvl[2][2][N2] =
	  fabs(dvl2(M_PI, x[2][2][N2 - 1])) +
	  fabs(dvl2(x[2][2][N2], M_PI));
    }
  }

  return (0);
}

void init_compsave(void)
{
  int i, j, k;

  // DS[type][direction]...
  LOOPF {
    // .5,.5 (used with centered variables)
    DS[1][1][k][j][i] = g[1][2][i] * g[1][3][i] * dvl[1][2][j];
    DS[1][2][k][j][i] = g[2][3][i] * dx[1][1][i] * g[1][4][j];
    DS[1][3][k][j][i] = 0;

    // 0,.5(eg. vx)
    DS[2][1][k][j][i] = g[2][2][i] * g[2][3][i] * dvl[1][2][j];
    DS[2][2][k][j][i] = g[1][3][i] * dx[2][1][i] * g[1][4][j];
    DS[2][3][k][j][i] = 0;

    // .5,0(eg. vy)
    DS[3][1][k][j][i] = g[1][2][i] * g[1][3][i] * dvl[2][2][j];
    DS[3][2][k][j][i] = g[2][3][i] * dx[1][1][i] * g[2][4][j];
    DS[3][3][k][j][i] = 0;

    // 0,0
    DS[4][1][k][j][i] = g[2][2][i] * g[2][3][i] * dvl[2][2][j];
    DS[4][2][k][j][i] = g[1][3][i] * dx[2][1][i] * g[2][4][j];
    DS[4][3][k][j][i] = 0;

    // OARCL[type][dir]...

    // .5,.5 (used with centered variables)
    // length along dir direction, used for the given variable
    // differences
    OARCL[1][1][k][j][i] = 1.0 / dx[2][1][i];
    OARCL[1][2][k][j][i] = 1.0 / (g[2][2][i] * dx[2][2][j]);
    OARCL[1][3][k][j][i] = 0;

    // 0,.5(eg. vx)
    OARCL[2][1][k][j][i] = 1.0 / dx[1][1][i];
    if ((COORD == 3) && (reflectix1 == 1) && (L[1][1] < 1.0E-6)
	&& (i == 0)) {
      OARCL[2][2][k][j][i] = 1.0;	// set to this to avoid
      // singularity
      OARCL[2][3][k][j][i] = 0;
    } else {
      OARCL[2][2][k][j][i] = 1.0 / (g[1][2][i] * dx[2][2][j]);
      OARCL[2][3][k][j][i] = 0;
    }

    // .5,0(eg. vy)
    OARCL[3][1][k][j][i] = 1.0 / dx[2][1][i];
    OARCL[3][2][k][j][i] = 1.0 / (g[2][2][i] * dx[1][2][j]);
    OARCL[3][3][k][j][i] = 0;

    // 0,0(eg. sigma_{r,theta})
    OARCL[4][1][k][j][i] = 1.0 / dx[1][1][i];
    if ((COORD == 3) && (reflectix1 == 1) && (L[1][1] < 1.0E-6)
	&& (i == 0)) {
      OARCL[4][2][k][j][i] = 1.0;	// to avoid sing
      OARCL[4][3][k][j][i] = 0;
    } else {
      OARCL[4][2][k][j][i] = 1.0 / (g[1][2][i] * dx[1][2][j]);
      OARCL[4][3][k][j][i] = 0;
    }


    // OVOL[type]...

    // .5,.5(eg. rho)
    OVOL[1][k][j][i] = 1.0 / (dvl[1][1][i] * dvl[1][2][j]);
    // .5,0(eg. vy)
    OVOL[3][k][j][i] = 1.0 / (dvl[1][1][i] * dvl[2][2][j]);
    // 0,.5(eg. vx)
    OVOL[2][k][j][i] = 1.0 / (dvl[2][1][i] * dvl[1][2][j]);

    // 0,0
    OVOL[4][k][j][i] = 1.0 / (dvl[2][1][i] * dvl[2][2][j]);
  }


}



// 1/2/3-d valid
int init_data(void)
{

  int i, j, k, l, m;

  /* BEGIN: zero out(assign) data */

  /* I break this loop from the previous because sca/vev/metric may
     need diffs in x */
  /* That is to say I do not have to worry about order problems */
  /* Note this data would be dynamic if dynamic grid */

  if (analoutput != 0) {

    analsolve(0);
    tdep_compute();		// compute time dep stuff


    for (l = 1; l <= NUMSCA; l++) {
      LOOPF {
	s[l][k][j][i] = sanal[l][k][j][i];
      }
    }

    for (m = 1; m <= NUMVEC - 1; m++) {
      for (l = 1; l <= 3; l++) {
	LOOPF {
	  v[m][l][k][j][i] = vanal[m][l][k][j][i];
	}
      }
    }
  } else {
    // POTENTIAL
    LOOPF {
      s[3][k][j][i] = -GRAVC * MASSBH / (x[2][1][i] - rgp);
      // s[3][k][j][i] = 0.0 ;
    }


    // MASS DENSITY 
    LOOPF {
      s[1][k][j][i] = 1.0;
    }

    // INTERNAL ENERGY DENSITY
    LOOPF {
      if (press == 1) {
	if (wgam) {
	  cs = pow(s[1][k][j][i], 0.5 * (gam - 1.));	// usually 
	  // 
	  // not
	  // needed
	  // here
	  s[2][k][j][i] = pow(s[1][k][j][i], gam) / gam / (gam - 1.);
	} else {
	  s[2][k][j][i] = cs * cs * s[1][k][j][i];
	}
      } else {
	s[2][k][j][i] = 0;
      }
    }

    // VECTORS(v and B)
    for (m = 1; m <= 3; m++) {
      LOOPF {
	v[1][m][k][j][i] = 0.0;
      }
    }

    for (m = 1; m <= 3; m++) {
      LOOPF {
	v[2][m][k][j][i] = 0.0;
      }
    }
  }
  /* END: zero out(assign) data */
  return (0);
}
