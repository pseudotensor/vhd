#include "global.h"
#include "defs.h"



// General analytic dump for given analytic solution
void analsolve(int gopp)
{
  static int timeless = 0;	// =0 not computed stationary solution
  // yet 1=have computed, so skip

  if (analoutput == 5) {
    if ((timeless == 0) || (gopp)) {
      timeless = 1;
      tori1sol(gopp);
    }
  } else if (analoutput == 8) {
    if ((timeless == 0) || (gopp)) {
      timeless = 1;
      injectsol(gopp);
    }
  }
}





#define ATMONLYT 0		// 0: normal torus 1: just atmosphere
#define ADATM 0			// 1: adiabatic atm 0: fraction pot atm

// Fill analytic solution with solution to tori solution, mostly used
// for initial conditions and comparisons as with bondi

FILE *tori1sol(int gopp)
{
  char check[50];
  int i, j, k;
  SFTYPE Cconst;
  SFTYPE C2const, Beta, RBeta;
  char temps[100];
  char filename[100];

  SFTYPE Rpot0 = 0;
  SFTYPE theta0;
  SFTYPE Rin;
  SFTYPE Dconst;

  SFTYPE phi_eff;
  static FILE *analyticout;
  static int firstsolve = 1;
  SFTYPE totalmass[2];		// 0: tori 1: atmosphere
  SFTYPE totalmass_full[2];	// 0: tori 1: atmosphere


  if (1 || (gopp == 0)) {

    if (!wgam) {
      fprintf(fail_file,
	      "tori1sol: No solution setup for gam==1 yet\n");
      myexit(1);
    }
    // want to write some interesting data on solution
    if (firstsolve == 1) {
      if (myid <= 0) {
	strcpy(temps, DATADIR);

	sprintf(filename, "%s0_analdata%s", temps, DATEXT);
	if ((analyticout = fopen(filename, "wt")) == NULL) {
	  fprintf(fail_file, "Cannot open %s\n", temps);
	  myexit(1);
	}
      }				// end if cpu write
    }
    // setup tori conditions

#define EQATM 0
#define DELTAR .5
#define QCONST 2.0
    theta0 = M_PI / 2.0;

#define UNITS 2
    // 1: original torus units(R0=1,Omega0=1,M=R0^3*rho0)
    // 2: GM=rg/2=Omega0=1

    if (UNITS == 1) {
      R0 = 1.0;			// must be 1.0
      Omega0 = 1.0;		// must be 1.0
      rho0 = 1.0;		// must be 1.0
      GM = R0 * R0 * R0 * Omega0 * Omega0 *
	  pow((1. / sin(theta0) - rgp / R0), 2.0);
    } else if (UNITS == 2) {
      if (fabs(rgp - 2.0) > ERR) {
	fprintf(fail_file,
		"These units require rgp to be set to"
		"2.0 in init.c!\n");
	myexit(1);
      }
      rho0 = 1.0;
      GM = 1.0;
      R0 = (.825) * (L[1][1] + L[2][1]);
      // below forced
      Omega0 = sqrt(GM / (R0 * R0 * R0 * pow(1.0 - rgp / R0, 2.0)));
    } else {
      // whatever

    }

    // rest should be generally true for any units

#define PHI_EFF(R) (R0*R0*Omega0*Omega0*pow((R/R0),(2.-2.*QCONST))\
                   /(2.*QCONST-2.))

#define PHI_G(r) (-GM/((r)-rgp))
#define V_PHI(R) (R0*Omega0*pow((R/R0),1.-QCONST))
    // below line on theta=Pi/2 axis only and rgp=0
#define R_P0 (R0*pow(GM*(2.*QCONST-2.)/(Omega0*Omega0*R0*R0*R0)\
                     ,1./(3.-2.*QCONST)))

    Rin = .75 * R0;		// SPB d=1.125 for
    // rgp=0,theta0=Pi/2,q=2

    Cconst = PHI_EFF(Rin) + PHI_G(Rin / sin(M_PI / 2.0));
    Dconst =
	(Cconst - PHI_EFF(R0) - PHI_G(R0 / sin(theta0))) * (1. / gam);

    // atmosphere params
    C2const = -1.0;		// for equilibrium atmosphere
    Beta = 1.0;			// factor on tori KAPPA for any
    // atmosphere

    // used for non adiabatic atm
    IEFRACT = 0.2;
    VZFRACT = .95;

    DENSITYFLOOR = 1.E-5;




    if (myid <= 0) {
      if (!wpw) {
	Rpot0 = R_P0;
	fprintf(analyticout, "GravPotential==0 @ R=%15.10g\n", Rpot0);
	if (Rin < Rpot0)
	  fprintf(analyticout,
		  "Warning, your set inner torus radius is smaller"
		  "than the 0 of the grav pot, so high dumping will"
		  "occur at large radius");
      } else {
	fprintf(analyticout, "No check to make sure Rin>Rpot0\n");
      }
    }
    if (Rin > R0) {
      fprintf(fail_file, "Invalid to have Rin>R0\n");
      myexit(1);
    }
    if (rgp > R0) {
      fprintf(fail_file, "Invalid to have rgp>R0\n");
      myexit(1);
    }

    totalmass[0] = 0.0;
    totalmass[1] = 0.0;

    LOOPF {
      sanal[3][k][j][i] = PHI_G(x[2][1][i]);

      phi_eff = PHI_EFF(x[2][1][i] * sin(x[2][2][j]));

      if (ATMONLYT == 0) {
	sanal[1][k][j][i] =
	    rho0 * pow((Cconst - sanal[3][k][j][i] - phi_eff) /
		       (gam * Dconst), 1. / (gam - 1.));
      } else {
	sanal[1][k][j][i] = DENSITYFLOOR * .1;	// .1 so forces 
	// 
	// check
      }
      sprintf(check, "%3f\n", sanal[1][k][j][i]);
      // atmosphere
      if ((check[0] == 'n') || (sanal[1][k][j][i] < DENSITYFLOOR)) {
	RBeta = Beta;

#if(EQATM==1)
	// implement equilibrium atmosphere
	sanal[1][k][j][i] =
	    pow((C2const -
		 sanal[3][k][j][i]) / (Beta * gam * Dconst),
		1. / (gam - 1.));

	sprintf(check, "%3f\n", sanal[1][k][j][i]);
	if ((check[0] == 'n') || sanal[1][k][j][i] < DENSITYFLOOR)
	  sanal[1][k][j][i] = DENSITYFLOOR;
#else
	sanal[1][k][j][i] = DENSITYFLOOR;
#endif
	vanal[1][1][k][j][i] = 0.0;
	vanal[1][2][k][j][i] = 0.0;
	vanal[1][3][k][j][i] = 0.0;

	if ((i >= 0) && (i < N1) && (j >= 0) && (j < N2)
	    && (k >= 0)
	    && (k < N3)) {	// only add that which is on
	  // real grid
	  totalmass[1] +=
	      sanal[1][k][j][i] * dvl[1][1][i] * dvl[1][2][j];
	}
	// SPB
	// sanal[2][k][j][i] =
	// RBeta*sanal[1][k][j][i]/((gam-1.)*x[2][1][i]); //en

	if (ADATM) {
	  if (wgam) {
	    sanal[2][k][j][i] = RBeta * Dconst * pow(sanal[1][k][j][i], gam);	// en
	  } else
	    sanal[2][k][j][i] = RBeta * alpha * cs * cs * sanal[1][k][j][i];	// en
	} else {
	  sanal[2][k][j][i] = IEFRACT * fabs(sanal[1][k][j][i] * sanal[3][k][j][i]);	// en
	}

      } else {			// torus
	RBeta = 1.0;
	vanal[1][3][k][j][i] = V_PHI(x[2][1][i] * sin(x[2][2][j]));

	vanal[1][1][k][j][i] = 0.0;
	vanal[1][2][k][j][i] = 0.0;

	if ((i >= 0) && (i < N1) && (j >= 0) && (j < N2)
	    && (k >= 0)
	    && (k < N3)) {	// only add that which is on
	  // real grid
	  totalmass[0] +=
	      sanal[1][k][j][i] * dvl[1][1][i] * dvl[1][2][j];
	}
	if (wgam) {
	  sanal[2][k][j][i] = RBeta * Dconst * pow(sanal[1][k][j][i], gam);	// en
	} else
	  sanal[2][k][j][i] = RBeta * alpha * cs * cs * sanal[1][k][j][i];	// en
      }
    }


    if (ADATM) {
      IEFLOOR = Dconst * pow(DENSITYFLOOR, gam);
    } else {
      IEFLOOR = IEFRACT * fabs(DENSITYFLOOR * sanal[3][0][0][N1]);
    }
    // output some interesting data
    // and setup initial problem stuff
    if (firstsolve == 1) {
      firstsolve = 0;

      if (numprocs > 1) {
      } else {
	totalmass_full[0] = totalmass[0];
	totalmass_full[1] = totalmass[1];
      }
      if (myid <= 0) {
	// now output some interesting analytic data

	fprintf(analyticout, "Rpot0: %15.10g Rin: %15.10g\n",
		Rpot0, Rin);
	fprintf(analyticout,
		"q: %15.10g Cconst: %15.10g Dconst: %15.10g\n",
		QCONST, Cconst, Dconst);
	fprintf(analyticout,
		"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",
		DENSITYFLOOR, IEFLOOR);
	if (1.0 * .001 < DENSITYFLOOR) {
	  fprintf(analyticout,
		  "Warning, density floor is higher than"
		  ".001 times of density max!\n");
	}

	fprintf(analyticout,
		"total tori mass: %15.10g total atm mass: %15.10g\n",
		totalmass_full[0], totalmass_full[1]);
	if (gopp != 2) {
	  fclose(analyticout);
	}
      }				// end if cpu write
    }				// end if firstsolve

  }				// end if gopp==0

  // ////////////////// VISUALIZATION settings

  // use image() to tell you how to set this(TOTALMM==1, and make
  // sure
  // DYNAMICMM==0)
  // for dvr=.1 Rin=.85


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



  // ////////////// ZOOM 
  i = 1;			// zoom view (vectors)
  j = 0;			// 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0] = 0.0;
  mmv[i][j][1][0][1] = 10.0;
  // vx1
  mmv[i][j][1][1][0] = -3.0;
  mmv[i][j][1][1][1] = .1;
  // vx2
  mmv[i][j][1][2][0] = -1.5;
  mmv[i][j][1][2][1] = 1.5;
  // vx3
  mmv[i][j][1][3][0] = -1;
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


  if (gopp == 2) {
    return (analyticout);
  } else {
    return (NULL);
  }
}



// inject solution(starting with just atmosphere)

#define DENSITYINJFLOOR (1.0E-5)
#define ENINJFLOOR (1.0E-8)

#define ADATM 0			// 1: adiabatic atm 0: fraction pot atm

#define TORIINIT 0		// 1: use tori as initial data 0: Don't

#define PATTERN 2		// 1: gaussian 2: tori

// relevant parameters for gaussian size, and for both gaussian/tori
// loop restrictions
// 
// if this is too far out, most of mass will outflow through outer
// boundary, making understanding of what fraction of injection the
// accretion is difficult
// x1
#define DIX1 ((0.8)*(L[1][1]+L[2][1]))
#define DOX1 ((0.85)*(L[1][1]+L[2][1]))
// x2
#define DIX2 ((3.0/8.0)*(L[1][2]+L[2][2]))
#define DOX2 ((5.0/8.0)*(L[1][2]+L[2][2]))

void injectsol(int gopp)
{
  int i, j, k;
  SFTYPE Cconst;
  SFTYPE Beta, RBeta;
  SFTYPE ftemp;
  char temps[100];
  char filename[100];
  SFTYPE ftemp1, ftemp2, ftempx1, ftempx2;
  SFTYPE R0;
  SFTYPE Dconst;
  SFTYPE MASSDOT, IEDOT;

  static FILE *analyticout;
  static int firstsolve = 1;
  SFTYPE totalmass[2];		// 0: tori 1: atmosphere
  SFTYPE totalmass_full[2];	// 0: tori 1: atmosphere
  SFTYPE massdot, massdottot, massdottot_full = 0, massrealdottot,
      massrealdottot_full = 0, massreal2dottot, massreal2dottot_full =
      0;
  SFTYPE iedottot, iedottot_full = 0, ierealdottot, ierealdottot_full =
      0, iereal2dottot, iereal2dottot_full = 0;
  int stoptag;

  if (1 || (gopp == 0)) {	// need computed values!

    if (!wgam) {
      fprintf(fail_file, "No solution setup for gam==1 yet\n");
      myexit(1);
    }

    if ((PATTERN == 2) || (TORIINIT == 1)) {
      analyticout = tori1sol(2);	// submit for solution and file 
					// 
      // 
      // setup
    } else {			// want to write some interesting data
      // on
      // solution
      if (firstsolve == 1) {
	if (myid <= 0) {
	  strcpy(temps, DATADIR);

	  sprintf(filename, "%s0_analdata%s", temps, DATEXT);
	  if ((analyticout = fopen(filename, "wt")) == NULL) {
	    fprintf(fail_file, "Cannot open %s\n", temps);
	    myexit(1);
	  }
	}			// end if cpu write
      }
    }
    // setup tori conditions

    MASSDOT = 1.0;
    IEDOT = 1.0;		// not really used currently
#define QCONST 2.0
#define MASSDOTTOTAL MASSDOT	// total amount of mass per unit time
    // injected
#define IEDOTTOTAL IEDOT	// total amount of internal energy per
    // unit
    // time injected
    GM = 1;
    // rg=2.0; // make sure set in init.c for these units
    IEFRACT = 0.2;		// normal IGU thing
    VZFRACT = .95;

    // in case want to use SPB visc with injection
    rho0 = 1;
    R0 = DIX1;
    // compare
    // forced 
    Omega0 = sqrt(GM / (R0 * R0 * R0 * pow(1.0 - rgp / R0, 2.0)));

#define PHI_G(r) (-GM/((r)-rgp))

    Beta = 1.0;			// factor on tori KAPPA for any
    // atmosphere
    RBeta = Beta;
    Dconst = 1.0;
    Cconst = 1.0;
    DENSITYFLOOR = 1.E-10;
    // DENSITYFLOOR=1.E-5;
    // IEFLOOR=Dconst*pow(DENSITYFLOOR,gam);

    tagii = 0;
    tagfi = N1;
    tagij = 0;
    tagfj = N2;
    tagik = 0;
    tagfk = N3;

    // determine region of loop
    for (i = -N1BND; i < N1 + N1BND; i++) {
      if (x[2][1][i] > DIX1) {
	tagii = i;
	break;
      }
    }
    if (i == N1 + N1BND) {	// then never found it, so not in this
      // domain
      tagii = i;
      tagfi = i;
      stoptag = 1;
    } else {
      for (i = N1 + N1BND - 1; i >= -N1BND; i--) {
	if (x[2][1][i] < DOX1) {
	  tagfi = i;
	  break;
	}
      }
      if (i == -N1BND - 1) {	// then never found it, so not
	// in this domain
	tagii = i + 1;
	tagfi = i + 1;
	stoptag = 1;
      } else {
	tagfi++;		// since loop says i<tagfi, not <=
      }
    }

    stoptag = 0;
    // now x2-dir
    for (j = -N2BND; j < N2 + N2BND; j++) {
      if (x[2][2][j] > DIX2) {
	tagij = j;
	break;
      }
    }
    if (j == N2 + N2BND) {
      tagij = j;
      tagfj = j;
      stoptag = 1;
    } else {
      for (j = N2 + N2BND - 1; j >= -N2BND; j--) {
	if (x[2][2][j] < DOX2) {
	  tagfj = j;
	  break;
	}
      }
      if (j == -N2BND - 1) {
	tagij = j + 1;
	tagfj = j + 1;
	stoptag = 1;
      } else {
	tagfj++;
      }
    }
    // correct for unresolved regions, moves out 1 or 2 zones(1 if
    // on
    // boundary already, 2 if not)
    if (tagii == tagfi) {
      if (tagii > 0) {
	tagii--;
      }
      if (tagfi < N1) {
	tagfi++;
      }
    }
    if (tagij == tagfj) {
      if (tagij > 0) {
	tagij--;
      }
      if (tagfj < N2) {
	tagfj++;
      }
    }

    t2gii = tagii;
    t2gfi = tagfi;
    t2gij = tagij;
    t2gfj = tagfj;
    t2gik = tagik = 0;
    t2gfk = tagfk = N3;

    // now correct normal tag
    if (tagii < 0)
      tagii = 0;
    if (tagii > N1)
      tagii = N1;
    if (tagfi < 0)
      tagfi = 0;
    if (tagfi > N1)
      tagfi = N1;
    if (tagij < 0)
      tagij = 0;
    if (tagij > N2)
      tagij = N2;
    if (tagfj < 0)
      tagfj = 0;
    if (tagfj > N2)
      tagfj = N2;

    // spread for rhoi,rhof in injection velocity part(interp needs
    // rhoi/rhof there)(both since v needs -1,-1 normally, and
    // since v
    // really changed at (+1,+1) need rhoi/rhof there outside it
    // for interp
    if (t2gii != t2gfi) {
      t3gii = t2gii - 1;
      t3gfi = t2gfi + 1;
    }
    if (t2gij != t2gfj) {
      t3gij = t2gij - 1;
      t3gfj = t2gfj + 1;
    }
    t3gik = t2gik;
    t3gfk = t2gfk;

    // spread (-1,-1) & (1,1) correction
    if (t3gii < -N1BND)
      t3gii = -N1BND;
    if (t3gij < -N2BND)
      t3gij = -N2BND;
    if (t3gfi > N1 + N1BND)
      t3gfi = N1 + N1BND;
    if (t3gfj > N2 + N2BND)
      t3gfj = N2 + N2BND;

    // spread for v since v really changed by next outer layer of
    // rho
    // interpolated between injection and no injection
    if (t2gii != t2gfi) {
      t4gii = t2gii;
      t4gfi = t2gfi + 1;
    }
    if (t2gij != t2gfj) {
      t4gij = t2gij;
      t4gfj = t2gfj + 1;
    }
    t4gik = t2gik;
    t4gfk = t2gfk;

    // spread (1,1) correction
    if (t4gfi > N1 + N1BND)
      t4gfi = N1 + N1BND;
    if (t4gfj > N2 + N2BND)
      t4gfj = N2 + N2BND;

    // NOW set the injection rate

    LOOPF {			// all zones
      rhoinject[k][j][i] = 0;
      eninject[k][j][i] = 0;
    }

    // for
    // multiple cpus, injection still only really occurs on grid
    // includes boundary zones so don't have to transfer rhoinject
    // effects
    // to other cpuss
    LOOPFINJ {
      if (PATTERN == 1) {
	// half widths
	ftemp1 = 0.5 * (DOX1 - DIX1);
	ftemp2 = 0.5 * (DOX2 - DIX2);
	// position minus position of central peak
	ftempx1 = x[2][1][i] - 0.5 * (DOX1 + DIX1);
	ftempx2 = x[2][2][j] - 0.5 * (DOX2 + DIX2);
	// distribution
	ftemp =
	    -ftempx1 * ftempx1 / (ftemp1 * ftemp1) -
	    ftempx2 * ftempx2 / (ftemp2 * ftemp2);
	massdot = 1.0 * exp(ftemp);
	// unnormalized
	rhoinject[k][j][i] = massdot * OVOL[1][k][j][i];
	eninject[k][j][i] =
	    IEFRACT * fabs(rhoinject[k][j][i] * s[3][k][j][i]);
      } else if (PATTERN == 2) {
	rhoinject[k][j][i] = sanal[1][k][j][i];	// from
	// tori1sol
	// eninject[k][j][i]=sanal[2][k][j][i]; // not really
	// well
	// defined since unnormalized
	eninject[k][j][i] = IEFRACT * fabs(rhoinject[k][j][i] * s[3][k][j][i]);	// to 
										// 
	// 
	// 
	// 
	// be 
	// similar
      }
      if (rhoinject[k][j][i] < DENSITYINJFLOOR)
	rhoinject[k][j][i] = DENSITYINJFLOOR;
      if (eninject[k][j][i] < ENINJFLOOR)
	eninject[k][j][i] = ENINJFLOOR;
    }
    // ?'s:
    // normalize eninject? Is en+en really tori en? No.

    // now only over normal grid
    massdottot = 0;
    iedottot = 0;
    // now only over real physical injection region
    LOOPINJ {			// if this were LOOP-> can't make this
      // LOOPINT because this is really determining 
      // the physical amount of mass injected,
      // regardless of accounting region
      // normalizer 
      massdottot += rhoinject[k][j][i] * dvl[1][1][i] * dvl[1][2][j];

      iedottot += eninject[k][j][i] * dvl[1][1][i] * dvl[1][2][j];
    }


    // unlike in most other places, the full value is needed by all 
    // 
    // CPUs
    if (numprocs > 1) {
    } else {
      massdottot_full = massdottot;
      iedottot_full = iedottot;
    }

    LOOPFINJ {			// LOOPF would work here too
      rhoinject[k][j][i] = rhoinject[k][j][i] * (MASSDOTTOTAL / massdottot_full);	// renormalize
      eninject[k][j][i] = eninject[k][j][i] * (MASSDOTTOTAL / massdottot_full);	// normalized 
										// 
      // 
      // 

      // same 
      // since 
      // eninject\propto 
      // rhoinject(same 
      // as 
      // computing 
      // eninject 
      // from 
      // above 
      // rhoinject 
      // instead)
      // *(IEDOTTOTAL/iedottot_full); // renormalize

      // fprintf(fail_file,"%d %d %d %15.10g
      // %15.10g\n",k,j,i,rhoinject[k][j][i],OVOL[1][k][j][i]);
    }

    // NOW FOR SOME CAUTIOUS COUNTING
    // 
    // now count up real mass per unit time to be injected
    massrealdottot = 0;
    ierealdottot = 0;
    // now only over real physical injection region
    LOOPINJ {			// if this were LOOP-> can't make this
      // LOOPINT because this is really determining 
      // the physical amount of mass injected,
      // regardless of accounting region
      massrealdottot += rhoinject[k][j][i]
	  * dvl[1][1][i] * dvl[1][2][j];
      ierealdottot += eninject[k][j][i] * dvl[1][1][i] * dvl[1][2][j];
    }


    // unlike in most other places, the full value is needed by all 
    // 
    // CPUs
    if (numprocs > 1) {
    } else {
      massrealdottot_full = massrealdottot;
      ierealdottot_full = ierealdottot;
    }

    // now count up real mass per unit time to be injected on full
    // active
    // grid
    massreal2dottot = 0;
    iereal2dottot = 0;
    // now only over real physical injection region
    LOOP {			// if this were LOOP-> can't make this
      // LOOPINT because this is really determining 
      // the physical amount of mass injected,
      // regardless of accounting region
      massreal2dottot += rhoinject[k][j][i]
	  * dvl[1][1][i] * dvl[1][2][j];
      iereal2dottot += eninject[k][j][i] * dvl[1][1][i] * dvl[1][2][j];
    }


    // unlike in most other places, the full value is needed by all 
    // 
    // CPUs
    if (numprocs > 1) {
    } else {
      massreal2dottot_full = massreal2dottot;
      iereal2dottot_full = iereal2dottot;
    }


    // setup initial atmosphere
    LOOPF {
      sanal[3][k][j][i] = PHI_G(x[2][1][i]);
    }

    if (ADATM) {
      IEFLOOR = Dconst * pow(DENSITYFLOOR, gam);
    } else {
      IEFLOOR = IEFRACT * fabs(DENSITYFLOOR * sanal[3][0][0][N1]);
    }

    totalmass[0] = 0.0;
    totalmass[1] = 0.0;
    LOOPF {
      if (TORIINIT == 0) {
	sanal[1][k][j][i] = DENSITYFLOOR;
	vanal[1][1][k][j][i] = 0.0;
	vanal[1][2][k][j][i] = 0.0;
	vanal[1][3][k][j][i] = 0.0;
      }
      if ((i >= 0) && (i < N1) && (j >= 0) && (j < N2) && (k >= 0)
	  && (k < N3)) {	// only add that which is on real grid
	totalmass[1] += sanal[1][k][j][i] * dvl[1][1][i] * dvl[1][2][j];
      }
      if (TORIINIT == 0) {
	if (ADATM) {
	  if (wgam) {
	    sanal[2][k][j][i] = RBeta * Dconst * pow(sanal[1][k][j][i], gam);	// en
	  } else
	    sanal[2][k][j][i] = RBeta * alpha * cs * cs * sanal[1][k][j][i];	// en
	} else {
	  sanal[2][k][j][i] = IEFRACT * fabs(sanal[1][k][j][i] * sanal[3][k][j][i]);	// en
	}
      }
    }

    // output some interesting data
    // and setup initial problem stuff
    if (firstsolve == 1) {
      firstsolve = 0;

      fprintf(log_file, "%7s %7s %7s %7s %7s %7s\n", "tagix1",
	      "tagox1", "tagix2", "tagox2", "tagix3", "tagox3");
      fprintf(log_file, "%7d %7d %7d %7d %7d %7d\n", tagii, tagfi,
	      tagij, tagfj, tagik, tagfk);
      fprintf(log_file, "%7s %7s %7s %7s %7s %7s\n", "t2gix1",
	      "t2gox1", "t2gix2", "t2gox2", "t2gix3", "t2gox3");
      fprintf(log_file, "%7d %7d %7d %7d %7d %7d\n", t2gii, t2gfi,
	      t2gij, t2gfj, t2gik, t2gfk);
      fprintf(log_file, "%7s %7s %7s %7s %7s %7s\n", "t3gix1",
	      "t3gox1", "t3gix2", "t3gox2", "t3gix3", "t3gox3");
      fprintf(log_file, "%7d %7d %7d %7d %7d %7d\n", t3gii, t3gfi,
	      t3gij, t3gfj, t3gik, t3gfk);
      fprintf(log_file, "%7s %7s %7s %7s %7s %7s\n", "t4gix1",
	      "t4gox1", "t4gix2", "t4gox2", "t4gix3", "t4gox3");
      fprintf(log_file, "%7d %7d %7d %7d %7d %7d\n", t4gii, t4gfi,
	      t4gij, t4gfj, t4gik, t4gfk);
      LOOPINJ {
	fprintf(log_file, "%d %d %d %21.15g\n", k, j, i,
		rhoinject[k][j][i]);
      }
      fflush(log_file);

      if (numprocs > 1) {
      } else {
	totalmass_full[0] = totalmass[0];
	totalmass_full[1] = totalmass[1];
      }
      if (myid <= 0) {
	// now output some interesting global analytic data

	fprintf(analyticout,
		"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",
		DENSITYFLOOR, IEFLOOR);
	fprintf(analyticout, "IEFRACT: %15.10g VZFRACT: %15.10g\n",
		IEFRACT, VZFRACT);
	fprintf(analyticout,
		"normalizer: %15.10g massdottotal: %15.10g\n",
		massdottot_full, MASSDOTTOTAL);
	fprintf(analyticout,
		"massdottotal: %21.15g massdottotal2: %21.15g\n",
		massrealdottot_full, massreal2dottot_full);
	fprintf(analyticout,
		"normalizer: %15.10g iedottotal: %15.10g\n",
		iedottot_full, IEDOTTOTAL);
	fprintf(analyticout,
		"iedottotal: %21.15g iedottotal2: %21.15g\n",
		ierealdottot_full, iereal2dottot_full);
	if (1.0 * .001 < DENSITYFLOOR) {
	  fprintf(analyticout,
		  "Warning, density floor is higher"
		  "than .001 times of density max!\n");
	}

	fprintf(analyticout,
		"total tori mass: %15.10g total atm mass: %15.10g\n",
		totalmass_full[0], totalmass_full[1]);
	fclose(analyticout);
      }				// end if cpu write
    }				// end if firstsolve

  }				// end if gopp==0


  // ////////////////// VISUALIZATION settings
  // 

  // scalars

  j = 0;			// normal output
  i = 0;			// outtype

  mms[i][j][1][0] = DENSITYFLOOR;
  mms[i][j][1][1] = .05;

  mms[i][j][2][0] = IEFLOOR;
  mms[i][j][2][1] = .006;

  mms[i][j][3][0] = sanal[3][0][0][0];
  mms[i][j][3][1] = sanal[3][0][0][N1];

  j = 1;			// second type of comp
  i = 0;			// view large


  mms[i][j][1][0] = DENSITYFLOOR;
  mms[i][j][1][1] = .05;

  mms[i][j][2][0] = IEFLOOR;
  mms[i][j][2][1] = .006;

  mms[i][j][3][0] = sanal[3][0][0][0];
  mms[i][j][3][1] = sanal[3][0][0][N1];


  j = 0;			// normal comp
  i = 1;			// view zoom


  mms[i][j][1][0] = DENSITYFLOOR;
  mms[i][j][1][1] = .05;

  mms[i][j][2][0] = IEFLOOR;
  mms[i][j][2][1] = .006;

  mms[i][j][3][0] = sanal[3][0][0][0];
  mms[i][j][3][1] = sanal[3][0][0][N1];

  j = 1;			// 2nd comp
  i = 1;			// view zoom


  mms[i][j][1][0] = DENSITYFLOOR;
  mms[i][j][1][1] = .05;

  mms[i][j][2][0] = IEFLOOR;
  mms[i][j][2][1] = .006;

  mms[i][j][3][0] = sanal[3][0][0][0];
  mms[i][j][3][1] = sanal[3][0][0][N1];

  // vectors

  i = 0;			// normal view
  j = 0;			// 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0] = 0.0;
  mmv[i][j][1][0][1] = 1.0;
  // vx1
  mmv[i][j][1][1][0] = -1.0;
  mmv[i][j][1][1][1] = .1;
  // vx2
  mmv[i][j][1][2][0] = -.3;
  mmv[i][j][1][2][1] = .3;
  // vx3
  mmv[i][j][1][3][0] = 0.0;
  mmv[i][j][1][3][1] = 1.26;

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
  mmv[i][j][1][0][1] = .1;
  // rho*vx1
  mmv[i][j][1][1][0] = -0.012;
  mmv[i][j][1][1][1] = 0.001;
  // rho*vx2
  mmv[i][j][1][2][0] = -0.002;
  mmv[i][j][1][2][1] = 0.002;
  // rho*vx3
  mmv[i][j][1][3][0] = 0.0;
  mmv[i][j][1][3][1] = 0.02;

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
  mmv[i][j][1][0][1] = 1.0;
  // vx1
  mmv[i][j][1][1][0] = -1.0;
  mmv[i][j][1][1][1] = .1;
  // vx2
  mmv[i][j][1][2][0] = -.3;
  mmv[i][j][1][2][1] = .3;
  // vx3
  mmv[i][j][1][3][0] = 0.0;
  mmv[i][j][1][3][1] = 1.26;

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
  mmv[i][j][1][0][1] = .1;
  // rho*vx1
  mmv[i][j][1][1][0] = -0.012;
  mmv[i][j][1][1][1] = 0.001;
  // rho*vx2
  mmv[i][j][1][2][0] = -0.002;
  mmv[i][j][1][2][1] = 0.002;
  // rho*vx3
  mmv[i][j][1][3][0] = 0.0;
  mmv[i][j][1][3][1] = 0.02;

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

void accountstoreset(void)
{
  int i, j, k;

  LOOPF {
    if ((i >= intix1) && (i < intox1) && (j >= intix2) && (j < intox2))
      accountstore[k][j][i] = 1;
    else
      accountstore[k][j][i] = 0;
  }

}
