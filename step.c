
#include "global.h"
#include "defs.h"

void stepvar(void)
{
  // source steps

  if (mdotin) {
    injection();
  }

  if (press) {
    step_pgc();			// includes curvature terms from
    // div(rho*v*v)
  }

  if (visc_art == 1) {
    step_visc();
  }

  if (ie) {
    if (wgam)
      step_ie();
  }

  if (visc_real == 1) {
    nu_compute();
    step_visc_real();
  }
  // transport steps
  if (trans == 1) {
    step_trans();
  }

  t += dt;
}

// compute things that are functions of time besides the analytic
// solution

#define OTAU 100.0		// 1 over tau, the time that alpha
				// becomes
				// more like alpha_real0

void tdep_compute(void)
{
  static int firsttime = 1;

  // ////MODIFY alpha_real
  // 
  if (t >= 1000.0)
    alpha_real = alpha_real0;
  if (runtype > 0) {
    if (t >= 0.0)
      alpha_real = alpha_real0;	// go ahead and turn on when
    // using 
    // reentrant data as initial data
  }
  // 
  // ///////////////////////////

  firsttime = 0;
}


void sp_compute(void)
{
  int i, j, k;
  static int firsttime = 1;
  static FTYPE Mach2x1l, Mach2x2l, Mach2totl;
  static FTYPE Mach2x1l_full, Mach2x2l_full, Mach2totl_full;
  static FTYPE Mach2x1h, Mach2x2h, Mach2toth;
  static FTYPE Mach2x1h_full, Mach2x2h_full, Mach2toth_full;
  FTYPE Mach2x1, Mach2x2, Mach2tot;
  char temps[MAXFILENAME];
  FTYPE cs2, vx1, vx2, vx12, vx22, vtot2;
  FTYPE signx1, signx2, signtot;
  static FTYPE tdumpsp;
  int dumpflag;
  static int dumpspc = 0;


  // not setup for seemless reentrance
  if (firsttime == 1) {
    dumpspc = 0;
    tdumpsp = t - 1.0E-6;

    Mach2x1l = Mach2x2l = Mach2totl = 100000.0;
    Mach2x1h = Mach2x2h = Mach2toth = -100000.0;

    if (myid <= 0) {
      Mach2x1l_full = Mach2x2l_full = Mach2totl_full = 100000.0;
      Mach2x1h_full = Mach2x2h_full = Mach2toth_full = -100000.0;

      sprintf(temps, "%s0_logsp%s", DATADIR, ".dat");
      if ((logsp_file = fopen(temps, WRITETYPE)) == NULL) {	// just 
								// 
	// 
	// naively 
	// append 
	// if 
	// appendold==1
	fprintf(fail_file, "error opening sp log output file %s\n",
		temps);
	myexit(1);
      }
      if (appendold == 0) {
	fprintf(logsp_file, "#%10s\n%10d %10d\n", "SPVER", SPVER,
		SPTYPE);
	fprintf(logsp_file,
		"#%16s %16s %16s %16s %16s %16s %16s\n", "time",
		"Mach2x1l", "Mach2x1h", "Mach2x2l", "Mach2x2h",
		"Mach2totl", "Mach2toth");
	fflush(logsp_file);
      }
    }

  }

  if (t >= tdumpsp) {
    dumpflag = 1;
  } else
    dumpflag = 0;


  // ////// Check on sonic point on inner edge
  // 
  k = 0;
  i = 0;			// don't choose intix1 since really
  // want
  // what inner grid is doing, which is
  // really what matters
  for (j = 0; j < N2; j++) {
    if (wgam)
      cs2 = gam * (gam - 1.) * s[2][k][j][i] / s[1][k][j][i];
    else
      cs2 = cs * cs;
    vx1 = v[1][1][k][j][i];
    if (vx1 < 0)
      signx1 = -1.0;
    else
      signx1 = 1.0;
    vx2 = v[1][2][k][j][i];
    if (vx2 < 0)
      signx2 = -1.0;
    else
      signx2 = 1.0;
    signtot = 1;
    vx12 = vx1 * vx1;
    vx22 = vx2 * vx2;
    vtot2 = vx12 + vx22;
    Mach2x1 = signx1 * vx12 / cs2;	// sign tells sign of sqrt of
    // this 
    // value actually
    Mach2x2 = signx2 * vx22 / cs2;	// sign tells sign of sqrt of
    // this 
    // value actually
    Mach2tot = signtot * vtot2 / cs2;


    if (Mach2x1 < Mach2x1l) {
      Mach2x1l = Mach2x1;
    }
    if (Mach2x2 < Mach2x2l) {
      Mach2x2l = Mach2x2;
    }
    if (Mach2tot < Mach2totl) {
      Mach2totl = Mach2tot;
    }

    if (Mach2x1 > Mach2x1h) {
      Mach2x1h = Mach2x1;
    }
    if (Mach2x2 > Mach2x2h) {
      Mach2x2h = Mach2x2;
    }
    if (Mach2tot > Mach2toth) {
      Mach2toth = Mach2tot;
    }

  }

  if (dumpflag == 1) {
    if (numprocs > 1) {
    } else {
      Mach2x1l_full = Mach2x1l;
      Mach2x2l_full = Mach2x2l;
      Mach2totl_full = Mach2totl;

      Mach2x1h_full = Mach2x1h;
      Mach2x2h_full = Mach2x2h;
      Mach2toth_full = Mach2toth;
    }

    if (myid <= 0) {
      fprintf(logsp_file,
	      "%16.10g %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g\n",
	      t, Mach2x1l_full, Mach2x1h_full, Mach2x2l_full,
	      Mach2x2h_full, Mach2totl_full, Mach2toth_full);
      fflush(logsp_file);
    }
    // clean up for next round of check
    if (myid <= 0) {
      Mach2x1l_full = Mach2x2l_full = Mach2totl_full = 100000.0;
      Mach2x1h_full = Mach2x2h_full = Mach2toth_full = -100000.0;
    }
    Mach2x1l = Mach2x2l = Mach2totl = 100000.0;
    Mach2x1h = Mach2x2h = Mach2toth = -100000.0;

    dumpspc++;
    tdumpsp = tstart + (FTYPE) (dumpspc) * DTsp;

  }
  // 
  // //////

  // 
  firsttime = 0;
}



// mass/etc injection routine

// inflows[] only takes what injected on accountable grid, not outside.
// So this doesn't match physical injection rate, which is the true
// reference
void injection(void)
{
  int i, j, k;
  static int firsttime = 1;
  static FTYPE(*rhoi)[N2M][N1M];
  static FTYPE(*rhof)[N2M][N1M];
  FTYPE vxa, vya, vza;
  FTYPE rhoiv[3 + 1], rhofv[3 + 1], massi, massf, die, den, dpot;
  FTYPE drhov[3 + 1], kei[3 + 1], kef[3 + 1], dpd[3 + 1];
  FTYPE dxdyc, dxdy1, dxdy2;
  FTYPE drho, dmass, volume;
  short storeit;

  rhoi = work1;
  rhof = work2;
  // //////////////////////
  // compute mass inflow as function of time over grid
  LOOPRINJ {			// so no need to bound (LOOPF would
    // work
    // too)
    // can't use LOOPFINJ because need rhoi,rhof -1,-1 from
    // original
    // LOOPFINJ modulo the BC
    storeit = accountstore[k][j][i];

    // below 2 only defined once, needed over whole grid due to
    // interpolation for velocities
    rhoi[k][j][i] = s[1][k][j][i];
    drho = rhoinject[k][j][i] * dt;
    rhof[k][j][i] = rhoi[k][j][i] + drho;	// PRECISION

    volume = dvl[1][1][i] * dvl[1][2][j];
    dmass = drho * volume;


    // mass
    if (storeit)
      inflows[1] += dmass;
    s[1][k][j][i] += drho;	// PRECISION
    // ie
    den = eninject[k][j][i] * dt;
    die = den * volume;
    if (storeit)
      inflows[2] += die;
    s[2][k][j][i] += den;

    // pot energy
    dpot = dmass * s[3][k][j][i];
    if (storeit)
      inflows[3] += dpot;


  }
  LOOPVINJ {			// LOOP would work here too if use
    // LOOPF/LOOPH above

    storeit = accountstore[k][j][i];

    dxdyc = dvl[1][1][i] * dvl[1][2][j];
    dxdy1 = dvl[2][1][i] * dvl[1][2][j];
    dxdy2 = dvl[1][1][i] * dvl[2][2][j];



    // for v's
    drhov[3] = rhoinject[k][j][i] * dt;
    drhov[2] = z2e_2(rhoinject[k], j, i) * dt;
    drhov[1] = z2e_1(rhoinject[k], j, i) * dt;

    // zone centered for ke only
    massi = rhoi[k][j][i] * dxdyc;
    massf = massi + drhov[3] * dxdyc;	// since v3 zone centered

    // for v's
    rhoiv[3] = rhoi[k][j][i];
    rhofv[3] = rhof[k][j][i];
    rhoiv[2] = z2e_2(rhoi[k], j, i);
    rhofv[2] = z2e_2(rhof[k], j, i);
    rhoiv[1] = z2e_1(rhoi[k], j, i);
    rhofv[1] = z2e_1(rhof[k], j, i);

    // initial ke
    vxa = e2z_1(v[1][1][k], j, i);
    vya = e2z_2(v[1][2][k], j, i);
    vza = v[1][3][k][j][i];	// When fake-3d

    kei[1] = 0.5 * massi * vxa * vxa;
    kei[2] = 0.5 * massi * vya * vya;
    kei[3] = 0.5 * massi * vza * vza;

    // vz
    dpd[3] = VZFRACT * drhov[3] * pow(g[2][3][i], -0.5);	// s3=r*v
    if (storeit)
      inflows[7] += g[2][3][i] * g[2][4][j] * dpd[3] * dxdyc;	// r*sin(theta)*m*v(angmom3)
    v[1][3][k][j][i] = (rhoiv[3] * v[1][3][k][j][i] + dpd[3]) / rhofv[3];	// cons 
										// 
    // 
    // of 
    // L, 
    // but 
    // r*sint 
    // same

    // vy
    dpd[2] = 0.0;		// s2=r*v (interpolated to vx2s edge)
    if (storeit)
      inflows[8] += g[2][3][i] * dpd[2] * dxdy2;	// r*m*v
    // (angmom2)
    v[1][2][k][j][i] = (rhoiv[2] * v[1][2][k][j][i] + dpd[2]) / rhofv[2];	// cons 
										// 
    // 
    // of 
    // L, 
    // but 
    // r 
    // same


    // vx
    dpd[1] = 0.0;		// s1=r*v (interpolated to vx1s edge)
    if (storeit)
      inflows[9] += dpd[1] * dxdy1;	// m*v(angmom1)
    v[1][1][k][j][i] = (rhoiv[1] * v[1][1][k][j][i] + dpd[1]) / rhofv[1];	// cons 
										// 
    // 
    // of 
    // L, 
    // but 
    // r 
    // same

    // fprintf(fail_file,"%d %d %d %15.10g
    // %15.10g\n",k,j,i,rhoiv[2]/rhofv[2],rhoiv[1]/rhofv[1]);

    // ke
    // compute final values
    vxa = e2z_1(v[1][1][k], j, i);
    vya = e2z_2(v[1][2][k], j, i);
    vza = v[1][3][k][j][i];	// When fake-3d

    kef[1] = 0.5 * massf * vxa * vxa;
    kef[2] = 0.5 * massf * vya * vya;
    kef[3] = 0.5 * massf * vza * vza;

    if (storeit)
      inflows[NUMSCA + 1] +=
	  (kef[3] - kei[3]) + (kef[2] - kei[2]) + (kef[1] - kei[1]);

  }
  // have to bound vectors if multiple cpus, since data changed
  // across
  // boundaries
  if (numprocs > 1) {		// could even isolate bound to
    // necessary
    // transfers only
    bound(NULL, NULL, 0, 1);	// just velocity for now
  }
  // fflush(fail_file);

  // 
  // /////////////////////


  firsttime = 0;
}



// must be consistent with nu_fact calc in init.c
void nu_compute(void)
{
  int i, j, k;
  FTYPE ftemp;
  FTYPE ftemp2;
  FTYPE cs2;

  ftemp = sqrt(gam * (gam - 1.));
  ftemp2 = gam * (gam - 1.);

  // no need to bound since can do LOOPF since nu_real only depends
  // on
  // zone centered quantities

  if (vreal == 1) {
    LOOPF {
      // alpha-param
      if (wgam)
	cs = ftemp * sqrt(s[2][k][j][i] / s[1][k][j][i]);
      nu_real[k][j][i] = alpha_real * nu_fact[k][j][i] * cs;
    }
  } else if (vreal == 2) {
    LOOPF {
      if (wgam)
	cs2 = ftemp2 * s[2][k][j][i] / s[1][k][j][i];
      else
	cs2 = cs * cs;
      // Igumenshchev alpha param: $\nu=\alpha*cs^2/Omega_{k}$
      nu_real[k][j][i] = alpha_real * nu_fact[k][j][i] * cs2;
    }
  } else if (vreal == 3) {
    LOOPF {
      // Stone et al.
      // make consistent with nu_fact in init.c
      // no need for nu_fact since just constant for now

      // run A through I'
      nu_real[k][j][i] = alpha_real * s[1][k][j][i] * nu_fact[k][j][i];

      // run J
      // nu_real[k][j][i]=alpha_real*s[1][k][j][i]*nu_fact[k][j][i];

      // run K
      // nu_real[k][j][i]=alpha_real*nu_fact[k][j][i];
    }
  } else if (vreal == 4) {
    LOOPF {
      if (wgam)
	cs2 = ftemp2 * s[2][k][j][i] / s[1][k][j][i];
      else
	cs2 = cs * cs;
      // Igumenshchev alpha param: $\nu=\alpha*cs^2/Omega_{k}$
      nu_real[k][j][i] = alpha_real * nu_fact[k][j][i] * cs2;
    }
  }
  if (analoutput == 6) {
    LOOPF {
      // for test of this code
      nu_real[k][j][i] = nu_fact[k][j][i];
    }
  }

}


// pressure and gravity step 

// Non-B part checked by Jon
void step_pgc(void)
{

  FTYPE rhoa, bxa, bxp, bxm, bya, byp, bym, bza, bzp, bzm;
  int i, j, k;
  FTYPE tempfx1 = 0, tempfx2 = 0, tempfs = 0;

  // first pressure gradient
  LOOPV1 {
    // x1
#if(PDEN)
#if(RELIE)
    // might be problem with division of linear interp
    tempfx1 =
	-(s[1][k][j][i] * s[2][k][j][i] -
	  s[1][k][j][i - 1] * s[2][k][j][i -
					 1]) * OARCL[1][1][k][j][i] /
	z2e_1(s[1][k], j, i);
#else
    if (wgam) {
      tempfx1 =
	  -(gam - 1.) * (s[2][k][j][i] -
			 s[2][k][j][i -
				    1]) * OARCL[1][1][k][j][i] /
	  z2e_1(s[1][k], j, i);
    } else {
      tempfx1 =
	  -cs * cs * (s[1][k][j][i] -
		      s[1][k][j][i -
				 1]) * OARCL[1][1][k][j][i] /
	  z2e_1(s[1][k], j, i);
    }
#endif
#else
    tempfx1 = 0;
#endif

#if(PGRAV)
    tempfx1 += -(s[3][k][j][i] - s[3][k][j][i - 1]) * OARCL[1][1][k][j][i];	// grav 
										// 
    // 
    // from 
    // pot

#endif

#if(COORD==3)
#if(CURVE==1)
    tempfs = e2e_v2(v[1][2][k], j, i);
    tempfs *= tempfs;
    tempfx1 += tempfs / g[1][2][i];	// * dg[1][2][i]=1
    tempfs = z2e_1(v[1][3][k], j, i);
    tempfs *= tempfs;
    tempfx1 += tempfs / g[1][3][i];	// * dg[1][3][i]=1
#endif
#endif

    v[1][1][k][j][i] += dt * tempfx1;

  }


  if (N2M > 1) {

    LOOPV2 {
      // x2
#if(PDEN==1)
#if(RELIE)
      // might be problem with division of linear interp
      tempfx2 =
	  -(s[1][k][j][i] * s[2][k][j][i] -
	    s[1][k][j - 1][i] * s[2][k][j -
					1][i]) *
	  OARCL[1][2][k][j][i] / z2e_2(s[1][k], j, i);
#else
      if (wgam) {
	tempfx2 =
	    -(gam - 1.) * (s[2][k][j][i] -
			   s[2][k][j -
				   1][i]) * OARCL[1][2][k][j][i] /
	    z2e_2(s[1][k], j, i);
      } else {
	tempfx2 =
	    -cs * cs * (s[1][k][j][i] -
			s[1][k][j -
				1][i]) * OARCL[1][2][k][j][i] /
	    z2e_2(s[1][k], j, i);
      }
#endif

#else
      tempfx2 = 0;
#endif

#if(PGRAV)
      tempfx2 +=
	  -(s[3][k][j][i] - s[3][k][j - 1][i]) * OARCL[1][2][k][j][i];
#endif

#if(CURVE==1)
#if( (COORD==2)||(COORD==3))
      tempfs = z2e_2(v[1][3][k], j, i);
      tempfs *= tempfs;
      tempfx2 += tempfs / (g[2][2][i] * g[1][4][j]) * dg[1][4][j];
#endif
#endif

      v[1][2][k][j][i] += dt * tempfx2;

    }
  }
  // below not correct for reflecting/N2-size conditions
  if (pmag == 1) {
    // magnetic pressure gradient 
    LOOP {
      rhoa = 0.5 * (s[1][k][j][i] + s[1][k][j][i - 1]);

      bya =
	  0.25 * (v[2][2][k][j][i] + v[2][2][k][j][i - 1] +
		  v[2][2][k][j + 1][i - 1] + v[2][2][k][j + 1][i]);
      byp = 0.5 * (v[2][2][k][j][i] + v[2][2][k][j + 1][i]);
      bym = 0.5 * (v[2][2][k][j][i - 1] + v[2][2][k][j + 1][i - 1]);

      bza = 0.5 * (v[2][3][k][j][i] + v[2][3][k][j][i - 1]);
      bzp = v[2][3][k][j][i];
      bzm = v[2][3][k][j][i - 1];

      v[1][1][k][j][i] +=
	  -(dt / (dx[2][1][i] * rhoa)) * (bya * (byp - bym) +
					  bza * (bzp - bzm)
	  );

#if(N2M>1)
      rhoa = 0.5 * (s[1][k][j][i] + s[1][k][j - 1][i]);

      bxa =
	  0.25 * (v[2][1][k][j][i] + v[2][1][k][j - 1][i] +
		  v[2][1][k][j - 1][i + 1] + v[2][1][k][j][i + 1]);
      bxp = 0.5 * (v[2][1][k][j][i] + v[2][1][k][j][i + 1]);
      bxm = 0.5 * (v[2][1][k][j - 1][i] + v[2][1][k][j - 1][i + 1]);

      bza = 0.5 * (v[2][3][k][j][i] + v[2][3][k][j - 1][i]);
      bzp = v[2][3][k][j][i];
      bzm = v[2][3][k][j - 1][i];


      v[1][2][k][j][i] +=
	  -(dt / (dx[2][2][j] * rhoa)) * (bxa * (bxp - bxm) +
					  bza * (bzp - bzm)
	  );
#endif
    }
  }
  // only needs to be done once 
  bound(NULL, NULL, 0, 1);	// step_pgc()


}

#define CRAPOLA 1.0
#define CRAPOLA2 1.0

// Checked by Jon
void step_visc(void)
{
  static FTYPE(*visc)[N3M][N2M][N1M];
  static FTYPE(*delv)[N2M][N1M];
  static FTYPE(*gradv)[N3M][N2M][N1M];
  static FTYPE(*l2_ten)[N2M][N1M];
  FTYPE dvx, dvy, qlx, qvnrx, qly, qvnry;
  int i, j, k;
  FTYPE ftemp;
  static FTYPE(*dv)[N3M][N2M][N1M];

  visc = workv1;
  dv = workv2;

  gradv = workv3;
  delv = work1;
  l2_ten = work2;

  // find x1-dir viscous stresses 
  LOOPVISC {			// no need to bound visc, all good

    // del v, at zone center
    dvx = dv[1][k][j][i] = v[1][1][k][j][i + 1] - v[1][1][k][j][i];
    // del v, at zone center 
    dvy = dv[2][k][j][i] = v[1][2][k][j + 1][i] - v[1][2][k][j][i];

    // linear viscosity
#if(VISC_LINEAR)
    if (wgam) {
      cs = sqrt(gam * (gam - 1.) * s[2][k][j][i] / s[1][k][j][i]);
    }
    qlx = -nu_l * s[1][k][j][i] * cs * dvx;
    qly = -nu_l * s[1][k][j][i] * cs * dvy;
#else
    qlx = qly = 0;
#endif

    // von neumann,richtmyer viscosity 
    if (dvx < 0) {
      qvnrx = nu_vnr * s[1][k][j][i] * dvx * dvx;
    } else
      qvnrx = 0.;

    // von neumann,richtmyer viscosity 
    if (dvy < 0) {
      qvnry = nu_vnr * s[1][k][j][i] * dvy * dvy;
    } else
      qvnry = 0.;

    visc[1][k][j][i] = qlx + qvnrx;
    visc[2][k][j][i] = qly + qvnry;

  }

  // update velocity, internal energy 
  if (ie) {
    if (wgam) {
      LOOP {

	s[2][k][j][i] +=
	    -dt * (visc[1][k][j][i] * dv[1][k][j][i] *
		   OARCL[2][1][k][j][i] +
		   visc[2][k][j][i] * dv[2][k][j][i] *
		   OARCL[3][2][k][j][i]);

#if(FORCEIE)
	if (s[2][k][j][i] < IEFLOOR) {
	  ftemp = (IEFLOOR - s[2][k][j][i]);
	  if ((i >= intix1) && (i < intox1) && (j >= intix2)
	      && (j < intox2))
	    floors[2] += ftemp * dvl[1][1][i] * dvl[1][2][j];
#if(FLOORDUMPFLAG==1)
	  floorvars[2][k][j][i] += ftemp;
#endif
#if(DOFLOORDIAG==1)
	  floorcnt[0][2]++;
	  if (s[2][k][j][i] < floorlowest[2]) {
	    floorlowest[2] = s[2][k][j][i];
	    wherelowest[2] = 0;
	  }
#endif
#if(DOFLOORD2==1)
	  fprintf(logfl_file,
		  "corrected en in step_visc: t: %15.10g %d %d %d %15.10g\n",
		  t, k, j, i, (IEFLOOR - s[2][k][j][i]));
#endif
	  s[2][k][j][i] = IEFLOOR;
	}
#endif
      }
    }
  }

  LOOPV1 {


    v[1][1][k][j][i] +=
	-dt * (visc[1][k][j][i] -
	       visc[1][k][j][i -
			     1]) * OARCL[1][1][k][j][i] /
	z2e_1(s[1][k], j, i);


  }

  if (N2M > 1) {
    LOOPV2 {



      v[1][2][k][j][i] +=
	  -dt * (visc[2][k][j][i] -
		 visc[2][k][j -
			    1][i]) * OARCL[1][2][k][j][i] /
	  z2e_2(s[1][k], j, i);

    }
  }
  // bound everything
  if ((ie) && (wgam)) {
    bound(NULL, NULL, 2, 1);	// both velocities and ie
    // (step_visc)
  } else
    bound(NULL, NULL, 0, 1);	// otherwise just velocities
  // (step_visc)
}


#define DIFFTYPE 0
#define VREALHARDCODE 0

/* 
   I don't compute momentum boundary flux due to v_r and v_theta
   components since they involve volume terms for spc */

void step_visc_real(void)
{
  static FTYPE(*sigma)[3][N3M][N2M][N1M];	// -2*rho*nu*e_{ij}=
  // sigma_{ij}
  static FTYPE(*rost)[3][N3M][N2M][N1M];	// e_{ij}
  static FTYPE(*rostnu)[3][N3M][N2M][N1M];	// e_{ij}*nu
  static FTYPE(*nurho_real)[N2M][N1M];

  static FTYPE(*delv)[N2M][N1M];	// deldotv
  int i, j, k;
  FTYPE ftemp;
  FTYPE subftemp;
  FTYPE flux;

  sigma = workt1;
  rost = workt2;
  rostnu = workt3;
  delv = work1;
  nurho_real = work2;

  compute_sigma(sigma, rost, rostnu, nurho_real, delv);


  // update velocity, internal energy 
  if (vischeat) {
    if ((ie) || (analoutput == 6)) {
      if (wgam) {
	LOOP {
	  ftemp = 0;
	  // this is some-what expensive, but I optimized for
	  // velocity terms and sigma, not this internal
	  // energy
	  // term!

#if(VISCE11)
	  ftemp += rost[1][1][k][j][i] * rost[1][1][k][j][i];
#endif

#if(VISCE22)
	  ftemp += rost[2][2][k][j][i] * rost[2][2][k][j][i];
#endif

#if(VISCE33)
	  ftemp += rost[3][3][k][j][i] * rost[3][3][k][j][i];
#endif

#if(VISCE13)
	  subftemp = e2z_1(rost[1][3][k], j, i);
	  // #if(ANALOUTPUT==6)
	  // v[1][2][k][j][i]=sigma[1][3][k][j][i];
	  // #endif 
	  ftemp += 2.0 * (subftemp * subftemp);	// 2
	  // accounts 
	  // for
	  // rost[3][1]
#endif

#if(VISCE23)
	  subftemp = e2z_2(rost[2][3][k], j, i);
	  ftemp += 2.0 * (subftemp * subftemp);	// 2
	  // accounts 
	  // for
	  // rost[3][2]
#endif

#if(VISCE12)
	  subftemp = c2z(rost[1][2][k], j, i);
	  ftemp += 2.0 * (subftemp * subftemp);	// 2
	  // accounts 
	  // for
	  // rost[2][1]
#endif

	  s[2][k][j][i] += dt * nurho_real[k][j][i] * ftemp;

#if(FORCEIE)
	  if (s[2][k][j][i] < IEFLOOR) {
	    ftemp = (IEFLOOR - s[2][k][j][i]);
	    if ((i >= intix1) && (i < intox1) && (j >= intix2)
		&& (j < intox2))
	      floors[2] += ftemp * dvl[1][1][i] * dvl[1][2][j];
#if(FLOORDUMPFLAG==1)
	    floorvars[2][k][j][i] += ftemp;
#endif
#if(DOFLOORDIAG==1)
	    floorcnt[5][2]++;
	    if (s[2][k][j][i] < floorlowest[2]) {
	      floorlowest[2] = s[2][k][j][i];
	      wherelowest[2] = 5;
	    }
#endif
#if(DOFLOORD2==1)
	    fprintf(logfl_file,
		    "corrected en in step_visc_real: t: %15.10g %d %d %d %15.10g\n",
		    t, k, j, i, (IEFLOOR - s[2][k][j][i]));
#endif
	    s[2][k][j][i] = IEFLOOR;
	  }
#endif
	}
      }
    }				// endif doing ie
  }				// endif vischeat==1


  // now right before going to advection, compute fluxes of energy
  // and
  // momentum due to viscosity

  if (DOLOSSDIAG) {

    if (vischeat) {
      // capture flux of internal energy

      if (VISCE31) {
	// only flux is radial!
	if ((reflectix1 == 0) || 1) {
	  k = intix3;
	  i = intix1;
	  for (j = intix2; j < intox2; j++) {
	    flux =
		dt * sigma[3][1][k][j][i] * z2e_1(v[1][3][k],
						  j,
						  i) *
		DS[1][1][k][j][i];
	    lossvisc[1][1][0][j] += -flux;
	  }
	}
	if ((reflectox1 == 0) || 1) {
	  k = intix3;
	  i = intox1;
	  for (j = intix2; j < intox2; j++) {
	    flux =
		dt * sigma[3][1][k][j][i] * z2e_1(v[1][3][k],
						  j,
						  i) *
		DS[1][1][k][j][i];
	    lossvisc[1][1][1][j] += flux;
	  }
	}
      }

      if (VISCE32) {
	if ((reflectix2 == 0) || 1) {
	  k = intix3;
	  j = intix2;
	  for (i = intix1; i < intox1; i++) {
	    flux =
		dt * sigma[3][2][k][j][i] * z2e_2(v[1][3][k],
						  j,
						  i) *
		DS[1][2][k][j][i];
	    lossvisc[1][2][0][i] += -flux;
	  }
	}
	if ((reflectox2 == 0) || 1) {
	  k = intix3;
	  j = intox2;
	  for (i = intix1; i < intox1; i++) {
	    flux =
		dt * sigma[3][2][k][j][i] * z2e_2(v[1][3][k],
						  j,
						  i) *
		DS[1][2][k][j][i];
	    lossvisc[1][2][1][i] += flux;
	  }
	}
      }

      if (VISCE11) {
	if ((reflectix1 == 0) || 1) {
	  k = intix3;
	  i = intix1;
	  for (j = intix2; j < intox2; j++) {
	    flux =
		dt * z2e_1(sigma[1][1][k], j,
			   i) * v[1][1][k][j][i] * DS[1][1][k][j][i];
	    lossvisc[1][1][0][j] += -flux;
	  }
	}
	if ((reflectox1 == 0) || 1) {
	  k = intix3;
	  i = intox1;
	  for (j = intix2; j < intox2; j++) {
	    flux =
		dt * z2e_1(sigma[1][1][k], j,
			   i) * v[1][1][k][j][i] * DS[1][1][k][j][i];
	    lossvisc[1][1][1][j] += flux;
	  }
	}
      }

      if (VISCE22) {
	if ((reflectix2 == 0) || 1) {
	  k = intix3;
	  j = intix2;
	  for (i = intix1; i < intox1; i++) {
	    flux =
		dt * z2e_2(sigma[2][2][k], j,
			   i) * v[1][2][k][j][i] * DS[1][2][k][j][i];
	    lossvisc[1][2][0][i] += -flux;
	  }
	}
	if ((reflectox2 == 0) || 1) {
	  k = intix3;
	  j = intox2;
	  for (i = intix1; i < intox1; i++) {
	    flux =
		dt * z2e_2(sigma[2][2][k], j,
			   i) * v[1][2][k][j][i] * DS[1][2][k][j][i];
	    lossvisc[1][2][1][i] += flux;
	  }
	}
      }

      if (VISCE21) {
	// only flux is radial!
	if ((reflectix1 == 0) || 1) {
	  k = intix3;
	  i = intix1;
	  for (j = intix2; j < intox2; j++) {
	    flux =
		dt * sigma[2][1][k][j][i] * z2e_1(v[1][2][k],
						  j,
						  i) *
		DS[3][1][k][j][i];
	    lossvisc[1][1][0][j] += -flux;
	  }
	}
	if ((reflectox1 == 0) || 1) {
	  k = intix3;
	  i = intox1;
	  for (j = intix2; j < intox2; j++) {
	    flux =
		dt * sigma[2][1][k][j][i] * z2e_1(v[1][2][k],
						  j,
						  i) *
		DS[3][1][k][j][i];
	    lossvisc[1][1][1][j] += flux;
	  }
	}
      }

      if (VISCE12) {
	if ((reflectix2 == 0) || 1) {
	  k = intix3;
	  j = intix2;
	  for (i = intix1; i < intox1; i++) {
	    flux =
		dt * sigma[1][2][k][j][i] * z2e_2(v[1][1][k],
						  j,
						  i) *
		DS[2][2][k][j][i];
	    lossvisc[1][2][0][i] += -flux;
	  }
	}
	if ((reflectox2 == 0) || 1) {
	  k = intix3;
	  j = intox2;
	  for (i = intix1; i < intox1; i++) {
	    flux =
		dt * sigma[1][2][k][j][i] * z2e_2(v[1][1][k],
						  j,
						  i) *
		DS[2][2][k][j][i];
	    lossvisc[1][2][1][i] += flux;
	  }
	}
      }
    }
    // capture flux of angular momentum(only l-phi for now)
    // don't do other momentum dirs since they have non-surface
    // type
    // integrals
    if (VISCE13) {
      // radial direction
      if ((reflectix1 == 0) || 1) {
	k = intix3;
	i = intix1;
	for (j = intix2; j < intox2; j++) {
	  ftemp = g[1][3][i] * g[2][4][j];
	  flux = dt * ftemp * sigma[1][3][k][j][i] * DS[1][1][k][j][i];
	  lossv[1][3][1][0][j] += -flux;
	}
      }
      if ((reflectox1 == 0) || 1) {
	k = intix3;
	i = intox1;
	for (j = intix2; j < intox2; j++) {
	  ftemp = g[1][3][i] * g[2][4][j];
	  flux = dt * ftemp * sigma[1][3][k][j][i] * DS[1][1][k][j][i];
	  lossv[1][3][1][1][j] += flux;
	}
      }
    }

    if (VISCE23) {
      // theta direction

      if ((reflectix2 == 0) || 1) {
	k = intix3;
	j = intix2;
	for (i = intix1; i < intox1; i++) {
	  ftemp = g[2][3][i] * g[1][4][j];
	  flux = dt * ftemp * sigma[2][3][k][j][i] * DS[1][2][k][j][i];
	  lossv[1][3][2][0][i] += -flux;
	}
      }
      if ((reflectox2 == 0) || 1) {
	k = intix3;
	j = intox2;
	for (i = intix1; i < intox1; i++) {
	  ftemp = g[2][3][i] * g[1][4][j];
	  flux = dt * ftemp * sigma[2][3][k][j][i] * DS[1][2][k][j][i];
	  lossv[1][3][2][1][i] += flux;
	}
      }
      // else leave at 0 and don't add anything
    }
  }

  if (VISCE11 || VISCE12 || VISCE22 || VISCE33) {
    LOOPV1 {
#if(COORD==3)

      v[1][1][k][j][i] +=
	  -dt / (z2e_1(s[1][k], j, i)) *
	  ((x[2][1][i] * x[2][1][i] * sigma[1][1][k][j][i] -
	    x[2][1][i - 1] * x[2][1][i - 1] * sigma[1][1][k][j][i -
								1]) /
	   dvl[2][1][i] +
	   (g[1][4][j + 1] * sigma[1][2][k][j + 1][i] -
	    g[1][4][j] * sigma[1][2][k][j][i]) / (g[1][3][i] *
						  dvl[1][2][j]) -
	   (z2e_1(sigma[2][2][k], j, i) +
	    z2e_1(sigma[3][3][k], j, i)) / x[1][1][i]);

#endif
    }
  }

  if (VISCE12 || VISCE22 || VISCE33) {
    if (N2M > 1) {
      LOOPV2 {
#if(COORD==3)

	v[1][2][k][j][i] +=
	    -dt / (z2e_2(s[1][k], j, i)) *
	    ((x[1][1][i + 1] * x[1][1][i + 1] * x[1][1][i + 1] *
	      sigma[1][2][k][j][i + 1] -
	      x[1][1][i] * x[1][1][i] * x[1][1][i] *
	      sigma[1][2][k][j][i]) / (x[2][1][i] * dvl[1][1][i]) +
	     (g[2][4][j] * sigma[2][2][k][j][i] -
	      g[2][4][j - 1] * sigma[2][2][k][j -
					      1][i]) /
	     (g[2][3][i] * dvl[2][2][j]) - z2e_2(sigma[3][3][k], j,
						 i) * cotan[1][j] /
	     x[2][1][i]);
#endif
      }
    }
  }

  if (VISCE13 || VISCE23) {
    LOOPV3 {			// absorbed sigma[1][3] terms together
#if(COORD==3)
      v[1][3][k][j][i] +=
#if(DIFFTYPE==0)
	  -dt / (s[1][k][j][i]) *
	  ((x[1][1][i + 1] * x[1][1][i + 1] * x[1][1][i + 1] *
	    sigma[1][3][k][j][i + 1] -
	    x[1][1][i] * x[1][1][i] * x[1][1][i] *
	    sigma[1][3][k][j][i]) / (x[2][1][i] * dvl[1][1][i]) +
	   (g[1][4][j + 1] * g[1][4][j + 1] *
	    sigma[2][3][k][j + 1][i] -
	    g[1][4][j] * g[1][4][j] * sigma[2][3][k][j][i]) /
	   (g[2][3][i] * g[2][4][j] * dvl[1][2][j]));
#endif
#if(DIFFTYPE==1)
      -dt * (-2.0 *
	     (x[1][1][i + 1] * x[1][1][i + 1] * x[1][1][i + 1] *
	      rostnu[1][3][k][j][i + 1] -
	      x[1][1][i] * x[1][1][i] * x[1][1][i] *
	      rostnu[1][3][k][j][i]) / (x[2][1][i] * dvl[1][1][i])
#if(VREALHARDCODE==3)
	     - 2.0 * alpha_real * e2z_1(rost[1][3][k], j,
					i) * (z2e_1(s[1][k], j,
						    i + 1) -
					      z2e_1(s[1][k], j,
						    i)) / dx[1][1][i]
#else
	     - 2.0 * e2z_1(rostnu[1][3][k], j,
			   i) / s[1][k][j][i] * (z2e_1(s[1][k], j,
						       i + 1) -
						 z2e_1(s[1][k], j,
						       i)) / dx[1][1][i]
#endif
	     +
	     (g[1][4][j + 1] * g[1][4][j + 1] *
	      sigma[2][3][k][j + 1][i] -
	      g[1][4][j] * g[1][4][j] * sigma[2][3][k][j][i]) /
	     (g[2][3][i] * g[2][4][j] * dvl[1][2][j]));
#endif
#endif
    }
  }
  // bound everything
  if ((ie) && (wgam)) {
    bound(NULL, NULL, 2, 1);	// both velocities and ie
    // (step_visc)
  } else
    bound(NULL, NULL, 0, 1);	// otherwise just velocities
  // (step_visc)

  if (analoutput == 6) {
    bound(NULL, NULL, 2, 1);	// for real visc test bound all
  }
}


// assume nu_real, s, v, and geom terms have already been computed as
// per
// timestep or initially
void compute_sigma(FTYPE(*sigma)[3][N3M][N2M][N1M],
		   FTYPE(*rost)[3][N3M][N2M][N1M],
		   FTYPE(*rostnu)[3][N3M][N2M][N1M],
		   FTYPE(*nurho_real)[N2M][N1M], FTYPE(*delv)[N2M][N1M])
{
  int i, j, k;
  FTYPE ftemp;

  // get important coefficient(nu*rho)
  LOOPF {			// needs full loop for interp for eij
    // on
    // half-full loop 
    nurho_real[k][j][i] = 2.0 * s[1][k][j][i] * nu_real[k][j][i];	// 2*rho*nu 
									// 
    // 
    // really
  }

  // find sigma=-2*rho*nu*e_{ij}
  LOOPH {
    delv[k][j][i] = deldotv(1, k, j, i);	// no need to interp
    // since 
    // this centered and used
    // to compute only
    // centered sigma's

    // only 2D versions
#if(COORD==1)
    fprintf(fail_file, "no cart real visc yet: %d %d\n", visc_real,
	    vreal);
    myexit(1);
#elif(COORD==2)
    fprintf(fail_file, "no cylindrical real visc yet: %d %d\n",
	    visc_real, vreal);
    myexit(1);
#elif(COORD==3)


#if(VISCE11)
    ftemp =
	((v[1][1][k][j][i + 1] -
	  v[1][1][k][j][i]) * OARCL[2][1][k][j][i] -
	 THIRD * delv[k][j][i]);
    rost[1][1][k][j][i] = ftemp;
    rostnu[1][1][k][j][i] = ftemp * nu_real[k][j][i];
    sigma[1][1][k][j][i] = -nurho_real[k][j][i] * ftemp;
#endif

#if(VISCE22)
    ftemp =
	((v[1][2][k][j + 1][i] -
	  v[1][2][k][j][i]) * OARCL[3][2][k][j][i] + e2z_1(v[1][1][k],
							   j,
							   i) /
	 g[2][2][i] - THIRD * delv[k][j][i]);
    rost[2][2][k][j][i] = ftemp;
    rostnu[2][2][k][j][i] = ftemp * nu_real[k][j][i];
    sigma[2][2][k][j][i] = -nurho_real[k][j][i] * ftemp;
#endif

#if(VISCE33)
    ftemp =
	((e2z_1(v[1][1][k], j, i) +
	  e2z_2(v[1][2][k], j,
		i) * cotan[2][j]) / x[2][1][i] - THIRD * delv[k][j][i]);
    rost[3][3][k][j][i] = ftemp;
    rostnu[3][3][k][j][i] = ftemp * nu_real[k][j][i];
    sigma[3][3][k][j][i] = -nurho_real[k][j][i] * ftemp;
#endif

#if(VISCE12)
    ftemp =
	0.5 * (x[1][1][i] *
	       (v[1][2][k][j][i] / x[2][1][i] -
		v[1][2][k][j][i - 1] / x[2][1][i -
					       1]) *
	       OARCL[3][1][k][j][i] + (v[1][1][k][j][i] -
				       v[1][1][k][j -
						  1][i]) *
	       OARCL[2][2][k][j][i]);
    rost[1][2][k][j][i] = rost[2][1][k][j][i] = ftemp;
    rostnu[1][2][k][j][i] = rostnu[2][1][k][j][i] =
	ftemp * z2c(nu_real[k], j, i);
    sigma[1][2][k][j][i] = sigma[2][1][k][j][i] =
	-z2c(nurho_real[k], j, i) * ftemp;
#endif

#if(VISCE13)
    ftemp =
	0.5 * (x[1][1][i] *
	       (v[1][3][k][j][i] / x[2][1][i] -
		v[1][3][k][j][i - 1] / x[2][1][i -
					       1]) *
	       OARCL[1][1][k][j][i]);
    rost[1][3][k][j][i] = rost[3][1][k][j][i] = ftemp;
    rostnu[1][3][k][j][i] = rostnu[3][1][k][j][i] =
	ftemp * z2e_1(nu_real[k], j, i);
    sigma[1][3][k][j][i] = sigma[3][1][k][j][i] =
	-z2e_1(nurho_real[k], j, i) * ftemp;
#endif

#if(VISCE23)
    ftemp =
	(0.5 * g[1][4][j] * OARCL[1][2][k][j][i]) * (v[1][3][k][j][i] /
						     g[2][4][j] -
						     v[1][3][k][j -
								1][i] /
						     g[2][4][j - 1]);
    rost[2][3][k][j][i] = rost[3][2][k][j][i] = ftemp;
    rostnu[2][3][k][j][i] = rostnu[3][2][k][j][i] =
	ftemp * z2e_2(nu_real[k], j, i);
    sigma[2][3][k][j][i] = sigma[3][2][k][j][i] =
	-z2e_2(nurho_real[k], j, i) * ftemp;
#endif

    // #if(ANALOUTPUT==6)
    // v[1][1][k][j][i]=sigma[1][3][k][j][i]; // for testing real
    // visc
    // #endif

#endif				// endif coord==3


  }				// end e_{ij} generation

}



// Checked by Jon
void step_ie(void)
{
  FTYPE dv, dvg;
  int i, j, k;
  FTYPE ftemp;


  // update internal energy 
  LOOP {
    dv = deldotv(1, k, j, i);
    dvg = 0.5 * dt * (gam - 1.0) * dv;
    s[2][k][j][i] *= (1. - dvg) / (1. + dvg);

    // s[2][k][j][i]+=dt*(gam-1.)*s[2][k][j][i]*dv;
#if(FORCEIE)
    if (s[2][k][j][i] < IEFLOOR) {
      ftemp = (IEFLOOR - s[2][k][j][i]);
      if ((i >= intix1) && (i < intox1) && (j >= intix2)
	  && (j < intox2))
	floors[2] += ftemp * dvl[1][1][i] * dvl[1][2][j];
#if(FLOORDUMPFLAG==1)
      floorvars[2][k][j][i] += ftemp;
#endif
#if(DOFLOORDIAG==1)
      floorcnt[1][2]++;
      if (s[2][k][j][i] < floorlowest[2]) {
	floorlowest[2] = s[2][k][j][i];
	wherelowest[2] = 1;
      }
#endif
#if(DOFLOORD2==1)
      fprintf(logfl_file,
	      "corrected en in step_ie: t: %15.10g %d %d %d %15.10g\n",
	      t, k, j, i, (IEFLOOR - s[2][k][j][i]));
#endif
      s[2][k][j][i] = IEFLOOR;
    }
#endif
  }

  bound(NULL, NULL, 2, 0);

}



// Do transport step.  Order of sweep direction is
// varied from step to step. 


void step_trans(void)
{
  static int nstep = 0;


  if (nstep % 2 == 0) {
    if (transx1) {
      sweepx();
    }
    if (transx2) {
      sweepy();
    }
  } else {
    if (transx2) {
      sweepy();
    }
    if (transx1) {
      sweepx();
    }
  }
  nstep++;
}
