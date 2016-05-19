// Checked by Jon

#include "global.h"
#include "defs.h"


/* Van Leer formula on staggered grid(scalars on zone centers, vectors
   on zone faces):

   scalar and non-directional vectors:

   S[i]* = S[i-1] + 0.5*dxb[i]*dq[i-1]-0.5*v[i]*dt*dq[i-1] if v[i]>0
   S[i]* = S[i] - 0.5*dxb[i]*dq[i]-0.5*v[i]*dt*dq[i] if v[i]<0

   directional vectors:

   V[i]* = V[i] + 0.5*dxa[i]*dq[i]-0.5*((v[i+1]+v[i])/2.0)*dt*dq[i] if
   ((v[i+1]+v[i])/2.0) > 0 V[i]* = V[i+1] -
   0.5*dxa[i]*dq[i+1]+0.5*((v[i+1]+v[i])/2.0)*dt*dq[i+1] if
   ((v[i+1]+v[i])/2.0) < 0

   Note: dq has 2* in it which cancels with one .5 above

 */

/* 
   Notes:

   Rate of convergence for cart is at best == 2 for nearly constant
   velocity

   Rate of convergence for spc is at best == 1.5 for nearly constant
   velocity.  Why not 2? Both have ROC of 1 for general velocity
   Perhaps could do double velocity loop in advection to see if can get 
   better ROC when vel not constant

   The floor may be implemented better somehow

 */

/* 
   Notes on what needs what: x1-dir: mdot= [s2-1,N2-1],[s1-1,N1]
   en/rho/vz/Bz= fl=[0,N2-1],[0,N1] , mdot=fl, dq=[fl][-1,N1] vx =
   fl=[0,N2-1],[s1-1,N1-1], mdot=[0,N2-1],[s1-1,N1], dq=[fl][-1,N1] vy
   = fl=[s2,N2-1],[0,N1], mdot=[s2-1,N2-1],[0,N1], dq=[fl][-1,N1]

   x2-dir: mdot= [s2-1,N2],[s1-1,N1-1] en/rho/vz/Bz=
   fl=[0,N2],[0,N1-1], mdot=fl, dq=[-1,N2][fl] vy =
   fl=[s2-1,N2-1],[0,N1-1], mdot=[s2-1,N2],[0,N1-1], dq=[-1,N2][fl] vx
   = fl=[0,N2-1][s1,N1-1], mdot=[0,N2],[s1-1,N1-1], dq=[-1,N2][fl]

   So when s1=0, need to at least bound mdot1, when s2=0, need to at
   least bound mdot2 due to dq[fl][-2,N1] needed for mdot1 and
   dq[-2,N2][fl] needed for mdot2


 */

void sweepx(void)
{
  // BEGIN variables
  FTYPE vpp, mdotp;
  static FTYPE(*dq)[N2M][N1M], (*dqv)[N3M][N2M][N1M],
      (*vp)[N3M][N2M][N1M];
  static FTYPE(*mdot)[N3M][N2M][N1M];
  static FTYPE(*p)[N3M][N2M][N1M];
  static FTYPE(*u)[N2M][N1M], (*Bzr)[N2M][N1M];
  static FTYPE(*fl)[N3M][N2M][N1M];


  register int i;
  int j, k;
  FTYPE ftemp;
  FTYPE ftemp0, ftemp1, ftemp2;
  FTYPE ftempv;
  FTYPE v1xa, v1ya, v1za;
  // END variables

  // pointer assignments
  p = workv1;
  vp = workv2;
  fl = workv3;
  mdot = workv4;
  dqv = workv5;
  dq = work6;
  u = work7;
  Bzr = work8;

  /* transform to momenta */
  LOOP {			// no need for LOOPV here since
    // multiplying and that's ok. 
    p[1][k][j][i] = z2e_1(s[1][k], j, i) * v[1][1][k][j][i];
    p[2][k][j][i] =
	z2e_2(s[1][k], j, i) * v[1][2][k][j][i] * g[2][2][i];
    p[3][k][j][i] =
	s[1][k][j][i] * v[1][3][k][j][i] * g[2][3][i] * g[2][4][j];
  }				// no need to bound

  LOOPH {			// only need half-full loop for sweep,
    // no
    // need to bound
    /* create spacial tranport variable */
    vp[1][k][j][i] = (v[1][1][k][j][i] - vg[1]) * dt;
  }



  /* first do mass: rho */
  dqx_calc(s[1], dq);

  // mdot is really mdot=Mdot*dt
  LOOPT0i {
    if (vp[1][k][j][i] > 0.) {
      mdot[1][k][j][i] =
	  (s[1][k][j][i - 1] +
	   (dx[2][1][i] - vp[1][k][j][i]) * dq[k][j][i -
						     1]) *
	  vp[1][k][j][i];
    } else {
      mdot[1][k][j][i] =
	  (s[1][k][j][i] -
	   (dx[2][1][i] +
	    vp[1][k][j][i]) * dq[k][j][i]) * vp[1][k][j][i];
    }
  }
  // only need to bound if doing periodic boundary conditions.  This
  // is
  // identically when you don't skip the first zone

  if ((skipix1 == 0) || (numprocs > 1)) {
    // bound(NULL,mdot,0,-3);
  }

  if (transiex1) {
    if (wgam) {
      /* then specific internal energy */
      LOOPH {			// only need half-full loop, no bound
	// needed with LOOPH
	u[k][j][i] = s[2][k][j][i] / s[1][k][j][i];
      }

      dqx_calc(u, dq);

      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1i {			// same as density w.r.t. loop/bound
	if (vp[1][k][j][i] > 0.)
	  fl[1][k][j][i] =
	      (u[k][j][i - 1] +
	       (dx[2][1][i] - vp[1][k][j][i]) * dq[k][j][i -
							 1]) *
	      mdot[1][k][j][i] * DS[1][1][k][j][i];
	else
	  fl[1][k][j][i] =
	      (u[k][j][i] -
	       (dx[2][1][i] +
		vp[1][k][j][i]) * dq[k][j][i]) *
	      mdot[1][k][j][i] * DS[1][1][k][j][i];
      }

      // assign internal energy density from surface fluxes
      LOOP {
	s[2][k][j][i] +=
	    (fl[1][k][j][i] - fl[1][k][j][i + 1]) * OVOL[1][k][j][i];
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
	  floorcnt[3][2]++;
	  if (s[2][k][j][i] < floorlowest[2]) {
	    floorlowest[2] = s[2][k][j][i];
	    wherelowest[2] = 3;
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
  } else {
    if (DOLOSSDIAG) {
      LOOPT1i {
	fl[1][k][j][i] = z2e_1(s[2][k], j, i) / z2e_1(s[1][k], j, i) * mdot[1][k][j][i] * DS[1][1][k][j][i];	// division 
														// 
	// 
	// by 
	// interp 
	// may 
	// be 
	// problem
      }
    }
  }


  if (DOLOSSDIAG) {
    // capture flux of enthalpy(only true for ideal gas EOS)
    if (reflectix1 == 0) {
      k = intix3;
      i = intix1;
      for (j = intix2; j < intox2; j++) {
	losss[2][1][0][j] += -gam * fl[1][k][j][i];
      }
    }
    if (reflectox1 == 0) {
      k = intix3;		// for fake 3d
      i = intox1;
      for (j = intix2; j < intox2; j++) {
	losss[2][1][1][j] += gam * fl[1][k][j][i];
      }
    }
  }





  if (transv1x1) {

    /* vx */
    dqvx_calc(DIRVEC, 1, v[1], dqv);

    LOOPT2i {			// only need to loop over -1 to N-1(1
      // bnd
      // zone) and should not bound fl!
      vpp = e2z_1(vp[1][k], j, i);
      mdotp = e2z_1(mdot[1][k], j, i);
      if (vpp > 0)
	fl[1][k][j][i] =
	    (v[1][1][k][j][i] +
	     (dx[1][1][i] -
	      vpp) * dqv[1][k][j][i]) * mdotp * DS[2][1][k][j][i];
      else
	fl[1][k][j][i] =
	    (v[1][1][k][j][i + 1] -
	     (dx[1][1][i] + vpp) * dqv[1][k][j][i +
						1]) * mdotp *
	    DS[2][1][k][j][i];
    }


    LOOPV1 {
      p[1][k][j][i] +=
	  (fl[1][k][j][i - 1] - fl[1][k][j][i]) * OVOL[2][k][j][i];
    }
  } else {
    if (DOLOSSDIAG) {
      LOOPT2i {
	fl[1][k][j][i] =
	    v[1][1][k][j][i] * mdot[1][k][j][i] * DS[1][1][k][j][i];
      }
    }
  }


  if (DOLOSSDIAG) {
    if (reflectix1 == 0) {
      k = intix3;
      i = intix1;
      for (j = intix2; j < intox2; j++) {
	ftemp = fl[1][k][j][i];	// not really on boundary!
	lossv[1][1][1][0][j] += -ftemp;
	if (kever == 0) {
	  lossv[1][0][1][0][j] += -ftemp * (v[1][1][k][j][i] * 0.5);
	}
      }
    }
    if (reflectox1 == 0) {
      k = intix3;
      i = intox1;
      for (j = intix2; j < intox2; j++) {
	ftemp = fl[1][k][j][i];	// not really on boundary!
	lossv[1][1][1][1][j] += ftemp;
	if (kever == 0) {
	  lossv[1][0][1][1][j] += ftemp * (v[1][1][k][j][i] * 0.5);
	}
      }
    }
  }



  if (transv2x1) {
    /* vy */
    dqvx_calc(!DIRVEC, 2, v[1], dqv);

    LOOPT3i {			// flux here lives on lower left corner
      mdotp = z2e_2(mdot[1][k], j, i);
      vpp = z2e_2(vp[1][k], j, i);

      if (vpp > 0.)
	fl[1][k][j][i] =
	    (v[1][2][k][j][i - 1] +
	     (dx[2][1][i] - vpp) * dqv[2][k][j][i -
						1]) * g[1][2][i] *
	    mdotp * DS[3][1][k][j][i];
      else
	fl[1][k][j][i] =
	    (v[1][2][k][j][i] -
	     (dx[2][1][i] +
	      vpp) * dqv[2][k][j][i]) * g[1][2][i] * mdotp *
	    DS[3][1][k][j][i];
    }

    LOOPV2 {
      p[2][k][j][i] +=
	  (fl[1][k][j][i] - fl[1][k][j][i + 1]) * OVOL[3][k][j][i];
    }
  } else {
    if (DOLOSSDIAG) {
      LOOPT3i {
	fl[1][k][j][i] =
	    e2e_v2(v[1][2][k], j,
		   i) * g[1][2][i] * mdot[1][k][j][i] *
	    DS[1][1][k][j][i];
      }
    }
  }


  if (DOLOSSDIAG) {
    if (reflectix1 == 0) {
      k = intix3;
      i = intix1;
      for (j = skipintix2; j < intox2; j++) {
	lossv[1][2][1][0][j] += -fl[1][k][j][i];
	// add in kinetic energy loss
	if (kever == 0) {
	  lossv[1][0][1][0][j] +=
	      -fl[1][k][j][i] * (z2e_1(v[1][2][k], j, i) * 0.5 /
				 g[1][2][i]);
	}
      }
    }
    if (reflectox1 == 0) {
      k = intix3;
      i = intox1;
      for (j = skipintix2; j < intox2; j++) {
	lossv[1][2][1][1][j] += fl[1][k][j][i];
	// add in kinetic energy loss
	if (kever == 0) {
	  lossv[1][0][1][1][j] +=
	      fl[1][k][j][i] * (z2e_1(v[1][2][k], j, i) * 0.5 /
				g[1][2][i]);
	}
      }
    }
  }


  if (transv3x1) {

    /* vz */
    dqvx_calc(!DIRVEC, 3, v[1], dqv);

    LOOPT1i {
      if (vp[1][k][j][i] > 0)
	fl[1][k][j][i] =
	    (v[1][3][k][j][i - 1] +
	     (dx[2][1][i] - vp[1][k][j][i]) * dqv[3][k][j][i -
							   1]) *
	    g[1][3][i] * g[2][4][j] * mdot[1][k][j][i] *
	    DS[1][1][k][j][i];
      else
	fl[1][k][j][i] =
	    (v[1][3][k][j][i] -
	     (dx[2][1][i] +
	      vp[1][k][j][i]) * dqv[3][k][j][i]) * g[1][3][i] *
	    g[2][4][j] * mdot[1][k][j][i] * DS[1][1][k][j][i];

    }

    LOOP {
      p[3][k][j][i] +=
	  (fl[1][k][j][i] - fl[1][k][j][i + 1]) * OVOL[1][k][j][i];
    }				// do not bound!
  } else {
    if (DOLOSSDIAG) {
      LOOPT1i {
	fl[1][k][j][i] =
	    z2e_1(v[1][3][k], j,
		  i) * g[1][3][i] * g[2][4][j] * mdot[1][k][j][i] *
	    DS[1][1][k][j][i];
      }
    }
  }


  if (DOLOSSDIAG) {
    // capture ang mom(s3) boundary loss
    if (reflectix1 == 0) {
      k = intix3;
      i = intix1;
      for (j = intix2; j < intox2; j++) {
	lossv[1][3][1][0][j] += -fl[1][k][j][i];
	// add in kinetic energy
	if (kever == 0) {
	  lossv[1][0][1][0][j] +=
	      -fl[1][k][j][i] * (z2e_1(v[1][3][k], j, i) * 0.5 /
				 (g[1][3][i] * g[2][4][j]));
	}
      }
    }
    if (reflectox1 == 0) {
      k = intix3;
      i = intox1;
      for (j = intix2; j < intox2; j++) {
	lossv[1][3][1][1][j] += fl[1][k][j][i];
	// add in kinetic energy
	if (kever == 0) {
	  lossv[1][0][1][1][j] +=
	      fl[1][k][j][i] * (z2e_1(v[1][3][k], j, i) * 0.5 /
				(g[1][3][i] * g[2][4][j]));
	}
      }
    }
  }



  if (transrhox1) {
    // fl is really fl=Flux*dt
    LOOPT1i {			// loops so gets only needed flux
      // values
      // for density calculation below, no BC
      // application should be applied
      fl[1][k][j][i] = mdot[1][k][j][i] * DS[1][1][k][j][i];
    }

    // correction and calculation of rho

    // compute local density from surface fluxes
    LOOP {
      s[1][k][j][i] +=
	  (fl[1][k][j][i] - fl[1][k][j][i + 1]) * OVOL[1][k][j][i];

#if(FORCERHO)
      if (s[1][k][j][i] < DENSITYFLOOR) {

	ftemp1 = (DENSITYFLOOR - s[1][k][j][i]);
	ftemp2 = ftemp1 * s[3][k][j][i];
	if ((i >= intix1) && (i < intox1) && (j >= intix2)
	    && (j < intox2))
	  floors[1] += ftemp1 * dvl[1][1][i] * dvl[1][2][j];
	if ((i >= intix1) && (i < intox1) && (j >= intix2)
	    && (j < intox2))
	  floors[3] += ftemp2 * dvl[1][1][i] * dvl[1][2][j];


	v1xa = e2z_1(v[1][1][k], j, i);
	v1ya = e2z_2(v[1][2][k], j, i);
	v1za = v[1][3][k][j][i];
	// floor ke located at zone center
	ftemp0 =
	    0.5 * (DENSITYFLOOR - s[1][k][j][i]) * (v1xa * v1xa +
						    v1ya * v1ya +
						    v1za * v1za);
	if ((i >= intix1) && (i < intox1) && (j >= intix2)
	    && (j < intox2))
	  floors[NUMSCA + 1] += ftemp0 * dvl[1][1][i] * dvl[1][2][j];

#if(FLOORDUMPFLAG==1)
	floorvars[1][k][j][i] += ftemp1;
	floorvars[3][k][j][i] += ftemp2;
	floorvar0[1][k][j][i] += ftemp0;
#endif
#if(DOFLOORDIAG==1)
	floorcnt[3][1]++;
	if (s[1][k][j][i] < floorlowest[1]) {
	  floorlowest[1] = s[1][k][j][i];
	  wherelowest[1] = 3;
	}
#endif
#if(DOFLOORD2==1)
	fprintf(logfl_file,
		"corrected rho in sweepx: t: %15.10g %d %d %d %15.10g\n",
		t, k, j, i, (DENSITYFLOOR - s[1][k][j][i]));
#endif
	s[1][k][j][i] = DENSITYFLOOR;
      }
#endif
#if(DOFLOORD2==1)
      if (s[1][k][j][i] < MIN) {
	fprintf(logfl_file, "small density: %d %g\n", i, s[1][k][j][i]);
      }
#endif

    }
    bound(NULL, NULL, 1, 0);
  } else {
    if (DOLOSSDIAG) {
      LOOPT1i {
	fl[1][k][j][i] = mdot[1][k][j][i] * DS[1][1][k][j][i];
      }
    }
  }



  if (DOLOSSDIAG) {
    // capture radial flux of mass
    if (reflectix1 == 0) {
      k = intix3;
      i = intix1;
      for (j = intix2; j < intox2; j++) {
	ftemp = fl[1][k][j][i];
	// mass
	losss[1][1][0][j] += -ftemp;
	// grav pot energy
	losss[3][1][0][j] += -ftemp * z2e_1(s[3][k], j, i);


	// enthalpy (computed above instead)
	// losss[2][1][0][j]+=-gam*z2e_1(s[2][k],j,i)/z2e_1(s[1][k],j,i)*ftemp;

	// adding up loss like this intead of as commented out
	// in
	// velocities gives better results in bondi runs

	if (kever == 1) {
	  // vx1-inner-ke
	  ftempv = v[1][1][k][j][i];
	  lossv[1][0][1][0][j] += -0.5 * ftemp * ftempv * ftempv;
	  // probably should do j=0 if skipix2=0
	  ftempv = e2e_v2(v[1][2][k], j, i);
	  lossv[1][0][1][0][j] += -0.5 * ftemp * ftempv * ftempv;	// vx2-inner-ke
	  ftempv = z2e_1(v[1][3][k], j, i);
	  lossv[1][0][1][0][j] += -0.5 * ftemp * ftempv * ftempv;	// vx3-inner-ke
	}
      }
    }
    if (reflectox1 == 0) {
      k = intix3;
      i = intox1;
      for (j = intix2; j < intox2; j++) {
	ftemp = fl[1][k][j][i];
	// mass
	losss[1][1][1][j] += ftemp;
	// grav pot energy
	losss[3][1][1][j] += ftemp * z2e_1(s[3][k], j, i);


	if (kever == 1) {
	  // vx1-outer-ke
	  ftempv = v[1][1][k][j][i];
	  lossv[1][0][1][1][j] += 0.5 * ftemp * ftempv * ftempv;
	  // probably should do j=0 if skipix2=0
	  ftempv = e2e_v2(v[1][2][k], j, i);
	  lossv[1][0][1][1][j] += 0.5 * ftemp * ftempv * ftempv;	// vx2-outer-ke
	  ftempv = z2e_1(v[1][3][k], j, i);
	  lossv[1][0][1][1][j] += 0.5 * ftemp * ftempv * ftempv;	// vx3-outer-ke
	}
      }
    }
  }


  /* return to velocities */

  if (transv1x1) {
    LOOPV1 {
      v[1][1][k][j][i] = p[1][k][j][i] / z2e_1(s[1][k], j, i);	// division 
								// 
      // 
      // by 
      // interp 
      // may 
      // be 
      // a 
      // problem

    }
  }
  if (transv2x1) {
    if (N2M > 1) {
      LOOPV2 {
	v[1][2][k][j][i] = p[2][k][j][i] / (z2e_2(s[1][k], j, i) * g[2][2][i]);	// division 
										// 
	// 
	// by 
	// interp 
	// may 
	// be 
	// a 
	// problem
      }
    }
  }
  if (transv3x1) {
    LOOP {
      v[1][3][k][j][i] =
	  p[3][k][j][i] / (s[1][k][j][i] * g[2][3][i] * g[2][4][j]);
    }
  }
  if ((transv1x1) || (transv2x1) || (transv3x1)) {
    bound(NULL, NULL, 0, 1);	// bound all grid velocity-vector
    // components
  }

}

void dqx_calc(FTYPE(*var)[N2M][N1M], FTYPE(*dq)[N2M][N1M])
{
  register int i, j, k;
  static FTYPE Dqp, Dqm, pr;

  if (advint == 1) {
    LOOPH {
      Dqp =
	  (var[k][j][i + 1] - var[k][j][i]) * OARCL[1][1][k][j][i + 1];
      Dqm = (var[k][j][i] - var[k][j][i - 1]) * OARCL[1][1][k][j][i];
      pr = Dqp * Dqm;
      dq[k][j][i] = (pr > 0.) ? pr / (Dqm + Dqp) : 0.;	// kill
      // factor
      // of 2
      // here
      // since q 
      // has 1/2
    }
  } else if (advint == 0) {
    LOOPH {
      dq[k][j][i] = 0;
    }
  }
}

void dqvx_calc(int wtype, int wcom, FTYPE(*var)[N3M][N2M][N1M],
	       FTYPE(*dqv)[N3M][N2M][N1M])
{
  register int i;
  int j, k;
  static FTYPE Dqp, Dqm, pr;

  // wtype=1 -> scalar type difference wtype=2 -> vector type
  // differencing

  if (advint == 1) {
    if (wtype != DIRVEC) {	// same as scalar above but do one comp 
				// 
      // of 
      // vector

      LOOPH {
	Dqp =
	    (var[wcom][k][j][i + 1] -
	     var[wcom][k][j][i]) * OARCL[3][1][k][j][i + 1];
	Dqm =
	    (var[wcom][k][j][i] -
	     var[wcom][k][j][i - 1]) * OARCL[3][1][k][j][i];
	pr = Dqp * Dqm;
	dqv[wcom][k][j][i] = (pr > 0.) ? pr / (Dqm + Dqp) : 0.;	// kill 
								// 
	// 
	// factor 
	// of 
	// 2 
	// here 
	// since 
	// q 
	// has 
	// 1/2
      }
    } else if (wtype == DIRVEC) {
      LOOPH {
	Dqp =
	    (var[wcom][k][j][i + 1] -
	     var[wcom][k][j][i]) * OARCL[2][1][k][j][i];
	Dqm =
	    (var[wcom][k][j][i] -
	     var[wcom][k][j][i - 1]) * OARCL[2][1][k][j][i - 1];
	pr = Dqp * Dqm;
	dqv[wcom][k][j][i] = (pr > 0.) ? pr / (Dqm + Dqp) : 0.;	// kill 
								// 
	// 
	// factor 
	// of 
	// 2 
	// here 
	// since 
	// q 
	// has 
	// 1/2
      }
    }
  } else if (advint == 0) {
    LOOPH {
      dqv[wcom][k][j][i] = 0;
    }
  }
}


void sweepy(void)
{
  FTYPE vpp, mdotp;
  static FTYPE(*dq)[N2M][N1M], (*dqv)[N3M][N2M][N1M],
      (*vp)[N3M][N2M][N1M];
  static FTYPE(*mdot)[N3M][N2M][N1M];
  static FTYPE(*p)[N3M][N2M][N1M];
  static FTYPE(*u)[N2M][N1M], (*Bzr)[N2M][N1M];
  static FTYPE(*fl)[N3M][N2M][N1M];
  register int i, j, k;

  FTYPE v1xa, v1ya, v1za, ftemp;
  FTYPE ftempv;
  FTYPE ftemp0, ftemp1, ftemp2;
  /* transform to momenta */

  p = workv1;
  vp = workv2;
  fl = workv3;
  mdot = workv4;
  dqv = workv5;
  dq = work6;
  u = work7;
  Bzr = work8;

  LOOP {			// no need for LOOPV here since
    // multiplying and that's ok. 
    p[1][k][j][i] = z2e_1(s[1][k], j, i) * v[1][1][k][j][i];
    p[2][k][j][i] =
	z2e_2(s[1][k], j, i) * v[1][2][k][j][i] * g[2][2][i];
    p[3][k][j][i] =
	s[1][k][j][i] * v[1][3][k][j][i] * g[2][3][i] * g[2][4][j];
  }				// no need to bound

  LOOPH {			// only need half full loop, no need to
    // bound
    /* create spacial tranport variable */
    vp[2][k][j][i] = (v[1][2][k][j][i] - vg[2]) * dt;
  }

  /* first do mass: rho */
  dqy_calc(s[1], dq);


  // mdot/fl is really fl=Flux*dt and mdot=Mdot*dt
  LOOPT0j {
    if (vp[2][k][j][i] > 0.) {
      mdot[2][k][j][i] =
	  (s[1][k][j - 1][i] +
	   (dx[2][2][j] * g[2][2][i] - vp[2][k][j][i]) * dq[k][j -
							       1][i])
	  * vp[2][k][j][i];
    } else {
      mdot[2][k][j][i] =
	  (s[1][k][j][i] -
	   (dx[2][2][j] * g[2][2][i] +
	    vp[2][k][j][i]) * dq[k][j][i]) * vp[2][k][j][i];
    }
  }
  if ((skipix2 == 0) || (numprocs > 1)) {
    bound(NULL, mdot, 0, -3);
  }


  if (transiex2) {
    if (wgam) {

      /* then specific internal energy */
      LOOPH {			// only need half-full loop
	u[k][j][i] = s[2][k][j][i] / s[1][k][j][i];
      }

      dqy_calc(u, dq);

      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1j {
	if (vp[2][k][j][i] > 0.)
	  fl[2][k][j][i] =
	      (u[k][j - 1][i] +
	       (dx[2][2][j] * g[2][2][i] -
		vp[2][k][j][i]) * dq[k][j -
					1][i]) *
	      mdot[2][k][j][i] * DS[1][2][k][j][i];
	else
	  fl[2][k][j][i] =
	      (u[k][j][i] -
	       (dx[2][2][j] * g[2][2][i] +
		vp[2][k][j][i]) * dq[k][j][i]) *
	      mdot[2][k][j][i] * DS[1][2][k][j][i];
      }

      LOOP {
	s[2][k][j][i] +=
	    (fl[2][k][j][i] - fl[2][k][j + 1][i]) * OVOL[1][k][j][i];

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
	  floorcnt[4][2]++;
	  if (s[2][k][j][i] < floorlowest[2]) {
	    floorlowest[2] = s[2][k][j][i];
	    wherelowest[2] = 4;
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
  } else {
    if (DOLOSSDIAG) {
      LOOPT1j {
	fl[2][k][j][i] = z2e_2(s[2][k], j, i) / z2e_2(s[1][k], j, i) * mdot[2][k][j][i] * DS[1][2][k][j][i];	// division 
														// 
	// 
	// by 
	// interp 
	// may 
	// be 
	// a 
	// problem
      }
    }
  }

  if (DOLOSSDIAG) {
    // capture flux of enthalpy
    if (reflectix2 == 0) {
      k = intix3;
      j = intix2;
      for (i = intix1; i < intox1; i++) {
	losss[2][2][0][i] += -gam * fl[2][k][j][i];
      }
    }
    if (reflectox2 == 0) {
      k = intix3;
      j = intox2;
      for (i = intix1; i < intox1; i++) {
	losss[2][2][1][i] += gam * fl[2][k][j][i];
      }
    }
  }

  // for sweepy, assume kinetic energy is added first in vy!
  if (transv2x2) {

    /* vy */
    dqvy_calc(DIRVEC, 2, v[1], dqv);

    LOOPT2j {
      vpp = e2z_2(vp[2][k], j, i);
      mdotp = e2z_2(mdot[2][k], j, i);
      if (vpp > 0)
	fl[2][k][j][i] =
	    (v[1][2][k][j][i] +
	     (dx[1][2][j] * g[2][2][i] -
	      vpp) * dqv[2][k][j][i]) * g[2][2][i] * mdotp *
	    DS[3][2][k][j][i];
      else
	fl[2][k][j][i] =
	    (v[1][2][k][j + 1][i] -
	     (dx[1][2][j] * g[2][2][i] + vpp) * dqv[2][k][j +
							  1][i]) *
	    g[2][2][i] * mdotp * DS[3][2][k][j][i];
    }

    LOOPV2 {
      p[2][k][j][i] +=
	  (fl[2][k][j - 1][i] - fl[2][k][j][i]) * OVOL[3][k][j][i];
    }
  } else {
    if (DOLOSSDIAG) {
      LOOPT2j {
	fl[2][k][j][i] =
	    v[1][2][k][j][i] * g[2][2][i] * mdot[2][k][j][i] *
	    DS[1][2][k][j][i];
      }
    }
  }

  if (DOLOSSDIAG) {
    if (reflectix2 == 0) {
      k = intix3;
      j = intix2;
      for (i = intix1; i < intox1; i++) {
	ftemp = fl[2][k][j][i];	// not on boundary!
	lossv[1][2][2][0][i] += -ftemp;
	if (kever == 0) {
	  lossv[1][0][2][0][i] +=
	      -ftemp * (v[1][2][k][j][i] * 0.5 / g[2][2][i]);
	}
      }
    }
    if (reflectox2 == 0) {
      k = intix3;
      j = intox2;
      for (i = intix1; i < intox1; i++) {
	ftemp = fl[2][k][j][i];
	lossv[1][2][2][1][i] += ftemp;
	if (kever == 0) {
	  lossv[1][0][2][1][i] +=
	      ftemp * (v[1][2][k][j][i] * 0.5 / g[2][2][i]);
	}
      }
    }
  }


  if (transv1x2) {

    /* vx */
    dqvy_calc(!DIRVEC, 1, v[1], dqv);

    LOOPT3j {
      // below avgs reason for loopt1j over other bzones
      mdotp = z2e_1(mdot[2][k], j, i);
      vpp = z2e_1(vp[2][k], j, i);
      if (vpp > 0.)
	fl[2][k][j][i] =
	    (v[1][1][k][j - 1][i] +
	     (dx[2][2][j] * g[1][2][i] - vpp) * dqv[1][k][j -
							  1][i]) *
	    mdotp * DS[2][2][k][j][i];
      else
	fl[2][k][j][i] =
	    (v[1][1][k][j][i] -
	     (dx[2][2][j] * g[1][2][i] +
	      vpp) * dqv[1][k][j][i]) * mdotp * DS[2][2][k][j][i];
    }

    LOOPV1 {
      p[1][k][j][i] +=
	  (fl[2][k][j][i] - fl[2][k][j + 1][i]) * OVOL[2][k][j][i];
    }
  } else {
    if (DOLOSSDIAG) {
      LOOPT3j {
	fl[2][k][j][i] =
	    e2e_v1(v[1][1][k], j,
		   i) * mdot[2][k][j][i] * DS[1][2][k][j][i];
      }
    }
  }

  if (DOLOSSDIAG) {
    if (reflectix2 == 0) {
      k = intix3;
      j = intix2;
      for (i = skipintix1; i < intox1; i++) {
	lossv[1][1][2][0][i] += -fl[2][k][j][i];
	// ke loss
	if (kever == 0) {
	  lossv[1][0][2][0][i] +=
	      -fl[2][k][j][i] * (z2e_2(v[1][1][k], j, i) * 0.5);
	}
      }
    }
    if (reflectox2 == 0) {
      k = intix3;
      j = intox2;
      for (i = skipintix1; i < intox1; i++) {
	lossv[1][1][2][1][i] += fl[2][k][j][i];
	// ke loss
	if (kever == 0) {
	  lossv[1][0][2][1][i] +=
	      fl[2][k][j][i] * (z2e_2(v[1][1][k], j, i) * 0.5);
	}
      }
    }
  }


  if (transv3x2) {

    /* vz */
    dqvy_calc(!DIRVEC, 3, v[1], dqv);

    LOOPT1j {
      if (vp[2][k][j][i] > 0)
	fl[2][k][j][i] =
	    (v[1][3][k][j - 1][i] +
	     (dx[2][2][j] * g[2][2][i] -
	      vp[2][k][j][i]) * dqv[3][k][j -
					  1][i]) * g[2][3][i] *
	    g[1][4][j] * mdot[2][k][j][i] * DS[1][2][k][j][i];
      else
	fl[2][k][j][i] =
	    (v[1][3][k][j][i] -
	     (dx[2][2][j] * g[2][2][i] +
	      vp[2][k][j][i]) * dqv[3][k][j][i]) * g[2][3][i] *
	    g[1][4][j] * mdot[2][k][j][i] * DS[1][2][k][j][i];
    }

    LOOP {
      p[3][k][j][i] +=
	  (fl[2][k][j][i] - fl[2][k][j + 1][i]) * OVOL[1][k][j][i];
    }				// do not bound!
  } else {
    if (DOLOSSDIAG) {
      LOOPT1j {
	fl[2][k][j][i] =
	    z2e_2(v[1][3][k], j,
		  i) * g[2][3][i] * g[1][4][j] * mdot[2][k][j][i] *
	    DS[1][2][k][j][i];
      }
    }
  }

  if (DOLOSSDIAG) {
    if (reflectix2 == 0) {
      k = intix3;
      j = intix2;
      for (i = intix1; i < intox1; i++) {
	lossv[1][3][2][0][i] += -fl[2][k][j][i];
	// ke loss
	if (kever == 0) {
	  lossv[1][0][2][0][i] +=
	      -fl[2][k][j][i] * (z2e_2(v[1][3][k], j, i) * 0.5 /
				 (g[2][3][i] * g[1][4][j]));
	}
      }
    }
    if (reflectox2 == 0) {
      k = intix3;
      j = intox2;
      for (i = intix1; i < intox1; i++) {
	lossv[1][3][2][1][i] += fl[2][k][j][i];
	// ke loss
	if (kever == 0) {
	  lossv[1][0][2][1][i] +=
	      fl[2][k][j][i] * (z2e_2(v[1][3][k], j, i) * 0.5 /
				(g[2][3][i] * g[1][4][j]));
	}
      }
    }
  }


  if (transrhox2) {

    // fl is really fl=Flux*dt
    LOOPT1j {
      fl[2][k][j][i] = mdot[2][k][j][i] * DS[1][2][k][j][i];
    }


    LOOP {
      s[1][k][j][i] +=
	  (fl[2][k][j][i] - fl[2][k][j + 1][i]) * OVOL[1][k][j][i];
#if(FORCERHO)
      if (s[1][k][j][i] < DENSITYFLOOR) {

	// floor

	ftemp1 = (DENSITYFLOOR - s[1][k][j][i]);
	ftemp2 = ftemp1 * s[3][k][j][i];
	// mass
	if ((i >= intix1) && (i < intox1) && (j >= intix2)
	    && (j < intox2))
	  floors[1] += ftemp1 * dvl[1][1][i] * dvl[1][2][j];
	// grav pot energy
	if ((i >= intix1) && (i < intox1) && (j >= intix2)
	    && (j < intox2))
	  floors[3] += ftemp2 * dvl[1][1][i] * dvl[1][2][j];

	v1xa = e2z_1(v[1][1][k], j, i);
	v1ya = e2z_2(v[1][2][k], j, i);
	v1za = v[1][3][k][j][i];
	// ke located at zone center
	ftemp0 =
	    0.5 * (DENSITYFLOOR - s[1][k][j][i]) * (v1xa * v1xa +
						    v1ya * v1ya +
						    v1za * v1za);
	if ((i >= intix1) && (i < intox1) && (j >= intix2)
	    && (j < intox2))
	  floors[NUMSCA + 1] += ftemp0 * dvl[1][1][i] * dvl[1][2][j];

#if(FLOORDUMPFLAG==1)
	floorvars[1][k][j][i] += ftemp1;
	floorvars[3][k][j][i] += ftemp2;
	floorvar0[1][k][j][i] += ftemp0;
#endif
#if(DOFLOORDIAG==1)
	floorcnt[4][1]++;
	if (s[1][k][j][i] < floorlowest[1]) {
	  floorlowest[1] = s[1][k][j][i];
	  wherelowest[1] = 4;
	}
#endif
#if(DOFLOORD2==1)
	fprintf(logfl_file,
		"corrected rho in sweepy: t: %15.10g %d %d %d %15.10g\n",
		t, k, j, i, (DENSITYFLOOR - s[1][k][j][i]));
#endif
	s[1][k][j][i] = DENSITYFLOOR;
      }
#endif

#if(DOFLOORD2==1)
      if (s[1][k][j][i] < MIN) {
	fprintf(logfl_file, "small density: %d %g\n", i, s[1][k][j][i]);
      }
#endif

    }
    bound(NULL, NULL, 1, 0);
  } else {
    if (DOLOSSDIAG) {
      LOOPT1j {
	fl[2][k][j][i] = mdot[2][k][j][i] * DS[1][2][k][j][i];
      }
    }
  }

  if (DOLOSSDIAG) {
    if (reflectix2 == 0) {
      k = intix3;
      j = intix2;
      for (i = intix1; i < intox1; i++) {
	ftemp = fl[2][k][j][i];
	// mass
	losss[1][2][0][i] += -ftemp;
	// losss[2][2][0][i]+=0; // don't use

	// grav pot energy
	losss[3][2][0][j] += -ftemp * z2e_2(s[3][k], j, i);

	if (kever == 1) {
	  // see sweepx for comment on why this done here
	  // instead of above commented ke's
	  ftempv = v[1][2][k][j][i];
	  lossv[1][0][2][0][i] += -0.5 * ftemp * ftempv * ftempv;	// vy
	  ftempv = e2e_v1(v[1][1][k], j, i);
	  lossv[1][0][2][0][i] += -0.5 * ftemp * ftempv * ftempv;	// vx
	  ftempv = z2e_2(v[1][3][k], j, i);
	  lossv[1][0][2][0][i] += -0.5 * ftemp * ftempv * ftempv;	// vz
	}
      }
    }
    if (reflectox2 == 0) {
      k = intix3;
      j = intox2;
      for (i = intix1; i < intox1; i++) {
	ftemp = fl[2][k][j][i];
	// mass
	losss[1][2][1][i] += ftemp;
	// losss[2][2][1][i]+=0; // don't use
	// grav pot energy
	losss[3][2][1][j] += ftemp * z2e_2(s[3][k], j, i);

	if (kever == 1) {
	  ftempv = v[1][2][k][j][i];
	  lossv[1][0][2][1][i] += 0.5 * ftemp * ftempv * ftempv;	// vy
	  ftempv = e2e_v1(v[1][1][k], j, i);
	  lossv[1][0][2][1][i] += 0.5 * ftemp * ftempv * ftempv;	// vx
	  ftempv = z2e_2(v[1][3][k], j, i);
	  lossv[1][0][2][1][i] += 0.5 * ftemp * ftempv * ftempv;	// vz
	}
      }
    }
  }

  /* return to velocities */

  if (transv1x2) {
    LOOPV1 {
      v[1][1][k][j][i] = p[1][k][j][i] / z2e_1(s[1][k], j, i);	// division 
								// 
      // 
      // by 
      // interp 
      // may 
      // be 
      // a 
      // problem
    }
  }
  if (transv2x2) {
    if (N2M > 1) {
      LOOPV2 {
	v[1][2][k][j][i] = p[2][k][j][i] / (z2e_2(s[1][k], j, i) * g[2][2][i]);	// division 
										// 
	// 
	// by 
	// interp 
	// may 
	// be 
	// a 
	// problem
      }
    }
  }
  if (transv3x2) {
    LOOP {
      v[1][3][k][j][i] =
	  p[3][k][j][i] / (s[1][k][j][i] * g[2][3][i] * g[2][4][j]);
    }
  }

  if ((transv1x2) || (transv2x2) || (transv3x2)) {
    bound(NULL, NULL, 0, 1);	// bound all grid velocity-vector
    // components
  }
}

void dqy_calc(FTYPE(*var)[N2M][N1M], FTYPE(*dq)[N2M][N1M])
{
  register int i, j, k;
  static FTYPE Dqp, Dqm, pr;

  if (advint == 1) {
    LOOPH {
      Dqp =
	  (var[k][j + 1][i] - var[k][j][i]) * OARCL[1][2][k][j + 1][i];
      Dqm = (var[k][j][i] - var[k][j - 1][i]) * OARCL[1][2][k][j][i];
      pr = Dqp * Dqm;
      dq[k][j][i] = (pr > 0.) ? pr / (Dqm + Dqp) : 0.;	// kill
      // factor
      // of 2
      // here
      // since q 
      // has 1/2
    }
  } else if (advint == 0) {
    LOOPH {
      dq[k][j][i] = 0;
    }
  }
}

// not sure if dq is right in spc, need grad?->curvature terms?
void dqvy_calc(int wtype, int wcom, FTYPE(*var)[N3M][N2M][N1M],
	       FTYPE(*dqv)[N3M][N2M][N1M])
{
  register int i, j, k;
  static FTYPE Dqp, Dqm, pr;

  if (advint == 1) {
    if (wtype != DIRVEC) {
      LOOPH {
	Dqp =
	    (var[wcom][k][j + 1][i] -
	     var[wcom][k][j][i]) * OARCL[2][2][k][j + 1][i];
	Dqm =
	    (var[wcom][k][j][i] -
	     var[wcom][k][j - 1][i]) * OARCL[2][2][k][j][i];
	pr = Dqp * Dqm;
	dqv[wcom][k][j][i] = (pr > 0.) ? pr / (Dqm + Dqp) : 0.;	// kill 
								// 
	// 
	// factor 
	// of 
	// 2 
	// here 
	// since 
	// q 
	// has 
	// 1/2
      }
    } else if (wtype == DIRVEC) {
      LOOPH {
	Dqp =
	    (var[wcom][k][j + 1][i] -
	     var[wcom][k][j][i]) * OARCL[3][2][k][j][i];
	Dqm =
	    (var[wcom][k][j][i] -
	     var[wcom][k][j - 1][i]) * OARCL[3][2][k][j - 1][i];
	pr = Dqp * Dqm;
	dqv[wcom][k][j][i] = (pr > 0.) ? pr / (Dqm + Dqp) : 0.;	// kill 
								// 
	// 
	// factor 
	// of 
	// 2 
	// here 
	// since 
	// q 
	// has 
	// 1/2
      }
    }
  } else if (advint == 0) {
    LOOPH {
      dqv[wcom][k][j][i] = 0;
    }
  }
}
