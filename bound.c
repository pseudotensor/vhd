/* enforce periodic boundary conditions */
#include "global.h"
#include "defs.h"

void legacyblock1(int bct, int bcdim, int bcdir, int k, int j, int i,
		  int *bcd1p, int *bcd2p, int *kkp, int *jjp, int *iip,
		  int *kkkp, int *jjjp, int *iiip);


void bound(FTYPE(*vars)[N2M][N1M],
	   FTYPE(*varv)[N3M][N2M][N1M], int wsca, int wvec)
{
  int jumpfactor;
  int god, godly = 0;
  int i, j, k, l, m;
  int ii, jj, kk, ll;
  int iii, jjj, kkk;
  int bct, bcdim, bcdir, bcd1, bcd2;
  int bcdir1 = 0, bcdir2 = 0, bcdir3 = 0;
  int bi, bj;
  FTYPE(*works)[N3M][N2M][N1M];
  FTYPE(*workv)[3][N3M][N2M][N1M];
  static int numhit = 0;
  int wbound = 0;
  int tempi = 0;
  int forces = 0, forcev = 0, forcecheck;

  // might not want to bound vectors if not advecting them
  // if(wvec==-1) return;
  // if(wvec==1) return;

  numhit++;
  if (wsca <= -2) {
    works = (FTYPE(*)[N3M][N2M][N1M]) (&vars[0][0][0]);
  } else
    works = s;
  if (wvec <= -2) {
    workv = (FTYPE(*)[3][N3M][N2M][N1M]) (&varv[0][0][0][0]);
  } else
    workv = v;

  // check to see if user wants to force outflow on some variable
  forces = forcev = 0;
  // force outflow for arbitrary incoming variable
  if (wsca == -2)
    forces++;			// force arbitrary to be outflow

  if (forces || forcev)
    forcecheck = 1;
  else
    forcecheck = 0;


#if(DOBOUNDSCA)
  // //////////// SCALARS
  if (wsca != 0) {
    for (l = 1; l <= NUMSCA; l++) {
      /* if not to do all, pick */
      if (wsca != -1) {
	if (wsca <= -2)
	  ll = 0;
	else
	  ll = wsca;
      } else
	ll = l;

      // determine how things are bounded for other vars
      if (ll == 0) {
	if (wsca == -2) {	// other
	  // bound like rho
	  wbound = 1;
	} else {
	  fprintf(fail_file,
		  "No definition for ll==0 for scalar case: wsca=%d\n",
		  wsca);
	  myexit(1);
	}
      } else
	wbound = ll;


      for (god = 1; god <= 5; god++) {
	if (god == 1)
	  godly = 3;
	else if (god == 2)
	  godly = 5;
	else if (god == 3)
	  godly = 4;
	else if (god == 4)
	  godly = 1;
	else if (god == 5)
	  godly = 2;



	LOOPBOUND {		// if using 2 boundary zones, still
	  // only
	  // loop over first layer, assuming 2nd
	  // layer is calculatable
#if( (!INNERBC)&&(N3BND!=0) )
	  if (((i >= 0) && (i < N1)) && ((j >= 0) && (j < N2))
	      && ((k >= 0) && (k < N3)))
	    i = N1;		// skip known non-boundary zones-3d
	  // version
#endif
#if( (!INNERBC)&&(N3BND==0) )
	  if (((i >= 0) && (i < N1)) && ((j >= 0) && (j < N2)))
	    i = N1;		// skip known non-boundary zones
#endif

	  bct = bcs[wbound][1][k][j][i];

	  if (forcecheck) {
	    if (bct == 3)
	      bct = 4;		// force dq or mdot to be
	    // outflowed(or other) if no
	    // assignments put in
	    else
	      bct = bcs[wbound][1][k][j][i];
	  }
	  if ((bct > 0) && (bct == godly)) {	/* skip zone if
						   computational zone
						   or null zone(corner
						   zones) */

	    bcdim = bcs[wbound][2][k][j][i];
	    bcdir = bcs[wbound][3][k][j][i];

	    if (bcdim == 1) {
	      bcdir1 = bcdir;
	      bcdir2 = 0;
	      bcdir3 = 0;
	    } else if (bcdim == 2) {
	      bcdir1 = 0;
	      bcdir2 = bcdir;
	      bcdir3 = 0;
	    } else if (bcdim == 3) {
	      bcdir1 = 0;
	      bcdir2 = 0;
	      bcdir3 = bcdir;
	    }
	    // below segment is legacy code for sanity
	    // checks
	    legacyblock1(bct, bcdim, bcdir, k, j, i, &bcd1,
			 &bcd2, &kk, &jj, &ii, &kkk, &jjj, &iii);

	    switch (bct) {

	    case 1:
	    case 2:
	      /* reflect/AOS same for scalar */
#if(LINEXT==0)
	      works[ll][k][j][i] = works[ll][kk][jj][ii];
#if(NBIGBND==2)
	      // copy 2nd layer of scalar
	      works[ll][k - bcdir3][j - bcdir2][i - bcdir1] =
		  works[ll][kkk][jjj][iii];
#endif
#else
	      if (bct == 4)
		tempi = 1;
	      else
		tempi = 0;
	      ex_s(bcdim, bcdir, works[ll][k], j, i, tempi);
#endif
	      break;
	    case 4:
	      /* outflow for scalar */
#if(LINEXT==0)
	      works[ll][k][j][i] = works[ll][kk][jj][ii];
#if(NBIGBND==2)
	      // copy 2nd layer of scalar
	      works[ll][k - bcdir3][j - bcdir2][i - bcdir1] =
		  works[ll][kk][jj][ii];
#endif
#else
#if(NBIGBND==2)
	      // fprintf(fail_file,"no 2bz for ex_s\n");
	      // myexit(1);
#endif
	      if (bct == 4)
		tempi = 1;
	      else
		tempi = 0;
	      ex_s(bcdim, bcdir, works[ll][k], j, i, tempi);
#endif
	      break;
	    case 3:
	      if ((wsca > 0) || (wsca == -1)) {
		works[ll][k][j][i] = sanal[ll][k][j][i];
#if(NBIGBND==2)
		works[ll][k - bcdir3][j - bcdir2][i -
						  bcdir1] =
		    sanal[ll][k - bcdir3][j - bcdir2][i - bcdir1];
#endif
	      } else if (wsca == -2) {
		fprintf(fail_file,
			"Case 3 bound with wsca==-2 has no definition\n");
		fprintf(fail_file, "%d %d %d %d %d\n",
			forcecheck, forces, forcev, numhit, bct);
		myexit(1);
	      }
	      break;
	    case 5:
	      /* periodic */
#if(LINEXT==0)
	      works[ll][k][j][i] = works[ll][kk][jj][ii];
#if(NBIGBND==2)
	      works[ll][k - bcdir3][j - bcdir2][i - bcdir1] =
		  works[ll][kk - bcdir3][jj - bcdir2][ii - bcdir1];
#endif
#else
#if(NBIGBND==2)
	      // fprintf(fail_file,"no 2bz for
	      // ex_s_p\n");
	      // myexit(1);
#endif
	      ex_s_p(bcdim, bcdir, works[ll][k], j, i);
#endif
	      break;
	    case 99:
	    case 98:
	      // do nothing
	      break;
	    default:
	      fprintf(fail_file,
		      "bound.c: switch(bct) error, case given is: %d\n\r",
		      bct);
	      myexit(1);
	    }			// end switch
	  }			// end if bct>0
	}			// end over current scalar boundary
	// zones(minus corners)
      }				// end over ordered conditions

#if(SCABCCORNER)
      // now do corner zones(assumes other boundary zones already
      // done!)

      // i=-1,j=-1 quadrant
      i = -1;
      j = -1;
      k = 0;

      bi = bcs[wbound][1][k][j + 1][i];	// i-boundary type
      bj = bcs[wbound][1][k][j][i + 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"s[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {

	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {

	  works[ll][k][j][i] = sanal[ll][k][j][i];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = sanal[ll][k][j][i - 1];
	  works[ll][k][j - 1][i] = sanal[ll][k][j - 1][i];
	  works[ll][k][j - 1][i - 1] = sanal[ll][k][j - 1][i - 1];
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {

	  works[ll][k][j][i] = works[ll][k][N2 - 1][N1 - 1];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = works[ll][k][N2 - 1][N1 - 2];
	  works[ll][k][j - 1][i] = works[ll][k][N2 - 2][N1 - 1];
	  works[ll][k][j - 1][i - 1] = works[ll][k][N2 - 2][N1 - 2];
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {

	  works[ll][k][j][i] = works[ll][k][j + 1][i];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = works[ll][k][j + 1][i - 1];
	  works[ll][k][j - 1][i] = works[ll][k][j + 1][i];
	  works[ll][k][j - 1][i - 1] = works[ll][k][j + 1][i - 1];
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {

	  works[ll][k][j][i] = works[ll][k][j][i + 1];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = works[ll][k][j][i + 1];
	  works[ll][k][j - 1][i] = works[ll][k][j - 1][i + 1];
	  works[ll][k][j - 1][i - 1] = works[ll][k][j - 1][i + 1];
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {

	  works[ll][k][j][i] = works[ll][k][j + 1][i];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i - 1] = works[ll][k][j + 1][i - 1];
	  works[ll][k][j - 1][i] = works[ll][k][j + jumpfactor][i];
	  works[ll][k][j - 1][i - 1] =
	      works[ll][k][j + jumpfactor][i - 1];
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {

	  works[ll][k][j][i] = works[ll][k][j][i + 1];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i - 1] = works[ll][k][j][i + jumpfactor];
	  works[ll][k][j - 1][i] = works[ll][k][j - 1][i + 1];
	  works[ll][k][j - 1][i - 1] =
	      works[ll][k][j - 1][i + jumpfactor];
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "s[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
      // i=N1,j=-1 quadrant
      i = N1;
      j = -1;
      k = 0;

      bi = bcs[wbound][1][k][j + 1][i];	// i-boundary type
      bj = bcs[wbound][1][k][j][i - 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"s[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {

	  works[ll][k][j][i] = sanal[ll][k][j][i];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = sanal[ll][k][j][i + 1];
	  works[ll][k][j - 1][i] = sanal[ll][k][j - 1][i];
	  works[ll][k][j - 1][i + 1] = sanal[ll][k][j - 1][i + 1];
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {

	  works[ll][k][j][i] = works[ll][k][N2 - 1][0];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = works[ll][k][N2 - 1][1];
	  works[ll][k][j - 1][i] = works[ll][k][N2 - 2][0];
	  works[ll][k][j - 1][i + 1] = works[ll][k][N2 - 2][1];
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {

	  works[ll][k][j][i] = works[ll][k][j + 1][i];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = works[ll][k][j + 1][i + 1];
	  works[ll][k][j - 1][i] = works[ll][k][j + 1][i];
	  works[ll][k][j - 1][i + 1] = works[ll][k][j + 1][i + 1];
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {

	  works[ll][k][j][i] = works[ll][k][j][i - 1];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = works[ll][k][j][i - 1];
	  works[ll][k][j - 1][i] = works[ll][k][j - 1][i - 1];
	  works[ll][k][j - 1][i + 1] = works[ll][k][j - 1][i - 1];
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {

	  works[ll][k][j][i] = works[ll][k][j + 1][i];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i + 1] = works[ll][k][j + 1][i + 1];
	  works[ll][k][j - 1][i] = works[ll][k][j + jumpfactor][i];
	  works[ll][k][j - 1][i + 1] =
	      works[ll][k][j + jumpfactor][i + 1];
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {

	  works[ll][k][j][i] = works[ll][k][j][i - 1];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i + 1] = works[ll][k][j][i - jumpfactor];
	  works[ll][k][j - 1][i] = works[ll][k][j - 1][i - 1];
	  works[ll][k][j - 1][i + 1] =
	      works[ll][k][j - 1][i - jumpfactor];
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "s[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
      // i=N1,j=N2 quadrant
      i = N1;
      j = N2;
      k = 0;

      bi = bcs[wbound][1][k][j - 1][i];	// i-boundary type
      bj = bcs[wbound][1][k][j][i - 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"s[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {

	  works[ll][k][j][i] = sanal[ll][k][j][i];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = sanal[ll][k][j][i + 1];
	  works[ll][k][j + 1][i] = sanal[ll][k][j + 1][i];
	  works[ll][k][j + 1][i + 1] = sanal[ll][k][j + 1][i + 1];
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {

	  works[ll][k][j][i] = works[ll][k][0][0];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = works[ll][k][0][1];
	  works[ll][k][j + 1][i] = works[ll][k][1][0];
	  works[ll][k][j + 1][i + 1] = works[ll][k][1][1];
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {

	  works[ll][k][j][i] = works[ll][k][j - 1][i];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = works[ll][k][j - 1][i + 1];
	  works[ll][k][j + 1][i] = works[ll][k][j - 1][i];
	  works[ll][k][j + 1][i + 1] = works[ll][k][j - 1][i + 1];
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {

	  works[ll][k][j][i] = works[ll][k][j][i - 1];
#if(NBIGBND==2)
	  works[ll][k][j][i + 1] = works[ll][k][j][i - 1];
	  works[ll][k][j + 1][i] = works[ll][k][j + 1][i - 1];
	  works[ll][k][j + 1][i + 1] = works[ll][k][j + 1][i - 1];
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {

	  works[ll][k][j][i] = works[ll][k][j - 1][i];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i + 1] = works[ll][k][j - 1][i + 1];
	  works[ll][k][j + 1][i] = works[ll][k][j - jumpfactor][i];
	  works[ll][k][j + 1][i + 1] =
	      works[ll][k][j - jumpfactor][i + 1];
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {

	  works[ll][k][j][i] = works[ll][k][j][i - 1];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i + 1] = works[ll][k][j][i - jumpfactor];
	  works[ll][k][j + 1][i] = works[ll][k][j + 1][i - 1];
	  works[ll][k][j + 1][i + 1] =
	      works[ll][k][j + 1][i - jumpfactor];
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "s[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }

      // i=-1,j=N2 quadrant
      i = -1;
      j = N2;
      k = 0;

      bi = bcs[wbound][1][k][j - 1][i];	// i-boundary type
      bj = bcs[wbound][1][k][j][i + 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"s[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {

	  works[ll][k][j][i] = sanal[ll][k][j][i];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = sanal[ll][k][j][i - 1];
	  works[ll][k][j + 1][i] = sanal[ll][k][j + 1][i];
	  works[ll][k][j + 1][i - 1] = sanal[ll][k][j + 1][i - 1];
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {

	  works[ll][k][j][i] = works[ll][k][0][N1 - 1];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = works[ll][k][0][N1 - 2];
	  works[ll][k][j + 1][i] = works[ll][k][1][N1 - 1];
	  works[ll][k][j + 1][i - 1] = works[ll][k][1][N1 - 2];
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {

	  works[ll][k][j][i] = works[ll][k][j - 1][i];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = works[ll][k][j - 1][i - 1];
	  works[ll][k][j + 1][i] = works[ll][k][j - 1][i];
	  works[ll][k][j + 1][i - 1] = works[ll][k][j - 1][i - 1];
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {

	  works[ll][k][j][i] = works[ll][k][j][i + 1];
#if(NBIGBND==2)
	  works[ll][k][j][i - 1] = works[ll][k][j][i + 1];
	  works[ll][k][j + 1][i] = works[ll][k][j + 1][i + 1];
	  works[ll][k][j + 1][i - 1] = works[ll][k][j + 1][i + 1];
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {

	  works[ll][k][j][i] = works[ll][k][j - 1][i];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i - 1] = works[ll][k][j - 1][i - 1];
	  works[ll][k][j + 1][i] = works[ll][k][j - jumpfactor][i];
	  works[ll][k][j + 1][i - 1] =
	      works[ll][k][j - jumpfactor][i - 1];
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {

	  works[ll][k][j][i] = works[ll][k][j][i + 1];
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
	  works[ll][k][j][i - 1] = works[ll][k][j][i + jumpfactor];
	  works[ll][k][j + 1][i] = works[ll][k][j + 1][i + 1];
	  works[ll][k][j + 1][i - 1] =
	      works[ll][k][j + 1][i + jumpfactor];
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "s[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
#endif



      /* cut short loop if only to do one */
      if (wsca != -1)
	l = NUMSCA;
    }				// end over scalars
  }				// end if any scalars

#endif








#if(DOBOUNDVEC)
  // //////////// VECTORS



  // Now do vectors if any
  if (wvec != 0) {
    for (l = 1; l <= NUMVEC - 1; l++) {
      /* if not to do all, pick */
      if (wvec != -1) {
	if (wvec <= -2)
	  ll = 0;
	else
	  ll = wvec;
      } else
	ll = l;

      if (ll == 0) {
	if (wvec == -3) {	// mdot
	  // bound like v
	  wbound = 1;
	} else {
	  fprintf(fail_file,
		  "No definition for ll==0, wvec=%d\n", wvec);
	  myexit(1);
	}
      } else
	wbound = ll;





      for (god = 1; god <= 5; god++) {
	if (god == 1)
	  godly = 3;
	else if (god == 2)
	  godly = 5;
	else if (god == 3)
	  godly = 4;
	else if (god == 4)
	  godly = 1;
	else if (god == 5)
	  godly = 2;

	LOOPBOUND {
#if( (!INNERBC)&&(N3BND!=0) )
	  if (((i >= 0) && (i < N1)) && ((j >= 0) && (j < N2))
	      && ((k >= 0) && (k < N3)))
	    i = N1;		// skip known non-boundary zones-3d
	  // version
#endif
#if( (!INNERBC)&&(N3BND==0) )
	  if (((i >= 0) && (i < N1)) && ((j >= 0) && (j < N2)))
	    i = N1;		// skip known non-boundary zones
#endif

	  bct = bcv[wbound][1][k][j][i];

	  if (forcecheck) {
	    if (bct == 3)
	      bct = 4;		// force dq or mdot to be
	    // outflowed(or other) if no
	    // assignments put in
	    else
	      bct = bcv[wbound][1][k][j][i];
	  }

	  if ((bct > 0) && (bct == godly)) {	/* skip zone if
						   computational zone
						   or null zone(corner
						   zones) */

	    bcdim = bcv[wbound][2][k][j][i];
	    bcdir = bcv[wbound][3][k][j][i];


	    if (bcdim == 1) {
	      bcdir1 = bcdir;
	      bcdir2 = 0;
	      bcdir3 = 0;
	    } else if (bcdim == 2) {
	      bcdir1 = 0;
	      bcdir2 = bcdir;
	      bcdir3 = 0;
	    } else if (bcdim == 3) {
	      bcdir1 = 0;
	      bcdir2 = 0;
	      bcdir3 = bcdir;
	    }
	    // below segment is legacy code for sanity
	    // checks
	    legacyblock1(bct, bcdim, bcdir, k, j, i, &bcd1,
			 &bcd2, &kk, &jj, &ii, &kkk, &jjj, &iii);

	    switch (bct) {

	    case 1:
	    case 2:
	      /* reflect and AOS */
	      if (bct == 1)
		tempi = 1;
	      else if (bct == 2)
		tempi = -1;
#if(LINEXT==0)
	      if (bcdim == 1) {
		workv[ll][1][k][j][i + bcd1] = 0;	/* reflect
							   inner
							   vector */
		if (i + bcd2 < N1 + N1BND)
		  workv[ll][1][k][j][i + bcd2] = -workv[ll][1][kk][jj][ii + bcd1];	/* reflect 
											   outer 
											   vector 
		 */// gets "2nd" velocity layer on right
		// of
		// grid

		workv[ll][2][k][j][i] = workv[ll][2][kk][jj][ii];
		workv[ll][2][k][j + 1][i] =
		    workv[ll][2][kk][jj + 1][ii];
		workv[ll][3][k][j][i] = tempi * workv[ll][3][kk][jj][ii];	/* same 
										   for 
										   DIM=2 
										   or 
										   3 
										 */
#if(NBIGBND==2)
		// correct for singular case
		if (N1 == 1)
		  jumpfactor = 0;
		else
		  jumpfactor = 1;
		if (i + bcd2 - bcdir1 < N1 + N1BND)
		  workv[ll][1][k][j][i + bcd2 - bcdir1] =
		      -workv[ll][1][kk][jj][ii + bcd1 +
					    bcdir1 * jumpfactor];
		workv[ll][2][k][j][i - bcdir1] =
		    workv[ll][2][kk][jj][ii + bcdir1 * jumpfactor];
		workv[ll][2][k][j + 1][i - bcdir1] =
		    workv[ll][2][kk][jj + 1][ii + bcdir1 * jumpfactor];
		workv[ll][3][k][j][i - bcdir1] = tempi * workv[ll][3][kk][jj][ii + bcdir1 * jumpfactor];	/* same 
														   for 
														   DIM=2 
														   or 
														   3 
														 */
#endif
	      }
	      if (bcdim == 2) {
		workv[ll][2][k][j + bcd1][i] = 0;	/* reflect
							   inner
							   vector */
		if (j + bcd2 < N2 + N2BND)
		  workv[ll][2][k][j + bcd2][i] = -workv[ll][2][kk][jj + bcd1][ii];	/* reflect 
											   outer 
											   vector 
											 */
		workv[ll][1][k][j][i] = workv[ll][1][kk][jj][ii];
		workv[ll][1][k][j][i + 1] =
		    workv[ll][1][kk][jj][ii + 1];
		workv[ll][3][k][j][i] = tempi * workv[ll][3][kk][jj][ii];	/* same 
										   for 
										   DIM=2 
										   or 
										   3 
										 */
#if(NBIGBND==2)
		// correct for singular case
		if (N2 == 1)
		  jumpfactor = 0;
		else
		  jumpfactor = 1;

		if (j + bcd2 - bcdir2 < N2 + N2BND)
		  workv[ll][2][k][j + bcd2 - bcdir2][i] = -workv[ll][2][kk][jj + bcd1 + bcdir2 * jumpfactor][ii];	// reflect 
															// 
		// 
		// outer 
		// vector
		workv[ll][1][k][j - bcdir2][i] =
		    workv[ll][1][kk][jj + bcdir2 * jumpfactor][ii];
		workv[ll][1][k][j - bcdir2][i + 1] =
		    workv[ll][1][kk][jj + bcdir2 * jumpfactor][ii + 1];
		workv[ll][3][k][j - bcdir2][i] = tempi * workv[ll][3][kk][jj + bcdir2 * jumpfactor][ii];	// same 
														// 
		// 
		// for 
		// DIM=2 
		// or 
		// 3
#endif
	      }
	      if (bcdim == 3) {	/* never reached unless really in 3-d */
		workv[ll][3][k + bcd1][j][i] = 0;	/* reflect
							   inner
							   vector */
		if (k + bcd2 < N3 + N3BND)
		  workv[ll][3][k + bcd2][j][i] = -workv[ll][3][kk + bcd1][jj][ii];	/* reflect 
											   outer 
											   vector 
											 */
		workv[ll][1][k][j][i] =
		    tempi * workv[ll][1][kk][jj][ii];
		workv[ll][1][k][j][i + 1] =
		    tempi * workv[ll][1][kk][jj][ii + 1];
		workv[ll][2][k][j][i] =
		    tempi * workv[ll][2][kk][jj][ii];
		workv[ll][2][k][j + 1][i] =
		    tempi * workv[ll][2][kk][jj + 1][ii];
#if(NBIGBND==2)
		// correct for singular case
		if (N3 == 1)
		  jumpfactor = 0;
		else
		  jumpfactor = 1;

		if (k + bcd2 - bcdir3 < N3 + N3BND)
		  workv[ll][3][k + bcd2 - bcdir3][j][i] = -workv[ll][3][kk + bcd1 + bcdir3 * jumpfactor][jj][ii];	/* reflect 
															   outer 
															   vector 
															 */
		workv[ll][1][k - bcdir3][j][i] =
		    tempi * workv[ll][1][kk + bcdir3 * jumpfactor][jj]
		    [ii];
		workv[ll][1][k - bcdir3][j][i + 1] =
		    tempi * workv[ll][1][kk +
					 bcdir3 *
					 jumpfactor][jj][ii + 1];
		workv[ll][2][k - bcdir3][j][i] =
		    tempi * workv[ll][2][kk + bcdir3 * jumpfactor][jj]
		    [ii];
		workv[ll][2][k - bcdir3][j + 1][i] =
		    tempi * workv[ll][2][kk +
					 bcdir3 * jumpfactor][jj + 1]
		    [ii];
#endif
	      }
#else
#if(NBIGBND==2)
	      // fprintf(fail_file,"no 2bz for ex_v\n");
	      // myexit(1);
#endif
	      ex_v(bcdim, bcdir, workv[ll][bcdim][k], j, i, 1);
	      ex_s(bcdim, bcdir,
		   workv[ll][3 - (4 - bcdim) % 3][k], j, i, 1);
	      ex_s(bcdim, bcdir, workv[ll][bcdim % 3 + 1][k], j, i, 1);
#endif
	      break;
	    case 3:
	      /* Fix/Time vary: Dirichlet */
	      if ((wvec > 0) || (wvec == -1)) {
		for (m = 1; m <= 3; m++) {
		  workv[ll][m][k][j][i] = vanal[ll][m][k][j][i];	// v[-1] 
									// 
		  // 
		  // or 
		  // v[N]
#if(NBIGBND==2)
		  workv[ll][m][k - bcdir3][j - bcdir2][i - bcdir1] = vanal[ll][m][k - bcdir3][j - bcdir2][i - bcdir1];	// v[-2] 
															// 
		  // 
		  // or 
		  // v[N+1]
#endif
		  // if on inner edge, also get
		  // element
		  // 0(due to staggered grid)
		  if (m == bcdim) {
		    if (bcdir == 1) {
		      workv[ll][m][k + bcdir3][j + bcdir2][i + bcdir1] = vanal[ll][m][k + bcdir3][j + bcdir2][i + bcdir1];	// v[0] 
																// 
		      // 
		      // 
		    }
		  }
		}
	      } else if (wvec == -2) {
		fprintf(fail_file,
			"Case 3 bound with wvec==-2 has no definition\n");
		myexit(1);
	      }
	      break;
	    case 4:
	      if (wvec == -3)
		break;		// this outflow is causing
	      // problems with wvec==-3!!
	      /* outflow */
#if(LINEXT==0)
	      // deal with asymmetry in velocity
	      // components
	      // on grid w.r.t. inner/outer edges
	      // deals also with inflow checking
	      if (bcdim == 1) {
		if (bcdir == 1) {
		  if (INFLOWCHECKIX1 && (wbound == 1)
		      && (workv[ll][1][k][j][ii + 1] > 0.0)) {
		    workv[ll][1][k][j][ii] = 0.0;
		    workv[ll][1][k][j][i] = 0.0;
#if(NBIGBND==2)
		    workv[ll][1][k][j][i - bcdir1] = 0.0;
#endif
		  } else {
		    workv[ll][1][k][j][ii] = workv[ll][1][k][j][ii + 1];
		    workv[ll][1][k][j][i] = workv[ll][1][k][j][ii];
#if(NBIGBND==2)
		    workv[ll][1][k][j][i - bcdir1] =
			workv[ll][1][k][j][i];
#endif
		  }
		} else {
		  if (INFLOWCHECKOX1 && (wbound == 1)
		      && (workv[ll][1][k][j][ii] < 0.0)) {
		    workv[ll][1][k][j][i] = 0.0;
#if(NBIGBND==2)
		    workv[ll][1][k][j][i + 1] = 0.0;
#endif
		  } else {
		    workv[ll][1][k][j][i] = workv[ll][1][k][j][ii];
#if(NBIGBND==2)
		    workv[ll][1][k][j][i + 1] = workv[ll][1][k][j][i];
#endif
		  }
		}
		workv[ll][2][k][j][i] = workv[ll][2][k][j][ii];
		workv[ll][2][k][j + 1][i] = workv[ll][2][k][j + 1][ii];
		workv[ll][3][k][j][i] = workv[ll][3][k][j][ii];
#if(NBIGBND==2)
		workv[ll][2][k][j][i - bcdir1] = workv[ll][2][k][j][i];
		workv[ll][2][k][j + 1][i - bcdir1] =
		    workv[ll][2][k][j + 1][i];
		workv[ll][3][k][j][i - bcdir1] = workv[ll][3][k][j][i];
#endif
	      } else if (bcdim == 2) {
		if (bcdir == 1) {
		  if (INFLOWCHECKIX2 && (wbound == 1)
		      && (workv[ll][2][k][jj + 1][i] > 0.0)) {
		    workv[ll][2][k][jj][i] = 0.0;
		    workv[ll][2][k][j][i] = 0.0;
#if(NBIGBND==2)
		    workv[ll][2][k][j - bcdir2][i] = 0.0;
#endif
		  } else {
		    workv[ll][2][k][jj][i] = workv[ll][2][k][jj + 1][i];
		    workv[ll][2][k][j][i] = workv[ll][2][k][jj][i];
#if(NBIGBND==2)
		    workv[ll][2][k][j - bcdir2][i] =
			workv[ll][2][k][j][i];
#endif
		  }
		} else {
		  if (INFLOWCHECKOX2 && (wbound == 1)
		      && (workv[ll][2][k][jj][i] < 0.0)) {
		    workv[ll][2][k][j][i] = 0.0;
#if(NBIGBND==2)
		    workv[ll][2][k][j + 1][i] = 0.0;
#endif
		  } else {
		    workv[ll][2][k][j][i] = workv[ll][2][k][jj][i];
#if(NBIGBND==2)
		    workv[ll][2][k][j + 1][i] = workv[ll][2][k][j][i];
#endif
		  }
		}
		workv[ll][1][k][j][i] = workv[ll][1][k][jj][i];
		workv[ll][1][k][j][i + 1] = workv[ll][1][k][jj][i + 1];
		workv[ll][3][k][j][i] = workv[ll][3][k][jj][i];
#if(NBIGBND==2)
		workv[ll][1][k][j - bcdir2][i] = workv[ll][1][k][j][i];
		workv[ll][1][k][j - bcdir2][i + 1] =
		    workv[ll][1][k][j][i + 1];
		workv[ll][3][k][j - bcdir2][i] = workv[ll][3][k][j][i];
#endif
	      } else if (bcdim == 3) {
		if (bcdir == 1) {
		  if (INFLOWCHECKIX3 && (wbound == 1)
		      && (workv[ll][3][kk + 1][j][i] > 0.0)) {
		    workv[ll][3][kk][j][i] = 0.0;
		    workv[ll][3][k][j][i] = 0.0;
#if(NBIGBND==2)
		    workv[ll][3][k - bcdir3][j][i] = 0.0;
#endif
		  } else {
		    workv[ll][3][kk][j][i] = workv[ll][3][kk + 1][j][i];
		    workv[ll][3][k][j][i] = workv[ll][3][kk][j][i];
#if(NBIGBND==2)
		    workv[ll][3][k - bcdir3][j][i] =
			workv[ll][3][k][j][i];
#endif
		  }
		} else {
		  if (INFLOWCHECKOX3 && (wbound == 1)
		      && (workv[ll][3][kk][j][i] < 0.0)) {
		    workv[ll][3][k][j][i] = 0.0;
#if(NBIGBND==2)
		    workv[ll][3][k + 1][j][i] = 0.0;
#endif
		  } else {
		    workv[ll][3][k][j][i] = workv[ll][3][kk][j][i];
#if(NBIGBND==2)
		    workv[ll][3][k + 1][j][i] = workv[ll][3][k][j][i];
#endif
		  }
		}
		workv[ll][1][k][j][i] = workv[ll][1][kk][j][i];
		workv[ll][1][k][j][i + 1] = workv[ll][1][kk][j][i + 1];
		workv[ll][2][k][j][i] = workv[ll][2][kk][j][i];
		workv[ll][2][k][j + 1][i] = workv[ll][2][kk][j + 1][i];
#if(NBIGBND==2)
		workv[ll][1][k - bcdir3][j][i] = workv[ll][1][k][j][i];
		workv[ll][1][k - bcdir3][j][i + 1] =
		    workv[ll][1][k][j][i + 1];
		workv[ll][2][k - bcdir3][j][i] = workv[ll][2][k][j][i];
		workv[ll][2][k - bcdir3][j + 1][i] =
		    workv[ll][2][k][j + 1][i];
#endif
	      }
#else
#if(NBIGBND==2)
	      // fprintf(fail_file,"no 2bz for ex_v\n");
	      // myexit(1);
#endif
	      for (m = 1; m <= 3; m++) {
		ex_v(bcdim, bcdir, workv[ll][m][k], j, i, 0);
	      }
#endif
	      break;
	    case 5:
	      /* periodic */
#if(LINEXT==0)
	      for (m = 1; m <= 3; m++) {
		workv[ll][m][k][j][i] = workv[ll][m][kk][jj][ii];
#if(NBIGBND==2)
		workv[ll][m][k - bcdir3][j - bcdir2][i - bcdir1]
		    =
		    workv[ll][m][kk - bcdir3][jj - bcdir3][ii - bcdir3];
#endif
	      }
#else
#if(NBIGBND==2)
	      // fprintf(fail_file,"no 2bz for ex_v\n");
	      // myexit(1);
#endif
	      ex_v_p(bcdim, bcdir, workv[ll][bcdim][k], j, i);
	      ex_s_p(bcdim, bcdir,
		     workv[ll][3 - (4 - bcdim) % 3][k], j, i);
	      ex_s_p(bcdim, bcdir, workv[ll][bcdim % 3 + 1][k], j, i);
#endif
	      break;
	    case 99:
	    case 98:
	      // no nothing
	      break;
	    default:
	      fprintf(fail_file,
		      "bound.c: switch(bct) error, case given is: %d\n\r",
		      bct);
	      myexit(1);
	    }			// end switch
	  }			// end if bct>0 for vectors
	}			// end over current vector
      }				// end over ordered condition

#if(VECBCCORNER)

      // now do corner zones(assumes other boundary zones already
      // done!)
      // normal boundary code takes care of 0 zone issue(due to
      // staggered grid)

      // i=-1,j=-1 quadrant
      i = -1;
      j = -1;
      k = 0;

      bi = bcv[wbound][1][k][j + 1][i];	// i-boundary type
      bj = bcv[wbound][1][k][j][i + 1];	// j-boundary type

      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"v[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {
	  for (m = 1; m <= 3; m++) {

	    workv[ll][m][k][j][i] = vanal[ll][m][k][j][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = vanal[ll][m][k][j][i - 1];
	    workv[ll][m][k][j][i - 2] = vanal[ll][m][k][j][i - 2];
	    workv[ll][m][k][j - 1][i] = vanal[ll][m][k][j - 1][i];
	    workv[ll][m][k][j - 1][i - 1] =
		vanal[ll][m][k][j - 1][i - 1];
	  }
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {

	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][N2 - 1][N1 - 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][N2 - 1][N1 - 2];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][N2 - 2][N1 - 1];
	    workv[ll][m][k][j - 1][i - 1] =
		workv[ll][m][k][N2 - 2][N1 - 2];
	  }
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j + 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][j + 1][i - 1];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][j + 1][i];
	    workv[ll][m][k][j - 1][i - 1] =
		workv[ll][m][k][j + 1][i - 1];
	  }
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {
	  for (m = 1; m <= 3; m++) {

	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i + 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][j][i + 1];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][j - 1][i + 1];
	    workv[ll][m][k][j - 1][i - 1] =
		workv[ll][m][k][j - 1][i + 1];
	  }
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j + 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][j + 1][i - 1];
	    workv[ll][m][k][j - 1][i] =
		workv[ll][m][k][j + jumpfactor][i];
	    workv[ll][m][k][j - 1][i - 1] =
		workv[ll][m][k][j + jumpfactor][i - 1];
	  }
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i + 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] =
		workv[ll][m][k][j][i + jumpfactor];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][j - 1][i + 1];
	    workv[ll][m][k][j - 1][i - 1] =
		workv[ll][m][k][j - 1][i + jumpfactor];
	  }
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "v[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
      // i=N1,j=-1 quadrant
      i = N1;
      j = -1;
      k = 0;

      bi = bcv[wbound][1][k][j + 1][i];	// i-boundary type
      bj = bcv[wbound][1][k][j][i - 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"v[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = vanal[ll][m][k][j][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = vanal[ll][m][k][j][i + 1];
	    workv[ll][m][k][j - 1][i] = vanal[ll][m][k][j - 1][i];
	    workv[ll][m][k][j - 1][i + 1] =
		vanal[ll][m][k][j - 1][i + 1];
	  }
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][N2 - 1][0];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][N2 - 1][1];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][N2 - 2][0];
	    workv[ll][m][k][j - 1][i + 1] = workv[ll][m][k][N2 - 2][1];
	  }
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j + 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][j + 1][i + 1];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][j + 1][i];
	    workv[ll][m][k][j - 1][i + 1] =
		workv[ll][m][k][j + 1][i + 1];
	  }
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i - 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][j][i - 1];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][j - 1][i - 1];
	    workv[ll][m][k][j - 1][i + 1] =
		workv[ll][m][k][j - 1][i - 1];
	  }
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j + 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][j + 1][i + 1];
	    workv[ll][m][k][j - 1][i] =
		workv[ll][m][k][j + jumpfactor][i];
	    workv[ll][m][k][j - 1][i + 1] =
		workv[ll][m][k][j + jumpfactor][i + 1];
	  }
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i - 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] =
		workv[ll][m][k][j][i - jumpfactor];
	    workv[ll][m][k][j - 1][i] = workv[ll][m][k][j - 1][i - 1];
	    workv[ll][m][k][j - 1][i + 1] =
		workv[ll][m][k][j - 1][i - jumpfactor];
	  }
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "v[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
      // i=N1,j=N2 quadrant
      i = N1;
      j = N2;
      k = 0;

      bi = bcv[wbound][1][k][j - 1][i];	// i-boundary type
      bj = bcv[wbound][1][k][j][i - 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"v[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = vanal[ll][m][k][j][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = vanal[ll][m][k][j][i + 1];
	    workv[ll][m][k][j + 1][i] = vanal[ll][m][k][j + 1][i];
	    workv[ll][m][k][j + 1][i + 1] =
		vanal[ll][m][k][j + 1][i + 1];
	  }
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][0][0];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][0][1];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][1][0];
	    workv[ll][m][k][j + 1][i + 1] = workv[ll][m][k][1][1];
	  }
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j - 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][j - 1][i + 1];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][j - 1][i];
	    workv[ll][m][k][j + 1][i + 1] =
		workv[ll][m][k][j - 1][i + 1];
	  }
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i - 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][j][i - 1];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][j + 1][i - 1];
	    workv[ll][m][k][j + 1][i + 1] =
		workv[ll][m][k][j + 1][i - 1];
	  }
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j - 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] = workv[ll][m][k][j - 1][i + 1];
	    workv[ll][m][k][j + 1][i] =
		workv[ll][m][k][j - jumpfactor][i];
	    workv[ll][m][k][j + 1][i + 1] =
		workv[ll][m][k][j - jumpfactor][i + 1];
	  }
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N1 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i - 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i + 1] =
		workv[ll][m][k][j][i - jumpfactor];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][j + 1][i - 1];
	    workv[ll][m][k][j + 1][i + 1] =
		workv[ll][m][k][j + 1][i - jumpfactor];
	  }
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "v[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
      // i=-1,j=N2 quadrant
      i = -1;
      j = N2;
      k = 0;

      bi = bcv[wbound][1][k][j - 1][i];	// i-boundary type
      bj = bcv[wbound][1][k][j][i + 1];	// j-boundary type
      if (((bi == 0) && (bj != 0)) || ((bj == 0) && (bi != 0))) {
	fprintf(fail_file,
		"v[%d]: j=%d,i=%d quad: undefined boundary bi: %d bj: %d\n",
		ll, j, i, bi, bj);
	myexit(1);
      } else if (!((bi == 0) && (bj == 0))) {
	// i-inflow AND j-inflow or if p-i or i-p 
	if (((bi == 3) && (bj == 3)) || ((bi == 5) && (bj == 3))
	    || ((bi == 3) && (bj == 5))) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = vanal[ll][m][k][j][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = vanal[ll][m][k][j][i - 1];
	    workv[ll][m][k][j + 1][i] = vanal[ll][m][k][j + 1][i];
	    workv[ll][m][k][j + 1][i - 1] =
		vanal[ll][m][k][j + 1][i - 1];
	  }
#endif
	}
	// i-periodic AND j-periodic
	else if ((bi == 5) && (bj == 5)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][0][N1 - 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][0][N1 - 2];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][1][N1 - 1];
	    workv[ll][m][k][j + 1][i - 1] = workv[ll][m][k][1][N1 - 2];
	  }
#endif
	}
	// i-outflow/inflow/periodic AND j-outflow (if o-o,
	// same
	// effect as next case)
	else if (((bi == 4) || (bi == 3) || (bi == 5))
		 && (bj == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j - 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][j - 1][i - 1];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][j - 1][i];
	    workv[ll][m][k][j + 1][i - 1] =
		workv[ll][m][k][j - 1][i - 1];
	  }
#endif
	}
	// j-outflow/inflow/periodic AND i-outflow
	else if (((bj == 4) || (bj == 3) || (bj == 5))
		 && (bi == 4)) {
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i + 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][j][i + 1];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][j + 1][i + 1];
	    workv[ll][m][k][j + 1][i - 1] =
		workv[ll][m][k][j + 1][i + 1];
	  }
#endif
	}
	// i-inflow/outflow/periodic/reflect(or AOS) AND
	// j-reflect(or AOS)
	else if (((bi == 3) || (bi == 4) || (bi == 5) || (bi == 1)
		  || (bi == 2)) && ((bj == 1) || (bj == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j - 1][i];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] = workv[ll][m][k][j - 1][i - 1];
	    workv[ll][m][k][j + 1][i] =
		workv[ll][m][k][j - jumpfactor][i];
	    workv[ll][m][k][j + 1][i - 1] =
		workv[ll][m][k][j - jumpfactor][i - 1];
	  }
#endif
	}
	// i-reflect(or AOS) AND
	// j-inflow/outflow/periodic/reflect(or AOS)
	else if (((bj == 3) || (bj == 4) || (bj == 5) || (bj == 1)
		  || (bj == 2)) && ((bi == 1) || (bi == 2))) {
#if(NBIGBND==2)
	  // correct for singular case
	  if (N2 == 1)
	    jumpfactor = 1;
	  else
	    jumpfactor = 2;
#endif
	  for (m = 1; m <= 3; m++) {
	    workv[ll][m][k][j][i] = workv[ll][m][k][j][i + 1];
#if(NBIGBND==2)
	    workv[ll][m][k][j][i - 1] =
		workv[ll][m][k][j][i + jumpfactor];
	    workv[ll][m][k][j + 1][i] = workv[ll][m][k][j + 1][i + 1];
	    workv[ll][m][k][j + 1][i - 1] =
		workv[ll][m][k][j + 1][i + jumpfactor];
	  }
#endif
	} else if ((bj == 99) || (bi == 99)) {
	} else {
	  fprintf(fail_file,
		  "v[%d]: j=%d,i=%d quad: No boundary case setup for bi: %d bj: %d\n",
		  ll, j, i, bi, bj);
	  myexit(1);
	}
      }
#endif





      /* cut short loop if only to do one */
      if (wvec != -1)
	l = NUMVEC;
    }				// end over vectors


#if(LINEXT==1)			// LINEXT==0 taken care of in above
    // automagically, should absorb this into
    // extrap routines
    // first setup special outflow conditions to avoid inflow
#if(INFLOWCHECKIX1)
    if ((wvec == 1) || (wvec == -3)) {
      if (wvec > 0)
	ll = wvec;
      else if (wvec == -3)
	ll = 0;

      for (k = -N3BND; k < N3 + N3BND; k++)
	for (j = -N2BND; j < N2 + N2BND; j++) {
	  if (workv[ll][1][k][j][1] > 0.0) {
	    workv[ll][1][k][j][0] = 0.0;
	    workv[ll][1][k][j][-1] = 0.0;
	  }
	}

    }
#endif
#if(INFLOWCHECKOX1)
    if ((wvec == 1) || (wvec == -3)) {
      if (wvec > 0)
	ll = wvec;
      else if (wvec == -3)
	ll = 0;

      for (k = -N3BND; k < N3 + N3BND; k++)
	for (j = -N2BND; j < N2 + N2BND; j++) {
	  if (workv[ll][1][k][j][N1 - 1] < 0.0) {
	    workv[ll][1][k][j][N1] = 0.0;
	  }
	}

    }
#endif

#if(INFLOWCHECKIX2)
    if ((wvec == 1) || (wvec == -3)) {
      if (wvec > 0)
	ll = wvec;
      else if (wvec == -3)
	ll = 0;

      for (k = -N3BND; k < N3 + N3BND; k++)
	for (i = -N1BND; i < N1 + N1BND; i++) {
	  if (workv[ll][2][k][1][i] > 0.0) {
	    workv[ll][2][k][0][i] = 0.0;
	    workv[ll][2][k][-1][i] = 0.0;
	  }
	}

    }
#endif
#if(INFLOWCHECKOX2)
    if ((wvec == 1) || (wvec == -3)) {
      if (wvec > 0)
	ll = wvec;
      else if (wvec == -3)
	ll = 0;

      for (k = -N3BND; k < N3 + N3BND; k++)
	for (i = -N1BND; i < N1 + N1BND; i++) {
	  if (workv[ll][2][k][N2 - 1][i] < 0.0) {
	    workv[ll][2][k][N2][i] = 0.0;
	  }
	}

    }
#endif

#if(INFLOWCHECKIX3)
    if ((wvec == 1) || (wvec == -3)) {
      if (wvec > 0)
	ll = wvec;
      else if (wvec == -3)
	ll = 0;

      for (j = -N2BND; j < N2 + N2BND; j++)
	for (i = -N1BND; i < N1 + N1BND; i++) {
	  if (workv[ll][3][1][j][i] > 0.0) {
	    workv[ll][3][0][j][i] = 0.0;
	    workv[ll][3][-1][j][i] = 0.0;
	  }
	}

    }
#endif
#if(INFLOWCHECKOX3)
    if ((wvec == 1) || (wvec == -3)) {
      if (wvec > 0)
	ll = wvec;
      else if (wvec == -3)
	ll = 0;

      for (j = -N2BND; j < N2 + N2BND; j++)
	for (i = -N1BND; i < N1 + N1BND; i++) {
	  if (workv[ll][3][N3 - 1][j][i] < 0.0) {
	    workv[ll][3][N3][j][i] = 0.0;
	  }
	}

    }
#endif
#endif				// end if outflowv==2

  }				// endif vectors to be done

#endif

}				// end function








void legacyblock1(int bct, int bcdim, int bcdir, int k, int j, int i,
		  int *bcd1p, int *bcd2p, int *kkp, int *jjp, int *iip,
		  int *kkkp, int *jjjp, int *iiip)
{

  if (bcdir == 1) {
    (*bcd1p) = 1;
    (*bcd2p) = 0;
  } else if (bcdir == (-1)) {
    (*bcd1p) = 0;
    (*bcd2p) = 1;
  } else {
    fprintf(fail_file,
	    "error: %d %d %d bound.c: bcdir out of bounds: %d\n\r", k,
	    j, i, bcdir);
    myexit(1);
  }
  if (bct != 5) {
    /* Determine which local cell is relevant */
    if (bcdim == 1) {
      *iip = i + bcdir;
      *jjp = j;
      *kkp = k;
      if (N1BND == 1) {
	*iiip = i + bcdir;
	*jjjp = j;
	*kkkp = k;
      } else {
	*iiip = i + 2 * bcdir;
	*jjjp = j;
	*kkkp = k;
      }
    } else if (bcdim == 2) {
      *iip = i;
      *jjp = j + bcdir;
      *kkp = k;
      if (N2BND == 1) {
	*iiip = i;
	*jjjp = j + bcdir;
	*kkkp = k;
      } else {
	*iiip = i;
	*jjjp = j + 2 * bcdir;
	*kkkp = k;
      }
    } else if (bcdim == 3) {	/* only has impact with 3-d problems */
      *iip = i;
      *jjp = j;
      *kkp = k + bcdir;
      if (N3BND == 1) {
	*iiip = i;
	*jjjp = j;
	*kkkp = k + bcdir;
      } else {
	*iiip = i;
	*jjjp = j;
	*kkkp = k + 2 * bcdir;
      }
    } else {
      fprintf(fail_file,
	      "error: bound.c: bcdim out of bounds on [%d][%d][%d]: %d\n\r",
	      k, j, i, bcdim);
      fprintf(fail_file, "bct: %d  bcdir: %d\n\r", bct, bcdir);
      myexit(1);
    }
  } else {

    /* Determine which zone to copy in period conditions */

    /* These non-local conditions are not general to any problem */
    /* The below assumes edges of rect-domain are periodic--standard */

    if (bcdim == 1) {
      *iip = i + bcdir * N1;
      *jjp = j;
      *kkp = k;
    } else if (bcdim == 2) {
      *iip = i;
      *jjp = j + bcdir * N2;
      *kkp = k;
    } else if (bcdim == 3) {
      *iip = i;
      *jjp = j;
      *kkp = k + bcdir * N3;
    } /* only has impact with 3-d problems */
    else {
      fprintf(fail_file,
	      "error: bound.c: bcdim out of bounds on [%d][%d][%d]: %d\n\r",
	      k, j, i, bcdim);
      fprintf(fail_file, "bct: %d  bcdir: %d\n\r", bct, bcdir);
      myexit(1);
    }
  }
}				// end function
