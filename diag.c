#include "global.h"
#include "defs.h"
/* diagnostics subroutine */

#define NUMINPUTS 10		// see how used! (overestimated)

#if(SENSITIVE==1)
#define INPUT3 "%lf"
#define INPUT4 " %lf %lf %lf"
#define INPUT5 " %lf %lf %lf %lf %lf %lf %lf"
#define INPUT6 " %lf %lf %lf %lf %lf"
#else
#define INPUT3 "%f"
#define INPUT4 " %f %f %f"
#define INPUT5 " %f %f %f %f %f %f %f"
#define INPUT6 " %f %f %f %f %f"
#endif




void diag(int call_code)
		 // call_code
		 // -1: do every call to diag(usually each time step as 
		 // called by
		 // main.c)
		 // 0: first call setup/dump
		 // 1: normal diag dump of variables each DTl (DT for
		 // logging)
		 // 2: final diag call after done
{
  SFTYPE tempf;
  SFTYPE tempfns;
  SFTYPE masstempc, masstemp1, masstemp2;
  char dfnam[MAXFILENAME];
  char dfnamtemp[MAXFILENAME];
  char dfnamback[MAXFILENAME];

  SFTYPE mass_full = 0, eth_full = 0, ek_full = 0, eg_full =
      0, eb_full = 0, se_full = 0, angmom_full[3 + 1] =
      { 0 }, cmode_amp_full = 0, smode_amp_full =
      0, floors_full[NUMLOSSVAR + 1] = {
  0}, inflows_full[NUMLOSSVAR + 1] = {
  0}, radiations_full[NUMLOSSVAR + 1] = {
  0};				// 0: 
  // etot
  SFTYPE ek, eth, eg, eb, se, cmode_amp, smode_amp, mass, angmom[3 + 1];	// 0: 
										// 
  // 
  // etot

  SFTYPE eki, ebi, ethi, egi;
  SFTYPE vxa, vya, vza;
  SFTYPE dxdyc, dxdy1, dxdy2;
  FILE *efboth = 0;
  FILE *dump_file;
  FILE *image_file;
  int i, j, k, l, m = 0;
  int loopdido;
  int dumi[10];
  SFTYPE tcheck;		// used to check time for loss append

  SFTYPE realloss;		// need for enthalpy and sm

  SFTYPE floattemp[NUMLOSSVAR + 1][NUMINPUTS];	// NUMINPUTS types 
  // of output for
  // 0_ener.dat per
  // variable
  static SFTYPE varinit[NUMLOSSVAR + 1];	// 0: etot
  static SFTYPE varfinal[NUMLOSSVAR + 1];	// 0: etot
  static char losstext[NUMLOSSVAR + 1][30];	// 0: etot

  static FTYPE tdump, tpdump, timage, tloss, tener;
  static FILE *ener_file;
  static FILE *ener_file_temp;
  static int dump_cnt, pdump_cnt, im_cnt;	// used to enumerate
						// files
  static int adump_cnt, npdump_cnt, floordump_cnt;
  static int dumpc, pdumpc, imagec, lossc;	// time-based counts
						// for
  // this run
  static int enerc;
  SFTYPE tempsumloss[NUMLOSSVAR + 1][4];
  static SFTYPE sumloss[NUMLOSSVAR + 1][4];
  static SFTYPE totloss[NUMLOSSVAR + 1];	// 0: etot
  static SFTYPE sumloss_full[NUMLOSSVAR + 1][4];
  static SFTYPE totloss_full[NUMLOSSVAR + 1];	// 0: etot
  static FILE *loss_file;
  static FILE *loss_file_temp;
  static FILE *loss_file2;
  static FILE *final_output;
  static int firsttime0 = 1;
  static int firsttimem1 = 1;

  char temps[MAXFILENAME];
  char filename[MAXFILENAME];
  char filenametemp[MAXFILENAME];
  char filenameback[MAXFILENAME];
  long fpos0 = 0;
  int gotit;


  if (DOLOSSDIAG) {

    // if every dt or initial or final, do special diags not
    // otherwise 
    // done
    if (call_code == -1) {	// don't do at end(2) since already
      // done
      // for last timestep
      if (firsttimem1 == 1) {	// the check for if first time
	// here
	for (i = 1; i <= NUMLOSSVAR; i++) {
	  totloss[i] = 0.0;
	}
	if (myid <= 0) {

	  if (DETAILMLOSS >= 0) {

	    sprintf(temps, DATADIR);

	    sprintf(dfnam, "%s0_loss%s", temps, DATEXT);
	    if ((loss_file = fopen(dfnam, WRITETYPE)) == NULL) {
	      fprintf(fail_file,
		      "error opening loss output file %s\n", dfnam);
	      myexit(1);
	    }
	    if (appendold == 0) {
	      // version header
	      fprintf(loss_file, "#%10s\n%10d %10d\n",
		      "LOSSVER", LOSSVER, LOSSTYPE);
	      if (DETAILMLOSS == 0) {
		fprintf(loss_file,
			"#%21s %21s %21s %21s %21s\n",
			"time", "Total Mass Loss",
			"Total Energy Loss",
			"Total AngMom Loss", "Total ViscE Loss");
	      } else if (DETAILMLOSS >= 1) {
		fprintf(loss_file, "#%21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s"
			" %21s %21s %21s %21s %21s\n",
			"time", "Total Mass Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total IEnergy Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total PEnergy Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total KEnergy Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total x1-Mom Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total x2-Mom Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total x3-Mom Loss",
			"Inner x1", "Outer x1", "Inner x2",
			"Outer x2", "Total visc-e Loss",
			"Inner x1", "Outer x1", "Inner x2", "Outer x2");
	      }
	    }			// if appendold==0
	    else {		// if appendold==1
	      if (DETAILMLOSS < 1) {
		if (myid <= 0) {
		  fprintf(logfull_file,
			  "appendold==1 with detailmassloss<1 won't allow correct append reentrance to loss data file\n");
		  fflush(logfull_file);
		}
	      } else {		// if detailmassloss>=1
		if (myid <= 0) {
		  fprintf(logfull_file,
			  "Start setup of loss file append\n");
		  fflush(logfull_file);
		}
		// need to read in totloss and sumloss
		// and 
		// stick on right cpu
		// need to make sure getting right
		// time.
		// If DTloss<DTd||DTi, then should be
		// able 
		// to get easily.  You should generally
		// have this true anyways

		rewind(loss_file);	// go to start
		// check version info
		while (fgetc(loss_file) != '\n');	// skip 
		// 
		// comment 
		// line
		fscanf(loss_file, "%d %d\n", &dumi[0], &dumi[1]);
		if ((dumi[0] != LOSSVER)
		    || (dumi[1] != LOSSTYPE)) {
		  fprintf(fail_file,
			  "Expected lossver/losstype: %d %d got %d %d\n",
			  LOSSVER, LOSSTYPE, dumi[0], dumi[1]);
		  myexit(6);
		}
		while (fgetc(loss_file) != '\n');	// skip 
		// 
		// comment 
		// line 
		// 
		gotit = 0;
		while ((!feof(loss_file)) && (gotit == 0)) {

		  fpos0 = ftell(loss_file);	// position 
		  // 
		  // to
		  // continue 
		  // writting 
		  // at if
		  // successful 
		  // get
		  fscanf(loss_file, INPUT3, &tcheck);
		  // fprintf(stderr,"%15.10g %15.10g
		  // %d\n",t,tcheck,gotit);
		  if (fabs(tcheck - t) < 1.0E-8 * t + 1.0E-6) {
		    gotit = 1;
		    for (l = 1; l <= NUMLOSSVAR; l++) {
		      if (l == NUMSCA + 1 + 3 + 1) {	// skip 
			// 
			// B 
			// field 
			// for 
			// now, 
			// go 
			// directly 
			// to 
			// visc 
			// energy
			l = NUMSCA + 1 + NUMVEC * 3 + 1;
		      }
		      fscanf(loss_file, INPUT6,
			     &totloss_full[l],
			     &sumloss_full[l][0],
			     &sumloss_full[l][1],
			     &sumloss_full[l][2], &sumloss_full[l][3]);
		    }
		  } else {
		    while ((fgetc(loss_file) != '\n') && (!feof(loss_file)));	// skip 
										// 
		    // 
		    // this 
		    // bad 
		    // line
		  }
		}
		if (gotit == 0) {
		  fprintf(fail_file,
			  "Never found right time in loss file when appending\n");
		  myexit(1);
		} else {
		  sprintf(temps, DATADIR);
		  sprintf(dfnamtemp, "%s0_loss%s.temp", temps, DATEXT);
		  sprintf(dfnam, "%s0_loss%s", temps, DATEXT);
		  sprintf(dfnamback, "%s0_loss%s.back", temps, DATEXT);

		  // now that done, fix up file
		  if ((loss_file_temp = fopen(dfnamtemp, "wt")) == NULL) {
		    fprintf(fail_file,
			    "Cannot open temp loss file for appending: %s\n",
			    dfnamtemp);
		    myexit(1);
		  } else {
		    rewind(loss_file);
		    while (ftell(loss_file) < fpos0) {
		      fputc(fgetc(loss_file), loss_file_temp);
		    }
		    fclose(loss_file_temp);
		    fclose(loss_file);
		    rename(dfnam, dfnamback);	// move 
		    // 
		    // old 
		    // to 
		    // backup 
		    // location
		    rename(dfnamtemp, dfnam);	// move 
		    // 
		    // new 
		    // to 
		    // old 
		    // name(normal 
		    // name)
		    // reopen loss_file (now normal
		    // name)
		    if ((loss_file = fopen(dfnam, "at")) == NULL) {
		      fprintf(fail_file,
			      "2: error opening loss output file %s\n",
			      dfnam);
		      myexit(1);
		    }
		    if (myid <= 0) {
		      fprintf(logfull_file,
			      "End setup of loss file append\n");
		      fflush(logfull_file);
		    }
		  }
		}
	      }			// end else if detailmassloss>=1
	    }			// end else if appendold==1
	  }			// end if doing normal mloss output

	  if (DETAILMLOSS == 2) {


	    sprintf(temps, DATADIR);

	    sprintf(dfnam, "%s0_lossd%s", temps, DATEXT);
	    if ((loss_file2 = fopen(dfnam, WRITETYPE)) == NULL) {
	      fprintf(fail_file,
		      "error opening loss detail output file %s\n",
		      dfnam);
	      myexit(1);
	    }

	    if (appendold == 0) {
	      fprintf(loss_file2, "#%21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s"
		      " %7s %7s %21s\n", "time", "i-coord",
		      "j-coord", "Mass Loss", "i-coord",
		      "j-coord", "IEnergy Loss", "i-coord",
		      "j-coord", "PEnergy Loss", "i-coord",
		      "j-coord", "KEnergy Loss", "i-coord",
		      "j-coord", "s1 loss", "i-coord",
		      "j-coord", "s2 loss", "i-coord",
		      "j-coord", "s3 loss", "i-coord",
		      "j-coord", "visc-e loss");
	    }
	  }
	}			// end if write cpu
	if (appendold == 1) {
	  // now must distribute loss data to appropriate cpu
	  // from root==0
	  // given only x2-dir splitting
	  if (myid <= 0) {
	    fprintf(logfull_file,
		    "Begin transfer setup of loss file append\n");
	    fflush(logfull_file);
	  }
	  for (l = 1; l <= NUMLOSSVAR; l++) {
	    if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	      // for now, go
	      // directly to
	      // visc energy
	      l = NUMSCA + 1 + NUMVEC * 3 + 1;
	    }
	    if (myid <= 0) {
	      // cpu=0 can keep total losses
	      totloss[l] = totloss_full[l];
	      // cpu=0 can keep radial fluxes(no way to
	      // recapture that seperate data anyways),
	      // and
	      // cpu=0 needs to keep theta=0 surface flux
	      sumloss[l][0] = sumloss_full[l][0];
	      sumloss[l][1] = sumloss_full[l][1];
	      sumloss[l][2] = sumloss_full[l][2];
	    }
	    if (numprocs > 1) {
	    } else {
	      sumloss[l][3] = sumloss_full[l][3];
	    }
	  }			// over loss loop
	  if (myid <= 0) {
	    fprintf(logfull_file,
		    "End transfer setup of loss file append\n");
	    fflush(logfull_file);
	  }
	}			// if appendold==1
	lossc = (int) ((t - tstart) / DTloss) + ireenter;
	tloss = tstart + (FTYPE) (lossc) * DTloss - 1.E-12;	// next 
								// 
	// 
	// loss 
	// time

      }				// end if first time in here 


      // continue with call_code==-1 for every call
      // compute total mass lost in this timestep
      for (l = 1; l <= NUMLOSSVAR; l++) {	// skip B field for now 
						// 
	// 

	if (l <= NUMSCA)
	  m = -1;
	else if (l == NUMSCA + 1)
	  m = 0;
	else if (l == NUMSCA + 1 + 1)
	  m = 1;
	else if (l == NUMSCA + 1 + 2)
	  m = 2;
	else if (l == NUMSCA + 1 + 3)
	  m = 3;
	else if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	  // for now, go
	  // directly to
	  // visc energy
	  l = NUMSCA + 1 + NUMVEC * 3 + 1;
	  m = -2;
	}

	j = 1;			// x1-direction
	for (k = 0; k < 2; k++) {	// i/o
	  // add up stuff that just went through boundary
	  tempsumloss[l][k] = 0;

	  for (i = 0; i < N2; i++) {	// list
	    if (m == -1)
	      tempf = losss[l][j][k][i];
	    else if (m == -2)
	      tempf = lossvisc[1][j][k][i];
	    else
	      tempf = lossv[1][m][j][k][i];
	    tempsumloss[l][k] += tempf;
	  }
	  // add new quantity passing through boundary
	  sumloss[l][k] += tempsumloss[l][k];
	}

	j = 2;			// x2-direction
	for (k = 0; k < 2; k++) {	// i/o
	  tempsumloss[l][k + 2] = 0;

	  if (((k == 0) && (myid == 0))
	      || ((k == 1) && (myid == numprocs - 1))) {
	    for (i = 0; i < N1; i++) {	// list
	      if (m == -1)
		tempf = losss[l][j][k][i];
	      else if (m == -2)
		tempf = lossvisc[1][j][k][i];
	      else
		tempf = lossv[1][m][j][k][i];
	      tempsumloss[l][k + 2] += tempf;
	    }
	    sumloss[l][k + 2] += tempsumloss[l][k + 2];
	  }
	}

	// only add new stuff to total sum
	for (k = 0; k < 4; k++) {
	  totloss[l] += tempsumloss[l][k];
	}
      }
    }				// end sum up every time step


    if ((t >= tener) || (t >= tloss) || (call_code == 2)) {	// must
      // compute 
      // sum
      // over
      // cpus
      // everytime 
      // since
      // needed
      // for
      // diagnostics 
      // on DTl
      // period, 
      // or for
      // loss
      // diagnostics
      for (l = 1; l <= NUMLOSSVAR; l++) {	// skip B field for now 
						// 
	// 

	if (l <= NUMSCA)
	  m = -1;
	else if (l == NUMSCA + 1)
	  m = 0;
	else if (l == NUMSCA + 1 + 1)
	  m = 1;
	else if (l == NUMSCA + 1 + 2)
	  m = 2;
	else if (l == NUMSCA + 1 + 3)
	  m = 3;
	else if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	  // for now, go
	  // directly to
	  // visc energy
	  l = NUMSCA + 1 + NUMVEC * 3 + 1;
	  m = -2;
	}

	if (numprocs > 1) {
	} else {
	  totloss_full[l] = totloss[l];
	  for (k = 0; k < 4; k++) {
	    sumloss_full[l][k] = sumloss[l][k];
	  }
	}
      }
    }
    // only dump loss data when wanted
    if ((t >= tloss) || (call_code == 2)) {
      if (myid <= 0) {
	// now output to file
	if (DETAILMLOSS == 0) {
	  fprintf(loss_file, " %21.15g", t);
	  for (l = 1; l <= NUMLOSSVAR; l++) {
	    if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	      // for now, go
	      // directly to
	      // visc energy
	      l = NUMSCA + 1 + NUMVEC * 3 + 1;
	    }
	    // totloss_full[l] holds total loss so far of
	    // variable l on entire grid through
	    // boundaries.
	    fprintf(loss_file, " %21.15g", totloss_full[l]);
	  }
	  fprintf(loss_file, "\n");
	  fflush(loss_file);
	} else if (DETAILMLOSS >= 1) {
	  fprintf(loss_file, " %21.15g", t);
	  for (l = 1; l <= NUMLOSSVAR; l++) {
	    if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	      // for now, go
	      // directly to
	      // visc energy
	      l = NUMSCA + 1 + NUMVEC * 3 + 1;
	    }
	    // sumloss_full[l][k] holds total loss so far
	    // of
	    // variable l on boundary k where k=0 is inner
	    // r-edge, k=1 outer r-edge, k=2 inner theta
	    // edge, 
	    // k=3 outer theta edge
	    fprintf(loss_file,
		    " %21.15g %21.15g %21.15g %21.15g %21.15g",
		    totloss_full[l], sumloss_full[l][0],
		    sumloss_full[l][1], sumloss_full[l][2],
		    sumloss_full[l][3]);
	  }
	  fprintf(loss_file, "\n");
	  fflush(loss_file);
	}

	if (DETAILMLOSS == 2) {
	  j = 1;
	  for (k = 0; k < 2; k++) {	// i/o
	    for (i = 0; i < N2; i++) {	// list
	      fprintf(loss_file2, " %21.15g", t);
	      for (l = 1; l <= NUMLOSSVAR; l++) {

		if (l <= NUMSCA)
		  m = -1;
		else if (l == NUMSCA + 1)
		  m = 0;
		else if (l == NUMSCA + 1 + 1)
		  m = 1;
		else if (l == NUMSCA + 1 + 2)
		  m = 2;
		else if (l == NUMSCA + 1 + 3)
		  m = 3;
		else if (l == NUMSCA + 1 + 3 + 1) {	// skip 
		  // 
		  // B 
		  // field 
		  // for 
		  // now, 
		  // go 
		  // directly 
		  // to 
		  // visc 
		  // energy
		  l = NUMSCA + 1 + NUMVEC * 3 + 1;
		  m = -2;
		}
		// loss holds this dt transported
		// variable 
		// l across boundary j=direction(1=x1,
		// 2=x2, 3=x3) k(=0 inner =1, outer)
		if (m == -1)
		  tempf = losss[l][j][k][i];
		else if (m == -2)
		  tempf = lossvisc[1][j][k][i];
		else
		  tempf = lossv[1][m][j][k][i];
		fprintf(loss_file2, " %7d %7d %21.15g",
			k * N1, i, tempf);
	      }
	      fprintf(loss_file2, "\n");
	    }
	  }

	  j = 2;
	  for (k = 0; k < 2; k++) {	// i/o
	    for (i = 0; i < N1; i++) {	// list
	      fprintf(loss_file2, " %21.15g", t);
	      for (l = 1; l <= NUMLOSSVAR; l++) {

		if (l <= NUMSCA)
		  m = -1;
		else if (l == NUMSCA + 1)
		  m = 0;
		else if (l == NUMSCA + 1 + 1)
		  m = 1;
		else if (l == NUMSCA + 1 + 2)
		  m = 2;
		else if (l == NUMSCA + 1 + 3)
		  m = 3;
		else if (l == NUMSCA + 1 + 3 + 1) {	// skip 
		  // 
		  // B 
		  // field 
		  // for 
		  // now, 
		  // go 
		  // directly 
		  // to 
		  // visc 
		  // energy
		  l = NUMSCA + 1 + NUMVEC * 3 + 1;
		  m = -2;
		}

		if (m == -1)
		  tempf = losss[l][j][k][i];
		else if (m == -2)
		  tempf = lossvisc[1][j][k][i];
		else
		  tempf = lossv[1][m][j][k][i];
		fprintf(loss_file2, " %7d %7d %21.15g", i,
			k * N2, tempf);
	      }
	      fprintf(loss_file2, "\n");
	    }
	  }
	  fflush(loss_file2);
	}
      }				// end if cpu that writes
      lossc++;
      tloss = tstart + (FTYPE) (lossc) * DTloss;
    }				// end if time to output or final
    // output
    if (call_code == -1) {
      // initialize loss data, so don't have to worry about how
      // added in routines!
      // this is pretty inexpensive, esp. for 1-D problems, can
      // fix
      // by knowing what adds first
      init_loss();
    }
    firsttimem1 = 0;
  }				// end if want loss diag


  if (call_code == 0) {
    if (firsttime0 == 1) {
      enerc = (int) ((t - tstart) / DTener) + ireenter;
      tener = tstart + (FTYPE) (enerc) * DTener - 1.E-12;	// next
      // ener
      // time
    }
    if (myid <= 0) {
      // setup labels for each varinit[] and varfinal[]
      strcpy(losstext[0], "ETotal");
      strcpy(losstext[1], "Mass");
      strcpy(losstext[2], "IE");
      strcpy(losstext[3], "Pot");
      strcpy(losstext[NUMSCA + 1], "KE");
      strcpy(losstext[NUMSCA + 1 + 1], "AngMom1");
      strcpy(losstext[NUMSCA + 1 + 2], "AngMom2");
      strcpy(losstext[NUMSCA + 1 + 3], "AngMom3");
      strcpy(losstext[NUMSCA + 1 + NUMVEC * 3 + 1], "ViscEn");

      sprintf(temps, DATADIR);
      sprintf(filename, "%s0_ener%s", temps, DATEXT);
      if ((ener_file = fopen(filename, WRITETYPE)) == NULL) {
	fprintf(fail_file, "error opening energy output file %s\n",
		filename);
	myexit(1);
      }

      sprintf(temps, DATADIR);
      sprintf(filename, "%s0_final%s", temps, DATEXT);
      if ((final_output = fopen(filename, "w")) == NULL) {
	fprintf(fail_file, "error opening final output file %s\n",
		filename);
	myexit(1);
      }
    }

    pdump_cnt = pdump_start;
    dump_cnt = dump_start;
    npdump_cnt = npdump_start;
    adump_cnt = adump_start;
    floordump_cnt = floor_start;
    im_cnt = image_start;

    pdumpc = (int) ((t - tstart) / DTpd) + ireenter;
    tpdump = tstart + ((FTYPE) (pdumpc) * DTpd) - 1.E-12;
    dumpc = (int) ((t - tstart) / DTd) + ireenter;
    tdump = tstart + ((FTYPE) (dumpc) * DTd) - 1.E-12;
    imagec = (int) ((t - tstart) / DTi) + ireenter;
    timage = tstart + ((FTYPE) (imagec) * DTi) - 1.E-12;

    if (myid <= 0) {
      if (appendold == 0) {
	// version header
	fprintf(ener_file, "#%10s\n%10d %10d\n", "ENERVER",
		ENERVER, ENERTYPE);
	fprintf(ener_file, "#%21s %21s %21s", "t", "cmode_amp",
		"smode_amp");
	for (l = 0; l <= NUMLOSSVAR; l++) {
	  if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	    // for now, go
	    // directly to
	    // visc energy
	    l = NUMSCA + 1 + NUMVEC * 3 + 1;
	  }

	  fprintf(ener_file,
		  " %20s%1d %20s%1d %20s%1d %20s%1d %20s%1d %20s%1d %20s%1d",
		  losstext[l], l, "Totalchange", l, "b_lost", l,
		  "fl_add", l, "inject_add", l, "radiate", l,
		  "totdel+bloss-fladd", l);
	}
	fprintf(ener_file, "\n");
      } else {			// if appendold==1
	if (myid <= 0) {
	  fprintf(logfull_file, "Start setup of energy file append\n");
	  fflush(logfull_file);
	}
	// need to read in first(or any) data line of
	// 0_ener.dat
	// file for varinit[l] values

	rewind(ener_file);	// go to start
	// check on version info
	while (fgetc(ener_file) != '\n');	// skip comment
	// line 
	fscanf(ener_file, "%d %d\n", &dumi[0], &dumi[1]);
	if ((dumi[0] != ENERVER) || (dumi[1] != ENERTYPE)) {
	  fprintf(fail_file,
		  "Expected enerver/enertype: %d %d got %d %d\n",
		  ENERVER, ENERTYPE, dumi[0], dumi[1]);
	  myexit(6);
	}

	while (fgetc(ener_file) != '\n');	// skip comment
	// line 
	gotit = 0;
	while ((!feof(ener_file)) && (gotit == 0)) {

	  fpos0 = ftell(ener_file);	// position to continue
	  // writting at if
	  // successful get
	  fscanf(ener_file, INPUT4, &tcheck, &tempfns, &tempfns);
	  // fprintf(stderr,"%21.15g %21.15g\n",tcheck,t);
	  if (fabs(tcheck - t) < 1.0E-8 * t + 1.0E-6) {
	    gotit = 1;
	    // read in init values
	    for (l = 0; l <= NUMLOSSVAR; l++) {
	      if (l == NUMSCA + 1 + 3 + 1) {	// skip B
		// field
		// for
		// now, go 
		// directly 
		// to visc 
		// energy
		l = NUMSCA + 1 + NUMVEC * 3 + 1;
	      }
	      fscanf(ener_file, INPUT5, &floattemp[l][0],
		     &floattemp[l][1], &floattemp[l][2],
		     &floattemp[l][3], &floattemp[l][4],
		     &floattemp[l][5], &floattemp[l][6]);
	      // now compute init values from last value
	      // and 
	      // difference
	      // varfinal computed already above, and
	      // correct
	      varinit[l] = floattemp[l][0] - floattemp[l][1];
	      // floors/inflows isn't really just on this
	      // CPU, but the final counting doesn't care
	      // and _full version is correct(and
	      // overrides
	      // collective above), and needed for later
	      // calcs
	      floors[l] = floors_full[l] = floattemp[l][3];
	      inflows[l] = inflows_full[l] = floattemp[l][4];
	      radiations[l] = radiations_full[l] = floattemp[l][5];
	    }
	  }			// if good time
	  else {
	    while ((fgetc(ener_file) != '\n') && (!feof(ener_file)));	// skip 
									// 
	    // 
	    // this 
	    // bad 
	    // line
	  }
	}			// while finding good time or hitting
	// end
	// of file
	if (gotit == 0) {
	  fprintf(fail_file,
		  "Never found right time in energy file when appending\n");
	  myexit(1);
	} else {
	  sprintf(filename, "%s0_ener%s", DATADIR, DATEXT);
	  sprintf(filenametemp, "%s0_ener%s.temp", DATADIR, DATEXT);
	  sprintf(filenameback, "%s0_ener%s.back", DATADIR, DATEXT);

	  if ((ener_file_temp = fopen(filenametemp, "wt")) == NULL) {
	    fprintf(fail_file,
		    "error opening temp energy output file %s\n",
		    filenametemp);
	    myexit(1);
	  } else {
	    rewind(ener_file);
	    while (ftell(ener_file) < fpos0) {
	      fputc(fgetc(ener_file), ener_file_temp);
	    }
	    fclose(ener_file_temp);
	    fclose(ener_file);
	    rename(filename, filenameback);	// move old to
	    // backup location
	    rename(filenametemp, filename);	// move new to
	    // old(normal)
	    // reopen ener_file
	    if ((ener_file = fopen(filename, "at")) == NULL) {
	      fprintf(fail_file,
		      "error opening energy output file %s\n",
		      filename);
	      myexit(1);
	    }
	    if (myid <= 0) {
	      fprintf(logfull_file,
		      "End setup of energy file append\n");
	      fflush(logfull_file);
	    }
	  }
	}			// end else if gotit==1
      }				// end else if appendold==1
    }				// end if myid<=0
  }				// end if call_code==0

  // add up, combine, and output integrated variables
  if ((call_code >= 0) && ((t >= tener) || (call_code == 2))) {

    // initialize variables for this dt step
    mass = 0;
    angmom[1] = 0;
    angmom[2] = 0;
    angmom[3] = 0;

    /* calculate energies */
    se = 0.;
    eb = 0.;
    ebi = 0;
    ek = 0.;
    eth = 0.;
    eg = 0.;
    cmode_amp = 0.;
    smode_amp = 0.;

    LOOPINT {			// begin loop over real domain to get
      // diagnostics

      vxa = e2z_1(v[1][1][k], j, i);
      vya = e2z_2(v[1][2][k], j, i);
      vza = v[1][3][k][j][i];	// When fake-3d


      eki = 0.5 * s[1][k][j][i] * (vxa * vxa + vya * vya + vza * vza);

      egi = s[1][k][j][i] * s[3][k][j][i];

      /* Equation of state */
      if (press == 1) {
	if (wgam)
	  ethi = s[2][k][j][i];
	else
	  ethi = cs * cs * s[1][k][j][i];
      } else {
	ethi = 0;
      }

      dxdyc = dvl[1][1][i] * dvl[1][2][j];
      dxdy1 = dvl[2][1][i] * dvl[1][2][j];
      dxdy2 = dvl[1][1][i] * dvl[2][2][j];
      masstempc = s[1][k][j][i] * dxdyc;	// mass in center of a
      // zone
      masstemp1 = z2e_1(s[1][k], j, i) * dxdy1;	// mass in zone
      // around edge 1
      masstemp2 = z2e_2(s[1][k], j, i) * dxdy2;	// mass in zone
      // around edge 2


      // sum up the zone values to totals
      ek += eki * dxdyc;
      // eb += ebi*dxdyc ;
      eth += ethi * dxdyc;
      eg += egi * dxdyc;

      se += s[2][k][j][i] / s[1][k][j][i];

      cmode_amp +=
	  cos(3. * 2. * M_PI * x[2][1][i] / L[2][1]) * masstempc;
      smode_amp +=
	  sin(3. * 2. * M_PI * x[2][1][i] / L[2][2]) * masstempc;

      mass += masstempc;
      angmom[1] += masstemp1 * v[1][1][k][j][i];
      angmom[2] += masstemp2 * v[1][2][k][j][i] * g[2][2][i];
      angmom[3] +=
	  masstempc * v[1][3][k][j][i] * g[2][3][i] * g[2][4][j];

    }

    if (numprocs > 1) {
    } else {
      cmode_amp_full = cmode_amp;
      smode_amp_full = smode_amp;
      mass_full = mass;
      eth_full = eth;
      eg_full = eg;
      ek_full = ek;
      eb_full = eb;
      se_full = se;
      angmom_full[1] = angmom[1];
      angmom_full[2] = angmom[2];
      angmom_full[3] = angmom[3];
      for (l = 1; l <= NUMLOSSVAR; l++) {
	floors_full[l] = floors[l];
	inflows_full[l] = inflows[l];
	radiations_full[l] = radiations[l];
      }
    }

    // now write results
    if (myid <= 0) {
      if (call_code == 0) {
	if (appendold == 0) {
	  // for call_code==0, varfinal final so far!
	  varinit[1] = mass_full + SSMALL;
	  varinit[2] = eth_full + SSMALL;
	  varinit[3] = eg_full + SSMALL;
	  varinit[NUMSCA + 1] = ek_full + SSMALL;
	  varinit[NUMSCA + 1 + 1] = angmom_full[1] + SSMALL;
	  varinit[NUMSCA + 1 + 2] = angmom_full[2] + SSMALL;
	  varinit[NUMSCA + 1 + 3] = angmom_full[3] + SSMALL;
	  // sum up all energy terms from totals
	  varinit[0] =
	      varinit[2] + varinit[3] + varinit[NUMSCA + 1] + eb_full;
	}			// end if callcode==0
      }
      // final so far, in any call_code
      varfinal[1] = mass_full;
      varfinal[2] = eth_full;
      varfinal[3] = eg_full;
      varfinal[NUMSCA + 1] = ek_full;
      varfinal[NUMSCA + 1 + 1] = angmom_full[1];
      varfinal[NUMSCA + 1 + 2] = angmom_full[2];
      varfinal[NUMSCA + 1 + 3] = angmom_full[3];
      // sum up all energy terms from totals
      varfinal[0] =
	  varfinal[2] + varfinal[3] + varfinal[NUMSCA + 1] + eb_full;

      // need to do both final out and ener out when call_code==2
      if (call_code == 2) {
	loopdido = 2;
      } else
	loopdido = 1;

      for (i = 1; i <= loopdido; i++) {
	if (i == 1)
	  efboth = ener_file;
	else if (i == 2)
	  efboth = final_output;

	if (i == 2)
	  fprintf(efboth, "#%21s %21s %21s\n", "t", "cmode_amp",
		  "smode_amp");
	fprintf(efboth, " %21.15g %21.15g %21.15g", t,
		cmode_amp_full, smode_amp_full);
	if (i == 2)
	  fprintf(efboth, "\n");
	if (i == 2)
	  fprintf(efboth,
		  "#%21s %21s %21s %21s %21s %21s %21s %21s\n",
		  "what", "value", "Totalchange", "b_lost",
		  "fl_add", "minject", "radiate", "totdel+bloss-fladd");

	// compute total energy variables from all data
	totloss_full[0] = totloss_full[2] + totloss_full[3] + totloss_full[NUMSCA + 1] + totloss_full[NUMSCA + 1 + NUMVEC * 3 + 1];	// [2] 
																	// 
	// 
	// is 
	// really 
	// enthalpy 
	// as 
	// it 
	// should 
	// be
	floors_full[0] =
	    floors_full[2] + floors_full[3] + floors_full[NUMSCA + 1];
	inflows_full[0] =
	    inflows_full[2] + inflows_full[3] +
	    inflows_full[NUMSCA + 1];
	radiations_full[0] =
	    radiations_full[2] + radiations_full[3] +
	    radiations_full[NUMSCA + 1];

	for (l = 0; l <= NUMLOSSVAR; l++) {
	  if (l == NUMSCA + 1 + 3 + 1) {	// skip B field
	    // for now, go
	    // directly to
	    // visc energy
	    l = NUMSCA + 1 + NUMVEC * 3 + 1;
	  }

	  if (i == 2)
	    fprintf(efboth, " %21s", losstext[l]);
	  // below line only true for ideal EOS
	  if (l == 2)
	    realloss = totloss_full[l] / gam;	// totloss_full 
	  // 
	  // holds
	  // enthalpy 
	  // across
	  // boundary, 
	  // so
	  // change
	  // to ie
	  else
	    realloss = totloss_full[l];
	  if (fabs(realloss) < 1.0E-99)
	    realloss = SSMALL;	// needed for sm to not
	  // die when using float
	  // version of sm
	  fprintf(efboth,
		  " %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g",
		  varfinal[l]
		  , (varfinal[l] - varinit[l])
		  // ,totloss_full[l]
		  , realloss, floors_full[l]
		  , inflows_full[l]
		  , radiations_full[l]
		  // ,(varfinal[l]+totloss_full[l]+radiations_full[l]-floors_full[l]-inflows_full[l]-varinit[l]));
		  ,
		  (varfinal[l] + realloss + radiations_full[l] -
		   floors_full[l] - inflows_full[l] - varinit[l]));
	  if (i == 2)
	    fprintf(efboth, "\n");
	}
	fprintf(efboth, "\n");
	fflush(efboth);

      }				// end over ener or final output of
      // integrated variables 
    }				// end if cpu to write to file
    // not combining dumps!
    /* dump when dt is too low, at first, or during regular intervals,
       or at end */
    enerc++;
    tener = tstart + (FTYPE) (enerc) * DTener;
  }				// end if doing integrated variable
  // output

  // dump
  if (call_code >= 0) {
    if (((dt < DTLOWDUMP) && (CHECKDTLOW == 1)) || (t >= tpdump)
	|| (call_code == 2)) {
      if (PDUMPFLAG == 1) {
	// the -1 means 0, but different other type
	dump(NULL, pdump_cnt, 0, -1);	// primitive variable
	// dumps
	pdump_cnt++;

	if (myid <= 0) {
	  sprintf(dfnam, "%s0_numpdumps%s", DATADIR, DATEXT);
	  if ((dump_file = fopen(dfnam, "w")) == NULL) {
	    fprintf(fail_file,
		    "error opening dump output file %s\n", dfnam);
	    myexit(1);
	  }
	  fprintf(dump_file, "#%30s\n", "Number of pdumps");
	  fprintf(dump_file, "%30d\n", pdump_cnt);
	  fclose(dump_file);
	}
      }
      // time-based period
      pdumpc++;
      tpdump = tstart + (FTYPE) (pdumpc) * DTpd;
    }
    if (((dt < DTLOWDUMP) && (CHECKDTLOW == 1)) || (t >= tdump)
	|| (call_code == 2)) {

      fprintf(stderr, "\nproc: %2d dump: %5d, cc: %5d\n", myid,
	      dumpc, call_code);

      if (DUMPFLAG == 1) {
	dump(NULL, dump_cnt, 0, 0);	// normal variable dumps
	dump_cnt++;
      }
      if (NPDUMPFLAG == 1) {
	dump(NULL, npdump_cnt, 3, 0);	// non-primitive variable
	// dumps(those things
	// complicated to compute)
	npdump_cnt++;
      }
      if (FLOORDUMPFLAG == 1) {
	dump(NULL, floordump_cnt, 2, 0);	// floor dump
	floordump_cnt++;
      }
      if (analoutput > 0) {
	if ((ADUMPFLAG == 1)
	    || ((call_code == 0) && (ADUMPFLAG == -1))) {
	  dump(NULL, adump_cnt, 1, 0);
	  adump_cnt++;
	}
      }
      // now write number of dumps out to file for pp control
      if (myid <= 0) {
	// normal dump
	sprintf(dfnam, "%s0_numdumps%s", DATADIR, DATEXT);
	if ((dump_file = fopen(dfnam, "w")) == NULL) {
	  fprintf(fail_file,
		  "error opening dump output file %s\n", dfnam);
	  myexit(1);
	}
	fprintf(dump_file, "#%30s\n", "Number of dumps");
	fprintf(dump_file, "%30d\n", dump_cnt);
	fclose(dump_file);

	// np dumps
	sprintf(dfnam, "%s0_numnpdumps%s", DATADIR, DATEXT);
	if ((dump_file = fopen(dfnam, "w")) == NULL) {
	  fprintf(fail_file,
		  "error opening dump output file %s\n", dfnam);
	  myexit(1);
	}
	fprintf(dump_file, "#%30s\n", "Number of npdumps");
	fprintf(dump_file, "%30d\n", npdump_cnt);
	fclose(dump_file);

	// floor dumps
	sprintf(dfnam, "%s0_numfloordumps%s", DATADIR, DATEXT);
	if ((dump_file = fopen(dfnam, "w")) == NULL) {
	  fprintf(fail_file,
		  "error opening dump output file %s\n", dfnam);
	  myexit(1);
	}
	fprintf(dump_file, "#%30s\n", "Number of floordumps");
	fprintf(dump_file, "%30d\n", floordump_cnt);
	fclose(dump_file);

	// adump
	sprintf(dfnam, "%s0_numadumps%s", DATADIR, DATEXT);
	if ((dump_file = fopen(dfnam, "w")) == NULL) {
	  fprintf(fail_file,
		  "error opening dump output file %s\n", dfnam);
	  myexit(1);
	}
	fprintf(dump_file, "#%30s\n", "Number of adumps");
	fprintf(dump_file, "%30d\n", adump_cnt);
	fclose(dump_file);

      }
      // time-based period
      dumpc++;
      tdump = tstart + (FTYPE) (dumpc) * DTd;
    }

    /* make image of variables at regular intervals */
    if ((t >= timage) || (call_code == 2)) {
      fprintf(stderr, "\nproc: %2d image: %5d, cc: %5d\n", myid,
	      imagec, call_code);
      if (IMAGEFLAG == 1) {
	image(im_cnt, -1, -1, call_code, 0);	// set which
	// images to
	// produce here
	im_cnt++;

	if (myid <= 0) {
	  sprintf(temps, "%s%s", DATADIR, "i/");

	  sprintf(dfnam, "%s0_numimages%s", temps, DATEXT);
	  if ((image_file = fopen(dfnam, "w")) == NULL) {
	    fprintf(fail_file,
		    "error opening dump output file %s\n", dfnam);
	    myexit(1);
	  }
	  fprintf(image_file, "#%30s\n", "Number of images");
	  fprintf(image_file, "%30d\n", im_cnt);
	  fclose(image_file);
	}
      }
      // time-based period
      imagec++;
      timage = tstart + (FTYPE) (imagec) * DTi;
    }
    // cleanup if final call
    if (call_code == 2) {
      if (myid <= 0) {
	fclose(ener_file);
	fclose(final_output);

	if (DOLOSSDIAG) {
	  fclose(loss_file);
	  if (DETAILMLOSS == 2) {
	    fclose(loss_file2);
	  }
	}
      }				// end over write cpu
    }				// end if final call
  }				// end if call_code>=0

  if (call_code >= 0) {
    firsttime0 = 0;
  }
}



// if change any headers, change init.c's input stuff!!

// which sign used to determine sample type, and so size/sample/zonec
// output
// currently used by all except rwhich==11 (avg1d)
void dump_header(FILE * fp, int rwhich)
{
  int realsample, realzone;
  int which;


  if (rwhich < 0) {
    which = -rwhich - 1;
    realsample = 0;
    realzone = 0;
  } else {
    which = rwhich;
    realsample = 0;
    realzone = DUMPSM;
  }
  // version header
  if ((which == 0) || (which == 1)) {
    fprintf(fp, "#%10s\n%10d %10d\n", "DVER", DVER, DTYPE);
  } else if (which == 2) {
    fprintf(fp, "#%10s\n%10d %10d\n", "FLVER", FLVER, FLTYPE);
  } else if (which == 3) {
    fprintf(fp, "#%10s\n%10d %10d\n", "NPVER", NPVER, NPTYPE);
  }
  if (which >= 0) {
    fprintf(fp, "#%21s %6s %6s\n", "t", "SAMPLE", "ZONEC");
    fprintf(fp, " %21.15g %6d %6d\n", t, realsample, realzone);
    fprintf(fp, "#%6s %6s %6s\n", "N1", "N2", "N3");
    if (rwhich == -1) {
      if (FULLOUTPUT == 0) {
	fprintf(fp, " %6d %6d %6d\n", N1, N2, N3);
      } else {
	fprintf(fp, " %6d %6d %6d\n", N1M, N2M, N3M);
      }
    } else {
      fprintf(fp, " %6d %6d %6d\n", DUMN1, DUMN2, DUMN3);
    }
  }
  if ((which == 10) || (which == 11)) {	// extra needed info
    fprintf(fp, "# Averagefromto %21.15g %21.15g\n", tavgstart,
	    tavgfinal);
    fprintf(fp, "# SAMPLENUM %d\n", avgcount);
  }

}

// which sign used to determine sample type, and so size/sample/zonec
// output
// currently only used for avg1d (rwhich==11)
void dump_header_full(FILE * fp, int rwhich)
{
  int realsample, realzone;
  int which;

  if (rwhich < 0) {
    which = -rwhich - 1;
    realsample = 0;
    realzone = 0;
  } else {
    which = rwhich;
    realsample = 0;
    realzone = DUMPSM;
  }
  // version header
  if ((which == 0) || (which == 1)) {
    fprintf(fp, "#%10s\n%10d %10d\n", "DVER", DVER, DTYPE);
  } else if (which == 2) {
    fprintf(fp, "#%10s\n%10d %10d\n", "FLVER", FLVER, FLTYPE);
  } else if (which == 3) {
    fprintf(fp, "#%10s\n%10d %10d\n", "NPVER", NPVER, NPTYPE);
  }
  if (which >= 0) {
    fprintf(fp, "#%21s %6s %6s\n", "t", "SAMPLE", "ZONEC");
    fprintf(fp, " %21.15g %6d %6d\n", t, realsample, realzone);
    fprintf(fp, "#%6s %6s %6s\n", "N1", "N2", "N3");
    if (rwhich == -1) {
      if (FULLOUTPUT == 0) {
	fprintf(fp, " %6d %6d %6d\n", totalsize[1], totalsize[2],
		totalsize[3]);
      } else {
	fprintf(fp, " %6d %6d %6d\n", totalsize[1] + 2 * N1BND,
		totalsize[2] + 2 * N2BND, totalsize[3] + 2 * N3BND);
      }
    } else {
      fprintf(fp, " %6d %6d %6d\n", itotalsize[1], itotalsize[2],
	      itotalsize[3]);
    }
  }
  if (which == 11) {		// extra needed info
    fprintf(fp, "# Averagefromto %21.15g %21.15g\n", tavgstart,
	    tavgfinal);
    fprintf(fp, "# SAMPLETCNT %d SAMPLE1DNUM %d %d\n", avgcount,
	    num1d_31, num1d_32);
  }
}


void dump_header2(FILE * fp, int rwhich)
{
  int realsample, realzone;
  int which;


  if (rwhich < 0) {
    which = -rwhich - 1;
    realsample = 0;
    realzone = 0;
  } else {
    which = rwhich;
    realsample = 0;
    realzone = DUMPSM;
  }


  if ((which == 0) || (which == 1) || (which == 2)) {
    fprintf(fp, "#%21s %21s %21s ", "rho", "u", "pot");
  }
  if (which == 2) {		// add in ke for floor
    fprintf(fp, "%21s ", "ke");
  }
  if ((which == 0) || (which == 1) || (which == 2)) {
    fprintf(fp, "%21s %21s %21s\n", "vx1", "vx2", "vx3");
  }
  if (which == 3) {
    fprintf(fp, "#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n",
	    "sigma11", "sigma12", "sigma13", "sigma21", "sigma22",
	    "sigma23", "sigma31", "sigma32", "sigma33", "nuvisc");
  }
  if (which == 10) {
    fprintf(fp, "#%15s %15s %15s %15s %15s %15s %15s %15s\n", "2D_rho",
	    "2D_en", "2D_Be", "2D_cs2", "2D_Exp[S]", "2D_vx1",
	    "2D_vx2", "2D_vx3");
  }
  if (which == 11) {
    fprintf(fp,
	    "#%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
	    "1D1_rho", "1D1_en", "1D1_Be", "1D1_cs2", "1D1_S",
	    "1D1_absvx1", "1D1_absvx2", "1D1_absvx3", "1D2_rho",
	    "1D2_en", "1D2_Be", "1D2_cs2", "1D2_S", "1D2_absvx1",
	    "1D2_absvx2", "1D2_absvx3");
  }

}

// which: -2=floor on native grid -1=primitive on native grid
// 0=primitive
// 1=analytic 2=floor 3=non primitive 9=nonprim data(special post
// process) 
// 10=average dump 8=average dump(post process)
// outtype: zoom or not for interp
void dump(FILE * fp, int dump_cnt, int which, int outtype)
{
  static FTYPE(*iq)[INTN2][INTN1];
  static FTYPE(*viq)[3][INTN2][INTN1];
  int i, j, k, l, m;
  FTYPE tempf = 0;

  FTYPE(*scain)[N3M][N2M][N1M];
  FTYPE(*vecin)[3][N3M][N2M][N1M];
  char dfheader[MAXFILENAME];
  char temps[MAXFILENAME];
  char dfnam[MAXFILENAME];

  static FTYPE(*sigma)[3][N3M][N2M][N1M];	// -2*rho*nu*e_{ij}=
  // sigma_{ij}
  static FTYPE(*rost)[3][N3M][N2M][N1M];	// e_{ij}
  static FTYPE(*rostnu)[3][N3M][N2M][N1M];	// e_{ij}*nu
  static FTYPE(*nurho_real)[N2M][N1M];

  static FTYPE(*delv)[N2M][N1M];	// deldotv

  // to speed up things use matrix for const grids
  // assumes never interpolate when postproc=0
  int fpnull;
  int nogridchoice;
  int it;
  int gopp;


  // grab work space for interp(note workspaced used below too)
  iq = workiq;
  viq = workviq;

  // use outtype sign as indication of what grid version of data
  // incoming
  if (outtype < 0) {
    outtype = -outtype - 1;
    nogridchoice = 1;
  } else
    nogridchoice = 0;

  // if post process call do special things
  if (which == 9) {
    gopp = 1;
    which = 3;
  } else
    gopp = 0;

  if (fp == NULL) {
    fpnull = 1;

    if (nogridchoice == 0) {
      if (which == 0) {
	strcpy(dfheader, "dump");
      }
      if (which == 1) {
	strcpy(dfheader, "adump");
      }
      if (which == 2) {
	strcpy(dfheader, "floor");
      }
      if (which == 3) {
	strcpy(dfheader, "npdump");
      }
    } else {
      if (which == 0) {
	strcpy(dfheader, "pdump");
      }
      if (which == 1) {
	strcpy(dfheader, "apdump");
      }
      if (which == 2) {
	strcpy(dfheader, "pfloor");
      }
      if (which == 3) {
	strcpy(dfheader, "nppdump");
      }

    }
  } else
    fpnull = 0;

  if (nogridchoice == 1)
    it = -which - 1;
  else
    it = which;

  dump_header(fp, which);
  dump_header2(fp, which);

  scain = NULL;
  vecin = NULL;
  if (which == 0) {
    scain = s;
    vecin = v;
  } else if (which == 1) {
    scain = sanal;
    vecin = vanal;
  } else if (which == 2) {
    scain = floorvars;
    vecin = floorvarv;
    // floorvar0 directly
  } else if (which == 3) {
    // compute sigma
    sigma = workt1;		// assume workt1 holds sigma data in
    // all
    // cases
    if (gopp == 0) {
      rost = workt2;
      rostnu = workt3;
      delv = work1;
      nurho_real = work2;
      if (visc_real == 1) {
	compute_sigma(sigma, rost, rostnu, nurho_real, delv);
      }
    }
  } else {
    fprintf(fail_file, "not setup for which=%d in dump\n", which);
    myexit(1);
  }


  if (nogridchoice == 1) {
#if(FULLOUTPUT)
    LOOPF
#else
    LOOP
#endif
    {
      fprintf(fp, " ");
      if ((which == 0) || (which == 1) || (which == 2)) {
	for (l = 1; l <= NUMSCA; l++) {
	  fprintf(fp, "%21.15g ", scain[l][k][j][i]);
	}

	for (l = 1; l <= NUMVEC - 1; l++) {
	  if (which == 2) {
	    // already in middle of zone
	    fprintf(fp, "%21.15g ", floorvar0[l][k][j][i]);
	  }
	  for (m = 1; m <= 3; m++) {
	    if ((DUMPSM == 1) && (nogridchoice == 0)) {
	      if (m == 1) {
		if (i == N1 + N1BND - 1) {
		  tempf = vecin[l][m][k][j][i];
		} else
		  tempf = e2z_1(vecin[l][m][k], j, i);
	      } else if (m == 2) {
		if (j == N2 + N2BND - 1) {
		  tempf = vecin[l][m][k][j][i];
		} else
		  tempf = e2z_2(vecin[l][m][k], j, i);
	      } else if (m == 3) {
		if (k == N3 + N3BND - 1) {
		  tempf = vecin[l][m][k][j][i];
		} else
		  tempf = vecin[l][m][k][j][i];
	      }
	    } else
	      tempf = vecin[l][m][k][j][i];
	    fprintf(fp, "%21.15g ", tempf);
	  }
	}
      } else if (which == 3) {
	for (l = 1; l <= 3; l++) {
	  for (m = 1; m <= 3; m++) {
	    if ((DUMPSM == 1) && (nogridchoice == 0)) {
	      if (((l == 1) && (m == 3))
		  || ((l == 3) && (m == 1))) {
		if (i == N1 + N1BND - 1) {
		  tempf = sigma[l][m][k][j][i];
		} else
		  tempf = e2z_1(sigma[l][m][k], j, i);
	      } else if (((l == 2) && (m == 3))
			 || ((l == 3) && (m == 2))) {
		if (j == N2 + N2BND - 1) {
		  tempf = sigma[l][m][k][j][i];
		} else
		  tempf = e2z_2(sigma[l][m][k], j, i);
	      } else if (((l == 2) && (m == 1))
			 || ((l == 1) && (m == 2))) {
		if ((i == N1 + N1BND - 1)
		    || (j == N2 + N2BND - 1)) {
		  tempf = sigma[l][m][k][j][i];
		} else
		  tempf = c2z(sigma[l][m][k], j, i);
	      } else if (((l == 1) && (m == 1))
			 || ((l == 2) && (m == 2))
			 || ((l == 3) && (m == 3))) {
		tempf = sigma[l][m][k][j][i];
	      }
	    } else
	      tempf = sigma[l][m][k][j][i];
	    fprintf(fp, "%21.15g ", tempf);
	  }			// over m
	}			// over l
	// now output nuvisc
	fprintf(fp, "%21.15g ", nu_real[k][j][i]);
      }				// if sigma dump
      fprintf(fp, "\n");
    }
  }

  fflush(fp);			// always flush
  if (fpnull == 1) {
    fclose(fp);			// otherwise assume caller will close
    // it
    fp = NULL;			// return as null since started
    // null(fclose should do this)
  }
}







// mapping from function to image space in value(color)
#define FMAPSCA1(x) log10(x)
#define FMAPSCA2(x) log10(x)
#define FMAPSCA3(x) (x)
#define FMAPSCAGEN(x) (x)

#define FMAPVEC1(x) (x)
#define FMAPVEC2(x) (x)

#define DYNAMICMM 3
// 0 : use static values as defined in initial array in image func(see
// analsol.c/defs.h/init.c)
// 1 : when creating image, use first min/max values as basis for all
// sequence of images
// 2 : when creating image, use latest min/max values as basis for this
// image
// 3 : when creating image, use latest *peak* min/max values as basis
// for
// this image, upto predetermined time: timagescale

#define ERASETILLT 1
// 1: erase latest min/max until the given time, so don't use
// crazy-high
// values that exist in early transient evolution
#define ERASETILLFACTOR 5.0	// timagescale/thisfactor is time at
																// which
				// we stop erasing

#define TOTALMM 1
// 1: find total time-max and time-min of all values and output to file
// when t=tf

void image(int im_cnt, int wsca, int wvec, int call_code, int outtype)
{
  // CUT PASTE START 0
  char ifheader[MAXFILENAME];
  FTYPE liq, a = 0, b = 0, c, rholiq;
  int iii, dualgo;
  int i, j, l, m, ll;
  int q;
  FTYPE ftempfix = 0;
  int floop;			// tells loop when first entered

  static FTYPE(*iq)[INTN2][INTN1];
  static FTYPE(*viq)[3][INTN2][INTN1];
  FILE *im_file;
  static FILE *ipartot;
  static FILE *iparuse;
  static char ifnam[MAXFILENAME], temp[MAXFILENAME];

  static int firstfirsttime = 1;
  static int lastlasttime = 0;
  static int firsttimes[ITYPES][NUMSCA + 1];
  static int firsttimev[ITYPES][NUMVEC + 1];
  // static int pal[3][256];
  static unsigned char pal[3][256];
  FILE *pal_file;
  FILE *size_file;

  int im_cnts[ITYPES][NUMSCA + 1];
  int im_cntv[ITYPES][NUMVEC + 1];

  FTYPE ftemp[2];

  char temps[MAXFILENAME];
  char command[MAXFILENAME];

  int qstart;
  static int dynamicmm3outs[ITYPES][CTYPES][NUMSCA][2];
  static int dynamicmm3outv[ITYPES][CTYPES][NUMVEC][3 + 1][2];


  // ////// BEGIN
  // 
  if (call_code == 2)
    lastlasttime++;

  // CUT PASTE END 0

  if (outtype > (ITYPES - 1)) {
    fprintf(fail_file, "outtype: %d > number of types allocated\n",
	    outtype);
    myexit(1);
  }
  // assign work space
  iq = workiq;
  viq = workviq;

  // CUT PASTE START 1

  // for now assume all im_cnt same
  for (i = 1; i <= NUMSCA; i++) {
    im_cnts[outtype][i] = im_cnt;
  }
  for (i = 1; i <= NUMVEC; i++) {
    im_cntv[outtype][i] = im_cnt;
  }


  // open image parameter file
  if (firstfirsttime) {

    for (i = 0; i < ITYPES; i++) {
      for (ll = 1; ll <= NUMSCA; ll++) {
	for (iii = 0; iii < CTYPES; iii++) {
	  for (j = 0; j <= 1; j++) {
	    dynamicmm3outs[i][iii][ll][j] = 0;
	  }
	}
      }
    }
    for (q = 0; q <= 3; q++) {
      for (i = 0; i < ITYPES; i++) {
	for (ll = 1; ll <= NUMVEC; ll++) {
	  for (iii = 0; iii < CTYPES; iii++) {
	    for (j = 0; j <= 1; j++) {
	      dynamicmm3outv[i][iii][ll][q][j] = 0;
	    }
	  }
	}
      }
    }

    for (i = 0; i <= 1; i++) {
      for (j = 1; j <= NUMSCA; j++) {
	firsttimes[i][j] = 1;
      }
      for (j = 1; j <= NUMVEC; j++) {
	firsttimev[i][j] = 1;
      }
    }
    if (myid <= 0) {
      // setup file output
      sprintf(temps, "%s%s", DATADIR, "i/");

      strcpy(ifnam, "");
      sprintf(temps, "%s%simsize.out", temps, ifnam);
      if ((size_file = fopen(temps, "wt")) == NULL) {
	fprintf(fail_file, "cannot open: %s\n", temps);
	myexit(1);
      } else {
	fprintf(size_file, "%5d %5d\n", IMGN1, IMGN2);
	fclose(size_file);
      }
    }
    // CUT PASTE HOLD 1 BEGIN
    if (myid <= 0) {

      sprintf(temps, "%s%s", DATADIR, "i/");

      for (i = 0; i < ITYPES; i++) {	// outtype
	for (ll = 1; ll <= NUMSCA; ll++) {
	  for (iii = 0; iii < CTYPES; iii++) {
	    sprintf(ifnam, "%simx%01d-%01d-s%01d/", temps, i, iii, ll);
	    sprintf(command, "mkdir %s", ifnam);
	    system(command);
	    if (deleteoldimg) {
	      sprintf(command, "rm %s*.dat*", ifnam);
	      system(command);
	    }
	  }
	}
      }
      for (i = 0; i < ITYPES; i++) {
	for (ll = 1; ll <= NUMVEC - 1; ll++) {
	  for (iii = 0; iii < CTYPES; iii++) {
	    for (q = 0; q <= 3; q++) {
	      sprintf(ifnam, "%simx%01d-%01d-v%01d-%01d/",
		      temps, i, iii, ll, q);
	      sprintf(command, "mkdir %s", ifnam);
	      system(command);
	      if (deleteoldimg) {
		sprintf(command, "rm %s*.dat*", ifnam);
		system(command);
	      }
	    }
	  }
	}
      }


    }
    // CUT PASTE HOLD 1 END 

    if (myid <= 0) {
      sprintf(temps, "%s%s", DATADIR, "i/");

      sprintf(ifnam, "%s0_ipartot%s", temps, PAREXT);
      if ((ipartot = fopen(ifnam, "w")) == NULL) {
	fprintf(fail_file,
		"error opening ipartot output file %s\n", ifnam);
	myexit(1);
      }

      sprintf(temps, "%s%s", DATADIR, "i/");

      sprintf(ifnam, "%s0_iparuse%s", temps, PAREXT);
      if ((iparuse = fopen(ifnam, "w")) == NULL) {
	fprintf(fail_file,
		"error opening iparuse output file %s\n", ifnam);
	myexit(1);
      }
      // read in, but make sure no overlap on read(i.e. read
      // sequentially)

      sprintf(temps, "%s%s", DATADIR, "i/");

      sprintf(ifnam, "%sgenimages%s", temps, ".pal");
      if ((pal_file = fopen(ifnam, "rb")) == NULL) {
	fprintf(fail_file, "error opening %s\n", ifnam);
	myexit(2);
      }
      for (j = 0; j < 3; j++) {
	for (i = 0; i < 256; i++) {
	  fread(&pal[j][i], sizeof(unsigned char), 1, pal_file);
	  // fprintf(stderr,"%d %d %d\n",j,i,pal[j][i]);
	}
      }
      fclose(pal_file);
    }
    if (numprocs > 1) {
    }
  }				// end if firstfirsttime
  // CUT PASTE END 1

  // ///////// SCALARS
  if (wsca != 0) {
    for (l = 1; l <= NUMSCA; l++) {
      /* if not to do all, pick */
      if (wsca != -1) {
	ll = wsca;
      } else
	ll = l;


      // write current min/max to temp values
      if (((DYNAMICMM == 1) && (firsttimes[outtype][ll] == 1))
	  || (TOTALMM == 1) || (DYNAMICMM == 2) || ((DYNAMICMM == 3)
						    && (t <
							timagescale))) {
	if (1) {
	  if ((firsttimes[outtype][ll] == 1)
	      || (ERASETILLT && (t < timagescale / ERASETILLFACTOR))) {
	    for (m = 0; m < CTYPES; m++) {
	      mmst[outtype][m][ll][0] = s[ll][0][IMGN2 - 1][0];
	      mmst[outtype][m][ll][1] = s[ll][0][IMGN2 - 1][0];
	    }
	  }
	  LOOPI {
	    for (m = 0; m < CTYPES; m++) {
	      if (mmst[outtype][m][ll][0] > s[ll][0][j][i])
		mmst[outtype][m][ll][0] = s[ll][0][j][i];
	      if (mmst[outtype][m][ll][1] < s[ll][0][j][i])
		mmst[outtype][m][ll][1] = s[ll][0][j][i];
	    }
	  }
	}
	// set used values from current/temp values
	if (DYNAMICMM > 0) {
	  if ((DYNAMICMM == 2)
	      || ((firsttimes[outtype][ll] == 1)
		  && (DYNAMICMM == 1)) || ((DYNAMICMM == 3)
					   && (t < timagescale))) {
	    for (m = 0; m < CTYPES; m++) {
	      mms[outtype][m][ll][0] = mmst[outtype][m][ll][0];
	      mms[outtype][m][ll][1] = mmst[outtype][m][ll][1];
	    }
	    if (numprocs > 1) {
	    }
	  }
	}
	if (myid <= 0) {
	  if (TOTALMM == 1) {
	    if ((call_code == 2) && (lastlasttime == 1)) {
	      for (i = 0; i <= 1; i++) {
		for (m = 0; m < CTYPES; m++) {
		  fprintf(ipartot,
			  "mms[%d][%d][%d][%d]= %15.10g ;\n",
			  outtype, m, ll, i, mmst[outtype][m][ll][i]);
		}
	      }
	    }
	  }
	  if (DYNAMICMM == 3) {
	    if (((call_code == 2) && (lastlasttime == 1))
		|| (t >= timagescale)) {
	      for (i = 0; i <= 1; i++) {	// min and max
		for (m = 0; m < CTYPES; m++) {
		  if (dynamicmm3outs[outtype][m][ll][i]
		      == 0) {
		    fprintf(iparuse,
			    "mms[%d][%d][%d][%d]= %15.10g ;\n",
			    outtype, m, ll, i, mms[outtype][m][ll][i]);
		    dynamicmm3outs[outtype][m][ll][i] = 1;
		  }
		}
	      }
	    }
	  }
	}
	firsttimes[outtype][ll] = 0;	// done first time with
	// this scalar
      }
      dualgo = 1;		// force so can use linear for linear
      // and
      // log for log, and linear+log to
      // recapture data stream at both levels

      for (iii = 0; iii <= dualgo; iii++) {
	// setup file output
	sprintf(temps, "%s%s", DATADIR, "i/");
	strcpy(ifheader, "imx");
	sprintf(temps, "%simx%01d-%01d-s%01d/", temps, outtype,
		iii, ll);
	sprintf(ifnam, "%s%s%01d-%01d-s%01d-%04d%s%s", temps,
		ifheader, outtype, iii, ll, im_cnts[outtype][ll],
		DATEXT, myidtxt);


	if (IMAGEFORMAT == 0) {
	  strcat(ifnam, ".r8");
	  if (GZIPIMAGE != 3)
	    im_file = fopen(ifnam, "wb");
	  else {
	    sprintf(temp, "gzip > %s.gz", ifnam);
	    strcpy(ifnam, temp);	// for below fprintf
	    im_file = popen(ifnam, "w");
	  }
	  if (im_file == NULL) {
	    fprintf(fail_file, "error opening image file: %s\n", ifnam);
	    myexit(2);
	  }
	  fprintf(im_file, "RAW\n");
	  SCAHEADER(im_file, ll, outtype, iii,
		    im_cnts[outtype][ll], t, IMGN1, IMGN2);

	}

	if (IMAGEFORMAT == 1) {
	  strcat(ifnam, ".ppm");
	  if (GZIPIMAGE != 3)
	    im_file = fopen(ifnam, "wt");
	  else {
	    sprintf(temp, "gzip > %s.gz", ifnam);
	    strcpy(ifnam, temp);	// for below fprintf
	    im_file = popen(ifnam, "w");
	  }
	  if (im_file == NULL) {
	    fprintf(fail_file, "error opening image file: %s\n", ifnam);
	    myexit(2);
	  }
	  fprintf(im_file, "P6\n");
	  SCAHEADER(im_file, ll, outtype, iii,
		    im_cnts[outtype][ll], t, IMGN1, IMGN2);
	  if (GZIPIMAGE != 3) {
	    fclose(im_file);
	    // reopen in binary append mode
	    if ((im_file = fopen(ifnam, "ab")) == NULL) {
	      fprintf(fail_file,
		      "error opening image file binary append: %s\n",
		      ifnam);
	      myexit(2);
	    }
	  }
	}
	// ///////////////////////////////////////////////////////////////////
	// now set image map using min/max
	if (ll == 1) {
	  if (iii == 0) {
	    b = FMAPSCA1(mms[outtype][iii][ll][0]);
	    if (fabs(b - FMAPSCA1(mms[outtype][iii][ll][1])) > 1E-10) {
	      a = 255. / (FMAPSCA1(mms[outtype][iii][ll][1]) - b);
	    } else {
	      a = 0;		// since then undefined
	    }
	  }
	  if (iii == 1) {
	    b = FMAPSCAGEN(mms[outtype][iii][ll][0]);
	    if (fabs(b - FMAPSCAGEN(mms[outtype][iii][ll][1]))
		> 1E-10) {
	      a = 255. / (FMAPSCAGEN(mms[outtype][iii][ll][1]) - b);
	    } else {
	      a = 0;		// since then undefined
	    }
	  }
	} else if (ll == 2) {
	  if (iii == 0) {
	    b = FMAPSCA2(mms[outtype][iii][ll][0]);
	    if (fabs(b - FMAPSCA2(mms[outtype][iii][ll][1])) > 1E-10) {
	      a = 255. / (FMAPSCA2(mms[outtype][iii][ll][1]) - b);
	    } else {
	      a = 0;		// since then undefined
	    }
	  }
	  if (iii == 1) {
	    b = FMAPSCAGEN(mms[outtype][iii][ll][0]);
	    if (fabs(b - FMAPSCAGEN(mms[outtype][iii][ll][1]))
		> 1E-10) {
	      a = 255. / (FMAPSCAGEN(mms[outtype][iii][ll][1]) - b);
	    } else {
	      a = 0;		// since then undefined
	    }
	  }
	} else if (ll == 3) {
	  b = FMAPSCA3(mms[outtype][iii][ll][0]);
	  if (fabs(b - FMAPSCA3(mms[outtype][iii][ll][1])) > 1E-10) {
	    a = 255. / (FMAPSCA3(mms[outtype][iii][ll][1]) - b);
	  } else {
	    a = 0;		// since then undefined
	  }
	}
	c = -a * b;

	// output image using map
	LOOPI {
	  if (1) {
	    liq = s[ll][0][j][i];
	  }
	  if (ll == 1) {
	    if (iii == 0) {
	      ftempfix = FMAPSCA1(liq);
	    }
	    if (iii == 1) {
	      ftempfix = FMAPSCAGEN(liq);
	    }
	  } else if (ll == 2) {
	    if (iii == 0) {
	      ftempfix = FMAPSCA2(liq);
	    }
	    if (iii == 1) {
	      ftempfix = FMAPSCAGEN(liq);
	    }
	  } else if (ll == 3) {
	    ftempfix = FMAPSCA3(liq);
	  }

	  if (ftempfix >= b) {
	    liq = a * ftempfix + c;
	  } else
	    liq = 0.0;

	  if (liq > 255.)
	    liq = 255.;
	  if (liq < 0.)
	    liq = 0.;

	  if (IMAGEFORMAT == 0) {
	    fputc((int) liq, im_file);	/* write value */
	  }
	  if (IMAGEFORMAT == 1) {
	    fputc((int) pal[0][(int) liq], im_file);	/* write red */
	    fputc((int) pal[1][(int) liq], im_file);	/* write green 
							 */
	    fputc((int) pal[2][(int) liq], im_file);	/* write blue */
	  }
	}


	if (GZIPIMAGE == 0) {
	  fclose(im_file);
	}
	if (GZIPIMAGE == 1) {
	  fclose(im_file);
	  strcpy(temp, "gzip ");
	  strcat(temp, ifnam);
	  system(temp);
	} else if (GZIPIMAGE == 2) {
	  fclose(im_file);
	  mysys("gzip", ifnam);
	} else if (GZIPIMAGE == 3) {
	  pclose(im_file);
	}
      }				// over dualgo
      /* cut short loop if only to do one */
      if (wsca != -1)
	l = NUMSCA;
    }
  }

  // ////////// VECTORS
  if (wvec != 0) {
    qstart = 1;

    for (l = 1; l <= NUMVEC - 1; l++) {
      /* if not to do all, pick */
      if (wvec != -1) {
	ll = wvec;
      } else
	ll = l;

      if (((DYNAMICMM == 1) && (firsttimev[outtype][ll] == 1))
	  || (TOTALMM == 1) || (DYNAMICMM == 2) || ((DYNAMICMM == 3)
						    && (t <
							timagescale))) {
	for (q = qstart; q <= 3; q++) {	// 0 is mag

	  if (q > 0) {
	    if ((firsttimev[outtype][ll] == 1)
		|| (ERASETILLT && (t < timagescale / ERASETILLFACTOR)))
	      floop = 1;
	    LOOPI {
	      for (m = 0; m < CTYPES; m++) {
		if ((mmvt[outtype][m][ll][q][0] >
		     v[ll][q][0][j][i]) || (floop == 1))
		  mmvt[outtype][m][ll][q][0] = v[ll][q][0][j][i];
		if ((mmvt[outtype][m][ll][q][1] <
		     v[ll][q][0][j][i]) || (floop == 1))
		  mmvt[outtype][m][ll][q][1] = v[ll][q][0][j][i];
	      }
	      floop = 0;
	    }
	  } else {
	    if ((firsttimev[outtype][ll] == 1)
		|| (ERASETILLT && (t < timagescale / ERASETILLFACTOR)))
	      floop = 1;
	    LOOPI {
	      ftemp[0] = 0.0;
	      for (m = 1; m <= 3; m++) {
		ftemp[0] += v[ll][m][0][j][i] * v[ll][m][0][j][i];
	      }
	      ftemp[0] = sqrt(ftemp[0]);
	      for (m = 0; m < CTYPES; m++) {
		if ((mmvt[outtype][m][ll][q][0] > ftemp[0])
		    || (floop == 1))
		  mmvt[outtype][m][ll][q][0] = ftemp[0];
		if ((mmvt[outtype][m][ll][q][1] < ftemp[0])
		    || (floop == 1))
		  mmvt[outtype][m][ll][q][1] = ftemp[0];
	      }
	      floop = 0;
	    }
	  }
	}			// over mag+dir


	if (DYNAMICMM > 0) {
	  // mmv is what's actually used.
	  if ((DYNAMICMM == 2)
	      || ((firsttimev[outtype][ll] == 1)
		  && (DYNAMICMM == 1)) || ((DYNAMICMM == 3)
					   && (t < timagescale))) {
	    for (q = qstart; q <= 3; q++) {
	      for (m = 0; m < CTYPES; m++) {
		mmv[outtype][m][ll][q][0] = mmvt[outtype][m][ll][q][0];
		mmv[outtype][m][ll][q][1] = mmvt[outtype][m][ll][q][1];
	      }
	    }
	    if (numprocs > 1) {
	    }
	  }
	}

	if (myid <= 0) {
	  // need mmvt diff from mmv so can find all-time
	  // min/max here
	  if (TOTALMM == 1) {
	    if ((call_code == 2) && (lastlasttime == 1)) {
	      for (q = qstart; q <= 3; q++) {
		for (i = 0; i < ITYPES; i++) {
		  for (m = 0; m < CTYPES; m++) {
		    fprintf(ipartot,
			    "mmv[%d][%d][%d][%d][%d]= %15.10g ;\n",
			    outtype, m, ll, q, i,
			    mmvt[outtype][m][ll][q]
			    [i]);
		  }
		}
	      }
	    }
	  }
	  if (DYNAMICMM == 3) {
	    if (((call_code == 2) && (lastlasttime == 1))
		|| (t >= timagescale)) {
	      for (q = qstart; q <= 3; q++) {
		for (i = 0; i <= 1; i++) {	// min and 
		  // max
		  for (m = 0; m < CTYPES; m++) {
		    if (dynamicmm3outv[outtype][m][ll]
			[q][i] == 0) {
		      fprintf(iparuse,
			      "mmv[%d][%d][%d][%d][%d]= %15.10g ;\n",
			      outtype, m, ll, q, i,
			      mmvt[outtype][m][ll][q]
			      [i]);
		      dynamicmm3outv[outtype][m][ll]
			  [q][i] = 1;
		    }
		  }
		}
	      }
	    }
	  }
	}
	firsttimev[outtype][ll] = 0;
      }				// if firsttime[outtype][ll]==1 or
      // dynamicmm==2 or totalmm==1


      for (q = qstart; q <= 3; q++) {	// over components





	dualgo = 0;

	for (iii = 0; iii <= dualgo; iii++) {
	  // setup file output
	  sprintf(temps, "%s%s", DATADIR, "i/");

	  strcpy(ifheader, "imx");
	  sprintf(temps, "%simx%01d-%01d-v%01d-%01d/", temps,
		  outtype, iii, ll, q);
	  sprintf(ifnam, "%s%s%01d-%01d-v%01d-%01d-%04d%s%s",
		  temps, ifheader, outtype, iii, ll, q,
		  im_cnts[outtype][ll], DATEXT, myidtxt);


	  if (IMAGEFORMAT == 0) {
	    strcat(ifnam, ".r8");
	    if (GZIPIMAGE != 3)
	      im_file = fopen(ifnam, "wb");
	    else {
	      sprintf(temp, "gzip > %s.gz", ifnam);
	      strcpy(ifnam, temp);	// for below
	      // fprintf
	      im_file = popen(ifnam, "w");
	    }
	    if (im_file == NULL) {
	      fprintf(fail_file,
		      "error opening image file: %s\n", ifnam);
	      myexit(2);
	    }
	    fprintf(im_file, "RAW\n");
	    VECHEADER(im_file, ll, outtype, iii, q,
		      im_cnts[outtype][ll], t, IMGN1, IMGN2);
	  }
	  if (IMAGEFORMAT == 1) {
	    strcat(ifnam, ".ppm");
	    if (GZIPIMAGE != 3)
	      im_file = fopen(ifnam, "wt");
	    else {
	      sprintf(temp, "gzip > %s.gz", ifnam);
	      strcpy(ifnam, temp);	// for below
	      // fprintf
	      im_file = popen(ifnam, "w");
	    }
	    if (im_file == NULL) {
	      fprintf(fail_file,
		      "error opening image file: %s\n", ifnam);
	      myexit(2);
	    }
	    fprintf(im_file, "P6\n");
	    VECHEADER(im_file, ll, outtype, iii, q,
		      im_cnts[outtype][ll], t, IMGN1, IMGN2);
	    if (GZIPIMAGE != 3) {
	      fclose(im_file);
	      // reopen in binary append mode
	      if ((im_file = fopen(ifnam, "ab")) == NULL) {
		fprintf(fail_file,
			"error opening image file binary append: %s\n",
			ifnam);
		myexit(2);
	      }
	    }
	  }
	  // setup map
	  if (ll == 1) {
	    if (iii == 0) {
	      // warning: silent fail if 0,1 are both 0
	      // for
	      // q=0 ?
	      b = FMAPVEC1(mmv[outtype][iii][ll][q][0]);
	      if (fabs
		  (b - FMAPVEC1(mmv[outtype][iii][ll][q][1])) > 1E-10) {
		a = 255. / (FMAPVEC1(mmv[outtype][iii][ll][q][1])
			    - b);
	      } else {
		a = 0;
	      }
	    }
	    if (iii == 1) {	// rho*v // only good if v
	      // range
	      // is - to + values
	      b = FMAPVEC1(mmv[outtype][iii][ll][q][0]);
	      if (fabs
		  (b - FMAPVEC1(mmv[outtype][iii][ll][q][1])) > 1E-10) {
		a = 255. / (FMAPVEC1(mmv[outtype][iii][ll][q][1])
			    - b);
	      } else {
		a = 0;
	      }
	    }
	  }
	  if (ll == 2) {
	    if (iii == 0) {	// B
	      b = FMAPVEC2(mmv[outtype][iii][ll][q][0]);
	      if (fabs
		  (b - FMAPVEC2(mmv[outtype][iii][ll][q][1])) > 1E-10) {
		a = 255. / (FMAPVEC2(mmv[outtype][iii][ll][q][1])
			    - b);
	      } else {
		a = 0;
	      }
	    }
	    if (iii == 1) {	// B*v // only good if v range
	      // is
	      // - to + values
	      b = FMAPVEC2(mmv[outtype][iii][ll][q][0]);
	      if (fabs
		  (b - FMAPVEC2(mmv[outtype][iii][ll][q][1])) > 1E-10) {
		a = 255. / (FMAPVEC2(mmv[outtype][iii][ll][q][1])
			    - b);
	      } else {
		a = 0;
	      }
	    }
	  }
	  c = -a * b;

	  LOOPI {
	    liq = 0;
	    if (q == 0) {
	      for (m = 1; m <= 3; m++) {
		if (1) {
		  liq = v[ll][m][0][j][i] * v[ll][m][0][j][i] + liq;
		  rholiq = s[1][0][j][i];
		}
	      }
	      liq = sqrt(liq);	// magnitude of vector at
	      // zone
	    } else {
	      if (1) {
		liq = v[ll][q][0][j][i];
		rholiq = s[1][0][j][i];
	      }
	    }
	    if (iii == 1) {
	      liq = rholiq * liq;
	    }


	    if (ll == 1) {
	      ftempfix = FMAPVEC1(liq);
	    } else if (ll == 2) {
	      ftempfix = FMAPVEC2(liq);
	    }
	    if (ftempfix >= b) {
	      liq = a * ftempfix + c;
	    } else
	      liq = 0.0;

	    if (liq > 255.)
	      liq = 255.;
	    if (liq < 0.)
	      liq = 0.;

	    if (IMAGEFORMAT == 0) {
	      fputc((int) liq, im_file);	/* write value */
	    }
	    if (IMAGEFORMAT == 1) {
	      fputc((int) pal[0][(int) liq], im_file);	/* write red */
	      fputc((int) pal[1][(int) liq], im_file);	/* write green 
							 */
	      fputc((int) pal[2][(int) liq], im_file);	/* write blue */
	    }
	  }
	  if (GZIPIMAGE == 0) {
	    fclose(im_file);
	  }
	  if (GZIPIMAGE == 1) {
	    fclose(im_file);
	    strcpy(temp, "gzip ");
	    strcat(temp, ifnam);
	    system(temp);
	  } else if (GZIPIMAGE == 2) {
	    fclose(im_file);
	    mysys("gzip", ifnam);
	  } else if (GZIPIMAGE == 3) {
	    pclose(im_file);
	  }
	}			// dualgo
      }				// over components of the vector
      /* cut short loop if only to do one */
      if (wvec != -1)
	l = NUMVEC - 1;
    }				// over vectors
  }				// if any vectors
  firstfirsttime = 0;

  if ((call_code == 2) && (lastlasttime == 1)) {
    // should be done now
    fclose(ipartot);		// assumes 1 outtype per run
    fclose(iparuse);		// assumes 1 outtype per run
  }

  if (call_code == 2)
    lastlasttime++;
}
