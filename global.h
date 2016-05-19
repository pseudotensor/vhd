#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>



#define N1 64			// must be even for totally general
				// consistency
#define N2 40			// must be even for totally general
				// consistency
#define N3 1

/* NBIG is bigger of N1 and N2 and N3 */
#define NBIG1 ((N1>N2) ? N1 : N2)
#define NBIG  ((NBIG1>N3) ? NBIG1 : N3)


// make sure has trailing /
#define DATADIR "./"

// size of data type used for all floats
// you must change all % f's to % lf's in global.h when using doubles
// you must change all % lf's to % f's in global.h when using floats!
#define FLOATTYPE 0		// 0: float 1: double (normal
				// non-sensitive or performance
				// critical datatypes)
#define SENSITIVE 1		// 0: float 1: double (non-perf
				// critical or sensitive data types)
// need not change below datatype stuff
#if(FLOATTYPE==0)
#define FTYPE float
#else
#define FTYPE double
#endif
#if(SENSITIVE==0)		// for sensitive counters
#define SFTYPE float
#else
#define SFTYPE double
#endif


// how often in steps to output step/dt/t data
#define NDTCCHECK 500
// how often in steps to check the go.go file to see if to continue
// running
#define NGOCHECK 500
// how often in steps to check speed in zonecycles/sec
#define NZCCHECK 1000


#define COMPDIM 2		// how many computational
				// dimensions(bounded directions)
#define MOVINGGRID 0		// whether using moving grid or not
#define DEBUG 0			// (temp stuff, none right now, moved
				// to real seperate defines)
// 0: no debug statements: real run
// 1: normal run with some checks/statements
// 2: as 1 but add output on changes to data

// force flush of fail_file data each dt
#define FLUSHFAILDT 1

// force floor on varibles
#define FORCERHO 1		// 1: force DENSITYFLOOR on rho
#define FORCEIE 1		// 1: force IEFLOOR on ie(s2 really)
#define TS0CHECK 0		// 1: check if rho or ie<0 in
				// timestep.c not needed if forcing
				// floor

// number of dt checks for timestep.c
#define NUMDTCHECKS 10


// Make sure boundary zones do not allow inflow on outflow boundary(see 
// 
// 
// 
// bound.c)
// only applies if outflow condition
#define INFLOWCHECKIX1 1	// inner x1-edge
#define INFLOWCHECKOX1 1	// outer x1-edge
#define INFLOWCHECKIX2 0	// inner x2-edge // only use in 2D+
#define INFLOWCHECKOX2 0	// outer x2-edge
#define INFLOWCHECKIX3 0	// inner x3-edge // only use in 3D+
#define INFLOWCHECKOX3 0	// outer x3-edge


// should either use 1 for all non-zero or 2 for all non-zero, don't
// mix 1 and 2.
#define N1BND 2
#define N2BND 2
#define N3BND 0
#define NBIGBND1 ((N1BND>N2BND) ? N1BND : N2BND)
#define NBIGBND  ((NBIGBND1>N3BND) ? NBIGBND1 : N3BND)

#define DOBOUNDSCA 1
#define DOBOUNDVEC 1
#define SCABCCORNER 1
#define VECBCCORNER 1

// whether to bound potential(generally should be kept as analytic
// solution)
#define NOBOUNDPOT 1		// 0: bound potential like other
				// scalars, 1: no bound of pot
#define BOUNDN1 1		// 0: don't do N1 layer bound, 1: do
// can use below if make N2==1 and don't care about theta structure.
// Assumes many things, so could be nasty, but much faster for 1-D
// problems
#if(N2==1)
#define BOUNDN2 0		// 0: don't do N2 layer bound, 1: do (0
				// for speed on 1D problems)
#else
#define BOUNDN2 1		// generally this should be 1
#endif

#define BOUNDN3 0		// 0: don't do N3 layer bound, 1: do

#define IOBOUNDARY 0		// 0: normal boundary, 1:
				// Inflow/outflow split outer
				// x2-boundary
#define NUMBTRANS 4		// number of boundary transitions(see
				// defs.h/init.c)

/* allocated memory uses this for active zones 0-N1-1 and bc beyond
   that */
#define N1M (N1+N1BND*2)
#define N2M (N2+N2BND*2)
#define N3M (N3+N3BND*2)

/* NBIGM is bigger of N1M and N2M and N3M */
#define NBIG1M ((N1M>N2M) ? N1M : N2M)
#define NBIGM  ((NBIG1M>N3M) ? NBIG1M : N3M)

#define INNERBC 0		// whether to include inner zones in
				// bound loop


#define CHECKDTLOW 1		// whether to check if dt went below
				// critical value
#define DTLOWEST (1.E-10)
#define IDTLOWEST (1.0/DTLOWEST)
#define SQIDTLOWEST (1.0/(DTLOWEST*DTLOWEST))	// actually dt over
						// cour then squared
#define DTLOWDUMP (DTLOWEST)	// lowest dt to initiate dumping


#define PDEN       1		// 0->no pressure 1->pressure from
				// density/energy
#define PGRAV      1		// 0->no ext pot 1->ext pot
#define CURVE      1		// 0->no curvature terms 1-> do


#define VISC_LINEAR 0		// 0-> linear viscosity off/1=on

// turn on/off different terms (1 or 0)
#define VISCE11 1
#define VISCE22 1
#define VISCE33 1
#define VISCE12 1
#define VISCE21 VISCE12
#define VISCE13 1
#define VISCE31 VISCE13
#define VISCE23 1
#define VISCE32 VISCE23

// 0: don't generate parameter output 1: do
#define DOPARDIAG 1		// generally do want to overwrite when
				// pp=0
// 0: don't output general diag 1: do output
#define DOGENDIAG 1
 // 0: don't output loss diag 1: do output/counting(required in
 // general)
#define DOLOSSDIAG 1
 // 0: don't output floor diag 1: do output level 1 2: do full floor
 // diag output (floor counters always active, this just controls
 // diagnostics of counters)
#define DOFLOORDIAG 1		// output details on what algorithm
				// caused most lows and how low values
				// got
#define DOFLOORD2 0		// VERY DETAILED: output where and what
				// corrected in logfl_file
 // 0: don't do dt diag 1: do
#define DODTDIAG 1		// output of dt constraints
#define DOTSTEPDIAG 1		// output of timescales(true global
				// data of just above)
#define DOSPDIAG 1		// sonic point location check log
				// output

#define DOLOGSTEP 1		// log the time step as per NDTCCHECK
#define DOLOGPERF 1		// log the performance as per NZCCHECK

#define PDUMPFLAG 1		// 0: don't dump primitive on own grid
				// for reentrance 1: do
#define DUMPFLAG 1		// 0: don't create dump files 1: do
				// create
#define NPDUMPFLAG 1		// 0: don't create np dump files 1: do
				// create
#define FLOORDUMPFLAG 1		// 0: don't create floor dump files 1:
				// do create
#define ADUMPFLAG -1		// 0: don't create analytic dump files
				// 1: do create -1: only create first
				// one

#define OLDSCHOOL 0		// 0: assume as current data sets 1:
				// use older format for data sets

#define IMAGEFLAG 1		// 0: don't create images 1: do create
// only need first data set for image dump before pp
#define NUMOUTTYPE 1		// number of outtypes to do(see diag.c
				// for how chosen)
#define OLDSCHOOL2 0		// 0: assume as current r8's 1: as old
				// r8's

// generally, should use sm to grid data correctly, but if don't care
// about reentering at all, can do here.
#define DUMPSM 1		// 0: dump data where gridded 1: dump
				// in center of zone for all
#define FULLOUTPUT 0		// 0: output data for active grid, 1:
				// full grid output
// can only take fullinput on fulloutput data
#define FULLINPUT 0		// 0: input data for active grid, 1:
				// full grid input

#define IMGN1 N1
#define IMGN2 N2
#define IMGN3 N3

#define DUMN1 N1
#define DUMN2 N2
#define DUMN3 N3


#define ITYPES 2		// number of types of image ranges, 0,
				// 1
#define CTYPES 2		// number of types of computed image
				// ranges, 0, 1

// determine largest interpolated grid sizes for memory allocation of
// working space
#define INTN1 ((DUMN1>IMGN1) ? DUMN1 : IMGN1)
#define INTN2 ((DUMN2>IMGN2) ? DUMN2 : IMGN2)
#define INTNBIG ((INTN1>INTN2) ? INTN1 : INTN2)

#define IMAGEFORMAT 0		// should always be 0 for post
				// processing ability
#define IMAGEFORMATINPUT 0	// dummy
// 0: r8(best in general except for immediate viewing)
// 1: ppm (best, esp. when used with gzip below)(can't be used during
// post process since can't reverse lookup easily, so don't use for
// now)

#define GZIPIMAGE 3		// best to always use this when post
				// proc
#define GZIPIMAGEINPUT 3	// choose was input type for post
				// processing

// 0: don't gzip(best if need not zipped)(necessary for gm over mpich)
// 1: do gzip using system call(best for compatibility)
// 2: do gzip using fork call (best for small files)
// 3: do gzip using popen (best for large files)


// 1 for cartesian coords
// 2 for cyl
// 3 for sph
// see init.c init_diffs() for usage
#define COORD 3

// keep below at 1 or 2, not 0 for all version stuff to work
#define DETAILMLOSS 1		// 0=only show total massloss each time
// 1= show mass loss total and through each boundary(4 in 2D, 8 in 3D)
// 2= Do 1 but in a seperate file output each grid points output each
// time

#define LINEARINTERP 1
// 0: just choose local value (for major speed, but large diffusion)
// 1: simple average or truely correct for uniform grid (for speed or
// uniform grid) (sorta corresponds to how differencing is done for 1st 
// 
// 
// 
// derivative in sweep.c)
// 2: correct linear average for generaly nonuniform grid (slower)

#define LINEXT 0
// 0: just copy

// major constants below
#define NUMSCA 3
#define NUMVEC 2
#define NUMSV (NUMSCA+NUMVEC)
#define NUMGRID 2
#define NUMMETRIC 4
#define NUMVOL 2
#define DIRVEC 1
#define NUMFLOORVAR 2		// 1: mass density 2: ie density
#define NUMFLOOROUT 7		// number of routines checked for floor
// 0 through 6 with 1 extra space


#define NUMLOSSVAR (NUMSCA+1+NUMVEC*3+1)
// 1: mass
// 2: enthalpy
// 3: grav pot energy

// 4: ke

// 5: s1
// 6: s2
// 7: s3
// 8: B1
// 9: B2
// 10: B3

// 11: etot(visc)

// extention for data files
#define DATEXT ".dat"
#define PAREXT ".par"
#define INEXT ".in"
#define OUTEXT ".out"
#define PPEXT ".pp"

#define INPUT2 "%f "
#define INPUT3I "%f"
#define INPUTIMGT "%f"
// INPUT1/1old/4/5/6/7/avgh1 defines in
// i/i/diag.c/diag.c/diag.c/init.c/init.c respectivly because sensitive
#define INPUT7 "%f"		// this not sensitive for now
#define INPUTRAD "%f"

#define ARGS 5
#define PARMTYPEARGS "%d %d %d %f %f"
#define MAXFILENAME 200

// file versions numbers(use sm for backwards compat)
#define PVER 6
#define GRIDVER 2
#define DVER 1			// dumps same as for pdumps, adumps
#define FLVER 2
#define NPVER 1
#define ENERVER 5
#define LOSSVER 5
#define SPVER   1
#define TSVER   1
#define LOGDTVER 1
#define STEPVER 1
#define PERFVER 3

// type designations for sm automagical read in correct format for
// similar things
#define PTYPE     1		// global par file
#define GRIDTYPE  2
#define DTYPE     3		// always same as adump and pdump
#define FLTYPE    4		// floor
#define NPTYPE    5		// np
#define ENERTYPE  8
#define LOSSTYPE  9
#define SPTYPE    10
#define TSTYPE    11
#define LOGDTTYPE 12
#define STEPTYPE  13
#define PERFTYPE  14

#define ZER 0x30
#define DIGILEN 10

// if you change these, make sure you know where they are used!!
#define ERR     1.e-6
#define SMALL	1.e-10
#define GMIN 	1.e-10
#define MIN	1.e-20		/* minimum density */
#define SSMALL  1.E-30

#define THIRD (0.3333333333333333333333333333333333333333333333333333333333333333333)

#define SCAHEADER(fp,ll,outtype,iii,cnt,t,nx,ny) \
	  fprintf(fp,"#scalar#: %d ot: %d mt: %d image#: %d t=%15.10g\n",ll,outtype,iii,im_cnts[outtype][ll],t); \
	  fprintf(fp,"%i %i\n", nx, ny); \
	  fprintf(fp,"255\n")

#define VECHEADER(fp,ll,outtype,iii,q,cnt,t,nx,ny) \
	    fprintf(fp,"#vector#: %d outtype#: %d maptype#: %d comp#: %d image#: %d t=%15.10g\n",ll,outtype,iii,q,cnt,t); \
	    fprintf(fp,"%i %i\n", nx, ny); \
	    fprintf(fp,"255\n")

#define HEADER3_S0 "%d %d\n"
#define HEADER3_S1 "%d %d %d\n"

// version of file format
// N1, N2, N3
 // L11, L12, L13, L21, L22, L23
 // rg rgp cs coolfact gam alpha_real0 n_real
 // nu_vnr nu_l nu_ten cour cour2
 // gravc massbh
 // tstart tf tavgi tavgf numavg
 // dtl dtd dti dtloss dtfloor dttimestep dtpd dtener
 // resist nu_sh
 // vg1 vg2 vg3
#define HEADER3_P  "# %4s %4s\n"\
                   "  %4d %4d\n"\
                   "# %4s %4s %4s\n"\
                   "  %4d %4d %d\n"\
                   "# %21s %21s %21s %21s %21s %21s\n"\
                   "  %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
                   "# %21s %21s %21s %21s %21s %21s %21s\n"\
                   "  %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s\n"\
		   "  %21.15g %21.15g\n"\
		   "# %21s %21s %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g %21d %21.15g\n"\
		   "# %21s %21s %21s %21s %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s\n"\
		   "  %21.15g %21.15g\n"\
		   "# %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g\n"


#define HEADER4_P  " %4d %4d %4d %4d "\
                   "%21.15g %21.15g %21.15g "\
                   "%21.15g %21.15g %21.15g "\
                   "%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "\
                   "%21.15g %21.15g %21.15g "

#define HEADER4_S  "%d %d %d %d "\
                   "%f %f %f "\
                   "%f %f %f "\
                   "%f %f %f %f %f %f %f %f %f "\
                   "%f %f %f "



// 2d only
#define deldotv(wvec,k,j,i) \
( (g[1][2][i+1]*g[1][3][i+1]*v[wvec][1][k][j][i+1]-g[1][2][i]*g[1][3][i]*v[wvec][1][k][j][i])/dvl[1][1][i]+\
(g[1][4][j+1]*v[wvec][2][k][j+1][i]-g[1][4][j]*v[wvec][2][k][j][i])/(g[2][2][i]*dvl[1][2][j]) )



#define gradv11(wvec,k,j,i) \
( (v[wvec][1][k][j][i+1]-v[wvec][1][k][j][i])/dx[1][1][i] )

#define gradv22(wvec,k,j,i) \
( (v[wvec][2][k][j+1][i]-v[wvec][2][k][j][i])/(g[2][2][i]*dx[1][2][j])+\
0.5*(v[wvec][1][k][j][i]+v[wvec][1][k][j][i+1])/(g[2][2][i])*dg[2][2][i] )

#define gradv33(wvec,k,j,i) \
(0.5*( (v[wvec][1][k][j][i]+v[wvec][1][k][j][i+1])/(g[2][3][i])*dg[2][3][i]+\
(v[wvec][2][k][j][i]+v[wvec][2][k][j+1][i])/(g[2][2][i]*g[2][4][j])*dg[2][4][j] ) )


// emulate functions as macros for speed, they aren't inlined for some
// stupid reason
// see numerics.c for details

#if(LINEARINTERP==2)
#define z2e_1(var,j,i) (0.5*(dx[1][1][i]*var[j][i-1]+dx[1][1][i-1]*var[j][i])/dx[2][1][i])
#define z2e_2(var,j,i) (0.5*(dx[1][2][j]*var[j-1][i]+dx[1][2][j-1]*var[j][i])/dx[2][2][j])
#define e2z_1(var,j,i) (0.5*(var[j][i] + var[j][i + 1]))
#define e2z_2(var,j,i) (0.5*(var[j][i] + var[j + 1][i]))
#define e2e_v2(var,j,i) (0.25* ((var[j][i] + var[1 + j][i])*dx[1][1][-1 + i] + (var[j][-1 + i] + var[1 + j][-1 + i])* dx[1][1][i])/(dx[2][1][i]))
#define e2e_v1(var,j,i) (0.25* ((var[j][i] + var[j][1 + i])* dx[1][2][-1 + j] + (var[-1 + j][i] + var[-1 + j][1 + i])* dx[1][2][j])/dx[2][2][j])
#define z2c(var,j,i) (0.25*((var[j][i]*dx[1][1][-1 + i] + var[j][-1 + i]*dx[1][1][i])*dx[1][2][-1 + j] + (var[-1 + j][i]*dx[1][1][-1 + i] + var[-1 + j][-1 + i]*dx[1][1][i])*dx[1][2][j])/(dx[2][1][i]*dx[2][2][j]))
#define c2z(var,j,i) (0.25*(var[j][i] + var[j][1 + i] + var[1 + j][i] + var[1 + j][1 + i]))

#elif(LINEARINTERP==1)

#define z2e_1(var,j,i) (0.5*(var[j][i-1]+var[j][i]))
#define z2e_2(var,j,i) (0.5*(var[j-1][i]+var[j][i]))
#define e2z_1(var,j,i) (0.5*(var[j][i] + var[j][i + 1]))
#define e2z_2(var,j,i) (0.5*(var[j][i] + var[j + 1][i]))
#define e2e_v2(var,j,i) (0.25* (var[j][i] + var[1 + j][i] + var[j][-1 + i] + var[1 + j][-1 + i] ))
#define e2e_v1(var,j,i) (0.25* (var[j][i] + var[j][1 + i] + var[-1 + j][i] + var[-1 + j][1 + i] ))
#define z2c(var,j,i) (0.25*(var[j][i] + var[j][-1 + i] + var[-1 + j][i]+ var[-1 + j][-1 + i] ))
#define c2z(var,j,i) (0.25*(var[j][i] + var[j][1 + i] + var[1 + j][i] + var[1 + j][1 + i]))

#elif(LINEARINTERP==0)

#define z2e_1(var,j,i) (var[j][i])
#define z2e_2(var,j,i) (var[j][i])
#define e2z_1(var,j,i) (var[j][i])
#define e2z_2(var,j,i) (var[j][i])
#define e2e_v2(var,j,i) (var[j][i])
#define e2e_v1(var,j,i) (var[j][i])
#define z2c(var,j,i) (var[j][i])
#define c2z(var,j,i) (var[j][i])

#endif




/* Function defines */
int main(int argc, char *argv[], char *envp[]
    );
int myexit(int call_code);

extern int init(int argc, char *argv[], char *envp[]);
extern void init_genfiles(int gopp);
extern void init_general(void);
extern void init_loss(void);
extern void init_compsave(void);
extern void init_visc(void);
extern void init_floor(void);
extern void init_inflows(void);
extern void init_radiations(void);
extern int init_mem(void);
extern int init_pointers(void);
extern int init_runpar(int gopp);
extern int init_rundat(int dump_cnt, int which);
extern int init_runimage(int im_cnt, int wsca, int wvec, int call_code,
			 int outtype);
extern int init_params(int seed, FTYPE beta, FTYPE nj);
extern int init_reentrance(void);
extern int init_reentrance2(SFTYPE time);
extern int init_otherparams(void);
extern int init_dx(int n1, int n2, int n3, int n1b, int n2b, int n3b,
		   SFTYPE * sx, int which, int outtype);
extern int init_x(int n1, int n2, int n3, int n1b, int n2b, int n3b,
		  SFTYPE * sx, int which, int outtype);
extern int init_diffs(void);
extern int init_data(void);
extern int init_bc(int simple, int ix1, int ox1, int ix2, int ox2,
		   int ix3, int ox3);
extern int init_mainbc(int px1, int six1, int rix1, int rox1, int px2,
		       int six2, int rix2, int rox2, int px3, int six3,
		       int rix3, int rox3);
extern int init_outgparm(int which);

extern FTYPE ranc(int iseed);
extern void bound(FTYPE(*vars)[N2M][N1M], FTYPE(*varv)[N3M][N2M][N1M], int wsca,	// which 
											// 
		  // 
		  // 
		  // scalars 
		  // to 
		  // update: 
		  // -1 
		  // for 
		  // all 
		  // or 
		  // 1-NUMSCA 
		  // for 
		  // individual
		  int wvec);	// which vectors to update: -1 for all
				// or 1-NUMVEC for individual, and -2
				// for any scalar/vector

extern void diag(int call_code);
extern void diagavg(int call_code);
extern void dump_header(FILE * fp, int which);
extern void dump_header_full(FILE * fp, int which);
extern void dump_header2(FILE * fp, int which);
extern void dump(FILE * fp, int dump_cnt, int which, int outtype);
extern void analdump(FILE * fp);
extern void analsolve(int gopp);
extern void sodsol(int gopp);
extern void advsol(int gopp);
extern void gausssol(int gopp);
extern void bondisol(int gopp);
extern FILE *tori1sol(int gopp);
extern void injectsol(int gopp);
extern void accountstoreset(void);
extern void visctermtest(int gopp);
extern void pulsesol(int gopp);
extern void test1sol(int gopp);
extern void test2sol(int gopp);

extern void image(int im_cnt, int wsca, int wvec, int call_code,
		  int outtype);
extern void spec_diag(void);

extern void stepvar(void);
extern void compute_sigma(FTYPE(*sigma)[3][N3M][N2M][N1M],
			  FTYPE(*rost)[3][N3M][N2M][N1M],
			  FTYPE(*rostnu)[3][N3M][N2M][N1M],
			  FTYPE(*nurho_real)[N2M][N1M],
			  FTYPE(*delv)[N2M][N1M]);
extern void timestep(void);
extern void idtcreate(FTYPE * idt, int k, int j, int i);
extern void timescale(void);
extern void timecheck(int failmode, FTYPE * idt, int k, int j, int i,
		      int reall);
extern void moc_ct(void);

extern void step_bz(void);
extern void bzsweepx(void);
extern void bzsweepy(void);
extern void step_pgc(void);
extern void step_visc(void);
extern void step_visc_real(void);
extern void tdep_compute(void);
extern void sp_compute(void);
extern void injection(void);
extern void cooling(void);
extern void compute_funcool(FTYPE * fun, FTYPE * thetai, int num,
			    FTYPE thetas, FTYPE thetaf,
			    FTYPE * dthetap);
extern void compute_funrelie(FTYPE * fun, FTYPE * thetai, int num,
			     FTYPE thetas, FTYPE thetaf,
			     FTYPE * dthetap);
extern void nu_compute(void);
extern void step_ie(void);
extern void step_relie(void);
extern void step_res(void);

extern void step_trans(void);
extern void sweepx(void);
extern void sweepy(void);
extern void dqx_calc(FTYPE(*var)[N2M][N1M], FTYPE(*dq)[N2M][N1M]);
extern void dqvx_calc(int wtype, int wcom, FTYPE(*var)[N3M][N2M][N1M],
		      FTYPE(*dqv)[N3M][N2M][N1M]);
extern void dqy_calc(FTYPE(*var)[N2M][N1M], FTYPE(*dq)[N2M][N1M]);
extern void dqvy_calc(int wtype, int wcom, FTYPE(*var)[N3M][N2M][N1M],
		      FTYPE(*dqv)[N3M][N2M][N1M]);


extern void ex_v(int dim, int dir, FTYPE(*var)[N1M], int j, int i,
		 int which);
extern void ex_s(int dim, int dir, FTYPE(*var)[N1M], int j, int i,
		 int which);
extern void ex_s_p(int dim, int dir, FTYPE(*var)[N1M], int j, int i);
extern void ex_v_p(int dim, int dir, FTYPE(*var)[N1M], int j, int i);

extern void interpolate(int wtype, int wsv, int stype, int whichg,
			int n1, int n2, FTYPE in[N2M][N1M],
			FTYPE out[INTN2][INTN1], FTYPE outerdef,
			int outtype);

extern FTYPE alnfact(FTYPE N);
extern void itoa(int x, char *p);
extern int mysys(char *com1, char *com2);
extern int fork2(void);
extern void ptraddr(int nstep);

// Macros

// note: LOOPF goes over corner zones.  I init corners to 1 or
// correctly in analsol.c so if I divide this useless data it does no
// harm.  Just for speed.  Could put ifs in to avoid corner zones, but
// slower.

// for compatibility with 1 bnd zone force no such H vs F to exist
#if(NBIGBND==2)
#define LOOPH3 for(k=-N3BND/2;k<N3+N3BND/2;k++)
#define LOOPH2 for(j=-N2BND/2;j<N2+N2BND/2;j++)
#define LOOPH1 for(i=-N1BND/2;i<N1+N1BND/2;i++)
#else
#define LOOPH3 for(k=-N3BND;k<N3+N3BND;k++)
#define LOOPH2 for(j=-N2BND;j<N2+N2BND;j++)
#define LOOPH1 for(i=-N1BND;i<N1+N1BND;i++)
#endif

#define LOOPC3 for(k=0;k<N3;k++)
#define LOOPC2 for(j=0;j<N2;j++)
#define LOOPC1 for(i=0;i<N1;i++)

#define LOOPF3 for(k=-N3BND;k<N3+N3BND;k++)
#define LOOPF2 for(j=-N2BND;j<N2+N2BND;j++)
#define LOOPF1 for(i=-N1BND;i<N1+N1BND;i++)

#define LOOPH LOOPH3 LOOPH2 LOOPH1
#define LOOPF LOOPF3 LOOPF2 LOOPF1

#define LOOPINT3 for(k=intix3;k<intox3;k++)
#define LOOPINT2 for(j=intix2;j<intox2;j++)
#define LOOPINT1 for(i=intix1;i<intox1;i++)

#define LOOPINT LOOPINT3 LOOPINT2 LOOPINT1

// boundary loops over only 1st layer of bc, bound.c takes care of rest
#if( (BOUNDN2==1)&&(BOUNDN1==1) )
#define LOOPBOUND LOOPH
#elif( (BOUNDN2==1)&&(BOUNDN1==0) )
#define LOOPBOUND k=0; i=0; LOOPH2
#elif( (BOUNDN2==0)&&(BOUNDN1==1) )
#define LOOPBOUND k=0; j=0; LOOPH1
#endif

#define LOOP LOOPC3 LOOPC2 LOOPC1

#define LOOPSK3 LOOPC3
#define LOOPSK2 for(j=skipix2;j<N2;j++)
#define LOOPSK1 for(i=skipix1;i<N1;i++)

// used for velocity calculations in step.c since when reflecting
// r/theta your metric components are 0 and floating point exception
// will occur, so just let bound do the work
// can either use below 2 loops to seperate x1/x2-directions or use
// LOOP above then have conditions using defines/ifs in loop
#define LOOPV1  LOOPC3 LOOPC2 LOOPSK1
#define LOOPV2  LOOPC3 LOOPSK2 LOOPC1
#define LOOPV3  LOOPSK3 LOOPC2 LOOPC1

#define LOOP2  LOOPC2 LOOPC1


#define LOOPT1i for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1+1;i++)

#if(NBIGBND==2)
// used for fluxes which never have i=-1 accessed
#define LOOPT0i for(k=0;k<N3;k++) for(j=skipix2-1;j<N2;j++) for(i=skipix1-1;i<N1+1;i++)
#else
#define LOOPT0i for(k=0;k<N3;k++) for(j=-1;j<N2;j++) for(i=0;i<N1+1;i++)
#endif

// used for fluxes which never have i=N1-1+1 accessed
// unless SKIPIX1==1, so that vx(0) is boundary zone so don't need
// fl(-1) for vx advection in x1-dir
#define LOOPT2i for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=skipix1-1;i<N1;i++)

// used for fluxes which never have i=-1 accessed, and don't need j=-1
// or 0 accessed for non-periodic bcs
#define LOOPT3i for(k=0;k<N3;k++) for(j=skipix2;j<N2;j++) for(i=0;i<N1+1;i++)


// used for fluxes which never have i=-1 accessed
#define LOOPT1j for(k=0;k<N3;k++) for(j=0;j<N2+1;j++) for(i=0;i<N1;i++)

#if(NBIGBND==2)
#define LOOPT0j for(k=0;k<N3;k++) for(j=-1+skipix2;j<N2+1;j++) for(i=-1;i<N1;i++)
#else
#define LOOPT0j for(k=0;k<N3;k++) for(j=0;j<N2+1;j++) for(i=-1;i<N1;i++)
#endif

// used for fluxes which never have i=N2-1+1 accessed
#define LOOPT2j for(k=0;k<N3;k++) for(j=skipix2-1;j<N2;j++) for(i=0;i<N1;i++)

// used for fluxes which never have i=-1 accessed
#define LOOPT3j for(k=0;k<N3;k++) for(j=0;j<N2+1;j++) for(i=skipix1;i<N1;i++)

// loops over more i for j and more j for i than needed, but faster to
// do all in one lump
#define LOOPVISC for(k=0;k<N3;k++) for(j=-1;j<N2;j++) for(i=-1;i<N1;i++)

// use below only if start grid in lower left corner.
// #define LOOPI for(j=IMGN2-1;j>=0;j--) for(i=0;i<IMGN1;i++)

// use for starting grid in upper left corner
#define LOOPI for(j=0;j<IMGN2;j++) for(i=0;i<IMGN1;i++)
#define LOOPINI for(j=0;j<N2;j++) for(i=0;i<N1;i++)	// used for
							// post process
#define LOOPD for(j=0;j<DUMN2;j++) for(i=0;i<DUMN1;i++)

#define LOOPINJ for(k=tagik;k<tagfk;k++) for(j=tagij;j<tagfj;j++) for(i=tagii;i<tagfi;i++)
#define LOOPFINJ for(k=t2gik;k<t2gfk;k++) for(j=t2gij;j<t2gfj;j++) for(i=t2gii;i<t2gfi;i++)
#define LOOPRINJ for(k=t3gik;k<t3gfk;k++) for(j=t3gij;j<t3gfj;j++) for(i=t3gii;i<t3gfi;i++)
#define LOOPVINJ for(k=t4gik;k<t4gfk;k++) for(j=t4gij;j<t4gfj;j++) for(i=t4gii;i<t4gfi;i++)

#define DEBUGP1 \
printf("\n\n\n\n"); \
  LOOPF{ \
printf("j: %2d i: %2d x1a: %21.15g x1b: %21.15g x2a: %21.15g x2b: %21.15g rho: %20.15g e: %21.15g pot: %21.15g Bx: %21.15g By: %21.15g Bz: %21.15g vx: %20.15g vy: %21.15g vz: %21.15g\n",j,i,x[1][1][i],x[2][1][i],x[1][2][j],x[2][2][j],s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i]); \
  }


#define BOUNDBUG printf("k: %d j: %d i: %d ll: %d kk: %d jj: %d ii: %d\n",k,j,i,ll,kk,jj,ii); printf("bcd1: %d bcd2: %d bcdir: %d\n",bcd1,bcd2,bcdir);


/* Notes: When memory is allocated for variables of size N I add the
   pointer address by NBND so the memory the pointer points to can be
   addressed with the index list: -NBND..-1 , N+1..N+NBND for Boundary
   zones and 0..N for active zones.

   BE CAREFUL with macros, always add parethesis around non-singular
   defines in case used later with multiplations, etc.

   rm nohup.out ; nohup time -v sh -c './twod > 0_o.out' &


   Use program "nm" to list objects and symbols, or use objdump -axfhp

   A // PRECISION in the code means there could be a precision
   problem(addition machine precision error)

 */

/* 
   Generally, switches that are performance related are here, while
   others are in init with variable complements for flexibility during
   runtime.

   Not everything is in here so can compile quicker on simple changes.


 */
