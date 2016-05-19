#include "global.h"

/* 1: a: + dx in 1/2/3-direction from cell face: x(i)=x(i-1)+dx(i-1) */
/* 2: b: + dx in 1/2/3-direction from cell center: x(i)=x(i-1)+dx(i) */
FTYPE dxa[NUMGRID][3][NBIGM];	// [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]
FTYPE(*dx)[3][NBIGM];

/* Positions and sizes of cell */
/* cell face position 1: a */
/* cell center position 2: b */
FTYPE xa[NUMGRID][3][NBIGM];	// [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]
FTYPE(*x)[3][NBIGM];

FTYPE idxa[NUMGRID][3][INTNBIG];
FTYPE(*idx)[3][INTNBIG];

/* Positions and sizes of cell */
/* cell face position 1: a */
/* cell center position 2: b */
FTYPE ixa[NUMGRID][3][INTNBIG];
FTYPE(*ix)[3][INTNBIG];


/* at cell center 1: d 2: e 3: phi */
FTYPE sa[NUMSCA][N3M][N2M][N1M];	// [rho=1,en=2,pot=3][k=-N3BND..N3+B3BND][j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
FTYPE(*s)[N3M][N2M][N1M];	/* 3 scalars */

/* at cell face, where position of cell-face var is less than cell
   center var at same index */
/* 1: v 2: B */
FTYPE va[NUMVEC][3][N3M][N2M][N1M];	// [velocity=1,B-field=2][1=x1dir,2=x2dir,3=x3dir][k=-N3BND..N3+B3BND][j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
FTYPE(*v)[3][N3M][N2M][N1M];

// //// analytic storage
/* at cell center 1: d 2: e 3: phi */
FTYPE sanala[NUMSCA][N3M][N2M][N1M];	// same as s
FTYPE(*sanal)[N3M][N2M][N1M];

/* at cell face, where position of cell-face var is less than cell
   center var at same index */
/* 1: v 2: B */
FTYPE vanala[NUMVEC][3][N3M][N2M][N1M];	// same as v
FTYPE(*vanal)[3][N3M][N2M][N1M];

// used for floor on grid for scalar values.
// 1: rho->rho
// 2: u->u
// 3: phi->rho*phi (gravitational potential energy density)
FTYPE floorvarsa[NUMSCA][N3M][N2M][N1M];
FTYPE(*floorvars)[N3M][N2M][N1M];

// used for floor on grid for vectors(note: need 0 element here)
// all due to rho changes
// 1: momentum density in radial direction
// 2: theta-angular momentum density
// 3: phi-angular momentum density
FTYPE floorvarva[NUMVEC][3][N3M][N2M][N1M];
FTYPE(*floorvarv)[3][N3M][N2M][N1M];

// 0: KE
FTYPE floorvar0a[NUMVEC][N3M][N2M][N1M];
FTYPE(*floorvar0)[N3M][N2M][N1M];

// counts at different 5 different places in code(see step.c and
// sweep.c)
// these 2 vars are only sources for floor changes in above arrays
int floorcnt[NUMFLOOROUT + 1][NUMFLOORVAR + 1];	// [0,1,2,3,4,5,6][1=rho,2=en]
// [0][] : step_visc
// [5][] : step_real_visc
// [1][] : step_ie
// [2][] : step_res
// [3][] : sweepx
// [4][] : sweepy


// stores lowest value of rho and en before floor fixed, over entire
// run
FTYPE floorlowest[NUMFLOORVAR + 1];	// [1=rho,2=en]

// what routine has the above lowest value generated
// same routine numbers as above
int wherelowest[NUMFLOORVAR + 1];	// [1=rho,2=en]



// comp savers (first[] is location of variable used on)
// direction specifies direction of surface normal vector
FTYPE DSa[4][3][N3M][N2M][N1M];	// [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][dir=1,2,3][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]
FTYPE(*DS)[3][N3M][N2M][N1M];

FTYPE OARCLa[4][3][N3M][N2M][N1M];	// [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][dir=1,2,3][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]
FTYPE(*OARCL)[3][N3M][N2M][N1M];

FTYPE OVOLa[4][N3M][N2M][N1M];	// [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]
FTYPE(*OVOL)[N3M][N2M][N1M];



  /* metric scales and their dervs--only applies to limited coord
     systems */
  /* 1: a-mesh */
  /* 2: b-mesh */
  /* 1: g1 2: g2 3: g31 4: g32 */
FTYPE ga[NUMGRID][NUMMETRIC][NBIGM];	// [1,2][1,2,3,4][...]
FTYPE(*g)[NUMMETRIC][NBIGM];
FTYPE dga[NUMGRID][NUMMETRIC][NBIGM];	// [1,2][1,2,3,4][...]
FTYPE(*dg)[NUMMETRIC][NBIGM];

FTYPE cotana[NUMGRID][NBIGM];	// [1,2][...]
FTYPE(*cotan)[NBIGM];


  /* Volume factors: dvl[a/b][dir][location] */
FTYPE dvla[NUMGRID][3][NBIGM];	// [1,2][1,2,3][...]
FTYPE(*dvl)[3][NBIGM];


short bcsa[NUMSCA][3][N3M][N2M][N1M];	// [which
					// scalar][1=type:assign=(1,2,3,4,5),2=dim:assign=(1,2,3),3=direction:assign=(-1,+1)][...][...][...]
short (*bcs)[3][N3M][N2M][N1M];

// bcv assumes all components bounded with same *type* of boundary
// condition
short bcva[NUMVEC][3][N3M][N2M][N1M];	// [which
					// vector][1=type:assign=(1,2,3,4,5),2=dim:assign=(1,2,3),3=direction:assign=(-1,+1)][...][...][...]
short (*bcv)[3][N3M][N2M][N1M];

  /* type=1-4 are local, type=5 is nonlocal copy */
  /* see bound.c for types etc. */
// -1: corner
// 0: comp
// 1: reflect
// 2: AOS
// 3: Dirichlet
// 4: outflow
// 5: periodic

/* Global parameters */

// int N[3+1] ; /* Number of cells in each direction */
SFTYPE L[2 + 1][3 + 1];		/* L[1] holds position of 1,1,1 cell
				   outer edges, L[2] holds size of grid 
				   in units */
SFTYPE iL[2 + 1][3 + 1];

// outerdefs for images
SFTYPE outerdefs[ITYPES][CTYPES][NUMSCA + 1];
SFTYPE outerdefv[ITYPES][CTYPES][NUMVEC + 1][3 + 1];	// [][0] here
							// is magnitude

// outerdef for dumps for interpolation
SFTYPE douterdefs[NUMSCA + 1];
SFTYPE douterdefv[NUMVEC + 1][3 + 1];	// [][0] here is magnitude
SFTYPE douterdef0[NUMVEC + 1];	// KE for floor(not set yet)


SFTYPE IOBound[NUMBTRANS];	// used to have inflow/outflow
				// transitions on outer theta edge

short accountstorea[N3M][N2M][N1M];
short (*accountstore)[N2M][N1M];

// average stuff needed elsewhere too
SFTYPE tavgstart, tavgfinal;	// real start and final
int avgcount;			// number of averages in time
int num1d_31, num1d_32;

SFTYPE floors[NUMLOSSVAR + 1];	// to keep track of floor input of
				// stuff
SFTYPE inflows[NUMLOSSVAR + 1];	// to keep track of injection of stuff
SFTYPE radiations[NUMLOSSVAR + 1];	// to keep track of radiation
					// of stuff
// 1: mass
// 2: enthalpy
// 3: grav pot energy

// 4: ke

// 5: s1 (m*v_r in COORD==3)
// 6: s2 (r*m*v_theta in COORD==3)
// 7: s3 (r*sin(theta)*m*v_phi in COORD==3)
// 8: B1
// 9: B2
// 10: B3

// 11: etot(visc)

// next block of changable global parameters(global.h just sets initial 
// 
// value)
int kever;
int mdotin, press, visc_art, ie, mag, pmag, res, visc_real, trans;
int vischeat;
int vreal;
int transx1, transx2;
int transiex1, transmagx1, transv1x1, transv2x1, transv3x1, transrhox1;
int transiex2, transmagx2, transv1x2, transv2x2, transv3x2, transrhox2;
int advint;

SFTYPE coolfact, resist;	// not used

int ireenter;

SFTYPE rho0, R0, Omega0;	// for tori problem and for SPB rho
				// visc
char WRITETYPE[10];
int reallaststep;
int tagik, tagfk, tagij, tagfj, tagii, tagfi;
int t2gik, t2gfk, t2gij, t2gfj, t2gii, t2gfi;
int t3gik, t3gfk, t3gij, t3gfj, t3gii, t3gfi;
int t4gik, t4gfk, t4gij, t4gfj, t4gii, t4gfi;
char myidtxt[MAXFILENAME];
int totalzones;
int itotalzones;
int sizes[3 + 1][1];
int isizes[3 + 1][1];
int totalsize[3 + 1];
int itotalsize[3 + 1];
int startj;
int endj;			// startj and endj are where this CPU
				// located on full grid
long nstep;
int runtype, directinput, appendold, deleteolddat, deleteoldpar,
    deleteoldimg;
int simplebc, bcix1, bcox1, bcix2, bcox2, bcix3, bcox3;
int nonunigridx1, nonunigridx2, nonunigridx3;
int analoutput;
SFTYPE x1in, x1out, x2in, x2out, x3in, x3out;
SFTYPE timereenter;
int pdump_start, dump_start, npdump_start, adump_start, floor_start,
    image_start;
FILE *fail_file;
FILE *logfull_file;
FILE *log_file;
FILE *logstep_file;
FILE *logperf_file;
FILE *logdt_file;
FILE *logfl_file;
FILE *logsp_file;
// Equation of state coefficients
SFTYPE tscycleto;		// time to subcycle to
SFTYPE tscyclefrom;		// time to subcycle from
SFTYPE dtlastscycle;		// dt used as basis for subcycle
int subcyclen;			// number of subycycles for viscosity
int nthsubcycle;		// number of subcycles so far
SFTYPE rg, rgp;
SFTYPE cour, invcour;		// invcour used in timestep1.h
SFTYPE cour2, invcour2;		// invcour2 used in timestep1.h
SFTYPE cs;
SFTYPE alpha;
SFTYPE gam;
int wgam;
int wgam1;
int wgam53;
int wpw;
int periodicx1, periodicx2, periodicx3;
int skipix1, reflectix1, reflectox1;
int skipix2, reflectix2, reflectox2;
int skipix3, reflectix3, reflectox3;
int intix1, intox1, intix2, intox2, intix3, intox3;
int skipintix1, skipintix2, skipintix3;
SFTYPE dt;
SFTYPE tstart, TSTART;
SFTYPE t, tf, timagescale;
SFTYPE tavgi, tavgf;
int numavg;
SFTYPE nu_vnr;
SFTYPE nu_l;
SFTYPE nu_ten;
SFTYPE alpha_real;
SFTYPE alpha_real0;
FTYPE nu_facta[N3M][N2M][N1M];
FTYPE(*nu_fact)[N2M][N1M];
FTYPE nu_reala[N3M][N2M][N1M];
FTYPE(*nu_real)[N2M][N1M];
SFTYPE n_real;
SFTYPE GRAVC;
SFTYPE MASSBH;
SFTYPE GM;
SFTYPE DTd;
SFTYPE DTpd;
SFTYPE DTavg;
SFTYPE DTi;
SFTYPE DTl;
SFTYPE DTener;
SFTYPE DTloss;
SFTYPE DTfloor;
SFTYPE DTtimestep;
SFTYPE DTsp;
SFTYPE nu_sh;
// int numbc[3+1];
SFTYPE startx[3 + 1];		// starting location for each cpu for
				// each direction
SFTYPE startix[3 + 1];		// starting location for each cpu for
				// each direction
/* grid velocity */
SFTYPE vg[3 + 1];
SFTYPE RHO0, V0;		// for advection
SFTYPE DENSITYFLOOR, IEFLOOR;
SFTYPE IEFRACT, VZFRACT;
SFTYPE massdot;


int numprocs, myid, procnamelen;	// not used

// only over active grid
SFTYPE losssa[NUMSCA][3][2][NBIG];	// [which=1...][dir=1,2,3][i/o=0,1][list=0...]
SFTYPE(*losss)[3][2][NBIG];
SFTYPE lossva[NUMVEC][3 + 1][3][2][NBIG];	// [which][v-dir][dir][i/o][list]
SFTYPE(*lossv)[3 + 1][3][2][NBIG];
SFTYPE lossvisca[1][3][2][NBIG];	// [which][dir][i/o][list]
SFTYPE(*lossvisc)[3][2][NBIG];
// v-dir: vector direction // 0=energy, 1,2,3
// dir: 1=x1 2=x2 3=x3
// i/o: 0=inner boundary 1=outer boundary
// list: var along that side: on active grid: 0..N-1(i.e. don't shift
// pointer!)



int globalinterpmod;


FTYPE rhoinjecta[N3M][N2M][N1M];
FTYPE(*rhoinject)[N2M][N1M];
FTYPE eninjecta[N3M][N2M][N1M];
FTYPE(*eninject)[N2M][N1M];

// used to hold boundary data(packed)
FTYPE worksbca[2][2][N2BND * N1M];	// scalar:
					// [out/in][-j/+j][datawidth][data]
FTYPE(*worksbc)[2][N2BND * N1M];
FTYPE workvbca[2][2][3 * N2BND * N1M];	// vector:
					// [out/in][[-j/+][datawidth][data]
FTYPE(*workvbc)[2][3 * N2BND * N1M];


// workXa is just like s[i] (i.e. one particular scalar)
FTYPE work1a[N3M][N2M][N1M];
FTYPE work2a[N3M][N2M][N1M];
FTYPE work3a[N3M][N2M][N1M];
FTYPE work4a[N3M][N2M][N1M];
FTYPE work5a[N3M][N2M][N1M];
FTYPE work6a[N3M][N2M][N1M];
FTYPE work7a[N3M][N2M][N1M];
FTYPE work8a[N3M][N2M][N1M];
FTYPE work9a[N3M][N2M][N1M];
FTYPE work10a[N3M][N2M][N1M];
FTYPE(*work1)[N2M][N1M];
FTYPE(*work2)[N2M][N1M];
FTYPE(*work3)[N2M][N1M];
FTYPE(*work4)[N2M][N1M];
FTYPE(*work5)[N2M][N1M];
FTYPE(*work6)[N2M][N1M];
FTYPE(*work7)[N2M][N1M];
FTYPE(*work8)[N2M][N1M];
FTYPE(*work9)[N2M][N1M];
FTYPE(*work10)[N2M][N1M];

// sensitives
SFTYPE works1a[N3M][N2M][N1M];
SFTYPE(*works1)[N2M][N1M];

// workvXa is just like v[i] (i.e. like one particular scalar)
FTYPE workv1a[3][N3M][N2M][N1M];
FTYPE workv2a[3][N3M][N2M][N1M];
FTYPE workv3a[3][N3M][N2M][N1M];
FTYPE workv4a[3][N3M][N2M][N1M];
FTYPE workv5a[3][N3M][N2M][N1M];
FTYPE(*workv1)[N3M][N2M][N1M];
FTYPE(*workv2)[N3M][N2M][N1M];
FTYPE(*workv3)[N3M][N2M][N1M];
FTYPE(*workv4)[N3M][N2M][N1M];
FTYPE(*workv5)[N3M][N2M][N1M];

// [1,2,3][1,2,3][...][...][...]
// below 3 tensors only used if real viscosity is present
FTYPE workt1a[3][3][N3M][N2M][N1M];
FTYPE(*workt1)[3][N3M][N2M][N1M];

FTYPE workt2a[3][3][N3M][N2M][N1M];
FTYPE(*workt2)[3][N3M][N2M][N1M];

FTYPE workt3a[3][3][N3M][N2M][N1M];
FTYPE(*workt3)[3][N3M][N2M][N1M];

// workiqa is like s but in 2D and only on active grid
// [0..n2-1][0..n1-1]
FTYPE workiqa[NUMSCA][INTN2][INTN1];
FTYPE(*workiq)[INTN2][INTN1];

// workviqa is like v but in 2D and only on active grid
FTYPE workviqa[NUMVEC][3][INTN2][INTN1];
FTYPE(*workviq)[3][INTN2][INTN1];

FTYPE work0iqa[NUMVEC][INTN2][INTN1];
FTYPE(*work0iq)[INTN2][INTN1];

FTYPE mms[ITYPES][CTYPES][NUMSCA + 1][2];
FTYPE mmv[ITYPES][CTYPES][NUMVEC + 1][3 + 1][2];	// 0 is
							// magnitude
FTYPE mmst[ITYPES][CTYPES][NUMSCA + 1][2];
FTYPE mmvt[ITYPES][CTYPES][NUMVEC + 1][3 + 1][2];
