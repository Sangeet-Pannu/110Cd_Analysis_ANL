
#ifndef _GEBSort_h
#define _GEBSort_h


#include "TFile.h"
#include "gdecomp.h"
#include "pull_tubes2.h"

#define TRACK2 1

#define LINPOL_BEAM 1
#define LINPOL_SOURCE 2

#define TAPE 0
#define NET  1
#define DISK 2
#define GEB  3

#define MEDIUMLEN 2048
#define LONGLEN 20000
#define SHORTLEN 1024
#define RBUFSIZE 500000
#define MAXPAYLOADSIZE 500000
#define STRLEN 256

#define MAXDETPOS 180
#define MAXCRYSTALNO 3
#define MAXDETNO 200
#define MAXMODNO 60
#define NUMAGATAPOS 180
#define NGATE_SPE 20000

/* veto cube */
#define MAXGTMODNO 30

#define RATELEN 60*24
#define DTBTEVLEN 1000
#define MAXNOSEG 9
#define NUMSEGS 36
#define MAXSEGNO MAXDETNO*NUMSEGS
#define RMIN 10
#define RMAX 35

#define MAXGEBS 1000
#define MAXLONG (long long)1<<62

#define MAX_GAMMA_RAYS 1000
#define GEB_BITS GEB_HEADER_BYTES*8

/* max values for # of bits */

#define M14BITS 0x3fff 
#define M13BITS 0x1fff
#define M12BITS 0x0fff
#define M11BITS 0x07ff
#define M10BITS 0x03ff

/* basic spectra lengths */

#define L14BITS  M14BITS+1
#define L13BITS  M13BITS+1
#define L12BITS  M12BITS+1
#define L11BITS  M11BITS+1
#define L10BITS  M10BITS+1

typedef struct EXCHANGE
{

  /* bin_dgs */

  /* bin_dfma */

  /* bin_mode3 */

  /* bin_mode2 */

  /* bin_mode1 */

  int ngates;

} EXCHANGE;

typedef struct GEB_event
{
  int mult;
  CRYS_INTPTS *ptinp[MAXGEBS];
  GEBDATA *ptgd[MAXGEBS];
} GEB_EVENT;

typedef struct PARS
{
  char ROOTFile[STRLEN];
  int nEvents;
  char ROOTFileOption[STRLEN];
  char GTSortInputFile[STRLEN];
  int UseShareMemFile;
  unsigned int StartMapAddress;
  char ShareMemFile[STRLEN];
  int InputSrc;
  int HaveRootFileName;
  int WeWereSignalled;
  int UseRootFile, SizeShareMemFile;
  int UpdateRootFile;
  char spname[STRLEN];
  int firstEvent;
  int GSudpPort;
  int NumToPrint;
  int DumpEvery;
  TFile *f1;
  TList *wlist;
  long long int curTS;
  long long int maxTS, tmpmaxTS;
  long long int dTS;
  long long int nbytes;
  int CurEvNo;
  char pHost[16];
  int Run_number;
  int grouping;
  int type;
  int enabled[MAXDETNO+1];
  float CCcal_gain[MAXDETNO+1];
  float CCcal_offset[MAXDETNO+1];
  float SEGcal_gain[MAXDETPOS + 1][MAXCRYSTALNO + 1];
  float SEGcal_offset[MAXDETPOS + 1][MAXCRYSTALNO + 1];
  float timeout;
  float crmat[MAXDETPOS + 1][MAXCRYSTALNO + 1][4][4]; 
  float detpolang[MAXDETPOS+1];
  float beta;
  float vc_xa;
  int GGMAX;
  int modwrite;
  int tsnumwrites;
  float fomlo[MAXNOSEG];
  float fomhi[MAXNOSEG];
  int ndetlimlo;
  int ndetlimhi;
  float beamdir[3];
  int nocrystaltoworldrot;
  int multlo;
  int multhi;
  int requiretracked;
  float  modCCang[MAXMODNO+1][MAXCRYSTALNO+1];
  float  modCCdopfac[MAXDETNO];
  int AGATA_data;
  double TrX[NUMAGATAPOS], TrY[NUMAGATAPOS], TrZ[NUMAGATAPOS];
  double rotxx[NUMAGATAPOS], rotxy[NUMAGATAPOS], rotxz[NUMAGATAPOS];
  double rotyx[NUMAGATAPOS], rotyy[NUMAGATAPOS], rotyz[NUMAGATAPOS];
  double rotzx[NUMAGATAPOS], rotzy[NUMAGATAPOS], rotzz[NUMAGATAPOS];
  int numgggates;
  int gg_gate_pos[100];
  int gg_gate_width[100];
  float target_x, target_y, target_z;
  char AGATA_data_fn[512];
  int minnumHitArray;
  int maxnumHitArray;
  int minNumGammas;
  int maxNumGammas;
  int nisomers;
  int minNumCC;
  int maxNumCC;
  float minCCe;
  int crystalID3D;
  float crystalID3D_elo;
  float crystalID3D_ehi;
  float crystalID3D_fomlo;
  float crystalID3D_fomhi;
  int do_bin_mode3;
  int do_bin_mode2;
  int do_bin_ndc;
  int do_bin_mode1;
  int do_bin_sicd;
  int do_bin_dfma;
  int do_bin_dgs;
  int do_bin_dub;
  int do_bin_XA;
  int do_bin_gtcal;
  int do_bin_template;  
  int do_bin_ft;  
  int do_bin_s800;  
  int do_bin_final;  
  int do_bin_angcor_GT;  
  int do_bin_DCO_GT;  
  int do_bin_angcor_DGS;  
  int do_bin_angdis;
  int do_bin_linpol;
  int echo_data;
  char echo_data_fn[512];
  off_t echo_data_pipe;
  unsigned int waitfordataseconds;
  float Hresolution;

  float linpol_rrlo;
  float linpol_rrhi;
  float linpol_polmin;
  float linpol_polmax;
  float linpol_scatmin;
  float linpol_scatmax;
  float linpol_source_ee;
  float linpol_source_de;
  int   linpol_usesource;
  int   linpol_nlo;
  int   linpol_nhi;
  int   linpol_ngmin;
  int   linpol_random;
  int   angcor_useplaneang;

  int   do_bin_g4sim;

  int   dgs_algo;
  double dgs_SZ_t1;
  double dgs_SZ_t2;
  float dgs_MM;
  float dgs_KK;
  char  dgs_PZfn[256];
  char  dgs_ecalfn[256];  
  char  dgs_factorfn[256];

  int   xa_algo;
  double xa_SZ_t1;
  double xa_SZ_t2;
  float xa_MM;
  float xa_KK;
  char  xa_PZfn[256];
  char  xa_ecalfn[256];  
  char  xa_factorfn[256];

  float dub_MM;
  float dub_PP;

  char  xyzoffsetfn[256];
  float xyzoffset[MAXDETNO+1][3];
  int   havexyzoffset;

  float DCO_ep1, DCO_ep2, DCO_pde, DCO_bde, DCO_cang; 
  float DCO_eb1, DCO_eb2;

  float crystalxx[MAXDETPOS + 1];
  float crystalyy[MAXDETPOS + 1];
  float crystalzz[MAXDETPOS + 1];

  float crystalrmin;
  float crystalrmax;

  float smapmine;
  float smapmaxe;

  float gb_dt_lo, gb_dt_hi;
  float decay_lo1, decay_hi1;
  float decay_lo2, decay_hi2;
  float decay_lo3, decay_hi3;
  float ggdt;
  int   have_decay_data;

  int HKintubes;

 int S800_PID_NX;
 float S800_PID_XLO;
 float S800_PID_XHI;
 int S800_PID_NY;
 float S800_PID_YLO;
 float S800_PID_YHI;
 float ytascale;
 float greta_ata;
 float greta_bta;
 int s800_status;
 int s800_inwindow;
 int havePIDwin;
 char PIDwinfile[STRLEN];
 float s800_csoff;
 float shellE;
 int justFilterPID;

} PARS;

/* structure for the tracked gamma rays */
/* written to the output with geb ID GEB_TYPE_TRACK==3 */

#if (TRACK2==1)

/* new format where we pass on the first hit crystal ID */

typedef struct TRACKED_GAMMA_RAY {
  float esum;                   /* gamma ray energy */
  int ndet;                     /* number of interactions */
  float fom;                    /* figure of merit */
  short int tracked;            /* 1==if tracked */
  long long int timestamp;      /* timestap of first interaction point */
  float x0, y0, z0, e0;         /* first interaction point */
  float x1, y1, z1, e1;         /* second interaction point */
  short int fhcrID;             /* first hit crystal ID */
  } TRACKED_GAMMA_RAY;

#endif

#if (TRACK2==0)

typedef struct TRACKED_GAMMA_RAY {
  float esum;                   /* gamma ray energy */
  int ndet;                     /* number of interactions */
  float fom;                    /* figure of merit */
  int tracked;            /* 1==if tracked */
  long long int timestamp;      /* timestap of first interaction point */
  float x0, y0, z0, e0;         /* first interaction point */
  float x1, y1, z1, e1;         /* second interaction point */
  } TRACKED_GAMMA_RAY;

#endif


typedef struct TRACKED_GAMMA_HIT {
      int ngam;
      int pad;
      TRACKED_GAMMA_RAY gr[MAX_GAMMA_RAYS];
    } TRACKED_GAMMA_HIT;

/*
 * S800 pre-processed (currently poor Dirk's version)
 *
 * PID plot made from tof_obje1/xfpe1/rfe1 versus ic _de
 *
 * ata/bta/dta/yta used for Doppler reconstruction
 *
 * gap = 1073mm, zfp = 0.
 * fp_afp = atan((crdc2_x - crdc1_x) / gap)
 * fp_bfp = atan((crdc2_y - crdc1_y) / gap)
 * fp_xfp = crdc1_x / 1000 + zfp * tan(fp_afp);
 * fp_yfp = crdc1_y / 1000 + zfp * tan(fp_bfp);
 *
 *
 */
#define S800PHYSDATA_NOTVALIDDATA -9999.999
#define S800PHYSDATA_TYPETAG  0xABCD1234
struct S800_physicsdata {
  int32_t type;    /* defined abcd1234 for indicating this version */
  float crdc1_x;   /* Crdc x/y positions in mm */
  float crdc1_y;
  float crdc2_x;
  float crdc2_y;
  float ic_sum;    /* ion chamber energy loss         */
  float tof_xfp;   /* TOF scintillator after A1900    */
  float tof_obj;   /* TOF scintillator in object box  */
  float rf;        /* Cyclotron RF for TOF            */ 
  int32_t trigger; /* Trigger register bit pattern    */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /* from here corrected values extracted from data above */ 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  float ic_de;
  /* TOF values with TOF correction applied (from afp/crdc x) */
  float tof_xfpe1;
  float tof_obje1;
  float tof_rfe1;
  /* Trajectory information at target position calculated from 
     a map and afp/bfp/xfp/yfp. New map and you need to re-calc */
  float ata; /* dispersive angle        */
  float bta; /* non-dispersive angle    */
  float dta; /* dT/T T:kinetic energy   */
  float yta; /* non-dispersive position */
};
typedef struct S800_physicsdata S800_PHYSICSDATA;

/* for Dirk's stuff */

typedef struct COOR {
  float x;
  float y;
  float z;
  float theta;
  float phi;
} COOR;



/* macros */

#define WRITEALLHISTS  \
  gDirectory->cd("root:/"); \
  wlist = gDirectory->GetList(); \
  if(!Pars.UpdateRootFile) f1= new TFile(Pars.RootFile, Pars.ROOTFileOption); \
  if (ComPressLevel>NOTDEF) f1->SetCompressionLevel(ComPressLevel); \
  printf("writing all spectra to [%s]\n", Pars.RootFile); \
  printf("be patient... "); \
  fflush(stdout); \
  t1 = time(NULL); \
  wlist->Write(0,TObject::kOverwrite); \
  t2 = time(NULL); \
  printf("DONE! on "); \
  time_stamp(stderr); \
  printf("file size: %i, ",f1->GetSize()); \
  printf("compression level: %i, ",f1->GetCompressionLevel()); \
  printf("and factor: %f\n",f1->GetCompressionFactor()); \
  printf("uncompressed root file size: %f\n",f1->GetSize()*f1->GetCompressionFactor()); \
  printf("writeout time: %i seconds\n", t2 - t1); \
  printf("at %7.2f Mbytes/sec\n", (float) f1->GetSize() / (t2 - t1) / 1000000.0); \
  printf("on "); \
  time_stamp(stderr); \
  fflush(stdout);

#define UPDSSHMEM \
  t1 = time(NULL); \
  mfile->Update(); \
  t2 = time(NULL); \
  printf("done! "); \
  printf("shared memory size: %i\n", mfile->GetSize()); \
  printf("update time: %i seconds ", t2 - t1); \
  printf("at %7.2f Mbytes/sec\n", (float) mfile->GetSize() / (t2 - t1) / 1000000.0); \
  printf("to mapfile [%s] on ",Pars.ShareMemFile); \
  time_stamp(stderr); \
  fflush(stdout);


  /* prototypes */

//  TH2F *mkTH2F (char *, char *, int , double , double , int , double , double );
//  TH1D *mkTH1D (char *, char *, int , double , double );

#endif	/* _GEBSort_h */
