//#include "GEBSort.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObject.h"
#include "TCutG.h"
#include "TMath.h"
#include "TString.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLatex.h"
#include "TTree.h"
//#include "TVector3.h"
using namespace std;

#include "GEBSort.h"
#include "GTMerge.h"
#include "veto_pos.h"
//#include "GTmode3.h"
#define Mvalue 128


int good1=0, bad1=0, good2=0, bad2=0, good3=0, bad3=0, type1=0, type3=0, paallonly=0, parinonly=0, paseconly=0, paallgam=0, paringam=0, pasecgam=0, gronly=0;
bool pippo=false;




const int si_ethreshold = 1; // in channel
const double si_calethreshold[2] = {1.,300.}; //Range of energy for good events
//const double GT_ethreshold = 150; //keV
// detector structure paramters
/*const double r_min = 11.532;
const double r_max = 35.00;
const double r_dead = 0.2;
const double n_rings = 24;
const double z = 500; // in mmm distance to target 5 is place holder value
*/
double ringID[24];
double Theta_Pb[24];
double Theta_C[24];
double Theta_com[24];
double Beta_Pb[24];
double Beta_C[24];
char namekinefile_C[100] = "kine_4MeVex.txt";
char namekinefile_Pb[100] = "kine_2MeVex.txt";
const unsigned long binGamE=4000, binGamE2D=4000; const double minGamE=50, maxGamE=4050.;
const int binSiE=380; const double minSiE=0, maxSiE=380; // Raw Energy Ch
const int binSiRawE=2000; const double minSiRawE=0, maxSiRawE=2000; // Raw Energy Ch
const int binDiffTS=1000; const double minDiffTS=-499.5, maxDiffTS=499.5;
const long timegate[2]={-20,240}, offtimegate[2]={100,120};
///
// Energy calibaration for SiCDs
double Ecal_ring_Gain[24] = {7.2843, 7.1341, 7.2843, 7.1341, 6.076, 7.2764, 7.257, 7.1797, 7.0558, 7.181, 7.1058, 6.93910, 7.295,
                              7.1633, 6.9159, 6.9285, 7.2465, 7.3192, 7.145, 7.214, 7.0, 7.251, 6.65, 5.709};

                       //ID 23 and 24 have almost no hits though.

double Ecal_sector[32] = {0.016595, 0.0164496, 0.0163965, 0.0164803, 0.01656, 0.0163513, 0.0164831, 0.0165203, 0, 0.0164597,
                               0.0168349, 0.0172811, 0.0167323, 0.0170783, 0.0169721, 0.0168064, 0.0168493, 0.0166977, 0.0167142, 0.0171135,
                               0.0169216, 0.0167393, 0.0166912, 0.0168298, 0, 0.0169825, 0.0172367, 0.0167626, 0.016903, 0.0166209, 0.0166205, 0.0166496};
                               // For sec22, a fit result for low E peak is manually obtained

double Ecal_sector_Gain[32] = {7.22177,7.22177,7.22177,7.22177,7.22177,7.22177,7.22177,7.22177,1.0,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.341,7.1369,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734,7.1734};
double Ecal_sector_Offset[32] = {987.225,987.225,987.225,987.225,987.225,987.225,987.225,987.225,0,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81,881.069,944.31,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81,964.81};

double Ecal_ring[24] = {0.0164025, 0.0164205, 0.0162874, 0.016611, 0.0192222, 0.0163881, 0.0171403, 0.0163095, 0.0165847, 0.0163237,
                       0.0165935, 0.0169492, 0.0161604, 0.016165, 0.016738, 0.0166998, 0.0164376, 0.0162197, 0.0163854, 0.0161243,
                       0.0164972, 0.015859, 0.0178253, 0.0197221};
                       //ID 23 and 24 have almost no hits though.

double Ecal_ring_Offset[24] = {851.7145938,857.71,851.7145938,857.71,902.76,917.6011,952.31,919.46,950.54,904.71,987.513,924.866,1032.43,
                        983.16,1018.03,1016.90,1030.19,1053.499,1057.678,1012.86,1162.07,1090.155,1077.42,1093.25};

static ofstream output1 ("XYZHitPatterns.txt");
static ofstream output2 ("SolidAngleCorr.txt");

#define MaxGTNum 28
#define MaxChicoNum 10


typedef struct SiCD_struct
{
  Double_t  LEDts_ring;
  Double_t  LEDts_sector;
  int id_ring;
  int id_sector;
  int detid_ring;
  int detid_sector;
  float rawE_ring;
  float rawE_sector;
  float calE_ring;
  float calE_sector;
  unsigned int status;
  double theta_C;
  double theta_Pb;
  double phi_C;
  double phi_Pb;
  double beta_C;
  double beta_Pb;
  //int			  RF;
  //bool			  SINGLE;
};

typedef struct part_evt {
  long ring_time;
  long sector_time;
  float ring_energy;
  float sector_energy;
  float phi_C;
  float theta_C;
  float phi_Pb;
  float theta_Pb;
  float beta_C;
  float beta_Pb;
  int ring_id;
  int sector_id;
};

typedef struct gamma_evt {
  long time;
  float gamma_energy;
  float gamma_theta;
  float gamma_phi;
};

extern PARS Pars;
extern EXCHANGE exchange;

extern int nCCenergies;
extern float CCenergies[MAX_GAMMA_RAYS];
extern unsigned long long int CCtimestamps[MAX_GAMMA_RAYS];
extern int CCids[MAX_GAMMA_RAYS];


///prototype
  int twoscomp_to_int_24 (unsigned int );


  TH2F *GidMul2;
// ProcessEventMode1
 TH2F *GT_TS;
  TH2F *SiCd_TS_R;
   TH2F *SiCd_TS_S;


TH1D * sectVsPhiAngle;

TH1D * EnergySpecVANG[24];

TH1D *gMultiplicity;
TH1D *gMultiplicityR;

TH2F *ParticleAngle;

TH2F *sicd_e_detid;

TH2F *SiCdAngleMap;

TH1D *LabEnoTG_C;

TH1D *Calculated_S;
TH2F *Qplot;

TH2F *PEVE,*GtVE;

TH2F *ggE, *ggEF;
TH2F *ggEBK, *ggEBKF;
TH2F *GTENCAL;
//TH2F *BetaCal;
TH1D *GTdiff;

TH2F *SiCD_Correlation;

TH1D *h1_ene88, *LabEnoTG, *LabEinTG, *LabEoffTG, *h1_ene882,*LabEoffTG_C,*LabEinTG_C,*SiCD_ESec_IdRing;
TH1D *LabEnoTGF, *LabEinTGF, *LabEoffTGF;
TH2F *h2_ene88;
TH2F *GTRVGTS;
TH2F *GT_diffTS_sector;
TH2F *GT_diffTS_ring;
//TLatex *Latex_BeamDirection;
TH2F *SiCD_XY;
TH1D *diff_SiCDTS, *Si_sect_diffTS, *Si_ring_diffTS;
TH2F *SiCDERatio_diffTS_sector, *SiCDERatio_diffTS_ring, *SiCDE_diffTS, *DoppEPb_diffTS, *labE_diffTS, *SiCD_ESec_ChRing, *SiCD_ERing_IdSec;
TH1D *doppE_C, *doppE_Pb, *doppEnoTG_C, *doppEnoTG_Pb, *doppEoffTG_C, *doppEoffTG_Pb, *labE, *labEoffTG;
TH2F *dcfac_corr, *doppE_angleC, *doppE_anglePb, *corr_angle_CPb;
TH2F *doppE_thetaC, *doppE_thetaPb, *doppE_phiC, *doppE_phiPb;
TH2F *labE_angleC, *labE_anglePb;
TH2F *labE_thetaC, *labE_thetaPb, *labE_phiC, *labE_phiPb;


//**********************


//////////////////////////////////////////////////////////////////////

void load_kinematics(){
/////// Load kinematic calculations
  fstream kinefile_C;
  kinefile_C.open(namekinefile_C);
  cout<<"Carbon kinematics file: "<<namekinefile_C<<endl;
  fstream kinefile_Pb;
  kinefile_Pb.open(namekinefile_Pb);
  cout<<"Lead kinematics file: "<<namekinefile_Pb<<endl<<"Carbon "<<endl;
  int num = 0;
  //int ch =0;
  int RingID;
  double Theta_Beam,Theta_Target,Theta_COM,Beta_Beam,Beta_Target;
  while (!kinefile_C.eof()) {
    if(num>=24) break;
    kinefile_C>>RingID>>Theta_Beam>>Theta_Target>>Theta_COM>>Beta_Beam>>Beta_Target;
    ringID[num] = RingID + 1; // Detector IDs are defined to start from 1
    Theta_Pb[num] = Theta_Beam;
    Theta_C[num] = Theta_Target;
    Theta_com[num] = Theta_COM;
    Beta_Pb[num] = Beta_Beam;
    Beta_C[num] = Beta_Target;
    cout<<num<<" "<<RingID<<" "<<Theta_Beam<<" "<<Theta_Target<<" "<<Theta_COM<<" "<<Beta_Beam<<" "<<Beta_Target<<endl;
    num++;
  }
}


int map_SiCD(SiCD_struct *SiCD){
  int sector_channel[32] = {10,11,12,13,14,15,16,17,20,21,22,23,24,25,26,27,30,31,32,33,34,35,36,37,0,1,2,3,4,5,6,7}; // DAQ Channel
  
 // int sector_channel[32] = {10,11,12,13,14,15,16,17,20,21,22,23,24,25,26,27,0,1,2,3,4,5,6,7,30,31,32,33,34,35,36,37}; // DAQ Channel
  int sector_number[32] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32}; // Physical ID
  double sector_phi[32] = {270,258.75,247.5,236.25,225,213.75,202.5,191.25,180,168.75,157.5,146.25,135,123.75,112.5,101.25,90,78.75,67.5,56.25,45,33.75,22.5,11.25,0,-11.25,-22.5,-33.75,-45,-56.25,-67.5,-78.75};
  
  for(int i=0; i<=32; i++){
    if(SiCD->id_sector != sector_channel[i]) continue;

    SiCD->detid_sector = sector_number[i];
    //Ecal
    SiCD->calE_sector =  (SiCD->rawE_sector * Ecal_sector[SiCD->detid_sector-1])*Ecal_sector_Gain[SiCD->detid_sector-1]+Ecal_sector_Offset[SiCD->detid_sector-1];
    //
    SiCD->phi_Pb = sector_phi[i] * TMath::DegToRad();

    if(SiCD->phi_Pb < 0)SiCD->phi_Pb+= 2. * TMath::Pi();

   //  SiCD->phi_C = SiCD->phi_Pb;

    //================[Flips xy coordinates to negative axis.]==================================================|
    if(SiCD->phi_Pb > TMath::Pi()){
      SiCD->phi_C = SiCD->phi_Pb - TMath::Pi();
    }else{
      SiCD->phi_C = SiCD->phi_Pb + TMath::Pi();
    }  

    //================[Flips xy coordinates to negative axis.]==================================================|

    break;
  }
  

 // double ring_ch[24] = {25,24,23,22,21,20,0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17,27,26}; // Mo double ring_ch[24] = {17,16,15,14,13,12,11,10,7,6,5,4,3,2,1,0,20,21,22,23,24,25,26,27};
  double ring_ch[24] = {27,26,25,24,23,22,21,20,10,11,12,13,14,15,16,17,0,1,2,3,4,5,6,7};
  double chan[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
  
  for(int i=0; i<=24; i++){

    if(SiCD->id_ring != ring_ch[i]) continue;
    SiCD->detid_ring = chan[i];
    //Ecal
    SiCD->calE_ring =  ((double)SiCD->rawE_ring * Ecal_ring[SiCD->detid_ring-1])*Ecal_ring_Gain[SiCD->detid_ring-1]+Ecal_ring_Offset[SiCD->detid_ring-1];
   
    SiCD->theta_C = Theta_C[i] * TMath::DegToRad();
    SiCD->theta_Pb = Theta_Pb[i] * TMath::DegToRad();

    SiCD->beta_C = Beta_C[i];
    SiCD->beta_Pb = Beta_Pb[i];
    break;
  }    

  return 0;
}



//////////////////////////////////////////////////////////////////////
float cos_Gamma_Recoil(float ptheta, float pphi, float gtheta, float gphi){
  float gccos;
  gccos=sinf(ptheta)*sinf(gtheta)*cosf(pphi-gphi)+cosf(ptheta)*cosf(gtheta);
  return gccos;
}


//////////
int sup_sicd(){
  char str1[STRLEN], str2[STRLEN];
  char str[127];
  char fn[127];

  int i,j,k;
  string OneLine;


  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);
  TH3F *mkTH3F (char *, char *, int, double, double, int, double, double,int, double, double);

  h1_ene88 = mkTH1D("ene88Sector","Calibrated energy of SiCD; Energy (keV)",700,0,70);
  h1_ene882 = mkTH1D("ene88Ring","Calibrated energy of SiCD; Energy (keV)",700,0,70);

  h2_ene88 = mkTH2F("ene88vsPhiAng","Energy vs. phi angle of detection",100,130,180,50000,0,50000);

 // bank88_multiplicity = mkTH1D("bank88_multiplicity","bank88_multiplicity",100,0,100);
  SiCd_TS_R = mkTH2F("SiCd_TS_R", "SiCD RING Time Stamp (ns) vs. Events",10000,0,1000000,100,7700,9730);
  SiCd_TS_S = mkTH2F("SiCd_TS_S", "SiCD SECTOR Time Stamp (ns) vs. Events",10000,0,1000000,100,7700,9730);

  GTENCAL = mkTH2F("SGTENCAL", "GRETINA Energy (Doppler CORRECTED) vs. Crystal number",100,40,140,3000,0,3000);

  gMultiplicity = mkTH1D("GMultiplicity","Gamma Mulitplicity",50,0,50); 
  gMultiplicityR = mkTH1D("GMultiplicityR","Gamma Mulitplicity",50,0,50); 

  GidMul2 = mkTH2F("ggMult2", "Hit 1 vs Hit 2 (ids)",90,40,130,90,40,130);

  PEVE = mkTH2F("ParticleEnergyVGammaEnergy", "Particle Energy Vs. Gamma Energy",700,0,70,3000,0,3000);
  GtVE = mkTH2F("GammaTimeVGammaEnergy", "Gamma time Vs. Gamma Energy",5000, -2500, 2500,3000,0,3000);

  GTdiff = mkTH1D("GTdiff","Gretina Time difference ns ",1000,-500,500);
  sectVsPhiAngle = mkTH1D("sectVsPhiAngle","Sector Phi angle vs. Counts",64,-360,360);
  sectVsPhiAngle->SetYTitle("Counts");
  sectVsPhiAngle->SetXTitle("Phi Angle");

  SiCdAngleMap = mkTH2F("SiCd Angle map","phi vs Theta",64,-360,360,50, 130, 180);

  Qplot = mkTH2F("QPlot", "Q-value (expected particle energy - detected) Vs Theta_p",50,130,180,100, -50, 50);

  ggE = mkTH2F("ggE", Form("#gamma-#gamma, t_{#gamma-#gamma} < 50 && t_{#gamma-#gamma} > 0"), 3000, 0, 3000,3000,0,3000);
  ggEBK = mkTH2F("ggEBK", Form("#gamma-#gamma, t_{#gamma-#gamma} < 150 && t_{#gamma-#gamma} > 100"), 3000, 0, 3000,3000,0,3000);

  ggEF = mkTH2F("ggE_fine", Form("#gamma-#gamma (fine-grain), t_{#gamma-#gamma} < 50 && t_{#gamma-#gamma} > 0"), 6000, 0, 3000,6000,0,3000);
  ggEBKF = mkTH2F("ggEBK_fine", Form("#gamma-#gamma (fine-grain), t_{#gamma-#gamma} < 150 && t_{#gamma-#gamma} > 100"), 6000, 0, 3000,6000,0,3000);

  SiCD_Correlation = mkTH2F("SICD_Cor","Correlation between Ring and Sector Energy", 100, 0, 100,100,0,100);
  SiCD_Correlation->SetYTitle("Ring Energy");
  SiCD_Correlation->SetXTitle("Sector Energy");

  SiCD_ESec_IdRing= mkTH1D("SiCD_ESec_IdRing","Id Ring vs. Correlation between Ring and Sector Energy",10,0,10);
  SiCD_ESec_IdRing->SetYTitle("Ratio of Ring to Sector Energy");
  SiCD_ESec_IdRing->SetXTitle("Id Ring");

  SiCD_ERing_IdSec= mkTH2F("SiCD_ERing_IdSec","Id Sector vs. Correlation between Ring and Sector Energy", 32, 0, 32,10,0,10);
  SiCD_ERing_IdSec->SetYTitle("Ratio of Ring to Sector Energy");
  SiCD_ERing_IdSec->SetXTitle("Id Sector");
///////
  //DoubleHitSector = mkTH2F("DoubleHitSector","Sector Hit 1 ID vs Hit 2 ID", 40, 0, 40,40, 0, 40);
 // DoubleHitSectorENG = mkTH2F("DoubleHitSectorENG","Sector Hit 1 Energy vs Hit 2 Energy", 5000, 0, 5000, 5000, 0, 5000);
  sicd_e_detid = mkTH2F("sicd_e_detid","Energy of SiCD vs Det ID; Detector ID (ring:1-24, sector:31-62; Energy", 65, 0, 65, 1000,0,100);//
 
 // SiCD_RawERing_ESec = mkTH2F("SiCD_RawERing_ESec", "Gain match: Raw Energy Ring / Cal Energy Sector; Detector ID Ring (ch); RawE_Ring/CalE_Sector",33,-0.5,32.5,binSiRawE,0,200);
  //
  SiCD_XY = mkTH2F("SiCD_XY", "XY image (randomised for each pads) in SiCD plane (Beam direction: #odot); X in GT coordinate (down) /mm; Y in GT coordinate (left when looking downstream) /mm", 1000,-50,50,1000,-50,50);

  diff_SiCDTS = mkTH1D("diff_SiCDTS","Time-Stamp difference btw ring and sector; LED time (ring - sector) / ns; counts", binDiffTS, minDiffTS, maxDiffTS);
  SiCDERatio_diffTS_sector = mkTH2F("SiCDERatio_diffTS_sector", "Energy of sector vs time difference of ring and sector; LED time (ring -sector) /ns; ", binDiffTS, minDiffTS, maxDiffTS, 5000, 0.5, 5000.5);
  SiCDERatio_diffTS_ring = mkTH2F("SiCDERatio_diffTS_ring", "Energy of ring  vs time difference of ring and sector; LED time (ring -sector) /ns;", binDiffTS, minDiffTS, maxDiffTS, 5000, 0.5, 5000.5);

  GT_diffTS_sector = mkTH2F("GT_diffTS_sector", "Energy of GRETINA vs time difference of GT and sector; LED time (GT -sector) /ns; ", binDiffTS, minDiffTS, maxDiffTS, 5000, 0.5, 5000.5);
  GT_diffTS_ring = mkTH2F("GT_diffTS_ring", "Energy of GRETINA  vs time difference of GT and RING; LED time (GT - ring) /ns;", binDiffTS, minDiffTS, maxDiffTS, 5000, 0.5, 5000.5);


  Si_sect_diffTS = mkTH1D("Si_sect_diffTS", "Time-Stamp difference btw sector and GT; TS diff (SiCD_sector - GT)/ns; counts", 5000, -2500, 2500);
  Si_ring_diffTS = mkTH1D("Si_ring_diffTS", "Time-Stamp difference btw ring and GT; TS diff (SiCD_ring - GT)/ns; counts", 5000, -2500, 2500);
  
  //Si_sect_diffTS2 = mkTH1D("Si_sect_diffTS_mode2", "Time-Stamp difference btw sector and GT; TS diff (SiCD_sector - GT)/ns; counts", 1000, -500, 500);
  //Si_ring_diffTS2 = mkTH1D("Si_ring_diffTS_mode2", "Time-Stamp difference btw ring and GT; TS diff (SiCD_ring - GT)/ns; counts", binDiffTS, minDiffTS, maxDiffTS);
 
  DoppEPb_diffTS = mkTH2F("DoppEPb_diffTS", "Doppler corrected Pb Energy vs Time-Stamp difference between GT and SiCD; TS diff (SiCD_sector - GT) / ns; Doppler corrected energy of Pb (keV)", binDiffTS, minDiffTS, maxDiffTS, binGamE2D, minGamE, maxGamE);

  ParticleAngle = mkTH2F("ParticleVTheta", "110Cd Particle Energy (KeV) Vs. Theta Angle (Degrees)",50, 130, 180,20000,0,20000);


  Calculated_S = mkTH1D("ParticleNRGSpec", "Particle Energy Spectrum Calculated", 700, 0, 70);

  LabEnoTG_C = mkTH1D("LabEnoTG_C", Form("energy for Lead Lines; Energy (keV); Counts / %.1f keV",(maxGamE-minGamE)/binGamE),4000, 0, 4000);
  LabEinTG_C = mkTH1D("LabEinTG_C", Form("energy for Lead Lines; Energy (keV); Counts / %.1f keV",(maxGamE-minGamE)/binGamE),binGamE, minGamE, maxGamE);
  LabEoffTG_C = mkTH1D("LabEoffTG_C", Form("energy for Lead Lines; Energy (keV); Counts / %.1f keV",(maxGamE-minGamE)/binGamE),binGamE, minGamE, maxGamE);

  LabEnoTG = mkTH1D("LabEnoTG", Form("energy for 110Cd; Energy (keV); Counts / %.1f keV",(maxGamE-minGamE)/binGamE),binGamE, minGamE, maxGamE);
  LabEinTG = mkTH1D("LabEinTG", Form("energy for 110Cd; Energy (keV); Counts / %.1f keV",(maxGamE-minGamE)/binGamE),binGamE, minGamE, maxGamE);
  LabEoffTG = mkTH1D("LabEoffTG", Form("energy for 110Cd; Energy (keV); Counts / %.1f keV",(maxGamE-minGamE)/binGamE),binGamE, minGamE, maxGamE);


  LabEnoTGF = mkTH1D("LabEnoTG_fine", Form("energy for 110Cd (fine_grain); Energy (keV); Counts / %.1f keV",0.5),6000, 0, 6000);
  LabEinTGF = mkTH1D("LabEinTG_fine", Form("energy for 110Cd (fine_grain); Energy (keV); Counts / %.1f keV",0.5),6000, 0, 6000);
  LabEoffTGF = mkTH1D("LabEoffTG_fine", Form("energy for 110Cd (fine_grain); Energy (keV); Counts / %.1f keV",0.5),6000, 0, 6000);


  doppE_anglePb = mkTH2F("doppE_anglePb","doppE_angle Pb;opening angle cos(theta);enegry (keV)",100,-1,1,binGamE, minGamE, maxGamE);
  doppE_thetaPb = mkTH2F("doppE_thetaPb", "Theta txt evaluated vs. calculated Theta ", 100,130,180,binGamE, minGamE, maxGamE);
  doppE_phiPb = mkTH2F("doppE_phiPb", "doppE_phiPb;phi of C (deg);energy (keV)", 32,0,360,binGamE, minGamE, maxGamE);


  load_kinematics();

    for (int i = 0; i < 24; ++i)
  {
    EnergySpecVANG[i] = mkTH1D(Form("EngVRing_%i",i+1), Form("Energy Spectrum Vs Scattering angle (#theta = %.2f)",Theta_Pb[i]),3000, 0, 3000); ;//mkTH1D(Form("EngVTheta_%d",Theta_Pb[i]), Form("Energy Spectrum Vs Scattering angle (#theta = %d)",Theta_Pb[i]),3000, 0, 3000); 
  }


}

int exit_sicd(void){
	printf("Bye S3....\n");

h1_ene88->Delete();
h1_ene882->Delete();
h2_ene88->Delete();
SiCd_TS_R->Delete();
SiCd_TS_S->Delete();
GTENCAL->Delete();
gMultiplicity->Delete();
GidMul2->Delete();
GTdiff->Delete();
sectVsPhiAngle->Delete();
Qplot->Delete();
ggE->Delete();
ggEBK->Delete();
SiCD_XY->Delete();
diff_SiCDTS->Delete();
SiCDERatio_diffTS_sector->Delete();
SiCDERatio_diffTS_ring->Delete();
GT_diffTS_sector->Delete();
GT_diffTS_ring->Delete();
Si_sect_diffTS->Delete();
Si_ring_diffTS->Delete();
DoppEPb_diffTS->Delete();
ParticleAngle->Delete();
LabEnoTG_C->Delete();
LabEinTG_C->Delete();
LabEoffTG_C->Delete();
LabEnoTG->Delete();
LabEinTG->Delete();
LabEoffTG->Delete();
doppE_anglePb->Delete();
doppE_thetaPb->Delete();
doppE_phiPb->Delete();

  return (0);
}



