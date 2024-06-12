#include "bin_sicd.h"
//////////////////////////////////////////////////////////////////////////////
int bin_sicd(GEB_EVENT * GEB_event)
{
//  printf("HERE IN BIN_SICD");
 int GebTypeStr (int type, char strg[]);

 int i, k;
  int nsubev;
  char strg[128];
int SiCdCount=0;

int HitChanRing[8]= {0,0,0,0,0,0,0,0};
float HitEnergyRing[8] = {0,0,0,0,0,0,0,0};

int HitChanSect[8]= {0,0,0,0,0,0,0,0};
float HitEnergySect[8] = {0,0,0,0,0,0,0,0};
 // float GT_offset[150];

  unsigned long long int tsEarly=0;
  unsigned long long int tsLate=0;

  unsigned long long int tsGT[MaxGTNum]={0};
  unsigned long long int tsChico[MaxChicoNum]={0};
  unsigned long long int tsTK=0;
  double Bank29dT;
  int j,g;
  int nGTRaw,nChicoRaw;
  int counter = 0;
  FILE *fp2;
  int i1;
  float r1,r2;

  ////for SiCD detector
  GTEVENT Event;
  unsigned int t1, t2, t3, t4;
  int ii, jj, nBadTestPat = 0, ncrystal;
  volatile unsigned int *bit32Pointer;
  unsigned short int *bit16Pointer;
  char str[128];
  int pos=0;
  unsigned int tempE;
  int rawE;
  int gtheaders=0, newcrystal=0;
  SiCD_struct SiCD_event;

  SiCD_event.LEDts_ring = 1;
  SiCD_event.LEDts_sector = 1;
  SiCD_event.detid_ring = -1;
  SiCD_event.detid_sector = -1;
  SiCD_event.id_ring = -9; // In case some wrong assingnemt, showing hits innner part of the ring;
  SiCD_event.id_sector = -1;
  SiCD_event.rawE_ring = -1;
  SiCD_event.rawE_sector = -1;
  SiCD_event.calE_ring = 1;
  SiCD_event.calE_sector = 1;
  SiCD_event.status = 0;
  SiCD_event.theta_C=1;
  SiCD_event.theta_Pb=1;
  SiCD_event.phi_C=1;
  SiCD_event.phi_Pb=1;
  SiCD_event.beta_C=1;
  SiCD_event.beta_Pb=1;
  
  ////For Doppler reconstruction
  TRACKED_GAMMA_HIT *grh;
  int nTrackedGammas=0, numHitArray1=0, numSiCDHit=0;
  float doppler_factor[MAX_GAMMA_RAYS], doppler_factor_C[MAX_GAMMA_RAYS], doppler_factor_Pb[MAX_GAMMA_RAYS];

  int status;
  unsigned long EventNum=0;
  unsigned long ChicoEventNum=0;
  unsigned long GTEventNum=0;
  unsigned long TKEventNum=0;
  unsigned long long int LEDts1 = 0;
  unsigned long long int LEDts2 = 0;

  long long int tsLastChico, lldt; 
  float dt, dt1, dt2;
  long long int bank88ts;

  nGTRaw = 0;
  nChicoRaw = 0;
  //InitCoinEvent(&CoinEvent);

  nsubev=0;

  bank88ts=0;
  int ch;
 // CoinEvent.nCoinTK=0;
  int nCC = 0;
  double nCCenergies=0;
  int nMode1 = 0;
  int firsttime = 0, t0 = 0, d1;
   bool ringStatus = false;
   bool sectorStatus = false;
  bool doubleHitSector = false;
  bool doubleHitring = false;
  int First_hit_ring=0;
  int First_hit_sector=0;
   SiCD_struct SiCD_Secevent;
    SiCD_Secevent.detid_ring = -90;
    SiCD_Secevent.detid_sector = -90;

  float BigEnergyRing = 0;
  float BigEnergySect = 0;
  double SecHit2Eng;
  int SecHit2Channel;
    //if(GEB_event->mult!=2)return (0);
    for (i = 0; i < GEB_event->mult; i++) {
    


    // BANK88
    /////// for Si detector
    if(GEB_event->ptgd[i]->type != GEB_TYPE_GT_MOD29) continue;
    firsttime++;
    SiCdCount++;
    pos=0;
    ncrystal=0;
    bank88ts=GEB_event->ptgd[i]->timestamp;
    


    ////////// Getting data for Bank 88 for SiCD detector 2022 ///////
    int dummy=1;

    if (Pars.CurEvNo <= Pars.NumToPrint)
      {
      GebTypeStr (GEB_event->ptgd[i]->type, str);
      printf ("\nbin_sicd: %2i> %2i, %s, TS=%lli, 0x%llx; ", i, GEB_event->ptgd[i]->type, str,
        GEB_event->ptgd[i]->timestamp, GEB_event->ptgd[i]->timestamp);
      printf ("payload length: %i bytes", GEB_event->ptgd[i]->length);
      printf ("\n");
      }
    /* byteswap the entire payload */

    bit32Pointer = (unsigned int *) GEB_event->ptinp[i];
    for (int p = 0; p < GEB_event->ptgd[i]->length / 4; p++)
      {
      t1 = (*(bit32Pointer + p) & 0x000000ff) << 24;
      t2 = (*(bit32Pointer + p) & 0x0000ff00) << 8;
      t3 = (*(bit32Pointer + p) & 0x00ff0000) >> 8;
      t4 = (*(bit32Pointer + p) & 0xff000000) >> 24;
      *(bit32Pointer + p) = t1 + t2 + t3 + t4;
      };

    /* inside the payload we have the a number */
    /* of header/trace data sets. These are the  */
    /* crystals that were in coincidence. Here  */
    /* we loop over these header/trace data sets. */
    /* We may have data from more than one crystal  */
    /* in this payload. Thus, there will be at  */
    /* least 40 traces, but there can also be  */
    /* 80, 120... etc. Be sure MAXPAYLOADSIZE  */
    /* is big enough to handle GT mode3 payloads. */

    /* loop over the header/traces of the mode 3 data */
    
    while (pos < GEB_event->ptgd[i]->length)
    {
      /* start of event (Event.len known from last event) */
      if (pos == 0)
        bit32Pointer = (unsigned int *) GEB_event->ptinp[i];
      else
        bit32Pointer += (Event.len / 4);

      /* check the EOE situation */

      if (*bit32Pointer != EOE)
        {
        nBadTestPat++;

        if (nBadTestPat == 10)
          printf ("** suspending warnings about bad EOE markers...\n");

        if (nBadTestPat < 10)
          {
          printf ("ooops: bit32Pointer=%8.8x after event # %i trace # %i\n", *bit32Pointer, Pars.CurEvNo,
              ncrystal);
          fflush (stdout);
          exit (1);
          };

        };

      if ( (gtheaders%40) == 0)
        {
        ncrystal++;
        newcrystal=1;
        }
      gtheaders++;
      //             printf("%i %i %i\n", gtheaders, ncrystal, newcrystal);

      /* fill event header, skip EOE */

      bit16Pointer = (unsigned short int *) (bit32Pointer + 1);
      for (j = 0; j < HDRLENWORDS; j++)
        {
        Event.hdr[j] = *(bit16Pointer + j);
        };

      /* debug list */
      /*
      if (Pars.CurEvNo <= Pars.NumToPrint)
        {
        printf("\n---------\nnew GTheader at pos=%i or %i words (%i):\n",pos,pos/2, pos/2+8);
        for (j = 0; j < HDRLENWORDS; j++)
          printf ("Event.hdr[%2i]=0x%4.4x, %6i\n", j, Event.hdr[j], Event.hdr[j]);
          
        };*/
      /* eventlength/tracelength in bytes */
      /* the +4 comes from the EOE 4 bit word */

      Event.len = 4 * (Event.hdr[1] & 0x7ff) + 4;
      Event.traceLen = Event.len - HDRLENBYTES - sizeof (unsigned int);

      if (Pars.CurEvNo <=0)
        {
        printf ("Event.len=%4i, Event.traceLen=%4i (in Bytes)\n", Event.len, Event.traceLen);
        printf ("Event.len=%4i, Event.traceLen=%4i (in 16 bit words)\n", Event.len / sizeof (short int),
        Event.traceLen / sizeof (short int));
        printf ("Event.len=%4i, Event.traceLen=%4i (in 32 bit words)\n", Event.len / sizeof (int),
        Event.traceLen / sizeof (int));
        };


      /* check the next EOE is there */
      /* before we go on (unless last set) */

      pos += Event.len;

      if (Pars.CurEvNo <= Pars.NumToPrint && pos < GEB_event->ptgd[i]->length)
        printf ("next start: 0x%8.8x, pos= %i\n", *(bit32Pointer + (Event.len / 4)), pos);

      /* potential debug info */

      if (*(bit32Pointer + (Event.len / 4)) != EOE && pos < GEB_event->ptgd[i]->length)
        {
        printf ("next EOE not found at pointer %p, seek manually\n", bit32Pointer + (Event.len / 4));
        printf ("Event.len=%i, Event.traceLen=%i\n", Event.len, Event.traceLen);
        int length = 1;
        while (length < (Event.len + 10))
          {
          printf ("%3i has 0x%8.8x\n", 4 * length, *(bit32Pointer + length));
          fflush (stdout);
          if (*(bit32Pointer + length) == EOE)
            {
            printf ("found it\n");
            break;
            };
          length++;
          };
        length--;
        printf ("%3i has 0x%8.8x\n", length, *(bit32Pointer + length));
        };

      /* extract the Board IDs etc */
      /* chan_id is 0-9 on digitizer */

      Event.chan_id = (Event.hdr[0] & 0x000f);
      Event.digitizer_id = (Event.hdr[0] >> 4) & 0x0003;
      Event.crystal_id = (Event.hdr[0] >> 6) & 0x0003;
      Event.module_id = (Event.hdr[0] >> 8) & 0x00ff;  //was 1f  /* zhu use 0x00ff */

      /* construct detector/channel/segment ID */

      Event.ge_id = Event.crystal_id;      

      Event.channel = Event.digitizer_id * 10 + Event.chan_id;

      
      Event.seg_id = Event.ge_id * 50 + Event.channel;

     
      

        /* count the crystals we have processed */
        if (Pars.CurEvNo <= Pars.NumToPrint)
        if (newcrystal==1)
          {
            newcrystal=0;
          };

      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf ("ncrystal=%i\n", ncrystal);


      tempE = (((unsigned int) Event.hdr[6]) & 0x00ff) << 16;
      tempE += (unsigned int) Event.hdr[5] & 0xffff;
      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf("two's comp tempE  :: ", tempE);

      rawE = twoscomp_to_int_24 (tempE);
      //if (Pars.CurEvNo <= Pars.NumToPrint)
      //  printf ("rawE              :: ", rawE);

      /* change the sign as the the Digs are neg */

      if (rawE <= 0)
        rawE = -rawE;
      if (Pars.CurEvNo <= Pars.NumToPrint)
       // printf ("rawE (sign change):: ", rawE);

      /* downscale */
      if(Event.ge_id==0) 
      {
        if(BigEnergyRing<((float) rawE )) 
        {
          BigEnergyRing = ((float) rawE );
        }
      }

      if(Event.ge_id==3) 
      {
        if(BigEnergySect<((float) rawE )) 
        {
          BigEnergySect = ((float) rawE );
        }
      }

//      Event.LEDts = (unsigned long long int) Event.hdr[2] + ((unsigned long long int) Event.hdr[3] << 16) + ((unsigned long long int) Event.hdr[4] << 32);
                         // printf("\nIT IS Event.ge_id: %i",Event.ge_id);

          numSiCDHit++;

      if((float) rawE < 1000.0) continue;

      if(Event.ge_id==0){
        HitChanRing[i] = Event.channel;
        HitEnergyRing[i] = ((float) rawE );
      }

       if(Event.ge_id==3){
        HitChanSect[i] = Event.channel;
        HitEnergySect[i] = ((float) rawE );
      }

      if( BigEnergyRing <=((float) rawE ) && Event.ge_id==0){

          if(Event.ge_id == 0) // For the case of ring
          { 

              if(First_hit_ring<i && ringStatus == true && doubleHitring ==false )
              {
                SiCD_Secevent.LEDts_ring = static_cast<Double_t>((GEB_event->ptgd[i]->timestamp/2)+gRandom->Uniform())*10;
                 SiCD_Secevent.id_ring = Event.channel;
                  SiCD_Secevent.rawE_ring = ((float) rawE);
                  doubleHitring = true;
              }

              if(ringStatus == false) First_hit_ring = i;

              SiCD_event.LEDts_ring = static_cast<Double_t>((GEB_event->ptgd[i]->timestamp/2)+gRandom->Uniform())*10;
              SiCD_event.id_ring = Event.channel;
              SiCD_event.rawE_ring = ((float) rawE);

              if(ringStatus == false)ringStatus = true;
          } 

            if (Pars.CurEvNo <= 0)
            {
              printf("bin_sicd RING PARAMETERS: \n");
              printf ("0x%4.4x ", Event.hdr[0]);
              printf ("\nmodule_id = %i ", Event.module_id);
              printf ("\ncrystal_id = %i ", Event.crystal_id);
              printf ("\ndigitizer_id = %i ", Event.digitizer_id);
              printf ("\nchan_id = %i ", Event.chan_id);
              printf ("\nge_id = %i ", Event.ge_id);
              printf ("\nchannel = %i ", Event.channel);
              printf("\nEnergy of SiCd is: %f\n",((float) rawE));
              printf ("\n Time Stamp of hit is = %f ", static_cast<Double_t>((GEB_event->ptgd[i]->timestamp/2)));
              printf ("\n");
            };


        }

      if( BigEnergySect <=((float) rawE ) && Event.ge_id==3){

         if(Event.ge_id == 3 ) // redefined sectors 
          { 


              if(First_hit_sector<i && sectorStatus == true && doubleHitSector ==false)
              {
                SiCD_Secevent.LEDts_sector = static_cast<Double_t>((GEB_event->ptgd[i]->timestamp/2)+gRandom->Uniform())*10;
                 SiCD_Secevent.id_sector = Event.channel;
                  SiCD_Secevent.rawE_sector = ((float) rawE);
                  doubleHitSector = true;

              }

                if(sectorStatus == false) First_hit_sector = i;

                SiCD_event.LEDts_sector = static_cast<Double_t>((GEB_event->ptgd[i]->timestamp/2)+gRandom->Uniform())*10;
                SiCD_event.id_sector = Event.channel;
                SiCD_event.rawE_sector = ((float) rawE);
                if(sectorStatus == false)sectorStatus = true;

          }


       
        } 
          
        }


      
      }  // End for Si hit reconstruction

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  int value_ink = map_SiCD(&SiCD_event); // Convert digitizer channel to id and theta, phi, beta
 

  
  if(SiCD_event.detid_ring == -9 || SiCD_event.detid_sector == -1 || SiCD_event.detid_ring == -1 || SiCD_event.detid_ring == 1 || SiCD_event.detid_ring == 2) {
      
    return (0);
  }

  int truechanid;
  
  if(SiCD_event.id_sector == 0 || SiCD_event.id_sector == 20) return (0); //There is no 0 or 20 but they fill up with anamolous hits.

  //---------------[SiCd Particle Energy Cut]--------------------------|
  int gammaCount=0;
  for (i = 0; i < GEB_event->mult; i++)
    {
          
          if (GEB_event->ptgd[i]->type != GEB_TYPE_TRACK) continue;
          gammaCount++;

    }


  //---------------[SiCd Particle Energy Cut]--------------------------|
  int id_ring = SiCD_event.detid_ring;
  int id_sector = SiCD_event.detid_sector;



  //---------------[Target Like Particle]--------------------------|

  double theta_C=SiCD_event.theta_C, phi_C=SiCD_event.phi_C;

  double degtheta_C=TMath::RadToDeg()*SiCD_event.theta_C;

  double degphi_C=TMath::RadToDeg()*SiCD_event.phi_C;

  double beta_C=SiCD_event.beta_C;
  
  //---------------[Target Like Particle]--------------------------|


  double SectorTotalEnergy=0;
  double RingTotalEnergy=0;

  if (((SiCD_event.calE_sector>=5000 && SiCD_event.calE_sector<=13000)&&(SiCD_event.calE_ring>=6000 && SiCD_event.calE_ring<=13000)))
  {
   for (int i = 0; i < 8; ++i)
    {
      SectorTotalEnergy += HitEnergySect[i];
      RingTotalEnergy += HitEnergyRing[i]; 
    }

     for (int i = 0; i < 8; ++i)
    {
      if(SiCdCount<=i) continue;

      RingPlot->Fill(SiCdCount,HitChanRing[0]-HitChanRing[i],(HitEnergyRing[i]/(RingTotalEnergy))*100.0);
      RingPlot->Fill(SiCdCount,HitChanSect[0]-HitChanSect[i],(HitEnergySect[i]/(SectorTotalEnergy))*100.0);

      if(SiCdCount==2){
      SiCdMul2->Fill(HitChanRing[0]-HitChanRing[i],(HitEnergyRing[i]/(RingTotalEnergy))*100.0);
      SiCdMul2->Fill(HitChanSect[0]-HitChanSect[i],(HitEnergySect[i]/(SectorTotalEnergy))*100.0);
    }

    if(SiCdCount==3){
      SiCdMul3->Fill(HitChanRing[0]-HitChanRing[i],(HitEnergyRing[i]/(RingTotalEnergy))*100.0);
      SiCdMul3->Fill(HitChanSect[0]-HitChanSect[i],(HitEnergySect[i]/(SectorTotalEnergy))*100.0);
    }

    if(SiCdCount==4){
      SiCdMul4->Fill(HitChanRing[0]-HitChanRing[i],(HitEnergyRing[i]/(RingTotalEnergy))*100.0);
      SiCdMul4->Fill(HitChanSect[0]-HitChanSect[i],(HitEnergySect[i]/(SectorTotalEnergy))*100.0);
    }

    }

  }
  //---------------[Beam Like Particle]--------------------------|

  double beta_Pb=SiCD_event.beta_Pb;

  double theta_Pb=SiCD_event.theta_Pb, phi_Pb=SiCD_event.phi_Pb;

  double degtheta_Pb=TMath::RadToDeg()*SiCD_event.theta_Pb;

  double degphi_Pb=TMath::RadToDeg()*SiCD_event.phi_Pb;

  //---------------[Beam Like Particle]--------------------------|

  //---------------[SiCd Particle Energy and TIMING Histograms]--------------------------|


  double Si_tdiff = SiCD_event.LEDts_ring - SiCD_event.LEDts_sector; 
  diff_SiCDTS->Fill(Si_tdiff);


  SiCDERatio_diffTS_sector->Fill(Si_tdiff, SiCD_event.calE_sector);
  SiCDERatio_diffTS_ring->Fill(Si_tdiff, SiCD_event.calE_ring);
  
 
  // To confirm ring mapping;
  sicd_e_detid -> Fill(id_ring, SiCD_event.calE_ring);
  sicd_e_detid -> Fill(id_sector + 30, SiCD_event.calE_sector);

  

  // Ecal rings
   SiCD_RawERing_ESec->Fill(id_ring, (double)SiCD_event.rawE_ring/SiCD_event.calE_sector);


  
   if(gammaCount<1) return 0;

   if(SiCdCount == 2)
  {
    SectorMultiEng->Fill((SiCD_event.calE_sector/1000.0)*4.937);
  }

  //---------------[SiCd Particle Energy and TIMING Histograms]--------------------------|
   h1_ene88->Fill((SiCD_event.calE_sector/1000.0)*4.937);
   h1_ene882->Fill((SiCD_event.calE_ring/1000.0)*4.937);
  //---------------[Histogram Build for X,Y plot of SiCD Hits]--------------------------|

   if((SiCD_event.calE_ring/1000.0)*4.937>=25 && (SiCD_event.calE_ring/1000.0)*4.937 <=35){
    gMultiplicity->Fill(SiCdCount);
   }

    double R_Si = 11. + ((double)id_ring - 0.5)*1.0 + gRandom->Uniform(-0.50,0.50); //11 mm is the smallest inner radius of disk.
                    
    double Phi_Si = ((degphi_Pb - gRandom->Uniform(-1.0,1.0) * 5.625)) * TMath::DegToRad();  
    Phi_Si = -Phi_Si; //by Sici

    double X_Si = R_Si * TMath::Cos(Phi_Si)+0.274;
    double Y_Si = R_Si * TMath::Sin(Phi_Si)-4.608;
                   
    Phi_Si = atan2f(Y_Si,X_Si);
    double r_sph = sqrt(pow(X_Si,2)+pow(Y_Si,2)+pow(40.0,2));
                                    
    double degtheta_Pb_calculated=TMath::Pi()-(acosf(40.0/(r_sph)));
      
    SiCdAngleMap->Fill(Phi_Si*TMath::RadToDeg(),degtheta_Pb_calculated*TMath::RadToDeg());

  //---------------[For Phi Correction from off-centered beam]--------------------------|

    TVector3 v1(R_Si * TMath::Cos(Phi_Si+TMath::Pi()),R_Si * TMath::Sin(Phi_Si)+TMath::Pi(),40.0); //Original Position 2 180 degrees from 1.
    TVector3 v2(R_Si * TMath::Cos(Phi_Si),R_Si * TMath::Sin(Phi_Si),40.0); //Original Position 1
    TVector3 v3(-0.274,+4.608,0); //Shift of Beam Center

    TVector3 v4 = v3 - v2;
    TVector3 v5 = v3 - v1;

    double R1 = sqrt(pow(v4.X(),2)+pow(v4.Y(),2)+pow(v4.Z(),2)); // the counts to be corrected..
    double R2 = sqrt(pow(v5.X(),2)+pow(v5.Y(),2)+pow(v5.Z(),2));

    double ShiftC = pow(R1,2)/pow(R2,2);

    if(Pars.CurEvNo <= 1000 && id_ring == 4) output2 << degphi_Pb << "\t" << R1 << "\t" << R2 << "\t" <<degtheta_Pb_calculated*TMath::RadToDeg()<< "\t" << ShiftC << endl;

   h2_ene88->Fill(Phi_Si*TMath::RadToDeg(),SiCD_event.calE_sector);
  //---------------[Histogram Build for X,Y plot of SiCD Hits]--------------------------|

   //---------------[Kinematics Fit Result Function Parameters]--------------------------|
   // Fit function and txt file for fit here: /home/sangeetpannu/Phd_110Cd_Project/KinematicsFits/
   double A = 3.8432e-9;
   double B = -3.14677e-6;
   double C = 0.000992028;
   double D = -0.141381;
   double E = 8.02178;

   //EBeam is in MeV/A so we multiply by the beam's atomic number.
   double Ebeam = (A*pow((degtheta_Pb_calculated)*TMath::RadToDeg(),4) + B*pow((degtheta_Pb_calculated)*TMath::RadToDeg(),3) + C*pow((degtheta_Pb_calculated)*TMath::RadToDeg(),2)+ D*pow((degtheta_Pb_calculated*TMath::RadToDeg()),1) + E )*109.9030021; 
    
   //---------------[Kinematics Fit Result Function Parameters]--------------------------|

   //---------------[SRIM Fit Calculations]--------------------------|

   double a1 = 0.3004;
   double a2 = -0.0026203;
   double a3 = 1.09e-05;
   double a4 = 1.5719;

   double SRIMCalENG = a4 + a1*pow(Ebeam,1)+a2*pow(Ebeam,2)+a3*pow(Ebeam,3);
   double truebeamENG = Ebeam - SRIMCalENG*(0.25/cosf(TMath::Pi()-degtheta_Pb_calculated));

    ParticleAngle->Fill(degtheta_Pb_calculated*TMath::RadToDeg(),truebeamENG*1000);
    Qplot->Fill(degtheta_Pb_calculated*TMath::RadToDeg(),(truebeamENG)-(SiCD_event.calE_ring/1000.0)*4.937);

    Calculated_S->Fill((truebeamENG));

    if(!((abs(Si_tdiff)<20))) return 0;

    
    SiCD_XY -> Fill(X_Si,Y_Si);
    sectVsPhiAngle->Fill(-Phi_Si*TMath::RadToDeg());

   double TS_adj = 0;
   if(Pars.Run_number>=43)TS_adj = 70.0;
   
   //---------------[SRIM Fit Calculations]--------------------------|

  ////////////////////////////////////////////////////////////////////
  //////////             GRETINA reconstructions             /////////
  ////////////////////////////////////////////////////////////////////

  for (i = 0; i < GEB_event->mult; i++)
    {
          
          if (GEB_event->ptgd[i]->type != GEB_TYPE_TRACK) continue;

          grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];
            
          if (grh->ngam < Pars.multlo || grh->ngam > Pars.multhi) return (0);
         
          if(grh->ngam==2)
          {
            GidMul2->Fill(grh->gr[0].fhcrID,grh->gr[1].fhcrID);
          }

            for (j = 0; j < grh->ngam; j++)
            {
             
              if (grh->gr[j].tracked)
                {  

                if((grh->gr[j].fhcrID==109)) continue;
                 if(grh->gr[j].fhcrID==117) continue;


                  //-----------------[GRETINA-SiCD TIMING Variables]---------------------|
                  /*


                    X The divison by two was an adjustment to the time stamps which were seen to be spaced by
                    units of 2. 

                    X The addition of gRandom is for the "smoothing" of the timestamps from "long long int" to "double". The
                    addition is completed between -1 to 1 randomly and this does not skew the data or add anomolous effects 
                    as the DIGITIZERS have a resolution of 10ns (i.e the hit timing is taken every 10ns intervals).

                  */
                  

                  double gtsi_diff = static_cast<Double_t>((grh->gr[j].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_sector-TS_adj;             
                  double gtri_diff = static_cast<Double_t>((grh->gr[j].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_ring-TS_adj;

                  Si_sect_diffTS->Fill(gtsi_diff);
                  Si_ring_diffTS->Fill(gtri_diff);

                  //-----------------[GRETINA-SiCD TIMING Variables]---------------------|
                 

                  //==================================================[Doppler Correction Parameters for 110Cd Particles]==================================================|

                  float gamma_det_radius = sqrt(pow(grh->gr[j].x0-Pars.target_x,2.)+pow(grh->gr[j].y0-Pars.target_y,2.)+pow(grh->gr[j].z0-Pars.target_z,2.));
                  float theta_G = acosf((grh->gr[j].z0-Pars.target_z)/ gamma_det_radius);
                  float phi_G = atan2f((grh->gr[j].y0),(grh->gr[j].x0));
                
                  float cos_GPb =  cos_Gamma_Recoil(degtheta_Pb_calculated,  Phi_Si,  theta_G,  phi_G); // MAYBE Fix the particle theta and phi.

                  double betaCalculated = sqrt((2*(SiCD_event.calE_ring/1000.0)*4.937)/(109.9030021*931.49410242));

                  float gamma_sq_Pb = 1.0 - (betaCalculated * betaCalculated);

                  doppler_factor_Pb[j] = (sqrt (gamma_sq_Pb))/((1.0 - (betaCalculated * cos_GPb)));

                  //==================================================[Doppler Correction Parameters for 110Cd Particles]==================================================|


                  //==================================================[Doppler Correction Applications to Corrections to Data]==================================================|

/*
                  double delE=4.0;
                  double delx=40000;

                
                  while(delE<=4.5 && (id_ring == 5 || id_ring == 4))
                      {

                        double betaCalculated_adj = sqrt((2*((SiCD_event.calE_ring/1000.0)*delE))/(109.9030021*931.49410242));

                        float gamma_sq_Pb_adj = 1.0 - (betaCalculated_adj * betaCalculated_adj);

                        float doppler_factor_adj = (sqrt (gamma_sq_Pb_adj))/((1.0 - (betaCalculated_adj * cos_GPb)));
                        //if(gamma_sq_Pb_adj>=0 && grh->gr[j].fhcrID==119 && id_ring == 5)OptimalE5->Fill(delE,((grh->gr[j].esum/doppler_factor_adj)*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                        if(gamma_sq_Pb_adj>=0 && grh->gr[j].fhcrID==119 && id_ring == 4)OptimalE4->Fill(delE,((grh->gr[j].esum/doppler_factor_adj)*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                        

                        ++delx;
                        delE = delx/10000.0;
                      }
                  */

          
                //==================================================[Doppler Correction Applications to Corrections to Data]==================================================|

               //===================================================[208-Pb Spectrum]===================================================|        

                  float gamma_sq_C = 1.0 - (beta_C * beta_C);

                  float cos_GC =  cos_Gamma_Recoil(SiCD_event.theta_C, SiCD_event.phi_C,  theta_G,  phi_G); // MAYBE Fix the particle theta and phi.

                  doppler_factor_C[j] = (sqrt (gamma_sq_C))/((1.0 - (beta_C * cos_GC)));

                  LabEnoTG_C->Fill(((grh->gr[j].esum/doppler_factor_C[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));

                  if((gtsi_diff>-20 && gtsi_diff<320)&& (gtri_diff>-20 && gtri_diff<320)) LabEinTG_C->Fill(((grh->gr[j].esum/doppler_factor_C[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                  if((abs(gtsi_diff)>400 && abs(gtsi_diff)<550)&& (abs(gtri_diff)>400 && abs(gtri_diff<550))) LabEoffTG_C->Fill(((grh->gr[j].esum/doppler_factor_C[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));


                //===================================================[208-Pb Spectrum]===================================================| 


                //===================================================[110Cd Particle Energy Gated Spectrum]===================================================|
                 

                  if(((SiCD_event.calE_sector>=5000 && SiCD_event.calE_sector<=13000)&&(SiCD_event.calE_ring>=6000 && SiCD_event.calE_ring<=13000)))
                  {
                    GTENCAL->Fill(grh->gr[j].fhcrID,((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    
                    GT_diffTS_sector->Fill(gtsi_diff,((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    GT_diffTS_ring->Fill(gtri_diff,((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    
                    //===================================================[110Cd Particle Time Gated Spectrum]===================================================|

                    LabEnoTG->Fill(((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    if((gtsi_diff>-20 && gtsi_diff<320)&& (gtri_diff>-20 && gtri_diff<320)) LabEinTG->Fill(((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID])); 
                    if((abs(gtsi_diff)>400 && abs(gtsi_diff)<550)&& (abs(gtri_diff)>400 && abs(gtri_diff<550))) LabEoffTG->Fill((((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID])));
                    
                    LabEnoTGF->Fill(((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    if((gtsi_diff>-20 && gtsi_diff<320)&& (gtri_diff>-20 && gtri_diff<320)) LabEinTGF->Fill(((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID])); 
                    if((abs(gtsi_diff)>400 && abs(gtsi_diff)<550)&& (abs(gtri_diff)>400 && abs(gtri_diff<550))) LabEoffTGF->Fill((((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID])));
                    

                   //===================================================[110Cd Particle Time Gated Spectrum]===================================================|
                    
                    doppE_anglePb->Fill(cos_GPb, ((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    doppE_thetaPb->Fill(TMath::RadToDeg()*(degtheta_Pb_calculated),((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                    doppE_phiPb->Fill(-Phi_Si*TMath::RadToDeg(),((grh->gr[j].esum/doppler_factor_Pb[j])*Pars.CCcal_gain[grh->gr[j].fhcrID]+Pars.CCcal_offset[grh->gr[j].fhcrID]));
                  }


                //===================================================[110Cd Particle Energy Gated Spectrum]===================================================|
                  
                  
                }
            }

            //===================================================[GRETINA gamma-gamma Energy Spectrum GATE ON]===================================================|
            for (int k = 0; k < grh->ngam; k++)
            {  
          
               double gtsi_diff = static_cast<Double_t>((grh->gr[k].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_sector-TS_adj;             
               double gtri_diff = static_cast<Double_t>((grh->gr[k].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_ring-TS_adj;

               if(!((gtsi_diff>-20 && gtsi_diff<320)&&(gtri_diff>-20 && gtri_diff<320))) continue;

                          if (grh->gr[k].tracked)
                            {
                              for (int l = k + 1; l < grh->ngam; l++)
                              {
                             
                                double gtsi_diff2 = static_cast<Double_t>((grh->gr[l].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_sector-TS_adj;             
                                double gtri_diff2 = static_cast<Double_t>((grh->gr[l].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_ring-TS_adj;
                                if(!((gtsi_diff2>-20 && gtsi_diff2<320)&&(gtri_diff2>-20 && gtri_diff2<320))) continue;

                                if (grh->gr[l].tracked)
                                    {
                                              
                                            if(grh->gr[k].fhcrID==109) continue; //THIS CRYSTAL WAS BAD!
                                            if(grh->gr[k].fhcrID==117) continue;
            
                                            if(grh->gr[l].fhcrID==109) continue; //THIS CRYSTAL WAS BAD!
                                            if(grh->gr[l].fhcrID==117) continue;
                                      
                                            //===================================================[110Cd Particle Energy Gated Spectrum]===================================================|
              
                                            if(((SiCD_event.calE_sector>=5000 && SiCD_event.calE_sector<=13000)&&(SiCD_event.calE_ring>=6000 && SiCD_event.calE_ring<=13000)))
                                            {

                                              ggE->Fill(((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]),\
                                                ((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]));
                                              ggE->Fill(((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]),\
                                                ((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]));
                                              

                                              ggEF->Fill(((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]),\
                                                ((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]));
                                              ggEF->Fill(((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]),\
                                                ((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]));
        
                                            }

                                           //===================================================[110Cd Particle Energy Gated Spectrum]===================================================|
             
                                          
                                              
            
                                    }
                                }
            
                            }
            
            }

            //===================================================[GRETINA gamma-gamma Energy Spectrum GATE ON]===================================================|




            //===================================================[GRETINA gamma-gamma Energy Spectrum GATE OFF]===================================================|

            for (int k = 0; k < grh->ngam; k++)
            {  
             
               double gtsi_diff = abs(static_cast<Double_t>((grh->gr[k].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_sector-TS_adj);             
               double gtri_diff = abs(static_cast<Double_t>((grh->gr[k].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_ring-TS_adj);

               if(!(gtsi_diff>400)&&(gtsi_diff<550) && (gtri_diff<400 && gtri_diff<550)) continue;
                          if (grh->gr[k].tracked)
                            {
                              for (int l = k + 1; l < grh->ngam; l++)
                              {
                               
                                double gtsi_diff2 = abs(static_cast<Double_t>((grh->gr[l].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_sector-TS_adj);             
                                double gtri_diff2 = abs(static_cast<Double_t>((grh->gr[l].timestamp/2)+gRandom->Uniform())*10 - SiCD_event.LEDts_ring-TS_adj);

                                if(!(gtsi_diff2>400)&&(gtsi_diff2<550) && (gtri_diff2<400 && gtri_diff2<550)) continue;

                                if (grh->gr[l].tracked)
                                    {
                                            
                                      if(grh->gr[k].fhcrID==109) continue; //THIS CRYSTAL WAS BAD!
                                      if(grh->gr[k].fhcrID==117) continue; //THIS CRYSTAL WAS BAD!
                                            
                                      if(grh->gr[l].fhcrID==109) continue; //THIS CRYSTAL WAS BAD!
                                      if(grh->gr[l].fhcrID==117) continue; //THIS CRYSTAL WAS BAD!
                                            
                                         
                                      //===================================================[110Cd Particle Energy Gated Spectrum]===================================================|
              
                                      if(((SiCD_event.calE_sector>=5000 && SiCD_event.calE_sector<=13000)&&(SiCD_event.calE_ring>=6000 && SiCD_event.calE_ring<=13000)))
                                      {
                                        ggEBK->Fill(((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]),\
                                          ((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]));
                                        ggEBK->Fill(((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]),\
                                          ((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]));    
                                        

                                        ggEBKF->Fill(((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]),\
                                          ((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]));
                                        ggEBKF->Fill(((grh->gr[l].esum/doppler_factor_Pb[l])*Pars.CCcal_gain[grh->gr[l].fhcrID]+Pars.CCcal_offset[grh->gr[l].fhcrID]),\
                                          ((grh->gr[k].esum/doppler_factor_Pb[k])*Pars.CCcal_gain[grh->gr[k].fhcrID]+Pars.CCcal_offset[grh->gr[k].fhcrID]));    
                                          
                                      }
                                      //===================================================[110Cd Particle Energy Gated Spectrum]===================================================|

                                    }
                                }
            
                            }
            
            }


            //===================================================[GRETINA gamma-gamma Energy Spectrum GATE OFF]===================================================|
            

    }


  return (0);
  
}
