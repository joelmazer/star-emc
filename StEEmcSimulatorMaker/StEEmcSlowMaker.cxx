// *-- Author : Hal Spinka
// 
// $Id: StEEmcSlowMaker.cxx,v 1.1 2004/12/15 17:02:56 balewski Exp $

#include <TFile.h>
#include <TH2.h>
#include <TRandom.h>

#include "StEEmcSlowMaker.h"

#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"

#include <StMessMgr.h>

#include "StEEmcDbMaker/EEmcDbItem.h"
#include "StEEmcDbMaker/StEEmcDbMaker.h"


ClassImp(StEEmcSlowMaker)

//________________________________________________
//________________________________________________
StEEmcSlowMaker::StEEmcSlowMaker( const char* self ,const char* muDstMakerName) : StMaker(self){
  mMuDstMaker = (StMuDstMaker*)GetMaker(muDstMakerName);
  assert(mMuDstMaker);
  eeDb=0;
  mHList=0;
  nInpEve=0; 
  //  printf("constr of calib =%s=\n",GetName());
  memset(hA,0,sizeof(hA));

  /// By default, we do not add a pedestal offset
  mAddPed   = false;
  /// By default, we do not smear the pedestal
  mSmearPed = false;
  /// By default, all channels are kept
  mDropBad  = false;
  /// By default, overwrite ADC in muDst
  mOverwrite = true;

}


//___________________ _____________________________
//________________________________________________
StEEmcSlowMaker::~StEEmcSlowMaker(){
}

 
//________________________________________________
//________________________________________________
Int_t StEEmcSlowMaker::Init(){
  eeDb=(StEEmcDbMaker*)GetMaker("eemcDb");
  assert( eeDb);
  if ( mHList ) InitHisto();
  return StMaker::Init();
}

//________________________________________________
//________________________________________________
void StEEmcSlowMaker::InitHisto(){

  char tt1[100], tt2[500];
  TH1F *h;
  TH2F *h2;

  int sectID=5;

  sprintf(tt1,"mm");
  sprintf(tt2,"freq vs. sector ID; sector ID");
  h=new TH1F(tt1,tt2,20,-0.5,19.5);
  hA[12]=h;
  
  sprintf(tt1,"PreADC");
  sprintf(tt2,"Revised ADC for Pre/Post Scint.");
  h=new TH1F(tt1,tt2,100,0,200.);
  hA[0]=h;
  
  sprintf(tt1,"Pre2ADC");
  sprintf(tt2,"Initial ADC for Pre/Post Scint.");
  h=new TH1F(tt1,tt2,100,0,200.);
  hA[1]=h;
  
  sprintf(tt1,"SMDADC");
  sprintf(tt2,"Revised ADC for SMD Strips");
  h=new TH1F(tt1,tt2,100,0,200.);
  hA[2]=h;

  sprintf(tt1,"SMD2ADC");
  sprintf(tt2,"Initial ADC for SMD Strips");
  h=new TH1F(tt1,tt2,100,0,200.);
  hA[3]=h;


  //..................
  sprintf(tt1,"adc%02d",sectID);
  sprintf(tt2," Adc vs, strip ID, sector=%02d; strip ID; ADC ",sectID);
  h2=new TH2F(tt1,tt2,30,0,300,100,0,200.);
  hA[20]=(TH1F*)h2;

  sprintf(tt1,"ADCvsstr");
  sprintf(tt2," Adc vs, strip ID; strip ID; ADC ");
  h2=new TH2F(tt1,tt2,60,0,300,100,0,200.);
  hA[19]=(TH1F*)h2;

  // add histos to the list (if provided)
  int i;
  if(mHList) {
    for(i=0;i<mxH;i++) {
      if(hA[i]==0) continue;
      mHList->Add(hA[i]);
    }
  }
}

//________________________________________________
//________________________________________________
Int_t StEEmcSlowMaker::InitRun(int runNo){
  if(runNo==0) {
    gMessMgr->Message("","W")<<GetName()<<"::InitRun("<<runNo<<") ??? changed to 555, it s OK for M-C - perhaps, JB"<<endm;
    runNo=555;
  }
  return kStOK;
}



//________________________________________________
//________________________________________________
Int_t StEEmcSlowMaker::Finish(){
  return kStOK;
}


//________________________________________________
//________________________________________________
Int_t StEEmcSlowMaker::Make(){
  nInpEve++;

  gMessMgr->Message("","D") <<GetName()<<"::Make() is called , iEve"<<nInpEve<<endm;
  
  /// Access to muDst .......................
  StMuEmcCollection* emc = mMuDstMaker->muDst()->muEmcCollection();
  if (!emc) { 
    gMessMgr->Message("","D") <<"No EMC data for this event"<<endm;    
    return kStOK;
  }
  
  mKSigma = eeDb -> KsigOverPed;

  /// Run slow simulator on pre/postshower
  MakePre(emc);
  /// Run slow simulator on smd
  MakeSMD(emc);
 
  return  kStOK;

}

// ----------------------------------------------------------------------------
void StEEmcSlowMaker::MakePre( StMuEmcCollection *emc )
{


  for ( Int_t i = 0; i < emc->getNEndcapPrsHits(); i++ ) {

    Int_t pre,sec,eta,sub;

    /// muDst ranges: sec:1-12, sub:1-5, eta:1-12 ,pre:1-3==>pre1/pre2/post
    StMuEmcHit *hit=emc->getEndcapPrsHit(i,sec,sub,eta,pre);
    
    /// tmp, for fasted analysis use only hits from sectors init in DB
    if (sec<eeDb->mfirstSecID || sec>eeDb->mlastSecID) continue;
     
    /// Db ranges: sec=1-12,sub=A-E,eta=1-12,type=T,P-R ; slow method
    const EEmcDbItem *x=eeDb-> getTile(sec,sub-1+'A', eta, pre-1+'P'); 
    if(x==0) continue;
  
    /// Drop broken channels
    if(mDropBad) if(x->fail ) continue; 
      
    Float_t adc   = hit -> getAdc();
    Float_t energy= adc / x->gain; // (GeV) Geant energy deposit 
    
    if(adc>0) {

      /// 
      if ( mHList ) 
	hA[1]->Fill(adc);

      /// Addition of Poisson fluctuations and Gaussian 
      /// 1 p.e. resolution
      Float_t mipval = energy / Pmip2ene;
      Float_t avgpe  = mipval * Pmip2pe;
      Float_t Npe    = gRandom->Poisson(avgpe);

      /// Lookup pedestal in database (possibly zero)
      Float_t ped    = (mAddPed)?x -> ped:0;
      /// Smear the pedestal
      Float_t sigped = ( x->thr - ped ) / mKSigma;
      if ( mSmearPed ) ped = gRandom -> Gaus(ped,sigped);
      

      if (Npe>0) {

	/// Determine number of photoelectrons
	Float_t sigmape   = sqrt(Npe) * sig1pe;
	Float_t smearedpe = gRandom -> Gaus(Npe,sigmape);

	/// Determine new ADC value
	Float_t newadc    = smearedpe * mip2ene*(x->gain)/Pmip2pe;
	Int_t   NUadc     = (Int_t)(newadc + 0.5 + ped );

	if ( mHList ) 
	  hA[0]->Fill(NUadc);

	///
	/// If we've made it here, overwrite the muDst
	///
	if ( mOverwrite ) hit -> setAdc( NUadc );

      }

    }

  }

}

// ----------------------------------------------------------------------------
void StEEmcSlowMaker::MakeSMD( StMuEmcCollection *emc ) 
{

  Char_t uv='U';
  for ( uv='U'; uv <= 'V'; uv++ ) {

    Int_t sec,strip;
    for (Int_t i = 0; i < emc->getNEndcapSmdHits(uv); i++ ) {
      
      StMuEmcHit *hit=emc->getEndcapSmdHit(uv,i,sec,strip);
      assert(sec>0 && sec<=MaxSectors);// total corruption of muDst
      
      // tmp, for fasted analysis use only hits from sectors init in DB
      if(sec<eeDb->mfirstSecID || sec>eeDb->mlastSecID) continue;
      
      const EEmcDbItem *x=eeDb->getByStrip(sec,uv,strip);
      assert(x); // it should never happened for muDst

      /// Drop broken channels
      if(mDropBad) {
        if(x->fail ) continue;
      }
    
      Float_t adc    = hit->getAdc();
      Float_t energy = adc/ x->gain; // (GeV) Geant energy deposit 
   
      if(adc>0) {
	
	if ( mHList ) { 
	  hA[12]->Fill(x->sec);
	  if(x->sec==5) ((TH2F*) hA[20])->Fill(x->strip,adc); 
	  hA[3]->Fill(adc);
	}
	
	// Addition of Poisson fluctuations and 
	// Gaussian 1 p.e. resolution
	Float_t mipval = energy / mip2ene;
	Float_t avgpe  = mipval * mip2pe[strip];
	Float_t Npe    = gRandom->Poisson(avgpe);

	/// Lookup pedestal in database (possibly zero)
	Float_t ped    = (mAddPed)?x -> ped:0;
	/// Smear the pedestal
	Float_t sigped = ( x->thr - ped ) / mKSigma;
	if ( mSmearPed ) ped = gRandom -> Gaus(ped,sigped);

	if (Npe>0) {

	  /// Determine number of photoelectrons
	  Float_t sigmape   = sqrt(Npe) * sig1pe;
	  Float_t smearedpe = gRandom -> Gaus(Npe,sigmape);

	  /// Determine new ADC value
	  Float_t newadc    = smearedpe * mip2ene * (x->gain) / mip2pe[strip];
	  Int_t NUadc       = (Int_t)(newadc + 0.5 + ped );

	  if ( mHList ) {
	    ((TH2F*) hA[19]) ->Fill(x->strip,NUadc);
	    hA[2]->Fill(NUadc);
	  }
	  
	  ///
	  /// If we've made it here, overwrite the muDst
	  ///
	  if ( mOverwrite ) hit -> setAdc( NUadc );
	  	  
	}

      }
      
      
    }
  }

}

// $Log: StEEmcSlowMaker.cxx,v $
// Revision 1.1  2004/12/15 17:02:56  balewski
// try 2
//