// *-- Author : J.Balewski, A.Ogawa, P.Zolnierczuk
// 
// $Id: StEEmcFastMaker.cxx,v 1.1 2003/01/28 23:12:59 balewski Exp $
// $Log: StEEmcFastMaker.cxx,v $
// Revision 1.1  2003/01/28 23:12:59  balewski
// star
//

#include "StChain.h"
#include "St_DataSetIter.h"
#include "StEventTypes.h"

#include "StEEmcFastMaker.h"


#include "StEEmcUtil/EEevent/EEevent.h"
#include "StEEmcUtil/EEevent/EEsectorDst.h"
#include "StEEmcUtil/EEevent/EEmcMCData.h"
#include "StEEmcUtil/EEevent/EEtwHitDst.h"
#include "StEEmcUtil/EEevent/EEsmdHitDst.h"


ClassImp(StEEmcFastMaker)

//--------------------------------------------
StEEmcFastMaker::StEEmcFastMaker(const char *name):StMaker(name){
  mlocalStEvent=0;
  mdbg=1;
  mevIN= new EEmcMCData;
  meeve=new EEevent;
  msamplingFraction=0.05;
  // towers are gain matched to fixed E_T
  const int maxAdc=4095;
  const int maxEtot=60;  // in GeV
  const float feta[kEEmcNumEtas]= {1.95,1.855,1.765,1.675,1.59,1.51,1.435,1.365,1.3,1.235,1.17,1.115}; 
  
  int i;
  
  mfixTgain=new float [kEEmcNumEtas];
  for (i=0;i<kEEmcNumEtas;i++) {
    mfixTgain[i]=maxAdc/maxEtot/cosh(feta[i-1])/msamplingFraction;
  }


  mfixSMDgain=23000;
  mfixPgain=23000;


}

//--------------------------------------------
StEEmcFastMaker::~StEEmcFastMaker(){
 delete  mevIN;
 delete  meeve;
 delete [] mfixTgain;
}

//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
Int_t StEEmcFastMaker::Init(){
  printf("\n\n%s::Init() \n\n",GetName());

  // Create tables
  // Create Histograms    
   return StMaker::Init();
}

//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
Int_t StEEmcFastMaker::Make(){
 
  static int first=1;
  printf("\n\n%s::Make()\n\n",GetName());
  meeve->clear();
  
  int nh=-1;
  if ( (nh = mevIN->readEventFromChain(this)) >0) {
    printf("%s  actual RAW geant EEMC hits =%d nh\n",GetName(),nh);
    if(first && mdbg>1)mevIN->print();
    first=0;
  } else {
    printf("%s  no geant EEMC hits found\n",GetName());
    return kStOK;
  }

  EEevent eeveRaw;    // raw M-C hits 
  
  // generation of TTree
  mevIN->write(&eeveRaw); // Clear & Store RAW EEevent
  if(mdbg) {printf("%s::  raw eeveRaw.print():\n",GetName());eeveRaw.print();}
  
  eeveRaw.sumRawMC(meeve); //sum hits with any detector
  if(mdbg){  printf("%s::  summed eeve.print():\n",GetName());meeve->print();}
   
  StEvent *stevent = mlocalStEvent;
  if(stevent==0) {
    printf("Access full StEvent ...\n");
    stevent =   (StEvent *) (StEvent *) GetInputDS("StEvent");
    assert(stevent); // do sth to provide StEvent first
    printf("check existence of emcCollection... StEvent=%p\n",stevent);
    assert(stevent->emcCollection()); 
    
  }  

  //  SetDumEE(eeve);

  mEE2ST(meeve, stevent);
  
  if(mdbg>2) { // test copying back 
    EEevent eeve2;   // after EE2St and ST2EE
    mST2EE(&eeve2, stevent);  
    printf( "***** before\n");
    meeve->print();
    printf( "***** after\n");
    eeve2.print();
  }

 return kStOK;
}


//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
void  StEEmcFastMaker::SetLocalStEvent(){
  mlocalStEvent = new StEvent;
  StEmcCollection* stemc = new StEmcCollection;
  mlocalStEvent->setEmcCollection(stemc);
}


//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
void  StEEmcFastMaker::mEE2ST(EEevent* eevt, StEvent* stevt){
  int mxSector = kEEmcNumSectors;
  assert(stevt); // fix dumm input
  if(mdbg)printf("EE2ST() start %p\n",stevt);

  StEmcCollection* emcC =(StEmcCollection*)stevt->emcCollection();
   if(mdbg)printf("EE2ST got emcCollection\n");
  
  for(int det = kEndcapEmcTowerId; det<= kEndcapSmdVStripId; det++){
      
    StDetectorId id = StDetectorId(det);
    StEmcDetector* d = new StEmcDetector(id,mxSector);
    emcC->setDetector(d);
    TClonesArray* tca;    
    if(mdbg) printf("EE2ST() copy hits from %d EEMC sectors, det=%d\n",eevt->getNSectors(),det);

    for(int isec=0; isec<mxSector; isec++){ // over used sectors
      int secID=isec+1;
      EEsectorDst* EEsec = (EEsectorDst*)eevt->getSec(secID);
      if(EEsec==0) continue;
      if(mdbg) printf("EE2ST() isec=%d sec_add=%p  secID=%d det=%d\n",isec,EEsec,secID,det);

      switch (det){
      case kEndcapEmcTowerId: //..............................     
	tca = EEsec->getTwHits();      
	for(int j=0; j<=tca->GetLast(); j++){
	  EEtwHitDst* t = (EEtwHitDst *) tca->At(j);
	  int jeta=t->eta()-1;
	  int jsub=t->sub()-'A';

	  // FAST SIMU:
	  int adc=(int) (t->energy() * mfixTgain[jeta]);
	  
	  StEmcRawHit* h = new StEmcRawHit(id,isec,jeta,jsub,adc,t->energy());
	  d->addHit(h);
	  //	  printf("yyy secID=%d, id2=%d\n",isec,h->module());
	  printf("Tw   %c  %d  %f  %d \n",t->sub(),t->eta(),t->energy(),adc);
	} break;
	
      case kEndcapEmcPreShowerId: {//............................
	tca = EEsec->getPre1Hits();      
	for(int j=0; j<=tca->GetLast(); j++){
	  EEtwHitDst* t = (EEtwHitDst *)(* tca)[j];
	  int jeta=t->eta()-1;
	  int jsub=t->sub()-'A';
	  int adc= (int) (t->energy()* mfixPgain);
	  StEmcRawHit* h = new StEmcRawHit(id,isec,jeta,jsub,adc,t->energy());
	  d->addHit(h);
	  printf("Pr1   %c  %d  %d %f\n",t->sub(),t->eta(),adc,t->energy());
	}

	tca = EEsec->getPre2Hits();      
	for(int j=0; j<=tca->GetLast(); j++){
	  EEtwHitDst* t = (EEtwHitDst *)(* tca)[j];
	  int jeta=t->eta()-1;
	  int jsub=t->sub()-'A'+5;
	  int adc= (int) (t->energy()* mfixPgain);
	  StEmcRawHit* h = new StEmcRawHit(id,isec,jeta,jsub,adc,t->energy());
	  d->addHit(h);
	  printf("Pr2   %c  %d  %d %f\n",t->sub(),t->eta(),adc,t->energy());
	}

	tca = EEsec->getPostHits();      
	for(int j=0; j<=tca->GetLast(); j++){
	  EEtwHitDst* t = (EEtwHitDst *)(* tca)[j];
	  int jeta=t->eta()-1;
	  int jsub=t->sub()-'A'+10;
	  int adc= (int) (t->energy()* mfixPgain);
	  StEmcRawHit* h = new StEmcRawHit(id,isec,jeta,jsub,adc,t->energy());
	  d->addHit(h);
	  printf("Post   %c  %d  %d %f\n",t->sub(),t->eta(),adc,t->energy());
	}
      } break;

      case kEndcapSmdUStripId: { //............................
	tca = EEsec->getSmdUHits(); 
	
	for(int j=0; j<=tca->GetLast(); j++){
	  EEsmdHitDst* t = (EEsmdHitDst *)(* tca)[j];
	  int jeta=t->strip()-1;
	  int adc= (int) (t->energy()* mfixPgain);
	  StEmcRawHit* h = new StEmcRawHit(id,isec,jeta,0,adc,t->energy());
	  d->addHit(h);
	  printf("SMDU     %d  %d %f\n",t->strip(),adc,t->energy());
	}
      } break;

      case kEndcapSmdVStripId:  {//............................

	tca = EEsec->getSmdVHits();
	for(int j=0; j<=tca->GetLast(); j++){
	  EEsmdHitDst* t = (EEsmdHitDst *)(* tca)[j];
	  int jeta=t->strip()-1;
	  int adc= (int) (t->energy()* mfixSMDgain);
	  StEmcRawHit* h = new StEmcRawHit(id,isec,jeta,0,adc,t->energy());
	  printf("SMDV    %d  %d  %f\n",t->strip(),adc,t->energy());
	  d->addHit(h);
	} 
      }break;
      default:
	assert(1==2); // the det is out of range, this code gaves up
      }   
    }
  }
}


//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
void  StEEmcFastMaker::mST2EE(EEevent* evt, StEvent* stevt){
  printf("ST2EE: started\n");
  printf("ST2EE: is not matched to EE2ST, fix the code first (J.B.)\n");
  assert(1==2); // not working method, fix it if you need it
  assert(stevt);
  printf("ST2EE:found StEvent\n");
  //  StTpcHitCollection* tpch = (StTpcHitCollection*)stevt->tpcHitCollection();
  StEmcCollection* emcC =(StEmcCollection*)stevt->emcCollection();
  assert(emcC);
  printf("ST2EE:found EmcCollection\n");

  evt->clear();
  for(int det = kEndcapEmcTowerId; det<= kEndcapSmdVStripId; det++){
    printf("ST2EE() det=%d \n",det);

    StDetectorId id = StDetectorId(det);
    StEmcDetector* d = emcC->detector(id);
    if(d==0) {
      printf("ST2EE() Found no detector collection, skipping id=%d\n",id);
      continue;
    }
    printf("ST2EE() det=%d add_d=%p\n",det,d);

    if(d->numberOfModules() < 1) {
      printf("ST2EE() Found no modules in the detector collection, skipping id=%d\n",id);
      continue;
    }
    for(unsigned int isec=0; isec<d->numberOfModules(); isec++){
      StEmcModule* stmod =  d->module(isec);
      if(stmod==0) { printf("ST2EE() couldn't get sector from StEvent, sector=%d\n",isec); continue;}
      StSPtrVecEmcRawHit & h = stmod->hits();
      if(h.size()>0){
	EEsectorDst* sec = evt->getSec((int)isec);
	if(sec==0) { 
	  printf("ST2EE() couldn't find a sector from EEevent, Adding sector=%d\n",isec);
	  sec = evt->addSectorDst(isec);
	}
	switch (det){
	case kEndcapEmcTowerId:     
	  for(unsigned int j=0; j<h.size() ;j++){
	    printf("Tw  %c %d %f\n",h[j]->sub()+'A', h[j]->eta(),h[j]->energy());
	    sec->addTwHit(h[j]->sub()+'A', h[j]->eta(),h[j]->energy());
	  } break;
	case kEndcapEmcPreShowerId: 
	  for(unsigned int j=0; j<h.size() ;j++){
	  int k = h[j]->sub();
	  printf("Pre %d %d %f\n", k, h[j]->eta(),h[j]->energy());
	  if(k<5)       { sec->addPre1Hit(k+'A',    h[j]->eta(),h[j]->energy()); } 
	  else if(k<10) { sec->addPre2Hit(k-5+'A',  h[j]->eta(),h[j]->energy()); } 
	  else          { sec->addPostHit(k-10+'A', h[j]->eta(),h[j]->energy()); } 
	  } break;
	case kEndcapSmdUStripId:
	  for(unsigned int j=0; j<h.size() ;j++){
	    printf("SmU %d %d %f\n",h[j]->eta(),h[j]->sub(),h[j]->energy());
	    sec->addSmdUHit(h[j]->eta(), h[j]->energy());
	  } break;
	case kEndcapSmdVStripId:
	  for(unsigned int j=0; j<h.size() ;j++){
	    printf("SmV %d %d %f\n",h[j]->eta(),h[j]->sub(),h[j]->energy());
	    sec->addSmdVHit(h[j]->eta(), h[j]->energy());
	  } break;
	default:
	  assert(1==2);

	}   
      }
    }
  }
}





