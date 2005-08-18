// *-- Author :Jan Balewski
// 
// $Id: StTrigOnlyPanitkinMaker.cxx,v 1.1 2005/06/17 17:41:13 balewski Exp $

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <StMessMgr.h>

#include "StTrigOnlyPanitkinMaker.h"
#include "StEEmcPool/HanksTriggerDataReader/StTriggerDataReader.h"
#include "StTriggerData2005.h" 
#include "EEdsmAna.h"


ClassImp(StTrigOnlyPanitkinMaker)

//________________________________________________
//________________________________________________
StTrigOnlyPanitkinMaker::StTrigOnlyPanitkinMaker( const char* self ,const char* muDstMakerName) : StMaker(self){
  nTrigEve=nInpEve=0;
  HList=0;
}


//___________________ _____________________________
//________________________________________________
StTrigOnlyPanitkinMaker::~StTrigOnlyPanitkinMaker(){
  delete dsm;
}

//___________________ _____________________________
//________________________________________________
void StTrigOnlyPanitkinMaker::saveHisto(TString fname){
  TString outName=fname+".hist.root";
  TFile f( outName,"recreate");
  assert(f.IsOpen());
  printf("%d histos are written  to '%s' ...\n",HList->GetEntries(),outName.Data());
  HList->Write();
  f.Close();

}
 
//________________________________________________
//________________________________________________
Int_t StTrigOnlyPanitkinMaker::Init(){

  assert(HList);
  dsm=new EEdsmAna(HList,"allTrig");
  dsm->initHisto();
  return StMaker::Init();
}


//________________________________________________
//________________________________________________
Int_t StTrigOnlyPanitkinMaker::Finish(){
  gMessMgr->Message("","I") <<GetName()<<"::Finish() inputEve="<<nInpEve<<" trigFilterEve="<<nTrigEve<<endm;
  return kStOK;
}

//________________________________________________
//________________________________________________
void StTrigOnlyPanitkinMaker::Clear(const Option_t*){
  dsm->clear();
}


//________________________________________________
//________________________________________________
Int_t StTrigOnlyPanitkinMaker::Make(){
  nInpEve++;
  gMessMgr->Message("","D") <<GetName()<<"::Make() is called "<<endm;
  
  StTriggerDataReader* trgreader = (StTriggerDataReader*)GetMaker("TRGRD");
  assert(trgreader);
  void *blob=trgreader->GetTriggerData();

  StTriggerData2005 trgAkio5( (const TrgDataType2005 *)blob);

  //........ filter trigger -none
  nTrigEve++;
  
  const unsigned char * dsm0inp= trgAkio5.getDsm0_EEMC();
  unsigned short int  * dsm1inp= trgAkio5.getDsm1_EEMC();
  unsigned short int  * dsm2inp= trgAkio5.getDsm2_EMC();
  unsigned short int  * dsm3inp= trgAkio5.getDsm3();
  //  trgAkio5->dump();
  dsm->sort( dsm0inp, dsm1inp, dsm2inp, dsm3inp);

  return kStOK;
} 

//---------------------------------------------------
// $Log: StTrigOnlyPanitkinMaker.cxx,v $
// Revision 1.1  2005/06/17 17:41:13  balewski
// *** empty log message ***
//
//
 