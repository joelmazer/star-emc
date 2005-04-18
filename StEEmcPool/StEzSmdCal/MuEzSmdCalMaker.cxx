// *-- Author : Jan Balewski
// 
// $Id: MuEzSmdCalMaker.cxx,v 1.1 2005/03/11 15:44:25 balewski Exp $

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <StMessMgr.h>

#include "MuEzSmdCalMaker.h"

#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"

#include "StMuDSTMaker/EZTREE/EztEventHeader.h"
#include "StMuDSTMaker/EZTREE/EztTrigBlob.h"
#include "StMuDSTMaker/EZTREE/EztEmcRawData.h"
#include "StMuDSTMaker/EZTREE/StTriggerDataMother.h"
//tmp
#include "StTriggerData2005.h" // tmp

#include "StEEmcDbMaker/StEEmcDbMaker.h"
#include "StEEmcDbMaker/EEmcDbItem.h"

ClassImp(MuEzSmdCalMaker)

//________________________________________________
//________________________________________________
MuEzSmdCalMaker::MuEzSmdCalMaker( const char* self ,const char* muDstMakerName) : StMaker(self){
  mMuDstMaker = (StMuDstMaker*)GetMaker(muDstMakerName);
  assert(mMuDstMaker);

  trgAkio=0;
  nAcceptEve=nTrigEve=nInpEve=nCorrEve=0;
  SetHList(0);
  SetTrigIdFilter(0);
  SetMaxCtbSum(0);

}

//________________________________________________
//________________________________________________
void MuEzSmdCalMaker::SetSector(int sec){
  setSector(sec);
  TString name=GetName();
  name+="-";
  name+=sec;
  //  printf("change name to %s sec=%d\n", name.Data(),sec);
  SetName(name);
}

//___________________ _____________________________
//________________________________________________
MuEzSmdCalMaker::~MuEzSmdCalMaker(){
  delete trgAkio;
}

//___________________ _____________________________
//________________________________________________
void 
MuEzSmdCalMaker::saveHisto(TString fname){
  TString outName=fname+".hist.root";
  TFile f( outName,"recreate");
  assert(f.IsOpen());
  printf("%d histos are written  to '%s' ...\n",HList->GetEntries(),outName.Data());

  HList->Write();
  f.Close();

}
 
//________________________________________________
//________________________________________________
Int_t 
MuEzSmdCalMaker::Init(){

  assert(HList);
  eeDb=(StEEmcDbMaker*)GetMaker("eemcDb");
  assert(eeDb);  
  EEsmdCal::init();
  
  gMessMgr->Message("","I") <<GetName()<<"::Init() filter trigID="<<trigID<<"  maxCtbSum="<<maxCtbSum<<endm;  
  return StMaker::Init();
}

//________________________________________________
//________________________________________________
Int_t 
MuEzSmdCalMaker::InitRun(int runNo){
  if(runNo==0) {
    gMessMgr->Message("","W")<<GetName()<<"::InitRun("<<runNo<<") ??? changed to 555, it s OK for M-C - perhaps, JB"<<endm;
    runNo=555;
  }

  initRun(runNo);
 
  return kStOK;
}

//________________________________________________
//________________________________________________
Int_t 
MuEzSmdCalMaker::Finish(){
  finish(1);// do not draw
  gMessMgr->Message("","I") <<GetName()<<"::Finish()\n    inputEve="<<nInpEve<<" trigFilterEve="<<nTrigEve<<" nCorrEve="<<nCorrEve<<" nAcceptEve="<<nAcceptEve<<endm;
  return kStOK;
}

//________________________________________________
//________________________________________________
void 
MuEzSmdCalMaker::Clear(const Option_t*){
  eHead=0;
  eETow=0;
  eESmd=0;
  eTrig=0;
  //  delete trgAkio; //JAN: perhaps is a memory leak? But crashes
}



//________________________________________________
//________________________________________________
Int_t 
MuEzSmdCalMaker::Make(){
  clear();
  nInpEve++;
  gMessMgr->Message("","D") <<GetName()<<"::Make() is called "<<endm;
  
  //..........  acquire EztHeader
  eHead= mMuDstMaker->muDst()->eztHeader();
  if(eHead==0) {
    gMessMgr->Message("","E") <<GetName()<<"::Make() no  EztEventHeader, skip event  "<<endm;    return kStOK;
  }

  if(nInpEve==1) eHead->print();

  if(trigID ) {// filter by triggerID on demand
    if (! mMuDstMaker->muDst()->event()->triggerIdCollection().nominal().isTrigger(trigID)) return kStOK;
  }
  nTrigEve++;
  //.... get data .....
  eETow=mMuDstMaker->muDst()->eztETow();
  eESmd=mMuDstMaker->muDst()->eztESmd();
  eTrig=mMuDstMaker->muDst()->eztTrig(); 
  // printf("pp %p %p %p\n",eETow, eESmd, eTrig);
  if(!eETow || !eESmd || !eTrig)  return kStOK;

  //  trgAkio=new StTriggerDataMother(eTrig);
  // trgAkio->dump();
  // eETow->print();
  // eESmd->print();
  //  eHead->print();
  
  // .... process adata ......
  void *blob=eTrig->trgd->GetArray();

  StTriggerData2005 trgAkio5( (const TrgDataType2005 *)blob);
  if(eETow->doTowerHeadCorruptionTest(trgAkio5.token())||
     eESmd->doMapmtHeadCorruptionTest(trgAkio5.token())
     ) {
    nCorrEve++;
    return kStOK;
  }
  
  int ctbSum=trgAkio5.ctbSum();
  if(maxCtbSum>0 && (ctbSum>maxCtbSum || ctbSum<maxCtbSum/2.))  return kStOK;
  nAcceptEve++;
  
  //  printf("\n ctbSum=%d \n",ctbSum);
  //x  hA[7]->Fill(ctbSum);


  unpackMuEzt(eETow);
  unpackMuEzt(eESmd);
    
  findSectorMip();// do real analysis
 
  return kStOK;
} 

//________________________________________________
//________________________________________________
void
MuEzSmdCalMaker::unpackMuEzt(EztEmcRawData  *eRaw){
  int  n1=0,n2=0,n3=0;

  if(eRaw==0) return ; // no data block
  int icr;
  for(icr=0;icr<eRaw->getNBlocks();icr++) {
    if(eRaw->isCrateVoid(icr)) continue;
    assert(!eRaw ->getCorruption(icr)); // zero-tolerance

    int crateID=eRaw->getCrateID(icr);
    int chan;
    const UShort_t* data=eRaw->data(icr);
    for(chan=0;chan<eRaw->sizeData(icr);chan++) {
      const  EEmcDbItem  *x=eeDb->getByCrate(crateID,chan);
      if(x==0) continue; 
      if(x->sec!=sectID && crateID>6 ) break;// assumes crates do not cross sectors for non-towers, faster code
      if(x->fail ) continue; // drop broken channels
      // accept this hit
      float adc=-100, rawAdc=-101;
 
      rawAdc=data[chan];  
      adc=rawAdc-x->ped; 

      if(x->isSMD()) {
	//........................  SMD  U or V .......
	if(rawAdc>x->thr)  n3++;
	int iuv=x->plane-'U';
	int istr=x->strip-1;
	assert(iuv>=0 && iuv<MaxSmdPlains);
	assert(istr>=0 && istr<MaxSmdStrips);
	smdAdc[iuv][istr]=adc;      
	if(x->gain<=0)continue; // drop channels w/o gains
	smdEne[iuv][istr]=adc/x->gain; 
	// if(rawAdc>x->thr) 	printf("%s %f %f \n",x->name,smdEne[iuv][istr],x->gain);
      } else { 
	//............................... tower/pre/post crates
	int iT=-1;// store T,P,Q,R depending on 'iT'
	if(x->name[2]=='T'){
	  iT=0;
	} else{
	  iT=x->name[2]-'P'+1;
	}
	assert(iT>=0 && iT<mxTile);

	int iphi=(x->sec-1)*MaxSubSec+(x->sub-'A');
	int ieta=x->eta-1;
	//	x->print();
	assert(iphi>=0 && iphi<MaxPhiBins);
	assert(ieta>=0 && ieta<MaxEtaBins);   
	tileAdc[iT][ieta][iphi]=adc;// store towers
	tileThr[iT][ieta][iphi]=rawAdc>x->thr;
	killT[iT][ieta][iphi]=false; // it is alive
	if(rawAdc>x->thr) { 
	  if(iT==0)  
	    n1++;
	  else
	    n2++;
	}
	
	if(x->gain<=0) continue;// drop channels w/o gains
	tileEne[iT][ieta][iphi]=adc/x->gain;
      }
    } // end of loop over one data block
  }// end of loop over blocks
  // printf("%s-->nTow=%d nPQR=%d nSMD=%d\n",GetName(),n1,n2,n3);
}

//---------------------------------------------------
// $Log: MuEzSmdCalMaker.cxx,v $
// Revision 1.1  2005/03/11 15:44:25  balewski
// works with muEzt, cucu200
//
  