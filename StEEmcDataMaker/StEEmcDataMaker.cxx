// $Id: StEEmcDataMaker.cxx,v 1.14 2004/04/03 06:32:45 balewski Exp $

#include <Stiostream.h>
#include <math.h>
#include <assert.h>

#include "StEventTypes.h"
#include "StEvent.h"

#include "StEEmcDataMaker.h"

#include "StDAQMaker/StDAQReader.h"
#include "StDAQMaker/StEEMCReader.h"

#include "StEEmcDbMaker/StEEmcDbMaker.h"
#include "StEEmcDbMaker/StEEmcDbIndexItem1.h"

#include "StEEmcDbMaker/EEmcDbCrate.h"

#include "StEEmcUtil/EEfeeRaw/EEfeeDataBlock.h" 

//#include "EEmcHealth.h"

ClassImp(StEEmcDataMaker)

//_____________________________________________________________________
StEEmcDataMaker::StEEmcDataMaker(const char *name):StMaker(name){
  mDb=0;

}

//___________________________________________________________
StEEmcDataMaker::~StEEmcDataMaker(){

}

//__________________________________________________
//__________________________________________________
//__________________________________________________

Int_t StEEmcDataMaker::Init(){
  if(mDb==0){ 
    mDb= (StEEmcDbMaker*) GetMaker("eeDb");
  } 
  
  if(mDb==0){  
    printf("\n\nWARN %s::Init() did not found \"eeDb-maker\", all EEMC data will be ignored\n\n", GetName());
  } 

  hs[0]= new TH1F("health","raw data health; X: 0=nEve, 1=raw, 2=OKhead , 3=tower(No ghost/n256)",9,-1.5,7.5);

  hs[1]= new TH1F("n256","No. of n256/eve, all header OK",100, -1.5,98.5);
  hs[2]= new TH1F("nGhost","No. of tower nGhost/eve, all header OK, chan>119",100,-1.5,98.5);
  
  return StMaker::Init();
}


//___________________________________________________
//___________________________________________________
//___________________________________________________

Int_t StEEmcDataMaker::InitRun  (int runNumber){
  printf("\n%s::InitRun(%d) list  DB content \n",GetName(),runNumber);
 if(mDb==0){  
    printf("\n\nWARN %s::InitRun() did not found \"eeDb-maker\", all EEMC data will be ignored\n\n", GetName());
  } else if ( mDb->valid()==0 ) {
    printf("\n\nWARN %s::InitRun()  found \"eeDb-maker\", but without any DB data, all EEMC data will be ignored\n\n", GetName());
    mDb=0;
  } else {
    mDb->print(0);
  }
  return kStOK;
}


//____________________________________________________
//____________________________________________________
//____________________________________________________
Int_t StEEmcDataMaker::Make(){
  if(mDb==0){  
    printf("WARN %s::Make() did not found \"eeDb-maker\" or no DB data for EEMC, all EEMC data will be ignored\n", GetName());
    return kStOK;
  } 
  
  StEvent* mEvent = (StEvent*)GetInputDS("StEvent");
  assert(mEvent);// fix your chain or open the right event file
  // printf("\n%s  accesing StEvent ID=%d\n",GetName(),mEvent->id());

  hs[0]->Fill(0);
 //::::::::::::::::: copy raw data to StEvent ::::::::::::: 
  if(!copyRawData(mEvent)) return kStOK;
  hs[0]->Fill(1);

 //::::::::::::::: assure raw data are sane  :::::::::::::::: 

 if(headersAreSick(mEvent))   return kStOK;
 hs[0]->Fill(2);

 if(towerDataAreSick(mEvent))   return kStOK;
 hs[0]->Fill(3);


 //:::::::::::::: store tower/pre/post/smd hits in StEvent
 raw2pixels(mEvent);
   
 return kStOK;
  }

//____________________________________________________
//____________________________________________________
//____________________________________________________
int   StEEmcDataMaker::copyRawData(StEvent* mEvent) {
  
  St_DataSet *daq = GetDataSet("StDAQReader");                 assert(daq);
  StDAQReader *fromVictor = (StDAQReader*) (daq->GetObject()); assert(fromVictor);
  StEEMCReader *eeReader  = fromVictor->getEEMCReader();  
  if(!eeReader) return false ;

  int icr; 

  StEmcCollection* emcC =(StEmcCollection*)mEvent->emcCollection();

  if(emcC==0) { // create this collection if not existing
    emcC=new StEmcCollection();
    mEvent->setEmcCollection(emcC);
    printf("%s::Make() has added a non existing StEmcCollection()\n", GetName());
  }
  
  StEmcRawData *raw=new  StEmcRawData;  
  emcC->setEemcRawData(raw);
  printf("%s::copy %d EEMC raw data blocks eveID=%d\n",GetName(),mDb->getNCrate(),mEvent->id());

  for(icr=0;icr<mDb->getNCrate();icr++) {
    const EEmcDbCrate *crate=mDb-> getCrate(icr);
    // printf("copy EEMC raw: ");crate->print();
    raw->createBank(icr,crate->nHead,crate->nCh);
    //  printf("aa=%p bb=%p\n",eeReader->getEemcHeadBlock(crate->fiber,crate->type),eeReader->getEemcDataBlock(crate->fiber,crate->type));

    raw->setHeader(icr,eeReader->getEemcHeadBlock(crate->fiber,crate->type));
    raw->setData(icr,eeReader->getEemcDataBlock(crate->fiber,crate->type));
  }
  return true;
}

//____________________________________________________
//____________________________________________________
//____________________________________________________
int  StEEmcDataMaker::headersAreSick(StEvent* mEvent) {
  
  StEmcCollection* emcC =(StEmcCollection*)mEvent->emcCollection();

  if(emcC==0) {
    printf("%s::headersAreSick() no emc collection, skip\n",GetName());
    return true;
  }

  StEmcRawData* raw=emcC->eemcRawData();
  assert(raw);
  
  StL0Trigger* trg=mEvent->l0Trigger(); assert(trg);
  int token=trg->triggerToken();

  EEfeeDataBlock block; // use utility class as the work horse

  int sick=false;
  int icr;
  for(icr=0;icr<mDb->getNCrate();icr++) {
    const EEmcDbCrate *crate=mDb-> getCrate(icr);
    if(!crate->useIt) continue; // drop masked out crates
    //printf(" EEMC raw-->pix crID=%d type=%c \n",crate->crID,crate->type);
    const  UShort_t* head=raw->header(icr);
    assert(head);
    block.clear();
    block.setHead(raw->header(icr));

    int lenCount=crate->nCh+crate->nHead;
    int errFlag=0;
    switch(crate->type) {
    case 'T': lenCount+=32; break; // one more board exist in harware
    case 'S':  errFlag=0x28; break; // sth is set wrong in the boxes
    default:;
    }

    int trigCommand=4; // physics, 9=laser/LED, 8=??
    int valid=block.isHeadValid(token,crate->crIDswitch,lenCount,trigCommand,errFlag);
 
    printf("valid=%d ",valid); crate->print();

    if(valid) continue;
    sick=true;
  }

  printf("%s::headersAreSick()==%d\n",GetName(),sick);
  return sick;
}



//____________________________________________________
//____________________________________________________
//____________________________________________________
int  StEEmcDataMaker::towerDataAreSick(StEvent* mEvent) {


  StEmcCollection* emcC =(StEmcCollection*)mEvent->emcCollection();

  assert(emcC);
  
  StEmcRawData* raw=emcC->eemcRawData();
  assert(raw);

  int nGhost=0, n256=0;
  int icr;
  for(icr=0;icr<mDb->getNCrate();icr++) {
    const EEmcDbCrate *crate=mDb-> getCrate(icr);
    if(!crate->useIt) continue; // drop masked out crates
    if(crate->type!='T') continue;
    const  UShort_t* data=raw->data(icr);
    assert(data);
    int i;
    for(i=0;i<raw->sizeData(icr);i++) {
      if((data[i] &0xff)==0 ) n256++;
      if(i>=121 && data[i]>40) nGhost++;
    }
  }

  hs[1]->Fill(n256);
  hs[2]->Fill(nGhost);

  printf("%s::towerDataAreSick() , n256=%d, nGhost=%d\n",GetName(),n256, nGhost);
  if(nGhost>0) return true;
  if(n256>20) return true;
  return false;
}



//____________________________________________________
//____________________________________________________
//____________________________________________________
void  StEEmcDataMaker::raw2pixels(StEvent* mEvent) {

  StEmcCollection* emcC =(StEmcCollection*)mEvent->emcCollection();

  assert(emcC);
  
  StEmcRawData* raw=emcC->eemcRawData();
  assert(raw);
  printf("%s::raw2pixels()\n",GetName());

  int mxSector = 12;

  // initialize tower/prePost /U/V in StEvent
  
  StEmcDetector* emcDet[kEndcapSmdVStripId+1]; // what a crap,JB
  StDetectorId emcId[kEndcapSmdVStripId+1];    // more crap
  memset(emcDet,0,sizeof(emcDet));

  int det;
  for(det = kEndcapEmcTowerId; det<= kEndcapSmdVStripId; det++){
    emcId[det] = StDetectorId(det);
    emcDet[det] = new StEmcDetector(emcId[det],mxSector);
    emcC->setDetector(emcDet[det]);
  }

  // store data from raw blocks to StEvent
  int nDrop=0;
  int nMap=0;
  int nTow=0;
  int icr;
  for(icr=0;icr<mDb->getNCrate();icr++) {
    const EEmcDbCrate *crate=mDb-> getCrate(icr);
    if(!crate->useIt) continue; // drop masked out crates
    
    const  UShort_t* data=raw->data(icr);
    assert(data);

    for(int chan=0;chan<raw->sizeData(icr);chan++) {
      const  StEEmcDbIndexItem1  *x=mDb->get(crate->crID,chan);
      if(x==0) {
	//printf("No EEMC mapping for crate=%3d chan=%3d\n",crate,chan);
	nDrop++;
	continue;
      }

      int isec=x->sec-1;   //range 0-11 ???
      int rawAdc=data[chan];
      float energy=123.456; // dumm value
      int det = 0;
      int ieta=0, isub=0;
      char type=x->name[2];

      switch(type) { // tw/pre/post/SMD
      case 'T': // towers
	isub=x->sub-'A'; //range 0-4
	ieta=x->eta-1;   //range 0-11
	det = kEndcapEmcTowerId;
	nTow++;
	break;
      case 'P': // pres1
      case 'Q': // pres2
      case 'R': // post
	isub=x->sub-'A'; //range 0-4
	ieta=x->eta-1;   //range 0-11
	isub+=5* (type-'P'); // this is even more sily
	det = kEndcapEmcPreShowerId;
	break;
      case 'U': //SMD
      case 'V': 
	isub=0; // not used for SMD
	ieta=x->strip-1;   //range 0-287
	det = kEndcapSmdUStripId;
	if(type=='V') det = kEndcapSmdVStripId;
	break;
      default:
	continue;
      }

      if(det==0) continue; // tmp
      //printf("EEMC crate=%3d chan=%3d  ADC: raw=%4d pedSub=%+.1f energy=%+10g  -->   %2.2dT%c%2.2d\n",crate,chan,rawAdc,adc,energy,x->sec,x->sub,x->eta);
      
      StEmcRawHit* h = new StEmcRawHit(emcId[det],isec,ieta,isub,rawAdc,energy);
      emcDet[det]->addHit(h);
      nMap++;

    } // end of loop over channels
  }// end of loop over crates 
  
    
  printf("%s event finished nDrop=%d nMap=%d nTow=%d\n",GetName(), nDrop,nMap,nTow);


  for(det = kEndcapEmcTowerId; det<= kEndcapSmdVStripId; det++){
   StEmcDetector* emcDetX= emcC->detector( StDetectorId(det));
   assert(emcDetX);
   printf("det=%d ID=%d nHits=%d\n",det,StDetectorId(det),emcDetX->numberOfHits());
  }  
}



 

// $Log: StEEmcDataMaker.cxx,v $
// Revision 1.14  2004/04/03 06:32:45  balewski
// firts attempt to store EEMC hits in StEvent & muDst,
// Implemented useIt for fibers in Db
// problems: - tower sector 12 is missing
// - no pres & msd in smd
//
// Revision 1.13  2004/04/02 06:38:43  balewski
// abort on any error in any header
//
// Revision 1.12  2004/03/30 04:44:52  balewski
// raw EEMC in StEvent done
//
// Revision 1.11  2004/03/28 04:08:12  balewski
// store raw EEMC data
//
// Revision 1.9  2004/03/20 20:25:40  balewski
// fix for empty ETOW/ESMD
//
// Revision 1.8  2004/03/19 21:29:28  balewski
// new EEMC daq reader
//
// Revision 1.7  2003/09/02 17:57:54  perev
// gcc 3.2 updates + WarnOff
//
// Revision 1.6  2003/08/18 21:31:42  balewski
// new adc--> energy formula, still no absolute gain
//
// Revision 1.5  2003/07/18 18:31:45  perev
// test for nonexistance of XXXReader added
//
// Revision 1.4  2003/05/01 02:19:19  balewski
// fixed empry event problem
//
// Revision 1.3  2003/04/29 16:18:39  balewski
// some ideas added
//
// Revision 1.2  2003/04/27 23:08:02  balewski
// clean up of daq-reader
//
// Revision 1.1  2003/04/25 14:15:59  jeromel
// Reshaped Jan's code
//
// Revision 1.4  2003/04/16 20:33:44  balewski
// small fixes in eemc daq reader
//
// Revision 1.3  2003/03/26 22:02:20  balewski
// fix auto db-find
//
// Revision 1.2  2003/03/26 21:27:44  balewski
// fix
//
// Revision 1.1  2003/03/25 18:30:14  balewski
// towards EEMC daq reader
//


  //old .......................



#if 0 // test of new access method  

  for(icr=0;icr<mDb->getNCrate();icr++) {
    const EEmcDbCrate *crate=mDb-> getCrate(icr);
    printf("geting data for fiber: ");crate->print();

    printf("---- HEAD ----\n");
    for(ch=0;ch<crate->nHead;ch++) {
      int val=-1;
      val=steemcreader->getEemcHead(crate->fiber,ch,crate->type); 
      printf("cr=%d ch=%d val=0x%04x\n",crate->crID,ch,val);
    }

    printf("---- DATA ----\n");
    for(ch=0;ch<crate->nch;ch++) {
      int val=-1;
      val=steemcreader->getEemcData(crate->fiber,ch,type); 
      printf("cr=%d ch=%d val=0x%04x\n",crate->crID,ch,val);
    }
    //   break;
  }
#endif


