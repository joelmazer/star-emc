// $Id: EEmcMCData.cxx,v 1.11 2004/04/08 21:34:55 perev Exp $

#include "StEventTypes.h"

//#include "emc_def.h"          /* FIXME - move it to pams */


#include "TError.h"
//
#include "tables/St_g2t_emc_hit_Table.h"
#include "tables/St_g2t_event_Table.h"
//
#include "StBFChain.h"
//
#include "EEmcException.h"
#include "EEmcMCData.h"

#include "StEEmcUtil/EEmcGeom/EEmcGeomDefs.h"
#include "StEEmcUtil/EEevent/EEeventDst.h"
#include "StEEmcUtil/EEevent/EEsectorDst.h"

#include <iostream>
using namespace std;


ClassImp(EEmcMCData)


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
EEmcMCData::EEmcMCData()
{ 
  mEventID   = -1;  
  mLastHit   = 0;
  mEthr      = kEEmcDefaultEnergyThreshold; 
  mSize      = kEEmcDefaultMCHitSize;
  mHit       = new struct EEmcMCHit[mSize];   // for now a static array
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
EEmcMCData::EEmcMCData( EEmcMCData &e )
{ 
  mEventID   = e.getEventID();
  mLastHit   = e.getLastHit();
  mEthr      = e.getEnergyThreshold();
  mSize      = e.getSize();
  mHit       = new struct EEmcMCHit[mSize];
  mSize      = e.getHitArray(mHit,mSize);
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
EEmcMCData::~EEmcMCData() {
  mSize  = 0;
  if(mHit) delete [] mHit;
}


//-------------------------------------------------------------------------
// decode Endcap Emc data
//-------------------------------------------------------------------------
Int_t             
EEmcMCData::readEventFromChain(StMaker *myMk)
{
  St_g2t_event   *g2t_event  = NULL; // g2t event
  St_g2t_emc_hit *emc_hit    = NULL; // tower hits
  St_g2t_emc_hit *smd_hit    = NULL; // smd hits
  g2t_event_st   *event_head = NULL; // g2t event header
  
  int err=0;

  mLastHit  = 0;
  memset(mHit,0x0,sizeof(struct EEmcMCHit)*mSize); 

  
  // get info from event header
  if( (g2t_event=(St_g2t_event *) myMk->GetDataSet("g2t_event")) == NULL ) 
    throw EEmcException1(kEEmcMCMissingEventHeader,"missing MC event header");

  if( (event_head= g2t_event->GetTable()) == NULL )
    throw EEmcException1(kEEmcMCMissingEventHeader,"missing MC event header table");

  mEventID=event_head->n_event;

  // get tower data
  if( (emc_hit = (St_g2t_emc_hit *) myMk->GetDataSet("g2t_eem_hit")) !=NULL ) {
    Int_t     nhits      = emc_hit->GetNRows();
    if(nhits<=0) { 
      Warning("readEventFromChain","no tower hits (%d)",nhits);
      goto skipTower;
    }

    g2t_emc_hit_st *hit  = emc_hit->GetTable();
    for(Int_t ihit=0; ihit<nhits; ihit++,hit++) {

      Int_t   ivid  = hit->volume_id;
      Short_t sec   = 0;
      Short_t ssec  = 0;
      // Short_t half  = ivid/kEEmcTowerHalfId;
      ivid %= kEEmcTowerHalfId;
      Short_t phi   = ivid/kEEmcTowerPhiId;  ivid %= kEEmcTowerPhiId;
      Short_t eta   = ivid/kEEmcTowerEtaId;  ivid %= kEEmcTowerEtaId;
      Short_t depth = ivid/kEEmcTowerDepId;  ivid %= kEEmcTowerDepId;  
      //printf("Tw hit->volume_id=%d hit->de=%g phi=%d\n",hit->volume_id,hit->de,phi);

      if(!ivid==0){err=1; goto crash;};
       
      ssec = (phi-1)%5 + 1; // sub-sector
      
      sec  = (phi-1)/5 + 1;

      if(!( 0<sec   && sec<=kEEmcNumSectors   )){err=2; goto crash;};
      if(!( 0<ssec  && ssec<=kEEmcNumSubSectors)){err=3; goto crash;};
      if(!( 0<eta   && eta<=kEEmcNumEtas)){err=4; goto crash;};
      if(!( 0<depth && depth<=kEEmcNumDepths)){err=5; goto crash;};
      
      switch(depth) {
      case kPreShower1Depth:
	mHit[mLastHit].detector = kEEmcMCPreShower1Id; 
	break; 
      case kPreShower2Depth: 
	mHit[mLastHit].detector = kEEmcMCPreShower2Id; 
	break;
      case kTower1Depth:    
      case kTower2Depth:    
	mHit[mLastHit].detector = kEEmcMCTowerId;     
	break;
      case kPostShowerDepth:
	mHit[mLastHit].detector = kEEmcMCPostShowerId; 
	break;
      default:
	Warning("readEventFromChain","unknown depth %d",depth);
	throw EEmcException1(kEEmcMCErr1,"invalid depth for tower tails");
	break;
      }
      
      mHit[mLastHit].sector      = sec;
      mHit[mLastHit].tower.ssec  = ssec;
      mHit[mLastHit].tower.eta   = eta;
      mHit[mLastHit].de          = hit->de;
      mLastHit++;
      // printf("depth=%d nH=%d\n",depth,mLastHit);

      if(mLastHit>=mSize && !expandMemory() ) 
	throw EEmcException1(kEEmcMCErr1,"failed expandMemory() for tower tails");

    } // end of tower hits
  } 

 skipTower:

  // get smd data
  if( (smd_hit = (St_g2t_emc_hit *) myMk->GetDataSet("g2t_esm_hit")) != NULL ) {
    Int_t nhits         = smd_hit->GetNRows();

    if(nhits<=0) { 
      Warning("readEventFromChain","no smd hits (%d)",nhits);
      goto done;
    }
    g2t_emc_hit_st *hit = smd_hit->GetTable();

    for(Int_t ihit=0; ihit<nhits ; ihit++,hit++) { 
      //printf("Smd hit->volume_id=%d hit->de=%g\n",hit->volume_id,hit->de);

      Int_t   ivid  = hit->volume_id;
      Short_t det   = 0;
      Short_t sec   = 0;
      Short_t half  = ivid/kEEmcSmdHalfId; ivid %= kEEmcSmdHalfId;
      Short_t phi   = ivid/kEEmcSmdPhiId;  ivid %= kEEmcSmdPhiId;
      Short_t plane = ivid/kEEmcSmdPlaneId;ivid %= kEEmcSmdPlaneId;
      Short_t strip = ivid/kEEmcSmdStripId;ivid %= kEEmcSmdStripId;  
      
      if(!ivid==0){err=10; goto crash;};

      switch(phi) { /* FIXME ONE DAY */
      case 1:
      case 4:
      case 7:
      case 10:
	switch(plane) {
	case 1:  det = kEEmcMCSmdVStripId; break; 
	case 3:  det = kEEmcMCSmdUStripId; break; 
	default: det = kUnknownId;         break;
	}
	break;
      case 2:
      case 5:
      case 8:
      case 11:
	switch(plane) {
	case 2:  det = kEEmcMCSmdVStripId; break;
	case 1:  det = kEEmcMCSmdUStripId; break;
	default: det = kUnknownId;         break;
	}
	break;
      case 3:
      case 6:
      case 9:
      case 12:
	switch(plane) {
	case 3:  det = kEEmcMCSmdVStripId; break;
	case 2:  det = kEEmcMCSmdUStripId; break;
	default: det = kUnknownId;         break;
	}
	break;
      default:
	det = kUnknownId;         
	break;
      }
      if(det!=kEEmcMCSmdVStripId && det!=kEEmcMCSmdUStripId ) { 
	Warning("readEventFromChain","unknown smd layer %d %d-%d-%d-%d",det,half,phi,plane,strip);
	throw EEmcException1(kEEmcMCErr1,"invalid depth for SMD");
      }

      sec = phi;

      if(! ( 0<sec   && sec  <=kEEmcNumSectors )){err=12; goto crash;};
      if(!( 0<strip && strip<=kEEmcNumStrips  )){err=13; goto crash;};

      // printf("Smd hit->volume_id=%d hit->de=%g sec=%d  U/V=%d, strip=%d\n",hit->volume_id,hit->de,sec, det,strip );
      // fill in
      mHit[mLastHit].detector = det;
      mHit[mLastHit].sector   = sec;
      mHit[mLastHit].strip    = strip;
      mHit[mLastHit].de       = hit->de;

      mLastHit++;
      if(mLastHit>=mSize && !expandMemory() ) 
	throw EEmcException1(kEEmcMCErr1,"failed expandMemory() for SMD strips");
    }
  }

 done:
  return mLastHit;

 crash: 
  printf("\n\n==============================\nEEmcMCData::readEventFromChain() Fatal error while decoding .fzd hits for EEMC err=%d,\n all EEMC data erased from StEvent\n=====================================\n",err);
  mLastHit  = 0;
  memset(mHit,0x0,sizeof(struct EEmcMCHit)*mSize); 
  return mLastHit;

}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
Int_t      
EEmcMCData::getHitArray(EEmcMCHit *h, Int_t size)
{
  if(size<=0) {
    cerr << "invalid size: " << size << endl;
    return 0;
  }
  if(size<mSize) 
    cerr << "truncating to " << size << " hits  (out of " << mSize << ")" << endl;
  int n = size; if (n > mSize) n = mSize;
  memcpy(h,mHit,size*sizeof(EEmcMCHit));
  return size;
}

Int_t      
EEmcMCData::setHitArray(EEmcMCHit *h, Int_t size)
{
  if(size<=0) {
    cerr << "invalid size: " << size << endl;
    return 0;
  }
  if(size>mSize) 
    cerr << "truncating to " << mSize << " hits  (out of " << size << ")" << endl;
  int n = size; if (n > mSize) n = mSize;
  memcpy(h,mHit,n*sizeof(EEmcMCHit));
  return mSize;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void
EEmcMCData::print()
{
  TString detName;
  Int_t   detId;
  EEmcMCHit *h  = mHit;  
  printf("EndcapEmc MC event #%d\n",mEventID);
  for(Int_t i=0; i<mLastHit; i++,h++) {
    detId=h->detector;
    switch(detId) {
    case kEEmcMCTowerId:      detName = "EndcapEmcTower      "; break;
    case kEEmcMCPreShower1Id: detName = "EndcapEmcPreShower1 "; break;
    case kEEmcMCPreShower2Id: detName = "EndcapEmcPreShower2 "; break;
    case kEEmcMCPostShowerId: detName = "EndcapEmcPostShower "; break;
    case kEEmcMCSmdUStripId:  detName = "EndcapEmcSmdUStrip  "; break;
    case kEEmcMCSmdVStripId:  detName = "EndcapEmcSmdVStrip  "; break;
    default:                  detName = "EndcapEmcUnknown    "; break;
    }
    

    switch(detId) {
    case kEEmcMCTowerId:
    case kEEmcMCPreShower1Id:
    case kEEmcMCPreShower2Id:
    case kEEmcMCPostShowerId:
      printf("%s sec=%2d sub=%c eta=%d de=%g\n",detName.Data(), h->sector,h->tower.ssec-1+'A',h->tower.eta,h->de); 
      break;
    case kEEmcMCSmdUStripId:
    case kEEmcMCSmdVStripId:
      printf("%s sec=%2d strip=%d de=%g\n",detName.Data(), h->sector,h->strip,h->de); 
      break;
    default:
      cout << "**** WARNING **** detectorId=" << detId << " is unknown";
      assert(1==2);
      break;
    }

  }
}

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
Int_t             
EEmcMCData::expandMemory()
{
  Int_t      newSize =  mSize + kEEmcDefaultMCHitSize;
  EEmcMCHit* newHit  =  new EEmcMCHit[newSize];
  Assert(newHit);

  if(mHit) {
    memcpy(newHit,mHit,mSize*sizeof(EEmcMCHit));
    delete [] mHit;
    mHit = 0;
  }
  Info("expandMemory"," MCHit memory expanded from 0x%04x to 0x%04x",mSize,newSize);
  mSize = newSize;
  mHit  = newHit;
  return kTRUE;
}


//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
Int_t             
EEmcMCData::read   (void *hit, int size)
{ 

  const int HitSize = sizeof(struct EEmcMCHit);
  const int MaxSize = mSize * HitSize;
  if(size>MaxSize) size=MaxSize;
  memcpy(mHit,hit,size);
  mLastHit =size/HitSize;
  return size;
}
  
//-------------------------------------------------------------------------
// writes binary hits
//-------------------------------------------------------------------------
Int_t             
EEmcMCData::write  (void *hit, int size)
{
  Int_t nbytes = mLastHit * sizeof(struct EEmcMCHit);
  if(size>nbytes) size=nbytes;
  memcpy(hit,mHit,size);
  return size;
}


//-------------------------------------------------------
//-------------------------------------------------------
Int_t EEmcMCData::write(EEeventDst *EEeve) {
  //printf("Clear & Export EEeventDst to TTree\n");
  EEeve->clear();

  EEeve->setType(EEeventDst::kRawMC);
  EEeve->setID(mEventID);

  EEmcMCHit *h  = mHit;
  for(Int_t i=0; i<mLastHit; i++,h++) {
    // Info("EEmcMCData().write(TTree)","depth=%2d sector=%2d sub=%2d eta=%2d energy=%f", h->depth,h->sec,h->ssec,h->eta,h->de);

    int secID=h->sector;
    EEsectorDst *EEsec= (EEsectorDst *)EEeve->getSec(secID,1);
    assert(EEsec);

    // temporary projection,  JB

    int subSec='A'+h->tower.ssec-1;
    switch( h->detector) {
    case kEEmcMCTowerId:
      EEsec->addTwHit(subSec,h->tower.eta,h->de);    break;
    case kEEmcMCPreShower1Id:
      EEsec->addPre1Hit(subSec,h->tower.eta,h->de);  break;
    case kEEmcMCPreShower2Id:
      EEsec->addPre2Hit(subSec,h->tower.eta,h->de);  break;
    case kEEmcMCPostShowerId:
      EEsec->addPostHit(subSec,h->tower.eta,h->de); break;
    case kEEmcMCSmdUStripId:
      EEsec->addSmdUHit(h->strip,h->de);  break;
    case kEEmcMCSmdVStripId:
      EEsec->addSmdVHit(h->strip,h->de);  break;
    default:
      throw EEmcException1(kEEmcMCInvalidDepth,"invalid MC depth");
      break;
    }
      
  } // end of loop over all hits

  return 0;
}

// $Log: EEmcMCData.cxx,v $
// Revision 1.11  2004/04/08 21:34:55  perev
// Leak off
//
// Revision 1.10  2003/11/12 19:59:13  balewski
// *** empty log message ***
//
// Revision 1.9  2003/09/17 22:05:36  zolnie
// delete mumbo-jumbo
//
// Revision 1.8  2003/09/11 19:41:08  zolnie
// updates for gcc3.2
//
// Revision 1.7  2003/09/02 17:57:56  perev
// gcc 3.2 updates + WarnOff
//
// Revision 1.6  2003/03/26 21:16:59  balewski
// swap U & V on Wei-Ming request
//
// Revision 1.5  2003/03/22 19:37:33  balewski
// *** empty log message ***
//
// Revision 1.4  2003/02/21 15:31:38  balewski
// do not kill the chain (it is against my will, JB)
//
// Revision 1.3  2003/02/20 21:27:06  zolnie
// added simple geometry class
//
// Revision 1.2  2003/02/20 20:13:20  balewski
// fixxy
// xy
//
// Revision 1.1  2003/02/20 05:14:07  balewski
// reorganization
//
// Revision 1.2  2003/02/14 00:04:39  balewski
// remove few printouts
//
// Revision 1.1  2003/01/28 23:16:07  balewski
// start
//
// Revision 1.21  2002/12/17 19:43:31  balewski
// remove some dependecies, choose better name for maker
//
// Revision 1.20  2002/12/10 19:02:33  balewski
// fix sector ID according to SN #299
// tower: volID contains absolute  subsectorID in the range 1-60, subSec #1 is the first subSec of sector 1,  and advances clockwise
// SMD :  volID contains absolute sector ID from 1-12,
// it does not matter what geometry has been used (i.e. 4 sectors, lower half or both halfs)
// The half ID info in volID is redundant now
//
// Revision 1.19  2002/12/05 14:21:58  balewski
// cleanup, time stamp corrected
//
// Revision 1.18  2002/12/04 15:03:36  balewski
// fix: read also just tower or just SMD data
//
// Revision 1.17  2002/11/13 21:53:51  zolnie
// restored "private" DetectorID in EEmcMCData
//
// Revision 1.16  2002/11/11 21:22:48  balewski
// EEMC added to StEvent
//
// Revision 1.15  2002/10/03 15:52:25  zolnie
// updates reflecting changes in *.g files
//
// Revision 1.14  2002/10/03 00:30:45  balewski
// tof taken away
//
// Revision 1.13  2002/10/01 14:41:54  balewski
// SMD added
//
// Revision 1.12  2002/10/01 06:03:15  balewski
// added smd & pre2 to TTree, tof removed
//
// Revision 1.11  2002/09/30 21:58:27  zolnie
// do we understand Oleg? (the depth problem)
//
// Revision 1.10  2002/09/30 21:04:02  zolnie
// SMD updates
//
// Revision 1.9  2002/09/30 20:15:55  zolnie
// Oleg's geometry updates
//
// Revision 1.8  2002/09/27 19:05:13  zolnie
// EEmcMCData updates
//
// Revision 1.5  2002/09/25 01:36:13  balewski
// fixed TOF in geant
//
// Revision 1.4  2002/09/24 22:47:35  zolnie
// major rewrite: SMD incorporated, use constants rather hard coded numbers
// 	introducing exceptions (rather assert)
//
// Revision 1.3  2002/09/20 21:58:13  balewski
// sum of MC hits over activ detectors
// produce total tower energy with weight 1 1 1 1
//
// Revision 1.2  2002/09/20 15:49:05  balewski
// add event ID
//
// Revision 1.1.1.1  2002/09/19 18:58:54  zolnie
// Imported sources
//
// Revision 1.1.1.1  2002/08/29 19:32:01  zolnie
// imported sources
// Revision 1.2  2002/08/28 01:42:59  zolnie
// version alpha ready
// Revision 1.1  2002/08/27 16:39:10  zolnie
// Initial revision
//
//
