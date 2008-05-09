//*-- Author : Jan Balewski
// 
// $Id: StEEmcFilterMaker.cxx,v 1.2 2008/05/09 22:14:35 balewski Exp $

#include <StThreeVector.hh>
#include <StPrimaryVertex.h>
#include "StEEmcPool/StEEmcA2EMaker/StEEmcA2EMaker.h"

#include "StEvent.h"
#include "StMessMgr.h"

#include "StEEmcFilterMaker.h"

ClassImp(StEEmcFilterMaker)

//_____________________________________________________________
//_____________________________________________________________
StEEmcFilterMaker::StEEmcFilterMaker(const char *name):StMaker(name){
  par_Et_thres=-1;
  par_Z0_vert=-9999;
  par_delZ_vert=-2;
  //
}


//_____________________________________________________________
//_____________________________________________________________
/// This is TLA destructor
StEEmcFilterMaker::~StEEmcFilterMaker(){
  
  //
}


//_____________________________________________________________
//_____________________________________________________________

Int_t StEEmcFilterMaker::Init(){
  SetAttr(".Privilege",1); //this maker so BFC 'listens' to it and skip events
  mGeomE = new EEmcGeomSimple(); //for trackXEndCap
  
  mH0=new TH1F("mH0","Event counter",5,0.5,5.5);
  nInpEve=nRecVert=nZverOK=nAccEve=0;
  LOG_INFO << Form("Init ee-filter cuts: ET>%.2f  Zvert=%.2f +/-%.2f (cm)", par_Et_thres, par_Z0_vert, par_delZ_vert)<<endm;

  assert(par_delZ_vert>0); // probably forgot to set all params for thsi maker
  return StMaker::Init();
}

//_____________________________________________________________
//_____________________________________________________________
Int_t StEEmcFilterMaker::FinishRun(int runumber){
  LOG_INFO << Form("Finish cuts: ET>%.2f  Zvert=%.2f +/-%.2f (cm)", par_Et_thres, par_Z0_vert, par_delZ_vert)<<endm;
  LOG_INFO << Form("Finish run=%d nInp=%d,nRecVer=%d, nZverOK=%d nAcc=%d",runumber,nInpEve,nRecVert,nZverOK,nAccEve) << endm;
  return kStOK;
}; 


//_____________________________________________________________
//_____________________________________________________________
Int_t StEEmcFilterMaker::Make(){
  LOG_INFO << Form("in::Make inp=%d,nVer=%d, nZverOK=%d nAcc=%d",nInpEve,nRecVert,nZverOK,nAccEve) << endm;
  nInpEve++;
  mH0->Fill(1);
  StEvent* mEvent = (StEvent*)GetInputDS("StEvent");
  assert(mEvent);// fix your chain or open the right event file
   
  int nV=mEvent->numberOfPrimaryVertices();
  if(nV<=0) return  kStSKIP;
 
  int iv=0; // pick always the first vertex, may be wrong w/ pileup, fix it,JB
  nRecVert++;
  mH0->Fill(2); 
  StPrimaryVertex*V=mEvent->primaryVertex(iv);
  assert(V);
  StThreeVectorF r=V->position();
  Float_t vertexPosZ = r.z();
  if(fabs(vertexPosZ - par_Z0_vert)> par_delZ_vert)  return  kStSKIP;

  nZverOK++;
  mH0->Fill(3); 
  StEEmcA2EMaker *a2eMK=(StEEmcA2EMaker *)GetMaker("EE_A2E");
  assert(a2eMK); //for endcap tower stuff
  StEEmcTower highTow = a2eMK->hightower(0);//find high tower

  Float_t highTEt = highTow.et();
  LOG_INFO << Form("eveID=%d  nPrimVert=%d  zVert=%.2f highTEt=%.2f\n", mEvent->id(),nV, vertexPosZ,highTEt);

  Float_t triggerPatchEt;
  Int_t triggerConditionReturn = triggerCondition(vertexPosZ,&highTow,triggerPatchEt);


  if(triggerConditionReturn != 1)    {
    printf("janKill-2 Me### eve=%d nv=%d zVer=%.2f highT_et=%.3f TP_ET=%.3f\n",nInpEve,nV,vertexPosZ, highTEt,triggerPatchEt);
    return kStSKIP;
  }

  nAccEve++;
  mH0->Fill(4); 
  //process event only if et above threshold
  printf("janKeep-2 Me### eve=%d  nv=%d zVer=%.2f highT_et=%.3f TP_ET=%.3f\n",nInpEve,nV,vertexPosZ, highTEt,triggerPatchEt);
  

  return kStOK;
}
//_____________________________________________________________
//_____________________________________________________________
Int_t StEEmcFilterMaker::triggerCondition( Float_t vertexPosZ, StEEmcTower *highTow,  Float_t &patchEt)
{

  //Float_t eemcPatchEt = highTow.et();//use my et function
  Float_t eemcPatchEt = transverseNRG(vertexPosZ,highTow);
  Int_t Nneigh = highTow->numberOfNeighbors();
  for(int Nn=0; Nn < Nneigh; Nn++)
    {
      StEEmcTower neighHighTow = highTow->neighbor(Nn);
      Float_t neighEemcPatchEt = transverseNRG(vertexPosZ,&neighHighTow);
      eemcPatchEt += neighEemcPatchEt;
    }

  patchEt = eemcPatchEt;

  if(eemcPatchEt <  par_Et_thres) return 0;

  return 1;
  
}
//_____________________________________________________________
//_____________________________________________________________
Float_t StEEmcFilterMaker::transverseNRG(Float_t vertexPosZ, StEEmcTower *tower)
{
  //find the transverse energy for a tower when vertex != 0
  //now we do everything with respect to the Smd plane

  TVector3 towVector = mGeomE->getTowerCenter(tower->sector(),tower->subsector(),tower->etabin());
  //get x,y coords at Smd plane
  Float_t xNew = towVector.x() - (288.2 - 279.54)*(towVector.x()/288.2);
  Float_t yNew = towVector.y() - (288.2 - 279.54)*(towVector.y()/288.2);

  TVector3 *eventTowVector = new TVector3(0,0,0);//make vector from vertex to smd plane
  eventTowVector->SetXYZ(xNew,yNew,(279.54 - vertexPosZ));

  //Float_t newEta = eventTowVector->Eta();
  //Float_t newPhi = eventTowVector->Phi();
  Float_t newTheta = eventTowVector->Theta();

  Float_t factor = sin(newTheta);//multiply energy by this to get event ET
  if(factor > 1.0 || factor < 0.0)
    {
      return -1.0;
    }

  Float_t energy = tower->energy();

  Float_t transverseEnergy = energy*factor;

  return transverseEnergy;
}

// $Log: StEEmcFilterMaker.cxx,v $
// Revision 1.2  2008/05/09 22:14:35  balewski
// new there are 2 filters
//
// Revision 1.1  2008/04/21 15:47:09  balewski
// star
//
