/**********************************************************************
* StSmdPedMaker
* Author: Alexandre A. P. Suaide 
*
* This maker does calibration on the EMC detector
***********************************************************************/

#include "StSmdPedMaker.h"
#include "StEvent/StEvent.h"
#include "StEvent/StEventTypes.h"
#include "TStopwatch.h"
#include <fstream.h>
#include "TFile.h"
#include "StEmcUtil/StEmcGeom.h"
/*#include "StEmcUtil/StEmcFilter.h"
#include "StEmcUtil/StEmcPosition.h"
#include "StMessMgr.h"
#include "StEmcUtil/emcDetectorName.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "tables/St_MagFactor_Table.h"
#include "stdlib.h"
#include "StThreeVectorD.hh"*/
#include "time.h"
#include "TDatime.h"

#include "StDbLib/StDbManager.hh"
#include "StDbLib/StDbTable.h"
#include "StDbLib/StDbConfigNode.hh"

#include "tables/St_smdPed_Table.h"

ClassImp(StSmdPedMaker);

//_____________________________________________________________________________
StSmdPedMaker::StSmdPedMaker(const char *name):StMaker(name)
{
  mPedDate = 0;
  mPedTime = 0;
  mPedInterval = 6; //in hours
  mMinEvents = 10000;
  mSaveToDB = kFALSE;
  mGeo[0] = StEmcGeom::instance("bsmde");
  mGeo[1] = StEmcGeom::instance("bsmdp");
}
//_____________________________________________________________________________
StSmdPedMaker::~StSmdPedMaker()
{
}
//_____________________________________________________________________________
Int_t StSmdPedMaker::Init()
{
  ZeroAll();
	return StMaker::Init();
}
//_____________________________________________________________________________
Int_t StSmdPedMaker::Make()
{
  cout <<"*** SMD Pedestal Maker\n";
  cout <<"Event Date = "<<GetDate()<<"  Time = "<<GetTime()<<endl;
	TStopwatch clock;
  clock.Start();
		
	if(!GetEvent()) return kStWarn;
	Int_t date = GetDate();
  Int_t time = GetTime();
  
  if(mPedStatus==-1)
  {
    Float_t dt = GetDeltaTime(date,time,mPedDate,mPedTime);
    if(dt>mPedInterval) 
    { 
      mPedDate = date; mPedTime = time; mPedStatus = 0;
      cout <<"New pedestal round is starting for SMD...\n"; 
    }
    else 
    {
      cout <<"last pedestal time for SMD = "<<mPedDate<<"  "<<mPedTime<<endl;
      cout <<"time left for new pedestal round for SMD = "<<dt<<" hours"<<endl;
    }
  }
  if(mPedStatus>=0)
  {
    FillPedestal();
    if(mNEvents>mMinEvents) mPedStatus=1;
    if(mPedStatus==1)
    {
      CalculatePedestals();
      if(mSaveToDB) SavePedestals(mPedDate,mPedTime);
      ZeroAll();
    }
  }

  clock.Stop();
  cout <<"Time to run StSmdPedMaker::Make() real = "<<clock.RealTime()<<"  cpu = "<<clock.CpuTime()<<" \n";
	
  return kStOK;
}
//_____________________________________________________________________________
Int_t StSmdPedMaker::Finish()
{
  return kStOk;
}
//_____________________________________________________________________________
Bool_t StSmdPedMaker::GetEvent()
{  
  StEvent *Event=(StEvent*)GetInputDS("StEvent");
  if(!Event) return kFALSE;
  mEmc=Event->emcCollection();
  if(!mEmc) return kFALSE;    	
	StL0Trigger* trg = Event->l0Trigger();
	Int_t trigger=0;
	if(trg) trigger = trg->triggerWord();
	if(trigger!=8192 && trigger!=4096) return kFALSE;
	return kTRUE;
}
//_____________________________________________________________________________
void StSmdPedMaker::FillPedestal()
{ 
  for(Int_t i=0; i<2; i++)
  {  
    StDetectorId id = static_cast<StDetectorId>(i+2+kBarrelEmcTowerId);
    StEmcDetector* detector=mEmc->detector(id);
    if(detector) for(UInt_t j=1;j<121;j++)
    {
      StEmcModule* module = detector->module(j);
      StSPtrVecEmcRawHit& rawHit=module->hits();
      for(Int_t k=0;k<rawHit.size();k++)
      {
        Int_t m = rawHit[k]->module();
        Int_t e = rawHit[k]->eta();
        Int_t s = rawHit[k]->sub();
        Float_t adc = (Float_t)rawHit[k]->adc();
        Int_t cap = rawHit[k]->calibrationType();
        Int_t id;
        mGeo[i]->getId(m,e,s,id);
        Int_t c = 0;
        if(cap==124) c=1;
        if(cap==125) c=2;
        if(adc>0)
        {
          mSmdPedX[i][c][id-1]+=adc;
          mSmdPedX2[i][c][id-1]+=adc*adc;
          mSmdPedSum[i][c][id-1]++;
        }
      }
    }
  }
  mNEvents++;
} 
//_____________________________________________________________________________
void StSmdPedMaker::CalculatePedestals()
{  
  for(Int_t i=0;i<2;i++)
    for(Int_t j=0;j<3;j++)
      for(Int_t k=0;k<18000;k++)
      if(mSmdPedSum[i][j][k]>0)
      {
        mSmdPed[i][j][k] = mSmdPedX[i][j][k]/mSmdPedSum[i][j][k];
        mSmdRMS[i][j][k] = sqrt(mSmdPedX2[i][j][k]/mSmdPedSum[i][j][k]-mSmdPed[i][j][k]*mSmdPed[i][j][k]);
      }
}
//_____________________________________________________________________________
void StSmdPedMaker::ZeroAll()
{  
  for(Int_t i=0;i<2;i++)
    for(Int_t j=0;j<3;j++)
      for(Int_t k=0;k<18000;k++)
      {
        mSmdPedX[i][j][k] = 0;
        mSmdPedX2[i][j][k] = 0;
        mSmdPedSum[i][j][k] = 0;
        mSmdPed[i][j][k] = 0;
        mSmdRMS[i][j][k] = 0;
      }
  mNEvents = 0;
 	mPedStatus=-1;
}
//_____________________________________________________________________________
Float_t StSmdPedMaker::GetDeltaTime(Int_t date, Int_t time, Int_t date0, Int_t time0)
{
  struct tm t0;
  struct tm t1;
    
  Int_t year  = (Int_t)(date/10000);
  t0.tm_year=year-1900;
  Int_t month = (Int_t)(date-year*10000)/100; 
  t0.tm_mon=month-1;
  Int_t day   = (Int_t)(date-year*10000-month*100);
  t0.tm_mday=day;
  Int_t hour  = (Int_t)(time/10000);
  t0.tm_hour=hour;
  Int_t minute= (Int_t)(time-hour*10000)/100;
  t0.tm_min=minute;
  Int_t second= (Int_t)(time-hour*10000-minute*100);
  t0.tm_sec=second;
  t0.tm_isdst = -1;
  
  year  = (Int_t)(date0/10000);
  t1.tm_year=year-1900;
  month = (Int_t)(date0-year*10000)/100; 
  t1.tm_mon=month-1;
  day   = (Int_t)(date0-year*10000-month*100);
  t1.tm_mday=day;
  hour  = (Int_t)(time0/10000);
  t1.tm_hour=hour;
  minute= (Int_t)(time0-hour*10000)/100;
  t1.tm_min=minute;
  second= (Int_t)(time0-hour*10000-minute*100);
  t1.tm_sec=second;
  t1.tm_isdst = -1;

  time_t ttime0=mktime(&t0);
  time_t ttime1=mktime(&t1);
  double diff=difftime(ttime0,ttime1);
  
  //cout <<"startdate = "<<startDate<<" now = "<<d<<"  dt = "<<diff/60<<endl;
  //cout <<diff<<"  "<<diff/60.<<endl;
  Float_t ddd=diff/3600.;
  return ddd;

}
//_____________________________________________________________________________
void StSmdPedMaker::SavePedestals(Int_t date, Int_t time)
{
	TDatime *tt = new TDatime(date,time);
  TString timestamp = tt->AsSQLString();
  delete tt;
          
  StDbManager* mgr=StDbManager::Instance();
	StDbConfigNode* node=mgr->initConfig(dbCalibrations,dbEmc);
		  
  TString tn[2];
  tn[0]="bsmdePed";
  tn[1]="bsmdpPed";
  	
  for(Int_t i=0;i<2;i++)
  {
    smdPed_st tnew;
	  for(Int_t j=0;j<18000;j++)
	  {
      Int_t status = 1;
      for(Int_t k = 0;k<3;k++)
      {
        short int ped = (short int)(100*mSmdPed[i][k][j]);
        short int rms = (short int)(100*mSmdRMS[i][k][j]);
        if(ped==0 && rms==0) status =0;
        tnew.AdcPedestal[j][k] = ped;
        tnew.AdcPedestalRMS[j][k] = rms;
      }
      tnew.Status[j] = (char) status;
	  }		
		/////////////////////
	  StDbTable* tab=node->addDbTable(tn[i].Data());
	  tab->SetTable((char*)&tnew,1);
	  mgr->setStoreTime(timestamp.Data());
	  mgr->storeDbTable(tab);
  }
}





