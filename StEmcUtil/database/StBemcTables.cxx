// 
// $Id: StBemcTables.cxx,v 1.1 2004/12/21 12:54:30 suaide Exp $
// $Log: StBemcTables.cxx,v $
// Revision 1.1  2004/12/21 12:54:30  suaide
// moved StBemcTables to StEmcUtil/database
// added interface to trigger tables in offline DB
//
// Revision 1.2  2004/10/19 17:53:00  suaide
// code clean up
//
// Revision 1.1  2004/10/18 18:20:07  suaide
// New Maker. Will replace StEmcADCtoEMaker in production.
// It reads only DAQ structures. Output is StEvent.
//
#include "StBemcTables.h"
#include "Stiostream.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "TString.h"

ClassImp(StBemcTables)

//_____________________________________________________________________________
/* 
   Default constructor. Set Initial values for some variables
*/
StBemcTables::StBemcTables():TObject()
{  
  mBtowP = NULL;
  mBprsP = NULL;
  mSmdeP = NULL;
  mSmdpP = NULL;
  mBtowS = NULL;
  mBprsS = NULL;
  mSmdeS = NULL;
  mSmdpS = NULL;
  mBtowC = NULL;
  mBprsC = NULL;
  mSmdeC = NULL;
  mSmdpC = NULL;
  mBtowG = NULL;
  mBprsG = NULL;
  mSmdeG = NULL;
  mSmdpG = NULL;
  mTrigS = NULL;
  mTrigP = NULL;
}
//_____________________________________________________________________________
/*! 
   Default destructor
*/
StBemcTables::~StBemcTables()
{
}
//_____________________________________________________________________________
/*!
  Get pedestals and status tables.  
*/
void StBemcTables::loadTables(Int_t det, TDataSet *DB)
{  
  if(det==BTOW)       { mBtowP = NULL; mBtowS = NULL; mBtowC = NULL; mBtowG = NULL; }
  else if(det==BPRS)  { mBprsP = NULL; mBprsS = NULL; mBprsC = NULL; mBprsG = NULL; }
  else if(det==BSMDE) { mSmdeP = NULL; mSmdeS = NULL; mSmdeC = NULL; mSmdeG = NULL; }
  else if(det==BSMDP) { mSmdpP = NULL; mSmdpS = NULL; mSmdpC = NULL; mSmdpG = NULL; }
  
  TString TableName;
  if(DB)
  {
    /////////////////////////////////////////////////
    TableName = detname[det-1]+"Ped";
    if(det==BTOW || det==BPRS)
    {
      St_emcPed* ped=(St_emcPed*)DB->Find(TableName.Data());
      if(ped)
      {
        if(det==BTOW) mBtowP = ped->GetTable();
        else       mBprsP = ped->GetTable();
      } 
    }
    else
    {
      St_smdPed* ped=(St_smdPed*)DB->Find(TableName.Data());
      if(ped)
      {
        if(det==BSMDE) mSmdeP = ped->GetTable();
        else       mSmdpP = ped->GetTable();
      }
    }
     
    /////////////////////////////////////////////////
    TableName = detname[det-1]+"Status";
    if(det==BTOW || det==BPRS)
    {
      St_emcStatus* status=(St_emcStatus*)DB->Find(TableName.Data());
      if(status)
      {
        if(det==BTOW) mBtowS = status->GetTable();
        else       mBprsS = status->GetTable();
      }
    }
    else
    {
      St_smdStatus* status=(St_smdStatus*)DB->Find(TableName.Data());
      if(status)
      {
        if(det==BSMDE) mSmdeS = status->GetTable();
        else       mSmdpS = status->GetTable();
      }
    }
      
    /////////////////////////////////////////////////
    TableName = detname[det-1]+"Calib";
    if(det==BTOW || det==BPRS)
    {
      St_emcCalib* calib=(St_emcCalib*)DB->Find(TableName.Data());
      if(calib)
      {
        if(det==BTOW) mBtowC = calib->GetTable();
        else       mBprsC = calib->GetTable();
      }
    }
    else
    {
      St_smdCalib* calib=(St_smdCalib*)DB->Find(TableName.Data());
      if(calib)
      {
        if(det==BSMDE) mSmdeC = calib->GetTable();
        else       mSmdpC = calib->GetTable();
      }
    }
      
    /////////////////////////////////////////////////
    TableName = detname[det-1]+"Gain";
    if(det==BTOW || det==BPRS)
    {
      St_emcGain* gain=(St_emcGain*)DB->Find(TableName.Data());
      if(gain)
      {
        if(det==BTOW) mBtowG = gain->GetTable();
        else       mBprsG = gain->GetTable();
      }
    }
    else
    {
      St_smdGain* gain=(St_smdGain*)DB->Find(TableName.Data());
      if(gain)
      {
        if(det==BSMDE) mSmdeG = gain->GetTable();
        else       mSmdpG = gain->GetTable();
      }
    }
  }
  return;  
}
//_____________________________________________________________________________
/*!
  Get pedestals and status tables.  
*/
void StBemcTables::loadTables(StMaker* maker)
{
  for(Int_t det = 1; det<=MAXDETBARREL;det++)
  {	  
    TString DbName = "Calibrations/emc/y3"+detname[det-1];
    TDataSet *DB   = maker->GetInputDB(DbName.Data());
    if(DB) loadTables(det,DB);
  }
  
  TString DbName = "Calibrations/emc/trigger";
  TDataSet *DB   = maker->GetInputDB(DbName.Data());
  
  mTrigS = NULL;
  mTrigP = NULL;
  if(DB)
  {
    St_emcTriggerStatus* s = (St_emcTriggerStatus*)DB->Find("bemcTriggerStatus");
    if(s) mTrigS = s->GetTable();
    St_emcTriggerPed* p = (St_emcTriggerPed*)DB->Find("bemcTriggerPed");
    if(s) mTrigP = p->GetTable();
  }
  
}  
//_____________________________________________________________________________
/*!
  Get pedestal value  
*/
void StBemcTables::getPedestal(Int_t det, Int_t id, Int_t CAP,Float_t& P, Float_t& R)
{
  P = 0;
  R = 0;
  if(det==BTOW && mBtowP) 
  {
    P = ((Float_t)mBtowP[0].AdcPedestal[id-1])/100;
    R = ((Float_t)mBtowP[0].AdcPedestalRMS[id-1])/100;
    return;
  }
  if(det==BPRS && mBprsP) 
  {
    P = ((Float_t)mBprsP[0].AdcPedestal[id-1])/100;
    R = ((Float_t)mBprsP[0].AdcPedestalRMS[id-1])/100;
    return;
  }
  if(det==BSMDE && mSmdeP) 
  {
    Int_t C = 0;
    if(CAP==CAP1) C = 1;
    if(CAP==CAP2) C = 2;
    P = ((Float_t)mSmdeP[0].AdcPedestal[id-1][C])/100;
    R = ((Float_t)mSmdeP[0].AdcPedestalRMS[id-1][C])/100;
    return;
  }
  if(det==BSMDP && mSmdpP) 
  {
    Int_t C = 0;
    if(CAP==CAP1) C = 1;
    if(CAP==CAP2) C = 2;
    P = ((Float_t)mSmdpP[0].AdcPedestal[id-1][C])/100;
    R = ((Float_t)mSmdpP[0].AdcPedestalRMS[id-1][C])/100;
    return;
  }
  return;  
}
//_____________________________________________________________________________
/*!
  Get status  
*/
void StBemcTables::getStatus(Int_t det, Int_t id, Int_t& S)
{
  S = STATUS_OK;
  if(det==BTOW && mBtowS) { S = (Int_t)mBtowS[0].Status[id-1];return;}
  if(det==BPRS && mBprsS) { S = (Int_t)mBprsS[0].Status[id-1];return;}
  if(det==BSMDE && mSmdeS) { S = (Int_t)mSmdeS[0].Status[id-1];return;}
  if(det==BSMDP && mSmdpS) { S = (Int_t)mSmdpS[0].Status[id-1];return;}
  return;  
}
//_____________________________________________________________________________
/*!
  Get Gain 
*/
void StBemcTables::getGain(Int_t det, Int_t id, Float_t& G)
{
  G = 1;
  if(det==BTOW && mBtowG) { G = (Float_t)mBtowG[0].Gain[id-1];return;}
  if(det==BPRS && mBprsG) { G = (Float_t)mBprsG[0].Gain[id-1];return;}
  if(det==BSMDE && mSmdeG) { G = (Float_t)mSmdeG[0].Gain[id-1];return;}
  if(det==BSMDP && mSmdpG) { G = (Float_t)mSmdpG[0].Gain[id-1];return;}
  return;  
}
//_____________________________________________________________________________
/*!
  Get Calib 
*/
void StBemcTables::getCalib(Int_t det, Int_t id, Int_t power, Float_t& C)
{
  C = 0;
  if(det==BTOW && mBtowC) { C = (Float_t)mBtowC[0].AdcToE[id-1][power];return;}
  if(det==BPRS && mBprsC) { C = (Float_t)mBprsC[0].AdcToE[id-1][power];return;}
  if(det==BSMDE && mSmdeC) { C = (Float_t)mSmdeC[0].AdcToE[id-1][power];return;}
  if(det==BSMDP && mSmdpC) { C = (Float_t)mSmdpC[0].AdcToE[id-1][power];return;}
  return;  
}
//_____________________________________________________________________________
/*!
  Get Trigger patch staus 
*/
void StBemcTables::getTriggerPatchStatus(Int_t patch, Int_t& STATUS)
{
  STATUS = 0;
  if(mTrigS && patch>=0 && patch<NBEMCTRIGGERTOWER) 
    STATUS = (Int_t)mTrigS[0].PatchStatus[patch];
}
//_____________________________________________________________________________
/*!
  Get Trigger high tower staus 
*/
void StBemcTables::getTriggerHighTowerStatus(Int_t hightower, Int_t& STATUS)
{
  STATUS = 0;
  if(mTrigS && hightower>=0 && hightower<NBEMCTRIGGERTOWER) 
    STATUS = (Int_t)mTrigS[0].HighTowerStatus[hightower];
}

//_____________________________________________________________________________
/*!
  Get Trigger single tower staus 
*/
void StBemcTables::getTriggerTowerStatus(Int_t crate,Int_t index, Int_t& STATUS)
{
  STATUS = 0;
  if(mTrigS && crate>0 && crate<=MAXCRATES && index>=0 && index<NTOWERSPERCRATE) 
    STATUS = (Int_t)mTrigS[0].TowerStatus[crate-1][index];
}
//_____________________________________________________________________________
/*!
  Get Trigger pedestal value 
*/
void StBemcTables::getTriggerPedestal(Int_t crate,Int_t index, Float_t& PEDESTAL)
{
  PEDESTAL = 0;
  if(mTrigP && crate>0 && crate<=MAXCRATES && index>=0 && index<NTOWERSPERCRATE) 
    PEDESTAL = ((Float_t)mTrigP[0].Ped[crate-1][index])/100.0;
}
//_____________________________________________________________________________
/*!
  Get Trigger pedestal value 
*/
void StBemcTables::getTriggerBitConv(Int_t crate,Int_t patch, Int_t& BIT)
{
  BIT = 0;
  if(mTrigP && crate>0 && crate<=MAXCRATES && patch>=0 && patch<NPATCHESPERCRATE) 
    BIT = (Int_t)mTrigP[0].BitConversionMode[crate-1][patch];
}





