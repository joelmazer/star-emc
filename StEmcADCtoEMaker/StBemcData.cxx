#include "StBemcData.h"
#include "StDaqLib/EMC/StEmcDecoder.h"
#include "iostream.h"
 
ClassImp(StBemcData)
 
StBemcData::StBemcData(char* name):TDataSet(name)
{
  mDecoder = NULL;
}
StBemcData::~StBemcData()
{
  if(mDecoder) delete mDecoder;
} 
Bool_t StBemcData::getTDCStatus(Int_t c)
{
	if(c<0 || c>29) return kFALSE;
	if(c<3  && EventDate<20011105) return kTRUE;
	if(c<5  && EventDate>=20011105 && EventDate<20020301) return kTRUE;
	if(c<15 && EventDate>=20020301) return kTRUE;
	return kFALSE;
} 
Bool_t StBemcData::getSMDStatus(Int_t c)
{
	if(c<0 || c>7) return kFALSE;
	if(c<1  && EventDate<20011105) return kTRUE;
	if(c<2  && EventDate>=20011105 && EventDate<20020301) return kTRUE;
	if(c<4 && EventDate>=20020301) return kTRUE;
	return kFALSE;
} 
Bool_t StBemcData::checkTDC(Int_t i)
{
	Bool_t ok = kTRUE;
  if(!mDecoder) mDecoder = new StEmcDecoder(EventDate,EventTime);
	Int_t crate;
  mDecoder->GetTowerCrateFromTDC(i,crate);
  if(TDCError[i]!=0) ok = kFALSE;
	if(TDCCrateId[i]!=crate) ok = kFALSE;
	if(TDCCount[i]!=164) ok = kFALSE;
  return ok;
}
void StBemcData::validateData()
{
	//towers first
	Bool_t ok = kTRUE;
	if(!TowerPresent) ok = kFALSE;
	if(TDCErrorFlag==1) ok = kFALSE;
  Int_t nbad = 0;
	for(Int_t i=0;i<30;i++) if(getTDCStatus(i)) if(!checkTDC(i)) nbad++;
  if(nbad==0) ok = kTRUE;
	ValidTowerEvent = ok;
	
	//SMD
	ok = kTRUE;
	if(!SMDPresent) ok = kFALSE;
	if(SMDErrorFlag ==1) ok = kFALSE;
	ValidSMDEvent = ok;
}
void StBemcData::printTower()
{
	cout <<"TOWER DATA -----------------------------------------------------------\n";
	cout <<"RunNumber = "<<RunNumber<<endl;
	cout <<"TriggerWord = "<<TriggerWord<<endl;
	cout <<"EventDate = "<<EventDate<<endl;
	cout <<"EventTime = "<<EventTime<<endl;
	cout <<"ValidTowerEvent = "<<(Int_t)ValidTowerEvent<<endl;
	cout <<"EventNumber = "<<EventNumber<<endl;
	cout <<"TowerByteCount = "<<TowerByteCount<<endl;
	cout <<"NTowerHits = "<<NTowerHits<<endl;
	cout <<"TDCErrorFlag = "<<TDCErrorFlag<<endl;
	cout <<"NTDCChannels = "<<NTDCChannels<<endl;
	for(Int_t i=0;i<30;i++)
	  cout <<"  TDC = "<<i<<"  Error = "<<TDCError[i]<<"  Token = "<<TDCToken[i]<<"  Trigger = "<<TDCTrigger[i]
		     <<"  Crate = "<<TDCCrateId[i]<<"  Count = "<<TDCCount[i]<<endl;
	for(Int_t k=0;k<4800;k++)
	{
	  cout <<" ID = "<<k<<" STAT = "<<(Int_t)TowerStatus[k]<<" ADC = "<<TowerADC[k]<<" E = "<<TowerEnergy[k]<<" |";
		if((k+1)%3==0) cout <<endl;
	}
}
void StBemcData::printSMD()
{
}
