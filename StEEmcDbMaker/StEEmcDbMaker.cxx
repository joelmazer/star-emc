// *-- Author : Jan Balewski
// 
// $Id: StEEmcDbMaker.cxx,v 1.43 2005/01/24 05:08:26 balewski Exp $
 

#include <time.h>
#include <string.h>

#include <TDatime.h>

#include "StChain.h"
#include "St_DataSetIter.h"

#include "St_db_Maker/St_db_Maker.h" // tmp to ovveride time stamp

#include "StEventTypes.h"

#include "StEEmcDbMaker.h"

#include "StEEmcDbMaker/EEmcDbItem.h"
#include "StEEmcDbMaker/EEmcDbCrate.h"
#include "StEEmcUtil/EEfeeRaw/EEname2Index.h" 


#include "tables/St_eemcDbADCconf_Table.h"
#include "tables/St_eemcDbPMTcal_Table.h"
#include "tables/St_eemcDbPMTname_Table.h"
#include "tables/St_eemcDbPIXcal_Table.h"
#include "tables/St_eemcDbPMTped_Table.h"
#include "tables/St_eemcDbPMTstat_Table.h"
#include "tables/St_kretDbBlobS_Table.h"
#include "cstructs/eemcConstDB.hh"
#include "cstructs/kretConstDB.hh"

#include <StMessMgr.h>

/*
http://www.star.bnl.gov/STAR/comp/pkg/dev/StRoot/StUtilities/doc/StMessMgr.html
gMessMgr->Message("This is an error message.","E");

The message type is specified with a single letter:
"E" = "Error"
"I" = "Info"
"W" = "Warning"
"D" = "Debug"
"Q" = "QAInfo"

gMessMgr->Message("","W") << big_num << " seems too big." << endm;

 char* myText = "hello";
  gMessMgr->Message(myText);


*/

ClassImp(StEEmcDbMaker)

//________________________________________________________
//________________________________________________________
StEEmcDbMaker::StEEmcDbMaker(const char *name):StMaker(name){
  gMessMgr->Message("","D") <<GetName()<<endm;
  mfirstSecID=mlastSecID=mNSector=0;
  myTimeStampDay=0;
  myTimeStampUnix=0;

  //................ allocate memory for lookup tables
  byIndex=new  EEmcDbItem[EEindexMax];
  
  byCrate=new  EEmcDbItem ** [MaxAnyCrate];
  
  int i;
  for(i=0;i<MaxAnyCrate;i++){
    byCrate[i]=NULL;
    if(i==0 || (i>MaxTwCrateID && i<MinMapmtCrateID) ) continue; // to save memory for nonexisting crates
    byCrate[i]=new EEmcDbItem * [MaxAnyCh];
    memset(byCrate[i],0,sizeof(EEmcDbItem *)*MaxAnyCh);// clear all pointers
  }


  mDbFiber=0; nFiber=0;
  setDBname("Calibrations/eemc");

  mAsciiDbase = "";

  mxChGain=8;
  chGainL = new TString [mxChGain];
  nChGain=0;

  mxChMask=8;
  chMaskL = new TString [mxChMask];
  nChMask=0;

  setThreshold(5.0);  // defines threshold for ADCs
  // should be +2 or +3 sigma in the future
  //  setPreferedFlavor("onlped","eemcPMTped"); // tmp for tests,JB

}


//________________________________________________________
//________________________________________________________
//_______________________________________________________
StEEmcDbMaker::~StEEmcDbMaker(){

  if( mNSector) {
    delete [] mDbADCconf;
    delete [] mDbPMTcal;
    delete [] mDbPMTname;
    delete [] mDbPIXcal;
    delete [] mDbPMTped;
    delete [] mDbPMTstat;
    delete [] mDbsectorID;
  }

  delete mDbFiber;
  delete [] chGainL;
  delete [] chMaskL;
}

//________________________________________________________
//________________________________________________________
void  StEEmcDbMaker::setTimeStampDay( int tD) {

  if(myTimeStampDay) {
    gMessMgr->Message("","E") <<GetName()<<"::setTimeStampDay() Operator Fatal error :  redeclaration ofmyTimeStampDay "<<myTimeStampDay<<endm;
    exit(1);
  }
  myTimeStampDay=tD;

  int iy=tD/10000; tD%=10000;
  int im=tD/100;  tD%=100;
  int id=tD;
  //  printf("%d %d %d\n",iy,im,id);

  // build Unix time stamp 
  struct tm t;
  memset(&t,0,sizeof t);
  t.tm_sec=0;    /* seconds after the minute - [0,61] */
  t.tm_min=0;      /* minutes after the hour - [0,59] */
  t.tm_hour=0;      /* hours - [0,23] */
  t.tm_mday=id;     /* day of month - [1,31] */
  t.tm_mon=im-1;      /* month of year - [0,11] */
  t.tm_year=iy-1900;     /* years since 1900 */
  t.tm_isdst=0;    /* daylight savings time flag */  
  //int tm_wday;     /* days since Sunday - [0,6] */
  //int tm_yday;     /* days since January 1 - [0,365] */
  myTimeStampUnix =(unsigned int) mktime(&t);

  gMessMgr->Message("","W") <<GetName()<<"::setTimeStampDay() set to "<<myTimeStampDay<<" decoded as "<< ctime((const time_t *)&myTimeStampUnix)<< endm;

}

//------------------
//------------------
void StEEmcDbMaker::setThreshold(float x){
 KsigOverPed=x;
 gMessMgr->Message("","I") <<GetName()<<"::setThres KsigOverPed="<<KsigOverPed<<", threshold=ped+sig*KsigOverPed"<< endm;
}


//------------------
//------------------
void StEEmcDbMaker::setPreferredFlavor(const char *flavor, const char *nameMask){
  strncpy(dbFlavor.flavor,flavor,DbFlavor::mx);
  strncpy(dbFlavor.nameMask,nameMask,DbFlavor::mx);
  gMessMgr->Message("","I") << GetName()<<"::setPreferredFlavor(flav='"<<dbFlavor.flavor <<"', mask='"<< dbFlavor.nameMask<<"')" <<endm;
}



//________________________________________________________
//________________________________________________________
//________________________________________________________
Int_t StEEmcDbMaker::Init(){
  if( mNSector==0) setSectors(1,12);//default

  return StMaker::Init();
}


//________________________________________________________
//________________________________________________________
//________________________________________________________
void StEEmcDbMaker::setSectors(int sec1,int sec2){
  assert(mNSector==0) ; //you can do it just once, no memory realocation implemented

  mfirstSecID=sec1;
  mlastSecID=sec2;
  mNSector=mlastSecID - mfirstSecID+1;

  mDbADCconf=(eemcDbADCconf_st **) new void *[mNSector];
  mDbPMTcal= (eemcDbPMTcal_st  **) new void *[mNSector];
  mDbPMTname=(eemcDbPMTname_st **) new void *[mNSector];
  mDbPIXcal= (eemcDbPIXcal_st  **) new void *[mNSector];
  mDbPMTped= (eemcDbPMTped_st  **) new void *[mNSector];
  mDbPMTstat=(eemcDbPMTstat_st **) new void *[mNSector];
  mDbsectorID=  new int [mNSector];


  clearItemArray();

  gMessMgr->Message("","I") << GetName()<<":: Use sectors from "<< mfirstSecID<<" to "<< mlastSecID<<endm;

}

//__________________________________________________
//__________________________________________________
//__________________________________________________

const EEmcDbCrate* StEEmcDbMaker::getFiber(int icr) {
  assert(icr>=0);
  assert(icr<nFiber);
  return mDbFiber+icr;
}



//--------------------------------------------------
//--------------------------------------------------
void StEEmcDbMaker::clearItemArray(){
  //printf("%s::clearItemArray()\n",GetName());
  nFound=0;

  int i;

  for(i=0; i<EEindexMax; i++)
    byIndex[i].clear();

  int j;
  for(i=0;i<MaxAnyCrate;i++) {
    if(byCrate[i]==NULL) continue;
    for(j=0;j<MaxAnyCh;j++)
      byCrate[i][j]=0;
  }

  memset(byStrip,0,sizeof(byStrip));
  
  if(mDbFiber) delete [] mDbFiber;
  nFiber=0;
  mDbFiberConfBlob=0;


  nFound=0;
  mDbADCconf[0]=0;
  for(i=0; i<mNSector; i++) {// clear pointers old DB tables
     mDbADCconf [i]=0;
     mDbPMTcal  [i]=0;
     mDbPMTname [i]=0;
     mDbPIXcal  [i]=0;
     mDbPMTped  [i]=0;
     mDbPMTstat [i]=0;
    mDbsectorID[i]=-1;
  }

}


//__________________________________________________
//__________________________________________________
//__________________________________________________

Int_t  StEEmcDbMaker::InitRun  (int runNumber)
{
  // Reloads database for each occurence of a new run number.
  // If an ascii file has been loaded, via setAsciiDatabase(),
  // issue a warning and return.

  if ( mAsciiDbase.Length() > 0 ) {
    Warning("InitRun","Database not reloaded for run number %i, values taken from %s",runNumber,mAsciiDbase.Data());
    return kStOK;
  }
  gMessMgr->Message("","I") << GetName()<<"::InitRun  :use(flav='"<<dbFlavor.flavor <<"', mask='"<< dbFlavor.nameMask<<"')" <<endm;

 

  clearItemArray();
  mRequestDataBase();

 //............  reload all lookup tables ...............
  int is;
  for(is=0; is< mNSector; is++) {
    if (  mOptimizeMapping(is) ) goto end;// on fatal error
    mOptimizeOthers(is); 
  }

  mOptimizeFibers();

  // overload some of DB info
  for(is=0;is<nChGain;is++) changeGainsAction(chGainL[is].Data());


  // exportAscii(); //tmp

 gMessMgr->Message("","I") << GetName()<<"::InitRun()  Found "<< nFound<<" EEMC related tables "<<endm;

 end:
  return StMaker::InitRun(runNumber);
}  


//__________________________________________________
//__________________________________________________

void  StEEmcDbMaker::mRequestDataBase(){

  // printf("%s::reloadDb using TimeStamp from 'StarDb'=%p or 'db'=%p \n",GetName(),(void*)GetMaker("StarDb"),(void*)GetMaker("db"));
  

  St_db_Maker* mydb = (St_db_Maker*)GetMaker("StarDb");
  if(mydb==0) mydb = (St_db_Maker*)GetMaker("db");
  assert(mydb);

  if(myTimeStampDay==0) { // use oryginal timestamp of event   
    
#if 0
    StEvent *stEvent= (StEvent *) GetInputDS("StEvent"); 
    assert(stEvent);     
    printf("StEvent time=%d, ID=%d, runID=%d\n",(int)stEvent->time(),(int)stEvent->id(),(int)stEvent->runId());
#endif
    
    StEvtHddr* fEvtHddr = (StEvtHddr*)GetDataSet("EvtHddr");
    
    gMessMgr->Message("","I") << GetName()<<"::RequestDataBase(), use EvtHddr actual event time stamp="<< (int)fEvtHddr->GetUTime()<< " , yyyy/mm/dd="<< fEvtHddr->GetDate()<<" hh/mm/ss="<<fEvtHddr->GetTime()<<endm;  

    //  int time0; // (sec) GMT of the first event
    //  time0=fEvtHddr->GetUTime( ); //<<==== this is used by DB

  } else { // WARN only if you wish to overwrite the global time stamp 
    gMessMgr->Message("","W") << GetName()<<"::RequestDataBase() replace  TimeStampDay to "<< myTimeStampDay<<endm;
    mydb->SetDateTime(myTimeStampDay,0); // set ~day & ~hour by hand
  }
  // mydb->SetDateTime(20021201,0); // set ~day & ~hour by hand

  gMessMgr->Message("","I") << GetName()<<"::RequestDataBase(), access DB='"<<dbName<<"  use timeStamp="<< mydb->GetDateTime().AsString()<<endm;

  int ifl;
  TString mask="";
  for(ifl=0;ifl<2;ifl++) { // loop over flavors
    if(ifl==1) {
      if( dbFlavor.flavor[0]==0) continue; // drop flavor change
      gMessMgr->Message("","I") << GetName()<<"::RequestDataBase()-->ifl=%d try flavor='"<<dbFlavor.flavor <<"' for  mask='"<< dbFlavor.nameMask<<"')" <<endm;
      
      SetFlavor(dbFlavor.flavor,dbFlavor.nameMask);
      mask=dbFlavor.nameMask;
    }
    
    TDataSet *eedb=GetDataBase(dbName );
    if(eedb==0) {
      gMessMgr->Message("","E") << GetName()<<"::RequestDataBase()  Could not find dbName ="<<dbName <<endm;
      return ;
      // down-stream makers should check for presence of dataset
    }
    //eedb->ls(2);  
        
    int is;
    for(is=0; is< mNSector; is++) {
      int secID=is+mfirstSecID;
      
      mDbsectorID[is]=secID;
      getTable<St_eemcDbADCconf,eemcDbADCconf_st>(eedb,secID,"eemcADCconf",mask,mDbADCconf+is);
   
      getTable<St_eemcDbPMTcal,eemcDbPMTcal_st>(eedb,secID,"eemcPMTcal",mask,mDbPMTcal+is);   

      getTable<St_eemcDbPMTname,eemcDbPMTname_st>(eedb,secID,"eemcPMTname",mask,mDbPMTname+is);   

      getTable<St_eemcDbPIXcal,eemcDbPIXcal_st>(eedb,secID,"eemcPIXcal",mask,mDbPIXcal+is);   
      
      getTable<St_eemcDbPMTped,eemcDbPMTped_st>(eedb,secID,"eemcPMTped",mask,mDbPMTped+is);
      
      getTable<St_eemcDbPMTstat,eemcDbPMTstat_st>(eedb,secID,"eemcPMTstat",mask,mDbPMTstat+is);
      
    } // end of loop over sectors
    
    // misc tables 
    
    getTable<St_kretDbBlobS,kretDbBlobS_st>(eedb,13,"eemcCrateConf",mask,&mDbFiberConfBlob);
    if(mDbFiberConfBlob) {
      gMessMgr->Message("","D") << GetName()<<"::RequestDataBase()  eemcCrateConf dump "<<endm;
      //  gMessMgr->Message(mDbFiberConfBlob->dataS,"D"); 
    }

  }// end of loop over flavors
 
#if 0

  //tmp
  TDatime aa1=mydb->GetDateTime(); 
  if (aa1.GetDate()<20040101) xxx;

#endif

}

//--------------------------------------------------
//--------------------------------------------------
int  StEEmcDbMaker::mOptimizeMapping(int is){

  gMessMgr->Message("","I") <<"\n  EEDB  conf ADC map for sector="<< mDbsectorID[is] <<endm;
  assert(mDbsectorID[is]>0);
  
  eemcDbADCconf_st *t= mDbADCconf[is];
  
  if(t==0) return -2; // it is fatal error
  gMessMgr->Message("","I") <<"      EEDB chanMap="<< t->comment <<endm;
  
  int j;
  for(j=0;j<EEMCDbMaxAdc; j++) { // loop over channels
    char *name=t->name+j*EEMCDbMaxName;
    
    if(*name==EEMCDbStringDelim) continue;
    
    //printf("%d '%s'  %d %d\n",j,name,t->crate[j],t->channel[j]);
    // printf("%d   %d %d\n",j,t->crate[j],t->channel[j]);
    
    int key=EEname2Index(name);
    assert(key>=0 && key<EEindexMax);
    EEmcDbItem *x=&byIndex[key];
    if(!x->isEmpty()) {
      x->print();
      assert(x->isEmpty());
    }
    x->crate=t->crate[j];
    x->chan=t->channel[j];
    x->setName(name);
    x->key=key;
    x->setDefaultTube(MinMapmtCrateID);
    // x->print();
    
    assert(x->crate>=0 && x->crate<MaxAnyCrate);
    assert(x->chan>=0 && x->chan<MaxAnyCh);
    assert(byCrate[x->crate]);// ERROR: duplicated crate ID from DB
    if(byCrate[x->crate][x->chan]) {
      gMessMgr->Message("","E") << GetName()<<"::Fatal Error of eemc DB records: the same crate="<<x->crate<<", ch="<<x->chan<<" entered twice \n\n Whole EEMC DB is erased from memory,\n  all data will be ignored,  FIX IT !, JB\n\n"<<endm;
      byCrate[x->crate][x->chan]->print(); // first time
      x->print(); // second time
      clearItemArray();
      return -1;
    }
    byCrate[x->crate][x->chan]=x;
    if(x->isSMD()) byStrip[x->sec-1][x->plane-'U'][x->strip-1]=x;
  }

  return 0;
  
}


//--------------------------------------------------
//--------------------------------------------------
void StEEmcDbMaker::mOptimizeOthers(int is){

  int secID= mDbsectorID[is];
  //  printf("\n  optimizeDB for sector=%d\n",secID); //tmp
  int ix1,ix2;
  EEindexRange(secID,ix1,ix2);
  
  gMessMgr->Message("","D") <<"\n   EEDB EEindexRange("<<secID<<","<<ix1<<","<<ix2<<")"<<endm;
  
  //  if(dbg)printf(" Size: ped=%d cal=%d name=%d stat=%d \n",sizeof(ped->name)/EEMCDbMaxName,sizeof(cal->name)/EEMCDbMaxName,sizeof(tubeTw->name)/EEMCDbMaxName,sizeof(stat->name)/EEMCDbMaxName);
  
  assert(secID>0);

  eemcDbPMTcal_st  *calT=mDbPMTcal[is];  
  if(calT)  gMessMgr->Message("","I") <<"      EEDB calTw="<<calT->comment<<endm;
  eemcDbPMTname_st *tubeTw=mDbPMTname[is];
  if(tubeTw) gMessMgr->Message("","I") <<"      EEDB tubeTw="<<tubeTw->comment<<endm;

  eemcDbPIXcal_st  *calM= mDbPIXcal[is];
  if(calM) gMessMgr->Message("","I") <<"      EEDB calMAPMT="<<calM->comment<<endm;

  eemcDbPMTped_st  *ped=mDbPMTped[is];
  if(ped) gMessMgr->Message("","I") <<"      EEDB ped="<<ped->comment<<endm;

  eemcDbPMTstat_st *stat=mDbPMTstat[is];
  if(stat) gMessMgr->Message("","I") <<"      EEDB stat="<<stat->comment<<endm;
  
  int key; 
  for(key=ix1;key<ix2; key++) { // loop  in this sector
    EEmcDbItem *x=byIndex+key;
    if(x->isEmpty()) continue;
    char *name=x->name;

    if(ped) { // pedestals 
      int j;
      int mx=sizeof(ped->name)/EEMCDbMaxName;
      for(j=0;j<mx; j++) {
	char *name1=ped->name+j*EEMCDbMaxName;
	if(strncmp(name,name1,strlen(name))) continue;
	x->ped=ped->ped[j];
	x->thr=ped->ped[j]+KsigOverPed*ped->sig[j];
	//printf("%d found %s %d %d\n",j,name,strlen(name),strncmp(name,name1,strlen(name)));
	//x->print();
	break;
      }
    } // end of pedestals


    if(calT&& name[2]=='T') { // calibration for towers only
      int j;
      int mx=sizeof(calT->name)/EEMCDbMaxName;
      for(j=0;j<mx; j++) {
	char *name1=calT->name+j*EEMCDbMaxName;
	if(strncmp(name,name1,strlen(name))) continue;
	x->gain=calT->gain[j];
	break;
      }
    } // end of Tower gains


    if(calM && name[2]!='T') { // calibration for MAPMT
      int j;
      int mx=sizeof(calM->name)/EEMCDbMaxName;

      for(j=0;j<mx; j++) {
	char *name1=calM->name+j*EEMCDbMaxName;
	if(strncmp(name,name1,strlen(name))) continue;
	x->gain=calM->gain[j];
	break;
      }
    } // end of gains

    if(tubeTw && name[2]=='T') { // change tube for towers only
      int j;
      int mx=sizeof(tubeTw->name)/EEMCDbMaxName;
      for(j=0;j<mx; j++) {
	char *name1=tubeTw->name+j*EEMCDbMaxName;
	if(strncmp(name,name1,strlen(name))) continue;
	x->setTube(tubeTw->tubeName+j*EEMCDbMaxName);
	//x->print();
	break;
      }
    } // end of tube

    
    if(stat) { // status
      int j;
      int mx=sizeof(stat->name)/EEMCDbMaxName;
      for(j=0;j<mx; j++) {
	char *name1=stat->name+j*EEMCDbMaxName;
	if(strncmp(name,name1,strlen(name))) continue;
	x->stat=stat->stat[j];
	x->fail=stat->fail[j];
	//x->print();
	break;
      }
    } // end of status


    
  }// end of pixels in this sector


}

//--------------------------------------------------
//--------------------------------------------------
void StEEmcDbMaker::exportAscii(char *fname) const{
  printf("EEmcDb::exportAscii(\'%s') ...\n",fname);

  FILE * fd=fopen(fname,"w");
  assert(fd);
  // fd=stdout;

  int nTot=0;

  fprintf(fd,"# EEmcDb::exportAscii()\ttime stamp   : %d / %s",(int)myTimeStampUnix,
          ctime((const time_t *)&myTimeStampUnix));
  fprintf(fd,"# see StRoot/StEEmcDbMaker/EEmcDbItem::exportAscii()  for definition\n");

  int j;
  
  fprintf(fd,"%d  #fibers: {name,crID,crIDswitch,fiber,nCh,nHead,type,useIt}\n",nFiber);
  for(j=0;j<nFiber;j++) 
    mDbFiber[j].exportAscii(fd);

  fprintf(fd,"#tw/pre/post:  {name,crate,chan,sec,plane,strip,gain,ped,thr,stat,fail,tube,key}\n");
  fprintf(fd,"#or \n");
  fprintf(fd,"#smd: {name,crate,chan,sec,sub,eta,gain,ped,thr,stat,fail,tube,key}\n");
 

  for(j=0;j<EEindexMax; j++) { // loop over channels
    const  EEmcDbItem *x=byIndex+j;
    if(x->isEmpty())continue;
    x->exportAscii(fd);
    nTot++;
  }
  printf("        nTot=%d, done\n",nTot);
  fclose(fd);
}


 
//__________________________________________________
//__________________________________________________
//__________________________________________________

void  StEEmcDbMaker::mOptimizeFibers  (){
  int icr=0;

  if( mDbFiberConfBlob==0) {
    gMessMgr->Message("","E") <<" EEDB EEMC  FiberConf not found"<<endm;
    goto fatal;
  }
  
  assert(mDbFiberConfBlob);
  assert(nFiber==0);
  
  //  printf("dataS='%s'\n",mDbFiberConfBlob->dataS);

  /// if the contents of the mDbFiberConfBlob is the same
  //  for a new runs
  /// the data is not fetched you are left with the stale data
  /// And since strtok overwrites its argument the next call to
  /// mOptimizeFibers fails
  /// Cheers
  ///     Piotr

  {
  static char blobBuffer[KRETmxBlobSlen+1];
  char  *blob = blobBuffer;


  int   len = strlen(mDbFiberConfBlob->dataS); 
  assert(len<KRETmxBlobSlen ); // strnlen() is better but aborts under redhat72,JB
  strncpy(blob,mDbFiberConfBlob->dataS,len);
 
  blob=strtok(blob,";"); // init iterator
  if(strstr(blob,"<ver1>")==0) {
    gMessMgr->Message("","E") <<" EEDB EEMC  FiberConf  missing opening key "<<endm;
    goto fatal;
  }
  
  int i=0;
  while((blob=strtok(0,";"))) {  // advance by one nam{
    i++;
    if(strstr(blob,"<#>")) continue; // ignore some records
    if(strstr(blob,"</ver1>")) goto done; // end of record, tmp
    
    // printf("i=%d -->'%s' \n",i,blob);
    if(nFiber==0) {
      nFiber=atoi(blob);
      mDbFiber=new EEmcDbCrate[ nFiber];
      gMessMgr->Message("","D") <<"\n  EEDB mOptimizeFibers() map start: , nFiber="<<nFiber<<endm;
      icr=0;
      continue;
    }
    assert(icr<nFiber);
    mDbFiber[icr].setAll(blob);
    //  mDbFiber[icr].print();
    icr++;
  };

  }
  
  gMessMgr->Message("","E") <<"\n  EEDB mOptimizeFibers() map missing/wrong terminating key "<<endm;
  
 fatal:
    gMessMgr->Message("","E") <<"\n\n      EEDB EEMC  FiberConf  error,\n without it decoding of the _raw_ EEMC data from StEvent is not working\n(all data are dropped).\nHowever, this missing table is irreleveant for muDst analysis, JB\n\n"<<endm;
    return;

 done:
    gMessMgr->Message("","I") <<"\n  EEDB mOptimizeFibers() map found for nFiber="<<nFiber<<endm; 
   if(icr==nFiber)  return;

    gMessMgr->Message("","E") <<"\n  EEDB mOptimizeFibers() map nFiber missmatch is="<<icr<<" should be="<<nFiber<<endm;
goto fatal;
} 


//_________________________________________________________
//_________________________________________________________
//_________________________________________________________

Int_t StEEmcDbMaker::Make(){
  
  //  printf("\n\nMake :::::: %s\n\n\n",GetName());

  return kStOK;

}

//--------------------------------------------------
//--------------------------------------------------
const  EEmcDbItem*  
StEEmcDbMaker::getByIndex(int i){
  // Gets database entry by absolute index

  assert(i>=0);
  assert(i<EEindexMax);
  const  EEmcDbItem *x=byIndex+i;
  if(x->isEmpty()) return 0;
  return x;
}


//--------------------------------------------------
//--------------------------------------------------
const  EEmcDbItem*  StEEmcDbMaker::getByCrate(int crateID, int channel) {
  // crateID counts from 1, channel from 0 
  int type=0;
  int max=0;
  // printf("cr=%d ch=%d\n",crateID, channel);
  if(crateID>=MinTwCrateID && crateID<=MaxTwCrateID) {
    // Towers
    type =1;
    max=MaxTwCrateCh;
    if(channel>=max) return 0; // not all data blocks are used 

  } else if (crateID>=MinMapmtCrateID && crateID<=MaxMapmtCrateID ){
    //MAPMT
    type =2;
    max=MaxMapmtCrateCh;
  } else if (crateID>=17 && crateID<=46 ){ // only if working with ezTree
    return 0;  //BTOW
  }

  // printf("id=%d  type=%d ch=%d\n",crateID,type,channel);
  //printf(" p=%p \n",byCrate[crateID]);

  assert(type);
  assert( byCrate[crateID]);
  assert(channel>=0);
  assert(channel<max);
  return byCrate[crateID][channel];
}



//_________________________________________________________
//_________________________________________________________
//_________________________________________________________

template <class St_T, class T_st> void StEEmcDbMaker 
::getTable(TDataSet *eedb, int secID, TString tabName, TString mask,  T_st** outTab ){

  //  printf("\n\n%s ::TTT --> %s, size=%d\n",GetName(),tabName.Data(),sizeof(T_st));

  //  printf("%s ::TTT --> mask='%s' p=%p ss=%d\n",tabName.Data(),mask.Data(),*outTab,tabName.Contains(mask));

  if(!mask.IsNull() && !tabName.Contains(mask)) return ;
  char name[1000];
  if(secID<13) 
    sprintf(name,"sector%2.2d/%s",secID,tabName.Data());
  else
    sprintf(name,"misc/%s",tabName.Data());

  
  gMessMgr->Message("","D") <<"EEDB request="<< name <<endm;

  St_T *ds= (St_T *)eedb->Find(name);
  if(ds==0) {
    gMessMgr->Message("","W") <<"EEDB sector="<<secID<<" table='"<< name <<" not Found in DB, continue "<<endm;
    return ;
  }

  if(ds->GetNRows()!=1) {
    gMessMgr->Message("","W") <<"EEDB sector="<<secID<<" table='"<< name <<" no records, continue "<<endm; 
    return ;
  }
  
  T_st *tab=(T_st *) ds->GetArray();

  if(tab==0) {
    printf(" GetArray() failed\n");
    return  ;
  }

  *outTab=tab;
  gMessMgr->Message("","D") <<"  EEDB map="<< (*outTab)->comment<<"'"<<endm;

  nFound++;
  return ; // copy the whole s-struct to allow flavor change;
}

//__________________________________________________
//__________________________________________________
//__________________________________________________

const EEmcDbItem*  
StEEmcDbMaker::getByStrip0(int isec, int iuv, int istrip){
  //  printf("isec=%d iuv=%d istrip=%d \n",isec,iuv,istrip);
  assert(isec>=0 & isec<MaxSectors);
  assert(iuv>=0 && iuv<MaxSmdPlains);
  assert(istrip>=0 && istrip<MaxSmdStrips);
  return byStrip[isec][iuv][istrip];  
}



//__________________________________________________
//__________________________________________________
//__________________________________________________

const EEmcDbItem*  
StEEmcDbMaker::getTile(int sec, char sub, int eta, char type){
  char name[20];
  sprintf(name,"%2.2d%c%c%2.2d",sec,type,sub,eta);
  int key=EEname2Index(name);
  return byIndex+key;  
}

//__________________________________________________
//__________________________________________________
//__________________________________________________

const EEmcDbItem* 
StEEmcDbMaker::getStrip(int sec, char uv, int strip){
  char name[20];
  sprintf(name,"%2.2d%c%3.3d",sec,uv,strip);
  int key=EEname2Index(name);
  return byIndex+key;
}



//__________________________________________________
//__________________________________________________
//__________________________________________________

void StEEmcDbMaker::setAsciiDatabase( const Char_t *ascii )
{
  // This method allows the user to initialize database
  // values for each EEMC detector from an ascii file
  // rather from a database query.  InitRun() will be
  // blocked from reloading the database on subsequent
  // run numbers.
  
  mAsciiDbase = ascii;

  //-- Clear database values
  //$$$  clearItemArray();

  //-- Initialize timestamp to something obviously dumb.
  myTimeStampDay = 1; 

  //-- No provision (for now) for loading in partial info
  Int_t sec1  = 1;
  Int_t sec2  = 12;
  nFound      = 0;
  mfirstSecID = sec1;
  mlastSecID  = sec2;
  mNSector    = mlastSecID - mfirstSecID + 1;  

  //-- Open the specified file
  FILE *fd=fopen(ascii,"r"); 
  EEmcDbItem item;
  int nd=0;
  if(fd==0) goto crashIt;


  
  //-- Loop over all entries in the specified data file
  while (1) {

    //-- Read in each individual entry.  The return value is checked
    //-- for success or failure
    int ret = item.importAscii(fd);

    if(ret==0) break;
    if(ret==1) continue;
    if(ret <0) goto crashIt;

    //-- Index for the specified (named) item
    int key=EEname2Index(item.name);

    //-- Apply a sanity check on the item
    if(key!=item.key) {
      printf("WARN: name='%s' key=%d!=inpKey=%d, inpKey ignored\n",item.name,key,item.key);
      item.key=key;
    }
    assert(key>=0 && key<EEindexMax);


    int isec=item.sec-sec1;
    if (sec1>sec2 ) isec+=12;
    
    if( isec <0 ) continue;
    if( isec >=  mNSector) continue; 
    nd++;

    //-- Copy the database record read in from the file into
    //-- our local arrays.
    EEmcDbItem *x=&byIndex[key];
    if(!x->isEmpty()) goto crashIt;
    *x=item; // copy this DB record
    //    x->print();

    assert(byCrate[x->crate]); // ERROR: unsupported crate ID
    if(byCrate[x->crate][x->chan]) {
      printf("Fatal Error of eemc DB records: the same crate=%d / channel=%d entered twice for :\n",x->crate,x->chan);
      byCrate[x->crate][x->chan]->print(); // first time
      x->print(); // second time
      assert(1==2);
    }
    byCrate[x->crate][x->chan]=x;
    // assert(2==311);
  }
  
  //--
  //-- Initialize the byStrip lookup table
  //--
  for ( Int_t mySec = 1; mySec <= 12; mySec++ ) 
    for ( Char_t uv = 'U'; uv <= 'V'; uv++ ) 
      for ( Int_t myStrip = 1; myStrip <= 288; myStrip++ ) byStrip[mySec-1][uv-'U'][myStrip-1] = (EEmcDbItem*)getStrip(mySec,uv,myStrip);

      
  //--
  printf("setAsciiDataBase() done, found %d valid records\n",nd);

  return;

 crashIt:// any failure of reading of the data

  clearItemArray();
  //  assert(2==3);
  printf("EEmcDb - no/corrupted input file, all cleared, continue\n");

}

//________________________________________________________
//________________________________________________________
void  StEEmcDbMaker::changeGains(char *fname) {
  assert(nChGain+1< mxChGain);
  chGainL[nChGain]=fname;
  nChGain++;
}

//________________________________________________________
//________________________________________________________
void  StEEmcDbMaker::changeMask(char *fname) {
  assert(nChMask+1< mxChMask);
  chGainL[nChMask]=fname;
  nChMask++;
}

//________________________________________________________
//________________________________________________________
void  StEEmcDbMaker::changeGainsAction(const char *fname) {
    
  /* Replace gains only for channels already initialized from the DB 
     format : {name, gains anythingElse}
     lines starting with '#' are ignored
     empty lines are not permitted
  */

  gMessMgr->Message("","W") <<"  EEDB ::changeGains('"<<fname<<"')"<<endm;
  FILE *fd=fopen(fname,"r");
  int nd=0,nl=0;
  const int mx=1000;
  char buf[mx];

  if(fd==0) goto end;
  char cVal[100];
  float xVal;
    
  while(1) {
    char *ret=fgets(buf,mx,fd);
    if(ret==0) break;
    
    nl++;
    if(buf[0]=='#') continue;
    int n=sscanf(buf,"%s %f",cVal,&xVal);
    assert(n==2);
    int key=EEname2Index(cVal);
    EEmcDbItem *x=byIndex+key;
    // printf("%s %p\n",cVal,x);
    if(x->isEmpty()) continue;
    // replace only initialized channels
    if(xVal<=0)  printf("Warning ! buf=%s=\n",buf);
    //    assert(xVal>0);
     x->gain=xVal;
     nd++;
  }
  fclose(fd);
    
 end:
  gMessMgr->Message("","I") <<"  EEDB ::changeGains('"<<fname<<"') done inpLines="<<nl<<" nChanged="<<nd<<endm;
  return;
  
} 

//________________________________________________________
//________________________________________________________
void  StEEmcDbMaker::changeMaskAction(const char *fname) {

  /* Replace gains only for channels already initialized from the DB 
     format : {name, stat, fail enythingElse}
     lines starting with '#' are ignored
     empty lines are not permitted
  */

  gMessMgr->Message("","I") <<"  EEDB ::changeMask('"<<fname<<"') "<<endm;

  FILE *fd=fopen(fname,"r");

  int nd=0,nl=0;
  const int mx=1000;
  char buf[mx];

  if(fd==0) goto end;
  char cVal[100];
  int xStat,xFail ;

  while(1) {
    char *ret=fgets(buf,mx,fd);
    if(ret==0) break;

    nl++;
    if(buf[0]=='#') continue;
    int n=sscanf(buf,"%s %d %d",cVal,&xStat,&xFail);
    assert(n==3);
    int key=EEname2Index(cVal);
    EEmcDbItem *x=byIndex+key;
    if(x->isEmpty()) continue;
    assert(xStat>=0);
    assert(xFail>=0);
    x->stat=xStat;
    x->fail=xFail;
    nd++;
  }
  fclose(fd);

 end:
 gMessMgr->Message("","I") <<"  EEDB ::changeMask('"<<fname<<"') done inpLines="<<nl<<" nChanged="<<nd<<endm;
  return;

}



// $Log: StEEmcDbMaker.cxx,v $
// Revision 1.43  2005/01/24 05:08:26  balewski
// more get-methods
//
// Revision 1.42  2004/10/27 17:02:46  balewski
// move setKsig from Init() to constructor where it belongs
//
// Revision 1.41  2004/10/20 20:06:55  balewski
// only printouts
//
// Revision 1.40  2004/09/20 13:32:59  balewski
// to make Valgrind happy
//
// Revision 1.39  2004/08/09 20:17:12  balewski
// a bit more printout
//
// Revision 1.38  2004/08/07 02:46:51  perev
// WarnOff
//
// Revision 1.37  2004/07/27 22:00:19  balewski
// can overwrite gains & stat from DB
//
// Revision 1.36  2004/06/25 22:55:53  balewski
// now it survives missing fiberMap in DB , also gMessMgr is used
//
// Revision 1.35  2004/06/04 13:30:24  balewski
// use gMessMgr for most of output
//
// Revision 1.33  2004/05/20 16:40:14  balewski
// fix of strnlen --> strlen
//
// Revision 1.32  2004/05/14 20:55:34  balewski
// fix to process many runs, by Piotr
//
// Revision 1.31  2004/05/05 22:01:44  jwebb
// byStrip[] is now initialized when reading database from a file.
//
// Revision 1.30  2004/05/04 16:24:18  balewski
// ready for analysis of 62GeV AuAU production
//
// Revision 1.29  2004/04/28 20:38:10  jwebb
// Added StEEmcDbMaker::setAsciiDatabase().  Currently not working, since
// tube name missing for some towers, triggereing a "clear" of all EEmcDbItems.
//
// Revision 1.28  2004/04/12 16:19:51  balewski
// DB cleanup & update
//
// Revision 1.27  2004/04/09 18:38:10  balewski
// more access methods, not important for 63GeV production
//
// Revision 1.26  2004/04/08 16:28:06  balewski
// *** empty log message ***
//
// Revision 1.25  2004/04/04 06:10:37  balewski
// *** empty log message ***
//
// Revision 1.24  2004/03/30 04:44:57  balewski
// *** empty log message ***
//
// Revision 1.23  2004/03/28 04:09:08  balewski
// storage of EEMC raw data, not finished
//
// Revision 1.22  2004/03/19 21:31:53  balewski
// new EEMC data decoder
//
// Revision 1.21  2004/01/06 21:19:34  jwebb
// Added methods for accessing preshower, postshower and SMD info.
//
// Revision 1.20  2003/11/20 16:01:25  balewski
// towards run4
//
// Revision 1.19  2003/10/03 22:44:27  balewski
// fix '$' problem in db-entries name
//
// Revision 1.18  2003/09/11 05:49:17  perev
// ansi corrs
//
// Revision 1.17  2003/09/02 19:02:49  balewski
// fix for TMemeStat
//
// Revision 1.16  2003/08/27 03:26:45  balewski
// flavor option added:  myMk1->setPreferedFlavor("set-b","eemcPMTcal");
//
// Revision 1.15  2003/08/26 03:02:30  balewski
// fix of pix-stat and other
//
// Revision 1.14  2003/08/25 17:57:12  balewski
// use teplate to access DB-tables
//
// Revision 1.13  2003/08/22 20:52:20  balewski
// access to stat-table
//
// Revision 1.12  2003/08/02 01:02:17  perev
// change %d to %p int printf
//
// Revision 1.11  2003/07/18 18:31:46  perev
// test for nonexistance of XXXReader added
//
// Revision 1.10  2003/04/27 23:08:13  balewski
// clean up of daq-reader
//
// Revision 1.9  2003/04/25 14:42:00  jeromel
// Minor change in messaging
//
// Revision 1.8  2003/04/16 20:33:51  balewski
// small fixes in eemc daq reader
//
// Revision 1.7  2003/04/02 20:42:23  balewski
// tower-->tube mapping
//
// Revision 1.6  2003/03/26 21:28:02  balewski
// fix
//
// Revision 1.5  2003/03/26 15:26:23  balewski
// add print()
//
// Revision 1.4  2003/03/07 15:35:44  balewski
// towards EEMC daq reader
//
// Revision 1.3  2003/02/18 22:01:40  balewski
// fixes
//
// Revision 1.2  2003/02/18 19:55:53  balewski
// add pedestals
//
// Revision 1.1  2003/01/28 23:18:34  balewski
// start
//
// Revision 1.5  2003/01/06 17:09:21  balewski
// DB-fix
//
// Revision 1.4  2003/01/03 23:37:56  balewski
// move to robinson
//
// Revision 1.3  2002/12/05 14:22:24  balewski
// cleanup, time stamp fixed
//
// Revision 1.2  2002/12/04 13:39:04  balewski
// remove dependency on dbase/
//
// Revision 1.1  2002/11/30 20:01:26  balewski
// start DB interface for EEMC RELATED ROUTINES
//


#if 0  
  // test of flavor
  {
    St_db_Maker* mydb = (St_db_Maker*)GetMaker("StarDb");
    mydb->SetDateTime(20030814,0); // set ~day & ~hour by hand
    printf("\nJB:test of flavor, use time stamp=");
    (mydb->GetDateTime()).Print();

    TDataSet *eedb1=GetDataBase("TestScheme/emc");
    assert(eedb1);

    St_eemcDbPMTcal *ds= (St_eemcDbPMTcal *)eedb1->Find("sector04/eemcPMTcal");
    assert(ds);
    eemcDbPMTcal_st *tab1=(eemcDbPMTcal_st *) ds->GetArray();
    assert(tab1);
    
    printf("1) p=%p tab1->comment=%s\n",tab1,tab1->comment);    

    printf("\n change flavor\n");
    SetFlavor("setb","eemcPMTcal");
    TDataSet *eedb2=GetDataBase("TestScheme/emc");
    assert(eedb2);

    ds= (St_eemcDbPMTcal *)eedb2->Find("sector04/eemcPMTcal");
    assert(ds);
    
    eemcDbPMTcal_st *tab=(eemcDbPMTcal_st *) ds->GetArray();
    assert(tab);

    printf("2) p=%p tab->comment=%s\n",tab,tab->comment);
    printf("1)\' p=%p tab1->comment=%s\n",tab1,tab1->comment);    
 

  }
  // end of flavor test
  assert(1==5);

#endif
