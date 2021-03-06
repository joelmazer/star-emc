//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Macro for running chain with different inputs                        //
// owner:  Yuri Fisyak                                                  //
//                                                                      //
// $Id: bfc.C,v 1.3 2011/04/08 22:35:02 balewski Exp $
//////////////////////////////////////////////////////////////////////////
class StBFChain;        
class StMessMgr;
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Stiostream.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TApplication.h"
#include "TInterpreter.h"
#include "StBFChain.h"
#include "StMessMgr.h"
#else
#endif
#define UseLogger
StBFChain    *chain=0; 
//_____________________________________________________________________
//_________________ Prototypes _______________________________________________
void Usage();
void Load(const Char_t *options="");
TString defChain("y2005h,Test.default.ITTF");
void bfc(Int_t First, Int_t Last,const Char_t *Chain = defChain + ",Display",
	 const Char_t *infile=0, const Char_t *outfile=0, const Char_t *TreeFile=0);
//	 const Char_t *Chain="gstar,y2005h,MakeEvent,trs,sss,svt,ssd,fss,bbcSim,emcY2,tpcI,fcf,ftpc,SvtCL,svtDb,ssdDb,svtIT,ssdIT,ITTF,genvtx,Idst,event,analysis,EventQA,tags,Tree,EvOut,McEvOut,GeantOut,IdTruth,miniMcMk,StarMagField,FieldOn,McAna,Display",//,,NoSimuDb, display, //McQa, 
void bfc(Int_t Last, const Char_t *Chain = defChain,
	 const Char_t *infile=0, const Char_t *outfile=0, const Char_t *TreeFile=0);
	 //	 const Char_t *Chain="gstar,y2005h,tpcDb,trs,tpc,Physics,Cdst,Kalman,tags,Tree,EvOut,McEvOut,IdTruth,miniMcMk,StarMagField,FieldOn,McAna", // McQA
//_____________________________________________________________________
void Load(const Char_t *options){
  cout << "Load system libraries" << endl;
  if ( gClassTable->GetID("TGiant3") < 0) { // ! root4star
    cout << endl << "Load ";
    if (!TString(options).Contains("nodefault",TString::kIgnoreCase) || 
	 TString(options).Contains("pgf77",TString::kIgnoreCase)) {
      const Char_t *pgf77 = "libpgf77VMC";
      if (gSystem->DynamicPathName(pgf77,kTRUE) ) {
	gSystem->Load(pgf77); cout << " " << pgf77 << " + ";
      }
    }
    if (!TString(options).Contains("nodefault",TString::kIgnoreCase) || 
	TString(options).Contains("cern",TString::kIgnoreCase)) {
      gSystem->Load("libminicern"); cout << "libminicern";
    }
    if (!TString(options).Contains("nodefault",TString::kIgnoreCase) || 
	TString(options).Contains("mysql",TString::kIgnoreCase)) {
      Char_t *mysql = "libmysqlclient";
      //
      // ATTENTION: The below will FAIL for 64 bits systems (JL 2009/10/22)
      //
      Char_t *libs[]  = {"", "/usr/mysql/lib/","/usr/lib/", 0}; // "$ROOTSYS/mysql-4.1.20/lib/",
      //Char_t *libs[]  = {"/usr/lib/", 0};
      Int_t i = 0;
      while ((libs[i])) {
	TString lib(libs[i]);
	lib += mysql;
	lib = gSystem->ExpandPathName(lib.Data());
	if (gSystem->DynamicPathName(lib,kTRUE)) {
	  gSystem->Load(lib.Data()); cout << " + " << lib.Data() << endl;
	  break;
	}
	i++;
      }
    }
  }
  //  if (gClassTable->GetID("TMatrix") < 0) gSystem->Load("StarRoot");// moved to rootlogon.C  StMemStat::PrintMem("load StarRoot");
#ifdef UseLogger
  // Look up for the logger option
  Bool_t needLogger  = kFALSE;
  if (!TString(options).Contains("-logger",TString::kIgnoreCase)) {
    needLogger = gSystem->Load("liblog4cxx.so") <= 0;              //  StMemStat::PrintMem("load log4cxx");
  }
#endif
  gSystem->Load("libSt_base");                                        //  StMemStat::PrintMem("load St_base");
#ifdef UseLogger
  if (needLogger) {
    gSystem->Load("libStStarLogger.so");
    gROOT->ProcessLine("StLoggerManager::StarLoggerInit();");      //  StMemStat::PrintMem("load StStarLogger");
  }
#endif
  gSystem->Load("libHtml");
  gSystem->Load("libStChain");                                        //  StMemStat::PrintMem("load StChain");
  gSystem->Load("libStUtilities");                                    //  StMemStat::PrintMem("load StUtilities");
  gSystem->Load("libStBFChain");                                      //  StMemStat::PrintMem("load StBFChain");
  gSystem->Load("libStChallenger");                                   //  StMemStat::PrintMem("load StChallenger");
}
//_____________________________________________________________________
void bfc(Int_t First, Int_t Last,
	 const Char_t *Chain,
	 const Char_t *infile,
	 const Char_t *outfile,
	 const Char_t *TreeFile)
{ // Chain variable define the chain configuration 
  // All symbols are significant (regardless of case)
  // "-" sign before requiest means that this option is disallowed
  // Chain = "gstar" run GEANT on flight with 10 muons in range |eta| < 1 amd pT = 1GeV/c (default)
  // Dynamically link some shared libs
  if (gClassTable->GetID("StBFChain") < 0) Load(Chain);
  chain = new StBFChain(); cout << "Create chain " << chain->GetName() << endl;
  TString tChain(Chain);
  chain->cd();
  chain->SetDebug(1);
  if (Last < -3) return;
  chain->SetFlags(Chain);
  if (tChain == "" || ! tChain.CompareTo("ittf",TString::kIgnoreCase)) Usage();
  chain->Set_IO_Files(infile,outfile);
  if (TreeFile) chain->SetTFile(new TFile(TreeFile,"RECREATE"));
  gMessMgr->QAInfo() << Form("Process [First=%6i/Last=%6i/Total=%6i] Events",First,Last,Last-First+1) << endm;
  if (Last < -2) return;
  if (chain->Load() > kStOk) {
    gMessMgr->Error() << "Problems with loading of shared library(ies)" << endm;
    gSystem->Exit(1);
  }
  if (Last < -1) return;
  if (chain->Instantiate() > kStOk)  { 
    gMessMgr->Error() << "Problems with instantiation of Maker(s)" << endm;
    gSystem->Exit(1);
  }
  if (Last < 0) return;
  StMaker *dbMk = chain->GetMaker("db");
  if (dbMk) dbMk->SetDebug(1);
#if 0
  // Insert your maker before "tpc_hits"
  Char_t *myMaker = "St_TLA_Maker";
  if (gClassTable->GetID(myMaker) < 0) {
	  gSystem->Load(myMaker);//  TString ts("load "; ts+=myMaker; StMemStat::PrintMem(ts.Data());
  }
  StMaker *myMk = chain->GetMaker(myMaker);
  if (myMk) delete myMk;
  myMk = chain->New(myMaker,"before");
  if (myMk) {
    Char_t *before = "tpc_hits";
    StMaker *tclmk = chain->GetMaker(before);
    if (tclmk) chain->AddBefore(before,myMk);
  }
  // Insert your maker after "tpc_hits"
  myMk = chain->New(myMaker,"after");
  if (myMk) {
    Char_t *after = "tpc_hits";
    StMaker *tclmk = chain->GetMaker(after);
    if (tclmk) chain->AddAfter(after,myMk);
  }
  // this block is meant as an example ONLY
  // The default values are set in StRoot/StPass0CalibMaker/ StTpcT0Maker 
  // constructor and are suitable for production. You can change it here
  // for test purposes.
  if (chain->GetOption("TpcT0")) {
    StTpcT0Maker *t0mk = (StTpcT0Maker *) chain->GetMaker("TpcT0");
    if (t0mk) t0mk->SetDesiredEntries(10);
  }
#endif
  {
    TDatime t;
    gMessMgr->QAInfo() << Form("Run is started at Date/Time %i/%i",t.GetDate(),t.GetTime()) << endm;
  }
  gMessMgr->QAInfo() << Form("Run on %s in %s",gSystem->HostName(),gSystem->WorkingDirectory()) << endm;
  gMessMgr->QAInfo() << Form("with %s", chain->GetCVS()) << endm;
  // Init the chain and all its makers
  TAttr::SetDebug(0);
  chain->SetAttr(".Privilege",0,"*"                ); 	  //All  makers are NOT priviliged
  chain->SetAttr(".Privilege",1,"StIOInterFace::*" ); 	  //All IO makers are priviliged
  chain->SetAttr(".Privilege",1,"St_geant_Maker::*"); 	  //It is also IO maker
  chain->SetAttr(".Privilege",1,"StTpcDbMaker::*"); 	  //It is also TpcDb maker to catch trips
  chain->SetAttr(".Privilege",1,"*::tpc_hits"); //May be allowed to act upon excessive events
  chain->SetAttr(".Privilege",1,"*::tpx_hits"); //May be allowed to act upon excessive events
  chain->SetAttr(".Privilege",1,"StTpcHitMover::*"); //May be allowed to act upon corrupt events
  chain->SetAttr(".Privilege",1,"*::tpcChain"); //May pass on messages from sub-makers
  chain->SetAttr(".Privilege",1,"StTriggerDataMaker::*"); //TriggerData could reject event based on corrupt triggers

  // add FGT stuff - jan
  TString &fnameOut1=chain->GetFileOut();
  cout<<"AAA"<<fnameOut1<<endl;
  TString fnameOut2=fnameOut1;
  fnameOut2.ReplaceAll(".root",".fgt");
  //fnameOut2+=".fgt";
  cout<<"FGT output core="<<fnameOut2<<endl;

  gSystem->Load("libStEEmcPoolFgtSandbox1");
  StFgtSlowSimuMaker  *myMkSM=new StFgtSlowSimuMaker;
  StFgtClustFindMaker *myMkCL=new StFgtClustFindMaker;
  StFgtClustEvalMaker *myMkEV=new StFgtClustEvalMaker;
  // ... position FGT maker in the chain
  Char_t *after = "geant"; StMaker *xMk = chain->GetMaker(after);  assert(xMk);
  // note , the lines should be in reversed order
  chain->AddAfter(after,myMkEV);  //position #3
  chain->AddAfter(after,myMkCL); //position #2
  chain->AddAfter(after,myMkSM); //position #1

  //...... auxiliary initialization
  HList=new  TObjArray;
  myMkSM->setHList(HList);  myMkCL->setHList(HList);  myMkEV->setHList(HList);
  myMkSM->initFrankModel("StRoot/StEEmcPool/FgtSandbox1/macros/BichselELossProbHighBG.dat");
  // myMkSM->forcePerpTracks(); //only  for testing,all tracks will go perp through the GEM gas volume
//  myMkSM->useOnlyDisk(1); //only  for testing, disk # 1...6
  myMkSM->useOnlyDisk(0); //All disks
  myMkSM->setRadStripGain(1.00,0.03); //relative gain variation: mean=1 & sigma=0.03  
  myMkSM->setPhiStripGain(1.00,0.03); //relative gain variation: mean=1.5 & sigma=0.03  
  myMkSM->setStripThresh(1.); // a.u. 
  myMkSM->setHexGemLatice(0.014, 60.); // (cm) pitch=0.014, (deg) 1st axis angle/deg
  myMkSM->setTransDiffusion(0.017); // (cm per 1 cm of path), default=0.017 cm/cm

  myMkCL->setSeedStripThres(30.); // a.u. 
  myMkCL->setStripNoiseSigma(12.); // a.u. , default=12=10%

  // chain is ready
  chain->ls(3);


  Int_t iInit = chain->Init();




  if (iInit >=  kStEOF) {chain->FatalErr(iInit,"on init"); return;}
  if (Last == 0) return;
  StEvtHddr *hd = (StEvtHddr*)chain->GetDataSet("EvtHddr");
  if (hd) hd->SetRunNumber(-2); // to be sure that InitRun calls at least once
    // skip if any
  chain->EventLoop(First,Last,0);
  gMessMgr->QAInfo() << "Run completed " << endm;
  gSystem->Exec("date");


  //FGT end....
  {
    TDatime t;
    gMessMgr->QAInfo() << Form("Run is finished at Date/Time %i/%i",
                               t.GetDate(),t.GetTime()) << endm;
  }
  myMkSM->saveHisto(fnameOut2);// one is enough, all histos are merged




}
//_____________________________________________________________________
void bfc(Int_t Last, 
	 const Char_t *Chain,
	 const Char_t *infile, 
	 const Char_t *outfile, 
	 const Char_t *TreeFile)
{
  bfc(1,Last,Chain,infile,outfile,TreeFile);
}
//____________________________________________________________
void Usage() {
  printf ("============= \t U S A G E =============\n");
  printf ("bfc(Int_t First,Int_t Last,const Char_t *Chain,const Char_t *infile,const Char_t *outfile,const Char_t *TreeFile)\n");
  printf ("bfc(Int_t Last,const Char_t *Chain,const Char_t *infile,const Char_t *outfile,const Char_t *TreeFile)\n");
  printf ("bfc(const Char_t *ChainShort,Int_t Last,const Char_t *infile,const Char_t *outfile)\n");
  printf ("where\n");
  printf (" First     \t- First event to process\t(Default = 1)\n");
  printf (" Last      \t- Last  event to process\t(Default = 1)\n");
  printf (" Chain     \t- Chain specification   \t(without First &  Last: Default is \"\" which gives this message)\n");
  printf ("           \t                        \t with    First || Last: Default is \"gstar tfs\")\n");
  printf (" infile    \t- Name of Input file    \t(Default = 0, i.e. use preset file names depending on Chain)\n"); 
  printf (" outfile   \t- Name of Output file   \t(Default = 0, i.e. define Output file name from Input one)\n");
  printf (" outfile   \t- Name of Tree File     \t(Default = 0, i.e. define Output file name from Input one (tags TNtuple))\n");
  printf (" ChainShort\t- Short cut for chain   \t(Default = \"\" -> print out of this message)\n");
  gSystem->Exit(1);
}
//_____________________________________________________________________
void bfc(const Char_t *Chain="ittf") {
  bfc(-2,Chain);
}
