#ifndef STAR_StEmcSimulatorMaker
#define STAR_StEmcSimulatorMaker
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// StEmcSimulatorMaker
//
// It replaced St_ems_Maker. This is the clear C++ code.
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include "StEmcUtil/others/emcInternalDef.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "tables/St_emc_hits_Table.h"
#include "StMessMgr.h"

#define EMCSIMNEW

class StMcEmcHitCollection;
class StEmcCollection;
class StEmcVirtualSimulator;
class St_controlEmcSimulatorMaker;
class St_controlEmcPmtSimulator;
class controlEmcSimulatorMaker_st;
class controlEmcPmtSimulator_st;
class St_emcStatus;
class emcStatus_st;
class St_smdStatus;
class smdStatus_st;
class St_g2t_emc_hit;

class StEmcSimulatorMaker : public StMaker
{
private:

    UInt_t                  mBEMC;            // Switch for BEMC; 0 => off; >0 => on
    UInt_t                  mEEMC;            // Switch for EEMC; 0 => off; >0 => on
    UInt_t                  mHistControl;     // Do histogramms (1) or no (0)
    UInt_t                  mDB;              // =0, no DB; >=1, using DB

    StEmcGeom*              mGeom[MAXDET];    // Geometry

    Bool_t                  mCompare;
    TCanvas*                mC1;              //

    Float_t                 mGain[MAXDET][18000];

    TDataSet*               geaIn;
    TDataSet*               ems;
    St_g2t_emc_hit*         g2t_emc_hit;
    St_g2t_emc_hit*         g2t_smd_hit;

    St_controlEmcSimulatorMaker* controlMaker;
    controlEmcSimulatorMaker_st* controlTable;

    St_controlEmcPmtSimulator* pmtSimulator;
    controlEmcPmtSimulator_st* pmtTable;

    TDataSet*               DB;
    TDataSet*               status;
    TDataSet*               ped;
    TDataSet*               calib;
    St_emcStatus*           statusEmc;   // status for BEMC  or BPRS
    emcStatus_st*           statusEmcRec;
    St_smdStatus*           statusSmd;   // status for BSMDE or BSMDP
    smdStatus_st*           statusSmdRec;

    StMcEmcHitCollection*   mEmcMcHits[MAXDET];  // For convinience
    St_emc_hits*            mEmcRawHits[MAXDET]; // For convinience
    StEmcCollection*        mEmcCollection;      // As in StEvent

    StEmcVirtualSimulator*  mSimulator[MAXDET];  //

    TH2F*                   m_nhit;           //!
    TH2F*                   m_etot;           //!
    TH2F*                   m_hits[MAXDET];   //!
    TH2F*                   m_energy[MAXDET]; //!
    TH1F*                   m_adc[MAXDET];    //!
    TH1F*                   mEnergySum[MAXDET];    //!

    TH1F*                   mhModule[MAXDET];    //! For testing only
    TH1F*                   mhSub[MAXDET];
    TH1F*                   mhDiffNumHits[4];    //!
    TH1F*                   mhDiffDe[4];         //!

    Bool_t                  mEmbed;

    void                    addBemcAndBprsHit(Int_t,Int_t,Int_t,Int_t,Float_t);
    void                    addBsmdeAndBsmdpHit(Int_t,Int_t,Int_t,Int_t,Float_t);

public:
    StEmcSimulatorMaker(const char *name="EmcSimulator");
    virtual                 ~StEmcSimulatorMaker();
    virtual Int_t           Init();
    virtual Int_t           InitRun(Int_t);
    virtual Int_t           Make();
    virtual void            Clear(const char *);
    Int_t                   fillStEvent();
    void                    saveRunco();

    Int_t                   makeBemc();
    Int_t                   makeBemcAndBprsMcHits();
    Int_t                   makeBsmdeAndBsmdpMcHits();
    Int_t                   makeAllRawHitsForBemc();

    Int_t                   makeEemc();

    void                    bookHistograms(const Int_t);
    void                    makeHistograms(const Int_t); //! must be changed

    UInt_t                  getBEMC()
    {
        return mBEMC;
    }
    UInt_t                  getEEMC()
    {
        return mEEMC;
    }
    UInt_t                  getHistControl()
    {
        return mHistControl;
    }
    StMcEmcHitCollection*   getEmcMcHits(Int_t det)
    {
        return mEmcMcHits[det-1];
    }
    StMcEmcHitCollection*   getBemcMcHits()
    {
        return getEmcMcHits(BEMC);
    }
    StMcEmcHitCollection*   getBprsMcHits()
    {
        return getEmcMcHits(BPRS);
    }
    StMcEmcHitCollection*   getBsmdeMcHits()
    {
        return getEmcMcHits(BSMDE);
    }
    StMcEmcHitCollection*   getBsmdpMcHits()
    {
        return getEmcMcHits(BSMDP);
    }

    void                    clearStEventStaf();
    void                    Browse(TBrowser* b); // StEvent staf will be visible in browser
    St_controlEmcSimulatorMaker* getControlSimulator();
    St_controlEmcPmtSimulator*   getControlPmtSimulator();
    StEmcCollection*        getEmcCollection();

    void                    pictureAllDetectors(Int_t print=0);                          // *MENU*
    void                    pictureForDetector(Int_t det, Int_t logy=1, Int_t print=0);  // *MENU*
    void                    pictureCompareDe(Int_t print=0);                             // *MENU*
    void                    printSimulator(Int_t det=0);                                 // *MENU*
    void                    printStatusTable(Int_t det=1, Int_t hist=0);                 // *MENU*
    TDataSet*               getStatus(const Int_t ind, TDataSet* statusDB);

    void                    compareOldSimulator();

    void                    printmBEMC();
    void                    setBEMC(UInt_t  key)
    {
        mBEMC = key;
        if (Debug())
            printmBEMC();
    }
    void                    setPrint(Bool_t a)
    {
		LOG_INFO << "::setPrint() is obsolete.  Use logger config file to set verbosity instead." << endm;
	}///< Obsolete function; users can control messages with logger config file.
    void                    setHistControl(UInt_t key)
    {
        mHistControl = key;
    }
    virtual const char*     GetCVS() const
    {
        static const char cvs[]="Tag $Name:  $ $Id: StEmcSimulatorMaker.h,v 1.16 2007/01/22 19:13:40 kocolosk Exp $ built "__DATE__" "__TIME__ ;
        return cvs;
    }

    ClassDef(StEmcSimulatorMaker,0)  // Simulation maker for BEMC and EEMC
};

#endif
//////////////////////////////////////////////////////////////////////////////////
//
// $Id: StEmcSimulatorMaker.h,v 1.16 2007/01/22 19:13:40 kocolosk Exp $
// $Log: StEmcSimulatorMaker.h,v $
// Revision 1.16  2007/01/22 19:13:40  kocolosk
// use STAR logger for all output
//
// Revision 1.15  2005/03/21 21:36:39  suaide
// fixed problem with chain
//
// Revision 1.14  2004/08/09 19:43:28  suaide
// moved global variables to private members and
// made small modifications to run in embedding mode
//
// Revision 1.13  2004/08/06 13:24:48  suaide
// New features added and fixed some bugs in the database
//
// Revision 1.12  2004/04/08 21:35:45  perev
// Leak off
//
// Revision 1.11  2003/09/23 15:19:55  suaide
// fixed bugs and modifications for embedding
//
// Revision 1.10  2003/09/10 19:47:12  perev
// ansi corrs
//
// Revision 1.9  2003/01/23 03:09:02  jeromel
// Include modif
//
// Revision 1.8  2002/09/10 16:51:32  pavlinov
// Discard line with mDbMaker->SetDateTime
//
// Revision 1.7  2002/06/04 16:09:37  pavlinov
// added option with DB(pedestal ans calibration  coefficients
//
// Revision 1.6  2002/06/03 23:35:11  pavlinov
// Last correction without DB for ped and calib. coeff.
//
// Revision 1.5  2001/09/22 00:29:47  pavlinov
// No public constructor for StEmcGeom
//
// Revision 1.4  2001/03/22 22:04:45  pavlinov
// Clean up for mdc4
//
// Revision 1.3  2001/02/03 00:00:01  pavlinov
// New function Browse() and cleanup for new version of BFC
//
// Revision 1.2  2000/10/28 00:33:45  pavlinov
// added methods getEmcCollectin()
//
// Revision 1.1  2000/10/23 22:53:14  pavlinov
// First working C++ version
//
//
//////////////////////////////////////////////////////////////////////////////////
