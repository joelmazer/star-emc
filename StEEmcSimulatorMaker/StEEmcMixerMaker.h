#ifndef STAR_StEEmcMixerMaker
#define STAR_StEEmcMixerMaker

/*!\class StEEmcMixerMaker
\author Wei-Ming Zhang, KSU 3/16/2005 
Based on StEmcMixerMaker of Alexandre A. P. Suaide

  This is the basic embedding scheme for EEMC. This maker gets the output
of two StEmcCollections objects:
- primary built from StEmcRawData, data from daq 
- secondary build in fly & owned by StEEmcFastMaker, data from fzd/geant.root 

and merges ONLY ADC hits for all EEMC pixels. 
Hits of the second event are added to the first one.

  The maker has full conectivity to the EEMC database through StEEmcDbMaker, 
  so it will calculate ADC from geant DE using actual gains for hits from secondary event.
  If a gain does not exist this pixel from secondary file is dropped.  
  
  It outputs an embedded StEmcCollection object which will be fetched
by StEmcRawMaker to build embedded data of event.root and MuDst.root.
*/

#include "StMaker.h"

class StEmcCollection;
class StEmcRawHit;
class StEEmcDbMaker;
class EEmcDbItem;
                                               
class StEEmcMixerMaker : public  StMaker {
  private:
    Bool_t        mClear;
    StEEmcDbMaker *mEEDb;

//  Embedded object which can be used by other makers.  
    StEmcCollection *mMixerEmcCollection;
    
    void    mergeADCs(StEmcCollection*A,StEmcCollection*B);   


  public: 
    StEEmcMixerMaker(const char *name="EEmcMixer");
    virtual       ~StEEmcMixerMaker();
    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
    StEmcCollection *getMixerEmcCollection() { return mMixerEmcCollection; }


   virtual const char *GetCVS() const {static const char cvs[]="Tag $Name:  $ $Id: StEEmcMixerMaker.h,v 1.1 2006/12/12 20:29:14 balewski Exp $ built "__DATE__" "__TIME__ ; return cvs;}

  ClassDef(StEEmcMixerMaker,0) 
};
#endif

///////////////////////////////////////////////////////////////////////////
//
// $Id: StEEmcMixerMaker.h,v 1.1 2006/12/12 20:29:14 balewski Exp $
// $Log: StEEmcMixerMaker.h,v $
// Revision 1.1  2006/12/12 20:29:14  balewski
// added hooks for Endcap embedding
//
// Revision 1.1.1.1  2005/05/31 18:53:25  wzhang
// First version
//
//
///////////////////////////////////////////////////////////////////////////

