// Hey Emacs this is really -*-c++-*- ! 
// \class  EEmcMCData
// \author Piotr A. Zolnierczuk             
// \date   Aug 26, 2002
#ifndef EEmcMCData_h
#define EEmcMCData_h
/*********************************************************************
 * $Id: EEmcMCData.h,v 1.2 2003/02/20 20:13:20 balewski Exp $
 *********************************************************************
 * Description:
 * STAR Endcap Electromagnetic Calorimeter Monte Carlo Data
 *********************************************************************
 * $Log: EEmcMCData.h,v $
 * Revision 1.2  2003/02/20 20:13:20  balewski
 * fixxy
 * xy
 *
 * Revision 1.1  2003/02/20 05:14:07  balewski
 * reorganization
 *
 * Revision 1.1  2003/01/28 23:16:07  balewski
 * start
 *
 * Revision 1.9  2002/11/13 21:53:52  zolnie
 * restored "private" DetectorID in EEmcMCData
 *
 * Revision 1.8  2002/11/11 21:22:48  balewski
 * EEMC added to StEvent
 *
 * Revision 1.7  2002/10/03 15:52:25  zolnie
 * updates reflecting changes in *.g files
 *
 * Revision 1.6  2002/10/01 06:03:15  balewski
 * added smd & pre2 to TTree, tof removed
 *
 * Revision 1.5  2002/09/30 21:58:27  zolnie
 * do we understand Oleg? (the depth problem)
 *
 * Revision 1.4  2002/09/30 20:15:55  zolnie
 * Oleg's geometry updates
 *
 * Revision 1.3  2002/09/24 22:47:35  zolnie
 * major rewrite: SMD incorporated, use constants rather hard coded numbers
 * 	introducing exceptions (rather assert)
 *
 * Revision 1.2  2002/09/20 15:49:05  balewski
 * add event ID
 *
 * Revision 1.1.1.1  2002/09/19 18:58:54  zolnie
 * Imported sources
 *
 *********************************************************************/
#include "TObject.h"
#include "EEmcDefs.h"
//#include "EEmcException.h"

class   StMaker;

class   St_g2t_emc_hit;
class   St_g2t_event;

class   EEeventDst;
//class   EEmcException1;

const   Float_t kEEmcDefaultEnergyThreshold = 0.0005; // 0.5 MeV
const   Int_t   kEEmcDefaultMCHitSize       = 0x1000; // 4k hitow


struct EEmcMCHit {
  UChar_t  detector;  // endcap detector part  (prs,tower,post,smdu,smdv)
  UChar_t  sector;    // endcap phi sector 1:12
  union {
    struct EEmcMCHitTower { 
      UChar_t  ssec;  // endcap subsector  1:5 (A-E)
      UChar_t  eta;   // endcap eta        1:12
    } tower;
    UShort_t strip;   // smd's  strip numbers
  };
  Float_t  de;        // energy loss in the element (GeV)
};

  

class   EEmcMCData : public TObject {
public:
  enum MCDepth {          // FIXME LATER: depths encoded in g2t data
    kUnknownDepth    = 0,
    kPreShower1Depth = 1,
    kPreShower2Depth = 2,
    kTower1Depth     = 3,
    kTower2Depth     = 4,
    kPostShowerDepth = 5
  };


  enum MCDetectorId {
    kEEmcMCUnknownId    = 0,
    kEEmcMCTowerId      = 1,
    kEEmcMCPreShower1Id = 2,
    kEEmcMCPreShower2Id = 3,
    kEEmcMCSmdUStripId  = 4,
    kEEmcMCSmdVStripId  = 5,
    kEEmcMCPostShowerId = 6
  };


  EEmcMCData();                              // default constructor
  EEmcMCData( EEmcMCData& );                 // copy constructor
  virtual ~EEmcMCData();                     // the destructor

  Int_t      readEventFromChain(StMaker *mk); // reads g2t event from chain
  
  Int_t      getSize()            { return mSize;    };
  Int_t      getLastHit()         { return mLastHit; }; 
  Int_t      getEventID()         { return mEventID; }; 

  Float_t    getEnergyThreshold() { return mEthr; };
  void       setEnergyThreshold(Float_t e) { mEthr = e; };

  Int_t      getHitArray(EEmcMCHit *h, Int_t size);
  Int_t      setHitArray(EEmcMCHit *h, Int_t size);

  void       print();                                      // diagnostic print

  // obsolete functions
  Int_t    read     (void *d, int s);   // reads  in  s bytes of hits 
  Int_t    write    (void *d, int s);   // writes out s bytes of hits 
  Int_t    write    (EEeventDst *);        // export hits to TTree
  
protected:
  Int_t             mEventID;  // WARN: not added in to binary I/O, JB
  struct EEmcMCHit *mHit ;     // array to hold hit data
  Int_t             mSize;     // size of the above array
  Int_t             mLastHit;  // last hit index
  Float_t           mEthr;     // energy threshold
  
private:
  //Float_t          *mDE; // temporary array
  Int_t             expandMemory();
  
  ClassDef(EEmcMCData,2) // Endcap Emc event
};

#endif
