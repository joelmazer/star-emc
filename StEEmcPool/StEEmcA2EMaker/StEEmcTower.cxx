/**
 * \class StEEmcTower
 * \brief Base class for representing tower, preshower and postshower elements
 *
 * This class inherits the StEEmcElement properies.  In addition, 
 * each tower object stores its index in a 720 element array.  It
 * knows its sector, subsector and etabin.  It also has pointers 
 * to its 8 nearest neighbors on the endcap.
 *
 * \author Jason C. Webb
 * $Date: 2010/08/26 22:49:21 $
 * $Revision: 1.3 $
 *
 */


#include "StEEmcTower.h"
#include <iostream>
#include "StEEmcUtil/EEmcGeom/EEmcGeomSimple.h"

ClassImp(StEEmcTower);

// ----------------------------------------------------------------------------
StEEmcTower::StEEmcTower()
    : StEEmcElement()
{
  mIndex=-1;
  mLayer=-1;
  mSector=-1;
  mPhibin=-1;
  mSubsector=-1;
  mEtabin=-1;
}

// ----------------------------------------------------------------------------
StEEmcTower::StEEmcTower( const StEEmcTower &other )
    : StEEmcElement(other)
{
  mIndex=other.mIndex;
  mLayer=other.mLayer;
  mSector=other.mSector;
  mPhibin=other.mPhibin;
  mSubsector=other.mSubsector;
  mEtabin=other.mEtabin;
  mET=other.mET;
  mNeighbors=other.mNeighbors;

  mRaw=other.mRaw;
  mAdc=other.mAdc;
  mEnergy=other.mEnergy;
  mStat=other.mStat;
  mFail=other.mFail;
  mName=other.mName;
  mstRawHit=other.mstRawHit;
}

// ----------------------------------------------------------------------------
void StEEmcTower::print() const
{

    std::cout << "--------------------------- StEEmcTower::print()" << std::endl; 
    std::cout << "name   = " << mName << std::endl;
    std::cout << "index  = " << mIndex << std::endl; 
    std::cout << "energy = " << energy() << std::endl;
    std::cout << "raw    = " << raw() << std::endl;
    std::cout << "adc    = " << adc() << std::endl;

}

void StEEmcTower::printLine() const
{
  std::cout << mName << " ADC=" << raw() << " " << adc() << " E=" << energy();
}

// ----------------------------------------------------------------------------
void StEEmcTower::index( Int_t i )
{

  mIndex  = i;
  mPhibin = i / 12;
  mEtabin = i % 12;
  mSector    = mPhibin / 5;
  mSubsector = mPhibin % 5;

  const Char_t *secs[]={"01","02","03","04","05","06","07","08","09","10","11","12"};
  const Char_t *subs[]={"A","B","C","D","E"};
  const Char_t *lays[]={"T","P","Q","R","S"};

  mName="";
  mName += secs[mSector];
  mName += lays[mLayer];
  mName += subs[mSubsector];
  mName += secs[mEtabin];
  
}

// ----------------------------------------------------------------------------
Bool_t StEEmcTower::isNeighbor( const StEEmcTower &t ) const
{

  /// passed a copy of this tower, return true
  if ( t.index() == index() ) return true;

  /// outside of eta range
  if ( TMath::Abs(t.etabin()-etabin()) > 1 ) return false;

  Int_t phi1 = t.phibin();
  Int_t phi2 = phibin();
  
  Int_t dphi = phi1 - phi2;

  /// neighbors in phi
  if ( TMath::Abs(dphi) <= 1 ) return true;

  /// test for cyclic neighbor in phi
  if ( phi1==59 ) phi1=-1;
  if ( phi2==59 ) phi2=-1;
  dphi = phi1 - phi2;

  /// neighbors in phi
  if ( TMath::Abs(dphi) <= 1 ) return true;

  return false;

}
// ----------------------------------------------------------------------------
TVector3 StEEmcTower::direction( Float_t zvertex ) const 
{
  TVector3 vertex(0.,0.,zvertex);
  const EEmcGeomSimple &geom=EEmcGeomSimple::Instance();
  TVector3 position=geom.getTowerCenter( (UInt_t)sector(), (UInt_t)subsector(), (UInt_t)etabin() ); 
  position -= vertex;
  return position.Unit();
}
