#include "StEEmcGenericClusterMaker.h"

#include "TH1F.h"
#include "TH2F.h"

#include "StEEmcPool/StEEmcA2EMaker/StEEmcA2EMaker.h"

#include "StEvent.h"
#include "StEmcDetector.h"
#include "StEmcCollection.h"

#include "StMessMgr.h"

ClassImp(StEEmcGenericClusterMaker);

// ----------------------------------------------------------------------------
StEEmcGenericClusterMaker::StEEmcGenericClusterMaker( const Char_t *name, StEEmcA2EMaker *a2e ) :
  StMaker(name)
{

  mEEanalysis=a2e; /* adc to energy maker */

  //  mCollection=new StEEmcClusterCollection();

  mEEmcGeom=new EEmcGeomSimple();
  mESmdGeom=EEmcSmdGeom::instance();
  mESmdMap =EEmcSmdMap::instance();

  ///
  /// Initialize storage banks for tower and SMD clusters.
  /// Both will be stored sector-wise and layer-wise.
  ///
  /// Tower clusters are stored according to the sector
  /// in which we find their seed tower.  
  ///
  /// Layers: T=0, P=1, Q=2, R=3   /  Planes: U=0, V=1
  ///
  /// mTowerClusters[sector][layer] returns vector of tower clusters
  /// mSmdClusters[sector][plane]   returns vector of SMD clusters
  ///

  StEEmcClusterVec_t t;
  std::vector< StEEmcClusterVec_t > layers;
  for ( Int_t i = 0; i < 4; i++ ) layers.push_back(t);
  for ( Int_t i = 0; i < 12; i++ ) mTowerClusters.push_back(layers);

  StEEmcSmdClusterVec_t s;
  std::vector< StEEmcSmdClusterVec_t > planes;
  planes.push_back(s);
  planes.push_back(s);
  for ( Int_t i = 0; i < 12; i++ ) mSmdClusters.push_back(planes);

  mSmdMatchRange = 40;
  Clear();

}

// ----------------------------------------------------------------------------
Int_t StEEmcGenericClusterMaker::Init()
{

  // If it's missing from the constructor, obtain pointer to 
  // the "default" eemc adc 2 energy maker.  If it doesn't
  // exist, crash and burn.
  if (!mEEanalysis) mEEanalysis=(StEEmcA2EMaker*)GetMaker("EEmcAdc2Energy");
  assert(mEEanalysis);// no input calib maker

  // book QA histograms
  const Char_t *layers[]={"T","P","Q","R","U","V"};
  const Char_t *names[]={"tower","preshower 1","preshower 2","postshower","smdu","smdv"};
  const Float_t max[]={50.0,200.0,400.0,200.0,2500.0,2500.0};
  for ( Int_t i=0;i<6;i++ )
    {

      hNumberOfClusters[i] = 
	new TH1F(TString("hNumberOfClusters")+layers[i],
		 TString("Number of ")+names[i]+" clusters",
		 20,0.,20.);

      hClusterSpectra[i] = 
	new TH1F(TString("hClusterSpectra")+layers[i],
		 TString(names[i])+" cluster energy spectra",
		 50,0.,max[i]);
    }

  return StMaker::Init();
}

// ----------------------------------------------------------------------------
Int_t StEEmcGenericClusterMaker::Make() // runs after user's Make()
{
 
  makeClusterMap();
  makeStEvent();
  makeHistograms();

  // add clusters to collection

  /***  
  mCollection->setNumberOfClusters(numberOfClusters());
  for ( Int_t sec=0;sec<12;sec++ )
    {
      for ( Int_t lay=0;lay<4;lay++ )
	{
	  for ( UInt_t i=0;i<mTowerClusters[sec][lay].size();i++ )
	    {
	      mCollection->add( mTowerClusters[sec][lay][i] );
	    }
	}
      for ( Int_t pln=0;pln<2;pln++ )
	{
	  for ( UInt_t i=0;i<mSmdClusters[sec][pln].size();i++ )
	    {
	      mCollection->add( mSmdClusters[sec][pln][i] );
	    }
	}
    }
  ***/
    
  return kStOK;

}

void StEEmcGenericClusterMaker::makeHistograms()
{
  /*
   * Fill QA histograms
   */
  for ( Int_t layer=0;layer<4;layer++ )
    {

      hNumberOfClusters[layer] -> Fill( mNumberOfClusters[layer] );
      for ( Int_t sec=0;sec<12;sec++ ) 
	for ( UInt_t i=0;i<mTowerClusters[sec][layer].size(); i++ ) 
	  {
	    StEEmcCluster cluster=mTowerClusters[sec][layer][i];
	    hClusterSpectra[layer]->Fill( cluster.energy() );
	  }

    }
  


  
}

void StEEmcGenericClusterMaker::makeClusterMap()
{

  
  // loop over all sectors
  for ( UInt_t sec=0;sec<12;sec++ )
    {

      // get vectors of clusters in this sector
      StEEmcClusterVec_t    &Tclusters = mTowerClusters[sec][0];
      StEEmcClusterVec_t    &Pclusters = mTowerClusters[sec][1];
      StEEmcClusterVec_t    &Qclusters = mTowerClusters[sec][2];
      StEEmcClusterVec_t    &Rclusters = mTowerClusters[sec][3];
      StEEmcSmdClusterVec_t &Uclusters = mSmdClusters[sec][0];
      StEEmcSmdClusterVec_t &Vclusters = mSmdClusters[sec][1];

      // loop over all cluster pairs and build associative maps
      for ( UInt_t t=0;t<Tclusters.size(); t++ )// towers
	{

	  StEEmcCluster &c = Tclusters[t];
	  EEmatch *mymatch = &mClusterMap[c.key()];

	  for ( UInt_t p=0;p<Pclusters.size();p++ )	    
	    if ( match(c,Pclusters[p]) ) mymatch->pre1.push_back(Pclusters[p]);
	  for ( UInt_t q=0;q<Qclusters.size();q++ )	    
	    if ( match(c,Qclusters[q]) ) mymatch->pre2.push_back(Qclusters[q]);
	  for ( UInt_t r=0;r<Rclusters.size();r++ )	    
	    if ( match(c,Rclusters[r]) ) mymatch->post.push_back(Rclusters[r]);
	  for ( UInt_t u=0;u<Uclusters.size();u++ )
	    if ( match(c,Uclusters[u]) ) mymatch->smdu.push_back(Uclusters[u]);
	  for ( UInt_t v=0;v<Vclusters.size();v++ )
	    if ( match(c,Vclusters[v]) ) mymatch->smdv.push_back(Vclusters[v]);


	  for ( UInt_t p=0;p<Pclusters.size();p++ )	    
	    if ( match(c,Pclusters[p]) ) c.addMatch( Pclusters[p].key(), 1 );
	  for ( UInt_t q=0;q<Qclusters.size();q++ )	    
	    if ( match(c,Qclusters[q]) ) c.addMatch( Qclusters[q].key(), 2 );
	  for ( UInt_t r=0;r<Rclusters.size();r++ )	    
	    if ( match(c,Rclusters[r]) ) c.addMatch( Rclusters[r].key(), 3 );
	  for ( UInt_t u=0;u<Uclusters.size();u++ )
	    if ( match(c,Uclusters[u]) ) c.addMatch( Uclusters[u].key(), 4 );
	  for ( UInt_t v=0;v<Vclusters.size();v++ )
	    if ( match(c,Vclusters[v]) ) c.addMatch( Vclusters[v].key(), 5 );
	   
	}

#if 0
      // print matches
      for ( UInt_t i=0;i<Tclusters.size();i++ )
	{

	  std::cout << "sec=" << sec << " i=" << i << " ------------------------------------" << std::endl;
	  StEEmcCluster c=Tclusters[i];
	  c.print();
	  Int_t k=c.key();
	  EEmatch m=mClusterMap[k];
	  for ( UInt_t j=0;j<m.pre1.size();j++ )	    
	    m.pre1[j].print();
	  for ( UInt_t j=0;j<m.pre2.size();j++ )
	    m.pre2[j].print();
	  for ( UInt_t j=0;j<m.post.size();j++ )
	    m.post[j].print();
	  for ( UInt_t j=0;j<m.smdu.size();j++ )
	    m.smdu[j].printLine();
	  for ( UInt_t j=0;j<m.smdv.size();j++ )
	    m.smdv[j].printLine();
	   
	}
    
#endif
#if 0
      // print matches
      for ( UInt_t i=0;i<Tclusters.size();i++ )
	{
	  std::cout << "sec=" << sec << " i=" << i << " ------------------------------------" << std::endl;
	  StEEmcCluster c=Tclusters[i];
	  c.print();

	  std::cout << "pre1 matches: ";
	  for ( Int_t j=0;j<c.numberOfMatchingClusters(1);j++ ) 
	    std::cout << c.getMatch(j,1) << " ";
	  std::cout << std::endl;

	  std::cout << "pre2 matches: ";
	  for ( Int_t j=0;j<c.numberOfMatchingClusters(2);j++ ) 
	    std::cout << c.getMatch(j,2) << " ";
	  std::cout << std::endl;


	  std::cout << "post matches: ";
	  for ( Int_t j=0;j<c.numberOfMatchingClusters(3);j++ ) 
	    std::cout << c.getMatch(j,3) << " ";
	  std::cout << std::endl;


	}
#endif


    }





}

void StEEmcGenericClusterMaker::makeStEvent()
{

  /*
   * create StEmcClusters and fill StEvent (to be added)
   */

  StEvent *stevent=(StEvent*)GetInputDS("StEvent");
  if ( !stevent ) {
    // assume we're running on MuDst and no collections need populated
    return;
  }
  
  ///
  /// First the eemc tower clusters
  ///
  StEmcDetector *detector=stevent->emcCollection()->detector(kEndcapEmcTowerId);
  if ( !detector )
    {
      Warning("fillStEvent","detector == NULL, MAJOR StEvent problem, continuing");
      return;
    }
  ///
  /// Next the pre and postshower clusters
  ///
  detector=stevent->emcCollection()->detector(kEndcapEmcPreShowerId);
  if ( !detector )
    {
      Warning("fillStEvent","detector == NULL for pre/post, no clusters for you");
      return;
    }
  ///
  /// Finally the U&V smd clusters
  ///
  StDetectorId ids[]={ kEndcapSmdUStripId, kEndcapSmdVStripId };

  
  for ( Int_t iplane=0; iplane<2; iplane++ )
    {

      detector=stevent->emcCollection()->detector(ids[iplane]);
      if ( !detector )
        {
          Warning("fillStEvent","detector == NULL for smd plane, no clusters for you");
	  return;
        }

    }


  return;
}

// ----------------------------------------------------------------------------
void StEEmcGenericClusterMaker::Clear(Option_t *opts)
{

  // clears storage arrays
  for ( Int_t sector=0; sector<12; sector++ ) {
    for ( Int_t layer=0; layer<4; layer++ )
      mTowerClusters[sector][layer].clear();
    for ( Int_t plane=0; plane<2; plane++ )
      mSmdClusters[sector][plane].clear();
  }
  for ( Int_t i=0;i<6;i++ ) mNumberOfClusters[i]=0;

  // clear the cluster map
  mClusterMap.clear();
  
  // reset the cluster id
  mClusterId = 0;

  // clear the collection
  //  mCollection->Clear();
  

}

// ----------------------------------------------------------------------------
void StEEmcGenericClusterMaker::add( StEEmcCluster &c )
{

  Int_t key    = nextClusterId();     /* next available cluster id */
  Int_t sector = c.tower(0).sector(); /* sector of seed tower      */
  Int_t layer  = c.tower(0).layer();  /* layer of the cluster      */
  c.key(key);
  mTowerClusters[sector][layer].push_back(c);
  
  // add this cluster to the cluster map
  EEmatch match;
  match.tower.push_back(c);
  mClusterMap[ key ] = match;

  mNumberOfClusters[c.tower(0).layer()]++;

}

// ----------------------------------------------------------------------------
Bool_t StEEmcGenericClusterMaker::match( StEEmcCluster &c1, StEEmcCluster &c2 )
{  
  return c1.tower(0).isNeighbor( c2.tower(0) );
}

Bool_t StEEmcGenericClusterMaker::match( StEEmcCluster &c1, StEEmcSmdCluster &c2 )
{
  StEEmcTower seed=c1.tower(0);
  if ( seed.sector() != c2.sector() ) return false;

  Int_t min[2],max[2];
  mESmdMap -> getRangeU( seed.sector(), seed.subsector(), seed.etabin(), min[0], max[0] );
  mESmdMap -> getRangeV( seed.sector(), seed.subsector(), seed.etabin(), min[1], max[1] );

  Int_t   plane = c2.plane();
  Float_t mean  = c2.mean();
  Int_t   mid   = (max[plane]+min[plane])/2;
  Float_t del   = mean - mid;
   
  return TMath::Abs(del)<mSmdMatchRange;
}

Bool_t StEEmcGenericClusterMaker::match( StEEmcSmdCluster &c1, StEEmcSmdCluster &c2 )
{

  Bool_t myMatch = false;
  if ( ( c1.energy() > 0.8 * c2.energy()   &&
	 c1.energy() < 1.2 * c2.energy() ) ||
       ( c2.energy() > 0.8 * c1.energy()   &&
	 c2.energy() < 1.2 * c1.energy() ) ) myMatch = true;

  if ( !myMatch ) return false;

  for ( Int_t sec=0;sec<12;sec++ )
    for ( UInt_t i=0; i < mTowerClusters[sec][0].size(); i++ )
      {

	StEEmcCluster &c=mTowerClusters[sec][0][i];
	if ( match ( c, c1 ) && match( c, c2 ) ) return true;

      }

  return false;

}

// ----------------------------------------------------------------------------
