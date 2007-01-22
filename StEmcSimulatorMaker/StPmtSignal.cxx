#include <math.h>
#include <StMessMgr.h>
#include "StPmtSignal.h"
ClassImp(StPmtSignal)

// Define work global variables => See code ipmsig.F (author V.Rykov)
const Int_t maxDyn=15;  // Max number of dynodes
Float_t sec[maxDyn],dyngain[maxDyn+1],g1[maxDyn+1],gdyn[maxDyn+1],dnw[maxDyn+1];
Float_t dn, sne, r;
Int_t nel, j;
// Rykov's notation
Int_t ndynpar, ndyn;
Float_t pmtgain, cnoise, dnoise;
Float_t gain, adcped, gnoise;

Float_t gausMean=0.0, gausRms=1.0;
// _____________________________________________________________________
StPmtSignal::StPmtSignal()
{
    mPmt = 1;    // Default PMT is ONE
    initPmtOne();
}
// _____________________________________________________________________
StPmtSignal::StPmtSignal(Int_t iPmt)
{
    mPmt = iPmt;
    switch (iPmt)
    {
    case 1:
        initPmtOne();
        break;
    default:
        mPmt = 1;
        LOG_ERROR << "StPmtSignal::StPmtSignal : Bad value of type of PMT \n" <<iPmt
        <<" , change to default value => "<<mPmt;
    }
}
// _____________________________________________________________________
void StPmtSignal::initPmtOne()
{
    Float_t base[11]={2.,2.,1.,1.,1.,1.,1.,1., 2., 3., 4.};
    SetName("HAMAMATSU-R6427");
    SetTitle("Bases used in Beam test-98");
    mNumDynodes      = 11;
    mPmtGain         = 1.5e+6;
    mProbPhotoCatode = 0.0;
    mProbPhotoDynode = 0.0;
    mNodeVoltage.Set(mNumDynodes,&base[0]);
    init();
}
// _____________________________________________________________________
void StPmtSignal::init() // See impsig_init
{
    //...Secondary electron conversion coefficients
    sne = 1.;
    for(j=0; j<mNumDynodes; j++)
    {
        sne = sne * mNodeVoltage[j];
    }
    sne = ::pow(double(mPmtGain/sne),1./(mNumDynodes));
    for(j=0; j<mNumDynodes; j++)
    {
        sec[j] = sne*mNodeVoltage[j];
    }

    //...Inversed partial gains after each dynode and "g1=G-1" constants (see SN301)
    dyngain[0]      = 1.0;
    g1[mNumDynodes] = 0.0;
    for(j=0; j<mNumDynodes; j++)
    {
        dyngain[j+1] = dyngain[j]/sec[j];
        g1[mNumDynodes-j-1] = (1. + g1[mNumDynodes-j])/sec[mNumDynodes-j-1];
    }
    //...Constants for the dnoise contribution calculation in the "fast" version
    dnw[mNumDynodes]  = 0.0;
    gdyn[mNumDynodes] = 0.0;
    for(j=mNumDynodes; j>0; j--)
    {
        dnw[j-1]  = mProbPhotoDynode*(1.+g1[j])*dyngain[j]*dyngain[j] + dnw[j];
        gdyn[j-1] = g1[j-1]*dyngain[j-1];
    }
    //...Shifts "dnoise" to avoid POISSN calls at mean=0.
    dn = mProbPhotoDynode + 1.0e-18;
}
// _____________________________________________________________________
Int_t StPmtSignal::getAdc(Int_t nphe, Int_t iver) // See impsig
{
    Int_t adc;
    if(mPmt==0)
    {
        LOG_ERROR << Form("Define PMT characteristics !!!") << endm;
        return -9999;
    }
    // Rykov's notation
    if(cnoise<1.e-18)
        nel=nphe;
    else
        nel = nphe + mRandom.Poisson(cnoise);

    // ------ F A S T   V E R S I O N ----------------------
    if(iver==0)
    {
        j=0;
        while(j<ndyn && nel<100)
        {
            sne=sec[j]*float(nel) + dn;
            nel=mRandom.Poisson(sne);
            j=j+1;
        }
        //...Normalized signal with noise
        sne = gain*float(nel)*dyngain[j];
        r   = mRandom.Gaus(gausMean,gausRms);
        sne = sne + gain*dnoise*gdyn[j] +
              r*::sqrt(gain*(sne*gdyn[j]+gain*dnw[j])+gnoise*gnoise);
    }
    // -----  F U L L   V E R S I O N ----------------------
    else
    {
        for(j=0;j<ndyn; j++)
        {
            sne = sec[j]*float(nel) + dn;
            nel = mRandom.Poisson(sne);
        }
        sne = mRandom.Gaus(gausMean,gausRms);
        sne = gain*float(nel)*dyngain[ndyn] + gnoise*sne;
    }
    sne    = sne + adcped;
    adc    = int(sne);
    if(sne<0.0)
        adc = adc - 1;
    return adc;
}
// _____________________________________________________________________
void StPmtSignal::printParameters()
{
    LOG_INFO << Form("Coefficient for transition from number of phe to adc (gain) => %f",
           mFromPheToAdc) << endm;
    LOG_INFO << Form(" Mean value of ADC pedestal(adcped) (ADC-counts) %7.3f", mMeanAdcPed) << endm;
    LOG_INFO << Form(" RMS  value of ADC pedestal(gnoise) (ADC-counts) %7.3f", mRmsAdcPed) << endm; 
}
// _____________________________________________________________________
void StPmtSignal::setAllParameters(Float_t gainw, Float_t adcpedw, Float_t gnoisew)
{
    mFromPheToAdc = gainw;
    mMeanAdcPed   = adcpedw;
    mRmsAdcPed    = gnoisew;
    // Rykov notation => for convenience only
    gain          = gainw;
    adcped        = adcpedw;
    gnoise        = gnoisew;

    ndynpar       = maxDyn;
    ndyn          = mNumDynodes;
    pmtgain       = mPmtGain;
    cnoise        = mProbPhotoCatode;
    dnoise        = mProbPhotoDynode;
}
// _____________________________________________________________________
void StPmtSignal::print(Int_t iprint=1)
{
    Int_t i;
    LOG_INFO << Form(" Type of PMT   => %30s Number %i",GetName(),mPmt) << endm;
    LOG_INFO << Form(" Type of bases => %30s",GetTitle()) << endm;
    LOG_INFO << Form(" Number of dynodes %2i",mNumDynodes) << endm;
    LOG_INFO << Form(" PMT Gain          %12.5e",mPmtGain) << endm;
    LOG_INFO << Form(" Noise on photocathode(cnoise) %8.6f",mProbPhotoCatode) << endm;
    LOG_INFO << Form(" Noise on       dynode(dnoise) %8.6f",mProbPhotoCatode) << endm;
    LOG_INFO << Form(" Relative Voltage distribution (base)") << endm;
    for(i=0;i<mNumDynodes;i++)
        LOG_INFO << Form("     i %2i %4.1f",i+1,mNodeVoltage[i]) << endm;
    if(iprint>=10)
    {  // Print also work arrays
        LOG_INFO << Form(" sec ") << endm;
        for(i=0; i<mNumDynodes; i++) {
            LOG_INFO << Form(" i %2i %20.7e",i,sec[i]) << endm;
		}
        LOG_INFO << Form(" dyngain ") << endm;
        for(i=0; i<=mNumDynodes; i++) {
            LOG_INFO << Form(" i %2i %20.7e",i,dyngain[i]) << endm;
		}
        LOG_INFO << Form(" g1 ") << endm;
        for(i=0; i<=mNumDynodes; i++) {
            LOG_INFO << Form(" i %2i %20.7e",i,g1[i]) << endm;
		}
        LOG_INFO << Form(" gdyn ") << endm;
        for(i=0; i<=mNumDynodes; i++) {
            LOG_INFO << Form(" i %2i %20.7e",i,gdyn[i]) << endm;
		}
        LOG_INFO << Form(" dnw ") << endm;
        for(i=0; i<=mNumDynodes; i++) {
            LOG_INFO << Form(" i %2i %20.7e",i,dnw[i]) << endm;
		}
    }
}
