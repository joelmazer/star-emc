#include "StPi0ResultsUtil.h"

#include "TDataProcessorMemberInspector.h"

#include <fstream>
#include <list>
#include <algorithm>
using namespace std;

#include <TObject.h>
#include <TROOT.h>
#include <TFile.h>
#include <TMatrix.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TImage.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TGraphSmooth.h>
#include <TDatime.h>
#include <TStopwatch.h>
#include <TFolder.h>
#include <TBrowser.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TMath.h>

#include <StEmcPool/StPi0Common/StPi0CommonUtil.h>
#include <StEmcPool/StPi0Common/Logger.h>

#include <StEmcPool/StPi0Analysis/TCuts.h>
#include <StEmcPool/StPi0Analysis/TCandidateDataProcessor.h>
#include <StEmcPool/StPi0Analysis/TPointDataProcessor.h>
#include <StEmcPool/StPi0Analysis/TEventDataProcessor.h>
#include <StEmcPool/StPi0Analysis/TSimuDataProcessor.h>
#include <StEmcPool/StPi0Analysis/TMCGammaDataProcessor.h>
#include <StEmcPool/StPi0Analysis/TDataProcessorPool.h>
#include <StEmcPool/StPi0Analysis/TInvariantMassDistribution.h>
#include <StEmcPool/StPi0Analysis/TWeightCalculator.h>
#include <StEmcPool/StPi0Analysis/StPi0AnalysisUtil.h>
#include <StEmcPool/StPi0Analysis/StPi0AnalysisUtilDB.h>
#include <StEmcPool/StPi0Analysis/TBinStatistics.h>                                                                                                                                      

// PHENIX d+Au pi0
Float_t Px[] =   {1.25,     1.75,     2.25,     2.75,     3.25,     3.75,     4.25,     4.75,     5.25,     5.75,     6.25,     6.75,     7.25,     7.75,     8.25,     8.75,     9.25,     9.75};
Float_t Pxe[] =  {0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0};
Float_t Py[] =   {5.34e-02, 1.02e-02, 2.31e-03, 6.49e-04, 2.07e-04, 7.33e-05, 2.56e-05, 1.11e-05, 4.67e-06, 2.17e-06, 1.10e-06, 5.78e-07, 3.27e-07, 1.89e-07, 9.56e-08, 7.07e-08, 2.48e-08, 1.86e-08};
Float_t Pye[] =  {3.0e-03,  6.1e-04,  1.5e-04,  4.5e-05,  1.5e-05,  8.1e-06,  2.8e-06,  1.2e-06,  5.1e-07,  2.4e-07,  1.2e-07,  6.5e-08,  3.9e-08,  2.3e-08,  1.3e-08,  1.3e-08,  5.7e-09,  6.1e-09};
Bool_t showPhenixPi0 = true;
const Char_t *showPhenixPi0Legend = "PHENIX d+Au #pi^{0}";

// PHENIX p+p pi0
Float_t Pppx[] =   {1.215,     1.719,     2.223,     2.726,     3.228,     3.730,     4.232,     4.733,     5.234,     5.735,     6.236,     6.737,     7.452,     8.457,     9.460,     10.861,    13.250};
Float_t Pppxe[] =  {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0};
Float_t Pppy[] =   {3.733e-01, 6.052e-02, 1.221e-02, 3.308e-03, 9.978e-04, 3.385e-04, 1.187e-04, 4.726e-05, 2.206e-05, 1.113e-05, 4.999e-06, 3.003e-06, 1.080e-06, 4.853e-07, 1.643e-07, 5.227e-08, 1.190e-8};
//#define P(i,x,y) (Pppy[i] * sqrt(((x) * (x)) + ((y) * (y)) + ((9.6) * (9.6))) / 100.0) // stat + syst + norm
#define P(i,x,y) (Pppy[i] * x / 100.0) // stat only
Float_t Pppye[] =  {P(0,1.6,7.3), P(1,1.8,7.1), P(2,2.5,7.1), P(3,3.6,7.2), P(4,5.7,7.3), P(5,7.3,7.7), P(6,2.4,8.3), P(7,4.2,8.5), P(8,5.0,8.7), P(9,4.5,9.2), P(10,6.3,9.5), P(11,7.7,9.8), P(12,8.8,10.1), P(13,12.0,10.8), P(14,19.3,11.0), P(15,22.4,11.7), P(16,41.3,15.9)};
Bool_t showPhenixPi0PP = true;
const Char_t *showPhenixPi0PPLegend = "PHENIX p+p #pi^{0}";

// STAR (h+ + h-)/2
Float_t Sx[] =   {1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.50, 2.70, 2.90, 3.16, 3.56, 4.07, 4.71, 5.50, 6.44, 7.45, 8.45, 9.46, 10.84};
Float_t Sxe[] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t Sy[] =   {0.16, 0.11, 0.081, 0.058, 0.042, 0.030, 0.022, 0.016, 0.012, 0.0091, 0.0068, 0.0050, 0.0038, 0.0028, 0.0019, 0.0011, 0.00072, 0.00041, 0.00017, 0.000062, 0.000019, 0.0000052, 0.0000013, 0.00000040, 0.00000012, 0.000000050, 0.000000016};
Float_t Sye[] =  {0.01, 0.01, 0.004, 0.003, 0.002, 0.001, 0.001, 0.001, 0.001, 0.0004, 0.0003, 0.0003, 0.0002, 0.0002, 0.0001, 0.0001, 0.00005, 0.00002, 0.00001, 0.000004, 0.000002, 0.0000006, 0.0000002, 0.00000006, 0.00000002, 0.000000009, 0.000000003};
Bool_t showStarChargedHadrons = false;
const Char_t *showStarChargedHadronsLegend = "STAR d+Au h^{#pm}";

// STAR pi+ from TOFr and rdE/dx
Float_t SPpx[] =  {0.35,    0.45,    0.55,    0.65,     0.75,     0.85,     0.95,     1.05,     1.15,      1.3,       1.5,       1.7,       1.9,        2.25,       2.75,        3.25,        3.75,        4.25,        4.75,        5.5,         6.5,         7.5,         8.5,         9.5};
Float_t SPpxe[] = {0,       0,       0,       0,        0,        0,        0,        0,        0,         0,         0,         0,         0,          0,          0,           0,           0,           0,           0,           0,           0,           0,           0,           0};
Float_t SPpy[] =  {3.13955, 1.82226, 1.08348, 0.617851, 0.413013, 0.266037, 0.169814, 0.114,    0.0769641, 0.0442544, 0.0210009, 0.0104954, 0.00585984, 0.00277035, 0.000589763, 0.000191322, 6.35517e-05, 2.29327e-05, 9.66819e-06, 2.88435e-06, 8.65872e-07, 2.20034e-07, 9.55296e-08, 4.13572e-08};
Float_t SPpye[] = {0.026,   0.1349,  0.084,   0.048,    0.03248,  0.021,    0.01344,  0.0041,   0.0033,    0.0017,    1.1e-07,   0.00072,   0.00055,    0.00039,    9.8e-05,     2.57927e-06, 1.66775e-06, 1.15619e-06, 4.02157e-07, 1.63822e-07, 5.84802e-08, 1.27507e-08, 3.748e-08,   2.11504e-08};
// STAR pi- from TOFr and rdE/dx
Float_t SPmx[] =  {0.35,    0.45,    0.55,    0.65,     0.75,     0.85,     0.95,     1.05,     1.15,      1.3,       1.5,       1.7,       1.9,        2.25,       2.75,        3.25,        3.75,        4.25,        4.75,        5.5,         6.5,         7.5,         8.5,         9.5};
Float_t SPmxe[] = {0,       0,       0,       0,        0,        0,        0,        0,        0,         0,         0,         0,         0,          0,          0,           0,           0,           0,           0,           0,           0,           0,           0,           0};
Float_t SPmy[] =  {3.19879, 1.82017, 1.08039, 0.670021, 0.404875, 0.258436, 0.167793, 0.116018, 0.0762589, 0.0435504, 0.0203983, 0.0102947, 0.00515995, 0.00241967, 0.000500648, 0.000189875, 6.48901e-05, 2.29395e-05, 9.48129e-06, 3.14379e-06, 8.46819e-07, 2.39244e-07, 8.01216e-08, 3.86e-08};
Float_t SPmye[] = {0.037,   0.024,   0.017,   0.013,    0.009,    0.007,    0.0054,   0.0042,   0.0033,    0.0017,    1.1e-07,   0.00072,   0.0005,     0.00042,    8.8e-05,     2.38112e-06, 1.3687e-06,  2.03362e-06, 5.54402e-07, 1.58e-07,    7.29168e-08, 1.84825e-08, 1.24556e-08, 9.09118e-09};
Bool_t showStarChargedPions = false;
const Char_t *showStarChargedPionsPlusLegend = "STAR d+Au #pi^{+}";
const Char_t *showStarChargedPionsMinusLegend = "STAR d+Au #pi^{-}";

// pi0 from BEMC (Andre Mischke)
Float_t Ax[] =   {1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 6, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5};
Float_t Axe[] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t Ay[] =   {0.00500784, 0.00178272, 0.000527415, 0.000173036, 7.47168e-05, 2.65451e-05, 1.0569e-05, 4.08643e-06, 1.29798e-06, 7.02737e-07, 2.8363e-07, 1.22991e-07, 5.68615e-08, 2.77471e-08, 1.61125e-08, 9.93835e-09, 5.44627e-09, 3.15326e-09, 1.56401e-09};
Float_t Aye[] =  {0.000146705, 5.55874e-05, 2.38637e-05, 1.17154e-05, 1.44958e-05, 4.05473e-06, 2.9797e-06, 1.23376e-06, 4.30548e-07, 1.00577e-07, 5.20925e-08, 2.57462e-08, 6.45742e-09, 3.97804e-09, 4.34288e-09, 2.42936e-09, 3.26254e-09, 2.1746e-09, 1.13746e-09};
Bool_t showMischkePi0 = false;
const Char_t *showMischkePi0Legend = "STAR d+Au #pi^{0} (preliminary, A. Mischke)";

// pi0 from BEMC (Alex Stolpovsky d+Au)
Float_t ASxMB[] =   {0.666852, 1.17992, 1.68955, 2.19691, 2.7027, 3.20736, 3.7112, 4.21441, 4.71714, 5.21947, 5.7215, 6.22328};
Float_t ASyMB[] =   {0.367177, 0.0485134, 0.00847456, 0.0018273, 0.000447192, 0.000158307, 5.73226e-005, 2.45748e-005, 9.31934e-006, 3.22783e-006, 1.82639e-006, 9.72972e-007};
Float_t ASxeMB[] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t ASyeMB[] =  {0.050165, 0.00241783, 0.000301056, 5.7797e-005, 1.10759e-005, 7.44382e-006, 4.17976e-006, 2.23793e-006, 1.88963e-006, 1.308e-006, 4.24687e-007, 3.68492e-007};
Float_t ASxHT1[] =  {4.36722, 5.38446, 6.39782, 7.40845, 8.41711, 9.42429};
Float_t ASyHT1[] =  {1.84773e-005, 3.22546e-006, 8.7002e-007, 1.7933e-007, 7.78323e-008, 2.87567e-008};
Float_t ASxeHT1[] = {0, 0, 0, 0, 0, 0};
Float_t ASyeHT1[] = {1.48682e-006, 2.06413e-007, 7.23673e-008, 8.00081e-008, 3.48125e-008, 2.59426e-008};
Float_t ASxHT2[] =  {5.38446, 6.39782, 7.40845, 8.41711, 9.42429, 10.4303, 11.4355, 12.4399, 13.4438};
Float_t ASyHT2[] =  {4.37073e-006, 8.80739e-007, 2.26148e-007, 7.89908e-008, 3.36603e-008, 1.25895e-008, 4.68547e-009, 2.1466e-009, 5.31072e-010};
Float_t ASxeHT2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t ASyeHT2[] = {1.25666e-006, 9.57822e-008, 2.08191e-008, 3.4561e-008, 1.70188e-008, 1.01504e-008, 5.77644e-009, 3.57516e-009, 2.78541e-009};
Bool_t showStolpovskyPi0 = false;
const Char_t *showStolpovskyPi0Legend = "STAR d+Au #pi^{0} (A. Stolpovsky)";
// pi0 from BEMC (Alex Stolpovsky p+p)
Float_t ASppxMB[] =   {1.24381,     2.29837,     /*3.33504*/    3.33504}; // the last point is average of MB and HT1
Float_t ASppxeMB[] =  {0,           0,           /*0*/          0};
Float_t ASppyMB[] =   {0.00876258,  0.000235669, /*1.27265e-5*/ 1.65e-5};
Float_t ASppyeMB[] =  {0.000897222, 2.24085e-5,  /*2.92395e-6*/ 2.4e-6};
Float_t ASppxHT1[] =  {/*3.33504,*/ 4.36091,     5.37997,    6.39454,    7.40601,    8.41526,    9.42287};
Float_t ASppxeHT1[] = {/*0,      */ 0,           0,          0,          0,          0,          0};
Float_t ASppyHT1[] =  {/*2.04e-5,*/ 3.32231e-6,  5.47952e-7, 9.8459e-8,  2.67994e-8, 1.25404e-8, 6.42752e-9};
Float_t ASppyeHT1[] = {/*2.48e-6,*/ 2.84655e-7,  4.08232e-8, 1.15377e-8, 2.05424e-8, 8.69315e-9, 5.1027e-9};
Bool_t showStolpovskyPi0PP = false;
const Char_t *showStolpovskyPi0PPLegend = "STAR p+p #pi^{0} (A. Stolpovsky)";

Bool_t showpQCD = true;
const Char_t *showpQCDLegend = "NLO pQCD (W. Vogelsang)";
const Char_t *showpQCDErrorsLegend = "pQCD scale uncertainty";

// STAR (h+ + h-)/2 Rcp
Float_t SRcpx[] =  {0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.5, 2.7, 2.9, 3.175, 3.575, 4.1, 4.75, 5.55, 6.5};
Float_t SRcpxe[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
// 0 <= eta < 0.5
Float_t SRcp1y[] =  {0.717, 0.730, 0.753, 0.781, 0.811, 0.838, 0.873, 0.898, 0.924, 0.948, 0.970, 0.989, 1.006, 1.019, 1.040, 1.048, 1.066, 1.061, 1.074, 1.101, 1.058, 1.102, 1.082, 1.087, 1.084, 1.095, 1.051, 1.116, 1.042, 0.856, 0.739};
Float_t SRcp1ye[] = {0.022, 0.019, 0.018, 0.018, 0.018, 0.018, 0.019, 0.019, 0.020, 0.020, 0.021, 0.021, 0.022, 0.022, 0.023, 0.023, 0.024, 0.024, 0.025, 0.026, 0.026, 0.028, 0.026, 0.029, 0.034, 0.035, 0.037, 0.051, 0.055, 0.058, 0.080};
// 0.5 <= eta < 1.0
Float_t SRcp2y[] =  {0.737, 0.750, 0.769, 0.786, 0.809, 0.832, 0.860, 0.874, 0.896, 0.914, 0.929, 0.941, 0.962, 0.975, 0.970, 0.998, 1.014, 1.006, 1.024, 1.042, 1.008, 1.035, 1.015, 1.070, 1.049, 1.031, 1.050, 1.020, 1.072, 0.919, 0.721};
Float_t SRcp2ye[] = {0.029, 0.024, 0.021, 0.020, 0.019, 0.019, 0.019, 0.020, 0.020, 0.020, 0.021, 0.021, 0.021, 0.022, 0.022, 0.023, 0.023, 0.023, 0.024, 0.026, 0.026, 0.029, 0.025, 0.028, 0.029, 0.028, 0.032, 0.036, 0.053, 0.065, 0.086};
Bool_t showStarRcp = true;
const Char_t *showStarRcpLegend1 = "STAR h^{#pm} R^{Au}_{CP}, 0 #leq #eta < 0.5";
const Char_t *showStarRcpLegend2 = "STAR h^{#pm} R^{Au}_{CP}, 0.5 #leq #eta < 1.0";

// PHENIX pi0 RdA
Float_t PdAx[] =  {1.25,  1.75,  2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,  6.25, 6.75,  7.25, 7.75, 8.25,  8.75, 9.25,  9.75};
Float_t PdAxe[] = {0,     0,     0,    0,    0,    0,    0,    0,    0,    0,     0,    0,     0,    0,    0,     0,    0,     0};
Float_t PdAy[] =  {0.817, 0.927, 1.02, 1.03, 1.08, 1.12, 1.11, 1.20, 1.08, 0.991, 1.11, 0.972, 1.13, 1.32, 0.748, 1.07, 0.565, 0.797};
#define S(i, x, y) (sqrt(((x) * (x)) + ((y) * (y)) + (0.1183 * 0.1183 * (PdAy[i]) * (PdAy[i]))))
Float_t PdAye[] = {S(0,0.053,0.034), S(1,0.06,0.043), S(2,0.069,0.053), S(3,0.076,0.059), S(4,0.094,0.066), S(5,0.16,0.072), S(6,0.12,0.075), S(7,0.12,0.084), S(8,0.087,0.11), S(9,0.078,0.1), S(10,0.1,0.11), S(11,0.1,0.097), S(12,0.15,0.11), S(13,0.22,0.13), S(14,0.14,0.076), S(15,0.27,0.11), S(16,0.18,0.06), S(17,0.49,0.086)};
Bool_t showPhenixPi0dA = true;
const Char_t *showPhenixPi0dALegend = "PHENIX #pi^{0} R_{dA} (preliminary)";

// Martijn's pi0 in d+Au
Float_t MartijnxMB[] = {1.25, 1.75, 2.25, 2.75, 3.25, 3.75/*, 4.5*/};
Float_t MartijnxeMB[] = {  0,    0,    0,    0,    0,    0/*,   0*/};
Float_t MartijnyMB[] = {0.0632152, 0.00995222, 0.00242639, 0.000654909, 0.000230016, 8.16159e-05/*, 1.78216e-05*/};
Float_t MartijnyeMB[] = {0.0119605, 0.00107061, 0.000189758, 5.37487e-05, 1.79199e-05, 7.84226e-06/*, 1.94046e-06*/};
Float_t MartijnxHT1[] = {4.5, 5.5/*, 6.5, 7.5, 8.5, 9.5, 10.5*/};
Float_t MartijnxeHT1[] = { 0,   0/*,   0,   0,   0,   0,    0*/};
Float_t MartijnyHT1[] = {1.70626e-05, 3.662e-06/*, 8.90958e-07, 3.04447e-07, 9.95878e-08, 3.46795e-08, 2.24438e-08*/};
Float_t MartijnyeHT1[] = {1.51102e-06, 2.53059e-07/*, 6.49528e-08, 2.69047e-08, 1.30154e-08, 6.88022e-09, 5.34097e-09*/};
Float_t MartijnxHT2[] = {6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.75, 14.25};
Float_t MartijnxeHT2[] = { 0,   0,   0,   0,    0,    0,     0,     0};
Float_t MartijnyHT2[] = {8.66779e-07, 2.42148e-07, 8.04315e-08, 3.40399e-08, 1.88435e-08, 7.50758e-09, 3.99289e-09, 7.20376e-10};
Float_t MartijnyeHT2[] = {7.87554e-08, 1.60149e-08, 5.5827e-09, 2.7661e-09, 1.78814e-09, 1.03377e-09, 5.51807e-10, 2.32924e-10};
Bool_t showMartijnPi0dAu = true;
const Char_t *showMartijnPi0dAuLegend = "STAR d+Au #pi^{0} (M. Russcher)";

// pi0 crossection in p+p 2005 from Frank Simon
Float_t FSxs_x[] =  {1.233,      1.72958,     2.22885,     2.72693,     3.42412,     4.00131,     4.44584,     4.94266,     5.59561,     6.479,       7.47177,     8.90215,     10.9064,    12.9013,     15.3554};
Float_t FSxs_xe[] = {0,          0,           0,           0,           0,           0,           0,           0,           0,           0,           0,           0,           0,          0,           0};
Float_t FSxs_y[] =  {0.241395,   0.0383293,   0.00886822,  0.0027491,   0.000552523, 0.000151854, 7.55546e-05, 3.15377e-05, 1.12085e-05, 3.22457e-06, 9.87101e-07, 2.64197e-07, 5.7152e-08, 1.45261e-08, 3.08004e-09};
Float_t FSxs_ye[] = {0.00513025, 0.000993099, 0.000333997, 0.000152168, 4.10583e-05, 3.57657e-06, 1.53147e-06, 6.87439e-07, 1.60763e-07, 5.3272e-08,  2.35091e-08, 7.71482e-09, 3.2278e-09, 1.47134e-09, 4.71187e-10};
Bool_t showSimonPi0PPCrossection = true;
const Char_t *showSimonPi0PPCrossectionLegend = "#pi^{0} p+p 2005, Frank Simon";

const Char_t *fitFuncZero = "gaus(0) + gaus(3)";
const Char_t *fitFuncConst = "gaus(0) + gaus(3) + [6]";
const Char_t *fitFuncLinear = "gaus(0) + gaus(3) + [6] + ([7]*x)";

#define PQCDFUNC(usePPpQCD)   ((Double_t (*)(Double_t*, Double_t*))((usePPpQCD) ? (&getNLOpQCDPP) : (&getNLOpQCD)))
#define PQCDWEIGHT(usePPpQCD) ((usePPpQCD) ? (pQCDPPweight) : (pQCDweight))

/*
const Char_t *fitFuncZero = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncConst = "[3] + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncLinear = "[3] + [4]*x + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncLinearZero = "[3]*x + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncParabola = "[3] + [4]*x + [5]*x*x + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncZero1x = "([4]/x) + ([3]*0) + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncConst1x = "([4]/x) + [3] + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncLinear1x = "([5]/x) + [3] + [4]*x + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncParabola1x = "([6]/x) + [3] + [4]*x + [5]*x*x + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))";
const Char_t *fitFuncConstExp = "[4]*exp(-x*[5]) + [3] + [0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; 

const Char_t *fitFuncZeroTwoPeaks = "([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))) + (0*([3]+[4]+[5]+[6])) + ([7]*exp(-0.5*((x-[8])/[9])*((x-[8])/[9]))) + [10]";
const Char_t *fitFuncConstTwoPeaks = "([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))) + (0*([3]+[4]+[5]+[6])) + ([7]*exp(-0.5*((x-[8])/[9])*((x-[8])/[9]))) + [10] + [11]";
const Char_t *fitFuncLinearTwoPeaks = "([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))) + (0*([3]+[4]+[5]+[6])) + ([7]*exp(-0.5*((x-[8])/[9])*((x-[8])/[9]))) + [10] + [11] + ([12]*x)";
const Char_t *fitFuncConstOnlyTwoPeaks = "([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))) + (0*([3]+[4]+[5]+[6]+[7]+[8]+[9])) + [10] + [11]";
const Char_t *fitFuncLinearOnlyTwoPeaks = "([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))) + (0*([3]+[4]+[5]+[6]+[7]+[8]+[9])) + [10] + [11] + ([12]*x)";
*/
const Char_t *fitFuncTrackDist = "0*gaus(0) + 0*gaus(3) + ([6]*x*exp(-x*[7]))";
const Char_t *fitFuncTrackDistPeak = "gaus(0) + 0*gaus(3) + ([6]*x*exp(-x*[7]))";
//const Char_t *fitFuncTrackDist = "([3]*x*exp(-x*x*[4])) + ([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))";
//const Char_t *fitFuncTrackDist = "([3]*x*exp(-x*(exp(-x*x*[5]))*[4])) + ([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))";
//const Char_t *fitFuncTrackDist = "([3]*x*( ([4]*exp(-x*exp(-x*[5])*[6])) + ((1-[4])*exp(-x*[6])) )) + ([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))";
//const Char_t *fitFuncTrackDist = "([3]*x*( ([4]*exp(-x*[6])*exp(-x*x*[5])) + ((1-[4])*exp(-x*[6])) )) + ([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))";

#define GET_DATA(var, class, pool, name) class *var = pool ? dynamic_cast<class*>(pool->getProcessor(name)) : 0;// cout << name << " -> " << var << endl;

ClassImp(TAnalysisSettingsTrigger);
ClassImp(TAnalysisSettings);
ClassImp(TAllSettings);
ClassImp(TAnalysisRawResultsTrigger);
ClassImp(TAnalysisResults);
ClassImp(TAllResults);

TH1F *histToFit1 = 0;
TH1F *histToFit2 = 0;
TH1F *histToFit3 = 0;
TH1F *histToFit4 = 0;
TH1F *histToFit5 = 0;
TH1F *histToFit6 = 0;
Double_t getHistWeightedSum(Double_t *x, Double_t *p) {
	Double_t result = 0.0;
	if (x && p) {
	    if (histToFit1) result += (p[0] * p[1] * histToFit1->GetBinContent(histToFit1->GetXaxis()->FindBin(x[0])));
	    if (histToFit2) result += (p[0] * (1.0 - p[1]) * histToFit2->GetBinContent(histToFit2->GetXaxis()->FindBin(x[0])));
	    if (histToFit3) result += (p[2] * histToFit3->GetBinContent(histToFit3->GetXaxis()->FindBin(x[0])));
	    if (histToFit4) result += (p[3] * histToFit4->GetBinContent(histToFit4->GetXaxis()->FindBin(x[0])));
	    if (histToFit5) result += (p[4] * histToFit5->GetBinContent(histToFit5->GetXaxis()->FindBin(x[0])));
	    if (histToFit6) result += (p[4] * p[5] * histToFit6->GetBinContent(histToFit6->GetXaxis()->FindBin(x[0])));
	}
        return result;
}

const Float_t noLimitMagic = 1.01010101;
Float_t getParameterValue(TFitSettings param, const TBinParameters &binParameters, Float_t *limitLow = 0, Float_t *limitHigh = 0) {
    Float_t binCenter = binParameters.getCenter();
    Float_t result = 0;
    if (param.func && (binCenter >= param.fitUseLow) && (binCenter < param.fitUseHigh)) {
	result = param.func->Eval(binCenter);
    } else {
	result = param.value + (param.valueDrift * binCenter);
    }
    if (limitLow) *limitLow = (param.useFit && param.wasFit) ? result : param.limitLow;
    if (limitHigh) *limitHigh = (param.useFit && param.wasFit) ? result : param.limitHigh;
    return result;
}

Int_t getBestRebin(const TH1F *hist, Float_t numBinsWeight, Float_t meanErrorWeight, Float_t rangeLeft, Float_t rangeRight, Int_t minBins) {
	Int_t nrebin = 1;
//	Int_t nrebinWidth = 1;
	if (((numBinsWeight != 0) || (meanErrorWeight != 0)) && hist) {
	    Float_t sumRebinMin = -1.0;
//	    Float_t widthRebinMin = -1.0;
	    Int_t nrebinMax = (hist->GetXaxis()->FindBin(rangeRight) - hist->GetXaxis()->FindBin(rangeLeft)) / minBins;
	    //cout << "nrebinMax = " << nrebinMax << endl;
	    for (Int_t ntemp = 1;ntemp < nrebinMax;ntemp++) {
		TH1F invRebin(*hist);
    		if (ntemp > 1) invRebin.Rebin(ntemp);
		Int_t binLeft = invRebin.GetXaxis()->FindBin(rangeLeft);
		Int_t binRight = invRebin.GetXaxis()->FindBin(rangeRight);
		Float_t binWidth = 0;
		Float_t errMean = 0;
		Int_t binNum = 0;
		for (Int_t ibin = binLeft;ibin < binRight;ibin++) {
		    Float_t val = invRebin.GetBinContent(ibin);
		    Float_t err = invRebin.GetBinError(ibin);
		    Float_t width = invRebin.GetBinWidth(ibin);
		    if (val != 0) {
			Float_t errRel = err / val;
		        errMean += errRel;
		        binWidth += width;
		        binNum += 1;
		    }
		}
		if (binNum != 0) errMean /= binNum;
		if (binNum != 0) binWidth /= binNum;
		Float_t sum = -1.0;
		//if (binWidth != 0) {
		    binWidth = binWidth / 0.01;
    		    errMean = errMean / 0.10;
		    sum = (binWidth * numBinsWeight) + (errMean * meanErrorWeight);
		    if (((sum >= 0.0) && (sum < sumRebinMin)) || (sumRebinMin < 0)) {
			sumRebinMin = sum;
			nrebin = ntemp;
		    }
		//}
		//cout << "ntemp = " << ntemp << ", (binWidth/0.01) = " << binWidth << ", (errMean/0.10) = " << errMean << ", sum = " << sum << endl;
/*
		TH1F invRebinToy(*distribution->getDistribution());
		Stat_t stats[4];
		invRebinToy.GetStats(stats);
		Float_t numberOfEntries = (stats[1] != 0) ? (stats[0] * stats[0] / stats[1]) : 1.0;
		Float_t area = invRebinToy.Integral();
		Float_t mean = (0.5 + Int_t(analysisSettings.peakPositionTrue / invRebinToy.GetXaxis()->GetBinWidth(1))) * invRebinToy.GetXaxis()->GetBinWidth(1);
		Float_t sigma = 0.023;
		for (Int_t i = 1;i <= invRebinToy.GetXaxis()->GetNbins();i++) {
		    Float_t x = invRebinToy.GetXaxis()->GetBinCenter(i);
		    Float_t val = area/sqrt(TMath::TwoPi()*sigma)*exp(-(x - mean)*(x - mean)/(2.0*sigma*sigma));
		    Float_t err = val / sqrt(numberOfEntries / (invRebinToy.GetXaxis()->GetNbins() / ntemp));
		    //Float_t fluct = (val != 0) ? gRandom->Gaus(1.0, err * 0.00 / val) : 1.0;
		    //val *= fluct;
		    //err *= fluct;
		    invRebinToy.SetBinContent(i, val);
		    invRebinToy.SetBinError(i, err);
		    //if (val != 0) cout << "--- " << val << " " << err << endl;
		}
    		if (ntemp > 1) invRebinToy.Rebin(ntemp);
		TF1 fg("fg","gaus");
		fg.SetParameter(0, area);
		fg.SetParLimits(0, 0, 2.0 * area);
		fg.SetParameter(1, mean);
		fg.SetParLimits(1, mean * 0.7, mean * 1.3);
		fg.SetParameter(2, sigma);
		fg.SetParLimits(2, sigma * 0.7, sigma * 1.3);
		invRebinToy.Fit(&fg, "QN");
		Float_t width = fg.GetParameter(2);
		Float_t widthErr = fg.GetParError(2);
		Float_t widthErrRel = widthErr / width;
		//cout << "ntemp = " << ntemp << ":\tsum = " << sum << ",\twidthErrRel = " << width << ", " << widthErr << ", " << widthErrRel << endl;
		if (((widthErrRel >= 0.0) && (widthErrRel < widthRebinMin)) || (widthRebinMin < 0)) {
		    widthRebinMin = widthErrRel;
		    nrebinWidth = ntemp;
		}
*/
	    }
	    //cout << "nrebin = " << nrebin << ", nrebinWidth = " << nrebinWidth << endl;
	    //nrebin = nrebinWidth;
	}
	return nrebin;
}

void setErrorOnZeroBins(TH1F *hist) {
    if (hist && hist->GetXaxis()) {
	Float_t valuePerEntry = hist->GetEntries() ? (hist->Integral() / hist->GetEntries()) : 0;
	for (Int_t mbin = 1;mbin < hist->GetXaxis()->GetNbins();mbin++) if (hist->GetBinContent(mbin) == 0) hist->SetBinError(mbin, valuePerEntry);
    }
}

void extractSignal(const TH1F &invSig, const TF1 &funcPureBg, Float_t peakLeftYield, Float_t peakRightYield
    , Bool_t useProportionalBgErrors
    , Bool_t useIndependentSigBgErrors
    , Bool_t useAnticorrelatedSigBgErrors
    , Bool_t normalizeSigBgErrors
    , Bool_t normalizeSigBgSquareErrors
    , Float_t &areaSigBg
    , Float_t &areaSigBgErr
    , Float_t &area
    , Float_t &areaErr
    , Float_t &areaBg
    , Float_t &areaBgErr
    , Float_t &mean
    , Float_t &meanErr
    , Float_t &width
    , Float_t &widthErr
) {
	areaSigBg = 0;
	areaSigBgErr = 0;
	area = 0;
	areaErr = 0;
	areaBg = 0;
	areaBgErr = 0;
	mean = 0;
	meanErr = 0;
	width = 0;
	widthErr = 0;
	if (peakLeftYield < invSig.GetXaxis()->GetBinLowEdge(1)) {
	    Float_t tmp = invSig.GetXaxis()->GetBinLowEdge(1);
	    //cout << "Yield left is out of range, changing from " << peakLeftYield << " to " << tmp << endl;
	    peakLeftYield = tmp;
	}
	if (peakRightYield > invSig.GetXaxis()->GetBinUpEdge(invSig.GetXaxis()->GetNbins())) {
	    Float_t tmp = invSig.GetXaxis()->GetBinUpEdge(invSig.GetXaxis()->GetNbins());
	    //cout << "Yield right is out of range, changing from " << peakRightYield << " to " << tmp << endl;
	    peakRightYield = tmp;
	}
	Int_t peakLeftYieldBin = invSig.GetXaxis()->FindBin(peakLeftYield);
	Int_t peakRightYieldBin = invSig.GetXaxis()->FindBin(peakRightYield);
	if ((peakRightYield <= peakLeftYield) || (peakLeftYieldBin < 1) || (peakRightYieldBin > invSig.GetXaxis()->GetNbins())) {
	    Int_t tmp = peakLeftYieldBin - 10;
	    //cout << "Yield range inconsistent, changing right bin from " << peakRightYieldBin << " to " << tmp << endl;
	    peakRightYieldBin = tmp;
	}

	Float_t sumW = 0;
	Float_t sumWX = 0;
	Float_t sumWW = 0;
	Float_t sumWWX = 0;
	Float_t sumWXX = 0;
	Float_t sumE = 0;
	Float_t sumEE = 0;

	for (Int_t mbin = peakLeftYieldBin;mbin <= peakRightYieldBin;mbin++) {
		Float_t xCenter = invSig.GetXaxis()->GetBinCenter(mbin);
		Float_t xLowEdge = invSig.GetXaxis()->GetBinLowEdge(mbin);
		Float_t xUpEdge = invSig.GetXaxis()->GetBinUpEdge(mbin);
		Float_t xWidth = invSig.GetXaxis()->GetBinWidth(mbin);
		Float_t xMin = max(peakLeftYield, xLowEdge);
		Float_t xMax = min(peakRightYield, xUpEdge);
		Float_t w = (xMax - xMin) / xWidth;
		if ((w < -0.0001) || (w > 1.0001)) {
		//if (w < 0.99) {
		    cout << "Bad bad bin!!! " << w << ": " << xLowEdge << ":" << xUpEdge << " " << peakLeftYield << ":" << peakRightYield << " " << xMin << ":" << xMax << " " << xWidth << endl;
		    continue;
		}
		if (w < 0.0001) w = 0.0;
		if (w > 0.999) w = 1.0;

		Float_t valSigBg = invSig.GetBinContent(mbin) * w;
		Float_t valSigBgErr = invSig.GetBinError(mbin) * w;
		Float_t valBg = funcPureBg.Eval(xCenter) * w;
		Float_t valBgErr = 0 * w;
		Float_t valSig = valSigBg - valBg;
		Float_t valSigErr = 0;
		if (useProportionalBgErrors && (valSigBg != 0)) valBgErr = valSigBgErr * (valBg / valSigBg);
		if (useIndependentSigBgErrors) {
		    valSigErr = sqrt((valSigBgErr * valSigBgErr) + (valBgErr * valBgErr));
		} else if (useAnticorrelatedSigBgErrors) {
		    valSigErr = valSigBgErr - valBgErr;
		}
		if (normalizeSigBgErrors) {
		    Float_t sumErr = (valBgErr + valSigErr);
		    if (sumErr != 0) {
    			Float_t valBgErr_new = valSigBgErr * (valBgErr / sumErr);
    			Float_t valSigErr_new = valSigBgErr * (valSigErr / sumErr);
			valBgErr = valBgErr_new;
			valSigErr = valSigErr_new;
		    }
		} else if (normalizeSigBgSquareErrors) {
		    Float_t sumSqErr = ((valBgErr*valBgErr) + (valSigErr*valSigErr));
		    if (sumSqErr != 0) {
    			Float_t valBgErr_new = valSigBgErr * (valBgErr*valBgErr / sumSqErr);
    			Float_t valSigErr_new = valSigBgErr * (valSigErr*valSigErr / sumSqErr);
			valBgErr = valBgErr_new;
			valSigErr = valSigErr_new;
		    }
		}
		areaSigBg += valSigBg;
		areaSigBgErr += valSigBgErr * valSigBgErr;
		areaBg += valBg;
		areaBgErr += valBgErr * valBgErr;
		area += valSig;
		areaErr += valSigErr * valSigErr;
		sumW += valSig;
		sumWX += valSig * xCenter;
		sumWW += valSig * valSig;
		sumWWX += valSig * valSig * xCenter;
		sumWXX += valSig * xCenter * xCenter;
		sumE += valSigErr;
		sumEE += valSigErr * valSigErr;
	}
	areaSigBgErr = sqrt(areaSigBgErr);
	areaErr = sqrt(areaErr);
	areaBgErr = sqrt(areaBgErr);

	Float_t numEntries = (sumEE != 0) ? (sumW * sumW / sumEE) : 0;
	mean = (sumW != 0) ? sumWX / sumW : 0;
	width = (sumW != 0) ? sqrt((sumWXX / sumW) - (mean * mean)) : 0;
	meanErr = (numEntries != 0) ? (width / sqrt(numEntries)) : 0;

//	Stat_t stats[4];
//	invResSig.GetStats(stats);
/*
	stats[1] = 0;
	for (Int_t mbin = 0;mbin < invSigSig.GetXaxis()->GetNbins() + 2;mbin++) {
		Float_t err = invSigSig.GetBinError(mbin);
		stats[1] += err * err; // This was wrong in ROOT <= 4.00.08
	}
*/
//	Float_t numberOfEntries = (stats[1] != 0) ? (stats[0] * stats[0] / stats[1]) : 0;
//	Float_t meanData = (stats[0] != 0) ? (stats[2] / stats[0]) : 0;
//	Float_t widthData = (stats[0]/*numberOfEntries*/ != 0) ? TMath::Sqrt(TMath::Abs((stats[3] / stats[0]/*numberOfEntries*/) - (meanData * meanData))) : 0;
//	Float_t meanDataErr = (numberOfEntries != 0) ? (widthData / TMath::Sqrt(numberOfEntries)) : 0;
//	Float_t widthDataErr = (numberOfEntries != 0) ?  (widthData / TMath::Sqrt(2.0 * numberOfEntries)) : 0;
}

TCanvas *findPeak(const distribution_type *distribution
		, const distribution_type *distributionBg
		, const distribution_type *distributionBgRandom
		, const distribution_type *distributionBgLowMass
		, const distribution_type *distributionNeutralMult
		, const distribution_type *distributionPeakShape
		, const distribution_type *distributionPeakShapeLowMass
		, const distribution_type *distributionPeakShapeEta
		, const distribution_type *distributionPeakShapeEtaLowMass
		, const distribution_type *distributionPeakShapeEtabg
		, const distribution_type *distributionPeakShapeEtabgLowMass
                , Int_t numBinsToShow, Int_t binToShow, TCanvas *canvas
		, const TAnalysisSettingsTrigger &analysisSettings
		, TAnalysisRawResultsTrigger &analysisResults
) {
	if (!distribution) return canvas;

	TBinParameters binParameters = distribution->getBinParameters();

	Float_t fitLeft = getParameterValue(analysisSettings.fitLeft, binParameters);
        Float_t fitRight = getParameterValue(analysisSettings.fitRight, binParameters);

	Float_t peakLeft = getParameterValue(analysisSettings.peakLeft, binParameters);
	Float_t peakRight = getParameterValue(analysisSettings.peakRight, binParameters);
	Float_t peakPosition_limitLow, peakPosition_limitHigh;
	Float_t peakPosition = getParameterValue(analysisSettings.peakPosition, binParameters, &peakPosition_limitLow, &peakPosition_limitHigh);
	Float_t peakWidth_limitLow, peakWidth_limitHigh;
	Float_t peakWidth = getParameterValue(analysisSettings.peakWidth, binParameters, &peakWidth_limitLow, &peakWidth_limitHigh);

	Float_t peakLeftEta = getParameterValue(analysisSettings.peakLeftEta, binParameters);
	Float_t peakRightEta = getParameterValue(analysisSettings.peakRightEta, binParameters);
	Float_t peakPositionEta_limitLow, peakPositionEta_limitHigh;
	Float_t peakPositionEta = getParameterValue(analysisSettings.peakPositionEta, binParameters, &peakPositionEta_limitLow, &peakPositionEta_limitHigh);
	Float_t peakWidthEta_limitLow, peakWidthEta_limitHigh;
	Float_t peakWidthEta = getParameterValue(analysisSettings.peakWidthEta, binParameters, &peakWidthEta_limitLow, &peakWidthEta_limitHigh);

	Int_t nrebin = getBestRebin(distribution->getDistribution(), analysisSettings.rebinNumBinsWeight, analysisSettings.rebinMeanErrorWeight, fitLeft, fitRight, 10);

	TH1F invSigBg(*distribution->getDistribution());
	invSigBg.SetNameTitle(analysisSettings.name + "_invSigBg_" + distribution->GetName(), distribution->GetTitle());
	if (nrebin > 1) invSigBg.Rebin(nrebin);
	if (analysisSettings.setErrorOnZeroBins) setErrorOnZeroBins(&invSigBg);

	TH1F invSig(invSigBg);
	invSig.SetName(analysisSettings.name + "_invSig_" + distribution->GetName());

	Float_t scaleBg = 1.0;
	TH1F *invBg = distributionBg ? (new TH1F(*distributionBg->getDistribution())) : 0;
	if (invBg) {
		invBg->SetName(analysisSettings.name + "_invBg_" + invBg->GetName());
		if (nrebin > 1) invBg->Rebin(nrebin);
		if (analysisSettings.smoothBg > 0) invBg->Smooth(analysisSettings.smoothBg, -1, -1);
		if (analysisSettings.setErrorOnZeroBinsBg) setErrorOnZeroBins(invBg);
		if (analysisSettings.mixNormFixed > 0) {
		    Float_t corr = analysisSettings.mixNormFixed;
		    if (analysisSettings.mixNormFixed > 0) scaleBg *= corr;
		}
		if (analysisSettings.mixNormMassRange) {
		    Float_t intBg = invBg->Integral(invBg->GetXaxis()->FindBin(analysisSettings.mixNormLeft), invBg->GetXaxis()->FindBin(analysisSettings.mixNormRight));
		    Float_t intSig = invSigBg.Integral(invSigBg.GetXaxis()->FindBin(analysisSettings.mixNormLeft), invSigBg.GetXaxis()->FindBin(analysisSettings.mixNormRight));
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormMassRange) scaleBg *= corr;
		}
		if (analysisSettings.mixNormEntries) {
		    Float_t intBg = invBg->Integral();
		    Float_t intSig = invSigBg.Integral();
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormEntries) scaleBg *= corr;
		}
		invBg->Scale(scaleBg);
		if (analysisSettings.subtractBg) {
		    cout << invSig.GetName() << " subtracting bg" << endl;
		    invSig.Add(invBg, -1);
		}
	}

	Float_t scaleBgRandom = 1.0;
	TH1F *invBgRandom = distributionBgRandom ? (new TH1F(*distributionBgRandom->getDistribution())) : 0;
	if (invBgRandom) {
		invBgRandom->SetName(analysisSettings.name + "_invBgRandom_" + invBgRandom->GetName());
		if (nrebin > 1) invBgRandom->Rebin(nrebin);
		if (analysisSettings.smoothBgRandom > 0) invBgRandom->Smooth(analysisSettings.smoothBgRandom, -1, -1);
		if (analysisSettings.setErrorOnZeroBinsBgRandom) setErrorOnZeroBins(invBgRandom);
		if (analysisSettings.mixNormFixedRandom > 0) {
		    Float_t corr = analysisSettings.mixNormFixedRandom;
		    if (analysisSettings.mixNormFixedRandom > 0) scaleBgRandom *= corr;
		}
		if (analysisSettings.mixNormMassRangeRandom) {
		    Float_t intBg = invBgRandom->Integral(invBgRandom->GetXaxis()->FindBin(analysisSettings.mixNormLeftRandom), invBgRandom->GetXaxis()->FindBin(analysisSettings.mixNormRightRandom));
		    Float_t intSig = invSigBg.Integral(invSigBg.GetXaxis()->FindBin(analysisSettings.mixNormLeftRandom), invSigBg.GetXaxis()->FindBin(analysisSettings.mixNormRightRandom));
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormMassRangeRandom) scaleBgRandom *= corr;
		}
		if (analysisSettings.mixNormEntriesRandom) {
		    Float_t intBg = invBgRandom->Integral();
		    Float_t intSig = invSigBg.Integral();
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormEntriesRandom) scaleBgRandom *= corr;
		}
		invBgRandom->Scale(scaleBgRandom);
		if (analysisSettings.subtractBgRandom) {
		    //cout << invSig.GetName() << " subtracting bg random" << endl;
		    invSig.Add(invBgRandom, -1);
		}
	}

	Float_t scaleBgLowMass = 1.0;
	TH1F *invBgLowMass = distributionBgLowMass ? (new TH1F(*distributionBgLowMass->getDistribution())) : 0;
	if (invBgLowMass) {
		invBgLowMass->SetName(analysisSettings.name + "_invBgLowMass_" + invBgLowMass->GetName());
		//cout << invBgLowMass->GetName() << " " << invBgLowMass->Integral();
		if (nrebin > 1) invBgLowMass->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsBgLowMass) setErrorOnZeroBins(invBgLowMass);
		if (analysisSettings.lowNormFixed > 0) {
		    Float_t corr = analysisSettings.lowNormFixed;
		    if (analysisSettings.lowNormFixed > 0) scaleBgLowMass *= corr;
		}
		if (analysisSettings.lowNormMassRange) {
		    Float_t intBg = invBgLowMass->Integral(invBgLowMass->GetXaxis()->FindBin(analysisSettings.lowNormMassRangeLeft), invBgLowMass->GetXaxis()->FindBin(analysisSettings.lowNormMassRangeRight));
		    Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(analysisSettings.lowNormMassRangeLeft), invSig.GetXaxis()->FindBin(analysisSettings.lowNormMassRangeRight));
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.lowNormMassRange) scaleBgLowMass *= corr;
		}
		if (analysisSettings.lowNormPoints) {
		    Float_t corr = analysisResults.lowNormPoints;
		    if (analysisSettings.lowNormPoints) scaleBgLowMass *= corr;
		}
		if (analysisSettings.lowNormPointsPt) {
		    Float_t corr = analysisResults.lowNormPointsPt;
		    if (analysisSettings.lowNormPointsPt) scaleBgLowMass *= corr;
		}
		//cout << " " << invBgLowMass->Integral() << " *" << scaleBgLowMass;
		invBgLowMass->Scale(scaleBgLowMass);
		//cout << " " << invBgLowMass->Integral() << endl;
		if (analysisSettings.subtractBgLowMass) {
		    //cout << invSig.GetName() << " subtracting bg low mass " << invSig.Integral();
		    invSig.Add(invBgLowMass, -1);
		    //cout << " " << invSig.Integral() << endl;
		}
	}

	Float_t scalePeakShape = 1.0;
	TH1F *invPeakShape = distributionPeakShape ? (new TH1F(*distributionPeakShape->getDistribution())) : 0;
	if (invPeakShape) {
		invPeakShape->SetName(analysisSettings.name + "_invPeakShape_" + invPeakShape->GetName());
		if (nrebin > 1) invPeakShape->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsPeakShape) setErrorOnZeroBins(invPeakShape);
		Float_t intBg = invPeakShape->Integral(invPeakShape->GetXaxis()->FindBin(peakLeft), invPeakShape->GetXaxis()->FindBin(peakRight));
		Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(peakLeft), invSig.GetXaxis()->FindBin(peakRight));
		scalePeakShape = (intBg != 0) ? (intSig / intBg) : 1.0;
		invPeakShape->Scale(scalePeakShape);
	}
	TH1F *invPeakShapeLowMass = distributionPeakShapeLowMass ? (new TH1F(*distributionPeakShapeLowMass->getDistribution())) : 0;
	if (invPeakShapeLowMass) {
		invPeakShapeLowMass->SetName(analysisSettings.name + "_invPeakShapeLowMass_" + invPeakShapeLowMass->GetName());
		if (nrebin > 1) invPeakShapeLowMass->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsPeakShape) setErrorOnZeroBins(invPeakShapeLowMass);
		invPeakShapeLowMass->Scale(scalePeakShape);
		Float_t scalePeakShapeLowMass = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeLowMass.begin();iter != analysisResults.scalePeakShapeLowMass.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeLowMass = (*iter).getValue();
		}
		invPeakShapeLowMass->Scale(scalePeakShapeLowMass);
	}
	if (invPeakShape && invPeakShapeLowMass && analysisSettings.subtractBgLowMassPeakShape) {
	    invPeakShape->Add(invPeakShapeLowMass, -1);
	}

	Float_t scalePeakShapeEta = 1.0;
	TH1F *invPeakShapeEta = distributionPeakShapeEta ? (new TH1F(*distributionPeakShapeEta->getDistribution())) : 0;
	if (invPeakShapeEta) {
		invPeakShapeEta->SetName(analysisSettings.name + "_invPeakShapeEta_" + invPeakShapeEta->GetName());
		if (nrebin > 1) invPeakShapeEta->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEta);
		Float_t intBg = invPeakShapeEta->Integral(invPeakShapeEta->GetXaxis()->FindBin(peakLeftEta), invPeakShapeEta->GetXaxis()->FindBin(peakRightEta));
		Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(peakLeftEta), invSig.GetXaxis()->FindBin(peakRightEta));
		scalePeakShapeEta = (intBg != 0) ? (intSig / intBg) : 1.0;
		invPeakShapeEta->Scale(scalePeakShapeEta);
	}
	TH1F *invPeakShapeEtaLowMass = distributionPeakShapeEtaLowMass ? (new TH1F(*distributionPeakShapeEtaLowMass->getDistribution())) : 0;
	if (invPeakShapeEtaLowMass) {
		invPeakShapeEtaLowMass->SetName(analysisSettings.name + "_invPeakShapeEtaLowMass_" + invPeakShapeEtaLowMass->GetName());
		if (nrebin > 1) invPeakShapeEtaLowMass->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtaLowMass);
		invPeakShapeEtaLowMass->Scale(scalePeakShapeEta);
		Float_t scalePeakShapeEtaLowMass = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeEtaLowMass.begin();iter != analysisResults.scalePeakShapeEtaLowMass.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeEtaLowMass = (*iter).getValue();
		}
		invPeakShapeEtaLowMass->Scale(scalePeakShapeEtaLowMass);
	}
	if (invPeakShapeEta && invPeakShapeEtaLowMass && analysisSettings.subtractBgLowMassPeakShapeEta) {
	    invPeakShapeEta->Add(invPeakShapeEtaLowMass, -1);
	}

	TH1F *invPeakShapeEtabg = distributionPeakShapeEtabg ? (new TH1F(*distributionPeakShapeEtabg->getDistribution())) : 0;
	if (invPeakShapeEtabg) {
		invPeakShapeEtabg->SetName(analysisSettings.name + "_invPeakShapeEtabg_" + invPeakShapeEtabg->GetName());
		if (nrebin > 1) invPeakShapeEtabg->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtabg);
		invPeakShapeEtabg->Scale(scalePeakShapeEta);
		if (analysisResults.scalePeakShapeEtabgRel > 0) invPeakShapeEtabg->Scale(analysisResults.scalePeakShapeEtabgRel);
	}
	TH1F *invPeakShapeEtabgLowMass = distributionPeakShapeEtabgLowMass ? (new TH1F(*distributionPeakShapeEtabgLowMass->getDistribution())) : 0;
	if (invPeakShapeEtabgLowMass) {
		invPeakShapeEtabgLowMass->SetName(analysisSettings.name + "_invPeakShapeEtabgLowMass_" + invPeakShapeEtabgLowMass->GetName());
		if (nrebin > 1) invPeakShapeEtabgLowMass->Rebin(nrebin);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtabgLowMass);
		invPeakShapeEtabgLowMass->Scale(scalePeakShapeEta);
		Float_t scalePeakShapeEtabgLowMass = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeEtabgLowMass.begin();iter != analysisResults.scalePeakShapeEtabgLowMass.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeEtabgLowMass = (*iter).getValue();
		}
		invPeakShapeEtabgLowMass->Scale(scalePeakShapeEtabgLowMass);
	}
	if (invPeakShapeEtabg && invPeakShapeEtabgLowMass && analysisSettings.subtractBgLowMassPeakShapeEta) {
	    invPeakShapeEtabg->Add(invPeakShapeEtabgLowMass, -1);
	}

	Float_t jetBgFraction_limitLow, jetBgFraction_limitHigh;
	Float_t jetBgFraction = getParameterValue(analysisSettings.jetBgFraction, binParameters, &jetBgFraction_limitLow, &jetBgFraction_limitHigh);
	if (jetBgFraction < 0.0) jetBgFraction = 0.0;
	if (jetBgFraction > 1.0) jetBgFraction = 1.0;
	if (jetBgFraction < jetBgFraction_limitLow) jetBgFraction_limitLow = jetBgFraction;
	if (jetBgFraction > jetBgFraction_limitHigh) jetBgFraction_limitHigh = jetBgFraction;
	histToFit1 = analysisSettings.fitBgDistribution ? invBg : 0;
	histToFit2 = analysisSettings.fitBgRandomDistribution ? invBgRandom : 0;
	histToFit3 = analysisSettings.fitBgLowMassDistribution ? invBgLowMass : 0;
	histToFit4 = analysisSettings.fitPeakShapeDistribution ? invPeakShape : 0;
	histToFit5 = analysisSettings.fitPeakShapeEtaDistribution ? invPeakShapeEta : 0;
	histToFit6 = analysisSettings.fitPeakShapeEtabgDistribution ? invPeakShapeEtabg : 0;
	TF1 fitDistributionFunc(analysisSettings.name + "_fitDistributionFunc", &getHistWeightedSum, analysisSettings.fitDistributionLeft, analysisSettings.fitDistributionRight, 6);
	for (Int_t ipar = 0;ipar < 6;ipar++) fitDistributionFunc.SetParameter(ipar, 0.0);
	Float_t bgJetToTotalRatio = 0.0;
	Float_t bgJetToTotalRatioErr = 0.0;
	TH1F *invFitDistribution = 0;
	if (histToFit1 || histToFit2 || histToFit3 || histToFit4 || histToFit5 || histToFit6) {
	    fitDistributionFunc.SetParameter(0, 1.0);
	    fitDistributionFunc.SetParLimits(0, 0.0, 100.0);
	    fitDistributionFunc.SetParameter(1, jetBgFraction);
    	    if ((jetBgFraction_limitLow != noLimitMagic) && (jetBgFraction_limitHigh != noLimitMagic)) fitDistributionFunc.SetParLimits(1, jetBgFraction_limitLow, jetBgFraction_limitHigh);
	    if (histToFit1) {if (!histToFit2) fitDistributionFunc.FixParameter(1, 1.0);} else fitDistributionFunc.FixParameter(0, 0.0);
	    if (histToFit3) {fitDistributionFunc.SetParameter(2, 1.0); fitDistributionFunc.SetParLimits(2, 0.0, 1000.0);} else fitDistributionFunc.FixParameter(2, 0.0);
	    if (histToFit4) {fitDistributionFunc.SetParameter(3, 1.0); fitDistributionFunc.SetParLimits(3, 0.0, 10.0);} else fitDistributionFunc.FixParameter(3, 0.0);
	    if (histToFit5) {fitDistributionFunc.SetParameter(4, 1.0); fitDistributionFunc.SetParLimits(4, 0.0, 10.0);} else fitDistributionFunc.FixParameter(4, 0.0);
	    if (histToFit6) fitDistributionFunc.FixParameter(5, (1.0 - trueEtaBranchingRatio) / trueEtaBranchingRatio); else fitDistributionFunc.FixParameter(5, 0.0);
	    invSigBg.Fit(&fitDistributionFunc, "RQN");
	    bgJetToTotalRatio = fitDistributionFunc.GetParameter(1);
	    bgJetToTotalRatioErr = fitDistributionFunc.GetParError(1);
	    invFitDistribution = new TH1F(invSigBg);
	    if (invFitDistribution) {
		invFitDistribution->SetName(analysisSettings.name + "_invFitDistribution_" + invSigBg.GetName());
		invFitDistribution->Reset();
		if (histToFit1) invFitDistribution->Add(histToFit1, fitDistributionFunc.GetParameter(0) * bgJetToTotalRatio);
		if (histToFit2) invFitDistribution->Add(histToFit2, fitDistributionFunc.GetParameter(0) * (1.0 - bgJetToTotalRatio));
		if (histToFit3) invFitDistribution->Add(histToFit3, fitDistributionFunc.GetParameter(2));
		if (histToFit6) invFitDistribution->Add(histToFit6, fitDistributionFunc.GetParameter(4) * fitDistributionFunc.GetParameter(5));
		if (analysisSettings.subtractBgFitDistribution) {
		    //cout << invSig.GetName() << " subtracting bg fit distr " << invSig.Integral();
		    invSig.Add(invFitDistribution, -1);
		    //cout << " " << invSig.Integral() << endl;
		}
	    }
	}
	histToFit1 = 0; histToFit2 = 0; histToFit3 = 0; histToFit4 = 0; histToFit5 = 0; histToFit6 = 0;

	Float_t fitParameter6_limitLow, fitParameter6_limitHigh;
	Float_t fitParameter6 = getParameterValue(analysisSettings.fitParameter6, binParameters, &fitParameter6_limitLow, &fitParameter6_limitHigh);
	Float_t fitParameter7_limitLow, fitParameter7_limitHigh;
	Float_t fitParameter7 = getParameterValue(analysisSettings.fitParameter7, binParameters, &fitParameter7_limitLow, &fitParameter7_limitHigh);
	Float_t fitParameters_limitLow, fitParameters_limitHigh;
	Float_t fitParameters = getParameterValue(analysisSettings.fitParameters, binParameters, &fitParameters_limitLow, &fitParameters_limitHigh);

	TF1 fitFunc(analysisSettings.name + "_fitFunc", analysisSettings.fitFunctionStr);
	fitFunc.SetRange(fitLeft, fitRight);
	fitFunc.SetParameter(0, invSig.Integral());
	fitFunc.SetParLimits(0, 0, invSig.Integral() * 2.0);
	fitFunc.SetParameter(1, peakPosition);
	if ((peakPosition_limitLow != noLimitMagic) && (peakPosition_limitHigh != noLimitMagic)) fitFunc.SetParLimits(1, peakPosition_limitLow, peakPosition_limitHigh);
	fitFunc.SetParameter(2, peakWidth);
	if ((peakWidth_limitLow != noLimitMagic) && (peakWidth_limitHigh != noLimitMagic)) fitFunc.SetParLimits(2, peakWidth_limitLow, peakWidth_limitHigh);
	fitFunc.SetParameter(3, invSig.Integral());
	fitFunc.SetParLimits(3, 0, invSig.Integral() * 2.0);
	fitFunc.SetParameter(4, peakPositionEta);
	if ((peakPositionEta_limitLow != noLimitMagic) && (peakPositionEta_limitHigh != noLimitMagic)) fitFunc.SetParLimits(4, peakPositionEta_limitLow, peakPositionEta_limitHigh);
	fitFunc.SetParameter(5, peakWidthEta);
	if ((peakWidthEta_limitLow != noLimitMagic) && (peakWidthEta_limitHigh != noLimitMagic)) fitFunc.SetParLimits(5, peakWidthEta_limitLow, peakWidthEta_limitHigh);
	fitFunc.SetParameter(6, fitParameter6);
	if ((fitParameter6_limitLow != noLimitMagic) && (fitParameter6_limitHigh != noLimitMagic)) fitFunc.SetParLimits(6, fitParameter6_limitLow, fitParameter6_limitHigh);
	fitFunc.SetParameter(7, fitParameter7);
	if ((fitParameter7_limitLow != noLimitMagic) && (fitParameter7_limitHigh != noLimitMagic)) fitFunc.SetParLimits(7, fitParameter7_limitLow, fitParameter7_limitHigh);
    	for (Int_t ipar = 8;ipar < fitFunc.GetNpar();ipar++) {
	    fitFunc.SetParameter(ipar, fitParameters);
	    if ((fitParameters_limitLow != noLimitMagic) && (fitParameters_limitHigh != noLimitMagic)) fitFunc.SetParLimits(ipar, fitParameters_limitLow, fitParameters_limitHigh);
	}
	invSig.Fit(&fitFunc, analysisSettings.fitOption);

	TF1 funcPurePeaks(fitFunc);
	funcPurePeaks.SetName(analysisSettings.name + "_funcPurePeaks_" + distribution->GetName());
	for (Int_t ipar = 6;ipar < funcPurePeaks.GetNpar();ipar++) funcPurePeaks.SetParameter(ipar, 0);

	TF1 funcPurePeak(fitFunc);
	funcPurePeak.SetName(analysisSettings.name + "_funcPurePeak_" + distribution->GetName());
	for (Int_t ipar = 6;ipar < funcPurePeak.GetNpar();ipar++) funcPurePeak.SetParameter(ipar, 0);
	funcPurePeak.SetParameter(3, 0);

	TF1 funcPurePeakEta(fitFunc);
	funcPurePeakEta.SetName(analysisSettings.name + "_funcPurePeakEta_" + distribution->GetName());
	for (Int_t ipar = 6;ipar < funcPurePeakEta.GetNpar();ipar++) funcPurePeakEta.SetParameter(ipar, 0);
	funcPurePeakEta.SetParameter(0, 0);

	TF1 funcPureBg(fitFunc);
	funcPureBg.SetName(analysisSettings.name + "_funcPureBg_" + distribution->GetName());
	funcPureBg.SetParameter(0, 0);
	funcPureBg.SetParameter(3, 0);

	TF1 funcBgPeak(fitFunc);
	funcBgPeak.SetName(analysisSettings.name + "_funcBgPeak_" + distribution->GetName());
	funcBgPeak.SetParameter(3, 0);

	TF1 funcBgPeakEta(fitFunc);
	funcBgPeakEta.SetName(analysisSettings.name + "_funcBgPeakEta_" + distribution->GetName());
	funcBgPeakEta.SetParameter(0, 0);

	TH1F invRes(invSig);
	invRes.SetName(analysisSettings.name + "_invRes_" + distribution->GetName());
	invRes.Add(&fitFunc, -1);

	TH1F invResBg(invSig);
	invResBg.SetName(analysisSettings.name + "_invResBg_" + distribution->GetName());
	invResBg.Add(&funcPurePeaks, -1);

	TH1F invResSig(invSig);
	invResSig.SetName(analysisSettings.name + "_invResSig_" + distribution->GetName());
	invResSig.Add(&funcPureBg, -1);

	//Float_t height = fitFunc.GetParameter(0);
	//Float_t heightErr = fitFunc.GetParError(0);
	Float_t mean = fitFunc.GetParameter(1);
	Float_t meanErr = fitFunc.GetParError(1);
	Float_t width = fitFunc.GetParameter(2);
	Float_t widthErr = fitFunc.GetParError(2);
	//Float_t heightEta = fitFunc.GetParameter(3);
	//Float_t heightEtaErr = fitFunc.GetParError(3);
	Float_t meanEta = fitFunc.GetParameter(4);
	Float_t meanEtaErr = fitFunc.GetParError(4);
	Float_t widthEta = fitFunc.GetParameter(5);
	Float_t widthEtaErr = fitFunc.GetParError(5);
	Float_t bgPar6 = fitFunc.GetParameter(6);
	Float_t bgPar6Err = fitFunc.GetParError(6);
	Float_t bgPar7 = fitFunc.GetParameter(7);
	Float_t bgPar7Err = fitFunc.GetParError(7);
	Float_t fitChi2 = fitFunc.GetChisquare();
	Float_t fitNdf = fitFunc.GetNDF();

	Float_t peakLeftYield = peakLeft;
	Float_t peakRightYield = peakRight;
	if (analysisSettings.usePeakRangeYield) {
	    peakLeftYield = getParameterValue(analysisSettings.peakLeftYield, binParameters);
	    peakRightYield = getParameterValue(analysisSettings.peakRightYield, binParameters);
	} else if (analysisSettings.usePeakRangeSigmaYield) {
	    Float_t leftSigma = getParameterValue(analysisSettings.peakLeftSigmaYield, binParameters);
	    Float_t rightSigma = getParameterValue(analysisSettings.peakRightSigmaYield, binParameters);
	    peakLeftYield = mean + (leftSigma * width);
	    peakRightYield = mean + (rightSigma * width);
	}

	Float_t peakLeftYieldEta = peakLeftEta;
	Float_t peakRightYieldEta = peakRightEta;
	if (analysisSettings.usePeakRangeYieldEta) {
	    peakLeftYieldEta = getParameterValue(analysisSettings.peakLeftYieldEta, binParameters);
	    peakRightYieldEta = getParameterValue(analysisSettings.peakRightYieldEta, binParameters);
	} else if (analysisSettings.usePeakRangeSigmaYieldEta) {
	    Float_t leftSigma = getParameterValue(analysisSettings.peakLeftSigmaYieldEta, binParameters);
	    Float_t rightSigma = getParameterValue(analysisSettings.peakRightSigmaYieldEta, binParameters);
	    peakLeftYieldEta = meanEta + (leftSigma * widthEta);
	    peakRightYieldEta = meanEta + (rightSigma * widthEta);
	}

	Float_t peakAreaShareBorder = 0;
	if (analysisSettings.usePeakRangeAreaShare) {
	    Int_t binLeft = invResSig.GetXaxis()->FindBin(fitLeft);
	    Int_t binRight = invResSig.GetXaxis()->FindBin(fitRight);
	    Float_t diff_min = -1;
	    Float_t borderMin = -1;
//cout << "binLeft = " << binLeft << ", binRight = " << binRight << endl;
	    for (Int_t ibin = binLeft;ibin < binRight;ibin++) {
		Float_t border = invResSig.GetXaxis()->GetBinUpEdge(ibin);
		Float_t intLeftPeak = funcPurePeak.Integral(border, fitRight);
		Float_t intRightPeak = funcPurePeakEta.Integral(fitLeft, border);
		Float_t diff = TMath::Abs(intLeftPeak - intRightPeak);
//cout << border << ":" << intLeftPeak << ":" << intRightPeak << ":" << diff << endl;
		if ((diff < diff_min) || (diff_min < 0)) {
		    diff_min = diff;
		    borderMin = border;
		}
	    }
	    if ((borderMin > fitLeft) && (borderMin < fitRight)) {
		peakRightYield = borderMin;
		peakLeftYieldEta = borderMin;
		peakAreaShareBorder = borderMin;
//cout << "borderMin = " << borderMin << endl;
	    }
	}

	if (peakLeftYieldEta < peakRightYield) peakLeftYieldEta = peakRightYield;

	Float_t areaFractionBg = getParameterValue(analysisSettings.areaFractionBg, binParameters);
	Float_t areaFractionBgErr = 0;
	if (areaFractionBg < 0) {
	    Float_t areaIntegralBg = funcPureBg.Integral(peakLeft, peakRight);
	    Float_t areaIntegralBgTotal = funcPureBg.Integral(fitLeft, fitRight);
	    areaFractionBg = (areaIntegralBgTotal != 0) ? ((areaIntegralBgTotal - areaIntegralBg) / areaIntegralBgTotal) : 0.0;
	}

	Float_t areaSigBg = 0;
	Float_t areaSigBgErr = 0;
	Float_t area = 0;
	Float_t areaErr = 0;
	Float_t areaBg = 0;
	Float_t areaBgErr = 0;
	Float_t meanData = 0;
	Float_t meanDataErr = 0;
	Float_t widthData = 0;
	Float_t widthDataErr = 0;

	Float_t areaSigBgEta = 0;
	Float_t areaSigBgEtaErr = 0;
	Float_t areaEta = 0;
	Float_t areaEtaErr = 0;
	Float_t areaBgEta = 0;
	Float_t areaBgEtaErr = 0;
	Float_t meanDataEta = 0;
	Float_t meanDataEtaErr = 0;
	Float_t widthDataEta = 0;
	Float_t widthDataEtaErr = 0;

	extractSignal(invSig, funcPureBg, peakLeftYield, peakRightYield
	    , analysisSettings.useProportionalBgErrors
	    , analysisSettings.useIndependentSigBgErrors
	    , analysisSettings.useAnticorrelatedSigBgErrors
	    , analysisSettings.normalizeSigBgErrors
	    , analysisSettings.normalizeSigBgSquareErrors
	    , areaSigBg
	    , areaSigBgErr
	    , area
	    , areaErr
	    , areaBg
	    , areaBgErr
	    , meanData
	    , meanDataErr
	    , widthData
	    , widthDataErr
	);

	extractSignal(invSig, funcPureBg, peakLeftYieldEta, peakRightYieldEta
	    , analysisSettings.useProportionalBgErrors
	    , analysisSettings.useIndependentSigBgErrors
	    , analysisSettings.useAnticorrelatedSigBgErrors
	    , analysisSettings.normalizeSigBgErrors
	    , analysisSettings.normalizeSigBgSquareErrors
	    , areaSigBgEta
	    , areaSigBgEtaErr
	    , areaEta
	    , areaEtaErr
	    , areaBgEta
	    , areaBgEtaErr
	    , meanDataEta
	    , meanDataEtaErr
	    , widthDataEta
	    , widthDataEtaErr
	);

	Float_t sigBg = getParameterValue(analysisSettings.signalToBackgroundRatio, binParameters);
	if (sigBg > 0) {
		Float_t area1 = areaSigBg / (1.0 + (1.0 / sigBg));
		areaErr = (area != 0) ? (areaErr * area1 / area) : 0;
		area = area1;
		Float_t areaBg1 = areaSigBg / (1.0 + sigBg);
		areaBgErr = (areaBg != 0) ? (areaBgErr * areaBg1 / areaBg) : 0;
		areaBg = areaBg1;
	}

	Float_t sigBgEta = getParameterValue(analysisSettings.signalToBackgroundRatioEta, binParameters);
	if (sigBgEta > 0) {
		Float_t areaEta1 = areaSigBgEta / (1.0 + (1.0 / sigBgEta));
		areaEtaErr = (areaEta != 0) ? (areaEtaErr * areaEta1 / areaEta) : 0;
		areaEta = areaEta1;
		Float_t areaBgEta1 = areaSigBgEta / (1.0 + sigBgEta);
		areaBgEtaErr = (areaBgEta != 0) ? (areaBgEtaErr * areaBgEta1 / areaBgEta) : 0;
		areaBgEta = areaBgEta1;
	}

	if (analysisSettings.normPeakShapeYieldShow && invPeakShape) {
		Float_t integral = invPeakShape->Integral();
		Float_t scalePeakShapeYield = (integral != 0) ? (area / integral) : 1.0;
		invPeakShape->Scale(scalePeakShapeYield);
	}

	if (analysisSettings.normPeakShapeEtaYieldShow && invPeakShapeEta) {
		Float_t integral = invPeakShapeEta->Integral();
		Float_t scalePeakShapeEtaYield = (integral != 0) ? (areaEta / integral) : 1.0;
		invPeakShapeEta->Scale(scalePeakShapeEtaYield);
		if (invPeakShapeEtabg) invPeakShapeEtabg->Scale(scalePeakShapeEtaYield);
	}

	TH1F invResPeakShapeBg(invSig);
	invResPeakShapeBg.SetName(analysisSettings.name + "_invResPeakShapeBg_" + distribution->GetName());
	if (invPeakShape) invResPeakShapeBg.Add(invPeakShape, -1);
	if (invPeakShapeEta) invResPeakShapeBg.Add(invPeakShapeEta, -1);

	if (analysisSettings.showBins) {
		if (!canvas) {
			canvas = new TCanvas(analysisSettings.name, analysisSettings.title);
			if (canvas) {
			    Int_t nyp = Int_t(sqrt(Float_t(numBinsToShow)));
			    Int_t nxp = Int_t((Float_t(numBinsToShow) - 0.0001) / Float_t(nyp)) + 1;
			    canvas->Divide(nxp, nyp, 0.001, 0.001);
			}
		}
		if (canvas) canvas->cd(binToShow + 1);
		Float_t max = invSigBg.GetMaximum();
		Float_t min = invSigBg.GetMinimum();
		if (!analysisSettings.showSigBg) {
		    invSigBg.Reset();
		}
		/*if (analysisSettings.showSigBg) */{
			invSigBg.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invSigBg.GetYaxis()->SetRangeUser(min - (TMath::Abs(max - min) * 0.1), max + (TMath::Abs(max - min) * 0.1));
			//invSigBg.SetMinimum(0);
			invSigBg.SetLineColor(analysisSettings.colorSigBg);
			invSigBg.SetMarkerColor(analysisSettings.colorSigBg);
			//invSigBg.SetFillColor(analysisSettings.colorSigBg);
			invSigBg.DrawCopy("");
		}
		if (invBg && analysisSettings.showBg) {
			invBg->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBg->SetLineColor(analysisSettings.colorBg);
			invBg->SetMarkerColor(analysisSettings.colorBg);
			//invBg->SetFillColor(analysisSettings.colorBg);
			invBg->DrawCopy("HIST SAME");
		}
		if (invBgRandom && analysisSettings.showBgRandom) {
			invBgRandom->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBgRandom->SetLineColor(analysisSettings.colorBgRandom);
			invBgRandom->SetMarkerColor(analysisSettings.colorBgRandom);
			//invBgRandom->SetFillColor(analysisSettings.colorBgRandom);
			invBgRandom->DrawCopy("HIST SAME");
		}
		if (invFitDistribution && analysisSettings.showBgFitDistribution) {
			invFitDistribution->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invFitDistribution->SetLineColor(analysisSettings.colorBgFitDistribution);
			invFitDistribution->DrawCopy("HIST SAME");
		}
		if (invBgLowMass && analysisSettings.showBgLowMass) {
			invBgLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBgLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invBgLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			invBgLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invBgLowMass->DrawCopy("HIST SAME");
		}
		if (analysisSettings.showSig) {
			invSig.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invSig.SetLineColor(analysisSettings.color);
			invSig.SetMarkerColor(analysisSettings.color);
			invSig.DrawCopy("SAME");
		}
		if (invPeakShape && analysisSettings.showPeakShape) {
			invPeakShape->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShape->SetLineColor(analysisSettings.color);
			invPeakShape->SetMarkerColor(analysisSettings.color);
			invPeakShape->DrawCopy("HIST SAME");
		}
		if (invPeakShapeLowMass && invPeakShape && analysisSettings.showBgLowMass && analysisSettings.showPeakShape) {
			invPeakShapeLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invPeakShapeLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invPeakShapeLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			invPeakShapeLowMass->DrawCopy("HIST SAME");
		}
		if (invPeakShapeEta && analysisSettings.showPeakShapeEta) {
			invPeakShapeEta->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEta->SetLineColor(analysisSettings.color);
			invPeakShapeEta->SetMarkerColor(analysisSettings.color);
			invPeakShapeEta->DrawCopy("HIST SAME");
		}
		if (invPeakShapeEtaLowMass && invPeakShapeEta && analysisSettings.showBgLowMass && analysisSettings.showPeakShapeEta) {
			invPeakShapeEtaLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEtaLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtaLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtaLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtaLowMass->DrawCopy("HIST SAME");
		}
		if (invPeakShapeEtabg && analysisSettings.showPeakShapeEtabg) {
			invPeakShapeEtabg->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEtabg->SetLineColor(analysisSettings.color);
			invPeakShapeEtabg->SetMarkerColor(analysisSettings.color);
			invPeakShapeEtabg->DrawCopy("HIST SAME");
		}
		if (invPeakShapeEtabgLowMass && invPeakShapeEtabg && analysisSettings.showBgLowMass && analysisSettings.showPeakShapeEtabg) {
			invPeakShapeEtabgLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEtabgLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtabgLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtabgLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtabgLowMass->DrawCopy("HIST SAME");
		}
		if (analysisSettings.showResidual) {
			invRes.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invRes.SetLineColor(analysisSettings.color);
			invRes.DrawCopy("HIST SAME");
		}
		if (analysisSettings.showResidualBg) {
			invResBg.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invResBg.SetLineColor(analysisSettings.color);
			invResBg.DrawCopy("HIST SAME");
		}
		if (analysisSettings.showResidualPeakShapeBg) {
			invResPeakShapeBg.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invResPeakShapeBg.SetLineColor(analysisSettings.color);
			invResPeakShapeBg.SetMarkerColor(analysisSettings.color);
			invResPeakShapeBg.DrawCopy("H SAME");
		}
		if (analysisSettings.showFuncPeak) {
			TArrow *arFitLeftPeak = new TArrow(peakLeft, 100000000, peakLeft, 0, 0.02, "-");
			arFitLeftPeak->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetLineStyle(2); arFitLeftPeak->Draw();
			TArrow *arFitLeftPeakYield = new TArrow(peakLeftYield, 100000000, peakLeftYield, 0, 0.05, "-");
			arFitLeftPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->Draw();
			TArrow *arFitRightPeak = new TArrow(peakRight, 100000000, peakRight, 0, 0.02, "-");
			arFitRightPeak->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeak->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeak->SetLineStyle(2); arFitRightPeak->Draw();
			TArrow *arFitRightPeakYield = new TArrow(peakRightYield, 100000000, peakRightYield, 0, 0.05, "-");
			arFitRightPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeakYield->Draw();
			funcBgPeak.SetLineColor(analysisSettings.color);
			funcBgPeak.DrawCopy("SAME");
		}
		if (analysisSettings.showFuncPeakEta) {
			TArrow *arFitLeftPeak = new TArrow(peakLeftEta, 100000000, peakLeftEta, 0, 0.02, "-");
			arFitLeftPeak->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetLineStyle(2); arFitLeftPeak->Draw();
			TArrow *arFitLeftPeakYield = new TArrow(peakLeftYieldEta, 100000000, peakLeftYieldEta, 0, 0.05, "-");
			arFitLeftPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->Draw();
			TArrow *arFitRightPeak = new TArrow(peakRightEta, 100000000, peakRightEta, 0, 0.02, "-");
			arFitRightPeak->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeak->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeak->SetLineStyle(2); arFitRightPeak->Draw();
			TArrow *arFitRightPeakYield = new TArrow(peakRightYieldEta, 100000000, peakRightYieldEta, 0, 0.05, "-");
			arFitRightPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeakYield->Draw();
			funcBgPeakEta.SetLineColor(analysisSettings.color);
			funcBgPeakEta.DrawCopy("SAME");
		}
		if (analysisSettings.showFuncBg) {
			TArrow *arFitRegionLeft = new TArrow(fitLeft, 100000000, fitLeft, 0, 0.06, "-");
			arFitRegionLeft->SetLineColor(analysisSettings.fitLinesColor); arFitRegionLeft->SetFillColor(analysisSettings.fitLinesColor); arFitRegionLeft->SetLineStyle(3); arFitRegionLeft->Draw();
			TArrow *arFitRegionRight = new TArrow(fitRight, 100000000, fitRight, 0, 0.06, "-");
			arFitRegionRight->SetLineColor(analysisSettings.fitLinesColor); arFitRegionRight->SetFillColor(analysisSettings.fitLinesColor); arFitRegionRight->SetLineStyle(3); arFitRegionRight->Draw();
			funcPureBg.SetLineColor(analysisSettings.bgColor);
			funcPureBg.SetFillColor(analysisSettings.fillColor);
			funcPureBg.SetFillStyle(analysisSettings.bgFillStyle);
			funcPureBg.DrawCopy("SAME");
		}
		if (analysisSettings.showTruePeakPos) {
			TArrow *arFitMinvPt = new TArrow(analysisSettings.peakPositionTrue, 100000000, analysisSettings.peakPositionTrue, 0, 0.03, "-");
			arFitMinvPt->SetLineColor(analysisSettings.color); arFitMinvPt->SetFillColor(analysisSettings.fillColor); arFitMinvPt->Draw();
			TArrow *arFitMinvPtEta = new TArrow(analysisSettings.peakPositionTrueEta, 100000000, analysisSettings.peakPositionTrueEta, 0, 0.03, "-");
			arFitMinvPtEta->SetLineColor(analysisSettings.color); arFitMinvPtEta->SetFillColor(analysisSettings.fillColor); arFitMinvPtEta->Draw();
		}
		if (analysisSettings.showText) {
			TString str = Form("#splitline{#splitline{S+B = %.3g +/- %.3g}{S = %.3g +/- %.3g}}{B = %.3g +/- %.3g}", areaSigBg, areaSigBgErr, area, areaErr, areaBg, areaBgErr);
		        TLatex *text = new TLatex((analysisSettings.binRangeLeft * 0.5) + (analysisSettings.binRangeRight * 0.5), invSigBg.GetMaximum() * 1.0, str.Data());
	                text->SetTextSize(0.05);
                        text->Draw();
		}
	}
	analysisResults.peakPosition.push_back(TBinStatistics(binParameters, mean, meanErr));
	analysisResults.peakWidth.push_back(TBinStatistics(binParameters, width, widthErr));
	analysisResults.peakPositionData.push_back(TBinStatistics(binParameters, meanData, meanDataErr));
	analysisResults.peakWidthData.push_back(TBinStatistics(binParameters, widthData, widthDataErr));
	analysisResults.rawYield.push_back(TBinStatistics(binParameters, area, areaErr));
	analysisResults.peakBackground.push_back(TBinStatistics(binParameters, areaBg, areaBgErr));
	analysisResults.peakSigBackground.push_back(TBinStatistics(binParameters, areaSigBg, areaSigBgErr));

	analysisResults.peakPositionEta.push_back(TBinStatistics(binParameters, meanEta, meanEtaErr));
	analysisResults.peakWidthEta.push_back(TBinStatistics(binParameters, widthEta, widthEtaErr));
	analysisResults.peakPositionDataEta.push_back(TBinStatistics(binParameters, meanDataEta, meanDataEtaErr));
	analysisResults.peakWidthDataEta.push_back(TBinStatistics(binParameters, widthDataEta, widthDataEtaErr));
	analysisResults.rawYieldEta.push_back(TBinStatistics(binParameters, areaEta, areaEtaErr));
	analysisResults.peakBackgroundEta.push_back(TBinStatistics(binParameters, areaBgEta, areaBgEtaErr));
	analysisResults.peakSigBackgroundEta.push_back(TBinStatistics(binParameters, areaSigBgEta, areaSigBgEtaErr));

	analysisResults.bgPar6.push_back(TBinStatistics(binParameters, bgPar6, bgPar6Err));
	analysisResults.bgPar7.push_back(TBinStatistics(binParameters, bgPar7, bgPar7Err));
	analysisResults.fitChi2Ndf.push_back(TBinStatistics(binParameters, fitChi2, fitNdf));
	analysisResults.bgJetToTotal.push_back(TBinStatistics(binParameters, bgJetToTotalRatio, bgJetToTotalRatioErr));
	analysisResults.areaFractionBg.push_back(TBinStatistics(binParameters, areaFractionBg, areaFractionBgErr));

	analysisResults.fitLeft.push_back(TBinStatistics(binParameters, fitLeft, 0));
	analysisResults.fitRight.push_back(TBinStatistics(binParameters, fitRight, 0));

	analysisResults.peakLeft.push_back(TBinStatistics(binParameters, peakLeft, 0));
	analysisResults.peakRight.push_back(TBinStatistics(binParameters, peakRight, 0));
	analysisResults.peakPositionForFit.push_back(TBinStatistics(binParameters, peakPosition, 0));
	analysisResults.peakWidthForFit.push_back(TBinStatistics(binParameters, peakWidth, 0));
	analysisResults.peakLeftYield.push_back(TBinStatistics(binParameters, peakLeftYield, 0));
	analysisResults.peakRightYield.push_back(TBinStatistics(binParameters, peakRightYield, 0));

	analysisResults.peakLeftEta.push_back(TBinStatistics(binParameters, peakLeftEta, 0));
	analysisResults.peakRightEta.push_back(TBinStatistics(binParameters, peakRightEta, 0));
	analysisResults.peakPositionForFitEta.push_back(TBinStatistics(binParameters, peakPositionEta, 0));
	analysisResults.peakWidthForFitEta.push_back(TBinStatistics(binParameters, peakWidthEta, 0));
	analysisResults.peakLeftYieldEta.push_back(TBinStatistics(binParameters, peakLeftYieldEta, 0));
	analysisResults.peakRightYieldEta.push_back(TBinStatistics(binParameters, peakRightYieldEta, 0));

	analysisResults.nrebin.push_back(TBinStatistics(binParameters, nrebin, 0));

	analysisResults.scaleBg.push_back(TBinStatistics(binParameters, scaleBg, 0));
	analysisResults.scaleBgRandom.push_back(TBinStatistics(binParameters, scaleBgRandom, 0));
	analysisResults.scaleBgLowMass.push_back(TBinStatistics(binParameters, scaleBgLowMass, 0));
	analysisResults.scalePeakShape.push_back(TBinStatistics(binParameters, scalePeakShape, 0));
	analysisResults.scalePeakShapeEta.push_back(TBinStatistics(binParameters, scalePeakShapeEta, 0));

	analysisResults.jetBgFractionForFit.push_back(TBinStatistics(binParameters, jetBgFraction, 0));

	analysisResults.fitParameter6.push_back(TBinStatistics(binParameters, fitParameter6, 0));
	analysisResults.fitParameter7.push_back(TBinStatistics(binParameters, fitParameter7, 0));

	analysisResults.peakAreaShareBorder.push_back(TBinStatistics(binParameters, peakAreaShareBorder, 0));

	if (invBg) delete invBg;
	if (invBgRandom) delete invBgRandom;
	if (invBgLowMass) delete invBgLowMass;
	if (invPeakShape) delete invPeakShape;
	if (invPeakShapeEta) delete invPeakShapeEta;

	return canvas;
}

TCanvas *findPeaks(const TAnalysisSettingsTrigger &analysisSettings
		 , TAnalysisRawResultsTrigger &analysisResults
		 , const TAnalysisRawResultsTrigger *analysisResultsSim
		 , const TAnalysisRawResultsTrigger *analysisResultsSimEta
		 , const TAnalysisRawResultsTrigger *analysisResultsSimEtabg
) {
	analysisResults.rawYield.clear();
	analysisResults.peakPosition.clear();
	analysisResults.peakWidth.clear();
	analysisResults.peakPositionData.clear();
	analysisResults.peakWidthData.clear();
	analysisResults.peakBackground.clear();
	analysisResults.peakSigBackground.clear();

	analysisResults.rawYieldEta.clear();
	analysisResults.peakPositionEta.clear();
	analysisResults.peakWidthEta.clear();
	analysisResults.peakPositionDataEta.clear();
	analysisResults.peakWidthDataEta.clear();
	analysisResults.peakBackgroundEta.clear();
	analysisResults.peakSigBackgroundEta.clear();

	analysisResults.bgPar6.clear();
	analysisResults.bgPar7.clear();
	analysisResults.fitChi2Ndf.clear();
	analysisResults.bgJetToTotal.clear();
	analysisResults.areaFractionBg.clear();

	analysisResults.fitLeft.clear();
	analysisResults.fitRight.clear();

	analysisResults.peakLeft.clear();
	analysisResults.peakRight.clear();
	analysisResults.peakPositionForFit.clear();
	analysisResults.peakWidthForFit.clear();
	analysisResults.peakLeftYield.clear();
	analysisResults.peakRightYield.clear();

	analysisResults.peakLeftEta.clear();
	analysisResults.peakRightEta.clear();
	analysisResults.peakPositionForFitEta.clear();
	analysisResults.peakWidthForFitEta.clear();
	analysisResults.peakLeftYieldEta.clear();
	analysisResults.peakRightYieldEta.clear();

	analysisResults.nrebin.clear();

	analysisResults.scaleBg.clear();
	analysisResults.scaleBgRandom.clear();
	analysisResults.scaleBgLowMass.clear();
	analysisResults.scalePeakShape.clear();
	analysisResults.scalePeakShapeEta.clear();

	analysisResults.jetBgFractionForFit.clear();

	analysisResults.fitParameter6.clear();
	analysisResults.fitParameter7.clear();

	analysisResults.peakAreaShareBorder.clear();

        Int_t binToShow = 0;
        TCanvas *canvas = 0;
        for (distribution_list_type::const_iterator iter = analysisResults.invlist.begin();iter != analysisResults.invlist.end();++iter) {
                const distribution_type &distr = *iter;
                const distribution_type *distrBg = 0;
                const distribution_type *distrBgRandom = 0;
                const distribution_type *distrBgLowMass = 0;
                const distribution_type *distrNeutralMult = 0;
                const distribution_type *distrPeakShape = 0;
                const distribution_type *distrPeakShapeLowMass = 0;
                const distribution_type *distrPeakShapeEta = 0;
                const distribution_type *distrPeakShapeEtaLowMass = 0;
                const distribution_type *distrPeakShapeEtabg = 0;
                const distribution_type *distrPeakShapeEtabgLowMass = 0;
                for (distribution_list_type::const_iterator iterBg = analysisResults.invBglist.begin();iterBg != analysisResults.invBglist.end();++iterBg) {
                                const distribution_type &distrBgRef = *iterBg;
                                if (distr == distrBgRef) distrBg = &distrBgRef;
                }
                for (distribution_list_type::const_iterator iterBgRandom = analysisResults.invBgRandomlist.begin();iterBgRandom != analysisResults.invBgRandomlist.end();++iterBgRandom) {
                                const distribution_type &distrBgRandomRef = *iterBgRandom;
                                if (distr == distrBgRandomRef) distrBgRandom = &distrBgRandomRef;
                }
                for (distribution_list_type::const_iterator iterBgLowMass = analysisResults.inv1gammalist.begin();iterBgLowMass != analysisResults.inv1gammalist.end();++iterBgLowMass) {
                                const distribution_type &distrBgLowMassRef = *iterBgLowMass;
                                if (distr == distrBgLowMassRef) distrBgLowMass = &distrBgLowMassRef;
                }
                for (distribution_list_type::const_iterator iterNeutralMult = analysisResults.invNeutralMultlist.begin();iterNeutralMult != analysisResults.invNeutralMultlist.end();++iterNeutralMult) {
                                const distribution_type &distrNeutralMultRef = *iterNeutralMult;
                                if (distr == distrNeutralMultRef) distrNeutralMult = &distrNeutralMultRef;
                }
                if (analysisResultsSim) for (distribution_list_type::const_iterator iterPeakShape = analysisResultsSim->invlist.begin();iterPeakShape != analysisResultsSim->invlist.end();++iterPeakShape) {
                                const distribution_type &distrPeakShapeRef = *iterPeakShape;
                                if (distr == distrPeakShapeRef) distrPeakShape = &distrPeakShapeRef;
                }
                if (analysisResultsSim) for (distribution_list_type::const_iterator iterPeakShapeLowMass = analysisResultsSim->inv1gammalist.begin();iterPeakShapeLowMass != analysisResultsSim->inv1gammalist.end();++iterPeakShapeLowMass) {
                                const distribution_type &distrPeakShapeLowMassRef = *iterPeakShapeLowMass;
                                if (distr == distrPeakShapeLowMassRef) distrPeakShapeLowMass = &distrPeakShapeLowMassRef;
                }
                if (analysisResultsSimEta) for (distribution_list_type::const_iterator iterPeakShapeEta = analysisResultsSimEta->invlist.begin();iterPeakShapeEta != analysisResultsSimEta->invlist.end();++iterPeakShapeEta) {
                                const distribution_type &distrPeakShapeEtaRef = *iterPeakShapeEta;
                                if (distr == distrPeakShapeEtaRef) distrPeakShapeEta = &distrPeakShapeEtaRef;
                }
                if (analysisResultsSimEta) for (distribution_list_type::const_iterator iterPeakShapeEtaLowMass = analysisResultsSimEta->inv1gammalist.begin();iterPeakShapeEtaLowMass != analysisResultsSimEta->inv1gammalist.end();++iterPeakShapeEtaLowMass) {
                                const distribution_type &distrPeakShapeEtaLowMassRef = *iterPeakShapeEtaLowMass;
                                if (distr == distrPeakShapeEtaLowMassRef) distrPeakShapeEtaLowMass = &distrPeakShapeEtaLowMassRef;
                }
                if (analysisResultsSimEtabg) for (distribution_list_type::const_iterator iterPeakShapeEtabg = analysisResultsSimEtabg->invlist.begin();iterPeakShapeEtabg != analysisResultsSimEtabg->invlist.end();++iterPeakShapeEtabg) {
                                const distribution_type &distrPeakShapeEtabgRef = *iterPeakShapeEtabg;
                                if (distr == distrPeakShapeEtabgRef) distrPeakShapeEtabg = &distrPeakShapeEtabgRef;
                }
                if (analysisResultsSimEtabg) for (distribution_list_type::const_iterator iterPeakShapeEtabgLowMass = analysisResultsSimEtabg->inv1gammalist.begin();iterPeakShapeEtabgLowMass != analysisResultsSimEtabg->inv1gammalist.end();++iterPeakShapeEtabgLowMass) {
                                const distribution_type &distrPeakShapeEtabgLowMassRef = *iterPeakShapeEtabgLowMass;
                                if (distr == distrPeakShapeEtabgLowMassRef) distrPeakShapeEtabgLowMass = &distrPeakShapeEtabgLowMassRef;
                }
                for (bin_stat_list_type::const_iterator iterLow = analysisResults.lowNormPointsPtCorr.begin();iterLow != analysisResults.lowNormPointsPtCorr.end();++iterLow) {
                                const bin_stat_type &binStatRef = *iterLow;
                                if (distr.getBinParameters() == binStatRef.getParameters()) analysisResults.lowNormPointsPt = binStatRef.getValue();
                }
                canvas = findPeak(&distr
				, distrBg
				, distrBgRandom
				, distrBgLowMass
				, distrNeutralMult
				, distrPeakShape
				, distrPeakShapeLowMass
				, distrPeakShapeEta
				, distrPeakShapeEtaLowMass
				, distrPeakShapeEtabg
				, distrPeakShapeEtabgLowMass
				, analysisResults.invlist.size(), binToShow, canvas, analysisSettings, analysisResults);
                binToShow++;
        }
	analysisResults.totalRawYield = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.rawYield.begin();iter != analysisResults.rawYield.end();++iter) analysisResults.totalRawYield += (*iter).getValue();
	analysisResults.totalRawYieldEta = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.rawYieldEta.begin();iter != analysisResults.rawYieldEta.end();++iter) analysisResults.totalRawYieldEta += (*iter).getValue();
        return canvas;
}


Bool_t fitParameter(const bin_stat_list_type *values, const TH1F *hist, TFitSettings &param, TFitSettings &paramSave, TString name) {
    Bool_t result = false;
    if (param.useFit) {
	param.valuesForFit.clear();
	param.fitSigma = 0;
	param.fitMeanAbs = 0;
	Int_t sizeMax = (values ? values->size() : 1) + (hist ? hist->GetXaxis()->GetNbins() : 1);
	Float_t *X = new Float_t[sizeMax];
	Float_t *Xerr = new Float_t[sizeMax];
	Float_t *Y = new Float_t[sizeMax];
	Float_t *Yerr = new Float_t[sizeMax];
	Int_t size = 0;
	Float_t maxX = -1;
	Float_t minX = -1;
	if (values) for (bin_stat_list_type::const_iterator iter = values->begin();iter != values->end();++iter) {
		const TBinStatistics &value = *iter;
		Float_t center = value.getParameters().getCenter();
		if ((center >= param.fitLow) && (center <= param.fitHigh)) {
			param.valuesForFit.push_back(value);
			X[size] = center;
			Xerr[size] = 0;
			Y[size] = value.getValue();
			Yerr[size] = value.getError();
			size++;
		}
		if ((center > maxX) || (maxX == -1)) maxX = center;
		if ((center < minX) || (minX == -1)) minX = center;
	}
	if (hist) {
	    for (Int_t ibin = 1;ibin <= hist->GetXaxis()->GetNbins();++ibin) {
		Float_t center = hist->GetXaxis()->GetBinCenter(ibin);
		if ((center >= param.fitLow) && (center <= param.fitHigh) && !((hist->GetBinContent(ibin) == 0) && (hist->GetBinError(ibin) == 0))) {
			TBinParameters binPar;
			binPar.min = hist->GetXaxis()->GetBinLowEdge(ibin);
			binPar.max = hist->GetXaxis()->GetBinUpEdge(ibin);
			TBinStatistics value;
			value.setParameters(binPar);
			value.setValue(hist->GetBinContent(ibin));
			value.setError(hist->GetBinError(ibin));
			param.valuesForFit.push_back(value);
			X[size] = center;
			Xerr[size] = 0;
			Y[size] = hist->GetBinContent(ibin);
			Yerr[size] = hist->GetBinError(ibin);
			size++;
		}
		if ((center > maxX) || (maxX == -1)) maxX = center;
		if ((center < minX) || (minX == -1)) minX = center;

	    }
	}
	if (size > 0) {
	    TGraphErrors *hgr = new TGraphErrors(size, X, Y, Xerr, Yerr);
	    if (hgr) {
		if (!param.func) param.func = new TF1(name, param.fitFuncStr);
		if (param.func) {
		    param.func->SetRange(minX, maxX);
		    hgr->Fit(param.func, param.fitOption);
		    param.wasFit = true;
		    Float_t sum = 0;
		    Float_t sumAbs = 0;
		    Float_t sum2 = 0;
		    for (Int_t i = 0;i < size;i++) {
			Float_t diff = Y[i] - param.func->Eval(X[i]);
			sum += diff;
			sumAbs += TMath::Abs(diff);
			sum2 += (diff * diff);
//cout << i <<  ": X = " << X[i] << ", Y = " << Y[i] << ", F = " << param.func->Eval(X[i]) << ", diff = " << diff << ", sum = " << sum << ", sumAbs = " << sumAbs << ", sum2 = " << sum2 << endl;
		    }
		    param.fitSigma = TMath::Sqrt((sum2 / size) - ((sum / size) * (sum / size)));
		    param.fitMeanAbs = sumAbs / size;
		    result = true;
		}
	    }
	}
    }
    paramSave = param;
    return result;
}

TCanvas *fitSpectraBins(const TAnalysisSettingsTrigger &analysisSettings
		      , TAnalysisRawResultsTrigger &analysisResults
		      , const TAnalysisRawResultsTrigger *analysisResultsSim
		      , const TAnalysisRawResultsTrigger *analysisResultsSimEta
		      , const TAnalysisRawResultsTrigger *analysisResultsSimEtabg
) {
	TAnalysisSettingsTrigger analysisSettingsNew = analysisSettings;

	if (analysisSettingsNew.correctCpv) {
	    analysisResults.cpvCorrection.clear();
	    TAnalysisSettingsTrigger analysisSettingsTemp = analysisSettingsNew;
	    TAnalysisRawResultsTrigger analysisResultsTemp = analysisResults;
	    analysisResultsTemp.invlist = analysisResults.trackdistlist;
	    analysisResultsTemp.invBglist.clear();
	    analysisResultsTemp.invBgRandomlist.clear();
	    analysisResultsTemp.inv1gammalist.clear();
	    //analysisResultsTemp.invBglist = analysisResults.trackdistBglist;
	    //analysisResultsTemp.invBgRandomlist = analysisResults.trackdistBgRandomlist;

	    analysisSettingsTemp.name += "_trackdist";
	    analysisSettingsTemp.title += ", CPV correction";
    	    analysisSettingsTemp.rebinNumBinsWeight = analysisSettingsTemp.rebinNumBinsWeightCpv;
    	    analysisSettingsTemp.rebinMeanErrorWeight = analysisSettingsTemp.rebinMeanErrorWeightCpv;
	    analysisSettingsTemp.setErrorOnZeroBins = false;
	    analysisSettingsTemp.setErrorOnZeroBinsBg = false;
	    analysisSettingsTemp.setErrorOnZeroBinsBgRandom = false;
	    analysisSettingsTemp.setErrorOnZeroBinsBgLowMass = false;

	    analysisSettingsTemp.fitLeft = analysisSettingsNew.fitLeftCpv;
	    analysisSettingsTemp.fitRight = analysisSettingsNew.fitRightCpv;

	    analysisSettingsTemp.peakLeft = analysisSettingsNew.peakLeftCpv;
	    analysisSettingsTemp.peakRight = analysisSettingsNew.peakRightCpv;
	    analysisSettingsTemp.peakPosition = analysisSettingsNew.peakPositionCpv;
	    analysisSettingsTemp.peakWidth = analysisSettingsNew.peakWidthCpv;

	    analysisSettingsTemp.peakLeftEta = TFitSettings();
	    analysisSettingsTemp.peakRightEta = TFitSettings();
	    analysisSettingsTemp.peakPositionEta = TFitSettings();
	    analysisSettingsTemp.peakWidthEta = TFitSettings();

	    analysisSettingsTemp.peakPositionTrue = 0;

	    analysisSettingsTemp.fitParameter6 = analysisSettingsNew.fitParameter6Cpv;
	    analysisSettingsTemp.fitParameter7 = analysisSettingsNew.fitParameter7Cpv;
	    analysisSettingsTemp.fitParameters = analysisSettingsNew.fitParametersCpv;

	    analysisSettingsTemp.fitFunctionStr = analysisSettingsNew.fitFunctionStrCpv;
	    analysisSettingsTemp.fitOption = analysisSettingsNew.fitOptionCpv;

	    analysisSettingsTemp.areaFractionBg = TFitSettings(-1);

	    analysisSettingsTemp.usePeakRangeYield = false;
	    analysisSettingsTemp.usePeakRangeYieldEta = false;
	    analysisSettingsTemp.usePeakRangeSigmaYield = false;
	    analysisSettingsTemp.usePeakRangeSigmaYieldEta = false;
	    analysisSettingsTemp.usePeakRangeAreaShare = false;

	    analysisSettingsTemp.signalToBackgroundRatio = TFitSettings(-1);
	    analysisSettingsTemp.signalToBackgroundRatioEta = TFitSettings(-1);

            analysisSettingsTemp.subtractBg = false;
            analysisSettingsTemp.subtractBgRandom = false;
            analysisSettingsTemp.subtractBgLowMass = false;
            analysisSettingsTemp.subtractBgFitDistribution = false;
            analysisSettingsTemp.fitBgDistribution = false;
            analysisSettingsTemp.fitBgRandomDistribution = false;
            analysisSettingsTemp.fitBgLowMassDistribution = false;
            analysisSettingsTemp.fitPeakShapeDistribution = false;
            analysisSettingsTemp.fitPeakShapeEtaDistribution = false;
            analysisSettingsTemp.fitPeakShapeEtabgDistribution = false;

	    analysisSettingsTemp.showBins = true;
	    analysisSettingsTemp.showSigBg = true;
	    analysisSettingsTemp.showSig = false;
	    analysisSettingsTemp.showPeakShape = false;
	    analysisSettingsTemp.showPeakShapeEta = false;
	    analysisSettingsTemp.showPeakShapeEtabg = false;
	    analysisSettingsTemp.showBg = false;
	    analysisSettingsTemp.showBgRandom = false;
	    analysisSettingsTemp.showBgLowMass = false;
	    analysisSettingsTemp.showBgFitDistribution = false;
	    analysisSettingsTemp.showFuncPeak = true;
	    analysisSettingsTemp.showFuncPeakEta = false;
	    analysisSettingsTemp.showFuncBg = true;
	    analysisSettingsTemp.showResidual = true;
	    analysisSettingsTemp.showResidualBg = true;
	    analysisSettingsTemp.showResidualPeakShapeBg = false;
	    analysisSettingsTemp.showText = false;
	    analysisSettingsTemp.showTruePeakPos = false;

	    analysisSettingsTemp.colorSigBg = analysisSettingsTemp.color;
	    analysisSettingsTemp.colorBg = analysisSettingsTemp.color;
	    analysisSettingsTemp.colorBgRandom = analysisSettingsTemp.color;
	    analysisSettingsTemp.colorBgFitDistribution = analysisSettingsTemp.color;
	    analysisSettingsTemp.fitLinesColor = analysisSettingsTemp.color;
	    analysisSettingsTemp.bgColor = analysisSettingsTemp.color;
	    analysisSettingsTemp.bgFillStyle = 0;

	    analysisSettingsTemp.binRangeLeft = analysisSettingsNew.binRangeLeftCpv;
	    analysisSettingsTemp.binRangeRight = analysisSettingsNew.binRangeRightCpv;

	    findPeaks(analysisSettingsTemp, analysisResultsTemp, 0, 0, 0);

	    if (
		fitParameter(&analysisResultsTemp.areaFractionBg, 0, analysisSettingsTemp.areaFractionBg, analysisResultsTemp.areaFractionBgParamFit, analysisSettingsNew.name + "_cpvFitFunc")
	    ) {
		analysisSettingsTemp.showBins = false;
		findPeaks(analysisSettingsTemp, analysisResultsTemp, 0, 0, 0);
	    }

	    analysisResults.cpvCorrection = analysisResultsTemp.areaFractionBg;
	    analysisResults.cpvCorrectionParamFit = analysisResultsTemp.areaFractionBgParamFit;
	    analysisResults.cpvCorrectionParam6 = analysisResultsTemp.bgPar6;
	    analysisResults.cpvCorrectionParam7 = analysisResultsTemp.bgPar7;
	}

	TCanvas *canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	analysisResults.bgJetToTotalFirst = analysisResults.bgJetToTotal;
	if (
	    fitParameter(&analysisResults.bgJetToTotal, 0, analysisSettingsNew.jetBgFraction, analysisResults.jetBgFractionParamFit, analysisSettingsNew.name + "_jetBgFractionFitFunc")
	    ) {
	    analysisSettingsNew.name += "_jetBgFractionFit";
	    analysisSettingsNew.title += ", jet bg fraction fit";
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	if (
	    fitParameter(&analysisResults.peakPosition, 0, analysisSettingsNew.peakPosition, analysisResults.peakPositionParamFit, analysisSettingsNew.name + "_peakPositionFitFunc")
	    ||
	    fitParameter(&analysisResults.peakWidth, 0, analysisSettingsNew.peakWidth, analysisResults.peakWidthParamFit, analysisSettingsNew.name + "_peakWidthFitFunc")
	    ||
	    fitParameter(&analysisResults.peakPositionEta, 0, analysisSettingsNew.peakPositionEta, analysisResults.peakPositionEtaParamFit, analysisSettingsNew.name + "_peakPositionEtaFitFunc")
	    ||
	    fitParameter(&analysisResults.peakWidthEta, 0, analysisSettingsNew.peakWidthEta, analysisResults.peakWidthEtaParamFit, analysisSettingsNew.name + "_peakWidthEtaFitFunc")
	    ) {
	    analysisSettingsNew.name += "_peakFit";
	    analysisSettingsNew.title += ", peak range fit";
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	if (
	    fitParameter(&analysisResults.bgPar6, 0, analysisSettingsNew.fitParameter6, analysisResults.bgPar6ParamFit, analysisSettingsNew.name + "_bgPar6FitFunc")
	    ||
	    fitParameter(&analysisResults.bgPar7, 0, analysisSettingsNew.fitParameter7, analysisResults.bgPar7ParamFit, analysisSettingsNew.name + "_bgPar7FitFunc")
	    ) {
	    analysisSettingsNew.name += "_bgParFit";
	    analysisSettingsNew.title += ", bg parameters fit";
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	if (
	    fitParameter(&analysisResults.peakSigBackground, 0, analysisSettingsNew.signalToBackgroundRatio, analysisResults.peakSigBackgroundParamFit, analysisSettingsNew.name + "_sigBgFitFunc")
	    ||
	    fitParameter(&analysisResults.peakSigBackgroundEta, 0, analysisSettingsNew.signalToBackgroundRatioEta, analysisResults.peakSigBackgroundEtaParamFit, analysisSettingsNew.name + "_sigBgEtaFitFunc")
	    ) {
	    analysisSettingsNew.name += "_sigBgFit";
	    analysisSettingsNew.title += ", sig/bg fit";
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	return canvas;
}

Float_t getScaleFactorFromPt(TString name, TString title, const TPointDataProcessor *proc1, const TPointDataProcessor *proc2
			, Float_t lowPt, Float_t highPt, const distribution_list_type &invlist, bin_stat_list_type &lowNormPointsPtCorr, Bool_t show, Bool_t print
) {
	Float_t scale = 0;
	lowNormPointsPtCorr.clear();
	if (proc1 && proc2) {
		const TH1F *histPt1 = proc1->getPt();
		const TH1F *histPt2 = proc2->getPt();
		if (histPt1 && histPt2) {
			Float_t int1 = histPt1->Integral(histPt1->GetXaxis()->FindBin(lowPt), histPt1->GetXaxis()->FindBin(highPt));
			Float_t int2 = histPt2->Integral(histPt2->GetXaxis()->FindBin(lowPt), histPt2->GetXaxis()->FindBin(highPt));
			if (int2 != 0) scale = int1 / int2;
			if (print) cout << "Scale factor calculation from Pt spectra, " << title.Data() << ": " << scale << endl;
    			for (distribution_list_type::const_iterator iter = invlist.begin();iter != invlist.end();++iter) {
            		    const distribution_type &distr = *iter;
			    const TBinParameters &param = distr.getBinParameters();
			    Float_t int1Pt = histPt1->Integral(histPt1->GetXaxis()->FindBin(param.min), histPt1->GetXaxis()->FindBin(param.max));
			    Float_t int2Pt = histPt2->Integral(histPt2->GetXaxis()->FindBin(param.min), histPt2->GetXaxis()->FindBin(param.max));
			    Float_t scalePt = (int2Pt != 0) ? (int1Pt / int2Pt) : 0.0;
			    lowNormPointsPtCorr.push_back(TBinStatistics(param, scalePt, 0));
			    //cout << param.min << "-" << param.max << ": " << int1Pt << "/" << int2Pt << " = " << scalePt << endl;
			}
			if (show) {
			    TCanvas *c = new TCanvas(name.Data(), title.Data());
			    if (c) {
				TH1F *hist1 = new TH1F(*histPt1);
				if (hist1) {
				    hist1->SetLineColor(kBlue);
				    hist1->SetMarkerColor(hist1->GetLineColor());
				    hist1->Draw();
				}
				TH1F *hist2 = new TH1F(*histPt2);
				if (hist2) {
				    hist2->Scale(scale);
				    hist2->SetLineColor(kRed);
				    hist2->SetMarkerColor(hist2->GetLineColor());
				    hist2->Draw("HIST SAME");
				}
				if (hist1) {
				    hist1->Draw("SAME");
				}
				TString text = Form("scale factor = %g", scale);
				TLatex *lat = new TLatex(0.5, 0.5, text.Data());
				if (lat) {
				    lat->Draw();
				}
				c->SetLogy();
			    }
			}
		}
	}
	return scale;
}
/*
const TMatrix *currentMatrix;
Double_t getMatrixConv(Double_t *x, Double_t *p) {
	if (!p || !x || !currentMatrix || (currentMatrix->GetNrows() != currentMatrix->GetNcols())) {
		cout << "Bad getMatrixConv call!!!" << endl;
		return 0;
	}
	Int_t ireco = ((Int_t)x[0]);
	Int_t size = currentMatrix->GetNrows();
	Double_t result = 0;
	for (Int_t isim = 0;isim < size;isim++) result += p[isim] * ((*currentMatrix)[isim][ireco]);
	return result;
}
*/
void calculateBinPurityMatrix(const distribution_list_type &simulatedPtDistributions
			, const bin_stat_list_type &areaIn
			, const bin_stat_list_type &areaInSim
			, const bin_stat_list_type &areaInSimReco
			, bin_stat_list_type &areaOut
			, bin_stat_list_type &areaOutR
) {
	if (!simulatedPtDistributions.size() || !areaInSim.size()) return;
	areaOut.clear();
	areaOutR.clear();
	Int_t matrixSize = areaIn.size();
	TMatrix matrix(matrixSize, matrixSize);
	TMatrix matrixErr(matrixSize, matrixSize);
	TMatrix vectorR(1, matrixSize);
	TMatrix vectorRErr(1, matrixSize);
	TMatrix vectorS(1, matrixSize);
	TMatrix vectorSErr(1, matrixSize);
	TMatrix vectorSSim(1, matrixSize);
	TMatrix vectorRSim(1, matrixSize);
	bin_stat_list_type::const_iterator iter = areaIn.begin();
	bin_stat_list_type::const_iterator iterSim = areaInSim.begin();
	for (Int_t icol = 0;icol < matrixSize;icol++, ++iter, ++iterSim) {
		vectorR[0][icol] = 0;
		vectorRErr[0][icol] = 0;
		vectorS[0][icol] = 0;
		vectorSErr[0][icol] = 0;
		vectorSSim[0][icol] = 0;
		vectorRSim[0][icol] = 0;
		const TBinStatistics &binStat = *iter;
		const TBinParameters &binPar = binStat.getParameters();
		for (bin_stat_list_type::const_iterator iterTmp = areaIn.begin();iterTmp != areaIn.end();++iterTmp) {
			const TBinStatistics &binStatTmp = *iterTmp;
			const TBinParameters &binParTmp = binStatTmp.getParameters();
			if (binParTmp == binPar) {
				vectorR[0][icol] = binStatTmp.getValue();
				vectorRErr[0][icol] = binStatTmp.getError();
//cout << "R[" << icol << "] = " << vectorR[0][icol] << ": " << binParTmp.min << " - " << binParTmp.max << endl;
			}
		}
		for (bin_stat_list_type::const_iterator iterTmp = areaInSimReco.begin();iterTmp != areaInSimReco.end();++iterTmp) {
			const TBinStatistics &binStatTmp = *iterTmp;
			const TBinParameters &binParTmp = binStatTmp.getParameters();
			if (binParTmp == binPar) {
				vectorRSim[0][icol] = binStatTmp.getValue();
//cout << "Rsim[" << icol << "] = " << vectorRSim[0][icol] << ": " << binParTmp.min << " - " << binParTmp.max << endl;
			}
		}
		for (bin_stat_list_type::const_iterator iterTmp = areaInSim.begin();iterTmp != areaInSim.end();++iterTmp) {
			const TBinStatistics &binStatTmp = *iterTmp;
			const TBinParameters &binParTmp = binStatTmp.getParameters();
			if (binParTmp == binPar) {
				vectorSSim[0][icol] = binStatTmp.getValue();
//cout << "Ssim[" << icol << "] = " << vectorSSim[0][icol] << ": " << binParTmp.min << " - " << binParTmp.max << endl;
			}
		}
		const TH1F *simPtHist = 0;
		for (distribution_list_type::const_iterator simPtIter = simulatedPtDistributions.begin();(simPtIter != simulatedPtDistributions.end()) && !simPtHist;++simPtIter) {
			const TInvariantMassDistribution &simPt = *simPtIter;
			const TBinParameters &simPtBinPar = simPt.getBinParameters();
			if (binPar == simPtBinPar) {
				simPtHist = simPt.getDistribution();
//cout << "simPt[" << icol << "]: " << simPtBinPar.min << " - " << simPtBinPar.max << endl;
			}
		}
		bin_stat_list_type::const_iterator iterRow = areaInSim.begin();
		for (Int_t irow = 0;irow < matrixSize;irow++, ++iterRow) {
			matrix[irow][icol] = 0;
			matrixErr[irow][icol] = 0;
			if (simPtHist) {
				const TBinParameters &binParRow = (*iterRow).getParameters();
				Float_t integralBin = simPtHist->Integral(simPtHist->GetXaxis()->FindBin(binParRow.min), simPtHist->GetXaxis()->FindBin(binParRow.max));
				matrix[irow][icol] = integralBin;
//cout << "\t[" << irow << "]: " << binParRow.min << " - " << binParRow.max << endl;
			}
		}
	}
	for (Int_t irow = 0;irow < matrixSize;irow++) {
		Float_t generated = vectorSSim[0][irow];
		if (generated != 0) for (Int_t icol = 0;icol < matrixSize;icol++) matrix[irow][icol] = matrix[irow][icol] / generated;
	}
/*
	for (Int_t irow = 0;irow < matrixSize;irow++) {
		Float_t rowIntegral = 0;
		for (Int_t icol = 0;icol < matrixSize;icol++) rowIntegral += matrix[irow][icol];
cout << "row " << irow << ": rowIntegral = " << rowIntegral << endl;
	}
*/
/*
	for (Int_t icol = 0;icol < matrixSize;icol++) {
		Float_t colIntegral = 0;
		for (Int_t irow = 0;irow < matrixSize;irow++) colIntegral += matrix[irow][icol];
		Float_t generated = vectorSSim[0][icol];
		Float_t reconstructed = vectorRSim[0][icol];
		for (Int_t irow = 0;irow < matrixSize;irow++) matrix[irow][icol] = matrix[irow][icol] / (colIntegral / (reconstructed / generated));
	}
*/
//	Float_t matrixIntegral = 0;
//	for (Int_t icol = 0;icol < matrixSize;icol++) for (Int_t irow = 0;irow < matrixSize;irow++) matrixIntegral += matrix[irow][icol];
//	Float_t det = matrix.Determinant();
//cout << "matrix determinant = " << det << endl;
//matrix.Print();
//matrixErr.Print();
//vectorR.Print();
//vectorRErr.Print();
//vectorSSim.Print();
//vectorRSim.Print();
//TMatrix vectorRSimCheck = vectorSSim * matrix;
//vectorRSimCheck.Print();
//	TMatrix matrixInvert(matrix);
//	matrixInvert *= (matrixSize * matrixSize) / (matrixIntegral);
//matrixInvert.Print();
//return;
//	Double_t detTmp = 0;
//	matrixInvert.Invert(&detTmp);
//	matrixInvert *= (matrixSize * matrixSize) / (matrixIntegral);
//cout << "Determinant = " << detTmp << endl;
//matrixInvert.Print();
//	vectorS = vectorR * matrixInvert;
//vectorS.Print();
/*
	currentMatrix = &matrix;
	TF1 f("f", &getMatrixConv, -0.5, matrixSize - 0.5, matrixSize);
	TH1F h("h", "h", matrixSize, -0.5, matrixSize - 0.5);
	for (Int_t i = 0;i < matrixSize;i++) {
		h.SetBinContent(i + 1, vectorR[0][i]);
		h.SetBinError(i + 1, vectorRErr[0][i]);
		f.SetParameter(i, vectorR[0][i]);
		f.SetParLimits(i, 0, (100 * vectorR[0][i]) + 1000);
	}
	h.Fit(&f, "RQ");
	h.Fit(&f, "RQ");
	h.Fit(&f, "R");
	currentMatrix = 0;
	iterSim = areaInSim.begin();
	for (Int_t i = 0;i < matrixSize;i++) {
		vectorS[0][i] = f.GetParameter(i);
		vectorSErr[0][i] = f.GetParError(i);
	}
*/
//matrix.Print();
//vectorR.Print();
	for (Int_t i = 0;i < matrixSize;i++) {
		if (matrix[i][i] != 0) {
			vectorS[0][i] = (vectorR[0][i] / matrix[i][i]) * (3.2 / 10.0);
		} else {
			vectorS[0][i] = 0;
		}
	}
	TMatrix vectorSnext = vectorS;
/*
{
TMatrix d = vectorS;
d.Transpose(d);
d.Print();
}
*/
	Int_t averageOver = 10;
	Int_t iterationsNumber = 50;
	if (iterationsNumber < averageOver) iterationsNumber = averageOver;
	TMatrix averageS(averageOver, matrixSize);
	Int_t averageIndex = 0;
	for (Int_t iIter = 0;iIter < iterationsNumber;iIter++) {
		Float_t lastResUpRatio = 0;
		Float_t firstResDownRatio = 0;
		for (Int_t i = 0;(i < matrixSize) && (firstResDownRatio == 0);i++) {
			Float_t resDown = 0;
			for (Int_t j = 0;j < i;j++) {
				Float_t tmp = vectorS[0][j] * matrix[j][i];
				resDown += tmp;
			}
			if ((resDown != 0) && (vectorR[0][i] != 0)) firstResDownRatio = resDown / vectorR[0][i];
		}
		for (Int_t i = 0;i < matrixSize;i++) {
			Float_t resDown = 0;
//cout << i << ": ";
			for (Int_t j = 0;j < i;j++) {
				Float_t tmp = vectorSnext[0][j] * matrix[j][i];
				resDown += tmp;
//cout << (Int_t)tmp << " ";
			}
			//if ((resDown == 0) && (firstResDownRatio != 0)) {
			//	resDown = firstResDownRatio * vectorR[0][i];
			//}
//cout << "| ";
			Float_t resUp = 0;
			for (Int_t j = i + 1;j < matrixSize;j++) {
				Float_t tmp = vectorS[0][j] * matrix[j][i];
				resUp += tmp;
//cout << (Int_t)tmp << " ";
			}
			if ((resUp == 0) && (lastResUpRatio != 0)) {
				resUp = lastResUpRatio * vectorR[0][i];
			}
			if (vectorR[0][i] != 0) lastResUpRatio = resUp / vectorR[0][i];
			Float_t denom = matrix[i][i];
			if (denom == 0) {
				for (Int_t j = 0;j < matrixSize;j++) denom += matrix[i][j];
			}
			if (denom != 0) {
				vectorSnext[0][i] = (vectorR[0][i] - resDown - resUp) / denom;
			} else {
				vectorSnext[0][i] = 0;
			}
//cout << "> " << (Int_t)res << " (" << Int_t(100.0 * res / vectorR[0][i]) << "%)" << endl;
		}

		//Float_t randomContrib = 1 - exp(-iIter / 100);
		for (Int_t i = 0;i < matrixSize;i++) {
			Float_t step = vectorSnext[0][i] - vectorS[0][i];
			if (vectorS[0][i] != 0) {
				Float_t stepRel = step / vectorS[0][i];
				Float_t newStepRel = stepRel * exp(-TMath::Abs(stepRel));
				//newStepRel = newStepRel * (((1 - randomContrib) * 1.0) + (randomContrib * gRandom->Uniform(0, 1)));
				Float_t newStep = newStepRel * vectorS[0][i];
				vectorSnext[0][i] = vectorS[0][i] + newStep;
//cout << i << ": stepRel " << stepRel << " -> " << newStepRel << endl;
			}
		}

		for (Int_t i = 1;i < matrixSize - 1;i++) {
			Float_t average = (vectorSnext[0][i - 1] + vectorSnext[0][i + 1]) / 2.0;
			Float_t diff = vectorSnext[0][i] - average;
			if (average != 0) {
				Float_t diffRel = diff / average;
				Float_t newDiffRel = diffRel * exp(-3.0*TMath::Abs(diffRel));
				//newDiffRel = newDiffRel * (((1 - randomContrib) * 1.0) + (randomContrib * gRandom->Uniform(0, 1)));
				Float_t newDiff = newDiffRel * average;
				vectorSnext[0][i] = average + newDiff;
//cout << i << ": diffRel " << diffRel << " -> " << newDiffRel << endl;
			}
		}
		Float_t maxStep = -1;
		for (Int_t i = 0;i < matrixSize;i++) {
			if (vectorS[0][i] != 0) {
				Float_t stepAbs = TMath::Abs((vectorSnext[0][i] - vectorS[0][i]) / vectorS[0][i]);
				if ((stepAbs > maxStep) || (maxStep == -1)) maxStep = stepAbs;
			}
		}
		Float_t chi2 = 0;
		Int_t ndf = 0;
		Float_t ratio = 0;
		for (Int_t i = 0;i < matrixSize - 0;i++) {
			if ((vectorR[0][i] != 0) && (vectorRErr[0][i] != 0)) {
				Float_t Rnew = 0;
				for (Int_t j = 0;j < matrixSize;j++) Rnew += vectorSnext[0][i] * matrix[j][i];
				chi2 += (vectorR[0][i] - Rnew) * (vectorR[0][i] - Rnew) / (vectorRErr[0][i] * vectorRErr[0][i]);
				ndf++;
				ratio += Rnew / vectorR[0][i];
//cout << Rnew / vectorR[0][i] << ", ";
			}
		}
//cout << endl;
		if (ndf) ratio /= ndf;
//cout << "Iteration " << iIter << ": Max step " << maxStep << ", chi2/ndf = " << chi2 << "/" << ndf << " = " << (chi2 / ndf) << ", ratio = " << ratio << endl;
		vectorS = vectorSnext;
		for (Int_t i = 0;i < matrixSize;i++) averageS[averageIndex][i] = vectorS[0][i];
		averageIndex++;
		if (averageIndex >= averageOver) averageIndex = 0;
/*
{
TMatrix d = vectorS;
d.Transpose(d);
d.Print();
}
*/
	}
	for (Int_t icol = 0;icol < matrixSize;icol++) {
		Float_t sum = 0;
		for (Int_t irow = 0;irow < averageOver;irow++) sum += averageS[irow][icol];
		sum /= (Float_t)averageOver;
		vectorS[0][icol] = sum;
	}
/*
{
TMatrix d = vectorS;
d.Transpose(d);
d.Print();
}
*/
	Float_t ratio = 0;
	Int_t ratioNdf = 0;
	for (Int_t i = 2;i < matrixSize - 2;i++) {
		Float_t Rnew = 0;
		for (Int_t j = 0;j < matrixSize;j++) Rnew += vectorS[0][i] * matrix[j][i];
		if (vectorR[0][i] != 0) {
			ratio += Rnew / vectorR[0][i];
			ratioNdf++;
		}
	}
	if (ratioNdf) {
		ratio /= ratioNdf;
//cout << "ratio = " << ratio << endl;
		for (Int_t i = 0;i < matrixSize;i++) {
			vectorS[0][i] = vectorS[0][i] / ratio;
			vectorSErr[0][i] = vectorSErr[0][i] / ratio;
		}
	}
	iterSim = areaInSim.begin();
	for (Int_t i = 0;i < matrixSize;i++, ++iterSim) {
		const TBinStatistics &binStat = *iterSim;
		//const TBinParameters &binPar = binStat.getParameters();
		TBinStatistics newBinStat = binStat;
		newBinStat.setValue(vectorS[0][i]);
		newBinStat.setError(vectorSErr[0][i]);
		areaOut.push_back(newBinStat);
		Float_t Rnew = 0;
		for (Int_t j = 0;j < matrixSize;j++) Rnew += vectorS[0][i] * matrix[j][i];
		TBinStatistics newRBin = binStat;
		newRBin.setValue(Rnew);
		newRBin.setError(0);
		areaOutR.push_back(newRBin);
	}
}

void calculateBinPurity(const TCandidateDataProcessor *candidate
			, const list<TInvariantMassDistribution> &inv
			, const Char_t *shiftFuncName, const Char_t *shiftFunc
			, TF1* *fPtShift
			, const Char_t *shiftCanvasName, const Char_t *shiftCanvasTitle
			, bool show
) {
	if (!candidate || !fPtShift) return;
	const TH2F *ptRes = candidate->getPtResolution();
	if (!ptRes) return;
	if (ptRes->GetEntries() <= 0) return;
	Int_t nx = ptRes->GetXaxis()->GetNbins();
	Int_t ny = ptRes->GetYaxis()->GetNbins();
	TString nameCopy = ptRes->GetName(); nameCopy += "_copy";
	TString titleCopy = ptRes->GetTitle();
	TH2F *ptResCopy = new TH2F(nameCopy.Data(), titleCopy.Data(), ny, ptRes->GetYaxis()->GetXmin(), ptRes->GetYaxis()->GetXmax(), nx, ptRes->GetXaxis()->GetXmin(), ptRes->GetXaxis()->GetXmax());
	ptResCopy->GetXaxis()->SetTitle(ptRes->GetYaxis()->GetTitle());
	ptResCopy->GetYaxis()->SetTitle(ptRes->GetXaxis()->GetTitle());
	Float_t *purity = new Float_t[inv.size()];
	Float_t *purityErr = new Float_t[inv.size()];
	Float_t *shiftX = new Float_t[inv.size()];
	Float_t *shiftXerr = new Float_t[inv.size()];
	Float_t *shiftY = new Float_t[inv.size()];
	Float_t *shiftYerr = new Float_t[inv.size()];
	Float_t minX = 100, maxX = 0, minY = 100, maxY = 0;
	Int_t shiftSize = 0;
	for (list<TInvariantMassDistribution>::const_iterator iter = inv.begin();iter != inv.end();++iter) {
		const TInvariantMassDistribution &binInv = *iter;
		const TBinParameters &binParams = binInv.getBinParameters();
		Float_t lowPt = binParams.min;
		Float_t highPt = binParams.max;
		Float_t sumW = 0;
		Float_t sumWTruepT = 0;
		Float_t sumWDetpT = 0;
		Float_t sumWTruepT2 = 0;
		Float_t sumWDetpT2 = 0;
		for (Int_t iy = 1;iy <= ny;++iy) {
			Axis_t biny = ptRes->GetYaxis()->GetBinCenter(iy);
			if ((biny >= lowPt) && (biny < highPt)) {
				for (Int_t ix = 1;ix <= nx;++ix) {
					Axis_t binx = ptRes->GetXaxis()->GetBinCenter(ix);
					Float_t w = ptRes->GetBinContent(ptRes->GetBin(ix, iy));
					sumW += w;
					sumWTruepT += (w * binx);
					sumWDetpT += (w * biny);
					sumWTruepT2 += (w * binx * binx);
					sumWDetpT2 += (w * biny * biny);
				}
			}
		}
		if (sumW > 0) {
			Float_t meanTruePt = sumWTruepT / sumW;
			Float_t meanTruePt2 = sumWTruepT2 / sumW;
			Float_t rmsTruePt2 = meanTruePt2 - (meanTruePt * meanTruePt);
			Float_t meanDetPt = sumWDetpT / sumW;
			Float_t meanDetPt2 = sumWDetpT2 / sumW;
			Float_t rmsDetPt2 = meanDetPt2 - (meanDetPt * meanDetPt);

			if (false) cout << "shiftSize " << shiftSize << ": x = " << meanTruePt << " +/- sqrt " << rmsTruePt2 << ", y = " << meanDetPt << " +/- sqrt " << rmsDetPt2 << " (y2 = " << meanDetPt2 << ")" << endl;

			minX = min(minX, meanDetPt);
			maxX = max(maxX, meanDetPt);
			minY = min(minY, meanTruePt);
			maxY = max(maxY, meanTruePt);
			shiftX[shiftSize] = meanDetPt;
			shiftXerr[shiftSize] = 0;//binParams.getWidth() / 2.0;//rmsDetPt2;
			shiftY[shiftSize] = meanTruePt;
			shiftYerr[shiftSize] = sqrt(rmsTruePt2);

			for (Int_t iy = 1;iy <= ny;++iy) {
				Axis_t biny = ptRes->GetYaxis()->GetBinCenter(iy);
				if ((biny >= lowPt) && (biny < highPt)) {
					for (Int_t ix = 1;ix <= nx;++ix) {
						//Axis_t binx = ptRes->GetXaxis()->GetBinCenter(ix);
						Float_t w = ptRes->GetBinContent(ptRes->GetBin(ix, iy));
						ptResCopy->SetBinContent(ptResCopy->GetBin(iy, ix), w / sumW);
					}
				}
			}
			
			Float_t correctedPtLow = lowPt * (meanTruePt / meanDetPt);
			Float_t correctedPtHigh = highPt * (meanTruePt / meanDetPt);
			Float_t sumWpure = 0;
			for (Int_t iy = 1;iy <= ny;++iy) {
				Axis_t biny = ptRes->GetYaxis()->GetBinCenter(iy);
				if ((biny >= lowPt) && (biny < highPt)) {
					for (Int_t ix = 1;ix <= nx;++ix) {
						Axis_t binx = ptRes->GetXaxis()->GetBinCenter(ix);
						if ((binx >= correctedPtLow) && (biny < correctedPtHigh)) {
							Float_t w = ptRes->GetBinContent(ptRes->GetBin(ix, iy));
							sumWpure += w;
						}
					}
				}
			}
			Float_t binPurity = sumWpure / sumW;
			purity[shiftSize] = binPurity;
			purityErr[shiftSize] = 0;
			shiftSize++;
		}
	}
	maxY = max(maxX, maxY);
	maxX = max(maxX, maxY);
	minY = min(minX, minY);
	minX = min(minX, minY);
	TGraphErrors *hgrShift = shiftSize ? (new TGraphErrors(shiftSize, shiftX, shiftY, shiftXerr, shiftYerr)) : 0;
	*fPtShift = new TF1(shiftFuncName, shiftFunc);
	(*fPtShift)->SetRange(minX, maxX);
	if (hgrShift) hgrShift->Fit(*fPtShift, "RQN");	
	cout << shiftFuncName << " = " << shiftFunc << endl;
	for (Int_t ipar = 0;ipar < (*fPtShift)->GetNpar();ipar++) {
		cout << shiftFuncName << ": [" << ipar << "] = " << (*fPtShift)->GetParameter(ipar) << endl;
	}
	if (show) {
		TCanvas *canvas = new TCanvas(shiftCanvasName, shiftCanvasTitle);
		if (canvas) {
			canvas->Divide(2, 1, 0.001, 0.001);
			canvas->cd(1);
			ptResCopy->Draw("COLZ");
			gPad->SetLogz(true);
			hgrShift->SetHistogram(ptResCopy);
			hgrShift->SetMarkerStyle(kFullCircle);
			hgrShift->Draw("P SAME");
			(*fPtShift)->SetLineColor(kRed);
			(*fPtShift)->SetLineWidth(3);
			(*fPtShift)->DrawCopy("SAME");
			TLine *l = new TLine(minX, minX, maxX, maxX);
			l->SetLineColor(kBlue);
			l->Draw();
			Float_t par1 = (*fPtShift)->GetParameter(1);
			TString textStr = Form("Slope = %1.2f", par1);
			TLatex *text = new TLatex((minX + maxX) / 2.0, maxY, textStr.Data());
			text->Draw();
			canvas->cd(2);
			TString dummyName(shiftCanvasName); dummyName += "_dummy";
			TString dummyTitle(shiftCanvasTitle); dummyTitle = "Bin purity;Corrected p_{T}, GeV/c;Bin purity";
			TH2F *dummy = new TH2F(dummyName.Data(), dummyTitle.Data(), 10, minX, maxX, 10, 0, 1);
			dummy->SetStats(0);
			dummy->Draw();
			TGraphErrors *hgrPurity = shiftSize ? (new TGraphErrors(shiftSize, shiftX, purity, shiftXerr, purityErr)) : 0;
			if (hgrPurity) {
				hgrPurity->SetHistogram(dummy);
				hgrPurity->SetLineColor(kBlack);
				hgrPurity->SetLineWidth(3);
				hgrPurity->Draw("L SAME");
			}
		}
	}
}

TF2 *ptBinPosFuncFromQCD = 0;
TGraph2D *ptBinPosGraphFromQCD = 0;
Double_t getPtBinPosFromQCD(Double_t *x, Double_t *p) {
	if (!x) return 0;
	//if (!ptBinPosGraphFromQCD) createPtBinPosFuncFromQCD();
	Float_t binCenter = x[0];
	if (!ptBinPosGraphFromQCD) return binCenter;
	Float_t binWidth = x[1];
	Float_t value = ptBinPosGraphFromQCD->Interpolate(binCenter, binWidth);
	//cout << "binCenter = " << binCenter << ", binWidth = " << binWidth << ", diffFrac = " << (value / binWidth) << endl;
	return binCenter - value;
}
TF1 *pQCDfunc = 0;
void clearPtBinPosFuncFromQCD() {
	if (ptBinPosGraphFromQCD) delete ptBinPosGraphFromQCD;
	ptBinPosGraphFromQCD = 0;
	if (ptBinPosFuncFromQCD) delete ptBinPosFuncFromQCD;
	ptBinPosFuncFromQCD = 0;
	if (pQCDfunc) delete pQCDfunc;
	pQCDfunc = 0;
}
void createPtBinPosFuncFromQCD(Double_t (* getNLOFunc) (Double_t*, Double_t*)) {
	if (!getNLOFunc) return;
	if (ptBinPosGraphFromQCD) return;
	//if (getNLOFunc) getNLOFunc(0, 0);
	pQCDfunc = new TF1("pQCDfuncDiv", getNLOFunc, 0, 20, 1);
	//pQCDfunc->SetParameter(0, 0); // scale mu = pT, inv. yield
	pQCDfunc->SetParameter(0, 3); // scale mu = pT, inv. yield * pT
	pQCDfunc->SetNpx(4);
	Float_t lowPtRange = 0.5;
	Float_t highPtRange = 18.0;
	Float_t stepPt = 0.5;
	Int_t numStepsPt = (Int_t)((highPtRange - lowPtRange) / stepPt);
	Float_t lowWidthRange = 0.5;
	Float_t highWidthRange = 2.51;
	Float_t stepWidth = 0.5;
	Int_t numStepsWidth = (Int_t)((highWidthRange - lowWidthRange) / stepWidth);
	ptBinPosGraphFromQCD = new TGraph2D();
	Int_t pointIndex = 0;
//	TH2F *hist = new TH2F("ptBinPosFuncFromQCDHist", "ptBinPosFuncFromQCDHist", numStepsPt + 1, lowPtRange, highPtRange, numStepsWidth + 1, lowWidthRange, highWidthRange);
	for (Float_t binWidth = lowWidthRange;binWidth < highWidthRange;binWidth += stepWidth) {
//		cout << "binWidth " << binWidth << endl;
		for (Float_t binCenter = lowPtRange;binCenter < highPtRange;binCenter += stepPt) {
			Float_t binLow = binCenter - (binWidth / 2.0);
			Float_t binHigh = binCenter + (binWidth / 2.0);
			Float_t integral = pQCDfunc->Integral(binLow, binHigh, (Double_t*)0, (Double_t)1e-6);
			Float_t mean = integral / (binHigh - binLow);
			Float_t meanX = pQCDfunc->GetX(mean, binLow, binHigh);
			//meanX = binLow;
			Float_t diff = binCenter - meanX;
			/*
			Float_t valOld = pQCDfunc->Eval(binCenter);
			Float_t val = pQCDfunc->Eval(meanX);
			Float_t diffFrac = diff / binWidth;
			Float_t corrVal = val / valOld;
			cout << "binCenter = " << binCenter << ", binWidth = " << binWidth << ", meanX = " << meanX << ", mean = " << mean << ", val = " << val << ", diffFrac = " << diffFrac << ", corr = " << corrVal << endl;
			*/
			//hist->Fill(binCenter, binWidth, diff);
			ptBinPosGraphFromQCD->SetPoint(pointIndex, binCenter, binWidth, diff);
			pointIndex++;
		}
	}
	//new TCanvas("1", "1");
	//hist->Draw("COLZ");
	//new TCanvas("2", "2");
	ptBinPosGraphFromQCD->SetNpx(numStepsPt);
	ptBinPosGraphFromQCD->SetNpy(numStepsWidth);
	//ptBinPosGraphFromQCD->Draw("COLZ");
	if (ptBinPosFuncFromQCD) delete ptBinPosFuncFromQCD;
	ptBinPosFuncFromQCD = 0;
	ptBinPosFuncFromQCD = new TF2("PtBinPosFromQCD", &getPtBinPosFromQCD, 0, 20, 0);
/*
	for (Float_t binCenter = lowPtRange;binCenter < highPtRange;binCenter += stepPt) {
		for (Float_t binWidth = lowWidthRange;binWidth < highWidthRange;binWidth += stepWidth) {
			Float_t trueVal = hist->GetBinContent(hist->FindBin(binCenter, binWidth));
//			Float_t approxVal = ptBinPosGraphFromQCD->Interpolate(6.5, 1.1);
			//Float_t approxVal = ptBinPosGraphFromQCD->Interpolate(binCenter, binWidth);
			Float_t approxVal = ptBinPosFuncFromQCD->Eval(binCenter, binWidth);
			Float_t diffFrac = (approxVal - trueVal) / trueVal;
			Float_t diffPtFrac = diffFrac * binWidth / binCenter;
			cout << "binCenter = " << binCenter << ", binWidth = " << binWidth
			 << ", diffFrac = " << diffFrac << ", diffPtFrac = " << diffPtFrac
			 << ", trueVal = " << trueVal << ", approxVal = " << approxVal << endl;
		}
	}
*/
}

void calculatePtShiftCorrInv(const list<TInvariantMassDistribution> &inv, const TF1* fPtShift, const TF2* fBinCenter, bin_stat_list_type &corr, bin_stat_list_type &centerShift) {
    for (list<TInvariantMassDistribution>::const_iterator iter = inv.begin();iter != inv.end();++iter) {
	const TInvariantMassDistribution &binInv = *iter;
	const TBinParameters &binParams = binInv.getBinParameters();
	Float_t lowPt = binParams.min;
	Float_t highPt = binParams.max;
	Float_t nominalCenter = (lowPt + highPt) / 2.0;
	Float_t newLowPt = fPtShift ? fPtShift->Eval(lowPt) : lowPt;
	Float_t newHighPt = fPtShift ? fPtShift->Eval(highPt) : highPt;
	Float_t newNominalCenter = (newLowPt + newHighPt) / 2.0;
	Float_t trueCenter = fBinCenter ? fBinCenter->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : newNominalCenter;
	if (pQCDfunc) {
	    TBinParameters newBinParams(binParams);
	    newBinParams.min = newLowPt;
	    newBinParams.max = newHighPt;
	    newBinParams.trueCenter = trueCenter;
	    TBinStatistics binCorr;
	    binCorr.setParameters(newBinParams);
	    Float_t valNominalCenter = pQCDfunc->Eval(nominalCenter);
	    Float_t valTrueCenter = pQCDfunc->Eval(trueCenter);
	    Float_t ratio = 1.0;
	    ratio *= (valTrueCenter != 0) ? (valNominalCenter / valTrueCenter) : 1.0;
	    ratio *= (trueCenter != 0) ? (nominalCenter / trueCenter) : 1.0;
	    binCorr.setValue(ratio);
	    binCorr.setError(0);
	    corr.push_back(binCorr);
	    TBinStatistics binShift;
	    binShift.setParameters(newBinParams);
	    Float_t fracShift = (nominalCenter != 0) ? ((trueCenter - nominalCenter) / nominalCenter) : 0;
	    binShift.setValue(fracShift);
	    binShift.setError(0);
	    centerShift.push_back(binShift);
	}
    }
}

void correctPtShiftInv(list<TInvariantMassDistribution> &inv, const TF1* fPtShift, const TF2* fBinCenter) {
	for (list<TInvariantMassDistribution>::iterator iter = inv.begin();iter != inv.end();++iter) {
		TInvariantMassDistribution &binInv = *iter;
		const TBinParameters &binParams = binInv.getBinParameters();
		Float_t lowPt = binParams.min;
		Float_t highPt = binParams.max;
		Float_t newLowPt = fPtShift ? fPtShift->Eval(lowPt) : lowPt;
		Float_t newHighPt = fPtShift ? fPtShift->Eval(highPt) : highPt;
    		Float_t newNominalCenter = (newLowPt + newHighPt) / 2.0;
		Float_t trueCenter = fBinCenter ? fBinCenter->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : newNominalCenter;
		TBinParameters newBinParams(binParams);
		newBinParams.min = newLowPt;
		newBinParams.max = newHighPt;
		newBinParams.trueCenter = trueCenter;
		binInv.setBinParameters(newBinParams);
	}
}

void correctPtShiftBins(bin_stat_list_type &bins, const TF1* fPtShift, const TF2* fBinCenter) {
	for (bin_stat_list_type::iterator iter = bins.begin();iter != bins.end();++iter) {
		TBinStatistics &bin = *iter;
		const TBinParameters &binParams = bin.getParameters();
		Float_t lowPt = binParams.min;
		Float_t highPt = binParams.max;
		Float_t newLowPt = fPtShift ? fPtShift->Eval(lowPt) : lowPt;
		Float_t newHighPt = fPtShift ? fPtShift->Eval(highPt) : highPt;
		Float_t newNominalCenter = (newLowPt + newHighPt) / 2.0;
		Float_t trueCenter = fBinCenter ? fBinCenter->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : newNominalCenter;
		TBinParameters newBinParams(binParams);
		newBinParams.min = newLowPt;
		newBinParams.max = newHighPt;
		newBinParams.trueCenter = trueCenter;
		bin.setParameters(newBinParams);
	}
}

void enumerateCombinations(UInt_t size, UInt_t N, const vector<UInt_t> cur, list<vector<UInt_t> > &result) {
//cout << "enum size = " << size << ", N = " << N << ", cur.size() = " << cur.size() << ", result.size() = " << result.size() << endl;
    if (cur.size() < N) {
	for (UInt_t i = 0;i < size;i++) {
	    Bool_t alreadyExists = false; for (UInt_t ii = 0;(ii < cur.size()) && (!alreadyExists);alreadyExists |= (cur[ii++] == i));
	    if (!alreadyExists) {
//cout << "i = " << i << endl;
		vector<UInt_t> newCur = cur;
		newCur.push_back(i);
		enumerateCombinations(size, N, newCur, result);
	    }
	}
    } else {
//cout << "including cur" << endl;
	result.push_back(cur);
    }
}

void selectInv(const TCandidateDataProcessor *proc, Float_t binSize, Float_t pTLow, Float_t pTHigh, Int_t mult, list<TInvariantMassDistribution> &binlist, const TPointDataProcessor *pointDataProcessor, bool debug) {
	const list<TInvariantMassDistribution> *listToUse = 0;
	if (mult == 0) {
	    listToUse = proc ? &proc->invariantMassDistributions : 0;
	} else if (mult == 1) {
	    listToUse = pointDataProcessor ? &pointDataProcessor->multiplicityPrimaryDistributions : &proc->multiplicityPrimaryDistributions;
	} else if (mult == 2) {
	    listToUse = pointDataProcessor ? &pointDataProcessor->multiplicityPointsDistributions : &proc->multiplicityPointsDistributions;
	} else if (mult == 4) {
	    listToUse = proc ? &proc->pointTrackDistDistributions : 0;
	} else if (mult == 5) {
	    listToUse = proc ? &proc->simulatedPtDistributions : 0;
	}
	if (listToUse) for (Float_t pTLowBin = pTLow, pTHighBin = pTLow + binSize;pTHighBin <= pTHigh;pTLowBin += binSize, pTHighBin += binSize) {
//if (debug) cout << "pTLowBin = " << pTLowBin << ", pTHighBin = " << pTHighBin << endl;
			vector<const TInvariantMassDistribution*> preselectedInv;
			Bool_t minExists = false;
			Bool_t maxExists = false;
			if (listToUse) for (list<TInvariantMassDistribution>::const_iterator iterInv = listToUse->begin();iterInv != listToUse->end();++iterInv) {
				if (((*iterInv).getBinParameters().min >= pTLowBin) && ((*iterInv).getBinParameters().max <= pTHighBin) && ((*iterInv).getBinParameters().min <= (*iterInv).getBinParameters().max)) {
					preselectedInv.push_back((const TInvariantMassDistribution*)&(*iterInv));
					if ((*iterInv).getBinParameters().min == pTLowBin) minExists = true;
					if ((*iterInv).getBinParameters().max == pTHighBin) maxExists = true;
				}
			}
			if (minExists && maxExists) {
			    list<const TInvariantMassDistribution*> selectedInv;
			    for (UInt_t nInv = 1;(nInv <= preselectedInv.size()) && (selectedInv.size() == 0);nInv++) {
				list<vector<UInt_t> > selIndVariants;
			        vector<UInt_t> selVariation;
			        enumerateCombinations(preselectedInv.size(), nInv, vector<UInt_t>(), selIndVariants);
			        for (list<vector<UInt_t> >::const_iterator iterVar = selIndVariants.begin();(iterVar != selIndVariants.end()) && (selVariation.size() == 0);++iterVar) {
				    const vector<UInt_t> &variation = *iterVar;
				    Bool_t allConnected = true;
				    for (UInt_t iVar = 0;iVar < variation.size();iVar++) {
					const TInvariantMassDistribution *iInv = preselectedInv[variation[iVar]];
					Bool_t thisConnectedLeft = false;
					Bool_t thisConnectedRight = false;
					if (iInv->getBinParameters().min == pTLowBin) thisConnectedLeft = true;
					if (iInv->getBinParameters().max == pTHighBin) thisConnectedRight = true;
					for (UInt_t jVar = 0;(jVar < variation.size()) && !(thisConnectedLeft && thisConnectedRight);jVar++) {
					    if (jVar != iVar) {
						const TInvariantMassDistribution *jInv = preselectedInv[variation[jVar]];
						if (iInv->getBinParameters().min == jInv->getBinParameters().max) thisConnectedLeft = true;
						if (iInv->getBinParameters().max == jInv->getBinParameters().min) thisConnectedRight = true;
					    }
					}
					if ((!thisConnectedLeft) || (!thisConnectedRight)) allConnected = false;
				    }
				    if (allConnected) {
					selVariation = variation;
				    }
				}
				if (selVariation.size() != 0) {
				    for (UInt_t iVar = 0;iVar < selVariation.size();iVar++) {
					const TInvariantMassDistribution *iInv = preselectedInv[selVariation[iVar]];
					selectedInv.push_back(iInv);
				    }
				}
			    }
            		    TBinParameters binPar;
            		    binPar.variable = (TBinVariable)0;
            		    binPar.min = pTLowBin;
            		    binPar.max = pTHighBin;
            		    TInvariantMassDistribution inv;
			    TString name;
			    TString title = Form("%g", pTLowBin); title += " - "; title += Form("%g", pTHighBin);
			    if (selectedInv.size() != 0) {
				for (list<const TInvariantMassDistribution*>::const_iterator iterInv = selectedInv.begin();iterInv != selectedInv.end();++iterInv) {
				    const TInvariantMassDistribution *iInv = *iterInv;
				    if (iterInv == selectedInv.begin()) {
					binPar.variable = iInv->getBinParameters().variable;
				    }
				    inv.add(*iInv, false);
//if (debug) cout << "adding " << iInv.getBinParameters().min << " - " << iInv.getBinParameters().max << ", integral " << iInv.getDistribution()->Integral() << endl;
				    name += iInv->GetName();
				}
			    }
            		    inv.setBinParameters(binPar);
			    inv.SetNameTitle(name, title);
//if (debug) cout << "name = " << name << ", title = " << title << ", integral " << inv.getDistribution()->Integral() << endl;
			    binlist.push_back(inv);
			}
		}
}

void selectInvSim(const TSimuDataProcessor *proc, Float_t binSize, Float_t pTLow, Float_t pTHigh, bin_stat_list_type &binlist) {
	const list<TBinStatistics> *listToUse = 0;
	if (proc) {
		listToUse = &proc->binStatistics;
	}
		for (Float_t pTLowBin = pTLow, pTHighBin = pTLow + binSize;pTHighBin <= pTHigh;pTLowBin += binSize, pTHighBin += binSize) {
			vector<const TBinStatistics*> preselectedInv;
			Bool_t minExists = false;
			Bool_t maxExists = false;
			if (listToUse) for (list<TBinStatistics>::const_iterator iterInv = listToUse->begin();iterInv != listToUse->end();++iterInv) {
				if (((*iterInv).getParameters().min >= pTLowBin) && ((*iterInv).getParameters().max <= pTHighBin) && ((*iterInv).getParameters().min <= (*iterInv).getParameters().max)) {
					preselectedInv.push_back((const TBinStatistics*)&(*iterInv));
					if ((*iterInv).getParameters().min == pTLowBin) minExists = true;
					if ((*iterInv).getParameters().max == pTHighBin) maxExists = true;
				}
			}
			if (minExists && maxExists) {
			    list<const TBinStatistics*> selectedInv;
			    for (UInt_t nInv = 1;(nInv <= preselectedInv.size()) && (selectedInv.size() == 0);nInv++) {
				list<vector<UInt_t> > selIndVariants;
			        vector<UInt_t> selVariation;
			        enumerateCombinations(preselectedInv.size(), nInv, vector<UInt_t>(), selIndVariants);
			        for (list<vector<UInt_t> >::const_iterator iterVar = selIndVariants.begin();(iterVar != selIndVariants.end()) && (selVariation.size() == 0);++iterVar) {
				    const vector<UInt_t> &variation = *iterVar;
				    Bool_t allConnected = true;
				    for (UInt_t iVar = 0;iVar < variation.size();iVar++) {
					const TBinStatistics *iInv = preselectedInv[variation[iVar]];
					Bool_t thisConnectedLeft = false;
					Bool_t thisConnectedRight = false;
					if (iInv->getParameters().min == pTLowBin) thisConnectedLeft = true;
					if (iInv->getParameters().max == pTHighBin) thisConnectedRight = true;
					for (UInt_t jVar = 0;(jVar < variation.size()) && !(thisConnectedLeft && thisConnectedRight);jVar++) {
					    if (jVar != iVar) {
						const TBinStatistics *jInv = preselectedInv[variation[jVar]];
						if (iInv->getParameters().min == jInv->getParameters().max) thisConnectedLeft = true;
						if (iInv->getParameters().max == jInv->getParameters().min) thisConnectedRight = true;
					    }
					}
					if ((!thisConnectedLeft) || (!thisConnectedRight)) allConnected = false;
				    }
				    if (allConnected) {
					selVariation = variation;
				    }
				}
				if (selVariation.size() != 0) {
				    for (UInt_t iVar = 0;iVar < selVariation.size();iVar++) {
					const TBinStatistics *iInv = preselectedInv[selVariation[iVar]];
					selectedInv.push_back(iInv);
				    }
				}
			    }
            		    TBinParameters binPar;
            		    binPar.variable = (TBinVariable)0;
            		    binPar.min = pTLowBin;
            		    binPar.max = pTHighBin;
            		    TBinStatistics inv;
			    TString name;
			    TString title = Form("%g", pTLowBin); title += " - "; title += Form("%g", pTHighBin);
			    if (selectedInv.size() != 0) {
				for (list<const TBinStatistics*>::const_iterator iterInv = selectedInv.begin();iterInv != selectedInv.end();++iterInv) {
				    const TBinStatistics *iInv = *iterInv;
				    if (iterInv == selectedInv.begin()) {
					binPar.variable = iInv->getParameters().variable;
				    }
				    inv.add(*iInv, false);
				    name += iInv->GetName();
				}
			    }
            		    inv.setParameters(binPar);
			    inv.SetNameTitle(name, title);
			    binlist.push_back(inv);
			}
		}
}

void calculatePSFromPoints(const Char_t *name, const Char_t *title
	, const TDataProcessorPool *poolMB, const TDataProcessorPool *poolHT1, const TDataProcessorPool *poolHT2
	, const TDataProcessorPool *pool1gammaMB, const TDataProcessorPool *pool1gammaHT1, const TDataProcessorPool *pool1gammaHT2
	, Bool_t chargedPoints, Bool_t allPoints
	, Float_t psPointsHT1Low, Float_t psPointsHT1High, Float_t psPointsHT2Low, Float_t psPointsHT2High, Float_t binSize
	, Float_t psPointsHT1Lowi, Float_t psPointsHT1Highi, Float_t psPointsHT2Lowi, Float_t psPointsHT2Highi
	, Float_t corrMBHT1, Float_t corrMBHT2, Float_t corrHT1HT2
	, Bool_t show, Bool_t print
	, Float_t &PSHT1, Float_t &PSHT2
	, Float_t &PSHT1Integral, Float_t &PSHT2Integral
) {
	GET_DATA(pointDataProcessorMB, const TPointDataProcessor, poolMB, (allPoints ? "pointAllMB" : (chargedPoints ? "pointChargedMB" : "pointMB")));
	GET_DATA(pointDataProcessorHT1, const TPointDataProcessor, poolHT1, (allPoints ? "pointAllHT1" : (chargedPoints ? "pointChargedHT1" : "pointHT1")));
	GET_DATA(pointDataProcessorHT2, const TPointDataProcessor, poolHT2, (allPoints ? "pointAllHT2" : (chargedPoints ? "pointChargedHT2" : "pointHT2")));

	GET_DATA(eventDataProcessorMB, const TEventDataProcessor, poolMB, "eventMB");
	GET_DATA(eventDataProcessorHT1, const TEventDataProcessor, poolHT1, "eventHT1");
	GET_DATA(eventDataProcessorHT2, const TEventDataProcessor, poolHT2, "eventHT2");

	GET_DATA(pointDataProcessorMB1gamma, const TPointDataProcessor, pool1gammaMB, (allPoints ? "pointAllMB" : (chargedPoints ? "pointChargedMB" : "pointMB")));
	GET_DATA(pointDataProcessorHT11gamma, const TPointDataProcessor, pool1gammaHT1, (allPoints ? "pointAllHT1" : (chargedPoints ? "pointChargedHT1" : "pointHT1")));
	GET_DATA(pointDataProcessorHT21gamma, const TPointDataProcessor, pool1gammaHT2, (allPoints ? "pointAllHT2" : (chargedPoints ? "pointChargedHT2" : "pointHT2")));

	GET_DATA(gammaDataProcessorMB, const TMCGammaDataProcessor, pool1gammaMB, "mcGammaMB");
	GET_DATA(gammaDataProcessorHT1, const TMCGammaDataProcessor, pool1gammaHT1, "mcGammaHT1");
	GET_DATA(gammaDataProcessorHT2, const TMCGammaDataProcessor, pool1gammaHT2, "mcGammaHT2");

	if ((corrMBHT1 < 0) || (corrMBHT1 > 1.0)) corrMBHT1 = 0;
	if ((corrMBHT2 < 0) || (corrMBHT2 > 1.0)) corrMBHT2 = 0;
	if ((corrHT1HT2 < 0) || (corrHT1HT2 > 1.0)) corrHT1HT2 = 0;

	const TH1F *histEnergyMBc = pointDataProcessorMB ? pointDataProcessorMB->getPt() : 0;
	const TH1F *histEnergyHT1c = pointDataProcessorHT1 ? pointDataProcessorHT1->getPt() : 0;
	const TH1F *histEnergyHT2c = pointDataProcessorHT2 ? pointDataProcessorHT2->getPt() : 0;

	const TH1F *histEnergySimMBc = gammaDataProcessorMB ? gammaDataProcessorMB->getPt() : 0;
	const TH1F *histEnergySimHT1c = gammaDataProcessorHT1 ? gammaDataProcessorHT1->getPt() : 0;
	const TH1F *histEnergySimHT2c = gammaDataProcessorHT2 ? gammaDataProcessorHT2->getPt() : 0;

	const TH1F *histEnergyRecoMBc = pointDataProcessorMB1gamma ? pointDataProcessorMB1gamma->getPt() : 0;
	const TH1F *histEnergyRecoHT1c = pointDataProcessorHT11gamma ? pointDataProcessorHT11gamma->getPt() : 0;
	const TH1F *histEnergyRecoHT2c = pointDataProcessorHT21gamma ? pointDataProcessorHT21gamma->getPt() : 0;

	TH1F *histEnergyMB = histEnergyMBc ? (new TH1F(*histEnergyMBc)) : 0;
	TH1F *histEnergyHT1 = histEnergyHT1c ? (new TH1F(*histEnergyHT1c)) : 0;
	TH1F *histEnergyHT2 = histEnergyHT2c ? (new TH1F(*histEnergyHT2c)) : 0;

	if (histEnergyMB && histEnergySimMBc && histEnergyRecoMBc) {
		histEnergyMB->Multiply(histEnergySimMBc);
		histEnergyMB->Divide(histEnergyRecoMBc);
	}
	if (histEnergyHT1 && histEnergySimHT1c && histEnergyRecoHT1c) {
		histEnergyHT1->Multiply(histEnergySimHT1c);
		histEnergyHT1->Divide(histEnergyRecoHT1c);
	}
	if (histEnergyHT2 && histEnergySimHT2c && histEnergyRecoHT2c) {
		histEnergyHT2->Multiply(histEnergySimHT2c);
		histEnergyHT2->Divide(histEnergyRecoHT2c);
	}

	if (histEnergyMB) if (histEnergyMB->GetSumw2N() == 0) histEnergyMB->Sumw2();
	if (histEnergyHT1) if (histEnergyHT1->GetSumw2N() == 0) histEnergyHT1->Sumw2();
	if (histEnergyHT2) if (histEnergyHT2->GetSumw2N() == 0) histEnergyHT2->Sumw2();

	if (histEnergyMB) while (histEnergyMB->GetBinWidth(1) < binSize) histEnergyMB->Rebin(Int_t((binSize / histEnergyMB->GetBinWidth(1)) + 1));
	if (histEnergyHT1) while (histEnergyHT1->GetBinWidth(1) < binSize) histEnergyHT1->Rebin(Int_t((binSize / histEnergyHT1->GetBinWidth(1)) + 1));
	if (histEnergyHT2) while (histEnergyHT2->GetBinWidth(1) < binSize) histEnergyHT2->Rebin(Int_t((binSize / histEnergyHT2->GetBinWidth(1)) + 1));

	if (histEnergyMB && eventDataProcessorMB && (eventDataProcessorMB->numPassedAllCuts != 0)) histEnergyMB->Scale(1.0 / eventDataProcessorMB->numPassedAllCuts);
	if (histEnergyHT1 && eventDataProcessorHT1 && (eventDataProcessorHT1->numPassedAllCuts != 0)) histEnergyHT1->Scale(1.0 / eventDataProcessorHT1->numPassedAllCuts);
	if (histEnergyHT2 && eventDataProcessorHT2 && (eventDataProcessorHT2->numPassedAllCuts != 0)) histEnergyHT2->Scale(1.0 / eventDataProcessorHT2->numPassedAllCuts);

	TH1F *histHT1r = histEnergyHT1 ? (new TH1F(*histEnergyHT1)) : 0;
	TH1F *histHT2r = histEnergyHT2 ? (new TH1F(*histEnergyHT2)) : 0;
	TH1F *histHT1i = histEnergyHT1 ? (new TH1F(*histEnergyHT1)) : 0;
	TH1F *histHT2i = histEnergyHT2 ? (new TH1F(*histEnergyHT2)) : 0;
	if (histHT1r) histHT1r->SetTitle("#frac{HighTower-1}{MinBias}, 100 * #frac{HighTower-2}{HighTower-1}");
	if (histHT1i) histHT1i->SetTitle("#frac{Int(HighTower-1)}{Int(MinBias)}, 100 * #frac{Int(HighTower-2)}{Int(HighTower-1)}");
	Float_t PSHT1r = 1.0, PSHT2r = 1.0, PSHT1i = 1.0, PSHT2i = 1.0;
	//TF1 *funcPSHT1r = new TF1("funcPSHT1r", "([0]*x*x) + ([1]*x) + [2]");
	//TF1 *funcPSHT2r = new TF1("funcPSHT2r", "([0]*x*x) + ([1]*x) + [2]");
	//TF1 *funcPSHT1i = new TF1("funcPSHT1i", "([0]*x*x) + ([1]*x) + [2]");
	//TF1 *funcPSHT2i = new TF1("funcPSHT2i", "([0]*x*x) + ([1]*x) + [2]");
	TF1 *funcPSHT1r = new TF1("funcPSHT1r", "[0]");
	TF1 *funcPSHT2r = new TF1("funcPSHT2r", "[0]");
	TF1 *funcPSHT1i = new TF1("funcPSHT1i", "[0]");
	TF1 *funcPSHT2i = new TF1("funcPSHT2i", "[0]");

	if (histHT1r && histEnergyMB) histHT1r->Divide(histEnergyMB);
	if (histHT2r && histEnergyHT1) histHT2r->Divide(histEnergyHT1);
	if (histHT2r) histHT2r->Scale(100.0);
	if (histHT1r && funcPSHT1r) {
		funcPSHT1r->SetRange(psPointsHT1Low, psPointsHT1High);
		histHT1r->Fit(funcPSHT1r, "QRN");
		PSHT1r = funcPSHT1r->GetParameter(0) / (1.0 + corrMBHT1);
		if (print) cout << "Prescales from points pT spectrum, " << title << ": PSHT1 = " << PSHT1r << endl;
		PSHT1 = PSHT1r;
	}
	if (histHT2r && funcPSHT2r) {
		funcPSHT2r->SetRange(psPointsHT2Low, psPointsHT2High);
		histHT2r->Fit(funcPSHT2r, "QRN");
		PSHT2r = PSHT1r * funcPSHT2r->GetParameter(0) / (100.0 * (1.0 + corrHT1HT2));
		if (print) cout << "Prescales from points pT spectrum, " << title << ": PSHT2 = " << PSHT2r << endl;
		PSHT2 = PSHT2r;
	}

	if (histEnergyMB && histEnergyHT1 && histEnergyHT2) {
		Float_t sumMB = 0, sumHT1 = 0, sumHT2 = 0;
		Float_t errMB = 0, errHT1 = 0, errHT2 = 0;
		for (Int_t mbin = histEnergyMB->GetXaxis()->GetNbins();mbin > 0;mbin--) {
			if (histEnergyMB) sumMB += histEnergyMB->GetBinContent(mbin);
			if (histEnergyHT1) sumHT1 += histEnergyHT1->GetBinContent(mbin);
			if (histEnergyHT2) sumHT2 += histEnergyHT2->GetBinContent(mbin);
			if (histEnergyMB) errMB += histEnergyMB->GetBinError(mbin) * histEnergyMB->GetBinError(mbin);
			if (histEnergyHT1) errHT1 += histEnergyHT1->GetBinError(mbin) * histEnergyHT1->GetBinError(mbin);
			if (histEnergyHT2) errHT2 += histEnergyHT2->GetBinError(mbin) * histEnergyHT2->GetBinError(mbin);
			Float_t relErrMB = (sumMB != 0) ? (errMB/(sumMB * sumMB)) : 0;
			Float_t relErrHT1 = (sumHT1 != 0) ? (errHT1/(sumHT1 * sumHT1)) : 0;
			Float_t relErrHT2 = (sumHT2 != 0) ? (errHT2/(sumHT2 * sumHT2)) : 0;
			Float_t relErrHT1i = sqrt(relErrMB + relErrHT1);
			Float_t relErrHT2i = sqrt(relErrHT1 + relErrHT2);
			if (histHT1i) {
				if ((sumMB != 0) && (sumHT1 != 0)) {
					Float_t integral = sumHT1 / sumMB;
					Float_t integralErr = integral * relErrHT1i;
					histHT1i->SetBinContent(mbin, integral);
					histHT1i->SetBinError(mbin, integralErr);
//cout << "bin " << mbin << ": sumMB=" << sumMB << ", sumHT1=" << sumHT1 << ", relErrMB=" << relErrMB << ", relErrHT1=" << relErrHT1 << ", relErrHT1i=" << relErrHT1i << ", integral=" << integral << ", err=" << integralErr<< endl;
				} else {
					histHT1i->SetBinContent(mbin, 0);
					histHT1i->SetBinError(mbin, 0);
				}
			}
			if (histHT2i) {
				if ((sumHT1 != 0) && (sumHT2 != 0)) {
					histHT2i->SetBinContent(mbin, sumHT2 / sumHT1);
					histHT2i->SetBinError(mbin, (sumHT2 / sumHT1) * relErrHT2i);
				} else {
					histHT2i->SetBinContent(mbin, 0);
					histHT2i->SetBinError(mbin, 0);
				}
			}
		}
		if (histHT2i) histHT2i->Scale(100.0);
		if (funcPSHT1i) {
			funcPSHT1i->SetRange(psPointsHT1Lowi, psPointsHT1Highi);
			if (histHT1i) histHT1i->Fit(funcPSHT1i, "QRN");
			PSHT1i = funcPSHT1i->GetParameter(0) / (1.0 + corrMBHT1);
			if (print) cout << "Prescales from points pT integral, " << title << ": PSHT1 = " << PSHT1i << endl;
			PSHT1Integral = PSHT1i;
		}
		if (funcPSHT2i) {
			funcPSHT2i->SetRange(psPointsHT2Lowi, psPointsHT2Highi);
			if (histHT2i) histHT2i->Fit(funcPSHT2i, "QRN");
			PSHT2i = PSHT1i * funcPSHT2i->GetParameter(0) / (100.0 * (1.0 + corrHT1HT2));
			if (print) cout << "Prescales from points pT integral, " << title << ": PSHT2 = " << PSHT2i << endl;
			PSHT2Integral = PSHT2i;
		}
	}
	if (show) {
		TString PSname = TString(name) + "_prescales_points";
		TString PStitle = TString(title) + ", Prescales calculation from points";
		TCanvas *psCanvas = new TCanvas(PSname.Data(), PStitle.Data());
		if (psCanvas) {
			psCanvas->Divide(1, 3, 0.001, 0.001);
			psCanvas->cd(1);
			if (histEnergyMB) {
				histEnergyMB->SetLineColor(kBlack);
				histEnergyMB->Draw();
			}
			if (histEnergyHT1) {
				histEnergyHT1->SetLineColor(kBlue);
				histEnergyHT1->Draw("SAME");
			}
			if (histEnergyHT2) {
				histEnergyHT2->SetLineColor(kRed);
				histEnergyHT2->Draw("SAME");
			}
			gPad->SetLogy();
			psCanvas->cd(2);
			if (histHT1r) {
				histHT1r->SetLineColor(kBlue);
				histHT1r->Draw();
			}
			if (histHT2r) {
				histHT2r->SetLineColor(kRed);
				histHT2r->Draw("SAME");
			}
			if (funcPSHT1r) {
				funcPSHT1r->SetLineColor(kBlue);
				funcPSHT1r->DrawCopy("SAME");
			}
			if (funcPSHT2r) {
				funcPSHT2r->SetLineColor(kRed);
				funcPSHT2r->DrawCopy("SAME");
			}
			Float_t x = 1;//histHT1r ? histHT1r->GetXaxis()->GetBinCenter(Int_t(histHT1r->GetXaxis()->GetNbins() / 1.7)) : 10;
			Float_t y = histHT1r ? (histHT1r->GetMaximum() * 0.6) : 400;
			TString textStr;
			textStr = Form("#splitline{PS_{HT1} = %4.0f}{PS_{HT2} = %4.0f}", PSHT1r, PSHT2r);
			TLatex *text = new TLatex(x, y, textStr.Data());
			text->Draw();
			psCanvas->cd(3);
			if (histHT1i) {
				histHT1i->SetLineColor(kBlue);
				histHT1i->Draw();
			}
			if (histHT2i) {
				histHT2i->SetLineColor(kRed);
				histHT2i->Draw("SAME");
			}
			if (funcPSHT1i) {
				funcPSHT1i->SetLineColor(kBlue);
				funcPSHT1i->DrawCopy("SAME");
			}
			if (funcPSHT2i) {
				funcPSHT2i->SetLineColor(kRed);
				funcPSHT2i->DrawCopy("SAME");
			}
			x = 1;//histHT1i ? histHT1i->GetXaxis()->GetBinCenter(Int_t(histHT1i->GetXaxis()->GetNbins() / 1.7)) : 10;
			y = histHT1i ? (histHT1i->GetMaximum() * 0.6) : 400;
			textStr = Form("#splitline{PS_{HT1} = %4.0f}{PS_{HT2} = %4.0f}", PSHT1i, PSHT2i);
			text = new TLatex(x, y, textStr.Data());
			text->Draw();
		}
	}
}

void calculatePSFromSim(const Char_t *name, const Char_t *title, const TDataProcessorPool *poolMB, const TDataProcessorPool *poolHT1
	, Float_t corrMBHT1, Float_t corrMBHT2, Float_t corrHT1HT2
	, Bool_t showReal, Bool_t print
	, Float_t &PSHT1, Float_t &PSHT2
) {

	GET_DATA(eventDataProcessorPSMB, const TEventDataProcessor, poolMB, "eventPSMB");
	GET_DATA(eventDataProcessorPSMBHT1, const TEventDataProcessor, poolMB, "eventPSMBHT1");
	GET_DATA(eventDataProcessorPSMBHT2, const TEventDataProcessor, poolMB, "eventPSMBHT2");
	GET_DATA(eventDataProcessorPSHT1, const TEventDataProcessor, poolHT1, "eventPSHT1");
	GET_DATA(eventDataProcessorPSHT1HT2, const TEventDataProcessor, poolHT1, "eventPSHT1HT2");

	const TH1F *maxAdcMBc = eventDataProcessorPSMB ? eventDataProcessorPSMB->getHighestAdcFinal() : 0;
	const TH1F *maxAdcMBHT1c = eventDataProcessorPSMBHT1 ? eventDataProcessorPSMBHT1->getHighestAdcFinal() : 0;
	const TH1F *maxAdcMBHT2c = eventDataProcessorPSMBHT2 ? eventDataProcessorPSMBHT2->getHighestAdcFinal() : 0;
	const TH1F *maxAdcHT1c = eventDataProcessorPSHT1 ? eventDataProcessorPSHT1->getHighestAdcFinal() : 0;
	const TH1F *maxAdcHT1HT2c = eventDataProcessorPSHT1HT2 ? eventDataProcessorPSHT1HT2->getHighestAdcFinal() : 0;
	const TH1F *maxEtMBc = eventDataProcessorPSMB ? eventDataProcessorPSMB->highestEtHit.getEt() : 0;
	const TH1F *maxEtMBHT1c = eventDataProcessorPSMBHT1 ? eventDataProcessorPSMBHT1->highestEtHit.getEt() : 0;
	const TH1F *maxEtMBHT2c = eventDataProcessorPSMBHT2 ? eventDataProcessorPSMBHT2->highestEtHit.getEt() : 0;
	const TH1F *maxEtHT1c = eventDataProcessorPSHT1 ? eventDataProcessorPSHT1->highestEtHit.getEt() : 0;
	const TH1F *maxEtHT1HT2c = eventDataProcessorPSHT1HT2 ? eventDataProcessorPSHT1HT2->highestEtHit.getEt() : 0;
	//Int_t PSbinLow = 1;
	if (eventDataProcessorPSMB && eventDataProcessorPSMBHT1) {
		Float_t evNumPSMB = /*maxAdcMBc ? maxAdcMBc->Integral(PSbinLow, maxAdcMBc->GetXaxis()->GetNbins()) : 0;//*/eventDataProcessorPSMB->numPassedAllCuts;
		Float_t evNumPSMBHT1 = /*maxAdcMBHT1c ? maxAdcMBHT1c->Integral(PSbinLow, maxAdcMBHT1c->GetXaxis()->GetNbins()) : 0;//*/eventDataProcessorPSMBHT1->numPassedAllCuts;
		if (evNumPSMBHT1 != 0) PSHT1 = evNumPSMB / (evNumPSMBHT1 * (1.0 + corrMBHT1));
		if (print) {
			cout << "HighTower-1 enhancement calculation from simulation, " << title << ": PSHT1 = " << PSHT1 << endl;
			//cout << "\tMinBias events: " << evNumPSMB << ", simulated HighTower-1 of them: " << evNumPSMBHT1 << ", PSHT1 = " << PSHT1 << endl;
		}
	}
	if (eventDataProcessorPSMB && eventDataProcessorPSMBHT2) {
		Float_t evNumPSMB = /*maxAdcMBc ? maxAdcMBc->Integral(PSbinLow, maxAdcMBc->GetXaxis()->GetNbins()) : 0;//*/eventDataProcessorPSMB->numPassedAllCuts;
		Float_t evNumPSMBHT2 = /*maxAdcMBHT2c ? maxAdcMBHT2c->Integral(PSbinLow, maxAdcMBHT2c->GetXaxis()->GetNbins()) : 0;//*/eventDataProcessorPSMBHT2->numPassedAllCuts;
		if (evNumPSMBHT2 != 0) PSHT2 = evNumPSMB / (evNumPSMBHT2 * (1.0 + corrMBHT2));
		if (print) {
			cout << "HighTower-2 enhancement calculation from simulation, " << title << ": PSHT2 = " << PSHT2 << endl;
			//cout << "\tMinBias events: " << evNumPSMB << ", simulated HighTower-2 of them: " << evNumPSMBHT2 << ", PSHT2 = " << PSHT2 << endl;
		}
	}
	if (eventDataProcessorPSHT1 && eventDataProcessorPSHT1HT2) {
		Float_t PSHT1HT2 = PSHT2 / PSHT1;
		//if (print) {
		//	cout << "PSHT1HT2 = " << PSHT1HT2 << endl;
		//}
		Float_t evNumPSHT1 = /*maxAdcHT1c ? maxAdcHT1c->Integral(PSbinLow, maxAdcHT1c->GetXaxis()->GetNbins()) : 0;//*/eventDataProcessorPSHT1->numPassedAllCuts;
		Float_t evNumPSHT1HT2 = /*maxAdcHT1HT2c ? maxAdcHT1HT2c->Integral(PSbinLow, maxAdcHT1HT2c->GetXaxis()->GetNbins()) : 0;//*/eventDataProcessorPSHT1HT2->numPassedAllCuts;
		Float_t PSHT1HT2_2 = (evNumPSHT1HT2 != 0) ? (evNumPSHT1 / (evNumPSHT1HT2 * (1.0 + corrHT1HT2))) : 0;
		if (print) {
			cout << "HighTower-2 enhancement calculation - 2 from simulation, " << title << ": PSHT1HT2 = " << PSHT1HT2_2 << endl;
			//cout << "\tHighTower-1 events: " << evNumPSHT1 << ", simulated HighTower-2 of them: " << evNumPSHT1HT2 << ", PSHT1HT2 = " << PSHT1HT2_2 << endl;
		}
		if (PSHT1HT2_2 != 0) {
			PSHT2 = (PSHT2 + (PSHT1 * PSHT1HT2_2)) / 2.0;
			if (print) {
				cout << "\tPSHT1HT2 difference: " << (200.0 * fabs(PSHT1HT2 - PSHT1HT2_2) / (PSHT1HT2 + PSHT1HT2_2)) << "%" << endl;
				cout << "\tNew HighTower-2 factor: " << PSHT2 << endl;
			}
		}
	}
	TCanvas *psAdcCanvas = 0;
	TString PSAdctitle = TString(title) + ", Prescales calculation from trigger simulation (max ADC)";
	if (showReal) {
		TString PSname = TString(name) + "_prescales_sim_adc";
		psAdcCanvas = new TCanvas(PSname.Data(), PSAdctitle.Data());
	}
	if (psAdcCanvas) {
		psAdcCanvas->cd();
		TH1F *maxAdcMB = maxAdcMBc ? (new TH1F(*maxAdcMBc)) : 0;
		TH1F *maxAdcMBHT1 = maxAdcMBHT1c ? (new TH1F(*maxAdcMBHT1c)) : 0;
		TH1F *maxAdcMBHT2 = maxAdcMBHT2c ? (new TH1F(*maxAdcMBHT2c)) : 0;
		TH1F *maxAdcHT1 = maxAdcHT1c ? (new TH1F(*maxAdcHT1c)) : 0;
		TH1F *maxAdcHT1HT2 = maxAdcHT1HT2c ? (new TH1F(*maxAdcHT1HT2c)) : 0;
		if (maxAdcMB) {
			maxAdcMB->SetFillColor(10);
			maxAdcMB->SetTitle(PSAdctitle.Data());
			maxAdcMB->Draw("H");
			if (maxAdcHT1) {
				maxAdcHT1->SetFillColor(10);
				maxAdcHT1->Draw("HSAME");
			}
			if (maxAdcHT1HT2) {
				maxAdcHT1HT2->SetFillColor(11);
				maxAdcHT1HT2->Draw("HSAME");
			}
			if (maxAdcMBHT1) {
				maxAdcMBHT1->SetFillColor(11);
				maxAdcMBHT1->Draw("HSAME");
			}
			if (maxAdcMBHT2) {
				maxAdcMBHT2->SetFillColor(12);
				maxAdcMBHT2->Draw("HSAME");
			}
			gPad->SetLogy();
			Float_t x = maxAdcMB->GetXaxis()->GetBinCenter(maxAdcMB->GetXaxis()->GetNbins() / 2);
			Float_t y = maxAdcMB->GetMaximum();
			TString textStr;
			textStr = Form("#splitline{PS_{HT1} = %4.0f}{PS_{HT2} = %4.0f}", PSHT1, PSHT2);
			TLatex *text = new TLatex(x, y, textStr.Data());
			text->Draw();
		}
	}
	TCanvas *psEtCanvas = 0;
	TString PSEttitle = TString(title) + ", Prescales calculation from trigger simulation (max Et)";
	if (showReal) {
		TString PSname = TString(name) + "_prescales_sim_et";
		psEtCanvas = new TCanvas(PSname.Data(), PSEttitle.Data());
	}
	if (psEtCanvas) {
		psEtCanvas->cd();
		TH1F *maxEtMB = maxEtMBc ? (new TH1F(*maxEtMBc)) : 0;
		TH1F *maxEtMBHT1 = maxEtMBHT1c ? (new TH1F(*maxEtMBHT1c)) : 0;
		TH1F *maxEtMBHT2 = maxEtMBHT2c ? (new TH1F(*maxEtMBHT2c)) : 0;
		TH1F *maxEtHT1 = maxEtHT1c ? (new TH1F(*maxEtHT1c)) : 0;
		TH1F *maxEtHT1HT2 = maxEtHT1HT2c ? (new TH1F(*maxEtHT1HT2c)) : 0;
		if (maxEtMB) {
			maxEtMB->SetFillColor(10);
			maxEtMB->SetTitle(PSEttitle.Data());
			maxEtMB->Draw("H");
			if (maxEtHT1) {
				maxEtHT1->SetFillColor(10);
				maxEtHT1->Draw("HSAME");
			}
			if (maxEtHT1HT2) {
				maxEtHT1HT2->SetFillColor(11);
				maxEtHT1HT2->Draw("HSAME");
			}
			if (maxEtMBHT1) {
				maxEtMBHT1->SetFillColor(11);
				maxEtMBHT1->Draw("HSAME");
			}
			if (maxEtMBHT2) {
				maxEtMBHT2->SetFillColor(12);
				maxEtMBHT2->Draw("HSAME");
			}
			gPad->SetLogy();
			Float_t x = maxEtMB->GetXaxis()->GetBinCenter(maxEtMB->GetXaxis()->GetNbins() / 2);
			Float_t y = maxEtMB->GetMaximum();
			TString textStr;
			textStr = Form("#splitline{PS_{HT1} = %4.0f}{PS_{HT2} = %4.0f}", PSHT1, PSHT2);
			TLatex *text = new TLatex(x, y, textStr.Data());
			text->Draw();
			}
	}
}

void calculatePSFromDB(const Int_t runInd, const Int_t *runNums, const Int_t *evNumsMB, const Int_t *evNumsHT1, const Int_t *evNumsHT2, const Int_t *psMB, const Int_t *psHT1, const Int_t *psHT2, Float_t &PS_HT1MB, Float_t &PS_HT2MB, Float_t &PS_HT1MB_alexst, Float_t &PS_HT2MB_alexst, Bool_t print) {
    if (!runNums || !evNumsMB || !evNumsHT1 || !evNumsHT2 || !psMB || !psHT1 || !psHT2) {
        if (print) cout << "No all the data supplied!" << endl;
        return;
    }
    if (print) cout << "Calculating..." << endl;
    Int_t runId = -1;
    Int_t evNumMB = -1;
    Int_t evNumHT1 = -1;
    Int_t evNumHT2 = -1;
    Int_t runPsMB = -1;
    Int_t runPsHT1 = -1;
    Int_t runPsHT2 = -1;
    Double_t evNumTotMB = 0;
    Double_t evNumTotHT1 = 0;
    Double_t evNumTotHT2 = 0;
    Double_t evNumTotMBtrue = 0;
    Double_t evNumTotHT1true = 0;
    Double_t evNumTotHT1Prectrue = 0;
    Double_t evNumTotHT2true = 0;
    Double_t evNumTotHT2Prectrue = 0;
    Double_t evNumTotMBHT1true = 0;
    Double_t evNumTotMBHT2true = 0;
    Double_t evNumTotMBHT1 = 0;
    Double_t evNumTotMBHT2 = 0;
    Double_t evNumTotHT1MBtrue = 0;
    Double_t evNumTotHT1HT2true = 0;
    Double_t evNumTotHT2MBtrue = 0;
    Double_t evNumTotHT2HT1true = 0;
    Double_t evNumTotHT1MBPrectrue = 0;
    Double_t evNumTotHT2MBPrectrue = 0;
    Double_t AHT1Sum = 0;
    Double_t BHT1Sum = 0;
    Double_t BHT1precSum = 0;
    Double_t AHT2Sum = 0;
    Double_t BHT2Sum = 0;
    Double_t BHT2precSum = 0;
    Double_t PSHT1MBSum = 0;
    Double_t PSHT2MBSum = 0;
    for (Int_t i = 0;i < runInd;i++) {
	runId = runNums[i];
	evNumMB = evNumsMB[i];
	evNumHT1 = evNumsHT1[i];
	evNumHT2 = evNumsHT2[i];
	runPsMB = psMB[i];
	runPsHT1 = psHT1[i];
	runPsHT2 = psHT2[i];
	if (runId != -1) {
	    Double_t A = -1;
	    Double_t BHT1 = -1;
	    Double_t BHT2 = -1;
	    Double_t BHT1prec = -1;
	    Double_t BHT2prec = -1;
	    if (evNumMB > 0) {
		if (runPsMB != -1) {
		    if (runPsMB == 0) if (print) cout << "PS for run " << runId << " MB is zero!" << endl;
		    evNumTotMB += evNumMB;
		    A = evNumMB * runPsMB;
		    evNumTotMBtrue += A;
		    if (evNumHT1 > 0) {
			evNumTotMBHT1 += evNumMB;
			evNumTotMBHT1true += A;
			PSHT1MBSum += runPsHT1 * evNumHT1 / Double_t(runPsMB);
		    }
		    if (evNumHT2 > 0) {
			evNumTotMBHT2 += evNumMB;
			evNumTotMBHT2true += A;
			PSHT2MBSum += runPsHT2 * evNumHT2 / Double_t(runPsMB);
		    }
		} else {
		    if (print) cout << "PS for run " << runId << " MB missing!" << endl;
		}
	    }
	    if (evNumHT1 > 0) {
		if (runPsHT1 != -1) {
		    if (runPsHT1 == 0) if (print) cout << "PS for run " << runId << " HT1 is zero!" << endl;
		    evNumTotHT1 += evNumHT1;
		    BHT1 = evNumHT1 * runPsHT1;
		    BHT1prec = ((runPsMB != 0) && (runPsHT1 != 0)) ? (BHT1 / (1 + (runPsHT1 - 1.0) / runPsMB)) : BHT1;
		    evNumTotHT1true += BHT1;
		    evNumTotHT1Prectrue += BHT1;
		    if (evNumMB > 0) {
			evNumTotHT1MBtrue += BHT1;
			evNumTotHT1MBPrectrue += BHT1;
		    }
		    if (evNumHT2 > 0) {
			evNumTotHT1HT2true += BHT1;
		    }
		} else {
		    if (print) cout << "PS for run " << runId << " HT1 missing!" << endl;
		}
	    }
	    if (evNumHT2 > 0) {
		if (runPsHT2 != -1) {
		    if (runPsHT2 == 0) if (print) cout << "PS for run " << runId << " HT2 is zero!" << endl;
		    evNumTotHT2 += evNumHT2;
		    BHT2 = evNumHT2 * runPsHT2;
		    BHT2prec = ((runPsMB != 0) && (runPsHT2 != 0)) ? (BHT2 / (1 + (runPsHT2 - 1.0) / runPsMB)) : BHT2;
		    evNumTotHT2true += BHT2;
		    evNumTotHT2Prectrue += BHT2;
		    if (evNumMB > 0) {
			evNumTotHT2MBtrue += BHT2;
			evNumTotHT2MBPrectrue += BHT2;
		    }
		    if (evNumHT1 > 0) {
			evNumTotHT2HT1true += BHT2;
		    }
		} else {
		    if (print) cout << "PS for run " << runId << " HT2 missing!" << endl;
		}
	    }
	    if (print) cout << (i + 1) << ". A=" << A << ", BHT1=" << BHT1 << ", BHT1prec=" << BHT1prec << ", BHT2=" << BHT2 << ", BHT2prec=" << BHT2prec;
	    if ((A > 0) && (BHT1 > 0)) {
		AHT1Sum += A;
		BHT1Sum += BHT1;
		BHT1precSum += BHT1prec;
	    }
	    if ((A > 0) && (BHT2 > 0)) {
		AHT2Sum += A;
		BHT2Sum += BHT2;
		BHT2precSum += BHT2prec;
	    }
	    if (print) cout << endl;
	}
    }
    if (print) cout << endl;
    if (print) cout << "Total number of events: " << evNumTotMB << " MB, " << evNumTotHT1 << " HT1, " << evNumTotHT2 << " HT2" << endl;
    if (print) cout << "Total true number of events: " << evNumTotMBtrue << " MB, " << evNumTotHT1true << " HT1, " << evNumTotHT2true << " HT2" << endl;
    if (print) cout << "Total precise true number of events: " << evNumTotMBtrue << " MB, " << evNumTotHT1Prectrue << " HT1, " << evNumTotHT2Prectrue << " HT2" << endl;
    Double_t psMBtrue = (evNumTotMB != 0) ? (evNumTotMBtrue / evNumTotMB) : 0;
    if (print) cout << "Average MB prescale = " << psMBtrue << endl;
    Double_t psHT1true = (evNumTotHT1 != 0) ? (evNumTotHT1true / evNumTotHT1) : 0;
    if (print) cout << "Average HT1 prescale = " << psHT1true << endl;
    Double_t psHT2true = (evNumTotHT2 != 0) ? (evNumTotHT2true / evNumTotHT2) : 0;
    if (print) cout << "Average HT2 prescale = " << psHT2true << endl;
    Double_t psHT1Prectrue = (evNumTotHT1 != 0) ? (evNumTotHT1Prectrue / evNumTotHT1) : 0;
    if (print) cout << "Average precise HT1 prescale = " << psHT1Prectrue << endl;
    Double_t psHT2Prectrue = (evNumTotHT2 != 0) ? (evNumTotHT2Prectrue / evNumTotHT2) : 0;
    if (print) cout << "Average precise HT2 prescale = " << psHT2Prectrue << endl;
    if (print) cout << endl;
    Double_t HT1MB = (evNumTotHT1MBtrue != 0) ? (evNumTotMBHT1true / evNumTotHT1MBtrue) : 0;
    if (print) cout << "HT1/MB = " << HT1MB << endl;
    Double_t HT2HT1 = (evNumTotHT2HT1true != 0) ? (evNumTotHT1HT2true / evNumTotHT2HT1true) : 0;
    if (print) cout << "HT2/HT1 = " << HT2HT1 << endl;
    Double_t HT2MB = (evNumTotHT2MBtrue != 0) ? (evNumTotMBHT2true / evNumTotHT2MBtrue) : 0;
    if (print) cout << "HT2/MB = " << HT2MB << ", HT1/MB * HT2/HT1 = " << (HT1MB * HT2HT1) << endl;
    HT2MB = (HT2MB + (HT1MB * HT2HT1)) / 2.0;
    if (print) cout << "HT2/MB (updated) = " << HT2MB << endl;
    if (print) cout << endl;
    PS_HT1MB = HT1MB;
    PS_HT2MB = HT2MB;
    Double_t HT1MBprec = (BHT1precSum != 0) ? (AHT1Sum / BHT1precSum) : 0;
    if (print) cout << "Precise HT1/MB = " << HT1MBprec << endl;
    Double_t HT2MBprec = (BHT2precSum != 0) ? (AHT2Sum / BHT2precSum) : 0;
    if (print) cout << "Precise HT2/MB = " << HT2MBprec << endl;
    if (print) cout << endl;
    Double_t FHT1 = (PSHT1MBSum != 0) ? (Double_t(evNumTotHT1) / PSHT1MBSum) : 0;
    Double_t FHT2 = (PSHT2MBSum != 0) ? (Double_t(evNumTotHT2) / PSHT2MBSum) : 0;
    PS_HT1MB_alexst = (evNumTotHT1 != 0) ? (evNumTotMBHT1 * FHT1 / evNumTotHT1) : 0;
    PS_HT2MB_alexst = (evNumTotHT2 != 0) ? (evNumTotMBHT2 * FHT2 / evNumTotHT2) : 0;
    if (print) cout << "Number of MB events in runs with HT1 = " << evNumTotMBHT1 << endl;
    if (print) cout << "Number of MB events in runs with HT2 = " << evNumTotMBHT2 << endl;
    if (print) cout << "F_HT1 (Alex St.) = " << FHT1 << " (" << PS_HT1MB_alexst << ")" << endl;
    if (print) cout << "F_HT2 (Alex St.) = " << FHT2 << " (" << PS_HT2MB_alexst << ")" << endl;
    if (print) cout << "Check HT1: " << evNumTotMBHT1 << " * " << FHT1 << " = " << (evNumTotMBHT1 * FHT1) << ", " << evNumTotHT1 << " * " << HT1MB << " = " << (evNumTotHT1 * HT1MB) << " (equal?)" << endl;
    if (print) cout << "Check HT2: " << evNumTotMBHT2 << " * " << FHT2 << " = " << (evNumTotMBHT2 * FHT2) << ", " << evNumTotHT2 << " * " << HT2MB << " = " << (evNumTotHT2 * HT2MB) << " (equal?)" << endl;
    if (print) cout << "Done." << endl;
}

Int_t getRunYear(const TEventDataProcessor *events) {
	Int_t year = 0;
	const TH1F *histStatYear = events ? events->getStatYear() : 0;
	if (histStatYear) {
	    for (Int_t ibin = 1;(ibin <= histStatYear->GetXaxis()->GetNbins()) && year == 0;ibin++) {
		if (histStatYear->GetBinContent(ibin) > 0) year = (Int_t)histStatYear->GetXaxis()->GetBinCenter(ibin);
	    }
	}
	if (year <= 0) year = 4;
	return year;
}

Int_t *prescalesRuns = 0;
Int_t *prescalesMB = 0;
Int_t *prescalesHT1 = 0;
Int_t *prescalesHT2 = 0;
Int_t prescalesNum = 0;
void clearCachedPrescales() {
    if (prescalesRuns) delete [] prescalesRuns; prescalesRuns = 0;
    if (prescalesMB) delete [] prescalesMB; prescalesMB = 0;
    if (prescalesHT1) delete [] prescalesHT1; prescalesHT1 = 0;
    if (prescalesHT2) delete [] prescalesHT2; prescalesHT2 = 0;
    prescalesNum = 0;
}

void calculatePSFromDB(const Char_t *name, const Char_t *title, const Char_t *psFilename
    , const TEventDataProcessor *eventDataProcessorMB, const TEventDataProcessor *eventDataProcessorHT1, const TEventDataProcessor *eventDataProcessorHT2
    , Bool_t show, Bool_t print, Float_t &PSHT1, Float_t &PSHT2, Float_t &PSHT1_alexst, Float_t &PSHT2_alexst) {
	if (print) cout << "Prescales calculation from the DB data, " << title << endl;
	const TH2F *histMB = eventDataProcessorMB ? eventDataProcessorMB->getDayRun() : 0;
	const TH2F *histHT1 = eventDataProcessorHT1 ? eventDataProcessorHT1->getDayRun() : 0;
	const TH2F *histHT2 = eventDataProcessorHT2 ? eventDataProcessorHT2->getDayRun() : 0;
	if (!histMB || !histHT1 || !histHT2) {
		if (print) cout << "No data found for " << title << " (histMB = " << histMB << ", histHT1 = " << histHT1 << ", histHT2 = " << histHT2 << ")" << endl;
		return;
	}
	Int_t year = getRunYear(eventDataProcessorMB);
	
	Int_t buffsize = 0;
	for (Int_t iday = 1;iday <= histMB->GetXaxis()->GetNbins();iday++) {
		for (Int_t irun = 1;irun <= histMB->GetYaxis()->GetNbins();irun++) {
			Bool_t runUsed = false;
			Int_t ibin = histMB->GetBin(iday, irun);
			//cout << "iday = " << iday << ", irun = " << ", ibin" << ibin << endl;
			Int_t evNumMB = (Int_t)histMB->GetBinContent(ibin);
			Int_t evNumHT1 = (Int_t)histHT1->GetBinContent(ibin);
			Int_t evNumHT2 = (Int_t)histHT2->GetBinContent(ibin);
			if (evNumMB != 0) {
				runUsed = true;
			}
			if (evNumHT1 != 0) {
				runUsed = true;
			}
			if (evNumHT2 != 0) {
				runUsed = true;
			}
			if (runUsed) {
				buffsize++;
			}
		}
	}

	Int_t runInd = 0;
	Int_t *runNums = new Int_t[buffsize];
	Int_t *evNumsMB = new Int_t[buffsize];
	Int_t *evNumsHT1 = new Int_t[buffsize];
	Int_t *evNumsHT2 = new Int_t[buffsize];
	if (runNums) for (Int_t i = 0;i < buffsize;runNums[i++] = 0);
	if (evNumsMB) for (Int_t i = 0;i < buffsize;evNumsMB[i++] = 0);
	if (evNumsHT1) for (Int_t i = 0;i < buffsize;evNumsHT1[i++] = 0);
	if (evNumsHT2) for (Int_t i = 0;i < buffsize;evNumsHT2[i++] = 0);
	for (Int_t iday = 1;(iday <= histMB->GetXaxis()->GetNbins()) && (runInd < buffsize);iday++) {
		for (Int_t irun = 1;(irun <= histMB->GetYaxis()->GetNbins()) && (runInd < buffsize);irun++) {
			Bool_t runUsed = false;
			Int_t ibin = histMB->GetBin(iday, irun);
			//cout << "iday = " << iday << ", irun = " << ", ibin" << ibin << endl;
			Int_t evNumMB = (Int_t)histMB->GetBinContent(ibin);
			Int_t evNumHT1 = (Int_t)histHT1->GetBinContent(ibin);
			Int_t evNumHT2 = (Int_t)histHT2->GetBinContent(ibin);
			if (evNumMB != 0) {
				runUsed = true;
				if (evNumsMB) evNumsMB[runInd] += evNumMB;
			}
			if (evNumHT1 != 0) {
				runUsed = true;
				if (evNumsHT1) evNumsHT1[runInd] += evNumHT1;
			}
			if (evNumHT2 != 0) {
				runUsed = true;
				if (evNumsHT2) evNumsHT2[runInd] += evNumHT2;
			}
			if (runUsed) {
				Int_t runId = (1000000 * year) + (1000 * (Int_t)histMB->GetXaxis()->GetBinLowEdge(iday)) + (1 * (Int_t)histMB->GetXaxis()->GetBinLowEdge(irun));
				if (runNums) runNums[runInd] = runId;
				//if (print) cout << "Added run #" << runInd << ": " << runId;
				//if (print && evNumsMB) cout << ", MB " << evNumsMB[runInd];
				//if (print && evNumsHT1) cout << ", HT1 " << evNumsHT1[runInd];
				//if (print && evNumsHT2) cout << ", HT2 " << evNumsHT2[runInd];
				//if (print) cout << endl;
				runInd++;
			}
		}
	}
	if (runInd >= buffsize) runInd = buffsize - 1;
	if (!prescalesRuns || !prescalesMB || !prescalesHT1 || !prescalesHT2) {
		clearCachedPrescales();
		const Int_t buffsizePS = 32768;
		prescalesRuns = new Int_t[buffsizePS];
		prescalesMB = new Int_t[buffsizePS];
		prescalesHT1 = new Int_t[buffsizePS];
		prescalesHT2 = new Int_t[buffsizePS];
		if (prescalesRuns) for (Int_t i = 0;i < buffsizePS;prescalesRuns[i++] = 0);
		if (prescalesMB) for (Int_t i = 0;i < buffsizePS;prescalesMB[i++] = -1);
		if (prescalesHT1) for (Int_t i = 0;i < buffsizePS;prescalesHT1[i++] = -1);
		if (prescalesHT2) for (Int_t i = 0;i < buffsizePS;prescalesHT2[i++] = -1);
		TString psFilenameExact = findFile(psFilename);
		if (print) cout << "Reading prescales from " << psFilename << " (" << psFilenameExact << ")" << endl;
		ifstream ipsStr(psFilenameExact);
		if (!ipsStr.good()) {
			cerr << "Cannot open prescales file " << psFilename << " !" << endl;
		} else {
			Int_t runId, runPsMB, runPsHT1, runPsHT2;
			while (ipsStr.good()) {
				ipsStr >> runId;
				ipsStr >> runPsMB;
				ipsStr >> runPsHT1;
				ipsStr >> runPsHT2;
				if (ipsStr.good() && (prescalesNum < buffsizePS)) {
					//if (print) cout << (prescalesNum + 1) << ". Read run " << runId << ": " << runPsMB << " " << runPsHT1 << " " << runPsHT2 << endl;
					prescalesRuns[prescalesNum] = runId;
					prescalesMB[prescalesNum] = runPsMB;
					prescalesHT1[prescalesNum] = runPsHT1;
					prescalesHT2[prescalesNum] = runPsHT2;
					prescalesNum++;
				}
			}
		}
		ipsStr.close();
		if (prescalesNum >= buffsizePS) prescalesNum = buffsizePS - 1;
	}

	Int_t *psMB = new Int_t[buffsize];
	Int_t *psHT1 = new Int_t[buffsize];
	Int_t *psHT2 = new Int_t[buffsize];
	if (psMB) for (Int_t i = 0;i < buffsize;psMB[i++] = -1);
	if (psHT1) for (Int_t i = 0;i < buffsize;psHT1[i++] = -1);
	if (psHT2) for (Int_t i = 0;i < buffsize;psHT2[i++] = -1);
	for (Int_t i = 0;i < runInd;i++) {
		Bool_t found = false;
		for (Int_t ips = 0;(ips < prescalesNum) && !found;ips++) {
			if (runNums && prescalesRuns && (runNums[i] == prescalesRuns[ips])) {
				if (psMB && prescalesMB) psMB[i] = prescalesMB[ips];
				if (psHT1 && prescalesHT1) psHT1[i] = prescalesHT1[ips];
				if (psHT2 && prescalesHT2) psHT2[i] = prescalesHT2[ips];
				found = true;
			}
		}
		if (found) {
			//if (print) cout << i << ". Run " << runNums[i];
			//if (print && psMB) cout << ", psMB " << psMB[i];
			//if (print && psHT1) cout << ", psHT1 " << psHT1[i];
			//if (print && psHT2) cout << ", psHT2 " << psHT2[i];
			//if (print) cout << endl;
		} else {
			if (print) cout << "Prescales not found for run #" << i << endl;
		}
	}
	if (psMB) for (Int_t i = 0;i < runInd;i++) if (psMB[i] == -1) if (print) cout << "No prescale MB found for run #" << i << " " << runNums[i] << endl;
	if (psHT1) for (Int_t i = 0;i < runInd;i++) if (psHT1[i] == -1) if (print) cout << "No prescale HT1 found for run #" << i << " " << runNums[i] << endl;
	if (psHT2) for (Int_t i = 0;i < runInd;i++) if (psHT2[i] == -1) if (print) cout << "No prescale HT2 found for run #" << i << " " << runNums[i] << endl;

	if (runNums && evNumsMB && evNumsHT1 && evNumsHT2 && psMB && psHT1 && psHT2) {
		//if (print) cout << "Calling the main routine..." << endl;
		calculatePSFromDB(runInd, runNums, evNumsMB, evNumsHT1, evNumsHT2, psMB, psHT1, psHT2, PSHT1, PSHT2, PSHT1_alexst, PSHT2_alexst, false);
		if (print) cout << "PSHT1 = " << PSHT1 << ", PSHT2 = " << PSHT2 << ", PSHT1_alexst = " << PSHT1_alexst << ", PSHT2_alexst = " << PSHT2_alexst << endl;
	}

	if (runNums) delete [] runNums;
	if (evNumsMB) delete [] evNumsMB;
	if (evNumsHT1) delete [] evNumsHT1;
	if (evNumsHT2) delete [] evNumsHT2;
	if (psMB) delete [] psMB;
	if (psHT1) delete [] psHT1;
	if (psHT2) delete [] psHT2;
	//if (print) cout << "...done " << title << endl;
}

void showResultsOthers(TH1 *histSpectrum, TLegend *legendSpectrum
	, Bool_t showpQCDWeight, Bool_t usePPpQCD
) {
	if (!histSpectrum) return;
	const TWeightCalculator &pQCDWeightCalc = PQCDWEIGHT(usePPpQCD);
	histSpectrum->SetMinimum(1e-10);
	histSpectrum->SetMaximum(10);
	Float_t Xmin = histSpectrum->GetXaxis()->GetBinLowEdge(1);
	Float_t Xmax = histSpectrum->GetXaxis()->GetBinUpEdge(histSpectrum->GetXaxis()->GetNbins());
	gPad->SetLogy();
	if (showpQCD) {
		TF1 *pQCDfunc = new TF1("pQCDfunc", PQCDFUNC(usePPpQCD), TMath::Max(Xmin, (Float_t)1.5), TMath::Min(Xmax, (Float_t)20), 1);
		pQCDfunc->SetLineWidth(2);
		pQCDfunc->SetLineStyle(1);
		pQCDfunc->SetParameter(0, 0);
		TF1 *pQCDfuncCopy = pQCDfunc->DrawCopy("SAME");
		TF1 *pQCD1func = new TF1("pQCD1func", PQCDFUNC(usePPpQCD), TMath::Max(Xmin, (Float_t)1.5), TMath::Min(Xmax, (Float_t)20), 1);
		pQCD1func->SetLineWidth(2);
		pQCD1func->SetLineStyle(2);
		pQCD1func->SetParameter(0, 1);
		TF1 *pQCD1funcCopy = pQCD1func->DrawCopy("SAME");
		TF1 *pQCD2func = new TF1("pQCD2func", PQCDFUNC(usePPpQCD), TMath::Max(Xmin, (Float_t)1.5), TMath::Min(Xmax, (Float_t)20), 1);
		pQCD2func->SetLineWidth(2);
		pQCD2func->SetLineStyle(2);
		pQCD2func->SetParameter(0, 2);
		TF1 *pQCD2funcCopy = pQCD2func->DrawCopy("SAME");
		if (legendSpectrum) {
			if (pQCDfuncCopy) legendSpectrum->AddEntry(pQCDfuncCopy, showpQCDLegend, "L");
			if (pQCD1funcCopy || pQCD2funcCopy) legendSpectrum->AddEntry(pQCD1funcCopy ? pQCD1funcCopy : pQCD2funcCopy, showpQCDErrorsLegend, "L");
		}
	}
	if (showStarChargedHadrons) {
		TGraphErrors *hgrS = new TGraphErrors(sizeof(Sx)/sizeof(Sx[0]), Sx, Sy, Sxe, Sye);
		hgrS->SetLineColor(kGreen);
		hgrS->SetMarkerColor(hgrS->GetLineColor());
		hgrS->SetMarkerStyle(22);
		hgrS->SetMarkerSize(1.3);
		hgrS->SetFillColor(hgrS->GetLineColor());
		hgrS->SetHistogram(histSpectrum);
		hgrS->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrS, showStarChargedHadronsLegend, "P");			
	}
	if (showStarChargedPions) {
		TGraphErrors *hgrSPp = new TGraphErrors(sizeof(SPpx)/sizeof(SPpx[0]), SPpx, SPpy, SPpxe, SPpye);
		hgrSPp->SetLineColor(kBlue);
		hgrSPp->SetMarkerColor(hgrSPp->GetLineColor());
		hgrSPp->SetMarkerStyle(29);
		hgrSPp->SetMarkerSize(1.3);
		hgrSPp->SetFillColor(hgrSPp->GetLineColor());
		hgrSPp->SetHistogram(histSpectrum);
		hgrSPp->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrSPp, showStarChargedPionsPlusLegend, "P");			
		TGraphErrors *hgrSPm = new TGraphErrors(sizeof(SPmx)/sizeof(SPmx[0]), SPmx, SPmy, SPmxe, SPmye);
		hgrSPm->SetLineColor(kBlue);
		hgrSPm->SetMarkerColor(hgrSPm->GetLineColor());
		hgrSPm->SetMarkerStyle(30);
		hgrSPm->SetMarkerSize(3);
		hgrSPm->SetFillColor(hgrSPm->GetLineColor());
		hgrSPm->SetHistogram(histSpectrum);
		hgrSPm->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrSPm, showStarChargedPionsMinusLegend, "P");			
	}
	if (showStolpovskyPi0) {
		TGraphErrors *hgrASMB = new TGraphErrors(sizeof(ASxMB)/sizeof(ASxMB[0]), ASxMB, ASyMB, ASxeMB, ASyeMB);
		hgrASMB->SetLineColor(kBlack);
		hgrASMB->SetMarkerColor(hgrASMB->GetLineColor());
		hgrASMB->SetFillColor(hgrASMB->GetLineColor());
		hgrASMB->SetMarkerStyle(23);
		hgrASMB->SetMarkerSize(1.3);
		hgrASMB->SetHistogram(histSpectrum);
		hgrASMB->Draw("P");
		TGraphErrors *hgrASHT1 = new TGraphErrors(sizeof(ASxHT1)/sizeof(ASxHT1[0]), ASxHT1, ASyHT1, ASxeHT1, ASyeHT1);
		hgrASHT1->SetLineColor(kBlue);
		hgrASHT1->SetMarkerColor(hgrASHT1->GetLineColor());
		hgrASHT1->SetFillColor(hgrASHT1->GetLineColor());
		hgrASHT1->SetMarkerStyle(hgrASMB->GetMarkerStyle());
		hgrASHT1->SetMarkerSize(hgrASMB->GetMarkerSize());
		hgrASHT1->SetHistogram(histSpectrum);
		hgrASHT1->Draw("P");
		TGraphErrors *hgrASHT2 = new TGraphErrors(sizeof(ASxHT2)/sizeof(ASxHT2[0]), ASxHT2, ASyHT2, ASxeHT2, ASyeHT2);
		hgrASHT2->SetLineColor(kRed);
		hgrASHT2->SetMarkerColor(hgrASHT2->GetLineColor());
		hgrASHT2->SetFillColor(hgrASHT2->GetLineColor());
		hgrASHT2->SetMarkerStyle(hgrASMB->GetMarkerStyle());
		hgrASHT2->SetMarkerSize(hgrASMB->GetMarkerSize());
		hgrASHT2->SetHistogram(histSpectrum);
		hgrASHT2->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrASMB, showStolpovskyPi0Legend, "P");
	}
	if (showPhenixPi0) {
		TGraphErrors *hgrP = new TGraphErrors(sizeof(Px)/sizeof(Px[0]), Px, Py, Pxe, Pye);
		hgrP->SetLineColor(1);
		hgrP->SetMarkerColor(hgrP->GetLineColor());
		hgrP->SetMarkerStyle(25);
		hgrP->SetMarkerSize(1.1);
		hgrP->SetFillColor(hgrP->GetLineColor());
		hgrP->SetHistogram(histSpectrum);
		hgrP->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrP, showPhenixPi0Legend, "P");
	}
	if (showPhenixPi0PP) {
		TGraphErrors *hgrPpp = new TGraphErrors(sizeof(Pppx)/sizeof(Pppx[0]), Pppx, Pppy, Pppxe, Pppye);
		hgrPpp->SetLineColor(1);
		hgrPpp->SetMarkerColor(hgrPpp->GetLineColor());
		hgrPpp->SetMarkerStyle(25);
		hgrPpp->SetMarkerSize(1.1);
		hgrPpp->SetFillColor(hgrPpp->GetLineColor());
		hgrPpp->SetHistogram(histSpectrum);
		hgrPpp->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrPpp, showPhenixPi0PPLegend, "P");
	}
	if (showStolpovskyPi0PP) {
		TGraphErrors *hgrASppMB = new TGraphErrors(sizeof(ASppxMB)/sizeof(ASppxMB[0]), ASppxMB, ASppyMB, ASppxeMB, ASppyeMB);
		hgrASppMB->SetLineColor(kBlack);
		hgrASppMB->SetMarkerColor(hgrASppMB->GetLineColor());
		hgrASppMB->SetFillColor(hgrASppMB->GetLineColor());
		hgrASppMB->SetMarkerStyle(21);
		hgrASppMB->SetMarkerSize(1.3);
		hgrASppMB->SetHistogram(histSpectrum);
		hgrASppMB->Draw("P");
		TGraphErrors *hgrASppHT1 = new TGraphErrors(sizeof(ASppxHT1)/sizeof(ASppxHT1[0]), ASppxHT1, ASppyHT1, ASppxeHT1, ASppyeHT1);
		hgrASppHT1->SetLineColor(kBlue);
		hgrASppHT1->SetMarkerColor(hgrASppHT1->GetLineColor());
		hgrASppHT1->SetFillColor(hgrASppHT1->GetLineColor());
		hgrASppHT1->SetMarkerStyle(hgrASppMB->GetMarkerStyle());
		hgrASppHT1->SetMarkerSize(hgrASppMB->GetMarkerSize());
		hgrASppHT1->SetHistogram(histSpectrum);
		hgrASppHT1->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrASppMB, showStolpovskyPi0PPLegend, "P");
	}
	if (showMischkePi0) {
		TGraphErrors *hgrA = new TGraphErrors(sizeof(Ax)/sizeof(Ax[0]), Ax, Ay, Axe, Aye);
		hgrA->SetLineColor(kGreen);
		hgrA->SetMarkerColor(hgrA->GetLineColor());
		hgrA->SetMarkerStyle(22);
		hgrA->SetMarkerSize(1.3);
		hgrA->SetFillColor(hgrA->GetLineColor());
		hgrA->SetHistogram(histSpectrum);
		hgrA->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrA, showMischkePi0Legend, "P");
	}
	if (showMartijnPi0dAu) {
		TGraphErrors *hgrMartijnMB = new TGraphErrors(sizeof(MartijnxMB)/sizeof(MartijnxMB[0]), MartijnxMB, MartijnyMB, MartijnxeMB, MartijnyeMB);
		hgrMartijnMB->SetLineColor(6);
		hgrMartijnMB->SetMarkerColor(hgrMartijnMB->GetLineColor());
		hgrMartijnMB->SetFillColor(hgrMartijnMB->GetLineColor());
		hgrMartijnMB->SetMarkerStyle(23);
		hgrMartijnMB->SetMarkerSize(1.3);
		hgrMartijnMB->SetHistogram(histSpectrum);
		hgrMartijnMB->Draw("P");
		TGraphErrors *hgrMartijnHT1 = new TGraphErrors(sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0]), MartijnxHT1, MartijnyHT1, MartijnxeHT1, MartijnyeHT1);
		hgrMartijnHT1->SetLineColor(6);
		hgrMartijnHT1->SetMarkerColor(hgrMartijnHT1->GetLineColor());
		hgrMartijnHT1->SetFillColor(hgrMartijnHT1->GetLineColor());
		hgrMartijnHT1->SetMarkerStyle(hgrMartijnMB->GetMarkerStyle());
		hgrMartijnHT1->SetMarkerSize(hgrMartijnMB->GetMarkerSize());
		hgrMartijnHT1->SetHistogram(histSpectrum);
		hgrMartijnHT1->Draw("P");
		TGraphErrors *hgrMartijnHT2 = new TGraphErrors(sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0]), MartijnxHT2, MartijnyHT2, MartijnxeHT2, MartijnyeHT2);
		hgrMartijnHT2->SetLineColor(6);
		hgrMartijnHT2->SetMarkerColor(hgrMartijnHT2->GetLineColor());
		hgrMartijnHT2->SetFillColor(hgrMartijnHT2->GetLineColor());
		hgrMartijnHT2->SetMarkerStyle(hgrMartijnMB->GetMarkerStyle());
		hgrMartijnHT2->SetMarkerSize(hgrMartijnMB->GetMarkerSize());
		hgrMartijnHT2->SetHistogram(histSpectrum);
		hgrMartijnHT2->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrMartijnMB, showMartijnPi0dAuLegend, "P");
	}
	if (showSimonPi0PPCrossection) {
		TGraphErrors *hgrFSxs = new TGraphErrors(sizeof(FSxs_x)/sizeof(FSxs_x[0]), FSxs_x, FSxs_y, FSxs_xe, FSxs_ye);
		hgrFSxs->SetLineColor(kGreen);
		hgrFSxs->SetMarkerColor(hgrFSxs->GetLineColor());
		hgrFSxs->SetMarkerStyle(22);
		hgrFSxs->SetMarkerSize(1.3);
		hgrFSxs->SetFillColor(hgrFSxs->GetLineColor());
		hgrFSxs->SetHistogram(histSpectrum);
		hgrFSxs->Draw("P");
		if (legendSpectrum) legendSpectrum->AddEntry(hgrFSxs, showSimonPi0PPCrossectionLegend, "P");			
	}
	if (showpQCDWeight) {
		TF1 *funcpQCD = pQCDWeightCalc.createFunc(true);
		TF1 *funcpQCDCopy = 0;
		if (funcpQCD) {
			funcpQCD->SetLineWidth(1);
			funcpQCD->SetLineColor(kBlue);
			funcpQCD->SetLineStyle(1);
			funcpQCD->SetRange(0, funcpQCD->GetXmax());
			funcpQCDCopy = funcpQCD->DrawCopy("SAME");
		}
		if (legendSpectrum) legendSpectrum->AddEntry(funcpQCDCopy, pQCDWeightCalc.GetTitle(), "L");
		cout << "+--------------------------------------" << endl;
		pQCDWeightCalc.Print("| ");
		cout << "+--------------------------------------" << endl;
	}
	if (legendSpectrum) legendSpectrum->Draw();
}

void showResultsOthersEtaToPi0(TH1 *hist, TLegend *legend) {
	if (!hist) return;
	hist->SetMinimum(0);
	hist->SetMaximum(2);
/*
	Float_t Xmin = hist->GetXaxis()->GetBinLowEdge(1);
	Float_t Xmax = hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins());
*/
}

void showResultsOthersDiv(TH1 *histDiv, TLegend *legend, TF1 *weightDiv, Bool_t usePPpQCD) {
	if (!histDiv) return;
	histDiv->SetMinimum(0);
	histDiv->SetMaximum(5);
	Float_t Xmin = histDiv->GetXaxis()->GetBinLowEdge(1);
	Float_t Xmax = histDiv->GetXaxis()->GetBinUpEdge(histDiv->GetXaxis()->GetNbins());
	Int_t curInd;
	if (showpQCD) {
		TF1 *pQCDfunc = new TF1("pQCDfuncDiv", PQCDFUNC(usePPpQCD), TMath::Max(Xmin, (Float_t)1.5), TMath::Min(Xmax, (Float_t)20), 1);
		pQCDfunc->SetLineWidth(2);
		pQCDfunc->SetLineStyle(1);
		pQCDfunc->SetParameter(0, 10);
		TF1 *pQCDfuncCopy = pQCDfunc->DrawCopy("SAME");
		TF1 *pQCD1func = new TF1("pQCD1funcDiv", PQCDFUNC(usePPpQCD), TMath::Max(Xmin, (Float_t)1.5), TMath::Min(Xmax, (Float_t)20), 1);
		pQCD1func->SetLineWidth(2);
		pQCD1func->SetLineStyle(2);
		pQCD1func->SetParameter(0, 11);
		TF1 *pQCD1funcCopy = pQCD1func->DrawCopy("SAME");
		TF1 *pQCD2func = new TF1("pQCD2funcDiv", PQCDFUNC(usePPpQCD), TMath::Max(Xmin, (Float_t)1.5), TMath::Min(Xmax, (Float_t)20), 1);
		pQCD2func->SetLineWidth(2);
		pQCD2func->SetLineStyle(2);
		pQCD2func->SetParameter(0, 12);
		TF1 *pQCD2funcCopy = pQCD2func->DrawCopy("SAME");
		if (legend) {
			if (pQCDfuncCopy) legend->AddEntry(pQCDfuncCopy, showpQCDLegend, "L");
			if (pQCD1funcCopy || pQCD2funcCopy) legend->AddEntry(pQCD1funcCopy ? pQCD1funcCopy : pQCD2funcCopy, showpQCDErrorsLegend, "L");
		}
	}
	if (showStarChargedHadrons) {
		Float_t *spectrumSX = new Float_t[sizeof(Sx)/sizeof(Sx[0])];
		Float_t *spectrumSXerr = new Float_t[sizeof(Sx)/sizeof(Sx[0])];
		Float_t *spectrumSY = new Float_t[sizeof(Sx)/sizeof(Sx[0])];
		Float_t *spectrumSYerr = new Float_t[sizeof(Sx)/sizeof(Sx[0])];
		for (curInd = 0;curInd < Int_t(sizeof(Sx)/sizeof(Sx[0])); curInd++) {
			spectrumSX[curInd] = Sx[curInd];
			spectrumSXerr[curInd] = Sxe[curInd];
			spectrumSY[curInd] = Sy[curInd] / weightDiv->Eval(spectrumSX[curInd]);
			spectrumSYerr[curInd] = Sye[curInd] / weightDiv->Eval(spectrumSX[curInd]);
		}
		TGraphErrors *hgrSpectrumS = new TGraphErrors(sizeof(Sx)/sizeof(Sx[0]), spectrumSX, spectrumSY, spectrumSXerr, spectrumSYerr);
		hgrSpectrumS->SetLineColor(kGreen);
		hgrSpectrumS->SetFillColor(hgrSpectrumS->GetLineColor());
		hgrSpectrumS->SetMarkerColor(hgrSpectrumS->GetLineColor());
		hgrSpectrumS->SetMarkerStyle(22);
		hgrSpectrumS->SetMarkerSize(1.3);
		hgrSpectrumS->SetHistogram(histDiv);
		hgrSpectrumS->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumS, showStarChargedHadronsLegend, "P");
	}
	if (showStarChargedPions) {
		Float_t *spectrumSPpX = new Float_t[sizeof(SPpx)/sizeof(SPpx[0])];
		Float_t *spectrumSPpXerr = new Float_t[sizeof(SPpx)/sizeof(SPpx[0])];
		Float_t *spectrumSPpY = new Float_t[sizeof(SPpx)/sizeof(SPpx[0])];
		Float_t *spectrumSPpYerr = new Float_t[sizeof(SPpx)/sizeof(SPpx[0])];
		for (curInd = 0;curInd < Int_t(sizeof(SPpx)/sizeof(SPpx[0])); curInd++) {
			spectrumSPpX[curInd] = SPpx[curInd];
			spectrumSPpXerr[curInd] = SPpxe[curInd];
			spectrumSPpY[curInd] = SPpy[curInd] / weightDiv->Eval(spectrumSPpX[curInd]);
			spectrumSPpYerr[curInd] = SPpye[curInd] / weightDiv->Eval(spectrumSPpX[curInd]);
		}
		TGraphErrors *hgrSpectrumSPp = new TGraphErrors(sizeof(SPpx)/sizeof(SPpx[0]), spectrumSPpX, spectrumSPpY, spectrumSPpXerr, spectrumSPpYerr);
		hgrSpectrumSPp->SetLineColor(kBlue);
		hgrSpectrumSPp->SetFillColor(hgrSpectrumSPp->GetLineColor());
		hgrSpectrumSPp->SetMarkerColor(hgrSpectrumSPp->GetLineColor());
		hgrSpectrumSPp->SetMarkerStyle(29);
		hgrSpectrumSPp->SetMarkerSize(1.3);
		hgrSpectrumSPp->SetHistogram(histDiv);
		hgrSpectrumSPp->Draw("P");
		Float_t *spectrumSPmX = new Float_t[sizeof(SPmx)/sizeof(SPmx[0])];
		Float_t *spectrumSPmXerr = new Float_t[sizeof(SPmx)/sizeof(SPmx[0])];
		Float_t *spectrumSPmY = new Float_t[sizeof(SPmx)/sizeof(SPmx[0])];
		Float_t *spectrumSPmYerr = new Float_t[sizeof(SPmx)/sizeof(SPmx[0])];
		for (curInd = 0;curInd < Int_t(sizeof(SPmx)/sizeof(SPmx[0])); curInd++) {
			spectrumSPmX[curInd] = SPmx[curInd];
			spectrumSPmXerr[curInd] = SPmxe[curInd];
			spectrumSPmY[curInd] = SPmy[curInd] / weightDiv->Eval(spectrumSPmX[curInd]);
			spectrumSPmYerr[curInd] = SPmye[curInd] / weightDiv->Eval(spectrumSPmX[curInd]);
		}
		if (legend) legend->AddEntry(hgrSpectrumSPp, showStarChargedPionsPlusLegend, "P");
		TGraphErrors *hgrSpectrumSPm = new TGraphErrors(sizeof(SPmx)/sizeof(SPmx[0]), spectrumSPmX, spectrumSPmY, spectrumSPmXerr, spectrumSPmYerr);
		hgrSpectrumSPm->SetLineColor(kBlue);
		hgrSpectrumSPm->SetFillColor(hgrSpectrumSPm->GetLineColor());
		hgrSpectrumSPm->SetMarkerColor(hgrSpectrumSPm->GetLineColor());
		hgrSpectrumSPm->SetMarkerStyle(30);
		hgrSpectrumSPm->SetMarkerSize(1.3);
		hgrSpectrumSPm->SetHistogram(histDiv);
		hgrSpectrumSPm->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumSPm, showStarChargedPionsMinusLegend, "P");
	}
	if (showStolpovskyPi0) {
		Float_t *spectrumASXMB = new Float_t[sizeof(ASxMB)/sizeof(ASxMB[0])];
		Float_t *spectrumASXerrMB = new Float_t[sizeof(ASxMB)/sizeof(ASxMB[0])];
		Float_t *spectrumASYMB = new Float_t[sizeof(ASxMB)/sizeof(ASxMB[0])];
		Float_t *spectrumASYerrMB = new Float_t[sizeof(ASxMB)/sizeof(ASxMB[0])];
		Float_t *spectrumASXHT1 = new Float_t[sizeof(ASxHT1)/sizeof(ASxHT1[0])];
		Float_t *spectrumASXerrHT1 = new Float_t[sizeof(ASxHT1)/sizeof(ASxHT1[0])];
		Float_t *spectrumASYHT1 = new Float_t[sizeof(ASxHT1)/sizeof(ASxHT1[0])];
		Float_t *spectrumASYerrHT1 = new Float_t[sizeof(ASxHT1)/sizeof(ASxHT1[0])];
		Float_t *spectrumASXHT2 = new Float_t[sizeof(ASxHT2)/sizeof(ASxHT2[0])];
		Float_t *spectrumASXerrHT2 = new Float_t[sizeof(ASxHT2)/sizeof(ASxHT2[0])];
		Float_t *spectrumASYHT2 = new Float_t[sizeof(ASxHT2)/sizeof(ASxHT2[0])];
		Float_t *spectrumASYerrHT2 = new Float_t[sizeof(ASxHT2)/sizeof(ASxHT2[0])];
		for (curInd = 0;curInd < Int_t(sizeof(ASxMB)/sizeof(ASxMB[0])); curInd++) {
			spectrumASXMB[curInd] = ASxMB[curInd];
			spectrumASXerrMB[curInd] = ASxeMB[curInd];
			spectrumASYMB[curInd] = ASyMB[curInd] / weightDiv->Eval(spectrumASXMB[curInd]);
			spectrumASYerrMB[curInd] = ASyeMB[curInd] / weightDiv->Eval(spectrumASXMB[curInd]);
		}
		for (curInd = 0;curInd < Int_t(sizeof(ASxHT1)/sizeof(ASxHT1[0])); curInd++) {
			spectrumASXHT1[curInd] = ASxHT1[curInd];
			spectrumASXerrHT1[curInd] = ASxeHT1[curInd];
			spectrumASYHT1[curInd] = ASyHT1[curInd] / weightDiv->Eval(spectrumASXHT1[curInd]);
			spectrumASYerrHT1[curInd] = ASyeHT1[curInd] / weightDiv->Eval(spectrumASXHT1[curInd]);
		}
		for (curInd = 0;curInd < Int_t(sizeof(ASxHT2)/sizeof(ASxHT2[0])); curInd++) {
			spectrumASXHT2[curInd] = ASxHT2[curInd];
			spectrumASXerrHT2[curInd] = ASxeHT2[curInd];
			spectrumASYHT2[curInd] = ASyHT2[curInd] / weightDiv->Eval(spectrumASXHT2[curInd]);
			spectrumASYerrHT2[curInd] = ASyeHT2[curInd] / weightDiv->Eval(spectrumASXHT2[curInd]);
		}
		TGraphErrors *hgrSpectrumASMB = new TGraphErrors(sizeof(ASxMB)/sizeof(ASxMB[0]), spectrumASXMB, spectrumASYMB, spectrumASXerrMB, spectrumASYerrMB);
		hgrSpectrumASMB->SetLineColor(kBlack);
		hgrSpectrumASMB->SetMarkerColor(hgrSpectrumASMB->GetLineColor());
		hgrSpectrumASMB->SetFillColor(hgrSpectrumASMB->GetLineColor());
		hgrSpectrumASMB->SetMarkerStyle(23);
		hgrSpectrumASMB->SetMarkerSize(1.3);
		hgrSpectrumASMB->SetHistogram(histDiv);
		hgrSpectrumASMB->Draw("P");
		TGraphErrors *hgrSpectrumASHT1 = new TGraphErrors(sizeof(ASxHT1)/sizeof(ASxHT1[0]), spectrumASXHT1, spectrumASYHT1, spectrumASXerrHT1, spectrumASYerrHT1);
		hgrSpectrumASHT1->SetLineColor(kBlue);
		hgrSpectrumASHT1->SetMarkerColor(hgrSpectrumASHT1->GetLineColor());
		hgrSpectrumASHT1->SetFillColor(hgrSpectrumASHT1->GetLineColor());
		hgrSpectrumASHT1->SetMarkerStyle(hgrSpectrumASMB->GetMarkerStyle());
		hgrSpectrumASHT1->SetMarkerSize(hgrSpectrumASMB->GetMarkerSize());
		hgrSpectrumASHT1->SetHistogram(histDiv);
		hgrSpectrumASHT1->Draw("P");
		TGraphErrors *hgrSpectrumASHT2 = new TGraphErrors(sizeof(ASxHT2)/sizeof(ASxHT2[0]), spectrumASXHT2, spectrumASYHT2, spectrumASXerrHT2, spectrumASYerrHT2);
		hgrSpectrumASHT2->SetLineColor(kRed);
		hgrSpectrumASHT2->SetMarkerColor(hgrSpectrumASHT2->GetLineColor());
		hgrSpectrumASHT2->SetFillColor(hgrSpectrumASHT2->GetLineColor());
		hgrSpectrumASHT2->SetMarkerStyle(hgrSpectrumASMB->GetMarkerStyle());
		hgrSpectrumASHT2->SetMarkerSize(hgrSpectrumASMB->GetMarkerSize());
		hgrSpectrumASHT2->SetHistogram(histDiv);
		hgrSpectrumASHT2->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumASMB, showStolpovskyPi0Legend, "P");
	}
	if (showPhenixPi0) {
		Float_t *spectrumPX = new Float_t[sizeof(Px)/sizeof(Px[0])];
		Float_t *spectrumPXerr = new Float_t[sizeof(Px)/sizeof(Px[0])];
		Float_t *spectrumPY = new Float_t[sizeof(Px)/sizeof(Px[0])];
		Float_t *spectrumPYerr = new Float_t[sizeof(Px)/sizeof(Px[0])];
		for (curInd = 0;curInd < Int_t(sizeof(Px)/sizeof(Px[0])); curInd++) {
			spectrumPX[curInd] = Px[curInd];
			spectrumPXerr[curInd] = Pxe[curInd];
			spectrumPY[curInd] = Py[curInd] / weightDiv->Eval(spectrumPX[curInd]);
			spectrumPYerr[curInd] = Pye[curInd] / weightDiv->Eval(spectrumPX[curInd]);
		}
		TGraphErrors *hgrSpectrumP = new TGraphErrors(sizeof(Px)/sizeof(Px[0]), spectrumPX, spectrumPY, spectrumPXerr, spectrumPYerr);
		hgrSpectrumP->SetLineColor(1);
		hgrSpectrumP->SetMarkerColor(hgrSpectrumP->GetLineColor());
		hgrSpectrumP->SetMarkerStyle(25);
		hgrSpectrumP->SetMarkerSize(1.1);
		hgrSpectrumP->SetFillColor(hgrSpectrumP->GetLineColor());
		hgrSpectrumP->SetHistogram(histDiv);
		hgrSpectrumP->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumP, showPhenixPi0Legend, "P");
	}
	if (showPhenixPi0PP) {
		Float_t *spectrumPppX = new Float_t[sizeof(Pppx)/sizeof(Pppx[0])];
		Float_t *spectrumPppXerr = new Float_t[sizeof(Pppx)/sizeof(Pppx[0])];
		Float_t *spectrumPppY = new Float_t[sizeof(Pppx)/sizeof(Pppx[0])];
		Float_t *spectrumPppYerr = new Float_t[sizeof(Pppx)/sizeof(Pppx[0])];
		for (curInd = 0;curInd < Int_t(sizeof(Pppx)/sizeof(Pppx[0])); curInd++) {
			spectrumPppX[curInd] = Pppx[curInd];
			spectrumPppXerr[curInd] = Pppxe[curInd];
			spectrumPppY[curInd] = Pppy[curInd] / weightDiv->Eval(spectrumPppX[curInd]);
			spectrumPppYerr[curInd] = Pppye[curInd] / weightDiv->Eval(spectrumPppX[curInd]);
		}
		TGraphErrors *hgrSpectrumPpp = new TGraphErrors(sizeof(Pppx)/sizeof(Pppx[0]), spectrumPppX, spectrumPppY, spectrumPppXerr, spectrumPppYerr);
		hgrSpectrumPpp->SetLineColor(1);
		hgrSpectrumPpp->SetMarkerColor(hgrSpectrumPpp->GetLineColor());
		hgrSpectrumPpp->SetMarkerStyle(25);
		hgrSpectrumPpp->SetMarkerSize(1.1);
		hgrSpectrumPpp->SetFillColor(hgrSpectrumPpp->GetLineColor());
		hgrSpectrumPpp->SetHistogram(histDiv);
		hgrSpectrumPpp->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumPpp, showPhenixPi0PPLegend, "P");
	}
	if (showMischkePi0) {
		Float_t *spectrumAX = new Float_t[sizeof(Ax)/sizeof(Ax[0])];
		Float_t *spectrumAXerr = new Float_t[sizeof(Ax)/sizeof(Ax[0])];
		Float_t *spectrumAY = new Float_t[sizeof(Ax)/sizeof(Ax[0])];
		Float_t *spectrumAYerr = new Float_t[sizeof(Ax)/sizeof(Ax[0])];
		for (curInd = 0;curInd < Int_t(sizeof(Ax)/sizeof(Ax[0])); curInd++) {
			spectrumAX[curInd] = Ax[curInd];
			spectrumAXerr[curInd] = Axe[curInd];
			spectrumAY[curInd] = Ay[curInd] / weightDiv->Eval(spectrumAX[curInd]);
			spectrumAYerr[curInd] = Aye[curInd] / weightDiv->Eval(spectrumAX[curInd]);
		}
		TGraphErrors *hgrSpectrumA = new TGraphErrors(sizeof(Ax)/sizeof(Ax[0]), spectrumAX, spectrumAY, spectrumAXerr, spectrumAYerr);
		hgrSpectrumA->SetLineColor(kGreen);
		hgrSpectrumA->SetFillColor(hgrSpectrumA->GetLineColor());
		hgrSpectrumA->SetMarkerColor(hgrSpectrumA->GetLineColor());
		hgrSpectrumA->SetMarkerStyle(22);
		hgrSpectrumA->SetMarkerSize(1.3);
		hgrSpectrumA->SetHistogram(histDiv);
		hgrSpectrumA->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumA, showMischkePi0Legend, "P");
	}
	if (showMartijnPi0dAu) {
		Float_t *spectrumMartijnXMB = new Float_t[sizeof(MartijnxMB)/sizeof(MartijnxMB[0])];
		Float_t *spectrumMartijnXerrMB = new Float_t[sizeof(MartijnxMB)/sizeof(MartijnxMB[0])];
		Float_t *spectrumMartijnYMB = new Float_t[sizeof(MartijnxMB)/sizeof(MartijnxMB[0])];
		Float_t *spectrumMartijnYerrMB = new Float_t[sizeof(MartijnxMB)/sizeof(MartijnxMB[0])];
		Float_t *spectrumMartijnXHT1 = new Float_t[sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0])];
		Float_t *spectrumMartijnXerrHT1 = new Float_t[sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0])];
		Float_t *spectrumMartijnYHT1 = new Float_t[sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0])];
		Float_t *spectrumMartijnYerrHT1 = new Float_t[sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0])];
		Float_t *spectrumMartijnXHT2 = new Float_t[sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0])];
		Float_t *spectrumMartijnXerrHT2 = new Float_t[sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0])];
		Float_t *spectrumMartijnYHT2 = new Float_t[sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0])];
		Float_t *spectrumMartijnYerrHT2 = new Float_t[sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0])];
		for (curInd = 0;curInd < Int_t(sizeof(MartijnxMB)/sizeof(MartijnxMB[0])); curInd++) {
			spectrumMartijnXMB[curInd] = MartijnxMB[curInd];
			spectrumMartijnXerrMB[curInd] = MartijnxeMB[curInd];
			spectrumMartijnYMB[curInd] = MartijnyMB[curInd] / weightDiv->Eval(spectrumMartijnXMB[curInd]);
			spectrumMartijnYerrMB[curInd] = MartijnyeMB[curInd] / weightDiv->Eval(spectrumMartijnXMB[curInd]);
		}
		for (curInd = 0;curInd < Int_t(sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0])); curInd++) {
			spectrumMartijnXHT1[curInd] = MartijnxHT1[curInd];
			spectrumMartijnXerrHT1[curInd] = MartijnxeHT1[curInd];
			spectrumMartijnYHT1[curInd] = MartijnyHT1[curInd] / weightDiv->Eval(spectrumMartijnXHT1[curInd]);
			spectrumMartijnYerrHT1[curInd] = MartijnyeHT1[curInd] / weightDiv->Eval(spectrumMartijnXHT1[curInd]);
		}
		for (curInd = 0;curInd < Int_t(sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0])); curInd++) {
			spectrumMartijnXHT2[curInd] = MartijnxHT2[curInd];
			spectrumMartijnXerrHT2[curInd] = MartijnxeHT2[curInd];
			spectrumMartijnYHT2[curInd] = MartijnyHT2[curInd] / weightDiv->Eval(spectrumMartijnXHT2[curInd]);
			spectrumMartijnYerrHT2[curInd] = MartijnyeHT2[curInd] / weightDiv->Eval(spectrumMartijnXHT2[curInd]);
		}
		TGraphErrors *hgrSpectrumMartijnMB = new TGraphErrors(sizeof(MartijnxMB)/sizeof(MartijnxMB[0]), spectrumMartijnXMB, spectrumMartijnYMB, spectrumMartijnXerrMB, spectrumMartijnYerrMB);
		hgrSpectrumMartijnMB->SetLineColor(6);
		hgrSpectrumMartijnMB->SetMarkerColor(hgrSpectrumMartijnMB->GetLineColor());
		hgrSpectrumMartijnMB->SetFillColor(hgrSpectrumMartijnMB->GetLineColor());
		hgrSpectrumMartijnMB->SetMarkerStyle(23);
		hgrSpectrumMartijnMB->SetMarkerSize(1.3);
		hgrSpectrumMartijnMB->SetHistogram(histDiv);
		hgrSpectrumMartijnMB->Draw("P");
		TGraphErrors *hgrSpectrumMartijnHT1 = new TGraphErrors(sizeof(MartijnxHT1)/sizeof(MartijnxHT1[0]), spectrumMartijnXHT1, spectrumMartijnYHT1, spectrumMartijnXerrHT1, spectrumMartijnYerrHT1);
		hgrSpectrumMartijnHT1->SetLineColor(6);
		hgrSpectrumMartijnHT1->SetMarkerColor(hgrSpectrumMartijnHT1->GetLineColor());
		hgrSpectrumMartijnHT1->SetFillColor(hgrSpectrumMartijnHT1->GetLineColor());
		hgrSpectrumMartijnHT1->SetMarkerStyle(hgrSpectrumMartijnMB->GetMarkerStyle());
		hgrSpectrumMartijnHT1->SetMarkerSize(hgrSpectrumMartijnMB->GetMarkerSize());
		hgrSpectrumMartijnHT1->SetHistogram(histDiv);
		hgrSpectrumMartijnHT1->Draw("P");
		TGraphErrors *hgrSpectrumMartijnHT2 = new TGraphErrors(sizeof(MartijnxHT2)/sizeof(MartijnxHT2[0]), spectrumMartijnXHT2, spectrumMartijnYHT2, spectrumMartijnXerrHT2, spectrumMartijnYerrHT2);
		hgrSpectrumMartijnHT2->SetLineColor(6);
		hgrSpectrumMartijnHT2->SetMarkerColor(hgrSpectrumMartijnHT2->GetLineColor());
		hgrSpectrumMartijnHT2->SetFillColor(hgrSpectrumMartijnHT2->GetLineColor());
		hgrSpectrumMartijnHT2->SetMarkerStyle(hgrSpectrumMartijnMB->GetMarkerStyle());
		hgrSpectrumMartijnHT2->SetMarkerSize(hgrSpectrumMartijnMB->GetMarkerSize());
		hgrSpectrumMartijnHT2->SetHistogram(histDiv);
		hgrSpectrumMartijnHT2->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumMartijnMB, showMartijnPi0dAuLegend, "P");
	}
	if (showSimonPi0PPCrossection) {
		Float_t *spectrumFSxsX = new Float_t[sizeof(FSxs_x)/sizeof(FSxs_x[0])];
		Float_t *spectrumFSxsXerr = new Float_t[sizeof(FSxs_x)/sizeof(FSxs_x[0])];
		Float_t *spectrumFSxsY = new Float_t[sizeof(FSxs_x)/sizeof(FSxs_x[0])];
		Float_t *spectrumFSxsYerr = new Float_t[sizeof(FSxs_x)/sizeof(FSxs_x[0])];
		for (curInd = 0;curInd < Int_t(sizeof(FSxs_x)/sizeof(FSxs_x[0])); curInd++) {
			spectrumFSxsX[curInd] = FSxs_x[curInd];
			spectrumFSxsXerr[curInd] = FSxs_xe[curInd];
			spectrumFSxsY[curInd] = FSxs_y[curInd] / weightDiv->Eval(spectrumFSxsX[curInd]);
			spectrumFSxsYerr[curInd] = FSxs_ye[curInd] / weightDiv->Eval(spectrumFSxsX[curInd]);
		}
		TGraphErrors *hgrSpectrumFSxs = new TGraphErrors(sizeof(FSxs_x)/sizeof(FSxs_x[0]), spectrumFSxsX, spectrumFSxsY, spectrumFSxsXerr, spectrumFSxsYerr);
		hgrSpectrumFSxs->SetLineColor(kGreen);
		hgrSpectrumFSxs->SetFillColor(hgrSpectrumFSxs->GetLineColor());
		hgrSpectrumFSxs->SetMarkerColor(hgrSpectrumFSxs->GetLineColor());
		hgrSpectrumFSxs->SetMarkerStyle(22);
		hgrSpectrumFSxs->SetMarkerSize(1.3);
		hgrSpectrumFSxs->SetHistogram(histDiv);
		hgrSpectrumFSxs->Draw("P");
		if (legend) legend->AddEntry(hgrSpectrumFSxs, showSimonPi0PPCrossectionLegend, "P");
	}
	if (legend) legend->Draw();
}

void showResultsOthersDAuRcp(TH1 *histRcp, TLegend *legendRcp
        , const TAnalysisResults &resultsDAuMostCentral
        , const TAnalysisResults &resultsDAuMostPeripheral
) {
	if (!histRcp) return;
	histRcp->SetMinimum(0);
	histRcp->SetMaximum(1.5);
	Float_t Xmin = histRcp->GetXaxis()->GetBinLowEdge(1);
	Float_t Xmax = histRcp->GetXaxis()->GetBinUpEdge(histRcp->GetXaxis()->GetNbins());
	TF1 *fRcp1 = new TF1("fRcp1", "1 + [0]", Xmin, Xmax);
	fRcp1->SetLineWidth(1);
	fRcp1->SetParameter(0, 0);
	fRcp1->DrawCopy("SAME");
	Float_t err1 = sqrt(((resultsDAuMostCentral.NBinaryCollisionsErr / resultsDAuMostCentral.NBinaryCollisions) * (resultsDAuMostCentral.NBinaryCollisionsErr / resultsDAuMostCentral.NBinaryCollisions)) + ((resultsDAuMostPeripheral.NBinaryCollisionsErr / resultsDAuMostPeripheral.NBinaryCollisions) * (resultsDAuMostPeripheral.NBinaryCollisionsErr / resultsDAuMostPeripheral.NBinaryCollisions)));
	Float_t *systNbinX = new Float_t[1];
	Float_t *systNbinXerr = new Float_t[1];
	Float_t *systNbinY = new Float_t[1];
	Float_t *systNbinYerr = new Float_t[1];
	systNbinX[0] = 0.25;
	systNbinXerr[0] = 0.25;
	systNbinY[0] = 1.0;
	systNbinYerr[0] = err1;
	TGraphErrors *systNbinGr = new TGraphErrors(1, systNbinX, systNbinY, systNbinXerr, systNbinYerr);
	systNbinGr->SetFillColor(11);
	systNbinGr->SetLineColor(11);
	systNbinGr->SetMarkerColor(11);
	systNbinGr->SetHistogram(histRcp);
	systNbinGr->Draw("2");
	TString textStr;
	textStr = Form("#splitline{Glauber:}{#splitline{<N_{bin}>_{%3.0f-%3.0f%%} = %3.1f #pm %3.1f}{#splitline{<N_{bin}>_{%3.0f-%3.0f%%} = %3.1f #pm %3.1f}{%3.0f%% error on ratio}}}"
		, resultsDAuMostCentral.PercentCentralityLow, resultsDAuMostCentral.PercentCentralityHigh, resultsDAuMostCentral.NBinaryCollisions, resultsDAuMostCentral.NBinaryCollisionsErr
		, resultsDAuMostPeripheral.PercentCentralityLow, resultsDAuMostPeripheral.PercentCentralityHigh, resultsDAuMostPeripheral.NBinaryCollisions, resultsDAuMostPeripheral.NBinaryCollisionsErr
		, 100.0 * err1);
	TLatex *text = new TLatex(10.0, 0.8, textStr.Data());
	text->SetTextSize(0.03);
	text->Draw();
	if (showStarRcp) {
		TGraphErrors *hgrSRcp1 = new TGraphErrors(sizeof(SRcpx)/sizeof(SRcpx[0]), SRcpx, SRcp1y, SRcpxe, SRcp1ye);
		hgrSRcp1->SetLineColor(kBlue);
		hgrSRcp1->SetMarkerColor(hgrSRcp1->GetLineColor());
		hgrSRcp1->SetFillColor(hgrSRcp1->GetLineColor());
		hgrSRcp1->SetMarkerStyle(25);
		hgrSRcp1->SetMarkerSize(1.0);
		hgrSRcp1->SetHistogram(histRcp);
		hgrSRcp1->Draw("P");
		if (legendRcp) legendRcp->AddEntry(hgrSRcp1, showStarRcpLegend1, "P");
		TGraphErrors *hgrSRcp2 = new TGraphErrors(sizeof(SRcpx)/sizeof(SRcpx[0]), SRcpx, SRcp2y, SRcpxe, SRcp2ye);
		hgrSRcp2->SetLineColor(hgrSRcp1->GetLineColor());
		hgrSRcp2->SetMarkerColor(hgrSRcp1->GetLineColor());
		hgrSRcp2->SetFillColor(hgrSRcp1->GetLineColor());
		hgrSRcp2->SetMarkerStyle(26);
		hgrSRcp2->SetMarkerSize(1.0);
		hgrSRcp2->SetHistogram(histRcp);
		hgrSRcp2->Draw("P");
		if (legendRcp) legendRcp->AddEntry(hgrSRcp2, showStarRcpLegend2, "P");
	}
	if (showPhenixPi0dA) {
		TGraphErrors *hgrPdA = new TGraphErrors(sizeof(PdAx)/sizeof(PdAx[0]), PdAx, PdAy, PdAxe, PdAye);
		hgrPdA->SetLineColor(1);
		hgrPdA->SetMarkerColor(hgrPdA->GetLineColor());
		hgrPdA->SetFillColor(hgrPdA->GetLineColor());
		hgrPdA->SetMarkerStyle(25);
		hgrPdA->SetMarkerSize(1.0);
		hgrPdA->SetHistogram(histRcp);
		hgrPdA->Draw("P");
		if (legendRcp) legendRcp->AddEntry(hgrPdA, showPhenixPi0dALegend, "P");
	}
	if (legendRcp) legendRcp->Draw();
}

void showPointsRcp(const Char_t *name, const Char_t *title
	, const TDataProcessorPool *poolMostcentral, const TDataProcessorPool *poolMostperipheral
	, Bool_t chargedPoints, Bool_t allPoints
	, Bool_t showSpectrum, Bool_t show, Bool_t print
	, Float_t binSizeMB, Float_t binSizeHT1, Float_t binSizeHT2
	, Float_t minPtMB, Float_t minPtHT1, Float_t minPtHT2
	, Float_t maxPtMB, Float_t maxPtHT1, Float_t maxPtHT2
) {
/*
	GET_DATA(pointDataProcessorMBMostcentral, const TPointDataProcessor, poolMostcentral, (allPoints ? "pointAllMB" : (chargedPoints ? "pointChargedMB" : "pointMB")));
	GET_DATA(pointDataProcessorHT1Mostcentral, const TPointDataProcessor, poolMostcentral, (allPoints ? "pointAllHT1" : (chargedPoints ? "pointChargedHT1" : "pointHT1")));
	GET_DATA(pointDataProcessorHT2Mostcentral, const TPointDataProcessor, poolMostcentral, (allPoints ? "pointAllHT2" : (chargedPoints ? "pointChargedHT2" : "pointHT2")));
	GET_DATA(pointDataProcessorMBMostperipheral, const TPointDataProcessor, poolMostperipheral, (allPoints ? "pointAllMB" : (chargedPoints ? "pointChargedMB" : "pointMB")));
	GET_DATA(pointDataProcessorHT1Mostperipheral, const TPointDataProcessor, poolMostperipheral, (allPoints ? "pointAllHT1" : (chargedPoints ? "pointChargedHT1" : "pointHT1")));
	GET_DATA(pointDataProcessorHT2Mostperipheral, const TPointDataProcessor, poolMostperipheral, (allPoints ? "pointAllHT2" : (chargedPoints ? "pointChargedHT2" : "pointHT2")));
	GET_DATA(eventDataProcessorMBMostcentral, const TEventDataProcessor, poolMostcentral, "eventMB");
	GET_DATA(eventDataProcessorHT1Mostcentral, const TEventDataProcessor, poolMostcentral, "eventHT1");
	GET_DATA(eventDataProcessorHT2Mostcentral, const TEventDataProcessor, poolMostcentral, "eventHT2");
	GET_DATA(eventDataProcessorMBMostperipheral, const TEventDataProcessor, poolMostperipheral, "eventMB");
	GET_DATA(eventDataProcessorHT1Mostperipheral, const TEventDataProcessor, poolMostperipheral, "eventHT1");
	GET_DATA(eventDataProcessorHT2Mostperipheral, const TEventDataProcessor, poolMostperipheral, "eventHT2");
	const TH1F *pointsMBc = pointDataProcessorMBMostcentral ? pointDataProcessorMBMostcentral->getPt() : 0;
	const TH1F *pointsHT1c = pointDataProcessorHT1Mostcentral ? pointDataProcessorHT1Mostcentral->getPt() : 0;
	const TH1F *pointsHT2c = pointDataProcessorHT2Mostcentral ? pointDataProcessorHT2Mostcentral->getPt() : 0;
	const TH1F *pointsMBp = pointDataProcessorMBMostperipheral ? pointDataProcessorMBMostperipheral->getPt() : 0;
	const TH1F *pointsHT1p = pointDataProcessorHT1Mostperipheral ? pointDataProcessorHT1Mostperipheral->getPt() : 0;
	const TH1F *pointsHT2p = pointDataProcessorHT2Mostperipheral ? pointDataProcessorHT2Mostperipheral->getPt() : 0;
	TH1F *pMBc = pointsMBc ? (new TH1F(*pointsMBc)) : 0;
	TH1F *pHT1c = pointsHT1c ? (new TH1F(*pointsHT1c)) : 0;
	TH1F *pHT2c = pointsHT2c ? (new TH1F(*pointsHT2c)) : 0;
	TH1F *pMBp = pointsMBp ? (new TH1F(*pointsMBp)) : 0;
	TH1F *pHT1p = pointsHT1p ? (new TH1F(*pointsHT1p)) : 0;
	TH1F *pHT2p = pointsHT2p ? (new TH1F(*pointsHT2p)) : 0;
	if (pMBc) if (pMBc->GetSumw2N() == 0) pMBc->Sumw2();
	if (pHT1c) if (pHT1c->GetSumw2N() == 0) pHT1c->Sumw2();
	if (pHT2c) if (pHT2c->GetSumw2N() == 0) pHT2c->Sumw2();
	if (pMBp) if (pMBp->GetSumw2N() == 0) pMBp->Sumw2();
	if (pHT1p) if (pHT1p->GetSumw2N() == 0) pHT1p->Sumw2();
	if (pHT2p) if (pHT2p->GetSumw2N() == 0) pHT2p->Sumw2();
	Float_t PSHT1c = 1.0, PSHT2c = 1.0, PSHT1p = 1.0, PSHT2p = 1.0;

	Float_t PSHT1c_sim, PSHT2c_sim;
	calculatePSFromSim(name, title, poolMostcentral, false, false, PSHT1c_sim, PSHT2c_sim);
	//Float_t PSHT1c_point, PSHT2c_point;
	//calculatePSFromPoints(name, title, poolMostcentral, false, true, false, true, 4.0, 5.0, 5.0, 10.0, 0.3, 4.0, 5.0, 5.0, 10.0, false, false, PSHT1c_point, PSHT2c_point);
	//PSHT1c = (PSHT1c_sim + PSHT1c_point) / 2.0;
	//PSHT2c = (PSHT2c_sim + PSHT2c_point) / 2.0;
	PSHT1c = PSHT1c_sim;
	PSHT2c = PSHT2c_sim;

	Float_t PSHT1p_sim, PSHT2p_sim;
	calculatePSFromSim(name, title, poolMostperipheral, false, false, PSHT1p_sim, PSHT2p_sim);
	//Float_t PSHT1p_point, PSHT2p_point;
	//calculatePSFromPoints(name, title, poolMostperipheral, false, true, false, true, 4.0, 5.0, 5.0, 10.0, 0.3, 4.0, 5.0, 5.0, 10.0, false, false, PSHT1p_point, PSHT2p_point);
	//PSHT1p = (PSHT1p_sim + PSHT1p_point) / 2.0;
	//PSHT2p = (PSHT2p_sim + PSHT2p_point) / 2.0;
	PSHT1p = PSHT1p_sim;
	PSHT2p = PSHT2p_sim;

	if (pMBc && eventDataProcessorMBMostcentral) {
		if (eventDataProcessorMBMostcentral->numPassedAllCuts != 0) pMBc->Scale(1.0 / eventDataProcessorMBMostcentral->numPassedAllCuts);
	}
	if (pHT1c && eventDataProcessorHT1Mostcentral) {
		if (eventDataProcessorHT1Mostcentral->numPassedAllCuts != 0) pHT1c->Scale(1.0 / eventDataProcessorHT1Mostcentral->numPassedAllCuts);
	}
	if (pHT2c && eventDataProcessorHT2Mostcentral) {
		if (eventDataProcessorHT2Mostcentral->numPassedAllCuts != 0) pHT2c->Scale(1.0 / eventDataProcessorHT2Mostcentral->numPassedAllCuts);
	}
	if (pMBp && eventDataProcessorMBMostperipheral) {
		if (eventDataProcessorMBMostperipheral->numPassedAllCuts != 0) pMBp->Scale(1.0 / eventDataProcessorMBMostperipheral->numPassedAllCuts);
	}
	if (pHT1p && eventDataProcessorHT1Mostperipheral) {
		if (eventDataProcessorHT1Mostperipheral->numPassedAllCuts != 0) pHT1p->Scale(1.0 / eventDataProcessorHT1Mostperipheral->numPassedAllCuts);
	}
	if (pHT2p && eventDataProcessorHT2Mostperipheral) {
		if (eventDataProcessorHT2Mostperipheral->numPassedAllCuts != 0) pHT2p->Scale(1.0 / eventDataProcessorHT2Mostperipheral->numPassedAllCuts);
	}
	if (pHT1c && (PSHT1c != 0)) pHT1c->Scale(1.0 / PSHT1c);
	if (pHT2c && (PSHT2c != 0)) pHT2c->Scale(1.0 / PSHT2c);
	if (pHT1p && (PSHT1p != 0)) pHT1p->Scale(1.0 / PSHT1p);
	if (pHT2p && (PSHT2p != 0)) pHT2p->Scale(1.0 / PSHT2p);
	if (showSpectrum && (pMBc || pHT1c || pHT2c)) {
		TString nameStr = TString(name) + "_points_pt_mostcentral";
		TString titleStr = TString(title) + ", EMC points p_{T} distribution per event - most central events";
		TCanvas *c = new TCanvas(nameStr.Data(), titleStr.Data());
		if (c) c->cd();
		if (pMBc) {
			pMBc->SetMarkerColor(kBlack);
			pMBc->SetLineColor(pMBc->GetMarkerColor());
			pMBc->SetFillColor(pMBc->GetMarkerColor());
			pMBc->SetTitle(titleStr.Data());
			pMBc->Draw();
		}
		if (pHT1c) {
			pHT1c->SetMarkerColor(kBlue);
			pHT1c->SetLineColor(pTH1c->GetMarkerColor());
			pHT1c->SetFillColor(pHT1c->GetMarkerColor());
			pHT1c->Draw("SAME");
		}
		if (pHT2c) {
			pHT2c->SetMarkerColor(kRed);
			pHT2c->SetLineColor(pHT2c->GetMarkerColor());
			pHT2c->SetFillColor(pHT2c->GetMarkerColor());
			pHT2c->Draw("SAME");
		}
		//if (pMBc && pHT2c) pMBc->GetYaxis()->SetRangeUser(pHT2c->GetMinimum(), pMBc->GetMaximum());
		if (pMBc && pHT2c) pMBc->GetYaxis()->SetRangeUser(1.0, 1.0e-9);
		gPad->SetLogy();
	}
	if (showSpectrum && (pMBp || pHT1p || pHT2p)) {
		TString nameStr = TString(name) + "_points_pt_mostperipheral";
		TString titleStr = TString(title) + ", EMC points p_{T} distribution per event - most peripheral events";
		TCanvas *c = new TCanvas(nameStr.Data(), titleStr.Data());
		if (c) c->cd();
		if (pMBp) {
			pMBp->SetMarkerColor(kBlack);
			pMBp->SetLineColor(pMBp->GetMarkerColor());
			pMBp->SetFillColor(pMBp->GetMarkerColor());
			pMBp->SetTitle(titleStr.Data());
			pMBp->Draw();
		}
		if (pHT1p) {
			pHT1p->SetMarkerColor(kBlue);
			pHT1p->SetLineColor(pHT1p->GetMarkerColor());
			pHT1p->SetFillColor(pHT1p->GetMarkerColor());
			pHT1p->Draw("SAME");
		}
		if (pHT2p) {
			pHT2p->SetMarkerColor(kRed);
			pHT2p->SetLineColor(pHT2p->GetMarkerColor());
			pHT2p->SetFillColor(pHT2p->GetMarkerColor());
			pHT2p->Draw("SAME");
		}
		//if (pMBp && pHT2p) pMBp->GetYaxis()->SetRangeUser(pHT2p->GetMinimum(), pMBp->GetMaximum());
		if (pMBp && pHT2p) pMBp->GetYaxis()->SetRangeUser(1.0, 1.0e-9);
		gPad->SetLogy();
	}
	if (show && ((pMBc && pMBp) || (pHT1c && pHT1p) || (pHT2c && pHT2p))) {
		TString nameStr = TString(name) + "_points_pt_Rcp";
		TString titleCentral = Form(", %3.0f-%3.0f%%/%3.0f-%3.0f%%", PercentMostCentralLow, PercentMostCentralHigh, PercentMostPeripheralLow, PercentMostPeripheralHigh);
		TString titleStr = TString(title) + ", EMC points R_{CP}" + titleCentral;
		TCanvas *c = new TCanvas(nameStr.Data(), titleStr.Data());
		TH1F *pMBcc = pMBc ? (new TH1F(*pMBc)) : 0;
		TH1F *pHT1cc = pHT1c ? (new TH1F(*pHT1c)) : 0;
		TH1F *pHT2cc = pHT2c ? (new TH1F(*pHT2c)) : 0;
		TH1F *pMBpp = pMBp ? (new TH1F(*pMBp)) : 0;
		TH1F *pHT1pp = pHT1p ? (new TH1F(*pHT1p)) : 0;
		TH1F *pHT2pp = pHT2p ? (new TH1F(*pHT2p)) : 0;
		if (c) c->cd();
		if (pMBcc && pMBpp) {
			if (pMBcc->GetBinWidth(1) < binSizeMB) pMBcc->Rebin(Int_t((binSizeMB / pMBcc->GetBinWidth(1)) + 1));
			if (pMBpp->GetBinWidth(1) < binSizeMB) pMBpp->Rebin(Int_t((binSizeMB / pMBpp->GetBinWidth(1)) + 1));
			pMBcc->Divide(pMBpp);
			pMBcc->Scale(NBinaryDAuMostPeripheral / NBinaryDAuMostCentral);
			Int_t mbin = 0;
			while (pMBcc->GetBinCenter(mbin) < minPtMB) {pMBcc->SetBinContent(mbin, 0);pMBcc->SetBinError(mbin, 0);mbin++;}
			mbin = pMBcc->GetNbinsX();
			while (pMBcc->GetBinCenter(mbin) > maxPtMB) {pMBcc->SetBinContent(mbin, 0);pMBcc->SetBinError(mbin, 0);mbin--;}
			pMBcc->SetMarkerSize(1.2);
			pMBcc->SetMarkerStyle(8);
			pMBcc->SetMarkerColor(kBlack);
			pMBcc->SetLineColor(pMBcc->GetMarkerColor());
			pMBcc->SetFillColor(pMBcc->GetMarkerColor());
			pMBcc->SetTitle(titleStr.Data());
			pMBcc->Draw("P");
		}
		if (pHT1cc && pHT1pp) {
			if (pHT1cc->GetBinWidth(1) < binSizeHT1) pHT1cc->Rebin(Int_t((binSizeHT1 / pHT1cc->GetBinWidth(1)) + 1));
			if (pHT1pp->GetBinWidth(1) < binSizeHT1) pHT1pp->Rebin(Int_t((binSizeHT1 / pHT1pp->GetBinWidth(1)) + 1));
			pHT1cc->Divide(pHT1pp);
			pHT1cc->Scale(NBinaryDAuMostPeripheral / NBinaryDAuMostCentral);
			Int_t mbin = 0;
			while (pHT1cc->GetBinCenter(mbin) < minPtHT1) {pHT1cc->SetBinContent(mbin, 0);pHT1cc->SetBinError(mbin, 0);mbin++;}
			mbin = pHT1cc->GetNbinsX();
			while (pHT1cc->GetBinCenter(mbin) > maxPtHT1) {pHT1cc->SetBinContent(mbin, 0);pHT1cc->SetBinError(mbin, 0);mbin--;}
			pHT1cc->SetMarkerColor(kBlue);
			pHT1cc->SetLineColor(pHT1cc->GetMarkerColor());
			pHT1cc->SetFillColor(pHT1cc->GetMarkerColor());
			pHT1cc->SetMarkerStyle(8);
			pHT1cc->SetMarkerSize(1.2);
			pHT1cc->Draw("SAME P");
		}
		if (pHT2cc && pHT2pp) {
			if (pHT2cc->GetBinWidth(1) < binSizeHT2) pHT2cc->Rebin(Int_t((binSizeHT2 / pHT2cc->GetBinWidth(1)) + 1));
			if (pHT2pp->GetBinWidth(1) < binSizeHT2) pHT2pp->Rebin(Int_t((binSizeHT2 / pHT2pp->GetBinWidth(1)) + 1));
			pHT2cc->Divide(pHT2pp);
			pHT2cc->Scale(NBinaryDAuMostPeripheral / NBinaryDAuMostCentral);
			Int_t mbin = 0;
			while (pHT2cc->GetBinCenter(mbin) < minPtHT2) {pHT2cc->SetBinContent(mbin, 0);pHT2cc->SetBinError(mbin, 0);mbin++;}
			mbin = pHT2cc->GetNbinsX();
			while (pHT2cc->GetBinCenter(mbin) > maxPtHT2) {pHT2cc->SetBinContent(mbin, 0);pHT2cc->SetBinError(mbin, 0);mbin--;}
			pHT2cc->SetMarkerColor(kRed);
			pHT2cc->SetLineColor(pHT2cc->GetMarkerColor());
			pHT2cc->SetFillColor(pHT2cc->GetMarkerColor());
			pHT2cc->SetMarkerStyle(8);
			pHT2cc->SetMarkerSize(1.2);
			pHT2cc->Draw("SAME P");
		}
		TLegend *legend = new TLegend(0.6, 0.8, 1.0, 1.0, "Legend");
		if (legend && pMBcc) legend->AddEntry(pMBcc, "MinBias", "P");
		if (legend && pHT1cc) legend->AddEntry(pHT1cc, "HighTower-1", "P");
		if (legend && pHT2cc) legend->AddEntry(pHT2cc, "HighTower-2", "P");
		if (legend) legend->Draw();
		showResultsOthersRcp(pMBcc, legend, 0, 20);
		if (pMBcc) pMBcc->SetMinimum(0);
		if (pMBcc) pMBcc->SetMaximum(1.5);
	}
*/
}

TH1 *ProjectX3(const TH3 *h, Int_t xLowBin, Int_t xHighBin, Int_t yLowBin, Int_t yHighBin, Int_t zLowBin, Int_t zHighBin) {
	TH1F *r = 0;
	if (h) {
		TString name = h->GetName();
		name += "_px";
		TString title = h->GetTitle();
		title += "_px";
		r = new TH1F(name.Data(), title.Data(), xHighBin - xLowBin, h->GetXaxis()->GetBinLowEdge(xLowBin), h->GetXaxis()->GetBinLowEdge(xHighBin));
		if (r) {
			r->Sumw2();
			for (Int_t ix = xLowBin;ix < xHighBin;ix++) {
				for (Int_t iy = yLowBin;iy <= yHighBin;iy++) {
					for (Int_t iz = zLowBin;iz <= zHighBin;iz++) {
						Int_t xbin = ix - xLowBin + 1;
						Float_t prevC = r->GetBinContent(xbin);
						Float_t prevE = r->GetBinError(xbin);
						Float_t nextC = h->GetBinContent(ix, iy, iz);
						Float_t nextE = h->GetBinError(ix, iy, iz);
						Float_t newC = prevC + nextC;
						Float_t newE = sqrt((prevE * prevE) + (nextE * nextE));
						r->SetBinContent(xbin, newC);
						r->SetBinError(xbin, newE);
					}
				}
			}
		}
	}
	return r;
}

void showPointsRcpEtaPhiCoord(const Char_t *name, const Char_t *title
	, const TDataProcessorPool *poolMostcentral, const TDataProcessorPool *poolMostperipheral
	, Bool_t chargedPoints, Bool_t allPoints
	, Bool_t show, Bool_t print
	, Float_t binSizeMBPt, Float_t binSizeHT1Pt, Float_t binSizeHT2Pt
	, Float_t binSizeEta, Float_t binSizePhi
	, Float_t minPtMB, Float_t minPtHT1, Float_t minPtHT2
	, Float_t maxPtMB, Float_t maxPtHT1, Float_t maxPtHT2
	, Float_t integralLow, Float_t maxError
) {
/*
	GET_DATA(pointDataProcessorMBMostcentral, const TPointDataProcessor, poolMostcentral, (allPoints ? "pointAllMB" : (chargedPoints ? "pointChargedMB" : "pointMB")));
	GET_DATA(pointDataProcessorHT1Mostcentral, const TPointDataProcessor, poolMostcentral, (allPoints ? "pointAllHT1" : (chargedPoints ? "pointChargedHT1" : "pointHT1")));
	GET_DATA(pointDataProcessorHT2Mostcentral, const TPointDataProcessor, poolMostcentral, (allPoints ? "pointAllHT2" : (chargedPoints ? "pointChargedHT2" : "pointHT2")));
	GET_DATA(pointDataProcessorMBMostperipheral, const TPointDataProcessor, poolMostperipheral, (allPoints ? "pointAllMB" : (chargedPoints ? "pointChargedMB" : "pointMB")));
	GET_DATA(pointDataProcessorHT1Mostperipheral, const TPointDataProcessor, poolMostperipheral, (allPoints ? "pointAllHT1" : (chargedPoints ? "pointChargedHT1" : "pointHT1")));
	GET_DATA(pointDataProcessorHT2Mostperipheral, const TPointDataProcessor, poolMostperipheral, (allPoints ? "pointAllHT2" : (chargedPoints ? "pointChargedHT2" : "pointHT2")));
	GET_DATA(eventDataProcessorMBMostcentral, const TEventDataProcessor, poolMostcentral, "eventMB");
	GET_DATA(eventDataProcessorHT1Mostcentral, const TEventDataProcessor, poolMostcentral, "eventHT1");
	GET_DATA(eventDataProcessorHT2Mostcentral, const TEventDataProcessor, poolMostcentral, "eventHT2");
	GET_DATA(eventDataProcessorMBMostperipheral, const TEventDataProcessor, poolMostperipheral, "eventMB");
	GET_DATA(eventDataProcessorHT1Mostperipheral, const TEventDataProcessor, poolMostperipheral, "eventHT1");
	GET_DATA(eventDataProcessorHT2Mostperipheral, const TEventDataProcessor, poolMostperipheral, "eventHT2");
	const TH3F *pointsMBc = pointDataProcessorMBMostcentral ? pointDataProcessorMBMostcentral->getPtEtaPhiCoord() : 0;
	const TH3F *pointsHT1c = pointDataProcessorHT1Mostcentral ? pointDataProcessorHT1Mostcentral->getPtEtaPhiCoord() : 0;
	const TH3F *pointsHT2c = pointDataProcessorHT2Mostcentral ? pointDataProcessorHT2Mostcentral->getPtEtaPhiCoord() : 0;
	const TH3F *pointsMBp = pointDataProcessorMBMostperipheral ? pointDataProcessorMBMostperipheral->getPtEtaPhiCoord() : 0;
	const TH3F *pointsHT1p = pointDataProcessorHT1Mostperipheral ? pointDataProcessorHT1Mostperipheral->getPtEtaPhiCoord() : 0;
	const TH3F *pointsHT2p = pointDataProcessorHT2Mostperipheral ? pointDataProcessorHT2Mostperipheral->getPtEtaPhiCoord() : 0;
	TH3F *pMBc = pointsMBc ? (new TH3F(*pointsMBc)) : 0;
	TH3F *pHT1c = pointsHT1c ? (new TH3F(*pointsHT1c)) : 0;
	TH3F *pHT2c = pointsHT2c ? (new TH3F(*pointsHT2c)) : 0;
	TH3F *pMBp = pointsMBp ? (new TH3F(*pointsMBp)) : 0;
	TH3F *pHT1p = pointsHT1p ? (new TH3F(*pointsHT1p)) : 0;
	TH3F *pHT2p = pointsHT2p ? (new TH3F(*pointsHT2p)) : 0;
	if (pMBc) if (pMBc->GetSumw2N() == 0) pMBc->Sumw2();
	if (pHT1c) if (pHT1c->GetSumw2N() == 0) pHT1c->Sumw2();
	if (pHT2c) if (pHT2c->GetSumw2N() == 0) pHT2c->Sumw2();
	if (pMBp) if (pMBp->GetSumw2N() == 0) pMBp->Sumw2();
	if (pHT1p) if (pHT1p->GetSumw2N() == 0) pHT1p->Sumw2();
	if (pHT2p) if (pHT2p->GetSumw2N() == 0) pHT2p->Sumw2();
	Float_t PSHT1c = 1.0, PSHT2c = 1.0, PSHT1p = 1.0, PSHT2p = 1.0;
	Float_t PSHT1c_sim, PSHT2c_sim;
	calculatePSFromSim(name, title, poolMostcentral, false, false, PSHT1c_sim, PSHT2c_sim);
	//Float_t PSHT1c_point, PSHT2c_point;
	//calculatePSFromPoints(name, title, poolMostcentral, false, true, false, true, 4.0, 5.0, 5.0, 10.0, 0.3, 4.0, 5.0, 5.0, 10.0, false, false, PSHT1c_point, PSHT2c_point);
	//PSHT1c = (PSHT1c_sim + PSHT1c_point) / 2.0;
	//PSHT2c = (PSHT2c_sim + PSHT2c_point) / 2.0;
	PSHT1c = PSHT1c_sim;
	PSHT2c = PSHT2c_sim;
	Float_t PSHT1p_sim, PSHT2p_sim;
	calculatePSFromSim(name, title, poolMostperipheral, false, false, PSHT1p_sim, PSHT2p_sim);
	//Float_t PSHT1p_point, PSHT2p_point;
	//calculatePSFromPoints(name, title, poolMostperipheral, false, true, false, true, 4.0, 5.0, 5.0, 10.0, 0.3, 4.0, 5.0, 5.0, 10.0, false, false, PSHT1p_point, PSHT2p_point);
	//PSHT1p = (PSHT1p_sim + PSHT1p_point) / 2.0;
	//PSHT2p = (PSHT2p_sim + PSHT2p_point) / 2.0;
	PSHT1p = PSHT1p_sim;
	PSHT2p = PSHT2p_sim;

	if (pMBc && eventDataProcessorMBMostcentral) {
		if (eventDataProcessorMBMostcentral->numPassedAllCuts != 0) pMBc->Scale(1.0 / eventDataProcessorMBMostcentral->numPassedAllCuts);
	}
	if (pHT1c && eventDataProcessorHT1Mostcentral) {
		if (eventDataProcessorHT1Mostcentral->numPassedAllCuts != 0) pHT1c->Scale(1.0 / eventDataProcessorHT1Mostcentral->numPassedAllCuts);
	}
	if (pHT2c && eventDataProcessorHT2Mostcentral) {
		if (eventDataProcessorHT2Mostcentral->numPassedAllCuts != 0) pHT2c->Scale(1.0 / eventDataProcessorHT2Mostcentral->numPassedAllCuts);
	}
	if (pMBp && eventDataProcessorMBMostperipheral) {
		if (eventDataProcessorMBMostperipheral->numPassedAllCuts != 0) pMBp->Scale(1.0 / eventDataProcessorMBMostperipheral->numPassedAllCuts);
	}
	if (pHT1p && eventDataProcessorHT1Mostperipheral) {
		if (eventDataProcessorHT1Mostperipheral->numPassedAllCuts != 0) pHT1p->Scale(1.0 / eventDataProcessorHT1Mostperipheral->numPassedAllCuts);
	}
	if (pHT2p && eventDataProcessorHT2Mostperipheral) {
		if (eventDataProcessorHT2Mostperipheral->numPassedAllCuts != 0) pHT2p->Scale(1.0 / eventDataProcessorHT2Mostperipheral->numPassedAllCuts);
	}
	if (pHT1c && (PSHT1c != 0)) pHT1c->Scale(1.0 / PSHT1c);
	if (pHT2c && (PSHT2c != 0)) pHT2c->Scale(1.0 / PSHT2c);
	if (pHT1p && (PSHT1p != 0)) pHT1p->Scale(1.0 / PSHT1p);
	if (pHT2p && (PSHT2p != 0)) pHT2p->Scale(1.0 / PSHT2p);
	if (show && ((pMBc && pMBp) || (pHT1c && pHT1p) || (pHT2c && pHT2p))) {
		Int_t numEtaBins = Int_t(Float_t(pMBc->GetYaxis()->GetXmax() - pMBc->GetYaxis()->GetXmin() + 0.0001) / binSizeEta);
		Int_t numPhiBins = Int_t(Float_t(pMBc->GetZaxis()->GetXmax() - pMBc->GetZaxis()->GetXmin() + 0.0001) / binSizePhi);
		if (numEtaBins <= 0) numEtaBins = 1;
		if (numPhiBins <= 0) numPhiBins = 1;
		TString nameStr = TString(name) + "_points_pt_eta_phi_Rcp";
		TString titleStr = TString(title) + ", EMC points R_{CP} vs. #eta and #phi";
		TCanvas *c = new TCanvas(nameStr.Data(), titleStr.Data());
                Int_t nyp = Int_t(sqrt(Float_t(numEtaBins * numPhiBins)));
                Int_t nxp = Int_t((Float_t(numEtaBins * numPhiBins) - 0.0001) / Float_t(nyp)) + 1;
		if (c) c->Divide(nxp, nyp, 0.001, 0.001);
		TString nameEtaPhiStr = TString(name) + "_points_pt_eta_phi_Rcp_mean_histo";
		TString titleEtaPhiStr = Form("%s, EMC points R_{CP} (average over p_{T} > %1.1f) vs. #eta and #phi;#eta;#phi", title, integralLow);
		TH2F *hEtaPhi = new TH2F(nameEtaPhiStr.Data(), titleEtaPhiStr.Data(), numEtaBins, pMBc->GetYaxis()->GetXmin(), pMBc->GetYaxis()->GetXmax(), numPhiBins, pMBc->GetZaxis()->GetXmin(), pMBc->GetZaxis()->GetXmax());
		for (Int_t ibin = 1;ibin <= numEtaBins * numPhiBins;ibin++) {
		    if (c) c->cd(ibin);
		    Int_t etaBin = (ibin - 1) / numPhiBins;
		    Int_t phiBin = (ibin - 1) - (etaBin * numPhiBins);
		    Float_t etaLow = pMBc->GetYaxis()->GetXmin() + (etaBin * binSizeEta);
		    Float_t etaHigh = etaLow + binSizeEta;
		    Float_t phiLow = pMBc->GetZaxis()->GetXmin() + (phiBin * binSizePhi);
		    Float_t phiHigh = phiLow + binSizePhi;
		    Int_t etaLowBin = pMBc->GetYaxis()->FindBin(etaLow);
		    Int_t etaHighBin = pMBc->GetYaxis()->FindBin(etaHigh) - 1;
		    Int_t phiLowBin = pMBc->GetZaxis()->FindBin(phiLow);
		    Int_t phiHighBin = pMBc->GetZaxis()->FindBin(phiHigh) - 1;
		    Int_t pTLowBinMB = pMBc->GetXaxis()->FindBin(minPtMB);
		    Int_t pTHighBinMB = pMBc->GetXaxis()->FindBin(maxPtMB);
		    Int_t pTLowBinHT1 = pHT1c->GetXaxis()->FindBin(minPtHT1);
		    Int_t pTHighBinHT1 = pHT1c->GetXaxis()->FindBin(maxPtHT1);
		    Int_t pTLowBinHT2 = pHT2c->GetXaxis()->FindBin(minPtHT2);
		    Int_t pTHighBinHT2 = pHT2c->GetXaxis()->FindBin(maxPtHT2);
		    if (etaLowBin <= 0) etaLowBin = 1;
		    if (etaHighBin > pMBc->GetYaxis()->GetNbins()) etaHighBin = pMBc->GetYaxis()->GetNbins();
		    if (phiLowBin <= 0) phiLowBin = 1;
		    if (phiHighBin > pMBc->GetZaxis()->GetNbins()) phiHighBin = pMBc->GetZaxis()->GetNbins();
		    if (pTLowBinMB <= 0) pTLowBinMB = 1;
		    if (pTHighBinMB > pMBc->GetXaxis()->GetNbins()) pTHighBinMB = pMBc->GetXaxis()->GetNbins();
		    if (pTLowBinHT1 <= 0) pTLowBinHT1 = 1;
		    if (pTHighBinHT1 > pHT1c->GetXaxis()->GetNbins()) pTHighBinHT1 = pHT1c->GetXaxis()->GetNbins();
		    if (pTLowBinHT2 <= 0) pTLowBinHT2 = 1;
		    if (pTHighBinHT2 > pHT2c->GetXaxis()->GetNbins()) pTHighBinHT2 = pHT2c->GetXaxis()->GetNbins();
		    TH1 *hMBc = ProjectX3(pMBc, pTLowBinMB, pTHighBinMB, etaLowBin, etaHighBin, phiLowBin, phiHighBin);
		    TH1 *hHT1c = ProjectX3(pHT1c, pTLowBinHT1, pTHighBinHT1, etaLowBin, etaHighBin, phiLowBin, phiHighBin);
		    TH1 *hHT2c = ProjectX3(pHT2c, pTLowBinHT2, pTHighBinHT2, etaLowBin, etaHighBin, phiLowBin, phiHighBin);
		    TH1 *hMBp = ProjectX3(pMBp, pTLowBinMB, pTHighBinMB, etaLowBin, etaHighBin, phiLowBin, phiHighBin);
		    TH1 *hHT1p = ProjectX3(pHT1p, pTLowBinHT1, pTHighBinHT1, etaLowBin, etaHighBin, phiLowBin, phiHighBin);
		    TH1 *hHT2p = ProjectX3(pHT2p, pTLowBinHT2, pTHighBinHT2, etaLowBin, etaHighBin, phiLowBin, phiHighBin);
		    TString namehhStr, titlehhStr;
		    namehhStr = Form("%s_%i.%i", name, etaBin, phiBin);
		    titlehhStr = Form("#splitline{EMC points R_{CP}, %3.0f-%3.0f%%/%3.0f-%3.0f%%%s}{%1.2f < #eta < %1.2f, %1.2f < #phi < %1.2f};p_{T};R_{CP}"
			, PercentMostCentralLow, PercentMostCentralHigh, PercentMostPeripheralLow, PercentMostPeripheralHigh, title, etaLow, etaHigh, phiLow, phiHigh);
		    TH1F *hh = new TH1F(namehhStr.Data(), titlehhStr.Data(), 1, minPtMB, maxPtHT2);
		    if (hh) hh->Draw();
		    TLegend *legend = new TLegend(0.6, 0.8, 1.0, 1.0, "Legend");
		    Float_t sumW = 0;
		    Float_t sumWR = 0;
		    if (hMBc && hMBp) {
			if (hMBc->GetBinWidth(1) < binSizeMBPt) hMBc->Rebin(Int_t((binSizeMBPt / hMBc->GetBinWidth(1)) + 1));
			if (hMBp->GetBinWidth(1) < binSizeMBPt) hMBp->Rebin(Int_t((binSizeMBPt / hMBp->GetBinWidth(1)) + 1));
			hMBc->Divide(hMBp);
			hMBc->Scale(NBinaryDAuMostPeripheral / NBinaryDAuMostCentral);
			for (Int_t ibin = 1;ibin <= hMBc->GetXaxis()->GetNbins();ibin++) {
				if ((hMBc->GetBinError(ibin) > maxError) || (hMBc->GetBinContent(ibin) <= 0)) {
					hMBc->SetBinContent(ibin, -1);
					hMBc->SetBinError(ibin, 0);
				}
				if (hMBc->GetXaxis()->GetBinCenter(ibin) >= integralLow) {
					Float_t val = hMBc->GetBinContent(ibin);
					Float_t err = hMBc->GetBinError(ibin);
					if ((val > 0) && (err != 0)) {
						Float_t w = 1.0 / (err * err);
						sumW += w;
						sumWR += w * val;
					}
				}
			}
			TGraphErrors *gr = new TGraphErrors(hMBc);
			gr->SetMarkerColor(MBcolor);
			gr->SetLineColor(MBcolor);
			gr->SetFillColor(MBcolor);
			gr->SetMarkerStyle(8);
			gr->SetMarkerSize(1.2);
			gr->SetHistogram(hh);
			gr->Draw("SAME P");
			if (legend && gr) legend->AddEntry(gr, "MinBias", "P");
		    }
		    if (hHT1c && hHT1p) {
			if (hHT1c->GetBinWidth(1) < binSizeHT1Pt) hHT1c->Rebin(Int_t((binSizeHT1Pt / hHT1c->GetBinWidth(1)) + 1));
			if (hHT1p->GetBinWidth(1) < binSizeHT1Pt) hHT1p->Rebin(Int_t((binSizeHT1Pt / hHT1p->GetBinWidth(1)) + 1));
			hHT1c->Divide(hHT1p);
			hHT1c->Scale(NBinaryDAuMostPeripheral / NBinaryDAuMostCentral);
			for (Int_t ibin = 1;ibin <= hHT1c->GetXaxis()->GetNbins();ibin++) {
				if ((hHT1c->GetBinError(ibin) > maxError) || (hHT1c->GetBinContent(ibin) <= 0)) {
					hHT1c->SetBinContent(ibin, -1);
					hHT1c->SetBinError(ibin, 0);
				}
				if (hHT1c->GetXaxis()->GetBinCenter(ibin) >= integralLow) {
					Float_t val = hHT1c->GetBinContent(ibin);
					Float_t err = hHT1c->GetBinError(ibin);
					if ((val > 0) && (err != 0)) {
						Float_t w = 1.0 / (err * err);
						sumW += w;
						sumWR += w * val;
					}
				}
			}
			TGraphErrors *gr = new TGraphErrors(hHT1c);
			gr->SetMarkerColor(HT1color);
			gr->SetLineColor(HT1color);
			gr->SetFillColor(HT1color);
			gr->SetMarkerStyle(8);
			gr->SetMarkerSize(1.2);
			gr->SetHistogram(hh);
			gr->Draw("SAME P");
			if (legend && gr) legend->AddEntry(gr, "HighTower-1", "P");
		    }
		    if (hHT2c && hHT2p) {
			if (hHT2c->GetBinWidth(1) < binSizeHT2Pt) hHT2c->Rebin(Int_t((binSizeHT2Pt / hHT2c->GetBinWidth(1)) + 1));
			if (hHT2p->GetBinWidth(1) < binSizeHT2Pt) hHT2p->Rebin(Int_t((binSizeHT2Pt / hHT2p->GetBinWidth(1)) + 1));
			hHT2c->Divide(hHT2p);
			hHT2c->Scale(NBinaryDAuMostPeripheral / NBinaryDAuMostCentral);
			for (Int_t ibin = 1;ibin <= hHT2c->GetXaxis()->GetNbins();ibin++) {
				if ((hHT2c->GetBinError(ibin) > maxError) || (hHT2c->GetBinContent(ibin) <= 0)) {
					hHT2c->SetBinContent(ibin, -1);
					hHT2c->SetBinError(ibin, 0);
				}
				if (hHT2c->GetXaxis()->GetBinCenter(ibin) >= integralLow) {
					Float_t val = hHT2c->GetBinContent(ibin);
					Float_t err = hHT2c->GetBinError(ibin);
					if ((val > 0) && (err != 0)) {
						Float_t w = 1.0 / (err * err);
						sumW += w;
						sumWR += w * val;
					}
				}
			}
			TGraphErrors *gr = new TGraphErrors(hHT2c);
			gr->SetMarkerColor(HT2color);
			gr->SetLineColor(HT2color);
			gr->SetFillColor(HT2color);
			gr->SetMarkerStyle(8);
			gr->SetMarkerSize(1.2);
			gr->SetHistogram(hh);
			gr->Draw("SAME P");
			if (legend && gr) legend->AddEntry(gr, "HighTower-2", "P");
		    }
		    Float_t Rmean = -1;
		    Float_t RmeanErr = 0;
		    if (hEtaPhi && (sumW != 0)) {
			Rmean = sumWR / sumW;
			RmeanErr = sqrt(1.0 / sumW);
			hEtaPhi->SetBinContent(etaBin + 1, phiBin + 1, Rmean);
			hEtaPhi->SetBinError(etaBin + 1, phiBin + 1, RmeanErr);
		    }
		    if (Rmean > 0) {
			Float_t *systNbinX = new Float_t[1];
			Float_t *systNbinXerr = new Float_t[1];
			Float_t *systNbinY = new Float_t[1];
			Float_t *systNbinYerr = new Float_t[1];
			systNbinX[0] = (maxPtHT2 + integralLow) / 2.0;
			systNbinXerr[0] = (maxPtHT2 - integralLow) / 2.0;
			systNbinY[0] = Rmean;
			systNbinYerr[0] = RmeanErr;
			TGraphErrors *systNbinGr = new TGraphErrors(1, systNbinX, systNbinY, systNbinXerr, systNbinYerr);
			systNbinGr->SetFillColor(15);
			systNbinGr->SetLineColor(15);
			systNbinGr->SetMarkerColor(15);
			systNbinGr->SetHistogram(hh);
			systNbinGr->Draw("2");
			if (legend) legend->AddEntry(systNbinGr, "Mean", "P");
		    }
		    if (legend) legend->Draw();
		    showResultsOthersRcp(hMBc, legend, 0, 20);
		    if (hh) hh->SetMinimum(0);
		    if (hh) hh->SetMaximum(1.5);
		}
		TString nameEtaPhiCStr = TString(name) + "_points_pt_eta_phi_Rcp_mean";
		TString titleEtaPhiCStr = TString(title) + ", EMC points <R_{CP}> vs. #eta and #phi";
		TCanvas *cEtaPhi = new TCanvas(nameEtaPhiCStr.Data(), titleEtaPhiCStr.Data());
		if (cEtaPhi) cEtaPhi->cd();
		if (hEtaPhi) {
			if (hEtaPhi->GetXaxis()->GetNbins() == 1) {
				TH1F *hEtaPhi1 = new TH1F(hEtaPhi->GetName(), hEtaPhi->GetTitle(), hEtaPhi->GetYaxis()->GetNbins(), hEtaPhi->GetYaxis()->GetXmin(), hEtaPhi->GetYaxis()->GetXmax());
				hEtaPhi1->Sumw2();
				for (Int_t ibin = 1;ibin <= hEtaPhi->GetYaxis()->GetNbins();ibin++) {
					hEtaPhi1->SetBinContent(ibin, hEtaPhi->GetBinContent(1, ibin));
					hEtaPhi1->SetBinError(ibin, hEtaPhi->GetBinError(1, ibin));
				}
				hEtaPhi1->GetXaxis()->SetTitle(hEtaPhi->GetYaxis()->GetTitle());
				hEtaPhi1->GetYaxis()->SetTitle("<R_{CP}>");
				hEtaPhi1->Draw();
			} else if (hEtaPhi->GetYaxis()->GetNbins() == 1) {
				TH1F *hEtaPhi1 = new TH1F(hEtaPhi->GetName(), hEtaPhi->GetTitle(), hEtaPhi->GetXaxis()->GetNbins(), hEtaPhi->GetXaxis()->GetXmin(), hEtaPhi->GetXaxis()->GetXmax());
				hEtaPhi1->Sumw2();
				for (Int_t ibin = 1;ibin <= hEtaPhi->GetXaxis()->GetNbins();ibin++) {
					hEtaPhi1->SetBinContent(ibin, hEtaPhi->GetBinContent(ibin, 1));
					hEtaPhi1->SetBinError(ibin, hEtaPhi->GetBinError(ibin, 1));
				}
				hEtaPhi1->GetXaxis()->SetTitle(hEtaPhi->GetXaxis()->GetTitle());
				hEtaPhi1->GetYaxis()->SetTitle("<R_{CP}>");
				hEtaPhi1->Draw();
			} else {
				hEtaPhi->Draw("COLZ");
			}
		}
	}
*/
}

void showResultsDAuCentrality(const Char_t *dAuCentralityName, const Char_t *dAuCentralityTitle
        , Bool_t showRcp, Bool_t showRcpEta, Bool_t dAuCentralityShow, Bool_t dAuCentralityPrint, Bool_t showDAuCentralityTriggersSeparately
        , const TAnalysisSettings &settingsDAuNoCentral
        , const TAnalysisSettings &settingsDAuAllCentral
        , const TAnalysisSettings &settingsDAuMostCentral
        , const TAnalysisSettings &settingsDAuMidCentral
        , const TAnalysisSettings &settingsDAuMostPeripheral
        , const TAnalysisResults &resultsDAuNoCentral
        , const TAnalysisResults &resultsDAuMostCentral
        , const TAnalysisResults &resultsDAuMidCentral
        , const TAnalysisResults &resultsDAuMostPeripheral
        , bin_stat_list_type &dAuRcpMB
        , bin_stat_list_type &dAuRcpHT1
        , bin_stat_list_type &dAuRcpHT2
        , bin_stat_list_type &dAuRcp
        , bin_stat_list_type &dAuRcpMBEta
        , bin_stat_list_type &dAuRcpHT1Eta
        , bin_stat_list_type &dAuRcpHT2Eta
        , bin_stat_list_type &dAuRcpEta
) {
	clearPtBinPosFuncFromQCD();
	createPtBinPosFuncFromQCD(PQCDFUNC(settingsDAuNoCentral.usePPpQCD));
	Float_t denomNbinary = resultsDAuMostCentral.NBinaryCollisions / resultsDAuMostPeripheral.NBinaryCollisions;
	divideBins(resultsDAuMostCentral.spectrumMBDiv, resultsDAuMostPeripheral.spectrumMBDiv, dAuRcpMB, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT1Div, resultsDAuMostPeripheral.spectrumHT1Div, dAuRcpHT1, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT2Div, resultsDAuMostPeripheral.spectrumHT2Div, dAuRcpHT2, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumDiv, resultsDAuMostPeripheral.spectrumDiv, dAuRcp, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumMBDivEta, resultsDAuMostPeripheral.spectrumMBDivEta, dAuRcpMBEta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT1DivEta, resultsDAuMostPeripheral.spectrumHT1DivEta, dAuRcpHT1Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT2DivEta, resultsDAuMostPeripheral.spectrumHT2DivEta, dAuRcpHT2Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumDivEta, resultsDAuMostPeripheral.spectrumDivEta, dAuRcpEta, denomNbinary, 0.0, false, false);

	TLegend *legendCentral = 0;
	TString nameCentral = "spectrum_div_centralitybins"; nameCentral += dAuCentralityName;
	TString titleCentral = "Centrality dependence of p_{T} spectrum"; titleCentral += dAuCentralityTitle; titleCentral += ";p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}n}{d#etadp_{T}} #frac{1}{Smooth function}";
	TH1F *histCentral = 0;
	if (settingsDAuMidCentral.showTriggersSeparately) {
		histCentral = showLists(&resultsDAuMidCentral.spectrumMBDiv, &resultsDAuMidCentral.spectrumHT1Div, &resultsDAuMidCentral.spectrumHT2Div
		    , nameCentral.Data(), titleCentral.Data()
		    , histCentral, &legendCentral
		    , settingsDAuMidCentral.settingsMB, settingsDAuMidCentral.settingsHT1, settingsDAuMidCentral.settingsHT2
		    , dAuCentralityShow);
	} else {
		TAnalysisSettingsTrigger temp = settingsDAuMidCentral.settingsMB;
		temp.legendTitle = "STAR d+Au #pi^{0}";
		histCentral = showList(&resultsDAuMidCentral.spectrumDiv
		    , nameCentral.Data(), titleCentral.Data()
		    , histCentral, &legendCentral
		    , temp
		    , dAuCentralityShow);
	}
	if (histCentral) {
		histCentral->GetXaxis()->SetLimits(0, settingsDAuMostCentral.settingsHT2.highPt);
		if (settingsDAuMostCentral.showTriggersSeparately) {
		    histCentral = showLists(&resultsDAuMostCentral.spectrumMBDiv, &resultsDAuMostCentral.spectrumHT1Div, &resultsDAuMostCentral.spectrumHT2Div
			, nameCentral.Data(), titleCentral.Data()
    			, histCentral, &legendCentral
			, settingsDAuMostCentral.settingsMB, settingsDAuMostCentral.settingsHT1, settingsDAuMostCentral.settingsHT2
			, dAuCentralityShow);
		} else {
		    TAnalysisSettingsTrigger temp = settingsDAuMostCentral.settingsMB;
		    temp.legendTitle = "Most central";
		    histCentral = showList(&resultsDAuMostCentral.spectrumDiv
			, nameCentral.Data(), titleCentral.Data()
    			, histCentral, &legendCentral
		        , temp
		        , dAuCentralityShow);
		}
		if (settingsDAuMostPeripheral.showTriggersSeparately) {
		    histCentral = showLists(&resultsDAuMostPeripheral.spectrumMBDiv, &resultsDAuMostPeripheral.spectrumHT1Div, &resultsDAuMostPeripheral.spectrumHT2Div
			, nameCentral.Data(), titleCentral.Data()
    			, histCentral, &legendCentral
			, settingsDAuMostPeripheral.settingsMB, settingsDAuMostPeripheral.settingsHT1, settingsDAuMostPeripheral.settingsHT2
			, dAuCentralityShow);
		} else {
		    TAnalysisSettingsTrigger temp = settingsDAuMostPeripheral.settingsMB;
		    temp.legendTitle = "Most peripheral";
		    histCentral = showList(&resultsDAuMostPeripheral.spectrumDiv
			, nameCentral.Data(), titleCentral.Data()
    			, histCentral, &legendCentral
		        , temp
		        , dAuCentralityShow);
		}
		histCentral->SetMinimum(0);
		histCentral->SetMaximum(6);
	}
	showResultsOthers(histCentral, legendCentral, false, settingsDAuNoCentral.usePPpQCD);

	TString titleRcpCentrality = Form(", %3.0f-%3.0f%%/%3.0f-%3.0f%%", resultsDAuMostCentral.PercentCentralityLow, resultsDAuMostCentral.PercentCentralityHigh, resultsDAuMostPeripheral.PercentCentralityLow, resultsDAuMostPeripheral.PercentCentralityHigh);
	TLegend *legendRcp = 0;
	TString nameRcp = "RdAuCP"; nameRcp += dAuCentralityName;
	TString titleRcp = "R^{dAu}_{CP}"; titleRcp += titleRcpCentrality; titleRcp += dAuCentralityTitle; titleRcp += ";p_{T}, GeV/c;R^{dAu}_{CP}";
	TH1F *histRcp = 0;
	if (showRcp) {
	    if (showDAuCentralityTriggersSeparately) {
		histRcp = showLists(&dAuRcpMB, &dAuRcpHT1, &dAuRcpHT2
		    , nameRcp.Data(), titleRcp.Data()
		    , histRcp, &legendRcp
		    , settingsDAuNoCentral.settingsMB, settingsDAuNoCentral.settingsHT1, settingsDAuNoCentral.settingsHT2
		    , dAuCentralityShow);
	    } else {
		TAnalysisSettingsTrigger temp = settingsDAuNoCentral.settingsMB;
		temp.legendTitle = "STAR #pi^{0} R^{dAu}_{CP}";
		histRcp = showList(&dAuRcp
		    , nameRcp.Data(), titleRcp.Data()
		    , histRcp, &legendRcp
		    , temp
		    , dAuCentralityShow);
	    }
	}
	if (showRcpEta) {
	    if (showDAuCentralityTriggersSeparately) {
		histRcp = showLists(&dAuRcpMBEta, &dAuRcpHT1Eta, &dAuRcpHT2Eta
		    , nameRcp.Data(), titleRcp.Data()
		    , histRcp, &legendRcp
		    , settingsDAuNoCentral.settingsSimEtaMB, settingsDAuNoCentral.settingsSimEtaHT1, settingsDAuNoCentral.settingsSimEtaHT2
		    , dAuCentralityShow);
	    } else {
		TAnalysisSettingsTrigger temp = settingsDAuNoCentral.settingsMB;
		temp.legendTitle = "STAR #pi^{0} R^{dAu}_{CP} #eta";
		histRcp = showList(&dAuRcpEta
		    , nameRcp.Data(), titleRcp.Data()
		    , histRcp, &legendRcp
		    , temp
		    , dAuCentralityShow);
	    }
	}
	if (histRcp) {
	    histRcp->GetXaxis()->SetLimits(0, settingsDAuMostCentral.settingsHT2.highPt);
	}
	showResultsOthersDAuRcp(histRcp, legendRcp, resultsDAuMostCentral, resultsDAuMostPeripheral);

	if (dAuCentralityPrint) {
		Float_t mostCentralPercent = (resultsDAuNoCentral.resultsMB.evNum > 0) ? (resultsDAuMostCentral.resultsMB.evNum / resultsDAuNoCentral.resultsMB.evNum) : 0.0;
		Float_t midCentralPercent = (resultsDAuNoCentral.resultsMB.evNum > 0) ? (resultsDAuMidCentral.resultsMB.evNum / resultsDAuNoCentral.resultsMB.evNum) : 0.0;
		Float_t mostPeripheralPercent = (resultsDAuNoCentral.resultsMB.evNum > 0) ? (resultsDAuMostPeripheral.resultsMB.evNum / resultsDAuNoCentral.resultsMB.evNum) : 0.0;
		Float_t mostCentralToMostPeripheral = (resultsDAuMostPeripheral.resultsMB.evNum > 0) ? (resultsDAuMostCentral.resultsMB.evNum / resultsDAuMostPeripheral.resultsMB.evNum) : 0.0;
		cout << "==============================================" << endl;
		cout << "Current centrality bins" << dAuCentralityTitle << ":" << endl;
		cout << "Most central: " << (mostCentralPercent * 100.0) << "%" << endl;
		cout << "Mid central: " << (midCentralPercent * 100.0) << "%" << endl;
		cout << "Most peripheral: " << (mostPeripheralPercent * 100.0) << "%" << endl;
		cout << "Sum: " << ((mostCentralPercent + midCentralPercent + mostPeripheralPercent) * 100) << "%" << endl;
		cout << "Most central / most peripheral: " << (mostCentralToMostPeripheral) << endl;
	}
}

void showResultsRDA(const Char_t *RDAName, const Char_t *RDATitle
        , Bool_t showRDA, Bool_t showRDAEta, Bool_t RDAShow, Bool_t showRDATriggersSeparately
        , const TAnalysisSettings &settingsDAuNoCentral
        , const TAnalysisSettings &settingsPP
        , const TAnalysisResults &resultsDAuNoCentral
        , const TAnalysisResults &resultsPP
        , bin_stat_list_type &RDAMB
        , bin_stat_list_type &RDAHT1
        , bin_stat_list_type &RDAHT2
        , bin_stat_list_type &RDA
        , bin_stat_list_type &RDAMBEta
        , bin_stat_list_type &RDAHT1Eta
        , bin_stat_list_type &RDAHT2Eta
        , bin_stat_list_type &RDAEta
) {
	Float_t denomNbinary = resultsDAuNoCentral.NBinaryCollisions / resultsPP.NBinaryCollisions;
	denomNbinary *= (settingsDAuNoCentral.crossectionMB / settingsDAuNoCentral.crossectionMBTotalFraction) / (settingsPP.crossectionMB / settingsPP.crossectionMBTotalFraction);
	divideBins(resultsDAuNoCentral.crossectionMB, resultsPP.crossectionMB, RDAMB, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossectionHT1, resultsPP.crossectionHT1, RDAHT1, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossectionHT2, resultsPP.crossectionHT2, RDAHT2, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossection, resultsPP.crossection, RDA, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossectionMBEta, resultsPP.crossectionMBEta, RDAMBEta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossectionHT1Eta, resultsPP.crossectionHT1Eta, RDAHT1Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossectionHT2Eta, resultsPP.crossectionHT2Eta, RDAHT2Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.crossectionEta, resultsPP.crossectionEta, RDAEta, denomNbinary, 0.0, false, false);

	TLegend *legendRDA = 0;
	TString nameRDA = "RDA"; nameRDA += RDAName;
	TString titleRDA = "R_{dA}"; titleRDA += RDATitle; titleRDA += ";p_{T}, GeV/c;R_{dA}";
	TH1F *histRDA = 0;
	if (showRDA) {
	    if (showRDATriggersSeparately) {
		histRDA = showLists(&RDAMB, &RDAHT1, &RDAHT2
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , settingsDAuNoCentral.settingsMB, settingsDAuNoCentral.settingsHT1, settingsDAuNoCentral.settingsHT2
		    , RDAShow);
	    } else {
		TAnalysisSettingsTrigger temp = settingsDAuNoCentral.settingsMB;
		temp.legendTitle = "STAR #pi^{0} R_{pA}";
		histRDA = showList(&RDA
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , temp
		    , RDAShow);
	    }
	}
	if (showRDAEta) {
	    if (showRDATriggersSeparately) {
		histRDA = showLists(&RDAMBEta, &RDAHT1Eta, &RDAHT2Eta
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , settingsDAuNoCentral.settingsMB, settingsDAuNoCentral.settingsHT1, settingsDAuNoCentral.settingsHT2
		    , RDAShow);
	    } else {
		TAnalysisSettingsTrigger temp = settingsDAuNoCentral.settingsMB;
		temp.legendTitle = "STAR #eta R_{pA}";
		histRDA = showList(&RDAEta
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , temp
		    , RDAShow);
	    }
	}
	if (histRDA) {
	    histRDA->GetXaxis()->SetLimits(0, settingsDAuNoCentral.settingsHT2.highPt);
	}
	//showResultsOthersDAuRcp(histRDA, legendRDA, resultsDAuNoCentral, resultsPP);
}

Float_t dAuMinBiasNBinaryCollisions = 7.5;
Float_t dAuMinBiasNBinaryCollisionsErr = 0.38;
void getCentralityBin(Int_t RefMultLow, Int_t RefMultHigh, Float_t &percentLow, Float_t &percentHigh, Float_t &Nbin, Float_t &NbinErr) {
	if ((RefMultLow <= 0) && (RefMultHigh > 100)) {
		percentLow = 0; percentHigh = 100; Nbin = dAuMinBiasNBinaryCollisions; NbinErr = dAuMinBiasNBinaryCollisionsErr;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 100; Nbin = 4.0; NbinErr = 0.3;
	} else if ((RefMultLow == 10) && (RefMultHigh == 17)) {
		percentLow = 20; percentHigh = 40; Nbin = 10.0; NbinErr = 0.6;
	} else if ((RefMultLow == 17) && (RefMultHigh > 100)) {
		percentLow = 0; percentHigh = 20; Nbin = 15.0; NbinErr = 1.1;
	} else if ((RefMultLow == 3) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 80; Nbin = 4.9; NbinErr = 0.5;
	} else if ((RefMultLow == 4) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 75; Nbin = 5.34; NbinErr = 0.65;
	} else if ((RefMultLow == 5) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 70; Nbin = 5.77; NbinErr = 0.8;
	} else if ((RefMultLow == 6) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 65; Nbin = 6.2; NbinErr = 0.9;
	} else if ((RefMultLow == 7) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 60; Nbin = 6.64; NbinErr = 0.98;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 3)) {
		percentLow = 80; percentHigh = 100; Nbin = 1.0; NbinErr = 0.01;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 4)) {
		percentLow = 75; percentHigh = 100; Nbin = 1.0; NbinErr = 0.01;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 5)) {
		percentLow = 70; percentHigh = 100; Nbin = 1.0; NbinErr = 0.01;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 6)) {
		percentLow = 65; percentHigh = 100; Nbin = 1.0; NbinErr = 0.01;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 7)) {
		percentLow = 60; percentHigh = 100; Nbin = 1.0; NbinErr = 0.01;
	} else {
		percentLow = 0; percentHigh = 100; Nbin = dAuMinBiasNBinaryCollisions; NbinErr = dAuMinBiasNBinaryCollisionsErr;
	}
//cout << "getCentralityBin: RefMultLow = " << RefMultLow << ", RefMultHigh = " << RefMultHigh;
//cout << ", percentLow = " << percentLow << ", percentHigh = " << percentHigh;
//cout << ", Nbin = " << Nbin << ", NbinErr = " << NbinErr << endl;
}

void calculateBunchCrossingId7bitOffset(const TEventDataProcessor *eventDataProcessorAllValid
	, bool show, bool print, const Char_t *outputFilename
	, const Char_t *outputHighbgFilename, const Char_t *highBgArrayname
	, Float_t bgFractionCut, Bool_t checkAbortGaps, Bool_t checkEmptyBunches)
{
	const TH3F *bcdistr = eventDataProcessorAllValid ? ((const TH3F*)eventDataProcessorAllValid->getDayRunBunchCrossing()) : 0;
	if (!bcdistr) {
		if (print) cout << "No data!" << endl;
		return;
	}
	Int_t year = getRunYear(eventDataProcessorAllValid);

	Int_t digitLow = 0;
	Int_t digitHigh = 120;
	Int_t digitLowBin = bcdistr->GetZaxis()->FindBin(digitLow);
	Int_t digitHighBin = bcdistr->GetZaxis()->FindBin(digitHigh);
	Int_t digitDiff = digitHigh - digitLow;

	TH1F *bunchCrossingOffset = new TH1F("bunchCrossingOffset", "BunchCrossing7bit + Offset;(Bunch crossing + offset) % 120", bcdistr->GetZaxis()->GetNbins(), bcdistr->GetZaxis()->GetBinLowEdge(1), bcdistr->GetZaxis()->GetBinUpEdge(bcdistr->GetZaxis()->GetNbins()));
	TH1F *bunchCrossingOffsetTrue = new TH1F("bunchCrossingOffsetTrue", "BunchCrossing7bit + Offset;(Bunch crossing + offset) % 120", bcdistr->GetZaxis()->GetNbins(), bcdistr->GetZaxis()->GetBinLowEdge(1), bcdistr->GetZaxis()->GetBinUpEdge(bcdistr->GetZaxis()->GetNbins()));
	TH1F *bunchCrossingNoOffset = new TH1F("bunchCrossingNoOffset", "Bunch crossing (7bit);Bunch crossing (7bit)", bcdistr->GetZaxis()->GetNbins(), bcdistr->GetZaxis()->GetBinLowEdge(1), bcdistr->GetZaxis()->GetBinUpEdge(bcdistr->GetZaxis()->GetNbins()));
	TH1F *bunchCrossingChi2 = new TH1F("bunchCrossingChi2", "#chi^{2};Offset - best offset;#chi^{2}", digitDiff, digitLow, digitHigh);
	TH1F *bunchCrossingChi2tmp = new TH1F("bunchCrossingChi2tmp", "bunchCrossing - chi2 (tmp)", digitDiff, digitLow, digitHigh);
	TH1F *offsets = new TH1F("offsets", "Offsets;Offset;Number of runs", digitDiff, 0, digitDiff);
	TH1F *background = new TH1F("background", "Background fraction;Background / total events per run;Number of runs", 100, 0, 2*18.0/120.0);
	Int_t ndays = bcdistr->GetXaxis()->GetNbins();
	Int_t nruns = bcdistr->GetYaxis()->GetNbins();
	Int_t *highbgRuns = new Int_t[1024];
	Float_t *highbgRatios = new Float_t[1024];
	Int_t numHighbgRuns = 0;
	Int_t numTotalRuns = 0;
	Float_t integralHighbgRuns = 0;
	Float_t integralTotalRuns = 0;
	ofstream ofstr(outputFilename);
	for (Int_t iday = 1;iday <= ndays;iday++) {
		for (Int_t irun = 1;irun <= nruns;irun++) {
			Int_t runDay = (Int_t)bcdistr->GetYaxis()->GetBinLowEdge(irun);
			Int_t day = (Int_t)bcdistr->GetXaxis()->GetBinLowEdge(iday);
			Int_t runId = (year * 1000000) + (day * 1000) + (runDay * 1);
//if (runId != 4060094) continue;
//if ((runId == 4034033) || (runId == 4038034) || (runId == 4041001) || (runId == 4041002) || (runId == 4055012) || (runId == 4055013)) continue;
//if (runId < 4079057) continue;
//if (runId > 4060000) continue;
//if (day != 39) continue;
			Float_t integral = 0;
			Int_t integralNum = 0;
			for (Int_t ibunch = digitLowBin;ibunch < digitHighBin;ibunch++) {
				Int_t bunch = (Int_t)bcdistr->GetZaxis()->GetBinLowEdge(ibunch);
				Float_t realNum = bcdistr->GetBinContent(bcdistr->GetBin(iday, irun, ibunch));
				if (bunchCrossingNoOffset) bunchCrossingNoOffset->Fill(bunch, realNum);
				integral += realNum;
				integralNum++;
			}
			Float_t average = (integralNum != 0) ? (integral / integralNum) : 0;
			integral = 0;
			integralNum = 0;
			Float_t integralErr = 0;
			for (Int_t ibunch = digitLowBin;ibunch < digitHighBin;ibunch++) {
				Float_t realNum = bcdistr->GetBinContent(bcdistr->GetBin(iday, irun, ibunch));
				Float_t realNumErr = bcdistr->GetBinError(bcdistr->GetBin(iday, irun, ibunch));
				if (realNum >= average) {
					integral += realNum;
					integralErr += realNumErr * realNumErr;
					integralNum++;
				}
			}
			integralErr = sqrt(integralErr);
			average = (integralNum != 0) ? (integral / integralNum) : 0;
			Float_t averageErr = (integralNum != 0) ? (integralErr / integralNum) : 1.0;
			Float_t chi2min = -1;
			Int_t bestOffset = -1;
			if (integral != 0) {
				numTotalRuns++;
				integralTotalRuns += integral;
				Int_t nbunches = 0;
				for (Int_t ibunch = digitLowBin;ibunch < digitHighBin;ibunch++) {
					Float_t realNum = bcdistr->GetBinContent(bcdistr->GetBin(iday, irun, ibunch));
					if (realNum >= (average / 4.0)) {
						nbunches++;
					}
				}
//cout << runId << " " << nbunches << endl;
				if (bunchCrossingChi2tmp) for (Int_t i = 1;i <= bunchCrossingChi2tmp->GetXaxis()->GetNbins();i++) bunchCrossingChi2tmp->SetBinContent(i, 0);
				for (Int_t offset = 0;offset < digitDiff;offset++) {
					Float_t chi2 = 0;
					for (Int_t ibunch = digitLowBin;ibunch < digitHighBin;ibunch++) {
						Int_t bunch = (Int_t)bcdistr->GetZaxis()->GetBinLowEdge(ibunch) + offset;
						while (bunch < digitLow) bunch += digitDiff;
						while (bunch >= digitHigh) bunch -= digitDiff;
						Float_t realNum = bcdistr->GetBinContent(bcdistr->GetBin(iday, irun, ibunch));
						Float_t realNumErr = bcdistr->GetBinError(bcdistr->GetBin(iday, irun, ibunch));
						if (realNumErr == 0) realNumErr = averageErr;
						Bool_t isGood = isGoodBunchCrossingId7bitPlusOffset(runId, year, day, runDay, bunch - digitLow, true, true);
						Float_t trueNum = isGood ? average : 0;
						Float_t trueNumErr = averageErr;
						Float_t reducedErr = ((realNumErr != 0) && (trueNumErr != 0)) ? ((1 / (realNumErr * realNumErr)) + (1 / (trueNumErr * trueNumErr))) : 0;
						reducedErr = (reducedErr > 0) ? (1.0 / sqrt(reducedErr)) : 0;
						Float_t diff = (realNum - trueNum);
						Float_t chi2tmp = (diff * diff) / (2 * reducedErr * reducedErr);
						if (realNum > (average / 2.0)) {
							chi2 += chi2tmp;
						}
					}
					chi2 = chi2 / 119.0;
					//if (print) cout << "offset " << offset << ": chi2 = " << chi2 << endl;
					if ((chi2 < chi2min) || (chi2min == -1)) {
						chi2min = chi2;
						bestOffset = offset;
						//cout << "this is the best offset " << bestOffset << endl;
					}
//cout << offset << " " << chi2 << endl;
					if (bunchCrossingChi2tmp) bunchCrossingChi2tmp->Fill(offset, chi2);
				}
			}
//bestOffset = 6; // for pp2005
			if ((chi2min != -1) && (bestOffset != -1) && (integral != 0)) {
				//if (print) cout << "Run " << runId << ": offset " << bestOffset << endl;
				ofstr << runId << " " << bestOffset << endl;
				if (offsets) offsets->Fill(bestOffset);
				Float_t integralGap = 0;
				for (Int_t ibunch = digitLowBin;ibunch < digitHighBin;ibunch++) {
					Int_t bunch = (Int_t)bcdistr->GetZaxis()->GetBinLowEdge(ibunch);
					Int_t bunchOffset = bunch + bestOffset;
					while (bunchOffset < digitLow) bunchOffset += digitDiff;
					while (bunchOffset >= digitHigh) bunchOffset -= digitDiff;
					Float_t realNumOffset = bcdistr->GetBinContent(bcdistr->GetBin(iday, irun, ibunch));
					Bool_t isGood = isGoodBunchCrossingId7bitPlusOffset(runId, year, day, runDay, bunchOffset - digitLow, true, true);
					Float_t trueNumOffset = isGood ? average : 0;
					if (bunchCrossingOffset) bunchCrossingOffset->Fill(bunchOffset, realNumOffset);
					if (bunchCrossingOffsetTrue) bunchCrossingOffsetTrue->Fill(bunchOffset, trueNumOffset);
					Bool_t isGoodBg = isGoodBunchCrossingId7bitPlusOffset(runId, year, day, runDay, bunchOffset - digitLow, checkAbortGaps, checkEmptyBunches);
					if (!isGoodBg) integralGap += realNumOffset;
				}
				Float_t gapFraction = integralGap / integral;
				if (background) background->Fill(gapFraction);
				if (gapFraction > bgFractionCut) {
					if (numHighbgRuns < 1024) {
						if (highbgRuns) highbgRuns[numHighbgRuns] = runId;
						if (highbgRatios) highbgRatios[numHighbgRuns] = gapFraction;
						numHighbgRuns++;
						integralHighbgRuns += integral;
					}
				}
				for (Int_t i = digitLow;i < digitHigh;i++) {
					Int_t bunchOffset = i + bestOffset;
					while (bunchOffset < digitLow) bunchOffset += digitDiff;
					while (bunchOffset >= digitHigh) bunchOffset -= digitDiff;
					Float_t c = bunchCrossingChi2tmp->GetBinContent(bunchCrossingChi2tmp->GetXaxis()->FindBin(bunchOffset));
					if (bunchCrossingChi2 && bunchCrossingChi2tmp) bunchCrossingChi2->Fill(i, c);
				}
			}
		}
	}
	ofstr.close();
	for (Int_t i = 0;i < numHighbgRuns - 1;i++) {
		for (Int_t j = i + 1;j < numHighbgRuns;j++) {
			if (highbgRatios && highbgRuns && (highbgRatios[i] < highbgRatios[j])) {
				Int_t tmpRun = highbgRuns[i];
				highbgRuns[i] = highbgRuns[j];
				highbgRuns[j] = tmpRun;
				Float_t tmpRatio = highbgRatios[i];
				highbgRatios[i] = highbgRatios[j];
				highbgRatios[j] = tmpRatio;
			}
		}
	}
	ofstream ofBgstr(outputHighbgFilename);
	ofBgstr << "// Runs with high background in the empty time buckets (cut at " << bgFractionCut << ")" << endl;
	ofBgstr << "const Int_t " << highBgArrayname << "[] = {" << endl;
	for (Int_t i = 0;i < numHighbgRuns;i++) {
		if (highbgRuns && highbgRatios) {
			if (print) cout << (i + 1) << ". " << highbgRuns[i] << " " << highbgRatios[i] << endl;
			ofBgstr << highbgRuns[i] << ", // " << highbgRatios[i] << endl;
		}
	}
	ofBgstr << "-1 // The end" << endl;
	ofBgstr << "};" << endl;
	Float_t rej = (integralTotalRuns != 0) ? (integralHighbgRuns / integralTotalRuns) : 0;
	ofBgstr << "// Rejected " << numHighbgRuns << " of " << numTotalRuns << " runs (" << (Int_t(rej * 1000)/10.0) << "% of events)" << endl;
	ofBgstr.close();
	if (bunchCrossingChi2tmp) delete bunchCrossingChi2tmp;
	if (highbgRuns) delete [] highbgRuns;
	if (highbgRatios) delete [] highbgRatios;
	if (show) {
		TCanvas *c = new TCanvas("bunchCrossingId7bitOffset", "bunchCrossingId7bit Offset");
		if (c) {
			c->Divide(2, 2, 0.001, 0.001);
			c->cd(1);
			if (bunchCrossingNoOffset) bunchCrossingNoOffset->Draw();
			c->cd(2);
			if (bunchCrossingChi2) bunchCrossingChi2->Draw();
			c->cd(4);
			if (background) {
				background->Draw();
				TLine *l = new TLine(bgFractionCut, 0, bgFractionCut, background->GetMaximum());
				if (l) {
					l->SetLineColor(kRed);
					l->Draw();
				}
				TString s = Form("#splitline{Rejected:}{#splitline{%.0f of %.0f runs}{%3.1f%% of events}}", Float_t(numHighbgRuns), Float_t(numTotalRuns), Float_t(Int_t(rej * 1000)/10.0));
				TLatex *t = new TLatex(bgFractionCut + 0.1, background->GetMaximum() * 3.0 / 4.0, s.Data());
				if (t) t->Draw();
			}
			c->cd(3);
			if (bunchCrossingOffset) {
				bunchCrossingOffset->SetLineColor(kBlack);
				bunchCrossingOffset->SetFillColor(16);
				bunchCrossingOffset->Draw();
			}
			if (bunchCrossingOffsetTrue) {
				bunchCrossingOffsetTrue->SetLineColor(kBlue);
				bunchCrossingOffsetTrue->SetLineWidth(2);
				bunchCrossingOffsetTrue->Draw(bunchCrossingOffset ? "SAME" : "");
			}
			//c->cd(5);
			//if (offsets) offsets->Draw();
			c->cd(0);
		}
	}
}

void getRunTimes(const TEventDataProcessor *events, const Char_t *runTimesFilename, bool print, const Char_t *dbStr) {
    if (!events || !runTimesFilename) return;
    const TH2F *dayRuns = events->getDayRun();
    if (!dayRuns) return;
    if (dbStr && dbStr[0]) {
	DBConnectString = dbStr;
	if (print) cout << "Using database " << DBConnectString << endl;
    }
    if (print) cout << "Saving timestamps from \"" << events->GetTitle() << "\" into " << runTimesFilename << endl;
    Int_t year = getRunYear(events);
    ofstream runTimesStr(runTimesFilename);
    for (Int_t ix = 1;ix <= dayRuns->GetXaxis()->GetNbins();ix++) {
	for (Int_t iy = 1;iy <= dayRuns->GetYaxis()->GetNbins();iy++) {
	    Float_t evNum = (Float_t)dayRuns->GetBinContent(dayRuns->GetBin(ix, iy));
	    if (evNum > 0) {
		Int_t day = (Int_t)dayRuns->GetXaxis()->GetBinCenter(ix);
		Int_t run = (Int_t)dayRuns->GetYaxis()->GetBinCenter(iy);
		Int_t runId = (year * 1000000) + (day * 1000) + run;
		Int_t runDate = -1;
		Int_t runTime = -1;
		getRunTimesFromDB(runId, runDate, runTime, print);
		if (runTimesStr.good() && (runDate != -1) && (runTime != -1)) {
		    runTimesStr << runId << "\t" << runDate << "\t" << runTime << "\t" << evNum << endl;
		    runTimesStr.flush();
		}
	    }
	}
    }
    runTimesStr.close();
}

void getPrescales(const TEventDataProcessor *eventsMB, const TEventDataProcessor *eventsHT1, const TEventDataProcessor *eventsHT2, const Char_t *prescalesFilename, bool print, const Char_t *dbStr, const Int_t *ourMinBiasTriggers, const Int_t *ourHt1Triggers, const Int_t *ourHt2Triggers) {
    const TEventDataProcessor *events = eventsMB ? eventsMB : (eventsHT1 ? eventsHT1 : eventsHT2);
    if (!events || !prescalesFilename) return;
    const TH2F *dayRunsMB = eventsMB ? eventsMB->getDayRun() : 0;
    const TH2F *dayRunsHT1 = eventsHT1 ? eventsHT1->getDayRun() : 0;
    const TH2F *dayRunsHT2 = eventsHT2 ? eventsHT2->getDayRun() : 0;
    const TH2F *dayRuns = dayRunsMB ? dayRunsMB : (dayRunsHT1 ? dayRunsHT1 : dayRunsHT2);
    if (!dayRuns) return;
    if (dbStr && dbStr[0]) {
	DBConnectString = dbStr;
	if (print) cout << "Using database " << DBConnectString << endl;
    }
    if (print) cout << "Saving prescales into " << prescalesFilename << ";" << endl;
    Int_t year = getRunYear(events);
    ofstream prescalesStr(prescalesFilename);
    for (Int_t ix = 1;ix <= dayRuns->GetXaxis()->GetNbins();ix++) {
	for (Int_t iy = 1;iy <= dayRuns->GetYaxis()->GetNbins();iy++) {
	    Float_t evNumMB = dayRunsMB ? (Float_t)dayRunsMB->GetBinContent(dayRuns->GetBin(ix, iy)) : 0;
	    Float_t evNumHT1 = dayRunsHT1 ? (Float_t)dayRunsHT1->GetBinContent(dayRuns->GetBin(ix, iy)) : 0;
	    Float_t evNumHT2 = dayRunsHT2 ? (Float_t)dayRunsHT2->GetBinContent(dayRuns->GetBin(ix, iy)) : 0;
	    if (((evNumMB > 0) || (evNumHT1 > 0) || (evNumHT2 > 0)) && (true)) {
		Int_t day = (Int_t)dayRuns->GetXaxis()->GetBinCenter(ix);
		Int_t run = (Int_t)dayRuns->GetYaxis()->GetBinCenter(iy);
		Int_t runId = (year * 1000000) + (day * 1000) + run;
		//if (runId < 6133071) continue;
		Int_t psMB = 0;
		Int_t psHT1 = 0;
		Int_t psHT2 = 0;
		getPrescalesFromDB(runId, psMB, psHT1, psHT2, print, ourMinBiasTriggers, ourHt1Triggers, ourHt2Triggers);
		if (prescalesStr.good()) {
		    prescalesStr << runId << "\t" << psMB << "\t" << psHT1 << "\t" << psHT2 << endl;
		    prescalesStr.flush();
		}
	    }
	}
    }
    prescalesStr.close();
}

TH1F *showList(const list<TBinStatistics> *points, const Char_t *name, const Char_t *title
	, TH1F *oldHist, TLegend* *legend
	, const TAnalysisSettingsTrigger &settings
	, Bool_t show
) {
	const Int_t bufferSize = (points ? points->size() : 0) + 1;
	Float_t *bufferX = new Float_t[bufferSize];
	Float_t *bufferY = new Float_t[bufferSize];
	Float_t *bufferXErr = new Float_t[bufferSize];
	Float_t *bufferYErr = new Float_t[bufferSize];
	Float_t minX = +1000;
	Float_t maxX = -1000;
	Float_t minY = +1000;
	Float_t maxY = -1000;
	Int_t numPoints = 0;
	if (points) for (list<TBinStatistics>::const_iterator iter = points->begin();iter != points->end();++iter) {
		const TBinStatistics &bin = *iter;
		Float_t x = bin.getParameters().getCenter();
		Float_t xErr = settings.useXerr ? (bin.getParameters().getWidth() / 2.0) : 0;
		Float_t v = bin.getValue();
		Float_t e = bin.getError();
		bufferX[numPoints] = x;
		bufferXErr[numPoints] = xErr;
		bufferY[numPoints] = v;
		bufferYErr[numPoints] = e;
		minX = min(minX, bin.getParameters().min - xErr);
		maxX = max(maxX, bin.getParameters().max + xErr);
		minY = min(minY, v - e);
		maxY = max(maxY, v + e);
		if ((numPoints + 1) < bufferSize) numPoints++;
	}
	TH1F *histRatio = oldHist;
	TLegend *legendL = legend ? *legend : 0;
	if (show && settings.showBins) {
	    if ((numPoints + 0) > 0) {
		if (histRatio) {
			minX = min(minX, (Float_t)histRatio->GetXaxis()->GetBinLowEdge(1));
			maxX = max(maxX, (Float_t)histRatio->GetXaxis()->GetBinUpEdge(histRatio->GetXaxis()->GetNbins()));
			minY = min(minY, (Float_t)histRatio->GetYaxis()->GetBinLowEdge(1));
			maxY = max(maxY, (Float_t)histRatio->GetYaxis()->GetBinUpEdge(histRatio->GetYaxis()->GetNbins()));
    		} else {
			TCanvas *clist = new TCanvas(name, title);
			if (clist) {
				clist->cd();
				clist->SetFillColor(kWhite);
				TString histName = TString(name) + "_hist";
				histRatio = new TH1F(histName.Data(), title, 10, 0, 1);
				if (histRatio) {
					histRatio->SetFillColor(kWhite);
					histRatio->SetMinimum(minY);
					histRatio->SetMaximum(maxY);
					histRatio->SetStats(0);
					histRatio->Draw();
				}
			}
		}
		if (histRatio) {
		    histRatio->GetXaxis()->SetLimits(minX, maxX);
		    histRatio->GetYaxis()->SetLimits(minY, maxY);
		}
		if (!legendL) {
		    legendL = new TLegend(0.6, 0.6, 0.89, 0.89);
		    if (legendL) {
			legendL->SetFillColor(kWhite);
		    }
		}
	    }
	    if (numPoints > 0) {
		TGraphErrors *ratio = new TGraphErrors(numPoints, bufferX, bufferY, bufferXErr, bufferYErr);
		if (ratio) {
		    ratio->SetTitle(title);
		    ratio->SetMarkerStyle(settings.markerStyle);
		    ratio->SetMarkerSize(settings.markerSize);
		    ratio->SetLineStyle(settings.lineStyle);
		    ratio->SetLineWidth(settings.lineSize);
		    ratio->SetLineColor(settings.color);
		    ratio->SetMarkerColor(settings.color);
		    ratio->SetFillColor(settings.fillColor);
		    ratio->SetHistogram(histRatio);
		    ratio->Draw(settings.drawOption);
		    if (legendL && (settings.legendTitle != "")) legendL->AddEntry(ratio, settings.legendTitle.Data(), settings.legendOption);
		}
	    }
	    if (legendL) legendL->Draw("SAME");
	}
	if (legend) *legend = legendL;
	delete [] bufferX;
	delete [] bufferXErr;
	delete [] bufferY;
	delete [] bufferYErr;
	return histRatio;
}

TH1F *showLists(const list<TBinStatistics> *pointsMB, const list<TBinStatistics> *pointsHT1, const list<TBinStatistics> *pointsHT2
        , const Char_t *name, const Char_t *title
        , TH1F *oldHist, TLegend* *legend
	, const TAnalysisSettingsTrigger &settingsMB, const TAnalysisSettingsTrigger &settingsHT1, const TAnalysisSettingsTrigger &settingsHT2
        , Bool_t show
) {
    TH1F *hist = oldHist;
    hist = showList(pointsMB, name, title, hist, legend, settingsMB, show);
    hist = showList(pointsHT1, name, title, hist, legend, settingsHT1, show);
    hist = showList(pointsHT2, name, title, hist, legend, settingsHT2, show);
    return hist;
}

void calculateVertexFindingEff(const Char_t *title, Bool_t calculate, const TEventDataProcessor *event, Float_t &eff) {
    if (calculate) {
	if (title) cout << title << ", vertex finding eff: ";
	if (!event) {
	    eff = 1.0;
    	    cout << "No event data processor!";
	} else {
	    const TCuts &cuts = event->getCuts();
	    if (cuts.EVENT_ALL_CUTS | EVENT_TPC_VERTEX_CUT) {
		Int_t tot = event->numTotal;
	        Int_t passed = 0;
		if (cuts.EVENT_passedCuts_separate.find(EVENT_TPC_VERTEX_CUT) != cuts.EVENT_passedCuts_separate.end()) {
		    passed = (*(cuts.EVENT_passedCuts_separate.find(EVENT_TPC_VERTEX_CUT))).second;
		}
	        eff = (tot != 0) ? (Float_t(passed) / Float_t(tot)) : 1.0;
	        cout << "tot " << tot << ", passed " << passed;
	    } else {
	        eff = 1.0;
	        cout << "No TPC vertex cut" << endl;
	    }
	}
        cout << ": eff = " << eff;
        cout << endl;
    } else {
	eff = 1.0;
    }
}

void calculateMeanAcceptance(const TEventDataProcessor *eventData, Float_t &acceptanceBTOW, Float_t &acceptanceBPRS, Float_t &acceptanceBSMDE, Float_t &acceptanceBSMDP) {
    acceptanceBTOW = 0.0;
    acceptanceBPRS = 0.0;
    acceptanceBSMDE = 0.0;
    acceptanceBSMDP = 0.0;
    if (eventData) {
	const TH2F *histDayRun = eventData->getDayRun();
	const TH2F *histDayRunAcceptanceBTOW = eventData->getDayRunAcceptanceBTOW();
	const TH2F *histDayRunAcceptanceBPRS = eventData->getDayRunAcceptanceBPRS();
	const TH2F *histDayRunAcceptanceBSMDE = eventData->getDayRunAcceptanceBSMDE();
	const TH2F *histDayRunAcceptanceBSMDP = eventData->getDayRunAcceptanceBSMDP();
	Float_t intDayRun = histDayRun ? histDayRun->Integral(1, histDayRun->GetXaxis()->GetNbins(), 1, histDayRun->GetYaxis()->GetNbins()) : 0;
	Float_t intDayRunAcceptanceBTOW = histDayRunAcceptanceBTOW ? histDayRunAcceptanceBTOW->Integral(1, histDayRunAcceptanceBTOW->GetXaxis()->GetNbins(), 1, histDayRunAcceptanceBTOW->GetYaxis()->GetNbins()) : 0;
	Float_t intDayRunAcceptanceBPRS = histDayRunAcceptanceBPRS ? histDayRunAcceptanceBPRS->Integral(1, histDayRunAcceptanceBPRS->GetXaxis()->GetNbins(), 1, histDayRunAcceptanceBPRS->GetYaxis()->GetNbins()) : 0;
	Float_t intDayRunAcceptanceBSMDE = histDayRunAcceptanceBSMDE ? histDayRunAcceptanceBSMDE->Integral(1, histDayRunAcceptanceBSMDE->GetXaxis()->GetNbins(), 1, histDayRunAcceptanceBSMDE->GetYaxis()->GetNbins()) : 0;
	Float_t intDayRunAcceptanceBSMDP = histDayRunAcceptanceBSMDP ? histDayRunAcceptanceBSMDP->Integral(1, histDayRunAcceptanceBSMDP->GetXaxis()->GetNbins(), 1, histDayRunAcceptanceBSMDP->GetYaxis()->GetNbins()) : 0;
	if (intDayRun > 0) acceptanceBTOW = intDayRunAcceptanceBTOW / intDayRun;
	if (intDayRun > 0) acceptanceBPRS = intDayRunAcceptanceBPRS / intDayRun;
	if (intDayRun > 0) acceptanceBSMDE = intDayRunAcceptanceBSMDE / intDayRun;
	if (intDayRun > 0) acceptanceBSMDP = intDayRunAcceptanceBSMDP / intDayRun;
    }
}

void calculatePrescaleTriggerBias(const TDataProcessorPool *pool, const TAnalysisSettings &settings, TAnalysisResults &results, Int_t step) {
    if (settings.correctPrescaleTriggerBias) {
	if (step == 0) {
	    if (results.prescaleTriggerBiasMBHT1 < 0) {
		cout << "Prescale trigger bias MBHT1 " << settings.title << " calculation..." << endl;
		results.prescaleTriggerBiasMBHT1 = -1.0E6;
	    }
	    if (results.prescaleTriggerBiasMBHT2 < 0) {
		cout << "Prescale trigger bias MBHT2 " << settings.title << " calculation..." << endl;
		results.prescaleTriggerBiasMBHT2 = -1.0E6;
	    }
	    if (results.prescaleTriggerBiasHT1HT2 < 0) {
		cout << "Prescale trigger bias MBHT1HT2 " << settings.title << " calculation..." << endl;
		results.prescaleTriggerBiasHT1HT2 = -1.0E6;
	    }
        } else if (step == 1) {
	    GET_DATA(eventDataProcessorPSMB, const TEventDataProcessor, pool, "eventPSMB");
	    GET_DATA(eventDataProcessorPSMBHT1, const TEventDataProcessor, pool, "eventPSMBHT1");
	    GET_DATA(eventDataProcessorPSMBHT2, const TEventDataProcessor, pool, "eventPSMBHT2");
	    if ((results.prescaleTriggerBiasMBHT1 <= -1.0E6) && eventDataProcessorPSMB && eventDataProcessorPSMBHT1) {
		Float_t evNumPSMB = eventDataProcessorPSMB->numPassedAllCuts;
		Float_t evNumPSMBHT1 = eventDataProcessorPSMBHT1->numPassedAllCuts;
		if (evNumPSMBHT1 != 0) results.prescaleTriggerBiasMBHT1 = -(evNumPSMB / evNumPSMBHT1) - 2.0E6;
		cout << "    Unbiased N_MB = " << evNumPSMB << ", N_MBHT1 = " << evNumPSMBHT1 << endl;
	    }
	    if ((results.prescaleTriggerBiasMBHT2 <= -1.0E6) && eventDataProcessorPSMB && eventDataProcessorPSMBHT2) {
		Float_t evNumPSMB = eventDataProcessorPSMB->numPassedAllCuts;
		Float_t evNumPSMBHT2 = eventDataProcessorPSMBHT2->numPassedAllCuts;
		if (evNumPSMBHT2 != 0) results.prescaleTriggerBiasMBHT2 = -(evNumPSMB / evNumPSMBHT2) - 2.0E6;
		cout << "    Unbiased N_MB = " << evNumPSMB << ", N_MBHT2 = " << evNumPSMBHT2 << endl;
	    }
	    if ((results.prescaleTriggerBiasHT1HT2 <= -1.0E6) && eventDataProcessorPSMBHT1 && eventDataProcessorPSMBHT2) {
		Float_t evNumPSMBHT1 = eventDataProcessorPSMBHT1->numPassedAllCuts;
		Float_t evNumPSMBHT2 = eventDataProcessorPSMBHT2->numPassedAllCuts;
		if (evNumPSMBHT2 != 0) results.prescaleTriggerBiasHT1HT2 = -(evNumPSMBHT1 / evNumPSMBHT2) - 2.0E6;
		cout << "    Unbiased N_MBHT1 = " << evNumPSMBHT1 << ", N_MBHT2 = " << evNumPSMBHT2 << endl;
	    }
	} else if (step == 2) {
	    GET_DATA(eventDataProcessorPSMB, const TEventDataProcessor, pool, "eventPSMB");
	    GET_DATA(eventDataProcessorPSMBHT1, const TEventDataProcessor, pool, "eventPSMBHT1");
	    GET_DATA(eventDataProcessorPSMBHT2, const TEventDataProcessor, pool, "eventPSMBHT2");
	    if ((results.prescaleTriggerBiasMBHT1 <= -2.0E6) && eventDataProcessorPSMB && eventDataProcessorPSMBHT1) {
		Float_t evNumPSMB = eventDataProcessorPSMB->numPassedAllCuts;
		Float_t evNumPSMBHT1 = eventDataProcessorPSMBHT1->numPassedAllCuts;
		if (evNumPSMBHT1 != 0) results.prescaleTriggerBiasMBHT1 = ((evNumPSMB / evNumPSMBHT1) / -(results.prescaleTriggerBiasMBHT1 + 2.0E6)) - 1.0;
		cout << "    Biased N_MB = " << evNumPSMB << ", N_MBHT1 = " << evNumPSMBHT1 << ", BIAS MBHT1 = " << results.prescaleTriggerBiasMBHT1 << endl;
	    }
	    if ((results.prescaleTriggerBiasMBHT2 <= -2.0E6) && eventDataProcessorPSMB && eventDataProcessorPSMBHT2) {
		Float_t evNumPSMB = eventDataProcessorPSMB->numPassedAllCuts;
		Float_t evNumPSMBHT2 = eventDataProcessorPSMBHT2->numPassedAllCuts;
		if (evNumPSMBHT2 != 0) results.prescaleTriggerBiasMBHT2 = ((evNumPSMB / evNumPSMBHT2) / -(results.prescaleTriggerBiasMBHT2 + 2.0E6)) - 1.0;
		cout << "    Biased N_MB = " << evNumPSMB << ", N_MBHT2 = " << evNumPSMBHT2 << ", BIAS MBHT2 = " << results.prescaleTriggerBiasMBHT2 << endl;
	    }
	    if (results.prescaleTriggerBiasHT1HT2 <= -2.0E6) {
		results.prescaleTriggerBiasHT1HT2 = (results.prescaleTriggerBiasHT1HT2 + 2.0E6) - 3.0E6;
	    }
        } else if (step == 3) {
	    GET_DATA(eventDataProcessorPSHT1, const TEventDataProcessor, pool, "eventPSHT1");
	    GET_DATA(eventDataProcessorPSHT1HT2, const TEventDataProcessor, pool, "eventPSHT1HT2");
	    if ((results.prescaleTriggerBiasHT1HT2 <= -3.0E6) && eventDataProcessorPSHT1 && eventDataProcessorPSHT1HT2) {
		Float_t evNumPSHT1 = eventDataProcessorPSHT1->numPassedAllCuts;
		Float_t evNumPSHT1HT2 = eventDataProcessorPSHT1HT2->numPassedAllCuts;
		if (evNumPSHT1HT2 != 0) results.prescaleTriggerBiasHT1HT2 = ((evNumPSHT1 / evNumPSHT1HT2) / -(results.prescaleTriggerBiasHT1HT2 + 3.0E6)) - 1.0;
		cout << "    Biased N_HT1 = " << evNumPSHT1 << ", N_HT1HT2 = " << evNumPSHT1HT2 << ", BIAS HT1HT2 = " << results.prescaleTriggerBiasHT1HT2 << endl;
	    }
        } else if (step == 4) {
	}
    } else {
	results.prescaleTriggerBiasMBHT1 = settings.prescaleTriggerBiasMBHT1;
	results.prescaleTriggerBiasMBHT2 = settings.prescaleTriggerBiasMBHT2;
	results.prescaleTriggerBiasHT1HT2 = settings.prescaleTriggerBiasHT1HT2;
    }
}

void outputCArrays(ostream &ostr, const bin_stat_list_type &values, const Char_t *title, const Char_t *suffix) {
        ostr << "//-------------------------------------" << endl;
        ostr << "// " << title << endl;
        ostr << "//" << endl;
        Bool_t comma_need = false;

#define OUTPUT_ARRAY(NAME, VALUE) \
        ostr << "Float_t " << #NAME << suffix << "[] = { "; \
        comma_need = false; \
        for (bin_stat_list_type::const_iterator iter = values.begin();iter != values.end();iter++) { \
            const TBinStatistics &bin = *iter; \
            if (comma_need) {ostr << ", "; comma_need = false;} \
            ostr << (VALUE); \
            comma_need = true; \
        } \
        ostr << " };" << endl; \

	OUTPUT_ARRAY(x, bin.getParameters().getCenter())
	OUTPUT_ARRAY(x_err_true, bin.getParameters().getWidth() / 2.0)
	OUTPUT_ARRAY(x_err, bin.getParameters().getWidth() * 0.0)
	OUTPUT_ARRAY(y, bin.getValue())
	OUTPUT_ARRAY(y_err, bin.getError())
	
        ostr << "//" << endl;
        ostr << "//-------------------------------------" << endl;
	ostr << endl;
}

void showResults(const TDataProcessorPool *poolMB,  const TDataProcessorPool *poolHT1,  const TDataProcessorPool *poolHT2
	, const TDataProcessorPool *poolEmbedMB, const TDataProcessorPool *poolEmbedHT1, const TDataProcessorPool *poolEmbedHT2
	, const TDataProcessorPool *poolEmbedEtaMB, const TDataProcessorPool *poolEmbedEtaHT1, const TDataProcessorPool *poolEmbedEtaHT2
	, const TDataProcessorPool *poolEmbedEtabgMB, const TDataProcessorPool *poolEmbedEtabgHT1, const TDataProcessorPool *poolEmbedEtabgHT2
	, const TDataProcessorPool *pool1gammaMB, const TDataProcessorPool *pool1gammaHT1, const TDataProcessorPool *pool1gammaHT2
	, const TDataProcessorPool *pool1gammaSimMB, const TDataProcessorPool *pool1gammaSimHT1, const TDataProcessorPool *pool1gammaSimHT2
	, const TDataProcessorPool *pool1gammaSimEtaMB, const TDataProcessorPool *pool1gammaSimEtaHT1, const TDataProcessorPool *pool1gammaSimEtaHT2
	, const TDataProcessorPool *poolPtshiftMB, const TDataProcessorPool *poolPtshiftHT1, const TDataProcessorPool *poolPtshiftHT2
        , const TAnalysisSettings &analysisSettings
	, Bool_t deletePools
        , TAnalysisResults &analysisResults
) {
	TAnalysisSettings analysisSettingsNew = analysisSettings;

	if (analysisSettingsNew.print) {
		cout << "==============================================" << endl;
		cout << "Processing " << analysisSettingsNew.title << endl;
	}

	GET_DATA(candidateDataProcessorMB, const TCandidateDataProcessor, poolMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessorHT1, const TCandidateDataProcessor, poolHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessorHT2, const TCandidateDataProcessor, poolHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(candidateDataProcessorMBNocpv, const TCandidateDataProcessor, poolMB, "candidateMBNocpv");
	if (!candidateDataProcessorMBNocpv) candidateDataProcessorMBNocpv = candidateDataProcessorMB;
	GET_DATA(candidateDataProcessorHT1Nocpv, const TCandidateDataProcessor, poolHT1, "candidateHT1Nocpv");
	if (!candidateDataProcessorHT1Nocpv) candidateDataProcessorHT1Nocpv = candidateDataProcessorHT1;
	GET_DATA(candidateDataProcessorHT2Nocpv, const TCandidateDataProcessor, poolHT2, "candidateHT2Nocpv");
	if (!candidateDataProcessorHT2Nocpv) candidateDataProcessorHT2Nocpv = candidateDataProcessorHT2;

	GET_DATA(candidateDataProcessorMBMixRandom, const TCandidateDataProcessor, poolMB, \
	    analysisSettingsNew.settingsMB.useEventMixingBgRandom ? \
		"candidateMBMix" : \
		(analysisSettingsNew.settingsMB.useEventMixingNotmatchedBgRandom ? \
		    "candidateMBMixNotmatched" : \
		    "" \
		) \
	    );
	GET_DATA(candidateDataProcessorHT1MixRandom, const TCandidateDataProcessor, poolHT1, \
	    analysisSettingsNew.settingsHT1.useEventMixingBgRandom ? \
		"candidateHT1Mix" : \
		(analysisSettingsNew.settingsHT1.useEventMixingNotmatchedBgRandom ? \
		    "candidateHT1MixNotmatched" : \
		    "" \
		) \
	    );
	GET_DATA(candidateDataProcessorHT2MixRandom, const TCandidateDataProcessor, poolHT2, \
	    analysisSettingsNew.settingsHT2.useEventMixingBgRandom ? \
		"candidateHT2Mix" : \
		(analysisSettingsNew.settingsHT2.useEventMixingNotmatchedBgRandom ? \
		    "candidateHT2MixNotmatched" : \
		    "" \
		) \
	    );

	GET_DATA(candidateDataProcessorMBMix, const TCandidateDataProcessor, poolMB, \
	    analysisSettingsNew.settingsMB.useEventJetMixingBg ? \
		"candidateMBJetmix" : \
		(analysisSettingsNew.settingsMB.useEventJetNotmatchedMixingBg ? \
		    "candidateMBJetmixNotmatched" : \
		    (analysisSettingsNew.settingsMB.useEventJetBackMixingBg ? \
			"candidateMBJetmixBack" : \
			(analysisSettingsNew.settingsMB.useEventShufflingBg ? \
			    "candidateMBShuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );
	GET_DATA(candidateDataProcessorHT1Mix, const TCandidateDataProcessor, poolHT1, \
	    analysisSettingsNew.settingsHT1.useEventJetMixingBg ? \
		"candidateHT1Jetmix" : \
		(analysisSettingsNew.settingsHT1.useEventJetNotmatchedMixingBg ? \
		    "candidateHT1JetmixNotmatched" : \
		    (analysisSettingsNew.settingsHT1.useEventJetBackMixingBg ? \
			"candidateHT1JetmixBack" : ( \
			analysisSettingsNew.settingsHT1.useEventShufflingBg ? \
			    "candidateHT1Shuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );
	GET_DATA(candidateDataProcessorHT2Mix, const TCandidateDataProcessor, poolHT2, \
	    analysisSettingsNew.settingsHT2.useEventJetMixingBg ? \
		"candidateHT2Jetmix" : \
		(analysisSettingsNew.settingsHT2.useEventJetNotmatchedMixingBg ? \
		    "candidateHT2JetmixNotmatched" : \
		    (analysisSettingsNew.settingsHT2.useEventJetBackMixingBg ? \
			"candidateHT2JetmixBack" : \
			(analysisSettingsNew.settingsHT2.useEventShufflingBg ? \
			    "candidateHT2Shuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );

	GET_DATA(pointDataProcessorMB, const TPointDataProcessor, poolMB, "pointMB");
	GET_DATA(pointDataProcessorHT1, const TPointDataProcessor, poolHT1, "pointHT1");
	GET_DATA(pointDataProcessorHT2, const TPointDataProcessor, poolHT2, "pointHT2");
//	const TPointDataProcessor *pointDataProcessorMB = candidateDataProcessorMB ? &candidateDataProcessorMB->points : 0;
//	const TPointDataProcessor *pointDataProcessorHT1 = candidateDataProcessorHT1 ? &candidateDataProcessorHT1->points : 0;
//	const TPointDataProcessor *pointDataProcessorHT2 = candidateDataProcessorHT2 ? &candidateDataProcessorHT2->points : 0;

	GET_DATA(eventDataProcessorMB, const TEventDataProcessor, poolMB, "eventMB");
	GET_DATA(eventDataProcessorHT1, const TEventDataProcessor, poolHT1, "eventHT1");
	GET_DATA(eventDataProcessorHT2, const TEventDataProcessor, poolHT2, "eventHT2");

	GET_DATA(candidateDataProcessorEmbedMB, const TCandidateDataProcessor, poolEmbedMB, "candidateMB");
	GET_DATA(candidateDataProcessorEmbedHT1, const TCandidateDataProcessor, poolEmbedHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorEmbedHT2, const TCandidateDataProcessor, poolEmbedHT2, "candidateHT2");

	GET_DATA(candidateDataProcessorEmbedEtaMB, const TCandidateDataProcessor, poolEmbedEtaMB, "candidateMB");
	GET_DATA(candidateDataProcessorEmbedEtaHT1, const TCandidateDataProcessor, poolEmbedEtaHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorEmbedEtaHT2, const TCandidateDataProcessor, poolEmbedEtaHT2, "candidateHT2");

	GET_DATA(candidateDataProcessorEmbedEtabgMB, const TCandidateDataProcessor, poolEmbedEtabgMB, "candidateMB");
	GET_DATA(candidateDataProcessorEmbedEtabgHT1, const TCandidateDataProcessor, poolEmbedEtabgHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorEmbedEtabgHT2, const TCandidateDataProcessor, poolEmbedEtabgHT2, "candidateHT2");

	//GET_DATA(candidateDataProcessorEmbedMBMix, const TCandidateDataProcessor, poolEmbedMB, "candidateMBMix");
	//GET_DATA(candidateDataProcessorEmbedHT1Mix, const TCandidateDataProcessor, poolEmbedHT1, "candidateHT1Mix");
	//GET_DATA(candidateDataProcessorEmbedHT2Mix, const TCandidateDataProcessor, poolEmbedHT2, "candidateHT2Mix");

	GET_DATA(eventDataProcessorEmbedMB, const TEventDataProcessor, poolEmbedMB, "eventMB");
	GET_DATA(eventDataProcessorEmbedHT1, const TEventDataProcessor, poolEmbedHT1, "eventHT1");
	GET_DATA(eventDataProcessorEmbedHT2, const TEventDataProcessor, poolEmbedHT2, "eventHT2");

	GET_DATA(eventDataProcessorEmbedEtaMB, const TEventDataProcessor, poolEmbedEtaMB, "eventMB");
	GET_DATA(eventDataProcessorEmbedEtaHT1, const TEventDataProcessor, poolEmbedEtaHT1, "eventHT1");
	GET_DATA(eventDataProcessorEmbedEtaHT2, const TEventDataProcessor, poolEmbedEtaHT2, "eventHT2");

	GET_DATA(eventDataProcessorEmbedEtabgMB, const TEventDataProcessor, poolEmbedEtabgMB, "eventMB");
	GET_DATA(eventDataProcessorEmbedEtabgHT1, const TEventDataProcessor, poolEmbedEtabgHT1, "eventHT1");
	GET_DATA(eventDataProcessorEmbedEtabgHT2, const TEventDataProcessor, poolEmbedEtabgHT2, "eventHT2");

	//GET_DATA(eventDataProcessorEmbedPSMB, const TEventDataProcessor, poolEmbedMB, "eventPSMB");
	//GET_DATA(eventDataProcessorEmbedPSMBHT1, const TEventDataProcessor, poolEmbedMB, "eventPSMBHT1");
	//GET_DATA(eventDataProcessorEmbedPSMBHT2, const TEventDataProcessor, poolEmbedMB, "eventPSMBHT2");
	//GET_DATA(eventDataProcessorEmbedPSHT1, const TEventDataProcessor, poolEmbedHT1, "eventPSHT1");
	//GET_DATA(eventDataProcessorEmbedPSHT1HT2, const TEventDataProcessor, poolEmbedHT1, "eventPSHT1HT2");

	GET_DATA(simuDataProcessorMB, const TSimuDataProcessor, poolEmbedMB, "simulatedMB");
	GET_DATA(simuDataProcessorHT1, const TSimuDataProcessor, poolEmbedHT1, "simulatedHT1");
	GET_DATA(simuDataProcessorHT2, const TSimuDataProcessor, poolEmbedHT2, "simulatedHT2");

	GET_DATA(simuDataProcessorEtaMB, const TSimuDataProcessor, poolEmbedEtaMB, "simulatedMB");
	GET_DATA(simuDataProcessorEtaHT1, const TSimuDataProcessor, poolEmbedEtaHT1, "simulatedHT1");
	GET_DATA(simuDataProcessorEtaHT2, const TSimuDataProcessor, poolEmbedEtaHT2, "simulatedHT2");

	GET_DATA(simuDataProcessorEtabgMB, const TSimuDataProcessor, poolEmbedEtabgMB, "simulatedMB");
	GET_DATA(simuDataProcessorEtabgHT1, const TSimuDataProcessor, poolEmbedEtabgHT1, "simulatedHT1");
	GET_DATA(simuDataProcessorEtabgHT2, const TSimuDataProcessor, poolEmbedEtabgHT2, "simulatedHT2");

	GET_DATA(pointDataProcessorEmbedMB, const TPointDataProcessor, poolEmbedMB, "pointMB");
	GET_DATA(pointDataProcessorEmbedHT1, const TPointDataProcessor, poolEmbedHT1, "pointHT1");
	GET_DATA(pointDataProcessorEmbedHT2, const TPointDataProcessor, poolEmbedHT2, "pointHT2");

	GET_DATA(pointDataProcessorEmbedEtaMB, const TPointDataProcessor, poolEmbedEtaMB, "pointMB");
	GET_DATA(pointDataProcessorEmbedEtaHT1, const TPointDataProcessor, poolEmbedEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessorEmbedEtaHT2, const TPointDataProcessor, poolEmbedEtaHT2, "pointHT2");

	GET_DATA(pointDataProcessorEmbedEtabgMB, const TPointDataProcessor, poolEmbedEtabgMB, "pointMB");
	GET_DATA(pointDataProcessorEmbedEtabgHT1, const TPointDataProcessor, poolEmbedEtabgHT1, "pointHT1");
	GET_DATA(pointDataProcessorEmbedEtabgHT2, const TPointDataProcessor, poolEmbedEtabgHT2, "pointHT2");

	GET_DATA(candidateDataProcessor1gammaMB, const TCandidateDataProcessor, pool1gammaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1, const TCandidateDataProcessor, pool1gammaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2, const TCandidateDataProcessor, pool1gammaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(candidateDataProcessor1gammaMBSim, const TCandidateDataProcessor, pool1gammaSimMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1Sim, const TCandidateDataProcessor, pool1gammaSimHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2Sim, const TCandidateDataProcessor, pool1gammaSimHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(candidateDataProcessor1gammaMBSimEta, const TCandidateDataProcessor, pool1gammaSimEtaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1SimEta, const TCandidateDataProcessor, pool1gammaSimEtaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2SimEta, const TCandidateDataProcessor, pool1gammaSimEtaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(candidateDataProcessor1gammaMBSimEtabg, const TCandidateDataProcessor, pool1gammaSimEtaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1SimEtabg, const TCandidateDataProcessor, pool1gammaSimEtaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2SimEtabg, const TCandidateDataProcessor, pool1gammaSimEtaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessor1gammaMB, const TPointDataProcessor, pool1gammaMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1, const TPointDataProcessor, pool1gammaHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2, const TPointDataProcessor, pool1gammaHT2, "pointHT2");

	GET_DATA(pointDataProcessor1gammaMBSim, const TPointDataProcessor, pool1gammaSimMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1Sim, const TPointDataProcessor, pool1gammaSimHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2Sim, const TPointDataProcessor, pool1gammaSimHT2, "pointHT2");

	GET_DATA(pointDataProcessor1gammaMBSimEta, const TPointDataProcessor, pool1gammaSimEtaMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1SimEta, const TPointDataProcessor, pool1gammaSimEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2SimEta, const TPointDataProcessor, pool1gammaSimEtaHT2, "pointHT2");

	GET_DATA(pointDataProcessor1gammaMBSimEtabg, const TPointDataProcessor, pool1gammaSimEtaMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1SimEtabg, const TPointDataProcessor, pool1gammaSimEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2SimEtabg, const TPointDataProcessor, pool1gammaSimEtaHT2, "pointHT2");

	GET_DATA(candidateDataProcessorPtshiftMB, const TCandidateDataProcessor, poolPtshiftMB, "candidateMB");
	GET_DATA(candidateDataProcessorPtshiftHT1, const TCandidateDataProcessor, poolPtshiftHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorPtshiftHT2, const TCandidateDataProcessor, poolPtshiftHT2, "candidateHT2");

	GET_DATA(pointDataProcessorPtshiftMB, const TPointDataProcessor, poolPtshiftMB, "pointMB");
	GET_DATA(pointDataProcessorPtshiftHT1, const TPointDataProcessor, poolPtshiftHT1, "pointHT1");
	GET_DATA(pointDataProcessorPtshiftHT2, const TPointDataProcessor, poolPtshiftHT2, "pointHT2");

	GET_DATA(gammaDataProcessorMB, const TMCGammaDataProcessor, pool1gammaMB, "mcGammaMB");
	GET_DATA(gammaDataProcessorHT1, const TMCGammaDataProcessor, pool1gammaHT1, "mcGammaHT1");
	GET_DATA(gammaDataProcessorHT2, const TMCGammaDataProcessor, pool1gammaHT2, "mcGammaHT2");

	clearCachedPrescales();
	clearPtBinPosFuncFromQCD();
	createPtBinPosFuncFromQCD(PQCDFUNC(analysisSettingsNew.usePPpQCD));

	if (analysisSettingsNew.saveRunTimes) {
	    getRunTimes(eventDataProcessorMB, Form(analysisSettingsNew.runTimesFilenameFormat, "MB"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	    getRunTimes(eventDataProcessorHT1, Form(analysisSettingsNew.runTimesFilenameFormat, "HT1"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	    getRunTimes(eventDataProcessorHT2, Form(analysisSettingsNew.runTimesFilenameFormat, "HT2"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	}
	if (analysisSettingsNew.saveRunTimesSim) {
	    getRunTimes(eventDataProcessorEmbedMB, Form(analysisSettingsNew.runTimesSimFilenameFormat, "MB"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	    getRunTimes(eventDataProcessorEmbedHT1, Form(analysisSettingsNew.runTimesSimFilenameFormat, "HT1"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	    getRunTimes(eventDataProcessorEmbedHT2, Form(analysisSettingsNew.runTimesSimFilenameFormat, "HT2"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	}
	if (analysisSettingsNew.saveRunTimesSimEta) {
	    getRunTimes(eventDataProcessorEmbedEtaMB, Form(analysisSettingsNew.runTimesSimEtaFilenameFormat, "MB"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	    getRunTimes(eventDataProcessorEmbedEtaHT1, Form(analysisSettingsNew.runTimesSimEtaFilenameFormat, "HT1"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	    getRunTimes(eventDataProcessorEmbedEtaHT2, Form(analysisSettingsNew.runTimesSimEtaFilenameFormat, "HT2"), analysisSettingsNew.print, analysisSettingsNew.dbConnectString);
	}
	if (analysisSettingsNew.savePrescales) {
	    Int_t *trgMB = 0;
	    Int_t *trgHT1 = 0;
	    Int_t *trgHT2 = 0;
	    if (poolMB && poolHT1 && poolHT2) {
		const Int_t trgSize = sizeof(poolMB->settings.triggers) / sizeof(poolMB->settings.triggers[0]);
		trgMB = new Int_t[trgSize];
		trgHT1 = new Int_t[trgSize];
		trgHT2 = new Int_t[trgSize];
		for (Int_t i = 0;i < trgSize;i++) {
		    if (trgMB) trgMB[i] = 0;
		    if (trgHT1) trgHT1[i] = 0;
		    if (trgHT2) trgHT2[i] = 0;
		}
		for (Int_t i = 0, ind = 0;i < trgSize;i++) if ((poolMB->settings.triggersMB & (1 << i)) && trgMB) trgMB[ind++] = poolMB->settings.triggers[i];
		for (Int_t i = 0, ind = 0;i < trgSize;i++) if ((poolHT1->settings.triggersHT1 & (1 << i)) && trgHT1) trgHT1[ind++] = poolHT1->settings.triggers[i];
		for (Int_t i = 0, ind = 0;i < trgSize;i++) if ((poolHT2->settings.triggersHT2 & (1 << i)) && trgHT2) trgHT2[ind++] = poolHT2->settings.triggers[i];
	    }
	    getPrescales(eventDataProcessorMB, eventDataProcessorHT1, eventDataProcessorHT2, analysisSettingsNew.prescalesSaveFilename, analysisSettingsNew.print, analysisSettingsNew.dbConnectString, trgMB, trgHT1, trgHT2);
	    if (trgMB) delete [] trgMB;
	    if (trgHT1) delete [] trgHT1;
	    if (trgHT2) delete [] trgHT2;
	}

	if ((analysisSettingsNew.NBinaryCollisions < 0) || (analysisSettingsNew.NBinaryCollisionsErr < 0)) {
    	    if (eventDataProcessorMB) getCentralityBin(eventDataProcessorMB->getCuts().getParametersEvent().ftpcRefMultLow, eventDataProcessorMB->getCuts().getParametersEvent().ftpcRefMultHigh, analysisResults.PercentCentralityLow, analysisResults.PercentCentralityHigh, analysisResults.NBinaryCollisions, analysisResults.NBinaryCollisionsErr);
	} else {
	    analysisResults.NBinaryCollisions = analysisSettingsNew.NBinaryCollisions;
	    analysisResults.NBinaryCollisionsErr = analysisSettingsNew.NBinaryCollisionsErr;
	}

	{
	Float_t PSMB = 0, PSHT1 = 0, PSHT2 = 0;
	Float_t PSMB_num = 0, PSHT1_num = 0, PSHT2_num = 0;
	if (analysisSettingsNew.settingsMB.fixedPrescale > 0) {
	    PSMB += analysisSettingsNew.settingsMB.fixedPrescale;
	    PSMB_num += 1;
	}
	if (analysisSettingsNew.settingsHT1.fixedPrescale > 0) {
	    PSHT1 += analysisSettingsNew.settingsHT1.fixedPrescale;
	    PSHT1_num += 1;
	}
	if (analysisSettingsNew.settingsHT2.fixedPrescale > 0) {
	    PSHT2 += analysisSettingsNew.settingsHT2.fixedPrescale;
	    PSHT2_num += 1;
	}
	if (analysisSettingsNew.settingsMB.calculatePrescaleFromDB || analysisSettingsNew.settingsHT1.calculatePrescaleFromDB || analysisSettingsNew.settingsHT2.calculatePrescaleFromDB) {
		Float_t PSHT1_tmp, PSHT2_tmp;
		Float_t PSHT1_tmp_alexst, PSHT2_tmp_alexst;
		calculatePSFromDB(analysisSettingsNew.name + "_data", analysisSettingsNew.title + " (data)", analysisSettingsNew.prescalesDBFilename, eventDataProcessorMB, eventDataProcessorHT1, eventDataProcessorHT2, analysisSettingsNew.showReal, analysisSettingsNew.print, PSHT1_tmp, PSHT2_tmp, PSHT1_tmp_alexst, PSHT2_tmp_alexst);
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromDB) {
		    PSMB += 1.0;
		    PSMB_num += 1;
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromDB && (PSHT1_tmp > 0)) {
		    PSHT1 += PSHT1_tmp;
		    PSHT1_num += 1;
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromDB && (PSHT2_tmp > 0)) {
		    PSHT2 += PSHT2_tmp;
		    PSHT2_num += 1;
		}
	}
	if (analysisSettingsNew.settingsMB.calculatePrescaleFromSim || analysisSettingsNew.settingsHT1.calculatePrescaleFromSim || analysisSettingsNew.settingsHT2.calculatePrescaleFromSim) {
		Float_t PSHT1_tmp, PSHT2_tmp;
		calculatePSFromSim(analysisSettingsNew.name + "_data", analysisSettingsNew.title + " (data)", poolMB, poolHT1
				, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0
				, analysisSettingsNew.showReal, analysisSettingsNew.print, PSHT1_tmp, PSHT2_tmp);
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromSim) {
		    PSMB += 1.0;
		    PSMB_num += 1;
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromSim && (PSHT1_tmp > 0)) {
		    PSHT1 += PSHT1_tmp;
		    PSHT1_num += 1;
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromSim && (PSHT2_tmp > 0)) {
		    PSHT2 += PSHT2_tmp;
		    PSHT2_num += 1;
		}
	}
	if (analysisSettingsNew.settingsMB.calculatePrescaleFromPoints || analysisSettingsNew.settingsHT1.calculatePrescaleFromPoints || analysisSettingsNew.settingsHT2.calculatePrescaleFromPoints || analysisSettingsNew.settingsMB.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsHT1.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsHT2.calculatePrescaleFromPointsIntegral) {
		Float_t PSHT1_tmp, PSHT2_tmp, PSHT1Integral_tmp, PSHT2Integral_tmp;
		calculatePSFromPoints(analysisSettingsNew.name + "_data", analysisSettingsNew.title + " (data)"
		    , poolMB, poolHT1, poolHT2
		    , pool1gammaMB, pool1gammaHT1, pool1gammaHT2
		    , false, false
		    , analysisSettingsNew.settingsHT1.prescalePointsLowPt, analysisSettingsNew.settingsHT1.prescalePointsHighPt
		    , analysisSettingsNew.settingsHT2.prescalePointsLowPt, analysisSettingsNew.settingsHT2.prescalePointsHighPt
		    , 0.3
		    , analysisSettingsNew.settingsHT1.prescalePointsLowPtIntegral, analysisSettingsNew.settingsHT1.prescalePointsHighPtIntegral
		    , analysisSettingsNew.settingsHT2.prescalePointsLowPtIntegral, analysisSettingsNew.settingsHT2.prescalePointsHighPtIntegral
		    , analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0
		    , analysisSettingsNew.showReal, analysisSettingsNew.print
		    , PSHT1_tmp, PSHT2_tmp
		    , PSHT1Integral_tmp, PSHT2Integral_tmp);
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromPoints) {
		    PSMB += 1.0;
		    PSMB_num += 1;
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromPoints && (PSHT1_tmp > 0)) {
		    PSHT1 += PSHT1_tmp;
		    PSHT1_num += 1;
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromPoints && (PSHT2_tmp > 0)) {
		    PSHT2 += PSHT2_tmp;
		    PSHT2_num += 1;
		}
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromPointsIntegral) {
		    PSMB += 1.0;
		    PSMB_num += 1;
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromPointsIntegral && (PSHT1Integral_tmp > 0)) {
		    PSHT1 += PSHT1Integral_tmp;
		    PSHT1_num += 1;
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromPointsIntegral && (PSHT2Integral_tmp > 0)) {
		    PSHT2 += PSHT2Integral_tmp;
		    PSHT2_num += 1;
		}
	}
	if (PSMB_num != 0) PSMB = PSMB / PSMB_num;
	if (PSHT1_num != 0) PSHT1 = PSHT1 / PSHT1_num;
	if (PSHT2_num != 0) PSHT2 = PSHT2 / PSHT2_num;
	if (analysisSettingsNew.print) cout << "HighTower enhancement factors, " << analysisSettingsNew.title << ": MB = " << PSMB << ", HT1 = " << PSHT1 << ", HT2 = " << PSHT2 << endl;
	analysisResults.resultsMB.prescaleUsed = PSMB;
	analysisResults.resultsHT1.prescaleUsed = PSHT1;
	analysisResults.resultsHT2.prescaleUsed = PSHT2;
	}

	{
	Float_t PSSimMB = 0, PSSimHT1 = 0, PSSimHT2 = 0;
	Float_t PSSimMB_num = 0, PSSimHT1_num = 0, PSSimHT2_num = 0;
	if (analysisSettingsNew.settingsSimMB.fixedPrescale > 0) {
	    PSSimMB += analysisSettingsNew.settingsSimMB.fixedPrescale;
	    PSSimMB_num += 1;
	}
	if (analysisSettingsNew.settingsSimHT1.fixedPrescale > 0) {
	    PSSimHT1 += analysisSettingsNew.settingsSimHT1.fixedPrescale;
	    PSSimHT1_num += 1;
	}
	if (analysisSettingsNew.settingsSimHT2.fixedPrescale > 0) {
	    PSSimHT2 += analysisSettingsNew.settingsSimHT2.fixedPrescale;
	    PSSimHT2_num += 1;
	}
	if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromDB || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromDB || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromDB) {
		Float_t PSSimHT1_tmp, PSSimHT2_tmp;
		Float_t PSSimHT1_tmp_alexst, PSSimHT2_tmp_alexst;
		calculatePSFromDB(analysisSettingsNew.name + "_sim", analysisSettingsNew.title + " (sim)", analysisSettingsNew.prescalesDBFilename, eventDataProcessorEmbedMB, eventDataProcessorEmbedHT1, eventDataProcessorEmbedHT2, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimHT1_tmp, PSSimHT2_tmp, PSSimHT1_tmp_alexst, PSSimHT2_tmp_alexst);
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromDB) {
		    PSSimMB += 1.0;
		    PSSimMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromDB && (PSSimHT1_tmp > 0)) {
		    PSSimHT1 += PSSimHT1_tmp;
		    PSSimHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromDB && (PSSimHT2_tmp > 0)) {
		    PSSimHT2 += PSSimHT2_tmp;
		    PSSimHT2_num += 1;
		}
	}
	if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromSim || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromSim || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromSim) {
		Float_t PSSimHT1_tmp, PSSimHT2_tmp;
		calculatePSFromSim(analysisSettingsNew.name + "_sim", analysisSettingsNew.title + " (sim)", poolEmbedMB, poolEmbedHT1
				, 0, 0, 0/*, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0*/
				, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimHT1_tmp, PSSimHT2_tmp);
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromSim) {
		    PSSimMB += 1.0;
		    PSSimMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromSim && (PSSimHT1_tmp > 0)) {
		    PSSimHT1 += PSSimHT1_tmp;
		    PSSimHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromSim && (PSSimHT2_tmp > 0)) {
		    PSSimHT2 += PSSimHT2_tmp;
		    PSSimHT2_num += 1;
		}
	}
	if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimMB.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPointsIntegral) {
		Float_t PSSimHT1_tmp, PSSimHT2_tmp, PSSimHT1Integral_tmp, PSSimHT2Integral_tmp;
		calculatePSFromPoints(analysisSettingsNew.name + "_sim", analysisSettingsNew.title + " (sim)"
		    , poolEmbedMB, poolEmbedHT1, poolEmbedHT2
		    , pool1gammaSimMB, pool1gammaSimHT1, pool1gammaSimHT2
		    , false, false
		    , analysisSettingsNew.settingsSimHT1.prescalePointsLowPt, analysisSettingsNew.settingsSimHT1.prescalePointsHighPt
		    , analysisSettingsNew.settingsSimHT2.prescalePointsLowPt, analysisSettingsNew.settingsSimHT2.prescalePointsHighPt
		    , 0.3
		    , analysisSettingsNew.settingsSimHT1.prescalePointsLowPtIntegral, analysisSettingsNew.settingsSimHT1.prescalePointsHighPtIntegral
		    , analysisSettingsNew.settingsSimHT2.prescalePointsLowPtIntegral, analysisSettingsNew.settingsSimHT2.prescalePointsHighPtIntegral
		    , 0, 0, 0/*, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0*/
		    , analysisSettingsNew.showSim, analysisSettingsNew.print
		    , PSSimHT1_tmp, PSSimHT2_tmp
		    , PSSimHT1Integral_tmp, PSSimHT2Integral_tmp);
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromPoints) {
		    PSSimMB += 1.0;
		    PSSimMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPoints && (PSSimHT1_tmp > 0)) {
		    PSSimHT1 += PSSimHT1_tmp;
		    PSSimHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPoints && (PSSimHT2_tmp > 0)) {
		    PSSimHT2 += PSSimHT2_tmp;
		    PSSimHT2_num += 1;
		}
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromPointsIntegral) {
		    PSSimMB += 1.0;
		    PSSimMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPointsIntegral && (PSSimHT1Integral_tmp > 0)) {
		    PSSimHT1 += PSSimHT1Integral_tmp;
		    PSSimHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPointsIntegral && (PSSimHT2Integral_tmp > 0)) {
		    PSSimHT2 += PSSimHT2Integral_tmp;
		    PSSimHT2_num += 1;
		}
	}
	if (PSSimMB_num != 0) PSSimMB = PSSimMB / PSSimMB_num;
	if (PSSimHT1_num != 0) PSSimHT1 = PSSimHT1 / PSSimHT1_num;
	if (PSSimHT2_num != 0) PSSimHT2 = PSSimHT2 / PSSimHT2_num;
	if (analysisSettingsNew.print) cout << "HighTower enhancement factors (sim), " << analysisSettingsNew.title << ": MB = " << PSSimMB << ", HT1 = " << PSSimHT1 << ", HT2 = " << PSSimHT2 << endl;
	analysisResults.resultsSimMB.prescaleUsed = PSSimMB;
	analysisResults.resultsSimHT1.prescaleUsed = PSSimHT1;
	analysisResults.resultsSimHT2.prescaleUsed = PSSimHT2;
	}

	{
	Float_t PSSimEtaMB = 0, PSSimEtaHT1 = 0, PSSimEtaHT2 = 0;
	Float_t PSSimEtaMB_num = 0, PSSimEtaHT1_num = 0, PSSimEtaHT2_num = 0;
	if (analysisSettingsNew.settingsSimEtaMB.fixedPrescale > 0) {
	    PSSimEtaMB += analysisSettingsNew.settingsSimEtaMB.fixedPrescale;
	    PSSimEtaMB_num += 1;
	}
	if (analysisSettingsNew.settingsSimEtaHT1.fixedPrescale > 0) {
	    PSSimEtaHT1 += analysisSettingsNew.settingsSimEtaHT1.fixedPrescale;
	    PSSimEtaHT1_num += 1;
	}
	if (analysisSettingsNew.settingsSimEtaHT2.fixedPrescale > 0) {
	    PSSimEtaHT2 += analysisSettingsNew.settingsSimEtaHT2.fixedPrescale;
	    PSSimEtaHT2_num += 1;
	}
	if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromDB || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromDB || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromDB) {
		Float_t PSSimEtaHT1_tmp, PSSimEtaHT2_tmp;
		Float_t PSSimEtaHT1_tmp_alexst, PSSimEtaHT2_tmp_alexst;
		calculatePSFromDB(analysisSettingsNew.name + "_simEta", analysisSettingsNew.title + " (sim eta)", analysisSettingsNew.prescalesDBFilename, eventDataProcessorEmbedEtaMB, eventDataProcessorEmbedEtaHT1, eventDataProcessorEmbedEtaHT2, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimEtaHT1_tmp, PSSimEtaHT2_tmp, PSSimEtaHT1_tmp_alexst, PSSimEtaHT2_tmp_alexst);
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromDB) {
		    PSSimEtaMB += 1.0;
		    PSSimEtaMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromDB && (PSSimEtaHT1_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1_tmp;
		    PSSimEtaHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromDB && (PSSimEtaHT2_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2_tmp;
		    PSSimEtaHT2_num += 1;
		}
	}
	if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromSim || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromSim || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromSim) {
		Float_t PSSimEtaHT1_tmp, PSSimEtaHT2_tmp;
		calculatePSFromSim(analysisSettingsNew.name + "_simEta", analysisSettingsNew.title + " (sim eta)", poolEmbedEtaMB, poolEmbedEtaHT1
				, 0, 0, 0/*, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0*/
				, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimEtaHT1_tmp, PSSimEtaHT2_tmp);
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromSim) {
		    PSSimEtaMB += 1.0;
		    PSSimEtaMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromSim && (PSSimEtaHT1_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1_tmp;
		    PSSimEtaHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromSim && (PSSimEtaHT2_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2_tmp;
		    PSSimEtaHT2_num += 1;
		}
	}
	if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPointsIntegral) {
		Float_t PSSimEtaHT1_tmp, PSSimEtaHT2_tmp, PSSimEtaHT1Integral_tmp, PSSimEtaHT2Integral_tmp;
		calculatePSFromPoints(analysisSettingsNew.name + "_simEta", analysisSettingsNew.title + " (sim eta)"
		    , poolEmbedEtaMB, poolEmbedEtaHT1, poolEmbedEtaHT2
		    , pool1gammaSimEtaMB, pool1gammaSimEtaHT1, pool1gammaSimEtaHT2
		    , false, false
		    , analysisSettingsNew.settingsSimEtaHT1.prescalePointsLowPt, analysisSettingsNew.settingsSimEtaHT1.prescalePointsHighPt
		    , analysisSettingsNew.settingsSimEtaHT2.prescalePointsLowPt, analysisSettingsNew.settingsSimEtaHT2.prescalePointsHighPt
		    , 0.3
		    , analysisSettingsNew.settingsSimEtaHT1.prescalePointsLowPtIntegral, analysisSettingsNew.settingsSimEtaHT1.prescalePointsHighPtIntegral
		    , analysisSettingsNew.settingsSimEtaHT2.prescalePointsLowPtIntegral, analysisSettingsNew.settingsSimEtaHT2.prescalePointsHighPtIntegral
		    , 0, 0, 0/*, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0*/
		    , analysisSettingsNew.showSim, analysisSettingsNew.print
		    , PSSimEtaHT1_tmp, PSSimEtaHT2_tmp
		    , PSSimEtaHT1Integral_tmp, PSSimEtaHT2Integral_tmp);
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPoints) {
		    PSSimEtaMB += 1.0;
		    PSSimEtaMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPoints && (PSSimEtaHT1_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1_tmp;
		    PSSimEtaHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPoints && (PSSimEtaHT2_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2_tmp;
		    PSSimEtaHT2_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPointsIntegral) {
		    PSSimEtaMB += 1.0;
		    PSSimEtaMB_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPointsIntegral && (PSSimEtaHT1Integral_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1Integral_tmp;
		    PSSimEtaHT1_num += 1;
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPointsIntegral && (PSSimEtaHT2Integral_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2Integral_tmp;
		    PSSimEtaHT2_num += 1;
		}
	}
	if (PSSimEtaMB_num != 0) PSSimEtaMB = PSSimEtaMB / PSSimEtaMB_num;
	if (PSSimEtaHT1_num != 0) PSSimEtaHT1 = PSSimEtaHT1 / PSSimEtaHT1_num;
	if (PSSimEtaHT2_num != 0) PSSimEtaHT2 = PSSimEtaHT2 / PSSimEtaHT2_num;
	if (analysisSettingsNew.print) cout << "HighTower enhancement factors (sim eta), " << analysisSettingsNew.title << ": MB = " << PSSimEtaMB << ", HT1 = " << PSSimEtaHT1 << ", HT2 = " << PSSimEtaHT2 << endl;
	analysisResults.resultsSimEtaMB.prescaleUsed = PSSimEtaMB;
	analysisResults.resultsSimEtaHT1.prescaleUsed = PSSimEtaHT1;
	analysisResults.resultsSimEtaHT2.prescaleUsed = PSSimEtaHT2;
	}

	analysisSettingsNew.settingsMB.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsMB.name;
	analysisSettingsNew.settingsMB.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsMB.title;
	analysisSettingsNew.settingsMB.showBins &= analysisSettingsNew.showReal;

	analysisSettingsNew.settingsHT1.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsHT1.name;
	analysisSettingsNew.settingsHT1.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsHT1.title;
	analysisSettingsNew.settingsHT1.showBins &= analysisSettingsNew.showReal;

	analysisSettingsNew.settingsHT2.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsHT2.name;
	analysisSettingsNew.settingsHT2.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsHT2.title;
	analysisSettingsNew.settingsHT2.showBins &= analysisSettingsNew.showReal;

	analysisSettingsNew.settingsSimMB.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsSimMB.name;
	analysisSettingsNew.settingsSimMB.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsSimMB.title;
	analysisSettingsNew.settingsSimMB.showBins &= analysisSettingsNew.showSim;

	analysisSettingsNew.settingsSimHT1.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsSimHT1.name;
	analysisSettingsNew.settingsSimHT1.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsSimHT1.title;
	analysisSettingsNew.settingsSimHT1.showBins &= analysisSettingsNew.showSim;

	analysisSettingsNew.settingsSimHT2.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsSimHT2.name;
	analysisSettingsNew.settingsSimHT2.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsSimHT2.title;
	analysisSettingsNew.settingsSimHT2.showBins &= analysisSettingsNew.showSim;

	analysisSettingsNew.settingsSimEtaMB.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsSimEtaMB.name;
	analysisSettingsNew.settingsSimEtaMB.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsSimEtaMB.title;
	analysisSettingsNew.settingsSimEtaMB.showBins &= analysisSettingsNew.showSim;

	analysisSettingsNew.settingsSimEtaHT1.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsSimEtaHT1.name;
	analysisSettingsNew.settingsSimEtaHT1.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsSimEtaHT1.title;
	analysisSettingsNew.settingsSimEtaHT1.showBins &= analysisSettingsNew.showSim;

	analysisSettingsNew.settingsSimEtaHT2.name = analysisSettingsNew.name + "_" + analysisSettingsNew.settingsSimEtaHT2.name;
	analysisSettingsNew.settingsSimEtaHT2.title = analysisSettingsNew.title + " " + analysisSettingsNew.settingsSimEtaHT2.title;
	analysisSettingsNew.settingsSimEtaHT2.showBins &= analysisSettingsNew.showSim;

	if (analysisSettingsNew.settingsMB.correctVertexFindingEff && (analysisSettingsNew.settingsMB.vertexFindingEff > 0)) {
	    analysisResults.resultsMB.vertexFindingEff = analysisSettingsNew.settingsMB.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsMB.title, analysisSettingsNew.settingsMB.correctVertexFindingEff, eventDataProcessorMB, analysisResults.resultsMB.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsHT1.correctVertexFindingEff && (analysisSettingsNew.settingsHT1.vertexFindingEff > 0)) {
	    analysisResults.resultsHT1.vertexFindingEff = analysisSettingsNew.settingsHT1.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsHT1.title, analysisSettingsNew.settingsHT1.correctVertexFindingEff, eventDataProcessorHT1, analysisResults.resultsHT1.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsHT2.correctVertexFindingEff && (analysisSettingsNew.settingsHT2.vertexFindingEff > 0)) {
	    analysisResults.resultsHT2.vertexFindingEff = analysisSettingsNew.settingsHT2.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsHT2.title, analysisSettingsNew.settingsHT2.correctVertexFindingEff, eventDataProcessorHT2, analysisResults.resultsHT2.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsHT1.correctVertexFindingEffMB) analysisResults.resultsHT1.vertexFindingEff *= analysisResults.resultsMB.vertexFindingEff;
	if (analysisSettingsNew.settingsHT2.correctVertexFindingEffMB) analysisResults.resultsHT2.vertexFindingEff *= analysisResults.resultsMB.vertexFindingEff;

	if (analysisSettingsNew.settingsSimMB.correctVertexFindingEff && (analysisSettingsNew.settingsSimMB.vertexFindingEff > 0)) {
	    analysisResults.resultsSimMB.vertexFindingEff = analysisSettingsNew.settingsSimMB.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsSimMB.title, analysisSettingsNew.settingsSimMB.correctVertexFindingEff, eventDataProcessorEmbedMB, analysisResults.resultsSimMB.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsSimHT1.correctVertexFindingEff && (analysisSettingsNew.settingsSimHT1.vertexFindingEff > 0)) {
	    analysisResults.resultsSimHT1.vertexFindingEff = analysisSettingsNew.settingsSimHT1.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsSimHT1.title, analysisSettingsNew.settingsSimHT1.correctVertexFindingEff, eventDataProcessorEmbedHT1, analysisResults.resultsSimHT1.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsSimHT2.correctVertexFindingEff && (analysisSettingsNew.settingsSimHT2.vertexFindingEff > 0)) {
	    analysisResults.resultsSimHT2.vertexFindingEff = analysisSettingsNew.settingsSimHT2.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsSimHT2.title, analysisSettingsNew.settingsSimHT2.correctVertexFindingEff, eventDataProcessorEmbedHT2, analysisResults.resultsSimHT2.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsSimHT1.correctVertexFindingEffMB) analysisResults.resultsSimHT1.vertexFindingEff *= analysisResults.resultsSimMB.vertexFindingEff;
	if (analysisSettingsNew.settingsSimHT2.correctVertexFindingEffMB) analysisResults.resultsSimHT2.vertexFindingEff *= analysisResults.resultsSimMB.vertexFindingEff;

	if (analysisSettingsNew.settingsSimEtaMB.correctVertexFindingEff && (analysisSettingsNew.settingsSimEtaMB.vertexFindingEff > 0)) {
	    analysisResults.resultsSimEtaMB.vertexFindingEff = analysisSettingsNew.settingsSimEtaMB.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsSimEtaMB.title, analysisSettingsNew.settingsSimEtaMB.correctVertexFindingEff, eventDataProcessorEmbedEtaMB, analysisResults.resultsSimEtaMB.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsSimEtaHT1.correctVertexFindingEff && (analysisSettingsNew.settingsSimEtaHT1.vertexFindingEff > 0)) {
	    analysisResults.resultsSimEtaHT1.vertexFindingEff = analysisSettingsNew.settingsSimEtaHT1.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsSimEtaHT1.title, analysisSettingsNew.settingsSimEtaHT1.correctVertexFindingEff, eventDataProcessorEmbedEtaHT1, analysisResults.resultsSimEtaHT1.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsSimEtaHT2.correctVertexFindingEff && (analysisSettingsNew.settingsSimEtaHT2.vertexFindingEff > 0)) {
	    analysisResults.resultsSimEtaHT2.vertexFindingEff = analysisSettingsNew.settingsSimEtaHT2.vertexFindingEff;
	} else {
	    calculateVertexFindingEff(analysisSettingsNew.settingsSimEtaHT2.title, analysisSettingsNew.settingsSimEtaHT2.correctVertexFindingEff, eventDataProcessorEmbedEtaHT2, analysisResults.resultsSimEtaHT2.vertexFindingEff);
	}
	if (analysisSettingsNew.settingsSimEtaHT1.correctVertexFindingEffMB) analysisResults.resultsSimEtaHT1.vertexFindingEff *= analysisResults.resultsSimEtaMB.vertexFindingEff;
	if (analysisSettingsNew.settingsSimEtaHT2.correctVertexFindingEffMB) analysisResults.resultsSimEtaHT2.vertexFindingEff *= analysisResults.resultsSimEtaMB.vertexFindingEff;

	selectInv(candidateDataProcessorMB, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	selectInv(candidateDataProcessorMB, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	selectInv(candidateDataProcessorHT1, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	selectInv(candidateDataProcessorHT1, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	selectInv(candidateDataProcessorHT2, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0, true);
	selectInv(candidateDataProcessorHT2, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0, true);
	
	if (!analysisSettingsNew.multiplicityDistributions) {
	    selectInv(candidateDataProcessorMBMix, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorMBMix, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorHT1Mix, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT1Mix, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT2Mix, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0, true);
	    selectInv(candidateDataProcessorHT2Mix, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0, true);
	    selectInv(candidateDataProcessorMBMixRandom, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorMBMixRandom, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorHT1MixRandom, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT1MixRandom, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT2MixRandom, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0, true);
	    selectInv(candidateDataProcessorHT2MixRandom, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0, true);
	}

	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB.lowPtBinSize, analysisSettingsNew.settingsSimMB.lowPt, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invlist);
	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB.highPtBinSize, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.settingsSimMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invlist);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1.lowPtBinSize, analysisSettingsNew.settingsSimHT1.lowPt, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invlist);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1.highPtBinSize, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.settingsSimHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invlist);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2.lowPtBinSize, analysisSettingsNew.settingsSimHT2.lowPt, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invlist);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2.highPtBinSize, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.settingsSimHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invlist);

	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB.lowPtBinSize, analysisSettingsNew.settingsSimMB.lowPt, analysisSettingsNew.settingsSimMB.switchBinSizePt, 5, analysisResults.resultsSimMB.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB.highPtBinSize, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.settingsSimMB.highPt, 5, analysisResults.resultsSimMB.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1.lowPtBinSize, analysisSettingsNew.settingsSimHT1.lowPt, analysisSettingsNew.settingsSimHT1.switchBinSizePt, 5, analysisResults.resultsSimHT1.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1.highPtBinSize, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.settingsSimHT1.highPt, 5, analysisResults.resultsSimHT1.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2.lowPtBinSize, analysisSettingsNew.settingsSimHT2.lowPt, analysisSettingsNew.settingsSimHT2.switchBinSizePt, 5, analysisResults.resultsSimHT2.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2.highPtBinSize, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.settingsSimHT2.highPt, 5, analysisResults.resultsSimHT2.simulatedPtlist);

	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invlist);
	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invlist);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invlist);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invlist);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invlist);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invlist);

	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, 5, analysisResults.resultsSimEtaMB.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, 5, analysisResults.resultsSimEtaMB.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, 5, analysisResults.resultsSimEtaHT1.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, 5, analysisResults.resultsSimEtaHT1.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, 5, analysisResults.resultsSimEtaHT2.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, 5, analysisResults.resultsSimEtaHT2.simulatedPtlist);

	selectInv(candidateDataProcessorEmbedEtabgMB, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.invlist);
	selectInv(candidateDataProcessorEmbedEtabgMB, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.invlist);
	selectInv(candidateDataProcessorEmbedEtabgHT1, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.invlist);
	selectInv(candidateDataProcessorEmbedEtabgHT1, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.invlist);
	selectInv(candidateDataProcessorEmbedEtabgHT2, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.invlist);
	selectInv(candidateDataProcessorEmbedEtabgHT2, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.invlist);

	selectInv(candidateDataProcessor1gammaMB, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMB : 0);
	selectInv(candidateDataProcessor1gammaMB, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMB : 0);
	selectInv(candidateDataProcessor1gammaHT1, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1 : 0);
	selectInv(candidateDataProcessor1gammaHT1, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1 : 0);
	selectInv(candidateDataProcessor1gammaHT2, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2 : 0);
	selectInv(candidateDataProcessor1gammaHT2, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2 : 0);

	selectInv(candidateDataProcessor1gammaMBSim, analysisSettingsNew.settingsSimMB.lowPtBinSize, analysisSettingsNew.settingsSimMB.lowPt, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSim : 0);
	selectInv(candidateDataProcessor1gammaMBSim, analysisSettingsNew.settingsSimMB.highPtBinSize, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.settingsSimMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSim : 0);
	selectInv(candidateDataProcessor1gammaHT1Sim, analysisSettingsNew.settingsSimHT1.lowPtBinSize, analysisSettingsNew.settingsSimHT1.lowPt, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1Sim : 0);
	selectInv(candidateDataProcessor1gammaHT1Sim, analysisSettingsNew.settingsSimHT1.highPtBinSize, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.settingsSimHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1Sim : 0);
	selectInv(candidateDataProcessor1gammaHT2Sim, analysisSettingsNew.settingsSimHT2.lowPtBinSize, analysisSettingsNew.settingsSimHT2.lowPt, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2Sim : 0);
	selectInv(candidateDataProcessor1gammaHT2Sim, analysisSettingsNew.settingsSimHT2.highPtBinSize, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.settingsSimHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2Sim : 0);

	selectInv(candidateDataProcessor1gammaMBSimEta, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSimEta : 0);
	selectInv(candidateDataProcessor1gammaMBSimEta, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSimEta : 0);
	selectInv(candidateDataProcessor1gammaHT1SimEta, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1SimEta : 0);
	selectInv(candidateDataProcessor1gammaHT1SimEta, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1SimEta : 0);
	selectInv(candidateDataProcessor1gammaHT2SimEta, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2SimEta : 0);
	selectInv(candidateDataProcessor1gammaHT2SimEta, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2SimEta : 0);

	selectInv(candidateDataProcessor1gammaMBSimEtabg, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSimEtabg : 0);
	selectInv(candidateDataProcessor1gammaMBSimEtabg, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSimEtabg : 0);
	selectInv(candidateDataProcessor1gammaHT1SimEtabg, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1SimEtabg : 0);
	selectInv(candidateDataProcessor1gammaHT1SimEtabg, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1SimEtabg : 0);
	selectInv(candidateDataProcessor1gammaHT2SimEtabg, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2SimEtabg : 0);
	selectInv(candidateDataProcessor1gammaHT2SimEtabg, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2SimEtabg : 0);

	selectInv(candidateDataProcessorMB, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, 3, analysisResults.resultsMB.invNeutralMultlist, 0);
	selectInv(candidateDataProcessorMB, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, 3, analysisResults.resultsMB.invNeutralMultlist, 0);
	selectInv(candidateDataProcessorHT1, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, 3, analysisResults.resultsHT1.invNeutralMultlist, 0);
	selectInv(candidateDataProcessorHT1, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, 3, analysisResults.resultsHT1.invNeutralMultlist, 0);
	selectInv(candidateDataProcessorHT2, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, 3, analysisResults.resultsHT2.invNeutralMultlist, 0);
	selectInv(candidateDataProcessorHT2, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, 3, analysisResults.resultsHT2.invNeutralMultlist, 0);
	
	selectInv(candidateDataProcessorMBNocpv, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, 4, analysisResults.resultsMB.trackdistlist, 0);
	selectInv(candidateDataProcessorMBNocpv, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, 4, analysisResults.resultsMB.trackdistlist, 0);
	selectInv(candidateDataProcessorHT1Nocpv, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, 4, analysisResults.resultsHT1.trackdistlist, 0);
	selectInv(candidateDataProcessorHT1Nocpv, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, 4, analysisResults.resultsHT1.trackdistlist, 0);
	selectInv(candidateDataProcessorHT2Nocpv, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, 4, analysisResults.resultsHT2.trackdistlist, 0);
	selectInv(candidateDataProcessorHT2Nocpv, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, 4, analysisResults.resultsHT2.trackdistlist, 0);
	
	selectInv(candidateDataProcessorMBMix, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, 4, analysisResults.resultsMB.trackdistBglist, 0);
	selectInv(candidateDataProcessorMBMix, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, 4, analysisResults.resultsMB.trackdistBglist, 0);
	selectInv(candidateDataProcessorHT1Mix, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, 4, analysisResults.resultsHT1.trackdistBglist, 0);
	selectInv(candidateDataProcessorHT1Mix, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, 4, analysisResults.resultsHT1.trackdistBglist, 0);
	selectInv(candidateDataProcessorHT2Mix, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, 4, analysisResults.resultsHT2.trackdistBglist, 0);
	selectInv(candidateDataProcessorHT2Mix, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, 4, analysisResults.resultsHT2.trackdistBglist, 0);
	
	selectInv(candidateDataProcessorMBMixRandom, analysisSettingsNew.settingsMB.lowPtBinSize, analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.switchBinSizePt, 4, analysisResults.resultsMB.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorMBMixRandom, analysisSettingsNew.settingsMB.highPtBinSize, analysisSettingsNew.settingsMB.switchBinSizePt, analysisSettingsNew.settingsMB.highPt, 4, analysisResults.resultsMB.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorHT1MixRandom, analysisSettingsNew.settingsHT1.lowPtBinSize, analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.switchBinSizePt, 4, analysisResults.resultsHT1.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorHT1MixRandom, analysisSettingsNew.settingsHT1.highPtBinSize, analysisSettingsNew.settingsHT1.switchBinSizePt, analysisSettingsNew.settingsHT1.highPt, 4, analysisResults.resultsHT1.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorHT2MixRandom, analysisSettingsNew.settingsHT2.lowPtBinSize, analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.switchBinSizePt, 4, analysisResults.resultsHT2.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorHT2MixRandom, analysisSettingsNew.settingsHT2.highPtBinSize, analysisSettingsNew.settingsHT2.switchBinSizePt, analysisSettingsNew.settingsHT2.highPt, 4, analysisResults.resultsHT2.trackdistBgRandomlist, 0);
	
	selectInv(candidateDataProcessorPtshiftMB, analysisSettingsNew.settingsSimMB.lowPtBinSize, analysisSettingsNew.settingsSimMB.lowPt, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftMB : 0);
	selectInv(candidateDataProcessorPtshiftMB, analysisSettingsNew.settingsSimMB.highPtBinSize, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.settingsSimMB.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftMB : 0);
	selectInv(candidateDataProcessorPtshiftHT1, analysisSettingsNew.settingsSimHT1.lowPtBinSize, analysisSettingsNew.settingsSimHT1.lowPt, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftHT1 : 0);
	selectInv(candidateDataProcessorPtshiftHT1, analysisSettingsNew.settingsSimHT1.highPtBinSize, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.settingsSimHT1.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftHT1 : 0);
	selectInv(candidateDataProcessorPtshiftHT2, analysisSettingsNew.settingsSimHT2.lowPtBinSize, analysisSettingsNew.settingsSimHT2.lowPt, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftHT2 : 0);
	selectInv(candidateDataProcessorPtshiftHT2, analysisSettingsNew.settingsSimHT2.highPtBinSize, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.settingsSimHT2.highPt, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftHT2 : 0);
	
	selectInvSim(simuDataProcessorMB, analysisSettingsNew.settingsSimMB.lowPtBinSize, analysisSettingsNew.settingsSimMB.lowPt, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisResults.areaSimuMB);
	selectInvSim(simuDataProcessorMB, analysisSettingsNew.settingsSimMB.highPtBinSize, analysisSettingsNew.settingsSimMB.switchBinSizePt, analysisSettingsNew.settingsSimMB.highPt, analysisResults.areaSimuMB);
	selectInvSim(simuDataProcessorHT1, analysisSettingsNew.settingsSimHT1.lowPtBinSize, analysisSettingsNew.settingsSimHT1.lowPt, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisResults.areaSimuHT1);
	selectInvSim(simuDataProcessorHT1, analysisSettingsNew.settingsSimHT1.highPtBinSize, analysisSettingsNew.settingsSimHT1.switchBinSizePt, analysisSettingsNew.settingsSimHT1.highPt, analysisResults.areaSimuHT1);
	selectInvSim(simuDataProcessorHT2, analysisSettingsNew.settingsSimHT2.lowPtBinSize, analysisSettingsNew.settingsSimHT2.lowPt, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisResults.areaSimuHT2);
	selectInvSim(simuDataProcessorHT2, analysisSettingsNew.settingsSimHT2.highPtBinSize, analysisSettingsNew.settingsSimHT2.switchBinSizePt, analysisSettingsNew.settingsSimHT2.highPt, analysisResults.areaSimuHT2);

	selectInvSim(simuDataProcessorEtaMB, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisResults.areaSimuEtaMB);
	selectInvSim(simuDataProcessorEtaMB, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, analysisResults.areaSimuEtaMB);
	selectInvSim(simuDataProcessorEtaHT1, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisResults.areaSimuEtaHT1);
	selectInvSim(simuDataProcessorEtaHT1, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, analysisResults.areaSimuEtaHT1);
	selectInvSim(simuDataProcessorEtaHT2, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisResults.areaSimuEtaHT2);
	selectInvSim(simuDataProcessorEtaHT2, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, analysisResults.areaSimuEtaHT2);

	selectInvSim(simuDataProcessorEtabgMB, analysisSettingsNew.settingsSimEtaMB.lowPtBinSize, analysisSettingsNew.settingsSimEtaMB.lowPt, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisResults.areaSimuEtabgMB);
	selectInvSim(simuDataProcessorEtabgMB, analysisSettingsNew.settingsSimEtaMB.highPtBinSize, analysisSettingsNew.settingsSimEtaMB.switchBinSizePt, analysisSettingsNew.settingsSimEtaMB.highPt, analysisResults.areaSimuEtabgMB);
	selectInvSim(simuDataProcessorEtabgHT1, analysisSettingsNew.settingsSimEtaHT1.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT1.lowPt, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisResults.areaSimuEtabgHT1);
	selectInvSim(simuDataProcessorEtabgHT1, analysisSettingsNew.settingsSimEtaHT1.highPtBinSize, analysisSettingsNew.settingsSimEtaHT1.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT1.highPt, analysisResults.areaSimuEtabgHT1);
	selectInvSim(simuDataProcessorEtabgHT2, analysisSettingsNew.settingsSimEtaHT2.lowPtBinSize, analysisSettingsNew.settingsSimEtaHT2.lowPt, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisResults.areaSimuEtabgHT2);
	selectInvSim(simuDataProcessorEtabgHT2, analysisSettingsNew.settingsSimEtaHT2.highPtBinSize, analysisSettingsNew.settingsSimEtaHT2.switchBinSizePt, analysisSettingsNew.settingsSimEtaHT2.highPt, analysisResults.areaSimuEtabgHT2);

	analysisResults.resultsMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleMB", analysisSettingsNew.title + ", low-mass norm MB", pointDataProcessorMB, pointDataProcessor1gammaMB, analysisSettingsNew.settingsMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsMB.lowmassbgPointsHighPt, analysisResults.resultsMB.invlist, analysisResults.resultsMB.lowNormPointsPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	analysisResults.resultsHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleHT1", analysisSettingsNew.title + ", low-mass norm HT1", pointDataProcessorHT1, pointDataProcessor1gammaHT1, analysisSettingsNew.settingsHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT1.lowmassbgPointsHighPt, analysisResults.resultsHT1.invlist, analysisResults.resultsHT1.lowNormPointsPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	analysisResults.resultsHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleHT2", analysisSettingsNew.title + ", low-mass norm HT2", pointDataProcessorHT2, pointDataProcessor1gammaHT2, analysisSettingsNew.settingsHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT2.lowmassbgPointsHighPt, analysisResults.resultsHT2.invlist, analysisResults.resultsHT2.lowNormPointsPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	analysisResults.resultsSimMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimMB", analysisSettingsNew.title + ", low-mass norm MB MC", pointDataProcessorEmbedMB, pointDataProcessor1gammaMBSim, analysisSettingsNew.settingsSimMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimMB.lowmassbgPointsHighPt, analysisResults.resultsSimMB.invlist, analysisResults.resultsSimMB.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimHT1", analysisSettingsNew.title + ", low-mass norm HT1 MC", pointDataProcessorEmbedHT1, pointDataProcessor1gammaHT1Sim, analysisSettingsNew.settingsSimHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT1.lowmassbgPointsHighPt, analysisResults.resultsSimHT1.invlist, analysisResults.resultsSimHT1.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimHT2", analysisSettingsNew.title + ", low-mass norm HT2 MC", pointDataProcessorEmbedHT2, pointDataProcessor1gammaHT2Sim, analysisSettingsNew.settingsSimHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT2.lowmassbgPointsHighPt, analysisResults.resultsSimHT2.invlist, analysisResults.resultsSimHT2.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimEtaMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtaMB", analysisSettingsNew.title + ", low-mass norm MB eta MC", pointDataProcessorEmbedEtaMB, pointDataProcessor1gammaMBSimEta, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsHighPt, analysisResults.resultsSimEtaMB.invlist, analysisResults.resultsSimEtaMB.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimEtaHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtaHT1", analysisSettingsNew.title + ", low-mass norm HT1 eta MC", pointDataProcessorEmbedEtaHT1, pointDataProcessor1gammaHT1SimEta, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsHighPt, analysisResults.resultsSimEtaHT1.invlist, analysisResults.resultsSimEtaHT1.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimEtaHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtaHT2", analysisSettingsNew.title + ", low-mass norm HT2 eta MC", pointDataProcessorEmbedEtaHT2, pointDataProcessor1gammaHT2SimEta, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt, analysisResults.resultsSimEtaHT2.invlist, analysisResults.resultsSimEtaHT2.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimEtabgMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtabgMB", analysisSettingsNew.title + ", low-mass norm MB etabg MC", pointDataProcessorEmbedEtabgMB, pointDataProcessor1gammaMBSimEtabg, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsHighPt, analysisResults.resultsSimEtabgMB.invlist, analysisResults.resultsSimEtabgMB.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimEtabgHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtabgHT1", analysisSettingsNew.title + ", low-mass norm HT1 etabg MC", pointDataProcessorEmbedEtabgHT1, pointDataProcessor1gammaHT1SimEtabg, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsHighPt, analysisResults.resultsSimEtabgHT1.invlist, analysisResults.resultsSimEtabgHT1.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	analysisResults.resultsSimEtabgHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtabgHT2", analysisSettingsNew.title + ", low-mass norm HT2 etabg MC", pointDataProcessorEmbedEtabgHT2, pointDataProcessor1gammaHT2SimEtabg, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt, analysisResults.resultsSimEtabgHT2.invlist, analysisResults.resultsSimEtabgHT2.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);

	TF1 *shiftMBfunc = 0;
	TF1 *shiftHT1func = 0;
	TF1 *shiftHT2func = 0;
	if (analysisSettingsNew.correctForPtShift || analysisSettingsNew.correctForPtShiftEff) {
		TString shiftMBname = analysisSettingsNew.name + "_shiftMB";
		TString shiftMBfuncname = analysisSettingsNew.name + "_shiftMBfunc";
		TString shiftMBtitle = analysisSettingsNew.title + ", p_{T} shift MB";
		calculateBinPurity(candidateDataProcessorPtshiftMB
				, analysisResults.resultsSimMB.invPtshiftlist
				, shiftMBfuncname.Data(), analysisSettingsNew.shiftFuncStr
				, &shiftMBfunc
				, shiftMBname.Data(), shiftMBtitle.Data()
				, analysisSettingsNew.showSim
				);
		TString shiftHT1name = analysisSettingsNew.name + "_shiftHT1";
		TString shiftHT1funcname = analysisSettingsNew.name + "_shiftHT1func";
		TString shiftHT1title = analysisSettingsNew.title + ", p_{T} shift HT1";
		calculateBinPurity(candidateDataProcessorPtshiftHT1
				, analysisResults.resultsSimHT1.invPtshiftlist
				, shiftHT1funcname.Data(), analysisSettingsNew.shiftFuncStr
				, &shiftHT1func
				, shiftHT1name.Data(), shiftHT1title.Data()
				, analysisSettingsNew.showSim
				);
		TString shiftHT2name = analysisSettingsNew.name + "_shiftHT2";
		TString shiftHT2funcname = analysisSettingsNew.name + "_shiftHT2func";
		TString shiftHT2title = analysisSettingsNew.title + ", p_{T} shift HT2";
		calculateBinPurity(candidateDataProcessorPtshiftHT2
				, analysisResults.resultsSimHT2.invPtshiftlist
				, shiftHT2funcname.Data(), analysisSettingsNew.shiftFuncStr
				, &shiftHT2func
				, shiftHT2name.Data(), shiftHT2title.Data()
				, analysisSettingsNew.showSim
				);
		}
	calculatePtShiftCorrInv(analysisResults.resultsMB.invlist, analysisSettingsNew.correctForPtShift ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0, analysisResults.ptShiftCorrMB, analysisResults.ptCenterShiftMB);
	calculatePtShiftCorrInv(analysisResults.resultsHT1.invlist, analysisSettingsNew.correctForPtShift ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0, analysisResults.ptShiftCorrHT1, analysisResults.ptCenterShiftHT1);
	calculatePtShiftCorrInv(analysisResults.resultsHT2.invlist, analysisSettingsNew.correctForPtShift ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0, analysisResults.ptShiftCorrHT2, analysisResults.ptCenterShiftHT2);
	if (analysisSettingsNew.correctForPtShift || analysisSettingsNew.correctForPreciseBinPt) {
		correctPtShiftInv(analysisResults.resultsMB.invlist, analysisSettingsNew.correctForPtShift ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsMB.invBglist, analysisSettingsNew.correctForPtShift ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsMB.inv1gammalist, analysisSettingsNew.correctForPtShift ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsHT1.invlist, analysisSettingsNew.correctForPtShift ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsHT1.invBglist, analysisSettingsNew.correctForPtShift ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsHT1.inv1gammalist, analysisSettingsNew.correctForPtShift ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsHT2.invlist, analysisSettingsNew.correctForPtShift ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsHT2.invBglist, analysisSettingsNew.correctForPtShift ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsHT2.inv1gammalist, analysisSettingsNew.correctForPtShift ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
	}
	if (analysisSettingsNew.correctForPtShiftEff || analysisSettingsNew.correctForPreciseBinPt) {
		correctPtShiftInv(analysisResults.resultsSimMB.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimMB.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuMB, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimHT1.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimHT1.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuHT1, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimHT2.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimHT2.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuHT2, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtaMB.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtaMB.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuEtaMB, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtaHT1.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtaHT1.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuEtaHT1, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtaHT2.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtaHT2.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuEtaHT2, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtabgMB.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtabgMB.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuEtabgMB, analysisSettingsNew.correctForPtShiftEff ? shiftMBfunc : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtabgHT1.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtabgHT1.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuEtabgHT1, analysisSettingsNew.correctForPtShiftEff ? shiftHT1func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtabgHT2.invlist, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftInv(analysisResults.resultsSimEtabgHT2.inv1gammalist, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
		correctPtShiftBins(analysisResults.areaSimuEtabgHT2, analysisSettingsNew.correctForPtShiftEff ? shiftHT2func : 0, analysisSettingsNew.correctForPreciseBinPt ? ptBinPosFuncFromQCD : 0);
	}

	if ((analysisSettingsNew.settingsMB.fitLeftCpv.value < 0) && candidateDataProcessorMB) analysisSettingsNew.settingsMB.fitLeftCpv.value = candidateDataProcessorMB->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsMB.fitLeftCpv.value < 0) analysisSettingsNew.settingsMB.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsHT1.fitLeftCpv.value < 0) && candidateDataProcessorHT1) analysisSettingsNew.settingsHT1.fitLeftCpv.value = candidateDataProcessorHT1->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsHT1.fitLeftCpv.value < 0) analysisSettingsNew.settingsHT1.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsHT2.fitLeftCpv.value < 0) && candidateDataProcessorHT2) analysisSettingsNew.settingsHT2.fitLeftCpv.value = candidateDataProcessorHT2->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsHT2.fitLeftCpv.value < 0) analysisSettingsNew.settingsHT2.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsSimMB.fitLeftCpv.value < 0) && candidateDataProcessorEmbedMB) analysisSettingsNew.settingsSimMB.fitLeftCpv.value = candidateDataProcessorEmbedMB->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimMB.fitLeftCpv.value < 0) analysisSettingsNew.settingsSimMB.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsSimHT1.fitLeftCpv.value < 0) && candidateDataProcessorEmbedHT1) analysisSettingsNew.settingsSimHT1.fitLeftCpv.value = candidateDataProcessorEmbedHT1->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimHT1.fitLeftCpv.value < 0) analysisSettingsNew.settingsSimHT1.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsSimHT2.fitLeftCpv.value < 0) && candidateDataProcessorEmbedHT2) analysisSettingsNew.settingsSimHT2.fitLeftCpv.value = candidateDataProcessorEmbedHT2->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimHT2.fitLeftCpv.value < 0) analysisSettingsNew.settingsSimHT2.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsSimEtaMB.fitLeftCpv.value < 0) && candidateDataProcessorEmbedEtaMB) analysisSettingsNew.settingsSimEtaMB.fitLeftCpv.value = candidateDataProcessorEmbedEtaMB->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimEtaMB.fitLeftCpv.value < 0) analysisSettingsNew.settingsSimEtaMB.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsSimEtaHT1.fitLeftCpv.value < 0) && candidateDataProcessorEmbedEtaHT1) analysisSettingsNew.settingsSimEtaHT1.fitLeftCpv.value = candidateDataProcessorEmbedEtaHT1->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimEtaHT1.fitLeftCpv.value < 0) analysisSettingsNew.settingsSimEtaHT1.fitLeftCpv.value = 0;
	if ((analysisSettingsNew.settingsSimEtaHT2.fitLeftCpv.value < 0) && candidateDataProcessorEmbedEtaHT2) analysisSettingsNew.settingsSimEtaHT2.fitLeftCpv.value = candidateDataProcessorEmbedEtaHT2->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimEtaHT2.fitLeftCpv.value < 0) analysisSettingsNew.settingsSimEtaHT2.fitLeftCpv.value = 0;

	if ((analysisSettingsNew.settingsMB.peakRightCpv.value < 0) && candidateDataProcessorMB) analysisSettingsNew.settingsMB.peakRightCpv.value = candidateDataProcessorMB->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsMB.peakRightCpv.value < 0) analysisSettingsNew.settingsMB.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsHT1.peakRightCpv.value < 0) && candidateDataProcessorHT1) analysisSettingsNew.settingsHT1.peakRightCpv.value = candidateDataProcessorHT1->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsHT1.peakRightCpv.value < 0) analysisSettingsNew.settingsHT1.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsHT2.peakRightCpv.value < 0) && candidateDataProcessorHT2) analysisSettingsNew.settingsHT2.peakRightCpv.value = candidateDataProcessorHT2->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsHT2.peakRightCpv.value < 0) analysisSettingsNew.settingsHT2.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsSimMB.peakRightCpv.value < 0) && candidateDataProcessorEmbedMB) analysisSettingsNew.settingsSimMB.peakRightCpv.value = candidateDataProcessorEmbedMB->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimMB.peakRightCpv.value < 0) analysisSettingsNew.settingsSimMB.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsSimHT1.peakRightCpv.value < 0) && candidateDataProcessorEmbedHT1) analysisSettingsNew.settingsSimHT1.peakRightCpv.value = candidateDataProcessorEmbedHT1->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimHT1.peakRightCpv.value < 0) analysisSettingsNew.settingsSimHT1.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsSimHT2.peakRightCpv.value < 0) && candidateDataProcessorEmbedHT2) analysisSettingsNew.settingsSimHT2.peakRightCpv.value = candidateDataProcessorEmbedHT2->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimHT2.peakRightCpv.value < 0) analysisSettingsNew.settingsSimHT2.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsSimEtaMB.peakRightCpv.value < 0) && candidateDataProcessorEmbedEtaMB) analysisSettingsNew.settingsSimEtaMB.peakRightCpv.value = candidateDataProcessorEmbedEtaMB->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimEtaMB.peakRightCpv.value < 0) analysisSettingsNew.settingsSimEtaMB.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsSimEtaHT1.peakRightCpv.value < 0) && candidateDataProcessorEmbedEtaHT1) analysisSettingsNew.settingsSimEtaHT1.peakRightCpv.value = candidateDataProcessorEmbedEtaHT1->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimEtaHT1.peakRightCpv.value < 0) analysisSettingsNew.settingsSimEtaHT1.peakRightCpv.value = 0;
	if ((analysisSettingsNew.settingsSimEtaHT2.peakRightCpv.value < 0) && candidateDataProcessorEmbedEtaHT2) analysisSettingsNew.settingsSimEtaHT2.peakRightCpv.value = candidateDataProcessorEmbedEtaHT2->getCuts().getParametersPoint().trackDistCutLow;
	if (analysisSettingsNew.settingsSimEtaHT2.peakRightCpv.value < 0) analysisSettingsNew.settingsSimEtaHT2.peakRightCpv.value = 0;

	analysisResults.resultsMB.evNum = eventDataProcessorMB ? (*(eventDataProcessorMB->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsMB.evNumErr = eventDataProcessorMB ? (*(eventDataProcessorMB->binStatistics.begin())).getError() : 0;
	analysisResults.resultsMB.evNumSeenByMaker = eventDataProcessorMB ? (eventDataProcessorMB->numTotal) : 0;
	analysisResults.resultsMB.evNumPassedAllCuts = eventDataProcessorMB ? (eventDataProcessorMB->numPassedAllCuts) : 0;
	analysisResults.resultsHT1.evNum = eventDataProcessorHT1 ? (*(eventDataProcessorHT1->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsHT1.evNumErr = eventDataProcessorHT1 ? (*(eventDataProcessorHT1->binStatistics.begin())).getError() : 0;
	analysisResults.resultsHT1.evNumSeenByMaker = eventDataProcessorHT1 ? (eventDataProcessorHT1->numTotal) : 0;
	analysisResults.resultsHT1.evNumPassedAllCuts = eventDataProcessorHT1 ? (eventDataProcessorHT1->numPassedAllCuts) : 0;
	analysisResults.resultsHT2.evNum = eventDataProcessorHT2 ? (*(eventDataProcessorHT2->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsHT2.evNumErr = eventDataProcessorHT2 ? (*(eventDataProcessorHT2->binStatistics.begin())).getError() : 0;
	analysisResults.resultsHT2.evNumSeenByMaker = eventDataProcessorHT2 ? (eventDataProcessorHT2->numTotal) : 0;
	analysisResults.resultsHT2.evNumPassedAllCuts = eventDataProcessorHT2 ? (eventDataProcessorHT2->numPassedAllCuts) : 0;

	analysisResults.resultsSimMB.evNum = eventDataProcessorEmbedMB ? (*(eventDataProcessorEmbedMB->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimMB.evNumErr = eventDataProcessorEmbedMB ? (*(eventDataProcessorEmbedMB->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimMB.evNumSeenByMaker = eventDataProcessorEmbedMB ? (eventDataProcessorEmbedMB->numTotal) : 0;
	analysisResults.resultsSimMB.evNumPassedAllCuts = eventDataProcessorEmbedMB ? (eventDataProcessorEmbedMB->numPassedAllCuts) : 0;
	analysisResults.resultsSimHT1.evNum = eventDataProcessorEmbedHT1 ? (*(eventDataProcessorEmbedHT1->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimHT1.evNumErr = eventDataProcessorEmbedHT1 ? (*(eventDataProcessorEmbedHT1->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimHT1.evNumSeenByMaker = eventDataProcessorEmbedHT1 ? (eventDataProcessorEmbedHT1->numTotal) : 0;
	analysisResults.resultsSimHT1.evNumPassedAllCuts = eventDataProcessorEmbedHT1 ? (eventDataProcessorEmbedHT1->numPassedAllCuts) : 0;
	analysisResults.resultsSimHT2.evNum = eventDataProcessorEmbedHT2 ? (*(eventDataProcessorEmbedHT2->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimHT2.evNumErr = eventDataProcessorEmbedHT2 ? (*(eventDataProcessorEmbedHT2->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimHT2.evNumSeenByMaker = eventDataProcessorEmbedHT2 ? (eventDataProcessorEmbedHT2->numTotal) : 0;
	analysisResults.resultsSimHT2.evNumPassedAllCuts = eventDataProcessorEmbedHT2 ? (eventDataProcessorEmbedHT2->numPassedAllCuts) : 0;

	analysisResults.resultsSimEtaMB.evNum = eventDataProcessorEmbedEtaMB ? (*(eventDataProcessorEmbedEtaMB->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimEtaMB.evNumErr = eventDataProcessorEmbedEtaMB ? (*(eventDataProcessorEmbedEtaMB->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimEtaMB.evNumSeenByMaker = eventDataProcessorEmbedEtaMB ? (eventDataProcessorEmbedEtaMB->numTotal) : 0;
	analysisResults.resultsSimEtaMB.evNumPassedAllCuts = eventDataProcessorEmbedEtaMB ? (eventDataProcessorEmbedEtaMB->numPassedAllCuts) : 0;
	analysisResults.resultsSimEtaHT1.evNum = eventDataProcessorEmbedEtaHT1 ? (*(eventDataProcessorEmbedEtaHT1->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimEtaHT1.evNumErr = eventDataProcessorEmbedEtaHT1 ? (*(eventDataProcessorEmbedEtaHT1->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimEtaHT1.evNumSeenByMaker = eventDataProcessorEmbedEtaHT1 ? (eventDataProcessorEmbedEtaHT1->numTotal) : 0;
	analysisResults.resultsSimEtaHT1.evNumPassedAllCuts = eventDataProcessorEmbedEtaHT1 ? (eventDataProcessorEmbedEtaHT1->numPassedAllCuts) : 0;
	analysisResults.resultsSimEtaHT2.evNum = eventDataProcessorEmbedEtaHT2 ? (*(eventDataProcessorEmbedEtaHT2->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimEtaHT2.evNumErr = eventDataProcessorEmbedEtaHT2 ? (*(eventDataProcessorEmbedEtaHT2->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimEtaHT2.evNumSeenByMaker = eventDataProcessorEmbedEtaHT2 ? (eventDataProcessorEmbedEtaHT2->numTotal) : 0;
	analysisResults.resultsSimEtaHT2.evNumPassedAllCuts = eventDataProcessorEmbedEtaHT2 ? (eventDataProcessorEmbedEtaHT2->numPassedAllCuts) : 0;

	analysisResults.resultsSimEtabgMB.evNum = eventDataProcessorEmbedEtabgMB ? (*(eventDataProcessorEmbedEtabgMB->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimEtabgMB.evNumErr = eventDataProcessorEmbedEtabgMB ? (*(eventDataProcessorEmbedEtabgMB->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimEtabgMB.evNumSeenByMaker = eventDataProcessorEmbedEtabgMB ? (eventDataProcessorEmbedEtabgMB->numTotal) : 0;
	analysisResults.resultsSimEtabgMB.evNumPassedAllCuts = eventDataProcessorEmbedEtabgMB ? (eventDataProcessorEmbedEtabgMB->numPassedAllCuts) : 0;
	analysisResults.resultsSimEtabgHT1.evNum = eventDataProcessorEmbedEtabgHT1 ? (*(eventDataProcessorEmbedEtabgHT1->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimEtabgHT1.evNumErr = eventDataProcessorEmbedEtabgHT1 ? (*(eventDataProcessorEmbedEtabgHT1->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimEtabgHT1.evNumSeenByMaker = eventDataProcessorEmbedEtabgHT1 ? (eventDataProcessorEmbedEtabgHT1->numTotal) : 0;
	analysisResults.resultsSimEtabgHT1.evNumPassedAllCuts = eventDataProcessorEmbedEtabgHT1 ? (eventDataProcessorEmbedEtabgHT1->numPassedAllCuts) : 0;
	analysisResults.resultsSimEtabgHT2.evNum = eventDataProcessorEmbedEtabgHT2 ? (*(eventDataProcessorEmbedEtabgHT2->binStatistics.begin())).getValue() : 0;
	analysisResults.resultsSimEtabgHT2.evNumErr = eventDataProcessorEmbedEtabgHT2 ? (*(eventDataProcessorEmbedEtabgHT2->binStatistics.begin())).getError() : 0;
	analysisResults.resultsSimEtabgHT2.evNumSeenByMaker = eventDataProcessorEmbedEtabgHT2 ? (eventDataProcessorEmbedEtabgHT2->numTotal) : 0;
	analysisResults.resultsSimEtabgHT2.evNumPassedAllCuts = eventDataProcessorEmbedEtabgHT2 ? (eventDataProcessorEmbedEtabgHT2->numPassedAllCuts) : 0;

	analysisResults.resultsMB.deltaEta = candidateDataProcessorMB ? (candidateDataProcessorMB->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorMB->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsHT1.deltaEta = candidateDataProcessorHT1 ? (candidateDataProcessorHT1->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorHT1->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsHT2.deltaEta = candidateDataProcessorHT2 ? (candidateDataProcessorHT2->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorHT2->getCuts().getParametersCandidate().etaLow) : 1.0;

	analysisResults.resultsSimMB.deltaEta = candidateDataProcessorEmbedMB ? (candidateDataProcessorEmbedMB->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedMB->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsSimHT1.deltaEta = candidateDataProcessorEmbedHT1 ? (candidateDataProcessorEmbedHT1->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedHT1->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsSimHT2.deltaEta = candidateDataProcessorEmbedHT2 ? (candidateDataProcessorEmbedHT2->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedHT2->getCuts().getParametersCandidate().etaLow) : 1.0;

	analysisResults.resultsSimEtaMB.deltaEta = candidateDataProcessorEmbedEtaMB ? (candidateDataProcessorEmbedEtaMB->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedEtaMB->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsSimEtaHT1.deltaEta = candidateDataProcessorEmbedEtaHT1 ? (candidateDataProcessorEmbedEtaHT1->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedEtaHT1->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsSimEtaHT2.deltaEta = candidateDataProcessorEmbedEtaHT2 ? (candidateDataProcessorEmbedEtaHT2->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedEtaHT2->getCuts().getParametersCandidate().etaLow) : 1.0;

	analysisResults.resultsSimEtabgMB.deltaEta = candidateDataProcessorEmbedEtabgMB ? (candidateDataProcessorEmbedEtabgMB->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedEtabgMB->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsSimEtabgHT1.deltaEta = candidateDataProcessorEmbedEtabgHT1 ? (candidateDataProcessorEmbedEtabgHT1->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedEtabgHT1->getCuts().getParametersCandidate().etaLow) : 1.0;
	analysisResults.resultsSimEtabgHT2.deltaEta = candidateDataProcessorEmbedEtabgHT2 ? (candidateDataProcessorEmbedEtabgHT2->getCuts().getParametersCandidate().etaHigh - candidateDataProcessorEmbedEtabgHT2->getCuts().getParametersCandidate().etaLow) : 1.0;

	if (analysisSettingsNew.print) {
		cout << "Number of sim events seen by maker: " << analysisResults.resultsSimMB.evNumSeenByMaker << endl;
		cout << "Number of sim events passed all cuts: " << analysisResults.resultsSimMB.evNumPassedAllCuts << endl;
    		cout << "Number of sim MinBias events: " << analysisResults.resultsSimMB.evNum << endl;
		cout << "Number of sim HighTower-1 events: " << analysisResults.resultsSimHT1.evNum << endl;
		cout << "Number of sim HighTower-2 events: " << analysisResults.resultsSimHT2.evNum << endl;
		cout << "Number of sim eta events seen by maker: " << analysisResults.resultsSimEtaMB.evNumSeenByMaker << endl;
		cout << "Number of sim eta events passed all cuts: " << analysisResults.resultsSimEtaMB.evNumPassedAllCuts << endl;
    		cout << "Number of sim eta MinBias events: " << analysisResults.resultsSimEtaMB.evNum << endl;
		cout << "Number of sim eta HighTower-1 events: " << analysisResults.resultsSimEtaHT1.evNum << endl;
		cout << "Number of sim eta HighTower-2 events: " << analysisResults.resultsSimEtaHT2.evNum << endl;
		cout << "Number of sim etabg events seen by maker: " << analysisResults.resultsSimEtabgMB.evNumSeenByMaker << endl;
		cout << "Number of sim etabg events passed all cuts: " << analysisResults.resultsSimEtabgMB.evNumPassedAllCuts << endl;
    		cout << "Number of sim etabg MinBias events: " << analysisResults.resultsSimEtabgMB.evNum << endl;
		cout << "Number of sim etabg HighTower-1 events: " << analysisResults.resultsSimEtabgHT1.evNum << endl;
		cout << "Number of sim etabg HighTower-2 events: " << analysisResults.resultsSimEtabgHT2.evNum << endl;
		cout << "Number of real MinBias events seen by maker: " << analysisResults.resultsMB.evNumSeenByMaker << endl;
		cout << "Number of real MinBias events passed all cuts: " << analysisResults.resultsMB.evNumPassedAllCuts << endl;
		cout << "Number of real MinBias events in the analysis: " << analysisResults.resultsMB.evNum << endl;
		cout << "Number of real HighTower-1 events seen by maker: " << analysisResults.resultsHT1.evNumSeenByMaker << endl;
		cout << "Number of real HighTower-1 events passed all cuts: " << analysisResults.resultsHT1.evNumPassedAllCuts << endl;
		cout << "Number of real HighTower-1 events in the analysis: " << analysisResults.resultsHT1.evNum << endl;
		cout << "Number of real HighTower-2 events seen by maker: " << analysisResults.resultsHT2.evNumSeenByMaker << endl;
		cout << "Number of real HighTower-2 events passed all cuts: " << analysisResults.resultsHT2.evNumPassedAllCuts << endl;
		cout << "Number of real HighTower-2 events in the analysis: " << analysisResults.resultsHT2.evNum << endl;
	}

	calculateMeanAcceptance(eventDataProcessorMB, analysisResults.resultsMB.meanAcceptanceBTOW, analysisResults.resultsMB.meanAcceptanceBPRS, analysisResults.resultsMB.meanAcceptanceBSMDE, analysisResults.resultsMB.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorHT1, analysisResults.resultsHT1.meanAcceptanceBTOW, analysisResults.resultsHT1.meanAcceptanceBPRS, analysisResults.resultsHT1.meanAcceptanceBSMDE, analysisResults.resultsHT1.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorHT2, analysisResults.resultsHT2.meanAcceptanceBTOW, analysisResults.resultsHT2.meanAcceptanceBPRS, analysisResults.resultsHT2.meanAcceptanceBSMDE, analysisResults.resultsHT2.meanAcceptanceBSMDP);

	calculateMeanAcceptance(eventDataProcessorEmbedMB, analysisResults.resultsSimMB.meanAcceptanceBTOW, analysisResults.resultsSimMB.meanAcceptanceBPRS, analysisResults.resultsSimMB.meanAcceptanceBSMDE, analysisResults.resultsSimMB.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorEmbedHT1, analysisResults.resultsSimHT1.meanAcceptanceBTOW, analysisResults.resultsSimHT1.meanAcceptanceBPRS, analysisResults.resultsSimHT1.meanAcceptanceBSMDE, analysisResults.resultsSimHT1.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorEmbedHT2, analysisResults.resultsSimHT2.meanAcceptanceBTOW, analysisResults.resultsSimHT2.meanAcceptanceBPRS, analysisResults.resultsSimHT2.meanAcceptanceBSMDE, analysisResults.resultsSimHT2.meanAcceptanceBSMDP);

	calculateMeanAcceptance(eventDataProcessorEmbedEtaMB, analysisResults.resultsSimEtaMB.meanAcceptanceBTOW, analysisResults.resultsSimEtaMB.meanAcceptanceBPRS, analysisResults.resultsSimEtaMB.meanAcceptanceBSMDE, analysisResults.resultsSimEtaMB.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorEmbedEtaHT1, analysisResults.resultsSimEtaHT1.meanAcceptanceBTOW, analysisResults.resultsSimEtaHT1.meanAcceptanceBPRS, analysisResults.resultsSimEtaHT1.meanAcceptanceBSMDE, analysisResults.resultsSimEtaHT1.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorEmbedEtaHT2, analysisResults.resultsSimEtaHT2.meanAcceptanceBTOW, analysisResults.resultsSimEtaHT2.meanAcceptanceBPRS, analysisResults.resultsSimEtaHT2.meanAcceptanceBSMDE, analysisResults.resultsSimEtaHT2.meanAcceptanceBSMDP);

	calculateMeanAcceptance(eventDataProcessorEmbedEtabgMB, analysisResults.resultsSimEtabgMB.meanAcceptanceBTOW, analysisResults.resultsSimEtabgMB.meanAcceptanceBPRS, analysisResults.resultsSimEtabgMB.meanAcceptanceBSMDE, analysisResults.resultsSimEtabgMB.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorEmbedEtabgHT1, analysisResults.resultsSimEtabgHT1.meanAcceptanceBTOW, analysisResults.resultsSimEtabgHT1.meanAcceptanceBPRS, analysisResults.resultsSimEtabgHT1.meanAcceptanceBSMDE, analysisResults.resultsSimEtabgHT1.meanAcceptanceBSMDP);
	calculateMeanAcceptance(eventDataProcessorEmbedEtabgHT2, analysisResults.resultsSimEtabgHT2.meanAcceptanceBTOW, analysisResults.resultsSimEtabgHT2.meanAcceptanceBPRS, analysisResults.resultsSimEtabgHT2.meanAcceptanceBSMDE, analysisResults.resultsSimEtabgHT2.meanAcceptanceBSMDP);

	analysisResults.meanAcceptanceRatioMB = 1.0;
	analysisResults.meanAcceptanceRatioHT1 = 1.0;
	analysisResults.meanAcceptanceRatioHT2 = 1.0;
	analysisResults.meanAcceptanceRatioEtaMB = 1.0;
	analysisResults.meanAcceptanceRatioEtaHT1 = 1.0;
	analysisResults.meanAcceptanceRatioEtaHT2 = 1.0;
	if (analysisSettingsNew.correctMeanAcceptanceRatioMB && (analysisResults.resultsMB.meanAcceptanceBTOW != 0) && (analysisResults.resultsSimMB.meanAcceptanceBTOW != 0)) {
	    analysisResults.meanAcceptanceRatioMB = analysisResults.resultsSimMB.meanAcceptanceBTOW / analysisResults.resultsMB.meanAcceptanceBTOW;
	    cout << "Acceptance correction MB: " << analysisResults.meanAcceptanceRatioMB << ", simulated/real = " << analysisResults.resultsSimMB.meanAcceptanceBTOW << "/" << analysisResults.resultsMB.meanAcceptanceBTOW << endl;
	}
	if (analysisSettingsNew.correctMeanAcceptanceRatioHT1 && (analysisResults.resultsHT1.meanAcceptanceBTOW != 0) && (analysisResults.resultsSimHT1.meanAcceptanceBTOW != 0)) {
	    analysisResults.meanAcceptanceRatioHT1 = analysisResults.resultsSimHT1.meanAcceptanceBTOW / analysisResults.resultsHT1.meanAcceptanceBTOW;
	    cout << "Acceptance correction HT1: " << analysisResults.meanAcceptanceRatioHT1 << ", simulated/real = " << analysisResults.resultsSimHT1.meanAcceptanceBTOW << "/" << analysisResults.resultsHT1.meanAcceptanceBTOW << endl;
	}
	if (analysisSettingsNew.correctMeanAcceptanceRatioHT2 && (analysisResults.resultsHT2.meanAcceptanceBTOW != 0) && (analysisResults.resultsSimHT2.meanAcceptanceBTOW != 0)) {
	    analysisResults.meanAcceptanceRatioHT2 = analysisResults.resultsSimHT2.meanAcceptanceBTOW / analysisResults.resultsHT2.meanAcceptanceBTOW;
	    cout << "Acceptance correction HT2: " << analysisResults.meanAcceptanceRatioHT2 << ", simulated/real = " << analysisResults.resultsSimHT2.meanAcceptanceBTOW << "/" << analysisResults.resultsHT2.meanAcceptanceBTOW << endl;
	}
	if (analysisSettingsNew.correctMeanAcceptanceRatioEtaMB && (analysisResults.resultsMB.meanAcceptanceBTOW != 0) && (analysisResults.resultsSimEtaMB.meanAcceptanceBTOW != 0)) {
	    analysisResults.meanAcceptanceRatioEtaMB = analysisResults.resultsSimEtaMB.meanAcceptanceBTOW / analysisResults.resultsMB.meanAcceptanceBTOW;
	    cout << "Acceptance correction eta MB: " << analysisResults.meanAcceptanceRatioEtaMB << ", simulated/real = " << analysisResults.resultsSimEtaMB.meanAcceptanceBTOW << "/" << analysisResults.resultsMB.meanAcceptanceBTOW << endl;
	}
	if (analysisSettingsNew.correctMeanAcceptanceRatioEtaHT1 && (analysisResults.resultsHT1.meanAcceptanceBTOW != 0) && (analysisResults.resultsSimEtaHT1.meanAcceptanceBTOW != 0)) {
	    analysisResults.meanAcceptanceRatioEtaHT1 = analysisResults.resultsSimEtaHT1.meanAcceptanceBTOW / analysisResults.resultsHT1.meanAcceptanceBTOW;
	    cout << "Acceptance correction eta HT1: " << analysisResults.meanAcceptanceRatioEtaHT1 << ", simulated/real = " << analysisResults.resultsSimEtaHT1.meanAcceptanceBTOW << "/" << analysisResults.resultsHT1.meanAcceptanceBTOW << endl;
	}
	if (analysisSettingsNew.correctMeanAcceptanceRatioEtaHT2 && (analysisResults.resultsHT2.meanAcceptanceBTOW != 0) && (analysisResults.resultsSimEtaHT2.meanAcceptanceBTOW != 0)) {
	    analysisResults.meanAcceptanceRatioEtaHT2 = analysisResults.resultsSimEtaHT2.meanAcceptanceBTOW / analysisResults.resultsHT2.meanAcceptanceBTOW;
	    cout << "Acceptance correction eta HT2: " << analysisResults.meanAcceptanceRatioEtaHT2 << ", simulated/real = " << analysisResults.resultsSimEtaHT2.meanAcceptanceBTOW << "/" << analysisResults.resultsHT2.meanAcceptanceBTOW << endl;
	}

	TH1F *histTPCVertexIneffMB = 0;
	TH1F *histTPCVertexIneffHT1 = 0;
	TH1F *histTPCVertexIneffHT2 = 0;
	{
	    histTPCVertexIneffMB = new TH1F("zSimuSmearingPtMB", "TPC vertex finding inefficiency vs. #pi^{0} p_{T} (MB);p_{T};TPC vertex finding inefficiency", 500, 0, 50);
	    histTPCVertexIneffHT1 = new TH1F("zSimuSmearingPtHT1", "TPC vertex finding inefficiency vs. #pi^{0} p_{T} (HT1);p_{T};TPC vertex finding inefficiency", 500, 0, 50);
	    histTPCVertexIneffHT2 = new TH1F("zSimuSmearingPtHT2", "TPC vertex finding inefficiency vs. #pi^{0} p_{T} (HT2);p_{T};TPC vertex finding inefficiency", 500, 0, 50);
	    TH1F *ratioMB = 0;
	    TH1F *ratioHT1 = 0;
	    TH1F *ratioHT2 = 0;
	    TF1 *ratioMBfunc = 0;
	    TF1 *ratioHT1func = 0;
	    TF1 *ratioHT2func = 0;
	    const TH1F *histTPCYesMB = candidateDataProcessorMB ? candidateDataProcessorMB->getPt() : 0;
	    const TH1F *histTPCNoMB = candidateDataProcessorMB ? candidateDataProcessorMB->getTPCVertNo() : 0;
	    const TH1F *histTPCYesHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getPt() : 0;
	    const TH1F *histTPCNoHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getTPCVertNo() : 0;
	    const TH1F *histTPCYesHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getPt() : 0;
	    const TH1F *histTPCNoHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getTPCVertNo() : 0;
	    if (histTPCVertexIneffMB && histTPCYesMB && histTPCNoMB) {
		ratioMB = new TH1F(*histTPCNoMB);
		if (ratioMB) {
		    ratioMB->Divide(histTPCYesMB);
		    //setErrorOnZeroBins(ratioMB);
		    ratioMBfunc = new TF1("ratioMBfunc", analysisSettingsNew.tpcVertexEffFunc);
		    if (ratioMBfunc) {
			//ratioMBfunc->SetRange(histTPCVertexIneffMB->GetXaxis()->GetBinLowEdge(1), histTPCVertexIneffMB->GetXaxis()->GetBinUpEdge(histTPCVertexIneffMB->GetXaxis()->GetNbins()));
			ratioMBfunc->SetRange(analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.highPt);
			//ratioMBfunc->SetParameter(0, 0);
			//ratioMBfunc->SetParLimits(0, 0, 1);
			//ratioMBfunc->SetParameter(1, 1);
			//ratioMBfunc->SetParLimits(1, 0, 100);
			ratioMB->Fit(ratioMBfunc, "RQN");
			for (Int_t ibin = 1;ibin <= histTPCVertexIneffMB->GetXaxis()->GetNbins();ibin++) histTPCVertexIneffMB->SetBinContent(ibin, ratioMBfunc->Eval(histTPCVertexIneffMB->GetXaxis()->GetBinCenter(ibin)));
		    }
		}
	    }
	    if (histTPCVertexIneffHT1 && histTPCYesHT1 && histTPCNoHT1) {
		ratioHT1 = new TH1F(*histTPCNoHT1);
		if (ratioHT1) {
		    ratioHT1->Divide(histTPCYesHT1);
		    //setErrorOnZeroBins(ratioHT1);
		    ratioHT1func = new TF1("ratioHT1func", analysisSettingsNew.tpcVertexEffFunc);
		    if (ratioHT1func) {
			//ratioHT1func->SetRange(histTPCVertexIneffHT1->GetXaxis()->GetBinLowEdge(1), histTPCVertexIneffHT1->GetXaxis()->GetBinUpEdge(histTPCVertexIneffHT1->GetXaxis()->GetNbins()));
			ratioHT1func->SetRange(analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.highPt);
			//ratioHT1func->SetParameter(0, 0);
			//ratioHT1func->SetParLimits(0, 0, 1);
			//ratioHT1func->SetParameter(1, 1);
			//ratioHT1func->SetParLimits(1, 0, 100);
			ratioHT1->Fit(ratioHT1func, "RQN");
			for (Int_t ibin = 1;ibin <= histTPCVertexIneffHT1->GetXaxis()->GetNbins();ibin++) histTPCVertexIneffHT1->SetBinContent(ibin, ratioHT1func->Eval(histTPCVertexIneffHT1->GetXaxis()->GetBinCenter(ibin)));
		    }
		}
	    }
	    if (histTPCVertexIneffHT2 && histTPCYesHT2 && histTPCNoHT2) {
		ratioHT2 = new TH1F(*histTPCNoHT2);
		if (ratioHT2) {
		    ratioHT2->Divide(histTPCYesHT2);
		    //setErrorOnZeroBins(ratioHT2);
		    ratioHT2func = new TF1("ratioHT2func", analysisSettingsNew.tpcVertexEffFunc);
		    if (ratioHT2func) {
			//ratioHT2func->SetRange(histTPCVertexIneffHT2->GetXaxis()->GetBinLowEdge(1), histTPCVertexIneffHT2->GetXaxis()->GetBinUpEdge(histTPCVertexIneffHT2->GetXaxis()->GetNbins()));
			ratioHT2func->SetRange(analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.highPt);
			//ratioHT2func->SetParameter(0, 0);
			//ratioHT2func->SetParLimits(0, 0, 1);
			//ratioHT2func->SetParameter(1, 1);
			//ratioHT2func->SetParLimits(1, 0, 100);
			ratioHT2->Fit(ratioHT2func, "RQN");
			for (Int_t ibin = 1;ibin <= histTPCVertexIneffHT2->GetXaxis()->GetNbins();ibin++) histTPCVertexIneffHT2->SetBinContent(ibin, ratioHT2func->Eval(histTPCVertexIneffHT2->GetXaxis()->GetBinCenter(ibin)));
		    }
		}
	    }
	    if (analysisSettingsNew.showTPCVertexEff) {
		TString name = analysisSettingsNew.name + "_TPCVert";
		TString title = analysisSettingsNew.title + ", TPC vertex finding inefficiency";
		TCanvas *c = new TCanvas(name, title);
		if (c) c->cd();
		if (ratioMB) {
		    ratioMB->SetLineColor(analysisSettingsNew.settingsMB.color);
		    ratioMB->SetMarkerColor(analysisSettingsNew.settingsMB.color);
		    ratioMB->DrawCopy("");
		    ratioMB->SetMaximum(1);
		    ratioMB->SetMinimum(0);
		}
		if (ratioMBfunc) {
		    ratioMBfunc->SetLineColor(analysisSettingsNew.settingsMB.color);
		    ratioMBfunc->SetMarkerColor(analysisSettingsNew.settingsMB.color);
		    ratioMBfunc->DrawCopy("SAME");
		}
		if (ratioHT1) {
		    ratioHT1->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1->SetMarkerColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1->DrawCopy("SAME");
		}
		if (ratioHT1func) {
		    ratioHT1func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1func->SetMarkerColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1func->DrawCopy("SAME");
		}
		if (ratioHT2) {
		    ratioHT2->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2->SetMarkerColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2->DrawCopy("SAME");
		}
		if (ratioHT2func) {
		    ratioHT2func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2func->SetMarkerColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2func->DrawCopy("SAME");
		}
	    }
	    if (ratioMB) delete ratioMB;
	    if (ratioHT1) delete ratioHT1;
	    if (ratioHT2) delete ratioHT2;
	    if (ratioMBfunc) delete ratioMBfunc;
	    if (ratioHT1func) delete ratioHT1func;
	    if (ratioHT2func) delete ratioHT2func;
	}

	TH1F *histBBCVertexIneffMB = 0;
	TH1F *histBBCVertexIneffHT1 = 0;
	TH1F *histBBCVertexIneffHT2 = 0;
	{
	    histBBCVertexIneffMB = new TH1F("zSimuSmearingZeroPtMB", "BBC vertex finding inefficiency vs. #pi^{0} p_{T} (MB);p_{T};BBC vertex finding inefficiency", 500, 0, 50);
	    histBBCVertexIneffHT1 = new TH1F("zSimuSmearingZeroPtHT1", "BBC vertex finding inefficiency vs. #pi^{0} p_{T} (HT1);p_{T};BBC vertex finding inefficiency", 500, 0, 50);
	    histBBCVertexIneffHT2 = new TH1F("zSimuSmearingZeroPtHT2", "BBC vertex finding inefficiency vs. #pi^{0} p_{T} (HT2);p_{T};BBC vertex finding inefficiency", 500, 0, 50);
	    TH1F *ratioMB = 0;
	    TH1F *ratioHT1 = 0;
	    TH1F *ratioHT2 = 0;
	    TF1 *ratioMBfunc = 0;
	    TF1 *ratioHT1func = 0;
	    TF1 *ratioHT2func = 0;
	    const TH1F *histBBCYesMB = candidateDataProcessorMB ? candidateDataProcessorMB->getTPCVertNo() : 0;
	    const TH1F *histBBCNoMB = candidateDataProcessorMB ? candidateDataProcessorMB->getTPCVertNoBBCVertNo() : 0;
	    const TH1F *histBBCYesHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getTPCVertNo() : 0;
	    const TH1F *histBBCNoHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getTPCVertNoBBCVertNo() : 0;
	    const TH1F *histBBCYesHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getTPCVertNo() : 0;
	    const TH1F *histBBCNoHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getTPCVertNoBBCVertNo() : 0;
	    if (histBBCVertexIneffMB && histBBCYesMB && histBBCNoMB) {
		ratioMB = new TH1F(*histBBCNoMB);
		if (ratioMB) {
		    ratioMB->Divide(histBBCYesMB);
		    //setErrorOnZeroBins(ratioMB);
		    ratioMBfunc = new TF1("ratioMBfunc", analysisSettingsNew.bbcVertexEffFunc);
		    if (ratioMBfunc) {
			//ratioMBfunc->SetRange(histBBCVertexIneffMB->GetXaxis()->GetBinLowEdge(1), histBBCVertexIneffMB->GetXaxis()->GetBinUpEdge(histBBCVertexIneffMB->GetXaxis()->GetNbins()));
			ratioMBfunc->SetRange(analysisSettingsNew.settingsMB.lowPt, analysisSettingsNew.settingsMB.highPt);
			//ratioMBfunc->SetParameter(0, 0);
			//ratioMBfunc->SetParLimits(0, 0, 1);
			//ratioMBfunc->SetParameter(1, 1);
			//ratioMBfunc->SetParLimits(1, 0, 100);
			ratioMB->Fit(ratioMBfunc, "RQN");
			for (Int_t ibin = 1;ibin <= histBBCVertexIneffMB->GetXaxis()->GetNbins();ibin++) histBBCVertexIneffMB->SetBinContent(ibin, ratioMBfunc->Eval(histBBCVertexIneffMB->GetXaxis()->GetBinCenter(ibin)));
		    }
		}
	    }
	    if (histBBCVertexIneffHT1 && histBBCYesHT1 && histBBCNoHT1) {
		ratioHT1 = new TH1F(*histBBCNoHT1);
		if (ratioHT1) {
		    ratioHT1->Divide(histBBCYesHT1);
		    //setErrorOnZeroBins(ratioHT1);
		    ratioHT1func = new TF1("ratioHT1func", analysisSettingsNew.bbcVertexEffFunc);
		    if (ratioHT1func) {
			//ratioHT1func->SetRange(histBBCVertexIneffHT1->GetXaxis()->GetBinLowEdge(1), histBBCVertexIneffHT1->GetXaxis()->GetBinUpEdge(histBBCVertexIneffHT1->GetXaxis()->GetNbins()));
			ratioHT1func->SetRange(analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT1.highPt);
			//ratioHT1func->SetParameter(0, 0);
			//ratioHT1func->SetParLimits(0, 0, 1);
			//ratioHT1func->SetParameter(1, 1);
			//ratioHT1func->SetParLimits(1, 0, 100);
			ratioHT1->Fit(ratioHT1func, "RQN");
			for (Int_t ibin = 1;ibin <= histBBCVertexIneffHT1->GetXaxis()->GetNbins();ibin++) histBBCVertexIneffHT1->SetBinContent(ibin, ratioHT1func->Eval(histBBCVertexIneffHT1->GetXaxis()->GetBinCenter(ibin)));
		    }
		}
	    }
	    if (histBBCVertexIneffHT2 && histBBCYesHT2 && histBBCNoHT2) {
		ratioHT2 = new TH1F(*histBBCNoHT2);
		if (ratioHT2) {
		    ratioHT2->Divide(histBBCYesHT2);
		    //setErrorOnZeroBins(ratioHT2);
		    ratioHT2func = new TF1("ratioHT2func", analysisSettingsNew.bbcVertexEffFunc);
		    if (ratioHT2func) {
			//ratioHT2func->SetRange(histBBCVertexIneffHT2->GetXaxis()->GetBinLowEdge(1), histBBCVertexIneffHT2->GetXaxis()->GetBinUpEdge(histBBCVertexIneffHT2->GetXaxis()->GetNbins()));
			ratioHT2func->SetRange(analysisSettingsNew.settingsHT2.lowPt, analysisSettingsNew.settingsHT2.highPt);
			//ratioHT2func->SetParameter(0, 0);
			//ratioHT2func->SetParLimits(0, 0, 1);
			//ratioHT2func->SetParameter(1, 1);
			//ratioHT2func->SetParLimits(1, 0, 100);
			ratioHT2->Fit(ratioHT2func, "RQN");
			for (Int_t ibin = 1;ibin <= histBBCVertexIneffHT2->GetXaxis()->GetNbins();ibin++) histBBCVertexIneffHT2->SetBinContent(ibin, ratioHT2func->Eval(histBBCVertexIneffHT2->GetXaxis()->GetBinCenter(ibin)));
		    }
		}
	    }
	    if (analysisSettingsNew.showBBCVertexEff) {
		TString name = analysisSettingsNew.name + "_BBCVert";
		TString title = analysisSettingsNew.title + ", BBC vertex finding inefficiency";
		TCanvas *c = new TCanvas(name, title);
		if (c) c->cd();
		if (ratioMB) {
		    ratioMB->SetLineColor(analysisSettingsNew.settingsMB.color);
		    ratioMB->SetMarkerColor(analysisSettingsNew.settingsMB.color);
		    ratioMB->DrawCopy("");
		    ratioMB->SetMaximum(1);
		    ratioMB->SetMinimum(0);
		}
		if (ratioMBfunc) {
		    ratioMBfunc->SetLineColor(analysisSettingsNew.settingsMB.color);
		    ratioMBfunc->SetMarkerColor(analysisSettingsNew.settingsMB.color);
		    ratioMBfunc->DrawCopy("SAME");
		}
		if (ratioHT1) {
		    ratioHT1->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1->SetMarkerColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1->DrawCopy("SAME");
		}
		if (ratioHT1func) {
		    ratioHT1func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1func->SetMarkerColor(analysisSettingsNew.settingsHT1.color);
		    ratioHT1func->DrawCopy("SAME");
		}
		if (ratioHT2) {
		    ratioHT2->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2->SetMarkerColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2->DrawCopy("SAME");
		}
		if (ratioHT2func) {
		    ratioHT2func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2func->SetMarkerColor(analysisSettingsNew.settingsHT2.color);
		    ratioHT2func->DrawCopy("SAME");
		}
	    }
	    if (ratioMB) delete ratioMB;
	    if (ratioHT1) delete ratioHT1;
	    if (ratioHT2) delete ratioHT2;
	    if (ratioMBfunc) delete ratioMBfunc;
	    if (ratioHT1func) delete ratioHT1func;
	    if (ratioHT2func) delete ratioHT2func;
	}

	{
    	    const TH2F *histRefMultEastDay = eventDataProcessorMB ? eventDataProcessorMB->getNumFtpcPrimariesEastDay() : 0;
	    TAxis *xaxRefMultEastDay = histRefMultEastDay ? const_cast<TAxis*>(histRefMultEastDay->GetXaxis()) : 0;
	    const TH2F *histRefMultWestDay = eventDataProcessorMB ? eventDataProcessorMB->getNumFtpcPrimariesWestDay() : 0;
	    TAxis *xaxRefMultWestDay = histRefMultWestDay ? const_cast<TAxis*>(histRefMultWestDay->GetXaxis()) : 0;
	    if (histRefMultEastDay && xaxRefMultEastDay) {
		for (Int_t iday = 1;iday <= xaxRefMultEastDay->GetNbins();iday++) {
			TH1 *projectionY = (TH1*)histRefMultEastDay->ProjectionY("_py", iday, iday, "");
			Float_t mean = projectionY ? projectionY->GetMean() : 0;
			if (mean != 0) {
				TBinParameters par; par.min = xaxRefMultEastDay->GetBinCenter(iday); par.max = par.min;
				TBinStatistics bin; bin.setParameters(par);
				bin.setValue(mean); bin.setError(projectionY->GetRMS() / sqrt(projectionY->GetSumOfWeights()));
				analysisResults.refMultEast.push_back(bin);
			}
			if (projectionY) delete projectionY;
		}
	    }
	    if (histRefMultWestDay && xaxRefMultWestDay) {
		for (Int_t iday = 1;iday <= xaxRefMultWestDay->GetNbins();iday++) {
			TH1 *projectionY = (TH1*)histRefMultWestDay->ProjectionY("_py", iday, iday, "");
			Float_t mean = projectionY ? projectionY->GetMean() : 0;
			if (mean != 0) {
				TBinParameters par; par.min = xaxRefMultWestDay->GetBinCenter(iday); par.max = par.min;
				TBinStatistics bin; bin.setParameters(par);
				bin.setValue(mean); bin.setError(projectionY->GetRMS() / sqrt(projectionY->GetSumOfWeights()));
				analysisResults.refMultWest.push_back(bin);
			}
			if (projectionY) delete projectionY;
		}
	    }
	    if (false) {
		const TH1F *histRefMult = eventDataProcessorMB ? eventDataProcessorMB->getNumFtpcPrimaries() : 0;
		TAxis *xaxRefMult = histRefMult ? const_cast<TAxis*>(histRefMult->GetXaxis()) : 0;
		if (histRefMult && xaxRefMult) {
			Int_t mostPeripheralLow = 0;
			Int_t mostPeripheralHigh = 0;
			Int_t mostCentralLow = 0;
			Int_t mostCentralHigh = 1000;
			Float_t mostCentralPercentTrue = 0.2;
			Float_t midCentralPercentTrue = 0.2;
			Float_t mostPeripheralPercentTrue = 1.0 - mostCentralPercentTrue - midCentralPercentTrue;
			Float_t mostCentralPercentTrueApprox = 0;
			Float_t midCentralPercentTrueApprox = 0;
			Float_t mostPeripheralPercentTrueApprox = 0;
			do {
				mostPeripheralHigh++;
				mostPeripheralPercentTrueApprox = histRefMult->Integral(xaxRefMult->FindBin(mostPeripheralLow), xaxRefMult->FindBin(mostPeripheralHigh) - 1) / histRefMult->Integral(1, xaxRefMult->GetNbins());
			} while (mostPeripheralPercentTrueApprox < mostPeripheralPercentTrue);
			mostCentralLow = mostPeripheralHigh;
			do {
				mostCentralLow++;
				mostCentralPercentTrueApprox = histRefMult->Integral(xaxRefMult->FindBin(mostCentralLow), xaxRefMult->GetNbins()) / histRefMult->Integral(1, xaxRefMult->GetNbins());
			} while (mostCentralPercentTrueApprox > mostCentralPercentTrue);
			midCentralPercentTrueApprox = histRefMult->Integral(xaxRefMult->FindBin(mostPeripheralHigh), xaxRefMult->FindBin(mostCentralLow) - 1) / histRefMult->Integral(1, xaxRefMult->GetNbins());
			cout << "Official 20-20-60 centrality bins:" << endl;
			cout << "Most central 20%: " << mostCentralLow << " <= RefMult < " << mostCentralHigh << " (" << (mostCentralPercentTrueApprox*100) << "%)" << endl;
			cout << "Mid central 20%: " << mostPeripheralHigh << " <= RefMult < " << mostCentralLow << " (" << (midCentralPercentTrueApprox*100) << "%)" << endl;
			cout << "Most peripheral 60%: " << mostPeripheralLow << " <= RefMult < " << mostPeripheralHigh << " (" << (mostPeripheralPercentTrueApprox*100) << "%)" << endl;
			cout << "Check: sum is " << ((mostCentralPercentTrueApprox + midCentralPercentTrueApprox + mostPeripheralPercentTrueApprox) * 100) << "%" << endl;
		}
	    }
	}

	TH1F *histJetFracMB = 0;
	TH1F *histJetFracHT1 = 0;
	TH1F *histJetFracHT2 = 0;
	{
		const TH1F *histPtMB = candidateDataProcessorMB ? candidateDataProcessorMB->getPt() : 0;
		const TH1F *histJetYesMB = candidateDataProcessorMB ? candidateDataProcessorMB->getJetYes() : 0;
		if (histPtMB && histJetYesMB) {
		    histJetFracMB = new TH1F(*histJetYesMB);
		    if (histJetFracMB) {
			histJetFracMB->Divide(histPtMB);
		    }
		}
	}
	{
		const TH1F *histPtHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getPt() : 0;
		const TH1F *histJetYesHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getJetYes() : 0;
		if (histPtHT1 && histJetYesHT1) {
		    histJetFracHT1 = new TH1F(*histJetYesHT1);
		    if (histJetFracHT1) {
			histJetFracHT1->Divide(histPtHT1);
		    }
		}
	}
	{
		const TH1F *histPtHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getPt() : 0;
		const TH1F *histJetYesHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getJetYes() : 0;
		if (histPtHT2 && histJetYesHT2) {
		    histJetFracHT2 = new TH1F(*histJetYesHT2);
		    if (histJetFracHT2) {
			histJetFracHT2->Divide(histPtHT2);
		    }
		}
	}
	if (analysisSettingsNew.settingsMB.fitJetBgFractionFromCandidates && histJetFracMB) {
	    analysisSettingsNew.settingsMB.jetBgFraction.useFit = true;
	    fitParameter(0, histJetFracMB, analysisSettingsNew.settingsMB.jetBgFraction, analysisResults.resultsMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandMB");
	    if (analysisSettingsNew.settingsMB.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsMB.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsMB.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsMB.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsMB.jetBgFraction.wasFit = false;
	}
	if (analysisSettingsNew.settingsHT1.fitJetBgFractionFromCandidates && histJetFracHT1) {
	    analysisSettingsNew.settingsHT1.jetBgFraction.useFit = true;
	    fitParameter(0, histJetFracHT1, analysisSettingsNew.settingsHT1.jetBgFraction, analysisResults.resultsHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandHT1");
	    if (analysisSettingsNew.settingsHT1.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsHT1.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsHT1.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsHT1.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsHT1.jetBgFraction.wasFit = false;
	}
	if (analysisSettingsNew.settingsHT2.fitJetBgFractionFromCandidates && histJetFracHT2) {
	    analysisSettingsNew.settingsHT2.jetBgFraction.useFit = true;
	    fitParameter(0, histJetFracHT2, analysisSettingsNew.settingsHT2.jetBgFraction, analysisResults.resultsHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandHT2");
	    if (analysisSettingsNew.settingsHT2.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsHT2.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsHT2.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsHT2.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsHT2.jetBgFraction.wasFit = false;
	}

	TWeightCalculator wgammas("weight_1gamma_pi0", "Single photon weight");
	Bool_t wgammasGood = false;
	{
	    const TH1F *histGammasPt = simuDataProcessorMB ? simuDataProcessorMB->gammas.getPt() : 0;
	    if (histGammasPt) {
		TH1F gammas(*histGammasPt);
    		gammas.GetXaxis()->SetRangeUser(analysisSettingsNew.settingsSimMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT2.lowmassbgPointsHighPt);
		wgammas.Fit(&gammas);
		//wgammas.mult = wgammas.multDrift;
		//wgammas.multDrift = 0.0;
		wgammas.rangeLow = 0;
		wgammas.rangeHigh = 100;
		wgammasGood = true;
	    }
	}

	TWeightCalculator wgammasEta("weight_1gamma_eta", "Single photon weight - eta");
	Bool_t wgammasEtaGood = false;
	{
	    const TH1F *histGammasEtaPt = simuDataProcessorEtaMB ? simuDataProcessorEtaMB->gammas.getPt() : 0;
	    if (histGammasEtaPt) {
		TH1F gammas(*histGammasEtaPt);
		gammas.GetXaxis()->SetRangeUser(analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt);
		wgammasEta.Fit(&gammas);
		//wgammasEta.mult = wgammasEta.multDrift;
		//wgammasEta.multDrift = 0.0;
		wgammasEta.rangeLow = 0;
		wgammasEta.rangeHigh = 100;
		wgammasEtaGood = true;
	    }
	}

	TWeightCalculator wpointsMB("weight_1gamma_data_MB", "Single photon weight - data MB");
	Bool_t wpointsMBGood = false;
	{
	    const TH1F *histPointPtMB = pointDataProcessorMB ? pointDataProcessorMB->getPt() : 0;
	    if (histPointPtMB) {
		TH1F points(*histPointPtMB);
		const TH1F *histEnergySimMBc = gammaDataProcessorMB ? gammaDataProcessorMB->getPt() : 0;
		const TH1F *histEnergyRecoMBc = pointDataProcessor1gammaMB ? pointDataProcessor1gammaMB->getPt() : 0;
		if (histEnergySimMBc && histEnergyRecoMBc) {
		    points.Multiply(histEnergySimMBc);
		    points.Divide(histEnergyRecoMBc);
		}
		points.GetXaxis()->SetRangeUser(analysisSettingsNew.settingsMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsMB.lowmassbgPointsHighPt);
		wpointsMB.Fit(&points);
		//wpointsMB.mult = wpointsMB.multDrift;
		//wpointsMB.multDrift = 0.0;
		wpointsMB.rangeLow = 0;
		wpointsMB.rangeHigh = 100;
		wpointsMBGood = true;
	    }
	}

	TWeightCalculator wpointsHT1("weight_1gamma_data_HT1", "Single photon weight - data HT1");
	Bool_t wpointsHT1Good = false;
	{
	    const TH1F *histPointPtHT1 = pointDataProcessorHT1 ? pointDataProcessorHT1->getPt() : 0;
	    if (histPointPtHT1) {
		TH1F points(*histPointPtHT1);
		const TH1F *histEnergySimHT1c = gammaDataProcessorHT1 ? gammaDataProcessorHT1->getPt() : 0;
		const TH1F *histEnergyRecoHT1c = pointDataProcessor1gammaHT1 ? pointDataProcessor1gammaHT1->getPt() : 0;
		if (histEnergySimHT1c && histEnergyRecoHT1c) {
		    points.Multiply(histEnergySimHT1c);
		    points.Divide(histEnergyRecoHT1c);
		}
		points.GetXaxis()->SetRangeUser(analysisSettingsNew.settingsHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT1.lowmassbgPointsHighPt);
		wpointsHT1.Fit(&points);
		//wpointsHT1.mult = wpointsHT1.multDrift;
		//wpointsHT1.multDrift = 0.0;
		wpointsHT1.rangeLow = 0;
		wpointsHT1.rangeHigh = 100;
		wpointsHT1Good = true;
	    }
	}

	TWeightCalculator wpointsHT2("weight_1gamma_data_HT2", "Single photon weight - data HT2");
	Bool_t wpointsHT2Good = false;
	{
	    const TH1F *histPointPtHT2 = pointDataProcessorHT2 ? pointDataProcessorHT2->getPt() : 0;
	    if (histPointPtHT2) {
		TH1F points(*histPointPtHT2);
		const TH1F *histEnergySimHT2c = gammaDataProcessorHT2 ? gammaDataProcessorHT2->getPt() : 0;
		const TH1F *histEnergyRecoHT2c = pointDataProcessor1gammaHT2 ? pointDataProcessor1gammaHT2->getPt() : 0;
		if (histEnergySimHT2c && histEnergyRecoHT2c) {
		    points.Multiply(histEnergySimHT2c);
		    points.Divide(histEnergyRecoHT2c);
		}
		points.GetXaxis()->SetRangeUser(analysisSettingsNew.settingsHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT2.lowmassbgPointsHighPt);
		wpointsHT2.Fit(&points);
		//wpointsHT2.mult = wpointsHT2.multDrift;
		//wpointsHT2.multDrift = 0.0;
		wpointsHT2.rangeLow = 0;
		wpointsHT2.rangeHigh = 100;
		wpointsHT2Good = true;
	    }
	}

	TWeightCalculator weightUsed("weightUsed", "Used weight");
	if (candidateDataProcessorEmbedMB) weightUsed = candidateDataProcessorEmbedMB->getWeightCalculator();

	analysisResults.totSimuYield = simuDataProcessorMB ? simuDataProcessorMB->binStatistics.back().getValue() : 0;
	analysisResults.totSimuYieldEta = simuDataProcessorEtaMB ? simuDataProcessorEtaMB->binStatistics.back().getValue() : 0;
	analysisResults.totSimuYieldEtabg = simuDataProcessorEtabgMB ? simuDataProcessorEtabgMB->binStatistics.back().getValue() : 0;

	analysisResults.resultsMB.scalePeakShapeEtabgRel = ((analysisResults.resultsSimEtabgMB.evNumPassedAllCuts > 0) && (analysisResults.resultsSimEtaMB.evNumPassedAllCuts > 0)) ? ((Float_t)analysisResults.resultsSimEtaMB.evNumPassedAllCuts / (Float_t)analysisResults.resultsSimEtabgMB.evNumPassedAllCuts) : 1.0;
	analysisResults.resultsHT1.scalePeakShapeEtabgRel = analysisResults.resultsMB.scalePeakShapeEtabgRel;
	analysisResults.resultsHT2.scalePeakShapeEtabgRel = analysisResults.resultsMB.scalePeakShapeEtabgRel;

	if (deletePools) {
	    removePool(poolEmbedMB, true, true);
	    removePool(poolEmbedHT1, true, true);
	    removePool(poolEmbedHT2, true, true);
	    removePool(pool1gammaSimMB, true, true);
	    removePool(pool1gammaSimHT1, true, true);
	    removePool(pool1gammaSimHT2, true, true);

	    removePool(poolEmbedEtaMB, true, true);
	    removePool(poolEmbedEtaHT1, true, true);
	    removePool(poolEmbedEtaHT2, true, true);
	    removePool(poolEmbedEtabgMB, true, true);
	    removePool(poolEmbedEtabgHT1, true, true);
	    removePool(poolEmbedEtabgHT2, true, true);
	    removePool(pool1gammaSimEtaMB, true, true);
	    removePool(pool1gammaSimEtaHT1, true, true);
	    removePool(pool1gammaSimEtaHT2, true, true);

	    removePool(poolMB, true, true);
	    removePool(poolHT1, true, true);
	    removePool(poolHT2, true, true);
	    removePool(pool1gammaMB, true, true);
	    removePool(pool1gammaHT1, true, true);
	    removePool(pool1gammaHT2, true, true);

	    removePool(poolPtshiftMB, true, true);
	    removePool(poolPtshiftHT1, true, true);
	    removePool(poolPtshiftHT2, true, true);
	}

	TCanvas *peakFinderSimMBcanvas = fitSpectraBins(analysisSettingsNew.settingsSimMB, analysisResults.resultsSimMB, 0, 0, 0);
	TCanvas *peakFinderSimHT1canvas = fitSpectraBins(analysisSettingsNew.settingsSimHT1, analysisResults.resultsSimHT1, 0, 0, 0);
	TCanvas *peakFinderSimHT2canvas = fitSpectraBins(analysisSettingsNew.settingsSimHT2, analysisResults.resultsSimHT2, 0, 0, 0);

	TCanvas *peakFinderSimEtaMBcanvas = fitSpectraBins(analysisSettingsNew.settingsSimEtaMB, analysisResults.resultsSimEtaMB, 0, 0, 0);
	TCanvas *peakFinderSimEtaHT1canvas = fitSpectraBins(analysisSettingsNew.settingsSimEtaHT1, analysisResults.resultsSimEtaHT1, 0, 0, 0);
	TCanvas *peakFinderSimEtaHT2canvas = fitSpectraBins(analysisSettingsNew.settingsSimEtaHT2, analysisResults.resultsSimEtaHT2, 0, 0, 0);

	TAnalysisSettingsTrigger settingsSimEtabgMB = analysisSettingsNew.settingsSimEtaMB;
	TAnalysisSettingsTrigger settingsSimEtabgHT1 = analysisSettingsNew.settingsSimEtaHT1;
	TAnalysisSettingsTrigger settingsSimEtabgHT2 = analysisSettingsNew.settingsSimEtaHT2;
	settingsSimEtabgMB.name += "_etabg";
	settingsSimEtabgHT1.name += "_etabg";
	settingsSimEtabgHT2.name += "_etabg";
	settingsSimEtabgMB.title += " (etabg)";
	settingsSimEtabgHT1.title += " (etabg)";
	settingsSimEtabgHT2.title += " (etabg)";
	TCanvas *peakFinderSimEtabgMBcanvas = fitSpectraBins(settingsSimEtabgMB, analysisResults.resultsSimEtabgMB, 0, 0, 0);
	TCanvas *peakFinderSimEtabgHT1canvas = fitSpectraBins(settingsSimEtabgHT1, analysisResults.resultsSimEtabgHT1, 0, 0, 0);
	TCanvas *peakFinderSimEtabgHT2canvas = fitSpectraBins(settingsSimEtabgHT2, analysisResults.resultsSimEtabgHT2, 0, 0, 0);

	analysisResults.resultsMB.scalePeakShapeLowMass = analysisResults.resultsSimMB.scaleBgLowMass;
	analysisResults.resultsMB.scalePeakShapeEtaLowMass = analysisResults.resultsSimEtaMB.scaleBgLowMass;
	analysisResults.resultsMB.scalePeakShapeEtabgLowMass = analysisResults.resultsSimEtabgMB.scaleBgLowMass;
	analysisResults.resultsHT1.scalePeakShapeLowMass = analysisResults.resultsSimHT1.scaleBgLowMass;
	analysisResults.resultsHT1.scalePeakShapeEtaLowMass = analysisResults.resultsSimEtaHT1.scaleBgLowMass;
	analysisResults.resultsHT1.scalePeakShapeEtabgLowMass = analysisResults.resultsSimEtabgHT1.scaleBgLowMass;
	analysisResults.resultsHT2.scalePeakShapeLowMass = analysisResults.resultsSimHT2.scaleBgLowMass;
	analysisResults.resultsHT2.scalePeakShapeEtaLowMass = analysisResults.resultsSimEtaHT2.scaleBgLowMass;
	analysisResults.resultsHT2.scalePeakShapeEtabgLowMass = analysisResults.resultsSimEtabgHT2.scaleBgLowMass;

	TCanvas *peakFinderMBcanvas = fitSpectraBins(analysisSettingsNew.settingsMB, analysisResults.resultsMB, &analysisResults.resultsSimMB, &analysisResults.resultsSimEtaMB, &analysisResults.resultsSimEtabgMB);
	TCanvas *peakFinderHT1canvas = fitSpectraBins(analysisSettingsNew.settingsHT1, analysisResults.resultsHT1, &analysisResults.resultsSimHT1, &analysisResults.resultsSimEtaHT1, &analysisResults.resultsSimEtabgHT1);
	TCanvas *peakFinderHT2canvas = fitSpectraBins(analysisSettingsNew.settingsHT2, analysisResults.resultsHT2, &analysisResults.resultsSimHT2, &analysisResults.resultsSimEtaHT2, &analysisResults.resultsSimEtabgHT2);

	if (analysisSettingsNew.savePtSlicesTrig) {
	    savePtSlices(analysisSettingsNew.PtSlicesTrigFilename
		, analysisResults.resultsMB.invlist, analysisResults.resultsHT1.invlist, analysisResults.resultsHT2.invlist
		, analysisSettingsNew.settingsMB.triggerThreshold, analysisSettingsNew.settingsHT1.triggerThreshold, analysisSettingsNew.settingsHT2.triggerThreshold
		, peakFinderMBcanvas, peakFinderHT1canvas, peakFinderHT2canvas
	    );
	}
	if (analysisSettingsNew.savePtSlicesTrigSim) {
	    savePtSlices(analysisSettingsNew.PtSlicesTrigSimFilename
		, analysisResults.resultsSimMB.invlist, analysisResults.resultsSimHT1.invlist, analysisResults.resultsSimHT2.invlist
		, analysisSettingsNew.settingsSimMB.triggerThreshold, analysisSettingsNew.settingsSimHT1.triggerThreshold, analysisSettingsNew.settingsSimHT2.triggerThreshold
		, peakFinderSimMBcanvas, peakFinderSimHT1canvas, peakFinderSimHT2canvas
	    );
	}
	if (analysisSettingsNew.savePtSlicesTrigSimEta) {
	    savePtSlices(analysisSettingsNew.PtSlicesTrigSimEtaFilename
		, analysisResults.resultsSimEtaMB.invlist, analysisResults.resultsSimEtaHT1.invlist, analysisResults.resultsSimEtaHT2.invlist
		, analysisSettingsNew.settingsSimEtaMB.triggerThreshold, analysisSettingsNew.settingsSimEtaHT1.triggerThreshold, analysisSettingsNew.settingsSimEtaHT2.triggerThreshold
		, peakFinderSimEtaMBcanvas, peakFinderSimEtaHT1canvas, peakFinderSimEtaHT2canvas
	    );
	}
	if (analysisSettingsNew.savePtSlicesTrigSimEtabg) {
	    savePtSlices(analysisSettingsNew.PtSlicesTrigSimEtabgFilename
		, analysisResults.resultsSimEtabgMB.invlist, analysisResults.resultsSimEtabgHT1.invlist, analysisResults.resultsSimEtabgHT2.invlist
		, analysisSettingsNew.settingsSimEtaMB.triggerThreshold, analysisSettingsNew.settingsSimEtaHT1.triggerThreshold, analysisSettingsNew.settingsSimEtaHT2.triggerThreshold
		, peakFinderSimEtabgMBcanvas, peakFinderSimEtabgHT1canvas, peakFinderSimEtabgHT2canvas
	    );
	}

	analysisResults.totSimuYieldMB = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuMB.begin();iter != analysisResults.areaSimuMB.end();++iter) analysisResults.totSimuYieldMB += (*iter).getValue();
	analysisResults.totSimuYieldHT1 = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuHT1.begin();iter != analysisResults.areaSimuHT1.end();++iter) analysisResults.totSimuYieldHT1 += (*iter).getValue();
	analysisResults.totSimuYieldHT2 = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuHT2.begin();iter != analysisResults.areaSimuHT2.end();++iter) analysisResults.totSimuYieldHT2 += (*iter).getValue();
	analysisResults.totSimuYieldEtaMB = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuEtaMB.begin();iter != analysisResults.areaSimuEtaMB.end();++iter) analysisResults.totSimuYieldEtaMB += (*iter).getValue();
	analysisResults.totSimuYieldEtaHT1 = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuEtaHT1.begin();iter != analysisResults.areaSimuEtaHT1.end();++iter) analysisResults.totSimuYieldEtaHT1 += (*iter).getValue();
	analysisResults.totSimuYieldEtaHT2 = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuEtaHT2.begin();iter != analysisResults.areaSimuEtaHT2.end();++iter) analysisResults.totSimuYieldEtaHT2 += (*iter).getValue();
	analysisResults.totSimuYieldEtabgMB = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuEtabgMB.begin();iter != analysisResults.areaSimuEtabgMB.end();++iter) analysisResults.totSimuYieldEtabgMB += (*iter).getValue();
	analysisResults.totSimuYieldEtabgHT1 = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuEtabgHT1.begin();iter != analysisResults.areaSimuEtabgHT1.end();++iter) analysisResults.totSimuYieldEtabgHT1 += (*iter).getValue();
	analysisResults.totSimuYieldEtabgHT2 = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.areaSimuEtabgHT2.begin();iter != analysisResults.areaSimuEtabgHT2.end();++iter) analysisResults.totSimuYieldEtabgHT2 += (*iter).getValue();
	if (analysisSettingsNew.print) {
		cout << "Total real yield: MinBias " << analysisResults.resultsMB.totalRawYield << ", HighTower-1 " << analysisResults.resultsHT1.totalRawYield << ", HighTower-2 " << analysisResults.resultsHT2.totalRawYield << endl;
		cout << "Total real eta yield: MinBias " << analysisResults.resultsMB.totalRawYieldEta << ", HighTower-1 " << analysisResults.resultsHT1.totalRawYieldEta << ", HighTower-2 " << analysisResults.resultsHT2.totalRawYieldEta << endl;
		cout << "Total simulated yield: All " << analysisResults.totSimuYield << ", MinBias " << analysisResults.totSimuYieldMB << ", HighTower-1 " << analysisResults.totSimuYieldHT1 << ", HighTower-2 " << analysisResults.totSimuYieldHT2 << endl;
		cout << "Total reco sim yield: MinBias " << analysisResults.resultsSimMB.totalRawYield << ", HighTower-1 " << analysisResults.resultsSimHT1.totalRawYield << ", HighTower-2 " << analysisResults.resultsSimHT2.totalRawYield << endl;
		cout << "Total simulated eta yield: All " << analysisResults.totSimuYieldEta << ", MinBias " << analysisResults.totSimuYieldEtaMB << ", HighTower-1 " << analysisResults.totSimuYieldEtaHT1 << ", HighTower-2 " << analysisResults.totSimuYieldEtaHT2 << endl;
		cout << "Total reco sim eta yield: MinBias " << analysisResults.resultsSimEtaMB.totalRawYield << ", HighTower-1 " << analysisResults.resultsSimEtaHT1.totalRawYield << ", HighTower-2 " << analysisResults.resultsSimEtaHT2.totalRawYield << endl;
		cout << "Total simulated etabg yield: All " << analysisResults.totSimuYieldEtabg << ", MinBias " << analysisResults.totSimuYieldEtabgMB << ", HighTower-1 " << analysisResults.totSimuYieldEtabgHT1 << ", HighTower-2 " << analysisResults.totSimuYieldEtabgHT2 << endl;
		cout << "Total reco sim etabg yield: MinBias " << analysisResults.resultsSimEtabgMB.totalRawYield << ", HighTower-1 " << analysisResults.resultsSimEtabgHT1.totalRawYield << ", HighTower-2 " << analysisResults.resultsSimEtabgHT2.totalRawYield << endl;
	}
	TBinStatistics totSimuYieldstat; totSimuYieldstat.setValue(analysisResults.totSimuYield);
	TBinStatistics totEmbedYieldMBstat; totEmbedYieldMBstat.setValue(analysisResults.resultsSimMB.totalRawYield);
	TBinStatistics totEmbedYieldHT1stat; totEmbedYieldHT1stat.setValue(analysisResults.resultsSimHT1.totalRawYield);
	TBinStatistics totEmbedYieldHT2stat; totEmbedYieldHT2stat.setValue(analysisResults.resultsSimHT2.totalRawYield);
	TBinStatistics totSimuYieldstatEta; totSimuYieldstatEta.setValue(analysisResults.totSimuYieldEta);
	TBinStatistics totEmbedYieldEtaMBstat; totEmbedYieldEtaMBstat.setValue(analysisResults.resultsSimEtaMB.totalRawYield);
	TBinStatistics totEmbedYieldEtaHT1stat; totEmbedYieldEtaHT1stat.setValue(analysisResults.resultsSimEtaHT1.totalRawYield);
	TBinStatistics totEmbedYieldEtaHT2stat; totEmbedYieldEtaHT2stat.setValue(analysisResults.resultsSimEtaHT2.totalRawYield);

	if (analysisSettingsNew.correctEffMatrix) {
		calculateBinPurityMatrix(analysisResults.resultsSimMB.simulatedPtlist, analysisResults.resultsMB.rawYield, analysisResults.areaSimuMB, analysisResults.resultsSimMB.rawYield, analysisResults.correctedMB, analysisResults.correctedAreaMB);
		calculateBinPurityMatrix(analysisResults.resultsSimHT1.simulatedPtlist, analysisResults.resultsHT1.rawYield, analysisResults.areaSimuHT1, analysisResults.resultsSimHT1.rawYield, analysisResults.correctedHT1, analysisResults.correctedAreaHT1);
		calculateBinPurityMatrix(analysisResults.resultsSimHT2.simulatedPtlist, analysisResults.resultsHT2.rawYield, analysisResults.areaSimuHT2, analysisResults.resultsSimHT2.rawYield, analysisResults.correctedHT2, analysisResults.correctedAreaHT2);
		if (false) {
			divideBins(analysisResults.correctedAreaMB, analysisResults.correctedMB, analysisResults.effMB, analysisResults.meanAcceptanceRatioMB, 0.0, false, false);
			divideBins(analysisResults.correctedAreaHT1, analysisResults.correctedHT1, analysisResults.effHT1, analysisResults.meanAcceptanceRatioHT1, 0.0, false, false);
			divideBins(analysisResults.correctedAreaHT2, analysisResults.correctedHT2, analysisResults.effHT2, analysisResults.meanAcceptanceRatioHT2, 0.0, false, false);
		} else {
			divideBins(analysisResults.resultsMB.rawYield, analysisResults.correctedMB, analysisResults.effMB, analysisResults.meanAcceptanceRatioMB, 0.0, false, false);
			divideBins(analysisResults.resultsHT1.rawYield, analysisResults.correctedHT1, analysisResults.effHT1, analysisResults.meanAcceptanceRatioHT1, 0.0, false, false);
			divideBins(analysisResults.resultsHT2.rawYield, analysisResults.correctedHT2, analysisResults.effHT2, analysisResults.meanAcceptanceRatioHT2, 0.0, false, false);
		}
	} else {
		divideBins(analysisResults.resultsSimMB.rawYield, analysisResults.areaSimuMB, analysisResults.effMB, analysisResults.meanAcceptanceRatioMB, 0.0, false, false);
		divideBins(analysisResults.resultsSimHT1.rawYield, analysisResults.areaSimuHT1, analysisResults.effHT1, analysisResults.meanAcceptanceRatioHT1, 0.0, false, false);
		divideBins(analysisResults.resultsSimHT2.rawYield, analysisResults.areaSimuHT2, analysisResults.effHT2, analysisResults.meanAcceptanceRatioHT2, 0.0, false, false);
	}
	if (analysisSettingsNew.smoothEffMB) {
		smoothPoints(analysisResults.effMB, analysisResults.effMBSmooth, analysisResults.effMBSmoothPoint);
	} else {
		analysisResults.effMBSmooth = analysisResults.effMB;
		analysisResults.effMBSmoothPoint = analysisResults.effMB;
	}
	if (analysisSettingsNew.smoothEffHT1) {
		smoothPoints(analysisResults.effHT1, analysisResults.effHT1Smooth, analysisResults.effHT1SmoothPoint);
	} else {
		analysisResults.effHT1Smooth = analysisResults.effHT1;
		analysisResults.effHT1SmoothPoint = analysisResults.effHT1;
	}
	if (analysisSettingsNew.smoothEffHT2) {
		smoothPoints(analysisResults.effHT2, analysisResults.effHT2Smooth, analysisResults.effHT2SmoothPoint);
	} else {
		analysisResults.effHT2Smooth = analysisResults.effHT2;
		analysisResults.effHT2SmoothPoint = analysisResults.effHT2;
	}
	if (analysisSettingsNew.dontCorrectEff) {
	    bin_stat_list_type tempMB, tempHT1, tempHT2;
	    resetBinsError(analysisResults.effMBSmooth, tempMB);
	    analysisResults.effMBSmooth = tempMB;
	    resetBinsError(analysisResults.effHT1Smooth, tempHT1);
	    analysisResults.effHT1Smooth = tempHT1;
	    resetBinsError(analysisResults.effHT2Smooth, tempHT2);
	    analysisResults.effHT2Smooth = tempHT2;
	}

	if (analysisSettingsNew.correctEffMatrixEta) {
		calculateBinPurityMatrix(analysisResults.resultsSimEtaMB.simulatedPtlist, analysisResults.resultsMB.rawYieldEta, analysisResults.areaSimuEtaMB, analysisResults.resultsSimMB.rawYield, analysisResults.correctedMBEta, analysisResults.correctedAreaMBEta);
		calculateBinPurityMatrix(analysisResults.resultsSimEtaHT1.simulatedPtlist, analysisResults.resultsHT1.rawYieldEta, analysisResults.areaSimuEtaHT1, analysisResults.resultsSimHT1.rawYield, analysisResults.correctedHT1Eta, analysisResults.correctedAreaHT1Eta);
		calculateBinPurityMatrix(analysisResults.resultsSimEtaHT2.simulatedPtlist, analysisResults.resultsHT2.rawYieldEta, analysisResults.areaSimuEtaHT2, analysisResults.resultsSimHT2.rawYield, analysisResults.correctedHT2Eta, analysisResults.correctedAreaHT2Eta);
		if (false) {
			divideBins(analysisResults.correctedAreaMBEta, analysisResults.correctedMBEta, analysisResults.effMBEta, analysisResults.meanAcceptanceRatioEtaMB, 0.0, false, false);
			divideBins(analysisResults.correctedAreaHT1Eta, analysisResults.correctedHT1Eta, analysisResults.effHT1Eta, analysisResults.meanAcceptanceRatioEtaHT1, 0.0, false, false);
			divideBins(analysisResults.correctedAreaHT2Eta, analysisResults.correctedHT2Eta, analysisResults.effHT2Eta, analysisResults.meanAcceptanceRatioEtaHT2, 0.0, false, false);
		} else {
			divideBins(analysisResults.resultsMB.rawYieldEta, analysisResults.correctedMBEta, analysisResults.effMBEta, analysisResults.meanAcceptanceRatioEtaMB, 0.0, false, false);
			divideBins(analysisResults.resultsHT1.rawYieldEta, analysisResults.correctedHT1Eta, analysisResults.effHT1Eta, analysisResults.meanAcceptanceRatioEtaHT1, 0.0, false, false);
			divideBins(analysisResults.resultsHT2.rawYieldEta, analysisResults.correctedHT2Eta, analysisResults.effHT2Eta, analysisResults.meanAcceptanceRatioEtaHT2, 0.0, false, false);
		}
	} else {
		divideBins(analysisResults.resultsSimEtaMB.rawYieldEta, analysisResults.areaSimuEtaMB, analysisResults.effMBEta, analysisResults.meanAcceptanceRatioEtaMB, 0.0, false, false);
		divideBins(analysisResults.resultsSimEtaHT1.rawYieldEta, analysisResults.areaSimuEtaHT1, analysisResults.effHT1Eta, analysisResults.meanAcceptanceRatioEtaHT1, 0.0, false, false);
		divideBins(analysisResults.resultsSimEtaHT2.rawYieldEta, analysisResults.areaSimuEtaHT2, analysisResults.effHT2Eta, analysisResults.meanAcceptanceRatioEtaHT2, 0.0, false, false);
	}
	if (analysisSettingsNew.smoothEffMBEta) {
		smoothPoints(analysisResults.effMBEta, analysisResults.effMBSmoothEta, analysisResults.effMBSmoothPointEta);
	} else {
		analysisResults.effMBSmoothEta = analysisResults.effMBEta;
		analysisResults.effMBSmoothPointEta = analysisResults.effMBEta;
	}
	if (analysisSettingsNew.smoothEffHT1Eta) {
		smoothPoints(analysisResults.effHT1Eta, analysisResults.effHT1SmoothEta, analysisResults.effHT1SmoothPointEta);
	} else {
		analysisResults.effHT1SmoothEta = analysisResults.effHT1Eta;
		analysisResults.effHT1SmoothPointEta = analysisResults.effHT1Eta;
	}
	if (analysisSettingsNew.smoothEffHT2Eta) {
		smoothPoints(analysisResults.effHT2Eta, analysisResults.effHT2SmoothEta, analysisResults.effHT2SmoothPointEta);
	} else {
		analysisResults.effHT2SmoothEta = analysisResults.effHT2Eta;
		analysisResults.effHT2SmoothPointEta = analysisResults.effHT2Eta;
	}
	if (analysisSettingsNew.dontCorrectEffEta) {
	    bin_stat_list_type tempMB, tempHT1, tempHT2;
	    resetBinsError(analysisResults.effMBSmoothEta, tempMB);
	    analysisResults.effMBSmoothEta = tempMB;
	    resetBinsError(analysisResults.effHT1SmoothEta, tempHT1);
	    analysisResults.effHT1SmoothEta = tempHT1;
	    resetBinsError(analysisResults.effHT2SmoothEta, tempHT2);
	    analysisResults.effHT2SmoothEta = tempHT2;
	}

	TLegend *legendRefmults = 0;
	TString nameRefmults = analysisSettingsNew.name + "_refmults";
	TString titleRefmults = analysisSettingsNew.title + ", RefMult per MinBias event vs. day;Day;RefMult per MinBias event";
	TH1F *histRefmults = 0;
	{
	    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsMB;
	    temp.legendTitle = "East FTPC RefMult per MinBias event";
	    temp.color = kBlue;
	    temp.fillColor = kBlue;
	    histRefmults = showList(&analysisResults.refMultEast, nameRefmults.Data(), titleRefmults.Data()
	        , histRefmults, &legendRefmults
	        , temp
	    );
	    temp.legendTitle = "West FTPC RefMult per MinBias event";
	    temp.color = kRed;
	    temp.fillColor = kRed;
	    histRefmults = showList(&analysisResults.refMultWest, nameRefmults.Data(), titleRefmults.Data()
	        , histRefmults, &legendRefmults
	        , temp
	    );
	}
	if (histRefmults) {
		histRefmults->SetMinimum(0);
		histRefmults->SetMaximum(20);
	}

	TLegend *legendPtShiftCorr = 0;
	TString ptShiftCorrName = analysisSettingsNew.name + "_ptshiftcorr";
	TString ptShiftCorrTitle = analysisSettingsNew.title + ", p_{T} shift correction due to finite bin width;p_{T}, GeV/c;(nominal/true)*(f(nominal)/f(true))";
	TH1F *histPtShiftCorr = 0;
	histPtShiftCorr = showLists(&analysisResults.ptShiftCorrMB, &analysisResults.ptShiftCorrHT1, &analysisResults.ptShiftCorrHT2
	    , ptShiftCorrName.Data(), ptShiftCorrTitle.Data(), histPtShiftCorr, &legendPtShiftCorr
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	    );
	if (histPtShiftCorr) {
		histPtShiftCorr->SetMinimum(0.6);
		histPtShiftCorr->SetMaximum(1.0);
		histPtShiftCorr->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	TLegend *legendPtShift = 0;
	TString ptShiftName = analysisSettingsNew.name + "_ptshift";
	TString ptShiftTitle = analysisSettingsNew.title + ", p_{T} shift due to finite bin width;p_{T}, GeV/c;(true - nominal) / nominal";
	TH1F *histPtShift = 0;
	histPtShift = showLists(&analysisResults.ptCenterShiftMB, &analysisResults.ptCenterShiftHT1, &analysisResults.ptCenterShiftHT2
	    , ptShiftName.Data(), ptShiftTitle.Data(), histPtShift, &legendPtShift
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	    );
	if (histPtShift) {
		histPtShift->SetMinimum(-0.1);
		histPtShift->SetMaximum(0.0);
		histPtShift->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	TLegend *legendEff = 0;
	TString effName = analysisSettingsNew.name + "_eff";
	TString effTitle = analysisSettingsNew.title + ", Efficiency and acceptance correction;p_{T}, GeV/c;Correction factor";
	TH1F *histEff = 0;
	histEff = showLists(&analysisResults.effMB, &analysisResults.effHT1, &analysisResults.effHT2
	    , effName.Data(), effTitle.Data(), histEff, &legendEff
	    , analysisSettingsNew.settingsSimMB, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.settingsSimHT2
	    );
	histEff = showLists(&analysisResults.effMBEta, &analysisResults.effHT1Eta, &analysisResults.effHT2Eta
	    , effName.Data(), effTitle.Data(), histEff, &legendEff
	    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
	    );
	if (histEff) {
		histEff->SetMinimum(0);
		histEff->SetMaximum(0.4);
		histEff->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
		if (analysisSettingsNew.smoothEffMB) {
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimMB;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effMBSmoothPoint
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT1) {
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimHT1;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT1SmoothPoint
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT2) {
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimHT2;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT2SmoothPoint
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffMBEta) {
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimEtaMB;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effMBSmoothPointEta
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT1Eta) {
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimEtaHT1;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT1SmoothPointEta
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT2Eta) {
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimEtaHT2;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT2SmoothPointEta
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
	}
	TLegend *legendBgJetToTotal = 0;
	TString BgJetToTotalName = analysisSettingsNew.name + "_BgJetToTotal";
	TString BgJetToTotalTitle = analysisSettingsNew.title + ", Jet-Mixed Background Fraction;p_{T}, GeV/c;Jet-Mixed Background Fraction";
	TH1F *histBgJetToTotal = 0;
	histBgJetToTotal = showLists(&analysisResults.resultsMB.bgJetToTotalFirst, &analysisResults.resultsHT1.bgJetToTotalFirst, &analysisResults.resultsHT2.bgJetToTotalFirst
	    , BgJetToTotalName.Data(), BgJetToTotalTitle.Data(), histBgJetToTotal, &legendBgJetToTotal
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	    );
	if (histBgJetToTotal) {
		histBgJetToTotal->SetMinimum(0);
		histBgJetToTotal->SetMaximum(1.0);
		histBgJetToTotal->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		/*
		if (histJetFracMB) {
		    histJetFracMB->SetLineColor(analysisSettingsNew.settingsMB.color);
		    histJetFracMB->Draw("HIST SAME");
		}
		if (histJetFracHT1) {
		    histJetFracHT1->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    histJetFracHT1->Draw("HIST SAME");
		}
		if (histJetFracHT2) {
		    histJetFracHT2->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    histJetFracHT2->Draw("HIST SAME");
		}
		*/
		if (analysisResults.resultsMB.jetBgFractionParamFit.func) {
		    analysisResults.resultsMB.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsMB.color);
		    analysisResults.resultsMB.jetBgFractionParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT1.jetBgFractionParamFit.func) {
		    analysisResults.resultsHT1.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    analysisResults.resultsHT1.jetBgFractionParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT2.jetBgFractionParamFit.func) {
		    analysisResults.resultsHT2.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    analysisResults.resultsHT2.jetBgFractionParamFit.func->Draw("SAME");
		}
	}

	//analysisResults.cpvCorrection.clear();
	//mergePoints(analysisResults.resultsMB.cpvCorrection, analysisResults.cpvCorrection);
	//mergePoints(analysisResults.resultsHT1.cpvCorrection, analysisResults.cpvCorrection);
	//mergePoints(analysisResults.resultsHT2.cpvCorrection, analysisResults.cpvCorrection);
	fitParameter(&analysisResults.resultsMB.cpvCorrection, 0, analysisSettingsNew.settingsMB.cpvCorrection, analysisResults.resultsMB.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitMB");
	fitParameter(&analysisResults.resultsHT1.cpvCorrection, 0, analysisSettingsNew.settingsHT1.cpvCorrection, analysisResults.resultsHT1.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT1");
	fitParameter(&analysisResults.resultsHT2.cpvCorrection, 0, analysisSettingsNew.settingsHT2.cpvCorrection, analysisResults.resultsHT2.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT2");
	analysisResults.cpvCorrection.clear();
	mergePoints(analysisResults.resultsMB.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrection);
	mergePoints(analysisResults.resultsHT1.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrection);
	mergePoints(analysisResults.resultsHT2.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrection);
	fitParameter(&analysisResults.cpvCorrection, 0, analysisSettingsNew.cpvCorrection, analysisResults.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfit");

	//analysisResults.cpvCorrectionSim.clear();
	//mergePoints(analysisResults.resultsSimMB.cpvCorrection, analysisResults.cpvCorrectionSim);
	//mergePoints(analysisResults.resultsSimHT1.cpvCorrection, analysisResults.cpvCorrectionSim);
	//mergePoints(analysisResults.resultsSimHT2.cpvCorrection, analysisResults.cpvCorrectionSim);
	fitParameter(&analysisResults.resultsSimMB.cpvCorrection, 0, analysisSettingsNew.settingsSimMB.cpvCorrection, analysisResults.resultsSimMB.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitMB_sim");
	fitParameter(&analysisResults.resultsSimHT1.cpvCorrection, 0, analysisSettingsNew.settingsSimHT1.cpvCorrection, analysisResults.resultsSimHT1.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT1_sim");
	fitParameter(&analysisResults.resultsSimHT2.cpvCorrection, 0, analysisSettingsNew.settingsSimHT2.cpvCorrection, analysisResults.resultsSimHT2.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT2_sim");
	analysisResults.cpvCorrectionSim.clear();
	mergePoints(analysisResults.resultsSimMB.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSim);
	mergePoints(analysisResults.resultsSimHT1.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSim);
	mergePoints(analysisResults.resultsSimHT2.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSim);
	fitParameter(&analysisResults.cpvCorrectionSim, 0, analysisSettingsNew.cpvCorrectionSim, analysisResults.cpvCorrectionFitSim, analysisSettingsNew.name + "_cpvfit_sim");

	//analysisResults.cpvCorrectionSimEta.clear();
	//mergePoints(analysisResults.resultsSimEtaMB.cpvCorrection, analysisResults.cpvCorrectionSimEta);
	//mergePoints(analysisResults.resultsSimEtaHT1.cpvCorrection, analysisResults.cpvCorrectionSimEta);
	//mergePoints(analysisResults.resultsSimEtaHT2.cpvCorrection, analysisResults.cpvCorrectionSimEta);
	fitParameter(&analysisResults.resultsSimEtaMB.cpvCorrection, 0, analysisSettingsNew.settingsSimEtaMB.cpvCorrection, analysisResults.resultsSimEtaMB.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitMB_simeta");
	fitParameter(&analysisResults.resultsSimEtaHT1.cpvCorrection, 0, analysisSettingsNew.settingsSimEtaHT1.cpvCorrection, analysisResults.resultsSimEtaHT1.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT1_simeta");
	fitParameter(&analysisResults.resultsSimEtaHT2.cpvCorrection, 0, analysisSettingsNew.settingsSimEtaHT2.cpvCorrection, analysisResults.resultsSimEtaHT2.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT2_simeta");
	analysisResults.cpvCorrectionSimEta.clear();
	mergePoints(analysisResults.resultsSimEtaMB.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSimEta);
	mergePoints(analysisResults.resultsSimEtaHT1.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSimEta);
	mergePoints(analysisResults.resultsSimEtaHT2.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSimEta);
	fitParameter(&analysisResults.cpvCorrectionSimEta, 0, analysisSettingsNew.cpvCorrectionSimEta, analysisResults.cpvCorrectionFitSimEta, analysisSettingsNew.name + "_cpvfit_simeta");

	TLegend *legendCpvCorrection = 0;
	TString CpvCorrectionName = analysisSettingsNew.name + "_CpvCorrection";
	TString CpvCorrectionTitle = analysisSettingsNew.title + ", CPV correction;p_{T}, GeV/c;CPV correction";
	TH1F *histCpvCorrection = 0;
	histCpvCorrection = showLists(&analysisResults.resultsMB.cpvCorrection, &analysisResults.resultsHT1.cpvCorrection, &analysisResults.resultsHT2.cpvCorrection
	    , CpvCorrectionName.Data(), CpvCorrectionTitle.Data(), histCpvCorrection, &legendCpvCorrection
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	    );
	if (histCpvCorrection) {
		histCpvCorrection->SetMinimum(0.6);
		histCpvCorrection->SetMaximum(1.0);
		histCpvCorrection->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		if (analysisResults.cpvCorrectionFit.func) {
		    analysisResults.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrection) legendCpvCorrection->AddEntry(analysisResults.cpvCorrectionFit.func, "Fit", "L");
		    if (analysisResults.cpvCorrectionFit.fitSigma > 0) {
			Float_t x = analysisResults.cpvCorrectionFit.func->GetXaxis()->GetXmax();
			Float_t y = analysisResults.cpvCorrectionFit.func->Eval(x);
			TH1F *histCpvFitErr = new TH1F(analysisSettingsNew.name + "_cpvfiterr", "fit err", 1, x - 1.5, x - 1.0);
			if (histCpvFitErr) {
			    histCpvFitErr->SetBinContent(1, y);
			    histCpvFitErr->SetBinError(1, analysisResults.cpvCorrectionFit.func->GetParError(0));
			    histCpvFitErr->SetMarkerColor(17);
			    histCpvFitErr->SetLineColor(17);
			    histCpvFitErr->SetFillColor(17);
			    histCpvFitErr->Draw("SAME ][ E2");
			    if (legendCpvCorrection) legendCpvCorrection->AddEntry(histCpvFitErr, Form("Par0 err = %.4f", analysisResults.cpvCorrectionFit.func->GetParError(0)), "EPF");
			}
			TH1F *histCpvMeanAbs = new TH1F(analysisSettingsNew.name + "_cpvmeanabs", "mean abs", 1, x - 1.0, x - 0.5);
			if (histCpvMeanAbs) {
			    histCpvMeanAbs->SetBinContent(1, y);
			    histCpvMeanAbs->SetBinError(1, analysisResults.cpvCorrectionFit.fitMeanAbs);
			    histCpvMeanAbs->SetMarkerColor(15);
			    histCpvMeanAbs->SetLineColor(15);
			    histCpvMeanAbs->SetFillColor(15);
			    histCpvMeanAbs->Draw("SAME ][ E2");
			    if (legendCpvCorrection) legendCpvCorrection->AddEntry(histCpvMeanAbs, Form("Mean |#Delta| = %.4f", analysisResults.cpvCorrectionFit.fitMeanAbs), "EPF");
			}
			TH1F *histCpvSigma = new TH1F(analysisSettingsNew.name + "_cpvsigma", "#sigma", 1, x - 0.5, x);
			if (histCpvSigma) {
			    histCpvSigma->SetBinContent(1, y);
			    histCpvSigma->SetBinError(1, analysisResults.cpvCorrectionFit.fitSigma);
			    histCpvSigma->SetMarkerColor(13);
			    histCpvSigma->SetLineColor(13);
			    histCpvSigma->SetFillColor(13);
			    histCpvSigma->Draw("SAME ][ E2");
			    if (legendCpvCorrection) legendCpvCorrection->AddEntry(histCpvSigma, Form("RMS = %.4f", analysisResults.cpvCorrectionFit.fitSigma), "EPF");
			}
		    }
		}
		if (analysisResults.resultsMB.cpvCorrectionFit.func) {
		    analysisResults.resultsMB.cpvCorrectionFit.func->SetLineColor(analysisSettingsNew.settingsMB.color);
		    analysisResults.resultsMB.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrection) legendCpvCorrection->AddEntry(analysisResults.resultsMB.cpvCorrectionFit.func, "Fit MB", "L");
		}
		if (analysisResults.resultsHT1.cpvCorrectionFit.func) {
		    analysisResults.resultsHT1.cpvCorrectionFit.func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    analysisResults.resultsHT1.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrection) legendCpvCorrection->AddEntry(analysisResults.resultsHT1.cpvCorrectionFit.func, "Fit HT1", "L");
		}
		if (analysisResults.resultsHT2.cpvCorrectionFit.func) {
		    analysisResults.resultsHT2.cpvCorrectionFit.func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    analysisResults.resultsHT2.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrection) legendCpvCorrection->AddEntry(analysisResults.resultsHT2.cpvCorrectionFit.func, "Fit HT2", "L");
		}
	}
	TLegend *legendCpvCorrectionParam7 = 0;
	TString CpvCorrectionParam7Name = analysisSettingsNew.name + "_CpvCorrectionParam7";
	TString CpvCorrectionParam7Title = analysisSettingsNew.title + ", CPV correction param7;p_{T}, GeV/c;CPV correction Param7";
	TH1F *histCpvCorrectionParam7 = 0;
	histCpvCorrectionParam7 = showLists(&analysisResults.resultsMB.cpvCorrectionParam7, &analysisResults.resultsHT1.cpvCorrectionParam7, &analysisResults.resultsHT2.cpvCorrectionParam7
	    , CpvCorrectionParam7Name.Data(), CpvCorrectionParam7Title.Data(), histCpvCorrectionParam7, &legendCpvCorrectionParam7
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	    );
	if (histCpvCorrectionParam7) {
		histCpvCorrectionParam7->SetMinimum(0.0);
		histCpvCorrectionParam7->SetMaximum(25.0);
		histCpvCorrectionParam7->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}
	TLegend *legendPosEmbed = 0;
	TString posEmbedName = analysisSettingsNew.name + "_position_embed";
	TString posEmbedTitle = analysisSettingsNew.title + ", Peak position and width vs. p_{T}, #pi^{0};p_{T}, GeV/c;Peak position and width, GeV/c^{2}";
	TH1F *histPosEmbed = 0;
	histPosEmbed = showLists(&analysisResults.resultsSimMB.peakPosition, &analysisResults.resultsSimHT1.peakPosition, &analysisResults.resultsSimHT2.peakPosition
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsSimMB, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.settingsSimHT2
	);
	histPosEmbed = showLists(&analysisResults.resultsMB.peakPosition, &analysisResults.resultsHT1.peakPosition, &analysisResults.resultsHT2.peakPosition
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	);
	histPosEmbed = showLists(&analysisResults.resultsSimMB.peakWidth, &analysisResults.resultsSimHT1.peakWidth, &analysisResults.resultsSimHT2.peakWidth
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsSimMB, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.settingsSimHT2
	);
	histPosEmbed = showLists(&analysisResults.resultsMB.peakWidth, &analysisResults.resultsHT1.peakWidth, &analysisResults.resultsHT2.peakWidth
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	);
	if (histPosEmbed) {
		histPosEmbed->SetMinimum(0);
		histPosEmbed->SetMaximum(0.4);
		histPosEmbed->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		Float_t Xmin = histPosEmbed->GetXaxis()->GetBinLowEdge(1);
		Float_t Xmax = histPosEmbed->GetXaxis()->GetBinUpEdge(histPosEmbed->GetXaxis()->GetNbins());
		TArrow *arGrMinvPt = new TArrow(Xmin, truePionMass, Xmax, truePionMass, 0.03, "-");
		arGrMinvPt->SetLineColor(1);
		arGrMinvPt->SetFillColor(1);
		arGrMinvPt->Draw("SAME");
	}
	TLegend *legendPosEmbedEta = 0;
	TString posEmbedEtaName = analysisSettingsNew.name + "_position_embed_eta";
	TString posEmbedEtaTitle = analysisSettingsNew.title + ", Peak position and width vs. p_{T}, #eta;p_{T}, GeV/c;Peak position and width, GeV/c^{2}";
	TH1F *histPosEmbedEta = 0;
	histPosEmbedEta = showLists(&analysisResults.resultsSimEtaMB.peakPositionEta, &analysisResults.resultsSimEtaHT1.peakPositionEta, &analysisResults.resultsSimEtaHT2.peakPositionEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
	);
	histPosEmbedEta = showLists(&analysisResults.resultsMB.peakPositionEta, &analysisResults.resultsHT1.peakPositionEta, &analysisResults.resultsHT2.peakPositionEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	);
	histPosEmbedEta = showLists(&analysisResults.resultsSimEtaMB.peakWidthEta, &analysisResults.resultsSimEtaHT1.peakWidthEta, &analysisResults.resultsSimEtaHT2.peakWidthEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
	);
	histPosEmbedEta = showLists(&analysisResults.resultsMB.peakWidthEta, &analysisResults.resultsHT1.peakWidthEta, &analysisResults.resultsHT2.peakWidthEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	);
	if (histPosEmbedEta) {
		histPosEmbedEta->SetMinimum(0);
		histPosEmbedEta->SetMaximum(1.0);
		histPosEmbedEta->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		Float_t Xmin = histPosEmbedEta->GetXaxis()->GetBinLowEdge(1);
		Float_t Xmax = histPosEmbedEta->GetXaxis()->GetBinUpEdge(histPosEmbed->GetXaxis()->GetNbins());
		TArrow *arGrMinvPtEta = new TArrow(Xmin, trueEtaMass, Xmax, trueEtaMass, 0.03, "-");
		arGrMinvPtEta->SetLineColor(1);
		arGrMinvPtEta->SetFillColor(1);
		arGrMinvPtEta->Draw("SAME");
	}
/*
		if (analysisSettingsNew.settingsSimMB.useFitRangeFit && analysisSettingsNew.settingsSimMB.fitLeftFunc) {
			analysisSettingsNew.settingsSimMB.fitLeftFunc->SetLineColor(analysisSettingsNew.settingsSimMB.color);
			analysisSettingsNew.settingsSimMB.fitLeftFunc->SetLineWidth(1);
			analysisSettingsNew.settingsSimMB.fitLeftFunc->SetLineStyle(2);
			analysisSettingsNew.settingsSimMB.fitLeftFunc->SetFillColor(analysisSettingsNew.settingsSimMB.fillColor);
			analysisSettingsNew.settingsSimMB.fitLeftFunc->Draw("SAME");
		}
		if (analysisSettingsNew.settingsSimMB.useFitRangeFit && analysisSettingsNew.settingsSimMB.fitRightFunc) {
			analysisSettingsNew.settingsSimMB.fitRightFunc->SetLineColor(analysisSettingsNew.settingsSimMB.color);
			analysisSettingsNew.settingsSimMB.fitRightFunc->SetLineWidth(1);
			analysisSettingsNew.settingsSimMB.fitRightFunc->SetLineStyle(2);
			analysisSettingsNew.settingsSimMB.fitRightFunc->SetFillColor(analysisSettingsNew.settingsSimMB.fillColor);
			analysisSettingsNew.settingsSimMB.fitRightFunc->Draw("SAME");
		}
		if (analysisSettingsNew.settingsSimHT1.useFitRangeFit && analysisSettingsNew.settingsSimHT1.fitLeftFunc) {
			analysisSettingsNew.settingsSimHT1.fitLeftFunc->SetLineColor(analysisSettingsNew.settingsSimHT1.color);
			analysisSettingsNew.settingsSimHT1.fitLeftFunc->SetLineWidth(1);
			analysisSettingsNew.settingsSimHT1.fitLeftFunc->SetLineStyle(2);
			analysisSettingsNew.settingsSimHT1.fitLeftFunc->SetFillColor(analysisSettingsNew.settingsSimHT1.fillColor);
			analysisSettingsNew.settingsSimHT1.fitLeftFunc->Draw("SAME");
		}
		if (analysisSettingsNew.settingsSimHT1.useFitRangeFit && analysisSettingsNew.settingsSimHT1.fitRightFunc) {
			analysisSettingsNew.settingsSimHT1.fitRightFunc->SetLineColor(analysisSettingsNew.settingsSimHT1.color);
			analysisSettingsNew.settingsSimHT1.fitRightFunc->SetLineWidth(1);
			analysisSettingsNew.settingsSimHT1.fitRightFunc->SetLineStyle(2);
			analysisSettingsNew.settingsSimHT1.fitRightFunc->SetFillColor(analysisSettingsNew.settingsSimHT1.fillColor);
			analysisSettingsNew.settingsSimHT1.fitRightFunc->Draw("SAME");
		}
		if (analysisSettingsNew.settingsSimHT2.useFitRangeFit && analysisSettingsNew.settingsSimHT2.fitLeftFunc) {
			analysisSettingsNew.settingsSimHT2.fitLeftFunc->SetLineColor(analysisSettingsNew.settingsSimHT2.color);
			analysisSettingsNew.settingsSimHT2.fitLeftFunc->SetLineWidth(1);
			analysisSettingsNew.settingsSimHT2.fitLeftFunc->SetLineStyle(2);
			analysisSettingsNew.settingsSimHT2.fitLeftFunc->SetFillColor(analysisSettingsNew.settingsSimHT2.fillColor);
			analysisSettingsNew.settingsSimHT2.fitLeftFunc->Draw("SAME");
		}
		if (analysisSettingsNew.settingsSimHT2.useFitRangeFit && analysisSettingsNew.settingsSimHT2.fitRightFunc) {
			analysisSettingsNew.settingsSimHT2.fitRightFunc->SetLineColor(analysisSettingsNew.settingsSimHT2.color);
			analysisSettingsNew.settingsSimHT2.fitRightFunc->SetLineWidth(1);
			analysisSettingsNew.settingsSimHT2.fitRightFunc->SetLineStyle(2);
			analysisSettingsNew.settingsSimHT2.fitRightFunc->SetFillColor(analysisSettingsNew.settingsSimHT2.fillColor);
			analysisSettingsNew.settingsSimHT2.fitRightFunc->Draw("SAME");
		}
*/
		/*
		if (analysisSettingsNew.settingsSimMB.useFitRangeFit || analysisSettingsNew.settingsSimSimHT1.useFitRangeFit || analysisSettingsNew.settingsSimSimHT2.useFitRangeFit) {
		    plots.showList(fitFuncLowUsedEmbedMB, fitFuncLowUsedEmbedHT1, fitFuncLowUsedEmbedHT2, posEmbedName.Data(), posEmbedTitle.Data(), "L E X0", histPosEmbed);
		    plots.showList(fitFuncHighUsedEmbedMB, fitFuncHighUsedEmbedHT1, fitFuncHighUsedEmbedHT2, posEmbedName.Data(), posEmbedTitle.Data(), "L E X0", histPosEmbed);
		}
		*/
	//}
/*
	bin_stat_list_type positionDiffMB, positionDiffHT1, positionDiffHT2;
	substractBins(positionMB, positionEmbedMB, positionDiffMB);
	substractBins(positionHT1, positionEmbedHT1, positionDiffHT1);
	substractBins(positionHT2, positionEmbedHT2, positionDiffHT2);
	TString positionDiffname = "position_diff"; positionDiffname += name;
	TString positionDifftitle = "Peak position difference vs. p_{T}"; positionDifftitle += title; positionDifftitle += ";p_{T}, GeV/c;Peak position: real - embed, GeV/c^{2}";
	plots.show = showReal || showSim;
	TH1F *histPosDiff = plots.showList(positionDiffMB, positionDiffHT1, positionDiffHT2, positionDiffname.Data(), positionDifftitle.Data(), "PZ");
	if (histPosDiff) {
		//histPosDiff->SetMinimum(0);
		//histPosDiff->SetMaximum(0.4);
	}
	bin_stat_list_type positionDivMB, positionDivHT1, positionDivHT2;
	divideBins(positionMB, positionEmbedMB, positionDivMB, 1.0, 0.0, false, false);
	divideBins(positionHT1, positionEmbedHT1, positionDivHT1, 1.0, 0.0, false, false);
	divideBins(positionHT2, positionEmbedHT2, positionDivHT2, 1.0, 0.0, false, false);
	TString positionDivname = "position_div"; positionDivname += name;
	TString positionDivtitle = "Peak position ratio vs. p_{T}"; positionDivtitle += title; positionDivtitle += ";p_{T}, GeV/c;Peak position: real / embed";
	plots.show = showReal || showSim;
	TH1F *histPosDiv = plots.showList(positionDivMB, positionDivHT1, positionDivHT2, positionDivname.Data(), positionDivtitle.Data(), "PZ");
	if (histPosDiv) {
		//histPosDiv->SetMinimum(0);
		//histPosDiv->SetMaximum(2);
	}

        TLegend *legendPar3 = 0;^
                TH1F *histBgPar3 = plots.showList(bgPar3MBleft, bgPar3HT1left, bgPar3HT2left, name, title, opt, 0, &legendPar3);^
                        if (histBgPar3) {^
                                        plots.markerStyle = 4;^
                                                        plots.showList(bgPar3MBright, bgPar3HT1right, bgPar3HT2right, name, title, opt, histBgPar3);^
                                                                        plots.markerStyle = 8;^
                                                                                        fBgPar3MB->SetLineWidth(1);^
                                                                                                        fBgPar3MB->SetLineColor(MBcolor);^
                                                                                                                        TF1 *fBgPar3MBr = new TF1(*fBgPar3MB);^
                                                                                                                                        fBgPar3MBr->SetRange(BgPar3MBHigh, plots.Xmax);^
                                                                                                                                                        fBgPar3MBr->SetLineColor(MBFillColor);^
                                                                                                                                                                        fBgPar3MBr->SetLineStyle(2);^
                                                                                                                                                                                        fBgPar3HT1->SetLineWidth(1);^
                                                                                                                                                                                                        fBgPar3HT1->SetLineColor(HT1color);^
                                                                                                                                                                                                                        TF1 *fBgPar3HT1r = new TF1(*fBgPar3HT1);^
                                                                                                                                                                                                                                        fBgPar3HT1r->SetRange(BgPar3HT1High, plots.Xmax);^
                                                                                                                                                                                                                                                        fBgPar3HT1r->SetLineColor(HT1FillColor);^
                                                                                                                                                                                                                                                                        fBgPar3HT1r->SetLineStyle(2);^
                                                                                                                                                                                                                                                                                        fBgPar3HT2->SetLineWidth(1);^
                                                                                                                                                                                                                                                                                                        fBgPar3HT2->SetLineColor(HT2color);^
                                                                                                                                                                                                                                                                                                                        TF1 *fBgPar3HT2r = new TF1(*fBgPar3HT2);^
                                                                                                                                                                                                                                                                                                                                        fBgPar3HT2r->SetRange(BgPar3HT2High, plots.Xmax);^
                                                                                                                                                                                                                                                                                                                                                        fBgPar3HT2r->SetLineColor(HT2FillColor);^
                                                                                                                                                                                                                                                                                                                                                                        fBgPar3HT2r->SetLineStyle(2);^
                                                                                                                                                                                                                                                                                                                                                                                        fBgPar3MB->Draw("SAME");^
                                                                                                                                                                                                                                                                                                                                                                                                        fBgPar3HT1->Draw("SAME");^
                                                                                                                                                                                                                                                                                                                                                                                                                        fBgPar3HT2->Draw("SAME");               ^
                                                                                                                                                                                                                                                                                                                                                                                                                                        fBgPar3MBr->Draw("SAME");^
                                                                                                                                                                                                                                                                                                                                                                                                                                                        fBgPar3HT1r->Draw("SAME");^
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        fBgPar3HT2r->Draw("SAME");              ^
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        histBgPar3->SetMinimum(-20);^
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        histBgPar3->SetMaximum(50);^
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
	TString bgPar3name = "bgPar3"; bgPar3name += name;
	TString bgPar3title = "Bg parameter3 vs. p_{T}"; bgPar3title += title; bgPar3title += ";p_{T}, GeV/c;Bg parameter3";
	plots.show = showReal && (smoothSettings.smoothMB.useParametersExtrapolation || smoothSettings.smoothHT1.useParametersExtrapolation || smoothSettings.smoothHT2.useParametersExtrapolation);
	plotBgPar(plots, bgPar3MBleft, bgPar3MBright, bgPar3HT1left, bgPar3HT1right, bgPar3HT2left, bgPar3HT2right, bgPar3name.Data(), bgPar3title.Data(), "PZ", fBgPar3MB, fBgPar3HT1, fBgPar3HT2, 4.0, 7.0, 8.0);

	TString bgPar4name = "bgPar4"; bgPar4name += name;
	TString bgPar4title = "Bg parameter4 vs. p_{T}"; bgPar4title += title; bgPar4title += ";p_{T}, GeV/c;Bg parameter4";
	plots.show = showReal && (smoothSettings.smoothMB.useParametersExtrapolation || smoothSettings.smoothHT1.useParametersExtrapolation || smoothSettings.smoothHT2.useParametersExtrapolation);
	plotBgPar(plots, bgPar4MBleft, bgPar4MBright, bgPar4HT1left, bgPar4HT1right, bgPar4HT2left, bgPar4HT2right, bgPar4name.Data(), bgPar4title.Data(), "PZ", fBgPar4MB, fBgPar4HT1, fBgPar4HT2, 4.0, 7.0, 8.0);

	TString sigBgname = "sigBg"; sigBgname += name;
	TString sigBgtitle = "Signal/Bg vs. p_{T}"; sigBgtitle += title; sigBgtitle += ";p_{T}, GeV/c;Signal/Bg";
	plots.show = showReal && (smoothSettings.smoothMB.useSigBgExtrapolation || smoothSettings.smoothHT1.useSigBgExtrapolation || smoothSettings.smoothHT2.useSigBgExtrapolation);
	plotBgPar(plots, sigBgMBleft, sigBgMBright, sigBgHT1left, sigBgHT1right, sigBgHT2left, sigBgHT2right, sigBgname.Data(), sigBgtitle.Data(), "PZ", fSigBgMB, fSigBgHT1, fSigBgHT2, 4.0, 7.0, 8.0);


	TString sig_bgembedname = "sig_bg_embed"; sig_bgembedname += name;
	TString sig_bgembedtitle = "Signal / Background - Embedding"; sig_bgembedtitle += title; sig_bgembedtitle += ";p_{T}, GeV/c;signal/background";
	plots.show = showSim;
	bin_stat_list_type sig_bgEmbedMB, sig_bgEmbedHT1, sig_bgEmbedHT2;
	divideBins(areaEmbedMB, peakFinderEmbedMB.getBackground(), sig_bgEmbedMB, 1.0, 0.0, false, false);
	divideBins(areaEmbedHT1, peakFinderEmbedHT1.getBackground(), sig_bgEmbedHT1, 1.0, 0.0, false, false);
	divideBins(areaEmbedHT2, peakFinderEmbedHT2.getBackground(), sig_bgEmbedHT2, 1.0, 0.0, false, false);
	plots.showList(sig_bgEmbedMB, sig_bgEmbedHT1, sig_bgEmbedHT2, sig_bgembedname.Data(), sig_bgembedtitle.Data(), "P");

	TString sig_bgname = "sig_bg"; sig_bgname += name;
	TString sig_bgtitle = "Signal / Background"; sig_bgtitle += title; sig_bgtitle += ";p_{T}, GeV/c;signal/background";
	plots.show = showReal;
	bin_stat_list_type sig_bgMB, sig_bgHT1, sig_bgHT2;
	divideBins(areaMB, peakFinderMB.getBackground(), sig_bgMB, 1.0, 0.0, false, false);
	divideBins(areaHT1, peakFinderHT1.getBackground(), sig_bgHT1, 1.0, 0.0, false, false);
	divideBins(areaHT2, peakFinderHT2.getBackground(), sig_bgHT2, 1.0, 0.0, false, false);
	plots.showList(sig_bgMB, sig_bgHT1, sig_bgHT2, sig_bgname.Data(), sig_bgtitle.Data(), "P");
*/

	TLegend *legendYieldEmbed = 0;
	TString yieldEmbedName = analysisSettingsNew.name + "_yield_embed";
	TString yieldEmbedTitle = analysisSettingsNew.title + ", Raw yield vs. p_{T} - simulation;p_{T}, GeV/c;Yield";
	TH1F *histYieldEmbed = 0;
	histYieldEmbed = showLists(&analysisResults.areaSimuMB, &analysisResults.areaSimuHT1, &analysisResults.areaSimuHT2
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimMB, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.settingsSimHT2
	);
	histYieldEmbed = showLists(&analysisResults.resultsSimMB.rawYield, &analysisResults.resultsSimHT1.rawYield, &analysisResults.resultsSimHT2.rawYield
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimMB, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.settingsSimHT2
	);
	histYieldEmbed = showLists(&analysisResults.areaSimuEtaMB, &analysisResults.areaSimuEtaHT1, &analysisResults.areaSimuEtaHT2
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
	);
	histYieldEmbed = showLists(&analysisResults.resultsSimEtaMB.rawYieldEta, &analysisResults.resultsSimEtaHT1.rawYieldEta, &analysisResults.resultsSimEtaHT2.rawYieldEta
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
	);
	if (histYieldEmbed) {
		if (histYieldEmbed->GetMinimum() <= 0) histYieldEmbed->SetMinimum(1e-11);
		gPad->SetLogy();
	}

	TLegend *legendYield = 0;
	TString yieldName = analysisSettingsNew.name + "_yield";
	TString yieldTitle = analysisSettingsNew.title + ", Raw yield vs. p_{T} - data, #pi^{0} and #eta;p_{T}, GeV/c;Yield";
	TH1F *histYield = 0;
	histYield = showLists(&analysisResults.resultsMB.rawYield, &analysisResults.resultsHT1.rawYield, &analysisResults.resultsHT2.rawYield
	    , yieldName.Data(), yieldTitle.Data(), histYield, &legendYield
	    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
	);
	histYield = showLists(&analysisResults.resultsMB.rawYieldEta, &analysisResults.resultsHT1.rawYieldEta, &analysisResults.resultsHT2.rawYieldEta
	    , yieldName.Data(), yieldTitle.Data(), histYield, &legendYield
	    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
	);
	if (histYield) {
		if (histYield->GetMinimum() <= 0) histYield->SetMinimum(1e-11);
		gPad->SetLogy();
	}
/*
	if (histYield) {
		if (correctEffMatrix) plots.showList(correctedAreaMB, correctedAreaHT1, correctedAreaHT2, yieldName.Data(), yieldTitle.Data(), "L", histYield, &legendYield);
		if (correctEffMatrix) plots.showList(correctedMB, correctedHT1, correctedHT2, yieldName.Data(), yieldTitle.Data(), "L", histYield, &legendYield);
		if (histYield->GetMinimum() <= 0) histYield->SetMinimum(1e-11);
		gPad->SetLogy();
	}
*/
	divideBins(analysisResults.resultsSimMB.rawYield, analysisResults.effMBSmooth, analysisResults.spectrumMBembed, 0.0, TMath::TwoPi() * analysisResults.resultsSimMB.evNum * analysisResults.resultsSimMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsSimMB.deltaEta * (1.0 / analysisResults.resultsSimMB.vertexFindingEff), true, false);
	divideBins(analysisResults.resultsSimHT1.rawYield, analysisResults.effHT1Smooth, analysisResults.spectrumHT1embed, 0.0, TMath::TwoPi() * analysisResults.resultsSimHT1.evNum * analysisResults.resultsSimHT1.prescaleUsed * truePionBranchingRatio * analysisResults.resultsSimHT1.deltaEta * (1.0 / analysisResults.resultsSimHT1.vertexFindingEff), true, false);
	divideBins(analysisResults.resultsSimHT2.rawYield, analysisResults.effHT2Smooth, analysisResults.spectrumHT2embed, 0.0, TMath::TwoPi() * analysisResults.resultsSimHT2.evNum * analysisResults.resultsSimHT2.prescaleUsed * truePionBranchingRatio * analysisResults.resultsSimHT2.deltaEta * (1.0 / analysisResults.resultsSimHT2.vertexFindingEff), true, false);

	divideBins(analysisResults.resultsSimEtaMB.rawYieldEta, analysisResults.effMBSmoothEta, analysisResults.spectrumMBembedEta, 0.0, TMath::TwoPi() * analysisResults.resultsSimEtaMB.evNum * analysisResults.resultsSimEtaMB.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsSimEtaMB.deltaEta * (1.0 / analysisResults.resultsSimEtaMB.vertexFindingEff), true, false);
	divideBins(analysisResults.resultsSimEtaHT1.rawYieldEta, analysisResults.effHT1SmoothEta, analysisResults.spectrumHT1embedEta, 0.0, TMath::TwoPi() * analysisResults.resultsSimEtaHT1.evNum * analysisResults.resultsSimEtaHT1.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsSimEtaHT1.deltaEta * (1.0 / analysisResults.resultsSimEtaHT1.vertexFindingEff), true, false);
	divideBins(analysisResults.resultsSimEtaHT2.rawYieldEta, analysisResults.effHT2SmoothEta, analysisResults.spectrumHT2embedEta, 0.0, TMath::TwoPi() * analysisResults.resultsSimEtaHT2.evNum * analysisResults.resultsSimEtaHT2.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsSimEtaHT2.deltaEta * (1.0 / analysisResults.resultsSimEtaHT2.vertexFindingEff), true, false);

	if (analysisSettingsNew.settingsSimMB.correctCpv) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFitSim.func) {
		divideBinsFunc(analysisResults.spectrumMBembed, analysisResults.cpvCorrectionFitSim.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsSimMB.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumMBembed, analysisResults.resultsSimMB.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumMBembed, analysisResults.resultsSimMB.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumMBembed = temp;
	}
	if (analysisSettingsNew.settingsSimHT1.correctCpv) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFitSim.func) {
		divideBinsFunc(analysisResults.spectrumHT1embed, analysisResults.cpvCorrectionFitSim.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsSimHT1.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT1embed, analysisResults.resultsSimHT1.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumHT1embed, analysisResults.resultsSimHT1.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT1embed = temp;
	}
	if (analysisSettingsNew.settingsSimHT2.correctCpv) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFitSim.func) {
		divideBinsFunc(analysisResults.spectrumHT2embed, analysisResults.cpvCorrectionFitSim.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsSimHT2.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT2embed, analysisResults.resultsSimHT2.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumHT2embed, analysisResults.resultsSimHT2.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT2embed = temp;
	}

	if (analysisSettingsNew.settingsSimEtaMB.correctCpvEta) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFitSimEta.func) {
		divideBinsFunc(analysisResults.spectrumMBembedEta, analysisResults.cpvCorrectionFitSimEta.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsSimEtaMB.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumMBembedEta, analysisResults.resultsSimEtaMB.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumMBembedEta, analysisResults.resultsSimEtaMB.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumMBembedEta = temp;
	}
	if (analysisSettingsNew.settingsSimEtaHT1.correctCpvEta) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFitSimEta.func) {
		divideBinsFunc(analysisResults.spectrumHT1embedEta, analysisResults.cpvCorrectionFitSimEta.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsSimEtaHT1.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT1embedEta, analysisResults.resultsSimEtaHT1.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumHT1embedEta, analysisResults.resultsSimEtaHT1.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT1embedEta = temp;
	}
	if (analysisSettingsNew.settingsSimEtaHT2.correctCpvEta) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFitSimEta.func) {
		divideBinsFunc(analysisResults.spectrumHT2embedEta, analysisResults.cpvCorrectionFitSimEta.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsSimEtaHT2.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT2embedEta, analysisResults.resultsSimEtaHT2.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumHT2embedEta, analysisResults.resultsSimEtaHT2.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT2embedEta = temp;
	}

	mergePoints(analysisResults.spectrumMBembed, analysisResults.spectrumembed);
	mergePoints(analysisResults.spectrumHT2embed, analysisResults.spectrumembed);
	mergePoints(analysisResults.spectrumHT1embed, analysisResults.spectrumembed);

	mergePoints(analysisResults.spectrumMBembedEta, analysisResults.spectrumembedEta);
	mergePoints(analysisResults.spectrumHT2embedEta, analysisResults.spectrumembedEta);
	mergePoints(analysisResults.spectrumHT1embedEta, analysisResults.spectrumembedEta);

	TString spectrumEmbedName = analysisSettingsNew.name + "_spectrum_embed";
	TString spectrumEmbedTitle = analysisSettingsNew.title + ", #pi^{0} Inclusive p_{T} Spectrum (simulation);p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}n}{d#etadp_{T}}";
	TLegend *legendSpectrumEmbed = 0;
	TH1F *histSpectrumEmbed = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
	    histSpectrumEmbed = showLists(&analysisResults.spectrumMBembed, &analysisResults.spectrumHT1embed, &analysisResults.spectrumHT2embed
		, spectrumEmbedName.Data(), spectrumEmbedTitle.Data()
		, histSpectrumEmbed, &legendSpectrumEmbed
		, analysisSettingsNew.settingsSimMB, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.settingsSimHT2
	    );
	} else {
		TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimMB;
		temp.legendTitle = "STAR d+Au #pi^{0}";
		histSpectrumEmbed = showList(&analysisResults.spectrumembed
		    , spectrumEmbedName.Data(), spectrumEmbedTitle.Data(), histSpectrumEmbed, &legendSpectrumEmbed
		    , temp
		);
	}
	if (histSpectrumEmbed) {
		histSpectrumEmbed->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
		histSpectrumEmbed->GetYaxis()->SetLimits(1.0e-11, 1.0);
		gPad->SetLogy();
		if (analysisSettingsNew.showWeightUsed) {
			TF1 *funcw = weightUsed.createFunc(true);
			TF1 *funcwCopy = 0;
			if (funcw) {
				funcw->SetLineWidth(1);
				funcw->SetLineColor(11);
				funcw->SetLineStyle(1);
				funcwCopy = funcw->DrawCopy("SAME");
			}
			if (legendSpectrumEmbed) legendSpectrumEmbed->AddEntry(funcwCopy, "Weight", "L");
		}
	}
	showResultsOthers(histSpectrumEmbed, legendSpectrumEmbed, analysisSettingsNew.showpQCDWeight, analysisSettingsNew.usePPpQCD);

	if (analysisSettingsNew.multiplicityDistributions || false) {
		analysisResults.spectrumMB = analysisResults.resultsMB.rawYield;
		analysisResults.spectrumHT1 = analysisResults.resultsHT1.rawYield;
		analysisResults.spectrumHT2 = analysisResults.resultsHT2.rawYield;
	} else {
		if (analysisSettingsNew.dontCorrectEff && false) {
			TF1 funcOne("funcOne", "1");
			divideBinsFunc(analysisResults.resultsMB.rawYield, &funcOne, analysisResults.spectrumMB, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff), true, false, false);
			divideBinsFunc(analysisResults.resultsHT1.rawYield, &funcOne, analysisResults.spectrumHT1, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff), true, false, false);
			divideBinsFunc(analysisResults.resultsHT2.rawYield, &funcOne, analysisResults.spectrumHT2, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff), true, false, false);
		} else {
			divideBins(analysisResults.resultsMB.rawYield, analysisResults.effMBSmooth, analysisResults.spectrumMB, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff), true, false);
			divideBins(analysisResults.resultsHT1.rawYield, analysisResults.effHT1Smooth, analysisResults.spectrumHT1, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff), true, false);
			divideBins(analysisResults.resultsHT2.rawYield, analysisResults.effHT2Smooth, analysisResults.spectrumHT2, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff), true, false);
		}
	}

	if (analysisSettingsNew.multiplicityDistributions || false) {
		analysisResults.spectrumMBEta = analysisResults.resultsMB.rawYieldEta;
		analysisResults.spectrumHT1Eta = analysisResults.resultsHT1.rawYieldEta;
		analysisResults.spectrumHT2Eta = analysisResults.resultsHT2.rawYieldEta;
	} else {
		if (analysisSettingsNew.dontCorrectEffEta && false) {
			TF1 funcOne("funcOne", "1");
			divideBinsFunc(analysisResults.resultsMB.rawYieldEta, &funcOne, analysisResults.spectrumMBEta, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff), true, false, false);
			divideBinsFunc(analysisResults.resultsHT1.rawYieldEta, &funcOne, analysisResults.spectrumHT1Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff), true, false, false);
			divideBinsFunc(analysisResults.resultsHT2.rawYieldEta, &funcOne, analysisResults.spectrumHT2Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff), true, false, false);
		} else {
			divideBins(analysisResults.resultsMB.rawYieldEta, analysisResults.effMBSmoothEta, analysisResults.spectrumMBEta, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff), true, false);
			divideBins(analysisResults.resultsHT1.rawYieldEta, analysisResults.effHT1SmoothEta, analysisResults.spectrumHT1Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff), true, false);
			divideBins(analysisResults.resultsHT2.rawYieldEta, analysisResults.effHT2SmoothEta, analysisResults.spectrumHT2Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff), true, false);
		}
	}

	if (analysisSettingsNew.settingsMB.correctCpv) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumMB, analysisResults.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsMB.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumMB, analysisResults.resultsMB.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumMB, analysisResults.resultsMB.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumMB = temp;
	}
	if (analysisSettingsNew.settingsHT1.correctCpv) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT1, analysisResults.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsHT1.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT1, analysisResults.resultsHT1.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumHT1, analysisResults.resultsHT1.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT1 = temp;
	}
	if (analysisSettingsNew.settingsHT2.correctCpv) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT2, analysisResults.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsHT2.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT2, analysisResults.resultsHT2.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
    		divideBins(analysisResults.spectrumHT2, analysisResults.resultsHT2.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT2 = temp;
	}

	if (analysisSettingsNew.settingsMB.correctCpvEta) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumMBEta, analysisResults.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsMB.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumMBEta, analysisResults.resultsMB.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumMBEta, analysisResults.resultsMB.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumMBEta = temp;
	}
	if (analysisSettingsNew.settingsHT1.correctCpvEta) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT1Eta, analysisResults.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsHT1.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT1Eta, analysisResults.resultsHT1.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
	        divideBins(analysisResults.spectrumHT1Eta, analysisResults.resultsHT1.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT1Eta = temp;
	}
	if (analysisSettingsNew.settingsHT2.correctCpvEta) {
	    bin_stat_list_type temp;
	    if (analysisResults.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT2Eta, analysisResults.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else if (analysisResults.resultsHT2.cpvCorrectionFit.func) {
		divideBinsFunc(analysisResults.spectrumHT2Eta, analysisResults.resultsHT2.cpvCorrectionFit.func, temp, 1.0, 0.0, false, false, false);
	    } else {
		divideBins(analysisResults.spectrumHT2Eta, analysisResults.resultsHT2.cpvCorrection, temp, 1.0, 0.0, false, false);
	    }
	    analysisResults.spectrumHT2Eta = temp;
	}

	mergePoints(analysisResults.spectrumMB, analysisResults.spectrum);
	mergePoints(analysisResults.spectrumHT2, analysisResults.spectrum);
	mergePoints(analysisResults.spectrumHT1, analysisResults.spectrum);

	mergePoints(analysisResults.spectrumMBEta, analysisResults.spectrumEta);
	mergePoints(analysisResults.spectrumHT2Eta, analysisResults.spectrumEta);
	mergePoints(analysisResults.spectrumHT1Eta, analysisResults.spectrumEta);

	{
	    TF1 funcOne("funcOne", "1");
	    Float_t cx = (analysisSettingsNew.crossectionMB != 0) ? (analysisSettingsNew.crossectionMBTotalFraction / analysisSettingsNew.crossectionMB) : 1.0;

	    divideBinsFunc(analysisResults.spectrumMB, &funcOne, analysisResults.crossectionMB, cx, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT1, &funcOne, analysisResults.crossectionHT1, cx, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT2, &funcOne, analysisResults.crossectionHT2, cx, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrum, &funcOne, analysisResults.crossection, cx, 0.0, false, false, false);

	    divideBinsFunc(analysisResults.spectrumMBEta, &funcOne, analysisResults.crossectionMBEta, cx, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT1Eta, &funcOne, analysisResults.crossectionHT1Eta, cx, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT2Eta, &funcOne, analysisResults.crossectionHT2Eta, cx, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumEta, &funcOne, analysisResults.crossectionEta, cx, 0.0, false, false, false);
	}

	divideBins(analysisResults.spectrumEta, analysisResults.spectrum, analysisResults.etaToPi0Ratio, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumMBEta, analysisResults.spectrumMB, analysisResults.etaToPi0RatioMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT1Eta, analysisResults.spectrumHT1, analysisResults.etaToPi0RatioHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT2Eta, analysisResults.spectrumHT2, analysisResults.etaToPi0RatioHT2, 1.0, 0.0, false, false);

	TString spectrumname = analysisSettingsNew.name + (analysisSettingsNew.plotCrossection ? "_crossection" : "_spectrum");
	//TString spectrumtitle = analysisSettingsNew.title + ", #pi^{0} Inclusive p_{T} spectrum;p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}n}{d#etadp_{T}}";
	TString spectrumtitle = analysisSettingsNew.title + (analysisSettingsNew.plotCrossection ? ", #pi^{0} invariant cross section;p_{T}, GeV/c;E#frac{d^{3}#sigma}{d^{3}p}, mb*GeV^{-2}*c^{3}" : ", #pi^{0} Inclusive p_{T} spectrum;p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}n}{d#etadp_{T}}");
	TLegend *legendSpectrum = 0;
	TH1F *histSpectrum = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		//histSpectrum = showLists(&analysisResults.spectrumMB, &analysisResults.spectrumHT1, &analysisResults.spectrumHT2
		histSpectrum = showLists(analysisSettingsNew.plotCrossection ? &analysisResults.crossectionMB : &analysisResults.spectrumMB, analysisSettingsNew.plotCrossection ? &analysisResults.crossectionHT1 : &analysisResults.spectrumHT1, analysisSettingsNew.plotCrossection ? &analysisResults.crossectionHT2 : &analysisResults.spectrumHT2
		    , spectrumname.Data(), spectrumtitle.Data(), histSpectrum, &legendSpectrum
		    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
		);
	} else {
		TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsMB;
		temp.legendTitle = "STAR d+Au #pi^{0}";
		//histSpectrum = showList(&analysisResults.spectrum
		histSpectrum = showList(analysisSettingsNew.plotCrossection ? &analysisResults.crossection : &analysisResults.spectrum
		    , spectrumname.Data(), spectrumtitle.Data(), histSpectrum, &legendSpectrum
		    , temp
		);
	}
	if (histSpectrum) {
		histSpectrum->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		if (analysisSettingsNew.showWeightUsed) {
			TF1 *funcw = weightUsed.createFunc(true);
			TF1 *funcwCopy = 0;
			if (funcw) {
				funcw->SetLineWidth(1);
				funcw->SetLineColor(11);
				funcw->SetLineStyle(1);
				funcwCopy = funcw->DrawCopy("SAME");
			}
			if (legendSpectrum) legendSpectrum->AddEntry(funcwCopy, "Weight", "L");
		}
	}
	showResultsOthers(histSpectrum, legendSpectrum, analysisSettingsNew.showpQCDWeight, analysisSettingsNew.usePPpQCD);

	TWeightCalculator w("weight_pi0", "Next weight iteration");
	Bool_t wGood = false;
	if (analysisResults.spectrum.size()) {
	    w.Fit(analysisResults.spectrum);
	    w.rangeLow = 0;
	    w.rangeHigh = 100;
	    wGood = true;
	}

	TWeightCalculator wEta("weight_eta", "Next weight iteration - eta");
	Bool_t wEtaGood = false;
	if (analysisResults.spectrumEta.size()) {
	    wEta.Fit(analysisResults.spectrumEta);
	    wEta.rangeLow = 0;
	    wEta.rangeHigh = 100;
	    wEtaGood = true;
	}


	if (analysisSettingsNew.saveWeightNextIteration) {
	    TFile f(analysisSettingsNew.weightFile, "UPDATE");
	    if (f.IsOpen()) {
		cout << "Writing new weights to " << analysisSettingsNew.weightFile << endl;
		if (wGood) w.Write(0, TObject::kOverwrite);
		if (wEtaGood) wEta.Write(0, TObject::kOverwrite);
		if (wgammasGood) wgammas.Write(0, TObject::kOverwrite);
		if (wgammasEtaGood) wgammasEta.Write(0, TObject::kOverwrite);
		if (wpointsMBGood) wpointsMB.Write(0, TObject::kOverwrite);
		if (wpointsHT1Good) wpointsHT1.Write(0, TObject::kOverwrite);
		if (wpointsHT2Good) wpointsHT2.Write(0, TObject::kOverwrite);
	    } else {
		cout << "Cannot open file " << analysisSettingsNew.weightFile << " !" << endl;
	    }
	}

	if (analysisSettingsNew.print || analysisSettingsNew.showWeightNextIteration) {
		if (wGood) {
		    cout << "+--------------------------------------" << endl;
		    w.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
		if (wEtaGood) {
		    cout << "+--------------------------------------" << endl;
		    wEta.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
		if (wgammasGood) {
		    cout << "+--------------------------------------" << endl;
		    wgammas.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
		if (wgammasEtaGood) {
		    cout << "+--------------------------------------" << endl;
		    wgammasEta.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
		if (wpointsMBGood) {
		    cout << "+--------------------------------------" << endl;
		    wpointsMB.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
		if (wpointsHT1Good) {
		    cout << "+--------------------------------------" << endl;
		    wpointsHT1.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
		if (wpointsHT2Good) {
		    cout << "+--------------------------------------" << endl;
		    wpointsHT2.Print("| ");
		    cout << "+--------------------------------------" << endl;
		}
	}

	if (histSpectrum) {
		if (analysisSettingsNew.showWeightNextIteration) {
			TF1 *func = w.createFunc(true);
			TF1 *funcCopy = 0;
			if (func) {
				func->SetLineWidth(1);
				func->SetLineColor(11);
				func->SetLineStyle(2);
				func->SetRange(0, func->GetXmax());
				funcCopy = func->DrawCopy("SAME");
			}
			if (legendSpectrum) legendSpectrum->AddEntry(funcCopy, w.GetTitle(), "L");
		}
		histSpectrum->GetYaxis()->SetLabelOffset(0);
		histSpectrum->GetYaxis()->SetTitleOffset(1.2);
	}


	if (analysisSettingsNew.saveTPCVertexEff) {
	    TFile f(analysisSettingsNew.weightFile, "UPDATE");
	    if (f.IsOpen()) {
		cout << "Writing TPC vertex efficiency to " << analysisSettingsNew.weightFile << endl;
		if (histTPCVertexIneffMB) histTPCVertexIneffMB->Write(0, TObject::kOverwrite);
		if (histTPCVertexIneffHT1) histTPCVertexIneffHT1->Write(0, TObject::kOverwrite);
		if (histTPCVertexIneffHT2) histTPCVertexIneffHT2->Write(0, TObject::kOverwrite);
	    } else {
		cout << "Cannot open file " << analysisSettingsNew.weightFile << " !" << endl;
	    }
	}
	if (histTPCVertexIneffMB) delete histTPCVertexIneffMB;
	if (histTPCVertexIneffHT1) delete histTPCVertexIneffHT1;
	if (histTPCVertexIneffHT2) delete histTPCVertexIneffHT2;


	if (analysisSettingsNew.saveBBCVertexEff) {
	    TFile f(analysisSettingsNew.weightFile, "UPDATE");
	    if (f.IsOpen()) {
		cout << "Writing BBC vertex efficiency to " << analysisSettingsNew.weightFile << endl;
		if (histBBCVertexIneffMB) histBBCVertexIneffMB->Write(0, TObject::kOverwrite);
		if (histBBCVertexIneffHT1) histBBCVertexIneffHT1->Write(0, TObject::kOverwrite);
		if (histBBCVertexIneffHT2) histBBCVertexIneffHT2->Write(0, TObject::kOverwrite);
	    } else {
		cout << "Cannot open file " << analysisSettingsNew.weightFile << " !" << endl;
	    }
	}
	if (histBBCVertexIneffMB) delete histBBCVertexIneffMB;
	if (histBBCVertexIneffHT1) delete histBBCVertexIneffHT1;
	if (histBBCVertexIneffHT2) delete histBBCVertexIneffHT2;

	TString etaToPi0name = analysisSettingsNew.name + "_etatopi";
	TString etaToPi0title = analysisSettingsNew.title + ", #eta / #pi^{0} Ratio;p_{T}, GeV/c;#eta / #pi^{0}";
	TLegend *legendEtaToPi0 = 0;
	TH1F *histEtaToPi0 = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		histEtaToPi0 = showLists(&analysisResults.etaToPi0RatioMB, &analysisResults.etaToPi0RatioHT1, &analysisResults.etaToPi0RatioHT2
		    , etaToPi0name.Data(), etaToPi0title.Data(), histEtaToPi0, &legendEtaToPi0
		    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
		);
	} else {
		TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsMB;
		temp.legendTitle = "STAR d+Au #eta / #pi^{0}";
		histEtaToPi0 = showList(&analysisResults.etaToPi0Ratio
		    , etaToPi0name.Data(), etaToPi0title.Data(), histEtaToPi0, &legendEtaToPi0
		    , temp
		);
	}
	if (histEtaToPi0) {
		histEtaToPi0->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}
	showResultsOthersEtaToPi0(histEtaToPi0, legendEtaToPi0);

	TF1 *weightDiv = new TF1("weightDiv", PQCDFUNC(analysisSettingsNew.usePPpQCD), min(analysisSettingsNew.settingsMB.lowPt, min(analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT2.lowPt)), max(analysisSettingsNew.settingsMB.highPt, max(analysisSettingsNew.settingsHT1.highPt, analysisSettingsNew.settingsHT2.highPt)), 1);
	weightDiv->SetParameter(0, 0);
	divideBinsFunc(analysisResults.spectrumMB, weightDiv, analysisResults.spectrumMBDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT1, weightDiv, analysisResults.spectrumHT1Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT2, weightDiv, analysisResults.spectrumHT2Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrum, weightDiv, analysisResults.spectrumDiv, 1.0, 0.0, false, false, false);

	divideBinsFunc(analysisResults.crossectionMB, weightDiv, analysisResults.crossectionMBDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT1, weightDiv, analysisResults.crossectionHT1Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT2, weightDiv, analysisResults.crossectionHT2Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossection, weightDiv, analysisResults.crossectionDiv, 1.0, 0.0, false, false, false);

	TF1 *weightDivEta = new TF1("weightDivEta", PQCDFUNC(analysisSettingsNew.usePPpQCD), min(analysisSettingsNew.settingsMB.lowPt, min(analysisSettingsNew.settingsHT1.lowPt, analysisSettingsNew.settingsHT2.lowPt)), max(analysisSettingsNew.settingsMB.highPt, max(analysisSettingsNew.settingsHT1.highPt, analysisSettingsNew.settingsHT2.highPt)), 1);
	weightDivEta->SetParameter(0, 0);
	divideBinsFunc(analysisResults.spectrumMBEta, weightDivEta, analysisResults.spectrumMBDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT1Eta, weightDivEta, analysisResults.spectrumHT1DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT2Eta, weightDivEta, analysisResults.spectrumHT2DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumEta, weightDivEta, analysisResults.spectrumDivEta, 1.0, 0.0, false, false, false);

	divideBinsFunc(analysisResults.crossectionMBEta, weightDivEta, analysisResults.crossectionMBDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT1Eta, weightDivEta, analysisResults.crossectionHT1DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT2Eta, weightDivEta, analysisResults.crossectionHT2DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionEta, weightDivEta, analysisResults.crossectionDivEta, 1.0, 0.0, false, false, false);

	TLegend *legendSpectrumDiv = 0;
	TString spectrum_divname = analysisSettingsNew.name + (analysisSettingsNew.plotCrossection ? "_crossection_div" : "_spectrum_div");
	//TString spectrum_divtitle = analysisSettingsNew.title + ", Spectrum divided by pQCD;p_{T}, GeV/c;Data / pQCD";
	TString spectrum_divtitle = analysisSettingsNew.title + (analysisSettingsNew.plotCrossection ? ", Cross section divided by pQCD;p_{T}, GeV/c;Data / pQCD" : ", Spectrum divided by pQCD;p_{T}, GeV/c;Data / pQCD");
	TH1F *histSpectrumDiv = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		//histSpectrumDiv = showLists(&analysisResults.spectrumMBDiv, &analysisResults.spectrumHT1Div, &analysisResults.spectrumHT2Div
		histSpectrumDiv = showLists(analysisSettingsNew.plotCrossection ? &analysisResults.crossectionMBDiv : &analysisResults.spectrumMBDiv, analysisSettingsNew.plotCrossection ? &analysisResults.crossectionHT1Div : &analysisResults.spectrumHT1Div, analysisSettingsNew.plotCrossection ? &analysisResults.crossectionHT2Div : &analysisResults.spectrumHT2Div
		    , spectrum_divname.Data(), spectrum_divtitle.Data()
		    , histSpectrumDiv, &legendSpectrumDiv
		    , analysisSettingsNew.settingsMB, analysisSettingsNew.settingsHT1, analysisSettingsNew.settingsHT2
		);
		//histSpectrumDiv = showLists(&analysisResults.spectrumMBDivEta, &analysisResults.spectrumHT1DivEta, &analysisResults.spectrumHT2DivEta
		/*histSpectrumDiv = showLists(&analysisResults.crossectionMBDivEta, &analysisResults.crossectionHT1DivEta, &analysisResults.crossectionHT2DivEta
		    , spectrum_divname.Data(), spectrum_divtitle.Data()
		    , histSpectrumDiv, &legendSpectrumDiv
		    , analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.settingsSimEtaHT2
		);*/
	} else {
		{
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsMB;
		    temp.legendTitle = "STAR d+Au #pi^{0}";
		    //histSpectrumDiv = showList(&analysisResults.spectrumDiv
		    histSpectrumDiv = showList(analysisSettingsNew.plotCrossection ? &analysisResults.crossectionDiv : &analysisResults.spectrumDiv
			, spectrum_divname.Data(), spectrum_divtitle.Data(), histSpectrumDiv, &legendSpectrumDiv
			, temp
		    );
		}
		/*{
		    TAnalysisSettingsTrigger temp = analysisSettingsNew.settingsSimEtaMB;
		    temp.legendTitle = "STAR d+Au #eta";
		    //histSpectrumDiv = showList(&analysisResults.spectrumDivEta
		    histSpectrumDiv = showList(&analysisResults.crossectionDivEta
			, spectrum_divname.Data(), spectrum_divtitle.Data(), histSpectrumDiv, &legendSpectrumDiv
			, temp
		    );
		}*/
	}
	if (histSpectrumDiv) {
		histSpectrumDiv->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}
	showResultsOthersDiv(histSpectrumDiv, legendSpectrumDiv, weightDiv, analysisSettingsNew.usePPpQCD);
	if (histSpectrumDiv) {
		if (analysisSettingsNew.showWeightNextIteration) {
		
		}
	}

	if (analysisSettingsNew.saveDataArrays) {
		TString name = analysisSettingsNew.name;
		TString title = analysisSettingsNew.title;
		TString filename = Form(analysisSettingsNew.dataArraysFilenameFormat.Data(), name.Data());
		cout << "Saving data arrays to " << filename << endl;
		ofstream ostr(filename.Data());
		ostr << "// " << title << endl;
		TDatime startTime;
		ostr << "// Generated: " << startTime.AsSQLString() << endl;
		TString arrayName;
		arrayName = "_MBspectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumMB, "pi0 pT spectrum - MinBias", arrayName.Data());
		arrayName = "_HT1spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT1, "pi0 pT spectrum - HighTower-1", arrayName.Data());
		arrayName = "_HT2spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT2, "pi0 pT spectrum - HighTower-2", arrayName.Data());
		arrayName = "_MBdiv_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumMBDiv, "pi0 pT spectrum / pQCD - MinBias", arrayName.Data());
		arrayName = "_HT1div_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT1Div, "pi0 pT spectrum / pQCD - HighTower-1", arrayName.Data());
		arrayName = "_HT2div_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT2Div, "pi0 pT spectrum / pQCD - HighTower-2", arrayName.Data());
		arrayName = "_MBcrossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionMB, "pi0 pT cross section - MinBias", arrayName.Data());
		arrayName = "_HT1crossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT1, "pi0 pT cross section - HighTower-1", arrayName.Data());
		arrayName = "_HT2crossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT2, "pi0 pT cross section - HighTower-2", arrayName.Data());
		arrayName = "_MBcrossectiondiv_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionMBDiv, "pi0 pT cross section / pQCD - MinBias", arrayName.Data());
		arrayName = "_HT1crossectiondiv_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT1Div, "pi0 pT cross section / pQCD - HighTower-1", arrayName.Data());
		arrayName = "_HT2crossectiondiv_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT2Div, "pi0 pT cross section / pQCD - HighTower-2", arrayName.Data());
		arrayName = "_MByield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.rawYield, "Raw yield - MinBias", arrayName.Data());
		arrayName = "_HT1yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.rawYield, "Raw yield - HighTower-1", arrayName.Data());
		arrayName = "_HT2yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.rawYield, "Raw yield - HighTower-2", arrayName.Data());
		arrayName = "_MBbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakBackground, "Background - MinBias", arrayName.Data());
		arrayName = "_HT1bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakBackground, "Background - HighTower-1", arrayName.Data());
		arrayName = "_HT2bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakBackground, "Background - HighTower-2", arrayName.Data());
		arrayName = "_MBsigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakSigBackground, "Signal+Background - MinBias", arrayName.Data());
		arrayName = "_HT1sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakSigBackground, "Signal+Background - HighTower-1", arrayName.Data());
		arrayName = "_HT2sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakSigBackground, "Signal+Background - HighTower-2", arrayName.Data());
		arrayName = "_MBeff_"; arrayName += name; outputCArrays(ostr, analysisResults.effMBSmooth, "Correction factor - MinBias", arrayName.Data());
		arrayName = "_HT1eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT1Smooth, "Correction factor - HighTower-1", arrayName.Data());
		arrayName = "_HT2eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT2Smooth, "Correction factor - HighTower-2", arrayName.Data());
		arrayName = "_MBsim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuMB, "Simulated yield - MinBias", arrayName.Data());
		arrayName = "_HT1sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuHT1, "Simulated yield - HighTower-1", arrayName.Data());
		arrayName = "_HT2sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuHT2, "Simulated yield - HighTower-2", arrayName.Data());
		arrayName = "_MBsimreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimMB.rawYield, "Simulated reco yield - MinBias", arrayName.Data());
		arrayName = "_HT1simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimHT1.rawYield, "Simulated reco yield - HighTower-1", arrayName.Data());
		arrayName = "_HT2simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimHT2.rawYield, "Simulated reco yield - HighTower-2", arrayName.Data());

		arrayName = "_EtaMBspectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumMBEta, "eta pT spectrum - MinBias", arrayName.Data());
		arrayName = "_EtaHT1spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT1Eta, "eta pT spectrum - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT2Eta, "eta pT spectrum - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBdiv_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumMBDivEta, "eta pT spectrum / pQCD - MinBias", arrayName.Data());
		arrayName = "_EtaHT1div_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT1DivEta, "eta pT spectrum / pQCD - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2div_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT2DivEta, "eta pT spectrum / pQCD - HighTower-2", arrayName.Data());
		arrayName = "_EtaMByield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.rawYieldEta, "Raw eta yield - MinBias", arrayName.Data());
		arrayName = "_EtaHT1yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.rawYieldEta, "Raw eta yield - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.rawYieldEta, "Raw eta yield - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakBackgroundEta, "Background eta - MinBias", arrayName.Data());
		arrayName = "_EtaHT1bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakBackgroundEta, "Background eta - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakBackgroundEta, "Background eta - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBsigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakSigBackgroundEta, "Signal+Background eta - MinBias", arrayName.Data());
		arrayName = "_EtaHT1sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakSigBackgroundEta, "Signal+Background eta - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakSigBackgroundEta, "Signal+Background eta - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBeff_"; arrayName += name; outputCArrays(ostr, analysisResults.effMBSmoothEta, "Correction factor eta - MinBias", arrayName.Data());
		arrayName = "_EtaHT1eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT1SmoothEta, "Correction factor eta - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT2SmoothEta, "Correction factor eta - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBsim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuEtaMB, "Simulated eta yield - MinBias", arrayName.Data());
		arrayName = "_EtaHT1sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuEtaHT1, "Simulated eta yield - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuEtaHT2, "Simulated eta yield - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBsimreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimEtaMB.rawYield, "Simulated reco eta yield - MinBias", arrayName.Data());
		arrayName = "_EtaHT1simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimEtaHT1.rawYield, "Simulated reco eta yield - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimEtaHT2.rawYield, "Simulated reco eta yield - HighTower-2", arrayName.Data());
		ostr.close();
	}

}

void savePtSlices(const Char_t *filenameOut
    , const list<TInvariantMassDistribution> &invlistMB
    , const list<TInvariantMassDistribution> &invlistHT1
    , const list<TInvariantMassDistribution> &invlistHT2
    , Float_t triggerThresholdMB
    , Float_t triggerThresholdHT1
    , Float_t triggerThresholdHT2
    , TCanvas *peakFinderMBcanvas
    , TCanvas *peakFinderHT1canvas
    , TCanvas *peakFinderHT2canvas
) {
		list<Float_t> centers, widths;
		ofstream ostr(filenameOut);
		ostr << "<HTML>" << endl;
		ostr << "<HEAD><TITLE>Comparison of different triggers in the same pT bin</TITLE></HEAD>" << endl;
		ostr << "<BODY>" << endl;
		ostr << "<H1>Comparison of different triggers in the same p<SUB>T</SUB> bin</H1>" << endl;
		ostr << "<TABLE border=0>\n<TR align=center>";
		ostr << "<TD bgcolor=e0e0e0>p<SUB>T</SUB> range, GeV/c<BR/>relative to the<BR/>trigger threshold T</TD>";
		ostr << "<TD bgcolor=e0e0e0>MinBias<BR/>T = " << triggerThresholdMB << " GeV/c</TD>";
		ostr << "<TD bgcolor=e0e0e0>HighTower-1<BR/>T = " << triggerThresholdHT1 << " GeV/c</TD>";
		ostr << "<TD bgcolor=e0e0e0>HighTower-2<BR/>T = " << triggerThresholdHT2 << " GeV/c</TD>";
		ostr << "</TR>\n" << endl;
		if (peakFinderMBcanvas) peakFinderMBcanvas->SetWindowSize(1024, 768);
		if (peakFinderHT1canvas) peakFinderHT1canvas->SetWindowSize(1024, 768);
		if (peakFinderHT2canvas) peakFinderHT2canvas->SetWindowSize(1024, 768);
		gSystem->ProcessEvents();
			//list<TInvariantMassDistribution>::const_iterator iter = invlistMB.begin();
			for (list<TInvariantMassDistribution>::const_iterator iterMB = invlistMB.begin();iterMB != invlistMB.end();++iterMB) {
				const TInvariantMassDistribution &inv = *iterMB;
				centers.push_back(inv.getBinParameters().getCenter() - triggerThresholdMB);
				widths.push_back(inv.getBinParameters().getWidth());
			}
			for (list<TInvariantMassDistribution>::const_iterator iterHT1 = invlistHT1.begin();iterHT1 != invlistHT1.end();++iterHT1) {
				const TInvariantMassDistribution &inv = *iterHT1;
				centers.push_back(inv.getBinParameters().getCenter() - triggerThresholdHT1);
				widths.push_back(inv.getBinParameters().getWidth());
			}
			for (list<TInvariantMassDistribution>::const_iterator iterHT2 = invlistHT2.begin();iterHT2 != invlistHT2.end();++iterHT2) {
				const TInvariantMassDistribution &inv = *iterHT2;
				centers.push_back(inv.getBinParameters().getCenter() - triggerThresholdHT2);
				widths.push_back(inv.getBinParameters().getWidth());
			}
			centers.sort();
			for (list<Float_t>::const_iterator centerIter = centers.begin(), widthIter = widths.begin();(centerIter != centers.end()) && (widthIter != widths.end());++centerIter, ++widthIter) {
				Float_t center = *centerIter, width = *widthIter;
				Bool_t processedAlready = false;
				for (list<Float_t>::const_iterator centerIter1 = centers.begin(), widthIter1 = widths.begin();(centerIter1 != centerIter) && (widthIter1 != widthIter);++centerIter1, ++widthIter1) {
					Float_t center1 = *centerIter1, width1 = *widthIter1;
					processedAlready |= (center == center1) && (width == width1);
				}
				if (processedAlready) continue;
				TString nameMB, nameHT1, nameHT2;
				Int_t padMB = -1, padHT1 = -1, padHT2 = -1;
				Int_t numPadsFound = 0;
				Int_t padiMB = 1;
				for (list<TInvariantMassDistribution>::const_iterator iter = invlistMB.begin();iter != invlistMB.end();++iter) {
					const TInvariantMassDistribution &inv = *iter;
					if ((center == (inv.getBinParameters().getCenter() - triggerThresholdMB)) && (width == inv.getBinParameters().getWidth())) {
						nameMB = "MB_pTbin"; nameMB += padiMB; nameMB += ".gif";
						padMB = padiMB;
						numPadsFound++;
					}
					padiMB++;
				}
				Int_t padiHT1 = 1;
				for (list<TInvariantMassDistribution>::const_iterator iter = invlistHT1.begin();iter != invlistHT1.end();++iter) {
					const TInvariantMassDistribution &inv = *iter;
					if ((center == (inv.getBinParameters().getCenter() - triggerThresholdHT1)) && (width == inv.getBinParameters().getWidth())) {
						nameHT1 = "HT1_pTbin"; nameHT1 += padiHT1; nameHT1 += ".gif";
						padHT1 = padiHT1;
						numPadsFound++;
					}
					padiHT1++;
				}
				Int_t padiHT2 = 1;
				for (list<TInvariantMassDistribution>::const_iterator iter = invlistHT2.begin();iter != invlistHT2.end();++iter) {
					const TInvariantMassDistribution &inv = *iter;
					if ((center == (inv.getBinParameters().getCenter() - triggerThresholdHT2)) && (width == inv.getBinParameters().getWidth())) {
						nameHT2 = "HT2_pTbin"; nameHT2 += padiHT2; nameHT2 += ".gif";
						padHT2 = padiHT2;
						numPadsFound++;
					}
					padiHT2++;
				}
				if (numPadsFound > 1) {
					if ((padMB != -1) && (peakFinderMBcanvas != 0)) {
						peakFinderMBcanvas->cd(padMB);
						TVirtualPad *pad = peakFinderMBcanvas->GetSelectedPad();
						TImage *img = TImage::Create();
						if (pad && img) {img->FromPad(pad); img->Scale(250, 155); img->WriteImage(nameMB.Data());}
						if (img) delete img;
					}
					if ((padHT1 != -1) && (peakFinderHT1canvas != 0)) {
						peakFinderHT1canvas->cd(padHT1);
						TVirtualPad *pad = peakFinderHT1canvas->GetSelectedPad();
						TImage *img = TImage::Create();
						if (pad && img) {img->FromPad(pad); img->Scale(250, 155); img->WriteImage(nameHT1.Data());}
						if (img) delete img;
					}
					if ((padHT2 != -1) && (peakFinderHT2canvas != 0)) {
						peakFinderHT2canvas->cd(padHT2);
						TVirtualPad *pad = peakFinderHT2canvas->GetSelectedPad();
						TImage *img = TImage::Create();
						if (pad && img) {img->FromPad(pad); img->Scale(250, 155); img->WriteImage(nameHT2.Data());}
						if (img) delete img;
					}
					ostr << "<TR align=center>";
					ostr << "<TD bgcolor=e0e0e0>T " << ((center - width/2.0) >= 0 ? "+" : "-") << " " << fabs(center - width/2.0) << " &lt; p<SUB>T</SUB> &lt; T " << ((center + width/2.0) >= 0 ? "+" : "-") << " " << fabs(center + width/2.0) << "</TD>";
					ostr << "<TD>";
					if (padMB != -1) ostr << "<IMG src=\"" << nameMB << "\"/>";
					ostr << "</TD>";
					ostr << "<TD>";
					if (padHT1 != -1) ostr << "<IMG src=\"" << nameHT1 << "\"/>";
					ostr << "</TD>";
					ostr << "<TD>";
					if (padHT2 != -1) ostr << "<IMG src=\"" << nameHT2 << "\"/>";
					ostr << "</TD>";
					ostr << "</TR>";
					ostr << endl;
				}
			}
		ostr << "</TABLE></BODY></HTML>" << endl;
		ostr.close();
}

void showAnalysis(const Char_t *DATA_DIR, const TAllSettings &allSettings, TAllResults &allResults) {
	cout << "================== STARTED ===================" << endl;
	TDatime startTime;
	TStopwatch timer;
	cout << "Started: " << startTime.AsSQLString() << endl;
	timer.Start();

        TFolder *rootFolder = gROOT->GetRootFolder();
        TFolder *dataFolder = (rootFolder && (allSettings.browseAnalysis || allSettings.browseResults)) ? rootFolder->AddFolder("processedData", "Processed data") : 0;
        TDataProcessorMemberInspector memberInspector;

#define FILE_NAME(NAME, PATH1, PATH2) TString NAME##str = PATH1; NAME##str += PATH2; const Char_t *NAME = NAME##str.Data();

#define DATA_FILE_DAU(NAME, CENTRALITY) \
	FILE_NAME(filenameDataMB_0##NAME, DATA_DIR, "/data_dAu2003_MB_dAuMinBias_" CENTRALITY ".root") \
	FILE_NAME(filenameDataMB_1##NAME, DATA_DIR, "/data_dAu2003_MB_dAuCombined_" CENTRALITY ".root") \
	FILE_NAME(filenameDataMB_2##NAME, DATA_DIR, "/data_dAu2003_MB_UPCCombined_" CENTRALITY ".root") \
	/*FILE_NAME(filenameDataMB_3##NAME, DATA_DIR, "/data_dAu2003_MB_" CENTRALITY ".root")*/ \
	FILE_NAME(filenameDataHT1_0##NAME, DATA_DIR, "/data_dAu2003_HT1_dAuCombined_" CENTRALITY ".root") \
	FILE_NAME(filenameDataHT1_1##NAME, DATA_DIR, "/data_dAu2003_HT1_UPCCombined_" CENTRALITY ".root") \
	/*FILE_NAME(filenameDataHT1_2##NAME, DATA_DIR, "/data_dAu2003_HT1_" CENTRALITY ".root")*/ \
	FILE_NAME(filenameDataHT2_0##NAME, DATA_DIR, "/data_dAu2003_HT2_dAuCombined_" CENTRALITY ".root") \
	FILE_NAME(filenameDataHT2_1##NAME, DATA_DIR, "/data_dAu2003_HT2_UPCCombined_" CENTRALITY ".root") \
	/*FILE_NAME(filenameDataHT2_2##NAME, DATA_DIR, "/data_dAu2003_HT2_" CENTRALITY ".root")*/ \
	FILE_NAME(filenameEmbed##NAME, DATA_DIR, "/embedding_dAu2003_" CENTRALITY ".root") \
	FILE_NAME(filenameEmbedNoweight##NAME, DATA_DIR, "/embedding_dAu2003_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameEmbedEta##NAME, DATA_DIR, "/embedding_dAu2003_eta_" CENTRALITY ".root") \
	FILE_NAME(filenameEmbedNoweightEta##NAME, DATA_DIR, "/embedding_dAu2003_eta_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuMB##NAME, DATA_DIR, "/simulation_dAu2003_pi0_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuHT1##NAME, DATA_DIR, "/simulation_dAu2003_pi0_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuHT2##NAME, DATA_DIR, "/simulation_dAu2003_pi0_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightMB##NAME, DATA_DIR, "/simulation_dAu2003_pi0_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightHT1##NAME, DATA_DIR, "/simulation_dAu2003_pi0_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightHT2##NAME, DATA_DIR, "/simulation_dAu2003_pi0_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuEtaMB##NAME, DATA_DIR, "/simulation_dAu2003_eta_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuEtaHT1##NAME, DATA_DIR, "/simulation_dAu2003_eta_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuEtaHT2##NAME, DATA_DIR, "/simulation_dAu2003_eta_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightEtaMB##NAME, DATA_DIR, "/simulation_dAu2003_eta_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightEtaHT1##NAME, DATA_DIR, "/simulation_dAu2003_eta_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightEtaHT2##NAME, DATA_DIR, "/simulation_dAu2003_eta_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuEtabgMB##NAME, DATA_DIR, "/simulation_dAu2003_etabg_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuEtabgHT1##NAME, DATA_DIR, "/simulation_dAu2003_etabg_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuEtabgHT2##NAME, DATA_DIR, "/simulation_dAu2003_etabg_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightEtabgMB##NAME, DATA_DIR, "/simulation_dAu2003_etabg_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightEtabgHT1##NAME, DATA_DIR, "/simulation_dAu2003_etabg_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimuNoweightEtabgHT2##NAME, DATA_DIR, "/simulation_dAu2003_etabg_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaMB##NAME, DATA_DIR, "/simulation_dAu2003_gamma_data_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaHT1##NAME, DATA_DIR, "/simulation_dAu2003_gamma_data_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaHT2##NAME, DATA_DIR, "/simulation_dAu2003_gamma_data_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightMB##NAME, DATA_DIR, "/simulation_dAu2003_gamma_data_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightHT1##NAME, DATA_DIR, "/simulation_dAu2003_gamma_data_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightHT2##NAME, DATA_DIR, "/simulation_dAu2003_gamma_data_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaSimMB##NAME, DATA_DIR, "/simulation_dAu2003_gamma_pi0_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaSimHT1##NAME, DATA_DIR, "/simulation_dAu2003_gamma_pi0_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaSimHT2##NAME, DATA_DIR, "/simulation_dAu2003_gamma_pi0_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightSimMB##NAME, DATA_DIR, "/simulation_dAu2003_gamma_pi0_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightSimHT1##NAME, DATA_DIR, "/simulation_dAu2003_gamma_pi0_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightSimHT2##NAME, DATA_DIR, "/simulation_dAu2003_gamma_pi0_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaSimEtaMB##NAME, DATA_DIR, "/simulation_dAu2003_gamma_eta_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaSimEtaHT1##NAME, DATA_DIR, "/simulation_dAu2003_gamma_eta_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaSimEtaHT2##NAME, DATA_DIR, "/simulation_dAu2003_gamma_eta_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightSimEtaMB##NAME, DATA_DIR, "/simulation_dAu2003_gamma_eta_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightSimEtaHT1##NAME, DATA_DIR, "/simulation_dAu2003_gamma_eta_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimu1gammaNoweightSimEtaHT2##NAME, DATA_DIR, "/simulation_dAu2003_gamma_eta_HT2_noweight_" CENTRALITY ".root") \
	TDataProcessorPool *poolMB##NAME = new TDataProcessorPool("data_dAu2003_MB_" CENTRALITY, "Real data dAu2003 - MB - " CENTRALITY); \
	TDataProcessorPool *poolHT1##NAME = new TDataProcessorPool("data_dAu2003_HT1_" CENTRALITY, "Real data dAu2003 - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolHT2##NAME = new TDataProcessorPool("data_dAu2003_HT2_" CENTRALITY, "Real data dAu2003 - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolEmbed0##NAME = new TDataProcessorPool("embedding_dAu2003_" CENTRALITY, "Embedding - dAu2003 - " CENTRALITY); \
	TDataProcessorPool *poolEmbedMB##NAME = poolEmbed0##NAME; \
	TDataProcessorPool *poolEmbedHT1##NAME = poolEmbed0##NAME; \
	TDataProcessorPool *poolEmbedHT2##NAME = poolEmbed0##NAME; \
	TDataProcessorPool *poolEmbedNoweight0##NAME = new TDataProcessorPool("embedding_dAu2003_noweight_" CENTRALITY, "Embedding - dAu2003 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolEmbedNoweightMB##NAME = poolEmbedNoweight0##NAME; \
	TDataProcessorPool *poolEmbedNoweightHT1##NAME = poolEmbedNoweight0##NAME; \
	TDataProcessorPool *poolEmbedNoweightHT2##NAME = poolEmbedNoweight0##NAME; \
	TDataProcessorPool *poolEmbedEta0##NAME = new TDataProcessorPool("embedding_dAu2003_eta_" CENTRALITY, "Embedding eta - dAu2003 - " CENTRALITY); \
	TDataProcessorPool *poolEmbedEtaMB##NAME = poolEmbedEta0##NAME; \
	TDataProcessorPool *poolEmbedEtaHT1##NAME = poolEmbedEta0##NAME; \
	TDataProcessorPool *poolEmbedEtaHT2##NAME = poolEmbedEta0##NAME; \
	TDataProcessorPool *poolEmbedNoweightEta0##NAME = new TDataProcessorPool("embedding_dAu2003_eta_noweight_" CENTRALITY, "Embedding eta - dAu2003 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolEmbedNoweightEtaMB##NAME = poolEmbedNoweightEta0##NAME; \
	TDataProcessorPool *poolEmbedNoweightEtaHT1##NAME = poolEmbedNoweightEta0##NAME; \
	TDataProcessorPool *poolEmbedNoweightEtaHT2##NAME = poolEmbedNoweightEta0##NAME; \
	TDataProcessorPool *poolSimuMB##NAME = new TDataProcessorPool("simulation_dAu2003_pi0_MB_" CENTRALITY, "Simulation - dAu2003 - pi0 - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimuHT1##NAME = new TDataProcessorPool("simulation_dAu2003_pi0_HT1_" CENTRALITY, "Simulation - dAu2003 - pi0 - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimuHT2##NAME = new TDataProcessorPool("simulation_dAu2003_pi0_HT2_" CENTRALITY, "Simulation - dAu2003 - pi0 - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightMB##NAME = new TDataProcessorPool("simulation_dAu2003_pi0_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - pi0 - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightHT1##NAME = new TDataProcessorPool("simulation_dAu2003_pi0_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - pi0 - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightHT2##NAME = new TDataProcessorPool("simulation_dAu2003_pi0_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - pi0 - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuEtaMB##NAME = new TDataProcessorPool("simulation_dAu2003_eta_MB_" CENTRALITY, "Simulation - dAu2003 - eta - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimuEtaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_eta_HT1_" CENTRALITY, "Simulation - dAu2003 - eta - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimuEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_eta_HT2_" CENTRALITY, "Simulation - dAu2003 - eta - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightEtaMB##NAME = new TDataProcessorPool("simulation_dAu2003_eta_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - eta - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightEtaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_eta_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - eta - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_eta_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - eta - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuEtabgMB##NAME = new TDataProcessorPool("simulation_dAu2003_etabg_MB_" CENTRALITY, "Simulation - dAu2003 - etabg - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimuEtabgHT1##NAME = new TDataProcessorPool("simulation_dAu2003_etabg_HT1_" CENTRALITY, "Simulation - dAu2003 - etabg - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimuEtabgHT2##NAME = new TDataProcessorPool("simulation_dAu2003_etabg_HT2_" CENTRALITY, "Simulation - dAu2003 - etabg - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightEtabgMB##NAME = new TDataProcessorPool("simulation_dAu2003_etabg_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - etabg - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightEtabgHT1##NAME = new TDataProcessorPool("simulation_dAu2003_etabg_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - etabg - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimuNoweightEtabgHT2##NAME = new TDataProcessorPool("simulation_dAu2003_etabg_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - etabg - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaMB##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_MB_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_HT1_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_HT2_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightMB##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightHT1##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaSimMB##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_pi0_MB_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaSimHT1##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_pi0_HT1_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaSimHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_pi0_HT2_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightSimMB##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_pi0_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightSimHT1##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_pi0_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightSimHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_pi0_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaSimEtaMB##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_MB_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaSimEtaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_HT1_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaSimEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_HT2_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightSimEtaMB##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightSimEtaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimu1gammaNoweightSimEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT2 - no weight - " CENTRALITY);

#define DATA_FILE_LOAD_DAU(NAME, CENTRALITY, READ, SETTINGS, RESULTS) \
	{ \
	if (READ) { \
		calculatePrescaleTriggerBias(0, SETTINGS, RESULTS, 0); \
		poolMB##NAME->processFile(filenameDataMB_0##NAME); \
		calculatePrescaleTriggerBias(poolMB##NAME, SETTINGS, RESULTS, 1); \
		poolMB##NAME->processFile(filenameDataMB_1##NAME); \
		poolMB##NAME->processFile(filenameDataMB_2##NAME); \
		/*poolMB##NAME->processFile(filenameDataMB_3##NAME);*/ \
		calculatePrescaleTriggerBias(poolMB##NAME, SETTINGS, RESULTS, 2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolMB##NAME, dataFolder); \
		poolHT1##NAME->processFile(filenameDataHT1_0##NAME); \
		poolHT1##NAME->processFile(filenameDataHT1_1##NAME); \
		/*poolHT1##NAME->processFile(filenameDataHT1_2##NAME);*/ \
		calculatePrescaleTriggerBias(poolHT1##NAME, SETTINGS, RESULTS, 3); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolHT1##NAME, dataFolder); \
		poolHT2##NAME->processFile(filenameDataHT2_0##NAME); \
		poolHT2##NAME->processFile(filenameDataHT2_1##NAME); \
		/*poolHT2##NAME->processFile(filenameDataHT2_2##NAME);*/ \
		calculatePrescaleTriggerBias(poolHT2##NAME, SETTINGS, RESULTS, 4); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolHT2##NAME, dataFolder); \
	} \
	if (allSettings.useDAuEff && ((allSettings.useDAuEffCentralities && READ) || (!allSettings.useDAuEffCentralities && (CENTRALITY == "nocentral")))) { \
		if (allSettings.useDAuEmbedding) { \
			if (allSettings.useDAuEffWeight || SETTINGS.correctForPtShift || SETTINGS.correctForPtShiftEff) { \
			    poolEmbed0##NAME->processFile(filenameEmbed##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbed0##NAME, dataFolder); \
    			} \
			if (!allSettings.useDAuEffWeight) { \
			    poolEmbedNoweight0##NAME->processFile(filenameEmbedNoweight##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbedNoweight0##NAME, dataFolder); \
			} \
		} else { \
			if (allSettings.useDAuEffWeight || SETTINGS.correctForPtShift || SETTINGS.correctForPtShiftEff) { \
				poolSimuMB##NAME->processFile(filenameSimuMB##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuMB##NAME, dataFolder); \
				poolSimuHT1##NAME->processFile(filenameSimuHT1##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuHT1##NAME, dataFolder); \
				poolSimuHT2##NAME->processFile(filenameSimuHT2##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuHT2##NAME, dataFolder); \
			} \
			if (!allSettings.useDAuEffWeight) { \
				poolSimuNoweightMB##NAME->processFile(filenameSimuNoweightMB##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightMB##NAME, dataFolder); \
				poolSimuNoweightHT1##NAME->processFile(filenameSimuNoweightHT1##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightHT1##NAME, dataFolder); \
				poolSimuNoweightHT2##NAME->processFile(filenameSimuNoweightHT2##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightHT2##NAME, dataFolder); \
			} \
		} \
	} \
	if (allSettings.useDAuEffEta && ((allSettings.useDAuEffCentralitiesEta && READ) || (!allSettings.useDAuEffCentralitiesEta && (CENTRALITY == "nocentral")))) { \
		if (allSettings.useDAuEmbeddingEta) { \
			if (allSettings.useDAuEffWeightEta || SETTINGS.correctForPtShift || SETTINGS.correctForPtShiftEff) { \
			    poolEmbedEta0##NAME->processFile(filenameEmbedEta##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbedEta0##NAME, dataFolder); \
    			} \
			if (!allSettings.useDAuEffWeightEta) { \
			    poolEmbedNoweightEta0##NAME->processFile(filenameEmbedNoweightEta##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbedNoweightEta0##NAME, dataFolder); \
			} \
		} else { \
			if (allSettings.useDAuEffWeightEta || SETTINGS.correctForPtShift || SETTINGS.correctForPtShiftEff) { \
				poolSimuEtaMB##NAME->processFile(filenameSimuEtaMB##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuEtaMB##NAME, dataFolder); \
				poolSimuEtaHT1##NAME->processFile(filenameSimuEtaHT1##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuEtaHT1##NAME, dataFolder); \
				poolSimuEtaHT2##NAME->processFile(filenameSimuEtaHT2##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuEtaHT2##NAME, dataFolder); \
				if (allSettings.useDAuEffEtabg) { \
				    poolSimuEtabgMB##NAME->processFile(filenameSimuEtabgMB##NAME); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuEtabgMB##NAME, dataFolder); \
				    poolSimuEtabgHT1##NAME->processFile(filenameSimuEtabgHT1##NAME); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuEtabgHT1##NAME, dataFolder); \
				    poolSimuEtabgHT2##NAME->processFile(filenameSimuEtabgHT2##NAME); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuEtabgHT2##NAME, dataFolder); \
				} \
			} \
			if (!allSettings.useDAuEffWeightEta) { \
				poolSimuNoweightEtaMB##NAME->processFile(filenameSimuNoweightEtaMB##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightEtaMB##NAME, dataFolder); \
				poolSimuNoweightEtaHT1##NAME->processFile(filenameSimuNoweightEtaHT1##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightEtaHT1##NAME, dataFolder); \
				poolSimuNoweightEtaHT2##NAME->processFile(filenameSimuNoweightEtaHT2##NAME); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightEtaHT2##NAME, dataFolder); \
				if (allSettings.useDAuEffEtabg) { \
				    poolSimuNoweightEtabgMB##NAME->processFile(filenameSimuNoweightEtabgMB##NAME); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightEtabgMB##NAME, dataFolder); \
				    poolSimuNoweightEtabgHT1##NAME->processFile(filenameSimuNoweightEtabgHT1##NAME); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightEtabgHT1##NAME, dataFolder); \
				    poolSimuNoweightEtabgHT2##NAME->processFile(filenameSimuNoweightEtabgHT2##NAME); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimuNoweightEtabgHT2##NAME, dataFolder); \
				}\
			} \
		} \
	} \
	if (allSettings.useDAu1gamma && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (CENTRALITY == "nocentral")))) { \
		poolSimu1gammaMB##NAME->processFile(filenameSimu1gammaMB##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaMB##NAME, dataFolder); \
		poolSimu1gammaHT1##NAME->processFile(filenameSimu1gammaHT1##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaHT1##NAME, dataFolder); \
		poolSimu1gammaHT2##NAME->processFile(filenameSimu1gammaHT2##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaHT2##NAME, dataFolder); \
		if (!allSettings.useDAuEffWeight) { \
		    poolSimu1gammaNoweightMB##NAME->processFile(filenameSimu1gammaNoweightMB##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightMB##NAME, dataFolder); \
		    poolSimu1gammaNoweightHT1##NAME->processFile(filenameSimu1gammaNoweightHT1##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightHT1##NAME, dataFolder); \
		    poolSimu1gammaNoweightHT2##NAME->processFile(filenameSimu1gammaNoweightHT2##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightHT2##NAME, dataFolder); \
		} \
	} \
	if (allSettings.useDAu1gammaSim && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (CENTRALITY == "nocentral")))) { \
		poolSimu1gammaSimMB##NAME->processFile(filenameSimu1gammaSimMB##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaSimMB##NAME, dataFolder); \
		poolSimu1gammaSimHT1##NAME->processFile(filenameSimu1gammaSimHT1##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaSimHT1##NAME, dataFolder); \
		poolSimu1gammaSimHT2##NAME->processFile(filenameSimu1gammaSimHT2##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaSimHT2##NAME, dataFolder); \
		if (!allSettings.useDAuEffWeight) { \
		    poolSimu1gammaNoweightSimMB##NAME->processFile(filenameSimu1gammaNoweightSimMB##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightSimMB##NAME, dataFolder); \
		    poolSimu1gammaNoweightSimHT1##NAME->processFile(filenameSimu1gammaNoweightSimHT1##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightSimHT1##NAME, dataFolder); \
		    poolSimu1gammaNoweightSimHT2##NAME->processFile(filenameSimu1gammaNoweightSimHT2##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightSimHT2##NAME, dataFolder); \
		} \
	} \
	if (allSettings.useDAu1gammaSimEta && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (CENTRALITY == "nocentral")))) { \
		poolSimu1gammaSimEtaMB##NAME->processFile(filenameSimu1gammaSimEtaMB##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaSimEtaMB##NAME, dataFolder); \
		poolSimu1gammaSimEtaHT1##NAME->processFile(filenameSimu1gammaSimEtaHT1##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaSimEtaHT1##NAME, dataFolder); \
		poolSimu1gammaSimEtaHT2##NAME->processFile(filenameSimu1gammaSimEtaHT2##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaSimEtaHT2##NAME, dataFolder); \
		if (!allSettings.useDAuEffWeight) { \
		    poolSimu1gammaNoweightSimEtaMB##NAME->processFile(filenameSimu1gammaNoweightSimEtaMB##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightSimEtaMB##NAME, dataFolder); \
		    poolSimu1gammaNoweightSimEtaHT1##NAME->processFile(filenameSimu1gammaNoweightSimEtaHT1##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightSimEtaHT1##NAME, dataFolder); \
		    poolSimu1gammaNoweightSimEtaHT2##NAME->processFile(filenameSimu1gammaNoweightSimEtaHT2##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimu1gammaNoweightSimEtaHT2##NAME, dataFolder); \
		} \
	} \
	}

#define DATA_FILE_DELETE_DAU(NAME) \
	{ \
	removePool(poolMB##NAME, true, true); \
	removePool(poolHT1##NAME, true, true); \
	removePool(poolHT2##NAME, true, true); \
	removePool(poolEmbed0##NAME, true, true); \
	removePool(poolEmbedMB##NAME, true, true); \
	removePool(poolEmbedHT1##NAME, true, true); \
	removePool(poolEmbedHT2##NAME, true, true); \
	removePool(poolEmbedNoweight0##NAME, true, true); \
	removePool(poolEmbedNoweightMB##NAME, true, true); \
	removePool(poolEmbedNoweightHT1##NAME, true, true); \
	removePool(poolEmbedNoweightHT2##NAME, true, true); \
	removePool(poolEmbedEta0##NAME, true, true); \
	removePool(poolEmbedEtaMB##NAME, true, true); \
	removePool(poolEmbedEtaHT1##NAME, true, true); \
	removePool(poolEmbedEtaHT2##NAME, true, true); \
	removePool(poolEmbedNoweightEta0##NAME, true, true); \
	removePool(poolEmbedNoweightEtaMB##NAME, true, true); \
	removePool(poolEmbedNoweightEtaHT1##NAME, true, true); \
	removePool(poolEmbedNoweightEtaHT2##NAME, true, true); \
	removePool(poolSimuMB##NAME, true, true); \
	removePool(poolSimuHT1##NAME, true, true); \
	removePool(poolSimuHT2##NAME, true, true); \
	removePool(poolSimuNoweightMB##NAME, true, true); \
	removePool(poolSimuNoweightHT1##NAME, true, true); \
	removePool(poolSimuNoweightHT2##NAME, true, true); \
	removePool(poolSimuEtaMB##NAME, true, true); \
	removePool(poolSimuEtaHT1##NAME, true, true); \
	removePool(poolSimuEtaHT2##NAME, true, true); \
	removePool(poolSimuNoweightEtaMB##NAME, true, true); \
	removePool(poolSimuNoweightEtaHT1##NAME, true, true); \
	removePool(poolSimuNoweightEtaHT2##NAME, true, true); \
	removePool(poolSimuEtabgMB##NAME, true, true); \
	removePool(poolSimuEtabgHT1##NAME, true, true); \
	removePool(poolSimuEtabgHT2##NAME, true, true); \
	removePool(poolSimuNoweightEtabgMB##NAME, true, true); \
	removePool(poolSimuNoweightEtabgHT1##NAME, true, true); \
	removePool(poolSimuNoweightEtabgHT2##NAME, true, true); \
	removePool(poolSimu1gammaMB##NAME, true, true); \
	removePool(poolSimu1gammaHT1##NAME, true, true); \
	removePool(poolSimu1gammaHT2##NAME, true, true); \
	removePool(poolSimu1gammaNoweightMB##NAME, true, true); \
	removePool(poolSimu1gammaNoweightHT1##NAME, true, true); \
	removePool(poolSimu1gammaNoweightHT2##NAME, true, true); \
	removePool(poolSimu1gammaSimMB##NAME, true, true); \
	removePool(poolSimu1gammaSimHT1##NAME, true, true); \
	removePool(poolSimu1gammaSimHT2##NAME, true, true); \
	removePool(poolSimu1gammaNoweightSimMB##NAME, true, true); \
	removePool(poolSimu1gammaNoweightSimHT1##NAME, true, true); \
	removePool(poolSimu1gammaNoweightSimHT2##NAME, true, true); \
	removePool(poolSimu1gammaSimEtaMB##NAME, true, true); \
	removePool(poolSimu1gammaSimEtaHT1##NAME, true, true); \
	removePool(poolSimu1gammaSimEtaHT2##NAME, true, true); \
	removePool(poolSimu1gammaNoweightSimEtaMB##NAME, true, true); \
	removePool(poolSimu1gammaNoweightSimEtaHT1##NAME, true, true); \
	removePool(poolSimu1gammaNoweightSimEtaHT2##NAME, true, true); \
	}

#define DATA_FILE_PP \
	FILE_NAME(filenamePPDataMB_0, DATA_DIR, "/data_pp2005_MB_ppProductionMinBias.root") \
	FILE_NAME(filenamePPDataMB_1, DATA_DIR, "/data_pp2005_MB_ppProduction.root") \
	/*FILE_NAME(filenamePPDataMB_2, DATA_DIR, "/data_pp2005_MB.root")*/ \
	FILE_NAME(filenamePPDataHT1_0, DATA_DIR, "/data_pp2005_HT1_ppProduction.root") \
	/*FILE_NAME(filenamePPDataHT1_1, DATA_DIR, "/data_pp2005_HT1.root")*/ \
	FILE_NAME(filenamePPDataHT2_0, DATA_DIR, "/data_pp2005_HT2_ppProduction.root") \
	/*FILE_NAME(filenamePPDataHT2_1, DATA_DIR, "/data_pp2005_HT2.root")*/ \
	FILE_NAME(filenamePPSimuMB, DATA_DIR, "/simulation_pp2005_pi0_MB.root") \
	FILE_NAME(filenamePPSimuHT1, DATA_DIR, "/simulation_pp2005_pi0_HT1.root") \
	FILE_NAME(filenamePPSimuHT2, DATA_DIR, "/simulation_pp2005_pi0_HT2.root") \
	FILE_NAME(filenamePPSimuNoweightMB, DATA_DIR, "/simulation_pp2005_pi0_MB_noweight.root") \
	FILE_NAME(filenamePPSimuNoweightHT1, DATA_DIR, "/simulation_pp2005_pi0_HT1_noweight.root") \
	FILE_NAME(filenamePPSimuNoweightHT2, DATA_DIR, "/simulation_pp2005_pi0_HT2_noweight.root") \
	FILE_NAME(filenamePPSimuEtaMB, DATA_DIR, "/simulation_pp2005_eta_MB.root") \
	FILE_NAME(filenamePPSimuEtaHT1, DATA_DIR, "/simulation_pp2005_eta_HT1.root") \
	FILE_NAME(filenamePPSimuEtaHT2, DATA_DIR, "/simulation_pp2005_eta_HT2.root") \
	FILE_NAME(filenamePPSimuNoweightEtaMB, DATA_DIR, "/simulation_pp2005_eta_MB_noweight.root") \
	FILE_NAME(filenamePPSimuNoweightEtaHT1, DATA_DIR, "/simulation_pp2005_eta_HT1_noweight.root") \
	FILE_NAME(filenamePPSimuNoweightEtaHT2, DATA_DIR, "/simulation_pp2005_eta_HT2_noweight.root") \
	FILE_NAME(filenamePPSimuEtabgMB, DATA_DIR, "/simulation_pp2005_etabg_MB.root") \
	FILE_NAME(filenamePPSimuEtabgHT1, DATA_DIR, "/simulation_pp2005_etabg_HT1.root") \
	FILE_NAME(filenamePPSimuEtabgHT2, DATA_DIR, "/simulation_pp2005_etabg_HT2.root") \
	FILE_NAME(filenamePPSimuNoweightEtabgMB, DATA_DIR, "/simulation_pp2005_etabg_MB_noweight.root") \
	FILE_NAME(filenamePPSimuNoweightEtabgHT1, DATA_DIR, "/simulation_pp2005_etabg_HT1_noweight.root") \
	FILE_NAME(filenamePPSimuNoweightEtabgHT2, DATA_DIR, "/simulation_pp2005_etabg_HT2_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaMB, DATA_DIR, "/simulation_pp2005_gamma_data_MB.root") \
	FILE_NAME(filenamePPSimu1gammaHT1, DATA_DIR, "/simulation_pp2005_gamma_data_HT1.root") \
	FILE_NAME(filenamePPSimu1gammaHT2, DATA_DIR, "/simulation_pp2005_gamma_data_HT2.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightMB, DATA_DIR, "/simulation_pp2005_gamma_data_MB_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightHT1, DATA_DIR, "/simulation_pp2005_gamma_data_HT1_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightHT2, DATA_DIR, "/simulation_pp2005_gamma_data_HT2_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaSimMB, DATA_DIR, "/simulation_pp2005_gamma_pi0_MB.root") \
	FILE_NAME(filenamePPSimu1gammaSimHT1, DATA_DIR, "/simulation_pp2005_gamma_pi0_HT1.root") \
	FILE_NAME(filenamePPSimu1gammaSimHT2, DATA_DIR, "/simulation_pp2005_gamma_pi0_HT2.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightSimMB, DATA_DIR, "/simulation_pp2005_gamma_pi0_MB_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightSimHT1, DATA_DIR, "/simulation_pp2005_gamma_pi0_HT1_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightSimHT2, DATA_DIR, "/simulation_pp2005_gamma_pi0_HT2_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaSimEtaMB, DATA_DIR, "/simulation_pp2005_gamma_eta_MB.root") \
	FILE_NAME(filenamePPSimu1gammaSimEtaHT1, DATA_DIR, "/simulation_pp2005_gamma_eta_HT1.root") \
	FILE_NAME(filenamePPSimu1gammaSimEtaHT2, DATA_DIR, "/simulation_pp2005_gamma_eta_HT2.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightSimEtaMB, DATA_DIR, "/simulation_pp2005_gamma_eta_MB_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightSimEtaHT1, DATA_DIR, "/simulation_pp2005_gamma_eta_HT1_noweight.root") \
	FILE_NAME(filenamePPSimu1gammaNoweightSimEtaHT2, DATA_DIR, "/simulation_pp2005_gamma_eta_HT2_noweight.root") \
	TDataProcessorPool *poolPPMB = new TDataProcessorPool("data_pp2005_MB", "Real data - pp2005 - MB"); \
	TDataProcessorPool *poolPPHT1 = new TDataProcessorPool("data_pp2005_HT1", "Real data - pp2005 - HT1"); \
	TDataProcessorPool *poolPPHT2 = new TDataProcessorPool("data_pp2005_HT2", "Real data - pp2005 - HT2"); \
	TDataProcessorPool *poolPPSimuMB = new TDataProcessorPool("simulation_pp2005_pi0_MB", "Simulation - pp2005 - pi0 - MB"); \
	TDataProcessorPool *poolPPSimuHT1 = new TDataProcessorPool("simulation_pp2005_pi0_HT1", "Simulation - pp2005 - pi0 - HT1"); \
	TDataProcessorPool *poolPPSimuHT2 = new TDataProcessorPool("simulation_pp2005_pi0_HT2", "Simulation - pp2005 - pi0 - HT2"); \
	TDataProcessorPool *poolPPSimuNoweightMB = new TDataProcessorPool("simulation_pp2005_pi0_MB_noweight", "Simulation - pp2005 - pi0 - MB - no weight"); \
	TDataProcessorPool *poolPPSimuNoweightHT1 = new TDataProcessorPool("simulation_pp2005_pi0_HT1_noweight", "Simulation - pp2005 - pi0 - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimuNoweightHT2 = new TDataProcessorPool("simulation_pp2005_pi0_HT2_noweight", "Simulation - pp2005 - pi0 - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimuEtaMB = new TDataProcessorPool("simulation_pp2005_eta_MB", "Simulation - pp2005 - eta - MB"); \
	TDataProcessorPool *poolPPSimuEtaHT1 = new TDataProcessorPool("simulation_pp2005_eta_HT1", "Simulation - pp2005 - eta - HT1"); \
	TDataProcessorPool *poolPPSimuEtaHT2 = new TDataProcessorPool("simulation_pp2005_eta_HT2", "Simulation - pp2005 - eta - HT2"); \
	TDataProcessorPool *poolPPSimuNoweightEtaMB = new TDataProcessorPool("simulation_pp2005_eta_MB_noweight", "Simulation - pp2005 - eta - MB - no weight"); \
	TDataProcessorPool *poolPPSimuNoweightEtaHT1 = new TDataProcessorPool("simulation_pp2005_eta_HT1_noweight", "Simulation - pp2005 - eta - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimuNoweightEtaHT2 = new TDataProcessorPool("simulation_pp2005_eta_HT2_noweight", "Simulation - pp2005 - eta - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimuEtabgMB = new TDataProcessorPool("simulation_pp2005_etabg_MB", "Simulation - pp2005 - etabg - MB"); \
	TDataProcessorPool *poolPPSimuEtabgHT1 = new TDataProcessorPool("simulation_pp2005_etabg_HT1", "Simulation - pp2005 - etabg - HT1"); \
	TDataProcessorPool *poolPPSimuEtabgHT2 = new TDataProcessorPool("simulation_pp2005_etabg_HT2", "Simulation - pp2005 - etabg - HT2"); \
	TDataProcessorPool *poolPPSimuNoweightEtabgMB = new TDataProcessorPool("simulation_pp2005_etabg_MB_noweight", "Simulation - pp2005 - etabg - MB - no weight"); \
	TDataProcessorPool *poolPPSimuNoweightEtabgHT1 = new TDataProcessorPool("simulation_pp2005_etabg_HT1_noweight", "Simulation - pp2005 - etabg - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimuNoweightEtabgHT2 = new TDataProcessorPool("simulation_pp2005_etabg_HT2_noweight", "Simulation - pp2005 - etabg - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaMB = new TDataProcessorPool("simulation_pp2005_gamma_MB", "Simulation - pp2005 - 1 photon/event - MB"); \
	TDataProcessorPool *poolPPSimu1gammaHT1 = new TDataProcessorPool("simulation_pp2005_gamma_HT1", "Simulation - pp2005 - 1 photon/event - HT1"); \
	TDataProcessorPool *poolPPSimu1gammaHT2 = new TDataProcessorPool("simulation_pp2005_gamma_HT2", "Simulation - pp2005 - 1 photon/event - HT2"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightMB = new TDataProcessorPool("simulation_pp2005_gamma_MB_noweight", "Simulation - pp2005 - 1 photon/event - MB - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightHT1 = new TDataProcessorPool("simulation_pp2005_gamma_HT1_noweight", "Simulation - pp2005 - 1 photon/event - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightHT2 = new TDataProcessorPool("simulation_pp2005_gamma_HT2_noweight", "Simulation - pp2005 - 1 photon/event - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaSimMB = new TDataProcessorPool("simulation_pp2005_gamma_pi0_MB", "Simulation - pp2005 - 1 photon/event pi0 - MB"); \
	TDataProcessorPool *poolPPSimu1gammaSimHT1 = new TDataProcessorPool("simulation_pp2005_gamma_pi0_HT1", "Simulation - pp2005 - 1 photon/event pi0 - HT1"); \
	TDataProcessorPool *poolPPSimu1gammaSimHT2 = new TDataProcessorPool("simulation_pp2005_gamma_pi0_HT2", "Simulation - pp2005 - 1 photon/event pi0 - HT2"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightSimMB = new TDataProcessorPool("simulation_pp2005_gamma_pi0_MB_noweight", "Simulation - pp2005 - 1 photon/event pi0 - MB - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightSimHT1 = new TDataProcessorPool("simulation_pp2005_gamma_pi0_HT1_noweight", "Simulation - pp2005 - 1 photon/event pi0 - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightSimHT2 = new TDataProcessorPool("simulation_pp2005_gamma_pi0_HT2_noweight", "Simulation - pp2005 - 1 photon/event pi0 - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaSimEtaMB = new TDataProcessorPool("simulation_pp2005_gamma_eta_MB", "Simulation - pp2005 - 1 photon/event eta - MB"); \
	TDataProcessorPool *poolPPSimu1gammaSimEtaHT1 = new TDataProcessorPool("simulation_pp2005_gamma_eta_HT1", "Simulation - pp2005 - 1 photon/event eta - HT1"); \
	TDataProcessorPool *poolPPSimu1gammaSimEtaHT2 = new TDataProcessorPool("simulation_pp2005_gamma_eta_HT2", "Simulation - pp2005 - 1 photon/event eta - HT2"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightSimEtaMB = new TDataProcessorPool("simulation_pp2005_gamma_eta_MB_noweight", "Simulation - pp2005 - 1 photon/event eta - MB - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightSimEtaHT1 = new TDataProcessorPool("simulation_pp2005_gamma_eta_HT1_noweight", "Simulation - pp2005 - 1 photon/event eta - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimu1gammaNoweightSimEtaHT2 = new TDataProcessorPool("simulation_pp2005_gamma_eta_HT2_noweight", "Simulation - pp2005 - 1 photon/event eta - HT2 - no weight");

#define DATA_FILE_LOAD_PP(READ, SETTINGS, RESULTS) \
	{ \
	if (READ) { \
		calculatePrescaleTriggerBias(0, SETTINGS, RESULTS, 0); \
		poolPPMB->processFile(filenamePPDataMB_0); \
		calculatePrescaleTriggerBias(poolPPMB, SETTINGS, RESULTS, 1); \
		poolPPMB->processFile(filenamePPDataMB_1); \
		/*poolPPMB->processFile(filenamePPDataMB_2);*/ \
		calculatePrescaleTriggerBias(poolPPMB, SETTINGS, RESULTS, 2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPMB, dataFolder); \
		poolPPHT1->processFile(filenamePPDataHT1_0); \
		/*poolPPHT1->processFile(filenamePPDataHT1_1);*/ \
		calculatePrescaleTriggerBias(poolPPHT1, SETTINGS, RESULTS, 3); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPHT1, dataFolder); \
		poolPPHT2->processFile(filenamePPDataHT2_0); \
		/*poolPPHT2->processFile(filenamePPDataHT2_1);*/ \
		calculatePrescaleTriggerBias(poolPPHT2, SETTINGS, RESULTS, 4); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPHT2, dataFolder); \
	} \
	if (allSettings.usePPEff) { \
			if (allSettings.usePPEffWeight || SETTINGS.correctForPtShift || SETTINGS.correctForPtShiftEff) { \
				poolPPSimuMB->processFile(filenamePPSimuMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuMB, dataFolder); \
				poolPPSimuHT1->processFile(filenamePPSimuHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuHT1, dataFolder); \
				poolPPSimuHT2->processFile(filenamePPSimuHT2); \
   				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuHT2, dataFolder); \
			} \
			if (!allSettings.usePPEffWeight) { \
				poolPPSimuNoweightMB->processFile(filenamePPSimuNoweightMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightMB, dataFolder); \
				poolPPSimuNoweightHT1->processFile(filenamePPSimuNoweightHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightHT1, dataFolder); \
				poolPPSimuNoweightHT2->processFile(filenamePPSimuNoweightHT2); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightHT2, dataFolder); \
			} \
	} \
	if (allSettings.usePPEffEta) { \
			if (allSettings.usePPEffWeightEta || SETTINGS.correctForPtShift || SETTINGS.correctForPtShiftEff) { \
				poolPPSimuEtaMB->processFile(filenamePPSimuEtaMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuEtaMB, dataFolder); \
				poolPPSimuEtaHT1->processFile(filenamePPSimuEtaHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuEtaHT1, dataFolder); \
				poolPPSimuEtaHT2->processFile(filenamePPSimuEtaHT2); \
   				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuEtaHT2, dataFolder); \
				if (allSettings.usePPEffEtabg) { \
				    poolPPSimuEtabgMB->processFile(filenamePPSimuEtabgMB); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuEtabgMB, dataFolder); \
				    poolPPSimuEtabgHT1->processFile(filenamePPSimuEtabgHT1); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuEtabgHT1, dataFolder); \
				    poolPPSimuEtabgHT2->processFile(filenamePPSimuEtabgHT2); \
   				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuEtabgHT2, dataFolder); \
				} \
			} \
			if (!allSettings.usePPEffWeightEta) { \
				poolPPSimuNoweightEtaMB->processFile(filenamePPSimuNoweightEtaMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightEtaMB, dataFolder); \
				poolPPSimuNoweightEtaHT1->processFile(filenamePPSimuNoweightEtaHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightEtaHT1, dataFolder); \
				poolPPSimuNoweightEtaHT2->processFile(filenamePPSimuNoweightEtaHT2); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightEtaHT2, dataFolder); \
				if (allSettings.usePPEffEtabg) { \
				    poolPPSimuNoweightEtabgMB->processFile(filenamePPSimuNoweightEtabgMB); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightEtabgMB, dataFolder); \
				    poolPPSimuNoweightEtabgHT1->processFile(filenamePPSimuNoweightEtabgHT1); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightEtabgHT1, dataFolder); \
				    poolPPSimuNoweightEtabgHT2->processFile(filenamePPSimuNoweightEtabgHT2); \
    				    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuNoweightEtabgHT2, dataFolder); \
				} \
			} \
	} \
	if (allSettings.usePP1gamma) { \
		poolPPSimu1gammaMB->processFile(filenamePPSimu1gammaMB); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaMB, dataFolder); \
		poolPPSimu1gammaHT1->processFile(filenamePPSimu1gammaHT1); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaHT1, dataFolder); \
		poolPPSimu1gammaHT2->processFile(filenamePPSimu1gammaHT2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaHT2, dataFolder); \
		if (!allSettings.usePPEffWeight) { \
		    poolPPSimu1gammaNoweightMB->processFile(filenamePPSimu1gammaNoweightMB); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightMB, dataFolder); \
		    poolPPSimu1gammaNoweightHT1->processFile(filenamePPSimu1gammaNoweightHT1); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightHT1, dataFolder); \
		    poolPPSimu1gammaNoweightHT2->processFile(filenamePPSimu1gammaNoweightHT2); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightHT2, dataFolder); \
		} \
	} \
	if (allSettings.usePP1gammaSim) { \
		poolPPSimu1gammaSimMB->processFile(filenamePPSimu1gammaSimMB); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaSimMB, dataFolder); \
		poolPPSimu1gammaSimHT1->processFile(filenamePPSimu1gammaSimHT1); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaSimHT1, dataFolder); \
		poolPPSimu1gammaSimHT2->processFile(filenamePPSimu1gammaSimHT2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaSimHT2, dataFolder); \
		if (!allSettings.usePPEffWeight) { \
		    poolPPSimu1gammaNoweightSimMB->processFile(filenamePPSimu1gammaNoweightSimMB); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightSimMB, dataFolder); \
		    poolPPSimu1gammaNoweightSimHT1->processFile(filenamePPSimu1gammaNoweightSimHT1); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightSimHT1, dataFolder); \
		    poolPPSimu1gammaNoweightSimHT2->processFile(filenamePPSimu1gammaNoweightSimHT2); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightSimHT2, dataFolder); \
		} \
	} \
	if (allSettings.usePP1gammaSimEta) { \
		poolPPSimu1gammaSimEtaMB->processFile(filenamePPSimu1gammaSimEtaMB); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaSimEtaMB, dataFolder); \
		poolPPSimu1gammaSimEtaHT1->processFile(filenamePPSimu1gammaSimEtaHT1); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaSimEtaHT1, dataFolder); \
		poolPPSimu1gammaSimEtaHT2->processFile(filenamePPSimu1gammaSimEtaHT2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaSimEtaHT2, dataFolder); \
		if (!allSettings.usePPEffWeight) { \
		    poolPPSimu1gammaNoweightSimEtaMB->processFile(filenamePPSimu1gammaNoweightSimEtaMB); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightSimEtaMB, dataFolder); \
		    poolPPSimu1gammaNoweightSimEtaHT1->processFile(filenamePPSimu1gammaNoweightSimEtaHT1); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightSimEtaHT1, dataFolder); \
		    poolPPSimu1gammaNoweightSimEtaHT2->processFile(filenamePPSimu1gammaNoweightSimEtaHT2); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimu1gammaNoweightSimEtaHT2, dataFolder); \
		} \
	} \
	}

#define DATA_FILE_DELETE_PP \
	{ \
	removePool(poolPPMB, true, true);\
	removePool(poolPPHT1, true, true);\
	removePool(poolPPHT2, true, true);\
	removePool(poolPPSimuMB, true, true);\
	removePool(poolPPSimuHT1, true, true);\
	removePool(poolPPSimuHT2, true, true);\
	removePool(poolPPSimuNoweightMB, true, true);\
	removePool(poolPPSimuNoweightHT1, true, true);\
	removePool(poolPPSimuNoweightHT2, true, true);\
	removePool(poolPPSimuEtaMB, true, true);\
	removePool(poolPPSimuEtaHT1, true, true);\
	removePool(poolPPSimuEtaHT2, true, true);\
	removePool(poolPPSimuNoweightEtaMB, true, true);\
	removePool(poolPPSimuNoweightEtaHT1, true, true);\
	removePool(poolPPSimuNoweightEtaHT2, true, true);\
	removePool(poolPPSimuEtabgMB, true, true);\
	removePool(poolPPSimuEtabgHT1, true, true);\
	removePool(poolPPSimuEtabgHT2, true, true);\
	removePool(poolPPSimuNoweightEtabgMB, true, true);\
	removePool(poolPPSimuNoweightEtabgHT1, true, true);\
	removePool(poolPPSimuNoweightEtabgHT2, true, true);\
	removePool(poolPPSimu1gammaMB, true, true);\
	removePool(poolPPSimu1gammaHT1, true, true);\
	removePool(poolPPSimu1gammaHT2, true, true);\
	removePool(poolPPSimu1gammaNoweightMB, true, true);\
	removePool(poolPPSimu1gammaNoweightHT1, true, true);\
	removePool(poolPPSimu1gammaNoweightHT2, true, true);\
	removePool(poolPPSimu1gammaSimMB, true, true);\
	removePool(poolPPSimu1gammaSimHT1, true, true);\
	removePool(poolPPSimu1gammaSimHT2, true, true);\
	removePool(poolPPSimu1gammaNoweightSimMB, true, true);\
	removePool(poolPPSimu1gammaNoweightSimHT1, true, true);\
	removePool(poolPPSimu1gammaNoweightSimHT2, true, true);\
	removePool(poolPPSimu1gammaSimEtaMB, true, true);\
	removePool(poolPPSimu1gammaSimEtaHT1, true, true);\
	removePool(poolPPSimu1gammaSimEtaHT2, true, true);\
	removePool(poolPPSimu1gammaNoweightSimEtaMB, true, true);\
	removePool(poolPPSimu1gammaNoweightSimEtaHT1, true, true);\
	removePool(poolPPSimu1gammaNoweightSimEtaHT2, true, true); \
	}

	DATA_FILE_DAU(Nocentral, "nocentral");
	DATA_FILE_DAU(Allcentral, "0-100");
	DATA_FILE_DAU(Mostcentral, "0-20");
	DATA_FILE_DAU(Midcentral, "20-40");
	DATA_FILE_DAU(Mostperipheral, "40-100");

	TDataProcessorPool *poolEffMBNo = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? poolEmbedMBNocentral : poolSimuMBNocentral) : (allSettings.useDAuEmbedding ? poolEmbedNoweightMBNocentral : poolSimuNoweightMBNocentral);
	TDataProcessorPool *poolEffMBAll = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBAllcentral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBAllcentral : poolSimuMBNocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightMBAllcentral : poolEmbedNoweightMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightMBAllcentral : poolSimuNoweightMBNocentral));
	TDataProcessorPool *poolEffMBMostC = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBMostcentral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBMostcentral : poolSimuMBNocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightMBMostcentral : poolEmbedNoweightMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightMBMostcentral : poolSimuNoweightMBNocentral));
	TDataProcessorPool *poolEffMBMidC = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBMidcentral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBMidcentral : poolSimuMBNocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightMBMidcentral : poolEmbedNoweightMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightMBMidcentral : poolSimuNoweightMBNocentral));
	TDataProcessorPool *poolEffMBMostP = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBMostperipheral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBMostperipheral : poolSimuMBNocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightMBMostperipheral : poolEmbedNoweightMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightMBMostperipheral : poolSimuNoweightMBNocentral));

	TDataProcessorPool *poolEffHT1No = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? poolEmbedHT1Nocentral : poolSimuHT1Nocentral) : (allSettings.useDAuEmbedding ? poolEmbedNoweightHT1Nocentral : poolSimuNoweightHT1Nocentral);
	TDataProcessorPool *poolEffHT1All = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Allcentral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Allcentral : poolSimuHT1Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT1Allcentral : poolEmbedNoweightHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT1Allcentral : poolSimuNoweightHT1Nocentral));
	TDataProcessorPool *poolEffHT1MostC = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Mostcentral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Mostcentral : poolSimuHT1Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT1Mostcentral : poolEmbedNoweightHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT1Mostcentral : poolSimuNoweightHT1Nocentral));
	TDataProcessorPool *poolEffHT1MidC = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Midcentral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Midcentral : poolSimuHT1Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT1Midcentral : poolEmbedNoweightHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT1Midcentral : poolSimuNoweightHT1Nocentral));
	TDataProcessorPool *poolEffHT1MostP = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Mostperipheral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Mostperipheral : poolSimuHT1Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT1Mostperipheral : poolEmbedNoweightHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT1Mostperipheral : poolSimuNoweightHT1Nocentral));

	TDataProcessorPool *poolEffHT2No = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? poolEmbedHT2Nocentral : poolSimuHT2Nocentral) : (allSettings.useDAuEmbedding ? poolEmbedNoweightHT2Nocentral : poolSimuNoweightHT2Nocentral);
	TDataProcessorPool *poolEffHT2All = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Allcentral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Allcentral : poolSimuHT2Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT2Allcentral : poolEmbedNoweightHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT2Allcentral : poolSimuNoweightHT2Nocentral));
	TDataProcessorPool *poolEffHT2MostC = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Mostcentral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Mostcentral : poolSimuHT2Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT2Mostcentral : poolEmbedNoweightHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT2Mostcentral : poolSimuNoweightHT2Nocentral));
	TDataProcessorPool *poolEffHT2MidC = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Midcentral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Midcentral : poolSimuHT2Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT2Midcentral : poolEmbedNoweightHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT2Midcentral : poolSimuNoweightHT2Nocentral));
	TDataProcessorPool *poolEffHT2MostP = allSettings.useDAuEffWeight ? (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Mostperipheral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Mostperipheral : poolSimuHT2Nocentral)) : (allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedNoweightHT2Mostperipheral : poolEmbedNoweightHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuNoweightHT2Mostperipheral : poolSimuNoweightHT2Nocentral));

	TDataProcessorPool *poolEffEtaMBNo = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? poolEmbedEtaMBNocentral : poolSimuEtaMBNocentral) : (allSettings.useDAuEmbeddingEta ? poolEmbedNoweightEtaMBNocentral : poolSimuNoweightEtaMBNocentral);
	TDataProcessorPool *poolEffEtaMBAll = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaMBAllcentral : poolEmbedEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaMBAllcentral : poolSimuEtaMBNocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaMBAllcentral : poolEmbedNoweightEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaMBAllcentral : poolSimuNoweightEtaMBNocentral));
	TDataProcessorPool *poolEffEtaMBMostC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaMBMostcentral : poolEmbedEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaMBMostcentral : poolSimuEtaMBNocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaMBMostcentral : poolEmbedNoweightEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaMBMostcentral : poolSimuNoweightEtaMBNocentral));
	TDataProcessorPool *poolEffEtaMBMidC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaMBMidcentral : poolEmbedEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaMBMidcentral : poolSimuEtaMBNocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaMBMidcentral : poolEmbedNoweightEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaMBMidcentral : poolSimuNoweightEtaMBNocentral));
	TDataProcessorPool *poolEffEtaMBMostP = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaMBMostperipheral : poolEmbedEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaMBMostperipheral : poolSimuEtaMBNocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaMBMostperipheral : poolEmbedNoweightEtaMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaMBMostperipheral : poolSimuNoweightEtaMBNocentral));

	TDataProcessorPool *poolEffEtaHT1No = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? poolEmbedEtaHT1Nocentral : poolSimuEtaHT1Nocentral) : (allSettings.useDAuEmbeddingEta ? poolEmbedNoweightEtaHT1Nocentral : poolSimuNoweightEtaHT1Nocentral);
	TDataProcessorPool *poolEffEtaHT1All = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT1Allcentral : poolEmbedEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT1Allcentral : poolSimuEtaHT1Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT1Allcentral : poolEmbedNoweightEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT1Allcentral : poolSimuNoweightEtaHT1Nocentral));
	TDataProcessorPool *poolEffEtaHT1MostC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT1Mostcentral : poolEmbedEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT1Mostcentral : poolSimuEtaHT1Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT1Mostcentral : poolEmbedNoweightEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT1Mostcentral : poolSimuNoweightEtaHT1Nocentral));
	TDataProcessorPool *poolEffEtaHT1MidC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT1Midcentral : poolEmbedEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT1Midcentral : poolSimuEtaHT1Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT1Midcentral : poolEmbedNoweightEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT1Midcentral : poolSimuNoweightEtaHT1Nocentral));
	TDataProcessorPool *poolEffEtaHT1MostP = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT1Mostperipheral : poolEmbedEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT1Mostperipheral : poolSimuEtaHT1Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT1Mostperipheral : poolEmbedNoweightEtaHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT1Mostperipheral : poolSimuNoweightEtaHT1Nocentral));

	TDataProcessorPool *poolEffEtaHT2No = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? poolEmbedEtaHT2Nocentral : poolSimuEtaHT2Nocentral) : (allSettings.useDAuEmbeddingEta ? poolEmbedNoweightEtaHT2Nocentral : poolSimuNoweightEtaHT2Nocentral);
	TDataProcessorPool *poolEffEtaHT2All = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT2Allcentral : poolEmbedEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT2Allcentral : poolSimuEtaHT2Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT2Allcentral : poolEmbedNoweightEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT2Allcentral : poolSimuNoweightEtaHT2Nocentral));
	TDataProcessorPool *poolEffEtaHT2MostC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT2Mostcentral : poolEmbedEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT2Mostcentral : poolSimuEtaHT2Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT2Mostcentral : poolEmbedNoweightEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT2Mostcentral : poolSimuNoweightEtaHT2Nocentral));
	TDataProcessorPool *poolEffEtaHT2MidC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT2Midcentral : poolEmbedEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT2Midcentral : poolSimuEtaHT2Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT2Midcentral : poolEmbedNoweightEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT2Midcentral : poolSimuNoweightEtaHT2Nocentral));
	TDataProcessorPool *poolEffEtaHT2MostP = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedEtaHT2Mostperipheral : poolEmbedEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuEtaHT2Mostperipheral : poolSimuEtaHT2Nocentral)) : (allSettings.useDAuEmbeddingEta ? (allSettings.useDAuEffCentralitiesEta ? poolEmbedNoweightEtaHT2Mostperipheral : poolEmbedNoweightEtaHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtaHT2Mostperipheral : poolSimuNoweightEtaHT2Nocentral));

	TDataProcessorPool *poolEffEtabgMBNo = allSettings.useDAuEffWeightEta ? poolSimuEtabgMBNocentral : poolSimuNoweightEtabgMBNocentral;
	TDataProcessorPool *poolEffEtabgMBAll = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgMBAllcentral : poolSimuEtabgMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgMBAllcentral : poolSimuNoweightEtabgMBNocentral);
	TDataProcessorPool *poolEffEtabgMBMostC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgMBMostcentral : poolSimuEtabgMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgMBMostcentral : poolSimuNoweightEtabgMBNocentral);
	TDataProcessorPool *poolEffEtabgMBMidC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgMBMidcentral : poolSimuEtabgMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgMBMidcentral : poolSimuNoweightEtabgMBNocentral);
	TDataProcessorPool *poolEffEtabgMBMostP = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgMBMostperipheral : poolSimuEtabgMBNocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgMBMostperipheral : poolSimuNoweightEtabgMBNocentral);

	TDataProcessorPool *poolEffEtabgHT1No = allSettings.useDAuEffWeightEta ? poolSimuEtabgHT1Nocentral : poolSimuNoweightEtabgHT1Nocentral;
	TDataProcessorPool *poolEffEtabgHT1All = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT1Allcentral : poolSimuEtabgHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT1Allcentral : poolSimuNoweightEtabgHT1Nocentral);
	TDataProcessorPool *poolEffEtabgHT1MostC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT1Mostcentral : poolSimuEtabgHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT1Mostcentral : poolSimuNoweightEtabgHT1Nocentral);
	TDataProcessorPool *poolEffEtabgHT1MidC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT1Midcentral : poolSimuEtabgHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT1Midcentral : poolSimuNoweightEtabgHT1Nocentral);
	TDataProcessorPool *poolEffEtabgHT1MostP = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT1Mostperipheral : poolSimuEtabgHT1Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT1Mostperipheral : poolSimuNoweightEtabgHT1Nocentral);

	TDataProcessorPool *poolEffEtabgHT2No = allSettings.useDAuEffWeightEta ? poolSimuEtabgHT2Nocentral : poolSimuNoweightEtabgHT2Nocentral;
	TDataProcessorPool *poolEffEtabgHT2All = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT2Allcentral : poolSimuEtabgHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT2Allcentral : poolSimuNoweightEtabgHT2Nocentral);
	TDataProcessorPool *poolEffEtabgHT2MostC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT2Mostcentral : poolSimuEtabgHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT2Mostcentral : poolSimuNoweightEtabgHT2Nocentral);
	TDataProcessorPool *poolEffEtabgHT2MidC = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT2Midcentral : poolSimuEtabgHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT2Midcentral : poolSimuNoweightEtabgHT2Nocentral);
	TDataProcessorPool *poolEffEtabgHT2MostP = allSettings.useDAuEffWeightEta ? (allSettings.useDAuEffCentralitiesEta ? poolSimuEtabgHT2Mostperipheral : poolSimuEtabgHT2Nocentral) : (allSettings.useDAuEffCentralitiesEta ? poolSimuNoweightEtabgHT2Mostperipheral : poolSimuNoweightEtabgHT2Nocentral);

	TDataProcessorPool *pool1gammaSimMBNo = allSettings.useDAuEffWeight ? (poolSimu1gammaSimMBNocentral) : (poolSimu1gammaNoweightSimMBNocentral);
	TDataProcessorPool *pool1gammaSimMBAll = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimMBAllcentral : poolSimu1gammaSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimMBAllcentral : poolSimu1gammaNoweightSimMBNocentral));
	TDataProcessorPool *pool1gammaSimMBMostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimMBMostcentral : poolSimu1gammaSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimMBMostcentral : poolSimu1gammaNoweightSimMBNocentral));
	TDataProcessorPool *pool1gammaSimMBMidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimMBMidcentral : poolSimu1gammaSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimMBMidcentral : poolSimu1gammaNoweightSimMBNocentral));
	TDataProcessorPool *pool1gammaSimMBMostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimMBMostperipheral : poolSimu1gammaSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimMBMostperipheral : poolSimu1gammaNoweightSimMBNocentral));

	TDataProcessorPool *pool1gammaSimHT1No = allSettings.useDAuEffWeight ? (poolSimu1gammaSimHT1Nocentral) : (poolSimu1gammaNoweightSimHT1Nocentral);
	TDataProcessorPool *pool1gammaSimHT1All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT1Allcentral : poolSimu1gammaSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT1Allcentral : poolSimu1gammaNoweightSimHT1Nocentral));
	TDataProcessorPool *pool1gammaSimHT1MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT1Mostcentral : poolSimu1gammaSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT1Mostcentral : poolSimu1gammaNoweightSimHT1Nocentral));
	TDataProcessorPool *pool1gammaSimHT1MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT1Midcentral : poolSimu1gammaSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT1Midcentral : poolSimu1gammaNoweightSimHT1Nocentral));
	TDataProcessorPool *pool1gammaSimHT1MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT1Mostperipheral : poolSimu1gammaSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT1Mostperipheral : poolSimu1gammaNoweightSimHT1Nocentral));

	TDataProcessorPool *pool1gammaSimHT2No = allSettings.useDAuEffWeight ? (poolSimu1gammaSimHT2Nocentral) : (poolSimu1gammaNoweightSimHT2Nocentral);
	TDataProcessorPool *pool1gammaSimHT2All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT2Allcentral : poolSimu1gammaSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT2Allcentral : poolSimu1gammaNoweightSimHT2Nocentral));
	TDataProcessorPool *pool1gammaSimHT2MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT2Mostcentral : poolSimu1gammaSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT2Mostcentral : poolSimu1gammaNoweightSimHT2Nocentral));
	TDataProcessorPool *pool1gammaSimHT2MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT2Midcentral : poolSimu1gammaSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT2Midcentral : poolSimu1gammaNoweightSimHT2Nocentral));
	TDataProcessorPool *pool1gammaSimHT2MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimHT2Mostperipheral : poolSimu1gammaSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimHT2Mostperipheral : poolSimu1gammaNoweightSimHT2Nocentral));

	TDataProcessorPool *pool1gammaSimEtaMBNo = allSettings.useDAuEffWeight ? (poolSimu1gammaSimEtaMBNocentral) : (poolSimu1gammaNoweightSimEtaMBNocentral);
	TDataProcessorPool *pool1gammaSimEtaMBAll = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaMBAllcentral : poolSimu1gammaSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaMBAllcentral : poolSimu1gammaNoweightSimEtaMBNocentral));
	TDataProcessorPool *pool1gammaSimEtaMBMostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaMBMostcentral : poolSimu1gammaSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaMBMostcentral : poolSimu1gammaNoweightSimEtaMBNocentral));
	TDataProcessorPool *pool1gammaSimEtaMBMidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaMBMidcentral : poolSimu1gammaSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaMBMidcentral : poolSimu1gammaNoweightSimEtaMBNocentral));
	TDataProcessorPool *pool1gammaSimEtaMBMostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaMBMostperipheral : poolSimu1gammaSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaMBMostperipheral : poolSimu1gammaNoweightSimEtaMBNocentral));

	TDataProcessorPool *pool1gammaSimEtaHT1No = allSettings.useDAuEffWeight ? (poolSimu1gammaSimEtaHT1Nocentral) : (poolSimu1gammaNoweightSimEtaHT1Nocentral);
	TDataProcessorPool *pool1gammaSimEtaHT1All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT1Allcentral : poolSimu1gammaSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT1Allcentral : poolSimu1gammaNoweightSimEtaHT1Nocentral));
	TDataProcessorPool *pool1gammaSimEtaHT1MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT1Mostcentral : poolSimu1gammaSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT1Mostcentral : poolSimu1gammaNoweightSimEtaHT1Nocentral));
	TDataProcessorPool *pool1gammaSimEtaHT1MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT1Midcentral : poolSimu1gammaSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT1Midcentral : poolSimu1gammaNoweightSimEtaHT1Nocentral));
	TDataProcessorPool *pool1gammaSimEtaHT1MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT1Mostperipheral : poolSimu1gammaSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT1Mostperipheral : poolSimu1gammaNoweightSimEtaHT1Nocentral));

	TDataProcessorPool *pool1gammaSimEtaHT2No = allSettings.useDAuEffWeight ? (poolSimu1gammaSimEtaHT2Nocentral) : (poolSimu1gammaNoweightSimEtaHT2Nocentral);
	TDataProcessorPool *pool1gammaSimEtaHT2All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT2Allcentral : poolSimu1gammaSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT2Allcentral : poolSimu1gammaNoweightSimEtaHT2Nocentral));
	TDataProcessorPool *pool1gammaSimEtaHT2MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT2Mostcentral : poolSimu1gammaSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT2Mostcentral : poolSimu1gammaNoweightSimEtaHT2Nocentral));
	TDataProcessorPool *pool1gammaSimEtaHT2MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT2Midcentral : poolSimu1gammaSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT2Midcentral : poolSimu1gammaNoweightSimEtaHT2Nocentral));
	TDataProcessorPool *pool1gammaSimEtaHT2MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaSimEtaHT2Mostperipheral : poolSimu1gammaSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimu1gammaNoweightSimEtaHT2Mostperipheral : poolSimu1gammaNoweightSimEtaHT2Nocentral));

	TDataProcessorPool *pool1gammaMBNo = poolSimu1gammaMBNocentral;
	TDataProcessorPool *pool1gammaMBAll = allSettings.useDAu1gammaCentralities ? poolSimu1gammaMBAllcentral : poolSimu1gammaMBNocentral;
	TDataProcessorPool *pool1gammaMBMostC = allSettings.useDAu1gammaCentralities ? poolSimu1gammaMBMostcentral : poolSimu1gammaMBNocentral;
	TDataProcessorPool *pool1gammaMBMidC = allSettings.useDAu1gammaCentralities ? poolSimu1gammaMBMidcentral : poolSimu1gammaMBNocentral;
	TDataProcessorPool *pool1gammaMBMostP = allSettings.useDAu1gammaCentralities ? poolSimu1gammaMBMostperipheral : poolSimu1gammaMBNocentral;

	TDataProcessorPool *pool1gammaHT1No = poolSimu1gammaHT1Nocentral;
	TDataProcessorPool *pool1gammaHT1All = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT1Allcentral : poolSimu1gammaHT1Nocentral;
	TDataProcessorPool *pool1gammaHT1MostC = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT1Mostcentral : poolSimu1gammaHT1Nocentral;
	TDataProcessorPool *pool1gammaHT1MidC = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT1Midcentral : poolSimu1gammaHT1Nocentral;
	TDataProcessorPool *pool1gammaHT1MostP = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT1Mostperipheral : poolSimu1gammaHT1Nocentral;

	TDataProcessorPool *pool1gammaHT2No = poolSimu1gammaHT2Nocentral;
	TDataProcessorPool *pool1gammaHT2All = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT2Allcentral : poolSimu1gammaHT2Nocentral;
	TDataProcessorPool *pool1gammaHT2MostC = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT2Mostcentral : poolSimu1gammaHT2Nocentral;
	TDataProcessorPool *pool1gammaHT2MidC = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT2Midcentral : poolSimu1gammaHT2Nocentral;
	TDataProcessorPool *pool1gammaHT2MostP = allSettings.useDAu1gammaCentralities ? poolSimu1gammaHT2Mostperipheral : poolSimu1gammaHT2Nocentral;

	TDataProcessorPool *poolPtshiftMBNo = allSettings.useDAuEmbedding ? poolEmbedMBNocentral : poolSimuMBNocentral;
	TDataProcessorPool *poolPtshiftMBAll = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBAllcentral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBAllcentral : poolSimuMBNocentral);
	TDataProcessorPool *poolPtshiftMBMostC = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBMostcentral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBMostcentral : poolSimuMBNocentral);
	TDataProcessorPool *poolPtshiftMBMidC = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBMidcentral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBMidcentral : poolSimuMBNocentral);
	TDataProcessorPool *poolPtshiftMBMostP = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedMBMostperipheral : poolEmbedMBNocentral) : (allSettings.useDAuEffCentralities ? poolSimuMBMostperipheral : poolSimuMBNocentral);

	TDataProcessorPool *poolPtshiftHT1No = allSettings.useDAuEmbedding ? poolEmbedHT1Nocentral : poolSimuHT1Nocentral;
	TDataProcessorPool *poolPtshiftHT1All = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Allcentral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Allcentral : poolSimuHT1Nocentral);
	TDataProcessorPool *poolPtshiftHT1MostC = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Mostcentral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Mostcentral : poolSimuHT1Nocentral);
	TDataProcessorPool *poolPtshiftHT1MidC = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Midcentral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Midcentral : poolSimuHT1Nocentral);
	TDataProcessorPool *poolPtshiftHT1MostP = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT1Mostperipheral : poolEmbedHT1Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT1Mostperipheral : poolSimuHT1Nocentral);

	TDataProcessorPool *poolPtshiftHT2No = allSettings.useDAuEmbedding ? poolEmbedHT2Nocentral : poolSimuHT2Nocentral;
	TDataProcessorPool *poolPtshiftHT2All = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Allcentral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Allcentral : poolSimuHT2Nocentral);
	TDataProcessorPool *poolPtshiftHT2MostC = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Mostcentral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Mostcentral : poolSimuHT2Nocentral);
	TDataProcessorPool *poolPtshiftHT2MidC = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Midcentral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Midcentral : poolSimuHT2Nocentral);
	TDataProcessorPool *poolPtshiftHT2MostP = allSettings.useDAuEmbedding ? (allSettings.useDAuEffCentralities ? poolEmbedHT2Mostperipheral : poolEmbedHT2Nocentral) : (allSettings.useDAuEffCentralities ? poolSimuHT2Mostperipheral : poolSimuHT2Nocentral);
        
	DATA_FILE_PP;

	TDataProcessorPool *poolPPEffMB = allSettings.usePPEffWeight ? poolPPSimuMB : poolPPSimuNoweightMB;
	TDataProcessorPool *poolPPEffHT1 = allSettings.usePPEffWeight ? poolPPSimuHT1 : poolPPSimuNoweightHT1;
	TDataProcessorPool *poolPPEffHT2 = allSettings.usePPEffWeight ? poolPPSimuHT2 : poolPPSimuNoweightHT2;
	TDataProcessorPool *poolPPEffEtaMB = allSettings.usePPEffWeightEta ? poolPPSimuEtaMB : poolPPSimuNoweightEtaMB;
	TDataProcessorPool *poolPPEffEtaHT1 = allSettings.usePPEffWeightEta ? poolPPSimuEtaHT1 : poolPPSimuNoweightEtaHT1;
	TDataProcessorPool *poolPPEffEtaHT2 = allSettings.usePPEffWeightEta ? poolPPSimuEtaHT2 : poolPPSimuNoweightEtaHT2;
	TDataProcessorPool *poolPPEffEtabgMB = allSettings.usePPEffWeightEta ? poolPPSimuEtabgMB : poolPPSimuNoweightEtabgMB;
	TDataProcessorPool *poolPPEffEtabgHT1 = allSettings.usePPEffWeightEta ? poolPPSimuEtabgHT1 : poolPPSimuNoweightEtabgHT1;
	TDataProcessorPool *poolPPEffEtabgHT2 = allSettings.usePPEffWeightEta ? poolPPSimuEtabgHT2 : poolPPSimuNoweightEtabgHT2;
	TDataProcessorPool *poolPP1gammaMB = poolPPSimu1gammaMB;
	TDataProcessorPool *poolPP1gammaHT1 = poolPPSimu1gammaHT1;
	TDataProcessorPool *poolPP1gammaHT2 = poolPPSimu1gammaHT2;
	TDataProcessorPool *poolPP1gammaSimMB = allSettings.usePPEffWeight ? poolPPSimu1gammaSimMB : poolPPSimu1gammaNoweightSimMB;
	TDataProcessorPool *poolPP1gammaSimHT1 = allSettings.usePPEffWeight ? poolPPSimu1gammaSimHT1 : poolPPSimu1gammaNoweightSimHT1;
	TDataProcessorPool *poolPP1gammaSimHT2 = allSettings.usePPEffWeight ? poolPPSimu1gammaSimHT2 : poolPPSimu1gammaNoweightSimHT2;
	TDataProcessorPool *poolPP1gammaSimEtaMB = allSettings.usePPEffWeight ? poolPPSimu1gammaSimEtaMB : poolPPSimu1gammaNoweightSimEtaMB;
	TDataProcessorPool *poolPP1gammaSimEtaHT1 = allSettings.usePPEffWeight ? poolPPSimu1gammaSimEtaHT1 : poolPPSimu1gammaNoweightSimEtaHT1;
	TDataProcessorPool *poolPP1gammaSimEtaHT2 = allSettings.usePPEffWeight ? poolPPSimu1gammaSimEtaHT2 : poolPPSimu1gammaNoweightSimEtaHT2;
	TDataProcessorPool *poolPPPtshiftMB = poolPPSimuMB;
	TDataProcessorPool *poolPPPtshiftHT1 = poolPPSimuHT1;
	TDataProcessorPool *poolPPPtshiftHT2 = poolPPSimuHT2;

	if (allSettings.showSpectrumDAu || allSettings.showSpectrumDAuEta || allSettings.showRDA || allSettings.showRDAEta) {
	    DATA_FILE_LOAD_DAU(Nocentral, "nocentral", true, allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
	    showResults(poolMBNocentral, poolHT1Nocentral, poolHT2Nocentral
		, poolEffMBNo, poolEffHT1No, poolEffHT2No
		, poolEffEtaMBNo, poolEffEtaHT1No, poolEffEtaHT2No
		, poolEffEtabgMBNo, poolEffEtabgHT1No, poolEffEtabgHT2No
		, pool1gammaMBNo, pool1gammaHT1No, pool1gammaHT2No
		, pool1gammaSimMBNo, pool1gammaSimHT1No, pool1gammaSimHT2No
		, pool1gammaSimEtaMBNo, pool1gammaSimEtaHT1No, pool1gammaSimEtaHT2No
		, poolPtshiftMBNo, poolPtshiftHT1No, poolPtshiftHT2No
		, allSettings.settingsDAuNoCentral, !allSettings.browseAnalysis, allResults.resultsDAuNoCentral
	    );
	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU(Nocentral);
	}

	if (allSettings.showRcp || allSettings.showRcpEta) {
	    if (allSettings.showAllCentral) {    
		DATA_FILE_LOAD_DAU(Allcentral, "0-100", true, allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		showResults(poolMBAllcentral, poolHT1Allcentral, poolHT2Allcentral
    		    , poolEffMBAll, poolEffHT1All, poolEffHT2All
		    , poolEffEtaMBAll, poolEffEtaHT1All, poolEffEtaHT2All
		    , poolEffEtabgMBAll, poolEffEtabgHT1All, poolEffEtabgHT2All
		    , pool1gammaMBAll, pool1gammaHT1All, pool1gammaHT2All
		    , pool1gammaSimMBAll, pool1gammaSimHT1All, pool1gammaSimHT2All
		    , pool1gammaSimEtaMBAll, pool1gammaSimEtaHT1All, pool1gammaSimEtaHT2All
		    , poolPtshiftMBAll, poolPtshiftHT1All, poolPtshiftHT2All
		    , allSettings.settingsDAuAllCentral, !allSettings.browseAnalysis, allResults.resultsDAuAllCentral
		);
		if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU(Allcentral);
	    }

    	    DATA_FILE_LOAD_DAU(Mostcentral, "0-20", true, allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
	    showResults(poolMBMostcentral, poolHT1Mostcentral, poolHT2Mostcentral
		, poolEffMBMostC, poolEffHT1MostC, poolEffHT2MostC
		, poolEffEtaMBMostC, poolEffEtaHT1MostC, poolEffEtaHT2MostC
		, poolEffEtabgMBMostC, poolEffEtabgHT1MostC, poolEffEtabgHT2MostC
		, pool1gammaMBMostC, pool1gammaHT1MostC, pool1gammaHT2MostC
		, pool1gammaSimMBMostC, pool1gammaSimHT1MostC, pool1gammaSimHT2MostC
		, pool1gammaSimEtaMBMostC, pool1gammaSimEtaHT1MostC, pool1gammaSimEtaHT2MostC
		, poolPtshiftMBMostC, poolPtshiftHT1MostC, poolPtshiftHT2MostC
		, allSettings.settingsDAuMostCentral, !allSettings.browseAnalysis, allResults.resultsDAuMostCentral
	    );
    	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU(Mostcentral);

	    if (allSettings.showMidCentral) {
    		DATA_FILE_LOAD_DAU(Midcentral, "20-40", true, allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
		showResults(poolMBMidcentral, poolHT1Midcentral, poolHT2Midcentral
		    , poolEffMBMidC, poolEffHT1MidC, poolEffHT2MidC
		    , poolEffEtaMBMidC, poolEffEtaHT1MidC, poolEffEtaHT2MidC
		    , poolEffEtabgMBMidC, poolEffEtabgHT1MidC, poolEffEtabgHT2MidC
		    , pool1gammaMBMidC, pool1gammaHT1MidC, pool1gammaHT2MidC
		    , pool1gammaSimMBMidC, pool1gammaSimHT1MidC, pool1gammaSimHT2MidC
		    , pool1gammaSimEtaMBMidC, pool1gammaSimEtaHT1MidC, pool1gammaSimEtaHT2MidC
		    , poolPtshiftMBMidC, poolPtshiftHT1MidC, poolPtshiftHT2MidC
		    , allSettings.settingsDAuMidCentral, !allSettings.browseAnalysis, allResults.resultsDAuMidCentral
		);
        	if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU(Midcentral);
	    }

	    DATA_FILE_LOAD_DAU(Mostperipheral, "40-100", true, allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
	    showResults(poolMBMostperipheral, poolHT1Mostperipheral, poolHT2Mostperipheral
		, poolEffMBMostP, poolEffHT1MostP, poolEffHT2MostP
		, poolEffEtaMBMostP, poolEffEtaHT1MostP, poolEffEtaHT2MostP
		, poolEffEtabgMBMostP, poolEffEtabgHT1MostP, poolEffEtabgHT2MostP
		, pool1gammaMBMostP, pool1gammaHT1MostP, pool1gammaHT2MostP
		, pool1gammaSimMBMostP, pool1gammaSimHT1MostP, pool1gammaSimHT2MostP
		, pool1gammaSimEtaMBMostP, pool1gammaSimEtaHT1MostP, pool1gammaSimEtaHT2MostP
		, poolPtshiftMBMostP, poolPtshiftHT1MostP, poolPtshiftHT2MostP
		, allSettings.settingsDAuMostPeripheral, !allSettings.browseAnalysis, allResults.resultsDAuMostPeripheral
	    );
    	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU(Mostperipheral);

	    showResultsDAuCentrality(allSettings.dAuCentralityName, allSettings.dAuCentralityTitle
                , allSettings.showRcp, allSettings.showRcpEta, allSettings.dAuCentralityShow, allSettings.dAuCentralityPrint, allSettings.showDAuCentralityTriggersSeparately
                , allSettings.settingsDAuNoCentral
                , allSettings.settingsDAuAllCentral
                , allSettings.settingsDAuMostCentral
                , allSettings.settingsDAuMidCentral
                , allSettings.settingsDAuMostPeripheral
                , allResults.resultsDAuNoCentral
                , allResults.resultsDAuMostCentral
                , allResults.resultsDAuMidCentral
                , allResults.resultsDAuMostPeripheral
                , allResults.dAuRcpMB
                , allResults.dAuRcpHT1
                , allResults.dAuRcpHT2
                , allResults.dAuRcp
                , allResults.dAuRcpMBEta
                , allResults.dAuRcpHT1Eta
                , allResults.dAuRcpHT2Eta
                , allResults.dAuRcpEta
	    );
	}

	if (allSettings.showSpectrumPP || allSettings.showSpectrumPPEta || allSettings.showRDA || allSettings.showRDAEta) {
	    DATA_FILE_LOAD_PP(true, allSettings.settingsPP, allResults.resultsPP);
	    showResults(poolPPMB, poolPPHT1, poolPPHT2
		, poolPPEffMB, poolPPEffHT1, poolPPEffHT2
		, poolPPEffEtaMB, poolPPEffEtaHT1, poolPPEffEtaHT2
		, poolPPEffEtabgMB, poolPPEffEtabgHT1, poolPPEffEtabgHT2
		, poolPP1gammaMB, poolPP1gammaHT1, poolPP1gammaHT2
		, poolPP1gammaSimMB, poolPP1gammaSimHT1, poolPP1gammaSimHT2
		, poolPP1gammaSimEtaMB, poolPP1gammaSimEtaHT1, poolPP1gammaSimEtaHT2
		, poolPPPtshiftMB, poolPPPtshiftHT1, poolPPPtshiftHT2
		, allSettings.settingsPP, !allSettings.browseAnalysis, allResults.resultsPP
	    );
	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_PP;
	}

	if (allSettings.showRDA || allSettings.showRDAEta) {
	    showResultsRDA(allSettings.RDAName, allSettings.RDATitle
	        , allSettings.showRDA, allSettings.showRDAEta, allSettings.RDAShow, allSettings.showRDATriggersSeparately
	        , allSettings.settingsDAuNoCentral
	        , allSettings.settingsPP
	        , allResults.resultsDAuNoCentral
	        , allResults.resultsPP
	        , allResults.RDAMB
	        , allResults.RDAHT1
	        , allResults.RDAHT2
	        , allResults.RDA
	        , allResults.RDAMBEta
	        , allResults.RDAHT1Eta
	        , allResults.RDAHT2Eta
	        , allResults.RDAEta
	    );
	}

    	if (dataFolder && allSettings.browseResults) memberInspector.Inspect(allSettings, dataFolder, "allSettings", "Analysis settings");
    	if (dataFolder && allSettings.browseResults) memberInspector.Inspect(allResults, dataFolder, "allResults", "Analysis results");
	if (dataFolder && (allSettings.browseAnalysis || allSettings.browseResults)) {
    	    TBrowser *browser = new TBrowser("dataBrowser", dataFolder, "Processed data browser");
    	    if (browser) {
		TCanvas *browserCanvas = new TCanvas("Canvas", "Processed data plot");
		if (browserCanvas) browserCanvas->cd();
	    }
	}

#undef DATA_FILE_DAU
#undef DATA_FILE_PP
#undef DATA_FILE_LOAD_DAU
#undef DATA_FILE_LOAD_PP
#undef DATA_FILE_DELETE_DAU
#undef DATA_FILE_DELETE_PP

	TDatime stopPreTime;
	cout << "Pre-finished: " << stopPreTime.AsSQLString() << endl;
	gSystem->ProcessEvents();
	TDatime stopTime;
	cout << "Finished: " << stopTime.AsSQLString() << endl;
	timer.Stop();
	timer.Print();
	cout << "================== FINISHED ==================" << endl;
}
