#include <TROOT.h>
#include <TSystem.h>
#include <TNamed.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TFrame.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TFolder.h>
#include <TBrowser.h>

#include <StEmcPool/StPi0Common/StPi0CommonUtil.h>

#include <StEmcPool/StPi0Analysis/TWeightCalculator.h>
#include <StEmcPool/StPi0Analysis/TDataProcessor.h>
#include <StEmcPool/StPi0Analysis/TDataProcessorPool.h>

#include "StPi0ResultsUtil.h"
#include "TDataProcessorMemberInspector.h"
#include "StPi0ResultsFinal.h"

ClassImp(TAllSettings);
ClassImp(TAllResults);

//#pragma GCC diagnostic ignored "-Wunused-variable"

//----------------------------------------------------------------------
void show_analysis_final(const Char_t *DATA_DIR) {

	cout << "Reading configuration..." << endl;

	getNLOpQCD(0, 0);
	getNLOpQCDPP(0, 0);

	Bool_t isPythia = false;

	Bool_t addDirName = false;

	Bool_t extendMCpT = true;

	Float_t vertEffDAuMinBias = 0.93;//0.837;
	Float_t vertEffDAuMostPeripheralMinBias = 0.88;//0.837;
	Float_t vertEffPPMinBias = 1.0;//0.65; // 0.85
	TFitSettings fitAllLeft = TFitSettings(0);
	TFitSettings fitAllRight = TFitSettings(0.9);
	TFitSettings fitPi0PeakLeft = TFitSettings(0.075 + 0.015, 0.025/15.0);
	TFitSettings fitPi0PeakRight = TFitSettings(0.25 - 0.080, 0.05/15.0);
	TFitSettings fitEtaPeakLeft = TFitSettings(0.35, 0.05/15.0);
	TFitSettings fitEtaPeakRight = TFitSettings(0.75, 0.0);
	TFitSettings fitLeft = fitPi0PeakLeft;
	TFitSettings fitRight = fitPi0PeakRight;
	TFitSettings fitLeftEta = fitEtaPeakLeft;
	TFitSettings fitRightEta = fitEtaPeakRight;

const Float_t pp_MB = 0.0;
const Float_t pp_MB_HT1 = 4.0;
const Float_t pp_HT1_HT2 = 7.0;
const Float_t pp_HT2 = 17.0;

const Float_t dAu_MB = 0.0;
const Float_t dAu_MB_HT1 = 4.0;
const Float_t dAu_HT1_HT2 = 7.0;
const Float_t dAu_HT2 = 17.0;

const Float_t pp_etatopi0_MB = 0.0;
const Float_t pp_etatopi0_MB_HT1 = 4.0;
const Float_t pp_etatopi0_HT1_HT2 = 7.0;
const Float_t pp_etatopi0_HT2 = 14.0;

const Float_t dAu_etatopi0_MB = 1.5;
const Float_t dAu_etatopi0_MB_HT1 = 4.0;
const Float_t dAu_etatopi0_HT1_HT2 = 7.0;
const Float_t dAu_etatopi0_HT2 = 14.0;

const Float_t RdA_MB = 0.0;
const Float_t RdA_MB_HT1 = 4.0;
const Float_t RdA_HT1_HT2 = 7.0;
const Float_t RdA_HT2 = 17.0;

const Float_t RdA_eta_MB = 4.0;
const Float_t RdA_eta_MB_HT1 = 4.0;
const Float_t RdA_eta_HT1_HT2 = 7.0;
const Float_t RdA_eta_HT2 = 14.0;


	TAnalysisSettingsTrigger analysisSettingsSpectrumMB;
	analysisSettingsSpectrumMB.lowPt = 1.0;
	analysisSettingsSpectrumMB.highPt = 5.0;
	analysisSettingsSpectrumMB.lowPtUse = dAu_MB;
	analysisSettingsSpectrumMB.highPtUse = dAu_MB_HT1;
	analysisSettingsSpectrumMB.lowPtBinSize = 0.5;
	analysisSettingsSpectrumMB.midPtBinSize = 1.0;
	analysisSettingsSpectrumMB.highPtBinSize = 2.0;
	analysisSettingsSpectrumMB.switchBinSizePt1 = 3.0;
	analysisSettingsSpectrumMB.switchBinSizePt2 = 5.0;
	analysisSettingsSpectrumMB.rebinNumBinsWeight = 1.0;
	analysisSettingsSpectrumMB.rebinMeanErrorWeight = 2.0;
	analysisSettingsSpectrumMB.setErrorOnZeroBins = false;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsBg = false;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsBgRandom = false;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsBgBack = false;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsBgLowMass = false;//true;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsPeakShape = false;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsPeakShapeEta = false;
	analysisSettingsSpectrumMB.setErrorOnZeroBinsSigma = 0;
	analysisSettingsSpectrumMB.fitLeft = fitLeft;
	analysisSettingsSpectrumMB.fitRight = fitRight;
	analysisSettingsSpectrumMB.fitLeftEta = fitLeftEta;
	analysisSettingsSpectrumMB.fitRightEta = fitRightEta;
	analysisSettingsSpectrumMB.peakLeft = TFitSettings(0.09, 0.002);
	analysisSettingsSpectrumMB.peakRight = TFitSettings(0.16, 0.02);
	analysisSettingsSpectrumMB.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.286683e-01) + (8.539953e-03 * x) + (-1.100295e-03 * x * x) + (4.733651e-05 * x * x *x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumMB.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(3.061734e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumMB.peakLeftEta = TFitSettings(0.3, 0.02);
	analysisSettingsSpectrumMB.peakRightEta = TFitSettings(0.8, 0.02);
	analysisSettingsSpectrumMB.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0, 100, 0, 100, "((5.809154e-01) + (1.123266e-03 * x)) + ([0]*0)", "RQN");
	analysisSettingsSpectrumMB.peakWidthEta = TFitSettings(0.03, 0, 0.03 * 0.5, 0.03 * 2.0, false, 0, 100, 0, 100, "(9.201360e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumMB.peakPositionTrue = truePionMass;
	analysisSettingsSpectrumMB.peakPositionTrueEta = trueEtaMass;
	analysisSettingsSpectrumMB.fitParameter6 = TFitSettings(0, 0, 0, 10000);
	analysisSettingsSpectrumMB.fitParameter7 = TFitSettings();
	analysisSettingsSpectrumMB.fitParameters = TFitSettings();
	analysisSettingsSpectrumMB.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumMB.fitOption = "RQN";
	analysisSettingsSpectrumMB.areaFractionBg = TFitSettings(-1);
	analysisSettingsSpectrumMB.usePeakRangeYield = true;
	analysisSettingsSpectrumMB.peakLeftYield = TFitSettings(0.075, 0.025/15.0);
	analysisSettingsSpectrumMB.peakRightYield = TFitSettings(0.25, 0.05/15.0);
	analysisSettingsSpectrumMB.usePeakRangeSigmaYield = false;
	analysisSettingsSpectrumMB.peakLeftSigmaYield = TFitSettings(-3.0, 0);
	analysisSettingsSpectrumMB.peakRightSigmaYield = TFitSettings(+4.0, 0.05);
	analysisSettingsSpectrumMB.usePeakRangeYieldEta = true;
	analysisSettingsSpectrumMB.peakLeftYieldEta = TFitSettings(0.35, 0.05/15.0);
	analysisSettingsSpectrumMB.peakRightYieldEta = TFitSettings(0.75, 0.0);
	analysisSettingsSpectrumMB.usePeakRangeSigmaYieldEta = false;
	analysisSettingsSpectrumMB.peakLeftSigmaYieldEta = TFitSettings(-2.0, 0);
	analysisSettingsSpectrumMB.peakRightSigmaYieldEta = TFitSettings(+3.0, 0.05);
	analysisSettingsSpectrumMB.usePeakRangeAreaShare = false;
	analysisSettingsSpectrumMB.signalToBackgroundRatio = TFitSettings(-1);
	analysisSettingsSpectrumMB.signalToBackgroundRatioEta = TFitSettings(-1);
	analysisSettingsSpectrumMB.correctCpv = true;
	analysisSettingsSpectrumMB.correctCpvEta = true;
	analysisSettingsSpectrumMB.rebinNumBinsWeightCpv = 1.0;
	analysisSettingsSpectrumMB.rebinMeanErrorWeightCpv = 1.0;
	analysisSettingsSpectrumMB.fitLeftCpv = TFitSettings(0);
	analysisSettingsSpectrumMB.fitRightCpv = TFitSettings(0.7);
	analysisSettingsSpectrumMB.peakLeftCpv = TFitSettings(0.0);
	analysisSettingsSpectrumMB.peakRightCpv = TFitSettings(-0.08);
	analysisSettingsSpectrumMB.peakPositionCpv = TFitSettings(0.0, 0.0, 0.0, 0.0);
	analysisSettingsSpectrumMB.peakWidthCpv = TFitSettings(0.03, 0.0, 0.006, 0.07);
	analysisSettingsSpectrumMB.fitParameter6Cpv = TFitSettings(1000, 0, 0, 1000000);
	analysisSettingsSpectrumMB.fitParameter7Cpv = TFitSettings(9.0, 0.0, 3.0, 20.0);
	analysisSettingsSpectrumMB.fitParametersCpv = TFitSettings();
	analysisSettingsSpectrumMB.fitFunctionStrCpv = fitFuncTrackDistPeak;
	analysisSettingsSpectrumMB.fitOptionCpv = "RQN";
	analysisSettingsSpectrumMB.binRangeLeftCpv = 0.0;
	analysisSettingsSpectrumMB.binRangeRightCpv = 0.7;
	analysisSettingsSpectrumMB.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 1, 3, 0, 100, "pol0", "RQN");
	analysisSettingsSpectrumMB.cpvCorrectionMult = TFitSettings(0, 0, 0, 0, true, 0, 0, 0, 0, "[0]*y*TMath::Exp(-([1]+[2]*x)*y)", "RN");
	analysisSettingsSpectrumMB.fitJetBgFractionFromCandidates = false;
	analysisSettingsSpectrumMB.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsSpectrumMB.getJetBgFractionFromCandidatesPtFit = false;
	analysisSettingsSpectrumMB.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 10, "([0]*0) + (3.288739e-02 - 0.15) + (6.936733e-02 * 1.4 * x)", "RQN");
	analysisSettingsSpectrumMB.bgToSigbg = TFitSettings(0.5);
	analysisSettingsSpectrumMB.smearOutliersSigma = 0;
	analysisSettingsSpectrumMB.smoothFitOption = "QN";
	analysisSettingsSpectrumMB.smoothBg = 0;
	analysisSettingsSpectrumMB.mixNormFixed = 0.0;
	analysisSettingsSpectrumMB.mixNormMassRange = false;
	analysisSettingsSpectrumMB.mixNormLeft = 1.0;
	analysisSettingsSpectrumMB.mixNormRight = 4.0;
	analysisSettingsSpectrumMB.mixNormEntries = true;
	analysisSettingsSpectrumMB.smoothBgRandom = 0.03;
	analysisSettingsSpectrumMB.mixNormFixedRandom = 0.0;
	analysisSettingsSpectrumMB.mixNormMassRangeRandom = false;
	analysisSettingsSpectrumMB.mixNormLeftRandom = 0.9;
	analysisSettingsSpectrumMB.mixNormRightRandom = 1.2;
	analysisSettingsSpectrumMB.mixNormEntriesRandom = true;
	analysisSettingsSpectrumMB.smoothBgBack = 0;
	analysisSettingsSpectrumMB.mixNormFixedBack = 0.0;
	analysisSettingsSpectrumMB.mixNormMassRangeBack = true;
	analysisSettingsSpectrumMB.mixNormLeftBack = 2.0;
	analysisSettingsSpectrumMB.mixNormRightBack = 5.0;
	analysisSettingsSpectrumMB.mixNormEntriesBack = false;
	analysisSettingsSpectrumMB.subtractBg = false;
	analysisSettingsSpectrumMB.subtractBgRandom = false;
	analysisSettingsSpectrumMB.subtractBgBack = false;
	analysisSettingsSpectrumMB.subtractBgLowMass = true;
	analysisSettingsSpectrumMB.subtractBgLowMassPeakShape = true;
	analysisSettingsSpectrumMB.subtractBgLowMassPeakShapeEta = true;
	analysisSettingsSpectrumMB.subtractBgFitDistribution = true;
	analysisSettingsSpectrumMB.useEventMixingBgRandom = true;
	analysisSettingsSpectrumMB.useEventMixingNotmatchedBgRandom = false;
	analysisSettingsSpectrumMB.useEventJetMixingBg = true;
	analysisSettingsSpectrumMB.useEventJetNotmatchedMixingBg = false;
	analysisSettingsSpectrumMB.useEventJetBackMixingBg = false;
	analysisSettingsSpectrumMB.useEventShufflingBg = false;
	analysisSettingsSpectrumMB.fitBgDistribution = true;
	analysisSettingsSpectrumMB.fitBgRandomDistribution = true;
	analysisSettingsSpectrumMB.fitBgLowMassDistribution = false;
	analysisSettingsSpectrumMB.fitPeakShapeDistribution = false;
	analysisSettingsSpectrumMB.fitPeakShapeEtaDistribution = false;
	analysisSettingsSpectrumMB.fitPeakShapeEtabgDistribution = false;
	analysisSettingsSpectrumMB.fitDistributionLeft = 0.3;//0.3; // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsSpectrumMB.fitDistributionRight = 0.4;//0.4; // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsSpectrumMB.fitDistributionLeft2 = 0.9;
	analysisSettingsSpectrumMB.fitDistributionRight2 = 1.2;
	analysisSettingsSpectrumMB.fitDistributionOption = "RQN UNORM";
	analysisSettingsSpectrumMB.lowNormFixed = TFitSettings(1.12, -0.19); // towerMB_smd2HT1
	//analysisSettingsSpectrumMB.lowNormFixed = TFitSettings(1.12, -0.09); // towerMB_smd2HT1_gammaetacut
	//analysisSettingsSpectrumMB.lowNormFixed = TFitSettings(1);
	//analysisSettingsSpectrumMB.lowNormFixed = TFitSettings(0.18, 0, 0, 1, true, 0, 100, 0, 2.5, "(TMath::Exp((1.28365)-(1.07833*x))) + (0*[0])", "RQN");
	analysisSettingsSpectrumMB.lowNormMassRange = false;
	analysisSettingsSpectrumMB.lowNormMassRangeLeft = 0.04;
	analysisSettingsSpectrumMB.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumMB.lowNormPoints = false;
	analysisSettingsSpectrumMB.lowmassbgPointsLowPt = 1.0;
	analysisSettingsSpectrumMB.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumMB.lowNormPointsPt = true;
	analysisSettingsSpectrumMB.lowmassbgNbar = TFitSettings(1.0, 0, 0, 1, false, 0, 100, 0, 100, "(1.0/((1.0/( ( 0.3*x*TMath::Exp(-0.3233*x) / ((7.334427e-01)+(-9.284123e-02*x)+(8.582639e-03*x*x)+(-2.596326e-04*x*x*x)) ) * (1-0.2+(0.2*0.358/0.639)) / (1.0+(2.48210e-03*x*x)) ))+1.0)) + (0*[0])", "RQN");
	analysisSettingsSpectrumMB.fixedPrescale = 1.0;
	analysisSettingsSpectrumMB.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumMB.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumMB.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumMB.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumMB.prescalePointsLowPt = 0.0;
	analysisSettingsSpectrumMB.prescalePointsHighPt = 2.0;
	analysisSettingsSpectrumMB.prescalePointsLowPtIntegral = 0.0;
	analysisSettingsSpectrumMB.prescalePointsHighPtIntegral = 2.0;
	analysisSettingsSpectrumMB.useProportionalBgErrors = false;
	analysisSettingsSpectrumMB.useIndependentSigBgErrors = true;
	analysisSettingsSpectrumMB.normalizeSigBgErrors = false;
	analysisSettingsSpectrumMB.normalizeSigBgSquareErrors = false;
	analysisSettingsSpectrumMB.useAnticorrelatedSigBgErrors = false;
	analysisSettingsSpectrumMB.normPeakShapeYieldShow = true;
	analysisSettingsSpectrumMB.normPeakShapeEtaYieldShow = true;
	analysisSettingsSpectrumMB.correctVertexFindingEff = true;
	analysisSettingsSpectrumMB.vertexFindingEff = vertEffDAuMinBias;
	analysisSettingsSpectrumMB.correctVertexFindingEffMB = false;
	analysisSettingsSpectrumMB.showBins = true;
	analysisSettingsSpectrumMB.showSigBg = true;
	analysisSettingsSpectrumMB.showSig = true;
	analysisSettingsSpectrumMB.showPeakShape = false;
	analysisSettingsSpectrumMB.showPeakShapeEta = false;
	analysisSettingsSpectrumMB.showPeakShapeEtabg = false;
	analysisSettingsSpectrumMB.showBg = false;
	analysisSettingsSpectrumMB.showBgRandom = false;
	analysisSettingsSpectrumMB.showBgBack = false;
	analysisSettingsSpectrumMB.showBgLowMass = true;
	analysisSettingsSpectrumMB.showBgFitDistribution = true;
	analysisSettingsSpectrumMB.showFuncPeak = true;
	analysisSettingsSpectrumMB.showFuncPeakEta = true;
	analysisSettingsSpectrumMB.showFuncBg = false;
	analysisSettingsSpectrumMB.showResidual = false;
	analysisSettingsSpectrumMB.showResidualBg = false;
	analysisSettingsSpectrumMB.showResidualPeakShapeBg = false;
	analysisSettingsSpectrumMB.showText = true;
	analysisSettingsSpectrumMB.showTruePeakPos = false;
	analysisSettingsSpectrumMB.color = kBlack;
	analysisSettingsSpectrumMB.colorSigBg = kBlack;
	analysisSettingsSpectrumMB.colorBg = 6;//kGreen;
	analysisSettingsSpectrumMB.colorBgRandom = 15;
	analysisSettingsSpectrumMB.colorBgBack = 6;
	analysisSettingsSpectrumMB.colorBgLowMass = 6; // hot pink
	analysisSettingsSpectrumMB.colorBgFitDistribution = 11;
	analysisSettingsSpectrumMB.peakLinesColor = kGreen;
	analysisSettingsSpectrumMB.fitLinesColor = kGreen;
	analysisSettingsSpectrumMB.bgColor = kGreen;
	analysisSettingsSpectrumMB.fillColor = 17;
	analysisSettingsSpectrumMB.bgFillStyle = 1001;
	analysisSettingsSpectrumMB.binRangeLeft = 0.03;
	analysisSettingsSpectrumMB.binRangeRight = 1.17;
	analysisSettingsSpectrumMB.name = "MB";
	analysisSettingsSpectrumMB.title = "MinBias";
	analysisSettingsSpectrumMB.useXerr = false;
	analysisSettingsSpectrumMB.triggerThreshold = 0.0;
	analysisSettingsSpectrumMB.useFittedPeakPos = false;
	analysisSettingsSpectrumMB.showInvBinsFitsSeparately = false;
        analysisSettingsSpectrumMB.drawOptions    = TDrawOptions(TNamed("MB", "MinBias"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullCircle, 1.5), "P", "P", TString("MinBias"));
        analysisSettingsSpectrumMB.drawOptionsEta = TDrawOptions(TNamed("MB", "MinBias"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleUp, 1.5), "P", "P", TString("MinBias #eta"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumHT1;
	analysisSettingsSpectrumHT1.lowPt = 4.0;
	analysisSettingsSpectrumHT1.highPt = 12.0;
	analysisSettingsSpectrumHT1.lowPtUse = dAu_MB_HT1;
	analysisSettingsSpectrumHT1.highPtUse = dAu_HT1_HT2;
	analysisSettingsSpectrumHT1.lowPtBinSize = 1.0;
	analysisSettingsSpectrumHT1.midPtBinSize = 2.0;
	analysisSettingsSpectrumHT1.highPtBinSize = 2.0;
	analysisSettingsSpectrumHT1.switchBinSizePt1 = 8.0;
	analysisSettingsSpectrumHT1.switchBinSizePt2 = 12.0;
	analysisSettingsSpectrumHT1.rebinNumBinsWeight = 1.0;
	analysisSettingsSpectrumHT1.rebinMeanErrorWeight = 1.5;
	analysisSettingsSpectrumHT1.setErrorOnZeroBins = false;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsBg = false;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsBgRandom = false;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsBgBack = false;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsBgLowMass = false;//true;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsPeakShape = false;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsPeakShapeEta = false;
	analysisSettingsSpectrumHT1.setErrorOnZeroBinsSigma = -0.04;
	analysisSettingsSpectrumHT1.fitLeft = fitLeft;
	analysisSettingsSpectrumHT1.fitRight = fitRight;
	analysisSettingsSpectrumHT1.fitLeftEta = fitLeftEta;
	analysisSettingsSpectrumHT1.fitRightEta = fitRightEta;
	analysisSettingsSpectrumHT1.peakLeft = TFitSettings(0.05, 0.002);
	analysisSettingsSpectrumHT1.peakRight = TFitSettings(0.2, 0.02);
	analysisSettingsSpectrumHT1.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.286683e-01) + (8.539953e-03 * x) + (-1.100295e-03 * x * x) + (4.733651e-05 * x * x *x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT1.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(3.061734e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT1.peakLeftEta = TFitSettings(0.2, 0.02);
	analysisSettingsSpectrumHT1.peakRightEta = TFitSettings(0.8, 0.02);
	analysisSettingsSpectrumHT1.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0, 100, 0, 100, "((5.759024e-01) + (1.209403e-03 * x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT1.peakWidthEta = TFitSettings(0.10, 0, 0.10 * 0.5, 0.10 * 1.2, false, 0, 100, 0, 100, "(9.201360e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT1.peakPositionTrue = truePionMass;
	analysisSettingsSpectrumHT1.peakPositionTrueEta = trueEtaMass;
	analysisSettingsSpectrumHT1.fitParameter6 = TFitSettings(0.3, 0, 0, 1.0);
	analysisSettingsSpectrumHT1.fitParameter7 = TFitSettings();
	analysisSettingsSpectrumHT1.fitParameters = TFitSettings();
	analysisSettingsSpectrumHT1.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumHT1.fitOption = "RQN";
	analysisSettingsSpectrumHT1.areaFractionBg = TFitSettings(-1);
	analysisSettingsSpectrumHT1.usePeakRangeYield = true;
	analysisSettingsSpectrumHT1.peakLeftYield = TFitSettings(0.075, 0.025/15.0);
	analysisSettingsSpectrumHT1.peakRightYield = TFitSettings(0.25, 0.05/15.0);
	analysisSettingsSpectrumHT1.usePeakRangeSigmaYield = false;
	analysisSettingsSpectrumHT1.peakLeftSigmaYield = TFitSettings(-2.5, 0);
	analysisSettingsSpectrumHT1.peakRightSigmaYield = TFitSettings(+4.0, 0.05);
	analysisSettingsSpectrumHT1.usePeakRangeYieldEta = true;
	analysisSettingsSpectrumHT1.peakLeftYieldEta = TFitSettings(0.35, 0.05/15.0);
	analysisSettingsSpectrumHT1.peakRightYieldEta = TFitSettings(0.75, 0.0);
	analysisSettingsSpectrumHT1.usePeakRangeSigmaYieldEta = false;
	analysisSettingsSpectrumHT1.peakLeftSigmaYieldEta = TFitSettings(-2.0, 0);
	analysisSettingsSpectrumHT1.peakRightSigmaYieldEta = TFitSettings(+3.0, 0.05);
	analysisSettingsSpectrumHT1.usePeakRangeAreaShare = false;
	analysisSettingsSpectrumHT1.signalToBackgroundRatio = TFitSettings(-1);
	analysisSettingsSpectrumHT1.signalToBackgroundRatioEta = TFitSettings(-1);
	analysisSettingsSpectrumHT1.correctCpv = true;
	analysisSettingsSpectrumHT1.correctCpvEta = true;
	analysisSettingsSpectrumHT1.rebinNumBinsWeightCpv = 1.0;
	analysisSettingsSpectrumHT1.rebinMeanErrorWeightCpv = 1.0;
	analysisSettingsSpectrumHT1.fitLeftCpv = TFitSettings(0);
	analysisSettingsSpectrumHT1.fitRightCpv = TFitSettings(0.7);
	analysisSettingsSpectrumHT1.peakLeftCpv = TFitSettings(0.0);
	analysisSettingsSpectrumHT1.peakRightCpv = TFitSettings(-0.08);
	analysisSettingsSpectrumHT1.peakPositionCpv = TFitSettings(0.0, 0.0, 0.0, 0.0);
	analysisSettingsSpectrumHT1.peakWidthCpv = TFitSettings(0.03, 0.0, 0.006, 0.07);
	analysisSettingsSpectrumHT1.fitParameter6Cpv = TFitSettings(1000, 0, 0, 1000000);
	analysisSettingsSpectrumHT1.fitParameter7Cpv = TFitSettings(9.0, 0.0, 3.0, 20.0);
	analysisSettingsSpectrumHT1.fitParametersCpv = TFitSettings();
	analysisSettingsSpectrumHT1.fitFunctionStrCpv = fitFuncTrackDistPeak;
	analysisSettingsSpectrumHT1.fitOptionCpv = "RQN";
	analysisSettingsSpectrumHT1.binRangeLeftCpv = 0.0;
	analysisSettingsSpectrumHT1.binRangeRightCpv = 0.7;
	analysisSettingsSpectrumHT1.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 0, 10, 0, 100, "pol1", "RQNW");
	analysisSettingsSpectrumHT1.cpvCorrectionMult = TFitSettings(0, 0, 0, 0, true, 0, 0, 0, 0, "[0]*y*TMath::Exp(-([1]+[2]*x)*y)", "RN");
	analysisSettingsSpectrumHT1.fitJetBgFractionFromCandidates = false;
	analysisSettingsSpectrumHT1.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsSpectrumHT1.getJetBgFractionFromCandidatesPtFit = false;
	analysisSettingsSpectrumHT1.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 10, "([0]*0) + (3.288739e-02 - 0.15) + (6.936733e-02 * 1.4 * x)", "RQN");
	analysisSettingsSpectrumHT1.bgToSigbg = TFitSettings(0.5);
	analysisSettingsSpectrumHT1.smearOutliersSigma = 0;
	analysisSettingsSpectrumHT1.smoothFitOption = "QN";
	analysisSettingsSpectrumHT1.smoothBg = 0;
	analysisSettingsSpectrumHT1.mixNormFixed = 0.0;
	analysisSettingsSpectrumHT1.mixNormMassRange = false;
	analysisSettingsSpectrumHT1.mixNormLeft = 1.0;
	analysisSettingsSpectrumHT1.mixNormRight = 4.0;
	analysisSettingsSpectrumHT1.mixNormEntries = true;
	analysisSettingsSpectrumHT1.smoothBgRandom = 0;
	analysisSettingsSpectrumHT1.mixNormFixedRandom = 0.0;
	analysisSettingsSpectrumHT1.mixNormMassRangeRandom = false;
	analysisSettingsSpectrumHT1.mixNormLeftRandom = 0.9;
	analysisSettingsSpectrumHT1.mixNormRightRandom = 1.2;
	analysisSettingsSpectrumHT1.mixNormEntriesRandom = true;
	analysisSettingsSpectrumHT1.smoothBgBack = 0;
	analysisSettingsSpectrumHT1.mixNormFixedBack = 0.0;
	analysisSettingsSpectrumHT1.mixNormMassRangeBack = true;
	analysisSettingsSpectrumHT1.mixNormLeftBack = 2.0;
	analysisSettingsSpectrumHT1.mixNormRightBack = 5.0;
	analysisSettingsSpectrumHT1.mixNormEntriesBack = false;
	analysisSettingsSpectrumHT1.subtractBg = false;
	analysisSettingsSpectrumHT1.subtractBgRandom = false;
	analysisSettingsSpectrumHT1.subtractBgBack = false;
	analysisSettingsSpectrumHT1.subtractBgLowMass = true;
	analysisSettingsSpectrumHT1.subtractBgLowMassPeakShape = true;
	analysisSettingsSpectrumHT1.subtractBgLowMassPeakShapeEta = true;
	analysisSettingsSpectrumHT1.subtractBgFitDistribution = true;
	analysisSettingsSpectrumHT1.useEventMixingBgRandom = true;
	analysisSettingsSpectrumHT1.useEventMixingNotmatchedBgRandom = false;
	analysisSettingsSpectrumHT1.useEventJetMixingBg = true;
	analysisSettingsSpectrumHT1.useEventJetNotmatchedMixingBg = false;
	analysisSettingsSpectrumHT1.useEventJetBackMixingBg = false;
	analysisSettingsSpectrumHT1.useEventShufflingBg = false;
	analysisSettingsSpectrumHT1.fitBgDistribution = true;
	analysisSettingsSpectrumHT1.fitBgRandomDistribution = true;
	analysisSettingsSpectrumHT1.fitBgLowMassDistribution = false;
	analysisSettingsSpectrumHT1.fitPeakShapeDistribution = false;
	analysisSettingsSpectrumHT1.fitPeakShapeEtaDistribution = false;
	analysisSettingsSpectrumHT1.fitPeakShapeEtabgDistribution = false;
	analysisSettingsSpectrumHT1.fitDistributionLeft = 0.3;
	analysisSettingsSpectrumHT1.fitDistributionRight = 0.4;
	analysisSettingsSpectrumHT1.fitDistributionLeft2 = 0.9;
	analysisSettingsSpectrumHT1.fitDistributionRight2 = 1.2;
	analysisSettingsSpectrumHT1.fitDistributionOption = "RQN UNORM";
	analysisSettingsSpectrumHT1.lowNormFixed = TFitSettings(1.0);
	analysisSettingsSpectrumHT1.lowNormMassRange = false;
	analysisSettingsSpectrumHT1.lowNormMassRangeLeft = 0.05;
	analysisSettingsSpectrumHT1.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumHT1.lowNormPoints = false;
	analysisSettingsSpectrumHT1.lowmassbgPointsLowPt = 4.0;
	analysisSettingsSpectrumHT1.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumHT1.lowNormPointsPt = true;
	analysisSettingsSpectrumHT1.lowmassbgNbar = TFitSettings(0.3, 0, 0, 1, true, 0, 100, 0, 100, "(1.0/((1.0/((0.3*x*TMath::Exp(-0.3233*x) / ((7.334427e-01)+(-9.284123e-02*x)+(8.582639e-03*x*x)+(-2.596326e-04*x*x*x))) * (1-0.2+(0.2*0.358/0.639)) * (10.39*TMath::Landau(x,3.214,0.695)) / (1.0+(2.48210e-03*x*x)) ))+1.0)) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT1.fixedPrescale = /*1547.33;//*/2873.35; // towerMB_smd2HT1_gammaetacut / towerMB_smd2HT1
	analysisSettingsSpectrumHT1.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumHT1.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumHT1.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumHT1.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumHT1.prescalePointsLowPt = 4.0;
	analysisSettingsSpectrumHT1.prescalePointsHighPt = 10.0;
	analysisSettingsSpectrumHT1.prescalePointsLowPtIntegral = 5.2;
	analysisSettingsSpectrumHT1.prescalePointsHighPtIntegral = 6.5;
	analysisSettingsSpectrumHT1.useProportionalBgErrors = false;
	analysisSettingsSpectrumHT1.useIndependentSigBgErrors = true;
	analysisSettingsSpectrumHT1.normalizeSigBgErrors = false;
	analysisSettingsSpectrumHT1.normalizeSigBgSquareErrors = false;
	analysisSettingsSpectrumHT1.useAnticorrelatedSigBgErrors = false;
	analysisSettingsSpectrumHT1.normPeakShapeYieldShow = true;
	analysisSettingsSpectrumHT1.normPeakShapeEtaYieldShow = true;
	analysisSettingsSpectrumHT1.correctVertexFindingEff = false;
	analysisSettingsSpectrumHT1.vertexFindingEff = 1.0;
	analysisSettingsSpectrumHT1.correctVertexFindingEffMB = true;
	analysisSettingsSpectrumHT1.showBins = true;
	analysisSettingsSpectrumHT1.showSigBg = true;
	analysisSettingsSpectrumHT1.showSig = true;
	analysisSettingsSpectrumHT1.showPeakShape = false;
	analysisSettingsSpectrumHT1.showPeakShapeEta = false;
	analysisSettingsSpectrumHT1.showPeakShapeEtabg = false;
	analysisSettingsSpectrumHT1.showBg = false;
	analysisSettingsSpectrumHT1.showBgRandom = false;
	analysisSettingsSpectrumHT1.showBgBack = false;
	analysisSettingsSpectrumHT1.showBgLowMass = true;
	analysisSettingsSpectrumHT1.showBgFitDistribution = true;
	analysisSettingsSpectrumHT1.showFuncPeak = true;
	analysisSettingsSpectrumHT1.showFuncPeakEta = true;
	analysisSettingsSpectrumHT1.showFuncBg = false;
	analysisSettingsSpectrumHT1.showResidual = false;
	analysisSettingsSpectrumHT1.showResidualBg = false;
	analysisSettingsSpectrumHT1.showResidualPeakShapeBg = false;
	analysisSettingsSpectrumHT1.showText = true;
	analysisSettingsSpectrumHT1.showTruePeakPos = false;
	analysisSettingsSpectrumHT1.color = kBlack;//kBlue;
	analysisSettingsSpectrumHT1.colorSigBg = kBlack;//kBlue;
	analysisSettingsSpectrumHT1.colorBg = 6;//kGreen;
	analysisSettingsSpectrumHT1.colorBgRandom = 15;
	analysisSettingsSpectrumHT1.colorBgBack = 6;
	analysisSettingsSpectrumHT1.colorBgLowMass = 6; // hot pink
	analysisSettingsSpectrumHT1.colorBgFitDistribution = 11;
	analysisSettingsSpectrumHT1.peakLinesColor = kGreen;
	analysisSettingsSpectrumHT1.fitLinesColor = kGreen;
	analysisSettingsSpectrumHT1.bgColor = kGreen;
	analysisSettingsSpectrumHT1.fillColor = 17;
	analysisSettingsSpectrumHT1.bgFillStyle = 1001;
	analysisSettingsSpectrumHT1.binRangeLeft = 0.03;
	analysisSettingsSpectrumHT1.binRangeRight = 1.17;
	analysisSettingsSpectrumHT1.name = "HT1";
	analysisSettingsSpectrumHT1.title = "HighTower-1";
	analysisSettingsSpectrumHT1.useXerr = false;
	analysisSettingsSpectrumHT1.triggerThreshold = 2.5 + 1.0;
	analysisSettingsSpectrumHT1.useFittedPeakPos = false;
	analysisSettingsSpectrumHT1.showInvBinsFitsSeparately = false;
        analysisSettingsSpectrumHT1.drawOptions    = TDrawOptions(TNamed("HT1", "HighTower-1"), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullCircle, 1.5), "P", "P", TString("HighTower-1"));
        analysisSettingsSpectrumHT1.drawOptionsEta = TDrawOptions(TNamed("HT1", "HighTower-1"), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullTriangleUp, 1.5), "P", "P", TString("HighTower-1 #eta"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumHT2;
	analysisSettingsSpectrumHT2.lowPt = 6.0;
	analysisSettingsSpectrumHT2.highPt = 17.0;
	analysisSettingsSpectrumHT2.lowPtUse = dAu_HT1_HT2;
	analysisSettingsSpectrumHT2.highPtUse = dAu_HT2;
	analysisSettingsSpectrumHT2.lowPtBinSize = 1.0;
	analysisSettingsSpectrumHT2.midPtBinSize = 2.0;
	analysisSettingsSpectrumHT2.highPtBinSize = 3.0;
	analysisSettingsSpectrumHT2.switchBinSizePt1 = 8.0;
	analysisSettingsSpectrumHT2.switchBinSizePt2 = 14.0;
	analysisSettingsSpectrumHT2.rebinNumBinsWeight = 1.0;
	analysisSettingsSpectrumHT2.rebinMeanErrorWeight = 1.5;
	analysisSettingsSpectrumHT2.setErrorOnZeroBins = false;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsBg = false;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsBgRandom = false;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsBgBack = false;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsBgLowMass = false;//true;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsPeakShape = false;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsPeakShapeEta = false;
	analysisSettingsSpectrumHT2.setErrorOnZeroBinsSigma = -0.04;
	analysisSettingsSpectrumHT2.fitLeft = fitLeft;
	analysisSettingsSpectrumHT2.fitRight = fitRight;
	analysisSettingsSpectrumHT2.fitLeftEta = fitLeftEta;
	analysisSettingsSpectrumHT2.fitRightEta = fitRightEta;
	analysisSettingsSpectrumHT2.peakLeft = TFitSettings(0.05, 0.002);
	analysisSettingsSpectrumHT2.peakRight = TFitSettings(0.2, 0.02);
	analysisSettingsSpectrumHT2.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.286683e-01) + (8.539953e-03 * x) + (-1.100295e-03 * x * x) + (4.733651e-05 * x * x *x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT2.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(3.061734e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT2.peakLeftEta = TFitSettings(0.1, 0.01);
	analysisSettingsSpectrumHT2.peakRightEta = TFitSettings(0.8, 0.02);
	analysisSettingsSpectrumHT2.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0, 100, 0, 100, "((5.154388e-01) + (9.879051e-03 * x) + (-3.015698e-04 * x * x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT2.peakWidthEta = TFitSettings(0.10, 0, 0.10 * 0.5, 0.10 * 1.2, false, 0, 100, 0, 100, "(9.201360e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT2.peakPositionTrue = truePionMass;
	analysisSettingsSpectrumHT2.peakPositionTrueEta = trueEtaMass;
	analysisSettingsSpectrumHT2.fitParameter6 = TFitSettings(0, 0, 0, 10000);
	analysisSettingsSpectrumHT2.fitParameter7 = TFitSettings();
	analysisSettingsSpectrumHT2.fitParameters = TFitSettings();
	analysisSettingsSpectrumHT2.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumHT2.fitOption = "RQN";
	analysisSettingsSpectrumHT2.areaFractionBg = TFitSettings(-1);
	analysisSettingsSpectrumHT2.usePeakRangeYield = true;
	analysisSettingsSpectrumHT2.peakLeftYield = TFitSettings(0.075, 0.025/15.0);
	analysisSettingsSpectrumHT2.peakRightYield = TFitSettings(0.25, 0.05/15.0);
	analysisSettingsSpectrumHT2.usePeakRangeSigmaYield = false;
	analysisSettingsSpectrumHT2.peakLeftSigmaYield = TFitSettings(-2.0, 0);
	analysisSettingsSpectrumHT2.peakRightSigmaYield = TFitSettings(+4.0, 0.05);
	analysisSettingsSpectrumHT2.usePeakRangeYieldEta = true;
	analysisSettingsSpectrumHT2.peakLeftYieldEta = TFitSettings(0.35, 0.05/15.0);
	analysisSettingsSpectrumHT2.peakRightYieldEta = TFitSettings(0.75, 0.0);
	analysisSettingsSpectrumHT2.usePeakRangeSigmaYieldEta = false;
	analysisSettingsSpectrumHT2.peakLeftSigmaYieldEta = TFitSettings(-2.0, 0);
	analysisSettingsSpectrumHT2.peakRightSigmaYieldEta = TFitSettings(+3.0, 0.05);
	analysisSettingsSpectrumHT2.usePeakRangeAreaShare = false;
	analysisSettingsSpectrumHT2.signalToBackgroundRatio = TFitSettings(-1);
	analysisSettingsSpectrumHT2.signalToBackgroundRatioEta = TFitSettings(-1);
	analysisSettingsSpectrumHT2.correctCpv = true;
	analysisSettingsSpectrumHT2.correctCpvEta = true;
	analysisSettingsSpectrumHT2.rebinNumBinsWeightCpv = 1.0;
	analysisSettingsSpectrumHT2.rebinMeanErrorWeightCpv = 1.0;
	analysisSettingsSpectrumHT2.fitLeftCpv = TFitSettings(0);
	analysisSettingsSpectrumHT2.fitRightCpv = TFitSettings(0.7);
	analysisSettingsSpectrumHT2.peakLeftCpv = TFitSettings(0.0);
	analysisSettingsSpectrumHT2.peakRightCpv = TFitSettings(-0.08);
	analysisSettingsSpectrumHT2.peakPositionCpv = TFitSettings(0.0, 0.0, 0.0, 0.0);
	analysisSettingsSpectrumHT2.peakWidthCpv = TFitSettings(0.03, 0.0, 0.006, 0.07);
	analysisSettingsSpectrumHT2.fitParameter6Cpv = TFitSettings(1000, 0, 0, 1000000);
	analysisSettingsSpectrumHT2.fitParameter7Cpv = TFitSettings(9.0, 0.0, 3.0, 20.0);
	analysisSettingsSpectrumHT2.fitParametersCpv = TFitSettings();
	analysisSettingsSpectrumHT2.fitFunctionStrCpv = fitFuncTrackDistPeak;
	analysisSettingsSpectrumHT2.fitOptionCpv = "RQN";
	analysisSettingsSpectrumHT2.binRangeLeftCpv = 0.0;
	analysisSettingsSpectrumHT2.binRangeRightCpv = 0.7;
	analysisSettingsSpectrumHT2.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 0, 14, 0, 100, "pol1", "RQNW");
	analysisSettingsSpectrumHT2.cpvCorrectionMult = TFitSettings(0, 0, 0, 0, true, 0, 0, 0, 0, "[0]*y*TMath::Exp(-([1]+[2]*x)*y)", "RN");
	analysisSettingsSpectrumHT2.fitJetBgFractionFromCandidates = false;
	analysisSettingsSpectrumHT2.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsSpectrumHT2.getJetBgFractionFromCandidatesPtFit = false;
	analysisSettingsSpectrumHT2.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 10, "([0]*0) + (3.288739e-02 - 0.15) + (6.936733e-02 * 1.4 * x)", "RQN");
	analysisSettingsSpectrumHT2.bgToSigbg = TFitSettings(0.5);
	analysisSettingsSpectrumHT2.smearOutliersSigma = 0;
	analysisSettingsSpectrumHT2.smoothFitOption = "QN";
	analysisSettingsSpectrumHT2.smoothBg = 0;
	analysisSettingsSpectrumHT2.mixNormFixed = 0.0;
	analysisSettingsSpectrumHT2.mixNormMassRange = false;
	analysisSettingsSpectrumHT2.mixNormLeft = 1.0;
	analysisSettingsSpectrumHT2.mixNormRight = 4.0;
	analysisSettingsSpectrumHT2.mixNormEntries = true;
	analysisSettingsSpectrumHT2.smoothBgRandom = 0;
	analysisSettingsSpectrumHT2.mixNormFixedRandom = 0.0;
	analysisSettingsSpectrumHT2.mixNormMassRangeRandom = false;
	analysisSettingsSpectrumHT2.mixNormLeftRandom = 0.9;
	analysisSettingsSpectrumHT2.mixNormRightRandom = 1.2;
	analysisSettingsSpectrumHT2.mixNormEntriesRandom = true;
	analysisSettingsSpectrumHT2.smoothBgBack = 0;
	analysisSettingsSpectrumHT2.mixNormFixedBack = 0.0;
	analysisSettingsSpectrumHT2.mixNormMassRangeBack = true;
	analysisSettingsSpectrumHT2.mixNormLeftBack = 2.0;
	analysisSettingsSpectrumHT2.mixNormRightBack = 5.0;
	analysisSettingsSpectrumHT2.mixNormEntriesBack = false;
	analysisSettingsSpectrumHT2.subtractBg = false;
	analysisSettingsSpectrumHT2.subtractBgRandom = false;
	analysisSettingsSpectrumHT2.subtractBgBack = false;
	analysisSettingsSpectrumHT2.subtractBgLowMass = true;
	analysisSettingsSpectrumHT2.subtractBgLowMassPeakShape = true;
	analysisSettingsSpectrumHT2.subtractBgLowMassPeakShapeEta = true;
	analysisSettingsSpectrumHT2.subtractBgFitDistribution = true;
	analysisSettingsSpectrumHT2.useEventMixingBgRandom = true;
	analysisSettingsSpectrumHT2.useEventJetMixingBg = true;
	analysisSettingsSpectrumHT2.useEventJetNotmatchedMixingBg = false;
	analysisSettingsSpectrumHT2.useEventJetBackMixingBg = false;
	analysisSettingsSpectrumHT2.useEventShufflingBg = false;
	analysisSettingsSpectrumHT2.fitBgDistribution = true;
	analysisSettingsSpectrumHT2.fitBgRandomDistribution = true;
	analysisSettingsSpectrumHT2.fitBgLowMassDistribution = false;
	analysisSettingsSpectrumHT2.fitPeakShapeDistribution = false;
	analysisSettingsSpectrumHT2.fitPeakShapeEtaDistribution = false;
	analysisSettingsSpectrumHT2.fitPeakShapeEtabgDistribution = false;
	analysisSettingsSpectrumHT2.fitDistributionLeft = TFitSettings(0.3, 0, 0, 1, true, 0, 100, 0, 14, "(0.1/9.5) + (x*(0.3-(7.0*0.1/9.5)))", "RQN");
	analysisSettingsSpectrumHT2.fitDistributionRight = TFitSettings(0.4, 0, 0, 1, true, 0, 100, 0, 14, "(0.1/9.5) + (x*(0.4-(7.0*0.1/9.5)))", "RQN");
	analysisSettingsSpectrumHT2.fitDistributionLeft2 = 0.9;
	analysisSettingsSpectrumHT2.fitDistributionRight2 = 1.2; // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut = 1.2 / 1.0
	analysisSettingsSpectrumHT2.fitDistributionOption = "RQN UNORM";
	analysisSettingsSpectrumHT2.lowNormFixed = TFitSettings(1.0);
	analysisSettingsSpectrumHT2.lowNormMassRange = false;
	analysisSettingsSpectrumHT2.lowNormMassRangeLeft = 0.05;
	analysisSettingsSpectrumHT2.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumHT2.lowNormPoints = false;
	analysisSettingsSpectrumHT2.lowmassbgPointsLowPt = 6.0;
	analysisSettingsSpectrumHT2.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumHT2.lowNormPointsPt = true;
	analysisSettingsSpectrumHT2.lowmassbgNbar = TFitSettings(0.3, 0, 0, 1, true, 0, 100, 0, 100, "(1.0/((1.0/((0.3*x*TMath::Exp(-0.3233*x) / ((7.334427e-01)+(-9.284123e-02*x)+(8.582639e-03*x*x)+(-2.596326e-04*x*x*x))) * (1-0.2+(0.2*0.358/0.639)) * (2.886*TMath::Landau(x,5.595,0.9)) / (1.0+(2.48210e-03*x*x)) ))+1.0)) + (0*[0])", "RQN");
	analysisSettingsSpectrumHT2.fixedPrescale = (28032.5 + 29070.3)/2.0;
	analysisSettingsSpectrumHT2.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumHT2.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumHT2.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumHT2.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumHT2.prescalePointsLowPt = 7.5;
	analysisSettingsSpectrumHT2.prescalePointsHighPt = 19.0;
	analysisSettingsSpectrumHT2.prescalePointsLowPtIntegral = 10.0;
	analysisSettingsSpectrumHT2.prescalePointsHighPtIntegral = 12.0;
	analysisSettingsSpectrumHT2.useProportionalBgErrors = false;
	analysisSettingsSpectrumHT2.useIndependentSigBgErrors = true;
	analysisSettingsSpectrumHT2.normalizeSigBgErrors = false;
	analysisSettingsSpectrumHT2.normalizeSigBgSquareErrors = false;
	analysisSettingsSpectrumHT2.useAnticorrelatedSigBgErrors = false;
	analysisSettingsSpectrumHT2.normPeakShapeYieldShow = true;
	analysisSettingsSpectrumHT2.normPeakShapeEtaYieldShow = true;
	analysisSettingsSpectrumHT2.correctVertexFindingEff = false;
	analysisSettingsSpectrumHT2.vertexFindingEff = 1.0;
	analysisSettingsSpectrumHT2.correctVertexFindingEffMB = true;
	analysisSettingsSpectrumHT2.showBins = true;
	analysisSettingsSpectrumHT2.showSigBg = true;
	analysisSettingsSpectrumHT2.showSig = true;
	analysisSettingsSpectrumHT2.showPeakShape = false;
	analysisSettingsSpectrumHT2.showPeakShapeEta = false;
	analysisSettingsSpectrumHT2.showPeakShapeEtabg = false;
	analysisSettingsSpectrumHT2.showBg = false;
	analysisSettingsSpectrumHT2.showBgRandom = false;
	analysisSettingsSpectrumHT2.showBgBack = false;
	analysisSettingsSpectrumHT2.showBgLowMass = true;
	analysisSettingsSpectrumHT2.showBgFitDistribution = true;
	analysisSettingsSpectrumHT2.showFuncPeak = true;
	analysisSettingsSpectrumHT2.showFuncPeakEta = true;
	analysisSettingsSpectrumHT2.showFuncBg = false;
	analysisSettingsSpectrumHT2.showResidual = false;
	analysisSettingsSpectrumHT2.showResidualBg = false;
	analysisSettingsSpectrumHT2.showResidualPeakShapeBg = false;
	analysisSettingsSpectrumHT2.showText = true;
	analysisSettingsSpectrumHT2.showTruePeakPos = false;
	analysisSettingsSpectrumHT2.color = kBlack;//kRed;
	analysisSettingsSpectrumHT2.colorSigBg = kBlack;//kRed;
	analysisSettingsSpectrumHT2.colorBg = 6;//kGreen;
	analysisSettingsSpectrumHT2.colorBgRandom = 15;
	analysisSettingsSpectrumHT2.colorBgBack = 6;
	analysisSettingsSpectrumHT2.colorBgLowMass = 6; // hot pink
	analysisSettingsSpectrumHT2.colorBgFitDistribution = 11;
	analysisSettingsSpectrumHT2.peakLinesColor = kGreen;
	analysisSettingsSpectrumHT2.fitLinesColor = kGreen;
	analysisSettingsSpectrumHT2.bgColor = kGreen;
	analysisSettingsSpectrumHT2.fillColor = 17;
	analysisSettingsSpectrumHT2.bgFillStyle = 1001;
	analysisSettingsSpectrumHT2.binRangeLeft = 0.03;
	analysisSettingsSpectrumHT2.binRangeRight = 1.17;
	analysisSettingsSpectrumHT2.name = "HT2";
	analysisSettingsSpectrumHT2.title = "HighTower-2";
	analysisSettingsSpectrumHT2.useXerr = false;
	analysisSettingsSpectrumHT2.triggerThreshold = 4.5 + 0.5;
	analysisSettingsSpectrumHT2.useFittedPeakPos = false;
	analysisSettingsSpectrumHT2.showInvBinsFitsSeparately = false;
        analysisSettingsSpectrumHT2.drawOptions    = TDrawOptions(TNamed("HT2", "HighTower-2"), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullCircle, 1.5), "P", "P", TString("HighTower-2"));
        analysisSettingsSpectrumHT2.drawOptionsEta = TDrawOptions(TNamed("HT2", "HighTower-2"), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullTriangleUp, 1.5), "P", "P", TString("HighTower-2 #eta"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumSimMB = analysisSettingsSpectrumMB;
	if (extendMCpT) analysisSettingsSpectrumSimMB.highPt = analysisSettingsSpectrumHT2.highPt;
	analysisSettingsSpectrumSimMB.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.286683e-01) + (8.539953e-03 * x) + (-1.100295e-03 * x * x) + (4.733651e-05 * x * x *x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimMB.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(3.061734e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimMB.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumSimMB.fixedPrescale = analysisSettingsSpectrumMB.fixedPrescale / (107.8*1.163);
	analysisSettingsSpectrumSimMB.correctCpv = false;
	analysisSettingsSpectrumSimMB.correctCpvEta = false;
	analysisSettingsSpectrumSimMB.jetBgFraction.useFit = false;
	analysisSettingsSpectrumSimMB.bgToSigbg.useFit = false;
	analysisSettingsSpectrumSimMB.smearOutliersSigma = 0.05;
	analysisSettingsSpectrumSimMB.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsSpectrumSimMB.lowNormMassRange = false;
	analysisSettingsSpectrumSimMB.lowNormMassRangeLeft = 0.015;
	analysisSettingsSpectrumSimMB.lowNormMassRangeRight = 0.04;
	analysisSettingsSpectrumSimMB.lowNormPoints = false;
	analysisSettingsSpectrumSimMB.lowmassbgPointsLowPt = 1.0;
	analysisSettingsSpectrumSimMB.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumSimMB.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsSpectrumSimMB.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsSpectrumSimMB.name = "SimMB";
	analysisSettingsSpectrumSimMB.title = "MinBias MC #pi^{0}";
        analysisSettingsSpectrumSimMB.drawOptions    = TDrawOptions(TNamed("SimMB", "MinBias MC #pi^{0}"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.5), "P", "P", TString("MinBias #pi^{0} MC"));
        analysisSettingsSpectrumSimMB.drawOptionsEta = TDrawOptions(TNamed("SimMB", "MinBias MC #pi^{0}"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenTriangleUp, 1.5), "P", "P", TString("MinBias #eta (from #pi^{0}) MC"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumSimHT1 = analysisSettingsSpectrumHT1;
	if (extendMCpT) analysisSettingsSpectrumSimHT1.highPt = analysisSettingsSpectrumHT2.highPt;
	analysisSettingsSpectrumSimHT1.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.286683e-01) + (8.539953e-03 * x) + (-1.100295e-03 * x * x) + (4.733651e-05 * x * x *x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimHT1.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(3.061734e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimHT1.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumSimHT1.correctCpv = false;
	analysisSettingsSpectrumSimHT1.correctCpvEta = false;
	analysisSettingsSpectrumSimHT1.jetBgFraction.useFit = false;
	analysisSettingsSpectrumSimHT1.bgToSigbg.useFit = false;
	analysisSettingsSpectrumSimHT1.fixedPrescale = analysisSettingsSpectrumHT1.fixedPrescale / (0.6116*1.039);
	analysisSettingsSpectrumSimHT1.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumSimHT1.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumSimHT1.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumSimHT1.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumSimHT1.smearOutliersSigma = 0.05;
	analysisSettingsSpectrumSimHT1.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsSpectrumSimHT1.lowNormMassRange = false;
	analysisSettingsSpectrumSimHT1.lowNormMassRangeLeft = 0.05;
	analysisSettingsSpectrumSimHT1.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumSimHT1.lowNormPoints = false;
	analysisSettingsSpectrumSimHT1.lowmassbgPointsLowPt = 4.0;
	analysisSettingsSpectrumSimHT1.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumSimHT1.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsSpectrumSimHT1.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsSpectrumSimHT1.name = "SimHT1";
	analysisSettingsSpectrumSimHT1.title = "HighTower-1 MC #pi^{0}";
        analysisSettingsSpectrumSimHT1.drawOptions    = TDrawOptions(TNamed("SimHT1", "HighTower-1 MC #pi^{0}"), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenCircle, 1.5), "P", "P", TString("HighTower-1 #pi^{0} MC"));
        analysisSettingsSpectrumSimHT1.drawOptionsEta = TDrawOptions(TNamed("SimHT1", "HighTower-1 MC #pi^{0}"), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenTriangleUp, 1.5), "P", "P", TString("HighTower-1 #eta (from #pi^{0}) MC"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumSimHT2 = analysisSettingsSpectrumHT2;
	if (extendMCpT) analysisSettingsSpectrumSimHT2.highPt = analysisSettingsSpectrumHT2.highPt;
	analysisSettingsSpectrumSimHT2.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.286683e-01) + (8.539953e-03 * x) + (-1.100295e-03 * x * x) + (4.733651e-05 * x * x *x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimHT2.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(3.061734e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimHT2.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumSimHT2.correctCpv = false;
	analysisSettingsSpectrumSimHT2.correctCpvEta = false;
	analysisSettingsSpectrumSimHT2.jetBgFraction.useFit = false;
	analysisSettingsSpectrumSimHT2.bgToSigbg.useFit = false;
	analysisSettingsSpectrumSimHT2.fixedPrescale = analysisSettingsSpectrumHT2.fixedPrescale / (0.5404*1.153);
	analysisSettingsSpectrumSimHT2.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumSimHT2.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumSimHT2.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumSimHT2.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumSimHT2.smearOutliersSigma = 0.05;
	analysisSettingsSpectrumSimHT2.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsSpectrumSimHT2.lowNormMassRange = false;
	analysisSettingsSpectrumSimHT2.lowNormMassRangeLeft = 0.05;
	analysisSettingsSpectrumSimHT2.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumSimHT2.lowNormPoints = false;
	analysisSettingsSpectrumSimHT2.lowmassbgPointsLowPt = 6.0;
	analysisSettingsSpectrumSimHT2.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumSimHT2.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsSpectrumSimHT2.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsSpectrumSimHT2.name = "SimHT2";
	analysisSettingsSpectrumSimHT2.title = "HighTower-2 MC #pi^{0}";
        analysisSettingsSpectrumSimHT2.drawOptions    = TDrawOptions(TNamed("SimHT2", "HighTower-2 MC #pi^{0}"), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kOpenCircle, 1.5), "P", "P", TString("HighTower-2 #pi^{0} MC"));
        analysisSettingsSpectrumSimHT2.drawOptionsEta = TDrawOptions(TNamed("SimHT2", "HighTower-2 MC #pi^{0}"), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kOpenTriangleUp, 1.5), "P", "P", TString("HighTower-2 #eta (from #pi^{0}) MC"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumSimEtaMB = analysisSettingsSpectrumMB;
	if (extendMCpT) analysisSettingsSpectrumSimEtaMB.highPt = analysisSettingsSpectrumHT2.highPt;
	analysisSettingsSpectrumSimEtaMB.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.5, trueEtaMass * 1.2, false, 0, 100, 0, 100, "((5.809154e-01) + (1.123266e-03 * x)) + ([0]*0)", "RQN");
	analysisSettingsSpectrumSimEtaMB.peakWidthEta = TFitSettings(0.10, 0, 0.10 * 0.3, 0.10 * 1.2, false, 0, 100, 0, 100, "(9.201360e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimEtaMB.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumSimEtaMB.fixedPrescale = analysisSettingsSpectrumMB.fixedPrescale * 0.5 / (34.36*0.5405/0.5);
	analysisSettingsSpectrumSimEtaMB.correctCpv = false;
	analysisSettingsSpectrumSimEtaMB.correctCpvEta = false;
	analysisSettingsSpectrumSimEtaMB.jetBgFraction.useFit = false;
	analysisSettingsSpectrumSimEtaMB.bgToSigbg.useFit = false;
	analysisSettingsSpectrumSimEtaMB.smearOutliersSigma = 0.05;
	analysisSettingsSpectrumSimEtaMB.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsSpectrumSimEtaMB.lowNormMassRange = false;
	analysisSettingsSpectrumSimEtaMB.lowNormMassRangeLeft = 0.015;
	analysisSettingsSpectrumSimEtaMB.lowNormMassRangeRight = 0.04;
	analysisSettingsSpectrumSimEtaMB.lowNormPoints = false;
	analysisSettingsSpectrumSimEtaMB.lowmassbgPointsLowPt = 1.0;
	analysisSettingsSpectrumSimEtaMB.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumSimEtaMB.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsSpectrumSimEtaMB.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsSpectrumSimEtaMB.name = "SimEtaMB";
	analysisSettingsSpectrumSimEtaMB.title = "MinBias MC #eta";
        analysisSettingsSpectrumSimEtaMB.drawOptions    = TDrawOptions(TNamed("SimEtaMB", "MinBias MC #eta"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.5), "P", "P", TString("MinBias #pi^{0} (from #eta) MC"));
        analysisSettingsSpectrumSimEtaMB.drawOptionsEta = TDrawOptions(TNamed("SimEtaMB", "MinBias MC #eta"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenTriangleUp, 1.5), "P", "P", TString("MinBias #eta MC"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumSimEtaHT1 = analysisSettingsSpectrumHT1;
	if (extendMCpT) analysisSettingsSpectrumSimEtaHT1.highPt = analysisSettingsSpectrumHT2.highPt;
	analysisSettingsSpectrumSimEtaHT1.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.7, trueEtaMass * 1.2, false, 0, 100, 0, 100, "((5.759024e-01) + (1.209403e-03 * x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimEtaHT1.peakWidthEta = TFitSettings(0.10, 0, 0.10 * 0.2, 0.10 * 1.2, false, 0, 100, 0, 100, "(9.201360e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimEtaHT1.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumSimEtaHT1.correctCpv = false;
	analysisSettingsSpectrumSimEtaHT1.correctCpvEta = false;
	analysisSettingsSpectrumSimEtaHT1.jetBgFraction.useFit = false;
	analysisSettingsSpectrumSimEtaHT1.bgToSigbg.useFit = false;
	analysisSettingsSpectrumSimEtaHT1.fixedPrescale = analysisSettingsSpectrumHT1.fixedPrescale * 0.5 / (0.1605*0.5016/0.5);
	analysisSettingsSpectrumSimEtaHT1.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumSimEtaHT1.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumSimEtaHT1.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumSimEtaHT1.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumSimEtaHT1.smearOutliersSigma = 0.05;
	analysisSettingsSpectrumSimEtaHT1.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsSpectrumSimEtaHT1.lowNormMassRange = false;
	analysisSettingsSpectrumSimEtaHT1.lowNormMassRangeLeft = 0.05;
	analysisSettingsSpectrumSimEtaHT1.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumSimEtaHT1.lowNormPoints = false;
	analysisSettingsSpectrumSimEtaHT1.lowmassbgPointsLowPt = 4.0;
	analysisSettingsSpectrumSimEtaHT1.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumSimEtaHT1.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsSpectrumSimEtaHT1.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsSpectrumSimEtaHT1.name = "SimEtaHT1";
	analysisSettingsSpectrumSimEtaHT1.title = "HighTower-1 MC #eta";
        analysisSettingsSpectrumSimEtaHT1.drawOptions    = TDrawOptions(TNamed("SimEtaHT1", "HighTower-1 MC #eta"), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenCircle, 1.5), "P", "P", TString("HighTower-1 #pi^{0} (from #eta) MC"));
        analysisSettingsSpectrumSimEtaHT1.drawOptionsEta = TDrawOptions(TNamed("SimEtaHT1", "HighTower-1 MC #eta"), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenTriangleUp, 1.5), "P", "P", TString("HighTower-1 #eta MC"));

	TAnalysisSettingsTrigger analysisSettingsSpectrumSimEtaHT2 = analysisSettingsSpectrumHT2;
	if (extendMCpT) analysisSettingsSpectrumSimEtaHT2.highPt = analysisSettingsSpectrumHT2.highPt;
	analysisSettingsSpectrumSimEtaHT2.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.7, trueEtaMass * 1.2, false, 0, 100, 0, 100, "((5.154388e-01) + (9.879051e-03 * x) + (-3.015698e-04 * x * x)) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimEtaHT2.peakWidthEta = TFitSettings(0.10, 0, 0.10 * 0.2, 0.10 * 1.2, false, 0, 100, 0, 100, "(9.201360e-02) + (0*[0])", "RQN");
	analysisSettingsSpectrumSimEtaHT2.fitFunctionStr = fitFuncZero;
	analysisSettingsSpectrumSimEtaHT2.correctCpv = false;
	analysisSettingsSpectrumSimEtaHT2.correctCpvEta = false;
	analysisSettingsSpectrumSimEtaHT2.jetBgFraction.useFit = false;
	analysisSettingsSpectrumSimEtaHT2.bgToSigbg.useFit = false;
	analysisSettingsSpectrumSimEtaHT2.fixedPrescale = analysisSettingsSpectrumHT2.fixedPrescale * 0.5 / (0.07455*0.4642/0.5);
	analysisSettingsSpectrumSimEtaHT2.calculatePrescaleFromDB = false;
	analysisSettingsSpectrumSimEtaHT2.calculatePrescaleFromSim = false;
	analysisSettingsSpectrumSimEtaHT2.calculatePrescaleFromPoints = false;
	analysisSettingsSpectrumSimEtaHT2.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsSpectrumSimEtaHT2.smearOutliersSigma = 0.05;
	analysisSettingsSpectrumSimEtaHT2.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsSpectrumSimEtaHT2.lowNormMassRange = false;
	analysisSettingsSpectrumSimEtaHT2.lowNormMassRangeLeft = 0.05;
	analysisSettingsSpectrumSimEtaHT2.lowNormMassRangeRight = 0.07;
	analysisSettingsSpectrumSimEtaHT2.lowNormPoints = false;
	analysisSettingsSpectrumSimEtaHT2.lowmassbgPointsLowPt = 6.0;
	analysisSettingsSpectrumSimEtaHT2.lowmassbgPointsHighPt = 100.0;
	analysisSettingsSpectrumSimEtaHT2.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsSpectrumSimEtaHT2.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsSpectrumSimEtaHT2.name = "SimEtaHT2";
	analysisSettingsSpectrumSimEtaHT2.title = "HighTower-2 MC #eta";
        analysisSettingsSpectrumSimEtaHT2.drawOptions    = TDrawOptions(TNamed("SimEtaHT2", "HighTower-2 MC #eta"), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kOpenCircle, 1.5), "P", "P", TString("HighTower-2 #pi^{0} (from #eta) MC"));
        analysisSettingsSpectrumSimEtaHT2.drawOptionsEta = TDrawOptions(TNamed("SimEtaHT2", "HighTower-2 MC #eta"), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kOpenTriangleUp, 1.5), "P", "P", TString("HighTower-2 #eta MC"));

	TAnalysisSettings analysisSettingsSpectrum;
	analysisSettingsSpectrum.settingsMB = analysisSettingsSpectrumMB;
	analysisSettingsSpectrum.settingsHT1 = analysisSettingsSpectrumHT1;
	analysisSettingsSpectrum.settingsHT2 = analysisSettingsSpectrumHT2;
	analysisSettingsSpectrum.settingsSimMB = analysisSettingsSpectrumSimMB;
	analysisSettingsSpectrum.settingsSimHT1 = analysisSettingsSpectrumSimHT1;
	analysisSettingsSpectrum.settingsSimHT2 = analysisSettingsSpectrumSimHT2;
	analysisSettingsSpectrum.settingsSimEtaMB = analysisSettingsSpectrumSimEtaMB;
	analysisSettingsSpectrum.settingsSimEtaHT1 = analysisSettingsSpectrumSimEtaHT1;
	analysisSettingsSpectrum.settingsSimEtaHT2 = analysisSettingsSpectrumSimEtaHT2;
	analysisSettingsSpectrum.smoothEffMB = false;
	analysisSettingsSpectrum.smoothEffHT1 = false;
	analysisSettingsSpectrum.smoothEffHT2 = false;
	analysisSettingsSpectrum.effFitMB = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "[0] * (1 + TMath::Erf(([1] * x) + [2]))", "RQN");
	analysisSettingsSpectrum.effFitHT1 = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "[0] * (1 + TMath::Erf(([1] * x) + [2]))", "RQN");
	analysisSettingsSpectrum.effFitHT2 = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "[0] * (1 + TMath::Erf(([1] * x) + [2]))", "RQN");
	analysisSettingsSpectrum.dontCorrectEff = false;
	analysisSettingsSpectrum.correctEffMatrix = false;
	analysisSettingsSpectrum.smoothEffMBEta = false;
	analysisSettingsSpectrum.smoothEffHT1Eta = false;
	analysisSettingsSpectrum.smoothEffHT2Eta = false;
	analysisSettingsSpectrum.effFitMBEta = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "[0] * (1 + TMath::Erf(([1] * x) + [2]))", "RQN");
	analysisSettingsSpectrum.effFitHT1Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "[0] * (1 + TMath::Erf(([1] * x) + [2]))", "RQN");
	analysisSettingsSpectrum.effFitHT2Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "[0] * (1 + TMath::Erf(([1] * x) + [2]))", "RQN");
	analysisSettingsSpectrum.dontCorrectEffEta = false;
	analysisSettingsSpectrum.correctEffMatrixEta = false;
	analysisSettingsSpectrum.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, true, -100, +100, -100, +100, "pol1", "RQN");
	analysisSettingsSpectrum.cpvCorrectionSim = TFitSettings(1.0, 0, 0.0, 1.0, false, -100, +100, -100, +100, "[0] + (0*x*[1])", "RQN");
	analysisSettingsSpectrum.cpvCorrectionSimEta = TFitSettings(1.0, 0, 0.0, 1.0, false, -100, +100, -100, +100, "[0] + (0*x*[1])", "RQN");
	analysisSettingsSpectrum.shiftFuncStr = "0*[0] + ([1]*x) + (0*[2]*x*x)";
	analysisSettingsSpectrum.showReal = true;
	analysisSettingsSpectrum.showSim = true;
	analysisSettingsSpectrum.print = true;
	analysisSettingsSpectrum.multiplicityDistributions = 0;
	analysisSettingsSpectrum.pointMultiplicityDistributions = true;
	analysisSettingsSpectrum.showTriggersSeparately = true;
	analysisSettingsSpectrum.showWeightUsed = false;
	analysisSettingsSpectrum.correctForPreciseBinPt = true;
	analysisSettingsSpectrum.correctForPreciseBinPtFromSim = false;
	analysisSettingsSpectrum.correctForPreciseBinPtHorizontal = false;
	analysisSettingsSpectrum.correctForPreciseBinPtNFitPoints = 2;
	analysisSettingsSpectrum.correctForPreciseBinPtFitOption = "RQN LL";
	analysisSettingsSpectrum.correctForPreciseBinPtFitOptionWeight = "POWERLAW";
	analysisSettingsSpectrum.correctForJacobian = true;
	analysisSettingsSpectrum.showWeightNextIteration = false;
	analysisSettingsSpectrum.saveWeightNextIteration = false;
	analysisSettingsSpectrum.tpcVertexEffFunc = "[0] + exp(-[1]*x)";
	analysisSettingsSpectrum.showTPCVertexEff = false;
	analysisSettingsSpectrum.saveTPCVertexEff = false;
	analysisSettingsSpectrum.bbcVertexEffFunc = "[0] + exp(-[1]*x)";
	analysisSettingsSpectrum.showBBCVertexEff = false;
	analysisSettingsSpectrum.saveBBCVertexEff = false;
	analysisSettingsSpectrum.weightFile = "./weights_dAu2003_new.root";
	analysisSettingsSpectrum.correctMeanAcceptanceRatioMB = true;
	analysisSettingsSpectrum.correctMeanAcceptanceRatioHT1 = true;
	analysisSettingsSpectrum.correctMeanAcceptanceRatioHT2 = true;
	analysisSettingsSpectrum.correctMeanAcceptanceRatioEtaMB = true;
	analysisSettingsSpectrum.correctMeanAcceptanceRatioEtaHT1 = true;
	analysisSettingsSpectrum.correctMeanAcceptanceRatioEtaHT2 = true;
	analysisSettingsSpectrum.name = TString("dAu") + (addDirName ? TString("_") + DATA_DIR : TString());
	analysisSettingsSpectrum.title = TString("#font[12]{d}+Au");
        analysisSettingsSpectrum.drawOptions    = TDrawOptions(TNamed("dAu", "#font[12]{d}+Au"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullCircle, 1.5), "P", "P", TString("#font[12]{d}+Au #pi^{0}"));
        analysisSettingsSpectrum.drawOptionsEta = TDrawOptions(TNamed("dAuEta", "#font[12]{d}+Au #eta"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleUp, 1.5), "P", "P", TString("#font[12]{d}+Au #eta"));
	analysisSettingsSpectrum.savePtSlicesTrig = false;
	analysisSettingsSpectrum.PtSlicesTrigFilename = "./pTslices.html";
	analysisSettingsSpectrum.savePtSlicesTrigSim = false;
	analysisSettingsSpectrum.PtSlicesTrigSimFilename = "./pTslices_sim.html";
	analysisSettingsSpectrum.savePtSlicesTrigSimEta = false;
	analysisSettingsSpectrum.PtSlicesTrigSimEtaFilename = "./pTslices_sim_eta.html";
	analysisSettingsSpectrum.savePtSlicesTrigSimEtabg = false;
	analysisSettingsSpectrum.PtSlicesTrigSimEtabgFilename = "./pTslices_sim_etabg.html";
	analysisSettingsSpectrum.triggerThresholdAdd = 4.0;
	analysisSettingsSpectrum.prescalesDBFilename = "./dAu2003Prescales.txt";
	analysisSettingsSpectrum.saveDataArrays = true;
	analysisSettingsSpectrum.dataArraysFilenameFormat = "./data_arrays_%s.C";
	analysisSettingsSpectrum.saveRunTimes = false;
	analysisSettingsSpectrum.runTimesFilenameFormat = "./dAu2003Times_%s.txt";
	analysisSettingsSpectrum.saveRunTimesSim = false;
	analysisSettingsSpectrum.runTimesSimFilenameFormat = "./dAu2003Times_pi0_%s.txt";
	analysisSettingsSpectrum.saveRunTimesSimEta = false;
	analysisSettingsSpectrum.runTimesSimEtaFilenameFormat = "./dAu2003Times_eta_%s.txt";
	analysisSettingsSpectrum.savePrescales = false;
	analysisSettingsSpectrum.prescalesSaveFilename = "./dAu2003Prescales.txt";
	analysisSettingsSpectrum.saveBadSimPi0EventsMB = false;
	analysisSettingsSpectrum.badSimPi0EventsMBFilename = "./badEvents_pythiaPi0Pt_dAu2003_MB.txt";
	analysisSettingsSpectrum.saveBadSimPi0EventsHT1 = false;
	analysisSettingsSpectrum.badSimPi0EventsHT1Filename = "./badEvents_pythiaPi0Pt_dAu2003_HT1.txt";
	analysisSettingsSpectrum.saveBadSimPi0EventsHT2 = false;
	analysisSettingsSpectrum.badSimPi0EventsHT2Filename = "./badEvents_pythiaPi0Pt_dAu2003_HT2.txt";
	analysisSettingsSpectrum.dbConnectString = "mysql://dbx.star.bnl.gov:3316/RunLog_onl";
	analysisSettingsSpectrum.NBinaryCollisions = -1;
	analysisSettingsSpectrum.NBinaryCollisionsErr = -1;
	analysisSettingsSpectrum.materialCorrectionSim = 1.15;
	analysisSettingsSpectrum.materialCorrectionSimGamma = 1.08;
	analysisSettingsSpectrum.crossectionMB = 2.21e3 * 0.95; // Total d+Au cross section @ 200 GeV/c is 2.21 +/- 0.09 b
	analysisSettingsSpectrum.crossectionMBErr = 0.09e3 * 0.95;
	analysisSettingsSpectrum.crossectionMBTotalFraction = 0.95; // d+Au 2003 MinBias trigger efficiency is 0.95 +/- 0.03
	analysisSettingsSpectrum.crossectionMBTotalFractionErr = 0.03;
	analysisSettingsSpectrum.crossectionSimMB = (2.21e3 * 0.95);
	analysisSettingsSpectrum.crossectionSimMBErr = 0;
	analysisSettingsSpectrum.crossectionSimMBTotalFraction = 0.95;
	analysisSettingsSpectrum.crossectionSimMBTotalFractionErr = 0;
	analysisSettingsSpectrum.crossectionSimHT1 = ((2.21e3 * 0.95));
	analysisSettingsSpectrum.crossectionSimHT2 = ((2.21e3 * 0.95));
	analysisSettingsSpectrum.crossectionSimEtaMB = (2.21e3 * 0.95) / 0.5;
	analysisSettingsSpectrum.crossectionSimEtaMBErr = 0;
	analysisSettingsSpectrum.crossectionSimEtaMBTotalFraction = 0.95;
	analysisSettingsSpectrum.crossectionSimEtaMBTotalFractionErr = 0;
	analysisSettingsSpectrum.crossectionSimEtaHT1 = ((2.21e3 * 0.95)) / 0.5;
	analysisSettingsSpectrum.crossectionSimEtaHT2 = ((2.21e3 * 0.95)) / 0.5;
	analysisSettingsSpectrum.crossectionInel = 2.21e3; // Total d+Au cross section @ 200 GeV/c is 2.21 +/- 0.09 b
	analysisSettingsSpectrum.crossectionInelErr = 0;
	analysisSettingsSpectrum.energyScaleErr = 0.05;
	analysisSettingsSpectrum.energyScaleErrFitOption = "RQN LL";
	analysisSettingsSpectrum.energyScaleErrFitOptionWeight = "POWERLAW POLYNOMIAL POL5";
	analysisSettingsSpectrum.energyScaleErrFitOptionWeightRaw = "POWERLAW POLYNOMIAL POL3";
	analysisSettingsSpectrum.usePPpQCD = false;

	TAnalysisSettings analysisSettingsRcp = analysisSettingsSpectrum;
	analysisSettingsRcp.dontCorrectEff = true;
	analysisSettingsRcp.dontCorrectEffEta = true;
	analysisSettingsRcp.saveWeightNextIteration = false;
	analysisSettingsRcp.saveTPCVertexEff = false;
	analysisSettingsRcp.showTPCVertexEff = false;
	analysisSettingsRcp.saveBBCVertexEff = false;
	analysisSettingsRcp.showBBCVertexEff = false;
	analysisSettingsRcp.saveRunTimes = false;
	analysisSettingsRcp.saveRunTimesSim = false;
	analysisSettingsRcp.saveRunTimesSimEta = false;
	analysisSettingsRcp.savePrescales = false;

	TAnalysisSettingsTrigger analysisSettingsPPMB = analysisSettingsSpectrumMB;
	analysisSettingsPPMB.lowPt = 1.0;
	analysisSettingsPPMB.highPt = 5.0;
	analysisSettingsPPMB.lowPtUse = pp_MB;
	analysisSettingsPPMB.highPtUse = pp_MB_HT1;
	analysisSettingsPPMB.lowPtBinSize = 0.5;
	analysisSettingsPPMB.midPtBinSize = 1.0;
	analysisSettingsPPMB.highPtBinSize = 2.0;
	analysisSettingsPPMB.switchBinSizePt1 = 3.0;
	analysisSettingsPPMB.switchBinSizePt2 = 5.0;
	analysisSettingsPPMB.rebinNumBinsWeight = 1.0;
	analysisSettingsPPMB.rebinMeanErrorWeight = 1.0;
	analysisSettingsPPMB.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 12, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPMB.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.7, 0.025 * 1.5, false, 0, 100, 0, 100, "(2.638381e-02) + (0*[0])", "RQN");
	analysisSettingsPPMB.peakLeftEta = TFitSettings(0.3, 0.02);
	analysisSettingsPPMB.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPMB.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.8, 0.06 * 1.5, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	analysisSettingsPPMB.fitParameter6 = TFitSettings(0, 0, 0, 10000);
	analysisSettingsPPMB.fitParameter7 = TFitSettings();
	analysisSettingsPPMB.fitParameters = TFitSettings();
	analysisSettingsPPMB.fitFunctionStr = fitFuncZero;
	analysisSettingsPPMB.fitOption = "RQN";
	analysisSettingsPPMB.fixedPrescale = 1.0;
	analysisSettingsPPMB.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPMB.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPMB.getJetBgFractionFromCandidatesPtFit = false;
	analysisSettingsPPMB.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 10, "([0]*0) + (3.288739e-02 - 0.15) + (6.936733e-02 * 1.4 * x)", "RQN");
	analysisSettingsPPMB.bgToSigbg = TFitSettings(0.5);
	analysisSettingsPPMB.smearOutliersSigma = 0;
	analysisSettingsPPMB.smoothFitOption = "QN";
	analysisSettingsPPMB.smoothBg = 0;//0.04;
	analysisSettingsPPMB.smoothBgRandom = 0;//0.04;
	analysisSettingsPPMB.smoothBgBack = 0;
	analysisSettingsPPMB.mixNormFixed = 0.0;
	analysisSettingsPPMB.mixNormMassRange = false;
	analysisSettingsPPMB.mixNormLeft = 1;//0.25;
	analysisSettingsPPMB.mixNormRight = 4;//0.4;
	analysisSettingsPPMB.mixNormEntries = true;
	analysisSettingsPPMB.mixNormFixedRandom = 0.0;
	analysisSettingsPPMB.mixNormMassRangeRandom = false;
	analysisSettingsPPMB.mixNormLeftRandom = 0.9;
	analysisSettingsPPMB.mixNormRightRandom = 1.2;
	analysisSettingsPPMB.mixNormEntriesRandom = true;
	analysisSettingsPPMB.mixNormFixedBack = 0.0;
	analysisSettingsPPMB.mixNormMassRangeBack = true;
	analysisSettingsPPMB.mixNormLeftBack = 2.0;
	analysisSettingsPPMB.mixNormRightBack = 5.0;
	analysisSettingsPPMB.mixNormEntriesBack = false;
	analysisSettingsPPMB.subtractBg = false;
	analysisSettingsPPMB.subtractBgRandom = false;
	analysisSettingsPPMB.subtractBgBack = false;
	analysisSettingsPPMB.subtractBgLowMass = true;
	analysisSettingsPPMB.subtractBgLowMassPeakShape = true;
	analysisSettingsPPMB.subtractBgLowMassPeakShapeEta = true;
	analysisSettingsPPMB.subtractBgFitDistribution = true;
	analysisSettingsPPMB.useEventMixingBgRandom = true;
	analysisSettingsPPMB.useEventMixingNotmatchedBgRandom = true;
	analysisSettingsPPMB.useEventJetMixingBg = true;
	analysisSettingsPPMB.useEventJetNotmatchedMixingBg = false;
	analysisSettingsPPMB.useEventJetBackMixingBg = false;
	analysisSettingsPPMB.useEventShufflingBg = false;
	analysisSettingsPPMB.fitBgDistribution = true;
	analysisSettingsPPMB.fitBgRandomDistribution = true;
	analysisSettingsPPMB.fitBgLowMassDistribution = false;
	analysisSettingsPPMB.fitPeakShapeDistribution = false && !isPythia;
	analysisSettingsPPMB.fitPeakShapeEtaDistribution = false && !isPythia;
	analysisSettingsPPMB.fitPeakShapeEtabgDistribution = false && !isPythia;
	analysisSettingsPPMB.fitDistributionLeft = 0.3;
	analysisSettingsPPMB.fitDistributionRight = 0.9;//0.5; // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsPPMB.fitDistributionLeft2 = 0.9;//0.85; // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsPPMB.fitDistributionRight2 = 1.2;//0.95; // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsPPMB.fitDistributionOption = "RQN UNORM";
	//analysisSettingsPPMB.lowNormFixed = TFitSettings(2.0, -1.0/4);
	//analysisSettingsPPMB.lowNormFixed = TFitSettings(0.75, -0.45/5);
	analysisSettingsPPMB.lowNormFixed = TFitSettings(1.12, -0.19);
	//analysisSettingsPPMB.lowNormFixed = TFitSettings(1.0);
	analysisSettingsPPMB.lowNormMassRange = false;
	analysisSettingsPPMB.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPMB.lowNormMassRangeRight = 0.095;
	analysisSettingsPPMB.lowNormPoints = false;
	analysisSettingsPPMB.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPMB.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPMB.lowNormPointsPt = true;
	analysisSettingsPPMB.lowmassbgNbar = TFitSettings(1.0, 0, 0, 1, false, 0, 100, 0, 100, "(1.0/((1.0/( (1+TMath::Exp(-x/1.5)) * (0.2689*x*TMath::Exp(-0.3258*x) / ((7.334427e-01)+(-9.284123e-02*x)+(8.582639e-03*x*x)+(-2.596326e-04*x*x*x))) * (1-0.2+(0.2*0.358/0.639)) / (1.0+(3.01325e-03*x*x)) ))+1.0)) + (0*[0])", "RQN");
	analysisSettingsPPMB.correctCpv = true;
	analysisSettingsPPMB.correctCpvEta = true;
	analysisSettingsPPMB.rebinNumBinsWeightCpv = 1.0;
	analysisSettingsPPMB.rebinMeanErrorWeightCpv = 1.0;
	analysisSettingsPPMB.fitLeftCpv = TFitSettings(0);
	analysisSettingsPPMB.fitRightCpv = TFitSettings(0.7);
	analysisSettingsPPMB.peakLeftCpv = TFitSettings(0.0);
	analysisSettingsPPMB.peakRightCpv = TFitSettings(-0.08);
	analysisSettingsPPMB.peakPositionCpv = TFitSettings(0.0, 0.0, 0.0, 0.0);
	analysisSettingsPPMB.peakWidthCpv = TFitSettings(0.03, 0.0, 0.006, 0.07);
	analysisSettingsPPMB.fitParameter6Cpv = TFitSettings(1000, 0, 0, 1000000);
	analysisSettingsPPMB.fitParameter7Cpv = TFitSettings(9.0, 0.0, 3.0, 15.0);
	analysisSettingsPPMB.fitParametersCpv = TFitSettings();
	analysisSettingsPPMB.fitFunctionStrCpv = fitFuncTrackDistPeak;
	analysisSettingsPPMB.fitOptionCpv = "RQN";
	analysisSettingsPPMB.binRangeLeftCpv = 0.0;
	analysisSettingsPPMB.binRangeRightCpv = 0.7;
	analysisSettingsPPMB.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 1, 4, -100, +100, "pol0", "RQN");
	analysisSettingsPPMB.correctVertexFindingEff = false;
	analysisSettingsPPMB.vertexFindingEff = vertEffPPMinBias;
	analysisSettingsPPMB.correctVertexFindingEffMB = false;
	analysisSettingsPPMB.showInvBinsFitsSeparately = false;

	TAnalysisSettingsTrigger analysisSettingsPPHT1 = analysisSettingsSpectrumHT1;
	analysisSettingsPPHT1.lowPt = 4.0;
	analysisSettingsPPHT1.highPt = 10.0;
	analysisSettingsPPHT1.lowPtUse = pp_MB_HT1;
	analysisSettingsPPHT1.highPtUse = pp_HT1_HT2;
	analysisSettingsPPHT1.lowPtBinSize = 0.5;
	analysisSettingsPPHT1.midPtBinSize = 1.0;
	analysisSettingsPPHT1.highPtBinSize = 2.0;
	analysisSettingsPPHT1.switchBinSizePt1 = 4.0;
	analysisSettingsPPHT1.switchBinSizePt2 = 8.0;
	analysisSettingsPPHT1.rebinNumBinsWeight = 1.0;
	analysisSettingsPPHT1.rebinMeanErrorWeight = 1.0;
	analysisSettingsPPHT1.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPHT1.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(2.704745e-02) + (0*[0])", "RQN");
	analysisSettingsPPHT1.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPHT1.peakWidthEta = TFitSettings(0.06, 0, 0.06/*0.06 * 0.5*/, 0.07/*0.06 * 2.0*/, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	analysisSettingsPPHT1.peakLeftEta = TFitSettings(0.2, 0.02);
	analysisSettingsPPHT1.peakRightEta = TFitSettings(0.8, 0.02);
	analysisSettingsPPHT1.fitParameter6 = TFitSettings(2, 0, 0, 10);
	analysisSettingsPPHT1.fitParameter7 = TFitSettings();
	analysisSettingsPPHT1.fitParameters = TFitSettings();
	//analysisSettingsPPHT1.fitFunctionStr = "gaus(0) + gaus(3) + ([6]*[6]*x*x)";
	analysisSettingsPPHT1.fitFunctionStr = fitFuncZero;
	analysisSettingsPPHT1.fitOption = "RQN";
	analysisSettingsPPHT1.fixedPrescale = 0;
	analysisSettingsPPHT1.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPHT1.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPHT1.getJetBgFractionFromCandidatesPtFit = false;
	analysisSettingsPPHT1.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 10, "([0]*0) + (3.288739e-02 - 0.15) + (6.936733e-02 * 1.4 * x)", "RQN");
	//analysisSettingsPPHT1.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 12, 0.0, 100.0, "pol2", "RQN");
	analysisSettingsPPHT1.bgToSigbg = TFitSettings(0.5);
	analysisSettingsPPHT1.fixedPrescale = /*7340.76;//*/4670.96; // towerMB_smd2HT1_gammaetacut / towerMB_smd2HT1
	analysisSettingsPPHT1.calculatePrescaleFromDB = false;
	analysisSettingsPPHT1.calculatePrescaleFromSim = false;
	analysisSettingsPPHT1.calculatePrescaleFromPoints = false;
	analysisSettingsPPHT1.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsPPHT1.prescalePointsLowPt = 4.0;
	analysisSettingsPPHT1.prescalePointsHighPt = 8.0;
	analysisSettingsPPHT1.prescalePointsLowPtIntegral = 4.0;
	analysisSettingsPPHT1.prescalePointsHighPtIntegral = 8.0;
	analysisSettingsPPHT1.correctCpv = true;
	analysisSettingsPPHT1.correctCpvEta = true;
	analysisSettingsPPHT1.rebinNumBinsWeightCpv = 1.0;
	analysisSettingsPPHT1.rebinMeanErrorWeightCpv = 1.0;
	analysisSettingsPPHT1.fitLeftCpv = TFitSettings(0);
	analysisSettingsPPHT1.fitRightCpv = TFitSettings(0.7);
	analysisSettingsPPHT1.peakLeftCpv = TFitSettings(0.0);
	analysisSettingsPPHT1.peakRightCpv = TFitSettings(-0.08);
	analysisSettingsPPHT1.peakPositionCpv = TFitSettings(0.0, 0.0, 0.0, 0.0);
	analysisSettingsPPHT1.peakWidthCpv = TFitSettings(0.03, 0.0, 0.006, 0.07);
	analysisSettingsPPHT1.fitParameter6Cpv = TFitSettings(1000, 0, 0, 1000000);
	analysisSettingsPPHT1.fitParameter7Cpv = TFitSettings(10.0, 0.5, 7.0, 22.0);
	analysisSettingsPPHT1.fitParametersCpv = TFitSettings();
	analysisSettingsPPHT1.fitFunctionStrCpv = fitFuncTrackDistPeak;
	analysisSettingsPPHT1.fitOptionCpv = "RQN";
	analysisSettingsPPHT1.binRangeLeftCpv = 0.0;
	analysisSettingsPPHT1.binRangeRightCpv = 0.7;
	analysisSettingsPPHT1.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 3, 11, -100, +100, "pol0", "RQNW");
	analysisSettingsPPHT1.smearOutliersSigma = 0;
	analysisSettingsPPHT1.smoothFitOption = "QN";
	analysisSettingsPPHT1.smoothBg = 0;//0.04;
	analysisSettingsPPHT1.mixNormFixed = 0.0;
	analysisSettingsPPHT1.mixNormMassRange = false;
	analysisSettingsPPHT1.mixNormLeft = 1.0;
	analysisSettingsPPHT1.mixNormRight = 4.0;
	analysisSettingsPPHT1.mixNormEntries = true;
	analysisSettingsPPHT1.smoothBgRandom = 0;//0.04;
	analysisSettingsPPHT1.mixNormFixedRandom = 0.0;
	analysisSettingsPPHT1.mixNormMassRangeRandom = false;
	analysisSettingsPPHT1.mixNormLeftRandom = 0.9;
	analysisSettingsPPHT1.mixNormRightRandom = 1.2;
	analysisSettingsPPHT1.mixNormEntriesRandom = true;
	analysisSettingsPPHT1.smoothBgBack = 0;
	analysisSettingsPPHT1.mixNormFixedBack = 0.0;
	analysisSettingsPPHT1.mixNormMassRangeBack = true;
	analysisSettingsPPHT1.mixNormLeftBack = 2.0;
	analysisSettingsPPHT1.mixNormRightBack = 5.0;
	analysisSettingsPPHT1.mixNormEntriesBack = false;
	analysisSettingsPPHT1.subtractBg = false;
	analysisSettingsPPHT1.subtractBgRandom = false;
	analysisSettingsPPHT1.subtractBgBack = false;
	analysisSettingsPPHT1.subtractBgLowMass = true;
	analysisSettingsPPHT1.subtractBgLowMassPeakShape = true;
	analysisSettingsPPHT1.subtractBgLowMassPeakShapeEta = true;
	analysisSettingsPPHT1.subtractBgFitDistribution = true;
	analysisSettingsPPHT1.useEventMixingBgRandom = true;
	analysisSettingsPPHT1.useEventMixingNotmatchedBgRandom = false;
	analysisSettingsPPHT1.useEventJetMixingBg = true;
	analysisSettingsPPHT1.useEventJetNotmatchedMixingBg = false;
	analysisSettingsPPHT1.useEventJetBackMixingBg = false;
	analysisSettingsPPHT1.useEventShufflingBg = false;
	analysisSettingsPPHT1.fitBgDistribution = true;
	analysisSettingsPPHT1.fitBgRandomDistribution = true;
	analysisSettingsPPHT1.fitBgLowMassDistribution = false;
	analysisSettingsPPHT1.fitPeakShapeDistribution = false && !isPythia;
	analysisSettingsPPHT1.fitPeakShapeEtaDistribution = false && !isPythia;
	analysisSettingsPPHT1.fitPeakShapeEtabgDistribution = false && !isPythia;
	analysisSettingsPPHT1.fitDistributionLeft = 0.3;//0.28;
	analysisSettingsPPHT1.fitDistributionRight = 0.4;//2;//0.38;
	analysisSettingsPPHT1.fitDistributionLeft2 = 0.9;
	analysisSettingsPPHT1.fitDistributionRight2 = 1.2;
	analysisSettingsPPHT1.fitDistributionOption = "RQN UNORM";
	analysisSettingsPPHT1.lowNormFixed = TFitSettings(1.0);
	analysisSettingsPPHT1.lowNormMassRange = false;
	analysisSettingsPPHT1.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPHT1.lowNormMassRangeRight = 0.06;
	analysisSettingsPPHT1.lowNormPoints = false;
	analysisSettingsPPHT1.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPHT1.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPHT1.lowNormPointsPt = true;
	analysisSettingsPPHT1.lowmassbgNbar = TFitSettings(0.3, 0, 0, 1, true, 0, 100, 0, 100, "(1.0/((1.0/((0.2689*x*TMath::Exp(-0.3258*x) / ((7.334427e-01)+(-9.284123e-02*x)+(8.582639e-03*x*x)+(-2.596326e-04*x*x*x))) * (1-0.2+(0.2*0.358/0.639)) * (10.39*TMath::Landau(x,3.214,0.695)) / (1.0+(3.01325e-03*x*x)) ))+1.0)) + (0*[0])", "RQN");
	analysisSettingsPPHT1.correctVertexFindingEff = false;
	analysisSettingsPPHT1.vertexFindingEff = 1.0;
	analysisSettingsPPHT1.correctVertexFindingEffMB = true;
	analysisSettingsPPHT1.triggerThreshold = 2.6 + 0.0;
	analysisSettingsPPHT1.showInvBinsFitsSeparately = false;

	TAnalysisSettingsTrigger analysisSettingsPPHT2 = analysisSettingsSpectrumHT2;
	analysisSettingsPPHT2.lowPt = 6.0;
	analysisSettingsPPHT2.highPt = 17.0;
	analysisSettingsPPHT2.lowPtUse = pp_HT1_HT2;
	analysisSettingsPPHT2.highPtUse = pp_HT2;
	analysisSettingsPPHT2.lowPtBinSize = 1.0;
	analysisSettingsPPHT2.midPtBinSize = 2.0;
	analysisSettingsPPHT2.highPtBinSize = 3.0;
	analysisSettingsPPHT2.switchBinSizePt1 = 8.0;
	analysisSettingsPPHT2.switchBinSizePt2 = 14.0;
	analysisSettingsPPHT2.rebinNumBinsWeight = 1.0;
	analysisSettingsPPHT2.rebinMeanErrorWeight = 1.0;
	analysisSettingsPPHT2.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.7, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPHT2.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(2.914530e-02) + (0*[0])", "RQN");
	analysisSettingsPPHT2.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, 0.58/*trueEtaMass * 1.3*/, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPHT2.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.6, 0.06 * 1.5, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	analysisSettingsPPHT2.peakLeftEta = TFitSettings(0.1, 0.01);
	analysisSettingsPPHT2.peakRightEta = TFitSettings(0.8, 0.02);
	analysisSettingsPPHT2.fitParameter6 = TFitSettings(2, 0, 0, 10);
	analysisSettingsPPHT2.fitParameter7 = TFitSettings();
	analysisSettingsPPHT2.fitParameters = TFitSettings();
	//analysisSettingsPPHT2.fitFunctionStr = "gaus(0) + gaus(3) + ([6]*[6]*x*x)";
	analysisSettingsPPHT2.fitFunctionStr = fitFuncZero;
	analysisSettingsPPHT2.fitOption = "RQN";
	analysisSettingsPPHT2.fixedPrescale = 0;
	analysisSettingsPPHT2.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPHT2.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPHT2.getJetBgFractionFromCandidatesPtFit = false;
	analysisSettingsPPHT2.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 10, "([0]*0) + (3.288739e-02 - 0.15) + (6.936733e-02 * 1.4 * x)", "RQN");
	//analysisSettingsPPHT2.jetBgFraction = TFitSettings(1.0, 0.0, 0.0, 1.0, true, 0.0, 12, 0.0, 100.0, "pol2", "RQN");
	analysisSettingsPPHT2.bgToSigbg = TFitSettings(0.5);
	analysisSettingsPPHT2.fixedPrescale = /*33897.6;//*/19632.8; // towerMB_smd2HT1_gammaetacut / towerMB_smd2HT1
	analysisSettingsPPHT2.calculatePrescaleFromDB = false;
	analysisSettingsPPHT2.calculatePrescaleFromSim = false;
	analysisSettingsPPHT2.calculatePrescaleFromPoints = false;
	analysisSettingsPPHT2.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsPPHT2.prescalePointsLowPt = 7.0;
	analysisSettingsPPHT2.prescalePointsHighPt = 25.0;
	analysisSettingsPPHT2.prescalePointsLowPtIntegral = 7.0;
	analysisSettingsPPHT2.prescalePointsHighPtIntegral = 25.0;
	analysisSettingsPPHT2.correctCpv = true;
	analysisSettingsPPHT2.correctCpvEta = true;
	analysisSettingsPPHT2.rebinNumBinsWeightCpv = 1.0;
	analysisSettingsPPHT2.rebinMeanErrorWeightCpv = 1.0;
	analysisSettingsPPHT2.fitLeftCpv = TFitSettings(0);
	analysisSettingsPPHT2.fitRightCpv = TFitSettings(0.7);
	analysisSettingsPPHT2.peakLeftCpv = TFitSettings(0.0);
	analysisSettingsPPHT2.peakRightCpv = TFitSettings(-0.08);
	analysisSettingsPPHT2.peakPositionCpv = TFitSettings(0.0, 0.0, 0.0, 0.0);
	analysisSettingsPPHT2.peakWidthCpv = TFitSettings(0.03, 0.0, 0.006, 0.07);
	analysisSettingsPPHT2.fitParameter6Cpv = TFitSettings(1000, 0, 0, 1000000);
	analysisSettingsPPHT2.fitParameter7Cpv = TFitSettings(10.0, 0.5, 7.0, 22.0);
	analysisSettingsPPHT2.fitParametersCpv = TFitSettings();
	analysisSettingsPPHT2.fitFunctionStrCpv = fitFuncTrackDistPeak;
	analysisSettingsPPHT2.fitOptionCpv = "RQN";
	analysisSettingsPPHT2.binRangeLeftCpv = 0.0;
	analysisSettingsPPHT2.binRangeRightCpv = 0.7;
	analysisSettingsPPHT2.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 0, 12, -100, +100, "pol0", "RQNW");
	analysisSettingsPPHT2.smearOutliersSigma = 0;
	analysisSettingsPPHT2.smoothFitOption = "QN";
	analysisSettingsPPHT2.smoothBg = 0;//0.04;
	analysisSettingsPPHT2.mixNormFixed = 0.0;
	analysisSettingsPPHT2.mixNormMassRange = false;
	analysisSettingsPPHT2.mixNormLeft = 1.0;
	analysisSettingsPPHT2.mixNormRight = 4.0;
	analysisSettingsPPHT2.mixNormEntries = true;
	analysisSettingsPPHT2.smoothBgRandom = 0;//0.04;
	analysisSettingsPPHT2.mixNormFixedRandom = 0.0;
	analysisSettingsPPHT2.mixNormMassRangeRandom = false;
	analysisSettingsPPHT2.mixNormLeftRandom = 0.9;
	analysisSettingsPPHT2.mixNormRightRandom = 1.2;
	analysisSettingsPPHT2.mixNormEntriesRandom = true;
	analysisSettingsPPHT2.smoothBgBack = 0;
	analysisSettingsPPHT2.mixNormFixedBack = 0.0;
	analysisSettingsPPHT2.mixNormMassRangeBack = true;
	analysisSettingsPPHT2.mixNormLeftBack = 2.0;
	analysisSettingsPPHT2.mixNormRightBack = 5.0;
	analysisSettingsPPHT2.mixNormEntriesBack = false;
	analysisSettingsPPHT2.subtractBg = false;
	analysisSettingsPPHT2.subtractBgRandom = false;
	analysisSettingsPPHT2.subtractBgBack = false;
	analysisSettingsPPHT2.subtractBgLowMass = true;
	analysisSettingsPPHT2.subtractBgLowMassPeakShape = true;
	analysisSettingsPPHT2.subtractBgLowMassPeakShapeEta = true;
	analysisSettingsPPHT2.subtractBgFitDistribution = true;
	analysisSettingsPPHT2.useEventMixingBgRandom = true;
	analysisSettingsPPHT2.useEventMixingNotmatchedBgRandom = false;
	analysisSettingsPPHT2.useEventJetMixingBg = true;
	analysisSettingsPPHT2.useEventJetNotmatchedMixingBg = false;
	analysisSettingsPPHT2.useEventJetBackMixingBg = false;
	analysisSettingsPPHT2.useEventShufflingBg = false;
	analysisSettingsPPHT2.fitBgDistribution = true;
	analysisSettingsPPHT2.fitBgRandomDistribution = true;
	analysisSettingsPPHT2.fitBgLowMassDistribution = false;
	analysisSettingsPPHT2.fitPeakShapeDistribution = false && !isPythia;
	analysisSettingsPPHT2.fitPeakShapeEtaDistribution = false && !isPythia;
	analysisSettingsPPHT2.fitPeakShapeEtabgDistribution = false && !isPythia;
	analysisSettingsPPHT2.fitDistributionLeft = 0.30;
	analysisSettingsPPHT2.fitDistributionRight = 0.4;//0.45;  // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsPPHT2.fitDistributionLeft2 = 0.9;//0.75;  // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsPPHT2.fitDistributionRight2 = 1.2;
	analysisSettingsPPHT2.fitDistributionOption = /*"RQN UNORM"; //*/ "RQN";  // towerMB_smd2HT1 / towerMB_smd2HT1_gammaetacut
	analysisSettingsPPHT2.lowNormFixed = TFitSettings(1.0);
	analysisSettingsPPHT2.lowNormMassRange = false;
	analysisSettingsPPHT2.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPHT2.lowNormMassRangeRight = 0.06;
	analysisSettingsPPHT2.lowNormPoints = false;
	analysisSettingsPPHT2.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPHT2.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPHT2.lowNormPointsPt = true;
	analysisSettingsPPHT2.lowmassbgNbar = TFitSettings(0.0, 0, 0, 1, false, 0, 100, 0, 100, "(1.0/((1.0/((0.2689*x*TMath::Exp(-0.3258*x) / ((7.334427e-01)+(-9.284123e-02*x)+(8.582639e-03*x*x)+(-2.596326e-04*x*x*x))) * (1-0.2+(0.2*0.358/0.639)) * (2.886*TMath::Landau(x,5.595,0.9)) / (1.0+(3.01325e-03*x*x)) ))+1.0)) + (0*[0])", "RQN");
	analysisSettingsPPHT2.correctVertexFindingEff = false;
	analysisSettingsPPHT2.vertexFindingEff = 1.0;
	analysisSettingsPPHT2.correctVertexFindingEffMB = true;
	analysisSettingsPPHT2.triggerThreshold = 3.5 - 0.5;
	analysisSettingsPPHT2.showInvBinsFitsSeparately = false;

	TAnalysisSettingsTrigger analysisSettingsPPSimMB = analysisSettingsPPMB;
	if (extendMCpT) analysisSettingsPPSimMB.highPt = analysisSettingsPPHT2.highPt;
	analysisSettingsPPSimMB.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPSimMB.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.7, 0.025 * 1.5, false, 0, 100, 0, 100, "(2.116816e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimMB.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPSimMB.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.8, 0.06 * 1.5, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimMB.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPSimMB.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPSimMB.getJetBgFractionFromCandidatesPtFit = isPythia;
	if (!isPythia) analysisSettingsPPSimMB.fitFunctionStr = fitFuncZero;
	if (!isPythia) analysisSettingsPPSimMB.jetBgFraction.useFit = false;
	if (!isPythia) analysisSettingsPPSimMB.bgToSigbg.useFit = false;
	analysisSettingsPPSimMB.fixedPrescale = analysisSettingsPPMB.fixedPrescale * 0.25/(1.102*1.916);
	analysisSettingsPPSimMB.correctCpv = isPythia;
	analysisSettingsPPSimMB.correctCpvEta = isPythia;
	analysisSettingsPPSimMB.smearOutliersSigma = -100;
	analysisSettingsPPSimMB.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsPPSimMB.lowNormMassRange = false;
	analysisSettingsPPSimMB.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPSimMB.lowNormMassRangeRight = 0.06;
	analysisSettingsPPSimMB.lowNormPoints = false;
	analysisSettingsPPSimMB.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPSimMB.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPSimMB.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsPPSimMB.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsPPSimMB.name = analysisSettingsSpectrumSimMB.name;
	analysisSettingsPPSimMB.title = analysisSettingsSpectrumSimMB.title;
        analysisSettingsPPSimMB.drawOptions = analysisSettingsSpectrumSimMB.drawOptions;
        analysisSettingsPPSimMB.drawOptionsEta = analysisSettingsSpectrumSimMB.drawOptionsEta;

	TAnalysisSettingsTrigger analysisSettingsPPSimHT1 = analysisSettingsPPHT1;
	if (extendMCpT) analysisSettingsPPSimHT1.highPt = analysisSettingsPPHT2.highPt;
	analysisSettingsPPSimHT1.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPSimHT1.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(2.116816e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimHT1.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPSimHT1.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.5, 0.06 * 2.0, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	if (!isPythia) analysisSettingsPPSimHT1.fitFunctionStr = fitFuncZero;
	if (!isPythia) analysisSettingsPPSimHT1.jetBgFraction.useFit = false;
	if (!isPythia) analysisSettingsPPSimHT1.bgToSigbg.useFit = false;
	analysisSettingsPPSimHT1.fixedPrescale = analysisSettingsPPHT1.fixedPrescale * 2.5/(1.109*1.076);
	analysisSettingsPPSimHT1.calculatePrescaleFromDB = false;
	analysisSettingsPPSimHT1.calculatePrescaleFromSim = false;
	analysisSettingsPPSimHT1.calculatePrescaleFromPoints = false;
	analysisSettingsPPSimHT1.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsPPSimHT1.prescalePointsLowPt = 7.0;
	analysisSettingsPPSimHT1.prescalePointsHighPt = 27.0;
	analysisSettingsPPSimHT1.prescalePointsLowPtIntegral = 7.0;
	analysisSettingsPPSimHT1.prescalePointsHighPtIntegral = 27.0;
	analysisSettingsPPSimHT1.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPSimHT1.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPSimHT1.getJetBgFractionFromCandidatesPtFit = isPythia;
	analysisSettingsPPSimHT1.correctCpv = isPythia;
	analysisSettingsPPSimHT1.correctCpvEta = isPythia;
	analysisSettingsPPSimHT1.smearOutliersSigma = -100;
	analysisSettingsPPSimHT1.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsPPSimHT1.lowNormMassRange = false;
	analysisSettingsPPSimHT1.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPSimHT1.lowNormMassRangeRight = 0.06;
	analysisSettingsPPSimHT1.lowNormPoints = false;
	analysisSettingsPPSimHT1.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPSimHT1.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPSimHT1.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsPPSimHT1.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsPPSimHT1.name = analysisSettingsSpectrumSimHT1.name;
	analysisSettingsPPSimHT1.title = analysisSettingsSpectrumSimHT1.title;
        analysisSettingsPPSimHT1.drawOptions = analysisSettingsSpectrumSimHT1.drawOptions;
        analysisSettingsPPSimHT1.drawOptionsEta = analysisSettingsSpectrumSimHT1.drawOptionsEta;

	TAnalysisSettingsTrigger analysisSettingsPPSimHT2 = analysisSettingsPPHT2;
	if (extendMCpT) analysisSettingsPPSimHT2.highPt = analysisSettingsPPHT2.highPt;
	analysisSettingsPPSimHT2.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.7, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPSimHT2.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(2.116816e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimHT2.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPSimHT2.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.6, 0.06 * 1.2, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	if (!isPythia) analysisSettingsPPSimHT2.fitFunctionStr = fitFuncZero;
	if (!isPythia) analysisSettingsPPSimHT2.jetBgFraction.useFit = false;
	if (!isPythia) analysisSettingsPPSimHT2.bgToSigbg.useFit = false;
	analysisSettingsPPSimHT2.fixedPrescale = analysisSettingsPPHT2.fixedPrescale * 2.5/(0.954*1.325);
	analysisSettingsPPSimHT2.calculatePrescaleFromDB = false;
	analysisSettingsPPSimHT2.calculatePrescaleFromSim = false;
	analysisSettingsPPSimHT2.calculatePrescaleFromPoints = false;
	analysisSettingsPPSimHT2.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsPPSimHT2.prescalePointsLowPt = 10.0;
	analysisSettingsPPSimHT2.prescalePointsHighPt = 27.0;
	analysisSettingsPPSimHT2.prescalePointsLowPtIntegral = 10.0;
	analysisSettingsPPSimHT2.prescalePointsHighPtIntegral = 27.0;
	analysisSettingsPPSimHT2.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPSimHT2.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPSimHT2.getJetBgFractionFromCandidatesPtFit = isPythia;
	analysisSettingsPPSimHT2.correctCpv = isPythia;
	analysisSettingsPPSimHT2.correctCpvEta = isPythia;
	analysisSettingsPPSimHT2.smearOutliersSigma = -100;
	analysisSettingsPPSimHT2.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsPPSimHT2.lowNormMassRange = false;
	analysisSettingsPPSimHT2.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPSimHT2.lowNormMassRangeRight = 0.06;
	analysisSettingsPPSimHT2.lowNormPoints = false;
	analysisSettingsPPSimHT2.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPSimHT2.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPSimHT2.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsPPSimHT2.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsPPSimHT2.name = analysisSettingsSpectrumSimHT2.name;
	analysisSettingsPPSimHT2.title = analysisSettingsSpectrumSimHT2.title;
        analysisSettingsPPSimHT2.drawOptions = analysisSettingsSpectrumSimHT2.drawOptions;
        analysisSettingsPPSimHT2.drawOptionsEta = analysisSettingsSpectrumSimHT2.drawOptionsEta;

	TAnalysisSettingsTrigger analysisSettingsPPSimEtaMB = analysisSettingsPPMB;
	if (extendMCpT) analysisSettingsPPSimEtaMB.highPt = analysisSettingsPPHT2.highPt;
	analysisSettingsPPSimEtaMB.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 12, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaMB.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.7, 0.025 * 1.5, false, 0, 100, 0, 100, "(2.116816e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaMB.peakLeftEta = TFitSettings(0.3, 0.02);
	analysisSettingsPPSimEtaMB.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaMB.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.8, 0.06 * 1.5, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	if (!isPythia) analysisSettingsPPSimEtaMB.fitFunctionStr = fitFuncZero;
	if (!isPythia) analysisSettingsPPSimEtaMB.jetBgFraction.useFit = false;
	if (!isPythia) analysisSettingsPPSimEtaMB.bgToSigbg.useFit = false;
	analysisSettingsPPSimEtaMB.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPSimEtaMB.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPSimEtaMB.getJetBgFractionFromCandidatesPtFit = isPythia;
	analysisSettingsPPSimEtaMB.fixedPrescale = analysisSettingsPPMB.fixedPrescale * 0.25 * 0.5/(0.198*0.683/0.5);
	analysisSettingsPPSimEtaMB.correctCpv = isPythia;
	analysisSettingsPPSimEtaMB.correctCpvEta = isPythia;
	analysisSettingsPPSimEtaMB.smearOutliersSigma = -100;
	analysisSettingsPPSimEtaMB.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsPPSimEtaMB.lowNormMassRange = false;
	analysisSettingsPPSimEtaMB.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPSimEtaMB.lowNormMassRangeRight = 0.06;
	analysisSettingsPPSimEtaMB.lowNormPoints = false;
	analysisSettingsPPSimEtaMB.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPSimEtaMB.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPSimEtaMB.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsPPSimEtaMB.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsPPSimEtaMB.name = analysisSettingsSpectrumSimEtaMB.name;
	analysisSettingsPPSimEtaMB.title = analysisSettingsSpectrumSimEtaMB.title;
        analysisSettingsPPSimEtaMB.drawOptions = analysisSettingsSpectrumSimEtaMB.drawOptions;
        analysisSettingsPPSimEtaMB.drawOptionsEta = analysisSettingsSpectrumSimEtaMB.drawOptionsEta;

	TAnalysisSettingsTrigger analysisSettingsPPSimEtaHT1 = analysisSettingsPPHT1;
	if (extendMCpT) analysisSettingsPPSimEtaHT1.highPt = analysisSettingsPPHT2.highPt;
	analysisSettingsPPSimEtaHT1.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.8, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaHT1.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(2.116816e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaHT1.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaHT1.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.5, 0.06 * 2.0, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	if (!isPythia) analysisSettingsPPSimEtaHT1.fitFunctionStr = fitFuncZero;
	if (!isPythia) analysisSettingsPPSimEtaHT1.jetBgFraction.useFit = false;
	if (!isPythia) analysisSettingsPPSimEtaHT1.bgToSigbg.useFit = false;
	analysisSettingsPPSimEtaHT1.fixedPrescale = analysisSettingsPPHT1.fixedPrescale * 2.5 * 0.5/(0.246*0.5192/0.5);
	analysisSettingsPPSimEtaHT1.calculatePrescaleFromDB = false;
	analysisSettingsPPSimEtaHT1.calculatePrescaleFromSim = false;
	analysisSettingsPPSimEtaHT1.calculatePrescaleFromPoints = false;
	analysisSettingsPPSimEtaHT1.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsPPSimEtaHT1.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPSimEtaHT1.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPSimEtaHT1.getJetBgFractionFromCandidatesPtFit = isPythia;
	analysisSettingsPPSimEtaHT1.correctCpv = isPythia;
	analysisSettingsPPSimEtaHT1.correctCpvEta = isPythia;
	analysisSettingsPPSimEtaHT1.smearOutliersSigma = -100;
	analysisSettingsPPSimEtaHT1.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsPPSimEtaHT1.lowNormMassRange = false;
	analysisSettingsPPSimEtaHT1.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPSimEtaHT1.lowNormMassRangeRight = 0.06;
	analysisSettingsPPSimEtaHT1.lowNormPoints = false;
	analysisSettingsPPSimEtaHT1.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPSimEtaHT1.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPSimEtaHT1.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsPPSimEtaHT1.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsPPSimEtaHT1.name = analysisSettingsSpectrumSimEtaHT1.name;
	analysisSettingsPPSimEtaHT1.title = analysisSettingsSpectrumSimEtaHT1.title;
        analysisSettingsPPSimEtaHT1.drawOptions = analysisSettingsSpectrumSimEtaHT1.drawOptions;
        analysisSettingsPPSimEtaHT1.drawOptionsEta = analysisSettingsSpectrumSimEtaHT1.drawOptionsEta;

	TAnalysisSettingsTrigger analysisSettingsPPSimEtaHT2 = analysisSettingsPPHT2;
	if (extendMCpT) analysisSettingsPPSimEtaHT2.highPt = analysisSettingsPPHT2.highPt;
	analysisSettingsPPSimEtaHT2.peakPosition = TFitSettings(truePionMass, 0, truePionMass * 0.7, truePionMass * 1.5, false, 0, 100, 0, 100, "((1.377651e-01) + (-1.082575e-03 * x) + (1.665054e-04 *x * x)) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaHT2.peakWidth = TFitSettings(0.025, 0, 0.025 * 0.8, 0.025 * 2.0, false, 0, 100, 0, 100, "(2.116816e-02) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaHT2.peakPositionEta = TFitSettings(trueEtaMass, 0, trueEtaMass * 0.95, trueEtaMass * 1.3, false, 0.0, 100.0, 0.0, 100.0, "(5.493869e-01) + (0*[0])", "RQN");
	analysisSettingsPPSimEtaHT2.peakWidthEta = TFitSettings(0.06, 0, 0.06 * 0.6, 0.06 * 1.2, false, 0.0, 100.0, 0.0, 100.0, "(4.808745e-02) + (0*[0])", "RQN");
	if (!isPythia) analysisSettingsPPSimEtaHT2.fitFunctionStr = fitFuncZero;
	if (!isPythia) analysisSettingsPPSimEtaHT2.jetBgFraction.useFit = false;
	if (!isPythia) analysisSettingsPPSimEtaHT2.bgToSigbg.useFit = false;
	analysisSettingsPPSimEtaHT2.fixedPrescale = analysisSettingsPPHT2.fixedPrescale * 2.5*0.5/(0.170*0.5604/0.5);
	analysisSettingsPPSimEtaHT2.calculatePrescaleFromDB = false;
	analysisSettingsPPSimEtaHT2.calculatePrescaleFromSim = false;
	analysisSettingsPPSimEtaHT2.calculatePrescaleFromPoints = false;
	analysisSettingsPPSimEtaHT2.calculatePrescaleFromPointsIntegral = false;
	analysisSettingsPPSimEtaHT2.fitJetBgFractionFromCandidates = false;
	analysisSettingsPPSimEtaHT2.getJetBgFractionFromCandidatesPt = false;
	analysisSettingsPPSimEtaHT2.getJetBgFractionFromCandidatesPtFit = isPythia;
	analysisSettingsPPSimEtaHT2.correctCpv = isPythia;
	analysisSettingsPPSimEtaHT2.correctCpvEta = isPythia;
	analysisSettingsPPSimEtaHT2.smearOutliersSigma = -100;
	analysisSettingsPPSimEtaHT2.lowNormFixed = TFitSettings(-1.0);
	analysisSettingsPPSimEtaHT2.lowNormMassRange = false;
	analysisSettingsPPSimEtaHT2.lowNormMassRangeLeft = 0.00;
	analysisSettingsPPSimEtaHT2.lowNormMassRangeRight = 0.06;
	analysisSettingsPPSimEtaHT2.lowNormPoints = false;
	analysisSettingsPPSimEtaHT2.lowmassbgPointsLowPt = 1.0;
	analysisSettingsPPSimEtaHT2.lowmassbgPointsHighPt = 100.0;
	analysisSettingsPPSimEtaHT2.lowNormPointsPt = true;
	if (!isPythia) analysisSettingsPPSimEtaHT2.lowmassbgNbar = TFitSettings(0, 0);
	analysisSettingsPPSimEtaHT2.name = analysisSettingsSpectrumSimEtaHT2.name;
	analysisSettingsPPSimEtaHT2.title = analysisSettingsSpectrumSimEtaHT2.title;
        analysisSettingsPPSimEtaHT2.drawOptions = analysisSettingsSpectrumSimEtaHT2.drawOptions;
        analysisSettingsPPSimEtaHT2.drawOptionsEta = analysisSettingsSpectrumSimEtaHT2.drawOptionsEta;

	TAnalysisSettings analysisSettingsPP = analysisSettingsSpectrum;
	analysisSettingsPP.settingsMB = analysisSettingsPPMB;
	analysisSettingsPP.settingsHT1 = analysisSettingsPPHT1;
	analysisSettingsPP.settingsHT2 = analysisSettingsPPHT2;
	analysisSettingsPP.settingsSimMB = analysisSettingsPPSimMB;
	analysisSettingsPP.settingsSimHT1 = analysisSettingsPPSimHT1;
	analysisSettingsPP.settingsSimHT2 = analysisSettingsPPSimHT2;
	analysisSettingsPP.settingsSimEtaMB = analysisSettingsPPSimEtaMB;
	analysisSettingsPP.settingsSimEtaHT1 = analysisSettingsPPSimEtaHT1;
	analysisSettingsPP.settingsSimEtaHT2 = analysisSettingsPPSimEtaHT2;
	analysisSettingsPP.smoothEffMB = false;
	analysisSettingsPP.smoothEffHT1 = false;
	analysisSettingsPP.smoothEffHT2 = false;
	//analysisSettingsPP.effFitMB = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 10.0, 0.0, 100.0, "(0.145) * 0.5 * (1 + TMath::Erf(pol1))", "RQN");
	//analysisSettingsPP.effFitHT1 = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0, 12, 0.0, 100.0, "(0.145) * 0.5 * (1 + TMath::Erf(pol1))", "RQN");
	//analysisSettingsPP.effFitHT2 = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0, 14, 0.0, 100.0, "(0.145) * 0.5 * (1 + TMath::Erf(pol1))", "RQN");
	if (isPythia) {
	    analysisSettingsPP.effFitMB = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 100.0, "pol2(0) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT1 = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 100.0, "((5.77197e-02) + (4.71842e-03*x) + (-3.17052e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT2 = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 100.0, "((5.77197e-02) + (4.71842e-03*x) + (-3.17052e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	} else {
	    analysisSettingsPP.effFitMB = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "pol2(0) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT1 = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "((5.77197e-02) + (4.71842e-03*x) + (-3.17052e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT2 = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "((5.77197e-02) + (4.71842e-03*x) + (-3.17052e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	}
	//analysisSettingsPP.effFitHT1 = TFitSettings(-7.987947e-02, 2.693218e-02, 0.0, 1.0, true, 8, 13, 8, 100.0, "pol0", "RQN");
	//analysisSettingsPP.effFitHT2 = TFitSettings(-1.079464e-01, 2.540717e-02, 0.0, 1.0, true, 8, 18, 8, 100.0, "pol0", "RQN");
	analysisSettingsPP.dontCorrectEff = false;
	analysisSettingsPP.smoothEffMBEta = false;
	analysisSettingsPP.smoothEffHT1Eta = false;
	analysisSettingsPP.smoothEffHT2Eta = false;
	//analysisSettingsPP.effFitMBEta = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0, 14, 0.0, 100.0, "0.5*(0.2) * (1 + TMath::Erf(pol1))", "RQN");
	//analysisSettingsPP.effFitHT1Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0, 16, 0.0, 100.0, "0.5*(0.2) * (1 + TMath::Erf(pol1))", "RQN");
	//analysisSettingsPP.effFitHT2Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0, 18, 0.0, 100.0, "0.5*(0.2) * (1 + TMath::Erf(pol1))", "RQN");
	if (isPythia) {
	    analysisSettingsPP.effFitMBEta = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 100.0, "0.5*((1.022050e-01) + (1.764432e-02*x) + (-8.041476e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT1Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 100.0, "0.5*((1.022050e-01) + (1.764432e-02*x) + (-8.041476e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT2Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, true, 0.0, 100.0, 0.0, 100.0, "0.5*((1.022050e-01) + (1.764432e-02*x) + (-8.041476e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	} else {
	    analysisSettingsPP.effFitMBEta = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "0.5*((1.022050e-01) + (1.764432e-02*x) + (-8.041476e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT1Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "0.5*((1.022050e-01) + (1.764432e-02*x) + (-8.041476e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	    analysisSettingsPP.effFitHT2Eta = TFitSettings(0.0, 0.0, 0.0, 1.0, false, 0.0, 100.0, 0.0, 100.0, "0.5*((1.022050e-01) + (1.764432e-02*x) + (-8.041476e-04*x*x)) * (1 + TMath::Erf(pol1(3)))", "RQN");
	}
	analysisSettingsPP.dontCorrectEffEta = false;
	analysisSettingsPP.cpvCorrection = TFitSettings(1.0, 0.0, 0.0, 1.0, true, -100, +100, -100, +100, "pol1", "RQNW");
	analysisSettingsPP.prescalesDBFilename = "./pp2005Prescales.txt";
	analysisSettingsPP.showWeightUsed = true;
	analysisSettingsPP.showWeightNextIteration = true;
	analysisSettingsPP.saveWeightNextIteration = false;
	analysisSettingsPP.tpcVertexEffFunc = "[0] + ([1]*x)";
	analysisSettingsPP.showTPCVertexEff = true;
	analysisSettingsPP.saveTPCVertexEff = true;
	analysisSettingsPP.bbcVertexEffFunc = "[0]";
	analysisSettingsPP.showBBCVertexEff = false;
	analysisSettingsPP.saveBBCVertexEff = false;
	analysisSettingsPP.saveRunTimes = false;
	analysisSettingsPP.runTimesFilenameFormat = "./pp2005Times_%s.txt";
	analysisSettingsPP.saveRunTimesSim = false;
	analysisSettingsPP.runTimesSimFilenameFormat = "./pp2005Times_pi0_%s.txt";
	analysisSettingsPP.saveRunTimesSimEta = false;
	analysisSettingsPP.runTimesSimEtaFilenameFormat = "./pp2005Times_eta_%s.txt";
	analysisSettingsPP.savePrescales = false;
	analysisSettingsPP.prescalesSaveFilename = "./pp2005Prescales.txt";
	analysisSettingsPP.saveBadSimPi0EventsMB = false;
	analysisSettingsPP.badSimPi0EventsMBFilename = "./badEvents_pythiaPi0Pt_pp2005_MB.txt";
	analysisSettingsPP.saveBadSimPi0EventsHT1 = false;
	analysisSettingsPP.badSimPi0EventsHT1Filename = "./badEvents_pythiaPi0Pt_pp2005_HT1.txt";
	analysisSettingsPP.saveBadSimPi0EventsHT2 = false;
	analysisSettingsPP.badSimPi0EventsHT2Filename = "./badEvents_pythiaPi0Pt_pp2005_HT2.txt";
	analysisSettingsPP.dbConnectString = "mysql://dbx.star.bnl.gov:3316/RunLog_onl";
	analysisSettingsPP.weightFile = "./weights_pp2005_new.root";
	analysisSettingsPP.NBinaryCollisions = 1.0;
	analysisSettingsPP.NBinaryCollisionsErr = 0.0;
	analysisSettingsPP.materialCorrectionSim = 1.10;
	analysisSettingsPP.materialCorrectionSimGamma = 1.05;
	analysisSettingsPP.crossectionMB = 26.1; // BBC cross section is 26.1 +/- 0.2(stat) +/- 1.8(sys) mb
	analysisSettingsPP.crossectionMBErr = 0.2;
	analysisSettingsPP.crossectionMBTotalFraction = 0.87; // BBC captured 85% of non-singly diffractive (NSD) events
	analysisSettingsPP.crossectionMBTotalFractionErr = 0.08;
	analysisSettingsPP.crossectionSimMB = isPythia ? 28.12 : (30.0); // PYTHIA crossection: hard processes (partonic pT > 2 GeV/c) - 9.88 mb, minbias - 28.12 mb.
	analysisSettingsPP.crossectionSimMBErr = 0.0;
	analysisSettingsPP.crossectionSimMBTotalFraction = 1.0;
	analysisSettingsPP.crossectionSimMBTotalFractionErr = 0.0;
	analysisSettingsPP.crossectionSimHT1 = isPythia ? 28.12 : (30.0);
	analysisSettingsPP.crossectionSimHT2 = isPythia ? 28.12 : (30.0);
	analysisSettingsPP.crossectionSimEtaMB = (isPythia ? 28.12 : ((30.0)/0.5));
	analysisSettingsPP.crossectionSimEtaMBErr = 0.0;
	analysisSettingsPP.crossectionSimEtaMBTotalFraction = 1.0;
	analysisSettingsPP.crossectionSimEtaMBTotalFractionErr = 0.0;
	analysisSettingsPP.crossectionSimEtaHT1 = (isPythia ? 28.12 : ((30.0)/0.5));
	analysisSettingsPP.crossectionSimEtaHT2 = (isPythia ? 28.12 : ((30.0)/0.5));
	analysisSettingsPP.crossectionInel = 42.2; // Total p+p inelastic cross section is 42.2 +/- xx mb at 200 GeV
	analysisSettingsPP.crossectionInelErr = 0;
	analysisSettingsPP.energyScaleErr = 0.05;
	analysisSettingsPP.energyScaleErrFitOption = "RQN LL";
	analysisSettingsPP.energyScaleErrFitOptionWeight = "POWERLAW POLYNOMIAL POL5";
	analysisSettingsPP.energyScaleErrFitOptionWeightRaw = "POWERLAW POLYNOMIAL POL3";
	analysisSettingsPP.usePPpQCD = true;
	analysisSettingsPP.name = TString("pp") + (addDirName ? TString("_") + DATA_DIR : TString());
	analysisSettingsPP.title = TString("#font[12]{p}+#font[12]{p}");
        analysisSettingsPP.drawOptions    = TDrawOptions(TNamed("pp", "#font[12]{p}+#font[12]{p}"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullCircle, 1.5), "P", "P", TString("#font[12]{p}+#font[12]{p} #pi^{0}"));
        analysisSettingsPP.drawOptionsEta = TDrawOptions(TNamed("pp", "#font[12]{p}+#font[12]{p}"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleUp, 1.5), "P", "P", TString("#font[12]{p}+#font[12]{p} #eta"));

	TAllSettings *p_allSettings = new TAllSettings();
	TAllSettings &allSettings = *p_allSettings;
	allSettings.settingsDAuNoCentral = analysisSettingsSpectrum;
	////allSettings.settingsDAuNoCentral.settingsMB.bgToSigbg = TFitSettings(0.86, 0.0, 0.0, 100.0, true, 1.0, 3.0, 0.0, +100.0, "pol1", "RQN");
	//allSettings.settingsDAuNoCentral.settingsMB.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 1.0, 2.5, 0.0, +100.0, "pol1", "RQN");
	////allSettings.settingsDAuNoCentral.settingsHT1.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 0, 8, 0, 8, "pol2", "RQN");
	////allSettings.settingsDAuNoCentral.settingsHT2.bgToSigbg = TFitSettings(0.536, -0.006, 0.0, 1.0, true, 0, 14, 0, 100, "pol0", "RQN");
	//allSettings.settingsDAuNoCentral.settingsHT1.signalToBackgroundRatio = TFitSettings(1.0, 0.0, 0.0, 100.0, true, 0, 10, 0, 10, "pol2", "RQN");
	//allSettings.settingsDAuNoCentral.settingsHT1.signalToBackgroundRatioEta = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 0, 10, 0, 100, "pol2", "RQN");
	//allSettings.settingsDAuNoCentral.settingsHT2.signalToBackgroundRatio = TFitSettings(1.2, 0.0, 0.0, 100.0, true, 0, 12, 0, 12, "pol2", "RQN");
	//allSettings.settingsDAuNoCentral.settingsHT2.signalToBackgroundRatioEta = TFitSettings(0.27, 0.0, 0.0, 100.0, true, 0, 17, 0, 100, "pol0", "RQN");
	allSettings.settingsDAuNoCentral.correctPrescaleTriggerBias = true;
	allSettings.settingsDAuNoCentral.prescaleTriggerBiasMBHT1 = -0.02;
	allSettings.settingsDAuNoCentral.prescaleTriggerBiasMBHT2 = -0.20;
	allSettings.settingsDAuNoCentral.prescaleTriggerBiasHT1HT2 = -0.25;
	allSettings.settingsDAuNoCentral.showSim = true;
	allSettings.settingsDAuAllCentral = analysisSettingsRcp;
	allSettings.settingsDAuAllCentral.correctPrescaleTriggerBias = false;
	allSettings.settingsDAuAllCentral.prescaleTriggerBiasMBHT1 = -0.02;
	allSettings.settingsDAuAllCentral.prescaleTriggerBiasMBHT2 = -0.20;
	allSettings.settingsDAuAllCentral.prescaleTriggerBiasHT1HT2 = -0.25;
	allSettings.settingsDAuAllCentral.settingsMB.correctCpv = true;
	allSettings.settingsDAuAllCentral.settingsMB.correctCpvEta = true;
	allSettings.settingsDAuAllCentral.settingsHT1.correctCpv = true;
	allSettings.settingsDAuAllCentral.settingsHT1.correctCpvEta = true;
	allSettings.settingsDAuAllCentral.settingsHT2.correctCpv = true;
	allSettings.settingsDAuAllCentral.settingsHT2.correctCpvEta = true;
	allSettings.settingsDAuAllCentral.showSim = false;
	allSettings.settingsDAuAllCentral.name = TString("dAuNocentral") + (addDirName ? TString("_") + DATA_DIR : TString());
	allSettings.settingsDAuAllCentral.title = TString("#font[12]{d}+Au 0-100%");
        allSettings.settingsDAuAllCentral.drawOptions    = TDrawOptions(TNamed("dAuNocentral", "#font[12]{d}+Au 0-100%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.5), "P", "P", TString("#font[12]{d}+Au 0-100% #pi^{0}"));
        allSettings.settingsDAuAllCentral.drawOptionsEta = TDrawOptions(TNamed("dAuNocentral", "#font[12]{d}+Au 0-100%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenTriangleUp, 1.5), "P", "P", TString("#font[12]{d}+Au 0-100% #eta"));
	allSettings.settingsDAuMostCentral = analysisSettingsRcp;
	allSettings.settingsDAuMostCentral.correctPrescaleTriggerBias = false;
	allSettings.settingsDAuMostCentral.prescaleTriggerBiasMBHT1 = -0.02;
	allSettings.settingsDAuMostCentral.prescaleTriggerBiasMBHT2 = -0.20;
	allSettings.settingsDAuMostCentral.prescaleTriggerBiasHT1HT2 = -0.25;
	////allSettings.settingsDAuMostCentral.settingsMB.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 0, 5, 0.0, +100.0, "pol2", "RQN");
	allSettings.settingsDAuMostCentral.settingsMB.correctVertexFindingEff = false;
	allSettings.settingsDAuMostCentral.settingsMB.vertexFindingEff = 1.0;
	allSettings.settingsDAuMostCentral.settingsMB.correctCpv = true;
	allSettings.settingsDAuMostCentral.settingsMB.correctCpvEta = true;
	allSettings.settingsDAuMostCentral.settingsMB.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 1, 3, 0, 100, "pol0", "RQNW");
	////allSettings.settingsDAuMostCentral.settingsHT1.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 4, 8, 0.0, +100.0, "pol0", "RQN");
	////allSettings.settingsDAuMostCentral.settingsHT1.signalToBackgroundRatio = TFitSettings(4.0, 0.0, 0.0, 100.0, true, 0, 100, 0, 0, "pol1", "RQN");
	allSettings.settingsDAuMostCentral.settingsHT1.correctCpv = true;
	allSettings.settingsDAuMostCentral.settingsHT1.correctCpvEta = true;
	allSettings.settingsDAuMostCentral.settingsHT1.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 4, 8, 0, 100, "pol0", "RQNW");
	allSettings.settingsDAuMostCentral.settingsHT1.fixedPrescale = 1524.52;
	////allSettings.settingsDAuMostCentral.settingsHT2.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 6, 10, 0, 100, "pol0", "RQN");
	////allSettings.settingsDAuMostCentral.settingsHT2.signalToBackgroundRatio = TFitSettings(5.5, -0.3, 0.0, 100.0, true, 0, 100, 0, 0, "pol0", "RQN");
	allSettings.settingsDAuMostCentral.settingsHT2.correctCpv = true;
	allSettings.settingsDAuMostCentral.settingsHT2.correctCpvEta = true;
	allSettings.settingsDAuMostCentral.settingsHT2.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 6, 12, 0, 100, "pol0", "RQN");
	allSettings.settingsDAuMostCentral.settingsHT2.fixedPrescale = 16336.4;
	allSettings.settingsDAuMostCentral.showSim = false;
	allSettings.settingsDAuMostCentral.name = TString("dAuMostcentral") + (addDirName ? TString("_") + DATA_DIR : TString());
	allSettings.settingsDAuMostCentral.title = TString("#font[12]{d}+Au 0-20%");
        allSettings.settingsDAuMostCentral.drawOptions    = TDrawOptions(TNamed("dAuMostcentral", "#font[12]{d}+Au 0-20%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kPlus, 1.5), "P", "P", TString("#font[12]{d}+Au 0-20% #pi^{0}"));
        allSettings.settingsDAuMostCentral.drawOptionsEta = TDrawOptions(TNamed("dAuMostcentral", "#font[12]{d}+Au 0-20%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kPlus, 1.5), "P", "P", TString("#font[12]{d}+Au 0-20% #eta"));
	allSettings.settingsDAuMidCentral = analysisSettingsRcp;
	allSettings.settingsDAuMidCentral.correctPrescaleTriggerBias = false;
	allSettings.settingsDAuMidCentral.prescaleTriggerBiasMBHT1 = -0.02;
	allSettings.settingsDAuMidCentral.prescaleTriggerBiasMBHT2 = -0.20;
	allSettings.settingsDAuMidCentral.prescaleTriggerBiasHT1HT2 = -0.25;
	allSettings.settingsDAuMidCentral.settingsMB.correctVertexFindingEff = false;
	allSettings.settingsDAuMidCentral.settingsMB.vertexFindingEff = 1.0;
	allSettings.settingsDAuMidCentral.settingsMB.correctCpv = true;
	allSettings.settingsDAuMidCentral.settingsMB.correctCpvEta = true;
	allSettings.settingsDAuMidCentral.settingsHT1.correctCpv = true;
	allSettings.settingsDAuMidCentral.settingsHT1.correctCpvEta = true;
	allSettings.settingsDAuMidCentral.settingsHT2.correctCpv = true;
	allSettings.settingsDAuMidCentral.settingsHT2.correctCpvEta = true;
	allSettings.settingsDAuMidCentral.showSim = false;
	allSettings.settingsDAuMidCentral.name = TString("dAuMidcentral") + (addDirName ? TString("_") + DATA_DIR : TString());
	allSettings.settingsDAuMidCentral.title = TString("#font[12]{d}+Au 20-40%");
        allSettings.settingsDAuMidCentral.drawOptions    = TDrawOptions(TNamed("dAuMidcentral", "#font[12]{d}+Au 20-40%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.5), "P", "P", TString("#font[12]{d}+Au 20-40% #pi^{0}"));
        allSettings.settingsDAuMidCentral.drawOptionsEta = TDrawOptions(TNamed("dAuMidcentral", "#font[12]{d}+Au 20-40%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenTriangleUp, 1.5), "P", "P", TString("#font[12]{d}+Au 20-40% #eta"));
	allSettings.settingsDAuMostPeripheral = analysisSettingsRcp;
	allSettings.settingsDAuMostPeripheral.correctPrescaleTriggerBias = false;
	allSettings.settingsDAuMostPeripheral.prescaleTriggerBiasMBHT1 = -0.02;
	allSettings.settingsDAuMostPeripheral.prescaleTriggerBiasMBHT2 = -0.20;
	allSettings.settingsDAuMostPeripheral.prescaleTriggerBiasHT1HT2 = -0.25;
	//allSettings.settingsDAuMostPeripheral.settingsMB.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 0, 4, 0.0, +100.0, "pol1", "RQN");
	allSettings.settingsDAuMostPeripheral.settingsMB.correctVertexFindingEff = true;
	allSettings.settingsDAuMostPeripheral.settingsMB.vertexFindingEff = vertEffDAuMostPeripheralMinBias;
	allSettings.settingsDAuMostPeripheral.settingsMB.correctCpv = true;
	allSettings.settingsDAuMostPeripheral.settingsMB.correctCpvEta = true;
	allSettings.settingsDAuMostPeripheral.settingsMB.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 2, 4, 0, 100, "pol0", "RQN");
	////allSettings.settingsDAuMostPeripheral.settingsHT1.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 4, 8, 0.0, +100.0, "pol1", "RQNW");
	////allSettings.settingsDAuMostPeripheral.settingsHT1.signalToBackgroundRatio = TFitSettings(2.0, 0.0, 0.0, 100.0, true, 6, 12, 6, 12, "pol5", "RQN");
	allSettings.settingsDAuMostPeripheral.settingsHT1.correctCpv = true;
	allSettings.settingsDAuMostPeripheral.settingsHT1.correctCpvEta = true;
	allSettings.settingsDAuMostPeripheral.settingsHT1.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 4, 7, 0, 100, "pol0", "RQN");
	allSettings.settingsDAuMostPeripheral.settingsHT1.fixedPrescale = 2803.91;//4960.25;
	////allSettings.settingsDAuMostPeripheral.settingsHT2.bgToSigbg = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 0, 100, 0.0, +100.0, "pol0", "RQNW");
	////allSettings.settingsDAuMostPeripheral.settingsHT2.signalToBackgroundRatio = TFitSettings(1.0, 0.0, 0.0, 100.0, true, 6, 12, 0, 12, "gaus", "RQN");
	//allSettings.settingsDAuMostPeripheral.settingsHT2.signalToBackgroundRatioEta = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 6, 12, 6, 12, "pol2", "RQN");
	allSettings.settingsDAuMostPeripheral.settingsHT2.correctCpv = true;
	allSettings.settingsDAuMostPeripheral.settingsHT2.correctCpvEta = true;
	allSettings.settingsDAuMostPeripheral.settingsHT2.cpvCorrection = TFitSettings(1.0, 0, 0.0, 1.0, false, 6, 10, 0, 100, "pol0", "RQN");
	allSettings.settingsDAuMostPeripheral.settingsHT2.fixedPrescale = 29681.8;//45297.6;
	allSettings.settingsDAuMostPeripheral.showSim = false;
	allSettings.settingsDAuMostPeripheral.name = TString("dAuMostperipheral") + (addDirName ? TString("_") + DATA_DIR : TString());
	allSettings.settingsDAuMostPeripheral.title = TString("#font[12]{d}+Au 40-100%");
        allSettings.settingsDAuMostPeripheral.drawOptions    = TDrawOptions(TNamed("dAuMostPeripheral", "#font[12]{d}+Au 40-100%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kMultiply, 1.5), "P", "P", TString("#font[12]{d}+Au 40-100% #pi^{0}"));
        allSettings.settingsDAuMostPeripheral.drawOptionsEta = TDrawOptions(TNamed("dAuMostPeripheral", "#font[12]{d}+Au 40-100%"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kMultiply, 1.5), "P", "P", TString("#font[12]{d}+Au 40-100% #eta"));
	allSettings.settingsPP = analysisSettingsPP;
	allSettings.settingsPP.correctPrescaleTriggerBias = false;
	allSettings.settingsPP.prescaleTriggerBiasMBHT1 = 0.02;
	allSettings.settingsPP.prescaleTriggerBiasMBHT2 = 0.20;
	allSettings.settingsPP.prescaleTriggerBiasHT1HT2 = 0.25;
	////allSettings.settingsPP.settingsMB.bgToSigbg = TFitSettings(0.65/*0.4*/, 0.0, 0, 1, true, 1, 3, 1, 3, "pol3", "RQN");
	//allSettings.settingsPP.settingsMB.signalToBackgroundRatioEta = TFitSettings(0.5, 0.0, 0.0, 100.0, true, 0, 4, 0, 100, "pol1", "RQN");
	//allSettings.settingsPP.settingsHT1.bgToSigbg = TFitSettings(0.7, 0.0, 0, 1, true, 0, 10, 0, +100.0, "0.7 - (0.0137 * x) + (0*[0])", "RQN");
	//allSettings.settingsPP.settingsHT1.bgToSigbg = TFitSettings(0.7, 0.0, 0, 1, true, 4, 8, 0, +100.0, "pol1", "RQN");
	////allSettings.settingsPP.settingsHT1.bgToSigbg = TFitSettings(0.94, -0.04, 0, 1, true, 4, 8, 4, 8, "pol3", "RQN");
	//allSettings.settingsPP.settingsHT1.signalToBackgroundRatio = TFitSettings(2.0, 0.0, 0.0, 100.0, true, 3, 11, 0, 100, "pol2", "RQN");
	//allSettings.settingsPP.settingsHT2.bgToSigbg = TFitSettings(0.8, 0.00, 0, 1, true, 0, 14, 0, 100, "0.8 - (0.025 * x) + (0*[0])", "RQN");
	//allSettings.settingsPP.settingsHT2.bgToSigbg = TFitSettings(0.67, 0.00, 0, 1, true, 6, 14, 0, 100, "pol1", "RQN");
	////allSettings.settingsPP.settingsHT2.bgToSigbg = TFitSettings(0.84 + 0.05, -0.028 - 0.006, 0, 1, true, 0, 12, 0, 100, "pol1", "RQN");
	//allSettings.settingsPP.settingsHT2.bgToSigbg = TFitSettings(0.84 + 0.05, -0.028 - 0.006, 0, 1, true, 4, 10, 4, 12, "pol2", "RQN");
	//allSettings.settingsPP.settingsHT2.signalToBackgroundRatio = TFitSettings(2.0, 0.0, 0.0, 100.0, true, 4, 12, 0, 100, "pol2", "RQN");
	allSettings.saveRcpDataArrays = true;
	allSettings.dataArraysRcpFilenameFormat = "./data_arrays_%s.C";
	allSettings.saveRdADataArrays = true;
	allSettings.dataArraysRdAFilenameFormat = "./data_arrays_%s.C";
	allSettings.browseAnalysis = false;
	allSettings.browseResults = false;
	allSettings.dAuCentralityName = TString("dAuRcp") + (addDirName ? TString("_") + DATA_DIR : TString());
	allSettings.dAuCentralityTitle = TString("#font[12]{R}^{#font[12]{d}+Au}_{#font[12]{CP}}");
	allSettings.dAuCentralityShow = true;
	allSettings.dAuCentralityPrint = true;
	allSettings.showDAuCentralityTriggersSeparately = true;
	allSettings.showSpectrumDAu = false;
	allSettings.showRcp = false;
	allSettings.showAllCentral = false;
	allSettings.showMidCentral = false;
	allSettings.useDAuEff = true;
	allSettings.useDAuEmbedding = false;
	allSettings.useDAuEffCentralities = false;
	allSettings.useDAu1gamma = true;
	allSettings.useDAu1gammaSim = true;
	allSettings.useDAu1gammaSimEta = true;
	allSettings.useDAu1gammaCentralities = false;
	allSettings.useDAunbar = true;
	allSettings.useDAunbarSim = isPythia;
	allSettings.useDAunbarSimEta = isPythia;
	allSettings.useDAuEffWeight = true;
	allSettings.showSpectrumDAuEta = false;
	allSettings.showRcpEta = false;
	allSettings.useDAuEffEta = true;
	allSettings.useDAuEffEtabg = false;
	allSettings.useDAuEmbeddingEta = false;
	allSettings.useDAuEffCentralitiesEta = false;
	allSettings.useDAuEffWeightEta = true;
	allSettings.showSpectrumPP = false;
	allSettings.usePPEff = true;
	allSettings.usePP1gamma = true;
	allSettings.usePP1gammaSim = true;
	allSettings.usePP1gammaSimEta = true;
	allSettings.usePPnbar = true;
	allSettings.usePPnbarSim = isPythia;
	allSettings.usePPnbarSimEta = isPythia;
	allSettings.usePPEffWeight = true;
	allSettings.usePPEffPythia = isPythia;
	allSettings.showSpectrumPPEta = false;
	allSettings.usePPEffEta = true;
	allSettings.usePPEffEtabg = false;
	allSettings.usePPEffWeightEta = true;
	allSettings.showRDA = false;
	allSettings.showRDAEta = false;
	allSettings.showRDAGamma = false;
	allSettings.RDAName = TString("RdA") + (addDirName ? TString("_") + DATA_DIR : TString());
	allSettings.RDATitle = TString("#font[12]{R}_{#font[12]{dA}}");
	allSettings.RDAShow = true;
	allSettings.showRDATriggersSeparately = true;

// Data from other experiments, theory etc.
#include "macros/my_results/data_arrays_others.C"

	Bool_t showPPKKP = true; // KKP or DSS fragmentation functions

	// fitted functions used to divide data
	//const Char_t *LevyFunc = "[0]*TMath::Power(1 + ((TMath::Sqrt(x*x + 0.135*0.135) - 0.135) / [1]), -[2])";
	TF1 *fPPPi0Invyield = 0;
	TF1 *fDAuPi0Invyield = 0;
	TF1 *fPPPi0Crossection = 0;
	TF1 *fDAuPi0Crossection = 0;
        TF1 *fNLOpQCD = new TF1("NLOpQCD", &getNLOpQCD, 1, 20.0/*17.5*/, 1);
        fNLOpQCD->SetParameter(0, 0);
        TF1 *fNLOpQCD1 = new TF1("NLOpQCD1", &getNLOpQCD, 3, 20.0/*17.5*/, 1);
        fNLOpQCD1->SetParameter(0, 1);
        TF1 *fNLOpQCD2 = new TF1("NLOpQCD2", &getNLOpQCD, 2, 20.0/*17.5*/, 1);
        fNLOpQCD2->SetParameter(0, 2);
        TF1 *fNLOpQCDPP = new TF1("NLOpQCDPP", &getNLOpQCDPP, 1, 20.0/*19.4*/, 1);
        fNLOpQCDPP->SetParameter(0, 0);
        TF1 *fNLOpQCDPP1 = new TF1("NLOpQCDPP1", &getNLOpQCDPP, 3, 20.0/*19.4*/, 1);
        fNLOpQCDPP1->SetParameter(0, 1);
        TF1 *fNLOpQCDPP2 = new TF1("NLOpQCDPP2", &getNLOpQCDPP, 2, 20.0/*19.4*/, 1);
        fNLOpQCDPP2->SetParameter(0, 2);
        TF1 *fNLOpQCDPP_2 = new TF1("NLOpQCDPP_2", &getNLOpQCDPP_2, 1, 20.0/*19.4*/, 1);
        fNLOpQCDPP_2->SetParameter(0, 0);
        TF1 *fNLOpQCDPP1_2 = new TF1("NLOpQCDPP1_2", &getNLOpQCDPP_2, 3, 20.0/*19.4*/, 1);
        fNLOpQCDPP1_2->SetParameter(0, 1);
        TF1 *fNLOpQCDPP2_2 = new TF1("NLOpQCDPP2_2", &getNLOpQCDPP_2, 2, 20.0/*19.4*/, 1);
        fNLOpQCDPP2_2->SetParameter(0, 2);
        TF1 *fNLOpQCDPPK = new TF1("NLOpQCDPPK", &getNLOpQCDPP, 1.5, 20.0/*19.4*/, 1);
        fNLOpQCDPPK->SetParameter(0, 20);
	//TF1 *fPPPi0CrossectionLevy = 0;
	//TF1 *fDAuPi0CrossectionLevy = 0;
	TF1 *fPPEtaInvyield = 0;
	TF1 *fDAuEtaInvyield = 0;

        TF1 *fDIV = fNLOpQCD;
        TF1 *fDIVPP = showPPKKP ? fNLOpQCDPP : fNLOpQCDPP_2;

/*
#define X_NEW(fD, x0, dx) {cout << (x0) << ": " << (fD)->GetX((fD)->Integral((x0) - (dx), (x0) + (dx))/((dx) + (dx)), (x0) - (dx), (x0) + (dx)) << endl;}
X_NEW(fDIV, 1.5+0.25, 0.25)
X_NEW(fDIV, 2.0+0.25, 0.25)
X_NEW(fDIV, 2.5+0.25, 0.25)
X_NEW(fDIV, 3.0+0.25, 0.25)
X_NEW(fDIV, 3.5+0.25, 0.25)
X_NEW(fDIV, 4.0+0.25, 0.25)
X_NEW(fDIV, 4.5+0.25, 0.25)
X_NEW(fDIV, 5.0+0.25, 0.25)
X_NEW(fDIV, 5.5+0.25, 0.25)
X_NEW(fDIV, 6.0+0.25, 0.25)
X_NEW(fDIV, 6.5+0.25, 0.25)
X_NEW(fDIV, 7.0+0.25, 0.25)
X_NEW(fDIV, 7.5+0.25, 0.25)
X_NEW(fDIV, 8.0+0.25, 0.25)
X_NEW(fDIV, 8.5+0.25, 0.25)
X_NEW(fDIV, 9.0+0.25, 0.25)
X_NEW(fDIV, 9.5+0.25, 0.25)
X_NEW(fDIV, 10.0+0.25, 0.25)
X_NEW(fDIV, 10.5+0.25, 0.25)
X_NEW(fDIV, 11.0+0.25, 0.25)
X_NEW(fDIV, 11.5+0.25, 0.25)
X_NEW(fDIV, 12.0+0.25, 0.25)
X_NEW(fDIV, 12.5+0.25, 0.25)
X_NEW(fDIV, 13.0+0.25, 0.25)
X_NEW(fDIV, 13.5+0.25, 0.25)
X_NEW(fDIV, 14.0+0.25, 0.25)
X_NEW(fDIV, 14.5+0.25, 0.25)
X_NEW(fDIV, 15.0+0.25, 0.25)
X_NEW(fDIV, 15.5+0.25, 0.25)
*/

#define CREATE_ARRAYS_FUNC(NAME, F, F1, F2, LOW, HIGH) \
	Float_t x_##NAME[100]; \
	Float_t x_err_##NAME[100]; \
	Float_t y1_##NAME[100]; \
	Float_t y1_err_##NAME[100]; \
	Float_t y2_##NAME[100]; \
	Float_t y2_err_##NAME[100]; \
	for (Int_t i = 0;i < 100;i++) { \
	    Float_t xi = (LOW) + (Float_t(i)*((HIGH)-(LOW))/(100.0)); \
	    Float_t exi = (0.5*((HIGH)-(LOW))/(100.0)); \
	    exi = 0; \
	    Float_t yi = F->Eval(xi); \
	    Float_t y1i = F1->Eval(xi); \
	    Float_t y2i = F2->Eval(xi); \
	    y1_##NAME[i] = (yi + y1i) / 2.0; \
	    y1_err_##NAME[i] = TMath::Abs(yi - y1i) / 2.0; \
	    y2_##NAME[i] = (yi + y2i) / 2.0; \
	    y2_err_##NAME[i] = TMath::Abs(yi - y2i) / 2.0; \
	    x_##NAME[i] = xi; \
	    x_err_##NAME[i] = exi; \
	}

	CREATE_ARRAYS_FUNC(NLOfill, fNLOpQCD, fNLOpQCD1, fNLOpQCD2, 2.0, 20.0/*17.5*/)
	CREATE_ARRAYS_FUNC(NLOPPfill, fNLOpQCDPP, fNLOpQCDPP1, fNLOpQCDPP2, 2.0, 20.0/*19.4*/)
	CREATE_ARRAYS_FUNC(NLOPPfill_2, fNLOpQCDPP_2, fNLOpQCDPP1_2, fNLOpQCDPP2_2, 2.0, 20.0/*19.4*/)
	
	// PHENIX d+Au pi0 cross section
	TDataPoints PHENIX_dAu_pi0 = TDataPoints(TDrawOptions(TNamed("PHENIX_dAu_pi0", "PHENIX #pi^{0}"/*showPhenixPi0Legend*/), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenSquare, 1.0), "P", "P", TString("PHENIX #pi^{0}^{}"/*showPhenixPi0Legend*/)), Px, Py, Pxe, Pye, sizeof(Px)/sizeof(Px[0]));
	// PHENIX d+Au eta cross section
	TDataPoints PHENIX_dAu_eta = TDataPoints(TDrawOptions(TNamed("PHENIX_dAu_eta", showPhenixPi0EtaLegend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenSquare, 1.0), "P", "P", TString(showPhenixPi0EtaLegend)), Petax, Petay, Petaxe, Petaye, sizeof(Petax)/sizeof(Petax[0]));
	// PHENIX p+p pi0 cross section
	TDataPoints PHENIX_pp_pi0 = TDataPoints(TDrawOptions(TNamed("PHENIX_pp_pi0", "PHENIX #pi^{0}"/*showPhenixPi0PPLegend*/), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.0), "P", "P", TString("PHENIX #pi^{0}^{}"/*showPhenixPi0PPLegend*/)), Pppx, Pppy, Pppxe, Pppye, sizeof(Pppx)/sizeof(Pppx[0]));
	// PHENIX p+p eta cross section
	TDataPoints PHENIX_pp_eta = TDataPoints(TDrawOptions(TNamed("PHENIX_pp_eta", showPhenixPi0PPEtaLegend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenSquare, 1.0), "P", "P", TString(showPhenixPi0PPEtaLegend)), Pppetax, Pppetay, Pppetaxe, Pppetaye, sizeof(Pppetax)/sizeof(Pppetax[0]));
	// PHENIX pi0 RdA
	TDataPoints PHENIX_pi0_RdA = TDataPoints(TDrawOptions(TNamed("PHENIX_Rda_pi0", "PHENIX #pi^{0}"/*showPhenixPi0dALegend*/), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.0), "P", "P", TString("PHENIX #pi^{0}^{}"/*showPhenixPi0dALegend*/)), PdAx, PdAy, PdAxe, PdAye, sizeof(PdAx)/sizeof(PdAx[0]));
	// PHENIX eta RdA
	TDataPoints PHENIX_eta_RdA = TDataPoints(TDrawOptions(TNamed("PHENIX_Rda_eta", "PHENIX #eta"/*showPhenixEtadALegend*/), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenSquare, 1.0), "P", "P", TString("PHENIX #eta"/*showPhenixEtadALegend*/)), PEtadAx, PEtadAy, PEtadAxe, PEtadAye, sizeof(PEtadAx)/sizeof(PEtadAx[0]));
	// PHENIX eta/pi0 ratio in p+p
	TDataPoints PHENIX_etaToPi0_pp = TDataPoints(TDrawOptions(TNamed("PHENIX_etatopi_pp", "PHENIX #font[12]{p}_{ }+_{ }#font[12]{p}"/*showPhenixEtaToPi0Legend*/), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenCircle, 1.0), "P", "P", TString("PHENIX #font[12]{p}_{ }+_{ }#font[12]{p}"/*showPhenixEtaToPi0Legend*/)), PEtaToPi0x, PEtaToPi0y, PEtaToPi0xe, PEtaToPi0ye, sizeof(PEtaToPi0x)/sizeof(PEtaToPi0x[0]));
	// PHENIX eta/pi0 ratio in d+Au
	TDataPoints PHENIX_etaToPi0_dAu = TDataPoints(TDrawOptions(TNamed("PHENIX_etatopi_dAu", "PHENIX #font[12]{d}_{ }+_{ }Au"/*showPhenixEtaToPi0dAuLegend*/), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kOpenSquare, 1.0), "P", "P", TString("PHENIX #font[12]{d}_{ }+_{ }Au"/*showPhenixEtaToPi0dAuLegend*/)), PEtaToPi0dAux, PEtaToPi0dAuy, PEtaToPi0dAuxe, PEtaToPi0dAuye, sizeof(PEtaToPi0dAux)/sizeof(PEtaToPi0dAux[0]));
	// STAR d+Au pi+
	TDataPoints STAR_dAu_piplus = TDataPoints(TDrawOptions(TNamed("STAR_dAu_piplus", "STAR #pi^{+}, #pi^{-}"/*showStarChargedPionsPlusLegend*/), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenStar, 1.3), "P", "P", TString("STAR #pi^{+}^{}, #pi^{#font[122]{-}}^{}"/*showStarChargedPionsPlusLegend*/)), SPpx, SPpy_xsec, SPpxe, SPpye_xsec, sizeof(SPpx)/sizeof(SPpx[0]));
	// STAR d+Au pi-
	TDataPoints STAR_dAu_piminus = TDataPoints(TDrawOptions(TNamed("STAR_dAu_piminus", showStarChargedPionsMinusLegend), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullStar, 1.3), "P", "P", TString(showStarChargedPionsMinusLegend)), SPmx, SPmy_xsec, SPmxe, SPmye_xsec, sizeof(SPmx)/sizeof(SPmx[0]));
	TDataPoints STAR_dAu_piminus_unnamed = TDataPoints(TDrawOptions(TNamed("STAR_dAu_piminus", ""), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullStar, 1.3), "P", "P", TString("")), SPmx, SPmy_xsec, SPmxe, SPmye_xsec, sizeof(SPmx)/sizeof(SPmx[0]));
	// STAR p+p pi+
	TDataPoints STAR_pp_piplus = TDataPoints(TDrawOptions(TNamed("STAR_pp_piplus", "STAR #pi^{+}, #pi^{-}"/*showStarChargedPionsPlusPPLegend*/), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenStar, 1.3), "P", "P", TString("STAR #pi^{+}^{}, #pi^{#font[122]{-}}^{}"/*showStarChargedPionsPlusPPLegend*/)), SPPPpx, SPPPpy_xsec, SPPPpxe, SPPPpye_xsec, sizeof(SPPPpx)/sizeof(SPPPpx[0]));
	// STAR p+p pi-
	TDataPoints STAR_pp_piminus = TDataPoints(TDrawOptions(TNamed("STAR_pp_piminus", showStarChargedPionsMinusPPLegend), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullStar, 1.3), "P", "P", TString(showStarChargedPionsMinusPPLegend)), SPPPmx, SPPPmy_xsec, SPPPmxe, SPPPmye_xsec, sizeof(SPPPmx)/sizeof(SPPPmx[0]));
	TDataPoints STAR_pp_piminus_unnamed = TDataPoints(TDrawOptions(TNamed("STAR_pp_piminus", ""), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullStar, 1.3), "P", "P", TString("")), SPPPmx, SPPPmy_xsec, SPPPmxe, SPPPmye_xsec, sizeof(SPPPmx)/sizeof(SPPPmx[0]));
	// STAR p+p (pi+ + pi-)/2
	TDataPoints STAR_pp_pipm = TDataPoints(TDrawOptions(TNamed("STAR_pp_pm", "STAR (#pi^{+}+#pi^{-})/2"/*showStarChargedPionsPlusPPLegend*/), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.3), "P", "P", TString("STAR (#pi^{+}^{ }+ #pi^{#font[122]{-}}^{})/2"/*showStarChargedPionsPlusPPLegend*/)), SPPPpmx, SPPPpmy, SPPPpmxe, SPPPpmye, sizeof(SPPPpmx)/sizeof(SPPPpmx[0]));
	// STAR d+Au (h+ + h-)/2
	TDataPoints STAR_dAu_hpm = TDataPoints(TDrawOptions(TNamed("STAR_dAu_hpm", showStarChargedHadronsLegend), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showStarChargedHadronsLegend)), Sx, Sy_xsec, Sxe, Sye_xsec, sizeof(Sx)/sizeof(Sx[0]));
	// STAR p+p (h+ + h-)/2
	TDataPoints STAR_pp_hpm = TDataPoints(TDrawOptions(TNamed("STAR_pp_hpm", showStarChargedHadronsPPLegend), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showStarChargedHadronsPPLegend)), SPPx, SPPy_xsec, SPPxe, SPPye_xsec, sizeof(SPPx)/sizeof(SPPx[0]));
	// STAR d+Au (h+ + h-)/2 Rcp
	TDataPoints STAR_dAu_Rcp_hpm_1 = TDataPoints(TDrawOptions(TNamed("STAR_Rcp_hpm_1", showStarRcpLegend1), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showStarRcpLegend1)), SRcpx, SRcp1y, SRcpxe, SRcp1ye, sizeof(SRcpx)/sizeof(SRcpx[0]));
	TDataPoints STAR_dAu_Rcp_hpm_2 = TDataPoints(TDrawOptions(TNamed("STAR_Rcp_hpm_2", showStarRcpLegend2), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kMultiply, 1.3), "P", "P", TString(showStarRcpLegend2)), SRcpx, SRcp2y, SRcpxe, SRcp2ye, sizeof(SRcpx)/sizeof(SRcpx[0]));
	// STAR d+Au pi+ Rcp
	TDataPoints STAR_dAu_Rcp_piplus = TDataPoints(TDrawOptions(TNamed("STAR_Rcp_piplus", "STAR #pi^{+}, #pi^{-}"/*showStarChargedPionsPlusRcpLegend*/), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenStar, 1.3), "P", "P", TString("STAR #pi^{+}^{}, #pi^{#font[122]{-}}^{}"/*showStarChargedPionsPlusRcpLegend*/)), SPpRcpx, SPpRcpy, SPpRcpxe, SPpRcpye, sizeof(SPpRcpx)/sizeof(SPpRcpx[0]));
	// STAR d+Au pi- Rcp
	TDataPoints STAR_dAu_Rcp_piminus = TDataPoints(TDrawOptions(TNamed("STAR_Rcp_piminus", showStarChargedPionsMinusRcpLegend), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullStar, 1.3), "P", "P", TString(showStarChargedPionsMinusRcpLegend)), SPmRcpx, SPmRcpy, SPmRcpxe, SPmRcpye, sizeof(SPmRcpx)/sizeof(SPmRcpx[0]));
	TDataPoints STAR_dAu_Rcp_piminus_unnamed = TDataPoints(TDrawOptions(TNamed("STAR_Rcp_piminus", ""/*showStarChargedPionsMinusRcpLegend*/), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kFullStar, 1.3), "P", "P", TString(/*showStarChargedPionsMinusRcpLegend*/)), SPmRcpx, SPmRcpy, SPmRcpxe, SPmRcpye, sizeof(SPmRcpx)/sizeof(SPmRcpx[0]));
	// STAR (h+ + h-)/2 RdA
	TDataPoints STAR_dAu_RdA_hpm = TDataPoints(TDrawOptions(TNamed("STAR_RdA_hpm", showStarRdALegend), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showStarRdALegend)), SRdAx, SRdAy, SRdAxe, SRdAye, sizeof(SRdAx)/sizeof(SRdAx[0]));
	// STAR (pi+ + pi-)/2 RdA
	TDataPoints STAR_dAu_RdA_pipm = TDataPoints(TDrawOptions(TNamed("STAR_RdA_pipm", "STAR (#pi^{+}+ #pi^{-})/2"/*showStarRdApLegend*/), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenCross, 1.3), "P", "P", TString("STAR (#pi^{+}^{}+ #pi^{#font[122]{-}}^{})/2"/*showStarRdApLegend*/)), SRdApx, SRdApy, SRdApxe, SRdApye, sizeof(SRdApx)/sizeof(SRdApx[0]));

	const Char_t *Rinf_titlestr = "Asymptotic #font[12]{R}^{_{ }#infty}^{} = 0.5";
	TDataPoints THEORY_etaToPi0_Rinf = TDataPoints(TDrawOptions(TNamed("THEORY_Rinf" , Rinf_titlestr), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(Rinf_titlestr)), RinfFunc);
	TDataPoints THEORY_etaToPi0_Rinf_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_Rinf" , ""), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), RinfFunc);
	const Char_t *mTscaling_titlestr = "#font[12]{m}_{#font[12]{T}} scaling";
	TDataPoints THEORY_etaToPi0_mTscaling = TDataPoints(TDrawOptions(TNamed("THEORY_mTscaling" , mTscaling_titlestr), TAttLine(kRed, kDashed, 2), TAttFill(kRed, 0), TAttMarker(kRed, kOpenSquare, 1.0), "L", "L", TString(mTscaling_titlestr)), mTScalingFunc);
	TDataPoints THEORY_etaToPi0_mTscaling_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_mTscaling" , ""), TAttLine(kRed, kDashed, 2), TAttFill(kRed, 0), TAttMarker(kRed, kOpenSquare, 1.0), "L", "L", TString()), mTScalingFunc);
	TDataPoints THEORY_Runity = TDataPoints(TDrawOptions(TNamed("THEORY_Runity" , "1"), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), RunityFunc);

	TDataPoints THEORY_pi0_dAu = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD" , "NLO pQCD"/*showpQCDLegend*/), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString("NLO pQCD"/*showpQCDLegend*/)), fNLOpQCD);
	TDataPoints THEORY_pi0_dAu_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD" , ""), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), fNLOpQCD);
	TDataPoints THEORY_pi0_dAu_1 = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD_1" , showpQCDErrorsLegend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showpQCDErrorsLegend)), fNLOpQCD1);
	TDataPoints THEORY_pi0_dAu_1_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD_1" , ""), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), fNLOpQCD1);
	TDataPoints THEORY_pi0_dAu_2 = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD_2" , ""), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), fNLOpQCD2);
	TDataPoints THEORY_pi0_dAu_fill1_nodraw = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD_fill1_nodraw", showpQCDErrorsLegend), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "x", "f", TString(showpQCDErrorsLegend)), x_NLOfill, y1_NLOfill, x_err_NLOfill, y1_err_NLOfill, sizeof(x_NLOfill)/sizeof(x_NLOfill[0]));
	TDataPoints THEORY_pi0_dAu_fill1_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD_fill1_unnamed", ""), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "z3", "", TString()), x_NLOfill, y1_NLOfill, x_err_NLOfill, y1_err_NLOfill, sizeof(x_NLOfill)/sizeof(x_NLOfill[0]));
	TDataPoints THEORY_pi0_dAu_fill2 = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCD_fill2", ""), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "z3", "", TString()), x_NLOfill, y2_NLOfill, x_err_NLOfill, y2_err_NLOfill, sizeof(x_NLOfill)/sizeof(x_NLOfill[0]));

	TDataPoints THEORY_pi0_pp = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP" , "NLO pQCD"/*showpQCDPPLegend*/), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString("NLO pQCD"/*showpQCDPPLegend*/)), fNLOpQCDPP);
	TDataPoints THEORY_pi0_pp_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP" , ""), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), fNLOpQCDPP);
	TDataPoints THEORY_pi0_pp_1 = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP_1" , showpQCDPPErrorsLegend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showpQCDPPErrorsLegend)), fNLOpQCDPP1);
	TDataPoints THEORY_pi0_pp_2 = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP_2" , ""), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), fNLOpQCDPP2);
	TDataPoints THEORY_pi0_pp_fill1_nodraw = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP_fill1_nodraw", showpQCDPPErrorsLegend), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "x", "f", TString(showpQCDPPErrorsLegend)), x_NLOPPfill, y1_NLOPPfill, x_err_NLOPPfill, y1_err_NLOPPfill, sizeof(x_NLOPPfill)/sizeof(x_NLOPPfill[0]));
	TDataPoints THEORY_pi0_pp_fill1_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP_fill1_unnamed", ""), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "z3", "", TString("")), x_NLOPPfill, y1_NLOPPfill, x_err_NLOPPfill, y1_err_NLOPPfill, sizeof(x_NLOPPfill)/sizeof(x_NLOPPfill[0]));
	TDataPoints THEORY_pi0_pp_fill2 = TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP_fill2", ""), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "z3", "", TString()), x_NLOPPfill, y2_NLOPPfill, x_err_NLOPPfill, y2_err_NLOPPfill, sizeof(x_NLOPPfill)/sizeof(x_NLOPPfill[0]));

	TDataPoints THEORY2_pi0_pp = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP" , "NLO pQCD"/*showVogel2Legend*/), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString("NLO pQCD"/*showVogel2Legend*/)), fNLOpQCDPP_2);
	TDataPoints THEORY2_pi0_pp_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP" , ""), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString("")), fNLOpQCDPP_2);
	TDataPoints THEORY2_pi0_pp_1 = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP_1" , showVogel2ErrorsLegend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showVogel2ErrorsLegend)), fNLOpQCDPP1_2);
	TDataPoints THEORY2_pi0_pp_2 = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP_2" , ""), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString()), fNLOpQCDPP2_2);
	TDataPoints THEORY2_pi0_pp_fill1_nodraw = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP_fill1_nodraw", showVogel2ErrorsLegend), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "x", "f", TString(showVogel2ErrorsLegend)), x_NLOPPfill_2, y1_NLOPPfill_2, x_err_NLOPPfill_2, y1_err_NLOPPfill_2, sizeof(x_NLOPPfill_2)/sizeof(x_NLOPPfill_2[0]));
	TDataPoints THEORY2_pi0_pp_fill1_unnamed = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP_fill1_unnamed", ""), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "z3", "", TString("")), x_NLOPPfill_2, y1_NLOPPfill_2, x_err_NLOPPfill_2, y1_err_NLOPPfill_2, sizeof(x_NLOPPfill_2)/sizeof(x_NLOPPfill_2[0]));
	TDataPoints THEORY2_pi0_pp_fill2 = TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP_fill2", ""), TAttLine(17, kSolid, 1), TAttFill(17, 1001), TAttMarker(17, kOpenSquare, 1.1), "z3", "", TString()), x_NLOPPfill_2, y2_NLOPPfill_2, x_err_NLOPPfill_2, y2_err_NLOPPfill_2, sizeof(x_NLOPPfill_2)/sizeof(x_NLOPPfill_2[0]));

	data_points_list dataPoints_dAu_crossSection;
	data_points_list dataPoints_dAu_crossSection_div;
	data_points_list dataPoints_dAu_crossSection_div1;
	data_points_list dataPoints_dAu_crossSection_div1_2;
	data_points_list dataPoints_dAu_invYield;
	data_points_list dataPoints_pp_crossSection;
	data_points_list dataPoints_pp_crossSection_spinpaper;
	data_points_list dataPoints_pp_crossSection_div;
	data_points_list dataPoints_pp_crossSection_div1;
	data_points_list dataPoints_pp_crossSection_div1_2;
	data_points_list dataPoints_pp_invYield;
	data_points_list dataPoints_R;
	data_points_list dataPoints_Rcp_pi0;
	data_points_list dataPoints_RdA_pi0;
	data_points_list dataPoints_Rcp_eta;
	data_points_list dataPoints_RdA_eta;
	data_points_list dataPoints_etatopi;
	data_points_list dataPoints_pp_etatopi;
	data_points_list dataPoints_dAu_etatopi;

	// This analysis, final results
	// pi0 inv. yield in pp and dAu min. bias
	data_points_list dataPoints_pi0_invyield;
	// pi0 inv. yield in dAu centrality bins divided by fit to min. bias
	data_points_list dataPoints_pi0_invyield_centr;
	// pi0 inv. yield in pp divided by fit
	data_points_list dataPoints_pi0_pp_invyield;
	// pi0 inv. yield in dAu min. bias divided by fit
	data_points_list dataPoints_pi0_dAu_invyield;
	// eta inv. yield in pp and dAu min. bias
	data_points_list dataPoints_eta_invyield;
	// eta inv. yield in pp divided by fit
	data_points_list dataPoints_eta_pp_invyield;
	// eta inv. yield in dAu divided by fit
	data_points_list dataPoints_eta_dAu_invyield;
	// pi0 cross section in pp
	data_points_list dataPoints_pi0_pp_crossection;
	data_points_list dataPoints_pi0_pp_crossection_1;
	data_points_list dataPoints_pi0_pp_crossection_div;
	data_points_list dataPoints_pi0_pp_crossection_div_2;
	data_points_list dataPoints_pi0_pp_crossection_div1;
	data_points_list dataPoints_pi0_pp_crossection_div1_2;
	// pi0 cross section in dAu
	data_points_list dataPoints_pi0_dAu_crossection;
	data_points_list dataPoints_pi0_dAu_crossection_div;
	data_points_list dataPoints_pi0_dAu_crossection_div_2;
	data_points_list dataPoints_pi0_dAu_crossection_div1;
	data_points_list dataPoints_pi0_dAu_crossection_div1_2;
	// R_CP and R_dA
	data_points_list dataPoints_pi0_dAu_Rcp;
	data_points_list dataPoints_pi0_dAu_Rcp_2;
	data_points_list dataPoints_pi0_RdA;
	data_points_list dataPoints_pi0_RdA_2;
	data_points_list dataPoints_eta_dAu_Rcp;
	data_points_list dataPoints_eta_dAu_Rcp_2;
	data_points_list dataPoints_eta_RdA;
	data_points_list dataPoints_eta_RdA_2;
	// eta/pi0 ratios
	data_points_list dataPoints_etaToPi0_pp;
	data_points_list dataPoints_etaToPi0_pp1;
	data_points_list dataPoints_etaToPi0_dAu;
	data_points_list dataPoints_etaToPi0_dAu1;

	dataPoints_dAu_crossSection.push_back(THEORY_pi0_dAu);
	//dataPoints_dAu_crossSection_div.push_back(THEORY_pi0_dAu_fill1_unnamed);
	//dataPoints_dAu_crossSection_div.push_back(THEORY_pi0_dAu_fill2);
	dataPoints_dAu_crossSection_div.push_back(THEORY_pi0_dAu_unnamed);
	//dataPoints_dAu_crossSection_div.push_back(THEORY_pi0_dAu_fill1_nodraw);
	dataPoints_dAu_crossSection_div.push_back(THEORY_pi0_dAu_1_unnamed);
	dataPoints_dAu_crossSection_div.push_back(THEORY_pi0_dAu_2);

	if (showPPKKP) {
	    dataPoints_pp_crossSection.push_back(THEORY_pi0_pp_unnamed);
	    dataPoints_pp_crossSection_spinpaper.push_back(THEORY_pi0_pp);
	    //dataPoints_pp_crossSection_div.push_back(THEORY_pi0_pp_fill1_unnamed);
	    //dataPoints_pp_crossSection_div.push_back(THEORY_pi0_pp_fill2);
	    dataPoints_pp_crossSection_div.push_back(THEORY_pi0_pp_unnamed);
	    //dataPoints_pp_crossSection_div.push_back(THEORY_pi0_pp_fill1_nodraw);
	    dataPoints_pp_crossSection_div.push_back(THEORY_pi0_pp_1);
	    dataPoints_pp_crossSection_div.push_back(THEORY_pi0_pp_2);
	    //dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp);
	} else {
	    dataPoints_pp_crossSection.push_back(THEORY2_pi0_pp_unnamed);
	    dataPoints_pp_crossSection_spinpaper.push_back(STAR_pp_pipm);
	    dataPoints_pp_crossSection_spinpaper.push_back(THEORY2_pi0_pp);
	    //dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp_fill1_unnamed);
	    //dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp_fill2);
	    dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp_unnamed);
	    //dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp_fill1_nodraw);
	    dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp_1);
	    dataPoints_pp_crossSection_div.push_back(THEORY2_pi0_pp_2);
	}

	dataPoints_dAu_crossSection_div1.push_back(THEORY_pi0_dAu_unnamed);
	//dataPoints_dAu_crossSection.push_back(STAR_dAu_hpm);
	//dataPoints_dAu_crossSection_div.push_back(STAR_dAu_hpm);
	//dataPoints_dAu_crossSection.push_back(STAR_dAu_piplus);
	dataPoints_dAu_crossSection_div1_2.push_back(STAR_dAu_piplus);
	//dataPoints_dAu_crossSection.push_back(STAR_dAu_piminus);
	dataPoints_dAu_crossSection_div1_2.push_back(STAR_dAu_piminus_unnamed);
	//dataPoints_dAu_crossSection.push_back(PHENIX_dAu_pi0);
	dataPoints_dAu_crossSection_div1_2.push_back(PHENIX_dAu_pi0);

	if (showPPKKP) {
	    dataPoints_pp_crossSection_div1.push_back(THEORY_pi0_pp_unnamed);
	} else {
	    dataPoints_pp_crossSection_div1.push_back(THEORY2_pi0_pp_unnamed);
	}
	//dataPoints_pp_crossSection.push_back(STAR_pp_hpm);
	//dataPoints_pp_crossSection_div.push_back(STAR_pp_hpm);
	//dataPoints_pp_crossSection.push_back(STAR_pp_piplus);
	dataPoints_pp_crossSection_div1_2.push_back(STAR_pp_piplus);
	//dataPoints_pp_crossSection.push_back(STAR_pp_piminus);
	dataPoints_pp_crossSection_div1_2.push_back(STAR_pp_piminus_unnamed);
	//dataPoints_pp_crossSection.push_back(PHENIX_pp_pi0);
	dataPoints_pp_crossSection_div1_2.push_back(PHENIX_pp_pi0);

	//dataPoints_R.push_back(STAR_dAu_Rcp_hpm_1);
	//dataPoints_R.push_back(STAR_dAu_Rcp_hpm_2);
	//dataPoints_R.push_back(STAR_dAu_RdA_hpm);
	///dataPoints_R.push_back(STAR_dAu_RdA_pipm);
	dataPoints_R.push_back(STAR_dAu_Rcp_piplus);
	dataPoints_R.push_back(STAR_dAu_Rcp_piminus);
	//dataPoints_Rcp_pi0.push_back(STAR_dAu_Rcp_hpm_1);
	//dataPoints_Rcp_pi0.push_back(STAR_dAu_Rcp_hpm_2);
	dataPoints_Rcp_pi0.push_back(STAR_dAu_Rcp_piplus);
	dataPoints_Rcp_pi0.push_back(STAR_dAu_Rcp_piminus_unnamed);
	//dataPoints_Rcp_eta.push_back(STAR_dAu_Rcp_hpm_1);
	//dataPoints_Rcp_eta.push_back(STAR_dAu_Rcp_hpm_2);
	///dataPoints_Rcp_eta.push_back(STAR_dAu_Rcp_piplus);
	///dataPoints_Rcp_eta.push_back(STAR_dAu_Rcp_piminus);
	//dataPoints_RdA_pi0.push_back(STAR_dAu_RdA_hpm);
	//dataPoints_RdA_eta.push_back(STAR_dAu_RdA_hpm);
	dataPoints_RdA_pi0.push_back(STAR_dAu_RdA_pipm);
	//dataPoints_RdA_eta.push_back(STAR_dAu_RdA_pipm);

	//dataPoints_R.push_back(PHENIX_pi0_RdA);
	dataPoints_RdA_pi0.push_back(PHENIX_pi0_RdA);

	//dataPoints_R.push_back(PHENIX_eta_RdA);
	dataPoints_RdA_eta.push_back(PHENIX_eta_RdA);

	dataPoints_etatopi.push_back(PHENIX_etaToPi0_pp);
	dataPoints_pp_etatopi.push_back(PHENIX_etaToPi0_pp);

	dataPoints_etatopi.push_back(PHENIX_etaToPi0_dAu);
	dataPoints_dAu_etatopi.push_back(PHENIX_etaToPi0_dAu);

	if (false) {
	    // pi0 from BEMC (Andre Mischke)
	    dataPoints_dAu_invYield.push_back(TDataPoints(TDrawOptions(TNamed("STAR_dAu_pi0_amischke", showMischkePi0Legend), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showMischkePi0Legend)), Ax, Ay, Axe, Aye, sizeof(Ax)/sizeof(Ax[0])));
	}
	if (false) {
	    // pi0 from BEMC (Alex Stolpovsky d+Au)
	    dataPoints_dAu_invYield.push_back(TDataPoints(TDrawOptions(TNamed("STAR_dAu_pi0_alexst_MB" , showStolpovskyPi0Legend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleDown, 1.3), "P", "P", TString(showStolpovskyPi0Legend)), ASxMB, ASyMB, ASxeMB, ASyeMB, sizeof(ASxMB)/sizeof(ASxMB[0])));
	    dataPoints_dAu_invYield.push_back(TDataPoints(TDrawOptions(TNamed("STAR_dAu_pi0_alexst_HT1", showStolpovskyPi0Legend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleDown, 1.3), "P", "P", TString()), ASxHT1, ASyHT1, ASxeHT1, ASyeHT1, sizeof(ASxHT1)/sizeof(ASxHT1[0])));
	    dataPoints_dAu_invYield.push_back(TDataPoints(TDrawOptions(TNamed("STAR_dAu_pi0_alexst_HT2", showStolpovskyPi0Legend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleDown, 1.3), "P", "P", TString()), ASxHT2, ASyHT2, ASxeHT2, ASyeHT2, sizeof(ASxHT2)/sizeof(ASxHT2[0])));
	}
	if (false) {
	    // pi0 from BEMC (Alex Stolpovsky p+p 2004)
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("STAR_pp_pi0_alexst_MB" , showStolpovskyPi0PPLegend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleDown, 1.3), "P", "P", TString(showStolpovskyPi0PPLegend)), ASppxMB, ASppyMB, ASppxeMB, ASppyeMB, sizeof(ASppxMB)/sizeof(ASppxMB[0])));
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("STAR_pp_pi0_alexst_HT1", showStolpovskyPi0PPLegend), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 1001), TAttMarker(kBlack, kFullTriangleDown, 1.3), "P", "P", TString()), ASppxHT1, ASppyHT1, ASppxeHT1, ASppyeHT1, sizeof(ASppxHT1)/sizeof(ASppxHT1[0])));
	}
	if (false) {
	    // Martijn's pi0 in d+Au
	    dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("STAR_dAu_pi0_mrusscher" , showMartijnPi0dAuLegend), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullTriangleDown, 1.3), "P", "P", TString(showMartijnPi0dAuLegend)), MartijnDAux, MartijnDAuy, MartijnDAuxe, MartijnDAuye, sizeof(MartijnDAux)/sizeof(MartijnDAux[0])));
	    dataPoints_dAu_crossSection_div1.push_back(TDataPoints(TDrawOptions(TNamed("STAR_dAu_pi0_mrusscher" , showMartijnPi0dAuLegend), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullTriangleDown, 1.3), "P", "P", TString(showMartijnPi0dAuLegend)), MartijnDAux, MartijnDAuy, MartijnDAuxe, MartijnDAuye, sizeof(MartijnDAux)/sizeof(MartijnDAux[0])));
	}
	if (false) {
	    // Martijn's pi0 in p+p
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("STAR_pp_pi0_mrusscher" , showMartijnPi0PPLegend), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullTriangleDown, 1.3), "P", "P", TString(showMartijnPi0PPLegend)), MartijnPPx, MartijnPPy, MartijnPPxe, MartijnPPye, sizeof(MartijnPPx)/sizeof(MartijnPPx[0])));
	    dataPoints_pp_crossSection_div1.push_back(TDataPoints(TDrawOptions(TNamed("STAR_pp_pi0_mrusscher" , showMartijnPi0PPLegend), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullTriangleDown, 1.3), "P", "P", TString(showMartijnPi0PPLegend)), MartijnPPx, MartijnPPy, MartijnPPxe, MartijnPPye, sizeof(MartijnPPx)/sizeof(MartijnPPx[0])));
	}
	if (false) {
	    // Martijn's pi0 RdA
	    dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("STAR_RdA_pi0_mrusscher" , showMartijnPi0RdALegend), TAttLine(kRed, kSolid, 1), TAttFill(kRed, 1001), TAttMarker(kRed, kFullTriangleDown, 1.3), "P", "P", TString(showMartijnPi0RdALegend)), MartijnRdAx, MartijnRdAy, MartijnRdAxe, MartijnRdAye, sizeof(MartijnRdAx)/sizeof(MartijnRdAx[0])));
	}
	if (false) {
	    // pi0 cross section in p+p 2005 from Frank Simon
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("STAR_pp_pi0_fsimon" , showSimonPi0PPCrossectionLegend), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showSimonPi0PPCrossectionLegend)), FSxs_x, FSxs_y, FSxs_xe, FSxs_ye, sizeof(FSxs_x)/sizeof(FSxs_x[0])));
	    dataPoints_pp_crossSection_div1.push_back(TDataPoints(TDrawOptions(TNamed("STAR_pp_pi0_fsimon" , showSimonPi0PPCrossectionLegend), TAttLine(kGreen, kSolid, 1), TAttFill(kGreen, 1001), TAttMarker(kGreen, kFullTriangleUp, 1.3), "P", "P", TString(showSimonPi0PPCrossectionLegend)), FSxs_x, FSxs_y, FSxs_xe, FSxs_ye, sizeof(FSxs_x)/sizeof(FSxs_x[0])));
	}

	if (false) {
	    dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("WEIGHT_NLOpQCD" , showpQCDWeightLegend), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kOpenSquare, 1.0), "L", "L", TString(showpQCDWeightLegend)), fNLOpQCDWeight));
	}
	if (false) {
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP" , showpQCDPPLegend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showpQCDPPLegend)), fNLOpQCDPP));
	    dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPP" , showpQCDPPLegend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showpQCDPPLegend)), fNLOpQCDPP));
	}
	if (false) {
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPPK" , showpQCDPPKLegend), TAttLine(kBlack, kDotted, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showpQCDPPKLegend)), fNLOpQCDPPK));
	    dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_NLOpQCDPPK" , showpQCDPPKLegend), TAttLine(kBlack, kDotted, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showpQCDPPKLegend)), fNLOpQCDPPK));
	}
	if (false) {
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP" , showVogel2Legend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showVogel2Legend)), fNLOpQCDPP_2));
	    dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY2_NLOpQCDPP" , showVogel2Legend), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenSquare, 1.0), "L", "L", TString(showVogel2Legend)), fNLOpQCDPP_2));
	}
	if (false) {
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("WEIGHT_NLOpQCDPP" , showpQCDPPWeightLegend), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kOpenSquare, 1.0), "L", "L", TString(showpQCDPPWeightLegend)), fNLOpQCDPPWeight));
	}
	dataPoints_etatopi.push_back(THEORY_etaToPi0_Rinf);
	dataPoints_pp_etatopi.push_back(THEORY_etaToPi0_Rinf);
	dataPoints_dAu_etatopi.push_back(THEORY_etaToPi0_Rinf_unnamed);
	dataPoints_etatopi.push_back(THEORY_etaToPi0_mTscaling);
	dataPoints_pp_etatopi.push_back(THEORY_etaToPi0_mTscaling);
	dataPoints_dAu_etatopi.push_back(THEORY_etaToPi0_mTscaling_unnamed);
	dataPoints_R.push_back(THEORY_Runity);
	dataPoints_Rcp_pi0.push_back(THEORY_Runity);
	dataPoints_RdA_pi0.push_back(THEORY_Runity);
	dataPoints_Rcp_eta.push_back(THEORY_Runity);
	dataPoints_RdA_eta.push_back(THEORY_Runity);

	if (false) {
	    // NLO pQCD p+p pi0 mu=0.9pT
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_09", pp_pi0_mu09_title), TAttLine(kBlue, kSolid, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pp_pi0_mu09_title)), pp_pi0_mu09_x, pp_pi0_mu09_y, pp_pi0_mu09_ex, pp_pi0_mu09_ey, sizeof(pp_pi0_mu09_x)/sizeof(pp_pi0_mu09_x[0])));
	    // NLO pQCD p+p pi0 mu=1.2pT
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_12", pp_pi0_mu12_title), TAttLine(kBlue, kDashed, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pp_pi0_mu12_title)), pp_pi0_mu12_x, pp_pi0_mu12_y, pp_pi0_mu12_ex, pp_pi0_mu12_ey, sizeof(pp_pi0_mu12_x)/sizeof(pp_pi0_mu12_x[0])));
	    // NLO pQCD p+p pi0 mu=1.5pT
	    dataPoints_pp_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_15", pp_pi0_mu15_title), TAttLine(kBlue, kDotted, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pp_pi0_mu15_title)), pp_pi0_mu15_x, pp_pi0_mu15_y, pp_pi0_mu15_ex, pp_pi0_mu15_ey, sizeof(pp_pi0_mu15_x)/sizeof(pp_pi0_mu15_x[0])));
	}
	if (false) {
	    // NLO pQCD p+A pi0 mu=1.2pT
	    //dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA", pA_pi0_mu12_title), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_title)), pA_pi0_mu12_x, pA_pi0_mu12_y, pA_pi0_mu12_ex, pA_pi0_mu12_ey, sizeof(pA_pi0_mu12_x)/sizeof(pA_pi0_mu12_x[0])));
	    dataPoints_dAu_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA", pA_pi0_mu12_title), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_title)), pA_pi0_mu12_x, pA_pi0_mu12_y, pA_pi0_mu12_ex, pA_pi0_mu12_ey, sizeof(pA_pi0_mu12_x)/sizeof(pA_pi0_mu12_x[0])));
	    // NLO pQCD p+A pi0 mu=1.2pT + EKS98
	    //dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_EKS98", pA_pi0_mu12_EKS98_title), TAttLine(kBlue, kSolid, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_EKS98_title)), pA_pi0_mu12_EKS98_x, pA_pi0_mu12_EKS98_y, pA_pi0_mu12_EKS98_ex, pA_pi0_mu12_EKS98_ey, sizeof(pA_pi0_mu12_EKS98_x)/sizeof(pA_pi0_mu12_EKS98_x[0])));
	    dataPoints_dAu_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_EKS98", pA_pi0_mu12_EKS98_title), TAttLine(kBlue, kSolid, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_EKS98_title)), pA_pi0_mu12_EKS98_x, pA_pi0_mu12_EKS98_y, pA_pi0_mu12_EKS98_ex, pA_pi0_mu12_EKS98_ey, sizeof(pA_pi0_mu12_EKS98_x)/sizeof(pA_pi0_mu12_EKS98_x[0])));
	    // NLO pQCD p+A pi0 mu=1.2pT + nDS
	    //dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_nDS", pA_pi0_mu12_nDS_title), TAttLine(kBlue, kDashed, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_nDS_title)), pA_pi0_mu12_nDS_x, pA_pi0_mu12_nDS_y, pA_pi0_mu12_nDS_ex, pA_pi0_mu12_nDS_ey, sizeof(pA_pi0_mu12_nDS_x)/sizeof(pA_pi0_mu12_nDS_x[0])));
	    dataPoints_dAu_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_nDS", pA_pi0_mu12_nDS_title), TAttLine(kBlue, kDashed, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_nDS_title)), pA_pi0_mu12_nDS_x, pA_pi0_mu12_nDS_y, pA_pi0_mu12_nDS_ex, pA_pi0_mu12_nDS_ey, sizeof(pA_pi0_mu12_nDS_x)/sizeof(pA_pi0_mu12_nDS_x[0])));
	    // NLO pQCD p+A pi0 mu=1.2pT + nPDF
	    //dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_nPDF", pA_pi0_mu12_nPDF_title), TAttLine(kBlue, kDotted, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_nPDF_title)), pA_pi0_mu12_nPDF_x, pA_pi0_mu12_nPDF_y, pA_pi0_mu12_nPDF_ex, pA_pi0_mu12_nPDF_ey, sizeof(pA_pi0_mu12_nPDF_x)/sizeof(pA_pi0_mu12_nPDF_x[0])));
	    dataPoints_dAu_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_nPDF", pA_pi0_mu12_nPDF_title), TAttLine(kBlue, kDotted, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_nPDF_title)), pA_pi0_mu12_nPDF_x, pA_pi0_mu12_nPDF_y, pA_pi0_mu12_nPDF_ex, pA_pi0_mu12_nPDF_ey, sizeof(pA_pi0_mu12_nPDF_x)/sizeof(pA_pi0_mu12_nPDF_x[0])));
	    // NLO pQCD p+A pi0 mu=1.2pT + HIJING
	    //dataPoints_dAu_crossSection.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_HIJING", pA_pi0_mu12_HIJING_title), TAttLine(kBlue, kDotted, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_HIJING_title)), pA_pi0_mu12_HIJING_x, pA_pi0_mu12_HIJING_y, pA_pi0_mu12_HIJING_ex, pA_pi0_mu12_HIJING_ey, sizeof(pA_pi0_mu12_HIJING_x)/sizeof(pA_pi0_mu12_HIJING_x[0])));
	    dataPoints_dAu_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_pA_HIJING", pA_pi0_mu12_HIJING_title), TAttLine(kBlue, kDotted, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(pA_pi0_mu12_HIJING_title)), pA_pi0_mu12_HIJING_x, pA_pi0_mu12_HIJING_y, pA_pi0_mu12_HIJING_ex, pA_pi0_mu12_HIJING_ey, sizeof(pA_pi0_mu12_HIJING_x)/sizeof(pA_pi0_mu12_HIJING_x[0])));
	}
	if (false) {
	    // NLO pQCD R_pA pi0 mu=1.2pT + EKS98
	    dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_RpA_EKS98", RpA_pi0_EKS98_title), TAttLine(kBlue, kSolid, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(RpA_pi0_EKS98_title)), RpA_pi0_EKS98_x, RpA_pi0_EKS98_y, RpA_pi0_EKS98_ex, RpA_pi0_EKS98_ey, sizeof(RpA_pi0_EKS98_x)/sizeof(RpA_pi0_EKS98_x[0])));
	    // NLO pQCD R_pA pi0 mu=1.2pT + nDS
	    dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_RpA_nDS", RpA_pi0_nDS_title), TAttLine(kBlue, kDashed, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(RpA_pi0_nDS_title)), RpA_pi0_nDS_x, RpA_pi0_nDS_y, RpA_pi0_nDS_ex, RpA_pi0_nDS_ey, sizeof(RpA_pi0_nDS_x)/sizeof(RpA_pi0_nDS_x[0])));
	    // NLO pQCD R_pA pi0 mu=1.2pT + nPDF
	    //dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_RpA_nPDF", RpA_pi0_nPDF_title), TAttLine(kBlue, kDotted, 3), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(RpA_pi0_nPDF_title)), RpA_pi0_nPDF_x, RpA_pi0_nPDF_y, RpA_pi0_nPDF_ex, RpA_pi0_nPDF_ey, sizeof(RpA_pi0_nPDF_x)/sizeof(RpA_pi0_nPDF_x[0])));
	    // NLO pQCD R_pA pi0 mu=1.2pT + HIJING
	    dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("THEORY_RpA_HIJING", RpA_pi0_HIJING_title), TAttLine(kBlue, kDotted, 1), TAttFill(kBlue, 1001), TAttMarker(kBlue, kOpenSquare, 1.1), "L", "L", TString(RpA_pi0_HIJING_title)), RpA_pi0_HIJING_x, RpA_pi0_HIJING_y, RpA_pi0_HIJING_ex, RpA_pi0_HIJING_ey, sizeof(RpA_pi0_HIJING_x)/sizeof(RpA_pi0_HIJING_x[0])));
	}
/*
	dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("PP_CURRENT", "current"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleUp, 2.0), "LP", "LP", TString("SMD required")), x_MBx_pp, y_MBx_pp_current, x_err_MBx_pp, y_err_MBx_pp_current, sizeof(y_MBx_pp_current)/sizeof(y_MBx_pp_current[0])));
	dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("PP_SMD2",    "smd2"),    TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleDown, 2.0), "LP", "LP", TString("2-strip cut")), x_MBx_pp, y_MBx_pp_smd2, x_err_MBx_pp, y_err_MBx_pp_smd2, sizeof(y_MBx_pp_smd2)/sizeof(y_MBx_pp_smd2[0])));
	dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("PP_NOSMD",   "nosmd"),   TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 2.0), "LP", "LP", TString("SMD not required")), x_MBx_pp, y_MBx_pp_nosmd, x_err_MBx_pp, y_err_MBx_pp_nosmd, sizeof(y_MBx_pp_nosmd)/sizeof(y_MBx_pp_nosmd[0])));
	dataPoints_pp_crossSection_div.push_back(TDataPoints(TDrawOptions(TNamed("PP_TOWER",   "tower"),   TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 2.0), "LP", "LP", TString("Towers only")), x_MBx_pp, y_MBx_pp_tower, x_err_MBx_pp, y_err_MBx_pp_tower, sizeof(y_MBx_pp_tower)/sizeof(y_MBx_pp_tower[0])));

	dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("RDA_CURRENT", "current"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleUp, 2.0), "LP", "LP", TString("SMD required")), x_RDAMB_RdA, y_RDAMB_RdA_current, x_err_RDAMB_RdA, y_err_RDAMB_RdA_current, sizeof(y_RDAMB_RdA_current)/sizeof(y_RDAMB_RdA_current[0])));
	dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("RDA_SMD2",    "smd2"),    TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleDown, 2.0), "LP", "LP", TString("2-strip cut")), x_RDAMB_RdA, y_RDAMB_RdA_smd2, x_err_RDAMB_RdA, y_err_RDAMB_RdA_smd2, sizeof(y_RDAMB_RdA_smd2)/sizeof(y_RDAMB_RdA_smd2[0])));
	dataPoints_RdA_pi0.push_back(TDataPoints(TDrawOptions(TNamed("RDA_TOWER",   "tower"),   TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 2.0), "LP", "LP", TString("Towers only")), x_RDAMB_RdA, y_RDAMB_RdA_tower, x_err_RDAMB_RdA, y_err_RDAMB_RdA_tower, sizeof(y_RDAMB_RdA_tower)/sizeof(y_RDAMB_RdA_tower[0])));

	dataPoints_Rcp_pi0.push_back(TDataPoints(TDrawOptions(TNamed("RCP_CURRENT", "current"), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleUp, 2.0), "LP", "LP", TString("SMD required")), x_RcpMB_dAuRcp, y_RcpMB_dAuRcp_current, x_err_RcpMB_dAuRcp, y_err_RcpMB_dAuRcp_current, sizeof(y_RcpMB_dAuRcp_current)/sizeof(y_RcpMB_dAuRcp_current[0])));
	dataPoints_Rcp_pi0.push_back(TDataPoints(TDrawOptions(TNamed("RCP_SMD2",    "smd2"),    TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleDown, 2.0), "LP", "LP", TString("2-strip cut")), x_RcpMB_dAuRcp, y_RcpMB_dAuRcp_smd2, x_err_RcpMB_dAuRcp, y_err_RcpMB_dAuRcp_smd2, sizeof(y_RcpMB_dAuRcp_smd2)/sizeof(y_RcpMB_dAuRcp_smd2[0])));
	dataPoints_Rcp_pi0.push_back(TDataPoints(TDrawOptions(TNamed("RCP_TOWER",   "tower"),   TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 2.0), "LP", "LP", TString("Towers only")), x_RcpMB_dAuRcp, y_RcpMB_dAuRcp_tower, x_err_RcpMB_dAuRcp, y_err_RcpMB_dAuRcp_tower, sizeof(y_RcpMB_dAuRcp_tower)/sizeof(y_RcpMB_dAuRcp_tower[0])));
*/
	if (true) {

// This analysis: final results
#include "macros/my_results/data_arrays_dAu.C"
#include "macros/my_results/data_arrays_dAuMostcentral.C"
#include "macros/my_results/data_arrays_dAuMostperipheral.C"
#include "macros/my_results/data_arrays_dAuRcp.C"
#include "macros/my_results/data_arrays_pp.C"
#include "macros/my_results/data_arrays_RdA.C"

#define SET_X_ERR(NAME, EX) \
    for (UInt_t i = 0;i < (sizeof((x_err_##NAME))/sizeof((x_err_##NAME)[0]));i++) { \
	(x_err_##NAME)[i] = (EX); \
    }
    
#define SET_X_ERRORS(NAME_MB, NAME_HT1, NAME_HT2) \
    SET_X_ERR(NAME_MB, 0.2) \
    SET_X_ERR(NAME_HT1, 0.2) \
    SET_X_ERR(NAME_HT2, 0.2)

#define ADD_ERR(X, Y, EY, ADD_EY, ADD_EY_PT) \
    for (UInt_t i = 0;i < (sizeof((X))/sizeof((X)[0]));i++) { \
	Float_t e = ((ADD_EY) + ((ADD_EY_PT)*((X)[i]))) * ((Y)[i]); \
	if (e < 0) e = 0; \
	(EY)[i] = TMath::Sqrt((((EY)[i])*((EY)[i])) + (e*e)); \
    }

#define ADD_ERR_ARRAY(X, Y, EY, ADD_EY) \
    for (UInt_t i = 0;i < (sizeof((X))/sizeof((X)[0]));i++) { \
	Float_t e = (i < (sizeof((ADD_EY))/sizeof((ADD_EY)[0]))) ? (((ADD_EY)[i]) * ((Y)[i])) : 0; \
	if (e < 0) e = 0; \
	(EY)[i] = TMath::Sqrt((((EY)[i])*((EY)[i])) + (e*e)); \
    }

#define ADD_ERR_POINTS(NAME_MB, NAME_HT1, NAME_HT2, ADD_EY, ADD_EY_PT) \
    ADD_ERR(x_##NAME_MB,  y_##NAME_MB,  y_err_##NAME_MB,  ADD_EY, ADD_EY_PT) \
    ADD_ERR(x_##NAME_HT1, y_##NAME_HT1, y_err_##NAME_HT1, ADD_EY, ADD_EY_PT) \
    ADD_ERR(x_##NAME_HT2, y_##NAME_HT2, y_err_##NAME_HT2, ADD_EY, ADD_EY_PT)

#define ADD_ERR_POINTS_ARRAY(NAME_MB, NAME_HT1, NAME_HT2, ADD_EY_MB, ADD_EY_HT1, ADD_EY_HT2) \
    ADD_ERR_ARRAY(x_##NAME_MB,  y_##NAME_MB,  y_err_##NAME_MB,  ADD_EY_MB) \
    ADD_ERR_ARRAY(x_##NAME_HT1, y_##NAME_HT1, y_err_##NAME_HT1, ADD_EY_HT1) \
    ADD_ERR_ARRAY(x_##NAME_HT2, y_##NAME_HT2, y_err_##NAME_HT2, ADD_EY_HT2)

#define S2 1.414

#define COPY_ARR(A, PREFIX) \
    Float_t A##_##PREFIX[sizeof(A)/sizeof((A)[0])]; \
    for (UInt_t i = 0;i < (sizeof(A)/sizeof((A)[0]));i++) { \
	A##_##PREFIX[i] = A[i]; \
    }

#define COPY_ARRAYS(AMB, AHT1, AHT2, PREFIX) \
    COPY_ARR(x_##AMB, PREFIX); COPY_ARR(y_##AMB, PREFIX); COPY_ARR(x_err_##AMB, PREFIX); COPY_ARR(y_err_##AMB, PREFIX) \
    COPY_ARR(x_##AHT1, PREFIX); COPY_ARR(y_##AHT1, PREFIX); COPY_ARR(x_err_##AHT1, PREFIX); COPY_ARR(y_err_##AHT1, PREFIX) \
    COPY_ARR(x_##AHT2, PREFIX); COPY_ARR(y_##AHT2, PREFIX); COPY_ARR(x_err_##AHT2, PREFIX); COPY_ARR(y_err_##AHT2, PREFIX)

#define CHECK_ERR(X, Y, EY, MAX_EY) \
    for (UInt_t i = 0;i < (sizeof((X))/sizeof((X)[0]));i++) { \
	if (((EY)[i]/(Y)[i]) > (MAX_EY)) (EY)[i] = (MAX_EY)*(Y)[i]; \
    }

#define ADD_ERRORS_POINTS(NAMEMB, NAMEHT1, NAMEHT2, PREFIX_ADD) \
    ADD_ERROR_POINTS(y_err_##NAMEMB, y_err_##NAMEMB##_##PREFIX_ADD) \
    ADD_ERROR_POINTS(y_err_##NAMEHT1, y_err_##NAMEHT1##_##PREFIX_ADD) \
    ADD_ERROR_POINTS(y_err_##NAMEHT2, y_err_##NAMEHT2##_##PREFIX_ADD)

#define ADD_ERRORS_POINTS1(NAMEMB, NAMEHT1, NAMEHT2, PREFIX, PREFIX_ADD) \
    ADD_ERROR_POINTS(y_err_##NAMEMB##_##PREFIX, y_err_##NAMEMB##_##PREFIX_ADD) \
    ADD_ERROR_POINTS(y_err_##NAMEHT1##_##PREFIX, y_err_##NAMEHT1##_##PREFIX_ADD) \
    ADD_ERROR_POINTS(y_err_##NAMEHT2##_##PREFIX, y_err_##NAMEHT2##_##PREFIX_ADD)

#define FORMAT_XSEC "%4.3e"
#define FORMAT_R "%4.3f"
#define FORMAT_X "%.2f"

/*#define OUTSTR(COUT, N, FORMAT, SUP) \
    { \
	TString s = Form(FORMAT, N); \
	s.ReplaceAll("e+0", "e"); \
	s.ReplaceAll("e-0", "e-"); \
	Bool_t foundE = false; \
	for (Int_t i = -100;i < +100;i++) { \
	    TString is = Form("%i", i); \
	    TString ss1 = "e"; ss1 += is; \
	    TString ss2 = "&\\ee{"; ss2 += is; ss2 += "}"; \
	    foundE |= s.Contains(ss1); \
	    s.ReplaceAll(ss1, ss2); \
	} \
	if (SUP && false) { \
	    TString ss1 = "&\\ee{0}"; \
	    TString ss2 = ""; \
	    foundE |= s.Contains(ss1); \
	    s.ReplaceAll(ss1, ss2); \
	} \
	if (!foundE && false) { \
	    s.ReplaceAll(".", "&"); \
	} \
	COUT << s.Data(); \
    }*/

#define OUTSTR(COUT, N, FORMAT, SUP) \
    { \
	TString s = Form(FORMAT, N); \
	s.ReplaceAll("e+0", "e"); \
	s.ReplaceAll("e-0", "e-"); \
	Bool_t foundE = false; \
	for (Int_t i = -100;i < +100;i++) { \
	    TString is = Form("%i", i); \
	    TString ss1 = "e"; ss1 += is; \
	    TString ss2 = "&times;10<SUP>"; ss2 += is; ss2 += "</SUP>"; \
	    foundE |= s.Contains(ss1); \
	    s.ReplaceAll(ss1, ss2); \
	} \
	if (SUP && false) { \
	    TString ss1 = "&times;10<SUP>0</SUP>"; \
	    TString ss2 = ""; \
	    foundE |= s.Contains(ss1); \
	    s.ReplaceAll(ss1, ss2); \
	} \
	COUT << s.Data(); \
    }

#define OUTPUT_TABLE_5(COUT, NAME, PTLOW, PTHIGH, S1, S2, S3, S4, S5, FORMAT, SUP) \
    for (UInt_t i = 0;i < sizeof(x_##NAME)/sizeof(x_##NAME[0]);i++) { \
	if ((x_##NAME[i] >= PTLOW) && (x_##NAME[i] < PTHIGH)) { \
	    OUTSTR(COUT, x_##NAME[i], FORMAT_X, true) \
	    COUT << "    & "; \
	    OUTSTR(COUT, y_##NAME[i], FORMAT, SUP) \
	    COUT << "    & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S1[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S2[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S3[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S4[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S5[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " \\\\ "; \
	    COUT << endl; \
	} \
    }

#define OUTPUT_TABLE_4(COUT, NAME, PTLOW, PTHIGH, S1, S2, S3, S4, FORMAT, SUP) \
    for (UInt_t i = 0;i < sizeof(x_##NAME)/sizeof(x_##NAME[0]);i++) { \
	if ((x_##NAME[i] >= PTLOW) && (x_##NAME[i] < PTHIGH)) { \
	    OUTSTR(COUT, x_##NAME[i], FORMAT_X, true) \
	    COUT << "    & "; \
	    OUTSTR(COUT, y_##NAME[i], FORMAT, SUP) \
	    COUT << "    & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S1[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S2[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S3[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S4[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " \\\\ "; \
	    COUT << endl; \
	} \
    }

#define OUTPUT_TABLE_3(COUT, NAME, PTLOW, PTHIGH, S1, S2, S3, FORMAT, SUP) \
    for (UInt_t i = 0;i < sizeof(x_##NAME)/sizeof(x_##NAME[0]);i++) { \
	if ((x_##NAME[i] >= PTLOW) && (x_##NAME[i] < PTHIGH)) { \
	    OUTSTR(COUT, x_##NAME[i], FORMAT_X, true) \
	    COUT << "    & "; \
	    OUTSTR(COUT, y_##NAME[i], FORMAT, SUP) \
	    COUT << "    & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S1[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S2[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " & "; \
	    OUTSTR(COUT, y_err_##NAME##_##S3[i], FORMAT, SUP) \
	    COUT << "    "; \
	    COUT << " \\\\ "; \
	    COUT << endl; \
	} \
    }

#define OUTPUT_TABLE_2(COUT, NAME, PTLOW, PTHIGH, S1, S2, FORMAT, SUP) \
    for (UInt_t i = 0;i < sizeof(x_##NAME)/sizeof(x_##NAME[0]);i++) { \
	if ((x_##NAME[i] >= PTLOW) && (x_##NAME[i] < PTHIGH)) { \
	    COUT << "<TR><TD>"; \
	    OUTSTR(COUT, x_##NAME[i], FORMAT_X, true) \
	    COUT << "</TD><TD>"; \
	    OUTSTR(COUT, y_##NAME[i], FORMAT, SUP) \
	    COUT << "</TD><TD>"; \
	    OUTSTR(COUT, y_err_##NAME##_##S1[i], FORMAT, SUP) \
	    COUT << "</TD><TD>"; \
	    OUTSTR(COUT, y_err_##NAME##_##S2[i], FORMAT, SUP) \
	    COUT << "</TD></TR>"; \
	    COUT << endl; \
	} \
    }

CHECK_ERR(x_EtaMBspectrum_pp, y_EtaMBspectrum_pp, y_err_EtaMBspectrum_pp, 0.9)
CHECK_ERR(x_EtaHT1spectrum_pp, y_EtaHT1spectrum_pp, y_err_EtaHT1spectrum_pp, 0.9)
CHECK_ERR(x_EtaHT2spectrum_pp, y_EtaHT2spectrum_pp, y_err_EtaHT2spectrum_pp, 0.9)

CHECK_ERR(x_EtaMBspectrum_dAu, y_EtaMBspectrum_dAu, y_err_EtaMBspectrum_dAu, 0.9)
CHECK_ERR(x_EtaHT1spectrum_dAu, y_EtaHT1spectrum_dAu, y_err_EtaHT1spectrum_dAu, 0.9)
CHECK_ERR(x_EtaHT2spectrum_dAu, y_EtaHT2spectrum_dAu, y_err_EtaHT2spectrum_dAu, 0.9)

// Systematic errors on pi0 cross section in p+p and d+Au
COPY_ARR(x_err_MBcrossection_pp, syst) COPY_ARR(x_err_HT1crossection_pp, syst) COPY_ARR(x_err_HT2crossection_pp, syst)
SET_X_ERRORS(MBcrossection_pp_syst, HT1crossection_pp_syst, HT2crossection_pp_syst)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, stat)
ZERO_ARRAY(y_err_MBcrossection_pp) ZERO_ARRAY(y_err_HT1crossection_pp) ZERO_ARRAY(y_err_HT2crossection_pp)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systA)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systB)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systE)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systBE)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systC)
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systN)
COPY_ARRAYS(MBspectrum_pp, HT1spectrum_pp, HT2spectrum_pp, smdsyst)
ZERO_ARRAY(y_err_MBspectrum_pp_smdsyst) ZERO_ARRAY(y_err_HT1spectrum_pp_smdsyst) ZERO_ARRAY(y_err_HT2spectrum_pp_smdsyst)
COPY_ARR(x_err_MBcrossection_dAu, syst) COPY_ARR(x_err_HT1crossection_dAu, syst) COPY_ARR(x_err_HT2crossection_dAu, syst)
SET_X_ERRORS(MBcrossection_dAu_syst, HT1crossection_dAu_syst, HT2crossection_dAu_syst)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, stat)
ZERO_ARRAY(y_err_MBcrossection_dAu) ZERO_ARRAY(y_err_HT1crossection_dAu) ZERO_ARRAY(y_err_HT2crossection_dAu)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systA)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systB)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systE)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systBE)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systC)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systN)
COPY_ARRAYS(MBspectrum_dAu, HT1spectrum_dAu, HT2spectrum_dAu, smdsyst)
ZERO_ARRAY(y_err_MBspectrum_dAu_smdsyst) ZERO_ARRAY(y_err_HT1spectrum_dAu_smdsyst) ZERO_ARRAY(y_err_HT2spectrum_dAu_smdsyst)
// Combinatorial background normalization: 0.5-3% for pi0
ADD_ERR_POINTS(MBcrossection_pp_systA, HT1crossection_pp_systA, HT2crossection_pp_systA, 0.005, (0.03 - 0.005)/15)
ADD_ERR_POINTS(MBcrossection_dAu_systA, HT1crossection_dAu_systA, HT2crossection_dAu_systA, 0.005, (0.03 - 0.005)/15)
// Mixed-event background composition: 5% for pi0
ADD_ERR_POINTS(MBcrossection_pp_systC, HT1crossection_pp_systC, HT2crossection_pp_systC, 0.05, 0)
ADD_ERR_POINTS(MBcrossection_dAu_systC, HT1crossection_dAu_systC, HT2crossection_dAu_systC, 0.05, 0)
// Low mass background: 1-15% for pi0
ADD_ERR_POINTS(MBcrossection_pp_systC, HT1crossection_pp_systC, HT2crossection_pp_systC, 0*0.01, 0*(0.15-0.01)/15)
ADD_ERR_POINTS(MBcrossection_dAu_systC, HT1crossection_dAu_systC, HT2crossection_dAu_systC, 0*0.01, 0*(0.15-0.01)/15)
// CPV correction: 2%
ADD_ERR_POINTS(MBcrossection_pp_systN, HT1crossection_pp_systN, HT2crossection_pp_systN, 0.02, 0)
ADD_ERR_POINTS(MBcrossection_dAu_systN, HT1crossection_dAu_systN, HT2crossection_dAu_systN, 0.02, 0)
// HighTower scale factors: 3-5% in HighTower
ADD_ERR(x_HT1crossection_pp_systB, y_HT1crossection_pp_systB, y_err_HT1crossection_pp_systB, 0.03, 0)
ADD_ERR(x_HT2crossection_pp_systB, y_HT2crossection_pp_systB, y_err_HT2crossection_pp_systB, 0.05, 0)
ADD_ERR(x_HT1crossection_dAu_systB, y_HT1crossection_dAu_systB, y_err_HT1crossection_dAu_systB, 0.03, 0)
ADD_ERR(x_HT2crossection_dAu_systB, y_HT2crossection_dAu_systB, y_err_HT2crossection_dAu_systB, 0.05, 0)
// Analysis cuts: 5%
ADD_ERR_POINTS(MBcrossection_pp_systA, HT1crossection_pp_systA, HT2crossection_pp_systA, 0.05, 0)
ADD_ERR_POINTS(MBcrossection_dAu_systA, HT1crossection_dAu_systA, HT2crossection_dAu_systA, 0.05, 0)
// Missing material in simulation correction
ADD_ERR_POINTS(MBcrossection_pp_systB, HT1crossection_pp_systB, HT2crossection_pp_systB, 0.02/1.10, 0)
ADD_ERR_POINTS(MBcrossection_dAu_systB, HT1crossection_dAu_systB, HT2crossection_dAu_systB, 0.03/1.15, 0)
// Detector simulation: 80-0%
//ADD_ERR_POINTS(MBcrossection_pp_systC, HT1crossection_pp_systC, HT2crossection_pp_systC, 0.80, -0.80/6)
ADD_ERR(x_HT1crossection_pp_systC, y_HT1crossection_pp_systC, y_err_HT1crossection_pp_systC, 0.40, -0.40/7)
ADD_ERR(x_MBspectrum_pp_smdsyst,  y_MBspectrum_pp_smdsyst,  y_err_MBspectrum_pp_smdsyst, 0.80, -0.80/6)
//ADD_ERR_POINTS(MBcrossection_dAu_systC, HT1crossection_dAu_systC, HT2crossection_dAu_systC, 0.80, -0.80/6)
ADD_ERR(x_HT1crossection_dAu_systC, y_HT1crossection_dAu_systC, y_err_HT1crossection_dAu_systC, 0.40, -0.40/7)
ADD_ERR(x_MBspectrum_dAu_smdsyst,  y_MBspectrum_dAu_smdsyst,  y_err_MBspectrum_dAu_smdsyst, 0.80, -0.80/6)
// SMD scale uncertainty
ADD_ERR(x_HT1crossection_pp_systB, y_HT1crossection_pp_systB, y_err_HT1crossection_pp_systB, 0.04, -0.0025)
ADD_ERR(x_HT2crossection_pp_systB, y_HT2crossection_pp_systB, y_err_HT2crossection_pp_systB, 0.04, -0.0025)
ADD_ERR(x_HT1crossection_dAu_systB, y_HT1crossection_dAu_systB, y_err_HT1crossection_dAu_systB, 0.04, -0.0025)
ADD_ERR(x_HT2crossection_dAu_systB, y_HT2crossection_dAu_systB, y_err_HT2crossection_dAu_systB, 0.04, -0.0025)
// Bin centering: 1.5%
ADD_ERR_POINTS(MBcrossection_pp_systC, HT1crossection_pp_systC, HT2crossection_pp_systC, 0.015, 0)
ADD_ERR_POINTS(MBcrossection_dAu_systC, HT1crossection_dAu_systC, HT2crossection_dAu_systC, 0.015, 0)
// Vertex finding efficiency: 1.08% in d+Au min. bias
ADD_ERR_POINTS(MBcrossection_dAu_systN, HT1crossection_dAu_systN, HT2crossection_dAu_systN, 0.0108, 0)
// BBC cross section: 6.94% in p+p
ADD_ERR_POINTS(MBcrossection_pp_systN, HT1crossection_pp_systN, HT2crossection_pp_systN, 0.0694, 0)
// BBC efficiency: 9.20% in p+p
ADD_ERR_POINTS(MBcrossection_pp_systN, HT1crossection_pp_systN, HT2crossection_pp_systN, 0.0920, 0)
// MinBias cross section: 4.07% in d+Au
ADD_ERR_POINTS(MBcrossection_dAu_systN, HT1crossection_dAu_systN, HT2crossection_dAu_systN, 0.0407, 0)
// MinBias efficiency: 3.16% in d+Au
ADD_ERR_POINTS(MBcrossection_dAu_systN, HT1crossection_dAu_systN, HT2crossection_dAu_systN, 0.0316, 0)
// Energy scale: 37-39% in p+p, 40-45% in d+Au
//ADD_ERR_POINTS(MBcrossection_pp_systE, HT1crossection_pp_systE, HT2crossection_pp_systE, 0.37, 0.0014)
//ADD_ERR_POINTS(MBcrossection_dAu_systE, HT1crossection_dAu_systE, HT2crossection_dAu_systE, 0.36, 0.0073)
// Energy scale: arrays
Float_t escale_err_pp_MB[] = {0.125, 0.185, 0.231, 0.271, 0.305, 0.301};
Float_t escale_err_pp_HT1[] = {0.159, 0.05, 0.103, 0.161, 0.221, 0.30};
Float_t escale_err_pp_HT2[] = {0.095, 0.144, 0.200, 0.264, 0.320, 0.378};
ADD_ERR_POINTS_ARRAY(MBcrossection_pp_systE, HT1crossection_pp_systE, HT2crossection_pp_systE, escale_err_pp_MB, escale_err_pp_HT1, escale_err_pp_HT2)
Float_t escale_err_dAu_MB[] = {0.095, 0.155, 0.209, 0.264, 0.356, 0.529};
Float_t escale_err_dAu_HT1[] = {0.244, 0.05, 0.089, 0.134, 0.186, 0.283};
Float_t escale_err_dAu_HT2[] = {0.149, 0.05, 0.114, 0.191, 0.278, 0.377};
ADD_ERR_POINTS_ARRAY(MBcrossection_dAu_systE, HT1crossection_dAu_systE, HT2crossection_dAu_systE, escale_err_dAu_MB, escale_err_dAu_HT1, escale_err_dAu_HT2)
// Low mass bg normalization: arrays
//Float_t lowmassbg_err_pp_MB[] = {0.005, 0.048, 0.303, 0.567, 0.797, 0.403};
//Float_t lowmassbg_err_dAu_MB[] = {0.005, 0.048, 0.303, 0.567, 0.797, 0.403};
//ADD_ERR_ARRAY(x_MBcrossection_pp, y_MBcrossection_pp, y_err_MBcrossection_pp_systB, lowmassbg_err_pp_MB)
//ADD_ERR_ARRAY(x_MBcrossection_dAu, y_MBcrossection_dAu, y_err_MBcrossection_dAu_systB, lowmassbg_err_dAu_MB)
// Combine B+E->BE
ADD_ERRORS_POINTS1(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systBE, systB)
ADD_ERRORS_POINTS1(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, systBE, systE)
ADD_ERRORS_POINTS1(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systBE, systB)
ADD_ERRORS_POINTS1(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, systBE, systE)
// Norm. errors
Float_t x_MBcrossection_pp_norm[] = {16.75}; Float_t x_err_MBcrossection_pp_norm[] = {0.25}; Float_t y_MBcrossection_pp_norm[1]; y_MBcrossection_pp_norm[0] = fDIVPP ? fDIVPP->Eval(x_MBcrossection_pp_norm[0]) : 1.0; Float_t y_err_MBcrossection_pp_norm[1]; y_err_MBcrossection_pp_norm[0] = y_MBcrossection_pp_norm[0] * y_err_MBcrossection_pp_systN[0] / y_MBcrossection_pp_systN[0];
Float_t x_MBcrossection_pp_norm1[] = {19.75}; Float_t x_err_MBcrossection_pp_norm1[] = {0.25}; Float_t y_MBcrossection_pp_norm1[1]; y_MBcrossection_pp_norm1[0] = fDIVPP ? fDIVPP->Eval(x_MBcrossection_pp_norm1[0]) : 1.0; Float_t y_err_MBcrossection_pp_norm1[1]; y_err_MBcrossection_pp_norm1[0] = y_MBcrossection_pp_norm1[0] * y_err_MBcrossection_pp_systN[0] / y_MBcrossection_pp_systN[0];
Float_t x_MBcrossection_dAu_norm[] = {16.75}; Float_t x_err_MBcrossection_dAu_norm[] = {0.25}; Float_t y_MBcrossection_dAu_norm[1]; y_MBcrossection_dAu_norm[0] = fDIV ? fDIV->Eval(x_MBcrossection_dAu_norm[0]) : 1.0; Float_t y_err_MBcrossection_dAu_norm[1]; y_err_MBcrossection_dAu_norm[0] = y_MBcrossection_dAu_norm[0] * y_err_MBcrossection_dAu_systN[0] / y_MBcrossection_dAu_systN[0];
Float_t x_MBcrossection_dAu_norm1[] = {19.75}; Float_t x_err_MBcrossection_dAu_norm1[] = {0.25}; Float_t y_MBcrossection_dAu_norm1[1]; y_MBcrossection_dAu_norm1[0] = fDIV ? fDIV->Eval(x_MBcrossection_dAu_norm1[0]) : 1.0; Float_t y_err_MBcrossection_dAu_norm1[1]; y_err_MBcrossection_dAu_norm1[0] = y_MBcrossection_dAu_norm1[0] * y_err_MBcrossection_dAu_systN[0] / y_MBcrossection_dAu_systN[0];
// Add errors for display
COPY_ARRAYS(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, esyst)
ADD_ERRORS_POINTS1(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, stat, systA)
ADD_ERRORS_POINTS1(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, esyst, systB)
ADD_ERRORS_POINTS1(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, esyst, systC)
ADD_ERRORS_POINTS1(MBcrossection_pp, HT1crossection_pp, HT2crossection_pp, esyst, systE)
COPY_ARRAYS(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, esyst)
ADD_ERRORS_POINTS1(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, stat, systA)
ADD_ERRORS_POINTS1(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, esyst, systB)
ADD_ERRORS_POINTS1(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, esyst, systC)
ADD_ERRORS_POINTS1(MBcrossection_dAu, HT1crossection_dAu, HT2crossection_dAu, esyst, systE)
// Write tables
//{
//ofstream ofstr("table_crossection_pi0_pp.tex");
//OUTPUT_TABLE_4(ofstr, MBcrossection_pp,  pp_MB,      pp_MB_HT1,  stat, systA, systBE, systC, FORMAT_XSEC, true)
//OUTPUT_TABLE_4(ofstr, HT1crossection_pp, pp_MB_HT1,  pp_HT1_HT2, stat, systA, systBE, systC, FORMAT_XSEC, true)
//OUTPUT_TABLE_4(ofstr, HT2crossection_pp, pp_HT1_HT2, pp_HT2,     stat, systA, systBE, systC, FORMAT_XSEC, true)
//}
{
ofstream ofstr("table_crossection_pi0_pp.html");
OUTPUT_TABLE_2(ofstr, MBcrossection_pp,  pp_MB,      pp_MB_HT1,  stat, esyst, FORMAT_XSEC, true)
OUTPUT_TABLE_2(ofstr, HT1crossection_pp, pp_MB_HT1,  pp_HT1_HT2, stat, esyst, FORMAT_XSEC, true)
OUTPUT_TABLE_2(ofstr, HT2crossection_pp, pp_HT1_HT2, pp_HT2,     stat, esyst, FORMAT_XSEC, true)
}
//{
//ofstream ofstr("table_crossection_pi0_dAu.tex");
//OUTPUT_TABLE_4(ofstr, MBcrossection_dAu,  dAu_MB,      dAu_MB_HT1,  stat, systA, systBE, systC, FORMAT_XSEC, true)
//OUTPUT_TABLE_4(ofstr, HT1crossection_dAu, dAu_MB_HT1,  dAu_HT1_HT2, stat, systA, systBE, systC, FORMAT_XSEC, true)
//OUTPUT_TABLE_4(ofstr, HT2crossection_dAu, dAu_HT1_HT2, dAu_HT2,     stat, systA, systBE, systC, FORMAT_XSEC, true)
//}

// Systematic errors on eta/pi0 ratio in p+p and d+Au
COPY_ARR(x_err_EtaToPi0MB_pp, syst) COPY_ARR(x_err_EtaToPi0HT1_pp, syst) COPY_ARR(x_err_EtaToPi0HT2_pp, syst)
SET_X_ERRORS(EtaToPi0MB_pp_syst, EtaToPi0HT1_pp_syst, EtaToPi0HT2_pp_syst)
COPY_ARRAYS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, stat)
ZERO_ARRAY(y_err_EtaToPi0MB_pp) ZERO_ARRAY(y_err_EtaToPi0HT1_pp) ZERO_ARRAY(y_err_EtaToPi0HT2_pp)
COPY_ARRAYS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, systA)
COPY_ARRAYS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, systC)
COPY_ARR(x_err_EtaToPi0MB_dAu, syst) COPY_ARR(x_err_EtaToPi0HT1_dAu, syst) COPY_ARR(x_err_EtaToPi0HT2_dAu, syst)
SET_X_ERRORS(EtaToPi0MB_dAu_syst, EtaToPi0HT1_dAu_syst, EtaToPi0HT2_dAu_syst)
COPY_ARRAYS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, stat)
ZERO_ARRAY(y_err_EtaToPi0MB_dAu) ZERO_ARRAY(y_err_EtaToPi0HT1_dAu) ZERO_ARRAY(y_err_EtaToPi0HT2_dAu)
COPY_ARRAYS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, systA)
COPY_ARRAYS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, systC)
// Combinatorial background normalization: (10-0.5)-(50-3)% (eta-pi0)
ADD_ERR_POINTS(EtaToPi0MB_pp_systA, EtaToPi0HT1_pp_systA, EtaToPi0HT2_pp_systA, 0.10-0.005, ((0.50-0.03) - (0.10-0.005))/15)
ADD_ERR_POINTS(EtaToPi0MB_dAu_systA, EtaToPi0HT1_dAu_systA, EtaToPi0HT2_dAu_systA, 0.10-0.005, ((0.50-0.03) - (0.10-0.005))/15)
// Mixed-event background composition: (5-1)% (eta-pi0)
ADD_ERR_POINTS(EtaToPi0MB_pp_systC, EtaToPi0HT1_pp_systC, EtaToPi0HT2_pp_systC, (0.05-0.01), 0)
ADD_ERR_POINTS(EtaToPi0MB_dAu_systC, EtaToPi0HT1_dAu_systC, EtaToPi0HT2_dAu_systC, (0.05-0.01), 0)
// Low mass background: 1-15% (pi0)
ADD_ERR_POINTS(EtaToPi0MB_pp_systC, EtaToPi0HT1_pp_systC, EtaToPi0HT2_pp_systC, 0*0.01, 0*(0.15-0.01)/15)
ADD_ERR_POINTS(EtaToPi0MB_dAu_systC, EtaToPi0HT1_dAu_systC, EtaToPi0HT2_dAu_systC, 0*0.01, 0*(0.15-0.01)/15)
// Analysis cuts: 5% for pi0 and 5% for eta
ADD_ERR_POINTS(EtaToPi0MB_pp_systA, EtaToPi0HT1_pp_systA, EtaToPi0HT2_pp_systA, S2*0.05, 0)
ADD_ERR_POINTS(EtaToPi0MB_dAu_systA, EtaToPi0HT1_dAu_systA, EtaToPi0HT2_dAu_systA, S2*0.05, 0)
// Write tables
{
ofstream ofstr("table_etatopi0_pp.tex");
OUTPUT_TABLE_3(ofstr, EtaToPi0MB_pp,  pp_etatopi0_MB,      pp_etatopi0_MB_HT1,  stat, systA, systC, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, EtaToPi0HT1_pp, pp_etatopi0_MB_HT1,  pp_etatopi0_HT1_HT2, stat, systA, systC, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, EtaToPi0HT2_pp, pp_etatopi0_HT1_HT2, pp_etatopi0_HT2,     stat, systA, systC, FORMAT_R, true)
}
{
ofstream ofstr("table_etatopi0_dAu.tex");
OUTPUT_TABLE_3(ofstr, EtaToPi0MB_dAu,  dAu_etatopi0_MB,      dAu_etatopi0_MB_HT1,  stat, systA, systC, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, EtaToPi0HT1_dAu, dAu_etatopi0_MB_HT1,  dAu_etatopi0_HT1_HT2, stat, systA, systC, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, EtaToPi0HT2_dAu, dAu_etatopi0_HT1_HT2, dAu_etatopi0_HT2,     stat, systA, systC, FORMAT_R, true)
}
// Add errors for display
COPY_ARRAYS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, totsyst)
ADD_ERRORS_POINTS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, stat)
ADD_ERRORS_POINTS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, systA)
//ADD_ERRORS_POINTS(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, systC)
ADD_ERRORS_POINTS1(EtaToPi0MB_pp, EtaToPi0HT1_pp, EtaToPi0HT2_pp, totsyst, systC)
COPY_ARRAYS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, totsyst)
ADD_ERRORS_POINTS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, stat)
ADD_ERRORS_POINTS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, systA)
//ADD_ERRORS_POINTS(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, systC)
ADD_ERRORS_POINTS1(EtaToPi0MB_dAu, EtaToPi0HT1_dAu, EtaToPi0HT2_dAu, totsyst, systC)

// Systematic errors on pi0 and eta Rcp
COPY_ARR(x_err_RcpMB_dAuRcp, syst) COPY_ARR(x_err_RcpHT1_dAuRcp, syst) COPY_ARR(x_err_RcpHT2_dAuRcp, syst)
SET_X_ERRORS(RcpMB_dAuRcp_syst, RcpHT1_dAuRcp_syst, RcpHT2_dAuRcp_syst)
COPY_ARRAYS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, stat)
ZERO_ARRAY(y_err_RcpMB_dAuRcp) ZERO_ARRAY(y_err_RcpHT1_dAuRcp) ZERO_ARRAY(y_err_RcpHT2_dAuRcp)
COPY_ARRAYS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, systA)
COPY_ARRAYS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, systB)
COPY_ARRAYS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, systN)
COPY_ARR(x_err_RcpMBEta_dAuRcp, syst) COPY_ARR(x_err_RcpHT1Eta_dAuRcp, syst) COPY_ARR(x_err_RcpHT2Eta_dAuRcp, syst)
SET_X_ERRORS(RcpMBEta_dAuRcp_syst, RcpHT1Eta_dAuRcp_syst, RcpHT2Eta_dAuRcp_syst)
COPY_ARRAYS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, stat)
ZERO_ARRAY(y_err_RcpMBEta_dAuRcp) ZERO_ARRAY(y_err_RcpHT1Eta_dAuRcp) ZERO_ARRAY(y_err_RcpHT2Eta_dAuRcp)
COPY_ARRAYS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, systA)
COPY_ARRAYS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, systB)
COPY_ARRAYS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, systN)
// Combinatorial background normalization: add
ADD_ERR_POINTS(RcpMB_dAuRcp_systA, RcpHT1_dAuRcp_systA, RcpHT2_dAuRcp_systA, 0.005, S2*(0.03 - 0.005)/15)
ADD_ERR_POINTS(RcpMBEta_dAuRcp_systA, RcpHT1Eta_dAuRcp_systA, RcpHT2Eta_dAuRcp_systA, 0.10, S2*(0.50 - 0.10)/15)
// CPV correction: add
ADD_ERR_POINTS(RcpMB_dAuRcp_systN, RcpHT1_dAuRcp_systN, RcpHT2_dAuRcp_systN, S2*0.02, 0)
ADD_ERR_POINTS(RcpMBEta_dAuRcp_systN, RcpHT1Eta_dAuRcp_systN, RcpHT2Eta_dAuRcp_systN, S2*0.02, 0)
// Analysis cuts: add
ADD_ERR_POINTS(RcpMB_dAuRcp_systA, RcpHT1_dAuRcp_systA, RcpHT2_dAuRcp_systA, S2*0.05, 0)
ADD_ERR_POINTS(RcpMBEta_dAuRcp_systA, RcpHT1Eta_dAuRcp_systA, RcpHT2Eta_dAuRcp_systA, S2*0.05, 0)
// HighTower scale factors: add
ADD_ERR(x_RcpHT1_dAuRcp_systB, y_RcpHT1_dAuRcp_systB, y_err_RcpHT1_dAuRcp_systB, S2*0.03, 0)
ADD_ERR(x_RcpHT2_dAuRcp_systB, y_RcpHT2_dAuRcp_systB, y_err_RcpHT2_dAuRcp_systB, S2*0.05, 0)
ADD_ERR(x_RcpHT1Eta_dAuRcp_systB, y_RcpHT1Eta_dAuRcp_systB, y_err_RcpHT1Eta_dAuRcp_systB, S2*0.03, 0)
ADD_ERR(x_RcpHT2Eta_dAuRcp_systB, y_RcpHT2Eta_dAuRcp_systB, y_err_RcpHT2Eta_dAuRcp_systB, S2*0.05, 0)
// Vertex finding efficiency: 2.27%
ADD_ERR_POINTS(RcpMB_dAuRcp_systN, RcpHT1_dAuRcp_systN, RcpHT2_dAuRcp_systN, 0.0227, 0)
ADD_ERR_POINTS(RcpMBEta_dAuRcp_systN, RcpHT1Eta_dAuRcp_systN, RcpHT2Eta_dAuRcp_systN, 0.0227, 0)
// <Ncoll> in most central: 1.1/15.0
ADD_ERR_POINTS(RcpMB_dAuRcp_systN, RcpHT1_dAuRcp_systN, RcpHT2_dAuRcp_systN, 1.1/15.0, 0)
ADD_ERR_POINTS(RcpMBEta_dAuRcp_systN, RcpHT1Eta_dAuRcp_systN, RcpHT2Eta_dAuRcp_systN, 1.1/15.0, 0)
// <Ncoll> in most peripheral: 0.3/4.0
ADD_ERR_POINTS(RcpMB_dAuRcp_systN, RcpHT1_dAuRcp_systN, RcpHT2_dAuRcp_systN, 0.3/4.0, 0)
ADD_ERR_POINTS(RcpMBEta_dAuRcp_systN, RcpHT1Eta_dAuRcp_systN, RcpHT2Eta_dAuRcp_systN, 0.3/4.0, 0)
// Norm. errors
Float_t x_RcpMB_dAuRcp_norm[] = {16.75}; Float_t x_err_RcpMB_dAuRcp_norm[] = {0.25}; Float_t y_RcpMB_dAuRcp_norm[] = {1.0}; Float_t y_err_RcpMB_dAuRcp_norm[1]; y_err_RcpMB_dAuRcp_norm[0] = y_err_RcpMB_dAuRcp_systN[0] / y_RcpMB_dAuRcp_systN[0];
Float_t x_RcpMBEta_dAuRcp_norm[] = {16.75}; Float_t x_err_RcpMBEta_dAuRcp_norm[] = {0.25}; Float_t y_RcpMBEta_dAuRcp_norm[] = {1.0}; Float_t y_err_RcpMBEta_dAuRcp_norm[1]; y_err_RcpMBEta_dAuRcp_norm[0] = y_err_RcpMBEta_dAuRcp_systN[0] / y_RcpMBEta_dAuRcp_systN[0];
// Write tables
{
ofstream ofstr("table_pi0_Rcp.tex");
OUTPUT_TABLE_3(ofstr, RcpMB_dAuRcp,  pp_MB,      pp_MB_HT1,  stat, systA, systB, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RcpHT1_dAuRcp, pp_MB_HT1,  pp_HT1_HT2, stat, systA, systB, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RcpHT2_dAuRcp, pp_HT1_HT2, pp_HT2,     stat, systA, systB, FORMAT_R, true)
}
{
ofstream ofstr("table_eta_Rcp.tex");
OUTPUT_TABLE_3(ofstr, RcpMBEta_dAuRcp,  dAu_MB,      dAu_MB_HT1,  stat, systA, systB, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RcpHT1Eta_dAuRcp, dAu_MB_HT1,  dAu_HT1_HT2, stat, systA, systB, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RcpHT2Eta_dAuRcp, dAu_HT1_HT2, dAu_HT2,     stat, systA, systB, FORMAT_R, true)
}
// Add errors for display
COPY_ARRAYS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, totsyst)
ADD_ERRORS_POINTS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, stat)
ADD_ERRORS_POINTS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, systA)
//ADD_ERRORS_POINTS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, systB)
ADD_ERRORS_POINTS1(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, totsyst, systB)
//ADD_ERRORS_POINTS(RcpMB_dAuRcp, RcpHT1_dAuRcp, RcpHT2_dAuRcp, systN)
COPY_ARRAYS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, totsyst)
ADD_ERRORS_POINTS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, stat)
ADD_ERRORS_POINTS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, systA)
//ADD_ERRORS_POINTS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, systB)
ADD_ERRORS_POINTS1(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, totsyst, systB)
//ADD_ERRORS_POINTS(RcpMBEta_dAuRcp, RcpHT1Eta_dAuRcp, RcpHT2Eta_dAuRcp, systN)

// Systematic errors on pi0 and eta RdA
COPY_ARR(x_err_RDAMB_RdA, syst) COPY_ARR(x_err_RDAHT1_RdA, syst) COPY_ARR(x_err_RDAHT2_RdA, syst)
SET_X_ERRORS(RDAMB_RdA_syst, RDAHT1_RdA_syst, RDAHT2_RdA_syst)
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, stat)
ZERO_ARRAY(y_err_RDAMB_RdA) ZERO_ARRAY(y_err_RDAHT1_RdA) ZERO_ARRAY(y_err_RDAHT2_RdA)
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systA)
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systB)
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systE)
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systBE)
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systN)
COPY_ARR(x_err_RDAMBEta_RdA, syst) COPY_ARR(x_err_RDAHT1Eta_RdA, syst) COPY_ARR(x_err_RDAHT2Eta_RdA, syst)
SET_X_ERRORS(RDAMBEta_RdA_syst, RDAHT1Eta_RdA_syst, RDAHT2Eta_RdA_syst)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, stat)
ZERO_ARRAY(y_err_RDAMBEta_RdA) ZERO_ARRAY(y_err_RDAHT1Eta_RdA) ZERO_ARRAY(y_err_RDAHT2Eta_RdA)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systA)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systB)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systE)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systBE)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systN)
// Combinatorial background normalization: add
ADD_ERR_POINTS(RDAMB_RdA_systA, RDAHT1_RdA_systA, RDAHT2_RdA_systA, 0.005, S2*(0.03 - 0.005)/15)
ADD_ERR_POINTS(RDAMBEta_RdA_systA, RDAHT1Eta_RdA_systA, RDAHT2Eta_RdA_systA, 0.10, S2*(0.50 - 0.10)/15)
// CPV correction: add
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, S2*0.02, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, S2*0.02, 0)
// Missing material in simulation correction: add
ADD_ERR_POINTS(RDAMB_RdA_systB, RDAHT1_RdA_systB, RDAHT2_RdA_systB, 0.02/1.10, 0)
ADD_ERR_POINTS(RDAMB_RdA_systB, RDAHT1_RdA_systB, RDAHT2_RdA_systB, 0.03/1.15, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systB, RDAHT1Eta_RdA_systB, RDAHT2Eta_RdA_systB, 0.02/1.10, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systB, RDAHT1Eta_RdA_systB, RDAHT2Eta_RdA_systB, 0.03/1.15, 0)
// SMD scale uncertainty
ADD_ERR(x_RDAHT1_RdA_systB, y_RDAHT1_RdA_systB, y_err_RDAHT1_RdA_systB, S2*0.04, S2*(-0.0025))
ADD_ERR(x_RDAHT2_RdA_systB, y_RDAHT2_RdA_systB, y_err_RDAHT2_RdA_systB, S2*0.04, S2*(-0.0025))
ADD_ERR(x_RDAHT1Eta_RdA_systB, y_RDAHT1Eta_RdA_systB, y_err_RDAHT1Eta_RdA_systB, S2*0.03, 0)
ADD_ERR(x_RDAHT2Eta_RdA_systB, y_RDAHT2Eta_RdA_systB, y_err_RDAHT2Eta_RdA_systB, S2*0.05, 0)
// Analysis cuts: add
ADD_ERR_POINTS(RDAMB_RdA_systA, RDAHT1_RdA_systA, RDAHT2_RdA_systA, S2*0.05, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systA, RDAHT1Eta_RdA_systA, RDAHT2Eta_RdA_systA, S2*0.05, 0)
// HighTower scale factors: add
ADD_ERR(x_RDAHT1_RdA_systB, y_RDAHT1_RdA_systB, y_err_RDAHT1_RdA_systB, S2*0.03, 0)
ADD_ERR(x_RDAHT2_RdA_systB, y_RDAHT2_RdA_systB, y_err_RDAHT2_RdA_systB, S2*0.05, 0)
ADD_ERR(x_RDAHT1Eta_RdA_systB, y_RDAHT1Eta_RdA_systB, y_err_RDAHT1Eta_RdA_systB, S2*0.03, 0)
ADD_ERR(x_RDAHT2Eta_RdA_systB, y_RDAHT2Eta_RdA_systB, y_err_RDAHT2Eta_RdA_systB, S2*0.05, 0)
// BBC cross section in p+p: 6.94%
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, 0.0694, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, 0.0694, 0)
// BBC efficiency in p+p: 9.20%
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, 0.0920, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, 0.0920, 0)
// Vertex finding efficiency in d+Au: 1.08%
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, 0.0108, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, 0.0108, 0)
// MinBias cross section in d+Au: 4.07%
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, 0.0407, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, 0.0407, 0)
// MinBias efficiency in d+Au: 3.16%
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, 0.0316, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, 0.0316, 0)
// <Ncoll> in min. bias d+Au: 0.4/7.5
ADD_ERR_POINTS(RDAMB_RdA_systN, RDAHT1_RdA_systN, RDAHT2_RdA_systN, 0.4/7.5, 0)
ADD_ERR_POINTS(RDAMBEta_RdA_systN, RDAHT1Eta_RdA_systN, RDAHT2Eta_RdA_systN, 0.4/7.5, 0)
// Energy scale: add!
//ADD_ERR_POINTS(RDAMB_RdA_systE, RDAHT1_RdA_systE, RDAHT2_RdA_systE, S2*0.40, 0)
//ADD_ERR_POINTS(RDAMBEta_RdA_systE, RDAHT1Eta_RdA_systE, RDAHT2Eta_RdA_systE, S2*0.40, 0)
// Energy scale: add arrays
ADD_ERR_POINTS_ARRAY(RDAMB_RdA_systE, RDAHT1_RdA_systE, RDAHT2_RdA_systE, escale_err_pp_MB, escale_err_pp_HT1, escale_err_pp_HT2)
ADD_ERR_POINTS_ARRAY(RDAMB_RdA_systE, RDAHT1_RdA_systE, RDAHT2_RdA_systE, escale_err_dAu_MB, escale_err_dAu_HT1, escale_err_dAu_HT2)
ADD_ERR_POINTS_ARRAY(RDAMBEta_RdA_systE, RDAHT1Eta_RdA_systE, RDAHT2Eta_RdA_systE, escale_err_pp_MB, escale_err_pp_HT1, escale_err_pp_HT2)
ADD_ERR_POINTS_ARRAY(RDAMBEta_RdA_systE, RDAHT1Eta_RdA_systE, RDAHT2Eta_RdA_systE, escale_err_dAu_MB, escale_err_dAu_HT1, escale_err_dAu_HT2)
// Combine B+E->BE
ADD_ERRORS_POINTS1(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systBE, systB)
ADD_ERRORS_POINTS1(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systBE, systE)
ADD_ERRORS_POINTS1(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systBE, systB)
ADD_ERRORS_POINTS1(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systBE, systE)
// Norm. errors
Float_t x_RDAMB_RdA_norm[] = {16.75}; Float_t x_err_RDAMB_RdA_norm[] = {0.25}; Float_t y_RDAMB_RdA_norm[] = {1.0}; Float_t y_err_RDAMB_RdA_norm[1]; y_err_RDAMB_RdA_norm[0] = y_err_RDAMB_RdA_systN[0] / y_RDAMB_RdA_systN[0];
Float_t x_RDAMBEta_RdA_norm[] = {16.75}; Float_t x_err_RDAMBEta_RdA_norm[] = {0.25}; Float_t y_RDAMBEta_RdA_norm[] = {1.0}; Float_t y_err_RDAMBEta_RdA_norm[1]; y_err_RDAMBEta_RdA_norm[0] = y_err_RDAMBEta_RdA_systN[0] / y_RDAMBEta_RdA_systN[0];
// Write tables
{
ofstream ofstr("table_pi0_RdA.tex");
OUTPUT_TABLE_3(ofstr, RDAMB_RdA,  RdA_MB,      RdA_MB_HT1,  stat, systA, systBE, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RDAHT1_RdA, RdA_MB_HT1,  RdA_HT1_HT2, stat, systA, systBE, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RDAHT2_RdA, RdA_HT1_HT2, RdA_HT2,     stat, systA, systBE, FORMAT_R, true)
}
{
ofstream ofstr("table_eta_RdA.tex");
OUTPUT_TABLE_3(ofstr, RDAMBEta_RdA, RdA_eta_MB,       RdA_eta_MB_HT1,  stat, systA, systBE, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RDAHT1Eta_RdA, RdA_eta_MB_HT1,  RdA_eta_HT1_HT2, stat, systA, systBE, FORMAT_R, true)
OUTPUT_TABLE_3(ofstr, RDAHT2Eta_RdA, RdA_eta_HT1_HT2, RdA_eta_HT2,     stat, systA, systBE, FORMAT_R, true)
}
// Add errors for display
COPY_ARRAYS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, esyst)
ADD_ERRORS_POINTS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, stat)
ADD_ERRORS_POINTS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systA)
//ADD_ERRORS_POINTS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systB)
//ADD_ERRORS_POINTS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systE)
//ADD_ERRORS_POINTS(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, systN)
ADD_ERRORS_POINTS1(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, esyst, systB)
ADD_ERRORS_POINTS1(RDAMB_RdA, RDAHT1_RdA, RDAHT2_RdA, esyst, systE)
COPY_ARRAYS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, esyst)
ADD_ERRORS_POINTS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, stat)
ADD_ERRORS_POINTS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systA)
//ADD_ERRORS_POINTS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systB)
//ADD_ERRORS_POINTS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systE)
//ADD_ERRORS_POINTS(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, systN)
ADD_ERRORS_POINTS1(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, esyst, systB)
ADD_ERRORS_POINTS1(RDAMBEta_RdA, RDAHT1Eta_RdA, RDAHT2Eta_RdA, esyst, systE)

	    bin_stat_list_type tmp;

	    const Char_t *titlePPMBPi0 = "#font[12]{p}_{ }+_{ }#font[12]{p} MinBias";
	    TDataPoints dataPPMBPi0Invyield = TDataPoints(TDrawOptions(TNamed("pp_MB_pi0_invyield", titlePPMBPi0), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 1.0), "P E Z", "P E Z", TString(titlePPMBPi0)), x_MBspectrum_pp, y_MBspectrum_pp, x_err_MBspectrum_pp, y_err_MBspectrum_pp, sizeof(x_MBspectrum_pp)/sizeof(x_MBspectrum_pp[0]));
	    bin_stat_list_type bins_pp_pi0_invyield_MB = dataPPMBPi0Invyield.getBins(pp_MB, pp_MB_HT1);
	    const Char_t *titlePPHT1Pi0 = "#font[12]{p}_{ }+_{ }#font[12]{p} HighTower-1";
	    TDataPoints dataPPHT1Pi0Invyield = TDataPoints(TDrawOptions(TNamed("pp_HT1_pi0_invyield", titlePPHT1Pi0), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kOpenSquare, 1.0), "P E Z", "P E Z", TString(titlePPHT1Pi0)), x_HT1spectrum_pp, y_HT1spectrum_pp, x_err_HT1spectrum_pp, y_err_HT1spectrum_pp, sizeof(x_HT1spectrum_pp)/sizeof(x_HT1spectrum_pp[0]));
	    bin_stat_list_type bins_pp_pi0_invyield_HT1 = dataPPHT1Pi0Invyield.getBins(pp_MB_HT1, pp_HT1_HT2);
	    const Char_t *titlePPHT2Pi0 = "#font[12]{p}_{ }+_{ }#font[12]{p} HighTower-2";
	    TDataPoints dataPPHT2Pi0Invyield = TDataPoints(TDrawOptions(TNamed("pp_HT2_pi0_invyield", titlePPHT2Pi0), TAttLine(kRed, kSolid, 2), TAttFill(kRed, 0), TAttMarker(kRed, kOpenTriangleUp, 1.0), "P E Z", "P E Z", TString(titlePPHT2Pi0)), x_HT2spectrum_pp, y_HT2spectrum_pp, x_err_HT2spectrum_pp, y_err_HT2spectrum_pp, sizeof(x_HT2spectrum_pp)/sizeof(x_HT2spectrum_pp[0]));
	    bin_stat_list_type bins_pp_pi0_invyield_HT2 = dataPPHT2Pi0Invyield.getBins(pp_HT1_HT2, pp_HT2);

	    bin_stat_list_type bins_pp_pi0_invyield;
	    bins_pp_pi0_invyield.merge(bins_pp_pi0_invyield_MB);
	    bins_pp_pi0_invyield.merge(bins_pp_pi0_invyield_HT1);
	    bins_pp_pi0_invyield.merge(bins_pp_pi0_invyield_HT2);
	    TWeightCalculator wPPPi0Invyield("wPPPi0Invyield", "#font[12]{p}+#font[12]{p} #pi^{0} invariant yield");
	    wPPPi0Invyield.Fit(bins_pp_pi0_invyield, "RQN LL M", "POWERLAW", 0, 100);
	    fPPPi0Invyield = wPPPi0Invyield.createFunc(true);
	    fPPPi0Invyield->SetTitle("Fit to #font[12]{p}_{ }+_{ }#font[12]{p}");
	    fPPPi0Invyield->SetRange(0.0, 20.0);

	    for (UInt_t i = 0;i < sizeof(x_MBspectrum_pp)/sizeof(x_MBspectrum_pp[0]);i++) {
		y_err_MBspectrum_pp_smdsyst[i] *= fPPPi0Invyield->Eval(x_MBspectrum_pp[i]) / y_MBspectrum_pp_smdsyst[i];
		y_MBspectrum_pp_smdsyst[i] *= fPPPi0Invyield->Eval(x_MBspectrum_pp[i]) / y_MBspectrum_pp_smdsyst[i];
	    };
	    tmp = TDataPoints(x_MBspectrum_pp, y_MBspectrum_pp_smdsyst, x_err_MBspectrum_pp, y_err_MBspectrum_pp_smdsyst, sizeof(x_MBspectrum_pp)/sizeof(x_MBspectrum_pp[0])).getBins();
	    const Char_t *titlePPMBPi0_smdsyst = "Additional SMD uncertainty";
	    TDataPoints dataPPMBPi0Invyield_smdsyst_1 = TDataPoints(TDrawOptions(TNamed("pp_MB_pi0_invyield_smdsyst_1", titlePPMBPi0_smdsyst), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 1.0), "L", "L", TString(titlePPMBPi0_smdsyst)), tmp, "LOWEDGE");
	    TDataPoints dataPPMBPi0Invyield_smdsyst_2 = TDataPoints(TDrawOptions(TNamed("pp_MB_pi0_invyield_smdsyst_2", ""), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 1.0), "L", "L", TString()), tmp, "HIGHEDGE");

	    dataPoints_pi0_invyield.push_back(dataPPMBPi0Invyield);
	    dataPoints_pi0_pp_invyield.push_back(dataPPMBPi0Invyield);
	    dataPoints_pi0_invyield.push_back(dataPPHT1Pi0Invyield);
	    dataPoints_pi0_pp_invyield.push_back(dataPPHT1Pi0Invyield);
	    dataPoints_pi0_invyield.push_back(dataPPHT2Pi0Invyield);
	    dataPoints_pi0_pp_invyield.push_back(dataPPHT2Pi0Invyield);
	    dataPoints_pi0_invyield.push_back(TDataPoints(TDrawOptions(TNamed("pp_pi0_invyield_fit", fPPPi0Invyield->GetTitle()), TAttLine(kBlack, kDashed, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPPi0Invyield->GetTitle())), fPPPi0Invyield));
	    dataPoints_pi0_pp_invyield.push_back(TDataPoints(TDrawOptions(TNamed("pp_pi0_invyield_fit", fPPPi0Invyield->GetTitle()), TAttLine(kBlack, kDashed, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPPi0Invyield->GetTitle())), fPPPi0Invyield));
	    //dataPoints_pi0_pp_invyield.push_back(dataPPMBPi0Invyield_smdsyst_1);
	    //dataPoints_pi0_pp_invyield.push_back(dataPPMBPi0Invyield_smdsyst_2);

	    bin_stat_list_type bins_pp_pi0_crossection;
	    tmp = TDataPoints(x_MBcrossection_pp, y_MBcrossection_pp, x_err_MBcrossection_pp, y_err_MBcrossection_pp_stat, sizeof(x_MBcrossection_pp)/sizeof(x_MBcrossection_pp[0])).getBins(pp_MB, pp_MB_HT1);
	    bins_pp_pi0_crossection.merge(tmp);
	    tmp = TDataPoints(x_HT1crossection_pp, y_HT1crossection_pp, x_err_HT1crossection_pp, y_err_HT1crossection_pp_stat, sizeof(x_HT1crossection_pp)/sizeof(x_HT1crossection_pp[0])).getBins(pp_MB_HT1, pp_HT1_HT2);
	    bins_pp_pi0_crossection.merge(tmp);
	    tmp = TDataPoints(x_HT2crossection_pp, y_HT2crossection_pp, x_err_HT2crossection_pp, y_err_HT2crossection_pp_stat, sizeof(x_HT2crossection_pp)/sizeof(x_HT2crossection_pp[0])).getBins(pp_HT1_HT2, pp_HT2);
	    bins_pp_pi0_crossection.merge(tmp);
	    const Char_t *titlePPPi0 = "#font[12]{p}_{ }+_{ }#font[12]{p}";
	    TDataPoints dataPPPi0ThisCrossection = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection", titlePPPi0), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.2), "P E Z", "P E Z", TString(titlePPPi0)), bins_pp_pi0_crossection);
	    TDataPoints dataPPPi0ThisCrossection_unnamed = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection", ""), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.2), "P E Z", "P E Z", TString()), bins_pp_pi0_crossection);
	    const Char_t *titlePPPi0_1 = "#pi^{0}^{} (this analysis)";
	    TDataPoints dataPPPi0ThisCrossection_name1 = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection", titlePPPi0_1), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.2), "P E Z", "P E Z", TString(titlePPPi0_1)), bins_pp_pi0_crossection);

	    bin_stat_list_type bins_pp_pi0_crossection_esyst;
	    tmp = TDataPoints(x_MBcrossection_pp, y_MBcrossection_pp, x_err_MBcrossection_pp_syst, y_err_MBcrossection_pp_esyst, sizeof(x_MBcrossection_pp)/sizeof(x_MBcrossection_pp[0])).getBins(pp_MB, pp_MB_HT1);
	    bins_pp_pi0_crossection_esyst.merge(tmp);
	    tmp = TDataPoints(x_HT1crossection_pp, y_HT1crossection_pp, x_err_HT1crossection_pp_syst, y_err_HT1crossection_pp_esyst, sizeof(x_HT1crossection_pp)/sizeof(x_HT1crossection_pp[0])).getBins(pp_MB_HT1, pp_HT1_HT2);
	    bins_pp_pi0_crossection_esyst.merge(tmp);
	    tmp = TDataPoints(x_HT2crossection_pp, y_HT2crossection_pp, x_err_HT2crossection_pp_syst, y_err_HT2crossection_pp_esyst, sizeof(x_HT2crossection_pp)/sizeof(x_HT2crossection_pp[0])).getBins(pp_HT1_HT2, pp_HT2);
	    bins_pp_pi0_crossection_esyst.merge(tmp);
	    TDataPoints dataPPPi0ThisCrossectionEsyst1 = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_esyst1", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullCircle, 1.2), "E2 Z", "E2 Z", TString()), bins_pp_pi0_crossection_esyst, ""/*"LOWEDGE"*/);
	    //TDataPoints dataPPPi0ThisCrossectionEsyst2 = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_esyst2", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullCircle, 1.2), "L", "L", TString()), bins_pp_pi0_crossection_esyst, ""/*"HIGHEDGE"*/);
	    TDataPoints dataPPPi0ThisCrossectionNorm = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_norm", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullCircle, 1.2), "E2 Z", "PF", TString()), x_MBcrossection_pp_norm, y_MBcrossection_pp_norm, x_err_MBcrossection_pp_norm, y_err_MBcrossection_pp_norm, sizeof(x_MBcrossection_pp_norm)/sizeof(x_MBcrossection_pp_norm[0]));
	    TDataPoints dataPPPi0ThisCrossectionNorm1 = TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_norm1", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullCircle, 1.2), "E2 Z", "PF", TString()), x_MBcrossection_pp_norm1, y_MBcrossection_pp_norm1, x_err_MBcrossection_pp_norm1, y_err_MBcrossection_pp_norm1, sizeof(x_MBcrossection_pp_norm1)/sizeof(x_MBcrossection_pp_norm1[0]));

	    TWeightCalculator wPPPi0Crossection("wPPPi0Crossection", "#font[12]{p}+#font[12]{p} #pi^{0} cross section");
	    wPPPi0Crossection.Fit(bins_pp_pi0_crossection, "RQN LL M", "POWERLAW", 0, 100);
	    fPPPi0Crossection = wPPPi0Crossection.createFunc(true);
	    fPPPi0Crossection->SetTitle("Fit to #font[12]{p}_{ }+_{ }#font[12]{p}");
	    fPPPi0Crossection->SetRange(0.0, 20.0);
	    //fPPPi0Crossection->Print("");

	    dataPoints_pi0_pp_crossection.push_back(dataPPPi0ThisCrossection);
	    dataPoints_pi0_pp_crossection_1.push_back(dataPPPi0ThisCrossection_name1);
	    dataPoints_pi0_pp_crossection_div.push_back(dataPPPi0ThisCrossectionEsyst1);
	    //dataPoints_pi0_pp_crossection_div.push_back(dataPPPi0ThisCrossectionEsyst2);
	    dataPoints_pi0_pp_crossection_div_2.push_back(dataPPPi0ThisCrossection_unnamed);
	    dataPoints_pi0_pp_crossection_div_2.push_back(dataPPPi0ThisCrossectionNorm);
	    //dataPoints_pi0_pp_crossection.push_back(TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_fit", fPPPi0Crossection->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPPi0Crossection->GetTitle())), fPPPi0Crossection));
	    //dataPoints_pi0_pp_crossection_div_2.push_back(TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_fit", fPPPi0Crossection->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPPi0Crossection->GetTitle())), fPPPi0Crossection));
	    dataPoints_pi0_pp_crossection_div1.push_back(dataPPPi0ThisCrossectionEsyst1);
	    //dataPoints_pi0_pp_crossection_div1.push_back(dataPPPi0ThisCrossectionEsyst2);
	    dataPoints_pi0_pp_crossection_div1_2.push_back(dataPPPi0ThisCrossection_name1);
	    dataPoints_pi0_pp_crossection_div1_2.push_back(dataPPPi0ThisCrossectionNorm1);

	    //TGraphErrors *pi0_pp_cx_gr = createGraph(bins_pp_pi0_crossection);
	    //fPPPi0CrossectionLevy = new TF1("fPPPi0CrossectionLevy", LevyFunc);
	    //fPPPi0CrossectionLevy->SetTitle("Levy fit to #font[12]{p}+#font[12]{p}");
	    //fPPPi0CrossectionLevy->SetRange(0.0, 20.0);
	    //fPPPi0CrossectionLevy->SetParameter(0, 100);
	    //fPPPi0CrossectionLevy->SetParLimits(0, 1, 50000);
	    //fPPPi0CrossectionLevy->SetParameter(1, 1);
	    //fPPPi0CrossectionLevy->SetParLimits(1, 0, 10);
	    //fPPPi0CrossectionLevy->SetParameter(2, 10);
	    //fPPPi0CrossectionLevy->SetParLimits(2, 5, 20);
	    //pi0_pp_cx_gr->Fit(fPPPi0CrossectionLevy, "RN LL M");
	    //dataPoints_pi0_pp_crossection.push_back(TDataPoints(TDrawOptions(TNamed("pp_pi0_crossection_levyfit", fPPPi0CrossectionLevy->GetTitle()), TAttLine(kRed, kDotted, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPPi0CrossectionLevy->GetTitle())), fPPPi0CrossectionLevy));

	    const Char_t *titlePPMBEta = "#font[12]{p}_{ }+_{ }#font[12]{p} MinBias";
	    TDataPoints dataPPMBEtaInvyield = TDataPoints(TDrawOptions(TNamed("pp_MB_eta_invyield", titlePPMBEta), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kOpenCircle, 1.0), "P E Z", "P E Z", TString(titlePPMBEta)), x_EtaMBspectrum_pp, y_EtaMBspectrum_pp, x_err_EtaMBspectrum_pp, y_err_EtaMBspectrum_pp, sizeof(x_EtaMBspectrum_pp)/sizeof(x_EtaMBspectrum_pp[0]));
	    dataPoints_eta_invyield.push_back(dataPPMBEtaInvyield);
	    dataPoints_eta_pp_invyield.push_back(dataPPMBEtaInvyield);
	    bin_stat_list_type bins_pp_eta_invyield_MB = dataPPMBEtaInvyield.getBins(pp_etatopi0_MB, pp_etatopi0_MB_HT1);

	    const Char_t *titlePPHT1Eta = "#font[12]{p}_{ }+_{ }#font[12]{p} HighTower-1";
	    TDataPoints dataPPHT1EtaInvyield = TDataPoints(TDrawOptions(TNamed("pp_HT1_eta_invyield", titlePPHT1Eta), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kOpenSquare, 1.0), "P E Z", "P E Z", TString(titlePPHT1Eta)), x_EtaHT1spectrum_pp, y_EtaHT1spectrum_pp, x_err_EtaHT1spectrum_pp, y_err_EtaHT1spectrum_pp, sizeof(x_EtaHT1spectrum_pp)/sizeof(x_EtaHT1spectrum_pp[0]));
	    dataPoints_eta_invyield.push_back(dataPPHT1EtaInvyield);
	    dataPoints_eta_pp_invyield.push_back(dataPPHT1EtaInvyield);
	    bin_stat_list_type bins_pp_eta_invyield_HT1 = dataPPHT1EtaInvyield.getBins(pp_etatopi0_MB_HT1, pp_etatopi0_HT1_HT2);

	    const Char_t *titlePPHT2Eta = "#font[12]{p}_{ }+_{ }#font[12]{p} HighTower-2";
	    TDataPoints dataPPHT2EtaInvyield = TDataPoints(TDrawOptions(TNamed("pp_HT2_eta_invyield", titlePPHT2Eta), TAttLine(kRed, kSolid, 2), TAttFill(kRed, 0), TAttMarker(kRed, kOpenTriangleUp, 1.0), "P E Z", "P E Z", TString(titlePPHT2Eta)), x_EtaHT2spectrum_pp, y_EtaHT2spectrum_pp, x_err_EtaHT2spectrum_pp, y_err_EtaHT2spectrum_pp, sizeof(x_EtaHT2spectrum_pp)/sizeof(x_EtaHT2spectrum_pp[0]));
	    dataPoints_eta_invyield.push_back(dataPPHT2EtaInvyield);
	    dataPoints_eta_pp_invyield.push_back(dataPPHT2EtaInvyield);
	    bin_stat_list_type bins_pp_eta_invyield_HT2 = dataPPHT2EtaInvyield.getBins(pp_etatopi0_HT1_HT2, pp_etatopi0_HT2);

	    bin_stat_list_type bins_pp_eta_invyield;
	    bins_pp_eta_invyield.merge(bins_pp_eta_invyield_MB);
	    bins_pp_eta_invyield.merge(bins_pp_eta_invyield_HT1);
	    bins_pp_eta_invyield.merge(bins_pp_eta_invyield_HT2);
	    TWeightCalculator wPPEtaInvyield("wPPEtaInvyield", "#font[12]{p}+#font[12]{p} #eta invariant yield");
	    wPPEtaInvyield.Fit(bins_pp_eta_invyield, "RQN LL M", "POWERLAW", 0, 100);
	    fPPEtaInvyield = wPPEtaInvyield.createFunc(true);
	    fPPEtaInvyield->SetTitle("Fit to #font[12]{p}_{ }+_{ }#font[12]{p}");
	    fPPEtaInvyield->SetRange(0.0, 20.0);
	    dataPoints_eta_invyield.push_back(TDataPoints(TDrawOptions(TNamed("pp_eta_invyield_fit", fPPEtaInvyield->GetTitle()), TAttLine(kBlack, kDashed, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPEtaInvyield->GetTitle())), fPPEtaInvyield));
	    dataPoints_eta_pp_invyield.push_back(TDataPoints(TDrawOptions(TNamed("pp_eta_invyield_fit", fPPEtaInvyield->GetTitle()), TAttLine(kBlack, kDashed, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fPPEtaInvyield->GetTitle())), fPPEtaInvyield));

	    bin_stat_list_type bins_etatopi_pp;
	    tmp = TDataPoints(x_EtaToPi0MB_pp, y_EtaToPi0MB_pp, x_err_EtaToPi0MB_pp, y_err_EtaToPi0MB_pp, sizeof(x_EtaToPi0MB_pp)/sizeof(x_EtaToPi0MB_pp[0])).getBins(pp_etatopi0_MB, pp_etatopi0_MB_HT1);
	    bins_etatopi_pp.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT1_pp, y_EtaToPi0HT1_pp, x_err_EtaToPi0HT1_pp, y_err_EtaToPi0HT1_pp, sizeof(x_EtaToPi0HT1_pp)/sizeof(x_EtaToPi0HT1_pp[0])).getBins(pp_etatopi0_MB_HT1, pp_etatopi0_HT1_HT2);
	    bins_etatopi_pp.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT2_pp, y_EtaToPi0HT2_pp, x_err_EtaToPi0HT2_pp, y_err_EtaToPi0HT2_pp, sizeof(x_EtaToPi0HT2_pp)/sizeof(x_EtaToPi0HT2_pp[0])).getBins(pp_etatopi0_HT1_HT2, pp_etatopi0_HT2);
	    bins_etatopi_pp.merge(tmp);
	    const Char_t *titlePPEtaToPi0 = "#font[12]{p}_{ }+_{ }#font[12]{p} (this analysis)";
	    TDataPoints dataPPEtaToPi0 = TDataPoints(TDrawOptions(TNamed("pp_etaToPi0", titlePPEtaToPi0), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "P E Z", "P E Z", TString(titlePPEtaToPi0)), bins_etatopi_pp);
	    dataPoints_etaToPi0_pp.push_back(dataPPEtaToPi0);

	    bin_stat_list_type bins_etatopi_pp_syst;
	    tmp = TDataPoints(x_EtaToPi0MB_pp, y_EtaToPi0MB_pp, x_err_EtaToPi0MB_pp_syst, y_err_EtaToPi0MB_pp_totsyst, sizeof(x_EtaToPi0MB_pp)/sizeof(x_EtaToPi0MB_pp[0])).getBins(pp_etatopi0_MB, pp_etatopi0_MB_HT1);
	    bins_etatopi_pp_syst.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT1_pp, y_EtaToPi0HT1_pp, x_err_EtaToPi0HT1_pp_syst, y_err_EtaToPi0HT1_pp_totsyst, sizeof(x_EtaToPi0HT1_pp)/sizeof(x_EtaToPi0HT1_pp[0])).getBins(pp_etatopi0_MB_HT1, pp_etatopi0_HT1_HT2);
	    bins_etatopi_pp_syst.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT2_pp, y_EtaToPi0HT2_pp, x_err_EtaToPi0HT2_pp_syst, y_err_EtaToPi0HT2_pp_totsyst, sizeof(x_EtaToPi0HT2_pp)/sizeof(x_EtaToPi0HT2_pp[0])).getBins(pp_etatopi0_HT1_HT2, pp_etatopi0_HT2);
	    bins_etatopi_pp_syst.merge(tmp);
	    TDataPoints dataPPEtaToPi0_syst = TDataPoints(TDrawOptions(TNamed("pp_etaToPi0_syst", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullCircle, 1.0), "E2 Z", "E2 Z", TString()), bins_etatopi_pp_syst);
	    dataPoints_etaToPi0_pp1.push_back(dataPPEtaToPi0_syst);

	    const Char_t *titleDAuMBPi0 = "#font[12]{d}_{ }+_{ }Au MinBias";
	    TDataPoints dataDAuMBPi0Invyield = TDataPoints(TDrawOptions(TNamed("dAu_MB_pi0_invyield", titleDAuMBPi0), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "P E Z", "P E Z", TString(titleDAuMBPi0)), x_MBspectrum_dAu, y_MBspectrum_dAu, x_err_MBspectrum_dAu, y_err_MBspectrum_dAu, sizeof(x_MBspectrum_dAu)/sizeof(x_MBspectrum_dAu[0]));
	    const Char_t *titleDAuMBPi0Centr = "#font[12]{d}_{ }+_{ }Au MinBias";
	    TDataPoints dataDAuMBPi0InvyieldCentr = TDataPoints(TDrawOptions(TNamed("dAu_MB_pi0_invyield", titleDAuMBPi0Centr), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "PL E Z", "PL E Z", TString(titleDAuMBPi0Centr)), x_MBspectrum_dAu, y_MBspectrum_dAu, x_err_MBspectrum_dAu, y_err_MBspectrum_dAu, sizeof(x_MBspectrum_dAu)/sizeof(x_MBspectrum_dAu[0]));
	    bin_stat_list_type bins_dAu_pi0_invyield_MB = dataDAuMBPi0Invyield.getBins(dAu_MB, dAu_MB_HT1);
	    const Char_t *titleDAuHT1Pi0 = "#font[12]{d}_{ }+_{ }Au HighTower-1";
	    TDataPoints dataDAuHT1Pi0Invyield = TDataPoints(TDrawOptions(TNamed("dAu_HT1_pi0_invyield", titleDAuHT1Pi0), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kFullSquare, 1.0), "P E Z", "P E Z", TString(titleDAuHT1Pi0)), x_HT1spectrum_dAu, y_HT1spectrum_dAu, x_err_HT1spectrum_dAu, y_err_HT1spectrum_dAu, sizeof(x_HT1spectrum_dAu)/sizeof(x_HT1spectrum_dAu[0]));
	    const Char_t *titleDAuHT1Pi0Centr = "#font[12]{d}_{ }+_{ }Au HighTower-1";
	    TDataPoints dataDAuHT1Pi0InvyieldCentr = TDataPoints(TDrawOptions(TNamed("dAu_HT1_pi0_invyield", titleDAuHT1Pi0Centr), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kFullSquare, 1.0), "PL E Z", "PL E Z", TString(titleDAuHT1Pi0Centr)), x_HT1spectrum_dAu, y_HT1spectrum_dAu, x_err_HT1spectrum_dAu, y_err_HT1spectrum_dAu, sizeof(x_HT1spectrum_dAu)/sizeof(x_HT1spectrum_dAu[0]));
	    bin_stat_list_type bins_dAu_pi0_invyield_HT1 = dataDAuHT1Pi0Invyield.getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    const Char_t *titleDAuHT2Pi0 = "#font[12]{d}_{ }+_{ }Au HighTower-2";
	    TDataPoints dataDAuHT2Pi0Invyield = TDataPoints(TDrawOptions(TNamed("dAu_HT2_pi0_invyield", titleDAuHT2Pi0), TAttLine(kRed, kSolid, 2), TAttFill(kRed, 0), TAttMarker(kRed, kFullTriangleUp, 1.0), "P E Z", "P E Z", TString(titleDAuHT2Pi0)), x_HT2spectrum_dAu, y_HT2spectrum_dAu, x_err_HT2spectrum_dAu, y_err_HT2spectrum_dAu, sizeof(x_HT2spectrum_dAu)/sizeof(x_HT2spectrum_dAu[0]));
	    const Char_t *titleDAuHT2Pi0Centr = "#font[12]{d}_{ }+_{ }Au HighTower-2";
	    TDataPoints dataDAuHT2Pi0InvyieldCentr = TDataPoints(TDrawOptions(TNamed("dAu_HT2_pi0_invyield", titleDAuHT2Pi0Centr), TAttLine(kRed, kSolid, 2), TAttFill(kRed, 0), TAttMarker(kRed, kFullTriangleUp, 1.0), "PL E Z", "PL E Z", TString(titleDAuHT2Pi0Centr)), x_HT2spectrum_dAu, y_HT2spectrum_dAu, x_err_HT2spectrum_dAu, y_err_HT2spectrum_dAu, sizeof(x_HT2spectrum_dAu)/sizeof(x_HT2spectrum_dAu[0]));
	    bin_stat_list_type bins_dAu_pi0_invyield_HT2 = dataDAuHT2Pi0Invyield.getBins(dAu_HT1_HT2, dAu_HT2);

	    bin_stat_list_type bins_dAu_pi0_invyield;
	    bins_dAu_pi0_invyield.merge(bins_dAu_pi0_invyield_MB);
	    bins_dAu_pi0_invyield.merge(bins_dAu_pi0_invyield_HT1);
	    bins_dAu_pi0_invyield.merge(bins_dAu_pi0_invyield_HT2);
	    TWeightCalculator wDAuPi0Invyield("wDAuPi0Invyield", "#font[12]{d}+Au #pi^{0} invariant yield");
	    wDAuPi0Invyield.Fit(bins_dAu_pi0_invyield, "RQN LL M", "POWERLAW", 0, 100);
	    fDAuPi0Invyield = wDAuPi0Invyield.createFunc(true);
	    fDAuPi0Invyield->SetTitle("Fit to #font[12]{d}_{ }+_{ }Au");
	    fDAuPi0Invyield->SetRange(0.0, 20.0);

	    for (UInt_t i = 0;i < sizeof(x_MBspectrum_dAu)/sizeof(x_MBspectrum_dAu[0]);i++) {
		y_err_MBspectrum_dAu_smdsyst[i] *= fDAuPi0Invyield->Eval(x_MBspectrum_dAu[i]) / y_MBspectrum_dAu_smdsyst[i];
		y_MBspectrum_dAu_smdsyst[i] *= fDAuPi0Invyield->Eval(x_MBspectrum_dAu[i]) / y_MBspectrum_dAu_smdsyst[i];
	    };
	    tmp = TDataPoints(x_MBspectrum_dAu, y_MBspectrum_dAu_smdsyst, x_err_MBspectrum_dAu, y_err_MBspectrum_dAu_smdsyst, sizeof(x_MBspectrum_dAu)/sizeof(x_MBspectrum_dAu[0])).getBins();
	    const Char_t *titleDAuMBPi0_smdsyst = "Additional SMD uncertainty";
	    TDataPoints dataDAuMBPi0Invyield_smdsyst_1 = TDataPoints(TDrawOptions(TNamed("dAu_MB_pi0_invyield_smdsyst_1", titleDAuMBPi0_smdsyst), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(titleDAuMBPi0_smdsyst)), tmp, "LOWEDGE");
	    TDataPoints dataDAuMBPi0Invyield_smdsyst_2 = TDataPoints(TDrawOptions(TNamed("dAu_MB_pi0_invyield_smdsyst_2", ""), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString()), tmp, "HIGHEDGE");

	    dataPoints_pi0_invyield.push_back(dataDAuMBPi0Invyield);
	    dataPoints_pi0_dAu_invyield.push_back(dataDAuMBPi0Invyield);
	    dataPoints_pi0_invyield_centr.push_back(dataDAuMBPi0InvyieldCentr);
	    dataPoints_pi0_invyield.push_back(dataDAuHT1Pi0Invyield);
	    dataPoints_pi0_dAu_invyield.push_back(dataDAuHT1Pi0Invyield);
	    dataPoints_pi0_invyield_centr.push_back(dataDAuHT1Pi0InvyieldCentr);
	    dataPoints_pi0_invyield.push_back(dataDAuHT2Pi0Invyield);
	    dataPoints_pi0_dAu_invyield.push_back(dataDAuHT2Pi0Invyield);
	    dataPoints_pi0_invyield_centr.push_back(dataDAuHT2Pi0InvyieldCentr);
	    dataPoints_pi0_invyield.push_back(TDataPoints(TDrawOptions(TNamed("dAu_pi0_invyield_fit", fDAuPi0Invyield->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuPi0Invyield->GetTitle())), fDAuPi0Invyield));
	    dataPoints_pi0_dAu_invyield.push_back(TDataPoints(TDrawOptions(TNamed("dAu_pi0_invyield_fit", fDAuPi0Invyield->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuPi0Invyield->GetTitle())), fDAuPi0Invyield));
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAu_pi0_invyield_fit", fDAuPi0Invyield->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuPi0Invyield->GetTitle())), fDAuPi0Invyield));
	    //dataPoints_pi0_dAu_invyield.push_back(dataDAuMBPi0Invyield_smdsyst_1);
	    //dataPoints_pi0_dAu_invyield.push_back(dataDAuMBPi0Invyield_smdsyst_2);

	    bin_stat_list_type bins_dAu_pi0_crossection;
	    tmp = TDataPoints(x_MBcrossection_dAu, y_MBcrossection_dAu, x_err_MBcrossection_dAu, y_err_MBcrossection_dAu_stat, sizeof(x_MBcrossection_dAu)/sizeof(x_MBcrossection_dAu[0])).getBins(dAu_MB, dAu_MB_HT1);
	    bins_dAu_pi0_crossection.merge(tmp);
	    tmp = TDataPoints(x_HT1crossection_dAu, y_HT1crossection_dAu, x_err_HT1crossection_dAu, y_err_HT1crossection_dAu_stat, sizeof(x_HT1crossection_dAu)/sizeof(x_HT1crossection_dAu[0])).getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    bins_dAu_pi0_crossection.merge(tmp);
	    tmp = TDataPoints(x_HT2crossection_dAu, y_HT2crossection_dAu, x_err_HT2crossection_dAu, y_err_HT2crossection_dAu_stat, sizeof(x_HT2crossection_dAu)/sizeof(x_HT2crossection_dAu[0])).getBins(dAu_HT1_HT2, dAu_HT2);
	    bins_dAu_pi0_crossection.merge(tmp);
	    const Char_t *titleDAuPi0 = "#font[12]{d}_{ }+_{ }Au";
	    TDataPoints dataDAuPi0ThisCrossection = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection", titleDAuPi0), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullSquare, 1.0), "P E Z", "P E Z", TString(titleDAuPi0)), bins_dAu_pi0_crossection);
	    TDataPoints dataDAuPi0ThisCrossection_unnamed = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection", ""), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullSquare, 1.0), "P E Z", "P E Z", TString()), bins_dAu_pi0_crossection);
	    const Char_t *titleDAuPi0_1 = "#pi^{0}^{}";
	    TDataPoints dataDAuPi0ThisCrossection_name1 = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection", titleDAuPi0_1), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullSquare, 1.0), "P E Z", "P E Z", TString(titleDAuPi0_1)), bins_dAu_pi0_crossection);

	    bin_stat_list_type bins_dAu_pi0_crossection_esyst;
	    tmp = TDataPoints(x_MBcrossection_dAu, y_MBcrossection_dAu, x_err_MBcrossection_dAu_syst, y_err_MBcrossection_dAu_esyst, sizeof(x_MBcrossection_dAu)/sizeof(x_MBcrossection_dAu[0])).getBins(dAu_MB, dAu_MB_HT1);
	    bins_dAu_pi0_crossection_esyst.merge(tmp);
	    tmp = TDataPoints(x_HT1crossection_dAu, y_HT1crossection_dAu, x_err_HT1crossection_dAu_syst, y_err_HT1crossection_dAu_esyst, sizeof(x_HT1crossection_dAu)/sizeof(x_HT1crossection_dAu[0])).getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    bins_dAu_pi0_crossection_esyst.merge(tmp);
	    tmp = TDataPoints(x_HT2crossection_dAu, y_HT2crossection_dAu, x_err_HT2crossection_dAu_syst, y_err_HT2crossection_dAu_esyst, sizeof(x_HT2crossection_dAu)/sizeof(x_HT2crossection_dAu[0])).getBins(dAu_HT1_HT2, dAu_HT2);
	    bins_dAu_pi0_crossection_esyst.merge(tmp);
	    TDataPoints dataDAuPi0ThisCrossection_esyst1 = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_esyst1", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullSquare, 1.0), "E2 Z", "E2 Z", TString()), bins_dAu_pi0_crossection_esyst, ""/*"LOWEDGE"*/);
	    //TDataPoints dataDAuPi0ThisCrossection_esyst2 = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_esyst2", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullSquare, 1.0), "L", "L", TString()), bins_dAu_pi0_crossection_esyst, ""/*"HIGHEDGE"*/);
	    TDataPoints dataDAuPi0ThisCrossectionNorm = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_norm", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullSquare, 1.0), "E2 Z", "PF", TString()), x_MBcrossection_dAu_norm, y_MBcrossection_dAu_norm, x_err_MBcrossection_dAu_norm, y_err_MBcrossection_dAu_norm, sizeof(x_MBcrossection_dAu_norm)/sizeof(x_MBcrossection_dAu_norm[0]));
	    TDataPoints dataDAuPi0ThisCrossectionNorm1 = TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_norm1", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullSquare, 1.0), "E2 Z", "PF", TString()), x_MBcrossection_dAu_norm1, y_MBcrossection_dAu_norm1, x_err_MBcrossection_dAu_norm1, y_err_MBcrossection_dAu_norm1, sizeof(x_MBcrossection_dAu_norm1)/sizeof(x_MBcrossection_dAu_norm1[0]));

	    TWeightCalculator wDAuPi0Crossection("wDAuPi0Crossection", "#font[12]{d}+Au #pi^{0} cross section");
	    wDAuPi0Crossection.Fit(bins_dAu_pi0_crossection, "RQN LL M", "POWERLAW", 0, 100);
	    fDAuPi0Crossection = wDAuPi0Crossection.createFunc(true);
	    fDAuPi0Crossection->SetTitle("Fit to #font[12]{d}+Au");
	    fDAuPi0Crossection->SetRange(0.0, 20.0);
	    //fDAuPi0Crossection->Print("");
	    //TGraphErrors *pi0_dAu_cx_gr = createGraph(bins_dAu_pi0_crossection);
	    //fDAuPi0CrossectionLevy = new TF1("fDAuPi0CrossectionLevy", LevyFunc);
	    //fDAuPi0CrossectionLevy->SetTitle("Levy fit to #font[12]{d}+Au");
	    //fDAuPi0CrossectionLevy->SetRange(0.0, 20.0);
	    //fDAuPi0CrossectionLevy->SetParameter(0, 100);
	    //fDAuPi0CrossectionLevy->SetParLimits(0, 1, 50000);
	    //fDAuPi0CrossectionLevy->SetParameter(1, 1);
	    //fDAuPi0CrossectionLevy->SetParLimits(1, 0, 10);
	    //fDAuPi0CrossectionLevy->SetParameter(2, 10);
	    //fDAuPi0CrossectionLevy->SetParLimits(2, 5, 20);
	    //pi0_dAu_cx_gr->Fit(fDAuPi0CrossectionLevy, "RN LL M");

	    dataPoints_pi0_dAu_crossection.push_back(dataDAuPi0ThisCrossection);
	    dataPoints_pi0_dAu_crossection_div.push_back(dataDAuPi0ThisCrossection_esyst1);
	    //dataPoints_pi0_dAu_crossection_div.push_back(dataDAuPi0ThisCrossection_esyst2);
	    dataPoints_pi0_dAu_crossection_div_2.push_back(dataDAuPi0ThisCrossection_unnamed);
	    dataPoints_pi0_dAu_crossection_div_2.push_back(dataDAuPi0ThisCrossectionNorm);
	    //dataPoints_pi0_dAu_crossection.push_back(TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_fit", fDAuPi0Crossection->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuPi0Crossection->GetTitle())), fDAuPi0Crossection));
	    //dataPoints_pi0_dAu_crossection_div_2.push_back(TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_fit", fDAuPi0Crossection->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuPi0Crossection->GetTitle())), fDAuPi0Crossection));
	    //dataPoints_pi0_dAu_crossection.push_back(TDataPoints(TDrawOptions(TNamed("dAu_pi0_crossection_levyfit", fDAuPi0CrossectionLevy->GetTitle()), TAttLine(kRed, kDotted, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuPi0CrossectionLevy->GetTitle())), fDAuPi0CrossectionLevy));
	    dataPoints_pi0_dAu_crossection_div1.push_back(dataDAuPi0ThisCrossection_esyst1);
	    //dataPoints_pi0_dAu_crossection_div1.push_back(dataDAuPi0ThisCrossection_esyst2);
	    dataPoints_pi0_dAu_crossection_div1_2.push_back(dataDAuPi0ThisCrossection_name1);
	    dataPoints_pi0_dAu_crossection_div1_2.push_back(dataDAuPi0ThisCrossectionNorm1);

	    const Char_t *titleDAuMBEta = "#font[12]{d}_{ }+_{ }Au MinBias";
	    TDataPoints dataDAuMBEtaInvyield = TDataPoints(TDrawOptions(TNamed("dAu_MB_eta_invyield", titleDAuMBEta), TAttLine(kBlack, kSolid, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "P E Z", "P E Z", TString(titleDAuMBEta)), x_EtaMBspectrum_dAu, y_EtaMBspectrum_dAu, x_err_EtaMBspectrum_dAu, y_err_EtaMBspectrum_dAu, sizeof(x_EtaMBspectrum_dAu)/sizeof(x_EtaMBspectrum_dAu[0]));
	    dataPoints_eta_invyield.push_back(dataDAuMBEtaInvyield);
	    dataPoints_eta_dAu_invyield.push_back(dataDAuMBEtaInvyield);
	    bin_stat_list_type bins_dAu_eta_invyield_MB = dataDAuMBEtaInvyield.getBins(dAu_etatopi0_MB, dAu_etatopi0_MB_HT1);

	    const Char_t *titleDAuHT1Eta = "#font[12]{d}_{ }+_{ }Au HighTower-1";
	    TDataPoints dataDAuHT1EtaInvyield = TDataPoints(TDrawOptions(TNamed("dAu_HT1_eta_invyield", titleDAuHT1Eta), TAttLine(kBlue, kSolid, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kFullSquare, 1.0), "P E Z", "P E Z", TString(titleDAuHT1Eta)), x_EtaHT1spectrum_dAu, y_EtaHT1spectrum_dAu, x_err_EtaHT1spectrum_dAu, y_err_EtaHT1spectrum_dAu, sizeof(x_EtaHT1spectrum_dAu)/sizeof(x_EtaHT1spectrum_dAu[0]));
	    dataPoints_eta_invyield.push_back(dataDAuHT1EtaInvyield);
	    dataPoints_eta_dAu_invyield.push_back(dataDAuHT1EtaInvyield);
	    bin_stat_list_type bins_dAu_eta_invyield_HT1 = dataDAuHT1EtaInvyield.getBins(dAu_etatopi0_MB_HT1, dAu_etatopi0_HT1_HT2);

	    const Char_t *titleDAuHT2Eta = "#font[12]{d}_{ }+_{ }Au HighTower-2";
	    TDataPoints dataDAuHT2EtaInvyield = TDataPoints(TDrawOptions(TNamed("dAu_HT2_eta_invyield", titleDAuHT2Eta), TAttLine(kRed, kSolid, 2), TAttFill(kRed, 0), TAttMarker(kRed, kFullTriangleUp, 1.0), "P E Z", "P E Z", TString(titleDAuHT2Eta)), x_EtaHT2spectrum_dAu, y_EtaHT2spectrum_dAu, x_err_EtaHT2spectrum_dAu, y_err_EtaHT2spectrum_dAu, sizeof(x_EtaHT2spectrum_dAu)/sizeof(x_EtaHT2spectrum_dAu[0]));
	    dataPoints_eta_invyield.push_back(dataDAuHT2EtaInvyield);
	    dataPoints_eta_dAu_invyield.push_back(dataDAuHT2EtaInvyield);
	    bin_stat_list_type bins_dAu_eta_invyield_HT2 = dataDAuHT2EtaInvyield.getBins(dAu_etatopi0_HT1_HT2, dAu_etatopi0_HT2);

	    bin_stat_list_type bins_dAu_eta_invyield;
	    bins_dAu_eta_invyield.merge(bins_dAu_eta_invyield_MB);
	    bins_dAu_eta_invyield.merge(bins_dAu_eta_invyield_HT1);
	    bins_dAu_eta_invyield.merge(bins_dAu_eta_invyield_HT2);
	    TWeightCalculator wDAuEtaInvyield("wDAuEtaInvyield", "#font[12]{d}+Au #eta invariant yield");
	    wDAuEtaInvyield.Fit(bins_dAu_eta_invyield, "RQN LL M", "POWERLAW", 0, 100);
	    fDAuEtaInvyield = wDAuEtaInvyield.createFunc(false);
	    for (Int_t i = 2;i < fDAuEtaInvyield->GetNpar();i++) {
		fDAuEtaInvyield->FixParameter(i, fPPEtaInvyield->GetParameter(i));
	    }
	    fDAuEtaInvyield->SetParameter(0, 1);
	    fDAuEtaInvyield->FixParameter(1, 0);//fDAuEtaInvyield->GetParameter(1));
	    wDAuEtaInvyield.Fit(bins_dAu_eta_invyield, "RQN LL M", "POWERLAW", 0, 100, fDAuEtaInvyield);
	    if (fDAuEtaInvyield) delete fDAuEtaInvyield;
	    fDAuEtaInvyield = wDAuEtaInvyield.createFunc(true);
	    fDAuEtaInvyield->SetTitle("Fit to #font[12]{d}_{ }+_{ }Au");
	    fDAuEtaInvyield->SetRange(0.0, 20.0);
	    dataPoints_eta_invyield.push_back(TDataPoints(TDrawOptions(TNamed("dAu_eta_invyield_fit", fDAuEtaInvyield->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuEtaInvyield->GetTitle())), fDAuEtaInvyield));
	    dataPoints_eta_dAu_invyield.push_back(TDataPoints(TDrawOptions(TNamed("dAu_eta_invyield_fit", fDAuEtaInvyield->GetTitle()), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString(fDAuEtaInvyield->GetTitle())), fDAuEtaInvyield));

	    bin_stat_list_type bins_etatopi_dAu;
	    tmp = TDataPoints(x_EtaToPi0MB_dAu, y_EtaToPi0MB_dAu, x_err_EtaToPi0MB_dAu, y_err_EtaToPi0MB_dAu, sizeof(x_EtaToPi0MB_dAu)/sizeof(x_EtaToPi0MB_dAu[0])).getBins(dAu_etatopi0_MB, dAu_etatopi0_MB_HT1);
	    bins_etatopi_dAu.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT1_dAu, y_EtaToPi0HT1_dAu, x_err_EtaToPi0HT1_dAu, y_err_EtaToPi0HT1_dAu, sizeof(x_EtaToPi0HT1_dAu)/sizeof(x_EtaToPi0HT1_dAu[0])).getBins(dAu_etatopi0_MB_HT1, dAu_etatopi0_HT1_HT2);
	    bins_etatopi_dAu.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT2_dAu, y_EtaToPi0HT2_dAu, x_err_EtaToPi0HT2_dAu, y_err_EtaToPi0HT2_dAu, sizeof(x_EtaToPi0HT2_dAu)/sizeof(x_EtaToPi0HT2_dAu[0])).getBins(dAu_etatopi0_HT1_HT2, dAu_etatopi0_HT2);
	    bins_etatopi_dAu.merge(tmp);
	    const Char_t *titleDAuEtaToPi0 = "#font[12]{d}_{ }+_{ }Au (this analysis)";
	    TDataPoints dataDAuEtaToPi0 = TDataPoints(TDrawOptions(TNamed("dAu_etaToPi0", titleDAuEtaToPi0), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullSquare, 1.0), "P E Z", "P E Z", TString(titleDAuEtaToPi0)), bins_etatopi_dAu);
	    dataPoints_etaToPi0_dAu.push_back(dataDAuEtaToPi0);

	    bin_stat_list_type bins_etatopi_dAu_syst;
	    tmp = TDataPoints(x_EtaToPi0MB_dAu, y_EtaToPi0MB_dAu, x_err_EtaToPi0MB_dAu_syst, y_err_EtaToPi0MB_dAu_totsyst, sizeof(x_EtaToPi0MB_dAu)/sizeof(x_EtaToPi0MB_dAu[0])).getBins(dAu_etatopi0_MB, dAu_etatopi0_MB_HT1);
	    bins_etatopi_dAu_syst.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT1_dAu, y_EtaToPi0HT1_dAu, x_err_EtaToPi0HT1_dAu_syst, y_err_EtaToPi0HT1_dAu_totsyst, sizeof(x_EtaToPi0HT1_dAu)/sizeof(x_EtaToPi0HT1_dAu[0])).getBins(dAu_etatopi0_MB_HT1, dAu_etatopi0_HT1_HT2);
	    bins_etatopi_dAu_syst.merge(tmp);
	    tmp = TDataPoints(x_EtaToPi0HT2_dAu, y_EtaToPi0HT2_dAu, x_err_EtaToPi0HT2_dAu_syst, y_err_EtaToPi0HT2_dAu_totsyst, sizeof(x_EtaToPi0HT2_dAu)/sizeof(x_EtaToPi0HT2_dAu[0])).getBins(dAu_etatopi0_HT1_HT2, dAu_etatopi0_HT2);
	    bins_etatopi_dAu_syst.merge(tmp);
	    TDataPoints dataDAuEtaToPi0_syst = TDataPoints(TDrawOptions(TNamed("dAu_etaToPi0_syst", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullSquare, 1.0), "E2 Z", "E2 Z", TString()), bins_etatopi_dAu_syst);
	    dataPoints_etaToPi0_dAu1.push_back(dataDAuEtaToPi0_syst);

	    const Char_t *titleDAuMostcentralMBPi0 = "0-20% most central";
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAuMostcentral_MB_pi0_invyield", titleDAuMostcentralMBPi0), TAttLine(kBlack, kDashed, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "PL", "L", TString(titleDAuMostcentralMBPi0)), x_MBspectrum_dAuMostcentral, y_MBspectrum_dAuMostcentral, x_err_MBspectrum_dAuMostcentral, y_err_MBspectrum_dAuMostcentral, sizeof(x_MBspectrum_dAuMostcentral)/sizeof(x_MBspectrum_dAuMostcentral[0])));

	    const Char_t *titleDAuMostcentralHT1Pi0 = 0;//"#font[12]{d}+Au 0-20% HighTower-1 #pi^{0}";
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAuMostcentral_HT1_pi0_invyield", titleDAuMostcentralHT1Pi0), TAttLine(kBlue, kDashed, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kFullSquare, 1.0), "PL", "L", TString(titleDAuMostcentralHT1Pi0)), x_HT1spectrum_dAuMostcentral, y_HT1spectrum_dAuMostcentral, x_err_HT1spectrum_dAuMostcentral, y_err_HT1spectrum_dAuMostcentral, sizeof(x_HT1spectrum_dAuMostcentral)/sizeof(x_HT1spectrum_dAuMostcentral[0])));

	    const Char_t *titleDAuMostcentralHT2Pi0 = 0;//"#font[12]{d}+Au 0-20% HighTower-2 #pi^{0}";
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAuMostcentral_HT2_pi0_invyield", titleDAuMostcentralHT2Pi0), TAttLine(kRed, kDashed, 2), TAttFill(kRed, 0), TAttMarker(kRed, kFullTriangleUp, 1.0), "PL", "L", TString(titleDAuMostcentralHT2Pi0)), x_HT2spectrum_dAuMostcentral, y_HT2spectrum_dAuMostcentral, x_err_HT2spectrum_dAuMostcentral, y_err_HT2spectrum_dAuMostcentral, sizeof(x_HT2spectrum_dAuMostcentral)/sizeof(x_HT2spectrum_dAuMostcentral[0])));

	    const Char_t *titleDAuMostperipheralMBPi0 = "40-100% most peripheral";
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAuMostperipheral_MB_pi0_invyield", titleDAuMostperipheralMBPi0), TAttLine(kBlack, kDotted, 2), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "PL", "L", TString(titleDAuMostperipheralMBPi0)), x_MBspectrum_dAuMostperipheral, y_MBspectrum_dAuMostperipheral, x_err_MBspectrum_dAuMostperipheral, y_err_MBspectrum_dAuMostperipheral, sizeof(x_MBspectrum_dAuMostperipheral)/sizeof(x_MBspectrum_dAuMostperipheral[0])));

	    const Char_t *titleDAuMostperipheralHT1Pi0 = 0;//"#font[12]{d}+Au 40-100% HighTower-1 #pi^{0}";
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAuMostperipheral_HT1_pi0_invyield", titleDAuMostperipheralHT1Pi0), TAttLine(kBlue, kDotted, 2), TAttFill(kBlue, 0), TAttMarker(kBlue, kFullSquare, 1.0), "PL", "L", TString(titleDAuMostperipheralHT1Pi0)), x_HT1spectrum_dAuMostperipheral, y_HT1spectrum_dAuMostperipheral, x_err_HT1spectrum_dAuMostperipheral, y_err_HT1spectrum_dAuMostperipheral, sizeof(x_HT1spectrum_dAuMostperipheral)/sizeof(x_HT1spectrum_dAuMostperipheral[0])));

	    const Char_t *titleDAuMostperipheralHT2Pi0 = 0;//"#font[12]{d}+Au 40-100% HighTower-2 #pi^{0}";
	    dataPoints_pi0_invyield_centr.push_back(TDataPoints(TDrawOptions(TNamed("dAuMostperipheral_HT2_pi0_invyield", titleDAuMostperipheralHT2Pi0), TAttLine(kRed, kDotted, 2), TAttFill(kRed, 0), TAttMarker(kRed, kFullTriangleUp, 1.0), "PL", "L", TString(titleDAuMostperipheralHT2Pi0)), x_HT2spectrum_dAuMostperipheral, y_HT2spectrum_dAuMostperipheral, x_err_HT2spectrum_dAuMostperipheral, y_err_HT2spectrum_dAuMostperipheral, sizeof(x_HT2spectrum_dAuMostperipheral)/sizeof(x_HT2spectrum_dAuMostperipheral[0])));

	    bin_stat_list_type bins_dAu_pi0_Rcp;
	    tmp = TDataPoints(x_RcpMB_dAuRcp, y_RcpMB_dAuRcp, x_err_RcpMB_dAuRcp, y_err_RcpMB_dAuRcp, sizeof(x_RcpMB_dAuRcp)/sizeof(x_RcpMB_dAuRcp[0])).getBins(dAu_MB, dAu_MB_HT1);
	    bins_dAu_pi0_Rcp.merge(tmp);
	    tmp = TDataPoints(x_RcpHT1_dAuRcp, y_RcpHT1_dAuRcp, x_err_RcpHT1_dAuRcp, y_err_RcpHT1_dAuRcp, sizeof(x_RcpHT1_dAuRcp)/sizeof(x_RcpHT1_dAuRcp[0])).getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    bins_dAu_pi0_Rcp.merge(tmp);
	    tmp = TDataPoints(x_RcpHT2_dAuRcp, y_RcpHT2_dAuRcp, x_err_RcpHT2_dAuRcp, y_err_RcpHT2_dAuRcp, sizeof(x_RcpHT2_dAuRcp)/sizeof(x_RcpHT2_dAuRcp[0])).getBins(dAu_HT1_HT2, dAu_HT2);
	    bins_dAu_pi0_Rcp.merge(tmp);
	    const Char_t *titleDAuPi0Rcp = "#pi^{0}^{} (this analysis)";
	    TDataPoints dataDAuPi0Rcp = TDataPoints(TDrawOptions(TNamed("dAu_pi0_Rcp", titleDAuPi0Rcp), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "P E Z", "P E Z", TString(titleDAuPi0Rcp)), bins_dAu_pi0_Rcp);
	    TDataPoints dataDAuPi0RcpNorm = TDataPoints(TDrawOptions(TNamed("dAu_pi0_Rcp_norm", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullCircle, 1.0), "E2 Z", "PF", TString()), x_RcpMB_dAuRcp_norm, y_RcpMB_dAuRcp_norm, x_err_RcpMB_dAuRcp_norm, y_err_RcpMB_dAuRcp_norm, sizeof(x_RcpMB_dAuRcp_norm)/sizeof(x_RcpMB_dAuRcp_norm[0]));

	    bin_stat_list_type bins_dAu_pi0_Rcp_syst;
	    tmp = TDataPoints(x_RcpMB_dAuRcp, y_RcpMB_dAuRcp, x_err_RcpMB_dAuRcp_syst, y_err_RcpMB_dAuRcp_totsyst, sizeof(x_RcpMB_dAuRcp)/sizeof(x_RcpMB_dAuRcp[0])).getBins(dAu_MB, dAu_MB_HT1);
	    bins_dAu_pi0_Rcp_syst.merge(tmp);
	    tmp = TDataPoints(x_RcpHT1_dAuRcp, y_RcpHT1_dAuRcp, x_err_RcpHT1_dAuRcp_syst, y_err_RcpHT1_dAuRcp_totsyst, sizeof(x_RcpHT1_dAuRcp)/sizeof(x_RcpHT1_dAuRcp[0])).getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    bins_dAu_pi0_Rcp_syst.merge(tmp);
	    tmp = TDataPoints(x_RcpHT2_dAuRcp, y_RcpHT2_dAuRcp, x_err_RcpHT2_dAuRcp_syst, y_err_RcpHT2_dAuRcp_totsyst, sizeof(x_RcpHT2_dAuRcp)/sizeof(x_RcpHT2_dAuRcp[0])).getBins(dAu_HT1_HT2, dAu_HT2);
	    bins_dAu_pi0_Rcp_syst.merge(tmp);
	    TDataPoints dataDAuPi0Rcp_syst = TDataPoints(TDrawOptions(TNamed("dAu_pi0_Rcp_syst", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullCircle, 1.0), "E2 Z", "E2 Z", TString()), bins_dAu_pi0_Rcp_syst);

	    dataPoints_pi0_dAu_Rcp_2.push_back(dataDAuPi0Rcp_syst);
	    dataPoints_pi0_dAu_Rcp.push_back(dataDAuPi0Rcp);
	    dataPoints_pi0_dAu_Rcp.push_back(dataDAuPi0RcpNorm);

	    bin_stat_list_type bins_dAu_eta_Rcp;
	    tmp = TDataPoints(x_RcpMBEta_dAuRcp, y_RcpMBEta_dAuRcp, x_err_RcpMBEta_dAuRcp, y_err_RcpMBEta_dAuRcp, sizeof(x_RcpMBEta_dAuRcp)/sizeof(x_RcpMBEta_dAuRcp[0])).getBins(dAu_MB, dAu_MB_HT1);
	    bins_dAu_eta_Rcp.merge(tmp);
	    tmp = TDataPoints(x_RcpHT1Eta_dAuRcp, y_RcpHT1Eta_dAuRcp, x_err_RcpHT1Eta_dAuRcp, y_err_RcpHT1Eta_dAuRcp, sizeof(x_RcpHT1Eta_dAuRcp)/sizeof(x_RcpHT1Eta_dAuRcp[0])).getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    bins_dAu_eta_Rcp.merge(tmp);
	    tmp = TDataPoints(x_RcpHT2Eta_dAuRcp, y_RcpHT2Eta_dAuRcp, x_err_RcpHT2Eta_dAuRcp, y_err_RcpHT2Eta_dAuRcp, sizeof(x_RcpHT2Eta_dAuRcp)/sizeof(x_RcpHT2Eta_dAuRcp[0])).getBins(dAu_HT1_HT2, dAu_HT2);
	    bins_dAu_eta_Rcp.merge(tmp);
	    const Char_t *titleDAuEtaRcp = "#eta (this analysis)";
	    TDataPoints dataDAuEtaRcp = TDataPoints(TDrawOptions(TNamed("dAu_MB_eta_Rcp", titleDAuEtaRcp), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleUp, 1.5), "P E Z", "P E Z", TString(titleDAuEtaRcp)), bins_dAu_eta_Rcp);
	    TDataPoints dataDAuEtaRcpNorm = TDataPoints(TDrawOptions(TNamed("dAu_eta_Rcp_norm", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullCircle, 1.5), "E2 Z", "PF", TString()), x_RcpMBEta_dAuRcp_norm, y_RcpMBEta_dAuRcp_norm, x_err_RcpMBEta_dAuRcp_norm, y_err_RcpMBEta_dAuRcp_norm, sizeof(x_RcpMBEta_dAuRcp_norm)/sizeof(x_RcpMBEta_dAuRcp_norm[0]));

	    bin_stat_list_type bins_dAu_eta_Rcp_syst;
	    tmp = TDataPoints(x_RcpMBEta_dAuRcp, y_RcpMBEta_dAuRcp, x_err_RcpMBEta_dAuRcp_syst, y_err_RcpMBEta_dAuRcp_totsyst, sizeof(x_RcpMBEta_dAuRcp)/sizeof(x_RcpMBEta_dAuRcp[0])).getBins(dAu_MB, dAu_MB_HT1);
	    bins_dAu_eta_Rcp_syst.merge(tmp);
	    tmp = TDataPoints(x_RcpHT1Eta_dAuRcp, y_RcpHT1Eta_dAuRcp, x_err_RcpHT1Eta_dAuRcp_syst, y_err_RcpHT1Eta_dAuRcp_totsyst, sizeof(x_RcpHT1Eta_dAuRcp)/sizeof(x_RcpHT1Eta_dAuRcp[0])).getBins(dAu_MB_HT1, dAu_HT1_HT2);
	    bins_dAu_eta_Rcp_syst.merge(tmp);
	    tmp = TDataPoints(x_RcpHT2Eta_dAuRcp, y_RcpHT2Eta_dAuRcp, x_err_RcpHT2Eta_dAuRcp_syst, y_err_RcpHT2Eta_dAuRcp_totsyst, sizeof(x_RcpHT2Eta_dAuRcp)/sizeof(x_RcpHT2Eta_dAuRcp[0])).getBins(dAu_HT1_HT2, dAu_HT2);
	    bins_dAu_eta_Rcp_syst.merge(tmp);
	    TDataPoints dataDAuEtaRcp_syst = TDataPoints(TDrawOptions(TNamed("dAu_MB_eta_Rcp_syst", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullTriangleUp, 1.0), "E2 Z", "E2 Z", TString()), bins_dAu_eta_Rcp_syst);

	    dataPoints_eta_dAu_Rcp_2.push_back(dataDAuEtaRcp_syst);
	    dataPoints_eta_dAu_Rcp.push_back(dataDAuEtaRcp);
	    dataPoints_eta_dAu_Rcp.push_back(dataDAuEtaRcpNorm);

	    bin_stat_list_type bins_pi0_RdA;
	    tmp = TDataPoints(x_RDAMB_RdA, y_RDAMB_RdA, x_err_RDAMB_RdA, y_err_RDAMB_RdA, sizeof(x_RDAMB_RdA)/sizeof(x_RDAMB_RdA[0])).getBins(RdA_MB, RdA_MB_HT1);
	    bins_pi0_RdA.merge(tmp);
	    tmp = TDataPoints(x_RDAHT1_RdA, y_RDAHT1_RdA, x_err_RDAHT1_RdA, y_err_RDAHT1_RdA, sizeof(x_RDAHT1_RdA)/sizeof(x_RDAHT1_RdA[0])).getBins(RdA_MB_HT1, RdA_HT1_HT2);
	    bins_pi0_RdA.merge(tmp);
	    tmp = TDataPoints(x_RDAHT2_RdA, y_RDAHT2_RdA, x_err_RDAHT2_RdA, y_err_RDAHT2_RdA, sizeof(x_RDAHT2_RdA)/sizeof(x_RDAHT2_RdA[0])).getBins(RdA_HT1_HT2, RdA_HT2);
	    bins_pi0_RdA.merge(tmp);
	    const Char_t *titlePi0RdA = "#pi^{0}^{} (this analysis)";
	    TDataPoints dataPi0ThisRdA = TDataPoints(TDrawOptions(TNamed("pi0_RdA", titlePi0RdA), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "P E Z", "P E Z", TString(titlePi0RdA)), bins_pi0_RdA);

	    bin_stat_list_type bins_pi0_RdA_esyst;
	    tmp = TDataPoints(x_RDAMB_RdA, y_RDAMB_RdA, x_err_RDAMB_RdA_syst, y_err_RDAMB_RdA_esyst, sizeof(x_RDAMB_RdA)/sizeof(x_RDAMB_RdA[0])).getBins(RdA_MB, RdA_MB_HT1);
	    bins_pi0_RdA_esyst.merge(tmp);
	    tmp = TDataPoints(x_RDAHT1_RdA, y_RDAHT1_RdA, x_err_RDAHT1_RdA_syst, y_err_RDAHT1_RdA_esyst, sizeof(x_RDAHT1_RdA)/sizeof(x_RDAHT1_RdA[0])).getBins(RdA_MB_HT1, RdA_HT1_HT2);
	    bins_pi0_RdA_esyst.merge(tmp);
	    tmp = TDataPoints(x_RDAHT2_RdA, y_RDAHT2_RdA, x_err_RDAHT2_RdA_syst, y_err_RDAHT2_RdA_esyst, sizeof(x_RDAHT2_RdA)/sizeof(x_RDAHT2_RdA[0])).getBins(RdA_HT1_HT2, RdA_HT2);
	    bins_pi0_RdA_esyst.merge(tmp);
	    TDataPoints dataPi0ThisRdA_esyst1 = TDataPoints(TDrawOptions(TNamed("pi0_RdA_esyst1", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullCircle, 1.0), "E2 Z", "E2 Z", TString()), bins_pi0_RdA_esyst, ""/*"LOWEDGE"*/);
	    //TDataPoints dataPi0ThisRdA_esyst2 = TDataPoints(TDrawOptions(TNamed("pi0_RdA_esyst2", ""), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullCircle, 1.0), "L", "L", TString()), bins_pi0_RdA_esyst, ""/*"HIGHEDGE"*/);
	    TDataPoints dataPi0ThisRdANorm = TDataPoints(TDrawOptions(TNamed("pi0_RdA_norm", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullCircle, 1.0), "E2 Z", "PF", TString()), x_RDAMB_RdA_norm, y_RDAMB_RdA_norm, x_err_RDAMB_RdA_norm, y_err_RDAMB_RdA_norm, sizeof(x_RDAMB_RdA_norm)/sizeof(x_RDAMB_RdA_norm[0]));

	    dataPoints_pi0_RdA_2.push_back(dataPi0ThisRdA_esyst1);
	    //dataPoints_pi0_RdA.push_back(dataPi0ThisRdA_esyst2);
	    dataPoints_pi0_RdA.push_back(dataPi0ThisRdA);
	    dataPoints_pi0_RdA.push_back(dataPi0ThisRdANorm);

	    bin_stat_list_type bins_eta_RdA;
	    tmp = TDataPoints(x_RDAMBEta_RdA, y_RDAMBEta_RdA, x_err_RDAMBEta_RdA, y_err_RDAMBEta_RdA, sizeof(x_RDAMBEta_RdA)/sizeof(x_RDAMBEta_RdA[0])).getBins(RdA_eta_MB, RdA_eta_MB_HT1);
	    bins_eta_RdA.merge(tmp);
	    tmp = TDataPoints(x_RDAHT1Eta_RdA, y_RDAHT1Eta_RdA, x_err_RDAHT1Eta_RdA, y_err_RDAHT1Eta_RdA, sizeof(x_RDAHT1Eta_RdA)/sizeof(x_RDAHT1Eta_RdA[0])).getBins(RdA_eta_MB_HT1, RdA_eta_HT1_HT2);
	    bins_eta_RdA.merge(tmp);
	    tmp = TDataPoints(x_RDAHT2Eta_RdA, y_RDAHT2Eta_RdA, x_err_RDAHT2Eta_RdA, y_err_RDAHT2Eta_RdA, sizeof(x_RDAHT2Eta_RdA)/sizeof(x_RDAHT2Eta_RdA[0])).getBins(RdA_eta_HT1_HT2, RdA_eta_HT2);
	    bins_eta_RdA.merge(tmp);
	    const Char_t *titleEtaRdA = "#eta (this analysis)";
	    TDataPoints dataEtaThisRdA = TDataPoints(TDrawOptions(TNamed("eta_RdA", titleEtaRdA), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleUp, 1.0), "P E Z", "P E Z", TString(titleEtaRdA)), bins_eta_RdA);

	    bin_stat_list_type bins_eta_RdA_esyst;
	    tmp = TDataPoints(x_RDAMBEta_RdA, y_RDAMBEta_RdA, x_err_RDAMBEta_RdA_syst, y_err_RDAMBEta_RdA_esyst, sizeof(x_RDAMBEta_RdA)/sizeof(x_RDAMBEta_RdA[0])).getBins(RdA_eta_MB, RdA_eta_MB_HT1);
	    bins_eta_RdA_esyst.merge(tmp);
	    tmp = TDataPoints(x_RDAHT1Eta_RdA, y_RDAHT1Eta_RdA, x_err_RDAHT1Eta_RdA_syst, y_err_RDAHT1Eta_RdA_esyst, sizeof(x_RDAHT1Eta_RdA)/sizeof(x_RDAHT1Eta_RdA[0])).getBins(RdA_eta_MB_HT1, RdA_eta_HT1_HT2);
	    bins_eta_RdA_esyst.merge(tmp);
	    tmp = TDataPoints(x_RDAHT2Eta_RdA, y_RDAHT2Eta_RdA, x_err_RDAHT2Eta_RdA_syst, y_err_RDAHT2Eta_RdA_esyst, sizeof(x_RDAHT2Eta_RdA)/sizeof(x_RDAHT2Eta_RdA[0])).getBins(RdA_eta_HT1_HT2, RdA_eta_HT2);
	    bins_eta_RdA_esyst.merge(tmp);
	    TDataPoints dataEtaThisRdA_esyst1 = TDataPoints(TDrawOptions(TNamed("eta_RdA_esyst1", ""), TAttLine(kBlack, kSolid, 1), TAttFill(17, 1001), TAttMarker(kBlack, kFullTriangleUp, 1.0), "E2 Z", "E2 Z", TString()), bins_eta_RdA_esyst, ""/*"LOWEDGE"*/);
	    //TDataPoints dataEtaThisRdA_esyst2 = TDataPoints(TDrawOptions(TNamed("eta_RdA_esyst2", ""), TAttLine(kBlack, kSolid, 1), TAttFill(kBlack, 0), TAttMarker(kBlack, kFullTriangleUp, 1.0), "L", "L", TString()), bins_eta_RdA_esyst, ""/*"HIGHEDGE"*/);
	    TDataPoints dataEtaThisRdANorm = TDataPoints(TDrawOptions(TNamed("eta_RdA_norm", ""), TAttLine(13, kSolid, 1), TAttFill(13, 1001), TAttMarker(13, kFullTriangleUp, 1.0), "E2 Z", "PF", TString()), x_RDAMBEta_RdA_norm, y_RDAMBEta_RdA_norm, x_err_RDAMBEta_RdA_norm, y_err_RDAMBEta_RdA_norm, sizeof(x_RDAMBEta_RdA_norm)/sizeof(x_RDAMBEta_RdA_norm[0]));

	    dataPoints_eta_RdA_2.push_back(dataEtaThisRdA_esyst1);
	    //dataPoints_eta_RdA.push_back(dataEtaThisRdA_esyst2);
	    dataPoints_eta_RdA.push_back(dataEtaThisRdA);
	    dataPoints_eta_RdA.push_back(dataEtaThisRdANorm);
	}

	allSettings.settingsDAuNoCentral.dataPointsInvYield = dataPoints_dAu_invYield;
	allSettings.settingsDAuNoCentral.dataPointsCrossSection = dataPoints_dAu_crossSection;
	allSettings.settingsDAuNoCentral.dataPointsEtaToPi = dataPoints_etatopi;
	allSettings.settingsDAuAllCentral.dataPointsInvYield = dataPoints_dAu_invYield;
	allSettings.settingsDAuAllCentral.dataPointsCrossSection = dataPoints_dAu_crossSection;
	allSettings.settingsDAuAllCentral.dataPointsEtaToPi = dataPoints_etatopi;
	allSettings.settingsDAuMostCentral.dataPointsInvYield = dataPoints_dAu_invYield;
	allSettings.settingsDAuMostCentral.dataPointsCrossSection = dataPoints_dAu_crossSection;
	allSettings.settingsDAuMostCentral.dataPointsEtaToPi = dataPoints_etatopi;
	allSettings.settingsDAuMidCentral.dataPointsInvYield = dataPoints_dAu_invYield;
	allSettings.settingsDAuMidCentral.dataPointsCrossSection = dataPoints_dAu_crossSection;
	allSettings.settingsDAuMidCentral.dataPointsEtaToPi = dataPoints_etatopi;
	allSettings.settingsDAuMostPeripheral.dataPointsInvYield = dataPoints_dAu_invYield;
	allSettings.settingsDAuMostPeripheral.dataPointsCrossSection = dataPoints_dAu_crossSection;
	allSettings.settingsDAuMostPeripheral.dataPointsEtaToPi = dataPoints_etatopi;
	allSettings.settingsPP.dataPointsInvYield = dataPoints_pp_invYield;
	allSettings.settingsPP.dataPointsCrossSection = dataPoints_pp_crossSection;
	allSettings.settingsPP.dataPointsEtaToPi = dataPoints_etatopi;
	allSettings.dataPointsR = dataPoints_R;
	
	TAllResults *p_allResults = new TAllResults();
	TAllResults &allResults = *p_allResults;

	TString dataDir = "./data_output/";
	dataDir += DATA_DIR;

	if (true) {
	    cout << "Starting analysis..." << endl;
	    showAnalysis(dataDir.Data(), allSettings, allResults);
	    cout << "Finished analysis." << endl;
	}

	if (false) saveCanvases("./plots_save", false);
	
	if (true) {
	    // Invariant yield plots
	    
	    TString namePi0InvYield = "pi0_invyield";
	    TString titlePi0InvYield = "#pi^{0} invariant yield";
    	    TCanvas *cPi0InvYield = new TCanvas(namePi0InvYield.Data(), titlePi0InvYield.Data(), 0, 0, 700, 1000);
	    cPi0InvYield->Divide(1, 3, 0, 0);
	    cPi0InvYield->cd(1);
	    gPad->SetPad(0, 0.5, 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString hPi0InvYieldName = namePi0InvYield; hPi0InvYieldName += "_hist";
	    TString hPi0InvYieldTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];#frac{1}{_{ }2#pi#font[12]{p}_{#font[12]{T}}} #frac{_{ }d^{2}#font[12]{N}}{_{ }d#font[12]{p}_{#font[12]{T}}_{ }d#font[12]{y}_{ }}   [(GeV/#font[12]{c})^{-2}^{ }]";
	    TH1F *hPi0InvYield = new TH1F(hPi0InvYieldName.Data(), hPi0InvYieldTitle.Data(), 1000, 0, 17);
	    hPi0InvYield->GetXaxis()->SetTitleOffset(1.0);
	    hPi0InvYield->GetYaxis()->SetTitleOffset(1.23);
	    setHistFontSize(hPi0InvYield, 240.0, 8.0);
	    hPi0InvYield->Draw();
	    TLegend *lPi0InvYield = new TLegend(0.48, 0.45, 0.96, 0.93);
	    lPi0InvYield->SetLineColor(0);
	    lPi0InvYield->SetFillColor(0);
	    SetShadowColor(lPi0InvYield, 0);
	    setLegendFontSize(lPi0InvYield, 240.0, 8.0);
	    lPi0InvYield->Draw();
	    showResultsOthers(hPi0InvYield, lPi0InvYield, dataPoints_pi0_invyield);
	    hPi0InvYield->SetMinimum(1.001e-11);
	    hPi0InvYield->SetMaximum(0.9e+00);
	    gPad->SetLogy();
	    TLatex *label = new TLatex(4.0, 0.177, "#pi^{0}^{}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();

	    if (false) {
	    TString namePi0InvYieldCentr = "pi0_invyield_centr";
	    TString titlePi0InvYieldCentr = "#pi^{0} invariant yield in centrality bins";
	    TString hPi0InvYieldCentrName = namePi0InvYieldCentr; hPi0InvYieldCentrName += "_hist";
	    TString hPi0InvYieldCentrTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];Invariant yield / Fit";
	    TH1F *hPi0InvYieldCentr = new TH1F(hPi0InvYieldCentrName.Data(), hPi0InvYieldCentrTitle.Data(), 1000, 0, 17);
	    hPi0InvYieldCentr->GetXaxis()->SetTitleOffset(1.0);
	    hPi0InvYieldCentr->GetYaxis()->SetTitleOffset(1.3);
	    setHistFontSize(hPi0InvYieldCentr, 240.0, 8.0);
	    hPi0InvYieldCentr->Draw();
	    TLegend *lPi0InvYieldCentr = new TLegend(0.2, 0.65, 0.6, 0.9);
	    lPi0InvYieldCentr->SetLineColor(0);
	    lPi0InvYieldCentr->SetFillColor(0);
	    setLegendFontSize(lPi0InvYieldCentr, 240.0, 8.0);
	    lPi0InvYieldCentr->Draw();
	    data_points_list dataPointsCentrDiv = dataPoints_pi0_invyield_centr;
	    if (fDAuPi0Invyield) dataPointsCentrDiv /= *fDAuPi0Invyield;
	    showResultsOthers(hPi0InvYieldCentr, lPi0InvYieldCentr, dataPointsCentrDiv);
	    hPi0InvYieldCentr->SetMinimum(0);
	    hPi0InvYieldCentr->SetMaximum(6);
	    }

	    TString namePi0InvYieldDivPP = "pi0_invyield_pp_div";
	    TString titlePi0InvYieldDivPP = "p+p #pi^{0} invariant yield / fit";
	    cPi0InvYield->cd(2);
	    gPad->SetPad(0, 0.5*4.0/7.0, 1, 0.5);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0);
	    TString hPi0InvYieldPPDivName = namePi0InvYieldDivPP; hPi0InvYieldPPDivName += "_hist";
	    TString hPi0InvYieldPPDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];Data_{ }/_{ }fit        ";
	    TH1F *hPi0InvYieldPPDiv = new TH1F(hPi0InvYieldPPDivName.Data(), hPi0InvYieldPPDivTitle.Data(), 1000, 0, 17);
	    setHistFontSize(hPi0InvYieldPPDiv, 240.0, 8.0);
	    hPi0InvYieldPPDiv->GetXaxis()->SetTitleOffset(1.0);
	    hPi0InvYieldPPDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0InvYieldPPDiv->GetYaxis()->SetTitleOffset(0.54);
	    hPi0InvYieldPPDiv->GetYaxis()->SetNdivisions(505);
	    hPi0InvYieldPPDiv->Draw();
	    TLegend *lPi0InvYieldPPDiv = new TLegend(0.2, 0.65, 0.6, 0.9);
	    lPi0InvYieldPPDiv->SetLineColor(0);
	    lPi0InvYieldPPDiv->SetFillColor(0);
	    SetShadowColor(lPi0InvYieldPPDiv, 0);
	    setLegendFontSize(lPi0InvYieldPPDiv, 240.0, 8.0);
	    lPi0InvYieldPPDiv->Draw();
	    data_points_list dataPointsPPPi0Div = dataPoints_pi0_pp_invyield;
	    if (fPPPi0Invyield) dataPointsPPPi0Div /= *fPPPi0Invyield;
	    showResultsOthers(hPi0InvYieldPPDiv, lPi0InvYieldPPDiv, dataPointsPPPi0Div);
	    hPi0InvYieldPPDiv->SetMinimum(0.001);
	    hPi0InvYieldPPDiv->SetMaximum(5.0-0.001);
	    delete lPi0InvYieldPPDiv; lPi0InvYieldPPDiv = 0;
	    TLatex *labelPP = new TLatex(4.0, 3.5, "#font[12]{p}_{ }+_{ }#font[12]{p}");
	    setLatexFontSize(labelPP, 240.0, 8.0);
	    labelPP->Draw();

	    TString namePi0InvYieldDivDAu = "pi0_invyield_dAu_div";
	    TString titlePi0InvYieldDivDAu = "d+Au #pi^{0} invariant yield / fit";
	    cPi0InvYield->cd(3);
	    gPad->SetPad(0, 0, 1, 0.5*4.0/7.0);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.25);
	    TString hPi0InvYieldDAuDivName = namePi0InvYieldDivDAu; hPi0InvYieldDAuDivName += "_hist";
	    TString hPi0InvYieldDAuDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];Data_{ }/_{ }fit        ";
	    TH1F *hPi0InvYieldDAuDiv = new TH1F(hPi0InvYieldDAuDivName.Data(), hPi0InvYieldDAuDivTitle.Data(), 1000, 0, 17);
	    hPi0InvYieldDAuDiv->GetXaxis()->SetTitleOffset(1.05);
	    hPi0InvYieldDAuDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0InvYieldDAuDiv->GetYaxis()->SetTitleOffset(0.72);
	    hPi0InvYieldDAuDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0InvYieldDAuDiv, 240.0, 8.0);
	    hPi0InvYieldDAuDiv->Draw();
	    TLegend *lPi0InvYieldDAuDiv = new TLegend(0.2, 0.65, 0.6, 0.9);
	    lPi0InvYieldDAuDiv->SetLineColor(0);
	    lPi0InvYieldDAuDiv->SetFillColor(0);
	    SetShadowColor(lPi0InvYieldDAuDiv, 0);
	    setLegendFontSize(lPi0InvYieldDAuDiv, 240.0, 8.0);
	    lPi0InvYieldDAuDiv->Draw();
	    data_points_list dataPointsDAuPi0Div = dataPoints_pi0_dAu_invyield;
	    if (fDAuPi0Invyield) dataPointsDAuPi0Div /= *fDAuPi0Invyield;
	    showResultsOthers(hPi0InvYieldDAuDiv, lPi0InvYieldDAuDiv, dataPointsDAuPi0Div);
	    hPi0InvYieldDAuDiv->SetMinimum(0.001);
	    hPi0InvYieldDAuDiv->SetMaximum(5.0-0.001);
	    delete lPi0InvYieldDAuDiv; lPi0InvYieldDAuDiv = 0;
	    TLatex *labelDAu = new TLatex(4.0, 3.5, "#font[12]{d}_{ }+_{ }Au");
	    setLatexFontSize(labelDAu, 240.0, 8.0);
	    labelDAu->Draw();

	    TString nameEtaInvYield = "eta_invyield";
	    TString titleEtaInvYield = "#eta invariant yield";
    	    TCanvas *cEtaInvYield = new TCanvas(nameEtaInvYield.Data(), titleEtaInvYield.Data(), 0, 0, 700, 1000);
	    cEtaInvYield->Divide(1, 3, 0, 0);
	    cEtaInvYield->cd(1);
	    gPad->SetPad(0, 0.5, 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString hEtaInvYieldName = nameEtaInvYield; hEtaInvYieldName += "_hist";
	    TString hEtaInvYieldTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];#frac{1}{_{ }2#pi#font[12]{p}_{#font[12]{T}}} #frac{_{ }d^{2}#font[12]{N}}{_{ }d#font[12]{p}_{#font[12]{T}}_{ }d#font[12]{y}_{ }}   [(GeV/#font[12]{c})^{-2}^{ }]";
	    TH1F *hEtaInvYield = new TH1F(hEtaInvYieldName.Data(), hEtaInvYieldTitle.Data(), 1000, 0, 17);
	    hEtaInvYield->GetXaxis()->SetTitleOffset(1.0);
	    hEtaInvYield->GetYaxis()->SetTitleOffset(1.23);
	    setHistFontSize(hEtaInvYield, 240.0, 8.0);
	    hEtaInvYield->Draw();
	    TLegend *lEtaInvYield = new TLegend(0.48, 0.45, 0.96, 0.93);
	    lEtaInvYield->SetLineColor(0);
	    lEtaInvYield->SetFillColor(0);
	    SetShadowColor(lEtaInvYield, 0);
	    setLegendFontSize(lEtaInvYield, 240.0, 8.0);
	    lEtaInvYield->Draw();
	    showResultsOthers(hEtaInvYield, lEtaInvYield, dataPoints_eta_invyield);
	    hEtaInvYield->SetMinimum(1.001e-11);
	    hEtaInvYield->SetMaximum(0.9e+00);
	    gPad->SetLogy();
	    TLatex *labelEta = new TLatex(4.0, 0.177, "#eta");
	    setLatexFontSize(labelEta, 240.0, 8.0);
	    labelEta->Draw();

	    TString nameEtaInvYieldDivPP = "eta_invyield_pp_div";
	    TString titleEtaInvYieldDivPP = "p+p #eta invariant yield / fit";
	    cEtaInvYield->cd(2);
	    gPad->SetPad(0, 0.5*4.0/7.0, 1, 0.5);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0);
	    TString hEtaInvYieldPPDivName = nameEtaInvYieldDivPP; hEtaInvYieldPPDivName += "_hist";
	    TString hEtaInvYieldPPDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];Data_{ }/_{ }fit        ";
	    TH1F *hEtaInvYieldPPDiv = new TH1F(hEtaInvYieldPPDivName.Data(), hEtaInvYieldPPDivTitle.Data(), 1000, 0, 17);
	    setHistFontSize(hEtaInvYieldPPDiv, 240.0, 8.0);
	    hEtaInvYieldPPDiv->GetXaxis()->SetTitleOffset(1.0);
	    hEtaInvYieldPPDiv->GetXaxis()->SetTickLength(0.06);
	    hEtaInvYieldPPDiv->GetYaxis()->SetTitleOffset(0.54);
	    hEtaInvYieldPPDiv->GetYaxis()->SetNdivisions(505);
	    hEtaInvYieldPPDiv->Draw();
	    TLegend *lEtaInvYieldPPDiv = new TLegend(0.2, 0.65, 0.6, 0.9);
	    lEtaInvYieldPPDiv->SetLineColor(0);
	    lEtaInvYieldPPDiv->SetFillColor(0);
	    SetShadowColor(lEtaInvYieldPPDiv, 0);
	    setLegendFontSize(lEtaInvYieldPPDiv, 240.0, 8.0);
	    lEtaInvYieldPPDiv->Draw();
	    data_points_list dataPointsPPEtaDiv = dataPoints_eta_pp_invyield;
	    if (fPPEtaInvyield) dataPointsPPEtaDiv /= *fPPEtaInvyield;
	    showResultsOthers(hEtaInvYieldPPDiv, lEtaInvYieldPPDiv, dataPointsPPEtaDiv);
	    hEtaInvYieldPPDiv->SetMinimum(0.001);
	    hEtaInvYieldPPDiv->SetMaximum(5.0-0.001);
	    delete lEtaInvYieldPPDiv; lEtaInvYieldPPDiv = 0;
	    TLatex *labelPPEta = new TLatex(4.0, 3.5, "#font[12]{p}_{ }+_{ }#font[12]{p}");
	    setLatexFontSize(labelPPEta, 240.0, 8.0);
	    labelPPEta->Draw();

	    TString nameEtaInvYieldDivDAu = "eta_invyield_dAu_div";
	    TString titleEtaInvYieldDivDAu = "d+Au #eta invariant yield / fit";
	    cEtaInvYield->cd(3);
	    gPad->SetPad(0, 0, 1, 0.5*4.0/7.0);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.25);
	    TString hEtaInvYieldDAuDivName = nameEtaInvYieldDivDAu; hEtaInvYieldDAuDivName += "_hist";
	    TString hEtaInvYieldDAuDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];Data_{ }/_{ }fit        ";
	    TH1F *hEtaInvYieldDAuDiv = new TH1F(hEtaInvYieldDAuDivName.Data(), hEtaInvYieldDAuDivTitle.Data(), 1000, 0, 17);
	    hEtaInvYieldDAuDiv->GetXaxis()->SetTitleOffset(1.05);
	    hEtaInvYieldDAuDiv->GetXaxis()->SetTickLength(0.06);
	    hEtaInvYieldDAuDiv->GetYaxis()->SetTitleOffset(0.72);
	    hEtaInvYieldDAuDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hEtaInvYieldDAuDiv, 240.0, 8.0);
	    hEtaInvYieldDAuDiv->Draw();
	    TLegend *lEtaInvYieldDAuDiv = new TLegend(0.2, 0.65, 0.6, 0.9);
	    lEtaInvYieldDAuDiv->SetLineColor(0);
	    lEtaInvYieldDAuDiv->SetFillColor(0);
	    SetShadowColor(lEtaInvYieldDAuDiv, 0);
	    setLegendFontSize(lEtaInvYieldDAuDiv, 240.0, 8.0);
	    lEtaInvYieldDAuDiv->Draw();
	    data_points_list dataPointsDAuEtaDiv = dataPoints_eta_dAu_invyield;
	    if (fDAuEtaInvyield) dataPointsDAuEtaDiv /= *fDAuEtaInvyield;
	    showResultsOthers(hEtaInvYieldDAuDiv, lEtaInvYieldDAuDiv, dataPointsDAuEtaDiv);
	    hEtaInvYieldDAuDiv->SetMinimum(0.001);
	    hEtaInvYieldDAuDiv->SetMaximum(5.0-0.001);
	    delete lEtaInvYieldDAuDiv; lEtaInvYieldDAuDiv = 0;
	    TLatex *labelDAuEta = new TLatex(4.0, 3.5, "#font[12]{d}_{ }+_{ }Au");
	    setLatexFontSize(labelDAuEta, 240.0, 8.0);
	    labelDAuEta->Draw();
	}

	if (true) {
	    // p+p and d+Au cross section plots, comparison to theory
	    
	    TString namePi0CrossSection = "pp_dAu_pi0_crossection";
	    TString titlePi0CrossSection = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p} and #font[12]{d}+Au";
    	    TCanvas *cPi0CrossSection = new TCanvas(namePi0CrossSection.Data(), titlePi0CrossSection.Data(), 0, 0, 700, 1000);
	    cPi0CrossSection->Divide(1, 3, 0, 0);
	    cPi0CrossSection->cd(1);
	    gPad->SetPad(0, 0.5, 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString hPi0CrossSectionName = namePi0CrossSection; hPi0CrossSectionName += "_hist";
	    TString hPi0CrossSectionTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];#font[12]{E} #frac{d^{3}#sigma}{d#font[22]{p}^{3}}   [mb GeV^{_{ }-2}_{ }#font[12]{c}^{3_{#color[0]{-}}}]";
	    TH1F *hPi0CrossSection = new TH1F(hPi0CrossSectionName.Data(), hPi0CrossSectionTitle.Data(), 1000, 0, 20);
	    hPi0CrossSection->GetXaxis()->SetTitleOffset(1.0);
	    hPi0CrossSection->GetYaxis()->SetTitleOffset(1.23);
	    setHistFontSize(hPi0CrossSection, 240.0, 8.0);
	    hPi0CrossSection->Draw();
	    TLegend *lPi0CrossSection = new TLegend(0.511494, 0.60, 0.933908, 0.80);
	    lPi0CrossSection->SetLineColor(0);
	    lPi0CrossSection->SetFillColor(0);
	    SetShadowColor(lPi0CrossSection, 0);
	    setLegendFontSize(lPi0CrossSection, 240.0, 8.0);
	    lPi0CrossSection->Draw();
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pi0_pp_crossection);
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pp_crossSection);
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pi0_dAu_crossection);
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_dAu_crossSection);
	    data_points_list dataPoints_pi0_pp_crossection_noname = dataPoints_pi0_pp_crossection;
	    for (data_points_list::iterator iter = dataPoints_pi0_pp_crossection_noname.begin();iter != dataPoints_pi0_pp_crossection_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pi0_pp_crossection_noname);
	    data_points_list dataPoints_pi0_dAu_crossection_noname = dataPoints_pi0_dAu_crossection;
	    for (data_points_list::iterator iter = dataPoints_pi0_dAu_crossection_noname.begin();iter != dataPoints_pi0_dAu_crossection_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pi0_dAu_crossection_noname);
	    hPi0CrossSection->SetMinimum(2e-10);
	    hPi0CrossSection->SetMaximum(0.9e+03);
	    hPi0CrossSection->GetXaxis()->SetRangeUser(0, 17);
	    gPad->SetLogy();
	    TLatex *label = new TLatex(10.2, 20.0, "#pi^{0}^{}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();

	    {
	    TString namePi0CrossSectionDiv = "pp_pi0_crossection_div";
	    TString titlePi0CrossSectionDiv = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p} / theory";
	    cPi0CrossSection->cd(2);
	    gPad->SetPad(0, 0.5*(4.0/7.0), 1, 0.5);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0);
	    TString hPi0CrossSectionDivName = namePi0CrossSectionDiv; hPi0CrossSectionDivName += "_hist";
	    TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];Data_{ }/_{ }theory    ";
	    TH1F *hPi0CrossSectionDiv = new TH1F(hPi0CrossSectionDivName.Data(), hPi0CrossSectionDivTitle.Data(), 1000, 0, 20);
	    hPi0CrossSectionDiv->GetXaxis()->SetTitleOffset(1.0);
	    hPi0CrossSectionDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0CrossSectionDiv->GetYaxis()->SetTitleOffset(0.54);
	    hPi0CrossSectionDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0CrossSectionDiv, 240.0, 8.0);
	    hPi0CrossSectionDiv->Draw();
	    TLegend *lPi0CrossSectionDiv = new TLegend(0.228448, 0.855967, 0.79023, 0.95679);
	    lPi0CrossSectionDiv->SetLineColor(0);
	    lPi0CrossSectionDiv->SetFillColor(0);
	    SetShadowColor(lPi0CrossSectionDiv, 0);
	    setLegendFontSize(lPi0CrossSectionDiv, 240.0, 8.0);
	    lPi0CrossSectionDiv->Draw();
	    data_points_list dataPointsDiv = dataPoints_pi0_pp_crossection_div;
	    if (fDIVPP) dataPointsDiv /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    data_points_list dataPointsOthersDiv = dataPoints_pp_crossSection_div;
	    if (fDIVPP) dataPointsOthersDiv /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv);
	    data_points_list dataPointsDiv_2 = dataPoints_pi0_pp_crossection_div_2;
	    if (fDIVPP) dataPointsDiv_2 /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    data_points_list dataPointsDiv_noname = dataPointsDiv_2;
	    for (data_points_list::iterator iter = dataPointsDiv_noname.begin();iter != dataPointsDiv_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_noname);
	    TLatex *label = new TLatex(7.55, 2.0, "#font[12]{p}_{ }+_{ }#font[12]{p}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    hPi0CrossSectionDiv->SetMinimum(0.001);
	    hPi0CrossSectionDiv->SetMaximum(3-0.001);
	    hPi0CrossSectionDiv->GetXaxis()->SetRangeUser(0, 17);
	    }

	    {
	    TString namePi0CrossSectionDiv = "dAu_pi0_crossection_div";
	    TString titlePi0CrossSectionDiv = "#pi^{0} invariant cross section in #font[12]{d}+Au / theory";
	    cPi0CrossSection->cd(3);
	    gPad->SetPad(0, 0, 1, 0.5*(4.0/7.0));
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.25);
	    TString hPi0CrossSectionDivName = namePi0CrossSectionDiv; hPi0CrossSectionDivName += "_hist";
	    TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];Data_{ }/_{ }theory    ";
	    TH1F *hPi0CrossSectionDiv = new TH1F(hPi0CrossSectionDivName.Data(), hPi0CrossSectionDivTitle.Data(), 1000, 0, 20);
	    hPi0CrossSectionDiv->GetXaxis()->SetTitleOffset(1.05);
	    hPi0CrossSectionDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0CrossSectionDiv->GetYaxis()->SetTitleOffset(0.72);
	    hPi0CrossSectionDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0CrossSectionDiv, 240.0, 8.0);
	    hPi0CrossSectionDiv->Draw();
	    TLegend *lPi0CrossSectionDiv = new TLegend(0.5, 0.6, 0.85, 0.85);
	    lPi0CrossSectionDiv->SetLineColor(0);
	    lPi0CrossSectionDiv->SetFillColor(0);
	    SetShadowColor(lPi0CrossSectionDiv, 0);
	    setLegendFontSize(lPi0CrossSectionDiv, 240.0, 8.0);
	    lPi0CrossSectionDiv->Draw();
	    data_points_list dataPointsDiv = dataPoints_pi0_dAu_crossection_div;
	    if (fDIV) dataPointsDiv /= *fDIV;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    data_points_list dataPointsOthersDiv = dataPoints_dAu_crossSection_div;
	    if (fDIV) dataPointsOthersDiv /= *fDIV;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv);
	    data_points_list dataPointsDiv_2 = dataPoints_pi0_dAu_crossection_div_2;
	    if (fDIV) dataPointsDiv_2 /= *fDIV;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    data_points_list dataPointsDiv_noname = dataPointsDiv_2;
	    for (data_points_list::iterator iter = dataPointsDiv_noname.begin();iter != dataPointsDiv_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_noname);
	    TLatex *label = new TLatex(7.55, 2.0, "#font[12]{d}_{ }+_{ }Au");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    hPi0CrossSectionDiv->SetMinimum(0.001);
	    hPi0CrossSectionDiv->SetMaximum(3-0.001);
	    hPi0CrossSectionDiv->GetXaxis()->SetRangeUser(0, 17);
	    delete lPi0CrossSectionDiv;
	    }
	}

	if (true) {
	    // p+p cross section for the spin paper
	    
	    TString namePi0CrossSection = "pp_pi0_crossection";
	    TString titlePi0CrossSection = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p}";
    	    TCanvas *cPi0CrossSection = new TCanvas(namePi0CrossSection.Data(), titlePi0CrossSection.Data(), 0, 0, 700, 700);
	    cPi0CrossSection->Divide(1, 2, 0, 0);
	    cPi0CrossSection->cd(1);
	    gPad->SetPad(0, 0.4, 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString hPi0CrossSectionName = namePi0CrossSection; hPi0CrossSectionName += "_hist";
	    TString hPi0CrossSectionTitle = ";#font[12]{p}_{#shiftx[-0.3]{#font[12]{T}}} #shifty[-0.07]{[}Ge#shiftx[-0.12]{V}#shiftx[-0.4]{#shifty[-0.06]{/}}#shiftx[-0.3]{#font[12]{c}}#shifty[-0.07]{]};#font[12]{E} #frac{d^{3}#sigma}{d#font[12]{p}^{3}}   [mb GeV^{_{ }-2}_{ }#font[12]{c}^{3_{#color[0]{-}}}]";
	    TH1F *hPi0CrossSection = new TH1F(hPi0CrossSectionName.Data(), hPi0CrossSectionTitle.Data(), 1000, 0, 20);
	    hPi0CrossSection->GetXaxis()->SetTitleOffset(1.0);
	    hPi0CrossSection->GetYaxis()->SetTitleOffset(1.03);
	    setHistFontSize(hPi0CrossSection, 240.0, 8.0);
	    hPi0CrossSection->Draw();
	    TLegend *lPi0CrossSection = new TLegend(0.511494, 0.53621, 0.933908, 0.80);
	    lPi0CrossSection->SetLineColor(0);
	    lPi0CrossSection->SetFillColor(0);
	    SetShadowColor(lPi0CrossSection, 0);
	    setLegendFontSize(lPi0CrossSection, 240.0, 8.0);
	    lPi0CrossSection->Draw();
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pi0_pp_crossection_1);
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pp_crossSection_spinpaper);
	    data_points_list dataPoints_pi0_pp_crossection_noname = dataPoints_pi0_pp_crossection;
	    for (data_points_list::iterator iter = dataPoints_pi0_pp_crossection_noname.begin();iter != dataPoints_pi0_pp_crossection_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSection, lPi0CrossSection, dataPoints_pi0_pp_crossection_noname);
	    hPi0CrossSection->SetMinimum(2e-10);
	    hPi0CrossSection->SetMaximum(0.9e+01);
	    hPi0CrossSection->GetXaxis()->SetRangeUser(0, 17);
	    gPad->SetLogy();
	    TLatex *label = new TLatex(7.45, 0.426, "#font[12]{p}_{ }+_{ }#font[12]{p} _{_{ }}at #surd#font[12]{#shiftx[-0.2]{s}}#shiftx[-0.67]{#shifty[-1]{#font[122]{-}}}=_{ }200_{ }Ge#shiftx[-0.12]{V}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();

	    {
	    TString namePi0CrossSectionDiv = "pp_pi0_crossection_div";
	    TString titlePi0CrossSectionDiv = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p} / theory";
	    cPi0CrossSection->cd(2);
	    gPad->SetPad(0, 0, 1, 0.4);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.30);
	    TString hPi0CrossSectionDivName = namePi0CrossSectionDiv; hPi0CrossSectionDivName += "_hist";
	    TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#shiftx[-0.3]{#font[12]{T}}} #shifty[-0.07]{[}Ge#shiftx[-0.12]{V}#shiftx[-0.4]{#shifty[-0.06]{/}}#shiftx[-0.3]{#font[12]{c}}#shifty[-0.07]{]};Data_{ }/_{ }theory    ";
	    TH1F *hPi0CrossSectionDiv = new TH1F(hPi0CrossSectionDivName.Data(), hPi0CrossSectionDivTitle.Data(), 1000, 0, 20);
	    hPi0CrossSectionDiv->GetXaxis()->SetTitleOffset(1.16);
	    hPi0CrossSectionDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0CrossSectionDiv->GetYaxis()->SetTitleOffset(0.69);
	    hPi0CrossSectionDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0CrossSectionDiv, 240.0, 8.0);
	    hPi0CrossSectionDiv->Draw();
	    TLegend *lPi0CrossSectionDiv = new TLegend(0.228448, 0.855967, 0.79023, 0.95679);
	    lPi0CrossSectionDiv->SetLineColor(0);
	    lPi0CrossSectionDiv->SetFillColor(0);
	    SetShadowColor(lPi0CrossSectionDiv, 0);
	    setLegendFontSize(lPi0CrossSectionDiv, 240.0, 8.0);
	    lPi0CrossSectionDiv->Draw();
	    data_points_list dataPointsDiv = dataPoints_pi0_pp_crossection_div;
	    if (fDIVPP) dataPointsDiv /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    data_points_list dataPointsOthersDiv = dataPoints_pp_crossSection_div;
	    if (fDIVPP) dataPointsOthersDiv /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv);
	    data_points_list dataPointsDiv_2 = dataPoints_pi0_pp_crossection_div_2;
	    if (fDIVPP) dataPointsDiv_2 /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    data_points_list dataPointsDiv_noname = dataPointsDiv_2;
	    for (data_points_list::iterator iter = dataPointsDiv_noname.begin();iter != dataPointsDiv_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_noname);
	    //TLatex *label = new TLatex(5.0, 2.5, "#font[12]{p}_{ }+_{ }#font[12]{p}");
	    //setLatexFontSize(label, 240.0, 8.0);
	    //label->Draw();
	    hPi0CrossSectionDiv->SetMinimum(0.001);
	    hPi0CrossSectionDiv->SetMaximum(3.5-0.001);
	    hPi0CrossSectionDiv->GetXaxis()->SetRangeUser(0, 17);
	    }
	}

	if (true) {
	    // p+p and d+Au cross section plots, comparison to other results
	    
	    TString namePi0CrossSection = "pp_dAu_pi0_crossection_others";
	    TString titlePi0CrossSection = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p} and #font[12]{d}+Au";
    	    TCanvas *cPi0CrossSection = new TCanvas(namePi0CrossSection.Data(), titlePi0CrossSection.Data(), 0, 0, 700, 700);
	    cPi0CrossSection->Divide(1, 2, 0, 0);

	    {
	    TString namePi0CrossSectionDiv = "pp_pi0_crossection_others_div";
	    TString titlePi0CrossSectionDiv = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p} / theory";
	    cPi0CrossSection->cd(1);
	    gPad->SetPad(0, 5.0*(1.0/9.0), 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString hPi0CrossSectionDivName = namePi0CrossSectionDiv; hPi0CrossSectionDivName += "_hist";
	    TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];Data_{ }/_{ }theory      ";
	    TH1F *hPi0CrossSectionDiv = new TH1F(hPi0CrossSectionDivName.Data(), hPi0CrossSectionDivTitle.Data(), 1000, 0, 20);
	    hPi0CrossSectionDiv->GetXaxis()->SetTitleOffset(1.00);
	    hPi0CrossSectionDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0CrossSectionDiv->GetYaxis()->SetTitleOffset(0.76);
	    hPi0CrossSectionDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0CrossSectionDiv, 240.0, 8.0);
	    hPi0CrossSectionDiv->Draw();
	    TLegend *lPi0CrossSectionDiv = new TLegend(0.278736, 0.611111, 0.666667, 0.909722);
	    lPi0CrossSectionDiv->SetLineColor(0);
	    lPi0CrossSectionDiv->SetFillColor(0);
	    SetShadowColor(lPi0CrossSectionDiv, 0);
	    setLegendFontSize(lPi0CrossSectionDiv, 240.0, 8.0);
	    lPi0CrossSectionDiv->Draw();
	    data_points_list dataPointsDiv = dataPoints_pi0_pp_crossection_div1;
	    if (fDIVPP) dataPointsDiv /= *fDIVPP;
	    data_points_list dataPointsDiv_2 = dataPoints_pi0_pp_crossection_div1_2;
	    if (fDIVPP) dataPointsDiv_2 /= *fDIVPP;
	    data_points_list dataPointsOthersDiv = dataPoints_pp_crossSection_div1;
	    if (fDIVPP) dataPointsOthersDiv /= *fDIVPP;
	    data_points_list dataPointsOthersDiv_2 = dataPoints_pp_crossSection_div1_2;
	    if (fDIVPP) dataPointsOthersDiv_2 /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv_2);
	    //showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    data_points_list dataPointsDiv_noname = dataPointsDiv_2;
	    for (data_points_list::iterator iter = dataPointsDiv_noname.begin();iter != dataPointsDiv_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_noname);
	    TLatex *label = new TLatex(16.7, 2.64, "#font[12]{p}_{ }+_{ }#font[12]{p}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    hPi0CrossSectionDiv->SetMinimum(0.001);
	    hPi0CrossSectionDiv->SetMaximum(3-0.001);
	    hPi0CrossSectionDiv->GetXaxis()->SetRangeUser(0, 20);
	    TLine *l1 = new TLine(3.4, 2.7, 4.3, 2.7);
	    l1->SetLineWidth(15);
	    l1->SetLineColor(0);
	    l1->Draw();
	    TMarker *m_pi0_1 = new TMarker(4.2, 2.7, kFullSquare);
	    m_pi0_1->SetMarkerColor(kBlack);
	    m_pi0_1->SetMarkerSize(1.1);
	    m_pi0_1->Draw();
	    TMarker *m_pi0_2 = new TMarker(3.6, 2.7, kFullCircle);
	    m_pi0_2->SetMarkerColor(kBlack);
	    m_pi0_2->SetMarkerSize(1.2);
	    m_pi0_2->Draw();
	    TLine *l = new TLine(3.4, 2.4, 4.3, 2.4);
	    l->SetLineWidth(15);
	    l->SetLineColor(0);
	    l->Draw();
	    TMarker *m_star_1 = new TMarker(4.2, 2.4, 29);
	    m_star_1->SetMarkerColor(kBlue);
	    m_star_1->SetMarkerSize(1.5);
	    m_star_1->Draw();
	    TMarker *m_star_2 = new TMarker(3.6, 2.4, 30);
	    m_star_2->SetMarkerColor(kBlue);
	    m_star_2->SetMarkerSize(1.5);
	    m_star_2->Draw();
	    TLine *l2 = new TLine(3.4, 2.1, 4.3, 2.1);
	    l2->SetLineWidth(15);
	    l2->SetLineColor(0);
	    l2->Draw();
	    TMarker *m_phenix_1 = new TMarker(4.2, 2.1, kOpenSquare);
	    m_phenix_1->SetMarkerColor(kBlack);
	    m_phenix_1->SetMarkerSize(1.1);
	    m_phenix_1->Draw();
	    TMarker *m_phenix_2 = new TMarker(3.6, 2.1, kOpenCircle);
	    m_phenix_2->SetMarkerColor(kBlack);
	    m_phenix_2->SetMarkerSize(1.1);
	    m_phenix_2->Draw();
	    }

	    {
	    TString namePi0CrossSectionDiv = "dAu_pi0_crossection_div";
	    TString titlePi0CrossSectionDiv = "#pi^{0} invariant cross section in #font[12]{d}+Au / theory";
	    cPi0CrossSection->cd(2);
	    gPad->SetPad(0, 0, 1, 5.0*(1.0/9.0));
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.20);
	    TString hPi0CrossSectionDivName = namePi0CrossSectionDiv; hPi0CrossSectionDivName += "_hist";
	    TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];Data_{ }/_{ }theory      ";
	    TH1F *hPi0CrossSectionDiv = new TH1F(hPi0CrossSectionDivName.Data(), hPi0CrossSectionDivTitle.Data(), 1000, 0, 20);
	    hPi0CrossSectionDiv->GetXaxis()->SetTitleOffset(1.05);
	    hPi0CrossSectionDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0CrossSectionDiv->GetYaxis()->SetTitleOffset(0.97);
	    hPi0CrossSectionDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0CrossSectionDiv, 240.0, 8.0);
	    hPi0CrossSectionDiv->Draw();
	    TLegend *lPi0CrossSectionDiv = new TLegend(0.428161, 0.634052, 0.656609, 0.934515);
	    lPi0CrossSectionDiv->SetLineColor(0);
	    lPi0CrossSectionDiv->SetFillColor(0);
	    SetShadowColor(lPi0CrossSectionDiv, 0);
	    setLegendFontSize(lPi0CrossSectionDiv, 240.0, 8.0);
	    lPi0CrossSectionDiv->Draw();
	    data_points_list dataPointsDiv = dataPoints_pi0_dAu_crossection_div1;
	    if (fDIV) dataPointsDiv /= *fDIV;
	    data_points_list dataPointsDiv_2 = dataPoints_pi0_dAu_crossection_div1_2;
	    if (fDIV) dataPointsDiv_2 /= *fDIV;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    data_points_list dataPointsOthersDiv = dataPoints_dAu_crossSection_div1;
	    if (fDIV) dataPointsOthersDiv /= *fDIV;
	    data_points_list dataPointsOthersDiv_2 = dataPoints_dAu_crossSection_div1_2;
	    if (fDIV) dataPointsOthersDiv_2 /= *fDIV;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv_2);
	    //showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    data_points_list dataPointsDiv_noname = dataPointsDiv_2;
	    for (data_points_list::iterator iter = dataPointsDiv_noname.begin();iter != dataPointsDiv_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_noname);
	    TLatex *label = new TLatex(16.7, 2.64, "#font[12]{d}_{ }+_{ }Au");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    hPi0CrossSectionDiv->SetMinimum(0.001);
	    hPi0CrossSectionDiv->SetMaximum(3-0.001);
	    hPi0CrossSectionDiv->GetXaxis()->SetRangeUser(0, 20);
	    delete lPi0CrossSectionDiv; lPi0CrossSectionDiv = 0;
	    }
	}

	if (true) {
	    // p+p cross section plot for the spin paper, comparison to other results

	    TString namePi0CrossSection = "pp_pi0_crossection_others";
	    TString titlePi0CrossSection = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p}";
    	    TCanvas *cPi0CrossSection = new TCanvas(namePi0CrossSection.Data(), titlePi0CrossSection.Data(), 0, 0, 700, 400);

	    {
	    TString namePi0CrossSectionDiv = "pp_pi0_crossection_others_div";
	    TString titlePi0CrossSectionDiv = "#pi^{0} invariant cross section in #font[12]{p}+#font[12]{p} / theory";
	    cPi0CrossSection->cd(1);
	    //gPad->SetPad(0, 0.5-0.06, 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.01);
	    gPad->SetBottomMargin(0.22);
	    TString hPi0CrossSectionDivName = namePi0CrossSectionDiv; hPi0CrossSectionDivName += "_hist";
	    TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#shiftx[-0.3]{#font[12]{T}}} #shifty[-0.07]{[}Ge#shiftx[-0.12]{V}#shiftx[-0.4]{#shifty[-0.06]{/}}#shiftx[-0.3]{#font[12]{c}}#shifty[-0.07]{]};Data_{ }/_{ }theory";
	    //TString hPi0CrossSectionDivTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [GeV/#font[12]{c}];Data_{ }/_{ }theory";
	    TH1F *hPi0CrossSectionDiv = new TH1F(hPi0CrossSectionDivName.Data(), hPi0CrossSectionDivTitle.Data(), 1000, 0, 20);
	    hPi0CrossSectionDiv->GetXaxis()->SetTitleOffset(1.16);
	    hPi0CrossSectionDiv->GetXaxis()->SetTickLength(0.06);
	    hPi0CrossSectionDiv->GetYaxis()->SetTitleOffset(0.96);
	    hPi0CrossSectionDiv->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0CrossSectionDiv, 240.0, 8.0);
	    hPi0CrossSectionDiv->Draw();
	    TLegend *lPi0CrossSectionDiv = new TLegend(0.280172, 0.702796, 0.668103, 0.950108);
	    lPi0CrossSectionDiv->SetLineColor(0);
	    lPi0CrossSectionDiv->SetFillColor(0);
	    SetShadowColor(lPi0CrossSectionDiv, 0);
	    setLegendFontSize(lPi0CrossSectionDiv, 240.0, 8.0);
	    lPi0CrossSectionDiv->Draw();
	    data_points_list dataPointsDiv = dataPoints_pi0_pp_crossection_div1;
	    if (fDIVPP) dataPointsDiv /= *fDIVPP;
	    data_points_list dataPointsDiv_2 = dataPoints_pi0_pp_crossection_div1_2;
	    if (fDIVPP) dataPointsDiv_2 /= *fDIVPP;
	    data_points_list dataPointsOthersDiv = dataPoints_pp_crossSection_div1;
	    if (fDIVPP) dataPointsOthersDiv /= *fDIVPP;
	    data_points_list dataPointsOthersDiv_2 = dataPoints_pp_crossSection_div1_2;
	    if (fDIVPP) dataPointsOthersDiv_2 /= *fDIVPP;
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv);
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsOthersDiv_2);
	    //showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_2);
	    data_points_list dataPointsDiv_noname = dataPointsDiv_2;
	    for (data_points_list::iterator iter = dataPointsDiv_noname.begin();iter != dataPointsDiv_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0CrossSectionDiv, lPi0CrossSectionDiv, dataPointsDiv_noname);
	    TLatex *label = new TLatex(17.0*17.0/20.0/*17.0*/, 2.777, "#font[12]{p}_{ }+_{ }#font[12]{p}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    hPi0CrossSectionDiv->SetMinimum(0.001);
	    hPi0CrossSectionDiv->SetMaximum(3.2 - 0.001);
	    hPi0CrossSectionDiv->GetXaxis()->SetRangeUser(0, 19.99);
	    /*TLine *l1 = new TLine(5.0, 2.55, 5.9, 2.55);
	    l1->SetLineWidth(15);
	    l1->SetLineColor(0);
	    l1->Draw();
	    TMarker *m_pi0_1 = new TMarker(5.8, 2.55, kFullSquare);
	    m_pi0_1->SetMarkerColor(kBlack);
	    m_pi0_1->SetMarkerSize(1.1);
	    m_pi0_1->Draw();
	    TMarker *m_pi0_2 = new TMarker(5.2, 2.55, kFullCircle);
	    m_pi0_2->SetMarkerColor(kBlack);
	    m_pi0_2->SetMarkerSize(1.2);
	    m_pi0_2->Draw();*/
	    TLine *l = new TLine(4.1*17.0/20.0, 2.53, 5.2*17.0/20.0, 2.53);
	    l->SetLineWidth(15);
	    l->SetLineColor(0);
	    l->Draw();
	    TMarker *m_star_1 = new TMarker(4.4*17.0/20.0, 2.53, 29);
	    m_star_1->SetMarkerColor(kBlue);
	    m_star_1->SetMarkerSize(1.5);
	    m_star_1->Draw();
	    TMarker *m_star_2 = new TMarker(4.9*17.0/20.0, 2.53, 30);
	    m_star_2->SetMarkerColor(kBlue);
	    m_star_2->SetMarkerSize(1.5);
	    m_star_2->Draw();
	    /*TLine *l2 = new TLine(5.0, 1.79, 5.9, 1.79);
	    l2->SetLineWidth(15);
	    l2->SetLineColor(0);
	    l2->Draw();
	    TMarker *m_phenix_1 = new TMarker(5.8, 1.82, kOpenSquare);
	    m_phenix_1->SetMarkerColor(kBlack);
	    m_phenix_1->SetMarkerSize(1.1);
	    m_phenix_1->Draw();
	    TMarker *m_phenix_2 = new TMarker(5.2, 1.82, kOpenCircle);
	    m_phenix_2->SetMarkerColor(kBlack);
	    m_phenix_2->SetMarkerSize(1.1);
	    m_phenix_2->Draw();*/
	    }
	}

	if (true) {
	    // d+Au pi0 and eta R_CP, comparison to other results, theory etc.

	    TString namePi0Rcp = "dAu_pi0_Rcp";
	    TString titlePi0Rcp = "#pi^{0}^{} #font[12]{R}_{#font[12]{CP}} in #font[12]{d}+Au";
    	    TCanvas *cPi0Rcp = new TCanvas(namePi0Rcp.Data(), titlePi0Rcp.Data(), 0, 0, 700, 400);

	    {
	    cPi0Rcp->cd();
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.03);
	    gPad->SetBottomMargin(0.20);
	    TString hPi0RcpName = namePi0Rcp; hPi0RcpName += "_hist";
	    TString hPi0RcpTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];#font[12]{R}_{#font[12]{CP}}";
	    TH1F *hPi0Rcp = new TH1F(hPi0RcpName.Data(), hPi0RcpTitle.Data(), 1000, 0, 20);
	    hPi0Rcp->GetXaxis()->SetTitleOffset(1.15);
	    hPi0Rcp->GetXaxis()->SetTickLength(0.05);
	    hPi0Rcp->GetYaxis()->SetTitleOffset(0.94);
	    hPi0Rcp->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hPi0Rcp, 240.0, 8.0);
	    hPi0Rcp->Draw();
	    TLegend *lPi0Rcp = new TLegend(0.557471, 0.78, 0.952586, 0.95);
	    lPi0Rcp->SetLineColor(0);
	    lPi0Rcp->SetFillColor(0);
	    SetShadowColor(lPi0Rcp, 0);
	    setLegendFontSize(lPi0Rcp, 240.0, 8.0);
	    lPi0Rcp->Draw();
	    showResultsOthers(hPi0Rcp, lPi0Rcp, dataPoints_pi0_dAu_Rcp_2);
	    showResultsOthers(hPi0Rcp, lPi0Rcp, dataPoints_pi0_dAu_Rcp);
	    showResultsOthers(hPi0Rcp, lPi0Rcp, dataPoints_Rcp_pi0);
	    data_points_list dataPoints_pi0_dAu_Rcp_noname = dataPoints_pi0_dAu_Rcp;
	    for (data_points_list::iterator iter = dataPoints_pi0_dAu_Rcp_noname.begin();iter != dataPoints_pi0_dAu_Rcp_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0Rcp, lPi0Rcp, dataPoints_pi0_dAu_Rcp_noname);
	    TLatex *label = new TLatex(2.4, 2.214, "#font[12]{R}_{#font[12]{CP}}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    TLine *l = new TLine(9.0, 2.02, 9.65, 2.02);
	    l->SetLineWidth(15);
	    l->SetLineColor(0);
	    l->Draw();
	    TMarker *m_star_1 = new TMarker(9.55, 2.02, 29);
	    m_star_1->SetMarkerColor(kBlue);
	    m_star_1->SetMarkerSize(1.5);
	    m_star_1->Draw();
	    TMarker *m_star_2 = new TMarker(9.05, 2.02, 30);
	    m_star_2->SetMarkerColor(kBlue);
	    m_star_2->SetMarkerSize(1.5);
	    m_star_2->Draw();
	    hPi0Rcp->SetMinimum(0.0);
	    hPi0Rcp->SetMaximum(2.499);
	    hPi0Rcp->GetXaxis()->SetRangeUser(0, 17);
	    }
	    
	    if (false) {
	    TString nameEtaRcp = "dAu_eta_Rcp";
	    TString titleEtaRcp = "#eta #font[12]{R}_{#font[12]{CP}} in #font[12]{d}+Au";
    	    TCanvas *cEtaRcp = new TCanvas(nameEtaRcp.Data(), titleEtaRcp.Data(), 0, 0, 700, 375);
	    cEtaRcp->cd();
	    gPad->SetLeftMargin(0.1);
	    gPad->SetRightMargin(0.05);
	    gPad->SetTopMargin(0.02);
	    gPad->SetBottomMargin(0.17);
	    TString hEtaRcpName = nameEtaRcp; hEtaRcpName += "_hist";
	    TString hEtaRcpTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];#font[12]{R}_{#font[12]{CP}}";
	    TH1F *hEtaRcp = new TH1F(hEtaRcpName.Data(), hEtaRcpTitle.Data(), 1000, 0, 20);
	    hEtaRcp->GetXaxis()->SetTitleOffset(1.0);
	    hEtaRcp->GetXaxis()->SetTickLength(0.05);
	    hEtaRcp->GetYaxis()->SetTitleOffset(0.70);
	    hEtaRcp->GetYaxis()->SetNdivisions(510);
	    setHistFontSize(hEtaRcp, 240.0, 8.0);
	    hEtaRcp->Draw();
	    TLegend *lEtaRcp = new TLegend(0.52, 0.78, 0.87, 0.95);
	    lEtaRcp->SetLineColor(0);
	    lEtaRcp->SetFillColor(0);
	    SetShadowColor(lEtaRcp, 0);
	    setLegendFontSize(lEtaRcp, 240.0, 8.0);
	    lEtaRcp->Draw();
	    showResultsOthers(hEtaRcp, lEtaRcp, dataPoints_eta_dAu_Rcp_2);
	    showResultsOthers(hEtaRcp, lEtaRcp, dataPoints_eta_dAu_Rcp);
	    showResultsOthers(hEtaRcp, lEtaRcp, dataPoints_Rcp_eta);
	    data_points_list dataPoints_eta_dAu_Rcp_noname = dataPoints_eta_dAu_Rcp;
	    for (data_points_list::iterator iter = dataPoints_eta_dAu_Rcp_noname.begin();iter != dataPoints_eta_dAu_Rcp_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hEtaRcp, lEtaRcp, dataPoints_eta_dAu_Rcp_noname);
	    hEtaRcp->SetMinimum(0.0);
	    hEtaRcp->SetMaximum(2.999);
	    hEtaRcp->GetXaxis()->SetRangeUser(0, 17);
	    }
	}

	if (true) {
	    // d+Au pi0 and eta R_dA, comparison to other results, theory etc.


	    TString nameRdA = "RdA_pi0_eta";
	    TString titleRdA = "#pi^{0} and #eta #font[12]{R}_{#font[12]{dA}}";
    	    TCanvas *cRdA = new TCanvas(nameRdA.Data(), titleRdA.Data(), 0, 0, 700, 700);
	    cRdA->Divide(1, 2, 0, 0);

	    {
	    cRdA->cd(1);
	    gPad->SetPad(0, 5.0*(1.0/9.0), 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString namePi0RdA = "pi0_RdA";
	    TString titlePi0RdA = "#pi^{0}^{} #font[12]{R}_{#font[12]{dA}}";
	    TString hPi0RdAName = namePi0RdA; hPi0RdAName += "_hist";
	    TString hPi0RdATitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];#font[12]{R}_{#font[12]{dA}}";
	    TH1F *hPi0RdA = new TH1F(hPi0RdAName.Data(), hPi0RdATitle.Data(), 1000, 0, 20);
	    hPi0RdA->GetXaxis()->SetTitleOffset(1.0);
	    hPi0RdA->GetXaxis()->SetTickLength(0.05);
	    hPi0RdA->GetYaxis()->SetTitleOffset(0.76);
	    hPi0RdA->GetYaxis()->SetNdivisions(510);
	    setHistFontSize(hPi0RdA, 240.0, 8.0);
	    hPi0RdA->Draw();
	    TLegend *lPi0RdA = new TLegend(0.54454, 0.628348, 0.95977, 0.926339);
	    lPi0RdA->SetLineColor(0);
	    lPi0RdA->SetFillColor(0);
	    SetShadowColor(lPi0RdA, 0);
	    setLegendFontSize(lPi0RdA, 240.0, 8.0);
	    lPi0RdA->Draw();
	    showResultsOthers(hPi0RdA, lPi0RdA, dataPoints_pi0_RdA_2);
	    showResultsOthers(hPi0RdA, lPi0RdA, dataPoints_pi0_RdA);
	    showResultsOthers(hPi0RdA, lPi0RdA, dataPoints_RdA_pi0);
	    data_points_list dataPoints_pi0_RdA_noname = dataPoints_pi0_RdA;
	    for (data_points_list::iterator iter = dataPoints_pi0_RdA_noname.begin();iter != dataPoints_pi0_RdA_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hPi0RdA, lPi0RdA, dataPoints_pi0_RdA_noname);
	    TLatex *label = new TLatex(2.4, 2.724, "#font[12]{R}_{#font[12]{dA}}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    hPi0RdA->SetMinimum(0.001);
	    hPi0RdA->SetMaximum(2.999);
	    hPi0RdA->GetXaxis()->SetRangeUser(0, 17);
	    }
	    
	    {
	    cRdA->cd(2);
	    gPad->SetPad(0, 0, 1, 5.0*(1.0/9.0));
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.20);
	    TString nameEtaRdA = "eta_RdA";
	    TString titleEtaRdA = "#eta #font[12]{R}_{#font[12]{dA}}";
	    TString hEtaRdAName = nameEtaRdA; hEtaRdAName += "_hist";
	    TString hEtaRdATitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];#font[12]{R}_{#font[12]{dA}}";
	    TH1F *hEtaRdA = new TH1F(hEtaRdAName.Data(), hEtaRdATitle.Data(), 1000, 0, 20);
	    hEtaRdA->GetXaxis()->SetTitleOffset(1.15);
	    hEtaRdA->GetXaxis()->SetTickLength(0.05);
	    hEtaRdA->GetYaxis()->SetTitleOffset(0.97);
	    hEtaRdA->GetYaxis()->SetNdivisions(505);
	    setHistFontSize(hEtaRdA, 240.0, 8.0);
	    hEtaRdA->Draw();
	    TLegend *lEtaRdA = new TLegend(0.54454, 0.80625, 0.95977, 0.976339);
	    lEtaRdA->SetLineColor(0);
	    lEtaRdA->SetFillColor(0);
	    SetShadowColor(lEtaRdA, 0);
	    setLegendFontSize(lEtaRdA, 240.0, 8.0);
	    lEtaRdA->Draw();
	    showResultsOthers(hEtaRdA, lEtaRdA, dataPoints_eta_RdA_2);
	    showResultsOthers(hEtaRdA, lEtaRdA, dataPoints_eta_RdA);
	    showResultsOthers(hEtaRdA, lEtaRdA, dataPoints_RdA_eta);
	    data_points_list dataPoints_eta_RdA_noname = dataPoints_eta_RdA;
	    for (data_points_list::iterator iter = dataPoints_eta_RdA_noname.begin();iter != dataPoints_eta_RdA_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hEtaRdA, lEtaRdA, dataPoints_eta_RdA_noname);
	    hEtaRdA->SetMinimum(0.0);
	    hEtaRdA->SetMaximum(4.499);
	    hEtaRdA->GetXaxis()->SetRangeUser(0, 17);
	    }
	}

	if (true) {
	    // d+Au and p+p eta/pi0
	    
	    TString nameEtaToPi0 = "etaToPi0_pp_dAu";
	    TString titleEtaToPi0 = "#eta/#pi^{0} in #font[12]{p}+#font[12]{p} and #font[12]{d}+Au";
    	    TCanvas *cEtaToPi0 = new TCanvas(nameEtaToPi0.Data(), titleEtaToPi0.Data(), 0, 0, 700, 700);
	    cEtaToPi0->Divide(1, 2, 0, 0);

	    {
	    cEtaToPi0->cd(1);
	    gPad->SetPad(0, 5.0*(1.0/9.0), 1, 1);
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0.05);
	    gPad->SetBottomMargin(0);
	    TString hEtaToPi0PPName = nameEtaToPi0; hEtaToPi0PPName += "_pp_hist";
	    TString hEtaToPi0PPTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];#font[12]{R}_{#eta/#pi}";
	    TH1F *hEtaToPi0PP = new TH1F(hEtaToPi0PPName.Data(), hEtaToPi0PPTitle.Data(), 1000, 0, 20);
	    hEtaToPi0PP->GetXaxis()->SetTitleOffset(1.0);
	    hEtaToPi0PP->GetXaxis()->SetTickLength(0.05);
	    hEtaToPi0PP->GetYaxis()->SetTitleOffset(0.74);
	    hEtaToPi0PP->GetYaxis()->SetNdivisions(510);
	    setHistFontSize(hEtaToPi0PP, 240.0, 8.0);
	    hEtaToPi0PP->Draw();
	    TLegend *lEtaToPi0PP = new TLegend(0.452586, 0.53125, 0.95977, 0.933036);
	    lEtaToPi0PP->SetLineColor(0);
	    lEtaToPi0PP->SetFillColor(0);
	    SetShadowColor(lEtaToPi0PP, 0);
	    setLegendFontSize(lEtaToPi0PP, 240.0, 8.0);
	    lEtaToPi0PP->Draw();
	    showResultsOthers(hEtaToPi0PP, lEtaToPi0PP, dataPoints_etaToPi0_pp1);
	    showResultsOthers(hEtaToPi0PP, lEtaToPi0PP, dataPoints_etaToPi0_pp);
	    showResultsOthers(hEtaToPi0PP, lEtaToPi0PP, dataPoints_pp_etatopi);
	    data_points_list dataPoints_etaToPi0_pp_noname = dataPoints_etaToPi0_pp;
	    for (data_points_list::iterator iter = dataPoints_etaToPi0_pp_noname.begin();iter != dataPoints_etaToPi0_pp_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hEtaToPi0PP, lEtaToPi0PP, dataPoints_etaToPi0_pp_noname);
	    hEtaToPi0PP->GetXaxis()->SetRangeUser(0, 17);
	    hEtaToPi0PP->SetMinimum(0.001);
	    hEtaToPi0PP->SetMaximum(1.499);
	    TLatex *label = new TLatex(2.4, 1.3, "#eta/#pi^{0}^{}");
	    setLatexFontSize(label, 240.0, 8.0);
	    label->Draw();
	    }
	    
	    {
	    TString nameEtaToPi0DAu = "etaToPi0_dAu";
	    TString titleEtaToPi0DAu = "#eta/#pi^{0} in #font[12]{d}+Au";
	    cEtaToPi0->cd(2);
	    gPad->SetPad(0, 0, 1, 5.0*(1.0/9.0));
	    gPad->SetLeftMargin(0.17);
	    gPad->SetRightMargin(0.03);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0.20);
	    TString hEtaToPi0DAuName = nameEtaToPi0DAu; hEtaToPi0DAuName += "_dAu_hist";
	    TString hEtaToPi0DAuTitle = ";#font[12]{p}_{#font[12]{T}}_{ } [_{_{_{ }}}GeV/#font[12]{c}_{_{ }}];#font[12]{R}_{#eta/#pi}";
	    TH1F *hEtaToPi0DAu = new TH1F(hEtaToPi0DAuName.Data(), hEtaToPi0DAuTitle.Data(), 1000, 0, 20);
	    hEtaToPi0DAu->GetXaxis()->SetTitleOffset(1.15);
	    hEtaToPi0DAu->GetXaxis()->SetTickLength(0.05);
	    hEtaToPi0DAu->GetYaxis()->SetTitleOffset(0.92);
	    hEtaToPi0DAu->GetYaxis()->SetNdivisions(510);
	    setHistFontSize(hEtaToPi0DAu, 240.0, 8.0);
	    hEtaToPi0DAu->Draw();
	    TLegend *lEtaToPi0DAu = new TLegend(0.452586, 0.833036, 0.95977, 0.983036);
	    lEtaToPi0DAu->SetLineColor(0);
	    lEtaToPi0DAu->SetFillColor(0);
	    SetShadowColor(lEtaToPi0DAu, 0);
	    setLegendFontSize(lEtaToPi0DAu, 240.0, 8.0);
	    lEtaToPi0DAu->Draw();
	    showResultsOthers(hEtaToPi0DAu, lEtaToPi0DAu, dataPoints_etaToPi0_dAu1);
	    showResultsOthers(hEtaToPi0DAu, lEtaToPi0DAu, dataPoints_etaToPi0_dAu);
	    showResultsOthers(hEtaToPi0DAu, lEtaToPi0DAu, dataPoints_dAu_etatopi);
	    data_points_list dataPoints_etaToPi0_dAu_noname = dataPoints_etaToPi0_dAu;
	    for (data_points_list::iterator iter = dataPoints_etaToPi0_dAu_noname.begin();iter != dataPoints_etaToPi0_dAu_noname.end();++iter) (*iter).legendTitle = "";
	    showResultsOthers(hEtaToPi0DAu, lEtaToPi0DAu, dataPoints_etaToPi0_dAu_noname);
	    hEtaToPi0DAu->GetXaxis()->SetRangeUser(0, 17);
	    hEtaToPi0DAu->SetMinimum(0.0);
	    hEtaToPi0DAu->SetMaximum(1.499);
	    }
	}

	if (true) {
	    // p+p HighTower-1 (4 < pT < 5) bin plots: event mixing, jet-aligned mixing, low mass background
	    
	    const TH1F *sigBg = allResults.resultsPP.resultsHT1.invlist.size() ? (*find(allResults.resultsPP.resultsHT1.invlist.begin(), allResults.resultsPP.resultsHT1.invlist.end(), TInvariantMassDistribution(TBinParameters(pT, 4, 5)))).getDistribution() : 0;
	    const TH1F *mixBg = allResults.resultsPP.resultsHT1.invBgRandomlist.size() ? (*find(allResults.resultsPP.resultsHT1.invBgRandomlist.begin(), allResults.resultsPP.resultsHT1.invBgRandomlist.end(), TInvariantMassDistribution(TBinParameters(pT, 4, 5)))).getDistribution() : 0;
	    const TH1F *jetmixBg = allResults.resultsPP.resultsHT1.invBglist.size() ? (*find(allResults.resultsPP.resultsHT1.invBglist.begin(), allResults.resultsPP.resultsHT1.invBglist.end(), TInvariantMassDistribution(TBinParameters(pT, 4, 5)))).getDistribution() : 0;
	    const TH1F *gammaBg = allResults.resultsPP.resultsHT1.inv1gammalist.size() ? (*find(allResults.resultsPP.resultsHT1.inv1gammalist.begin(), allResults.resultsPP.resultsHT1.inv1gammalist.end(), TInvariantMassDistribution(TBinParameters(pT, 4, 5)))).getDistribution() : 0;
	    const TH1F *nbarBg = allResults.resultsPP.resultsHT1.invnbarlist.size() ? (*find(allResults.resultsPP.resultsHT1.invnbarlist.begin(), allResults.resultsPP.resultsHT1.invnbarlist.end(), TInvariantMassDistribution(TBinParameters(pT, 4, 5)))).getDistribution() : 0;
	    Int_t nrebin = allResults.resultsPP.resultsHT1.nrebin.size() ? (Int_t)(*find(allResults.resultsPP.resultsHT1.nrebin.begin(), allResults.resultsPP.resultsHT1.nrebin.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    Float_t scaleBg = allResults.resultsPP.resultsHT1.scaleBg.size() ? (*find(allResults.resultsPP.resultsHT1.scaleBg.begin(), allResults.resultsPP.resultsHT1.scaleBg.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    Float_t scaleBgRandom = allResults.resultsPP.resultsHT1.scaleBgRandom.size() ? (*find(allResults.resultsPP.resultsHT1.scaleBgRandom.begin(), allResults.resultsPP.resultsHT1.scaleBgRandom.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    Float_t scaleBgLowMass = allResults.resultsPP.resultsHT1.scaleBgLowMass.size() ? (*find(allResults.resultsPP.resultsHT1.scaleBgLowMass.begin(), allResults.resultsPP.resultsHT1.scaleBgLowMass.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    Float_t scaleBgLowMass2 = allResults.resultsPP.resultsHT1.scaleBgLowMass2.size() ? (*find(allResults.resultsPP.resultsHT1.scaleBgLowMass2.begin(), allResults.resultsPP.resultsHT1.scaleBgLowMass2.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    Float_t bgJetToTotal = allResults.resultsPP.resultsHT1.bgJetToTotal.size() ? (*find(allResults.resultsPP.resultsHT1.bgJetToTotal.begin(), allResults.resultsPP.resultsHT1.bgJetToTotal.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    Float_t bgToSigbg = allResults.resultsPP.resultsHT1.bgToSigbg.size() ? (*find(allResults.resultsPP.resultsHT1.bgToSigbg.begin(), allResults.resultsPP.resultsHT1.bgToSigbg.end(), TBinParameters(pT, 4, 5))).getValue() : 0;
	    if (sigBg && mixBg && jetmixBg && gammaBg && nbarBg) {
	    
	    TString nameMix = "pp_HT1_45_rndmix";
	    TString titleMix = "p+p HighTower-1 4 < pT < 5 GeV/c, event mixing";
	    TString labelText = "#splitline{#font[12]{p}_{ }+_{ }#font[12]{p} HighTower-1}{4 < #font[12]{p}_{#font[12]{T}}_{ }< 5 GeV/#font[12]{c}}";
	    TCanvas *cMix = new TCanvas(nameMix.Data(), titleMix.Data(), 0, 0, 600, 3*300);
	    //cMix->SetCanvasSize(600, 3*300);
	    gPad->SetLeftMargin(0);
	    gPad->SetRightMargin(0);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0);
	    cMix->Divide(1, 3, 0, 0);
	    cMix->cd(1);
	    gPad->SetPad(0, 2.25*(1.0/3.25), 1, 1);
	    gPad->SetLeftMargin(0.13);
	    gPad->SetRightMargin(0.01);
	    TH1F *hSigBg = new TH1F(*sigBg);
	    hSigBg->SetNameTitle("same", ";#font[12]{M}_{#gamma#gamma} [GeV/#font[12]{c}^{2}];Counts  ");
	    if (nrebin > 1) hSigBg->Rebin(nrebin);
	    hSigBg->GetXaxis()->SetLabelFont(132);
	    hSigBg->GetXaxis()->SetTitleFont(132);
	    hSigBg->GetYaxis()->SetLabelFont(132);
	    hSigBg->GetYaxis()->SetTitleFont(132);
	    hSigBg->GetXaxis()->SetRangeUser(0.0 + hSigBg->GetXaxis()->GetBinWidth(1), 1.2 - 3*hSigBg->GetXaxis()->GetBinWidth(1));
	    hSigBg->SetMinimum(-149);
	    hSigBg->SetMaximum(1049);
	    hSigBg->GetXaxis()->SetTitleOffset(1.1);
	    hSigBg->GetYaxis()->SetTitleOffset(0.68);
	    setHistFontSize(hSigBg, 240.0, 8.0);
	    hSigBg->DrawCopy();
	    TLine *zeroLineSigbgMix = new TLine(0.02, 0, 1.18, 0);
	    zeroLineSigbgMix->Draw();
	    TLatex *labelMix = new TLatex(0.4, 750, labelText.Data());
	    setLatexFontSize(labelMix, 240.0, 8.0);
	    labelMix->Draw();
	    TLatex *labelMix1 = new TLatex(0.4, 400, "(a) Event mixing");
	    setLatexFontSize(labelMix1, 240.0, 8.0);
	    labelMix1->Draw();
	    TH1F *hMixBg = new TH1F(*mixBg);
	    hMixBg->SetName("mix");
	    if (nrebin > 1) hMixBg->Rebin(nrebin);
	    hMixBg->Scale(hSigBg->Integral(hSigBg->GetXaxis()->FindBin(0.9), hSigBg->GetXaxis()->FindBin(1.2)) / hMixBg->Integral(hMixBg->GetXaxis()->FindBin(0.9), hMixBg->GetXaxis()->FindBin(1.2)));
	    hMixBg->Draw("HIST SAME");
	    TH1F *hMixBgNorm = new TH1F(*hMixBg);
	    hMixBgNorm->SetName("mix_norm");
	    hMixBgNorm->GetXaxis()->SetRangeUser(0.9, 1.2);
	    hMixBgNorm->SetFillColor(17);
	    hMixBgNorm->Draw("HIST SAME ][");
	    cMix->cd(2);
	    gPad->SetPad(0, 1.25*(1.0/3.25), 1, 2.25*(1.0/3.25));
	    gPad->SetLeftMargin(0.13);
	    gPad->SetRightMargin(0.01);
/*	    TH1F *hMixSub = new TH1F(*hSigBg);
	    hMixSub->SetName("mixsub");
	    hMixSub->Add(hMixBg, -1);
	    hMixSub->SetMinimum(-199);
	    hMixSub->SetMaximum(1249);
	    //hMixSub->GetYaxis()->LabelsOption("");
	    setHistFontSize(hMixSub, 240.0, 8.0);
	    hMixSub->Draw();
	    TLine *zeroLineMixSub = new TLine(0, 0, 1.2, 0);
	    zeroLineMixSub->Draw();
	    
	    TString nameJetmix = "pp_HT1_45_jetmix";
	    TString titleJetmix = "p+p HighTower-1 4 < pT < 5 GeV/c, jet-aligned event mixing";
	    TCanvas *cJetmix = new TCanvas(nameJetmix.Data(), titleJetmix.Data());
	    cJetmix->SetCanvasSize(600, 2*400);
	    cJetmix->Divide(1, 2, 0, 0);
	    cJetmix->cd(1);
*/	    setHistFontSize(hSigBg, 240.0, 8.0);
	    hSigBg->DrawCopy();
	    TLine *zeroLineSigbgJetmix = new TLine(0.02, 0, 1.18, 0);
	    zeroLineSigbgJetmix->Draw();
	    //TLatex *labelJetmix = new TLatex(0.5, 600, labelText.Data());
	    //labelJetmix->SetTextSize(0.055);
	    //setLatexFontSize(labelJetmix, 240.0, 8.0);
	    //labelJetmix->Draw();
	    TLatex *labelMix2 = new TLatex(0.4, 400, "(b) Jet-aligned mixing");
	    setLatexFontSize(labelMix2, 240.0, 8.0);
	    labelMix2->Draw();
	    TH1F *hMixBgRnd = new TH1F(*mixBg);
	    hMixBgRnd->SetName("mix_rnd");
	    if (nrebin > 1) hMixBgRnd->Rebin(nrebin);
	    hMixBgRnd->Scale(scaleBgRandom * bgToSigbg * (1.0 - bgJetToTotal));
	    TH1F *hMixBgJet = new TH1F(*jetmixBg);
	    hMixBgJet->SetName("mix_jet");
	    if (nrebin > 1) hMixBgJet->Rebin(nrebin);
	    hMixBgJet->Scale(scaleBg * bgToSigbg * bgJetToTotal);
	    TH1F *hJetmixBg = new TH1F(*hMixBgRnd);
	    hJetmixBg->SetName("jetmix");
	    hJetmixBg->Add(hMixBgJet);
	    hJetmixBg->Draw("HIST SAME");
	    TH1F *hJetmixBgNorm1 = new TH1F(*hJetmixBg);
	    hJetmixBgNorm1->SetName("jetmix_norm1");
	    hJetmixBgNorm1->GetXaxis()->SetRangeUser(0.3, 0.4);
	    hJetmixBgNorm1->SetLineColor(kWhite);
	    hJetmixBgNorm1->SetFillColor(17);
	    hJetmixBgNorm1->Draw("HIST SAME ][");
	    TH1F *hJetmixBgNorm2 = new TH1F(*hJetmixBg);
	    hJetmixBgNorm2->SetName("jetmix_norm2");
	    hJetmixBgNorm2->GetXaxis()->SetRangeUser(0.9, 1.2);
	    hJetmixBgNorm2->SetFillColor(17);
	    hJetmixBgNorm2->Draw("HIST SAME ][");
	    hJetmixBg->Draw("HIST SAME");
//	    cJetmix->cd(2);
	    cMix->cd(3);
	    gPad->SetPad(0, 0, 1, 1.25*(1.0/3.25));
	    gPad->SetLeftMargin(0.13);
	    gPad->SetRightMargin(0.01);
	    gPad->SetBottomMargin(0.20);
	    TH1F *hJetmixSub = new TH1F(*hSigBg);
	    hJetmixSub->SetName("jetmixsub");
	    hJetmixSub->Add(hJetmixBg, -1);
	    hJetmixSub->SetMinimum(-149);
	    hJetmixSub->SetMaximum(1049);
	    //hJetmixSub->GetYaxis()->LabelsOption("");
	    setHistFontSize(hJetmixSub, 240.0, 8.0);
	    hJetmixSub->GetYaxis()->SetTitleOffset(0.85);
	    hJetmixSub->GetXaxis()->SetTitleOffset(1.1);
	    hJetmixSub->Draw();
	    TLine *zeroLineJetmixSub = new TLine(0.02, 0, 1.18, 0);
	    zeroLineJetmixSub->Draw();
	    TLatex *labelMix3 = new TLatex(0.4, 400, "#splitline{(c) Combinatorial background}{     _{ }_{ }subtracted}");
	    setLatexFontSize(labelMix3, 240.0, 8.0);
	    labelMix3->Draw();

	    TString nameLowmass = "pp_HT1_45_lowmass";
	    TString titleLowmass = "p+p HighTower-1 4 < pT < 5 GeV/c, low mass background";
	    TCanvas *cLowmass = new TCanvas(nameLowmass.Data(), titleLowmass.Data(), 0, 0, 600, 2*300);
	    //cLowmass->SetCanvasSize(600, 2*300);
	    gPad->SetLeftMargin(0);
	    gPad->SetRightMargin(0);
	    gPad->SetTopMargin(0);
	    gPad->SetBottomMargin(0);
	    cLowmass->Divide(1, 2, 0, 0);
	    cLowmass->cd(1);
	    gPad->SetPad(0, 1.25*(1.0/2.25), 1, 1);
	    gPad->SetLeftMargin(0.13);
	    gPad->SetRightMargin(0.01);
	    TH1F *hSig = new TH1F(*hJetmixSub);
	    hSig->SetNameTitle("same", ";#font[12]{M}_{#gamma#gamma} [GeV/#font[12]{c}^{2}];Counts  ");
	    hSig->GetYaxis()->LabelsOption("v");
	    hSig->GetYaxis()->SetTitleOffset(0.68);
	    setHistFontSize(hSig, 240.0, 8.0);
	    hSig->Draw();
	    TLine *zeroLineSigLowmass = new TLine(0.02, 0, 1.18, 0);
	    zeroLineSigLowmass->Draw();
	    TLatex *labelLowmass = new TLatex(0.5, 600, labelText.Data());
	    labelLowmass->SetTextSize(0.055);
	    setLatexFontSize(labelLowmass, 240.0, 8.0);
	    labelLowmass->Draw();
	    TH1F *hBgLowmass = new TH1F(*gammaBg);
	    hBgLowmass->SetName("lowmass");
	    if (nrebin > 1) hBgLowmass->Rebin(nrebin);
	    hBgLowmass->Scale(scaleBgLowMass);
	    TH1F *hBgLowmass2 = new TH1F(*nbarBg);
	    hBgLowmass2->SetName("lowmass2");
	    if (nrebin > 1) hBgLowmass2->Rebin(nrebin);
	    hBgLowmass2->Scale(scaleBgLowMass2);
	    hBgLowmass->Add(hBgLowmass2);
	    hBgLowmass->Draw("HIST SAME");
	    cLowmass->cd(2);
	    gPad->SetPad(0, 0, 1, 1.25*(1.0/2.25));
	    gPad->SetLeftMargin(0.13);
	    gPad->SetRightMargin(0.01);
	    gPad->SetBottomMargin(0.20);
	    TH1F *hLowmassSub = new TH1F(*hSig);
	    hLowmassSub->SetName("lowmasssub");
	    hLowmassSub->Add(hBgLowmass, -1);
	    hLowmassSub->SetMinimum(-149);
	    hLowmassSub->SetMaximum(1049);
	    //hLowmassSub->GetYaxis()->LabelsOption("");
	    hLowmassSub->GetXaxis()->SetTitleOffset(1.1);
	    hLowmassSub->GetYaxis()->SetTitleOffset(0.85);
	    setHistFontSize(hLowmassSub, 240.0, 8.0);
	    hLowmassSub->Draw();
	    TLine *zeroLineLowmassSub = new TLine(0.02, 0, 1.18, 0);
	    zeroLineLowmassSub->Draw();
	    TLatex *labelBgSubtracted = new TLatex(0.4, 400, "#splitline{Combinatorial and low mass}{background subtracted}");
	    setLatexFontSize(labelBgSubtracted, 240.0, 8.0);
	    labelBgSubtracted->Draw();

	    TString nameLowmassScale = "pp_HT1_45_lowmass_scale";
	    TString titleLowmassScale = "p+p HighTower-1 4 < pT < 5 GeV/c, low mass background scale";
	    TCanvas *cLowmassScale = new TCanvas(nameLowmassScale.Data(), titleLowmassScale.Data(), 0, 0, 600, 300);
	    //cLowmassScale->SetCanvasSize(600, 400);
	    cLowmassScale->cd();
	    TH1F *hSigLowmassScale = new TH1F(*hSig);
	    hSigLowmassScale->SetNameTitle("same", ";#font[12]{M}_{#gamma#gamma} [GeV/#font[12]{c}^{2}];Counts  ");
	    hSigLowmassScale->GetXaxis()->SetRangeUser(0, 0.35);
	    setHistFontSize(hSigLowmassScale, 240.0, 8.0);
	    hSigLowmassScale->Draw();
	    TLine *zeroLineSigLowmassScale = new TLine(0, 0, 0.35, 0);
	    zeroLineSigLowmassScale->Draw();
	    hBgLowmass->DrawCopy("HIST SAME");
	    //TH1F *hBgLowmassScale = new TH1F(*hBgLowmass);
	    //hBgLowmassScale->SetName("lowmass_scale");
	    //hBgLowmassScale->Scale(1.0 / 1.5);
	    //hBgLowmassScale->SetFillColor(17);
	    //hBgLowmassScale->DrawCopy("HIST SAME");
	    TArrow *arrow = new TArrow(0.15,1100,0.0825,1100,0.03,"|>");
            arrow->SetFillColor(1);
            arrow->SetFillStyle(1001);
	    arrow->SetAngle(45);
	    arrow->Draw();
	    arrow = new TArrow(0.22,1100,0.265,1100,0.03,"|>");
	    arrow->SetFillColor(1);
	    arrow->SetFillStyle(1001);
	    arrow->SetAngle(45);
	    arrow->Draw();
	    arrow = new TArrow(0.0825,1000,0.0825,1200,0,"");
	    arrow->SetFillColor(1);
	    arrow->SetFillStyle(1001);
	    arrow->SetAngle(0);
	    arrow->Draw();
	    arrow = new TArrow(0.265,1000,0.265,1200,0,"");
	    arrow->SetFillColor(1);
	    arrow->SetFillStyle(1001);
	    arrow->SetAngle(0);
	    arrow->Draw();
	    TLatex *   tex = new TLatex(0.1537283,1078,"#pi^{0}^{} window");
    	    tex->SetTextFont(132);
    	    tex->SetTextSize(0.04700855);
	    tex->SetLineWidth(2);
	    setLatexFontSize(tex, 240.0, 8.0);
	    tex->Draw();
	    
	    }
	}
}

//----------------------------------------------------------------------
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
	FILE_NAME(filenameSimunbarMB##NAME, DATA_DIR, "/simulation_dAu2003_nbar_data_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarHT1##NAME, DATA_DIR, "/simulation_dAu2003_nbar_data_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarHT2##NAME, DATA_DIR, "/simulation_dAu2003_nbar_data_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightMB##NAME, DATA_DIR, "/simulation_dAu2003_nbar_data_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightHT1##NAME, DATA_DIR, "/simulation_dAu2003_nbar_data_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightHT2##NAME, DATA_DIR, "/simulation_dAu2003_nbar_data_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarSimMB##NAME, DATA_DIR, "/simulation_dAu2003_nbar_pi0_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarSimHT1##NAME, DATA_DIR, "/simulation_dAu2003_nbar_pi0_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarSimHT2##NAME, DATA_DIR, "/simulation_dAu2003_nbar_pi0_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightSimMB##NAME, DATA_DIR, "/simulation_dAu2003_nbar_pi0_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightSimHT1##NAME, DATA_DIR, "/simulation_dAu2003_nbar_pi0_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightSimHT2##NAME, DATA_DIR, "/simulation_dAu2003_nbar_pi0_HT2_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarSimEtaMB##NAME, DATA_DIR, "/simulation_dAu2003_nbar_eta_MB_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarSimEtaHT1##NAME, DATA_DIR, "/simulation_dAu2003_nbar_eta_HT1_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarSimEtaHT2##NAME, DATA_DIR, "/simulation_dAu2003_nbar_eta_HT2_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightSimEtaMB##NAME, DATA_DIR, "/simulation_dAu2003_nbar_eta_MB_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightSimEtaHT1##NAME, DATA_DIR, "/simulation_dAu2003_nbar_eta_HT1_noweight_" CENTRALITY ".root") \
	FILE_NAME(filenameSimunbarNoweightSimEtaHT2##NAME, DATA_DIR, "/simulation_dAu2003_nbar_eta_HT2_noweight_" CENTRALITY ".root") \
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
	TDataProcessorPool *poolSimu1gammaNoweightSimEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_gamma_eta_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarMB##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_MB_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarHT1##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_HT1_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarHT2##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_HT2_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightMB##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightHT1##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightHT2##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarSimMB##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_pi0_MB_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarSimHT1##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_pi0_HT1_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarSimHT2##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_pi0_HT2_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightSimMB##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_pi0_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightSimHT1##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_pi0_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightSimHT2##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_pi0_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event pi0 - HT2 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarSimEtaMB##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_eta_MB_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - MB - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarSimEtaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_eta_HT1_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT1 - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarSimEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_eta_HT2_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT2 - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightSimEtaMB##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_eta_MB_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - MB - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightSimEtaHT1##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_eta_HT1_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT1 - no weight - " CENTRALITY); \
	TDataProcessorPool *poolSimunbarNoweightSimEtaHT2##NAME = new TDataProcessorPool("simulation_dAu2003_nbar_eta_HT2_noweight_" CENTRALITY, "Simulation - dAu2003 - 1 photon/event eta - HT2 - no weight - " CENTRALITY);

#define DATA_FILE_LOAD_DAU_DATA(NAME, CENTRALITY, READ, SETTINGS, RESULTS) \
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
	}

#define DATA_FILE_LOAD_DAU_SIM(NAME, CENTRALITY, READ, SETTINGS, RESULTS) \
	{ \
	if (allSettings.useDAuEff && ((allSettings.useDAuEffCentralities && READ) || (!allSettings.useDAuEffCentralities && (TString("nocentral") == CENTRALITY)))) { \
		if (allSettings.useDAuEmbedding) { \
			if (allSettings.useDAuEffWeight) { \
			    poolEmbed0##NAME->processFile(filenameEmbed##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbed0##NAME, dataFolder); \
    			} \
			if (!allSettings.useDAuEffWeight) { \
			    poolEmbedNoweight0##NAME->processFile(filenameEmbedNoweight##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbedNoweight0##NAME, dataFolder); \
			} \
		} else { \
			if (allSettings.useDAuEffWeight) { \
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
	}

#define DATA_FILE_LOAD_DAU_SIMETA(NAME, CENTRALITY, READ, SETTINGS, RESULTS) \
	{ \
	if (allSettings.useDAuEffEta && ((allSettings.useDAuEffCentralitiesEta && READ) || (!allSettings.useDAuEffCentralitiesEta && (TString("nocentral") == CENTRALITY)))) { \
		if (allSettings.useDAuEmbeddingEta) { \
			if (allSettings.useDAuEffWeightEta) { \
			    poolEmbedEta0##NAME->processFile(filenameEmbedEta##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbedEta0##NAME, dataFolder); \
    			} \
			if (!allSettings.useDAuEffWeightEta) { \
			    poolEmbedNoweightEta0##NAME->processFile(filenameEmbedNoweightEta##NAME); \
    			    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolEmbedNoweightEta0##NAME, dataFolder); \
			} \
		} else { \
			if (allSettings.useDAuEffWeightEta) { \
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
	}

#define DATA_FILE_LOAD_DAU_SIM1GAMMA(NAME, CENTRALITY, READ, SETTINGS, RESULTS) \
	{ \
	if (allSettings.useDAu1gamma && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (TString("nocentral") == CENTRALITY)))) { \
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
	if (allSettings.useDAu1gammaSim && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (TString("nocentral") == CENTRALITY)))) { \
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
	if (allSettings.useDAu1gammaSimEta && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (TString("nocentral") == CENTRALITY)))) { \
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

#define DATA_FILE_LOAD_DAU_SIMNBAR(NAME, CENTRALITY, READ, SETTINGS, RESULTS) \
	{ \
	if (allSettings.useDAunbar && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (TString("nocentral") == CENTRALITY)))) { \
		poolSimunbarMB##NAME->processFile(filenameSimunbarMB##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarMB##NAME, dataFolder); \
		poolSimunbarHT1##NAME->processFile(filenameSimunbarHT1##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarHT1##NAME, dataFolder); \
		poolSimunbarHT2##NAME->processFile(filenameSimunbarHT2##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarHT2##NAME, dataFolder); \
		if (!allSettings.useDAuEffWeight) { \
		    poolSimunbarNoweightMB##NAME->processFile(filenameSimunbarNoweightMB##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightMB##NAME, dataFolder); \
		    poolSimunbarNoweightHT1##NAME->processFile(filenameSimunbarNoweightHT1##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightHT1##NAME, dataFolder); \
		    poolSimunbarNoweightHT2##NAME->processFile(filenameSimunbarNoweightHT2##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightHT2##NAME, dataFolder); \
		} \
	} \
	if (allSettings.useDAunbarSim && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (TString("nocentral") == CENTRALITY)))) { \
		poolSimunbarSimMB##NAME->processFile(filenameSimunbarSimMB##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarSimMB##NAME, dataFolder); \
		poolSimunbarSimHT1##NAME->processFile(filenameSimunbarSimHT1##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarSimHT1##NAME, dataFolder); \
		poolSimunbarSimHT2##NAME->processFile(filenameSimunbarSimHT2##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarSimHT2##NAME, dataFolder); \
		if (!allSettings.useDAuEffWeight) { \
		    poolSimunbarNoweightSimMB##NAME->processFile(filenameSimunbarNoweightSimMB##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightSimMB##NAME, dataFolder); \
		    poolSimunbarNoweightSimHT1##NAME->processFile(filenameSimunbarNoweightSimHT1##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightSimHT1##NAME, dataFolder); \
		    poolSimunbarNoweightSimHT2##NAME->processFile(filenameSimunbarNoweightSimHT2##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightSimHT2##NAME, dataFolder); \
		} \
	} \
	if (allSettings.useDAunbarSimEta && ((allSettings.useDAu1gammaCentralities && READ) || (!allSettings.useDAu1gammaCentralities && (TString("nocentral") == CENTRALITY)))) { \
		poolSimunbarSimEtaMB##NAME->processFile(filenameSimunbarSimEtaMB##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarSimEtaMB##NAME, dataFolder); \
		poolSimunbarSimEtaHT1##NAME->processFile(filenameSimunbarSimEtaHT1##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarSimEtaHT1##NAME, dataFolder); \
		poolSimunbarSimEtaHT2##NAME->processFile(filenameSimunbarSimEtaHT2##NAME); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarSimEtaHT2##NAME, dataFolder); \
		if (!allSettings.useDAuEffWeight) { \
		    poolSimunbarNoweightSimEtaMB##NAME->processFile(filenameSimunbarNoweightSimEtaMB##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightSimEtaMB##NAME, dataFolder); \
		    poolSimunbarNoweightSimEtaHT1##NAME->processFile(filenameSimunbarNoweightSimEtaHT1##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightSimEtaHT1##NAME, dataFolder); \
		    poolSimunbarNoweightSimEtaHT2##NAME->processFile(filenameSimunbarNoweightSimEtaHT2##NAME); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolSimunbarNoweightSimEtaHT2##NAME, dataFolder); \
		} \
	} \
	}

#define DATA_FILE_DELETE_DAU_DATA(NAME) \
	{ \
	removePool(poolMB##NAME, true, true); \
	removePool(poolHT1##NAME, true, true); \
	removePool(poolHT2##NAME, true, true); \
	}

#define DATA_FILE_DELETE_DAU_SIM(NAME) \
	{ \
	removePool(poolEmbed0##NAME, true, true); \
	removePool(poolEmbedMB##NAME, true, true); \
	removePool(poolEmbedHT1##NAME, true, true); \
	removePool(poolEmbedHT2##NAME, true, true); \
	removePool(poolEmbedNoweight0##NAME, true, true); \
	removePool(poolEmbedNoweightMB##NAME, true, true); \
	removePool(poolEmbedNoweightHT1##NAME, true, true); \
	removePool(poolEmbedNoweightHT2##NAME, true, true); \
	removePool(poolSimuMB##NAME, true, true); \
	removePool(poolSimuHT1##NAME, true, true); \
	removePool(poolSimuHT2##NAME, true, true); \
	removePool(poolSimuNoweightMB##NAME, true, true); \
	removePool(poolSimuNoweightHT1##NAME, true, true); \
	removePool(poolSimuNoweightHT2##NAME, true, true); \
	}

#define DATA_FILE_DELETE_DAU_SIMETA(NAME) \
	{ \
	removePool(poolEmbedEta0##NAME, true, true); \
	removePool(poolEmbedEtaMB##NAME, true, true); \
	removePool(poolEmbedEtaHT1##NAME, true, true); \
	removePool(poolEmbedEtaHT2##NAME, true, true); \
	removePool(poolEmbedNoweightEta0##NAME, true, true); \
	removePool(poolEmbedNoweightEtaMB##NAME, true, true); \
	removePool(poolEmbedNoweightEtaHT1##NAME, true, true); \
	removePool(poolEmbedNoweightEtaHT2##NAME, true, true); \
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
	}

#define DATA_FILE_DELETE_DAU_SIM1GAMMA(NAME) \
	{ \
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

#define DATA_FILE_DELETE_DAU_SIMNBAR(NAME) \
	{ \
	removePool(poolSimunbarMB##NAME, true, true); \
	removePool(poolSimunbarHT1##NAME, true, true); \
	removePool(poolSimunbarHT2##NAME, true, true); \
	removePool(poolSimunbarNoweightMB##NAME, true, true); \
	removePool(poolSimunbarNoweightHT1##NAME, true, true); \
	removePool(poolSimunbarNoweightHT2##NAME, true, true); \
	removePool(poolSimunbarSimMB##NAME, true, true); \
	removePool(poolSimunbarSimHT1##NAME, true, true); \
	removePool(poolSimunbarSimHT2##NAME, true, true); \
	removePool(poolSimunbarNoweightSimMB##NAME, true, true); \
	removePool(poolSimunbarNoweightSimHT1##NAME, true, true); \
	removePool(poolSimunbarNoweightSimHT2##NAME, true, true); \
	removePool(poolSimunbarSimEtaMB##NAME, true, true); \
	removePool(poolSimunbarSimEtaHT1##NAME, true, true); \
	removePool(poolSimunbarSimEtaHT2##NAME, true, true); \
	removePool(poolSimunbarNoweightSimEtaMB##NAME, true, true); \
	removePool(poolSimunbarNoweightSimEtaHT1##NAME, true, true); \
	removePool(poolSimunbarNoweightSimEtaHT2##NAME, true, true); \
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
	FILE_NAME(filenamePPSimunbarMB, DATA_DIR, "/simulation_pp2005_nbar_data_MB.root") \
	FILE_NAME(filenamePPSimunbarHT1, DATA_DIR, "/simulation_pp2005_nbar_data_HT1.root") \
	FILE_NAME(filenamePPSimunbarHT2, DATA_DIR, "/simulation_pp2005_nbar_data_HT2.root") \
	FILE_NAME(filenamePPSimunbarNoweightMB, DATA_DIR, "/simulation_pp2005_nbar_data_MB_noweight.root") \
	FILE_NAME(filenamePPSimunbarNoweightHT1, DATA_DIR, "/simulation_pp2005_nbar_data_HT1_noweight.root") \
	FILE_NAME(filenamePPSimunbarNoweightHT2, DATA_DIR, "/simulation_pp2005_nbar_data_HT2_noweight.root") \
	FILE_NAME(filenamePPSimunbarSimMB, DATA_DIR, "/simulation_pp2005_nbar_pi0_MB.root") \
	FILE_NAME(filenamePPSimunbarSimHT1, DATA_DIR, "/simulation_pp2005_nbar_pi0_HT1.root") \
	FILE_NAME(filenamePPSimunbarSimHT2, DATA_DIR, "/simulation_pp2005_nbar_pi0_HT2.root") \
	FILE_NAME(filenamePPSimunbarNoweightSimMB, DATA_DIR, "/simulation_pp2005_nbar_pi0_MB_noweight.root") \
	FILE_NAME(filenamePPSimunbarNoweightSimHT1, DATA_DIR, "/simulation_pp2005_nbar_pi0_HT1_noweight.root") \
	FILE_NAME(filenamePPSimunbarNoweightSimHT2, DATA_DIR, "/simulation_pp2005_nbar_pi0_HT2_noweight.root") \
	FILE_NAME(filenamePPSimunbarSimEtaMB, DATA_DIR, "/simulation_pp2005_nbar_eta_MB.root") \
	FILE_NAME(filenamePPSimunbarSimEtaHT1, DATA_DIR, "/simulation_pp2005_nbar_eta_HT1.root") \
	FILE_NAME(filenamePPSimunbarSimEtaHT2, DATA_DIR, "/simulation_pp2005_nbar_eta_HT2.root") \
	FILE_NAME(filenamePPSimunbarNoweightSimEtaMB, DATA_DIR, "/simulation_pp2005_nbar_eta_MB_noweight.root") \
	FILE_NAME(filenamePPSimunbarNoweightSimEtaHT1, DATA_DIR, "/simulation_pp2005_nbar_eta_HT1_noweight.root") \
	FILE_NAME(filenamePPSimunbarNoweightSimEtaHT2, DATA_DIR, "/simulation_pp2005_nbar_eta_HT2_noweight.root") \
	FILE_NAME(filenamePPPythiaMB, DATA_DIR, "/pythia_pp2005_pi0_MB.root") \
	FILE_NAME(filenamePPPythiaHT1, DATA_DIR, "/pythia_pp2005_pi0_HT1.root") \
	FILE_NAME(filenamePPPythiaHT2, DATA_DIR, "/pythia_pp2005_pi0_HT2.root") \
	FILE_NAME(filenamePPPythiaEtaMB, DATA_DIR, "/pythia_pp2005_eta_MB.root") \
	FILE_NAME(filenamePPPythiaEtaHT1, DATA_DIR, "/pythia_pp2005_eta_HT1.root") \
	FILE_NAME(filenamePPPythiaEtaHT2, DATA_DIR, "/pythia_pp2005_eta_HT2.root") \
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
	TDataProcessorPool *poolPPSimu1gammaNoweightSimEtaHT2 = new TDataProcessorPool("simulation_pp2005_gamma_eta_HT2_noweight", "Simulation - pp2005 - 1 photon/event eta - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimunbarMB = new TDataProcessorPool("simulation_pp2005_nbar_MB", "Simulation - pp2005 - 1 photon/event - MB"); \
	TDataProcessorPool *poolPPSimunbarHT1 = new TDataProcessorPool("simulation_pp2005_nbar_HT1", "Simulation - pp2005 - 1 photon/event - HT1"); \
	TDataProcessorPool *poolPPSimunbarHT2 = new TDataProcessorPool("simulation_pp2005_nbar_HT2", "Simulation - pp2005 - 1 photon/event - HT2"); \
	TDataProcessorPool *poolPPSimunbarNoweightMB = new TDataProcessorPool("simulation_pp2005_nbar_MB_noweight", "Simulation - pp2005 - 1 photon/event - MB - no weight"); \
	TDataProcessorPool *poolPPSimunbarNoweightHT1 = new TDataProcessorPool("simulation_pp2005_nbar_HT1_noweight", "Simulation - pp2005 - 1 photon/event - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimunbarNoweightHT2 = new TDataProcessorPool("simulation_pp2005_nbar_HT2_noweight", "Simulation - pp2005 - 1 photon/event - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimunbarSimMB = new TDataProcessorPool("simulation_pp2005_nbar_pi0_MB", "Simulation - pp2005 - 1 photon/event pi0 - MB"); \
	TDataProcessorPool *poolPPSimunbarSimHT1 = new TDataProcessorPool("simulation_pp2005_nbar_pi0_HT1", "Simulation - pp2005 - 1 photon/event pi0 - HT1"); \
	TDataProcessorPool *poolPPSimunbarSimHT2 = new TDataProcessorPool("simulation_pp2005_nbar_pi0_HT2", "Simulation - pp2005 - 1 photon/event pi0 - HT2"); \
	TDataProcessorPool *poolPPSimunbarNoweightSimMB = new TDataProcessorPool("simulation_pp2005_nbar_pi0_MB_noweight", "Simulation - pp2005 - 1 photon/event pi0 - MB - no weight"); \
	TDataProcessorPool *poolPPSimunbarNoweightSimHT1 = new TDataProcessorPool("simulation_pp2005_nbar_pi0_HT1_noweight", "Simulation - pp2005 - 1 photon/event pi0 - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimunbarNoweightSimHT2 = new TDataProcessorPool("simulation_pp2005_nbar_pi0_HT2_noweight", "Simulation - pp2005 - 1 photon/event pi0 - HT2 - no weight"); \
	TDataProcessorPool *poolPPSimunbarSimEtaMB = new TDataProcessorPool("simulation_pp2005_nbar_eta_MB", "Simulation - pp2005 - 1 photon/event eta - MB"); \
	TDataProcessorPool *poolPPSimunbarSimEtaHT1 = new TDataProcessorPool("simulation_pp2005_nbar_eta_HT1", "Simulation - pp2005 - 1 photon/event eta - HT1"); \
	TDataProcessorPool *poolPPSimunbarSimEtaHT2 = new TDataProcessorPool("simulation_pp2005_nbar_eta_HT2", "Simulation - pp2005 - 1 photon/event eta - HT2"); \
	TDataProcessorPool *poolPPSimunbarNoweightSimEtaMB = new TDataProcessorPool("simulation_pp2005_nbar_eta_MB_noweight", "Simulation - pp2005 - 1 photon/event eta - MB - no weight"); \
	TDataProcessorPool *poolPPSimunbarNoweightSimEtaHT1 = new TDataProcessorPool("simulation_pp2005_nbar_eta_HT1_noweight", "Simulation - pp2005 - 1 photon/event eta - HT1 - no weight"); \
	TDataProcessorPool *poolPPSimunbarNoweightSimEtaHT2 = new TDataProcessorPool("simulation_pp2005_nbar_eta_HT2_noweight", "Simulation - pp2005 - 1 photon/event eta - HT2 - no weight"); \
	TDataProcessorPool *poolPPPythiaMB = new TDataProcessorPool("pythia_pp2005_pi0_MB", "Pythia - pp2005 - pi0 - MB"); \
	TDataProcessorPool *poolPPPythiaHT1 = new TDataProcessorPool("pythia_pp2005_pi0_HT1", "Pythia - pp2005 - pi0 - HT1"); \
	TDataProcessorPool *poolPPPythiaHT2 = new TDataProcessorPool("pythia_pp2005_pi0_HT2", "Pythia - pp2005 - pi0 - HT2"); \
	TDataProcessorPool *poolPPPythiaEtaMB = new TDataProcessorPool("pythia_pp2005_eta_MB", "Pythia - pp2005 - eta - MB"); \
	TDataProcessorPool *poolPPPythiaEtaHT1 = new TDataProcessorPool("pythia_pp2005_eta_HT1", "Pythia - pp2005 - eta - HT1"); \
	TDataProcessorPool *poolPPPythiaEtaHT2 = new TDataProcessorPool("pythia_pp2005_eta_HT2", "Pythia - pp2005 - eta - HT2");

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
			if (allSettings.usePPEffWeight) { \
			    if (allSettings.usePPEffPythia) { \
				poolPPPythiaMB->processFile(filenamePPPythiaMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPPythiaMB, dataFolder); \
				poolPPPythiaHT1->processFile(filenamePPPythiaHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPPythiaHT1, dataFolder); \
				poolPPPythiaHT2->processFile(filenamePPPythiaHT2); \
   				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPPythiaHT2, dataFolder); \
			    } else { \
				poolPPSimuMB->processFile(filenamePPSimuMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuMB, dataFolder); \
				poolPPSimuHT1->processFile(filenamePPSimuHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuHT1, dataFolder); \
				poolPPSimuHT2->processFile(filenamePPSimuHT2); \
   				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimuHT2, dataFolder); \
			    } \
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
			if (allSettings.usePPEffWeightEta) { \
			    if (allSettings.usePPEffPythia) { \
				poolPPPythiaEtaMB->processFile(filenamePPPythiaEtaMB); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPPythiaEtaMB, dataFolder); \
				poolPPPythiaEtaHT1->processFile(filenamePPPythiaEtaHT1); \
    				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPPythiaEtaHT1, dataFolder); \
				poolPPPythiaEtaHT2->processFile(filenamePPPythiaEtaHT2); \
   				if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPPythiaEtaHT2, dataFolder); \
			    } else { \
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
	if (allSettings.usePPnbar) { \
		poolPPSimunbarMB->processFile(filenamePPSimunbarMB); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarMB, dataFolder); \
		poolPPSimunbarHT1->processFile(filenamePPSimunbarHT1); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarHT1, dataFolder); \
		poolPPSimunbarHT2->processFile(filenamePPSimunbarHT2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarHT2, dataFolder); \
		if (!allSettings.usePPEffWeight) { \
		    poolPPSimunbarNoweightMB->processFile(filenamePPSimunbarNoweightMB); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightMB, dataFolder); \
		    poolPPSimunbarNoweightHT1->processFile(filenamePPSimunbarNoweightHT1); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightHT1, dataFolder); \
		    poolPPSimunbarNoweightHT2->processFile(filenamePPSimunbarNoweightHT2); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightHT2, dataFolder); \
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
	if (allSettings.usePPnbarSim) { \
		poolPPSimunbarSimMB->processFile(filenamePPSimunbarSimMB); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarSimMB, dataFolder); \
		poolPPSimunbarSimHT1->processFile(filenamePPSimunbarSimHT1); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarSimHT1, dataFolder); \
		poolPPSimunbarSimHT2->processFile(filenamePPSimunbarSimHT2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarSimHT2, dataFolder); \
		if (!allSettings.usePPEffWeight) { \
		    poolPPSimunbarNoweightSimMB->processFile(filenamePPSimunbarNoweightSimMB); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightSimMB, dataFolder); \
		    poolPPSimunbarNoweightSimHT1->processFile(filenamePPSimunbarNoweightSimHT1); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightSimHT1, dataFolder); \
		    poolPPSimunbarNoweightSimHT2->processFile(filenamePPSimunbarNoweightSimHT2); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightSimHT2, dataFolder); \
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
	if (allSettings.usePPnbarSimEta) { \
		poolPPSimunbarSimEtaMB->processFile(filenamePPSimunbarSimEtaMB); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarSimEtaMB, dataFolder); \
		poolPPSimunbarSimEtaHT1->processFile(filenamePPSimunbarSimEtaHT1); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarSimEtaHT1, dataFolder); \
		poolPPSimunbarSimEtaHT2->processFile(filenamePPSimunbarSimEtaHT2); \
    		if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarSimEtaHT2, dataFolder); \
		if (!allSettings.usePPEffWeight) { \
		    poolPPSimunbarNoweightSimEtaMB->processFile(filenamePPSimunbarNoweightSimEtaMB); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightSimEtaMB, dataFolder); \
		    poolPPSimunbarNoweightSimEtaHT1->processFile(filenamePPSimunbarNoweightSimEtaHT1); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightSimEtaHT1, dataFolder); \
		    poolPPSimunbarNoweightSimEtaHT2->processFile(filenamePPSimunbarNoweightSimEtaHT2); \
    		    if (dataFolder && allSettings.browseAnalysis) memberInspector.Inspect(poolPPSimunbarNoweightSimEtaHT2, dataFolder); \
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
	removePool(poolPPSimunbarMB, true, true);\
	removePool(poolPPSimunbarHT1, true, true);\
	removePool(poolPPSimunbarHT2, true, true);\
	removePool(poolPPSimunbarNoweightMB, true, true);\
	removePool(poolPPSimunbarNoweightHT1, true, true);\
	removePool(poolPPSimunbarNoweightHT2, true, true);\
	removePool(poolPPSimunbarSimMB, true, true);\
	removePool(poolPPSimunbarSimHT1, true, true);\
	removePool(poolPPSimunbarSimHT2, true, true);\
	removePool(poolPPSimunbarNoweightSimMB, true, true);\
	removePool(poolPPSimunbarNoweightSimHT1, true, true);\
	removePool(poolPPSimunbarNoweightSimHT2, true, true);\
	removePool(poolPPSimunbarSimEtaMB, true, true);\
	removePool(poolPPSimunbarSimEtaHT1, true, true);\
	removePool(poolPPSimunbarSimEtaHT2, true, true);\
	removePool(poolPPSimunbarNoweightSimEtaMB, true, true);\
	removePool(poolPPSimunbarNoweightSimEtaHT1, true, true);\
	removePool(poolPPSimunbarNoweightSimEtaHT2, true, true); \
	removePool(poolPPPythiaMB, true, true);\
	removePool(poolPPPythiaHT1, true, true);\
	removePool(poolPPPythiaHT2, true, true);\
	removePool(poolPPPythiaEtaMB, true, true);\
	removePool(poolPPPythiaEtaHT1, true, true);\
	removePool(poolPPPythiaEtaHT2, true, true);\
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

	TDataProcessorPool *poolnbarSimMBNo = allSettings.useDAuEffWeight ? (poolSimunbarSimMBNocentral) : (poolSimunbarNoweightSimMBNocentral);
	TDataProcessorPool *poolnbarSimMBAll = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimMBAllcentral : poolSimunbarSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimMBAllcentral : poolSimunbarNoweightSimMBNocentral));
	TDataProcessorPool *poolnbarSimMBMostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimMBMostcentral : poolSimunbarSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimMBMostcentral : poolSimunbarNoweightSimMBNocentral));
	TDataProcessorPool *poolnbarSimMBMidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimMBMidcentral : poolSimunbarSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimMBMidcentral : poolSimunbarNoweightSimMBNocentral));
	TDataProcessorPool *poolnbarSimMBMostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimMBMostperipheral : poolSimunbarSimMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimMBMostperipheral : poolSimunbarNoweightSimMBNocentral));

	TDataProcessorPool *poolnbarSimHT1No = allSettings.useDAuEffWeight ? (poolSimunbarSimHT1Nocentral) : (poolSimunbarNoweightSimHT1Nocentral);
	TDataProcessorPool *poolnbarSimHT1All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT1Allcentral : poolSimunbarSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT1Allcentral : poolSimunbarNoweightSimHT1Nocentral));
	TDataProcessorPool *poolnbarSimHT1MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT1Mostcentral : poolSimunbarSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT1Mostcentral : poolSimunbarNoweightSimHT1Nocentral));
	TDataProcessorPool *poolnbarSimHT1MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT1Midcentral : poolSimunbarSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT1Midcentral : poolSimunbarNoweightSimHT1Nocentral));
	TDataProcessorPool *poolnbarSimHT1MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT1Mostperipheral : poolSimunbarSimHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT1Mostperipheral : poolSimunbarNoweightSimHT1Nocentral));

	TDataProcessorPool *poolnbarSimHT2No = allSettings.useDAuEffWeight ? (poolSimunbarSimHT2Nocentral) : (poolSimunbarNoweightSimHT2Nocentral);
	TDataProcessorPool *poolnbarSimHT2All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT2Allcentral : poolSimunbarSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT2Allcentral : poolSimunbarNoweightSimHT2Nocentral));
	TDataProcessorPool *poolnbarSimHT2MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT2Mostcentral : poolSimunbarSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT2Mostcentral : poolSimunbarNoweightSimHT2Nocentral));
	TDataProcessorPool *poolnbarSimHT2MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT2Midcentral : poolSimunbarSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT2Midcentral : poolSimunbarNoweightSimHT2Nocentral));
	TDataProcessorPool *poolnbarSimHT2MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimHT2Mostperipheral : poolSimunbarSimHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimHT2Mostperipheral : poolSimunbarNoweightSimHT2Nocentral));

	TDataProcessorPool *poolnbarSimEtaMBNo = allSettings.useDAuEffWeight ? (poolSimunbarSimEtaMBNocentral) : (poolSimunbarNoweightSimEtaMBNocentral);
	TDataProcessorPool *poolnbarSimEtaMBAll = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaMBAllcentral : poolSimunbarSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaMBAllcentral : poolSimunbarNoweightSimEtaMBNocentral));
	TDataProcessorPool *poolnbarSimEtaMBMostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaMBMostcentral : poolSimunbarSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaMBMostcentral : poolSimunbarNoweightSimEtaMBNocentral));
	TDataProcessorPool *poolnbarSimEtaMBMidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaMBMidcentral : poolSimunbarSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaMBMidcentral : poolSimunbarNoweightSimEtaMBNocentral));
	TDataProcessorPool *poolnbarSimEtaMBMostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaMBMostperipheral : poolSimunbarSimEtaMBNocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaMBMostperipheral : poolSimunbarNoweightSimEtaMBNocentral));

	TDataProcessorPool *poolnbarSimEtaHT1No = allSettings.useDAuEffWeight ? (poolSimunbarSimEtaHT1Nocentral) : (poolSimunbarNoweightSimEtaHT1Nocentral);
	TDataProcessorPool *poolnbarSimEtaHT1All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT1Allcentral : poolSimunbarSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT1Allcentral : poolSimunbarNoweightSimEtaHT1Nocentral));
	TDataProcessorPool *poolnbarSimEtaHT1MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT1Mostcentral : poolSimunbarSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT1Mostcentral : poolSimunbarNoweightSimEtaHT1Nocentral));
	TDataProcessorPool *poolnbarSimEtaHT1MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT1Midcentral : poolSimunbarSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT1Midcentral : poolSimunbarNoweightSimEtaHT1Nocentral));
	TDataProcessorPool *poolnbarSimEtaHT1MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT1Mostperipheral : poolSimunbarSimEtaHT1Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT1Mostperipheral : poolSimunbarNoweightSimEtaHT1Nocentral));

	TDataProcessorPool *poolnbarSimEtaHT2No = allSettings.useDAuEffWeight ? (poolSimunbarSimEtaHT2Nocentral) : (poolSimunbarNoweightSimEtaHT2Nocentral);
	TDataProcessorPool *poolnbarSimEtaHT2All = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT2Allcentral : poolSimunbarSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT2Allcentral : poolSimunbarNoweightSimEtaHT2Nocentral));
	TDataProcessorPool *poolnbarSimEtaHT2MostC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT2Mostcentral : poolSimunbarSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT2Mostcentral : poolSimunbarNoweightSimEtaHT2Nocentral));
	TDataProcessorPool *poolnbarSimEtaHT2MidC = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT2Midcentral : poolSimunbarSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT2Midcentral : poolSimunbarNoweightSimEtaHT2Nocentral));
	TDataProcessorPool *poolnbarSimEtaHT2MostP = allSettings.useDAuEffWeight ? ((allSettings.useDAu1gammaCentralities ? poolSimunbarSimEtaHT2Mostperipheral : poolSimunbarSimEtaHT2Nocentral)) : ((allSettings.useDAu1gammaCentralities ? poolSimunbarNoweightSimEtaHT2Mostperipheral : poolSimunbarNoweightSimEtaHT2Nocentral));

	TDataProcessorPool *poolnbarMBNo = poolSimunbarMBNocentral;
	TDataProcessorPool *poolnbarMBAll = allSettings.useDAu1gammaCentralities ? poolSimunbarMBAllcentral : poolSimunbarMBNocentral;
	TDataProcessorPool *poolnbarMBMostC = allSettings.useDAu1gammaCentralities ? poolSimunbarMBMostcentral : poolSimunbarMBNocentral;
	TDataProcessorPool *poolnbarMBMidC = allSettings.useDAu1gammaCentralities ? poolSimunbarMBMidcentral : poolSimunbarMBNocentral;
	TDataProcessorPool *poolnbarMBMostP = allSettings.useDAu1gammaCentralities ? poolSimunbarMBMostperipheral : poolSimunbarMBNocentral;

	TDataProcessorPool *poolnbarHT1No = poolSimunbarHT1Nocentral;
	TDataProcessorPool *poolnbarHT1All = allSettings.useDAu1gammaCentralities ? poolSimunbarHT1Allcentral : poolSimunbarHT1Nocentral;
	TDataProcessorPool *poolnbarHT1MostC = allSettings.useDAu1gammaCentralities ? poolSimunbarHT1Mostcentral : poolSimunbarHT1Nocentral;
	TDataProcessorPool *poolnbarHT1MidC = allSettings.useDAu1gammaCentralities ? poolSimunbarHT1Midcentral : poolSimunbarHT1Nocentral;
	TDataProcessorPool *poolnbarHT1MostP = allSettings.useDAu1gammaCentralities ? poolSimunbarHT1Mostperipheral : poolSimunbarHT1Nocentral;

	TDataProcessorPool *poolnbarHT2No = poolSimunbarHT2Nocentral;
	TDataProcessorPool *poolnbarHT2All = allSettings.useDAu1gammaCentralities ? poolSimunbarHT2Allcentral : poolSimunbarHT2Nocentral;
	TDataProcessorPool *poolnbarHT2MostC = allSettings.useDAu1gammaCentralities ? poolSimunbarHT2Mostcentral : poolSimunbarHT2Nocentral;
	TDataProcessorPool *poolnbarHT2MidC = allSettings.useDAu1gammaCentralities ? poolSimunbarHT2Midcentral : poolSimunbarHT2Nocentral;
	TDataProcessorPool *poolnbarHT2MostP = allSettings.useDAu1gammaCentralities ? poolSimunbarHT2Mostperipheral : poolSimunbarHT2Nocentral;

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

	TDataProcessorPool *poolPPEffMB = allSettings.usePPEffWeight ? (allSettings.usePPEffPythia ? poolPPPythiaMB : poolPPSimuMB) : poolPPSimuNoweightMB;
	TDataProcessorPool *poolPPEffHT1 = allSettings.usePPEffWeight ? (allSettings.usePPEffPythia ? poolPPPythiaHT1 : poolPPSimuHT1) : poolPPSimuNoweightHT1;
	TDataProcessorPool *poolPPEffHT2 = allSettings.usePPEffWeight ? (allSettings.usePPEffPythia ? poolPPPythiaHT2 : poolPPSimuHT2) : poolPPSimuNoweightHT2;
	TDataProcessorPool *poolPPEffEtaMB = allSettings.usePPEffWeightEta ? (allSettings.usePPEffPythia ? poolPPPythiaEtaMB : poolPPSimuEtaMB) : poolPPSimuNoweightEtaMB;
	TDataProcessorPool *poolPPEffEtaHT1 = allSettings.usePPEffWeightEta ? (allSettings.usePPEffPythia ? poolPPPythiaEtaHT1 : poolPPSimuEtaHT1) : poolPPSimuNoweightEtaHT1;
	TDataProcessorPool *poolPPEffEtaHT2 = allSettings.usePPEffWeightEta ? (allSettings.usePPEffPythia ? poolPPPythiaEtaHT2 : poolPPSimuEtaHT2) : poolPPSimuNoweightEtaHT2;
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
	TDataProcessorPool *poolPPnbarMB = poolPPSimunbarMB;
	TDataProcessorPool *poolPPnbarHT1 = poolPPSimunbarHT1;
	TDataProcessorPool *poolPPnbarHT2 = poolPPSimunbarHT2;
	TDataProcessorPool *poolPPnbarSimMB = allSettings.usePPEffWeight ? poolPPSimunbarSimMB : poolPPSimunbarNoweightSimMB;
	TDataProcessorPool *poolPPnbarSimHT1 = allSettings.usePPEffWeight ? poolPPSimunbarSimHT1 : poolPPSimunbarNoweightSimHT1;
	TDataProcessorPool *poolPPnbarSimHT2 = allSettings.usePPEffWeight ? poolPPSimunbarSimHT2 : poolPPSimunbarNoweightSimHT2;
	TDataProcessorPool *poolPPnbarSimEtaMB = allSettings.usePPEffWeight ? poolPPSimunbarSimEtaMB : poolPPSimunbarNoweightSimEtaMB;
	TDataProcessorPool *poolPPnbarSimEtaHT1 = allSettings.usePPEffWeight ? poolPPSimunbarSimEtaHT1 : poolPPSimunbarNoweightSimEtaHT1;
	TDataProcessorPool *poolPPnbarSimEtaHT2 = allSettings.usePPEffWeight ? poolPPSimunbarSimEtaHT2 : poolPPSimunbarNoweightSimEtaHT2;
	TDataProcessorPool *poolPPPtshiftMB = allSettings.usePPEffPythia ? poolPPPythiaMB : poolPPSimuMB;
	TDataProcessorPool *poolPPPtshiftHT1 = allSettings.usePPEffPythia ? poolPPPythiaHT1 : poolPPSimuHT1;
	TDataProcessorPool *poolPPPtshiftHT2 = allSettings.usePPEffPythia ? poolPPPythiaHT2 : poolPPSimuHT2;

	typedef const TDataProcessorPool *pcPool;

	if (allSettings.showSpectrumDAu || allSettings.showSpectrumDAuEta || allSettings.showRcp || allSettings.showRcpEta || allSettings.showRDA || allSettings.showRDAEta || allSettings.showRDAGamma) {
	    DATA_FILE_LOAD_DAU_SIM(Nocentral, "nocentral", true, allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
	    DATA_FILE_LOAD_DAU_SIMETA(Nocentral, "nocentral", true, allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
	    DATA_FILE_LOAD_DAU_SIM1GAMMA(Nocentral, "nocentral", true, allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
	    DATA_FILE_LOAD_DAU_SIMNBAR(Nocentral, "nocentral", true, allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
	    if (allSettings.showSpectrumDAu || allSettings.showSpectrumDAuEta || allSettings.showRDA || allSettings.showRDAEta || allSettings.showRDAGamma) {
		DATA_FILE_LOAD_DAU_DATA(Nocentral, "nocentral", true, allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
		pcPool poolsNo[] = {poolMBNocentral, poolHT1Nocentral, poolHT2Nocentral
		    , poolEffMBNo, poolEffHT1No, poolEffHT2No
		    , poolEffEtaMBNo, poolEffEtaHT1No, poolEffEtaHT2No
		    , poolEffEtabgMBNo, poolEffEtabgHT1No, poolEffEtabgHT2No
		    , pool1gammaMBNo, pool1gammaHT1No, pool1gammaHT2No
		    , pool1gammaSimMBNo, pool1gammaSimHT1No, pool1gammaSimHT2No
		    , pool1gammaSimEtaMBNo, pool1gammaSimEtaHT1No, pool1gammaSimEtaHT2No
		    , poolnbarMBNo, poolnbarHT1No, poolnbarHT2No
		    , poolnbarSimMBNo, poolnbarSimHT1No, poolnbarSimHT2No
		    , poolnbarSimEtaMBNo, poolnbarSimEtaHT1No, poolnbarSimEtaHT2No
		    , poolPtshiftMBNo, poolPtshiftHT1No, poolPtshiftHT2No};
		showResults(&poolsNo[0], allSettings.settingsDAuNoCentral, allResults.resultsDAuNoCentral);
		if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU_DATA(Nocentral);
		if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralities) DATA_FILE_DELETE_DAU_SIM(Nocentral);
		if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralitiesEta) DATA_FILE_DELETE_DAU_SIMETA(Nocentral);
		if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIM1GAMMA(Nocentral);
		if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIMNBAR(Nocentral);
	    }
	}

	if (allSettings.showRcp || allSettings.showRcpEta) {
	    if (allSettings.showAllCentral) {    
		DATA_FILE_LOAD_DAU_DATA(Allcentral, "0-100", true, allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		DATA_FILE_LOAD_DAU_SIM(Allcentral, "0-100", true, allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		DATA_FILE_LOAD_DAU_SIMETA(Allcentral, "0-100", true, allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		DATA_FILE_LOAD_DAU_SIM1GAMMA(Allcentral, "0-100", true, allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		DATA_FILE_LOAD_DAU_SIMNBAR(Allcentral, "0-100", true, allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		pcPool poolsAll[] = {poolMBAllcentral, poolHT1Allcentral, poolHT2Allcentral
    		    , poolEffMBAll, poolEffHT1All, poolEffHT2All
		    , poolEffEtaMBAll, poolEffEtaHT1All, poolEffEtaHT2All
		    , poolEffEtabgMBAll, poolEffEtabgHT1All, poolEffEtabgHT2All
		    , pool1gammaMBAll, pool1gammaHT1All, pool1gammaHT2All
		    , pool1gammaSimMBAll, pool1gammaSimHT1All, pool1gammaSimHT2All
		    , pool1gammaSimEtaMBAll, pool1gammaSimEtaHT1All, pool1gammaSimEtaHT2All
		    , poolnbarMBAll, poolnbarHT1All, poolnbarHT2All
		    , poolnbarSimMBAll, poolnbarSimHT1All, poolnbarSimHT2All
		    , poolnbarSimEtaMBAll, poolnbarSimEtaHT1All, poolnbarSimEtaHT2All
		    , poolPtshiftMBAll, poolPtshiftHT1All, poolPtshiftHT2All};
		showResults(&poolsAll[0], allSettings.settingsDAuAllCentral, allResults.resultsDAuAllCentral);
		if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU_DATA(Allcentral);
		if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralities) DATA_FILE_DELETE_DAU_SIM(Allcentral);
		if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralitiesEta) DATA_FILE_DELETE_DAU_SIMETA(Allcentral);
		if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIM1GAMMA(Allcentral);
		if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIMNBAR(Allcentral);
	    }

    	    DATA_FILE_LOAD_DAU_DATA(Mostcentral, "0-20", true, allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
    	    DATA_FILE_LOAD_DAU_SIM(Mostcentral, "0-20", true, allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
    	    DATA_FILE_LOAD_DAU_SIMETA(Mostcentral, "0-20", true, allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
    	    DATA_FILE_LOAD_DAU_SIM1GAMMA(Mostcentral, "0-20", true, allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
    	    DATA_FILE_LOAD_DAU_SIMNBAR(Mostcentral, "0-20", true, allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
	    pcPool poolsMostC[] = {poolMBMostcentral, poolHT1Mostcentral, poolHT2Mostcentral
		, poolEffMBMostC, poolEffHT1MostC, poolEffHT2MostC
		, poolEffEtaMBMostC, poolEffEtaHT1MostC, poolEffEtaHT2MostC
		, poolEffEtabgMBMostC, poolEffEtabgHT1MostC, poolEffEtabgHT2MostC
		, pool1gammaMBMostC, pool1gammaHT1MostC, pool1gammaHT2MostC
		, pool1gammaSimMBMostC, pool1gammaSimHT1MostC, pool1gammaSimHT2MostC
		, pool1gammaSimEtaMBMostC, pool1gammaSimEtaHT1MostC, pool1gammaSimEtaHT2MostC
		, poolnbarMBMostC, poolnbarHT1MostC, poolnbarHT2MostC
		, poolnbarSimMBMostC, poolnbarSimHT1MostC, poolnbarSimHT2MostC
		, poolnbarSimEtaMBMostC, poolnbarSimEtaHT1MostC, poolnbarSimEtaHT2MostC
		, poolPtshiftMBMostC, poolPtshiftHT1MostC, poolPtshiftHT2MostC};
	    showResults(&poolsMostC[0], allSettings.settingsDAuMostCentral, allResults.resultsDAuMostCentral);
    	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU_DATA(Mostcentral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralities) DATA_FILE_DELETE_DAU_SIM(Mostcentral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralitiesEta) DATA_FILE_DELETE_DAU_SIMETA(Mostcentral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIM1GAMMA(Mostcentral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIMNBAR(Mostcentral);

	    if (allSettings.showMidCentral) {
    		DATA_FILE_LOAD_DAU_DATA(Midcentral, "20-40", true, allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
    		DATA_FILE_LOAD_DAU_SIM(Midcentral, "20-40", true, allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
    		DATA_FILE_LOAD_DAU_SIMETA(Midcentral, "20-40", true, allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
    		DATA_FILE_LOAD_DAU_SIM1GAMMA(Midcentral, "20-40", true, allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
    		DATA_FILE_LOAD_DAU_SIMNBAR(Midcentral, "20-40", true, allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
		pcPool poolsMidC[] = {poolMBMidcentral, poolHT1Midcentral, poolHT2Midcentral
		    , poolEffMBMidC, poolEffHT1MidC, poolEffHT2MidC
		    , poolEffEtaMBMidC, poolEffEtaHT1MidC, poolEffEtaHT2MidC
		    , poolEffEtabgMBMidC, poolEffEtabgHT1MidC, poolEffEtabgHT2MidC
		    , pool1gammaMBMidC, pool1gammaHT1MidC, pool1gammaHT2MidC
		    , pool1gammaSimMBMidC, pool1gammaSimHT1MidC, pool1gammaSimHT2MidC
		    , pool1gammaSimEtaMBMidC, pool1gammaSimEtaHT1MidC, pool1gammaSimEtaHT2MidC
		    , poolnbarMBMidC, poolnbarHT1MidC, poolnbarHT2MidC
		    , poolnbarSimMBMidC, poolnbarSimHT1MidC, poolnbarSimHT2MidC
		    , poolnbarSimEtaMBMidC, poolnbarSimEtaHT1MidC, poolnbarSimEtaHT2MidC
		    , poolPtshiftMBMidC, poolPtshiftHT1MidC, poolPtshiftHT2MidC};
		showResults(&poolsMidC[0], allSettings.settingsDAuMidCentral, allResults.resultsDAuMidCentral);
        	if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU_DATA(Midcentral);
        	if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralities) DATA_FILE_DELETE_DAU_SIM(Midcentral);
        	if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralitiesEta) DATA_FILE_DELETE_DAU_SIMETA(Midcentral);
        	if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIM1GAMMA(Midcentral);
        	if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIMNBAR(Midcentral);
	    }

	    DATA_FILE_LOAD_DAU_DATA(Mostperipheral, "40-100", true, allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
	    DATA_FILE_LOAD_DAU_SIM(Mostperipheral, "40-100", true, allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
	    DATA_FILE_LOAD_DAU_SIMETA(Mostperipheral, "40-100", true, allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
	    DATA_FILE_LOAD_DAU_SIM1GAMMA(Mostperipheral, "40-100", true, allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
	    DATA_FILE_LOAD_DAU_SIMNBAR(Mostperipheral, "40-100", true, allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
	    pcPool poolsMostP[] = {poolMBMostperipheral, poolHT1Mostperipheral, poolHT2Mostperipheral
		, poolEffMBMostP, poolEffHT1MostP, poolEffHT2MostP
		, poolEffEtaMBMostP, poolEffEtaHT1MostP, poolEffEtaHT2MostP
		, poolEffEtabgMBMostP, poolEffEtabgHT1MostP, poolEffEtabgHT2MostP
		, pool1gammaMBMostP, pool1gammaHT1MostP, pool1gammaHT2MostP
		, pool1gammaSimMBMostP, pool1gammaSimHT1MostP, pool1gammaSimHT2MostP
		, pool1gammaSimEtaMBMostP, pool1gammaSimEtaHT1MostP, pool1gammaSimEtaHT2MostP
		, poolnbarMBMostP, poolnbarHT1MostP, poolnbarHT2MostP
		, poolnbarSimMBMostP, poolnbarSimHT1MostP, poolnbarSimHT2MostP
		, poolnbarSimEtaMBMostP, poolnbarSimEtaHT1MostP, poolnbarSimEtaHT2MostP
		, poolPtshiftMBMostP, poolPtshiftHT1MostP, poolPtshiftHT2MostP};
	    showResults(&poolsMostP[0], allSettings.settingsDAuMostPeripheral, allResults.resultsDAuMostPeripheral);
    	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_DAU_DATA(Mostperipheral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralities) DATA_FILE_DELETE_DAU_SIM(Mostperipheral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAuEffCentralitiesEta) DATA_FILE_DELETE_DAU_SIMETA(Mostperipheral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIM1GAMMA(Mostperipheral);
    	    if (!allSettings.browseAnalysis && allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIMNBAR(Mostperipheral);

	    showResultsDAuCentrality(allSettings.dAuCentralityName, allSettings.dAuCentralityTitle
                , allSettings.showRcp, allSettings.showRcpEta, allSettings.dAuCentralityShow, allSettings.dAuCentralityPrint, allSettings.showDAuCentralityTriggersSeparately
		, allSettings.saveRcpDataArrays, allSettings.dataArraysRcpFilenameFormat.Data()
                , allSettings.settingsDAuNoCentral
                , allSettings.settingsDAuAllCentral
                , allSettings.settingsDAuMostCentral
                , allSettings.settingsDAuMidCentral
                , allSettings.settingsDAuMostPeripheral
                , allResults.resultsDAuNoCentral
                , allResults.resultsDAuMostCentral
                , allResults.resultsDAuMidCentral
                , allResults.resultsDAuMostPeripheral
		, allSettings.dataPointsR
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
	if (!allSettings.browseAnalysis && !allSettings.useDAuEffCentralities) DATA_FILE_DELETE_DAU_SIM(Nocentral);
	if (!allSettings.browseAnalysis && !allSettings.useDAuEffCentralitiesEta) DATA_FILE_DELETE_DAU_SIMETA(Nocentral);
	if (!allSettings.browseAnalysis && !allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIM1GAMMA(Nocentral);
	if (!allSettings.browseAnalysis && !allSettings.useDAu1gammaCentralities) DATA_FILE_DELETE_DAU_SIMNBAR(Nocentral);

	if (allSettings.showSpectrumPP || allSettings.showSpectrumPPEta || allSettings.showRDA || allSettings.showRDAEta || allSettings.showRDAGamma) {
	    DATA_FILE_LOAD_PP(true, allSettings.settingsPP, allResults.resultsPP);
	    pcPool poolsPP[] = {poolPPMB, poolPPHT1, poolPPHT2
		, poolPPEffMB, poolPPEffHT1, poolPPEffHT2
		, poolPPEffEtaMB, poolPPEffEtaHT1, poolPPEffEtaHT2
		, poolPPEffEtabgMB, poolPPEffEtabgHT1, poolPPEffEtabgHT2
		, poolPP1gammaMB, poolPP1gammaHT1, poolPP1gammaHT2
		, poolPP1gammaSimMB, poolPP1gammaSimHT1, poolPP1gammaSimHT2
		, poolPP1gammaSimEtaMB, poolPP1gammaSimEtaHT1, poolPP1gammaSimEtaHT2
		, poolPPnbarMB, poolPPnbarHT1, poolPPnbarHT2
		, poolPPnbarSimMB, poolPPnbarSimHT1, poolPPnbarSimHT2
		, poolPPnbarSimEtaMB, poolPPnbarSimEtaHT1, poolPPnbarSimEtaHT2
		, poolPPPtshiftMB, poolPPPtshiftHT1, poolPPPtshiftHT2};
	    showResults(&poolsPP[0], allSettings.settingsPP, allResults.resultsPP);
	    if (!allSettings.browseAnalysis) DATA_FILE_DELETE_PP;
	}

	if (allSettings.showRDA || allSettings.showRDAEta || allSettings.showRDAGamma) {
	    showResultsRDA(allSettings.RDAName, allSettings.RDATitle
	        , allSettings.showRDA, allSettings.showRDAEta, allSettings.showRDAGamma, allSettings.RDAShow, allSettings.showRDATriggersSeparately
		, allSettings.saveRdADataArrays, allSettings.dataArraysRdAFilenameFormat.Data()
	        , allSettings.settingsDAuNoCentral
	        , allSettings.settingsPP
	        , allResults.resultsDAuNoCentral
	        , allResults.resultsPP
		, allSettings.dataPointsR
	        , allResults.RDAMB
	        , allResults.RDAHT1
	        , allResults.RDAHT2
	        , allResults.RDA
	        , allResults.RDAMBEta
	        , allResults.RDAHT1Eta
	        , allResults.RDAHT2Eta
	        , allResults.RDAEta
	        , allResults.RDAMBGamma
	        , allResults.RDAHT1Gamma
	        , allResults.RDAHT2Gamma
	        , allResults.RDAGamma
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

	TDatime stopPreTime;
	cout << "Finished calculations, drawing... " << stopPreTime.AsSQLString() << endl;
	gSystem->ProcessEvents();
	TDatime stopTime;
	cout << "Finished: " << stopTime.AsSQLString() << endl;
	timer.Stop();
	timer.Print();
	cout << "================== FINISHED ==================" << endl;
}
