#include "StPi0ResultsUtil.h"

#include "TDataProcessorMemberInspector.h"

#include <fstream>
#include <list>
#include <algorithm>
#include <map>
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
#include <TCanvas.h>
#include <TCollection.h>
#include <TLegend.h>
#include <TF1.h>
#include <TMath.h>
#include <TVirtualFitter.h>

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

const Char_t *fitFuncZero = "gaus(0) + gaus(3)";
const Char_t *fitFuncConst = "gaus(0) + gaus(3) + [6]";
const Char_t *fitFuncLinear = "gaus(0) + gaus(3) + [6] + ([7]*x)";
const Char_t *fitFuncZeroLinear = "gaus(0) + gaus(3) + ([6]*0) + ([7]*x)";

const Char_t *fitFuncTrackDist = "0*gaus(0) + 0*gaus(3) + ([6]*x*exp(-x*[7]))";
const Char_t *fitFuncTrackDistPeak = "gaus(0) + 0*gaus(3) + ([6]*x*exp(-x*[7]))";

TF1 funcOne("funcOne", "1");
TF1 funcZero("funcZero","0");

#define PQCDFUNC(usePPpQCD)   ((Double_t (*)(Double_t*, Double_t*))((usePPpQCD) ? (&getNLOpQCDPP) : (&getNLOpQCD)))
#define PQCDWEIGHT(usePPpQCD) ((usePPpQCD) ? (pQCDPPweight) : (pQCDweight))

#define GET_DATA(var, class, pool, name) \
    class *var = pool ? dynamic_cast<class*>(pool->getProcessor(name)) : 0;

ClassImp(TDrawOptions);
ClassImp(TDataPoints);
ClassImp(TAnalysisSettingsTrigger);
ClassImp(TAnalysisSettings);
ClassImp(TAnalysisRawResultsTrigger);
ClassImp(TAnalysisResults);

Double_t TDataPoints::operator() (Double_t *x, Double_t * /*p*/) const {
    Float_t val = funcDenom ? (funcDenom->Eval(x[0])) : 1;
    return (func && (val != 0)) ? (func->Eval(x[0]) / val) : 0;
}
TDataPoints &TDataPoints::operator=(const TDataPoints &p) {
        TDrawOptions::operator=(p);
        x = p.x;
        y = p.y;
        ex = p.ex;
        ey = p.ey;
        n = p.n;
        func = p.func;
        funcDenom = p.funcDenom;// ? (new TF1(*(p.funcDenom))) : 0;
        return *this;
}
TDataPoints &TDataPoints::operator/=(const TF1 &f) {
        funcDenom = &f;// ? (new TF1(*f)) : 0;
        if (funcDenom && x && n) {
            Float_t *newY = new Float_t[n];
            Float_t *newYe = new Float_t[n];
            for (UInt_t i = 0;i < n;i++) {
                Float_t val = funcDenom->Eval(x[i]);
                if (y && newY) newY[i] = ((val != 0) ? (y[i] / val) : 0);
                if (ey && newYe) newYe[i] = ((val != 0) ? (ey[i] / val) : 0);
            }
            y = newY;
            ey = newYe;
        }
	if (funcDenom) {
	    TString *newName = new TString(GetName());
	    *newName += "_div_";
	    *newName += funcDenom->GetName();
	    SetName(newName->Data()); 
	}
        return *this;
}
bin_stat_list_type TDataPoints::getBins(Float_t low, Float_t high) const {
    bin_stat_list_type bins;
    if (!this->isFunction() && (this->n > 0) && this->x) {
	for (UInt_t i = 0;i < this->n;i++) {
	    Float_t center = this->x[i];
	    if ((low <= center) && (center < high)) {
		TBinParameters binPar;
		binPar.variable = pT;
		binPar.min = center - (this->ex ? this->ex[i] : 0);
		binPar.max = center + (this->ex ? this->ex[i] : 0);
		TBinStatistics binStat;
		binStat.setParameters(binPar);
		binStat.setValue(this->y ? this->y[i] : 0);
		binStat.setError(this->ey ? this->ey[i] : 0);
		bins.push_back(binStat);
	    }
	}
    }
    return bins;
}
TDataPoints &TDataPoints::crop(Float_t lowPt, Float_t highPt) {
    if (!this->isFunction() && (this->n > 0) && this->x) {
	UInt_t j = 0;
	Float_t *x_new = new Float_t[this->n];
	Float_t *y_new = new Float_t[this->n];
	Float_t *ex_new = new Float_t[this->n];
	Float_t *ey_new = new Float_t[this->n];
	for (UInt_t i = 0;i < this->n;i++) {
	    Float_t center = this->x[i];
	    if ((lowPt <= center) && (center < highPt))	{
		x_new[j] = this->x ? this->x[i] : 0;
		y_new[j] = this->y ? this->y[i] : 0;
		ex_new[j] = this->ex ? this->ex[i] : 0;
		ey_new[j] = this->ey ? this->ey[i] : 0;
		j++;
	    }
	}
	this->n = j;
	this->x = x_new;
	this->y = y_new;
	this->ex = ex_new;
	this->ey = ey_new;
    }
    return *this;
}

TDataPoints::TDataPoints(const TDrawOptions &draw, const bin_stat_list_type &bins, Option_t *option)
       : TDrawOptions(draw), x(0), y(0), ex(0), ey(0), n(0), func(0), funcDenom(0) 
{
    Float_t *xx = new Float_t[bins.size()];
    Float_t *yy = new Float_t[bins.size()];
    Float_t *exx = new Float_t[bins.size()];
    Float_t *eyy = new Float_t[bins.size()];
    Int_t i = 0;
    for (bin_stat_list_type::const_iterator iter = bins.begin();iter != bins.end();++iter) {
	const TBinStatistics &bin = *iter;
	const TBinParameters &par = bin.getParameters();
	Float_t bx = par.getCenter();
	Float_t bex = (par.max - par.min) / 2.0;
	Float_t by = bin.getValue();
	Float_t bey = bin.getError();
	xx[i] = bx;
	exx[i] = bex;
	TString opt(option);
	if (opt.Contains("LOWEDGE")) {
	    yy[i] = by - bey;
	    eyy[i] = 0;
	} else if (opt.Contains("HIGHEDGE")) {
	    yy[i] = by + bey;
	    eyy[i] = 0;
	} else {
	    yy[i] = by;
	    eyy[i] = bey;
	}
	i++;
    }
    this->x = xx;
    this->y = yy;
    this->ex = exx;
    this->ey = eyy;
    this->n = i;
}

Bool_t operator==(const TDataPoints &/*p1*/, const TDataPoints &/*p2*/) {return false;}
Bool_t operator<(const TDataPoints &/*p1*/, const TDataPoints &/*p2*/) {return true;}
data_points_list &operator/=(data_points_list &l, const TF1 &f) {
    for (data_points_list::iterator it = l.begin();it != l.end();++it) {
        (*it) /= f;
    }
    return l;
}

void calculateDerivatives(const bin_stat_list_type &values, Float_t dx_relative, bin_stat_list_type &dy_relative, Option_t *fitOpt, Option_t *fitOptWeight) {
            Float_t minPt = -1;
            Float_t maxPt = -1;
            for (bin_stat_list_type::const_iterator it = values.begin();it != values.end();++it) {
                Float_t binCenter = (*it).getParameters().getCenter();
                if ((binCenter < minPt) || (minPt < 0)) minPt = binCenter;
                if ((binCenter > maxPt) || (maxPt < 0)) maxPt = binCenter;
            }
//cout << "DER: minPt = " << minPt << ", maxPt = " << maxPt << endl;
            TWeightCalculator w;
            w.Fit(values, fitOpt, fitOptWeight);
            TF1 *func = w.createFunc(true);
            if (func) for (bin_stat_list_type::const_iterator it = values.begin();it != values.end();++it) {
                const TBinStatistics &binStat = *it;
                const TBinParameters &binPar = binStat.getParameters();
                Float_t center = binPar.getCenter();
                Bool_t minPtBin = (center == minPt);
                Bool_t maxPtBin = (center == maxPt);
//if (minPtBin) cout << "DER: minPtBin center = " << center << endl;
//if (maxPtBin) cout << "DER: maxPtBin center = " << center << endl;
                Float_t rangeLeft = center * (1.0 - (minPtBin ? 0.0 : dx_relative));
                Float_t rangeRight = center * (1.0 + (maxPtBin ? 0.0 : dx_relative));
                Float_t val = func->Eval(center);
                Float_t valLeft = func->Eval(rangeLeft);
                Float_t valRight = func->Eval(rangeRight);
                Float_t systErr = ((val != 0) && (rangeLeft != rangeRight) && (dx_relative != 0)) ? (((TMath::Abs(valLeft - val) + TMath::Abs(val - valRight)) / (((rangeRight - rangeLeft)/center)/dx_relative)) / val) : 0;
                TBinStatistics binStatNew;
                binStatNew.setParameters(binPar);
                binStatNew.setValue(systErr);
                binStatNew.setError(0);
                dy_relative.push_back(binStatNew);
            }
            delete func; func = 0;
}

TH1F *histBeingFit = 0;
TH1F *histToFit1 = 0;
TH1F *histToFit2 = 0;
TH1F *histToFit3 = 0;
TH1F *histToFit4 = 0;
TH1F *histToFit5 = 0;
TH1F *histToFit6 = 0;
Float_t fitRangeLeft = 0;
Float_t fitRangeRight = 0;
Float_t fitRangeLeft2 = 0;
Float_t fitRangeRight2 = 0;
Double_t getHistWeightedSumErr(Double_t *x, Double_t *p, Double_t *err) {
	Double_t result = 0.0;
	Double_t error = 0.0;
	if (x && p) {
	    if (((x[0] >= fitRangeLeft) && (x[0] <= fitRangeRight)) || ((x[0] >= fitRangeLeft2) && (x[0] <= fitRangeRight2))) {
		Int_t ibin = 0;
		Float_t e = 0;
		if (histToFit1) {
		    ibin = histToFit1->GetXaxis()->FindBin(x[0]);
		    e = (p[0] * p[1] * histToFit1->GetBinError(ibin));
		    error += (e * e);
		    result += (p[0] * p[1] * histToFit1->GetBinContent(ibin));
		}
		if (histToFit2) {
		    ibin = histToFit2->GetXaxis()->FindBin(x[0]);
		    e = (p[0] * (1.0 - p[1]) * histToFit2->GetBinError(ibin));
		    error += (e * e);
		    result += (p[0] * (1.0 - p[1]) * histToFit2->GetBinContent(ibin));
		}
		if (histToFit3) {
		    ibin = histToFit3->GetXaxis()->FindBin(x[0]);
		    e = (p[2] * histToFit3->GetBinError(ibin));
		    error += (e * e);
		    result += (p[2] * histToFit3->GetBinContent(ibin));
		}
		if (histToFit4) {
		    ibin = histToFit4->GetXaxis()->FindBin(x[0]);
		    e = (p[3] * histToFit4->GetBinError(ibin));
		    error += (e * e);
		    result += (p[3] * histToFit4->GetBinContent(ibin));
		}
		if (histToFit5) {
		    ibin = histToFit5->GetXaxis()->FindBin(x[0]);
		    e = (p[4] * histToFit5->GetBinError(ibin));
		    error += (e * e);
		    result += (p[4] * histToFit5->GetBinContent(ibin));
		}
		if (histToFit6) {
		    ibin = histToFit6->GetXaxis()->FindBin(x[0]);
		    e = (p[4] * p[5] * histToFit6->GetBinError(ibin));
		    error += (e * e);
		    result += (p[4] * p[5] * histToFit6->GetBinContent(ibin));
		}
	    }
	}
	error = TMath::Sqrt(error);
	if (err) *err = error;
//cerr << "-"; cerr.flush();
        return result;
}
Double_t getHistWeightedSum(Double_t *x, Double_t *p) {
    Double_t err = 0;
    return getHistWeightedSumErr(x, p, &err);
}

void myHistHistFitChi2(Int_t& /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag*/) {
    fval = 0;
    //nPar = 6;
    if (!histBeingFit || !p) return;
//cerr << "<" << nPar << " " << iflag; cerr.flush();
    Int_t binMax = histBeingFit->GetXaxis()->GetNbins();
    Float_t chi2 = 0;
    for (Int_t ibin = 1;ibin <= binMax;ibin++) {
	Double_t binCenter = histBeingFit->GetXaxis()->GetBinCenter(ibin);
	if (((binCenter >= fitRangeLeft) && (binCenter <= fitRangeRight)) || ((binCenter >= fitRangeLeft2) && (binCenter <= fitRangeRight2))) {
	    Float_t valBeingFit = histBeingFit->GetBinContent(ibin);
    	    Float_t errBeingFit = histBeingFit->GetBinError(ibin);
    	    Double_t errFunc = 0;
	    Float_t valFunc = getHistWeightedSumErr(&binCenter, p, &errFunc);
	    Float_t valDiff = valBeingFit - valFunc;
	    Float_t errDiff2 = (errBeingFit * errBeingFit) + (errFunc * errFunc);
	    if (errDiff2 != 0) chi2 += (valDiff * valDiff / errDiff2);
	}
    }
    fval = chi2;
//cerr << " " << fval << "> "; cerr.flush();
}

const Float_t noLimitMagic = 1.01010101;
Float_t getParameterValue(const TFitSettings &param, const TBinParameters &binParameters, Float_t *limitLow = 0, Float_t *limitHigh = 0) {
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

void setBinsToFit(bin_stat_list_type &bins, const TFitSettings &param) {
    for (bin_stat_list_type::iterator iter = bins.begin();iter != bins.end();++iter) {
	TBinStatistics &binStat = *iter;
	Float_t val = getParameterValue(param, binStat.getParameters());
	binStat.setValue(val);
    }
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
		if (!param.func) param.func = new TF1(name, (const Char_t*)param.fitFuncStr);
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

Int_t getBestRebin(const TH1F *hist, Float_t numBinsWeight, Float_t meanErrorWeight, Float_t rangeLeft, Float_t rangeRight, Int_t minBins, Float_t minBinWidth) {
	Int_t nrebin = 1;
//	Int_t nrebinWidth = 1;
	if (((numBinsWeight != 0) || (meanErrorWeight != 0)) && hist) {
	    Float_t sumRebinMin = -1.0;
//	    Float_t widthRebinMin = -1.0;
	    Int_t nrebinMax = (hist->GetXaxis()->FindBin(rangeRight) - hist->GetXaxis()->FindBin(rangeLeft)) / minBins;
	    Int_t nrebinMin = Int_t(minBinWidth / hist->GetXaxis()->GetBinWidth(1)) + 1;
	    //cout << "nrebinMax = " << nrebinMax << endl;
	    for (Int_t ntemp = nrebinMin;ntemp < nrebinMax;ntemp++) {
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
			Float_t errRel = TMath::Abs(err / val);
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

void setErrorOnZeroBins(TH1F *hist, Float_t sigma) {
    if (hist && hist->GetXaxis()) {
	Float_t valuePerEntry = hist->GetEntries() ? (hist->Integral() / hist->GetEntries()) : 0;
	Float_t nsigma = 2.0;
	for (Int_t mbin = 1;mbin <= hist->GetXaxis()->GetNbins();mbin++) {
	    if (hist->GetBinContent(mbin) == 0) {
		Float_t valueFraction = 1.0;
		if (sigma > 0) {
		    Float_t binCenter = hist->GetXaxis()->GetBinCenter(mbin);
		    Float_t binWidth = hist->GetXaxis()->GetBinWidth(mbin);
		    //sigma = (binCenter * 0.07282) + 0.02;
		    Int_t nsigmaBins = Int_t((sigma * nsigma / binWidth) + 0.5);
		    Float_t sum = 0;
		    Float_t sumInt = 0;
		    for (Int_t sbin = mbin - nsigmaBins;sbin <= mbin + nsigmaBins;sbin++) {
			if (sbin != mbin) {
			    Float_t x = TMath::Abs(hist->GetXaxis()->GetBinCenter(sbin) - binCenter) - binWidth;
			    Float_t val = hist->GetBinContent(sbin);
			    sum += val * TMath::Exp(- (x*x) / (2.0*sigma*sigma));
			    sumInt += valuePerEntry;
			}
		    }
		    if (sumInt > 0) sum /= sumInt;
		    Float_t sumMin = 0.418 / (nsigmaBins + 1 + nsigmaBins);
		    if (sum < sumMin) sum = sumMin;
		    valueFraction = sum;
		}
		hist->SetBinError(mbin, valuePerEntry * valueFraction);
	    }
	}
    }
}

void smearOutliers(TH1F *hist, Float_t sigma, Float_t nsigmaCut, Bool_t zeroOutliers) {
    if (hist && hist->GetXaxis() && (sigma > 0) && (nsigmaCut > 0)) {
	//cout << "Integral " << hist->GetName() << ": " << hist->Integral();
	Int_t nOut = 0;
	TH1F smeared(*hist);
	for (Int_t mbin = 1;mbin <= smeared.GetXaxis()->GetNbins();mbin++) {
	    smeared.SetBinContent(mbin, 0);
	    smeared.SetBinError(mbin, 0);
	}
	for (Int_t mbin = 1;mbin <= hist->GetXaxis()->GetNbins();mbin++) {
	    Float_t val = hist->GetBinContent(mbin);
	    Float_t err = hist->GetBinError(mbin);
	    Float_t valPre = hist->GetBinContent(mbin - 1);
	    Float_t errPre = hist->GetBinError(mbin - 1);
	    Float_t valPost = hist->GetBinContent(mbin + 1);
	    Float_t errPost = hist->GetBinError(mbin + 1);
	    Float_t nent = (val * val) / (err * err);
	    Float_t nentPre = (valPre * valPre) / (errPre * errPre);
	    Float_t nentPost = (valPost * valPost) / (errPost * errPost);
	    // defining the outlier
	    Bool_t isOut =     (val > 0) && (err > 0) && (((nent > 1.001) && (nent < 1.999)) || ((nent > 2.001) && (nent < 2.999)));
	    Bool_t isOutPre =  (valPre > 0) && (errPre > 0) && (((nentPre > 1.001) && (nentPre < 1.999)) || ((nentPre > 2.001) && (nentPre < 2.999)));
	    Bool_t isOutPost = (valPost > 0) && (errPost > 0) && (((nentPost > 1.001) && (nentPost < 1.999)) || ((nentPost > 2.001) && (nentPost < 2.999)));
	    if (isOut) {
		Float_t valNear = 0;
		Int_t valNearNum = 0;
		if ((valPre > 0) && !isOutPre) {valNear += valPre; valNearNum++;}
		if ((valPost > 0) && !isOutPost) {valNear += valPost; valNearNum++;}
		if (valNearNum) valNear /= valNearNum;
		Float_t valGaus = val - valNear;
		if (valGaus > 0) {
		    if (zeroOutliers) {
			hist->SetBinContent(mbin, valNear);
			hist->SetBinError(mbin, hist->GetBinError(mbin) * valNear / valGaus);
		    } else {
			//cout << mbin << " ";
			//cout << "out = " << isOutPre << "," << isOut << "," << isOutPost << " ";
			//cout << ", val=" << val << ":" << err << ", valNear=" << valNear << ", valGaus=" << valGaus << " ";
			hist->SetBinContent(mbin, valNear);
			hist->SetBinError(mbin, 0);
			Float_t binCenter = hist->GetXaxis()->GetBinCenter(mbin);
		        Float_t binWidth = hist->GetXaxis()->GetBinWidth(mbin);
			Int_t nbins = (Int_t)(((nsigmaCut * sigma) / binWidth + 0.5));
			//cout << ", center=" << binCenter << ", width=" << binWidth << ", nbins=" << nbins << " ";
			TH1F smearedGaus(smeared);
			for (Int_t mmbin = 1;mmbin < smearedGaus.GetXaxis()->GetNbins();mmbin++) {
			    smearedGaus.SetBinContent(mmbin, 0);
			    smearedGaus.SetBinError(mmbin, 0);
			}
			//cout << ": ";
			Float_t sumErr2 = 0;
			for (Int_t sbin = mbin - nbins;sbin <= mbin + nbins;sbin++) {
			    Float_t x = hist->GetXaxis()->GetBinCenter(sbin) - binCenter;
			    Float_t curVal = valGaus * binWidth * (1.0 / (TMath::Sqrt(2.0 * 3.1415) * sigma)) * TMath::Exp(-(x * x) / (2.0 * sigma * sigma));
			    Float_t curErr = curVal;
			    //cout << curVal << " ";
			    smearedGaus.SetBinContent(sbin, curVal);
			    smearedGaus.SetBinError(sbin, curErr);
			    sumErr2 += (curErr * curErr);
			}
			Float_t errMult = (sumErr2 != 0) ? TMath::Sqrt(err * err / sumErr2) : 0.0;
			//cout << ", errMult=" << errMult;
			if (errMult > 0) {
			    for (Int_t sbin = mbin - nbins;sbin <= mbin + nbins;sbin++) {
				smearedGaus.SetBinError(sbin, smearedGaus.GetBinError(sbin) * errMult);
			    }
			}
			Float_t totGausVal = 0;
			Float_t totGausErr = 0;
		        for (Int_t sbin = mbin - nbins;sbin <= mbin + nbins;sbin++) {
			    totGausVal += smearedGaus.GetBinContent(sbin);
			    totGausErr += smearedGaus.GetBinError(sbin) * smearedGaus.GetBinError(sbin);
			    smeared.SetBinContent(sbin, smeared.GetBinContent(sbin) + smearedGaus.GetBinContent(sbin));
			    smeared.SetBinError(sbin, TMath::Sqrt((smeared.GetBinError(sbin) * smeared.GetBinError(sbin)) + (smearedGaus.GetBinError(sbin) * smearedGaus.GetBinError(sbin))));
			}
		        totGausErr = TMath::Sqrt(totGausErr);
		        //cout << "; totGaus=" << totGausVal << ":" << totGausErr;
		        //cout << endl;
		        nOut++;
		    }
		}
	    }
	}
	if (nOut) {
	    for (Int_t mbin = 1;mbin <= hist->GetXaxis()->GetNbins();mbin++) {
		hist->SetBinContent(mbin, smeared.GetBinContent(mbin) + hist->GetBinContent(mbin));
		hist->SetBinError(mbin, TMath::Sqrt((smeared.GetBinError(mbin) * smeared.GetBinError(mbin)) + (hist->GetBinError(mbin) * hist->GetBinError(mbin))));
	    }
	}
	//cout << " " << hist->Integral() << ", " << nOut << " outliers" << endl;
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

void smoothHistPol2(TH1F *hist, Int_t nbins, Option_t *opt) {
    if (!hist || (nbins <= 0)) return;
    TH1F histSmoothed(*hist);
    TF1 func("func", "pol2");
    Int_t maxBin = hist->GetXaxis()->GetNbins();
    for (Int_t i = 1;i <= maxBin;i++) {
	Int_t binLow = i - nbins;
	Int_t binHigh = i + nbins;
	if (binLow < 1) binLow = 1;
	if (binHigh > maxBin) binHigh = maxBin;
	Float_t x = hist->GetXaxis()->GetBinCenter(i);
	Float_t rangeLow = hist->GetXaxis()->GetBinLowEdge(binLow);
	Float_t rangeHigh = hist->GetXaxis()->GetBinUpEdge(binHigh);
	func.SetRange(rangeLow, rangeHigh);
	hist->GetXaxis()->SetRangeUser(rangeLow, rangeHigh);
	hist->Fit(&func, opt, "", rangeLow, rangeHigh);
	//Float_t oldY = hist->GetBinContent(i);
	Float_t oldYerr = hist->GetBinError(i);
	Float_t newY = func.Eval(x);
	Float_t newYerr = oldYerr;
	histSmoothed.SetBinContent(i, newY);
	histSmoothed.SetBinError(i, newYerr);
    }
    for (Int_t i = 1;i <= maxBin;i++) {
	hist->SetBinContent(i, histSmoothed.GetBinContent(i));
	hist->SetBinError(i, histSmoothed.GetBinError(i));
    }
}

void smoothHistGaus(TH1F *hist, Float_t sigma, Float_t nsigma) {
    if (!hist || (sigma <= 0) || (nsigma <= 0)) return;
    TH1F histSmoothed(*hist);
    Int_t maxBin = hist->GetXaxis()->GetNbins();
    for (Int_t i = 1;i <= maxBin;i++) {
	Float_t binCenter = hist->GetXaxis()->GetBinCenter(i);
	Float_t binWidth = hist->GetXaxis()->GetBinWidth(i);
	Int_t nsigmaBins = Int_t((sigma * nsigma / binWidth) + 0.5);
	Float_t sum = 0;
	Float_t sumW = 0;
	for (Int_t sbin = i - nsigmaBins;sbin <= i + nsigmaBins;sbin++) {
	    Float_t x = TMath::Abs(hist->GetXaxis()->GetBinCenter(sbin) - binCenter);
	    Float_t val = hist->GetBinContent(sbin);
	    Float_t w = TMath::Exp(- (x*x) / (2.0*sigma*sigma));
	    sum  += val * w;
	    sumW += 1.0 * w;
	}
	if (sumW != 0) sum /= sumW;
	Float_t oldYerr = hist->GetBinError(i);
	histSmoothed.SetBinContent(i, sum);
	histSmoothed.SetBinError(i, oldYerr);
    }
    for (Int_t i = 1;i <= maxBin;i++) {
	hist->SetBinContent(i, histSmoothed.GetBinContent(i));
	hist->SetBinError(i, histSmoothed.GetBinError(i));
    }
}

TCanvas *findPeak(const TFindPeakDistributions &s, Int_t numBinsToShow, Int_t binToShow, TCanvas *canvas, const TAnalysisSettingsTrigger &analysisSettings, TAnalysisRawResultsTrigger &analysisResults) {
	if (!s.distribution) return canvas;

	TBinParameters binParameters = s.distribution->getBinParameters();

	Float_t fitLeft = getParameterValue(analysisSettings.fitLeft, binParameters);
        Float_t fitRight = getParameterValue(analysisSettings.fitRight, binParameters);
	Float_t fitLeftEta = getParameterValue(analysisSettings.fitLeftEta, binParameters);
        Float_t fitRightEta = getParameterValue(analysisSettings.fitRightEta, binParameters);

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

	Float_t lowmassbgNbar = getParameterValue(analysisSettings.lowmassbgNbar, binParameters);
//cout << "pT = " << binParameters.getCenter() << ": nbar = " << lowmassbgNbar << ": " << analysisSettings.name << endl;
	
	TH1F invSigBg(*(s.distribution)->getDistribution());
	TString newTitle;
	//newTitle += "  ";
	//newTitle += analysisSettings.title;
	//newTitle += ", ";
	newTitle += s.distribution->GetTitle();
	//newTitle += "  ";
	invSigBg.SetNameTitle(analysisSettings.name + "_invSigBg_" + s.distribution->GetName(), newTitle.Data());
	Int_t nrebin = getBestRebin(&invSigBg, analysisSettings.rebinNumBinsWeight, analysisSettings.rebinMeanErrorWeight, 0.0, 1.0, 10, 0.01);
	if (nrebin > 1) invSigBg.Rebin(nrebin);
	Float_t smearOutliersSigma = analysisSettings.smearOutliersSigma;
	Float_t zeroOutliers = (analysisSettings.smearOutliersSigma > 10);
	Float_t smearOutliersNSigma = 3.0;
	Float_t smoothGausNSigma = 3.0;
	Bool_t doSmearOutliers = (smearOutliersSigma > 0) || zeroOutliers;
//cout << invSigBg.GetName() << " invSigBg.Integral() = " << invSigBg.Integral();
	if (doSmearOutliers) smearOutliers(&invSigBg, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
//cout << " " << invSigBg.Integral() << endl;
	if (analysisSettings.setErrorOnZeroBins) setErrorOnZeroBins(&invSigBg, analysisSettings.setErrorOnZeroBinsSigma);

	TH1F invSig(invSigBg);
	invSig.SetName(analysisSettings.name + "_invSig_" + s.distribution->GetName());

	Float_t scaleBgRandom = 1.0;
	TH1F *invBgRandom = s.distributionBgRandom ? (new TH1F(*(s.distributionBgRandom)->getDistribution())) : 0;
	if (invBgRandom) {
		invBgRandom->SetName(analysisSettings.name + "_invBgRandom_" + invBgRandom->GetName());
		if (nrebin > 1) invBgRandom->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invBgRandom, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.smoothBgRandom > 0) {
		    smoothHistGaus(invBgRandom, TMath::Abs(analysisSettings.smoothBgRandom), smoothGausNSigma);
		    //invBgRandom->Smooth(analysisSettings.smoothBgRandom, -1, -1);
		} else if (analysisSettings.smoothBgRandom < 0) {
		    smoothHistPol2(invBgRandom, (Int_t)TMath::Abs(analysisSettings.smoothBgRandom), analysisSettings.smoothFitOption);
		}
		if (analysisSettings.setErrorOnZeroBinsBgRandom) setErrorOnZeroBins(invBgRandom, analysisSettings.setErrorOnZeroBinsSigma);
		if (analysisSettings.mixNormFixedRandom > 0) {
		    Float_t corr = analysisSettings.mixNormFixedRandom;
		    if (analysisSettings.mixNormFixedRandom > 0) scaleBgRandom *= corr;
		}
		if (analysisSettings.mixNormMassRangeRandom) {
		    Float_t intBg = invBgRandom->Integral(invBgRandom->GetXaxis()->FindBin(analysisSettings.mixNormLeftRandom), invBgRandom->GetXaxis()->FindBin(analysisSettings.mixNormRightRandom));
		    Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(analysisSettings.mixNormLeftRandom), invSig.GetXaxis()->FindBin(analysisSettings.mixNormRightRandom));
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormMassRangeRandom) scaleBgRandom *= corr;
		}
		if (analysisSettings.mixNormEntriesRandom) {
		    Float_t intBg = invBgRandom->Integral();
		    Float_t intSig = invSig.Integral();
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormEntriesRandom) scaleBgRandom *= corr;
		}
		invBgRandom->Scale(scaleBgRandom);
		if (analysisSettings.subtractBgRandom) {
		    //cout << invSig.GetName() << " subtracting bg random" << endl;
		    invSig.Add(invBgRandom, -1);
		}
	}

	Float_t scaleBg = 1.0;
	TH1F *invBg = s.distributionBg ? (new TH1F(*(s.distributionBg)->getDistribution())) : 0;
	if (invBg) {
		invBg->SetName(analysisSettings.name + "_invBg_" + invBg->GetName());
		if (nrebin > 1) invBg->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invBg, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.smoothBg > 0) {
		    smoothHistGaus(invBg, TMath::Abs(analysisSettings.smoothBg), smoothGausNSigma);
		    //invBg->Smooth(analysisSettings.smoothBg, -1, -1);
		} else if (analysisSettings.smoothBg < 0) {
		    smoothHistPol2(invBg, (Int_t)TMath::Abs(analysisSettings.smoothBg), analysisSettings.smoothFitOption);
		}
		if (analysisSettings.setErrorOnZeroBinsBg) setErrorOnZeroBins(invBg, analysisSettings.setErrorOnZeroBinsSigma);
		if (analysisSettings.mixNormFixed > 0) {
		    Float_t corr = analysisSettings.mixNormFixed;
		    if (analysisSettings.mixNormFixed > 0) scaleBg *= corr;
		}
		if (analysisSettings.mixNormMassRange) {
		    Float_t intBg = invBg->Integral(invBg->GetXaxis()->FindBin(analysisSettings.mixNormLeft), invBg->GetXaxis()->FindBin(analysisSettings.mixNormRight));
		    Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(analysisSettings.mixNormLeft), invSig.GetXaxis()->FindBin(analysisSettings.mixNormRight));
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormMassRange) scaleBg *= corr;
		}
		if (analysisSettings.mixNormEntries) {
		    Float_t intBg = invBg->Integral();
		    Float_t intSig = invSig.Integral();
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.mixNormEntries) scaleBg *= corr;
		}
		invBg->Scale(scaleBg);
		if (analysisSettings.subtractBg) {
		    //cout << invSig.GetName() << " subtracting bg" << endl;
		    invSig.Add(invBg, -1);
		}
	}

	Float_t scaleBgLowMass = 1.0;
	Float_t scaleBgLowMass2 = 1.0;
	TH1F *invBgLowMass = s.distributionBgLowMass ? (new TH1F(*(s.distributionBgLowMass)->getDistribution())) : 0;
	TH1F *invBgLowMass2 = s.distributionBgLowMass2 ? (new TH1F(*(s.distributionBgLowMass2)->getDistribution())) : 0;
	if (invBgLowMass2) {
		invBgLowMass2->SetName(analysisSettings.name + "_invBgLowMass2_" + invBgLowMass2->GetName());
		if (nrebin > 1) invBgLowMass2->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invBgLowMass2, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsBgLowMass) setErrorOnZeroBins(invBgLowMass2, analysisSettings.setErrorOnZeroBinsSigma);
		Float_t lowNormFixed_limitLow, lowNormFixed_limitHigh;
		Float_t lowNormFixed = getParameterValue(analysisSettings.lowNormFixed, binParameters, &lowNormFixed_limitLow, &lowNormFixed_limitHigh);
		if (lowNormFixed > 0) {
		    Float_t corr = lowNormFixed;
		    if (lowNormFixed > 0) scaleBgLowMass2 *= corr;
		}
		if (analysisSettings.lowNormMassRange) {
		    Float_t intBg = invBgLowMass2->Integral(invBgLowMass2->GetXaxis()->FindBin(analysisSettings.lowNormMassRangeLeft), invBgLowMass2->GetXaxis()->FindBin(analysisSettings.lowNormMassRangeRight));
		    Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(analysisSettings.lowNormMassRangeLeft), invSig.GetXaxis()->FindBin(analysisSettings.lowNormMassRangeRight));
		    Float_t corr = (intBg != 0) ? (intSig / intBg) : 1.0;
		    if (analysisSettings.lowNormMassRange) scaleBgLowMass2 *= corr;
		}
		if (analysisSettings.lowNormPoints) {
		    Float_t corr = analysisResults.lowNormPointsNbar;
		    if (analysisSettings.lowNormPoints) scaleBgLowMass2 *= corr;
		}
		if (analysisSettings.lowNormPointsPt) {
		    Float_t corr = analysisResults.lowNormPointsNbarPt;
		    if (analysisSettings.lowNormPointsPt) scaleBgLowMass2 *= corr;
		}
		scaleBgLowMass2 *= lowmassbgNbar;
		invBgLowMass2->Scale(scaleBgLowMass2);
	}
	if (invBgLowMass) {
		invBgLowMass->SetName(analysisSettings.name + "_invBgLowMass_" + invBgLowMass->GetName());
		if (nrebin > 1) invBgLowMass->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invBgLowMass, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsBgLowMass) setErrorOnZeroBins(invBgLowMass, analysisSettings.setErrorOnZeroBinsSigma);
		Float_t lowNormFixed_limitLow, lowNormFixed_limitHigh;
		Float_t lowNormFixed = getParameterValue(analysisSettings.lowNormFixed, binParameters, &lowNormFixed_limitLow, &lowNormFixed_limitHigh);
		if (lowNormFixed > 0) {
		    Float_t corr = lowNormFixed;
		    if (lowNormFixed > 0) scaleBgLowMass *= corr;
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
		scaleBgLowMass *= (1.0 - lowmassbgNbar);
		invBgLowMass->Scale(scaleBgLowMass);
		if (invBgLowMass2) invBgLowMass->Add(invBgLowMass2);
		if (analysisSettings.subtractBgLowMass) {
		    invSig.Add(invBgLowMass, -1);
		}
	}

	Float_t scalePeakShape = 1.0;
	TH1F *invPeakShape = s.distributionPeakShape ? (new TH1F(*(s.distributionPeakShape)->getDistribution())) : 0;
	if (invPeakShape) {
		invPeakShape->SetName(analysisSettings.name + "_invPeakShape_" + invPeakShape->GetName());
		if (nrebin > 1) invPeakShape->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShape, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShape) setErrorOnZeroBins(invPeakShape, analysisSettings.setErrorOnZeroBinsSigma);
		Float_t intBg = invPeakShape->Integral(invPeakShape->GetXaxis()->FindBin(peakLeft), invPeakShape->GetXaxis()->FindBin(peakRight));
		Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(peakLeft), invSig.GetXaxis()->FindBin(peakRight));
		scalePeakShape = (intBg != 0) ? (intSig / intBg) : 1.0;
		invPeakShape->Scale(scalePeakShape);
	}
	TH1F *invPeakShapeLowMass = s.distributionPeakShapeLowMass ? (new TH1F(*(s.distributionPeakShapeLowMass)->getDistribution())) : 0;
	TH1F *invPeakShapeLowMass2 = s.distributionPeakShapeLowMass2 ? (new TH1F(*(s.distributionPeakShapeLowMass2)->getDistribution())) : 0;
	if (invPeakShapeLowMass) {
		invPeakShapeLowMass->SetName(analysisSettings.name + "_invPeakShapeLowMass_" + invPeakShapeLowMass->GetName());
		if (nrebin > 1) invPeakShapeLowMass->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeLowMass, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShape) setErrorOnZeroBins(invPeakShapeLowMass, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeLowMass->Scale(scalePeakShape);
		Float_t scalePeakShapeLowMass = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeLowMass.begin();iter != analysisResults.scalePeakShapeLowMass.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeLowMass = (*iter).getValue();
		}
		invPeakShapeLowMass->Scale(scalePeakShapeLowMass);
	}
	if (invPeakShapeLowMass2) {
		invPeakShapeLowMass2->SetName(analysisSettings.name + "_invPeakShapeLowMass2_" + invPeakShapeLowMass2->GetName());
		if (nrebin > 1) invPeakShapeLowMass2->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeLowMass2, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShape) setErrorOnZeroBins(invPeakShapeLowMass2, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeLowMass2->Scale(scalePeakShape);
		Float_t scalePeakShapeLowMass2 = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeLowMass2.begin();iter != analysisResults.scalePeakShapeLowMass2.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeLowMass2 = (*iter).getValue();
		}
		invPeakShapeLowMass2->Scale(scalePeakShapeLowMass2);
	}
	if (invPeakShape && invPeakShapeLowMass && analysisSettings.subtractBgLowMassPeakShape) {
	    invPeakShape->Add(invPeakShapeLowMass, -1);
	}
	if (invPeakShape && invPeakShapeLowMass2 && analysisSettings.subtractBgLowMassPeakShape) {
	    invPeakShape->Add(invPeakShapeLowMass2, -1);
	}

	Float_t scalePeakShapeEta = 1.0;
	TH1F *invPeakShapeEta = s.distributionPeakShapeEta ? (new TH1F(*(s.distributionPeakShapeEta)->getDistribution())) : 0;
	if (invPeakShapeEta) {
		invPeakShapeEta->SetName(analysisSettings.name + "_invPeakShapeEta_" + invPeakShapeEta->GetName());
		if (nrebin > 1) invPeakShapeEta->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeEta, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEta, analysisSettings.setErrorOnZeroBinsSigma);
		Float_t intBg = invPeakShapeEta->Integral(invPeakShapeEta->GetXaxis()->FindBin(peakLeftEta), invPeakShapeEta->GetXaxis()->FindBin(peakRightEta));
		Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(peakLeftEta), invSig.GetXaxis()->FindBin(peakRightEta));
		scalePeakShapeEta = (intBg != 0) ? (intSig / intBg) : 1.0;
		invPeakShapeEta->Scale(scalePeakShapeEta);
	}
	TH1F *invPeakShapeEtaLowMass = s.distributionPeakShapeEtaLowMass ? (new TH1F(*(s.distributionPeakShapeEtaLowMass)->getDistribution())) : 0;
	TH1F *invPeakShapeEtaLowMass2 = s.distributionPeakShapeEtaLowMass2 ? (new TH1F(*(s.distributionPeakShapeEtaLowMass2)->getDistribution())) : 0;
	if (invPeakShapeEtaLowMass) {
		invPeakShapeEtaLowMass->SetName(analysisSettings.name + "_invPeakShapeEtaLowMass_" + invPeakShapeEtaLowMass->GetName());
		if (nrebin > 1) invPeakShapeEtaLowMass->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeEtaLowMass, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtaLowMass, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeEtaLowMass->Scale(scalePeakShapeEta);
		Float_t scalePeakShapeEtaLowMass = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeEtaLowMass.begin();iter != analysisResults.scalePeakShapeEtaLowMass.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeEtaLowMass = (*iter).getValue();
		}
		invPeakShapeEtaLowMass->Scale(scalePeakShapeEtaLowMass);
	}
	if (invPeakShapeEtaLowMass2) {
		invPeakShapeEtaLowMass2->SetName(analysisSettings.name + "_invPeakShapeEtaLowMass2_" + invPeakShapeEtaLowMass2->GetName());
		if (nrebin > 1) invPeakShapeEtaLowMass2->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeEtaLowMass2, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtaLowMass2, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeEtaLowMass2->Scale(scalePeakShapeEta);
		Float_t scalePeakShapeEtaLowMass2 = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeEtaLowMass2.begin();iter != analysisResults.scalePeakShapeEtaLowMass2.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeEtaLowMass2 = (*iter).getValue();
		}
		invPeakShapeEtaLowMass2->Scale(scalePeakShapeEtaLowMass2);
	}
	if (invPeakShapeEta && invPeakShapeEtaLowMass && analysisSettings.subtractBgLowMassPeakShapeEta) {
	    invPeakShapeEta->Add(invPeakShapeEtaLowMass, -1);
	}
	if (invPeakShapeEta && invPeakShapeEtaLowMass2 && analysisSettings.subtractBgLowMassPeakShapeEta) {
	    invPeakShapeEta->Add(invPeakShapeEtaLowMass2, -1);
	}

	TH1F *invPeakShapeEtabg = s.distributionPeakShapeEtabg ? (new TH1F(*(s.distributionPeakShapeEtabg)->getDistribution())) : 0;
	if (invPeakShapeEtabg) {
		invPeakShapeEtabg->SetName(analysisSettings.name + "_invPeakShapeEtabg_" + invPeakShapeEtabg->GetName());
		if (nrebin > 1) invPeakShapeEtabg->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeEtabg, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtabg, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeEtabg->Scale(scalePeakShapeEta);
		if (analysisResults.scalePeakShapeEtabgRel > 0) invPeakShapeEtabg->Scale(analysisResults.scalePeakShapeEtabgRel);
	}
	TH1F *invPeakShapeEtabgLowMass = s.distributionPeakShapeEtabgLowMass ? (new TH1F(*(s.distributionPeakShapeEtabgLowMass)->getDistribution())) : 0;
	TH1F *invPeakShapeEtabgLowMass2 = s.distributionPeakShapeEtabgLowMass2 ? (new TH1F(*(s.distributionPeakShapeEtabgLowMass2)->getDistribution())) : 0;
	if (invPeakShapeEtabgLowMass) {
		invPeakShapeEtabgLowMass->SetName(analysisSettings.name + "_invPeakShapeEtabgLowMass_" + invPeakShapeEtabgLowMass->GetName());
		if (nrebin > 1) invPeakShapeEtabgLowMass->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeEtabgLowMass, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtabgLowMass, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeEtabgLowMass->Scale(scalePeakShapeEta);
		if (analysisResults.scalePeakShapeEtabgRel > 0) invPeakShapeEtabgLowMass->Scale(analysisResults.scalePeakShapeEtabgRel);
		Float_t scalePeakShapeEtabgLowMass = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeEtabgLowMass.begin();iter != analysisResults.scalePeakShapeEtabgLowMass.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeEtabgLowMass = (*iter).getValue();
		}
		invPeakShapeEtabgLowMass->Scale(scalePeakShapeEtabgLowMass);
	}
	if (invPeakShapeEtabgLowMass2) {
		invPeakShapeEtabgLowMass2->SetName(analysisSettings.name + "_invPeakShapeEtabgLowMass2_" + invPeakShapeEtabgLowMass2->GetName());
		if (nrebin > 1) invPeakShapeEtabgLowMass2->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invPeakShapeEtabgLowMass2, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.setErrorOnZeroBinsPeakShapeEta) setErrorOnZeroBins(invPeakShapeEtabgLowMass2, analysisSettings.setErrorOnZeroBinsSigma);
		invPeakShapeEtabgLowMass2->Scale(scalePeakShapeEta);
		if (analysisResults.scalePeakShapeEtabgRel > 0) invPeakShapeEtabgLowMass2->Scale(analysisResults.scalePeakShapeEtabgRel);
		Float_t scalePeakShapeEtabgLowMass2 = 1.0;
		for (bin_stat_list_type::const_iterator iter = analysisResults.scalePeakShapeEtabgLowMass2.begin();iter != analysisResults.scalePeakShapeEtabgLowMass2.end();++iter) {
		    if (binParameters == (*iter).getParameters()) scalePeakShapeEtabgLowMass2 = (*iter).getValue();
		}
		invPeakShapeEtabgLowMass2->Scale(scalePeakShapeEtabgLowMass2);
	}
	if (invPeakShapeEtabg && invPeakShapeEtabgLowMass && analysisSettings.subtractBgLowMassPeakShapeEta) {
	    invPeakShapeEtabg->Add(invPeakShapeEtabgLowMass, -1);
	}
	if (invPeakShapeEtabg && invPeakShapeEtabgLowMass2 && analysisSettings.subtractBgLowMassPeakShapeEta) {
	    invPeakShapeEtabg->Add(invPeakShapeEtabgLowMass2, -1);
	}

	Float_t jetBgFraction_limitLow, jetBgFraction_limitHigh;
	Float_t jetBgFraction = getParameterValue(analysisSettings.jetBgFraction, binParameters, &jetBgFraction_limitLow, &jetBgFraction_limitHigh);
	if (jetBgFraction < 0.0) jetBgFraction = 0.0;
	if (jetBgFraction > 1.0) jetBgFraction = 1.0;
	if (jetBgFraction < jetBgFraction_limitLow) jetBgFraction_limitLow = jetBgFraction;
	if (jetBgFraction > jetBgFraction_limitHigh) jetBgFraction_limitHigh = jetBgFraction;
	if (analysisSettings.getJetBgFractionFromCandidatesPt) {
	    jetBgFraction = analysisResults.bgJetToTotalFromCandidates;
	    jetBgFraction_limitLow = jetBgFraction;
	    jetBgFraction_limitHigh = jetBgFraction;
	}

	histBeingFit = &invSig;
	histToFit1 = analysisSettings.fitBgDistribution ? invBg : 0;
	histToFit2 = analysisSettings.fitBgRandomDistribution ? invBgRandom : 0;
	histToFit3 = analysisSettings.fitBgLowMassDistribution ? invBgLowMass : 0;
	histToFit4 = analysisSettings.fitPeakShapeDistribution ? invPeakShape : 0;
	histToFit5 = analysisSettings.fitPeakShapeEtaDistribution ? invPeakShapeEta : 0;
	histToFit6 = analysisSettings.fitPeakShapeEtabgDistribution ? invPeakShapeEtabg : 0;
	fitRangeLeft = getParameterValue(analysisSettings.fitDistributionLeft, binParameters);
	fitRangeRight = getParameterValue(analysisSettings.fitDistributionRight, binParameters);
	fitRangeLeft2 = getParameterValue(analysisSettings.fitDistributionLeft2, binParameters);
	fitRangeRight2 = getParameterValue(analysisSettings.fitDistributionRight2, binParameters);
	TF1 fitDistributionFunc(analysisSettings.name + "_fitDistributionFunc", &getHistWeightedSum, (fitRangeLeft2 > 0) ? TMath::Min(fitRangeLeft, fitRangeLeft2) : fitRangeLeft, (fitRangeRight2 > 0) ? TMath::Max(fitRangeRight, fitRangeRight2) : fitRangeRight, 6);
	for (Int_t ipar = 0;ipar < 6;ipar++) fitDistributionFunc.SetParameter(ipar, 0.0);
	Float_t bgJetToTotalRatio = 0.0;
	Float_t bgJetToTotalRatioErr = 0.0;
	Float_t bgToSigbg_limitLow, bgToSigbg_limitHigh;
	Float_t bgToSigbg = getParameterValue(analysisSettings.bgToSigbg, binParameters, &bgToSigbg_limitLow, &bgToSigbg_limitHigh);
	Float_t bgToSigbgErr = 0.0;
	Float_t bgLowMassFitNorm = 0.0;
	Float_t bgLowMassFitNormErr = 0.0;
	Float_t temp;
	TH1F *invFitDistribution = 0;
	if (histToFit1 || histToFit2 || histToFit3 || histToFit4 || histToFit5 || histToFit6) {
	    fitDistributionFunc.SetParameter(0, bgToSigbg);
    	    if ((bgToSigbg_limitLow != noLimitMagic) && (bgToSigbg_limitHigh != noLimitMagic)) fitDistributionFunc.SetParLimits(0, bgToSigbg_limitLow, bgToSigbg_limitHigh);
		else fitDistributionFunc.SetParLimits(0, 0.0, 10.0);
	    fitDistributionFunc.SetParameter(1, jetBgFraction);
    	    if ((jetBgFraction_limitLow != noLimitMagic) && (jetBgFraction_limitHigh != noLimitMagic)) fitDistributionFunc.SetParLimits(1, jetBgFraction_limitLow, jetBgFraction_limitHigh);
	    if (histToFit1) {if (!histToFit2) fitDistributionFunc.FixParameter(1, 1.0);} else fitDistributionFunc.FixParameter(1, 0.0);
	    if (histToFit3) {fitDistributionFunc.SetParameter(2, 1.0); fitDistributionFunc.SetParLimits(2, 0.0, 1000.0);} else fitDistributionFunc.FixParameter(2, 0.0);
	    if (histToFit4) {fitDistributionFunc.SetParameter(3, 1.0); fitDistributionFunc.SetParLimits(3, 0.0, 10.0);} else fitDistributionFunc.FixParameter(3, 0.0);
	    if (histToFit5) {fitDistributionFunc.SetParameter(4, 1.0); fitDistributionFunc.SetParLimits(4, 0.0, 10.0);} else fitDistributionFunc.FixParameter(4, 0.0);
	    if (histToFit6) fitDistributionFunc.FixParameter(5, (1.0 - trueEtaBranchingRatio) / trueEtaBranchingRatio); else fitDistributionFunc.FixParameter(5, 0.0);

	    if (TString(analysisSettings.fitDistributionOption).Contains("U")) {
		//TVirtualFitter::Fitter(&invSig)->SetFCN((void (*)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t))0);
		delete TVirtualFitter::GetFitter();
		TVirtualFitter::SetFitter(0);
		TVirtualFitter::Fitter(&invSig)->SetFCN(&myHistHistFitChi2);
	    }
//cout << "fitDistributionOption=" << analysisSettings.fitDistributionOption << endl;
	    if (TString(analysisSettings.fitDistributionOption).Contains("UNORM")) {

	    Float_t sigNorm = 0;
	    Float_t sigNormErr = 0;
	    Float_t sigNorm2 = 0;
	    Float_t sigNorm2Err = 0;
	    Float_t bgMixNorm = 0;
	    Float_t bgMixNormErr = 0;
	    Float_t bgMixNorm2 = 0;
	    Float_t bgMixNorm2Err = 0;
	    Float_t bgJetmixNorm = 0;
	    Float_t bgJetmixNormErr = 0;
	    Float_t bgJetmixNorm2 = 0;
	    Float_t bgJetmixNorm2Err = 0;
	    if (fitRangeLeft < fitRangeRight) {
		extractSignal(invSig, funcZero, fitRangeLeft, fitRangeRight
		    , false, true, false, false, false
		    , sigNorm, sigNormErr
		    , temp, temp, temp, temp, temp, temp, temp, temp
		);
		extractSignal(*invBgRandom, funcZero, fitRangeLeft, fitRangeRight
		    , false, true, false, false, false
		    , bgMixNorm, bgMixNormErr
		    , temp, temp, temp, temp, temp, temp, temp, temp
		);
		extractSignal(*invBg, funcZero, fitRangeLeft, fitRangeRight
		    , false, true, false, false, false
		    , bgJetmixNorm, bgJetmixNormErr
		    , temp, temp, temp, temp, temp, temp, temp, temp
		);
	    }
	    if (fitRangeLeft2 < fitRangeRight2) {
		extractSignal(invSig, funcZero, fitRangeLeft2, fitRangeRight2
		    , false, true, false, false, false
		    , sigNorm2, sigNorm2Err
		    , temp, temp, temp, temp, temp, temp, temp, temp
		);
		extractSignal(*invBgRandom, funcZero, fitRangeLeft2, fitRangeRight2
		    , false, true, false, false, false
		    , bgMixNorm2, bgMixNorm2Err
		    , temp, temp, temp, temp, temp, temp, temp, temp
		);
		extractSignal(*invBg, funcZero, fitRangeLeft2, fitRangeRight2
		    , false, true, false, false, false
		    , bgJetmixNorm2, bgJetmixNorm2Err
		    , temp, temp, temp, temp, temp, temp, temp, temp
		);
	    }
	    Float_t sigN = sigNorm+sigNorm2;
	    Float_t sigNe = TMath::Sqrt(sigNormErr*sigNormErr + sigNorm2Err*sigNorm2Err);
	    Float_t bgMixN = bgMixNorm+bgMixNorm2;
	    Float_t bgMixNe = TMath::Sqrt(bgMixNormErr*bgMixNormErr + bgMixNorm2Err*bgMixNorm2Err);
	    Float_t bgJetmixN = bgJetmixNorm+bgJetmixNorm2;
	    Float_t bgJetmixNe = TMath::Sqrt(bgJetmixNormErr*bgJetmixNormErr + bgJetmixNorm2Err*bgJetmixNorm2Err);
	    Float_t newAlpha = 1. / (1. - ((bgJetmixNorm*bgMixNorm2 - bgJetmixNorm2*bgMixNorm) / (bgMixNorm*sigNorm2 - bgMixNorm2*sigNorm)));	    
//cout << "jetBgFraction = " << jetBgFraction << ", new = " << newAlpha << endl;
    	    if ((jetBgFraction_limitLow != noLimitMagic) && (jetBgFraction_limitHigh != noLimitMagic)) {
		newAlpha = TMath::Max(newAlpha, jetBgFraction_limitLow);
		newAlpha = TMath::Min(newAlpha, jetBgFraction_limitHigh);
	    }
	    //newAlpha = jetBgFraction;
	    Float_t A = (sigN) / (newAlpha * bgJetmixN + (1 - newAlpha) * bgMixN);
	    Float_t Aerr = TMath::Sqrt((sigNe*sigNe/(sigN*sigN)) + (newAlpha*newAlpha*bgJetmixNe*bgJetmixNe/(bgJetmixN*bgJetmixN)) 
		+ ((1-newAlpha)*(1-newAlpha)*bgMixNe*bgMixNe/(bgMixN*bgMixN)));
    	    if ((bgToSigbg_limitLow != noLimitMagic) && (bgToSigbg_limitHigh != noLimitMagic)) {
		A = TMath::Max(A, bgToSigbg_limitLow);
		A = TMath::Min(A, bgToSigbg_limitHigh);
	    }
//cout << "pT " << binParameters.min << "-" << binParameters.max << ": sigNorm=" << sigNorm << ", sigNorm2=" << sigNorm2 << ", bgJetmixNorm=" << bgJetmixNorm << ", bgJetmixNorm2=" << bgJetmixNorm2 << ", A=" << A << ", alpha=" << newAlpha << endl;
	    fitDistributionFunc.SetParameter(0, A);
	    fitDistributionFunc.SetParError(0, Aerr);
	    fitDistributionFunc.SetParameter(1, newAlpha);
	    fitDistributionFunc.SetParError(1, 0);
	    
	    } else {
	    invSig.Fit(&fitDistributionFunc, analysisSettings.fitDistributionOption);
	    }
	    
	    if (TString(analysisSettings.fitDistributionOption).Contains("U")) {
		delete TVirtualFitter::GetFitter();
		TVirtualFitter::SetFitter(0);
	    }

	    bgJetToTotalRatio = fitDistributionFunc.GetParameter(1);
	    bgJetToTotalRatioErr = fitDistributionFunc.GetParError(1);
	    bgToSigbg = fitDistributionFunc.GetParameter(0);
	    bgToSigbgErr = fitDistributionFunc.GetParError(0);
	    bgLowMassFitNorm = fitDistributionFunc.GetParameter(2);
	    bgLowMassFitNormErr = fitDistributionFunc.GetParError(2);
	    invFitDistribution = new TH1F(invSig);
	    if (invFitDistribution) {
		invFitDistribution->SetName(analysisSettings.name + "_invFitDistribution_" + invSig.GetName());
		invFitDistribution->Reset();
		if (analysisSettings.fitBgDistribution && invBg) invFitDistribution->Add(invBg, bgToSigbg * bgJetToTotalRatio);
		if (analysisSettings.fitBgRandomDistribution && invBgRandom) invFitDistribution->Add(invBgRandom, bgToSigbg * (1.0 - bgJetToTotalRatio));
		if (analysisSettings.fitBgLowMassDistribution && invBgLowMass) invFitDistribution->Add(invBgLowMass, bgLowMassFitNorm);
		//if (histToFit5 && histToFit6) invFitDistribution->Add(histToFit6, fitDistributionFunc.GetParameter(4) * fitDistributionFunc.GetParameter(5));
		if (analysisSettings.subtractBgFitDistribution) {
		    //cout << invSig.GetName() << " subtracting bg fit distr " << invSig.Integral();
		    invSig.Add(invFitDistribution, -1);
		    //cout << " " << invSig.Integral() << endl;
		}
	    }
	}
	histBeingFit = 0; histToFit1 = 0; histToFit2 = 0; histToFit3 = 0; histToFit4 = 0; histToFit5 = 0; histToFit6 = 0;

	Float_t scaleBgBack = 1.0;
	TH1F *invBgBack = s.distributionBgBack ? (new TH1F(*(s.distributionBgBack)->getDistribution())) : 0;
	if (invBgBack) {
		invBgBack->SetName(analysisSettings.name + "_invBgBack_" + invBgBack->GetName());
		if (nrebin > 1) invBgBack->Rebin(nrebin);
		if (doSmearOutliers) smearOutliers(invBgBack, smearOutliersSigma, smearOutliersNSigma, zeroOutliers);
		if (analysisSettings.smoothBgBack > 0) {
		    smoothHistGaus(invBgBack, TMath::Abs(analysisSettings.smoothBgBack), smoothGausNSigma);
		    //invBgBack->Smooth(analysisSettings.smoothBgBack, -1, -1);
		} else if (analysisSettings.smoothBgBack < 0) {
		    smoothHistPol2(invBgBack, (Int_t)TMath::Abs(analysisSettings.smoothBgBack), analysisSettings.smoothFitOption);
		}
		if (analysisSettings.setErrorOnZeroBinsBgBack) setErrorOnZeroBins(invBgBack, analysisSettings.setErrorOnZeroBinsSigma);
		if (analysisSettings.mixNormFixedBack > 0) {
		    Float_t corr = analysisSettings.mixNormFixedBack;
		    if (analysisSettings.mixNormFixedBack > 0) scaleBgBack *= corr;
		}
		if (analysisSettings.mixNormMassRangeBack) {
		    Float_t intBgBack = invBgBack->Integral(invBgBack->GetXaxis()->FindBin(analysisSettings.mixNormLeftBack), invBgBack->GetXaxis()->FindBin(analysisSettings.mixNormRightBack));
		    Float_t intSig = invSig.Integral(invSig.GetXaxis()->FindBin(analysisSettings.mixNormLeftBack), invSig.GetXaxis()->FindBin(analysisSettings.mixNormRightBack));
		    Float_t corr = (intBgBack != 0) ? (intSig / intBgBack) : 1.0;
		    if (analysisSettings.mixNormMassRangeBack) scaleBgBack *= corr;
		}
		if (analysisSettings.mixNormEntriesBack) {
		    Float_t intBgBack = invBgBack->Integral();
		    Float_t intSig = invSig.Integral();
		    Float_t corr = (intBgBack != 0) ? (intSig / intBgBack) : 1.0;
		    if (analysisSettings.mixNormEntriesBack) scaleBgBack *= corr;
		}
		invBgBack->Scale(scaleBgBack);
		if (analysisSettings.subtractBgBack) {
		    //cout << invSig.GetName() << " subtracting bg back" << endl;
		    invSig.Add(invBgBack, -1);
		}
	}

	Float_t fitParameter6_limitLow, fitParameter6_limitHigh;
	Float_t fitParameter6 = getParameterValue(analysisSettings.fitParameter6, binParameters, &fitParameter6_limitLow, &fitParameter6_limitHigh);
	Float_t fitParameter7_limitLow, fitParameter7_limitHigh;
	Float_t fitParameter7 = getParameterValue(analysisSettings.fitParameter7, binParameters, &fitParameter7_limitLow, &fitParameter7_limitHigh);
	Float_t fitParameters_limitLow, fitParameters_limitHigh;
	Float_t fitParameters = getParameterValue(analysisSettings.fitParameters, binParameters, &fitParameters_limitLow, &fitParameters_limitHigh);

	TF1 fitFunc(analysisSettings.name + "_fitFunc", (const Char_t*)analysisSettings.fitFunctionStr);
	fitFunc.SetParameter(0, invSig.Integral());
	fitFunc.SetParLimits(0, 0, invSig.Integral() * 2.0);
	fitFunc.SetParameter(1, peakPosition);
	if ((peakPosition_limitLow != noLimitMagic) && (peakPosition_limitHigh != noLimitMagic)) fitFunc.SetParLimits(1, peakPosition_limitLow, peakPosition_limitHigh);
	fitFunc.SetParameter(2, peakWidth);
	if ((peakWidth_limitLow != noLimitMagic) && (peakWidth_limitHigh != noLimitMagic)) fitFunc.SetParLimits(2, peakWidth_limitLow, peakWidth_limitHigh);
	fitFunc.SetParameter(6, fitParameter6);
	if ((fitParameter6_limitLow != noLimitMagic) && (fitParameter6_limitHigh != noLimitMagic)) fitFunc.SetParLimits(6, fitParameter6_limitLow, fitParameter6_limitHigh);
	fitFunc.SetParameter(7, fitParameter7);
	if ((fitParameter7_limitLow != noLimitMagic) && (fitParameter7_limitHigh != noLimitMagic)) fitFunc.SetParLimits(7, fitParameter7_limitLow, fitParameter7_limitHigh);
    	for (Int_t ipar = 8;ipar < fitFunc.GetNpar();ipar++) {
	    fitFunc.SetParameter(ipar, fitParameters);
	    if ((fitParameters_limitLow != noLimitMagic) && (fitParameters_limitHigh != noLimitMagic)) fitFunc.SetParLimits(ipar, fitParameters_limitLow, fitParameters_limitHigh);
	}
	fitFunc.SetRange(fitLeft, fitRight);
	fitFunc.FixParameter(3, 0);
	fitFunc.FixParameter(4, 0);
	fitFunc.FixParameter(5, 0);
	invSig.Fit(&fitFunc, analysisSettings.fitOption);
	Float_t mean = fitFunc.GetParameter(1);
	Float_t meanErr = fitFunc.GetParError(1);
	Float_t width = fitFunc.GetParameter(2);
	Float_t widthErr = fitFunc.GetParError(2);
	fitFunc.SetRange(fitLeftEta, fitRightEta);
	fitFunc.FixParameter(0, fitFunc.GetParameter(0));
	fitFunc.FixParameter(1, fitFunc.GetParameter(1));
	fitFunc.FixParameter(2, fitFunc.GetParameter(2));
	fitFunc.SetParameter(3, invSig.Integral());
	fitFunc.SetParLimits(3, 0, invSig.Integral() * 2.0);
	fitFunc.SetParameter(4, peakPositionEta);
	if ((peakPositionEta_limitLow != noLimitMagic) && (peakPositionEta_limitHigh != noLimitMagic)) fitFunc.SetParLimits(4, peakPositionEta_limitLow, peakPositionEta_limitHigh);
	fitFunc.SetParameter(5, peakWidthEta);
	if ((peakWidthEta_limitLow != noLimitMagic) && (peakWidthEta_limitHigh != noLimitMagic)) fitFunc.SetParLimits(5, peakWidthEta_limitLow, peakWidthEta_limitHigh);
	invSig.Fit(&fitFunc, analysisSettings.fitOption);
	Float_t meanEta = fitFunc.GetParameter(4);
	Float_t meanEtaErr = fitFunc.GetParError(4);
	Float_t widthEta = fitFunc.GetParameter(5);
	Float_t widthEtaErr = fitFunc.GetParError(5);
	fitFunc.SetRange(fitLeft, fitRightEta);

	TF1 funcPurePeaks(fitFunc);
	funcPurePeaks.SetName(analysisSettings.name + "_funcPurePeaks_" + s.distribution->GetName());
	for (Int_t ipar = 6;ipar < funcPurePeaks.GetNpar();ipar++) funcPurePeaks.SetParameter(ipar, 0);

	TF1 funcPurePeak(fitFunc);
	funcPurePeak.SetName(analysisSettings.name + "_funcPurePeak_" + s.distribution->GetName());
	for (Int_t ipar = 6;ipar < funcPurePeak.GetNpar();ipar++) funcPurePeak.SetParameter(ipar, 0);
	funcPurePeak.SetParameter(3, 0);
	funcPurePeak.SetRange(fitLeft, fitRight);

	TF1 funcPurePeakEta(fitFunc);
	funcPurePeakEta.SetName(analysisSettings.name + "_funcPurePeakEta_" + s.distribution->GetName());
	for (Int_t ipar = 6;ipar < funcPurePeakEta.GetNpar();ipar++) funcPurePeakEta.SetParameter(ipar, 0);
	funcPurePeakEta.SetParameter(0, 0);
	funcPurePeakEta.SetRange(fitLeftEta, fitRightEta);

	TF1 funcPureBg(fitFunc);
	funcPureBg.SetName(analysisSettings.name + "_funcPureBg_" + s.distribution->GetName());
	funcPureBg.SetParameter(0, 0);
	funcPureBg.SetParameter(3, 0);

	TF1 funcBgPeak(fitFunc);
	funcBgPeak.SetName(analysisSettings.name + "_funcBgPeak_" + s.distribution->GetName());
	funcBgPeak.SetParameter(3, 0);

	TF1 funcBgPeakEta(fitFunc);
	funcBgPeakEta.SetName(analysisSettings.name + "_funcBgPeakEta_" + s.distribution->GetName());
	funcBgPeakEta.SetParameter(0, 0);

	TH1F invRes(invSig);
	invRes.SetName(analysisSettings.name + "_invRes_" + s.distribution->GetName());
	invRes.Add(&fitFunc, -1);

	TH1F invResBg(invSig);
	invResBg.SetName(analysisSettings.name + "_invResBg_" + s.distribution->GetName());
	invResBg.Add(&funcPurePeaks, -1);

	TH1F invResSig(invSig);
	invResSig.SetName(analysisSettings.name + "_invResSig_" + s.distribution->GetName());
	invResSig.Add(&funcPureBg, -1);

	//Float_t height = fitFunc.GetParameter(0);
	//Float_t heightErr = fitFunc.GetParError(0);
	//Float_t heightEta = fitFunc.GetParameter(3);
	//Float_t heightEtaErr = fitFunc.GetParError(3);
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
	if (true) {
	    Float_t a = bgPar7;
	    Float_t aErr = bgPar7Err;
	    Float_t b = peakRight;
	    areaFractionBg = TMath::Exp(-a*b) * ((a*b) + 1);
	    areaFractionBgErr = aErr * (a*b*b*TMath::Exp(-a*b));
	}


	Float_t areaSigBg = 0;
	Float_t areaSigBgErr = 0;
	Float_t areaSigBgFull = 0;
	Float_t areaSigBgFullErr = 0;
	Float_t areaSigBgFullRange = 0;
	Float_t areaSigBgFullRangeErr = 0;
	Float_t area = 0;
	Float_t areaErr = 0;
	Float_t areaBg = 0;
	Float_t areaBgErr = 0;
	Float_t meanData = 0;
	Float_t meanDataErr = 0;
	Float_t widthData = 0;
	Float_t widthDataErr = 0;
	Float_t lowmassYield = 0;
	Float_t lowmassYieldErr = 0;
	Float_t etabgYield = 0;
	Float_t etabgYieldErr = 0;
	Float_t areaPeakShape = 0;
	Float_t areaPeakShapeErr = 0;
	Float_t yieldBgJetmix = 0;
	Float_t yieldBgJetmixErr = 0;
	Float_t yieldBgJetmixNorm = 0;
	Float_t yieldBgJetmixNormErr = 0;
	Float_t yieldBgJetmixNorm2 = 0;
	Float_t yieldBgJetmixNorm2Err = 0;
	Float_t yieldBgMix = 0;
	Float_t yieldBgMixErr = 0;
	Float_t yieldBgMixNorm = 0;
	Float_t yieldBgMixNormErr = 0;
	Float_t yieldBgMixNorm2 = 0;
	Float_t yieldBgMixNorm2Err = 0;

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

	extractSignal(invSigBg, funcZero, peakLeftYield, peakRightYield
	    , false, true, false, false, false
	    , areaSigBgFull, areaSigBgFullErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	areaSigBgFullRange = invSigBg.Integral();
	areaSigBgFullRangeErr = 0;
	/*extractSignal(invSigBg, funcZero, invSigBg.GetXaxis()->GetBinLowEdge(1), invSigBg.GetXaxis()->GetBinUpEdge(invSigBg.GetXaxis()->GetNbins())
	    , false, true, false, false, false
	    , areaSigBgFullRange, areaSigBgFullRangeErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);*/
	if (invBg) extractSignal(*invBg, funcZero, peakLeftYield, peakRightYield
	    , false, true, false, false, false
	    , yieldBgJetmix, yieldBgJetmixErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);
	if (invBg && (fitRangeLeft < fitRangeRight)) extractSignal(*invBg, funcZero, fitRangeLeft, fitRangeRight
	    , false, true, false, false, false
	    , yieldBgJetmixNorm, yieldBgJetmixNormErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);
	if (invBg && (fitRangeLeft2 < fitRangeRight2)) extractSignal(*invBg, funcZero, fitRangeLeft2, fitRangeRight2
	    , false, true, false, false, false
	    , yieldBgJetmixNorm2, yieldBgJetmixNorm2Err
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invBgRandom) extractSignal(*invBgRandom, funcZero, peakLeftYield, peakRightYield
	    , false, true, false, false, false
	    , yieldBgMix, yieldBgMixErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);
	if (invBgRandom && (fitRangeLeft < fitRangeRight)) extractSignal(*invBgRandom, funcZero, fitRangeLeft, fitRangeRight
	    , false, true, false, false, false
	    , yieldBgMixNorm, yieldBgMixNormErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);
	if (invBgRandom && (fitRangeLeft2 < fitRangeRight2)) extractSignal(*invBgRandom, funcZero, fitRangeLeft2, fitRangeRight2
	    , false, true, false, false, false
	    , yieldBgMixNorm2, yieldBgMixNorm2Err
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invBgLowMass) extractSignal(*invBgLowMass, funcZero, peakLeftYield, peakRightYield
	    , false, true, false, false, false
	    , lowmassYield, lowmassYieldErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invPeakShape) extractSignal(*invPeakShape, funcZero, peakLeftYield, peakRightYield
	    , false, true, false, false, false
	    , areaPeakShape, areaPeakShapeErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invPeakShapeEtabg) extractSignal(*invPeakShapeEtabg, funcZero, peakLeftYield, peakRightYield
	    , false, true, false, false, false
	    , etabgYield, etabgYieldErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	Float_t sigBg = getParameterValue(analysisSettings.signalToBackgroundRatio, binParameters);
	if ((sigBg > 0) && (analysisSettings.signalToBackgroundRatio.wasFit)) {
		/*
		Float_t area1 = areaSigBgFull / (1.0 + (1.0 / sigBg));
		areaErr = (area != 0) ? (areaErr * area1 / area) : 0;
		area = area1;
		Float_t areaBg1 = areaSigBgFull / (1.0 + sigBg);
		areaBgErr = (areaBg != 0) ? (areaBgErr * areaBg1 / areaBg) : 0;
		areaBg = areaBg1;
		*/
		area = areaSigBgFull / (1.0 + (1.0 / sigBg));
		areaErr = areaSigBgFullErr / (1.0 + (1.0 / sigBg));
		//areaBg = areaSigBgFull / (1.0 + sigBg);
		//areaBgErr = areaSigBgFullErr / (1.0 + sigBg);		
	}

	Float_t yieldToAreaSigBgFull = (areaSigBgFull != 0) ? (area / areaSigBgFull) : 0;
	Float_t yieldToAreaSigBgFullErr = ((yieldToAreaSigBgFull != 0) && (area != 0) && (areaSigBgFull != 0)) ? (yieldToAreaSigBgFull * TMath::Sqrt((areaErr * areaErr / (area * area)) + (areaSigBgFullErr * areaSigBgFullErr / (areaSigBgFull * areaSigBgFull)))) : 0;

	Float_t yieldToAreaSigBgFullRange = (areaSigBgFullRange != 0) ? (area / areaSigBgFullRange) : 0;
	Float_t yieldToAreaSigBgFullRangeErr = ((yieldToAreaSigBgFullRange != 0) && (area != 0) && (areaSigBgFullRange != 0)) ? (yieldToAreaSigBgFullRange * TMath::Sqrt((areaErr*areaErr/(area*area)) + (areaSigBgFullRangeErr*areaSigBgFullRangeErr/(areaSigBgFullRange*areaSigBgFullRange)))) : 0;

	Float_t yieldToBgFull = (area != areaSigBgFull) ? (area / (areaSigBgFull - area)) : 0;
	Float_t yieldToBgFullErr = (yieldToAreaSigBgFull != 0) ? (yieldToAreaSigBgFullErr * yieldToBgFull / yieldToAreaSigBgFull) : 0;

	Float_t lowmassToYield = (area != 0) ? (lowmassYield / area) : 0;
	Float_t lowmassToYieldErr = ((area != 0) && (lowmassYield != 0)) ? (lowmassToYield * TMath::Sqrt(((areaErr/area)*(areaErr/area)) + ((lowmassYieldErr/lowmassYield)*(lowmassYieldErr/lowmassYield)))) : 0;

	Float_t areaSigBgEta = 0;
	Float_t areaSigBgEtaErr = 0;
	Float_t areaSigBgEtaFull = 0;
	Float_t areaSigBgEtaFullErr = 0;
	Float_t areaEta = 0;
	Float_t areaEtaErr = 0;
	Float_t areaBgEta = 0;
	Float_t areaBgEtaErr = 0;
	Float_t meanDataEta = 0;
	Float_t meanDataEtaErr = 0;
	Float_t widthDataEta = 0;
	Float_t widthDataEtaErr = 0;
	Float_t areaPeakShapeEta = 0;
	Float_t areaPeakShapeEtaErr = 0;
	Float_t yieldBgJetmixEta = 0;
	Float_t yieldBgJetmixEtaErr = 0;
	Float_t yieldBgMixEta = 0;
	Float_t yieldBgMixEtaErr = 0;

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

	extractSignal(invSigBg, funcZero, peakLeftYieldEta, peakRightYieldEta
	    , false, true, false, false, false
	    , areaSigBgEtaFull, areaSigBgEtaFullErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invBg) extractSignal(*invBg, funcZero, peakLeftYieldEta, peakRightYieldEta
	    , false, true, false, false, false
	    , yieldBgJetmixEta, yieldBgJetmixEtaErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invBgRandom) extractSignal(*invBgRandom, funcZero, peakLeftYieldEta, peakRightYieldEta
	    , false, true, false, false, false
	    , yieldBgMixEta, yieldBgMixEtaErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	if (invPeakShapeEta) extractSignal(*invPeakShapeEta, funcZero, peakLeftYieldEta, peakRightYieldEta
	    , false, true, false, false, false
	    , areaPeakShapeEta, areaPeakShapeEtaErr
	    , temp, temp, temp, temp, temp, temp, temp, temp
	);

	Float_t sigBgEta = getParameterValue(analysisSettings.signalToBackgroundRatioEta, binParameters);
	if ((sigBgEta > 0) && (analysisSettings.signalToBackgroundRatioEta.wasFit)) {
		/*
		Float_t areaEta1 = areaSigBgEtaFull / (1.0 + (1.0 / sigBgEta));
		areaEtaErr = (areaEta != 0) ? (areaEtaErr * areaEta1 / areaEta) : 0;
		areaEta = areaEta1;
		Float_t areaBgEta1 = areaSigBgEtaFull / (1.0 + sigBgEta);
		areaBgEtaErr = (areaBgEta != 0) ? (areaBgEtaErr * areaBgEta1 / areaBgEta) : 0;
		areaBgEta = areaBgEta1;
		*/
		areaEta = areaSigBgEtaFull / (1.0 + (1.0 / sigBgEta));
		areaEtaErr = areaSigBgEtaFullErr / (1.0 + (1.0 / sigBgEta));
		//areaBgEta = areaSigBgEtaFull / (1.0 + sigBgEta);
		//areaBgEtaErr = areaSigBgEtaFullErr / (1.0 + sigBgEta);
	}

	Float_t yieldToAreaSigBgEtaFull = (areaSigBgEtaFull != 0) ? (areaEta / areaSigBgEtaFull) : 0;
	Float_t yieldToAreaSigBgEtaFullErr = ((yieldToAreaSigBgEtaFull != 0) && (areaEta != 0) && (areaSigBgEtaFull != 0)) ? (yieldToAreaSigBgEtaFull * TMath::Sqrt((areaEtaErr * areaEtaErr / (areaEta * areaEta)) + (areaSigBgEtaFullErr * areaSigBgEtaFullErr / (areaSigBgEtaFull * areaSigBgEtaFull)))) : 0;

	Float_t yieldToBgEtaFull = (areaEta != areaSigBgEtaFull) ? (areaEta / (areaSigBgEtaFull - areaEta)) : 0;
	Float_t yieldToBgEtaFullErr = (yieldToAreaSigBgEtaFull != 0) ? (yieldToAreaSigBgEtaFullErr * yieldToBgEtaFull / yieldToAreaSigBgEtaFull) : 0;

	if (analysisSettings.normPeakShapeYieldShow && invPeakShape) {
		Float_t integral = areaPeakShape;
		Float_t scalePeakShapeYield = (integral != 0) ? (area / integral) : 1.0;
		invPeakShape->Scale(scalePeakShapeYield);
	}

	if (analysisSettings.normPeakShapeEtaYieldShow && invPeakShapeEta) {
		Float_t integral = areaPeakShapeEta;
		Float_t scalePeakShapeEtaYield = (integral != 0) ? (areaEta / integral) : 1.0;
		invPeakShapeEta->Scale(scalePeakShapeEtaYield);
		if (invPeakShapeEtabg) invPeakShapeEtabg->Scale(scalePeakShapeEtaYield);
	}

	Float_t etabgToYield = ((area != 0) && (areaPeakShapeEta != 0)) ? ((etabgYield * (areaEta / areaPeakShapeEta)) / area) : 0;
	Float_t etabgToYieldErr = etabgToYield * TMath::Sqrt(
	    ((area != 0) ? (areaErr/area)*(areaErr/area) : 0) + 
	    ((etabgYield != 0) ? (etabgYieldErr/etabgYield)*(etabgYieldErr/etabgYield) : 0) +
	    ((areaEta != 0) ? (areaEtaErr/areaEta)*(areaEtaErr/areaEta) : 0) +
	    ((areaPeakShapeEta != 0) ? (areaPeakShapeEtaErr/areaPeakShapeEta)*(areaPeakShapeEtaErr/areaPeakShapeEta) : 0)
	);

	yieldBgJetmix *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixErr *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixEta *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixEtaErr *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixNorm *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixNormErr *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixNorm2 *= bgToSigbg * bgJetToTotalRatio;
	yieldBgJetmixNorm2Err *= bgToSigbg * bgJetToTotalRatio;
	yieldBgMix *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixErr *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixEta *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixEtaErr *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixNorm *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixNormErr *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixNorm2 *= bgToSigbg * (1 - bgJetToTotalRatio);
	yieldBgMixNorm2Err *= bgToSigbg * (1 - bgJetToTotalRatio);
	Float_t mixJetmixDiff = 0;
	Float_t mixJetmixDiffErr = 0;
	Float_t mixJetmixDiffEta = 0;
	Float_t mixJetmixDiffEtaErr = 0;
	{
	Float_t jetMixOnly =  ((yieldBgJetmixNorm + yieldBgJetmixNorm2) != 0) ? yieldBgJetmix * (yieldBgJetmixNorm + yieldBgJetmixNorm2 + yieldBgMixNorm + yieldBgMixNorm2) / (yieldBgJetmixNorm + yieldBgJetmixNorm2) : 0;
	Float_t mixOnly =  ((yieldBgMixNorm + yieldBgMixNorm2) != 0) ? yieldBgMix * (yieldBgJetmixNorm + yieldBgJetmixNorm2 + yieldBgMixNorm + yieldBgMixNorm2) / (yieldBgMixNorm + yieldBgMixNorm2) : 0;
	mixJetmixDiff = ((yieldBgJetmix + yieldBgMix) != 0) ? ((TMath::Abs(jetMixOnly - mixOnly) / 2.0) / (yieldBgJetmix + yieldBgMix)) : 0;
	mixJetmixDiffErr = 0;
	}
	{
	Float_t jetMixOnlyEta =  ((yieldBgJetmixNorm + yieldBgJetmixNorm2) != 0) ? yieldBgJetmixEta * (yieldBgJetmixNorm + yieldBgJetmixNorm2 + yieldBgMixNorm + yieldBgMixNorm2) / (yieldBgJetmixNorm + yieldBgJetmixNorm2) : 0;
	Float_t mixOnlyEta =  ((yieldBgMixNorm + yieldBgMixNorm2) != 0) ? yieldBgMixEta * (yieldBgJetmixNorm + yieldBgJetmixNorm2 + yieldBgMixNorm + yieldBgMixNorm2) / (yieldBgMixNorm + yieldBgMixNorm2) : 0;
	mixJetmixDiffEta = ((yieldBgJetmixEta + yieldBgMixEta) != 0) ? ((TMath::Abs(jetMixOnlyEta - mixOnlyEta) / 2.0) / (yieldBgJetmixEta + yieldBgMixEta)) : 0;
	mixJetmixDiffEtaErr = 0;
	}

	TH1F invResPeakShapeBg(invSig);
	invResPeakShapeBg.SetName(analysisSettings.name + "_invResPeakShapeBg_" + s.distribution->GetName());
	if (invPeakShape) invResPeakShapeBg.Add(invPeakShape, -1);
	if (invPeakShapeEta) invResPeakShapeBg.Add(invPeakShapeEta, -1);

	if (analysisSettings.showBins) {
		Int_t nyp = Int_t(sqrt(Float_t(numBinsToShow)));
		Int_t nxp = Int_t((Float_t(numBinsToShow) - 0.0001) / Float_t(nyp)) + 1;
		if (!canvas) {
			canvas = new TCanvas(analysisSettings.name, analysisSettings.title);
			if (canvas) {
			    Int_t xsize = nxp * 400;
			    Int_t ysize = nyp * 400;
			    canvas->SetCanvasSize(xsize, ysize);
			    canvas->Divide(nxp, nyp, 0.00, 0.00);
			}
		}
		if (canvas) canvas->cd(binToShow + 1);
		Bool_t lastRow = (binToShow >= ((nyp - 1) * nxp));
		Bool_t firstCol = ((binToShow % nxp) == 0);

		Bool_t scaleVert = false;

		if (scaleVert) {
		    invSigBg.Scale(1.0 / invSig.GetMaximum());
		    invSig.Scale(1.0 / invSig.GetMaximum());
		}
		
		Float_t textSize = 0.055;

		Float_t max = invSigBg.GetMaximum();
		Float_t min = invSigBg.GetMinimum();
		Float_t maxSig = invSig.GetMaximum();
		Float_t minSig = invSig.GetMinimum();
		TLegend *legend = 0;
		if (binToShow == 0) {
		    legend = new TLegend(0.4, 0.4, 0.89, 0.89);//(analysisSettings.binRangeLeft * 0.5) + (analysisSettings.binRangeRight * 0.5), invSigBg.GetMaximum() * 0.65, analysisSettings.binRangeRight, invSigBg.GetMaximum() * 0.3);
		}
		if (!analysisSettings.showSigBg) {
		    invSigBg.Reset();
		}
		{
		invSigBg.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
		if (analysisSettings.showSigBg) {
			invSigBg.GetYaxis()->SetRangeUser(min - (TMath::Abs(max - min) * 0.18), max + (TMath::Abs(max - min) * 0.18));
		} else if (analysisSettings.showSig) {
			Float_t minSigVal = minSig;
			minSigVal = 0;
			invSigBg.GetYaxis()->SetRangeUser(minSigVal - (TMath::Abs(maxSig - minSigVal) * 0.18), maxSig + (TMath::Abs(maxSig - minSigVal) * 0.18));
		}
		invSigBg.GetXaxis()->SetTitle("#font[12]{M}_{#gamma#gamma} [GeV/#font[12]{c}^{2}]");
		if (scaleVert) {
		    invSigBg.GetYaxis()->SetTitle("");
		} else {
		    invSigBg.GetYaxis()->SetTitle("Counts");
		}
		invSigBg.GetXaxis()->SetTitleOffset(1.2);
		invSigBg.GetYaxis()->SetTitleOffset(1.2);
		if (!firstCol) {
		    invSigBg.GetYaxis()->SetTitle("");
		    invSigBg.GetYaxis()->LabelsOption("");
		}
		if (!lastRow) {
		    invSigBg.GetXaxis()->SetTitle("");
		    invSigBg.GetXaxis()->LabelsOption("");
		}
		//invSigBg.GetXaxis()->SetTitleSize(textSize);
		//invSigBg.GetXaxis()->SetLabelSize(textSize);
		//invSigBg.GetYaxis()->SetTitleSize(textSize);
		//invSigBg.GetYaxis()->SetLabelSize(textSize);
		invSigBg.SetTitle("");
			//invSigBg.SetMinimum(0);
			invSigBg.SetLineColor(analysisSettings.colorSigBg);
			invSigBg.SetMarkerColor(analysisSettings.colorSigBg);
			//invSigBg.SetFillColor(analysisSettings.colorSigBg);
			TH1 *newHist = invSigBg.DrawCopy("");
			if (newHist) {
			    newHist->GetXaxis()->SetLabelFont(gStyle->GetLabelFont("x"));
			    newHist->GetYaxis()->SetLabelFont(gStyle->GetLabelFont("y"));
			    newHist->GetXaxis()->SetTitleFont(gStyle->GetTitleFont("x"));
			    newHist->GetYaxis()->SetTitleFont(gStyle->GetTitleFont("y"));
			}
			if (analysisSettings.showSigBg) {
			    if (legend && newHist) legend->AddEntry(newHist, "#font[12]{M}_{#gamma#gamma} distribution", "L E");
			}
		}
		if (invBg && analysisSettings.showBg) {
			invBg->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBg->SetLineColor(analysisSettings.colorBg);
			invBg->SetMarkerColor(analysisSettings.colorBg);
			//invBg->SetFillColor(analysisSettings.colorBg);
			TH1 *newHist = invBg->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Jet mixing", "L");
		}
		if (invBgRandom && analysisSettings.showBgRandom) {
			invBgRandom->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBgRandom->SetLineColor(analysisSettings.colorBgRandom);
			invBgRandom->SetMarkerColor(analysisSettings.colorBgRandom);
			//invBgRandom->SetFillColor(analysisSettings.colorBgRandom);
			TH1 *newHist = invBgRandom->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Random mixing", "L");
		}
		if (invBgBack && analysisSettings.showBgBack) {
			invBgBack->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBgBack->SetLineColor(analysisSettings.colorBgBack);
			invBgBack->SetMarkerColor(analysisSettings.colorBgBack);
			//invBgBack->SetFillColor(analysisSettings.colorBgBack);
			TH1 *newHist = invBgBack->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Jet back-to-back", "L");
		}
		if (invFitDistribution && analysisSettings.showBgFitDistribution) {
			invFitDistribution->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invFitDistribution->SetLineColor(analysisSettings.colorBgFitDistribution);
			TH1 *newHist = invFitDistribution->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Combinatorial bg", "L");
		}
		if (analysisSettings.showSig) {
			TH1F invSigPi0Peak(invSig);
			TString name = invSigPi0Peak.GetName();
			name += "_pi0peak";
			invSigPi0Peak.SetName(name.Data());
			invSigPi0Peak.GetXaxis()->SetRangeUser(peakLeftYield, peakRightYield);
			invSigPi0Peak.SetLineColor(analysisSettings.fillColor);
			invSigPi0Peak.SetFillColor(analysisSettings.fillColor);
			invSigPi0Peak.SetMarkerColor(analysisSettings.fillColor);
			invSigPi0Peak.SetLineWidth(0);
			TH1 *newHist = invSigPi0Peak.DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "#pi^{0} peak", "L E");
		}
		if (analysisSettings.showSig) {
			TH1F invSigEtaPeak(invSig);
			TString name = invSigEtaPeak.GetName();
			name += "_etapeak";
			invSigEtaPeak.SetName(name.Data());
			invSigEtaPeak.GetXaxis()->SetRangeUser(peakLeftYieldEta, peakRightYieldEta);
			invSigEtaPeak.SetLineColor(analysisSettings.fillColor);
			invSigEtaPeak.SetFillColor(analysisSettings.fillColor);
			invSigEtaPeak.SetMarkerColor(analysisSettings.fillColor);
			invSigEtaPeak.SetLineWidth(0);
			TH1 *newHist = invSigEtaPeak.DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "#eta peak", "L E");
		}
		if (analysisSettings.showSig) {
			invSig.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invSig.SetLineColor(analysisSettings.color);
			invSig.SetMarkerColor(analysisSettings.color);
			//invSig.SetLineWidth(3);
			TH1 *newHist = invSig.DrawCopy("SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Background subtracted", "L E");
		}
		if (invBgLowMass && analysisSettings.showBgLowMass) {
			invBgLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invBgLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invBgLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			//invBgLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			TH1 *newHist = invBgLowMass->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Low mass bg", "L");
		}
		if (invPeakShape && analysisSettings.showPeakShape) {
			invPeakShape->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShape->SetLineColor(analysisSettings.color);
			invPeakShape->SetMarkerColor(analysisSettings.color);
			TH1 *newHist = invPeakShape->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Sig MC", "L");
		}
		if (invPeakShapeLowMass && invPeakShape && analysisSettings.showBgLowMass && analysisSettings.showPeakShape) {
			invPeakShapeLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invPeakShapeLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invPeakShapeLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			TH1 *newHist = invPeakShapeLowMass->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Low mass bg MC", "L");
		}
		if (invPeakShapeEta && analysisSettings.showPeakShapeEta) {
			invPeakShapeEta->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEta->SetLineColor(analysisSettings.color);
			invPeakShapeEta->SetMarkerColor(analysisSettings.color);
			TH1 *newHist = invPeakShapeEta->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Sig #eta MC", "L");
		}
		if (invPeakShapeEtaLowMass && invPeakShapeEta && analysisSettings.showBgLowMass && analysisSettings.showPeakShapeEta) {
			invPeakShapeEtaLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEtaLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtaLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtaLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			TH1 *newHist = invPeakShapeEtaLowMass->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Low mass bg #eta MC", "L");
		}
		if (invPeakShapeEtabg && analysisSettings.showPeakShapeEtabg) {
			invPeakShapeEtabg->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEtabg->SetLineColor(analysisSettings.color);
			invPeakShapeEtabg->SetMarkerColor(analysisSettings.color);
			TH1 *newHist = invPeakShapeEtabg->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Bg #eta MC", "L");
		}
		if (invPeakShapeEtabgLowMass && invPeakShapeEtabg && analysisSettings.showBgLowMass && analysisSettings.showPeakShapeEtabg) {
			invPeakShapeEtabgLowMass->GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invPeakShapeEtabgLowMass->SetLineColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtabgLowMass->SetFillColor(analysisSettings.colorBgLowMass);
			invPeakShapeEtabgLowMass->SetMarkerColor(analysisSettings.colorBgLowMass);
			TH1 *newHist = invPeakShapeEtabgLowMass->DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Low mass bg bg #eta MC", "L");
		}
		if (analysisSettings.showResidual) {
			invRes.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invRes.SetLineColor(analysisSettings.color);
			TH1 *newHist = invRes.DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Residual", "L");
		}
		if (analysisSettings.showResidualBg) {
			invResBg.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invResBg.SetLineColor(analysisSettings.color);
			TH1 *newHist = invResBg.DrawCopy("HIST SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Residual bg", "L");
		}
		if (analysisSettings.showResidualPeakShapeBg) {
			invResPeakShapeBg.GetXaxis()->SetRangeUser(analysisSettings.binRangeLeft, analysisSettings.binRangeRight);
			invResPeakShapeBg.SetLineColor(analysisSettings.color);
			invResPeakShapeBg.SetMarkerColor(analysisSettings.color);
			TH1 *newHist = invResPeakShapeBg.DrawCopy("H SAME");
			if (legend && newHist) legend->AddEntry(newHist, "Residual MC bg", "L");
		}
		if (analysisSettings.showFuncPeak) {
			//TArrow *arFitLeftPeak = new TArrow(peakLeft, invSigBg.GetYaxis()->GetBinUpEdge(invSigBg.GetYaxis()->GetLast()), peakLeft, 0, 0.02, "-");
			//arFitLeftPeak->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetLineStyle(2); arFitLeftPeak->Draw();
			TArrow *arFitLeftPeakYield = new TArrow(peakLeftYield, invSigBg.GetMaximum(), peakLeftYield, 0, 0.05, "-");
			arFitLeftPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->Draw();
			//TArrow *arFitRightPeak = new TArrow(peakRight, invSigBg.GetMaximum(), peakRight, 0, 0.02, "-");
			//arFitRightPeak->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeak->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeak->SetLineStyle(2); arFitRightPeak->Draw();
			TArrow *arFitRightPeakYield = new TArrow(peakRightYield, invSigBg.GetMaximum(), peakRightYield, 0, 0.05, "-");
			arFitRightPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeakYield->Draw();
			if (legend) legend->AddEntry(arFitRightPeakYield, "Peak area", "L");
			funcBgPeak.SetLineColor(analysisSettings.color);
			TF1 *newFunc = funcBgPeak.DrawCopy("SAME");
			if (legend && newFunc) legend->AddEntry(newFunc, "Peak fit", "L");
		}
		if (analysisSettings.showFuncPeakEta) {
			//TArrow *arFitLeftPeak = new TArrow(peakLeftEta, invSigBg.GetMaximum(), peakLeftEta, 0, 0.02, "-");
			//arFitLeftPeak->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeak->SetLineStyle(2); arFitLeftPeak->Draw();
			TArrow *arFitLeftPeakYield = new TArrow(peakLeftYieldEta, invSigBg.GetMaximum(), peakLeftYieldEta, 0, 0.05, "-");
			arFitLeftPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitLeftPeakYield->Draw();
			//TArrow *arFitRightPeak = new TArrow(peakRightEta, invSigBg.GetMaximum(), peakRightEta, 0, 0.02, "-");
			//arFitRightPeak->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeak->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeak->SetLineStyle(2); arFitRightPeak->Draw();
			TArrow *arFitRightPeakYield = new TArrow(peakRightYieldEta, invSigBg.GetMaximum(), peakRightYieldEta, 0, 0.05, "-");
			arFitRightPeakYield->SetLineColor(analysisSettings.peakLinesColor); arFitRightPeakYield->SetFillColor(analysisSettings.peakLinesColor); arFitRightPeakYield->Draw();
			funcBgPeakEta.SetLineColor(analysisSettings.color);
			funcBgPeakEta.DrawCopy("SAME");
		}
		if (analysisSettings.showFuncBg) {
			TArrow *arFitRegionLeft = new TArrow(fitLeft, invSigBg.GetMaximum(), fitLeft, 0, 0.06, "-");
			arFitRegionLeft->SetLineColor(analysisSettings.fitLinesColor); arFitRegionLeft->SetFillColor(analysisSettings.fitLinesColor); arFitRegionLeft->SetLineStyle(3); arFitRegionLeft->Draw();
			TArrow *arFitRegionRight = new TArrow(fitRight, invSigBg.GetMaximum(), fitRight, 0, 0.06, "-");
			arFitRegionRight->SetLineColor(analysisSettings.fitLinesColor); arFitRegionRight->SetFillColor(analysisSettings.fitLinesColor); arFitRegionRight->SetLineStyle(3); arFitRegionRight->Draw();
			funcPureBg.SetLineColor(analysisSettings.bgColor);
			funcPureBg.SetFillColor(analysisSettings.fillColor);
			funcPureBg.SetFillStyle(analysisSettings.bgFillStyle);
			TF1 *newFunc = funcPureBg.DrawCopy("SAME");
			if (legend && newFunc) legend->AddEntry(newFunc, "Pure bg func", "L");
		}
		if (analysisSettings.showTruePeakPos) {
			TArrow *arFitMinvPt = new TArrow(analysisSettings.peakPositionTrue, invSigBg.GetMaximum(), analysisSettings.peakPositionTrue, 0, 0.03, "-");
			arFitMinvPt->SetLineColor(analysisSettings.color); arFitMinvPt->SetFillColor(analysisSettings.fillColor); arFitMinvPt->Draw();
			TArrow *arFitMinvPtEta = new TArrow(analysisSettings.peakPositionTrueEta, invSigBg.GetMaximum(), analysisSettings.peakPositionTrueEta, 0, 0.03, "-");
			arFitMinvPtEta->SetLineColor(analysisSettings.color); arFitMinvPtEta->SetFillColor(analysisSettings.fillColor); arFitMinvPtEta->Draw();
			if (legend) legend->AddEntry(arFitMinvPtEta, "True mass", "L");
		}
		if (analysisSettings.showText) {
			TString strDataset = analysisSettings.title;
			TString strPt = Form("%g < #font[12]{p}_{#font[12]{T}} < %g GeV/#font[12]{c}", binParameters.min, binParameters.max);
			TString strPi0 = Form("#pi^{0} = %.1f #pm %.1f", area, areaErr);
			TString strEta = Form("#eta^{^{  }}  = %.1f #pm %.1f", areaEta, areaEtaErr);
			TString strYield = "#splitline{" + strPi0 + "}{" + strEta + "}";
			TString strPtYield = "#splitline{" + strPt + "}{" + strYield + "}";
			TString strDatasetPtYield = "#splitline{" + strDataset + "}{" + strPtYield + "}";
			TString str = strDatasetPtYield;
		        TLatex *text = new TLatex((analysisSettings.binRangeLeft * 0.5) + (analysisSettings.binRangeRight * 0.5), invSigBg.GetMaximum() * 0.65, str.Data());
                        if (!(legend && analysisSettings.showSigBg)) text->Draw();
	                text->SetTextSize(textSize);
		}
		if (legend && analysisSettings.showSigBg) {
		    legend->Draw();
		}
	}
	analysisResults.peakPosition.push_back(TBinStatistics(binParameters, mean, meanErr));
	analysisResults.peakWidth.push_back(TBinStatistics(binParameters, width, widthErr));
	analysisResults.peakPositionData.push_back(TBinStatistics(binParameters, meanData, meanDataErr));
	analysisResults.peakWidthData.push_back(TBinStatistics(binParameters, widthData, widthDataErr));
	analysisResults.rawYield.push_back(TBinStatistics(binParameters, area, areaErr));
	analysisResults.peakBackground.push_back(TBinStatistics(binParameters, areaBg, areaBgErr));
	analysisResults.peakSigBackground.push_back(TBinStatistics(binParameters, areaSigBgFull, areaSigBgFullErr));
	analysisResults.peakYieldToSigBackground.push_back(TBinStatistics(binParameters, yieldToAreaSigBgFull, yieldToAreaSigBgFullErr));
	analysisResults.peakYieldToSigBackgroundRange.push_back(TBinStatistics(binParameters, yieldToAreaSigBgFullRange, yieldToAreaSigBgFullRangeErr));
	analysisResults.peakYieldToBackground.push_back(TBinStatistics(binParameters, yieldToBgFull, yieldToBgFullErr));
	analysisResults.lowmassToYield.push_back(TBinStatistics(binParameters, lowmassToYield, lowmassToYieldErr));
	analysisResults.etabgToYield.push_back(TBinStatistics(binParameters, etabgToYield, etabgToYieldErr));
	analysisResults.mixJetmixDiff.push_back(TBinStatistics(binParameters, mixJetmixDiff, mixJetmixDiffErr));

	analysisResults.peakPositionEta.push_back(TBinStatistics(binParameters, meanEta, meanEtaErr));
	analysisResults.peakWidthEta.push_back(TBinStatistics(binParameters, widthEta, widthEtaErr));
	analysisResults.peakPositionDataEta.push_back(TBinStatistics(binParameters, meanDataEta, meanDataEtaErr));
	analysisResults.peakWidthDataEta.push_back(TBinStatistics(binParameters, widthDataEta, widthDataEtaErr));
	analysisResults.rawYieldEta.push_back(TBinStatistics(binParameters, areaEta, areaEtaErr));
	analysisResults.peakBackgroundEta.push_back(TBinStatistics(binParameters, areaBgEta, areaBgEtaErr));
	analysisResults.peakSigBackgroundEta.push_back(TBinStatistics(binParameters, areaSigBgEtaFull, areaSigBgEtaFullErr));
	analysisResults.peakYieldToSigBackgroundEta.push_back(TBinStatistics(binParameters, yieldToAreaSigBgEtaFull, yieldToAreaSigBgEtaFullErr));
	analysisResults.peakYieldToBackgroundEta.push_back(TBinStatistics(binParameters, yieldToBgEtaFull, yieldToBgEtaFullErr));
	analysisResults.mixJetmixDiffEta.push_back(TBinStatistics(binParameters, mixJetmixDiffEta, mixJetmixDiffEtaErr));

	analysisResults.bgPar6.push_back(TBinStatistics(binParameters, bgPar6, bgPar6Err));
	analysisResults.bgPar7.push_back(TBinStatistics(binParameters, bgPar7, bgPar7Err));
	analysisResults.fitChi2Ndf.push_back(TBinStatistics(binParameters, fitChi2, fitNdf));
	analysisResults.bgJetToTotal.push_back(TBinStatistics(binParameters, bgJetToTotalRatio, bgJetToTotalRatioErr));
	analysisResults.areaFractionBg.push_back(TBinStatistics(binParameters, areaFractionBg, areaFractionBgErr));
	analysisResults.bgToSigbg.push_back(TBinStatistics(binParameters, bgToSigbg, bgToSigbgErr));

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
	analysisResults.scaleBgBack.push_back(TBinStatistics(binParameters, scaleBgBack, 0));
	analysisResults.scaleBgLowMass.push_back(TBinStatistics(binParameters, scaleBgLowMass, 0));
	analysisResults.scaleBgLowMass2.push_back(TBinStatistics(binParameters, scaleBgLowMass2, 0));
	analysisResults.scalePeakShape.push_back(TBinStatistics(binParameters, scalePeakShape, 0));
	analysisResults.scalePeakShapeEta.push_back(TBinStatistics(binParameters, scalePeakShapeEta, 0));

	analysisResults.jetBgFractionForFit.push_back(TBinStatistics(binParameters, jetBgFraction, 0));

	analysisResults.fitParameter6.push_back(TBinStatistics(binParameters, fitParameter6, 0));
	analysisResults.fitParameter7.push_back(TBinStatistics(binParameters, fitParameter7, 0));

	analysisResults.peakAreaShareBorder.push_back(TBinStatistics(binParameters, peakAreaShareBorder, 0));

	if (invBg) delete invBg;
	if (invBgRandom) delete invBgRandom;
	if (invBgBack) delete invBgBack;
	if (invBgLowMass) delete invBgLowMass;
	if (invPeakShape) delete invPeakShape;
	if (invPeakShapeEta) delete invPeakShapeEta;
	if (invPeakShapeEtabg) delete invPeakShapeEtabg;
	if (invPeakShapeLowMass) delete invPeakShapeLowMass;
	if (invPeakShapeEtaLowMass) delete invPeakShapeEtaLowMass;
	if (invPeakShapeEtabgLowMass) delete invPeakShapeEtabgLowMass;

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
	analysisResults.peakYieldToSigBackground.clear();
	analysisResults.peakYieldToSigBackgroundRange.clear();
	analysisResults.peakYieldToBackground.clear();
	analysisResults.lowmassToYield.clear();
	analysisResults.etabgToYield.clear();
	analysisResults.mixJetmixDiff.clear();

	analysisResults.rawYieldEta.clear();
	analysisResults.peakPositionEta.clear();
	analysisResults.peakWidthEta.clear();
	analysisResults.peakPositionDataEta.clear();
	analysisResults.peakWidthDataEta.clear();
	analysisResults.peakBackgroundEta.clear();
	analysisResults.peakSigBackgroundEta.clear();
	analysisResults.peakYieldToSigBackgroundEta.clear();
	analysisResults.peakYieldToBackgroundEta.clear();
	analysisResults.mixJetmixDiffEta.clear();

	analysisResults.bgPar6.clear();
	analysisResults.bgPar7.clear();
	analysisResults.fitChi2Ndf.clear();
	analysisResults.bgJetToTotal.clear();
	analysisResults.areaFractionBg.clear();
	analysisResults.bgToSigbg.clear();

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
	analysisResults.scaleBgBack.clear();
	analysisResults.scaleBgLowMass.clear();
	analysisResults.scaleBgLowMass2.clear();
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
                const distribution_type *distrBgBack = 0;
                const distribution_type *distrBgLowMass = 0;
                const distribution_type *distrBgLowMassNbar = 0;
                const distribution_type *distrMult = 0;
                const distribution_type *distrPeakShape = 0;
                const distribution_type *distrPeakShapeLowMass = 0;
                const distribution_type *distrPeakShapeLowMassNbar = 0;
                const distribution_type *distrPeakShapeEta = 0;
                const distribution_type *distrPeakShapeEtaLowMass = 0;
                const distribution_type *distrPeakShapeEtaLowMassNbar = 0;
                const distribution_type *distrPeakShapeEtabg = 0;
                const distribution_type *distrPeakShapeEtabgLowMass = 0;
                const distribution_type *distrPeakShapeEtabgLowMassNbar = 0;
                for (distribution_list_type::const_iterator iterBg = analysisResults.invBglist.begin();iterBg != analysisResults.invBglist.end();++iterBg) {
                                const distribution_type &distrBgRef = *iterBg;
                                if (distr == distrBgRef) distrBg = &distrBgRef;
                }
                for (distribution_list_type::const_iterator iterBgRandom = analysisResults.invBgRandomlist.begin();iterBgRandom != analysisResults.invBgRandomlist.end();++iterBgRandom) {
                                const distribution_type &distrBgRandomRef = *iterBgRandom;
                                if (distr == distrBgRandomRef) distrBgRandom = &distrBgRandomRef;
                }
                for (distribution_list_type::const_iterator iterBgBack = analysisResults.invBgBacklist.begin();iterBgBack != analysisResults.invBgBacklist.end();++iterBgBack) {
                                const distribution_type &distrBgBackRef = *iterBgBack;
                                if (distr == distrBgBackRef) distrBgBack = &distrBgBackRef;
                }
                for (distribution_list_type::const_iterator iterBgLowMass = analysisResults.inv1gammalist.begin();iterBgLowMass != analysisResults.inv1gammalist.end();++iterBgLowMass) {
                                const distribution_type &distrBgLowMassRef = *iterBgLowMass;
                                if (distr == distrBgLowMassRef) distrBgLowMass = &distrBgLowMassRef;
                }
                for (distribution_list_type::const_iterator iterBgLowMassNbar = analysisResults.invnbarlist.begin();iterBgLowMassNbar != analysisResults.invnbarlist.end();++iterBgLowMassNbar) {
                                const distribution_type &distrBgLowMassNbarRef = *iterBgLowMassNbar;
                                if (distr == distrBgLowMassNbarRef) distrBgLowMassNbar = &distrBgLowMassNbarRef;
                }
                for (distribution_list_type::const_iterator iterMult = analysisResults.multlist.begin();iterMult != analysisResults.multlist.end();++iterMult) {
                                const distribution_type &distrMultRef = *iterMult;
                                if (distr == distrMultRef) distrMult = &distrMultRef;
                }
                if (analysisResultsSim) for (distribution_list_type::const_iterator iterPeakShape = analysisResultsSim->invlist.begin();iterPeakShape != analysisResultsSim->invlist.end();++iterPeakShape) {
                                const distribution_type &distrPeakShapeRef = *iterPeakShape;
                                if (distr == distrPeakShapeRef) distrPeakShape = &distrPeakShapeRef;
                }
                if (analysisResultsSim) for (distribution_list_type::const_iterator iterPeakShapeLowMass = analysisResultsSim->inv1gammalist.begin();iterPeakShapeLowMass != analysisResultsSim->inv1gammalist.end();++iterPeakShapeLowMass) {
                                const distribution_type &distrPeakShapeLowMassRef = *iterPeakShapeLowMass;
                                if (distr == distrPeakShapeLowMassRef) distrPeakShapeLowMass = &distrPeakShapeLowMassRef;
                }
                if (analysisResultsSim) for (distribution_list_type::const_iterator iterPeakShapeLowMassNbar = analysisResultsSim->invnbarlist.begin();iterPeakShapeLowMassNbar != analysisResultsSim->invnbarlist.end();++iterPeakShapeLowMassNbar) {
                                const distribution_type &distrPeakShapeLowMassNbarRef = *iterPeakShapeLowMassNbar;
                                if (distr == distrPeakShapeLowMassNbarRef) distrPeakShapeLowMassNbar = &distrPeakShapeLowMassNbarRef;
                }
                if (analysisResultsSimEta) for (distribution_list_type::const_iterator iterPeakShapeEta = analysisResultsSimEta->invlist.begin();iterPeakShapeEta != analysisResultsSimEta->invlist.end();++iterPeakShapeEta) {
                                const distribution_type &distrPeakShapeEtaRef = *iterPeakShapeEta;
                                if (distr == distrPeakShapeEtaRef) distrPeakShapeEta = &distrPeakShapeEtaRef;
                }
                if (analysisResultsSimEta) for (distribution_list_type::const_iterator iterPeakShapeEtaLowMass = analysisResultsSimEta->inv1gammalist.begin();iterPeakShapeEtaLowMass != analysisResultsSimEta->inv1gammalist.end();++iterPeakShapeEtaLowMass) {
                                const distribution_type &distrPeakShapeEtaLowMassRef = *iterPeakShapeEtaLowMass;
                                if (distr == distrPeakShapeEtaLowMassRef) distrPeakShapeEtaLowMass = &distrPeakShapeEtaLowMassRef;
                }
                if (analysisResultsSimEta) for (distribution_list_type::const_iterator iterPeakShapeEtaLowMassNbar = analysisResultsSimEta->invnbarlist.begin();iterPeakShapeEtaLowMassNbar != analysisResultsSimEta->invnbarlist.end();++iterPeakShapeEtaLowMassNbar) {
                                const distribution_type &distrPeakShapeEtaLowMassNbarRef = *iterPeakShapeEtaLowMassNbar;
                                if (distr == distrPeakShapeEtaLowMassNbarRef) distrPeakShapeEtaLowMassNbar = &distrPeakShapeEtaLowMassNbarRef;
                }
                if (analysisResultsSimEtabg) for (distribution_list_type::const_iterator iterPeakShapeEtabg = analysisResultsSimEtabg->invlist.begin();iterPeakShapeEtabg != analysisResultsSimEtabg->invlist.end();++iterPeakShapeEtabg) {
                                const distribution_type &distrPeakShapeEtabgRef = *iterPeakShapeEtabg;
                                if (distr == distrPeakShapeEtabgRef) distrPeakShapeEtabg = &distrPeakShapeEtabgRef;
                }
                if (analysisResultsSimEtabg) for (distribution_list_type::const_iterator iterPeakShapeEtabgLowMass = analysisResultsSimEtabg->inv1gammalist.begin();iterPeakShapeEtabgLowMass != analysisResultsSimEtabg->inv1gammalist.end();++iterPeakShapeEtabgLowMass) {
                                const distribution_type &distrPeakShapeEtabgLowMassRef = *iterPeakShapeEtabgLowMass;
                                if (distr == distrPeakShapeEtabgLowMassRef) distrPeakShapeEtabgLowMass = &distrPeakShapeEtabgLowMassRef;
                }
                if (analysisResultsSimEtabg) for (distribution_list_type::const_iterator iterPeakShapeEtabgLowMassNbar = analysisResultsSimEtabg->invnbarlist.begin();iterPeakShapeEtabgLowMassNbar != analysisResultsSimEtabg->invnbarlist.end();++iterPeakShapeEtabgLowMassNbar) {
                                const distribution_type &distrPeakShapeEtabgLowMassNbarRef = *iterPeakShapeEtabgLowMassNbar;
                                if (distr == distrPeakShapeEtabgLowMassNbarRef) distrPeakShapeEtabgLowMassNbar = &distrPeakShapeEtabgLowMassNbarRef;
                }
                for (bin_stat_list_type::const_iterator iterLow = analysisResults.lowNormPointsPtCorr.begin();iterLow != analysisResults.lowNormPointsPtCorr.end();++iterLow) {
                                const bin_stat_type &binStatRef = *iterLow;
                                if (distr.getBinParameters() == binStatRef.getParameters()) analysisResults.lowNormPointsPt = binStatRef.getValue();
                }
                for (bin_stat_list_type::const_iterator iterLow = analysisResults.lowNormPointsNbarPtCorr.begin();iterLow != analysisResults.lowNormPointsNbarPtCorr.end();++iterLow) {
                                const bin_stat_type &binStatRef = *iterLow;
                                if (distr.getBinParameters() == binStatRef.getParameters()) analysisResults.lowNormPointsNbarPt = binStatRef.getValue();
                }
		analysisResults.bgJetToTotalFromCandidates = 0.0;
                for (bin_stat_list_type::const_iterator iterLow = analysisResults.bgJetToTotalFromCandidatesCorr.begin();iterLow != analysisResults.bgJetToTotalFromCandidatesCorr.end();++iterLow) {
                                const bin_stat_type &binStatRef = *iterLow;
                                if (distr.getBinParameters() == binStatRef.getParameters()) analysisResults.bgJetToTotalFromCandidates = binStatRef.getValue();
                }
		TFindPeakDistributions s = {&distr
				, distrBg
				, distrBgRandom
				, distrBgBack
				, distrBgLowMass
				, distrBgLowMassNbar
				, distrMult
				, distrPeakShape
				, distrPeakShapeLowMass
				, distrPeakShapeLowMassNbar
				, distrPeakShapeEta
				, distrPeakShapeEtaLowMass
				, distrPeakShapeEtaLowMassNbar
				, distrPeakShapeEtabg
				, distrPeakShapeEtabgLowMass
				, distrPeakShapeEtabgLowMassNbar};
                canvas = findPeak(s, analysisResults.invlist.size(), binToShow, canvas, analysisSettings, analysisResults);
                binToShow++;
        }
	analysisResults.totalRawYield = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.rawYield.begin();iter != analysisResults.rawYield.end();++iter) analysisResults.totalRawYield += (*iter).getValue();
	analysisResults.totalRawYieldEta = 0.0;
	for (bin_stat_list_type::const_iterator iter = analysisResults.rawYieldEta.begin();iter != analysisResults.rawYieldEta.end();++iter) analysisResults.totalRawYieldEta += (*iter).getValue();
        return canvas;
}

TCanvas *fitSpectraBins(const TAnalysisSettingsTrigger &analysisSettings
		      , TAnalysisRawResultsTrigger &analysisResults
		      , const TAnalysisRawResultsTrigger *analysisResultsSim
		      , const TAnalysisRawResultsTrigger *analysisResultsSimEta
		      , const TAnalysisRawResultsTrigger *analysisResultsSimEtabg
) {
	TAnalysisSettingsTrigger analysisSettingsNew = analysisSettings;

	Bool_t fitDone = false;
	if (analysisSettingsNew.correctCpv) {
	    analysisResults.cpvCorrection.clear();
	    analysisResults.cpvCorrectionMult.clear();
	    TAnalysisSettingsTrigger analysisSettingsTemp = analysisSettingsNew;
	    TAnalysisRawResultsTrigger analysisResultsTemp = analysisResults;
	    analysisResultsTemp.invlist = analysisResults.trackdistlist;
	    analysisResultsTemp.invBglist.clear();
	    analysisResultsTemp.invBgRandomlist.clear();
	    //analysisResultsTemp.invBgRandomlist = analysisResults.trackdist2list;
	    analysisResultsTemp.invBgBacklist.clear();
	    analysisResultsTemp.inv1gammalist.clear();
	    analysisResultsTemp.invnbarlist.clear();
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
	    analysisSettingsTemp.setErrorOnZeroBinsSigma = 0;

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

    	    analysisSettingsTemp.mixNormFixedRandom = 0.0;
    	    analysisSettingsTemp.mixNormMassRangeRandom = false;
    	    analysisSettingsTemp.mixNormLeftRandom = 0.25;
    	    analysisSettingsTemp.mixNormRightRandom = 0.7;
    	    analysisSettingsTemp.mixNormEntriesRandom = true;
 
            analysisSettingsTemp.subtractBg = false;
            analysisSettingsTemp.subtractBgRandom = false;
            analysisSettingsTemp.subtractBgBack = false;
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
	    analysisSettingsTemp.showBgBack = false;
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
	    fitDone = false;
	    fitDone |= fitParameter(&analysisResultsTemp.areaFractionBg, 0, analysisSettingsTemp.areaFractionBg, analysisResultsTemp.areaFractionBgParamFit, analysisSettingsNew.name + "_cpvFitFunc");
	    if (fitDone) {
		analysisSettingsTemp.showBins = false;
		findPeaks(analysisSettingsTemp, analysisResultsTemp, 0, 0, 0);
	    }

	    analysisResults.cpvCorrection = analysisResultsTemp.areaFractionBg;
	    analysisResults.cpvCorrectionParamFit = analysisResultsTemp.areaFractionBgParamFit;
	    analysisResults.cpvCorrectionParam6 = analysisResultsTemp.bgPar6;
	    analysisResults.cpvCorrectionParam7 = analysisResultsTemp.bgPar7;

	    if (analysisResultsTemp.cpvCorrectionMultFunc.func) {
		TF2 *funcLow = new TF2(*((TF2*)(analysisResultsTemp.cpvCorrectionMultFunc.func)));
		if (funcLow) for (Int_t p = 0;p < analysisResultsTemp.cpvCorrectionMultFunc.func->GetNpar();p++) 
		    funcLow->SetParameter(p, TMath::Max(analysisResultsTemp.cpvCorrectionMultFunc.func->GetParameter(p) - analysisResultsTemp.cpvCorrectionMultFunc.func->GetParError(p), analysisResultsTemp.cpvCorrectionMultFunc.func->GetParameter(p) / 1.05));
		TF2 *funcHigh = new TF2(*((TF2*)(analysisResultsTemp.cpvCorrectionMultFunc.func)));
		if (funcHigh) for (Int_t p = 0;p < analysisResultsTemp.cpvCorrectionMultFunc.func->GetNpar();p++) 
		    funcHigh->SetParameter(p, TMath::Min(analysisResultsTemp.cpvCorrectionMultFunc.func->GetParameter(p) + analysisResultsTemp.cpvCorrectionMultFunc.func->GetParError(p), analysisResultsTemp.cpvCorrectionMultFunc.func->GetParameter(p) * 1.05));
    	    for (distribution_list_type::const_iterator iterMult = analysisResults.multlist.begin();iterMult != analysisResults.multlist.end();++iterMult) {
                const distribution_type &distrMult = *iterMult;
		const TH1F *h = distrMult.getDistribution();
		if (h) {
		    //Float_t integral = 0;
		    //Float_t integral2 = 0;
		    Float_t peakRightCpv = getParameterValue(analysisSettingsNew.peakRightCpv, distrMult.getBinParameters());
		    Float_t cpv = 0;
		    Float_t cpvLow = 0;
		    Float_t cpvHigh = 0;
		    Float_t aTot = 0, bTot = 0;
		    Float_t aTotLow = 0, bTotLow = 0;
		    Float_t aTotHigh = 0, bTotHigh = 0;
		    Float_t aTotErr = 0, bTotErr = 0;
		    //Float_t aTotErrLow = 0, bTotErrLow = 0;
		    //Float_t aTotErrHigh = 0, bTotErrHigh = 0;
		    Float_t binNum = 0;
		    for (Int_t i = 1;i <= h->GetXaxis()->GetNbins();i++) {
			Float_t m = h->GetXaxis()->GetBinCenter(i);
			Float_t mLow = h->GetXaxis()->GetBinLowEdge(i);
			Float_t mHigh = h->GetXaxis()->GetBinUpEdge(i);
			if (m >= 1) {
			    Float_t a = ((TF2*)(analysisResultsTemp.cpvCorrectionMultFunc.func))->Integral(mLow, mHigh, 0.0, peakRightCpv);
			    Float_t b = ((TF2*)(analysisResultsTemp.cpvCorrectionMultFunc.func))->Integral(mLow, mHigh, 0.0, 1.0);
			    if ((b != 0) && ((h->GetBinContent(i) != 0) || (h->GetBinError(i) != 0))) {
				aTot += a * h->GetBinContent(i);
				bTot += b * h->GetBinContent(i);
				aTotErr += (a * h->GetBinError(i)) * (a * h->GetBinError(i));
				bTotErr += (b * h->GetBinError(i)) * (b * h->GetBinError(i));
				//Float_t cpv_i = (b - a)/b;
				//cpv += cpv_i * h->GetBinContent(i);
				//integral += h->GetBinContent(i);
				//integral2 += h->GetBinContent(i) * h->GetBinContent(i);
				binNum += 1.0;
				if (funcLow) {
				    Float_t aLow = funcLow->Integral(mLow, mHigh, 0.0, peakRightCpv);
				    Float_t bLow = funcLow->Integral(mLow, mHigh, 0.0, 1.0);
				    if (bLow != 0) {
					aTotLow += aLow * h->GetBinContent(i);
					bTotLow += bLow * h->GetBinContent(i);
					//aTotErrLow += aLow * h->GetBinContent(i) * h->GetBinContent(i);
					//bTotErrLow += bLow * h->GetBinContent(i) * h->GetBinContent(i);
					//Float_t cpv_iLow = (bLow - aLow)/bLow;
					//cpvLow += cpv_iLow * h->GetBinContent(i);
				    }
				}
				if (funcHigh) {
				    Float_t aHigh = funcHigh->Integral(mLow, mHigh, 0.0, peakRightCpv);
				    Float_t bHigh = funcHigh->Integral(mLow, mHigh, 0.0, 1.0);
				    if (bHigh != 0) {
					aTotHigh += aHigh * h->GetBinContent(i);
					bTotHigh += bHigh * h->GetBinContent(i);
					//aTotErrHigh += aHigh * h->GetBinContent(i) * h->GetBinContent(i);
					//bTotErrHigh += bHigh * h->GetBinContent(i) * h->GetBinContent(i);
					//Float_t cpv_iHigh = (bHigh - aHigh)/bHigh;
					//cpvHigh += cpv_iHigh * h->GetBinContent(i);
				    }
				}
				
			    }
			}
		    }
		    aTotErr = TMath::Sqrt(aTotErr);
		    bTotErr = TMath::Sqrt(bTotErr);
		    //integral2 = TMath::Sqrt(integral2);
		    //aTotErrLow = TMath::Sqrt(aTotErrLow);
		    //bTotErrLow = TMath::Sqrt(bTotErrLow);
		    //aTotErrHigh = TMath::Sqrt(aTotErrHigh);
		    //bTotErrHigh = TMath::Sqrt(bTotErrHigh);
		    //aTot = (integral != 0) ? aTot/integral : 0;
		    //bTot = (integral != 0) ? bTot/integral : 0;
			//aTotLow /= integral;
			//bTotLow /= integral;
			//aTotHigh /= integral;
			//bTotHigh /= integral;
		    //aTotErr = (integral2 != 0) ? aTotErr/integral2 : 0;
		    //bTotErr = (integral2 != 0) ? aTotErr/integral2 : 0;
			//aTotErrLow /= integral;
			//bTotErrLow /= integral;
			//aTotErrHigh /= integral;
			//bTotErrHigh /= integral;
		    aTotErr = (aTot != 0) ? (aTotErr / aTot) : 0;
		    bTotErr = (bTot != 0) ? (bTotErr / bTot) : 0;
		    Float_t e = TMath::Sqrt((aTotErr*aTotErr) + (bTotErr*bTotErr) - ((binNum > 0) ? (2.0*aTotErr*bTotErr/binNum) : 0)); // 1/binNum correlated
		    e = (binNum > 0) ? (e/TMath::Sqrt(binNum)) : e;
		    //Float_t e = TMath::Abs(aTotErr - bTotErr); // fully correlated
//cout << "binNum = " << binNum << ", aTotErr = " << aTotErr << ", bTotErr = " << bTotErr << ", e = " << e << endl;
		    cpv = (bTot != 0) ? ((bTot - aTot) / bTot) : 0.0; if (cpv == 0.0) cpv = 1.0;
		    cpvLow = (bTotLow != 0) ? ((bTotLow - aTotLow) / bTotLow) : 0.0; if (cpvLow == 0.0) cpvLow = 1.0;
		    cpvHigh = (bTotHigh != 0) ? ((bTotHigh - aTotHigh) / bTotHigh) : 0.0; if (cpvHigh == 0.0) cpvHigh = 1.0;
		    Float_t e1 = (cpvLow + cpvHigh != 0) ? (TMath::Abs(cpvLow - cpvHigh) / (cpvLow + cpvHigh)) : 0;
		    TBinStatistics binStat;
		    binStat.setParameters(distrMult.getBinParameters());
		    binStat.setValue(cpv);
		        //binStat.setError(TMath::Abs(cpvHigh - cpvLow)/2.0);
		    binStat.setError(cpv * TMath::Sqrt(e*e + e1*e1));
		    analysisResults.cpvCorrectionMult.push_back(binStat);
		}
            }
		//if (funcLow) delete funcLow;
		//if (funcHigh) delete funcHigh;
		analysisResults.cpvCorrection = analysisResults.cpvCorrectionMult;
	    }
	}

	{
	    TH1F h("h","h",1,0,1); h.Fill(0.5, 1);
	    fitParameter(0, &h, analysisSettingsNew.fitDistributionLeft, analysisSettingsNew.fitDistributionLeft, analysisSettingsNew.name + "_fitDistributionLeft");
	    fitParameter(0, &h, analysisSettingsNew.fitDistributionRight, analysisSettingsNew.fitDistributionRight, analysisSettingsNew.name + "_fitDistributionRight");
	    fitParameter(0, &h, analysisSettingsNew.fitDistributionLeft2, analysisSettingsNew.fitDistributionLeft2, analysisSettingsNew.name + "_fitDistributionLeft2");
	    fitParameter(0, &h, analysisSettingsNew.fitDistributionRight2, analysisSettingsNew.fitDistributionRight2, analysisSettingsNew.name + "_fitDistributionRight2");
	}

	{
	    bin_stat_list_type dummy; dummy.push_back(TBinStatistics());
	    fitParameter(&dummy, 0, analysisSettingsNew.lowmassbgNbar, analysisResults.lowmassbgNbar, analysisSettingsNew.name + "_lowmassbgNbar");
	}
	{
	    bin_stat_list_type dummy; dummy.push_back(TBinStatistics());
	    fitParameter(&dummy, 0, analysisSettingsNew.lowNormFixed, analysisResults.lowNormFixed, analysisSettingsNew.name + "_lowNormFixed");
	}

	TCanvas *canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	analysisResults.bgJetToTotalFirst = analysisResults.bgJetToTotal;
	fitDone = false;
	fitDone |= fitParameter(&analysisResults.bgJetToTotal, 0, analysisSettingsNew.jetBgFraction, analysisResults.jetBgFractionParamFit, analysisSettingsNew.name + "_jetBgFractionFitFunc");
	if (fitDone) {
	    if (analysisSettingsNew.showInvBinsFitsSeparately) {
		analysisSettingsNew.name += "_jetBgFractionFit";
		analysisSettingsNew.title += ", jet bg fraction fit";
	    }
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	fitDone = false;
	fitDone |= fitParameter(&analysisResults.peakPosition, 0, analysisSettingsNew.peakPosition, analysisResults.peakPositionParamFit, analysisSettingsNew.name + "_peakPositionFitFunc");
	fitDone |= fitParameter(&analysisResults.peakWidth, 0, analysisSettingsNew.peakWidth, analysisResults.peakWidthParamFit, analysisSettingsNew.name + "_peakWidthFitFunc");
	fitDone |= fitParameter(&analysisResults.peakPositionEta, 0, analysisSettingsNew.peakPositionEta, analysisResults.peakPositionEtaParamFit, analysisSettingsNew.name + "_peakPositionEtaFitFunc");
	fitDone |= fitParameter(&analysisResults.peakWidthEta, 0, analysisSettingsNew.peakWidthEta, analysisResults.peakWidthEtaParamFit, analysisSettingsNew.name + "_peakWidthEtaFitFunc");
	if (fitDone) {
	    if (analysisSettingsNew.showInvBinsFitsSeparately) {
		analysisSettingsNew.name += "_peakFit";
		analysisSettingsNew.title += ", peak range fit";
	    }
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	fitDone = false;
	fitDone |= fitParameter(&analysisResults.bgToSigbg, 0, analysisSettingsNew.bgToSigbg, analysisResults.bgToSigbgParamFit, analysisSettingsNew.name + "_bgToSigbgFitFunc");
	analysisResults.bgToSigbgFirst = analysisResults.bgToSigbg;
	if (fitDone) {
	    if (analysisSettingsNew.showInvBinsFitsSeparately) {
		analysisSettingsNew.name += "_bgToSigbgFit";
	        analysisSettingsNew.title += ", bg rel fit";
	    }
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	fitDone = false;
	fitDone |= fitParameter(&analysisResults.bgPar6, 0, analysisSettingsNew.fitParameter6, analysisResults.bgPar6ParamFit, analysisSettingsNew.name + "_bgPar6FitFunc");
	fitDone |= fitParameter(&analysisResults.bgPar7, 0, analysisSettingsNew.fitParameter7, analysisResults.bgPar7ParamFit, analysisSettingsNew.name + "_bgPar7FitFunc");
	if (fitDone) {
	    if (analysisSettingsNew.showInvBinsFitsSeparately) {
		analysisSettingsNew.name += "_bgParFit";
		analysisSettingsNew.title += ", bg parameters fit";
	    }
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	fitDone = false;
	fitDone |= fitParameter(&analysisResults.peakYieldToBackground, 0, analysisSettingsNew.signalToBackgroundRatio, analysisResults.peakSigBackgroundParamFit, analysisSettingsNew.name + "_sigBgFitFunc");
	fitDone |= fitParameter(&analysisResults.peakYieldToBackgroundEta, 0, analysisSettingsNew.signalToBackgroundRatioEta, analysisResults.peakSigBackgroundEtaParamFit, analysisSettingsNew.name + "_sigBgEtaFitFunc");
	if (fitDone) {
	    if (analysisSettingsNew.showInvBinsFitsSeparately) {
		analysisSettingsNew.name += "_sigBgFit";
		analysisSettingsNew.title += ", sig/sigbg fit";
	    }
	    canvas = findPeaks(analysisSettingsNew, analysisResults, analysisResultsSim, analysisResultsSimEta, analysisResultsSimEtabg);
	}
	return canvas;
}

void fitCpvMult(const TH2F *hist, const TFitSettings &funcSet, TFitSettings &func, TString name) {
    if (!hist) return;
    TH2F *h = new TH2F(*hist);
    for (Int_t i = 1;i <= h->GetXaxis()->GetNbins();i++) {
	Float_t integral = 0;
	for (Int_t j = 1;j <= h->GetYaxis()->GetNbins();j++) integral += h->GetBinContent(h->GetBin(i, j));
	if (integral > 0) for (Int_t j = 1;j <= h->GetYaxis()->GetNbins();j++) {
	    h->SetBinContent(h->GetBin(i, j), h->GetBinContent(h->GetBin(i, j)) / integral);
	    h->SetBinError(h->GetBin(i, j), h->GetBinError(h->GetBin(i, j)) / integral);
	}
    }
    func = funcSet;
    h->GetYaxis()->SetRangeUser(0.05, 0.6);
    if (!func.func) {
	TString nm = name; nm += "_cpvmultfunc";
	func.func = new TF2(nm.Data(), func.fitFuncStr.Data());
	if (func.func) {
	    func.func->FixParameter(0, 1.0);
	    func.func->FixParameter(1, 11.07);
	    func.func->FixParameter(2, 0.0852);
	    //func.func->FixParameter(1, 8.59);
	    //func.func->FixParameter(2, 0.114);
	}
    }
    if (func.func) {
	h->Fit(func.func, func.fitOption);
    }
}

Float_t getScaleFactorFromPt(TString name, TString title, const TH1F *histPt1, const TH1F *histPt2
			, Float_t lowPt, Float_t highPt, const distribution_list_type &invlist, bin_stat_list_type &lowNormPointsPtCorr, Bool_t show, Bool_t print
) {
	Float_t scale = 0;
	lowNormPointsPtCorr.clear();
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
				    if (scale != 0) hist2->Scale(scale);
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
			
			Float_t correctedPtLow = lowPt/* * (meanTruePt / meanDetPt)*/;
			Float_t correctedPtHigh = highPt/* * (meanTruePt / meanDetPt)*/;
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
	if (fPtShift) *fPtShift = new TF1(shiftFuncName, shiftFunc);
	if (fPtShift && *fPtShift) (*fPtShift)->SetRange(minX, maxX);
	if (hgrShift && fPtShift && *fPtShift) hgrShift->Fit(*fPtShift, "RQN");	
	//cout << shiftFuncName << " = " << shiftFunc << endl;
	//for (Int_t ipar = 0;ipar < (*fPtShift)->GetNpar();ipar++) {
	//	cout << shiftFuncName << ": [" << ipar << "] = " << (*fPtShift)->GetParameter(ipar) << endl;
	//}
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
			if (fPtShift && *fPtShift) (*fPtShift)->SetLineColor(kRed);
			if (fPtShift && *fPtShift) (*fPtShift)->SetLineWidth(3);
			if (fPtShift && *fPtShift) (*fPtShift)->DrawCopy("SAME");
			TLine *l = new TLine(minX, minX, maxX, maxX);
			l->SetLineColor(kBlue);
			l->Draw();
			Float_t par1 = (fPtShift && *fPtShift) ? (*fPtShift)->GetParameter(1) : 0.0;
			Float_t par1Err = (fPtShift && *fPtShift) ? (*fPtShift)->GetParError(1) : 0.0;
			TString textStr = Form("Slope = %1.2f +/- %1.2f", par1, par1Err);
			TLatex *text = new TLatex((minX + maxX) / 2.0, maxY, textStr.Data());
			text->Draw();
			canvas->cd(2);
			TString dummyName(shiftCanvasName); dummyName += "_dummy";
			TString dummyTitle(shiftCanvasTitle); dummyTitle = "Bin purity;p_{T reco}, GeV/c;Bin purity";
			TH2F *dummy = new TH2F(dummyName.Data(), dummyTitle.Data(), 10, minX, maxX, 10, 0, 1);
			dummy->SetStats(0);
			dummy->Draw();
			TGraphErrors *hgrPurity = shiftSize ? (new TGraphErrors(shiftSize, shiftX, purity, shiftXerr, purityErr)) : 0;
			if (hgrPurity) {
				hgrPurity->SetHistogram(dummy);
				hgrPurity->SetLineColor(kBlack);
				hgrPurity->SetLineWidth(3);
				hgrPurity->Draw("PL SAME");
			}
		}
	}
}

TF2 *ptBinPosFuncFromQCD = 0;
TF2 *ptBinPosCorrFuncFromQCD = 0;
typedef pair<Double_t, Double_t> ptBinPosFromQCD_point_type;
typedef map<ptBinPosFromQCD_point_type, Double_t> ptBinPosFromQCD_map_type;
ptBinPosFromQCD_map_type ptBinPosMapFromQCD;
ptBinPosFromQCD_map_type ptBinPosCorrMapFromQCD;
TF1 *pQCDfunc = 0;
Double_t getPtBinPosFromQCD(Double_t *x, Double_t * /*p*/) {
	Double_t newBinCenter = 0;
	if (!x) return newBinCenter;
	Float_t binCenter = x[0];
	newBinCenter = binCenter;
	Float_t binWidth = x[1];
	ptBinPosFromQCD_point_type point(binCenter, binWidth);
	ptBinPosFromQCD_map_type::const_iterator it = ptBinPosMapFromQCD.find(point);
	if (it != ptBinPosMapFromQCD.end()) {
	    newBinCenter = (*it).second;
	    //cout << "Cached: nominal center " << binCenter << ", width " << binWidth << ", true center " << newBinCenter << endl;
	} else {
	    if (pQCDfunc) {
		Float_t binLow = binCenter - (binWidth / 2.0);
		Float_t binHigh = binCenter + (binWidth / 2.0);
		Float_t integral = pQCDfunc->Integral(binLow, binHigh, (Double_t*)0, (Double_t)1e-6);
		Float_t mean = integral / (binHigh - binLow);
		Float_t meanX = pQCDfunc->GetX(mean, binLow, binHigh);
		newBinCenter = meanX;
		ptBinPosMapFromQCD[point] = newBinCenter;
		Float_t valNominalCenter = pQCDfunc->Eval(binCenter);
	        Float_t valTrueCenter = pQCDfunc->Eval(newBinCenter);
	        Float_t ratio = (valTrueCenter != 0) ? (valNominalCenter / valTrueCenter) : 1.0;
		ptBinPosCorrMapFromQCD[point] = ratio;
		//cout << "Calculated: nominal center " << binCenter << ", width " << binWidth << ", true center " << newBinCenter << ", correction " << ratio << endl;
	    }
	}
	return newBinCenter;
}
Double_t getPtBinPosCorrFromQCD(Double_t *x, Double_t * /*p*/) {
	Double_t newBinCenterCorr = 0;
	if (!x) return newBinCenterCorr;
	Float_t binCenter = x[0];
	Float_t binWidth = x[1];
	ptBinPosFromQCD_point_type point(binCenter, binWidth);
	ptBinPosFromQCD_map_type::const_iterator it = ptBinPosCorrMapFromQCD.find(point);
	if (it != ptBinPosCorrMapFromQCD.end()) {
	    newBinCenterCorr = (*it).second;
	    //cout << "Cached: nominal center " << binCenter << ", width " << binWidth << ", correction " << newBinCenterCorr << endl;
	}
	return newBinCenterCorr;
}
void clearPtBinPosFuncFromQCD() {
	ptBinPosMapFromQCD.clear();
	ptBinPosCorrMapFromQCD.clear();
	if (ptBinPosFuncFromQCD) delete ptBinPosFuncFromQCD;
	ptBinPosFuncFromQCD = 0;
	if (ptBinPosCorrFuncFromQCD) delete ptBinPosCorrFuncFromQCD;
	ptBinPosCorrFuncFromQCD = 0;
	//if (pQCDfunc) delete pQCDfunc;
	pQCDfunc = 0;
}
void createPtBinPosFuncFromQCD(TF1 *getWeightFunc) {
	clearPtBinPosFuncFromQCD();
	//if (!getWeightFunc) return;
	pQCDfunc = getWeightFunc;
	ptBinPosFuncFromQCD = new TF2("PtBinPosFromQCD", &getPtBinPosFromQCD, 0, 100, 0);
	ptBinPosCorrFuncFromQCD = new TF2("PtBinPosCorrFromQCD", &getPtBinPosCorrFromQCD, 0, 100, 0);
}

void calculatePtShiftCorr(const bin_stat_list_type &bins, bin_stat_list_type &corr, bin_stat_list_type &centerShift) {
    for (bin_stat_list_type::const_iterator iter = bins.begin();iter != bins.end();++iter) {
	const TBinStatistics &binStat = *iter;
	const TBinParameters &binParams = binStat.getParameters();
	Float_t lowPt = binParams.min;
	Float_t highPt = binParams.max;
	Float_t nominalCenter = (lowPt + highPt) / 2.0;
	Float_t newLowPt = lowPt;
	Float_t newHighPt = highPt;
	Float_t newNominalCenter = (newLowPt + newHighPt) / 2.0;
	Float_t trueCenter = ptBinPosFuncFromQCD ? ptBinPosFuncFromQCD->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : newNominalCenter;
	Float_t trueCenterCorr = ptBinPosCorrFuncFromQCD ? ptBinPosCorrFuncFromQCD->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : 0;
	    TBinParameters newBinParams(binParams);
	    newBinParams.min = newLowPt;
	    newBinParams.max = newHighPt;
	    newBinParams.trueCenter = trueCenter;
	    TBinStatistics binCorr;
	    binCorr.setParameters(newBinParams);
	    binCorr.setValue(trueCenterCorr);
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

void correctPtShift(bin_stat_list_type &bins, Bool_t horizontal, const distribution_list_type * /*simdistrlist*/) {
	for (bin_stat_list_type::iterator iter = bins.begin();iter != bins.end();++iter) {
		TBinStatistics &bin = *iter;
		const TBinParameters &binParams = bin.getParameters();
		Float_t lowPt = binParams.min;
		Float_t highPt = binParams.max;
		Float_t newLowPt = lowPt;
		Float_t newHighPt = highPt;
		Float_t newNominalCenter = (newLowPt + newHighPt) / 2.0;
		Float_t trueCenter = ptBinPosFuncFromQCD ? ptBinPosFuncFromQCD->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : newNominalCenter;
		Float_t trueCenterCorr = ptBinPosCorrFuncFromQCD ? ptBinPosCorrFuncFromQCD->Eval((Double_t)((newLowPt + newHighPt) / 2.0), (Double_t)(newHighPt - newLowPt)) : 0;
		TBinParameters newBinParams(binParams);
		newBinParams.min = newLowPt;
		newBinParams.max = newHighPt;
		newBinParams.trueCenter = trueCenter;
		if (horizontal) {
		    bin.setParameters(newBinParams);
		} else {
		    bin.setValue(bin.getValue() * trueCenterCorr);
		    bin.setError(bin.getError() * trueCenterCorr);
		}
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

void selectInv(const TCandidateDataProcessor *proc, Float_t binSize, Float_t pTLow, Float_t pTHigh, Int_t mult, list<TInvariantMassDistribution> &binlist, const TPointDataProcessor *pointDataProcessor, bool /*debug*/) {
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
	} else if (mult == 6) {
	    listToUse = proc ? &proc->pointTrackDist2Distributions : 0;
	}
	if (listToUse) for (Float_t pTLowBin = pTLow, pTHighBin = pTLow + binSize;pTHighBin <= pTHigh;pTLowBin += binSize, pTHighBin += binSize) {
//cout << "pTLowBin = " << pTLowBin << ", pTHighBin = " << pTHighBin << endl;
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
			    if (true) {
				for (UInt_t i = 0;i < preselectedInv.size();i++) {
				    selectedInv.push_back(preselectedInv[i]);
				}
			    } else {
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
			    }
            		    TBinParameters binPar;
            		    binPar.variable = (TBinVariable)0;
            		    binPar.min = pTLowBin;
            		    binPar.max = pTHighBin;
            		    TInvariantMassDistribution inv;
			    TString name;
			    TString title = Form("%g #leq p_{T} < %g GeV/c", pTLowBin, pTHighBin);
			    if (selectedInv.size() != 0) {
				for (list<const TInvariantMassDistribution*>::const_iterator iterInv = selectedInv.begin();iterInv != selectedInv.end();++iterInv) {
				    const TInvariantMassDistribution *iInv = *iterInv;
				    if (iterInv == selectedInv.begin()) {
					binPar.variable = iInv->getBinParameters().variable;
				    }
				    inv.add(*iInv, false);
//cout << "adding " << iInv->getBinParameters().min << " - " << iInv->getBinParameters().max << ", integral " << iInv->getDistribution()->Integral() << endl;
				    name += iInv->GetName();
				}
			    }
            		    inv.setBinParameters(binPar);
			    inv.SetNameTitle(name, title);
//cout << "name = " << name << ", title = " << title << ", integral " << inv.getDistribution()->Integral() << endl;
			    binlist.push_back(inv);
			}
		}
}

void selectInv(const TCandidateDataProcessor *proc, const TAnalysisSettingsTrigger &s, Int_t mult, list<TInvariantMassDistribution> &binlist, const TPointDataProcessor *pointDataProcessor, bool debug) {
    selectInv(proc, s.lowPtBinSize, s.lowPt, s.switchBinSizePt1, mult, binlist, pointDataProcessor, debug);
    selectInv(proc, s.midPtBinSize, s.switchBinSizePt1, s.switchBinSizePt2, mult, binlist, pointDataProcessor, debug);
    selectInv(proc, s.highPtBinSize, s.switchBinSizePt2, s.highPt, mult, binlist, pointDataProcessor, debug);
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
			    if (true) {
				for (UInt_t i = 0;i < preselectedInv.size();i++) {
				    selectedInv.push_back(preselectedInv[i]);
				}
			    } else {
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
			    }
            		    TBinParameters binPar;
            		    binPar.variable = (TBinVariable)0;
            		    binPar.min = pTLowBin;
            		    binPar.max = pTHighBin;
            		    TBinStatistics inv;
			    TString name;
			    TString title = Form("%g #leq p_{T} < %g GeV/c", pTLowBin, pTHighBin);
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

void selectInvSim(const TSimuDataProcessor *proc, const TAnalysisSettingsTrigger &s, bin_stat_list_type &binlist) {
    selectInvSim(proc, s.lowPtBinSize, s.lowPt, s.switchBinSizePt1, binlist);
    selectInvSim(proc, s.midPtBinSize, s.switchBinSizePt1, s.switchBinSizePt2, binlist);
    selectInvSim(proc, s.highPtBinSize, s.switchBinSizePt2, s.highPt, binlist);
}

void selectInvPoint(const TPointDataProcessor *proc, Float_t binSize, Float_t pTLow, Float_t pTHigh, bin_stat_list_type &binlist) {
    if (!proc) return;
    const TH1F *histPt = proc->getPt();
    if (!histPt) return;
    for (Float_t pTLowBin = pTLow, pTHighBin = pTLow + binSize;pTHighBin <= pTHigh;pTLowBin += binSize, pTHighBin += binSize) {
	Float_t temp;
	Float_t count = 0, countErr = 0;
	extractSignal(*histPt, funcZero, pTLowBin, pTHighBin, false, true, false, false, false, count, countErr, temp, temp, temp, temp, temp, temp, temp, temp);
        TBinParameters binPar;
        binPar.variable = (TBinVariable)0;
        binPar.min = pTLowBin;
        binPar.max = pTHighBin;
        TBinStatistics inv;
	TString name;
	TString title = Form("%g #leq p_{T} < %g GeV/c", pTLowBin, pTHighBin);
	name = title;
	inv.setValue(count);
	inv.setError(countErr);
        inv.setParameters(binPar);
	inv.SetNameTitle(name, title);
	binlist.push_back(inv);
    }
}

void selectInvPoint(const TPointDataProcessor *proc, const TAnalysisSettingsTrigger &s, bin_stat_list_type &binlist) {
    selectInvPoint(proc, s.lowPtBinSize, s.lowPt, s.switchBinSizePt1, binlist);
    selectInvPoint(proc, s.midPtBinSize, s.switchBinSizePt1, s.switchBinSizePt2, binlist);
    selectInvPoint(proc, s.highPtBinSize, s.switchBinSizePt2, s.highPt, binlist);
}

void selectInvMCGamma(const TMCGammaDataProcessor *proc, Float_t binSize, Float_t pTLow, Float_t pTHigh, bin_stat_list_type &binlist) {
    if (!proc) return;
    const TH1F *histPt = proc->getPt();
    if (!histPt) return;
    for (Float_t pTLowBin = pTLow, pTHighBin = pTLow + binSize;pTHighBin <= pTHigh;pTLowBin += binSize, pTHighBin += binSize) {
	Float_t temp;
	Float_t count = 0, countErr = 0;
	extractSignal(*histPt, funcZero, pTLowBin, pTHighBin, false, true, false, false, false, count, countErr, temp, temp, temp, temp, temp, temp, temp, temp);
        TBinParameters binPar;
        binPar.variable = (TBinVariable)0;
        binPar.min = pTLowBin;
        binPar.max = pTHighBin;
        TBinStatistics inv;
	TString name;
	TString title = Form("%g #leq p_{T} < %g GeV/c", pTLowBin, pTHighBin);
	name = title;
	inv.setValue(count);
	inv.setError(countErr);
        inv.setParameters(binPar);
	inv.SetNameTitle(name, title);
	binlist.push_back(inv);
    }
}

void selectInvMCGamma(const TMCGammaDataProcessor *proc, const TAnalysisSettingsTrigger &s, bin_stat_list_type &binlist) {
    selectInvMCGamma(proc, s.lowPtBinSize, s.lowPt, s.switchBinSizePt1, binlist);
    selectInvMCGamma(proc, s.midPtBinSize, s.switchBinSizePt1, s.switchBinSizePt2, binlist);
    selectInvMCGamma(proc, s.highPtBinSize, s.switchBinSizePt2, s.highPt, binlist);
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
	, Float_t &PSHT1_err, Float_t &PSHT2_err
	, Float_t &PSHT1Integral, Float_t &PSHT2Integral
	, Float_t &PSHT1Integral_err, Float_t &PSHT2Integral_err
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
	Float_t PSHT1r_err = 1.0, PSHT2r_err = 1.0, PSHT1i_err = 1.0, PSHT2i_err = 1.0;
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
		PSHT1r_err = funcPSHT1r->GetParError(0) / (1.0 + corrMBHT1);
		if (print) cout << "Prescales from points pT spectrum, " << title << ": PSHT1 = " << PSHT1r << " +/- " << PSHT1r_err << endl;
		PSHT1 = PSHT1r;
		PSHT1_err = PSHT1r_err;
	}
	if (histHT2r && funcPSHT2r) {
		funcPSHT2r->SetRange(psPointsHT2Low, psPointsHT2High);
		histHT2r->Fit(funcPSHT2r, "QRN");
		PSHT2r = PSHT1r * funcPSHT2r->GetParameter(0) / (100.0 * (1.0 + corrHT1HT2));
		PSHT2r_err = PSHT1r * funcPSHT2r->GetParError(0) / (100.0 * (1.0 + corrHT1HT2));
		PSHT2r_err = TMath::Sqrt(((PSHT1r != 0) ? (PSHT2r*PSHT1r_err/PSHT1r)*(PSHT2r*PSHT1r_err/PSHT1r) : 0) + PSHT2r_err*PSHT2r_err);
		if (print) cout << "Prescales from points pT spectrum, " << title << ": PSHT2 = " << PSHT2r << " +/- " << PSHT2r_err << endl;
		PSHT2 = PSHT2r;
		PSHT2_err = PSHT2r_err;
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
			PSHT1i_err = funcPSHT1i->GetParError(0) / (1.0 + corrMBHT1);
			if (print) cout << "Prescales from points pT integral, " << title << ": PSHT1 = " << PSHT1i << " +/- " << PSHT1i_err << endl;
			PSHT1Integral = PSHT1i;
			PSHT1Integral_err = PSHT1i_err;
		}
		if (funcPSHT2i) {
			funcPSHT2i->SetRange(psPointsHT2Lowi, psPointsHT2Highi);
			if (histHT2i) histHT2i->Fit(funcPSHT2i, "QRN");
			PSHT2i = PSHT1i * funcPSHT2i->GetParameter(0) / (100.0 * (1.0 + corrHT1HT2));
			PSHT2i_err = PSHT1i * funcPSHT2i->GetParError(0) / (100.0 * (1.0 + corrHT1HT2));
			PSHT2i_err = TMath::Sqrt(((PSHT1i != 0) ? (PSHT2i*PSHT1i_err/PSHT1i)*(PSHT2i*PSHT1i_err/PSHT1i) : 0) + PSHT2i_err*PSHT2i_err);
			if (print) cout << "Prescales from points pT integral, " << title << ": PSHT2 = " << PSHT2i << " +/- " << PSHT2i_err << endl;
			PSHT2Integral = PSHT2i;
			PSHT2Integral_err = PSHT2i_err;
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
			textStr = Form("#splitline{PS_{HT1} = %5.1f #pm %5.1f}{PS_{HT2} = %5.1f #pm %5.1f}", PSHT1r, PSHT1r_err, PSHT2r, PSHT2r_err);
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
			textStr = Form("#splitline{PS_{HT1} = %5.1f #pm %5.1f}{PS_{HT2} = %5.1f #pm %5.1f}", PSHT1i, PSHT1i_err, PSHT2i, PSHT2i_err);
			text = new TLatex(x, y, textStr.Data());
			text->Draw();
		}
	}
}

void calculatePSFromSim(const Char_t *name, const Char_t *title, const TDataProcessorPool *poolMB, const TDataProcessorPool *poolHT1
	, Float_t corrMBHT1, Float_t corrMBHT2, Float_t corrHT1HT2
	, Bool_t showReal, Bool_t print
	, Float_t &PSHT1, Float_t &PSHT2
	, Float_t &PSHT1_err, Float_t &PSHT2_err
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
		Float_t evNumPSMB = eventDataProcessorPSMB->binStatistics.size() ? ((*(eventDataProcessorPSMB->binStatistics.begin())).getValue()) : (eventDataProcessorPSMB->getZused() ? eventDataProcessorPSMB->getZused()->Integral() : eventDataProcessorPSMB->numPassedAllCuts);
		Float_t evNumPSMB_err = 0;//eventDataProcessorPSMB->binStatistics.size() ? ((*(eventDataProcessorPSMB->binStatistics.begin())).getError()) : 0;
		Float_t evNumPSMBHT1 = eventDataProcessorPSMBHT1->binStatistics.size() ? ((*(eventDataProcessorPSMBHT1->binStatistics.begin())).getValue()) : (eventDataProcessorPSMBHT1->getZused() ? eventDataProcessorPSMBHT1->getZused()->Integral() : eventDataProcessorPSMBHT1->numPassedAllCuts);
		Float_t evNumPSMBHT1_err = eventDataProcessorPSMBHT1->binStatistics.size() ? ((*(eventDataProcessorPSMBHT1->binStatistics.begin())).getError()) : 0;
		if (evNumPSMBHT1 != 0) PSHT1 = evNumPSMB / (evNumPSMBHT1 * (1.0 + corrMBHT1));
		PSHT1_err = PSHT1 * TMath::Sqrt(((evNumPSMBHT1 != 0) ? (evNumPSMBHT1_err/evNumPSMBHT1)*(evNumPSMBHT1_err/evNumPSMBHT1) : 0) + ((evNumPSMB != 0) ? (evNumPSMB_err/evNumPSMB)*(evNumPSMB_err/evNumPSMB) : 0));
		if (print) {
			cout << "HighTower-1 enhancement calculation from simulation, " << title << ": PSHT1 = " << PSHT1 << " +/- " << PSHT1_err << endl;
		}
	}
	if (eventDataProcessorPSMB && eventDataProcessorPSMBHT2) {
		Float_t evNumPSMB = eventDataProcessorPSMB->binStatistics.size() ? ((*(eventDataProcessorPSMB->binStatistics.begin())).getValue()) : (eventDataProcessorPSMB->getZused() ? eventDataProcessorPSMB->getZused()->Integral() : eventDataProcessorPSMB->numPassedAllCuts);
		Float_t evNumPSMB_err = 0;//eventDataProcessorPSMB->binStatistics.size() ? ((*(eventDataProcessorPSMB->binStatistics.begin())).getError()) : 0;
		Float_t evNumPSMBHT2 = eventDataProcessorPSMBHT2->binStatistics.size() ? ((*(eventDataProcessorPSMBHT2->binStatistics.begin())).getValue()) : (eventDataProcessorPSMBHT2->getZused() ? eventDataProcessorPSMBHT2->getZused()->Integral() : eventDataProcessorPSMBHT2->numPassedAllCuts);
		Float_t evNumPSMBHT2_err = eventDataProcessorPSMBHT2->binStatistics.size() ? ((*(eventDataProcessorPSMBHT2->binStatistics.begin())).getError()) : 0;
		if (evNumPSMBHT2 != 0) PSHT2 = evNumPSMB / (evNumPSMBHT2 * (1.0 + corrMBHT2));
		PSHT2_err = PSHT2 * TMath::Sqrt(((evNumPSMBHT2 != 0) ? (evNumPSMBHT2_err/evNumPSMBHT2)*(evNumPSMBHT2_err/evNumPSMBHT2) : 0) + ((evNumPSMB != 0) ? (evNumPSMB_err/evNumPSMB)*(evNumPSMB_err/evNumPSMB) : 0));
		if (print) {
			cout << "HighTower-2 enhancement calculation from simulation, " << title << ": PSHT2 = " << PSHT2 << " +/- " << PSHT2_err << endl;
		}
	}
	if (eventDataProcessorPSHT1 && eventDataProcessorPSHT1HT2) {
		Float_t PSHT1HT2 = PSHT2 / PSHT1;
		Float_t PSHT1HT2_err = PSHT1HT2 * TMath::Sqrt(((PSHT1 != 0) ? (PSHT1_err/PSHT1)*(PSHT1_err/PSHT1) : 0) + ((PSHT2 != 0) ? (PSHT2_err/PSHT2)*(PSHT2_err/PSHT2) : 0));
		if (print) {
			cout << "PSHT1HT2 = " << PSHT1HT2 << " +/- " << PSHT1HT2_err << endl;
		}
		Float_t evNumPSHT1 = eventDataProcessorPSHT1->binStatistics.size() ? ((*(eventDataProcessorPSHT1->binStatistics.begin())).getValue()) : (eventDataProcessorPSHT1->getZused() ? eventDataProcessorPSHT1->getZused()->Integral() : eventDataProcessorPSHT1->numPassedAllCuts);
		Float_t evNumPSHT1_err = 0;//eventDataProcessorPSHT1->binStatistics.size() ? ((*(eventDataProcessorPSHT1->binStatistics.begin())).getError()) : 0;
		Float_t evNumPSHT1HT2 = eventDataProcessorPSHT1HT2->binStatistics.size() ? ((*(eventDataProcessorPSHT1HT2->binStatistics.begin())).getValue()) : (eventDataProcessorPSHT1HT2->getZused() ? eventDataProcessorPSHT1HT2->getZused()->Integral() : eventDataProcessorPSHT1HT2->numPassedAllCuts);
		Float_t evNumPSHT1HT2_err = eventDataProcessorPSHT1HT2->binStatistics.size() ? ((*(eventDataProcessorPSHT1HT2->binStatistics.begin())).getError()) : 0;
		Float_t PSHT1HT2_2 = (evNumPSHT1HT2 != 0) ? (evNumPSHT1 / (evNumPSHT1HT2 * (1.0 + corrHT1HT2))) : 0;
		Float_t PSHT1HT2_2_err = PSHT1HT2_2 * TMath::Sqrt(((evNumPSHT1HT2 != 0) ? (evNumPSHT1HT2_err/evNumPSHT1HT2)*(evNumPSHT1HT2_err/evNumPSHT1HT2) : 0) + ((evNumPSHT1 != 0) ? (evNumPSHT1_err/evNumPSHT1)*(evNumPSHT1_err/evNumPSHT1) : 0));
		if (print) {
			cout << "HighTower-2 enhancement calculation - 2 from simulation, " << title << ": PSHT1HT2 = " << PSHT1HT2_2 << " +/-" << PSHT1HT2_2_err << endl;
		}
		if (PSHT1HT2_2 != 0) {
			Float_t PSHT2_2 = PSHT1 * PSHT1HT2_2;
			Float_t PSHT2_2_err = PSHT2_2 * TMath::Sqrt(((PSHT1 != 0) ? (PSHT1_err/PSHT1)*(PSHT1_err/PSHT1) : 0) + ((PSHT1HT2_2 != 0) ? (PSHT1HT2_2_err/PSHT1HT2_2)*(PSHT1HT2_2_err/PSHT1HT2_2) : 0));
			Float_t w1 = 0.5, w2 = 0.5;
			if ((PSHT2_err != 0) && (PSHT2_2_err != 0)) {
			    w1 = 1.0 / (PSHT2_err * PSHT2_err);
			    w2 = 1.0 / (PSHT2_2_err * PSHT2_2_err);
			}
			PSHT2 = ((PSHT2 * w1) + (PSHT2_2 * w2)) / (w1 + w2);
			PSHT2_err = TMath::Sqrt(1.0 / (w1 + w2));
			if (print) {
				cout << "\tPSHT1HT2 difference: " << (200.0 * TMath::Abs(PSHT1HT2 - PSHT1HT2_2) / (PSHT1HT2 + PSHT1HT2_2)) << "%" << endl;
				cout << "\tNew HighTower-2 factor: " << PSHT2 << " +/- " << PSHT2_err << endl;
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

void calculatePSFromDB(const Int_t runInd, const Int_t *runNums, const Int_t *evNumsMB, const Int_t *evNumsHT1, const Int_t *evNumsHT2, const Int_t *psMB, const Int_t *psHT1, const Int_t *psHT2, 
    Float_t &PS_HT1MB, Float_t &PS_HT2MB, 
    Float_t &PS_HT1MB_err, Float_t &PS_HT2MB_err, 
    Float_t &PS_HT1MB_alexst, Float_t &PS_HT2MB_alexst, Bool_t print) {
//print = true;
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
    Double_t evNumTotMBHT1true_err = 0;
    Double_t evNumTotMBHT2true_err = 0;
    Double_t evNumTotMBHT1 = 0;
    Double_t evNumTotMBHT2 = 0;
    Double_t evNumTotHT1MBtrue = 0;
    Double_t evNumTotHT1MBtrue_err = 0;
    Double_t evNumTotHT1HT2true = 0;
    Double_t evNumTotHT1HT2true_err = 0;
    Double_t evNumTotHT2MBtrue = 0;
    Double_t evNumTotHT2MBtrue_err = 0;
    Double_t evNumTotHT2HT1true = 0;
    Double_t evNumTotHT2HT1true_err = 0;
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
			evNumTotMBHT1true_err += (1.0 * runPsMB) * (1.0 * runPsMB);
			PSHT1MBSum += runPsHT1 * evNumHT1 / Double_t(runPsMB);
		    }
		    if (evNumHT2 > 0) {
			evNumTotMBHT2 += evNumMB;
			evNumTotMBHT2true += A;
			evNumTotMBHT2true_err += (1.0 * runPsMB) * (1.0 * runPsMB);
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
			evNumTotHT1MBtrue_err += (1.0 * runPsHT1) * (1.0 * runPsHT1);
			evNumTotHT1MBPrectrue += BHT1;
		    }
		    if (evNumHT2 > 0) {
			evNumTotHT1HT2true += BHT1;
			evNumTotHT1HT2true_err += (1.0 * runPsHT1) * (1.0 * runPsHT1);
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
			evNumTotHT2MBtrue_err += (1.0 * runPsHT2) * (1.0 * runPsHT2);
			evNumTotHT2MBPrectrue += BHT2;
		    }
		    if (evNumHT1 > 0) {
			evNumTotHT2HT1true += BHT2;
			evNumTotHT2HT1true_err += (1.0 * runPsHT2) * (1.0 * runPsHT2);
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
    evNumTotMBHT1true_err = TMath::Sqrt(evNumTotMBHT1true_err);
    evNumTotMBHT2true_err = TMath::Sqrt(evNumTotMBHT2true_err);
    evNumTotHT1MBtrue_err = TMath::Sqrt(evNumTotHT1MBtrue_err);
    evNumTotHT2MBtrue_err = TMath::Sqrt(evNumTotHT2MBtrue_err);
    evNumTotHT1HT2true_err = TMath::Sqrt(evNumTotHT1HT2true_err);
    evNumTotHT2HT1true_err = TMath::Sqrt(evNumTotHT2HT1true_err);
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
    Double_t HT1MB_err = HT1MB * TMath::Sqrt(((evNumTotMBHT1true != 0) ? (evNumTotMBHT1true_err/evNumTotMBHT1true)*(evNumTotMBHT1true_err/evNumTotMBHT1true) : 0) + ((evNumTotHT1MBtrue) ? (evNumTotHT1MBtrue_err/evNumTotHT1MBtrue)*(evNumTotHT1MBtrue_err/evNumTotHT1MBtrue) : 0));
    if (print) cout << "HT1/MB = " << HT1MB << " +/- " << HT1MB_err << endl;
    Double_t HT2HT1 = (evNumTotHT2HT1true != 0) ? (evNumTotHT1HT2true / evNumTotHT2HT1true) : 0;
    Double_t HT2HT1_err = HT2HT1 * TMath::Sqrt(((evNumTotHT1HT2true != 0) ? (evNumTotHT1HT2true_err/evNumTotHT1HT2true)*(evNumTotHT1HT2true_err/evNumTotHT1HT2true) : 0) + ((evNumTotHT2HT1true) ? (evNumTotHT2HT1true_err/evNumTotHT2HT1true)*(evNumTotHT2HT1true_err/evNumTotHT2HT1true) : 0));
    if (print) cout << "HT2/HT1 = " << HT2HT1 << " +/- " << HT2HT1 << endl;
    Double_t HT2MB = (evNumTotHT2MBtrue != 0) ? (evNumTotMBHT2true / evNumTotHT2MBtrue) : 0;
    Double_t HT2MB_err = HT2MB * TMath::Sqrt(((evNumTotMBHT2true != 0) ? (evNumTotMBHT2true_err/evNumTotMBHT2true)*(evNumTotMBHT2true_err/evNumTotMBHT2true) : 0) + ((evNumTotHT2MBtrue) ? (evNumTotHT2MBtrue_err/evNumTotHT2MBtrue)*(evNumTotHT2MBtrue_err/evNumTotHT2MBtrue) : 0));
    if (print) cout << "HT2/MB = " << HT2MB << " +/- " << HT2MB_err << endl;
    Double_t HT2MB_2 = HT1MB * HT2HT1;
    Double_t HT2MB_2_err = HT2MB_2 * TMath::Sqrt(((HT1MB != 0) ? (HT1MB_err/HT1MB)*(HT1MB_err/HT1MB) : 0) + ((HT2HT1) ? (HT2HT1_err/HT2HT1)*(HT2HT1_err/HT2HT1) : 0));
    if (print) cout << "HT1/MB * HT2/HT1 = " << HT2MB_2 << " +/- " << HT2MB_2_err << endl;
    Double_t w1 = 0.5, w2 = 0.5;
    if ((HT2MB_err != 0) && (HT2MB_2_err != 0)) {
	w1 = 1.0 / (HT2MB_err * HT2MB_err);
	w2 = 1.0 / (HT2MB_2_err * HT2MB_2_err);
    }
    HT2MB = ((HT2MB * w1) + (HT2MB_2 * w2)) / (w1 + w2);
    HT2MB_err = TMath::Sqrt(1.0 / (w1 + w2));
    if (print) cout << "HT2/MB (updated) = " << HT2MB << " +/- " << HT2MB_err << endl;
    if (print) cout << endl;
    PS_HT1MB = HT1MB;
    PS_HT1MB_err = HT1MB_err;
    PS_HT2MB = HT2MB;
    PS_HT2MB_err = HT2MB_err;
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

void calculatePSFromDB(const Char_t * /*name*/, const Char_t *title, const Char_t *psFilename
    , const TEventDataProcessor *eventDataProcessorMB, const TEventDataProcessor *eventDataProcessorHT1, const TEventDataProcessor *eventDataProcessorHT2
    , Bool_t /*show*/, Bool_t print, Float_t &PSHT1, Float_t &PSHT2, Float_t &PSHT1_err, Float_t &PSHT2_err, Float_t &PSHT1_alexst, Float_t &PSHT2_alexst) {
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
		calculatePSFromDB(runInd, runNums, evNumsMB, evNumsHT1, evNumsHT2, psMB, psHT1, psHT2, PSHT1, PSHT2, PSHT1_err, PSHT2_err, PSHT1_alexst, PSHT2_alexst, print && false);
		if (print) cout << "PSHT1 = " << PSHT1 << " +/- " << PSHT1_err << ", PSHT2 = " << PSHT2 << " +/- " << PSHT2_err << ", PSHT1_alexst = " << PSHT1_alexst << ", PSHT2_alexst = " << PSHT2_alexst << endl;
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

void showResultsOthers(TH1F *histSpectrum, TLegend *legendSpectrum, const data_points_list &dataPointsList) {
	if (!histSpectrum) return;
	for (data_points_list::const_iterator it = dataPointsList.begin();it != dataPointsList.end();++it) {
	    const TDataPoints &dataPoints = *it;
	    if (dataPoints.isFunction()) {
		if (dataPoints.func) {
		    
		    TDataPoints *newPoints = new TDataPoints(dataPoints);
		    Double_t xmin = 0, xmax = 0;
		    newPoints->func->GetRange(xmin, xmax);
		    TString newName;
		    {
			Int_t i = 0;
			Bool_t goodNewName = false;
			TCollection *listF = gROOT->GetListOfFunctions();
			do {
			    newName = newPoints->GetName(); newName += "_"; newName += newPoints->func->GetName(); newName += "_draw"; if (i) {newName += "_"; newName += i;}
			    if (listF && listF->Contains(newName.Data())) i++; else goodNewName = true;
			} while (!goodNewName);
		    }
		    //TString fstr = newPoints->func->GetName();
		    //if (newPoints->funcDenom) fstr += TString("/") + newPoints->funcDenom->GetName();
		    //TF1 *f = new TF1(newName.Data(), fstr.Data(), xmin, xmax);
		    TF1 *f = 0;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,21,0) 
		    f = new TF1(newName.Data(), ROOT::Math::ParamFunctor(newPoints, &TDataPoints::operator()), xmin, xmax, newPoints->func->GetNpar());
#endif
		    if (f) {
			f->SetParameters(newPoints->func->GetParameters());
			newPoints->TAttLine::Copy(*f);
			newPoints->TAttFill::Copy(*f);
			newPoints->TAttMarker::Copy(*f);
			TString drawOpt = newPoints->drawOption;
			drawOpt += " SAME";
			f->Draw((Option_t*)drawOpt.Data());
			if (legendSpectrum && (newPoints->legendTitle != "")) legendSpectrum->AddEntry(f, newPoints->legendTitle, newPoints->legendOption);
		    }
		    
		    /*
		    TString newName;
		    {
			Int_t i = 0;
			Bool_t goodNewName = false;
			TCollection *listF = gROOT->GetListOfFunctions();
			do {
			    newName = dataPoints.GetName(); newName += "_"; newName += dataPoints.func->GetName(); newName += "_draw"; if (i) {newName += "_"; newName += i;}
			    if (listF && listF->Contains(newName.Data())) i++; else goodNewName = true;
			} while (!goodNewName);
		    }
		    TString drawOpt = dataPoints.drawOption;
		    drawOpt += " SAME";
		    TF1 *f = dataPoints.func->DrawCopy((Option_t*)drawOpt.Data());
		    if (f) {
			f->SetName(newName.Data());
			dataPoints.TAttLine::Copy(*f);
			dataPoints.TAttFill::Copy(*f);
			dataPoints.TAttMarker::Copy(*f);
			if (legendSpectrum && (dataPoints.legendTitle != "")) legendSpectrum->AddEntry(f, dataPoints.legendTitle, dataPoints.legendOption);
		    }
		    */
		}
	    } else {
		showArrays(dataPoints.n, dataPoints.x, dataPoints.y, dataPoints.ex, dataPoints.ey, dataPoints.GetName(), dataPoints.GetTitle(), histSpectrum, &legendSpectrum, dataPoints, true);
	    }
	}
	if (legendSpectrum) legendSpectrum->Draw();
}

/*
void showPointsRcp(const Char_t *name, const Char_t *title
	, const TDataProcessorPool *poolMostcentral, const TDataProcessorPool *poolMostperipheral
	, Bool_t chargedPoints, Bool_t allPoints
	, Bool_t showSpectrum, Bool_t show, Bool_t print
	, Float_t binSizeMB, Float_t binSizeHT1, Float_t binSizeHT2
	, Float_t minPtMB, Float_t minPtHT1, Float_t minPtHT2
	, Float_t maxPtMB, Float_t maxPtHT1, Float_t maxPtHT2
) {
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
}
*/

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

/*
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
}
*/

void showResultsDAuCentrality(const Char_t *dAuCentralityName, const Char_t *dAuCentralityTitle
        , Bool_t showRcp, Bool_t showRcpEta, Bool_t dAuCentralityShow, Bool_t dAuCentralityPrint, Bool_t showDAuCentralityTriggersSeparately
	, Bool_t saveDataArrays, const Char_t *dataArraysFilenameFormat
        , const TAnalysisSettings &settingsDAuNoCentral
        , const TAnalysisSettings &/*settingsDAuAllCentral*/
        , const TAnalysisSettings &settingsDAuMostCentral
        , const TAnalysisSettings &settingsDAuMidCentral
        , const TAnalysisSettings &settingsDAuMostPeripheral
        , const TAnalysisResults &resultsDAuNoCentral
        , const TAnalysisResults &resultsDAuMostCentral
        , const TAnalysisResults &resultsDAuMidCentral
        , const TAnalysisResults &resultsDAuMostPeripheral
	, const data_points_list &dataPointsRcp
        , bin_stat_list_type &dAuRcpMB
        , bin_stat_list_type &dAuRcpHT1
        , bin_stat_list_type &dAuRcpHT2
        , bin_stat_list_type &dAuRcp
        , bin_stat_list_type &dAuRcpMBEta
        , bin_stat_list_type &dAuRcpHT1Eta
        , bin_stat_list_type &dAuRcpHT2Eta
        , bin_stat_list_type &dAuRcpEta
) {
	Float_t denomNbinary = resultsDAuMostCentral.NBinaryCollisions / resultsDAuMostPeripheral.NBinaryCollisions;
	divideBins(resultsDAuMostCentral.spectrumMB, resultsDAuMostPeripheral.spectrumMB, dAuRcpMB, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT1, resultsDAuMostPeripheral.spectrumHT1, dAuRcpHT1, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT2, resultsDAuMostPeripheral.spectrumHT2, dAuRcpHT2, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrum, resultsDAuMostPeripheral.spectrum, dAuRcp, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumMBEta, resultsDAuMostPeripheral.spectrumMBEta, dAuRcpMBEta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT1Eta, resultsDAuMostPeripheral.spectrumHT1Eta, dAuRcpHT1Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumHT2Eta, resultsDAuMostPeripheral.spectrumHT2Eta, dAuRcpHT2Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuMostCentral.spectrumEta, resultsDAuMostPeripheral.spectrumEta, dAuRcpEta, denomNbinary, 0.0, false, false);

	TLegend *legendCentral = 0;
	TString nameCentral = "spectrum_centralitybins"; nameCentral += dAuCentralityName;
	TString titleCentral = "Centrality dependence of p_{T} inv. yield"; titleCentral += ";p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}n}{d#etadp_{T}}";
	TH1F *histCentral = 0;
	if (settingsDAuNoCentral.showTriggersSeparately) {
	    histCentral = showLists(&resultsDAuNoCentral.spectrumMB, &resultsDAuNoCentral.spectrumHT1, &resultsDAuNoCentral.spectrumHT2
		, nameCentral.Data(), titleCentral.Data()
    		, histCentral, &legendCentral
		, settingsDAuNoCentral.settingsMB.drawOptions, settingsDAuNoCentral.settingsHT1.drawOptions, settingsDAuNoCentral.settingsHT2.drawOptions
		, dAuCentralityShow);
	} else {
	    histCentral = showList(&resultsDAuNoCentral.spectrum
		, nameCentral.Data(), titleCentral.Data()
    		, histCentral, &legendCentral
	        , settingsDAuNoCentral.drawOptions
	        , dAuCentralityShow);
	}
	if (settingsDAuMostCentral.showTriggersSeparately) {
	    histCentral = showLists(&resultsDAuMostCentral.spectrumMB, &resultsDAuMostCentral.spectrumHT1, &resultsDAuMostCentral.spectrumHT2
		, nameCentral.Data(), titleCentral.Data()
    		, histCentral, &legendCentral
		, settingsDAuMostCentral.settingsMB.drawOptions, settingsDAuMostCentral.settingsHT1.drawOptions, settingsDAuMostCentral.settingsHT2.drawOptions
		, dAuCentralityShow);
	} else {
	    histCentral = showList(&resultsDAuMostCentral.spectrum
		, nameCentral.Data(), titleCentral.Data()
    		, histCentral, &legendCentral
	        , settingsDAuMostCentral.drawOptions
	        , dAuCentralityShow);
	}
	if (settingsDAuMidCentral.showTriggersSeparately) {
		histCentral = showLists(&resultsDAuMidCentral.spectrumMB, &resultsDAuMidCentral.spectrumHT1, &resultsDAuMidCentral.spectrumHT2
		    , nameCentral.Data(), titleCentral.Data()
		    , histCentral, &legendCentral
		    , settingsDAuMidCentral.settingsMB.drawOptions, settingsDAuMidCentral.settingsHT1.drawOptions, settingsDAuMidCentral.settingsHT2.drawOptions
		    , dAuCentralityShow);
	} else {
		histCentral = showList(&resultsDAuMidCentral.spectrum
		    , nameCentral.Data(), titleCentral.Data()
		    , histCentral, &legendCentral
		    , settingsDAuMidCentral.drawOptions
		    , dAuCentralityShow);
	}
	if (settingsDAuMostPeripheral.showTriggersSeparately) {
	    histCentral = showLists(&resultsDAuMostPeripheral.spectrumMB, &resultsDAuMostPeripheral.spectrumHT1, &resultsDAuMostPeripheral.spectrumHT2
		, nameCentral.Data(), titleCentral.Data()
    		, histCentral, &legendCentral
		, settingsDAuMostPeripheral.settingsMB.drawOptions, settingsDAuMostPeripheral.settingsHT1.drawOptions, settingsDAuMostPeripheral.settingsHT2.drawOptions
		, dAuCentralityShow);
	} else {
	    histCentral = showList(&resultsDAuMostPeripheral.spectrum
		, nameCentral.Data(), titleCentral.Data()
    		, histCentral, &legendCentral
	        , settingsDAuMostPeripheral.drawOptions
	        , dAuCentralityShow);
	}
	showResultsOthers(histCentral, legendCentral, settingsDAuNoCentral.dataPointsInvYield);
	if (histCentral) {
		histCentral->GetXaxis()->SetLimits(0, settingsDAuMostCentral.settingsHT2.highPt);
		histCentral->SetMinimum(1e-11);
		histCentral->SetMaximum(1e+01);
		gPad->SetLogy();
	}

	TString titleRcpCentrality = Form(", %g-%g%%/%g-%g%%", (Float_t)(Int_t)resultsDAuMostCentral.PercentCentralityLow, (Float_t)(Int_t)resultsDAuMostCentral.PercentCentralityHigh, (Float_t)(Int_t)resultsDAuMostPeripheral.PercentCentralityLow, (Float_t)(Int_t)resultsDAuMostPeripheral.PercentCentralityHigh);
	TLegend *legendRcp = 0;
	TString nameRcp = dAuCentralityName;
	TString titleRcp = dAuCentralityTitle; titleRcp += titleRcpCentrality; titleRcp += ";p_{T}, GeV/c;R^{dAu}_{CP}";
	TH1F *histRcp = 0;
	if (showRcp) {
	    if (showDAuCentralityTriggersSeparately) {
		histRcp = showLists(&dAuRcpMB, &dAuRcpHT1, &dAuRcpHT2
		    , nameRcp.Data(), titleRcp.Data()
		    , histRcp, &legendRcp
		    , settingsDAuNoCentral.settingsMB.drawOptions, settingsDAuNoCentral.settingsHT1.drawOptions, settingsDAuNoCentral.settingsHT2.drawOptions
		    , dAuCentralityShow);
	    } else {
		TDrawOptions temp = settingsDAuNoCentral.drawOptionsEta;
		temp.legendTitle = "STAR R^{dAu}_{CP} #pi^{0}";
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
		    , settingsDAuNoCentral.settingsMB.drawOptionsEta, settingsDAuNoCentral.settingsHT1.drawOptionsEta, settingsDAuNoCentral.settingsHT2.drawOptionsEta
		    , dAuCentralityShow);
	    } else {
		TDrawOptions temp = settingsDAuNoCentral.drawOptionsEta;
		temp.legendTitle = "STAR R^{dAu}_{CP} #eta";
		histRcp = showList(&dAuRcpEta
		    , nameRcp.Data(), titleRcp.Data()
		    , histRcp, &legendRcp
		    , temp
		    , dAuCentralityShow);
	    }
	}
	if (histRcp) {
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
	    textStr = Form("#splitline{Glauber:}{#splitline{<N_{bin}>_{%g-%g%%} = %3.1f #pm %3.1f}{#splitline{<N_{bin}>_{%g-%g%%} = %3.1f #pm %3.1f}{%g%% error on ratio}}}"
		, (Float_t)(Int_t)resultsDAuMostCentral.PercentCentralityLow, (Float_t)(Int_t)resultsDAuMostCentral.PercentCentralityHigh, resultsDAuMostCentral.NBinaryCollisions, resultsDAuMostCentral.NBinaryCollisionsErr
		, (Float_t)(Int_t)resultsDAuMostPeripheral.PercentCentralityLow, (Float_t)(Int_t)resultsDAuMostPeripheral.PercentCentralityHigh, resultsDAuMostPeripheral.NBinaryCollisions, resultsDAuMostPeripheral.NBinaryCollisionsErr
		, 100.0 * err1);
	    TLatex *text = new TLatex(10.0, 0.8, textStr.Data());
	    text->SetTextSize(0.03);
	    text->Draw();
	}
	showResultsOthers(histRcp, legendRcp, dataPointsRcp);
	if (histRcp) {
	    histRcp->GetXaxis()->SetLimits(0, settingsDAuMostCentral.settingsHT2.highPt);
	    histRcp->SetMinimum(0);
	    histRcp->SetMaximum(2.0);
	}

	if (dAuCentralityPrint) {
		Float_t mostCentralPercent = (resultsDAuNoCentral.resultsMB.evNum > 0) ? (resultsDAuMostCentral.resultsMB.evNum / resultsDAuNoCentral.resultsMB.evNum) : 0.0;
		Float_t midCentralPercent = (resultsDAuNoCentral.resultsMB.evNum > 0) ? (resultsDAuMidCentral.resultsMB.evNum / resultsDAuNoCentral.resultsMB.evNum) : 0.0;
		Float_t mostPeripheralPercent = (resultsDAuNoCentral.resultsMB.evNum > 0) ? (resultsDAuMostPeripheral.resultsMB.evNum / resultsDAuNoCentral.resultsMB.evNum) : 0.0;
		Float_t mostCentralToMostPeripheral = (resultsDAuMostPeripheral.resultsMB.evNum > 0) ? (resultsDAuMostCentral.resultsMB.evNum / resultsDAuMostPeripheral.resultsMB.evNum) : 0.0;
		cout << "==============================================" << endl;
		cout << "Current centrality bins " << dAuCentralityTitle << ":" << endl;
		cout << "Most central: " << (mostCentralPercent * 100.0) << "%" << endl;
		cout << "Mid central: " << (midCentralPercent * 100.0) << "%" << endl;
		cout << "Most peripheral: " << (mostPeripheralPercent * 100.0) << "%" << endl;
		cout << "Sum: " << ((mostCentralPercent + midCentralPercent + mostPeripheralPercent) * 100) << "%" << endl;
		cout << "Most central / most peripheral: " << (mostCentralToMostPeripheral) << endl;
	}
	if (saveDataArrays) {
		TString name = dAuCentralityName;
		TString title = dAuCentralityTitle;
		TString filename = Form(dataArraysFilenameFormat, name.Data());
		cout << "Saving data arrays to " << filename << endl;
		ofstream ostr(filename.Data());
		ostr << "#ifndef DATA_ARRAYS_" << name << endl;
		ostr << "#define DATA_ARRAYS_" << name << endl;
		ostr << "// " << title << endl;
		TDatime startTime;
		ostr << "// Generated: " << startTime.AsSQLString() << endl;
		TString arrayName;
		arrayName = "_RcpMB_"; arrayName += name; outputCArrays(ostr, dAuRcpMB, "R_CP - MinBias", arrayName.Data());
		arrayName = "_RcpHT1_"; arrayName += name; outputCArrays(ostr, dAuRcpHT1, "R_CP - HighTower-1", arrayName.Data());
		arrayName = "_RcpHT2_"; arrayName += name; outputCArrays(ostr, dAuRcpHT2, "R_CP - HighTower-2", arrayName.Data());
		arrayName = "_Rcp_"; arrayName += name; outputCArrays(ostr, dAuRcp, "R_CP - combined triggers", arrayName.Data());
		arrayName = "_RcpMBEta_"; arrayName += name; outputCArrays(ostr, dAuRcpMBEta, "R_CP - eta - MinBias", arrayName.Data());
		arrayName = "_RcpHT1Eta_"; arrayName += name; outputCArrays(ostr, dAuRcpHT1Eta, "R_CP - eta - HighTower-1", arrayName.Data());
		arrayName = "_RcpHT2Eta_"; arrayName += name; outputCArrays(ostr, dAuRcpHT2Eta, "R_CP - eta - HighTower-2", arrayName.Data());
		arrayName = "_RcpEta_"; arrayName += name; outputCArrays(ostr, dAuRcpEta, "R_CP - eta - combined triggers", arrayName.Data());
		ostr << "#endif" << endl;
		ostr << "" << endl;
		ostr.close();
	}
}

void showResultsRDA(const Char_t *RDAName, const Char_t *RDATitle
        , Bool_t showRDA, Bool_t showRDAEta, Bool_t showRDAGamma, Bool_t RDAShow, Bool_t showRDATriggersSeparately
	, Bool_t saveDataArrays, const Char_t *dataArraysFilenameFormat
        , const TAnalysisSettings &settingsDAuNoCentral
        , const TAnalysisSettings &settingsPP
        , const TAnalysisResults &resultsDAuNoCentral
        , const TAnalysisResults &resultsPP
	, const data_points_list &dataPointsRda
        , bin_stat_list_type &RDAMB
        , bin_stat_list_type &RDAHT1
        , bin_stat_list_type &RDAHT2
        , bin_stat_list_type &RDA
        , bin_stat_list_type &RDAMBEta
        , bin_stat_list_type &RDAHT1Eta
        , bin_stat_list_type &RDAHT2Eta
        , bin_stat_list_type &RDAEta
        , bin_stat_list_type &RDAMBGamma
        , bin_stat_list_type &RDAHT1Gamma
        , bin_stat_list_type &RDAHT2Gamma
        , bin_stat_list_type &RDAGamma
) {
	Float_t denomNbinary = (resultsDAuNoCentral.NBinaryCollisions / resultsPP.NBinaryCollisions) / settingsPP.crossectionInel;
	divideBins(resultsDAuNoCentral.spectrumMB, resultsPP.crossectionMB, RDAMB, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrumHT1, resultsPP.crossectionHT1, RDAHT1, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrumHT2, resultsPP.crossectionHT2, RDAHT2, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrum, resultsPP.crossection, RDA, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrumMBEta, resultsPP.crossectionMBEta, RDAMBEta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrumHT1Eta, resultsPP.crossectionHT1Eta, RDAHT1Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrumHT2Eta, resultsPP.crossectionHT2Eta, RDAHT2Eta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.spectrumEta, resultsPP.crossectionEta, RDAEta, denomNbinary, 0.0, false, false);
	divideBins(resultsDAuNoCentral.RGammaMB, resultsPP.RGammaMB, RDAMBGamma, 1.0, 0.0, false, false);
	divideBins(resultsDAuNoCentral.RGammaHT1, resultsPP.RGammaHT1, RDAHT1Gamma, 1.0, 0.0, false, false);
	divideBins(resultsDAuNoCentral.RGammaHT2, resultsPP.RGammaHT2, RDAHT2Gamma, 1.0, 0.0, false, false);
	divideBins(resultsDAuNoCentral.RGamma, resultsPP.RGamma, RDAGamma, 1.0, 0.0, false, false);

	TLegend *legendRDA = 0;
	TString nameRDA = RDAName;
	TString titleRDA = RDATitle; titleRDA += ";p_{T}, GeV/c;R_{dA}";
	TH1F *histRDA = 0;
	if (showRDA) {
	    if (showRDATriggersSeparately) {
		histRDA = showLists(&RDAMB, &RDAHT1, &RDAHT2
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , settingsDAuNoCentral.settingsMB.drawOptions, settingsDAuNoCentral.settingsHT1.drawOptions, settingsDAuNoCentral.settingsHT2.drawOptions
		    , RDAShow);
	    } else {
		TDrawOptions temp = settingsDAuNoCentral.drawOptionsEta;
		temp.legendTitle = "STAR R_{dA} #pi^{0}";
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
		    , settingsDAuNoCentral.settingsMB.drawOptionsEta, settingsDAuNoCentral.settingsHT1.drawOptionsEta, settingsDAuNoCentral.settingsHT2.drawOptionsEta
		    , RDAShow);
	    } else {
		TDrawOptions temp = settingsDAuNoCentral.drawOptionsEta;
		temp.legendTitle = "STAR R_{dA} #eta";
		histRDA = showList(&RDAEta
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , temp
		    , RDAShow);
	    }
	}
	if (showRDAGamma) {
		histRDA = showLists(&RDAMBGamma, &RDAHT1Gamma, &RDAHT2Gamma
		    , nameRDA.Data(), titleRDA.Data()
		    , histRDA, &legendRDA
		    , settingsDAuNoCentral.settingsMB.drawOptions, settingsDAuNoCentral.settingsHT1.drawOptions, settingsDAuNoCentral.settingsHT2.drawOptions
		    , RDAShow);
	}
	showResultsOthers(histRDA, legendRDA, dataPointsRda);
	if (histRDA) {
	    histRDA->GetXaxis()->SetLimits(0, settingsDAuNoCentral.settingsHT2.highPt);
	    histRDA->SetMinimum(0);
	    histRDA->SetMaximum(2.0);
	}
	if (saveDataArrays) {
		TString name = RDAName;
		TString title = RDATitle;
		TString filename = Form(dataArraysFilenameFormat, name.Data());
		cout << "Saving data arrays to " << filename << endl;
		ofstream ostr(filename.Data());
		ostr << "#ifndef DATA_ARRAYS_" << name << endl;
		ostr << "#define DATA_ARRAYS_" << name << endl;
		ostr << "// " << title << endl;
		TDatime startTime;
		ostr << "// Generated: " << startTime.AsSQLString() << endl;
		TString arrayName;
		arrayName = "_RDAMB_"; arrayName += name; outputCArrays(ostr, RDAMB, "R_dA - MinBias", arrayName.Data());
		arrayName = "_RDAHT1_"; arrayName += name; outputCArrays(ostr, RDAHT1, "R_dA - HighTower-1", arrayName.Data());
		arrayName = "_RDAHT2_"; arrayName += name; outputCArrays(ostr, RDAHT2, "R_dA - HighTower-2", arrayName.Data());
		arrayName = "_RDA_"; arrayName += name; outputCArrays(ostr, RDA, "R_dA - combined triggers", arrayName.Data());
		arrayName = "_RDAMBEta_"; arrayName += name; outputCArrays(ostr, RDAMBEta, "R_dA - eta MinBias", arrayName.Data());
		arrayName = "_RDAHT1Eta_"; arrayName += name; outputCArrays(ostr, RDAHT1Eta, "R_dA - eta - HighTower-1", arrayName.Data());
		arrayName = "_RDAHT2Eta_"; arrayName += name; outputCArrays(ostr, RDAHT2Eta, "R_dA - eta - HighTower-2", arrayName.Data());
		arrayName = "_RDAEta_"; arrayName += name; outputCArrays(ostr, RDAEta, "R_dA - eta - combined triggers", arrayName.Data());
		arrayName = "_RDAMBGamma_"; arrayName += name; outputCArrays(ostr, RDAMBGamma, "R_dA - gamma - MinBias", arrayName.Data());
		arrayName = "_RDAHT1Gamma_"; arrayName += name; outputCArrays(ostr, RDAHT1Gamma, "R_dA - gamma - HighTower-1", arrayName.Data());
		arrayName = "_RDAHT2Gamma_"; arrayName += name; outputCArrays(ostr, RDAHT2Gamma, "R_dA - gamma - HighTower-2", arrayName.Data());
		ostr << "#endif" << endl;
		ostr << "" << endl;
		ostr.close();
	}
}

Float_t dAuMinBiasNBinaryCollisions = 7.5;
Float_t dAuMinBiasNBinaryCollisionsErr = 0.38; // +/- 0.4 in M. Miller's thesis
void getCentralityBin(Int_t RefMultLow, Int_t RefMultHigh, Float_t &percentLow, Float_t &percentHigh, Float_t &Nbin, Float_t &NbinErr) {
	if ((RefMultLow <= 0) && (RefMultHigh > 100)) {
		percentLow = 0; percentHigh = 100; Nbin = dAuMinBiasNBinaryCollisions; NbinErr = dAuMinBiasNBinaryCollisionsErr;
	} else if ((RefMultLow <= 0) && (RefMultHigh == 10)) {
		percentLow = 40; percentHigh = 100; Nbin = 4.0; NbinErr = 0.3;
	} else if ((RefMultLow == 10) && (RefMultHigh == 17)) {
		percentLow = 20; percentHigh = 40; Nbin = 10.0; NbinErr = 0.6; // 10.2 +1.0/-0.8 in M. Miller's thesis
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

TH1F *showArrays(Int_t n, const Float_t *x, const Float_t *y, const Float_t *ex, const Float_t *ey
	, const Char_t *name, const Char_t *title
	, TH1F *oldHist, TLegend* *legend
	, const TDrawOptions &settings
	, Bool_t show
) {
	Float_t minX = +1000;
	Float_t maxX = -1000;
	Float_t minY = +1000;
	Float_t maxY = -1000;
	if (n && x && y && ex && ey) for (Int_t i = 0;i < n;i++) {
		Float_t xx = x[i];
		Float_t xxe = ex[i];
		Float_t yy = y[i];
		Float_t yye = ey[i];
		minX = min(minX, xx - xxe);
		maxX = max(maxX, xx + xxe);
		minY = min(minY, yy - yye);
		maxY = max(maxY, yy + yye);
	}
	TH1F *histRatio = oldHist;
	TLegend *legendL = legend ? *legend : 0;
	if (n && x && y && ex && ey && show) {
	    if ((n + 0) > 0) {
		if (histRatio) {
			minX = min(minX, (Float_t)histRatio->GetXaxis()->GetBinLowEdge(1));
			maxX = max(maxX, (Float_t)histRatio->GetXaxis()->GetBinUpEdge(histRatio->GetXaxis()->GetNbins()));
			minY = min(minY, (Float_t)histRatio->GetMinimum());
			maxY = max(maxY, (Float_t)histRatio->GetMaximum());
    		} else {
			TCanvas *clist = new TCanvas(name, title);
			if (clist) {
				clist->cd();
				TString histName = TString(name) + "_hist";
				histRatio = new TH1F(histName.Data(), title, 1000, 0, 1);
				if (histRatio) {
					histRatio->SetStats(0);
					histRatio->Draw();
				}
			}
		}
		if (histRatio) {
		    histRatio->GetXaxis()->SetLimits(minX, maxX);
		    histRatio->GetYaxis()->SetLimits(minY, maxY);
		    histRatio->GetXaxis()->SetRangeUser(minX, maxX);
		    histRatio->GetYaxis()->SetRangeUser(minY, maxY);
		    histRatio->SetMinimum(minY);
		    histRatio->SetMaximum(maxY);
		}
		if (!legendL) {
		    legendL = new TLegend(0.6, 0.6, 0.89, 0.89);
		}
	    }
	    //if (legendL) legendL->SetTextFont(gStyle->GetTextFont()); // bug in ROOT? TLegend text font is not always set to gStyle->GetTextFont()
	    if (n > 0) {
		TGraphErrors *ratio = new TGraphErrors(n, x, y, ex, ey);
		if (ratio) {
		    TString gname = name; gname += "_graph";
		    ratio->SetName(gname.Data());
		    ratio->SetTitle(title);
		    settings.TAttLine::Copy(*ratio);
		    settings.TAttFill::Copy(*ratio);
		    settings.TAttMarker::Copy(*ratio);
		    ratio->SetHistogram(histRatio);
		    ratio->Draw(settings.drawOption);
		    if (legendL && (settings.legendTitle != "")) legendL->AddEntry(ratio, settings.legendTitle.Data(), settings.legendOption);
		}
	    }
	    if (legendL) legendL->Draw();
	}
	if (legend) *legend = legendL;
	return histRatio;
}

TH1F *showList(const list<TBinStatistics> *points, const Char_t *name, const Char_t *title
	, TH1F *oldHist, TLegend* *legend
	, const TDrawOptions &settings
	, Bool_t show
) {
	const Int_t bufferSize = (points ? points->size() : 0) + 1;
	Float_t *bufferX = new Float_t[bufferSize];
	Float_t *bufferY = new Float_t[bufferSize];
	Float_t *bufferXErr = new Float_t[bufferSize];
	Float_t *bufferYErr = new Float_t[bufferSize];
	Int_t numPoints = 0;
	if (points) for (list<TBinStatistics>::const_iterator iter = points->begin();iter != points->end();++iter) {
		const TBinStatistics &bin = *iter;
		Float_t x = bin.getParameters().getCenter();
		Float_t xErr = /*settings.useXerr ? (bin.getParameters().getWidth() / 2.0) : */0;
		Float_t v = bin.getValue();
		Float_t e = bin.getError();
		bufferX[numPoints] = x;
		bufferXErr[numPoints] = xErr;
		bufferY[numPoints] = v;
		bufferYErr[numPoints] = e;
		if ((numPoints + 1) < bufferSize) numPoints++;
	}
	TH1F *histRatio = showArrays(numPoints, bufferX, bufferY, bufferXErr, bufferYErr, name, title, oldHist, legend, settings, show);
	delete [] bufferX;
	delete [] bufferXErr;
	delete [] bufferY;
	delete [] bufferYErr;
	return histRatio;
}

TH1F *showLists(const list<TBinStatistics> *pointsMB, const list<TBinStatistics> *pointsHT1, const list<TBinStatistics> *pointsHT2
        , const Char_t *name, const Char_t *title
        , TH1F *oldHist, TLegend* *legend
	, const TDrawOptions &settingsMB, const TDrawOptions &settingsHT1, const TDrawOptions &settingsHT2
        , Bool_t show
) {
    TH1F *hist = oldHist;
    hist = showList(pointsMB, name, title, hist, legend, settingsMB, show);
    hist = showList(pointsHT1, name, title, hist, legend, settingsHT1, show);
    hist = showList(pointsHT2, name, title, hist, legend, settingsHT2, show);
    return hist;
}

Float_t calculateCutEff(const TCuts::cuts_map_type &cutsmap, TCuts::cut_type all, TCuts::cut_type passed) {
    Int_t totNum = 0;
    Int_t passedNum = 0;
    for (TCuts::cuts_map_type::const_iterator it = cutsmap.begin();it != cutsmap.end();++it) {
	if (PASSED((*it).first, all)) totNum += (*it).second;
	if (PASSED((*it).first, passed)) passedNum += (*it).second;
    }
    Float_t eff = (totNum != 0) ? (Float_t(passedNum) / Float_t(totNum)) : 1.0;
    cout << "Calculating cut efficiency: total " << totNum << ", passed " << passedNum << ", eff = " << eff << endl;
    return eff;
}

void calculateVertexFindingEff(const Char_t *title, Bool_t calculate, const TEventDataProcessor *event, Float_t &eff) {
    if (calculate) {
	if (title) cout << title << ", vertex finding eff: ";
	if (!event) {
	    eff = 1.0;
    	    cout << "No event data processor!";
	} else {
	    const TCuts &cuts = event->getCuts();
	    if (PASSED(cuts.EVENT_ALL_CUTS, EVENT_TPC_VERTEX_CUT)) {
		{
		    Int_t tot = event->numTotal;
	    	    Int_t passed = 0;
		    TCuts::cuts_map_type::const_iterator it = cuts.EVENT_passedCuts_separate.find(EVENT_TPC_VERTEX_CUT);
		    if (it != cuts.EVENT_passedCuts_separate.end()) {
			passed = (*it).second;
		    }
	    	    eff = (tot != 0) ? (Float_t(passed) / Float_t(tot)) : 1.0;
	    	    cout << "(1: tot " << tot << ", passed " << passed << ", eff " << eff << ")";
		}
		{
		    eff = calculateCutEff(cuts.EVENT_passedCuts_all, cuts.EVENT_ALL_CUTS & ~(EVENT_TPC_VERTEX_CUT | EVENT_PT_VS_ET_CUT), cuts.EVENT_ALL_CUTS);
	    	    cout << "(2: eff " << eff << ")";
		}
	    } else {
	        eff = 1.0;
	        cout << "No TPC vertex cut";
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

void correctJacobian(bin_stat_list_type &bins, Float_t m, Float_t eta1, Float_t eta2) {
	for (bin_stat_list_type::iterator iter = bins.begin();iter != bins.end();++iter) {
		TBinStatistics &bin = *iter;
		const TBinParameters &binParams = bin.getParameters();
		Float_t pT = binParams.getCenter();
		Float_t y1 = getRapidity(m, pT, eta1);
		Float_t y2 = getRapidity(m, pT, eta2);
		Float_t ratio = (eta2 - eta1) / (y2 - y1);
		//cout << "Jacobian corr: m = " << m << ", pT = " << pT << ", eta = " << eta1 << " - " << eta2 << ", y = " << y1 << " - " << y2 << ", corr = " << ratio << endl;
		bin.setValue(bin.getValue() * ratio);
		bin.setError(bin.getError() * ratio);
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

void showResults(const TDataProcessorPool **pools, const TAnalysisSettings &analysisSettings, TAnalysisResults &analysisResults) {
	Int_t iii = 0;
	const TDataProcessorPool *poolMB = pools[iii++];
	const TDataProcessorPool *poolHT1 = pools[iii++];
	const TDataProcessorPool *poolHT2 = pools[iii++];
	const TDataProcessorPool *poolEmbedMB = pools[iii++];
	const TDataProcessorPool *poolEmbedHT1 = pools[iii++];
	const TDataProcessorPool *poolEmbedHT2 = pools[iii++];
	const TDataProcessorPool *poolEmbedEtaMB = pools[iii++];
	const TDataProcessorPool *poolEmbedEtaHT1 = pools[iii++];
	const TDataProcessorPool *poolEmbedEtaHT2 = pools[iii++];
	const TDataProcessorPool *poolEmbedEtabgMB = pools[iii++];
	const TDataProcessorPool *poolEmbedEtabgHT1 = pools[iii++];
	const TDataProcessorPool *poolEmbedEtabgHT2 = pools[iii++];
	const TDataProcessorPool *pool1gammaMB = pools[iii++];
	const TDataProcessorPool *pool1gammaHT1 = pools[iii++];
	const TDataProcessorPool *pool1gammaHT2 = pools[iii++];
	const TDataProcessorPool *pool1gammaSimMB = pools[iii++];
	const TDataProcessorPool *pool1gammaSimHT1 = pools[iii++];
	const TDataProcessorPool *pool1gammaSimHT2 = pools[iii++];
	const TDataProcessorPool *pool1gammaSimEtaMB = pools[iii++];
	const TDataProcessorPool *pool1gammaSimEtaHT1 = pools[iii++];
	const TDataProcessorPool *pool1gammaSimEtaHT2 = pools[iii++];
	const TDataProcessorPool *poolnbarMB = pools[iii++];
	const TDataProcessorPool *poolnbarHT1 = pools[iii++];
	const TDataProcessorPool *poolnbarHT2 = pools[iii++];
	const TDataProcessorPool *poolnbarSimMB = pools[iii++];
	const TDataProcessorPool *poolnbarSimHT1 = pools[iii++];
	const TDataProcessorPool *poolnbarSimHT2 = pools[iii++];
	const TDataProcessorPool *poolnbarSimEtaMB = pools[iii++];
	const TDataProcessorPool *poolnbarSimEtaHT1 = pools[iii++];
	const TDataProcessorPool *poolnbarSimEtaHT2 = pools[iii++];
	const TDataProcessorPool *poolPtshiftMB = pools[iii++];
	const TDataProcessorPool *poolPtshiftHT1 = pools[iii++];
	const TDataProcessorPool *poolPtshiftHT2 = pools[iii++];

	TAnalysisSettings analysisSettingsNew = analysisSettings;

	if (analysisSettingsNew.print) {
		cout << "==============================================" << endl;
		cout << "Processing " << analysisSettingsNew.title << endl;
	}

	// Real data processors

	GET_DATA(candidateDataProcessorMB, const TCandidateDataProcessor, poolMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessorHT1, const TCandidateDataProcessor, poolHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessorHT2, const TCandidateDataProcessor, poolHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(candidateDataProcessorMBNocpv, const TCandidateDataProcessor, poolMB, "candidateMBNocpv");
	if (!candidateDataProcessorMBNocpv) candidateDataProcessorMBNocpv = candidateDataProcessorMB;
	GET_DATA(candidateDataProcessorHT1Nocpv, const TCandidateDataProcessor, poolHT1, "candidateHT1Nocpv");
	if (!candidateDataProcessorHT1Nocpv) candidateDataProcessorHT1Nocpv = candidateDataProcessorHT1;
	GET_DATA(candidateDataProcessorHT2Nocpv, const TCandidateDataProcessor, poolHT2, "candidateHT2Nocpv");
	if (!candidateDataProcessorHT2Nocpv) candidateDataProcessorHT2Nocpv = candidateDataProcessorHT2;

	GET_DATA(candidateDataProcessorMBNotcpv, const TCandidateDataProcessor, poolMB, "candidateMBNotcpv");
	if (!candidateDataProcessorMBNotcpv) candidateDataProcessorMBNotcpv = candidateDataProcessorMB;
	GET_DATA(candidateDataProcessorHT1Notcpv, const TCandidateDataProcessor, poolHT1, "candidateHT1Notcpv");
	if (!candidateDataProcessorHT1Notcpv) candidateDataProcessorHT1Notcpv = candidateDataProcessorHT1;
	GET_DATA(candidateDataProcessorHT2Notcpv, const TCandidateDataProcessor, poolHT2, "candidateHT2Notcpv");
	if (!candidateDataProcessorHT2Notcpv) candidateDataProcessorHT2Notcpv = candidateDataProcessorHT2;

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

	GET_DATA(candidateDataProcessorMBMixBack, const TCandidateDataProcessor, poolMB, "candidateMBJetmixBack");
	GET_DATA(candidateDataProcessorHT1MixBack, const TCandidateDataProcessor, poolHT1, "candidateHT1JetmixBack");
	GET_DATA(candidateDataProcessorHT2MixBack, const TCandidateDataProcessor, poolHT2, "candidateHT2JetmixBack");

	GET_DATA(pointDataProcessorMB, const TPointDataProcessor, poolMB, "pointMB");
	GET_DATA(pointDataProcessorHT1, const TPointDataProcessor, poolHT1, "pointHT1");
	GET_DATA(pointDataProcessorHT2, const TPointDataProcessor, poolHT2, "pointHT2");

	GET_DATA(eventDataProcessorMB, const TEventDataProcessor, poolMB, "eventMB");
	GET_DATA(eventDataProcessorHT1, const TEventDataProcessor, poolHT1, "eventHT1");
	GET_DATA(eventDataProcessorHT2, const TEventDataProcessor, poolHT2, "eventHT2");

	// Simulated data processors

	GET_DATA(candidateDataProcessorEmbedMB, const TCandidateDataProcessor, poolEmbedMB, "candidateMB");
	GET_DATA(candidateDataProcessorEmbedHT1, const TCandidateDataProcessor, poolEmbedHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorEmbedHT2, const TCandidateDataProcessor, poolEmbedHT2, "candidateHT2");

	GET_DATA(candidateDataProcessorEmbedMBNocpv, const TCandidateDataProcessor, poolEmbedMB, "candidateMBNocpv");
	if (!candidateDataProcessorEmbedMBNocpv) candidateDataProcessorEmbedMBNocpv = candidateDataProcessorEmbedMB;
	GET_DATA(candidateDataProcessorEmbedHT1Nocpv, const TCandidateDataProcessor, poolEmbedHT1, "candidateHT1Nocpv");
	if (!candidateDataProcessorEmbedHT1Nocpv) candidateDataProcessorEmbedHT1Nocpv = candidateDataProcessorEmbedHT1;
	GET_DATA(candidateDataProcessorEmbedHT2Nocpv, const TCandidateDataProcessor, poolEmbedHT2, "candidateHT2Nocpv");
	if (!candidateDataProcessorEmbedHT2Nocpv) candidateDataProcessorEmbedHT2Nocpv = candidateDataProcessorEmbedHT2;

	GET_DATA(candidateDataProcessorEmbedMBNotcpv, const TCandidateDataProcessor, poolEmbedMB, "candidateMBNotcpv");
	if (!candidateDataProcessorEmbedMBNotcpv) candidateDataProcessorEmbedMBNotcpv = candidateDataProcessorEmbedMB;
	GET_DATA(candidateDataProcessorEmbedHT1Notcpv, const TCandidateDataProcessor, poolEmbedHT1, "candidateHT1Notcpv");
	if (!candidateDataProcessorEmbedHT1Notcpv) candidateDataProcessorEmbedHT1Notcpv = candidateDataProcessorEmbedHT1;
	GET_DATA(candidateDataProcessorEmbedHT2Notcpv, const TCandidateDataProcessor, poolEmbedHT2, "candidateHT2Notcpv");
	if (!candidateDataProcessorEmbedHT2Notcpv) candidateDataProcessorEmbedHT2Notcpv = candidateDataProcessorEmbedHT2;

	GET_DATA(candidateDataProcessorEmbedMBMixRandom, const TCandidateDataProcessor, poolEmbedMB, \
	    analysisSettingsNew.settingsSimMB.useEventMixingBgRandom ? \
		"candidateMBMix" : \
		(analysisSettingsNew.settingsSimMB.useEventMixingNotmatchedBgRandom ? \
		    "candidateMBMixNotmatched" : \
		    "" \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedHT1MixRandom, const TCandidateDataProcessor, poolEmbedHT1, \
	    analysisSettingsNew.settingsSimHT1.useEventMixingBgRandom ? \
		"candidateHT1Mix" : \
		(analysisSettingsNew.settingsSimHT1.useEventMixingNotmatchedBgRandom ? \
		    "candidateHT1MixNotmatched" : \
		    "" \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedHT2MixRandom, const TCandidateDataProcessor, poolEmbedHT2, \
	    analysisSettingsNew.settingsSimHT2.useEventMixingBgRandom ? \
		"candidateHT2Mix" : \
		(analysisSettingsNew.settingsSimHT2.useEventMixingNotmatchedBgRandom ? \
		    "candidateHT2MixNotmatched" : \
		    "" \
		) \
	    );

	GET_DATA(candidateDataProcessorEmbedMBMix, const TCandidateDataProcessor, poolEmbedMB, \
	    analysisSettingsNew.settingsSimMB.useEventJetMixingBg ? \
		"candidateMBJetmix" : \
		(analysisSettingsNew.settingsSimMB.useEventJetNotmatchedMixingBg ? \
		    "candidateMBJetmixNotmatched" : \
		    (analysisSettingsNew.settingsSimMB.useEventJetBackMixingBg ? \
			"candidateMBJetmixBack" : \
			(analysisSettingsNew.settingsSimMB.useEventShufflingBg ? \
			    "candidateMBShuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedHT1Mix, const TCandidateDataProcessor, poolEmbedHT1, \
	    analysisSettingsNew.settingsSimHT1.useEventJetMixingBg ? \
		"candidateHT1Jetmix" : \
		(analysisSettingsNew.settingsSimHT1.useEventJetNotmatchedMixingBg ? \
		    "candidateHT1JetmixNotmatched" : \
		    (analysisSettingsNew.settingsSimHT1.useEventJetBackMixingBg ? \
			"candidateHT1JetmixBack" : ( \
			analysisSettingsNew.settingsSimHT1.useEventShufflingBg ? \
			    "candidateHT1Shuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedHT2Mix, const TCandidateDataProcessor, poolEmbedHT2, \
	    analysisSettingsNew.settingsSimHT2.useEventJetMixingBg ? \
		"candidateHT2Jetmix" : \
		(analysisSettingsNew.settingsSimHT2.useEventJetNotmatchedMixingBg ? \
		    "candidateHT2JetmixNotmatched" : \
		    (analysisSettingsNew.settingsSimHT2.useEventJetBackMixingBg ? \
			"candidateHT2JetmixBack" : \
			(analysisSettingsNew.settingsSimHT2.useEventShufflingBg ? \
			    "candidateHT2Shuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );

	GET_DATA(candidateDataProcessorEmbedMBMixBack, const TCandidateDataProcessor, poolEmbedMB, "candidateMBJetmixBack");
	GET_DATA(candidateDataProcessorEmbedHT1MixBack, const TCandidateDataProcessor, poolEmbedHT1, "candidateHT1JetmixBack");
	GET_DATA(candidateDataProcessorEmbedHT2MixBack, const TCandidateDataProcessor, poolEmbedHT2, "candidateHT2JetmixBack");

	GET_DATA(pointDataProcessorEmbedMB, const TPointDataProcessor, poolEmbedMB, "pointMB");
	GET_DATA(pointDataProcessorEmbedHT1, const TPointDataProcessor, poolEmbedHT1, "pointHT1");
	GET_DATA(pointDataProcessorEmbedHT2, const TPointDataProcessor, poolEmbedHT2, "pointHT2");

	GET_DATA(eventDataProcessorEmbedMB, const TEventDataProcessor, poolEmbedMB, "eventMB");
	GET_DATA(eventDataProcessorEmbedHT1, const TEventDataProcessor, poolEmbedHT1, "eventHT1");
	GET_DATA(eventDataProcessorEmbedHT2, const TEventDataProcessor, poolEmbedHT2, "eventHT2");

	GET_DATA(simuDataProcessorMB, const TSimuDataProcessor, poolEmbedMB, "simulatedMB");
	GET_DATA(simuDataProcessorHT1, const TSimuDataProcessor, poolEmbedHT1, "simulatedHT1");
	GET_DATA(simuDataProcessorHT2, const TSimuDataProcessor, poolEmbedHT2, "simulatedHT2");

	GET_DATA(gammaPi0DataProcessorMB, const TMCGammaDataProcessor, poolEmbedMB, "mcGammaMB");
	GET_DATA(gammaPi0DataProcessorHT1, const TMCGammaDataProcessor, poolEmbedHT1, "mcGammaHT1");
	GET_DATA(gammaPi0DataProcessorHT2, const TMCGammaDataProcessor, poolEmbedHT2, "mcGammaHT2");

	// Simulated Eta data processors

	GET_DATA(candidateDataProcessorEmbedEtaMB, const TCandidateDataProcessor, poolEmbedEtaMB, "candidateMB");
	GET_DATA(candidateDataProcessorEmbedEtaHT1, const TCandidateDataProcessor, poolEmbedEtaHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorEmbedEtaHT2, const TCandidateDataProcessor, poolEmbedEtaHT2, "candidateHT2");

	GET_DATA(candidateDataProcessorEmbedEtaMBNocpv, const TCandidateDataProcessor, poolEmbedEtaMB, "candidateMBNocpv");
	if (!candidateDataProcessorEmbedEtaMBNocpv) candidateDataProcessorEmbedEtaMBNocpv = candidateDataProcessorEmbedEtaMB;
	GET_DATA(candidateDataProcessorEmbedEtaHT1Nocpv, const TCandidateDataProcessor, poolEmbedEtaHT1, "candidateHT1Nocpv");
	if (!candidateDataProcessorEmbedEtaHT1Nocpv) candidateDataProcessorEmbedEtaHT1Nocpv = candidateDataProcessorEmbedEtaHT1;
	GET_DATA(candidateDataProcessorEmbedEtaHT2Nocpv, const TCandidateDataProcessor, poolEmbedEtaHT2, "candidateHT2Nocpv");
	if (!candidateDataProcessorEmbedEtaHT2Nocpv) candidateDataProcessorEmbedEtaHT2Nocpv = candidateDataProcessorEmbedEtaHT2;

	GET_DATA(candidateDataProcessorEmbedEtaMBNotcpv, const TCandidateDataProcessor, poolEmbedEtaMB, "candidateMBNotcpv");
	if (!candidateDataProcessorEmbedEtaMBNotcpv) candidateDataProcessorEmbedEtaMBNotcpv = candidateDataProcessorEmbedEtaMB;
	GET_DATA(candidateDataProcessorEmbedEtaHT1Notcpv, const TCandidateDataProcessor, poolEmbedEtaHT1, "candidateHT1Notcpv");
	if (!candidateDataProcessorEmbedEtaHT1Notcpv) candidateDataProcessorEmbedEtaHT1Notcpv = candidateDataProcessorEmbedEtaHT1;
	GET_DATA(candidateDataProcessorEmbedEtaHT2Notcpv, const TCandidateDataProcessor, poolEmbedEtaHT2, "candidateHT2Notcpv");
	if (!candidateDataProcessorEmbedEtaHT2Notcpv) candidateDataProcessorEmbedEtaHT2Notcpv = candidateDataProcessorEmbedEtaHT2;

	GET_DATA(candidateDataProcessorEmbedEtaMBMixRandom, const TCandidateDataProcessor, poolEmbedEtaMB, \
	    analysisSettingsNew.settingsSimEtaMB.useEventMixingBgRandom ? \
		"candidateMBMix" : \
		(analysisSettingsNew.settingsSimEtaMB.useEventMixingNotmatchedBgRandom ? \
		    "candidateMBMixNotmatched" : \
		    "" \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedEtaHT1MixRandom, const TCandidateDataProcessor, poolEmbedEtaHT1, \
	    analysisSettingsNew.settingsSimEtaHT1.useEventMixingBgRandom ? \
		"candidateHT1Mix" : \
		(analysisSettingsNew.settingsSimEtaHT1.useEventMixingNotmatchedBgRandom ? \
		    "candidateHT1MixNotmatched" : \
		    "" \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedEtaHT2MixRandom, const TCandidateDataProcessor, poolEmbedEtaHT2, \
	    analysisSettingsNew.settingsSimEtaHT2.useEventMixingBgRandom ? \
		"candidateHT2Mix" : \
		(analysisSettingsNew.settingsSimEtaHT2.useEventMixingNotmatchedBgRandom ? \
		    "candidateHT2MixNotmatched" : \
		    "" \
		) \
	    );

	GET_DATA(candidateDataProcessorEmbedEtaMBMix, const TCandidateDataProcessor, poolEmbedEtaMB, \
	    analysisSettingsNew.settingsSimEtaMB.useEventJetMixingBg ? \
		"candidateMBJetmix" : \
		(analysisSettingsNew.settingsSimEtaMB.useEventJetNotmatchedMixingBg ? \
		    "candidateMBJetmixNotmatched" : \
		    (analysisSettingsNew.settingsSimEtaMB.useEventJetBackMixingBg ? \
			"candidateMBJetmixBack" : \
			(analysisSettingsNew.settingsSimEtaMB.useEventShufflingBg ? \
			    "candidateMBShuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedEtaHT1Mix, const TCandidateDataProcessor, poolEmbedEtaHT1, \
	    analysisSettingsNew.settingsSimEtaHT1.useEventJetMixingBg ? \
		"candidateHT1Jetmix" : \
		(analysisSettingsNew.settingsSimEtaHT1.useEventJetNotmatchedMixingBg ? \
		    "candidateHT1JetmixNotmatched" : \
		    (analysisSettingsNew.settingsSimEtaHT1.useEventJetBackMixingBg ? \
			"candidateHT1JetmixBack" : ( \
			analysisSettingsNew.settingsSimEtaHT1.useEventShufflingBg ? \
			    "candidateHT1Shuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );
	GET_DATA(candidateDataProcessorEmbedEtaHT2Mix, const TCandidateDataProcessor, poolEmbedEtaHT2, \
	    analysisSettingsNew.settingsSimEtaHT2.useEventJetMixingBg ? \
		"candidateHT2Jetmix" : \
		(analysisSettingsNew.settingsSimEtaHT2.useEventJetNotmatchedMixingBg ? \
		    "candidateHT2JetmixNotmatched" : \
		    (analysisSettingsNew.settingsSimEtaHT2.useEventJetBackMixingBg ? \
			"candidateHT2JetmixBack" : \
			(analysisSettingsNew.settingsSimEtaHT2.useEventShufflingBg ? \
			    "candidateHT2Shuffle" : \
			    "" \
			) \
		    ) \
		) \
	    );

	GET_DATA(candidateDataProcessorEmbedEtaMBMixBack, const TCandidateDataProcessor, poolEmbedEtaMB, "candidateMBJetmixBack");
	GET_DATA(candidateDataProcessorEmbedEtaHT1MixBack, const TCandidateDataProcessor, poolEmbedEtaHT1, "candidateHT1JetmixBack");
	GET_DATA(candidateDataProcessorEmbedEtaHT2MixBack, const TCandidateDataProcessor, poolEmbedEtaHT2, "candidateHT2JetmixBack");

	GET_DATA(pointDataProcessorEmbedEtaMB, const TPointDataProcessor, poolEmbedEtaMB, "pointMB");
	GET_DATA(pointDataProcessorEmbedEtaHT1, const TPointDataProcessor, poolEmbedEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessorEmbedEtaHT2, const TPointDataProcessor, poolEmbedEtaHT2, "pointHT2");

	GET_DATA(eventDataProcessorEmbedEtaMB, const TEventDataProcessor, poolEmbedEtaMB, "eventMB");
	GET_DATA(eventDataProcessorEmbedEtaHT1, const TEventDataProcessor, poolEmbedEtaHT1, "eventHT1");
	GET_DATA(eventDataProcessorEmbedEtaHT2, const TEventDataProcessor, poolEmbedEtaHT2, "eventHT2");

	GET_DATA(simuDataProcessorEtaMB, const TSimuDataProcessor, poolEmbedEtaMB, "simulatedMB");
	GET_DATA(simuDataProcessorEtaHT1, const TSimuDataProcessor, poolEmbedEtaHT1, "simulatedHT1");
	GET_DATA(simuDataProcessorEtaHT2, const TSimuDataProcessor, poolEmbedEtaHT2, "simulatedHT2");

	GET_DATA(gammaEtaDataProcessorMB, const TMCGammaDataProcessor, poolEmbedEtaMB, "mcGammaMB");
	GET_DATA(gammaEtaDataProcessorHT1, const TMCGammaDataProcessor, poolEmbedEtaHT1, "mcGammaHT1");
	GET_DATA(gammaEtaDataProcessorHT2, const TMCGammaDataProcessor, poolEmbedEtaHT2, "mcGammaHT2");

	// Simulated Etabg data processors

	GET_DATA(candidateDataProcessorEmbedEtabgMB, const TCandidateDataProcessor, poolEmbedEtabgMB, "candidateMB");
	GET_DATA(candidateDataProcessorEmbedEtabgHT1, const TCandidateDataProcessor, poolEmbedEtabgHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorEmbedEtabgHT2, const TCandidateDataProcessor, poolEmbedEtabgHT2, "candidateHT2");

	GET_DATA(pointDataProcessorEmbedEtabgMB, const TPointDataProcessor, poolEmbedEtabgMB, "pointMB");
	GET_DATA(pointDataProcessorEmbedEtabgHT1, const TPointDataProcessor, poolEmbedEtabgHT1, "pointHT1");
	GET_DATA(pointDataProcessorEmbedEtabgHT2, const TPointDataProcessor, poolEmbedEtabgHT2, "pointHT2");

	GET_DATA(eventDataProcessorEmbedEtabgMB, const TEventDataProcessor, poolEmbedEtabgMB, "eventMB");
	GET_DATA(eventDataProcessorEmbedEtabgHT1, const TEventDataProcessor, poolEmbedEtabgHT1, "eventHT1");
	GET_DATA(eventDataProcessorEmbedEtabgHT2, const TEventDataProcessor, poolEmbedEtabgHT2, "eventHT2");

	GET_DATA(simuDataProcessorEtabgMB, const TSimuDataProcessor, poolEmbedEtabgMB, "simulatedMB");
	GET_DATA(simuDataProcessorEtabgHT1, const TSimuDataProcessor, poolEmbedEtabgHT1, "simulatedHT1");
	GET_DATA(simuDataProcessorEtabgHT2, const TSimuDataProcessor, poolEmbedEtabgHT2, "simulatedHT2");

	// Simulated 1gamma data processors

	GET_DATA(candidateDataProcessor1gammaMB, const TCandidateDataProcessor, pool1gammaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1, const TCandidateDataProcessor, pool1gammaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2, const TCandidateDataProcessor, pool1gammaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessor1gammaMB, const TPointDataProcessor, pool1gammaMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1, const TPointDataProcessor, pool1gammaHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2, const TPointDataProcessor, pool1gammaHT2, "pointHT2");

	GET_DATA(gammaDataProcessorMB, const TMCGammaDataProcessor, pool1gammaMB, "mcGammaMB");
	GET_DATA(gammaDataProcessorHT1, const TMCGammaDataProcessor, pool1gammaHT1, "mcGammaHT1");
	GET_DATA(gammaDataProcessorHT2, const TMCGammaDataProcessor, pool1gammaHT2, "mcGammaHT2");

	// Simulated 1gamma pi0 data processors

	GET_DATA(candidateDataProcessor1gammaMBSim, const TCandidateDataProcessor, pool1gammaSimMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1Sim, const TCandidateDataProcessor, pool1gammaSimHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2Sim, const TCandidateDataProcessor, pool1gammaSimHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessor1gammaMBSim, const TPointDataProcessor, pool1gammaSimMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1Sim, const TPointDataProcessor, pool1gammaSimHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2Sim, const TPointDataProcessor, pool1gammaSimHT2, "pointHT2");

	// Simulated 1gamma eta data processors

	GET_DATA(candidateDataProcessor1gammaMBSimEta, const TCandidateDataProcessor, pool1gammaSimEtaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1SimEta, const TCandidateDataProcessor, pool1gammaSimEtaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2SimEta, const TCandidateDataProcessor, pool1gammaSimEtaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessor1gammaMBSimEta, const TPointDataProcessor, pool1gammaSimEtaMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1SimEta, const TPointDataProcessor, pool1gammaSimEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2SimEta, const TPointDataProcessor, pool1gammaSimEtaHT2, "pointHT2");

	// Simulated 1gamma etabg data processors

	GET_DATA(candidateDataProcessor1gammaMBSimEtabg, const TCandidateDataProcessor, pool1gammaSimEtaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessor1gammaHT1SimEtabg, const TCandidateDataProcessor, pool1gammaSimEtaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessor1gammaHT2SimEtabg, const TCandidateDataProcessor, pool1gammaSimEtaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessor1gammaMBSimEtabg, const TPointDataProcessor, pool1gammaSimEtaMB, "pointMB");
	GET_DATA(pointDataProcessor1gammaHT1SimEtabg, const TPointDataProcessor, pool1gammaSimEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessor1gammaHT2SimEtabg, const TPointDataProcessor, pool1gammaSimEtaHT2, "pointHT2");

	// Simulated antineutron data processors

	GET_DATA(candidateDataProcessornbarMB, const TCandidateDataProcessor, poolnbarMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessornbarHT1, const TCandidateDataProcessor, poolnbarHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessornbarHT2, const TCandidateDataProcessor, poolnbarHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessornbarMB, const TPointDataProcessor, poolnbarMB, "pointMB");
	GET_DATA(pointDataProcessornbarHT1, const TPointDataProcessor, poolnbarHT1, "pointHT1");
	GET_DATA(pointDataProcessornbarHT2, const TPointDataProcessor, poolnbarHT2, "pointHT2");

	GET_DATA(nbarDataProcessorMB, const TMCGammaDataProcessor, poolnbarMB, "mcGammaMB");
	GET_DATA(nbarDataProcessorHT1, const TMCGammaDataProcessor, poolnbarHT1, "mcGammaHT1");
	GET_DATA(nbarDataProcessorHT2, const TMCGammaDataProcessor, poolnbarHT2, "mcGammaHT2");

	// Simulated antineutron pi0 data processors

	GET_DATA(candidateDataProcessornbarMBSim, const TCandidateDataProcessor, poolnbarSimMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessornbarHT1Sim, const TCandidateDataProcessor, poolnbarSimHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessornbarHT2Sim, const TCandidateDataProcessor, poolnbarSimHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessornbarMBSim, const TPointDataProcessor, poolnbarSimMB, "pointMB");
	GET_DATA(pointDataProcessornbarHT1Sim, const TPointDataProcessor, poolnbarSimHT1, "pointHT1");
	GET_DATA(pointDataProcessornbarHT2Sim, const TPointDataProcessor, poolnbarSimHT2, "pointHT2");

	// Simulated antineutron eta data processors

	GET_DATA(candidateDataProcessornbarMBSimEta, const TCandidateDataProcessor, poolnbarSimEtaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessornbarHT1SimEta, const TCandidateDataProcessor, poolnbarSimEtaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessornbarHT2SimEta, const TCandidateDataProcessor, poolnbarSimEtaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessornbarMBSimEta, const TPointDataProcessor, poolnbarSimEtaMB, "pointMB");
	GET_DATA(pointDataProcessornbarHT1SimEta, const TPointDataProcessor, poolnbarSimEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessornbarHT2SimEta, const TPointDataProcessor, poolnbarSimEtaHT2, "pointHT2");

	// Simulated antineutron etabg data processors

	GET_DATA(candidateDataProcessornbarMBSimEtabg, const TCandidateDataProcessor, poolnbarSimEtaMB, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakMB" : "candidateMB"));
	GET_DATA(candidateDataProcessornbarHT1SimEtabg, const TCandidateDataProcessor, poolnbarSimEtaHT1, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT1" : "candidateHT1"));
	GET_DATA(candidateDataProcessornbarHT2SimEtabg, const TCandidateDataProcessor, poolnbarSimEtaHT2, (analysisSettingsNew.multiplicityDistributions ? "candidatePeakHT2" : "candidateHT2"));

	GET_DATA(pointDataProcessornbarMBSimEtabg, const TPointDataProcessor, poolnbarSimEtaMB, "pointMB");
	GET_DATA(pointDataProcessornbarHT1SimEtabg, const TPointDataProcessor, poolnbarSimEtaHT1, "pointHT1");
	GET_DATA(pointDataProcessornbarHT2SimEtabg, const TPointDataProcessor, poolnbarSimEtaHT2, "pointHT2");

	// Data processors for the pT shift calculation

	GET_DATA(candidateDataProcessorPtshiftMB, const TCandidateDataProcessor, poolPtshiftMB, "candidateMB");
	GET_DATA(candidateDataProcessorPtshiftHT1, const TCandidateDataProcessor, poolPtshiftHT1, "candidateHT1");
	GET_DATA(candidateDataProcessorPtshiftHT2, const TCandidateDataProcessor, poolPtshiftHT2, "candidateHT2");

	GET_DATA(pointDataProcessorPtshiftMB, const TPointDataProcessor, poolPtshiftMB, "pointMB");
	GET_DATA(pointDataProcessorPtshiftHT1, const TPointDataProcessor, poolPtshiftHT1, "pointHT1");
	GET_DATA(pointDataProcessorPtshiftHT2, const TPointDataProcessor, poolPtshiftHT2, "pointHT2");


	if (analysisSettings.saveBadSimPi0EventsMB && simuDataProcessorMB) {
	    writeEventListToFile(analysisSettings.badSimPi0EventsMBFilename, &simuDataProcessorMB->badEvents);
	}
	if (analysisSettings.saveBadSimPi0EventsHT1 && simuDataProcessorHT1) {
	    writeEventListToFile(analysisSettings.badSimPi0EventsHT1Filename, &simuDataProcessorHT1->badEvents);
	}
	if (analysisSettings.saveBadSimPi0EventsHT2 && simuDataProcessorHT2) {
	    writeEventListToFile(analysisSettings.badSimPi0EventsHT2Filename, &simuDataProcessorHT2->badEvents);
	}

	clearCachedPrescales();

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
	    PSMB += analysisSettingsNew.settingsMB.fixedPrescale / (1.0 * 1.0);
	    PSMB_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsHT1.fixedPrescale > 0) {
	    PSHT1 += analysisSettingsNew.settingsHT1.fixedPrescale / (1.0 * 1.0);
	    PSHT1_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsHT2.fixedPrescale > 0) {
	    PSHT2 += analysisSettingsNew.settingsHT2.fixedPrescale / (1.0 * 1.0);
	    PSHT2_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsMB.calculatePrescaleFromDB || analysisSettingsNew.settingsHT1.calculatePrescaleFromDB || analysisSettingsNew.settingsHT2.calculatePrescaleFromDB) {
		Float_t PSHT1_tmp, PSHT2_tmp;
		Float_t PSHT1_tmp_err, PSHT2_tmp_err;
		Float_t PSHT1_tmp_alexst, PSHT2_tmp_alexst;
		calculatePSFromDB(analysisSettingsNew.name + "_data", analysisSettingsNew.title + " (data)", analysisSettingsNew.prescalesDBFilename, eventDataProcessorMB, eventDataProcessorHT1, eventDataProcessorHT2, analysisSettingsNew.showReal, analysisSettingsNew.print, PSHT1_tmp, PSHT2_tmp, PSHT1_tmp_err, PSHT2_tmp_err, PSHT1_tmp_alexst, PSHT2_tmp_alexst);
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromDB) {
		    PSMB += 1.0 / (1.0 * 1.0);
		    PSMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromDB && (PSHT1_tmp > 0)) {
		    PSHT1 += PSHT1_tmp / (PSHT1_tmp_err * PSHT1_tmp_err);
		    PSHT1_num += 1.0 / (PSHT1_tmp_err * PSHT1_tmp_err);
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromDB && (PSHT2_tmp > 0)) {
		    PSHT2 += PSHT2_tmp / (PSHT2_tmp_err * PSHT2_tmp_err);
		    PSHT2_num += 1.0 / (PSHT2_tmp_err * PSHT2_tmp_err);
		}
	}
	if (analysisSettingsNew.settingsMB.calculatePrescaleFromSim || analysisSettingsNew.settingsHT1.calculatePrescaleFromSim || analysisSettingsNew.settingsHT2.calculatePrescaleFromSim) {
		Float_t PSHT1_tmp, PSHT2_tmp;
		Float_t PSHT1_err_tmp, PSHT2_err_tmp;
		calculatePSFromSim(analysisSettingsNew.name + "_data", analysisSettingsNew.title + " (data)", poolMB, poolHT1
				, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0
				, analysisSettingsNew.showReal, analysisSettingsNew.print, PSHT1_tmp, PSHT2_tmp, PSHT1_err_tmp, PSHT2_err_tmp);
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromSim) {
		    PSMB += 1.0 / (1.0 * 1.0);
		    PSMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromSim && (PSHT1_tmp > 0)) {
		    PSHT1 += PSHT1_tmp / (PSHT1_err_tmp * PSHT1_err_tmp);
		    PSHT1_num += 1.0 / (PSHT1_err_tmp * PSHT1_err_tmp);
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromSim && (PSHT2_tmp > 0)) {
		    PSHT2 += PSHT2_tmp / (PSHT2_err_tmp * PSHT2_err_tmp);
		    PSHT2_num += 1.0 / (PSHT2_err_tmp * PSHT2_err_tmp);
		}
	}
	if (analysisSettingsNew.settingsMB.calculatePrescaleFromPoints || analysisSettingsNew.settingsHT1.calculatePrescaleFromPoints || analysisSettingsNew.settingsHT2.calculatePrescaleFromPoints || analysisSettingsNew.settingsMB.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsHT1.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsHT2.calculatePrescaleFromPointsIntegral) {
		Float_t PSHT1_tmp, PSHT2_tmp, PSHT1Integral_tmp, PSHT2Integral_tmp;
		Float_t PSHT1_tmp_err, PSHT2_tmp_err, PSHT1Integral_tmp_err, PSHT2Integral_tmp_err;
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
		    , PSHT1_tmp_err, PSHT2_tmp_err
		    , PSHT1Integral_tmp, PSHT2Integral_tmp
		    , PSHT1Integral_tmp_err, PSHT2Integral_tmp_err);
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromPoints) {
		    PSMB += 1.0 / (1.0 * 1.0);
		    PSMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromPoints && (PSHT1_tmp > 0)) {
		    PSHT1 += PSHT1_tmp / (PSHT1_tmp_err * PSHT1_tmp_err);
		    PSHT1_num += 1.0 / (PSHT1_tmp_err * PSHT1_tmp_err);
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromPoints && (PSHT2_tmp > 0)) {
		    PSHT2 += PSHT2_tmp / (PSHT2_tmp_err * PSHT2_tmp_err);
		    PSHT2_num += 1.0 / (PSHT2_tmp_err * PSHT2_tmp_err);
		}
		if (analysisSettingsNew.settingsMB.calculatePrescaleFromPointsIntegral) {
		    PSMB += 1.0 / (1.0 * 1.0);
		    PSMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsHT1.calculatePrescaleFromPointsIntegral && (PSHT1Integral_tmp > 0)) {
		    PSHT1 += PSHT1Integral_tmp / (PSHT1Integral_tmp_err * PSHT1Integral_tmp_err);
		    PSHT1_num += 1.0 / (PSHT1Integral_tmp_err * PSHT1Integral_tmp_err);
		}
		if (analysisSettingsNew.settingsHT2.calculatePrescaleFromPointsIntegral && (PSHT2Integral_tmp > 0)) {
		    PSHT2 += PSHT2Integral_tmp / (PSHT2Integral_tmp_err * PSHT2Integral_tmp_err);
		    PSHT2_num += 1.0 / (PSHT2Integral_tmp_err * PSHT2Integral_tmp_err);
		}
	}
	if (PSMB_num != 0) {
	    PSMB = PSMB / PSMB_num;
	    PSMB_num = TMath::Sqrt(1.0 / PSMB_num);
	}
	if (PSHT1_num != 0) {
	    PSHT1 = PSHT1 / PSHT1_num;
	    PSHT1_num = TMath::Sqrt(1.0 / PSHT1_num);
	}
	if (PSHT2_num != 0) {
	    PSHT2 = PSHT2 / PSHT2_num;
	    PSHT2_num = TMath::Sqrt(1.0 / PSHT2_num);
	}
	if (analysisSettingsNew.print) cout << "HighTower enhancement factors, " << analysisSettingsNew.title << ": MB = " << PSMB << " +/- " << PSMB_num << ", HT1 = " << PSHT1 << " +/- " << PSHT1_num << ", HT2 = " << PSHT2 << " +/- " << PSHT2_num << endl;
	analysisResults.resultsMB.prescaleUsed = PSMB;
	analysisResults.resultsHT1.prescaleUsed = PSHT1;
	analysisResults.resultsHT2.prescaleUsed = PSHT2;
	}

	{
	Float_t PSSimMB = 0, PSSimHT1 = 0, PSSimHT2 = 0;
	Float_t PSSimMB_num = 0, PSSimHT1_num = 0, PSSimHT2_num = 0;
	if (analysisSettingsNew.settingsSimMB.fixedPrescale > 0) {
	    PSSimMB += analysisSettingsNew.settingsSimMB.fixedPrescale / (1.0 * 1.0);
	    PSSimMB_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsSimHT1.fixedPrescale > 0) {
	    PSSimHT1 += analysisSettingsNew.settingsSimHT1.fixedPrescale / (1.0 * 1.0);
	    PSSimHT1_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsSimHT2.fixedPrescale > 0) {
	    PSSimHT2 += analysisSettingsNew.settingsSimHT2.fixedPrescale / (1.0 * 1.0);
	    PSSimHT2_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromDB || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromDB || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromDB) {
		Float_t PSSimHT1_tmp, PSSimHT2_tmp;
		Float_t PSSimHT1_tmp_err, PSSimHT2_tmp_err;
		Float_t PSSimHT1_tmp_alexst, PSSimHT2_tmp_alexst;
		calculatePSFromDB(analysisSettingsNew.name + "_sim", analysisSettingsNew.title + " (sim)", analysisSettingsNew.prescalesDBFilename, eventDataProcessorEmbedMB, eventDataProcessorEmbedHT1, eventDataProcessorEmbedHT2, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimHT1_tmp, PSSimHT2_tmp, PSSimHT1_tmp_err, PSSimHT2_tmp_err, PSSimHT1_tmp_alexst, PSSimHT2_tmp_alexst);
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromDB) {
		    PSSimMB += 1.0 / (1.0 * 1.0);
		    PSSimMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromDB && (PSSimHT1_tmp > 0)) {
		    PSSimHT1 += PSSimHT1_tmp / (PSSimHT1_tmp_err * PSSimHT1_tmp_err);
		    PSSimHT1_num += 1.0 / (PSSimHT1_tmp_err * PSSimHT1_tmp_err);
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromDB && (PSSimHT2_tmp > 0)) {
		    PSSimHT2 += PSSimHT2_tmp / (PSSimHT2_tmp_err * PSSimHT2_tmp_err);
		    PSSimHT2_num += 1.0 / (PSSimHT2_tmp_err * PSSimHT2_tmp_err);
		}
	}
	if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromSim || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromSim || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromSim) {
		Float_t PSSimHT1_tmp, PSSimHT2_tmp;
		Float_t PSSimHT1_err_tmp, PSSimHT2_err_tmp;
		calculatePSFromSim(analysisSettingsNew.name + "_sim", analysisSettingsNew.title + " (sim)", poolEmbedMB, poolEmbedHT1
				, 0, 0, 0/*, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0*/
				, true || analysisSettingsNew.showSim, true || analysisSettingsNew.print, PSSimHT1_tmp, PSSimHT2_tmp, PSSimHT1_err_tmp, PSSimHT2_err_tmp);
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromSim) {
		    PSSimMB += 1.0 / (1.0 * 1.0);
		    PSSimMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromSim && (PSSimHT1_tmp > 0)) {
		    PSSimHT1 += PSSimHT1_tmp / (PSSimHT1_err_tmp * PSSimHT1_err_tmp);
		    PSSimHT1_num += 1.0 / (PSSimHT1_err_tmp * PSSimHT1_err_tmp);
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromSim && (PSSimHT2_tmp > 0)) {
		    PSSimHT2 += PSSimHT2_tmp / (PSSimHT2_err_tmp * PSSimHT2_err_tmp);
		    PSSimHT2_num += 1.0 / (PSSimHT2_err_tmp * PSSimHT2_err_tmp);
		}
	}
	if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimMB.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPointsIntegral) {
		Float_t PSSimHT1_tmp, PSSimHT2_tmp, PSSimHT1Integral_tmp, PSSimHT2Integral_tmp;
		Float_t PSSimHT1_tmp_err, PSSimHT2_tmp_err, PSSimHT1Integral_tmp_err, PSSimHT2Integral_tmp_err;
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
		    , PSSimHT1_tmp_err, PSSimHT2_tmp_err
		    , PSSimHT1Integral_tmp, PSSimHT2Integral_tmp
		    , PSSimHT1Integral_tmp_err, PSSimHT2Integral_tmp_err);
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromPoints) {
		    PSSimMB += 1.0 / (1.0 * 1.0);
		    PSSimMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPoints && (PSSimHT1_tmp > 0)) {
		    PSSimHT1 += PSSimHT1_tmp / (PSSimHT1_tmp_err * PSSimHT1_tmp_err);
		    PSSimHT1_num += 1.0 / (PSSimHT1_tmp_err * PSSimHT1_tmp_err);
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPoints && (PSSimHT2_tmp > 0)) {
		    PSSimHT2 += PSSimHT2_tmp / (PSSimHT2_tmp_err * PSSimHT2_tmp_err);
		    PSSimHT2_num += 1.0 / (PSSimHT2_tmp_err * PSSimHT2_tmp_err);
		}
		if (analysisSettingsNew.settingsSimMB.calculatePrescaleFromPointsIntegral) {
		    PSSimMB += 1.0 / (1.0 * 1.0);
		    PSSimMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimHT1.calculatePrescaleFromPointsIntegral && (PSSimHT1Integral_tmp > 0)) {
		    PSSimHT1 += PSSimHT1Integral_tmp / (PSSimHT1Integral_tmp_err * PSSimHT1Integral_tmp_err);
		    PSSimHT1_num += 1.0 / (PSSimHT1Integral_tmp_err * PSSimHT1Integral_tmp_err);
		}
		if (analysisSettingsNew.settingsSimHT2.calculatePrescaleFromPointsIntegral && (PSSimHT2Integral_tmp > 0)) {
		    PSSimHT2 += PSSimHT2Integral_tmp / (PSSimHT2Integral_tmp_err * PSSimHT2Integral_tmp_err);
		    PSSimHT2_num += 1.0 / (PSSimHT2Integral_tmp_err * PSSimHT2Integral_tmp_err);
		}
	}
	if (PSSimMB_num != 0) {
	    PSSimMB = PSSimMB / PSSimMB_num;
	    PSSimMB_num = TMath::Sqrt(1.0 / PSSimMB_num);
	}
	if (PSSimHT1_num != 0) {
	    PSSimHT1 = PSSimHT1 / PSSimHT1_num;
	    PSSimHT1_num = TMath::Sqrt(1.0 / PSSimHT1_num);
	}
	if (PSSimHT2_num != 0) {
	    PSSimHT2 = PSSimHT2 / PSSimHT2_num;
	    PSSimHT2_num = TMath::Sqrt(1.0 / PSSimHT2_num);
	}
	if (analysisSettingsNew.print) cout << "HighTower enhancement factors (sim), " << analysisSettingsNew.title << ": MB = " << PSSimMB << " +/- " << PSSimMB_num << ", HT1 = " << PSSimHT1 << " +/- " << PSSimHT1_num << ", HT2 = " << PSSimHT2 << " +/- " << PSSimHT2_num << endl;
	analysisResults.resultsSimMB.prescaleUsed = PSSimMB;
	analysisResults.resultsSimHT1.prescaleUsed = PSSimHT1;
	analysisResults.resultsSimHT2.prescaleUsed = PSSimHT2;
	}

	{
	Float_t PSSimEtaMB = 0, PSSimEtaHT1 = 0, PSSimEtaHT2 = 0;
	Float_t PSSimEtaMB_num = 0, PSSimEtaHT1_num = 0, PSSimEtaHT2_num = 0;
	if (analysisSettingsNew.settingsSimEtaMB.fixedPrescale > 0) {
	    PSSimEtaMB += analysisSettingsNew.settingsSimEtaMB.fixedPrescale / (1.0 * 1.0);
	    PSSimEtaMB_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsSimEtaHT1.fixedPrescale > 0) {
	    PSSimEtaHT1 += analysisSettingsNew.settingsSimEtaHT1.fixedPrescale / (1.0 * 1.0);
	    PSSimEtaHT1_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsSimEtaHT2.fixedPrescale > 0) {
	    PSSimEtaHT2 += analysisSettingsNew.settingsSimEtaHT2.fixedPrescale / (1.0 * 1.0);
	    PSSimEtaHT2_num += 1.0 / (1.0 * 1.0);
	}
	if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromDB || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromDB || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromDB) {
		Float_t PSSimEtaHT1_tmp, PSSimEtaHT2_tmp;
		Float_t PSSimEtaHT1_tmp_err, PSSimEtaHT2_tmp_err;
		Float_t PSSimEtaHT1_tmp_alexst, PSSimEtaHT2_tmp_alexst;
		calculatePSFromDB(analysisSettingsNew.name + "_simEta", analysisSettingsNew.title + " (sim eta)", analysisSettingsNew.prescalesDBFilename, eventDataProcessorEmbedEtaMB, eventDataProcessorEmbedEtaHT1, eventDataProcessorEmbedEtaHT2, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimEtaHT1_tmp, PSSimEtaHT2_tmp, PSSimEtaHT1_tmp_err, PSSimEtaHT2_tmp_err, PSSimEtaHT1_tmp_alexst, PSSimEtaHT2_tmp_alexst);
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromDB) {
		    PSSimEtaMB += 1.0 / (1.0 * 1.0);
		    PSSimEtaMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromDB && (PSSimEtaHT1_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1_tmp / (PSSimEtaHT1_tmp_err * PSSimEtaHT1_tmp_err);
		    PSSimEtaHT1_num += 1.0 / (PSSimEtaHT1_tmp_err * PSSimEtaHT1_tmp_err);
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromDB && (PSSimEtaHT2_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2_tmp / (PSSimEtaHT2_tmp_err * PSSimEtaHT2_tmp_err);
		    PSSimEtaHT2_num += 1.0 / (PSSimEtaHT2_tmp_err * PSSimEtaHT2_tmp_err);
		}
	}
	if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromSim || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromSim || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromSim) {
		Float_t PSSimEtaHT1_tmp, PSSimEtaHT2_tmp;
		Float_t PSSimEtaHT1_err_tmp, PSSimEtaHT2_err_tmp;
		calculatePSFromSim(analysisSettingsNew.name + "_simEta", analysisSettingsNew.title + " (sim eta)", poolEmbedEtaMB, poolEmbedEtaHT1
				, 0, 0, 0/*, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT1 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasMBHT2 : 0, analysisSettingsNew.correctPrescaleTriggerBias ? analysisResults.prescaleTriggerBiasHT1HT2 : 0*/
				, analysisSettingsNew.showSim, analysisSettingsNew.print, PSSimEtaHT1_tmp, PSSimEtaHT2_tmp, PSSimEtaHT1_err_tmp, PSSimEtaHT2_err_tmp);
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromSim) {
		    PSSimEtaMB += 1.0 / (1.0 * 1.0);
		    PSSimEtaMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromSim && (PSSimEtaHT1_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1_tmp / (PSSimEtaHT1_err_tmp * PSSimEtaHT1_err_tmp);
		    PSSimEtaHT1_num += 1.0 / (PSSimEtaHT1_err_tmp * PSSimEtaHT1_err_tmp);
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromSim && (PSSimEtaHT2_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2_tmp / (PSSimEtaHT2_err_tmp * PSSimEtaHT2_err_tmp);
		    PSSimEtaHT2_num += 1.0 / (PSSimEtaHT2_err_tmp * PSSimEtaHT2_err_tmp);
		}
	}
	if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPoints || analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPointsIntegral || analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPointsIntegral) {
		Float_t PSSimEtaHT1_tmp, PSSimEtaHT2_tmp, PSSimEtaHT1Integral_tmp, PSSimEtaHT2Integral_tmp;
		Float_t PSSimEtaHT1_tmp_err, PSSimEtaHT2_tmp_err, PSSimEtaHT1Integral_tmp_err, PSSimEtaHT2Integral_tmp_err;
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
		    , PSSimEtaHT1_tmp_err, PSSimEtaHT2_tmp_err
		    , PSSimEtaHT1Integral_tmp, PSSimEtaHT2Integral_tmp
		    , PSSimEtaHT1Integral_tmp_err, PSSimEtaHT2Integral_tmp_err);
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPoints) {
		    PSSimEtaMB += 1.0 / (1.0 * 1.0);
		    PSSimEtaMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPoints && (PSSimEtaHT1_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1_tmp / (PSSimEtaHT1_tmp_err * PSSimEtaHT1_tmp_err);
		    PSSimEtaHT1_num += 1.0 / (PSSimEtaHT1_tmp_err * PSSimEtaHT1_tmp_err);
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPoints && (PSSimEtaHT2_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2_tmp / (PSSimEtaHT2_tmp_err * PSSimEtaHT2_tmp_err);
		    PSSimEtaHT2_num += 1.0 / (PSSimEtaHT2_tmp_err * PSSimEtaHT2_tmp_err);
		}
		if (analysisSettingsNew.settingsSimEtaMB.calculatePrescaleFromPointsIntegral) {
		    PSSimEtaMB += 1.0 / (1.0 * 1.0);
		    PSSimEtaMB_num += 1.0 / (1.0 * 1.0);
		}
		if (analysisSettingsNew.settingsSimEtaHT1.calculatePrescaleFromPointsIntegral && (PSSimEtaHT1Integral_tmp > 0)) {
		    PSSimEtaHT1 += PSSimEtaHT1Integral_tmp / (PSSimEtaHT1Integral_tmp_err * PSSimEtaHT1Integral_tmp_err);
		    PSSimEtaHT1_num += 1.0 / (PSSimEtaHT1Integral_tmp_err * PSSimEtaHT1Integral_tmp_err);
		}
		if (analysisSettingsNew.settingsSimEtaHT2.calculatePrescaleFromPointsIntegral && (PSSimEtaHT2Integral_tmp > 0)) {
		    PSSimEtaHT2 += PSSimEtaHT2Integral_tmp / (PSSimEtaHT2Integral_tmp_err * PSSimEtaHT2Integral_tmp_err);
		    PSSimEtaHT2_num += 1.0 / (PSSimEtaHT2Integral_tmp_err * PSSimEtaHT2Integral_tmp_err);
		}
	}
	if (PSSimEtaMB_num != 0) {
	    PSSimEtaMB = PSSimEtaMB / PSSimEtaMB_num;
	    PSSimEtaMB_num = TMath::Sqrt(1.0 / PSSimEtaMB_num);
	}
	if (PSSimEtaHT1_num != 0) {
	    PSSimEtaHT1 = PSSimEtaHT1 / PSSimEtaHT1_num;
	    PSSimEtaHT1_num = TMath::Sqrt(1.0 / PSSimEtaHT1_num);
	}
	if (PSSimEtaHT2_num != 0) {
	    PSSimEtaHT2 = PSSimEtaHT2 / PSSimEtaHT2_num;
	    PSSimEtaHT2_num = TMath::Sqrt(1.0 / PSSimEtaHT2_num);
	}
	if (analysisSettingsNew.print) cout << "HighTower enhancement factors (sim eta), " << analysisSettingsNew.title << ": MB = " << PSSimEtaMB << " +/- " << PSSimEtaMB_num << ", HT1 = " << PSSimEtaHT1 << " +/- " << PSSimEtaHT1 << ", HT2 = " << PSSimEtaHT2 << " +/- " << PSSimEtaHT2_num << endl;
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

	selectInv(candidateDataProcessorMB, analysisSettingsNew.settingsMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	selectInv(candidateDataProcessorHT1, analysisSettingsNew.settingsHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	selectInv(candidateDataProcessorHT2, analysisSettingsNew.settingsHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0);
	
	if (!analysisSettingsNew.multiplicityDistributions) {
	    selectInv(candidateDataProcessorMBMix, analysisSettingsNew.settingsMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorHT1Mix, analysisSettingsNew.settingsHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT2Mix, analysisSettingsNew.settingsHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0);

	    selectInv(candidateDataProcessorMBMixRandom, analysisSettingsNew.settingsMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorHT1MixRandom, analysisSettingsNew.settingsHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT2MixRandom, analysisSettingsNew.settingsHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0);

	    selectInv(candidateDataProcessorMBMixBack, analysisSettingsNew.settingsMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorMB : 0);
	    selectInv(candidateDataProcessorHT1MixBack, analysisSettingsNew.settingsHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT1 : 0);
	    selectInv(candidateDataProcessorHT2MixBack, analysisSettingsNew.settingsHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorHT2 : 0);
	}

	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invlist);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invlist);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invlist);
	
	if (!analysisSettingsNew.multiplicityDistributions) {
	    selectInv(candidateDataProcessorEmbedMBMix, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedMB : 0);
	    selectInv(candidateDataProcessorEmbedHT1Mix, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedHT1 : 0);
	    selectInv(candidateDataProcessorEmbedHT2Mix, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedHT2 : 0);

	    selectInv(candidateDataProcessorEmbedMBMixRandom, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedMB : 0);
	    selectInv(candidateDataProcessorEmbedHT1MixRandom, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedHT1 : 0);
	    selectInv(candidateDataProcessorEmbedHT2MixRandom, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedHT2 : 0);

	    selectInv(candidateDataProcessorEmbedMBMixBack, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedMB : 0);
	    selectInv(candidateDataProcessorEmbedHT1MixBack, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedHT1 : 0);
	    selectInv(candidateDataProcessorEmbedHT2MixBack, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedHT2 : 0);
	}

	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB, 5, analysisResults.resultsSimMB.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1, 5, analysisResults.resultsSimHT1.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2, 5, analysisResults.resultsSimHT2.simulatedPtlist);

	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invlist);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invlist);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invlist);
	
	if (!analysisSettingsNew.multiplicityDistributions) {
	    selectInv(candidateDataProcessorEmbedEtaMBMix, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaMB : 0);
	    selectInv(candidateDataProcessorEmbedEtaHT1Mix, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaHT1 : 0);
	    selectInv(candidateDataProcessorEmbedEtaHT2Mix, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invBglist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaHT2 : 0);

	    selectInv(candidateDataProcessorEmbedEtaMBMixRandom, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaMB : 0);
	    selectInv(candidateDataProcessorEmbedEtaHT1MixRandom, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaHT1 : 0);
	    selectInv(candidateDataProcessorEmbedEtaHT2MixRandom, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invBgRandomlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaHT2 : 0);

	    selectInv(candidateDataProcessorEmbedEtaMBMixBack, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaMB : 0);
	    selectInv(candidateDataProcessorEmbedEtaHT1MixBack, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaHT1 : 0);
	    selectInv(candidateDataProcessorEmbedEtaHT2MixBack, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invBgBacklist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorEmbedEtaHT2 : 0);
	}

	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB, 5, analysisResults.resultsSimEtaMB.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1, 5, analysisResults.resultsSimEtaHT1.simulatedPtlist);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2, 5, analysisResults.resultsSimEtaHT2.simulatedPtlist);

	selectInv(candidateDataProcessorEmbedEtabgMB, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.invlist);
	selectInv(candidateDataProcessorEmbedEtabgHT1, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.invlist);
	selectInv(candidateDataProcessorEmbedEtabgHT2, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.invlist);

	selectInv(candidateDataProcessor1gammaMB, analysisSettingsNew.settingsMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMB : 0);
	selectInv(candidateDataProcessor1gammaHT1, analysisSettingsNew.settingsHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1 : 0);
	selectInv(candidateDataProcessor1gammaHT2, analysisSettingsNew.settingsHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2 : 0);

	selectInv(candidateDataProcessor1gammaMBSim, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSim : 0);
	selectInv(candidateDataProcessor1gammaHT1Sim, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1Sim : 0);
	selectInv(candidateDataProcessor1gammaHT2Sim, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2Sim : 0);

	selectInv(candidateDataProcessor1gammaMBSimEta, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSimEta : 0);
	selectInv(candidateDataProcessor1gammaHT1SimEta, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1SimEta : 0);
	selectInv(candidateDataProcessor1gammaHT2SimEta, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2SimEta : 0);

	selectInv(candidateDataProcessor1gammaMBSimEtabg, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaMBSimEtabg : 0);
	selectInv(candidateDataProcessor1gammaHT1SimEtabg, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT1SimEtabg : 0);
	selectInv(candidateDataProcessor1gammaHT2SimEtabg, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.inv1gammalist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessor1gammaHT2SimEtabg : 0);

	selectInv(candidateDataProcessornbarMB, analysisSettingsNew.settingsMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsMB.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarMB : 0);
	selectInv(candidateDataProcessornbarHT1, analysisSettingsNew.settingsHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT1.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT1 : 0);
	selectInv(candidateDataProcessornbarHT2, analysisSettingsNew.settingsHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsHT2.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT2 : 0);

	selectInv(candidateDataProcessornbarMBSim, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarMBSim : 0);
	selectInv(candidateDataProcessornbarHT1Sim, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT1Sim : 0);
	selectInv(candidateDataProcessornbarHT2Sim, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT2Sim : 0);

	selectInv(candidateDataProcessornbarMBSimEta, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaMB.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarMBSimEta : 0);
	selectInv(candidateDataProcessornbarHT1SimEta, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT1.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT1SimEta : 0);
	selectInv(candidateDataProcessornbarHT2SimEta, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtaHT2.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT2SimEta : 0);

	selectInv(candidateDataProcessornbarMBSimEtabg, analysisSettingsNew.settingsSimEtaMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgMB.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarMBSimEtabg : 0);
	selectInv(candidateDataProcessornbarHT1SimEtabg, analysisSettingsNew.settingsSimEtaHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT1.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT1SimEtabg : 0);
	selectInv(candidateDataProcessornbarHT2SimEtabg, analysisSettingsNew.settingsSimEtaHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimEtabgHT2.invnbarlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessornbarHT2SimEtabg : 0);

	selectInv(candidateDataProcessorMB, analysisSettingsNew.settingsMB, 1, analysisResults.resultsMB.multlist, 0);
	selectInv(candidateDataProcessorHT1, analysisSettingsNew.settingsHT1, 1, analysisResults.resultsHT1.multlist, 0);
	selectInv(candidateDataProcessorHT2, analysisSettingsNew.settingsHT2, 1, analysisResults.resultsHT2.multlist, 0);

	selectInv(candidateDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB, 1, analysisResults.resultsSimMB.multlist, 0);
	selectInv(candidateDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1, 1, analysisResults.resultsSimHT1.multlist, 0);
	selectInv(candidateDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2, 1, analysisResults.resultsSimHT2.multlist, 0);

	selectInv(candidateDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimEtaMB, 1, analysisResults.resultsSimEtaMB.multlist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimEtaHT1, 1, analysisResults.resultsSimEtaHT1.multlist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimEtaHT2, 1, analysisResults.resultsSimEtaHT2.multlist, 0);
	
	selectInv(candidateDataProcessorMBNocpv, analysisSettingsNew.settingsMB, 4, analysisResults.resultsMB.trackdistlist, 0);
	selectInv(candidateDataProcessorHT1Nocpv, analysisSettingsNew.settingsHT1, 4, analysisResults.resultsHT1.trackdistlist, 0);
	selectInv(candidateDataProcessorHT2Nocpv, analysisSettingsNew.settingsHT2, 4, analysisResults.resultsHT2.trackdistlist, 0);
	
	selectInv(candidateDataProcessorEmbedMBNocpv, analysisSettingsNew.settingsSimMB, 4, analysisResults.resultsSimMB.trackdistlist, 0);
	selectInv(candidateDataProcessorEmbedHT1Nocpv, analysisSettingsNew.settingsSimHT1, 4, analysisResults.resultsSimHT1.trackdistlist, 0);
	selectInv(candidateDataProcessorEmbedHT2Nocpv, analysisSettingsNew.settingsSimHT2, 4, analysisResults.resultsSimHT2.trackdistlist, 0);
	
	selectInv(candidateDataProcessorEmbedEtaMBNocpv, analysisSettingsNew.settingsSimEtaMB, 4, analysisResults.resultsSimEtaMB.trackdistlist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT1Nocpv, analysisSettingsNew.settingsSimEtaHT1, 4, analysisResults.resultsSimEtaHT1.trackdistlist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT2Nocpv, analysisSettingsNew.settingsSimEtaHT2, 4, analysisResults.resultsSimEtaHT2.trackdistlist, 0);
	/*
	selectInv(candidateDataProcessorMBNotcpv, analysisSettingsNew.settingsMB, 6, analysisResults.resultsMB.trackdist2list, 0);
	selectInv(candidateDataProcessorHT1Notcpv, analysisSettingsNew.settingsHT1, 6, analysisResults.resultsHT1.trackdist2list, 0);
	selectInv(candidateDataProcessorHT2Notcpv, analysisSettingsNew.settingsHT2, 6, analysisResults.resultsHT2.trackdist2list, 0);
	
	selectInv(candidateDataProcessorEmbedMBNotcpv, analysisSettingsNew.settingsSimMB, 6, analysisResults.resultsSimMB.trackdist2list, 0);
	selectInv(candidateDataProcessorEmbedHT1Notcpv, analysisSettingsNew.settingsSimHT1, 6, analysisResults.resultsSimHT1.trackdist2list, 0);
	selectInv(candidateDataProcessorEmbedHT2Notcpv, analysisSettingsNew.settingsSimHT2, 6, analysisResults.resultsSimHT2.trackdist2list, 0);
	
	selectInv(candidateDataProcessorEmbedEtaMBNotcpv, analysisSettingsNew.settingsSimEtaMB, 6, analysisResults.resultsSimEtaMB.trackdist2list, 0);
	selectInv(candidateDataProcessorEmbedEtaHT1Notcpv, analysisSettingsNew.settingsSimEtaHT1, 6, analysisResults.resultsSimEtaHT1.trackdist2list, 0);
	selectInv(candidateDataProcessorEmbedEtaHT2Notcpv, analysisSettingsNew.settingsSimEtaHT2, 6, analysisResults.resultsSimEtaHT2.trackdist2list, 0);
	*/
	/*
	selectInv(candidateDataProcessorMBMix, analysisSettingsNew.settingsMB, 4, analysisResults.resultsMB.trackdistBglist, 0);
	selectInv(candidateDataProcessorHT1Mix, analysisSettingsNew.settingsHT1, 4, analysisResults.resultsHT1.trackdistBglist, 0);
	selectInv(candidateDataProcessorHT2Mix, analysisSettingsNew.settingsHT2, 4, analysisResults.resultsHT2.trackdistBglist, 0);

	selectInv(candidateDataProcessorEmbedMBMix, analysisSettingsNew.settingsSimMB, 4, analysisResults.resultsSimMB.trackdistBglist, 0);
	selectInv(candidateDataProcessorEmbedHT1Mix, analysisSettingsNew.settingsSimHT1, 4, analysisResults.resultsSimHT1.trackdistBglist, 0);
	selectInv(candidateDataProcessorEmbedHT2Mix, analysisSettingsNew.settingsSimHT2, 4, analysisResults.resultsSimHT2.trackdistBglist, 0);
	
	selectInv(candidateDataProcessorEmbedEtaMBMix, analysisSettingsNew.settingsSimEtaMB, 4, analysisResults.resultsSimEtaMB.trackdistBglist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT1Mix, analysisSettingsNew.settingsSimEtaHT1, 4, analysisResults.resultsSimEtaHT1.trackdistBglist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT2Mix, analysisSettingsNew.settingsSimEtaHT2, 4, analysisResults.resultsSimEtaHT2.trackdistBglist, 0);
	*/
	/*
	selectInv(candidateDataProcessorMBMixRandom, analysisSettingsNew.settingsMB, 4, analysisResults.resultsMB.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorHT1MixRandom, analysisSettingsNew.settingsHT1, 4, analysisResults.resultsHT1.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorHT2MixRandom, analysisSettingsNew.settingsHT2, 4, analysisResults.resultsHT2.trackdistBgRandomlist, 0);
	
	selectInv(candidateDataProcessorEmbedMBMixRandom, analysisSettingsNew.settingsSimMB, 4, analysisResults.resultsSimMB.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorEmbedHT1MixRandom, analysisSettingsNew.settingsSimHT1, 4, analysisResults.resultsSimHT1.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorEmbedHT2MixRandom, analysisSettingsNew.settingsSimHT2, 4, analysisResults.resultsSimHT2.trackdistBgRandomlist, 0);
	
	selectInv(candidateDataProcessorEmbedEtaMBMixRandom, analysisSettingsNew.settingsSimEtaMB, 4, analysisResults.resultsSimEtaMB.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT1MixRandom, analysisSettingsNew.settingsSimEtaHT1, 4, analysisResults.resultsSimEtaHT1.trackdistBgRandomlist, 0);
	selectInv(candidateDataProcessorEmbedEtaHT2MixRandom, analysisSettingsNew.settingsSimEtaHT2, 4, analysisResults.resultsSimEtaHT2.trackdistBgRandomlist, 0);
	*/
	selectInv(candidateDataProcessorPtshiftMB, analysisSettingsNew.settingsSimMB, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimMB.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftMB : 0);
	selectInv(candidateDataProcessorPtshiftHT1, analysisSettingsNew.settingsSimHT1, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT1.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftHT1 : 0);
	selectInv(candidateDataProcessorPtshiftHT2, analysisSettingsNew.settingsSimHT2, analysisSettingsNew.multiplicityDistributions, analysisResults.resultsSimHT2.invPtshiftlist, analysisSettingsNew.pointMultiplicityDistributions ? pointDataProcessorPtshiftHT2 : 0);
	
	selectInvSim(simuDataProcessorMB, analysisSettingsNew.settingsSimMB, analysisResults.areaSimuMB);
	selectInvSim(simuDataProcessorHT1, analysisSettingsNew.settingsSimHT1, analysisResults.areaSimuHT1);
	selectInvSim(simuDataProcessorHT2, analysisSettingsNew.settingsSimHT2, analysisResults.areaSimuHT2);

	selectInvSim(simuDataProcessorEtaMB, analysisSettingsNew.settingsSimEtaMB, analysisResults.areaSimuEtaMB);
	selectInvSim(simuDataProcessorEtaHT1, analysisSettingsNew.settingsSimEtaHT1, analysisResults.areaSimuEtaHT1);
	selectInvSim(simuDataProcessorEtaHT2, analysisSettingsNew.settingsSimEtaHT2, analysisResults.areaSimuEtaHT2);

	selectInvSim(simuDataProcessorEtabgMB, analysisSettingsNew.settingsSimEtaMB, analysisResults.areaSimuEtabgMB);
	selectInvSim(simuDataProcessorEtabgHT1, analysisSettingsNew.settingsSimEtaHT1, analysisResults.areaSimuEtabgHT1);
	selectInvSim(simuDataProcessorEtabgHT2, analysisSettingsNew.settingsSimEtaHT2, analysisResults.areaSimuEtabgHT2);

	selectInvPoint(pointDataProcessorMB, analysisSettingsNew.settingsMB, analysisResults.resultsMB.rawYieldGamma);
	selectInvPoint(pointDataProcessorHT1, analysisSettingsNew.settingsHT1, analysisResults.resultsHT1.rawYieldGamma);
	selectInvPoint(pointDataProcessorHT2, analysisSettingsNew.settingsHT2, analysisResults.resultsHT2.rawYieldGamma);

	selectInvPoint(pointDataProcessorEmbedMB, analysisSettingsNew.settingsSimMB, analysisResults.resultsSimMB.rawYieldGamma);
	selectInvPoint(pointDataProcessorEmbedHT1, analysisSettingsNew.settingsSimHT1, analysisResults.resultsSimHT1.rawYieldGamma);
	selectInvPoint(pointDataProcessorEmbedHT2, analysisSettingsNew.settingsSimHT2, analysisResults.resultsSimHT2.rawYieldGamma);

	selectInvMCGamma(gammaPi0DataProcessorMB, analysisSettingsNew.settingsSimMB, analysisResults.resultsSimMB.inputYieldGamma);
	selectInvMCGamma(gammaPi0DataProcessorHT1, analysisSettingsNew.settingsSimHT1, analysisResults.resultsSimHT1.inputYieldGamma);
	selectInvMCGamma(gammaPi0DataProcessorHT2, analysisSettingsNew.settingsSimHT2, analysisResults.resultsSimHT2.inputYieldGamma);

	selectInvPoint(pointDataProcessorEmbedEtaMB, analysisSettingsNew.settingsSimMB, analysisResults.resultsSimEtaMB.rawYieldGamma);
	selectInvPoint(pointDataProcessorEmbedEtaHT1, analysisSettingsNew.settingsSimHT1, analysisResults.resultsSimEtaHT1.rawYieldGamma);
	selectInvPoint(pointDataProcessorEmbedEtaHT2, analysisSettingsNew.settingsSimHT2, analysisResults.resultsSimEtaHT2.rawYieldGamma);

	selectInvMCGamma(gammaEtaDataProcessorMB, analysisSettingsNew.settingsSimEtaMB, analysisResults.resultsSimEtaMB.inputYieldGamma);
	selectInvMCGamma(gammaEtaDataProcessorHT1, analysisSettingsNew.settingsSimEtaHT1, analysisResults.resultsSimEtaHT1.inputYieldGamma);
	selectInvMCGamma(gammaEtaDataProcessorHT2, analysisSettingsNew.settingsSimEtaHT2, analysisResults.resultsSimEtaHT2.inputYieldGamma);

	selectInvPoint(pointDataProcessor1gammaMB, analysisSettingsNew.settingsSimMB, analysisResults.rawYieldGammaMB);
	selectInvPoint(pointDataProcessor1gammaHT1, analysisSettingsNew.settingsSimHT1, analysisResults.rawYieldGammaHT1);
	selectInvPoint(pointDataProcessor1gammaHT2, analysisSettingsNew.settingsSimHT2, analysisResults.rawYieldGammaHT2);

	selectInvMCGamma(gammaDataProcessorMB, analysisSettingsNew.settingsSimMB, analysisResults.inputYieldGammaMB);
	selectInvMCGamma(gammaDataProcessorHT1, analysisSettingsNew.settingsSimHT1, analysisResults.inputYieldGammaHT1);
	selectInvMCGamma(gammaDataProcessorHT2, analysisSettingsNew.settingsSimHT2, analysisResults.inputYieldGammaHT2);

	selectInvPoint(pointDataProcessornbarMB, analysisSettingsNew.settingsSimMB, analysisResults.rawYieldNbarMB);
	selectInvPoint(pointDataProcessornbarHT1, analysisSettingsNew.settingsSimHT1, analysisResults.rawYieldNbarHT1);
	selectInvPoint(pointDataProcessornbarHT2, analysisSettingsNew.settingsSimHT2, analysisResults.rawYieldNbarHT2);

	selectInvMCGamma(nbarDataProcessorMB, analysisSettingsNew.settingsSimMB, analysisResults.inputYieldNbarMB);
	selectInvMCGamma(nbarDataProcessorHT1, analysisSettingsNew.settingsSimHT1, analysisResults.inputYieldNbarHT1);
	selectInvMCGamma(nbarDataProcessorHT2, analysisSettingsNew.settingsSimHT2, analysisResults.inputYieldNbarHT2);

	if (pointDataProcessorMB && pointDataProcessor1gammaMB) analysisResults.resultsMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleMB", analysisSettingsNew.title + ", low-mass norm MB"
	    , pointDataProcessorMB->getPt(), pointDataProcessor1gammaMB->getPt()
	    , analysisSettingsNew.settingsMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsMB.lowmassbgPointsHighPt
	    , analysisResults.resultsMB.invlist, analysisResults.resultsMB.lowNormPointsPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	if (pointDataProcessorHT1 && pointDataProcessor1gammaHT1) analysisResults.resultsHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleHT1", analysisSettingsNew.title + ", low-mass norm HT1"
	    , pointDataProcessorHT1->getPt(), pointDataProcessor1gammaHT1->getPt()
	    , analysisSettingsNew.settingsHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsHT1.invlist, analysisResults.resultsHT1.lowNormPointsPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	if (pointDataProcessorHT2 && pointDataProcessor1gammaHT2) analysisResults.resultsHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleHT2", analysisSettingsNew.title + ", low-mass norm HT2"
	    , pointDataProcessorHT2->getPt(), pointDataProcessor1gammaHT2->getPt()
	    , analysisSettingsNew.settingsHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsHT2.invlist, analysisResults.resultsHT2.lowNormPointsPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	if (pointDataProcessorEmbedMB && pointDataProcessor1gammaMBSim) analysisResults.resultsSimMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimMB", analysisSettingsNew.title + ", low-mass norm MB MC"
	    , pointDataProcessorEmbedMB->getPt(), pointDataProcessor1gammaMBSim->getPt()
	    , analysisSettingsNew.settingsSimMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimMB.lowmassbgPointsHighPt
	    , analysisResults.resultsSimMB.invlist, analysisResults.resultsSimMB.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedHT1 && pointDataProcessor1gammaHT1Sim) analysisResults.resultsSimHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimHT1", analysisSettingsNew.title + ", low-mass norm HT1 MC"
	    , pointDataProcessorEmbedHT1->getPt(), pointDataProcessor1gammaHT1Sim->getPt()
	    , analysisSettingsNew.settingsSimHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsSimHT1.invlist, analysisResults.resultsSimHT1.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedHT2 && pointDataProcessor1gammaHT2Sim) analysisResults.resultsSimHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimHT2", analysisSettingsNew.title + ", low-mass norm HT2 MC"
	    , pointDataProcessorEmbedHT2->getPt(), pointDataProcessor1gammaHT2Sim->getPt()
	    , analysisSettingsNew.settingsSimHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsSimHT2.invlist, analysisResults.resultsSimHT2.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtaMB && pointDataProcessor1gammaMBSimEta) analysisResults.resultsSimEtaMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtaMB", analysisSettingsNew.title + ", low-mass norm MB eta MC"
	    , pointDataProcessorEmbedEtaMB->getPt(), pointDataProcessor1gammaMBSimEta->getPt()
	    , analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtaMB.invlist, analysisResults.resultsSimEtaMB.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtaHT1 && pointDataProcessor1gammaHT1SimEta) analysisResults.resultsSimEtaHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtaHT1", analysisSettingsNew.title + ", low-mass norm HT1 eta MC"
	    , pointDataProcessorEmbedEtaHT1->getPt(), pointDataProcessor1gammaHT1SimEta->getPt()
	    , analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtaHT1.invlist, analysisResults.resultsSimEtaHT1.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtaHT2 && pointDataProcessor1gammaHT2SimEta) analysisResults.resultsSimEtaHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtaHT2", analysisSettingsNew.title + ", low-mass norm HT2 eta MC"
	    , pointDataProcessorEmbedEtaHT2->getPt(), pointDataProcessor1gammaHT2SimEta->getPt()
	    , analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtaHT2.invlist, analysisResults.resultsSimEtaHT2.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtabgMB && pointDataProcessor1gammaMBSimEtabg) analysisResults.resultsSimEtabgMB.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtabgMB", analysisSettingsNew.title + ", low-mass norm MB etabg MC"
	    , pointDataProcessorEmbedEtabgMB->getPt(), pointDataProcessor1gammaMBSimEtabg->getPt()
	    , analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtabgMB.invlist, analysisResults.resultsSimEtabgMB.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtabgHT1 && pointDataProcessor1gammaHT1SimEtabg) analysisResults.resultsSimEtabgHT1.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtabgHT1", analysisSettingsNew.title + ", low-mass norm HT1 etabg MC"
	    , pointDataProcessorEmbedEtabgHT1->getPt(), pointDataProcessor1gammaHT1SimEtabg->getPt()
	    , analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtabgHT1.invlist, analysisResults.resultsSimEtabgHT1.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtabgHT2 && pointDataProcessor1gammaHT2SimEtabg) analysisResults.resultsSimEtabgHT2.lowNormPoints = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscaleSimEtabgHT2", analysisSettingsNew.title + ", low-mass norm HT2 etabg MC"
	    , pointDataProcessorEmbedEtabgHT2->getPt(), pointDataProcessor1gammaHT2SimEtabg->getPt()
	    , analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtabgHT2.invlist, analysisResults.resultsSimEtabgHT2.lowNormPointsPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);

	if (pointDataProcessorMB && pointDataProcessornbarMB) analysisResults.resultsMB.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarMB", analysisSettingsNew.title + ", low-mass norm nbar MB"
	    , pointDataProcessorMB->getPt(), pointDataProcessornbarMB->getPt()
	    , analysisSettingsNew.settingsMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsMB.lowmassbgPointsHighPt
	    , analysisResults.resultsMB.invlist, analysisResults.resultsMB.lowNormPointsNbarPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	if (pointDataProcessorHT1 && pointDataProcessornbarHT1) analysisResults.resultsHT1.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarHT1", analysisSettingsNew.title + ", low-mass norm nbar HT1"
	    , pointDataProcessorHT1->getPt(), pointDataProcessornbarHT1->getPt()
	    , analysisSettingsNew.settingsHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsHT1.invlist, analysisResults.resultsHT1.lowNormPointsNbarPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	if (pointDataProcessorHT2 && pointDataProcessornbarHT2) analysisResults.resultsHT2.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarHT2", analysisSettingsNew.title + ", low-mass norm nbar HT2"
	    , pointDataProcessorHT2->getPt(), pointDataProcessornbarHT2->getPt()
	    , analysisSettingsNew.settingsHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsHT2.invlist, analysisResults.resultsHT2.lowNormPointsNbarPtCorr, analysisSettingsNew.showReal, analysisSettingsNew.print);
	if (pointDataProcessorEmbedMB && pointDataProcessornbarMBSim) analysisResults.resultsSimMB.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimMB", analysisSettingsNew.title + ", low-mass norm nbar MB MC"
	    , pointDataProcessorEmbedMB->getPt(), pointDataProcessornbarMBSim->getPt()
	    , analysisSettingsNew.settingsSimMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimMB.lowmassbgPointsHighPt
	    , analysisResults.resultsSimMB.invlist, analysisResults.resultsSimMB.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedHT1 && pointDataProcessornbarHT1Sim) analysisResults.resultsSimHT1.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimHT1", analysisSettingsNew.title + ", low-mass norm nbar HT1 MC"
	    , pointDataProcessorEmbedHT1->getPt(), pointDataProcessornbarHT1Sim->getPt()
	    , analysisSettingsNew.settingsSimHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsSimHT1.invlist, analysisResults.resultsSimHT1.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedHT2 && pointDataProcessornbarHT2Sim) analysisResults.resultsSimHT2.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimHT2", analysisSettingsNew.title + ", low-mass norm nbar HT2 MC"
	    , pointDataProcessorEmbedHT2->getPt(), pointDataProcessornbarHT2Sim->getPt()
	    , analysisSettingsNew.settingsSimHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsSimHT2.invlist, analysisResults.resultsSimHT2.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtaMB && pointDataProcessornbarMBSimEta) analysisResults.resultsSimEtaMB.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimEtaMB", analysisSettingsNew.title + ", low-mass norm nbar MB eta MC"
	    , pointDataProcessorEmbedEtaMB->getPt(), pointDataProcessornbarMBSimEta->getPt()
	    , analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtaMB.invlist, analysisResults.resultsSimEtaMB.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtaHT1 && pointDataProcessornbarHT1SimEta) analysisResults.resultsSimEtaHT1.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimEtaHT1", analysisSettingsNew.title + ", low-mass norm nbar HT1 eta MC"
	    , pointDataProcessorEmbedEtaHT1->getPt(), pointDataProcessornbarHT1SimEta->getPt()
	    , analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtaHT1.invlist, analysisResults.resultsSimEtaHT1.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtaHT2 && pointDataProcessornbarHT2SimEta) analysisResults.resultsSimEtaHT2.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimEtaHT2", analysisSettingsNew.title + ", low-mass norm nbar HT2 eta MC"
	    , pointDataProcessorEmbedEtaHT2->getPt(), pointDataProcessornbarHT2SimEta->getPt()
	    , analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtaHT2.invlist, analysisResults.resultsSimEtaHT2.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtabgMB && pointDataProcessornbarMBSimEtabg) analysisResults.resultsSimEtabgMB.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimEtabgMB", analysisSettingsNew.title + ", low-mass norm nbar MB etabg MC"
	    , pointDataProcessorEmbedEtabgMB->getPt(), pointDataProcessornbarMBSimEtabg->getPt()
	    , analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaMB.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtabgMB.invlist, analysisResults.resultsSimEtabgMB.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtabgHT1 && pointDataProcessornbarHT1SimEtabg) analysisResults.resultsSimEtabgHT1.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimEtabgHT1", analysisSettingsNew.title + ", low-mass norm nbar HT1 etabg MC"
	    , pointDataProcessorEmbedEtabgHT1->getPt(), pointDataProcessornbarHT1SimEtabg->getPt()
	    , analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT1.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtabgHT1.invlist, analysisResults.resultsSimEtabgHT1.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);
	if (pointDataProcessorEmbedEtabgHT2 && pointDataProcessornbarHT2SimEtabg) analysisResults.resultsSimEtabgHT2.lowNormPointsNbar = getScaleFactorFromPt(analysisSettingsNew.name + "_lowscalenbarSimEtabgHT2", analysisSettingsNew.title + ", low-mass norm nbar HT2 etabg MC"
	    , pointDataProcessorEmbedEtabgHT2->getPt(), pointDataProcessornbarHT2SimEtabg->getPt()
	    , analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimEtaHT2.lowmassbgPointsHighPt
	    , analysisResults.resultsSimEtabgHT2.invlist, analysisResults.resultsSimEtabgHT2.lowNormPointsNbarPtCorr, analysisSettingsNew.showSim, analysisSettingsNew.print);


	{
	        TF1 *shiftMBfunc = 0;
	        TF1 *shiftHT1func = 0;
	        TF1 *shiftHT2func = 0;
		TString shiftMBname = analysisSettingsNew.name + "_shiftMB";
		TString shiftMBfuncname = analysisSettingsNew.name + "_shiftMBfunc";
		TString shiftMBtitle = analysisSettingsNew.title + ", p_{T} shift MB";
		calculateBinPurity(candidateDataProcessorEmbedMB
				, analysisResults.resultsSimMB.invlist
				, shiftMBfuncname.Data(), analysisSettingsNew.shiftFuncStr
				, &shiftMBfunc
				, shiftMBname.Data(), shiftMBtitle.Data()
				, analysisSettingsNew.showSim
				);
		TString shiftHT1name = analysisSettingsNew.name + "_shiftHT1";
		TString shiftHT1funcname = analysisSettingsNew.name + "_shiftHT1func";
		TString shiftHT1title = analysisSettingsNew.title + ", p_{T} shift HT1";
		calculateBinPurity(candidateDataProcessorEmbedHT1
				, analysisResults.resultsSimHT1.invlist
				, shiftHT1funcname.Data(), analysisSettingsNew.shiftFuncStr
				, &shiftHT1func
				, shiftHT1name.Data(), shiftHT1title.Data()
				, analysisSettingsNew.showSim
				);
		TString shiftHT2name = analysisSettingsNew.name + "_shiftHT2";
		TString shiftHT2funcname = analysisSettingsNew.name + "_shiftHT2func";
		TString shiftHT2title = analysisSettingsNew.title + ", p_{T} shift HT2";
		calculateBinPurity(candidateDataProcessorEmbedHT2
				, analysisResults.resultsSimHT2.invlist
				, shiftHT2funcname.Data(), analysisSettingsNew.shiftFuncStr
				, &shiftHT2func
				, shiftHT2name.Data(), shiftHT2title.Data()
				, analysisSettingsNew.showSim
				);
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

	{
	const TH2F *h1 = candidateDataProcessorMB ? candidateDataProcessorMB->getTpcRefmultTrackDist() : 0;
	const TH2F *h2 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getTpcRefmultTrackDist() : 0;
	const TH2F *h3 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getTpcRefmultTrackDist() : 0;
	const TH2F *h = h1 ? h1 : (h2 ? h2 : h3);
	TH2F *hist = h ? (new TH2F(*h)) : 0;
	if (hist && (h != h1)) hist->Add(h1);
	if (hist && (h != h2)) hist->Add(h2);
	if (hist && (h != h3)) hist->Add(h3);
	fitCpvMult(hist, analysisSettingsNew.settingsMB.cpvCorrectionMult, analysisResults.resultsMB.cpvCorrectionMultFunc, analysisSettingsNew.settingsMB.name);
	fitCpvMult(hist, analysisSettingsNew.settingsHT1.cpvCorrectionMult, analysisResults.resultsHT1.cpvCorrectionMultFunc, analysisSettingsNew.settingsHT1.name);
	fitCpvMult(hist, analysisSettingsNew.settingsHT2.cpvCorrectionMult, analysisResults.resultsHT2.cpvCorrectionMultFunc, analysisSettingsNew.settingsHT2.name);
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
		    ratioMBfunc = new TF1("ratioMBfunc", (const Char_t*)analysisSettingsNew.tpcVertexEffFunc);
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
		    ratioHT1func = new TF1("ratioHT1func", (const Char_t*)analysisSettingsNew.tpcVertexEffFunc);
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
		    ratioHT2func = new TF1("ratioHT2func", (const Char_t*)analysisSettingsNew.tpcVertexEffFunc);
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
		    ratioMBfunc = new TF1("ratioMBfunc", (const Char_t*)analysisSettingsNew.bbcVertexEffFunc);
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
		    ratioHT1func = new TF1("ratioHT1func", (const Char_t*)analysisSettingsNew.bbcVertexEffFunc);
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
		    ratioHT2func = new TF1("ratioHT2func", (const Char_t*)analysisSettingsNew.bbcVertexEffFunc);
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
		const TH1F *histInJetYesMB = candidateDataProcessorMB ? candidateDataProcessorMB->getInJetPointYes() : 0;
		if (histInJetYesMB) histJetFracMB = new TH1F(*histInJetYesMB);
		if (histJetFracMB && histPtMB) histJetFracMB->Divide(histPtMB);
		if (histJetFracMB && histPtMB && histJetYesMB) {
		    histJetFracMB->Multiply(histJetYesMB);
		    histJetFracMB->Divide(histPtMB);
		}
	}
	{
		const TH1F *histPtHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getPt() : 0;
		const TH1F *histJetYesHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getJetYes() : 0;
		const TH1F *histInJetYesHT1 = candidateDataProcessorHT1 ? candidateDataProcessorHT1->getInJetPointYes() : 0;
		if (histInJetYesHT1) histJetFracHT1 = new TH1F(*histInJetYesHT1);
		if (histJetFracHT1 && histPtHT1) histJetFracHT1->Divide(histPtHT1);
		if (histJetFracHT1 && histPtHT1 && histJetYesHT1) {
			histJetFracHT1->Multiply(histJetYesHT1);
			histJetFracHT1->Divide(histPtHT1);
		}
	}
	{
		const TH1F *histPtHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getPt() : 0;
		const TH1F *histJetYesHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getJetYes() : 0;
		const TH1F *histInJetYesHT2 = candidateDataProcessorHT2 ? candidateDataProcessorHT2->getInJetPointYes() : 0;
		if (histInJetYesHT2) histJetFracHT2 = new TH1F(*histInJetYesHT2);
		if (histJetFracHT2 && histPtHT2) histJetFracHT2->Divide(histPtHT2);
		if (histJetFracHT2 && histPtHT2 && histJetYesHT2) {
		    histJetFracHT2->Multiply(histJetYesHT2);
		    histJetFracHT2->Divide(histPtHT2);
		}
	}
	if (candidateDataProcessorMB) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandMB", analysisSettingsNew.title + ", in jet cand fraction MB"
		, candidateDataProcessorMB->getInJetPointYes(), candidateDataProcessorMB->getPt(), 0, 100
		, analysisResults.resultsMB.invlist, temp1, analysisSettingsNew.showReal && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvMB", analysisSettingsNew.title + ", jet ev fraction MB"
		, candidateDataProcessorMB->getPt(), candidateDataProcessorMB->getJetYes(), 0, 100
		, analysisResults.resultsMB.invlist, temp2, analysisSettingsNew.showReal && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsMB.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	if (candidateDataProcessorHT1) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandHT1", analysisSettingsNew.title + ", in jet cand fraction HT1"
		, candidateDataProcessorHT1->getInJetPointYes(), candidateDataProcessorHT1->getPt(), 0, 100
		, analysisResults.resultsHT1.invlist, temp1, analysisSettingsNew.showReal && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvHT1", analysisSettingsNew.title + ", jet ev fraction HT1"
		, candidateDataProcessorHT1->getPt(), candidateDataProcessorHT1->getJetYes(), 0, 100
		, analysisResults.resultsHT1.invlist, temp2, analysisSettingsNew.showReal && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsHT1.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	if (candidateDataProcessorHT2) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandHT2", analysisSettingsNew.title + ", in jet cand fraction HT2"
		, candidateDataProcessorHT2->getInJetPointYes(), candidateDataProcessorHT2->getPt(), 0, 100
		, analysisResults.resultsHT2.invlist, temp1, analysisSettingsNew.showReal && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvHT2", analysisSettingsNew.title + ", jet ev fraction HT2"
		, candidateDataProcessorHT2->getPt(), candidateDataProcessorHT2->getJetYes(), 0, 100
		, analysisResults.resultsHT2.invlist, temp2, analysisSettingsNew.showReal && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsHT2.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	list<bin_stat_list_type> bgJetToTotalFromCandidatesCorrList;
	bgJetToTotalFromCandidatesCorrList.push_back(analysisResults.resultsMB.bgJetToTotalFromCandidatesCorr);
	bgJetToTotalFromCandidatesCorrList.push_back(analysisResults.resultsHT1.bgJetToTotalFromCandidatesCorr);
	bgJetToTotalFromCandidatesCorrList.push_back(analysisResults.resultsHT2.bgJetToTotalFromCandidatesCorr);
	bin_stat_list_type bgJetToTotalFromCandidatesCorr, temp;
	mergePoints(bgJetToTotalFromCandidatesCorrList, bgJetToTotalFromCandidatesCorr, temp);
	if ((analysisSettingsNew.settingsMB.fitJetBgFractionFromCandidates && histJetFracMB) || analysisSettingsNew.settingsMB.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsMB.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsMB.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracMB, analysisSettingsNew.settingsMB.jetBgFraction, analysisResults.resultsMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandMB");
	    } else if (analysisSettingsNew.settingsMB.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorr, 0, analysisSettingsNew.settingsMB.jetBgFraction, analysisResults.resultsMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandMB");
	    }
	    if (analysisSettingsNew.settingsMB.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsMB.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsMB.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsMB.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsMB.jetBgFraction.wasFit = false;
	}
	if ((analysisSettingsNew.settingsHT1.fitJetBgFractionFromCandidates && histJetFracHT1) || analysisSettingsNew.settingsHT1.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsHT1.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsHT1.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracHT1, analysisSettingsNew.settingsHT1.jetBgFraction, analysisResults.resultsHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandHT1");
	    } else if (analysisSettingsNew.settingsHT1.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorr, 0, analysisSettingsNew.settingsHT1.jetBgFraction, analysisResults.resultsHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandHT1");
	    }
	    if (analysisSettingsNew.settingsHT1.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsHT1.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsHT1.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsHT1.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsHT1.jetBgFraction.wasFit = false;
	}
	if ((analysisSettingsNew.settingsHT2.fitJetBgFractionFromCandidates && histJetFracHT2) || analysisSettingsNew.settingsHT2.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsHT2.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsHT2.fitJetBgFractionFromCandidates) { 
		fitParameter(0, histJetFracHT2, analysisSettingsNew.settingsHT2.jetBgFraction, analysisResults.resultsHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandHT2");
	    } else if (analysisSettingsNew.settingsHT2.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorr, 0, analysisSettingsNew.settingsHT2.jetBgFraction, analysisResults.resultsHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandHT2");
	    }
	    if (analysisSettingsNew.settingsHT2.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsHT2.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsHT2.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsHT2.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsHT2.jetBgFraction.wasFit = false;
	}

	TH1F *histJetFracSimMB = 0;
	TH1F *histJetFracSimHT1 = 0;
	TH1F *histJetFracSimHT2 = 0;
	{
		const TH1F *histPtSimMB = candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getPt() : 0;
		const TH1F *histJetYesSimMB = candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getJetYes() : 0;
		const TH1F *histInJetYesSimMB = candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getInJetPointYes() : 0;
		if (histInJetYesSimMB) histJetFracSimMB = new TH1F(*histInJetYesSimMB);
		if (histJetFracSimMB && histPtSimMB) histJetFracSimMB->Divide(histPtSimMB);
		if (histJetFracSimMB && histPtSimMB && histJetYesSimMB) {
		    histJetFracSimMB->Multiply(histJetYesSimMB);
		    histJetFracSimMB->Divide(histPtSimMB);
		}
	}
	{
		const TH1F *histPtSimHT1 = candidateDataProcessorEmbedHT1 ? candidateDataProcessorEmbedHT1->getPt() : 0;
		const TH1F *histJetYesSimHT1 = candidateDataProcessorEmbedHT1 ? candidateDataProcessorEmbedHT1->getJetYes() : 0;
		const TH1F *histInJetYesSimHT1 = candidateDataProcessorEmbedHT1 ? candidateDataProcessorEmbedHT1->getInJetPointYes() : 0;
		if (histInJetYesSimHT1) histJetFracSimHT1 = new TH1F(*histInJetYesSimHT1);
		if (histJetFracSimHT1 && histPtSimHT1) histJetFracSimHT1->Divide(histPtSimHT1);
		if (histJetFracSimHT1 && histPtSimHT1 && histJetYesSimHT1) {
		    histJetFracSimHT1->Multiply(histJetYesSimHT1);
		    histJetFracSimHT1->Divide(histPtSimHT1);
		}
	}
	{
		const TH1F *histPtSimHT2 = candidateDataProcessorEmbedHT2 ? candidateDataProcessorEmbedHT2->getPt() : 0;
		const TH1F *histJetYesSimHT2 = candidateDataProcessorEmbedHT2 ? candidateDataProcessorEmbedHT2->getJetYes() : 0;
		const TH1F *histInJetYesSimHT2 = candidateDataProcessorEmbedHT2 ? candidateDataProcessorEmbedHT2->getInJetPointYes() : 0;
		if (histInJetYesSimHT2) histJetFracSimHT2 = new TH1F(*histInJetYesSimHT2);
		if (histJetFracSimHT2 && histPtSimHT2) histJetFracSimHT2->Divide(histPtSimHT2);
		if (histJetFracSimHT2 && histPtSimHT2 && histJetYesSimHT2) {
		    histJetFracSimHT2->Multiply(histJetYesSimHT2);
		    histJetFracSimHT2->Divide(histPtSimHT2);
		}
	}
	if (candidateDataProcessorEmbedMB) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandSimMB", analysisSettingsNew.title + ", in jet cand fraction MB MC"
		, candidateDataProcessorEmbedMB->getInJetPointYes(), candidateDataProcessorEmbedMB->getPt(), 0, 100
		, analysisResults.resultsSimMB.invlist, temp1, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvSimMB", analysisSettingsNew.title + ", jet ev fraction MB MC"
		, candidateDataProcessorEmbedMB->getPt(), candidateDataProcessorEmbedMB->getJetYes(), 0, 100
		, analysisResults.resultsSimMB.invlist, temp2, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsSimMB.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	if (candidateDataProcessorEmbedHT1) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandSimHT1", analysisSettingsNew.title + ", in jet cand fraction HT1 MC"
		, candidateDataProcessorEmbedHT1->getInJetPointYes(), candidateDataProcessorEmbedHT1->getPt(), 0, 100
		, analysisResults.resultsSimHT1.invlist, temp1, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvSimHT1", analysisSettingsNew.title + ", jet ev fraction HT1 MC"
		, candidateDataProcessorEmbedHT1->getPt(), candidateDataProcessorEmbedHT1->getJetYes(), 0, 100
		, analysisResults.resultsSimHT1.invlist, temp2, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsSimHT1.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	if (candidateDataProcessorEmbedHT2) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandSimHT2", analysisSettingsNew.title + ", in jet cand fraction HT2 MC"
		, candidateDataProcessorEmbedHT2->getInJetPointYes(), candidateDataProcessorEmbedHT2->getPt(), 0, 100
		, analysisResults.resultsSimHT2.invlist, temp1, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvSimHT2", analysisSettingsNew.title + ", jet ev fraction HT2 MC"
		, candidateDataProcessorEmbedHT2->getPt(), candidateDataProcessorEmbedHT2->getJetYes(), 0, 100
		, analysisResults.resultsSimHT2.invlist, temp2, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsSimHT2.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
        }
	list<bin_stat_list_type> bgJetToTotalFromCandidatesCorrListSim;
	bgJetToTotalFromCandidatesCorrListSim.push_back(analysisResults.resultsSimMB.bgJetToTotalFromCandidatesCorr);
	bgJetToTotalFromCandidatesCorrListSim.push_back(analysisResults.resultsSimHT1.bgJetToTotalFromCandidatesCorr);
	bgJetToTotalFromCandidatesCorrListSim.push_back(analysisResults.resultsSimHT2.bgJetToTotalFromCandidatesCorr);
	bin_stat_list_type bgJetToTotalFromCandidatesCorrSim, tempSim;
	mergePoints(bgJetToTotalFromCandidatesCorrListSim, bgJetToTotalFromCandidatesCorrSim, temp);
	if ((analysisSettingsNew.settingsSimMB.fitJetBgFractionFromCandidates && histJetFracSimMB) || analysisSettingsNew.settingsSimMB.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsSimMB.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsSimMB.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracSimMB, analysisSettingsNew.settingsSimMB.jetBgFraction, analysisResults.resultsSimMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimMB");
	    } else if (analysisSettingsNew.settingsSimMB.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorrSim, 0, analysisSettingsNew.settingsSimMB.jetBgFraction, analysisResults.resultsSimMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimMB");
	    }
	    if (analysisSettingsNew.settingsSimMB.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsSimMB.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsSimMB.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsSimMB.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsSimMB.jetBgFraction.wasFit = false;
	}
	if ((analysisSettingsNew.settingsSimHT1.fitJetBgFractionFromCandidates && histJetFracSimHT1) || analysisSettingsNew.settingsSimHT1.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsSimHT1.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsSimHT1.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracSimHT1, analysisSettingsNew.settingsSimHT1.jetBgFraction, analysisResults.resultsSimHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimHT1");
	    } else if (analysisSettingsNew.settingsSimHT1.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorrSim, 0, analysisSettingsNew.settingsSimHT1.jetBgFraction, analysisResults.resultsSimHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimHT1");
	    }
	    if (analysisSettingsNew.settingsSimHT1.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsSimHT1.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsSimHT1.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsSimHT1.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsSimHT1.jetBgFraction.wasFit = false;
	}
	if ((analysisSettingsNew.settingsSimHT2.fitJetBgFractionFromCandidates && histJetFracSimHT2) || analysisSettingsNew.settingsSimHT2.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsSimHT2.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsSimHT2.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracSimHT2, analysisSettingsNew.settingsSimHT2.jetBgFraction, analysisResults.resultsSimHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimHT2");
	    } else if (analysisSettingsNew.settingsSimHT2.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorrSim, 0, analysisSettingsNew.settingsSimHT2.jetBgFraction, analysisResults.resultsSimHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimHT2");
	    }
	    if (analysisSettingsNew.settingsSimHT2.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsSimHT2.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsSimHT2.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsSimHT2.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsSimHT2.jetBgFraction.wasFit = false;
	}

	TH1F *histJetFracSimEtaMB = 0;
	TH1F *histJetFracSimEtaHT1 = 0;
	TH1F *histJetFracSimEtaHT2 = 0;
	{
		const TH1F *histPtSimEtaMB = candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getPt() : 0;
		const TH1F *histJetYesSimEtaMB = candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getJetYes() : 0;
		const TH1F *histInJetYesSimEtaMB = candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getInJetPointYes() : 0;
		if (histInJetYesSimEtaMB) histJetFracSimEtaMB = new TH1F(*histInJetYesSimEtaMB);
		if (histJetFracSimEtaMB && histPtSimEtaMB) histJetFracSimEtaMB->Divide(histPtSimEtaMB);
		if (histJetFracSimEtaMB && histPtSimEtaMB && histJetYesSimEtaMB) {
		    histJetFracSimEtaMB->Multiply(histJetYesSimEtaMB);
		    histJetFracSimEtaMB->Divide(histPtSimEtaMB);
		}
	}
	{
		const TH1F *histPtSimEtaHT1 = candidateDataProcessorEmbedEtaHT1 ? candidateDataProcessorEmbedEtaHT1->getPt() : 0;
		const TH1F *histJetYesSimEtaHT1 = candidateDataProcessorEmbedEtaHT1 ? candidateDataProcessorEmbedEtaHT1->getJetYes() : 0;
		const TH1F *histInJetYesSimEtaHT1 = candidateDataProcessorEmbedEtaHT1 ? candidateDataProcessorEmbedEtaHT1->getInJetPointYes() : 0;
		if (histInJetYesSimEtaHT1) histJetFracSimEtaHT1 = new TH1F(*histInJetYesSimEtaHT1);
		if (histJetFracSimEtaHT1 && histPtSimEtaHT1) histJetFracSimEtaHT1->Divide(histPtSimEtaHT1);
		if (histJetFracSimEtaHT1 && histPtSimEtaHT1 && histJetYesSimEtaHT1) {
		    histJetFracSimEtaHT1->Multiply(histJetYesSimEtaHT1);
		    histJetFracSimEtaHT1->Divide(histPtSimEtaHT1);
		}
	}
	{
		const TH1F *histPtSimEtaHT2 = candidateDataProcessorEmbedEtaHT2 ? candidateDataProcessorEmbedEtaHT2->getPt() : 0;
		const TH1F *histJetYesSimEtaHT2 = candidateDataProcessorEmbedEtaHT2 ? candidateDataProcessorEmbedEtaHT2->getJetYes() : 0;
		const TH1F *histInJetYesSimEtaHT2 = candidateDataProcessorEmbedEtaHT2 ? candidateDataProcessorEmbedEtaHT2->getInJetPointYes() : 0;
		if (histInJetYesSimEtaHT2) histJetFracSimEtaHT2 = new TH1F(*histInJetYesSimEtaHT2);
		if (histJetFracSimEtaHT2 && histPtSimEtaHT2) histJetFracSimEtaHT2->Divide(histPtSimEtaHT2);
		if (histJetFracSimEtaHT2 && histPtSimEtaHT2 && histJetYesSimEtaHT2) {
		    histJetFracSimEtaHT2->Multiply(histJetYesSimEtaHT2);
		    histJetFracSimEtaHT2->Divide(histPtSimEtaHT2);
		}
	}
	if (candidateDataProcessorEmbedEtaMB) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandSimEtaMB", analysisSettingsNew.title + ", in jet cand fraction MB eta MC"
		, candidateDataProcessorEmbedEtaMB->getInJetPointYes(), candidateDataProcessorEmbedEtaMB->getPt(), 0, 100
		, analysisResults.resultsSimEtaMB.invlist, temp1, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvSimEtaMB", analysisSettingsNew.title + ", jet ev fraction MB eta MC"
		, candidateDataProcessorEmbedEtaMB->getPt(), candidateDataProcessorEmbedEtaMB->getJetYes(), 0, 100
		, analysisResults.resultsSimEtaMB.invlist, temp2, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsSimEtaMB.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	if (candidateDataProcessorEmbedEtaHT1) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandSimEtaHT1", analysisSettingsNew.title + ", in jet cand fraction HT1 eta MC"
		, candidateDataProcessorEmbedEtaHT1->getInJetPointYes(), candidateDataProcessorEmbedEtaHT1->getPt(), 0, 100
		, analysisResults.resultsSimEtaHT1.invlist, temp1, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvSimEtaHT1", analysisSettingsNew.title + ", jet ev fraction HT1 eta MC"
		, candidateDataProcessorEmbedEtaHT1->getPt(), candidateDataProcessorEmbedEtaHT1->getJetYes(), 0, 100
		, analysisResults.resultsSimEtaHT1.invlist, temp2, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsSimEtaHT1.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	if (candidateDataProcessorEmbedEtaHT2) {
	    bin_stat_list_type temp1, temp2;
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgInJetCandSimEtaHT2", analysisSettingsNew.title + ", in jet cand fraction HT2 eta MC"
		, candidateDataProcessorEmbedEtaHT2->getInJetPointYes(), candidateDataProcessorEmbedEtaHT2->getPt(), 0, 100
		, analysisResults.resultsSimEtaHT2.invlist, temp1, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    getScaleFactorFromPt(analysisSettingsNew.name + "_bgJetEvSimEtaHT2", analysisSettingsNew.title + ", jet ev fraction HT2 eta MC"
		, candidateDataProcessorEmbedEtaHT2->getPt(), candidateDataProcessorEmbedEtaHT2->getJetYes(), 0, 100
		, analysisResults.resultsSimEtaHT2.invlist, temp2, analysisSettingsNew.showSim && false, analysisSettingsNew.print && false);
	    divideBins(temp1, temp2, analysisResults.resultsSimEtaHT2.bgJetToTotalFromCandidatesCorr, 1, 0, false, false);
	}
	list<bin_stat_list_type> bgJetToTotalFromCandidatesCorrListSimEta;
	bgJetToTotalFromCandidatesCorrListSimEta.push_back(analysisResults.resultsSimEtaMB.bgJetToTotalFromCandidatesCorr);
	bgJetToTotalFromCandidatesCorrListSimEta.push_back(analysisResults.resultsSimEtaHT1.bgJetToTotalFromCandidatesCorr);
	bgJetToTotalFromCandidatesCorrListSimEta.push_back(analysisResults.resultsSimEtaHT2.bgJetToTotalFromCandidatesCorr);
	bin_stat_list_type bgJetToTotalFromCandidatesCorrSimEta, tempSimEta;
	mergePoints(bgJetToTotalFromCandidatesCorrListSimEta, bgJetToTotalFromCandidatesCorrSimEta, tempSimEta);
	if ((analysisSettingsNew.settingsSimEtaMB.fitJetBgFractionFromCandidates && histJetFracSimEtaMB) || analysisSettingsNew.settingsSimEtaMB.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsSimEtaMB.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsSimEtaMB.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracSimEtaMB, analysisSettingsNew.settingsSimEtaMB.jetBgFraction, analysisResults.resultsSimEtaMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimEtaMB");
	    } else if (analysisSettingsNew.settingsSimEtaMB.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorrSimEta, 0, analysisSettingsNew.settingsSimEtaMB.jetBgFraction, analysisResults.resultsSimEtaMB.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimEtaMB");
	    }
	    if (analysisSettingsNew.settingsSimEtaMB.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsSimEtaMB.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsSimEtaMB.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsSimEtaMB.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsSimEtaMB.jetBgFraction.wasFit = false;
	}
	if ((analysisSettingsNew.settingsSimEtaHT1.fitJetBgFractionFromCandidates && histJetFracSimEtaHT1) || analysisSettingsNew.settingsSimEtaHT1.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsSimEtaHT1.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsSimEtaHT1.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracSimEtaHT1, analysisSettingsNew.settingsSimEtaHT1.jetBgFraction, analysisResults.resultsSimEtaHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimEtaHT1");
	    } else if (analysisSettingsNew.settingsSimEtaHT1.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorrSimEta, 0, analysisSettingsNew.settingsSimEtaHT1.jetBgFraction, analysisResults.resultsSimEtaHT1.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimEtaHT1");
	    }
	    if (analysisSettingsNew.settingsSimEtaHT1.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsSimEtaHT1.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsSimEtaHT1.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsSimEtaHT1.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsSimEtaHT1.jetBgFraction.wasFit = false;
	}
	if ((analysisSettingsNew.settingsSimEtaHT2.fitJetBgFractionFromCandidates && histJetFracSimEtaHT2) || analysisSettingsNew.settingsSimEtaHT2.getJetBgFractionFromCandidatesPtFit) {
	    analysisSettingsNew.settingsSimEtaHT2.jetBgFraction.useFit = true;
	    if (analysisSettingsNew.settingsSimEtaHT2.fitJetBgFractionFromCandidates) {
		fitParameter(0, histJetFracSimEtaHT2, analysisSettingsNew.settingsSimEtaHT2.jetBgFraction, analysisResults.resultsSimEtaHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimEtaHT2");
	    } else if (analysisSettingsNew.settingsSimEtaHT2.getJetBgFractionFromCandidatesPtFit) {
		fitParameter(&bgJetToTotalFromCandidatesCorrSimEta, 0, analysisSettingsNew.settingsSimEtaHT2.jetBgFraction, analysisResults.resultsSimEtaHT2.jetBgFractionCandidatesFit, analysisSettingsNew.name + "_jetbgfitCandSimEtaHT2");
	    }
	    if (analysisSettingsNew.settingsSimEtaHT2.jetBgFraction.func) for (Int_t ipar = 0;ipar < analysisSettingsNew.settingsSimEtaHT2.jetBgFraction.func->GetNpar();ipar++) analysisSettingsNew.settingsSimEtaHT2.jetBgFraction.func->FixParameter(ipar, analysisSettingsNew.settingsSimEtaHT2.jetBgFraction.func->GetParameter(ipar));
	    analysisSettingsNew.settingsSimEtaHT2.jetBgFraction.wasFit = false;
	}

	TWeightCalculator wgammas("weight_1gamma_pi0", "Single photon weight");
	Bool_t wgammasGood = false;
	{
	    const TH1F *histGammasPt = simuDataProcessorMB ? simuDataProcessorMB->gammas.getPt() : 0;
	    if (histGammasPt) {
		TH1F gammas(*histGammasPt);
    		gammas.GetXaxis()->SetRangeUser(analysisSettingsNew.settingsSimMB.lowmassbgPointsLowPt, analysisSettingsNew.settingsSimHT2.lowmassbgPointsHighPt);
		wgammas.Fit(&gammas);
		wgammas.mult = wgammas.multDrift;
		wgammas.multDrift = 0.0;
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
		wgammasEta.mult = wgammasEta.multDrift;
		wgammasEta.multDrift = 0.0;
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
		wpointsMB.mult = wpointsMB.multDrift;
		wpointsMB.multDrift = 0.0;
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
		wpointsHT1.mult = wpointsHT1.multDrift;
		wpointsHT1.multDrift = 0.0;
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
		wpointsHT2.mult = wpointsHT2.multDrift;
		wpointsHT2.multDrift = 0.0;
		wpointsHT2.rangeLow = 0;
		wpointsHT2.rangeHigh = 100;
		wpointsHT2Good = true;
	    }
	}

	TWeightCalculator weightUsed("weightUsed", "Used weight");
	if (candidateDataProcessorEmbedMB) weightUsed = candidateDataProcessorEmbedMB->getWeightCalculator();
	TWeightCalculator weightUsedEta("weightUsedEta", "Used weight #eta");
	if (candidateDataProcessorEmbedEtaMB) weightUsedEta = candidateDataProcessorEmbedEtaMB->getWeightCalculator();

	analysisResults.totSimuYield = simuDataProcessorMB ? simuDataProcessorMB->binStatistics.back().getValue() : 0;
	analysisResults.totSimuYieldEta = simuDataProcessorEtaMB ? simuDataProcessorEtaMB->binStatistics.back().getValue() : 0;
	analysisResults.totSimuYieldEtabg = simuDataProcessorEtabgMB ? simuDataProcessorEtabgMB->binStatistics.back().getValue() : 0;

	analysisResults.resultsMB.scalePeakShapeEtabgRel = 
	    ((analysisResults.resultsSimEtabgMB.evNumPassedAllCuts > 0) && (analysisResults.resultsSimEtaMB.evNumPassedAllCuts > 0)) ? 
		((Float_t)analysisResults.resultsSimEtaMB.evNumPassedAllCuts / (Float_t)analysisResults.resultsSimEtabgMB.evNumPassedAllCuts) : 
		1.0;
	analysisResults.resultsHT1.scalePeakShapeEtabgRel = analysisResults.resultsMB.scalePeakShapeEtabgRel;
	analysisResults.resultsHT2.scalePeakShapeEtabgRel = analysisResults.resultsMB.scalePeakShapeEtabgRel;

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
	analysisResults.resultsMB.scalePeakShapeLowMass2 = analysisResults.resultsSimMB.scaleBgLowMass2;
	analysisResults.resultsMB.scalePeakShapeEtaLowMass2 = analysisResults.resultsSimEtaMB.scaleBgLowMass2;
	analysisResults.resultsMB.scalePeakShapeEtabgLowMass2 = analysisResults.resultsSimEtabgMB.scaleBgLowMass2;
	analysisResults.resultsHT1.scalePeakShapeLowMass2 = analysisResults.resultsSimHT1.scaleBgLowMass2;
	analysisResults.resultsHT1.scalePeakShapeEtaLowMass2 = analysisResults.resultsSimEtaHT1.scaleBgLowMass2;
	analysisResults.resultsHT1.scalePeakShapeEtabgLowMass2 = analysisResults.resultsSimEtabgHT1.scaleBgLowMass2;
	analysisResults.resultsHT2.scalePeakShapeLowMass2 = analysisResults.resultsSimHT2.scaleBgLowMass2;
	analysisResults.resultsHT2.scalePeakShapeEtaLowMass2 = analysisResults.resultsSimEtaHT2.scaleBgLowMass2;
	analysisResults.resultsHT2.scalePeakShapeEtabgLowMass2 = analysisResults.resultsSimEtabgHT2.scaleBgLowMass2;

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
	fitParameter(&analysisResults.effMBSmooth, 0, analysisSettingsNew.effFitMB, analysisResults.effFitMB, analysisSettingsNew.name + "_effFitMB");
	if (analysisSettingsNew.effFitMB.useFit && analysisSettingsNew.effFitMB.wasFit) setBinsToFit(analysisResults.effMBSmooth, analysisSettingsNew.effFitMB);
	if (analysisSettingsNew.smoothEffHT1) {
		smoothPoints(analysisResults.effHT1, analysisResults.effHT1Smooth, analysisResults.effHT1SmoothPoint);
	} else {
		analysisResults.effHT1Smooth = analysisResults.effHT1;
		analysisResults.effHT1SmoothPoint = analysisResults.effHT1;
	}
	fitParameter(&analysisResults.effHT1Smooth, 0, analysisSettingsNew.effFitHT1, analysisResults.effFitHT1, analysisSettingsNew.name + "_effFitHT1");
	if (analysisSettingsNew.effFitHT1.useFit && analysisSettingsNew.effFitHT1.wasFit) setBinsToFit(analysisResults.effHT1Smooth, analysisSettingsNew.effFitHT1);
	if (analysisSettingsNew.smoothEffHT2) {
		smoothPoints(analysisResults.effHT2, analysisResults.effHT2Smooth, analysisResults.effHT2SmoothPoint);
	} else {
		analysisResults.effHT2Smooth = analysisResults.effHT2;
		analysisResults.effHT2SmoothPoint = analysisResults.effHT2;
	}
	fitParameter(&analysisResults.effHT2Smooth, 0, analysisSettingsNew.effFitHT2, analysisResults.effFitHT2, analysisSettingsNew.name + "_effFitHT2");
	if (analysisSettingsNew.effFitHT2.useFit && analysisSettingsNew.effFitHT2.wasFit) setBinsToFit(analysisResults.effHT2Smooth, analysisSettingsNew.effFitHT2);
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
	fitParameter(&analysisResults.effMBSmoothEta, 0, analysisSettingsNew.effFitMBEta, analysisResults.effFitMBEta, analysisSettingsNew.name + "_effFitMBEta");
	if (analysisSettingsNew.effFitMBEta.useFit && analysisSettingsNew.effFitMBEta.wasFit) setBinsToFit(analysisResults.effMBSmoothEta, analysisSettingsNew.effFitMBEta);
	if (analysisSettingsNew.smoothEffHT1Eta) {
		smoothPoints(analysisResults.effHT1Eta, analysisResults.effHT1SmoothEta, analysisResults.effHT1SmoothPointEta);
	} else {
		analysisResults.effHT1SmoothEta = analysisResults.effHT1Eta;
		analysisResults.effHT1SmoothPointEta = analysisResults.effHT1Eta;
	}
	fitParameter(&analysisResults.effHT1SmoothEta, 0, analysisSettingsNew.effFitHT1Eta, analysisResults.effFitHT1Eta, analysisSettingsNew.name + "_effFitHT1Eta");
	if (analysisSettingsNew.effFitHT1Eta.useFit && analysisSettingsNew.effFitHT1Eta.wasFit) setBinsToFit(analysisResults.effHT1SmoothEta, analysisSettingsNew.effFitHT1Eta);
	if (analysisSettingsNew.smoothEffHT2Eta) {
		smoothPoints(analysisResults.effHT2Eta, analysisResults.effHT2SmoothEta, analysisResults.effHT2SmoothPointEta);
	} else {
		analysisResults.effHT2SmoothEta = analysisResults.effHT2Eta;
		analysisResults.effHT2SmoothPointEta = analysisResults.effHT2Eta;
	}
	fitParameter(&analysisResults.effHT2SmoothEta, 0, analysisSettingsNew.effFitHT2Eta, analysisResults.effFitHT2Eta, analysisSettingsNew.name + "_effFitHT2Eta");
	if (analysisSettingsNew.effFitHT2Eta.useFit && analysisSettingsNew.effFitHT2Eta.wasFit) setBinsToFit(analysisResults.effHT2SmoothEta, analysisSettingsNew.effFitHT2Eta);
	if (analysisSettingsNew.dontCorrectEffEta) {
	    bin_stat_list_type tempMB, tempHT1, tempHT2;
	    resetBinsError(analysisResults.effMBSmoothEta, tempMB);
	    analysisResults.effMBSmoothEta = tempMB;
	    resetBinsError(analysisResults.effHT1SmoothEta, tempHT1);
	    analysisResults.effHT1SmoothEta = tempHT1;
	    resetBinsError(analysisResults.effHT2SmoothEta, tempHT2);
	    analysisResults.effHT2SmoothEta = tempHT2;
	}

	divideBins(analysisResults.resultsSimMB.inputYieldGamma, analysisResults.areaSimuMB, analysisResults.gammaToPi0DecayMB, 1.0/truePionBranchingRatio, 0.0, false, false);
	divideBins(analysisResults.resultsSimHT1.inputYieldGamma, analysisResults.areaSimuHT1, analysisResults.gammaToPi0DecayHT1, 1.0/truePionBranchingRatio, 0.0, false, false);
	divideBins(analysisResults.resultsSimHT2.inputYieldGamma, analysisResults.areaSimuHT2, analysisResults.gammaToPi0DecayHT2, 1.0/truePionBranchingRatio, 0.0, false, false);

	divideBins(analysisResults.resultsSimEtaMB.inputYieldGamma, analysisResults.areaSimuEtaMB, analysisResults.gammaToEtaDecayMB, 1.0/trueEtaBranchingRatio, 0.0, false, false);
	divideBins(analysisResults.resultsSimEtaHT1.inputYieldGamma, analysisResults.areaSimuEtaHT1, analysisResults.gammaToEtaDecayHT1, 1.0/trueEtaBranchingRatio, 0.0, false, false);
	divideBins(analysisResults.resultsSimEtaHT2.inputYieldGamma, analysisResults.areaSimuEtaHT2, analysisResults.gammaToEtaDecayHT2, 1.0/trueEtaBranchingRatio, 0.0, false, false);

	divideBins(analysisResults.resultsSimMB.rawYieldGamma, analysisResults.resultsSimMB.inputYieldGamma, analysisResults.effGammaPi0MB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsSimHT1.rawYieldGamma, analysisResults.resultsSimHT1.inputYieldGamma, analysisResults.effGammaPi0HT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsSimHT2.rawYieldGamma, analysisResults.resultsSimHT2.inputYieldGamma, analysisResults.effGammaPi0HT2, 1.0, 0.0, false, false);

	divideBins(analysisResults.resultsSimEtaMB.rawYieldGamma, analysisResults.resultsSimEtaMB.inputYieldGamma, analysisResults.effGammaEtaMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsSimEtaHT1.rawYieldGamma, analysisResults.resultsSimEtaHT1.inputYieldGamma, analysisResults.effGammaEtaHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsSimEtaHT2.rawYieldGamma, analysisResults.resultsSimEtaHT2.inputYieldGamma, analysisResults.effGammaEtaHT2, 1.0, 0.0, false, false);

	divideBins(analysisResults.rawYieldGammaMB, analysisResults.inputYieldGammaMB, analysisResults.effGammaMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.rawYieldGammaHT1, analysisResults.inputYieldGammaHT1, analysisResults.effGammaHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.rawYieldGammaHT2, analysisResults.inputYieldGammaHT2, analysisResults.effGammaHT2, 1.0, 0.0, false, false);

	divideBins(analysisResults.rawYieldNbarMB, analysisResults.inputYieldNbarMB, analysisResults.effNbarMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.rawYieldNbarHT1, analysisResults.inputYieldNbarHT1, analysisResults.effNbarHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.rawYieldNbarHT2, analysisResults.inputYieldNbarHT2, analysisResults.effNbarHT2, 1.0, 0.0, false, false);

	TLegend *legendRefmults = 0;
	TString nameRefmults = analysisSettingsNew.name + "_refmults";
	TString titleRefmults = analysisSettingsNew.title + ", RefMult per MinBias event vs. day;Day;RefMult per MinBias event";
	TH1F *histRefmults = 0;
	{
	    TDrawOptions temp = analysisSettingsNew.drawOptions;
	    temp.legendTitle = "East FTPC RefMult per MinBias event";
	    temp.setColor(kBlue);
	    histRefmults = showList(&analysisResults.refMultEast, nameRefmults.Data(), titleRefmults.Data()
	        , histRefmults, &legendRefmults
	        , temp
	    );
	    temp.legendTitle = "West FTPC RefMult per MinBias event";
	    temp.setColor(kRed);
	    histRefmults = showList(&analysisResults.refMultWest, nameRefmults.Data(), titleRefmults.Data()
	        , histRefmults, &legendRefmults
	        , temp
	    );
	}
	if (histRefmults) {
		histRefmults->SetMinimum(0);
		histRefmults->SetMaximum(20);
	}

	TLegend *legendEff = 0;
	TString effName = analysisSettingsNew.name + "_eff";
	TString effTitle = analysisSettingsNew.title + ", Efficiency and acceptance correction;p_{T}, GeV/c;Correction factor";
	TH1F *histEff = 0;
	histEff = showLists(&analysisResults.effMB, &analysisResults.effHT1, &analysisResults.effHT2
	    , effName.Data(), effTitle.Data(), histEff, &legendEff
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	histEff = showLists(&analysisResults.effMBEta, &analysisResults.effHT1Eta, &analysisResults.effHT2Eta
	    , effName.Data(), effTitle.Data(), histEff, &legendEff
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	    );
	if (histEff) {
		histEff->SetMinimum(0);
		histEff->SetMaximum(0.4);
		histEff->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
		if (analysisSettingsNew.smoothEffMB) {
		    TDrawOptions temp = analysisSettingsNew.settingsSimMB.drawOptions;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effMBSmoothPoint
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT1) {
		    TDrawOptions temp = analysisSettingsNew.settingsSimHT1.drawOptions;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT1SmoothPoint
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT2) {
		    TDrawOptions temp = analysisSettingsNew.settingsSimHT2.drawOptions;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT2SmoothPoint
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.effFitMB.func) {
		    analysisSettingsNew.effFitMB.func->SetLineColor(analysisSettingsNew.settingsMB.color);
		    analysisSettingsNew.effFitMB.func->Draw("SAME");
		}
		if (analysisSettingsNew.effFitHT1.func) {
		    analysisSettingsNew.effFitHT1.func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    analysisSettingsNew.effFitHT1.func->Draw("SAME");
		}
		if (analysisSettingsNew.effFitHT2.func) {
		    analysisSettingsNew.effFitHT2.func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    analysisSettingsNew.effFitHT2.func->Draw("SAME");
		}

		if (analysisSettingsNew.smoothEffMBEta) {
		    TDrawOptions temp = analysisSettingsNew.settingsSimEtaMB.drawOptionsEta;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effMBSmoothPointEta
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT1Eta) {
		    TDrawOptions temp = analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT1SmoothPointEta
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.smoothEffHT2Eta) {
		    TDrawOptions temp = analysisSettingsNew.settingsSimEtaHT2.drawOptions;
		    temp.drawOption = "C";
		    temp.legendOption = "C";
		    temp.legendTitle += " smoothed";
    		    histEff = showList(&analysisResults.effHT2SmoothPointEta
			, effName.Data(), effTitle.Data(), histEff, &legendEff
			, temp
		    );
		}
		if (analysisSettingsNew.effFitMBEta.func) {
		    analysisSettingsNew.effFitMBEta.func->SetLineColor(analysisSettingsNew.settingsSimEtaMB.color);
		    analysisSettingsNew.effFitMBEta.func->Draw("SAME");
		}
		if (analysisSettingsNew.effFitHT1Eta.func) {
		    analysisSettingsNew.effFitHT1Eta.func->SetLineColor(analysisSettingsNew.settingsSimEtaHT1.color);
		    analysisSettingsNew.effFitHT1Eta.func->Draw("SAME");
		}
		if (analysisSettingsNew.effFitHT2Eta.func) {
		    analysisSettingsNew.effFitHT2Eta.func->SetLineColor(analysisSettingsNew.settingsSimEtaHT2.color);
		    analysisSettingsNew.effFitHT2Eta.func->Draw("SAME");
		}
	}

	TLegend *legendEffGamma = 0;
	TString effGammaName = analysisSettingsNew.name + "_effGamma";
	TString effGammaTitle = analysisSettingsNew.title + ", Efficiency and acceptance correction for #gamma;p_{T}, GeV/c;Correction factor";
	TH1F *histEffGamma = 0;
	histEffGamma = showLists(&analysisResults.effGammaPi0MB, &analysisResults.effGammaPi0HT1, &analysisResults.effGammaPi0HT2
	    , effGammaName.Data(), effGammaTitle.Data(), histEffGamma, &legendEffGamma
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	histEffGamma = showLists(&analysisResults.effGammaEtaMB, &analysisResults.effGammaEtaHT1, &analysisResults.effGammaEtaHT2
	    , effGammaName.Data(), effGammaTitle.Data(), histEffGamma, &legendEffGamma
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	    );
	histEffGamma = showLists(&analysisResults.effGammaMB, &analysisResults.effGammaHT1, &analysisResults.effGammaHT2
	    , effGammaName.Data(), effGammaTitle.Data(), histEffGamma, &legendEffGamma
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histEffGamma = showLists(&analysisResults.effNbarMB, &analysisResults.effNbarHT1, &analysisResults.effNbarHT2
	    , effGammaName.Data(), effGammaTitle.Data(), histEffGamma, &legendEffGamma
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histEffGamma) {
		histEffGamma->SetMinimum(0);
		histEffGamma->SetMaximum(1.5);
		histEffGamma->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
	}
/*
	TLegend *legendEffGammaIncl = 0;
	TString effGammaInclName = analysisSettingsNew.name + "_effGammaIncl";
	TString effGammaInclTitle = analysisSettingsNew.title + ", Efficiency and acceptance correction for inclusive #gamma;p_{T}, GeV/c;Correction factor";
	TH1F *histEffGammaIncl = 0;
	histEffGammaIncl = showLists(&analysisResults.effGammaMBIncl, &analysisResults.effGammaHT1Incl, &analysisResults.effGammaHT2Incl
	    , effGammaInclName.Data(), effGammaInclTitle.Data(), histEffGammaIncl, &legendEffGammaIncl
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histEffGammaIncl) {
		histEffGammaIncl->SetMinimum(0);
		histEffGammaIncl->SetMaximum(1.5);
		histEffGammaIncl->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
	}
*/
/*
	TLegend *legendFractionGamma = 0;
	TString fractionGammaName = analysisSettingsNew.name + "_fractionGamma";
	TString fractionGammaTitle = analysisSettingsNew.title + ", Fraction of decay #gamma;p_{T}, GeV/c;Fraction";
	TH1F *histFractionGamma = 0;
	histFractionGamma = showLists(&analysisResults.resultsSimMB.rawYieldGammaFraction, &analysisResults.resultsSimHT1.rawYieldGammaFraction, &analysisResults.resultsSimHT2.rawYieldGammaFraction
	    , fractionGammaName.Data(), fractionGammaTitle.Data(), histFractionGamma, &legendFractionGamma
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	histFractionGamma = showLists(&analysisResults.resultsSimEtaMB.rawYieldGammaFraction, &analysisResults.resultsSimEtaHT1.rawYieldGammaFraction, &analysisResults.resultsSimEtaHT2.rawYieldGammaFraction
	    , fractionGammaName.Data(), fractionGammaTitle.Data(), histFractionGamma, &legendFractionGamma
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	    );
	if (histFractionGamma) {
		histFractionGamma->SetMinimum(0);
		histFractionGamma->SetMaximum(1.0);
		histFractionGamma->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
	}
*/
	TLegend *legendBgToSigbg = 0;
	TString BgToSigbgName = analysisSettingsNew.name + "_BgToSigbg";
	TString BgToSigbgTitle = analysisSettingsNew.title + ", Bg / Sig+Bg;p_{T}, GeV/c;Bg / Sig+Bg";
	TH1F *histBgToSigbg = 0;
	histBgToSigbg = showLists(&analysisResults.resultsMB.bgToSigbgFirst, &analysisResults.resultsHT1.bgToSigbgFirst, &analysisResults.resultsHT2.bgToSigbgFirst
	    , BgToSigbgName.Data(), BgToSigbgTitle.Data(), histBgToSigbg, &legendBgToSigbg
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histBgToSigbg = showLists(&analysisResults.resultsMB.bgToSigbg, &analysisResults.resultsHT1.bgToSigbg, &analysisResults.resultsHT2.bgToSigbg
	    , BgToSigbgName.Data(), BgToSigbgTitle.Data(), histBgToSigbg, &legendBgToSigbg
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histBgToSigbg) {
		histBgToSigbg->SetMinimum(0);
		histBgToSigbg->SetMaximum(1.0);
		histBgToSigbg->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		if (analysisResults.resultsMB.bgToSigbgParamFit.func) {
		    analysisResults.resultsMB.bgToSigbgParamFit.func->SetLineColor(analysisSettingsNew.settingsMB.color);
		    analysisResults.resultsMB.bgToSigbgParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT1.bgToSigbgParamFit.func) {
		    analysisResults.resultsHT1.bgToSigbgParamFit.func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    analysisResults.resultsHT1.bgToSigbgParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT2.bgToSigbgParamFit.func) {
		    analysisResults.resultsHT2.bgToSigbgParamFit.func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    analysisResults.resultsHT2.bgToSigbgParamFit.func->Draw("SAME");
		}
	}
	TLegend *legendSigbg = 0;
	TString SigbgName = analysisSettingsNew.name + "_Sigbg";
	TString SigbgTitle = analysisSettingsNew.title + ", Sig / Bg;p_{T}, GeV/c;Sig / Bg";
	TH1F *histSigbg = 0;
	histSigbg = showLists(&analysisResults.resultsMB.peakYieldToBackground, &analysisResults.resultsHT1.peakYieldToBackground, &analysisResults.resultsHT2.peakYieldToBackground
	    , SigbgName.Data(), SigbgTitle.Data(), histSigbg, &legendSigbg
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histSigbg = showLists(&analysisResults.resultsMB.peakYieldToBackgroundEta, &analysisResults.resultsHT1.peakYieldToBackgroundEta, &analysisResults.resultsHT2.peakYieldToBackgroundEta
	    , SigbgName.Data(), SigbgTitle.Data(), histSigbg, &legendSigbg
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	    );
	if (histSigbg) {
		histSigbg->SetMinimum(0);
		histSigbg->SetMaximum(12);
		histSigbg->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		if (analysisResults.resultsMB.peakSigBackgroundParamFit.wasFit || analysisResults.resultsHT1.peakSigBackgroundParamFit.wasFit || analysisResults.resultsHT2.peakSigBackgroundParamFit.wasFit) {
		    histSigbg = showLists(&analysisResults.resultsMB.peakSigBackgroundParamFit.valuesForFit, &analysisResults.resultsHT1.peakSigBackgroundParamFit.valuesForFit, &analysisResults.resultsHT2.peakSigBackgroundParamFit.valuesForFit
			, SigbgName.Data(), SigbgTitle.Data(), histSigbg, &legendSigbg
			, analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    	    );
		}
		if (analysisResults.resultsMB.peakSigBackgroundEtaParamFit.wasFit || analysisResults.resultsHT1.peakSigBackgroundEtaParamFit.wasFit || analysisResults.resultsHT2.peakSigBackgroundEtaParamFit.wasFit) {
		    histSigbg = showLists(&analysisResults.resultsMB.peakSigBackgroundEtaParamFit.valuesForFit, &analysisResults.resultsHT1.peakSigBackgroundEtaParamFit.valuesForFit, &analysisResults.resultsHT2.peakSigBackgroundEtaParamFit.valuesForFit
			, SigbgName.Data(), SigbgTitle.Data(), histSigbg, &legendSigbg
			, analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	    	    );
		}
		if (analysisResults.resultsMB.peakSigBackgroundParamFit.func) {
		    analysisResults.resultsMB.peakSigBackgroundParamFit.func->SetLineColor(analysisSettingsNew.settingsMB.color);
		    analysisResults.resultsMB.peakSigBackgroundParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT1.peakSigBackgroundParamFit.func) {
		    analysisResults.resultsHT1.peakSigBackgroundParamFit.func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    analysisResults.resultsHT1.peakSigBackgroundParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT2.peakSigBackgroundParamFit.func) {
		    analysisResults.resultsHT2.peakSigBackgroundParamFit.func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    analysisResults.resultsHT2.peakSigBackgroundParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsMB.peakSigBackgroundEtaParamFit.func) {
		    analysisResults.resultsMB.peakSigBackgroundEtaParamFit.func->SetLineColor(analysisSettingsNew.settingsMB.color);
		    analysisResults.resultsMB.peakSigBackgroundEtaParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT1.peakSigBackgroundEtaParamFit.func) {
		    analysisResults.resultsHT1.peakSigBackgroundEtaParamFit.func->SetLineColor(analysisSettingsNew.settingsHT1.color);
		    analysisResults.resultsHT1.peakSigBackgroundEtaParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsHT2.peakSigBackgroundEtaParamFit.func) {
		    analysisResults.resultsHT2.peakSigBackgroundEtaParamFit.func->SetLineColor(analysisSettingsNew.settingsHT2.color);
		    analysisResults.resultsHT2.peakSigBackgroundEtaParamFit.func->Draw("SAME");
		}
	}
	TLegend *legendYieldToSigbgRange = 0;
	TString YieldToSigbgRangeName = analysisSettingsNew.name + "_YieldToSigbgRange";
	TString YieldToSigbgRangeTitle = analysisSettingsNew.title + ", Sig / Sig+Bg full;p_{T}, GeV/c;Sig / Sig+Bg full";
	TH1F *histYieldToSigbgRange = 0;
	histYieldToSigbgRange = showLists(&analysisResults.resultsMB.peakYieldToSigBackgroundRange, &analysisResults.resultsHT1.peakYieldToSigBackgroundRange, &analysisResults.resultsHT2.peakYieldToSigBackgroundRange
	    , YieldToSigbgRangeName.Data(), YieldToSigbgRangeTitle.Data(), histYieldToSigbgRange, &legendYieldToSigbgRange
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histYieldToSigbgRange) {
		histYieldToSigbgRange->SetMinimum(0);
		histYieldToSigbgRange->SetMaximum(1);
	}
	TLegend *legendBgJetToTotal = 0;
	TString BgJetToTotalName = analysisSettingsNew.name + "_BgJetToTotal";
	TString BgJetToTotalTitle = analysisSettingsNew.title + ", Jet-Mixed Background Fraction;p_{T}, GeV/c;Jet-Mixed Background Fraction";
	TH1F *histBgJetToTotal = 0;
	histBgJetToTotal = showLists(&analysisResults.resultsMB.bgJetToTotalFirst, &analysisResults.resultsHT1.bgJetToTotalFirst, &analysisResults.resultsHT2.bgJetToTotalFirst
	    , BgJetToTotalName.Data(), BgJetToTotalTitle.Data(), histBgJetToTotal, &legendBgJetToTotal
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histBgJetToTotal) {
		histBgJetToTotal->SetMinimum(0);
		histBgJetToTotal->SetMaximum(1.0);
		histBgJetToTotal->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
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

	TLegend *legendBgJetToTotalSim = 0;
	TString BgJetToTotalSimName = analysisSettingsNew.name + "_BgJetToTotalSim";
	TString BgJetToTotalSimTitle = analysisSettingsNew.title + ", Jet-Mixed Background Fraction, MC pi0;p_{T}, GeV/c;Jet-Mixed Background Fraction";
	TH1F *histBgJetToTotalSim = 0;
	histBgJetToTotalSim = showLists(&analysisResults.resultsSimMB.bgJetToTotalFirst, &analysisResults.resultsSimHT1.bgJetToTotalFirst, &analysisResults.resultsSimHT2.bgJetToTotalFirst
	    , BgJetToTotalSimName.Data(), BgJetToTotalSimTitle.Data(), histBgJetToTotalSim, &legendBgJetToTotalSim
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	if (histBgJetToTotalSim) {
		histBgJetToTotalSim->SetMinimum(0);
		histBgJetToTotalSim->SetMaximum(1.0);
		histBgJetToTotalSim->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
		if (histJetFracSimMB) {
		    histJetFracSimMB->SetLineColor(analysisSettingsNew.settingsSimMB.color);
		    histJetFracSimMB->Draw("HIST SAME");
		}
		if (histJetFracSimHT1) {
		    histJetFracSimHT1->SetLineColor(analysisSettingsNew.settingsSimHT1.color);
		    histJetFracSimHT1->Draw("HIST SAME");
		}
		if (histJetFracSimHT2) {
		    histJetFracSimHT2->SetLineColor(analysisSettingsNew.settingsSimHT2.color);
		    histJetFracSimHT2->Draw("HIST SAME");
		}
		if (analysisResults.resultsSimMB.jetBgFractionParamFit.func) {
		    analysisResults.resultsSimMB.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsSimMB.color);
		    analysisResults.resultsSimMB.jetBgFractionParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsSimHT1.jetBgFractionParamFit.func) {
		    analysisResults.resultsSimHT1.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsSimHT1.color);
		    analysisResults.resultsSimHT1.jetBgFractionParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsSimHT2.jetBgFractionParamFit.func) {
		    analysisResults.resultsSimHT2.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsSimHT2.color);
		    analysisResults.resultsSimHT2.jetBgFractionParamFit.func->Draw("SAME");
		}
	    histBgJetToTotalSim = showLists(&analysisResults.resultsMB.bgJetToTotalFirst, &analysisResults.resultsHT1.bgJetToTotalFirst, &analysisResults.resultsHT2.bgJetToTotalFirst
		, BgJetToTotalSimName.Data(), BgJetToTotalSimTitle.Data(), histBgJetToTotalSim, &legendBgJetToTotalSim
		, analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
		);
	}
	TLegend *legendBgJetToTotalSimEta = 0;
	TString BgJetToTotalSimEtaName = analysisSettingsNew.name + "_BgJetToTotalSimEta";
	TString BgJetToTotalSimEtaTitle = analysisSettingsNew.title + ", Jet-Mixed Background Fraction, MC eta;p_{T}, GeV/c;Jet-Mixed Background Fraction";
	TH1F *histBgJetToTotalSimEta = 0;
	histBgJetToTotalSimEta = showLists(&analysisResults.resultsSimEtaMB.bgJetToTotalFirst, &analysisResults.resultsSimEtaHT1.bgJetToTotalFirst, &analysisResults.resultsSimEtaHT2.bgJetToTotalFirst
	    , BgJetToTotalSimEtaName.Data(), BgJetToTotalSimEtaTitle.Data(), histBgJetToTotalSimEta, &legendBgJetToTotalSimEta
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	    );
	if (histBgJetToTotalSimEta) {
		histBgJetToTotalSimEta->SetMinimum(0);
		histBgJetToTotalSimEta->SetMaximum(1.0);
		histBgJetToTotalSimEta->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimEtaHT2.highPt);
		if (histJetFracSimEtaMB) {
		    histJetFracSimEtaMB->SetLineColor(analysisSettingsNew.settingsSimEtaMB.color);
		    histJetFracSimEtaMB->Draw("HIST SAME");
		}
		if (histJetFracSimEtaHT1) {
		    histJetFracSimEtaHT1->SetLineColor(analysisSettingsNew.settingsSimEtaHT1.color);
		    histJetFracSimEtaHT1->Draw("HIST SAME");
		}
		if (histJetFracSimEtaHT2) {
		    histJetFracSimEtaHT2->SetLineColor(analysisSettingsNew.settingsSimEtaHT2.color);
		    histJetFracSimEtaHT2->Draw("HIST SAME");
		}
		if (analysisResults.resultsSimEtaMB.jetBgFractionParamFit.func) {
		    analysisResults.resultsSimEtaMB.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsSimEtaMB.color);
		    analysisResults.resultsSimEtaMB.jetBgFractionParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsSimEtaHT1.jetBgFractionParamFit.func) {
		    analysisResults.resultsSimEtaHT1.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsSimEtaHT1.color);
		    analysisResults.resultsSimEtaHT1.jetBgFractionParamFit.func->Draw("SAME");
		}
		if (analysisResults.resultsSimEtaHT2.jetBgFractionParamFit.func) {
		    analysisResults.resultsSimEtaHT2.jetBgFractionParamFit.func->SetLineColor(analysisSettingsNew.settingsSimEtaHT2.color);
		    analysisResults.resultsSimEtaHT2.jetBgFractionParamFit.func->Draw("SAME");
		}
	    histBgJetToTotalSimEta = showLists(&analysisResults.resultsMB.bgJetToTotalFirst, &analysisResults.resultsHT1.bgJetToTotalFirst, &analysisResults.resultsHT2.bgJetToTotalFirst
		, BgJetToTotalSimEtaName.Data(), BgJetToTotalSimEtaTitle.Data(), histBgJetToTotalSimEta, &legendBgJetToTotalSimEta
		, analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	        );
	}

	bin_stat_list_type lowmassToYieldExpectMB, lowmassToYieldExpectHT1, lowmassToYieldExpectHT2;
	divideBins(analysisResults.resultsSimMB.lowmassToYield, analysisResults.resultsMB.peakYieldToSigBackgroundRange, lowmassToYieldExpectMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsSimHT1.lowmassToYield, analysisResults.resultsHT1.peakYieldToSigBackgroundRange, lowmassToYieldExpectHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsSimHT2.lowmassToYield, analysisResults.resultsHT2.peakYieldToSigBackgroundRange, lowmassToYieldExpectHT2, 1.0, 0.0, false, false);

	TLegend *legendLowmassToYield = 0;
	TString lowmassToYieldName = analysisSettingsNew.name + "_lowmassToYield";
	TString lowmassToYieldTitle = analysisSettingsNew.title + ", low mass bg. / yield;p_{T}, GeV/c;Low mass bg. / yield";
	TH1F *histLowmassToYield = 0;
	histLowmassToYield = showLists(&analysisResults.resultsMB.lowmassToYield, &analysisResults.resultsHT1.lowmassToYield, &analysisResults.resultsHT2.lowmassToYield
	    , lowmassToYieldName.Data(), lowmassToYieldTitle.Data(), histLowmassToYield, &legendLowmassToYield
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histLowmassToYield = showLists(&analysisResults.resultsSimMB.lowmassToYield, &analysisResults.resultsSimHT1.lowmassToYield, &analysisResults.resultsSimHT2.lowmassToYield
	    , lowmassToYieldName.Data(), lowmassToYieldTitle.Data(), histLowmassToYield, &legendLowmassToYield
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	{
	TDrawOptions optMB = analysisSettingsNew.settingsMB.drawOptions;
	TDrawOptions optHT1 = analysisSettingsNew.settingsHT1.drawOptions;
	TDrawOptions optHT2 = analysisSettingsNew.settingsHT2.drawOptions;
	optMB.legendTitle += " (expected)";
	optHT1.legendTitle += " (expected)";
	optHT2.legendTitle += " (expected)";
	optMB.drawOption = "lp";
	optHT1.drawOption = "lp";
	optHT2.drawOption = "lp";
	optMB.legendOption = "lp";
	optHT1.legendOption = "lp";
	optHT2.legendOption = "lp";
	histLowmassToYield = showLists(&lowmassToYieldExpectMB, &lowmassToYieldExpectHT1, &lowmassToYieldExpectHT2
	    , lowmassToYieldName.Data(), lowmassToYieldTitle.Data(), histLowmassToYield, &legendLowmassToYield
	    , optMB, optHT1, optHT2
	    );
	}
	if (histLowmassToYield) {
		histLowmassToYield->SetMinimum(0);
		histLowmassToYield->SetMaximum(0.5);
		histLowmassToYield->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	TLegend *legendEtabgToYield = 0;
	TString etabgToYieldName = analysisSettingsNew.name + "_etabgToYield";
	TString etabgToYieldTitle = analysisSettingsNew.title + ", etabg / yield;p_{T}, GeV/c;etabg / yield";
	TH1F *histEtabgToYield = 0;
	histEtabgToYield = showLists(&analysisResults.resultsMB.etabgToYield, &analysisResults.resultsHT1.etabgToYield, &analysisResults.resultsHT2.etabgToYield
	    , etabgToYieldName.Data(), etabgToYieldTitle.Data(), histEtabgToYield, &legendEtabgToYield
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histEtabgToYield) {
	    histEtabgToYield->SetMinimum(0.0);
	    histEtabgToYield->SetMaximum(0.1);
	    histEtabgToYield->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	TLegend *legendMixJetmixDiff = 0;
	TString mixJetmixDiffName = analysisSettingsNew.name + "_mixJetmixDiff";
	TString mixJetmixDiffTitle = analysisSettingsNew.title + ", 0.5 * (mix - jetmix bg) bg/ ;p_{T}, GeV/c;0.5 * (mix-jetmix bg) / bg";
	TH1F *histMixJetmixDiff = 0;
	histMixJetmixDiff = showLists(&analysisResults.resultsMB.mixJetmixDiff, &analysisResults.resultsHT1.mixJetmixDiff, &analysisResults.resultsHT2.mixJetmixDiff
	    , mixJetmixDiffName.Data(), mixJetmixDiffTitle.Data(), histMixJetmixDiff, &legendMixJetmixDiff
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histMixJetmixDiff = showLists(&analysisResults.resultsMB.mixJetmixDiffEta, &analysisResults.resultsHT1.mixJetmixDiffEta, &analysisResults.resultsHT2.mixJetmixDiffEta
	    , mixJetmixDiffName.Data(), mixJetmixDiffTitle.Data(), histMixJetmixDiff, &legendMixJetmixDiff
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	    );
	if (histMixJetmixDiff) {
	    histMixJetmixDiff->SetMinimum(0.0);
	    histMixJetmixDiff->SetMaximum(1.0);
	    histMixJetmixDiff->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	fitParameter(&analysisResults.resultsMB.cpvCorrection, 0, analysisSettingsNew.settingsMB.cpvCorrection, analysisResults.resultsMB.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitMB");
	fitParameter(&analysisResults.resultsHT1.cpvCorrection, 0, analysisSettingsNew.settingsHT1.cpvCorrection, analysisResults.resultsHT1.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT1");
	fitParameter(&analysisResults.resultsHT2.cpvCorrection, 0, analysisSettingsNew.settingsHT2.cpvCorrection, analysisResults.resultsHT2.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT2");
	analysisResults.cpvCorrection.clear();
	mergePoints(analysisResults.resultsMB.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrection);
	mergePoints(analysisResults.resultsHT1.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrection);
	mergePoints(analysisResults.resultsHT2.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrection);
	fitParameter(&analysisResults.cpvCorrection, 0, analysisSettingsNew.cpvCorrection, analysisResults.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfit");

	fitParameter(&analysisResults.resultsSimMB.cpvCorrection, 0, analysisSettingsNew.settingsSimMB.cpvCorrection, analysisResults.resultsSimMB.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitMB_sim");
	fitParameter(&analysisResults.resultsSimHT1.cpvCorrection, 0, analysisSettingsNew.settingsSimHT1.cpvCorrection, analysisResults.resultsSimHT1.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT1_sim");
	fitParameter(&analysisResults.resultsSimHT2.cpvCorrection, 0, analysisSettingsNew.settingsSimHT2.cpvCorrection, analysisResults.resultsSimHT2.cpvCorrectionFit, analysisSettingsNew.name + "_cpvfitHT2_sim");
	analysisResults.cpvCorrectionSim.clear();
	mergePoints(analysisResults.resultsSimMB.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSim);
	mergePoints(analysisResults.resultsSimHT1.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSim);
	mergePoints(analysisResults.resultsSimHT2.cpvCorrectionFit.valuesForFit, analysisResults.cpvCorrectionSim);
	fitParameter(&analysisResults.cpvCorrectionSim, 0, analysisSettingsNew.cpvCorrectionSim, analysisResults.cpvCorrectionFitSim, analysisSettingsNew.name + "_cpvfit_sim");

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
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histCpvCorrection = showLists(&analysisResults.resultsMB.cpvCorrectionMult, &analysisResults.resultsHT1.cpvCorrectionMult, &analysisResults.resultsHT2.cpvCorrectionMult
	    , CpvCorrectionName.Data(), CpvCorrectionTitle.Data(), histCpvCorrection, &legendCpvCorrection
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
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
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histCpvCorrectionParam7) {
		histCpvCorrectionParam7->SetMinimum(0.0);
		histCpvCorrectionParam7->SetMaximum(25.0);
		histCpvCorrectionParam7->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}
	TLegend *legendCpvCorrectionSim = 0;
	TString CpvCorrectionSimName = analysisSettingsNew.name + "_CpvCorrectionSim";
	TString CpvCorrectionSimTitle = analysisSettingsNew.title + ", CPV correction MC pi0;p_{T}, GeV/c;CPV correction";
	TH1F *histCpvCorrectionSim = 0;
	histCpvCorrectionSim = showLists(&analysisResults.resultsSimMB.cpvCorrection, &analysisResults.resultsSimHT1.cpvCorrection, &analysisResults.resultsSimHT2.cpvCorrection
	    , CpvCorrectionSimName.Data(), CpvCorrectionSimTitle.Data(), histCpvCorrectionSim, &legendCpvCorrectionSim
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	if (histCpvCorrectionSim) {
		histCpvCorrectionSim->SetMinimum(0.6);
		histCpvCorrectionSim->SetMaximum(1.0);
		histCpvCorrectionSim->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
		if (analysisResults.cpvCorrectionFitSim.func) {
		    analysisResults.cpvCorrectionFitSim.func->DrawCopy("SAME");
		    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(analysisResults.cpvCorrectionFitSim.func, "Fit", "L");
		    if (analysisResults.cpvCorrectionFitSim.fitSigma > 0) {
			Float_t x = analysisResults.cpvCorrectionFitSim.func->GetXaxis()->GetXmax();
			Float_t y = analysisResults.cpvCorrectionFitSim.func->Eval(x);
			TH1F *histCpvFitErr = new TH1F(analysisSettingsNew.name + "_cpvfiterrSim", "fit err", 1, x - 1.5, x - 1.0);
			if (histCpvFitErr) {
			    histCpvFitErr->SetBinContent(1, y);
			    histCpvFitErr->SetBinError(1, analysisResults.cpvCorrectionFitSim.func->GetParError(0));
			    histCpvFitErr->SetMarkerColor(17);
			    histCpvFitErr->SetLineColor(17);
			    histCpvFitErr->SetFillColor(17);
			    histCpvFitErr->Draw("SAME ][ E2");
			    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(histCpvFitErr, Form("Par0 err = %.4f", analysisResults.cpvCorrectionFitSim.func->GetParError(0)), "EPF");
			}
			TH1F *histCpvMeanAbs = new TH1F(analysisSettingsNew.name + "_cpvmeanabsSim", "mean abs", 1, x - 1.0, x - 0.5);
			if (histCpvMeanAbs) {
			    histCpvMeanAbs->SetBinContent(1, y);
			    histCpvMeanAbs->SetBinError(1, analysisResults.cpvCorrectionFitSim.fitMeanAbs);
			    histCpvMeanAbs->SetMarkerColor(15);
			    histCpvMeanAbs->SetLineColor(15);
			    histCpvMeanAbs->SetFillColor(15);
			    histCpvMeanAbs->Draw("SAME ][ E2");
			    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(histCpvMeanAbs, Form("Mean |#Delta| = %.4f", analysisResults.cpvCorrectionFitSim.fitMeanAbs), "EPF");
			}
			TH1F *histCpvSigma = new TH1F(analysisSettingsNew.name + "_cpvsigmaSim", "#sigma", 1, x - 0.5, x);
			if (histCpvSigma) {
			    histCpvSigma->SetBinContent(1, y);
			    histCpvSigma->SetBinError(1, analysisResults.cpvCorrectionFitSim.fitSigma);
			    histCpvSigma->SetMarkerColor(13);
			    histCpvSigma->SetLineColor(13);
			    histCpvSigma->SetFillColor(13);
			    histCpvSigma->Draw("SAME ][ E2");
			    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(histCpvSigma, Form("RMS = %.4f", analysisResults.cpvCorrectionFitSim.fitSigma), "EPF");
			}
		    }
		}
		if (analysisResults.resultsSimMB.cpvCorrectionFit.func) {
		    analysisResults.resultsSimMB.cpvCorrectionFit.func->SetLineColor(analysisSettingsNew.settingsSimMB.color);
		    analysisResults.resultsSimMB.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(analysisResults.resultsSimMB.cpvCorrectionFit.func, "Fit MB", "L");
		}
		if (analysisResults.resultsSimHT1.cpvCorrectionFit.func) {
		    analysisResults.resultsSimHT1.cpvCorrectionFit.func->SetLineColor(analysisSettingsNew.settingsSimHT1.color);
		    analysisResults.resultsSimHT1.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(analysisResults.resultsSimHT1.cpvCorrectionFit.func, "Fit HT1", "L");
		}
		if (analysisResults.resultsSimHT2.cpvCorrectionFit.func) {
		    analysisResults.resultsSimHT2.cpvCorrectionFit.func->SetLineColor(analysisSettingsNew.settingsSimHT2.color);
		    analysisResults.resultsSimHT2.cpvCorrectionFit.func->DrawCopy("SAME");
		    if (legendCpvCorrectionSim) legendCpvCorrectionSim->AddEntry(analysisResults.resultsSimHT2.cpvCorrectionFit.func, "Fit HT2", "L");
		}
	}
	TLegend *legendCpvCorrectionParam7Sim = 0;
	TString CpvCorrectionParam7SimName = analysisSettingsNew.name + "_CpvCorrectionParam7Sim";
	TString CpvCorrectionParam7SimTitle = analysisSettingsNew.title + ", CPV correction param7 MC pi0;p_{T}, GeV/c;CPV correction Param7";
	TH1F *histCpvCorrectionParam7Sim = 0;
	histCpvCorrectionParam7Sim = showLists(&analysisResults.resultsSimMB.cpvCorrectionParam7, &analysisResults.resultsSimHT1.cpvCorrectionParam7, &analysisResults.resultsSimHT2.cpvCorrectionParam7
	    , CpvCorrectionParam7SimName.Data(), CpvCorrectionParam7SimTitle.Data(), histCpvCorrectionParam7Sim, &legendCpvCorrectionParam7Sim
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	if (histCpvCorrectionParam7Sim) {
		histCpvCorrectionParam7Sim->SetMinimum(0.0);
		histCpvCorrectionParam7Sim->SetMaximum(25.0);
		histCpvCorrectionParam7Sim->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsSimHT2.highPt);
	}

	divideBins(analysisResults.resultsMB.peakPosition, analysisResults.resultsSimMB.peakPosition, analysisResults.peakPositionDataToSimMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT1.peakPosition, analysisResults.resultsSimHT1.peakPosition, analysisResults.peakPositionDataToSimHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT2.peakPosition, analysisResults.resultsSimHT2.peakPosition, analysisResults.peakPositionDataToSimHT2, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsMB.peakWidth, analysisResults.resultsSimMB.peakWidth, analysisResults.peakWidthDataToSimMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT1.peakWidth, analysisResults.resultsSimHT1.peakWidth, analysisResults.peakWidthDataToSimHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT2.peakWidth, analysisResults.resultsSimHT2.peakWidth, analysisResults.peakWidthDataToSimHT2, 1.0, 0.0, false, false);

	divideBins(analysisResults.resultsMB.peakPositionEta, analysisResults.resultsSimEtaMB.peakPositionEta, analysisResults.peakPositionDataToSimEtaMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT1.peakPositionEta, analysisResults.resultsSimEtaHT1.peakPositionEta, analysisResults.peakPositionDataToSimEtaHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT2.peakPositionEta, analysisResults.resultsSimEtaHT2.peakPositionEta, analysisResults.peakPositionDataToSimEtaHT2, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsMB.peakWidth, analysisResults.resultsSimEtaMB.peakWidth, analysisResults.peakWidthDataToSimEtaMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT1.peakWidth, analysisResults.resultsSimEtaHT1.peakWidth, analysisResults.peakWidthDataToSimEtaHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT2.peakWidth, analysisResults.resultsSimEtaHT2.peakWidth, analysisResults.peakWidthDataToSimEtaHT2, 1.0, 0.0, false, false);

	TLegend *legendPosEmbed = 0;
	TString posEmbedName = analysisSettingsNew.name + "_position_embed";
	TString posEmbedTitle = analysisSettingsNew.title + ", Peak position and width vs. p_{T}, #pi^{0};p_{T}, GeV/c;Peak position and width, GeV/c^{2}";
	TH1F *histPosEmbed = 0;
	histPosEmbed = showLists(&analysisResults.resultsSimMB.peakPosition, &analysisResults.resultsSimHT1.peakPosition, &analysisResults.resultsSimHT2.peakPosition
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	);
	histPosEmbed = showLists(&analysisResults.resultsMB.peakPosition, &analysisResults.resultsHT1.peakPosition, &analysisResults.resultsHT2.peakPosition
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	);
	histPosEmbed = showLists(&analysisResults.resultsSimMB.peakWidth, &analysisResults.resultsSimHT1.peakWidth, &analysisResults.resultsSimHT2.peakWidth
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	);
	histPosEmbed = showLists(&analysisResults.resultsMB.peakWidth, &analysisResults.resultsHT1.peakWidth, &analysisResults.resultsHT2.peakWidth
	    , posEmbedName.Data(), posEmbedTitle.Data(), histPosEmbed, &legendPosEmbed
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
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
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	);
	histPosEmbedEta = showLists(&analysisResults.resultsMB.peakPositionEta, &analysisResults.resultsHT1.peakPositionEta, &analysisResults.resultsHT2.peakPositionEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	);
	histPosEmbedEta = showLists(&analysisResults.resultsSimEtaMB.peakWidthEta, &analysisResults.resultsSimEtaHT1.peakWidthEta, &analysisResults.resultsSimEtaHT2.peakWidthEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	);
	histPosEmbedEta = showLists(&analysisResults.resultsMB.peakWidthEta, &analysisResults.resultsHT1.peakWidthEta, &analysisResults.resultsHT2.peakWidthEta
	    , posEmbedEtaName.Data(), posEmbedEtaTitle.Data(), histPosEmbedEta, &legendPosEmbedEta
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
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

	TLegend *legendPosDiv = 0;
	TString posDivName = analysisSettingsNew.name + "_position_div";
	TString posDivTitle = analysisSettingsNew.title + ", Peak position data/sim;p_{T}, GeV/c;Peak position, data/sim";
	TH1F *histPosDiv = 0;
	histPosDiv = showLists(&analysisResults.peakPositionDataToSimMB, &analysisResults.peakPositionDataToSimHT1, &analysisResults.peakPositionDataToSimHT2
	    , posDivName.Data(), posDivTitle.Data(), histPosDiv, &legendPosDiv
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	);
	histPosDiv = showLists(&analysisResults.peakPositionDataToSimEtaMB, &analysisResults.peakPositionDataToSimEtaHT1, &analysisResults.peakPositionDataToSimEtaHT2
	    , posDivName.Data(), posDivTitle.Data(), histPosDiv, &legendPosDiv
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	);
	if (histPosDiv) {
		histPosDiv->SetMinimum(0);
		histPosDiv->SetMaximum(2.0);
		histPosDiv->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	{
	    bin_stat_list_type temp = analysisResults.peakPositionDataToSimMB;
	    mergePoints(analysisResults.peakPositionDataToSimHT1, temp);
	    mergePoints(analysisResults.peakPositionDataToSimHT2, temp);
	    bin_stat_list_type tempEta = analysisResults.peakPositionDataToSimEtaMB;
	    mergePoints(analysisResults.peakPositionDataToSimEtaHT1, tempEta);
	    mergePoints(analysisResults.peakPositionDataToSimEtaHT2, tempEta);
	    bin_stat_list_type tempPi0Eta = temp;
	    mergePoints(tempEta, tempPi0Eta);
	    TLegend *legendPosDivMerged = 0;
	    TString posDivMergedName = analysisSettingsNew.name + "_position_div_merged";
	    TString posDivMergedTitle = analysisSettingsNew.title + ", Peak position data/sim merged;p_{T}, GeV/c;Peak position, data/sim";
	    TH1F *histPosDivMerged = 0;
	    histPosDivMerged = showList(&temp
		, posDivMergedName.Data(), posDivMergedTitle.Data(), histPosDivMerged, &legendPosDivMerged
		, analysisSettingsNew.settingsMB.drawOptions
	    );
	    histPosDivMerged = showList(&tempEta
		, posDivMergedName.Data(), posDivMergedTitle.Data(), histPosDivMerged, &legendPosDivMerged
		, analysisSettingsNew.settingsMB.drawOptionsEta
	    );
	    histPosDivMerged = showList(&tempPi0Eta
		, posDivMergedName.Data(), posDivMergedTitle.Data(), histPosDivMerged, &legendPosDivMerged
		, analysisSettingsNew.settingsSimMB.drawOptions
	    );
	    if (histPosDivMerged) {
		    histPosDivMerged->SetMinimum(0);
		    histPosDivMerged->SetMaximum(2.0);
		    histPosDivMerged->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	    }
	}

	/*{
	    bin_stat_list_type temp = analysisResults.peakPositionDataToSimMB;
	    temp.insert(temp.end(), analysisResults.peakPositionDataToSimHT1.begin(), analysisResults.peakPositionDataToSimHT1.end());
	    temp.insert(temp.end(), analysisResults.peakPositionDataToSimHT2.begin(), analysisResults.peakPositionDataToSimHT2.end());
	    bin_stat_list_type tempEta = analysisResults.peakPositionDataToSimEtaMB;
	    tempEta.insert(tempEta.end(), analysisResults.peakPositionDataToSimEtaHT1.begin(), analysisResults.peakPositionDataToSimEtaHT1.end());
	    tempEta.insert(tempEta.end(), analysisResults.peakPositionDataToSimEtaHT2.begin(), analysisResults.peakPositionDataToSimEtaHT2.end());
	    bin_stat_list_type tempPi0Eta = temp;
	    tempPi0Eta.insert(tempPi0Eta.end(), tempEta.begin(), tempEta.end());
	    TLegend *legendPosDivAdded = 0;
	    TString posDivAddedName = analysisSettingsNew.name + "_position_added";
	    TString posDivAddedTitle = analysisSettingsNew.title + ", Peak position data/sim added;p_{T}, GeV/c;Peak position, data/sim";
	    TH1F *histPosDivAdded = 0;
	    histPosDivAdded = showList(&temp
		, posDivAddedName.Data(), posDivAddedTitle.Data(), histPosDivAdded, &legendPosDivAdded
		, analysisSettingsNew.settingsMB.drawOptions
	    );
	    histPosDivAdded = showList(&tempEta
		, posDivAddedName.Data(), posDivAddedTitle.Data(), histPosDivAdded, &legendPosDivAdded
		, analysisSettingsNew.settingsMB.drawOptionsEta
	    );
	    histPosDivAdded = showList(&tempPi0Eta
		, posDivAddedName.Data(), posDivAddedTitle.Data(), histPosDivAdded, &legendPosDivAdded
		, analysisSettingsNew.settingsSimMB.drawOptions
	    );
	    if (histPosDivAdded) {
		    histPosDivAdded->SetMinimum(0);
		    histPosDivAdded->SetMaximum(2.0);
		    histPosDivAdded->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	    }
	}*/

	TLegend *legendPosWDiv = 0;
	TString posWDivName = analysisSettingsNew.name + "_width_sim_div";
	TString posWDivTitle = analysisSettingsNew.title + ", Peak width data/sim;p_{T}, GeV/c;Peak width, data/sim";
	TH1F *histPosWDiv = 0;
	histPosWDiv = showLists(&analysisResults.peakWidthDataToSimMB, &analysisResults.peakWidthDataToSimHT1, &analysisResults.peakWidthDataToSimHT2
	    , posWDivName.Data(), posWDivTitle.Data(), histPosWDiv, &legendPosWDiv
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	);
	histPosWDiv = showLists(&analysisResults.peakWidthDataToSimEtaMB, &analysisResults.peakWidthDataToSimEtaHT1, &analysisResults.peakWidthDataToSimEtaHT2
	    , posWDivName.Data(), posWDivTitle.Data(), histPosWDiv, &legendPosWDiv
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	);
	if (histPosWDiv) {
		histPosWDiv->SetMinimum(0);
		histPosWDiv->SetMaximum(2.0);
		histPosWDiv->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	TLegend *legendYieldEmbed = 0;
	TString yieldEmbedName = analysisSettingsNew.name + "_yield_embed";
	TString yieldEmbedTitle = analysisSettingsNew.title + ", Raw yield vs. p_{T} - simulation;p_{T}, GeV/c;Yield";
	TH1F *histYieldEmbed = 0;
	histYieldEmbed = showLists(&analysisResults.areaSimuMB, &analysisResults.areaSimuHT1, &analysisResults.areaSimuHT2
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	);
	histYieldEmbed = showLists(&analysisResults.resultsSimMB.rawYield, &analysisResults.resultsSimHT1.rawYield, &analysisResults.resultsSimHT2.rawYield
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	);
	histYieldEmbed = showLists(&analysisResults.areaSimuEtaMB, &analysisResults.areaSimuEtaHT1, &analysisResults.areaSimuEtaHT2
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	);
	histYieldEmbed = showLists(&analysisResults.resultsSimEtaMB.rawYieldEta, &analysisResults.resultsSimEtaHT1.rawYieldEta, &analysisResults.resultsSimEtaHT2.rawYieldEta
	    , yieldEmbedName.Data(), yieldEmbedTitle.Data(), histYieldEmbed, &legendYieldEmbed
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
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
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	);
	histYield = showLists(&analysisResults.resultsMB.rawYieldEta, &analysisResults.resultsHT1.rawYieldEta, &analysisResults.resultsHT2.rawYieldEta
	    , yieldName.Data(), yieldTitle.Data(), histYield, &legendYield
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
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
	bin_stat_list_type yieldRatioMB, yieldRatioHT1, yieldRatioHT2, yieldRatioEtaMB, yieldRatioEtaHT1, yieldRatioEtaHT2;
	divideBins(analysisResults.resultsMB.rawYield, analysisResults.resultsSimMB.rawYield, yieldRatioMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT1.rawYield, analysisResults.resultsSimHT1.rawYield, yieldRatioHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT2.rawYield, analysisResults.resultsSimHT2.rawYield, yieldRatioHT2, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsMB.rawYieldEta, analysisResults.resultsSimEtaMB.rawYieldEta, yieldRatioEtaMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT1.rawYieldEta, analysisResults.resultsSimEtaHT1.rawYieldEta, yieldRatioEtaHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.resultsHT2.rawYieldEta, analysisResults.resultsSimEtaHT2.rawYieldEta, yieldRatioEtaHT2, 1.0, 0.0, false, false);
	TLegend *legendYieldRatio = 0;
	TString yieldRatioName = analysisSettingsNew.name + "_yieldRatio";
	TString yieldRatioTitle = analysisSettingsNew.title + ", Raw yield ratio vs. p_{T} - data/simulation, #pi^{0} and #eta;p_{T}, GeV/c;Yield ratio";
	TH1F *histYieldRatio = 0;
	histYieldRatio = showLists(&yieldRatioMB, &yieldRatioHT1, &yieldRatioHT2
	    , yieldRatioName.Data(), yieldRatioTitle.Data(), histYieldRatio, &legendYieldRatio
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	);
	histYieldRatio = showLists(&yieldRatioEtaMB, &yieldRatioEtaHT1, &yieldRatioEtaHT2
	    , yieldRatioName.Data(), yieldRatioTitle.Data(), histYieldRatio, &legendYieldRatio
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	);

	divideBins(analysisResults.resultsSimMB.rawYield, analysisResults.effMBSmooth, analysisResults.spectrumMBembed, 0.0, 
TMath::TwoPi() * analysisResults.resultsSimMB.evNum * analysisResults.resultsSimMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsSimMB.deltaEta 
* (1.0 / analysisResults.resultsSimMB.vertexFindingEff)
, true, false);
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

	if (analysisSettingsNew.multiplicityDistributions || false) {
		analysisResults.spectrumMB = analysisResults.resultsMB.rawYield;
		analysisResults.spectrumHT1 = analysisResults.resultsHT1.rawYield;
		analysisResults.spectrumHT2 = analysisResults.resultsHT2.rawYield;
	} else {
		if (analysisSettingsNew.dontCorrectEff && false) {
			divideBinsFunc(analysisResults.resultsMB.rawYield, &funcOne, analysisResults.spectrumMB, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
			divideBinsFunc(analysisResults.resultsHT1.rawYield, &funcOne, analysisResults.spectrumHT1, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
			divideBinsFunc(analysisResults.resultsHT2.rawYield, &funcOne, analysisResults.spectrumHT2, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		} else {
		    if (analysisResults.effMBSmooth.size()) {
			divideBins(analysisResults.resultsMB.rawYield, analysisResults.effMBSmooth, analysisResults.spectrumMB, 0.0, 
TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsMB.deltaEta
 * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim)
, true, false);
		    } else {
			divideBinsFunc(analysisResults.resultsMB.rawYield, &funcOne, analysisResults.spectrumMB, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * truePionBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		    }
		    if (analysisResults.effHT1Smooth.size()) {
			divideBins(analysisResults.resultsHT1.rawYield, analysisResults.effHT1Smooth, analysisResults.spectrumHT1, 0.0, 
TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT1.deltaEta
 * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false);
		    } else {
			divideBinsFunc(analysisResults.resultsHT1.rawYield, &funcOne, analysisResults.spectrumHT1, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		    }
		    if (analysisResults.effHT2Smooth.size()) {
			divideBins(analysisResults.resultsHT2.rawYield, analysisResults.effHT2Smooth, analysisResults.spectrumHT2, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false);
		    } else {
			divideBinsFunc(analysisResults.resultsHT2.rawYield, &funcOne, analysisResults.spectrumHT2, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * truePionBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		    }
		}
	}

	if (analysisSettingsNew.multiplicityDistributions || false) {
		analysisResults.spectrumMBEta = analysisResults.resultsMB.rawYieldEta;
		analysisResults.spectrumHT1Eta = analysisResults.resultsHT1.rawYieldEta;
		analysisResults.spectrumHT2Eta = analysisResults.resultsHT2.rawYieldEta;
	} else {
		if (analysisSettingsNew.dontCorrectEffEta && false) {
			divideBinsFunc(analysisResults.resultsMB.rawYieldEta, &funcOne, analysisResults.spectrumMBEta, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
			divideBinsFunc(analysisResults.resultsHT1.rawYieldEta, &funcOne, analysisResults.spectrumHT1Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
			divideBinsFunc(analysisResults.resultsHT2.rawYieldEta, &funcOne, analysisResults.spectrumHT2Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		} else {
		    if (analysisResults.effMBSmoothEta.size()) {
			divideBins(analysisResults.resultsMB.rawYieldEta, analysisResults.effMBSmoothEta, analysisResults.spectrumMBEta, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false);
		    } else {
			divideBinsFunc(analysisResults.resultsMB.rawYieldEta, &funcOne, analysisResults.spectrumMBEta, 0.0, TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsMB.deltaEta * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		    }
		    if (analysisResults.effHT1SmoothEta.size()) {
			divideBins(analysisResults.resultsHT1.rawYieldEta, analysisResults.effHT1SmoothEta, analysisResults.spectrumHT1Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false);
		    } else {
			divideBinsFunc(analysisResults.resultsHT1.rawYieldEta, &funcOne, analysisResults.spectrumHT1Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT1.deltaEta * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		    }
		    if (analysisResults.effHT2SmoothEta.size()) {
			divideBins(analysisResults.resultsHT2.rawYieldEta, analysisResults.effHT2SmoothEta, analysisResults.spectrumHT2Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false);
		    } else {
			divideBinsFunc(analysisResults.resultsHT2.rawYieldEta, &funcOne, analysisResults.spectrumHT2Eta, 0.0, TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * trueEtaBranchingRatio * analysisResults.resultsHT2.deltaEta * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim), true, false, false);
		    }
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

	{
	bin_stat_list_type temp;
	cropPoints(analysisResults.spectrumMB, temp, analysisSettingsNew.settingsMB.lowPtUse, analysisSettingsNew.settingsMB.highPtUse);
	mergePoints(temp, analysisResults.spectrum);
	}
	{
	bin_stat_list_type temp;
	cropPoints(analysisResults.spectrumHT1, temp, analysisSettingsNew.settingsHT1.lowPtUse, analysisSettingsNew.settingsHT1.highPtUse);
	mergePoints(temp, analysisResults.spectrum);
	}
	{
	bin_stat_list_type temp;
	cropPoints(analysisResults.spectrumHT2, temp, analysisSettingsNew.settingsHT2.lowPtUse, analysisSettingsNew.settingsHT2.highPtUse);
	mergePoints(temp, analysisResults.spectrum);
	}

	{
	bin_stat_list_type temp;
	cropPoints(analysisResults.spectrumMBEta, temp, analysisSettingsNew.settingsMB.lowPtUse, analysisSettingsNew.settingsMB.highPtUse);
	mergePoints(temp, analysisResults.spectrumEta);
	}
	{
	bin_stat_list_type temp;
	cropPoints(analysisResults.spectrumHT1Eta, temp, analysisSettingsNew.settingsHT1.lowPtUse, analysisSettingsNew.settingsHT1.highPtUse);
	mergePoints(temp, analysisResults.spectrumEta);
	}
	{
	bin_stat_list_type temp;
	cropPoints(analysisResults.spectrumHT2Eta, temp, analysisSettingsNew.settingsHT2.lowPtUse, analysisSettingsNew.settingsHT2.highPtUse);
	mergePoints(temp, analysisResults.spectrumEta);
	}


	bin_stat_list_type tempMB;
	divideBins(analysisResults.resultsMB.rawYieldGamma, analysisResults.effGammaMB, tempMB, 0.0, 
(TMath::TwoPi() * analysisResults.resultsMB.evNum * analysisResults.resultsMB.prescaleUsed * analysisResults.resultsMB.deltaEta
 * (1.0 / analysisResults.resultsMB.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim)), false, false);
	bin_stat_list_type tempHT1;
	divideBins(analysisResults.resultsHT1.rawYieldGamma, analysisResults.effGammaHT1, tempHT1, 0.0, 
(TMath::TwoPi() * analysisResults.resultsHT1.evNum * analysisResults.resultsHT1.prescaleUsed * analysisResults.resultsHT1.deltaEta
 * (1.0 / analysisResults.resultsHT1.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim)), false, false);
	bin_stat_list_type tempHT2;
	divideBins(analysisResults.resultsHT2.rawYieldGamma, analysisResults.effGammaHT2, tempHT2, 0.0, 
(TMath::TwoPi() * analysisResults.resultsHT2.evNum * analysisResults.resultsHT2.prescaleUsed * analysisResults.resultsHT2.deltaEta
 * (1.0 / analysisResults.resultsHT2.vertexFindingEff) * (1.0 / analysisSettingsNew.materialCorrectionSim)), false, false);
	analysisResults.invYieldGammaMB = tempMB
		- (analysisResults.spectrumMB * (
			(analysisResults.gammaToPi0DecayMB * ((analysisResults.effGammaPi0MB / analysisResults.effGammaMB) - 1.0)) +
			(0.5 * analysisResults.gammaToEtaDecayMB * ((analysisResults.effGammaEtaMB / analysisResults.effGammaMB) - 1.0))
// here nbar
		));
	analysisResults.invYieldGammaHT1 = tempHT1
		- (analysisResults.spectrumHT1 * (
			(analysisResults.gammaToPi0DecayHT1 * ((analysisResults.effGammaPi0HT1 / analysisResults.effGammaHT1) - 1.0)) +
			(0.5 * analysisResults.gammaToEtaDecayHT1 * ((analysisResults.effGammaEtaHT1 / analysisResults.effGammaHT1) - 1.0))
// here nbar
		));
	analysisResults.invYieldGammaHT2 = tempHT2
		- (analysisResults.spectrumHT2 * (
			(analysisResults.gammaToPi0DecayHT2 * ((analysisResults.effGammaPi0HT2 / analysisResults.effGammaHT2) - 1.0)) +
			(0.5 * analysisResults.gammaToEtaDecayHT2 * ((analysisResults.effGammaEtaHT2 / analysisResults.effGammaHT2) - 1.0))
// here nbar
		));
/*
	{
	bin_stat_list_type tempMB = ((analysisResults.effGammaPi0MB / analysisResults.effGammaMB) - 1.0);
	bin_stat_list_type tempHT1 = ((analysisResults.effGammaPi0HT1 / analysisResults.effGammaHT1) - 1.0);
	bin_stat_list_type tempHT2 = ((analysisResults.effGammaPi0HT2 / analysisResults.effGammaHT2) - 1.0);
	bin_stat_list_type tempMBEta = ((analysisResults.effGammaEtaMB / analysisResults.effGammaMB) - 1.0);
	bin_stat_list_type tempHT1Eta = ((analysisResults.effGammaEtaHT1 / analysisResults.effGammaHT1) - 1.0);
	bin_stat_list_type tempHT2Eta = ((analysisResults.effGammaEtaHT2 / analysisResults.effGammaHT2) - 1.0);
	TLegend *legend = 0;
	TString EnergyScaleSystErrName = analysisSettingsNew.name + "_temp";
	TString EnergyScaleSystErrTitle = analysisSettingsNew.title + ", TEMP";
	TH1F *histEnergyScaleSystErr = 0;
	histEnergyScaleSystErr = showLists(&tempMB, &tempHT1, &tempHT2
	    , EnergyScaleSystErrName.Data(), EnergyScaleSystErrTitle.Data(), histEnergyScaleSystErr, &legend
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histEnergyScaleSystErr = showLists(&tempMBEta, &tempHT1Eta, &tempHT2Eta
	    , EnergyScaleSystErrName.Data(), EnergyScaleSystErrTitle.Data(), histEnergyScaleSystErr, &legend
	    , analysisSettingsNew.settingsMB.drawOptionsEta, analysisSettingsNew.settingsHT1.drawOptionsEta, analysisSettingsNew.settingsHT2.drawOptionsEta
	    );
	}
*/
        mergePoints(analysisResults.invYieldGammaMB, analysisResults.invYieldGamma);
        mergePoints(analysisResults.invYieldGammaHT1, analysisResults.invYieldGamma);
        mergePoints(analysisResults.invYieldGammaHT2, analysisResults.invYieldGamma);

	divideBins(analysisResults.invYieldGammaMB, analysisResults.spectrumMB, analysisResults.gammaToPi0InclMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.invYieldGammaHT1, analysisResults.spectrumHT1, analysisResults.gammaToPi0InclHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.invYieldGammaHT2, analysisResults.spectrumHT2, analysisResults.gammaToPi0InclHT2, 1.0, 0.0, false, false);

	TF1 decayFuncMartijn("decayFuncMartijn", "(5.166342e-01) + (-4.843337e-02*x) + (4.257745e-03*x*x) + (-1.273812e-04*x*x*x)"); // Martijn's
	decayFuncMartijn.SetRange(0, 20);
	decayFuncMartijn.SetLineStyle(kDashed);
	//TF1 decayFunc("decayFunc", "(7.735118e-01) + (-6.636997e-02*x) + (5.317628e-03*x*x) + (-1.464556e-04*x*x*x)"); // mine
	TF1 decayFunc("decayFunc", "(4.907793e-01) + (-3.421260e-02*x) + (2.738541e-03*x*x) + (-7.777598e-05*x*x*x)"); // mine with eta cut on gammas
	decayFunc.SetRange(0, 20);

	divideBinsFunc(analysisResults.gammaToPi0InclMB, &decayFunc, analysisResults.RGammaMB, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.gammaToPi0InclHT1, &decayFunc, analysisResults.RGammaHT1, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.gammaToPi0InclHT2, &decayFunc, analysisResults.RGammaHT2, 1.0, 0.0, false, false, false);

	mergePoints(analysisResults.RGammaMB, analysisResults.RGamma);
	mergePoints(analysisResults.RGammaHT1, analysisResults.RGamma);
	mergePoints(analysisResults.RGammaHT2, analysisResults.RGamma);

	{
	    bin_stat_list_type bins = analysisResults.spectrum;
	    //TF1 *funcWPtShift = 0;
	    TWeightCalculator wPtShift;
	    for (Int_t i = 0;i < 5;i++) {
//cout << "QQQ: i = " << i << endl;
/*
		wPtShift.Fit(bins, analysisSettingsNew.correctForPreciseBinPtFitOption);
		//cout << "||   ----- pt shift iter " << i << endl;
		//wPtShift.Print("||   ");
		if (funcWPtShift) delete funcWPtShift;
		funcWPtShift = 0;
		funcWPtShift = wPtShift.createFunc(true);
		createPtBinPosFuncFromQCD(funcWPtShift);
		correctPtShift(bins, true, 0);
*/
                createPtBinPosFuncFromQCD(0);
                for (bin_stat_list_type::const_iterator it = bins.begin();it != bins.end();++it) {
                    bin_stat_list_type temp_bins;
                    bin_stat_list_type::const_iterator it_prev = it;
                    for (Int_t i = 0;i < analysisSettingsNew.correctForPreciseBinPtNFitPoints;i++) {
                        if (it_prev != bins.begin()) {--it_prev; temp_bins.push_front(*it_prev);}
                    }
                    temp_bins.push_back(*it);
                    bin_stat_list_type::const_iterator it_next = it; ++it_next;
                    for (Int_t i = 0;i < analysisSettingsNew.correctForPreciseBinPtNFitPoints;i++) {
                        if (it_next != bins.end()) {temp_bins.push_back(*it_next); ++it_next;}
                    }
//cout << "QQQ: " << temp_bins.size() << endl;
                    wPtShift.Fit(temp_bins, analysisSettingsNew.correctForPreciseBinPtFitOption, analysisSettingsNew.correctForPreciseBinPtFitOptionWeight);
                    pQCDfunc = wPtShift.createFunc(true);
                    Float_t binCenter = ((*it).getParameters().max + (*it).getParameters().min) / 2.0;
                    Float_t binWidth = ((*it).getParameters().max - (*it).getParameters().min);
                    if (ptBinPosFuncFromQCD) ptBinPosFuncFromQCD->Eval(binCenter, binWidth);
                    if (ptBinPosCorrFuncFromQCD) ptBinPosCorrFuncFromQCD->Eval(binCenter, binWidth);
                    delete pQCDfunc; pQCDfunc = 0;
                }
                correctPtShift(bins, true, 0);
	    }
//cout << "AAA: " << bins.size() << endl;
//for (bin_stat_list_type::const_iterator it = bins.begin();it != bins.end();++it) {
//cout << "bin " << (*it).getParameters().min << ": " << (*it).getValue() << endl;
//}

	    if (analysisSettingsNew.correctForPreciseBinPt) {
		correctPtShift(analysisResults.spectrumMB, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimMB.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT1, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimHT1.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT2, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimHT2.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrum, analysisSettingsNew.correctForPreciseBinPtHorizontal, 0);

		correctPtShift(analysisResults.spectrumMBEta, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimEtaMB.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT1Eta, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimEtaHT1.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT2Eta, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimEtaHT2.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumEta, analysisSettingsNew.correctForPreciseBinPtHorizontal, 0);

		correctPtShift(analysisResults.spectrumMBembed, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimMB.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT1embed, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimHT1.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT2embed, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimHT2.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumembed, analysisSettingsNew.correctForPreciseBinPtHorizontal, 0);

		correctPtShift(analysisResults.spectrumMBembedEta, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimEtaMB.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT1embedEta, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimEtaHT1.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumHT2embedEta, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimEtaHT2.simulatedPtlist : 0);
		correctPtShift(analysisResults.spectrumembedEta, analysisSettingsNew.correctForPreciseBinPtHorizontal, 0);

		correctPtShift(analysisResults.invYieldGammaMB, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimMB.simulatedPtlist : 0);
		correctPtShift(analysisResults.invYieldGammaHT1, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimHT1.simulatedPtlist : 0);
		correctPtShift(analysisResults.invYieldGammaHT2, analysisSettingsNew.correctForPreciseBinPtHorizontal, analysisSettingsNew.correctForPreciseBinPtFromSim ? &analysisResults.resultsSimHT2.simulatedPtlist : 0);
		correctPtShift(analysisResults.invYieldGamma, analysisSettingsNew.correctForPreciseBinPtHorizontal, 0);
	    }
	}
        {
            calculateDerivatives(analysisResults.spectrum, analysisSettingsNew.energyScaleErr, analysisResults.energyScaleSystErr, analysisSettingsNew.energyScaleErrFitOption, analysisSettingsNew.energyScaleErrFitOptionWeight);
            calculateDerivatives(analysisResults.resultsMB.rawYield, analysisSettingsNew.energyScaleErr, analysisResults.energyScaleSystErrRawMB, analysisSettingsNew.energyScaleErrFitOption, analysisSettingsNew.energyScaleErrFitOptionWeightRaw);
            calculateDerivatives(analysisResults.resultsHT1.rawYield, analysisSettingsNew.energyScaleErr, analysisResults.energyScaleSystErrRawHT1, analysisSettingsNew.energyScaleErrFitOption, analysisSettingsNew.energyScaleErrFitOptionWeightRaw);
            calculateDerivatives(analysisResults.resultsHT2.rawYield, analysisSettingsNew.energyScaleErr, analysisResults.energyScaleSystErrRawHT2, analysisSettingsNew.energyScaleErrFitOption, analysisSettingsNew.energyScaleErrFitOptionWeightRaw);
        }
	if (analysisSettingsNew.correctForJacobian) {
		correctJacobian(analysisResults.spectrumMB, truePionMass, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsMB.deltaEta);
		correctJacobian(analysisResults.spectrumHT1, truePionMass, candidateDataProcessorHT1 ? candidateDataProcessorHT1->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorHT1 ? candidateDataProcessorHT1->getCuts().getParametersCandidate().etaLow : analysisResults.resultsHT1.deltaEta);
		correctJacobian(analysisResults.spectrumHT2, truePionMass, candidateDataProcessorHT2 ? candidateDataProcessorHT2->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorHT2 ? candidateDataProcessorHT2->getCuts().getParametersCandidate().etaLow : analysisResults.resultsHT2.deltaEta);
		correctJacobian(analysisResults.spectrum, truePionMass, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsMB.deltaEta);

		correctJacobian(analysisResults.spectrumMBEta, trueEtaMass, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsMB.deltaEta);
		correctJacobian(analysisResults.spectrumHT1Eta, trueEtaMass, candidateDataProcessorHT1 ? candidateDataProcessorHT1->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorHT1 ? candidateDataProcessorHT1->getCuts().getParametersCandidate().etaLow : analysisResults.resultsHT1.deltaEta);
		correctJacobian(analysisResults.spectrumHT2Eta, trueEtaMass, candidateDataProcessorHT2 ? candidateDataProcessorHT2->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorHT2 ? candidateDataProcessorHT2->getCuts().getParametersCandidate().etaLow : analysisResults.resultsHT2.deltaEta);
		correctJacobian(analysisResults.spectrumEta, trueEtaMass, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorMB ? candidateDataProcessorMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsMB.deltaEta);

		correctJacobian(analysisResults.spectrumMBembed, truePionMass, candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimMB.deltaEta);
		correctJacobian(analysisResults.spectrumHT1embed, truePionMass, candidateDataProcessorEmbedHT1 ? candidateDataProcessorEmbedHT1->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedHT1 ? candidateDataProcessorEmbedHT1->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimHT1.deltaEta);
		correctJacobian(analysisResults.spectrumHT2embed, truePionMass, candidateDataProcessorEmbedHT2 ? candidateDataProcessorEmbedHT2->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedHT2 ? candidateDataProcessorEmbedHT2->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimHT2.deltaEta);
		correctJacobian(analysisResults.spectrumembed, truePionMass, candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedMB ? candidateDataProcessorEmbedMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimMB.deltaEta);

		correctJacobian(analysisResults.spectrumMBembedEta, trueEtaMass, candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimEtaMB.deltaEta);
		correctJacobian(analysisResults.spectrumHT1embedEta, trueEtaMass, candidateDataProcessorEmbedEtaHT1 ? candidateDataProcessorEmbedEtaHT1->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedEtaHT1 ? candidateDataProcessorEmbedEtaHT1->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimEtaHT1.deltaEta);
		correctJacobian(analysisResults.spectrumHT2embedEta, trueEtaMass, candidateDataProcessorEmbedEtaHT2 ? candidateDataProcessorEmbedEtaHT2->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedEtaHT2 ? candidateDataProcessorEmbedEtaHT2->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimEtaHT2.deltaEta);
		correctJacobian(analysisResults.spectrumembedEta, trueEtaMass, candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getCuts().getParametersCandidate().etaHigh : 0, candidateDataProcessorEmbedEtaMB ? candidateDataProcessorEmbedEtaMB->getCuts().getParametersCandidate().etaLow : analysisResults.resultsSimEtaMB.deltaEta);
	}

	TLegend *legendEnergyScaleSystErr = 0;
	TString EnergyScaleSystErrName = analysisSettingsNew.name + "_energyscalesysterr";
	TString EnergyScaleSystErrTitle = analysisSettingsNew.title + ", energy scale syst err;p_{T}, GeV/c;syst err";
	TH1F *histEnergyScaleSystErr = 0;
	histEnergyScaleSystErr = showLists(&analysisResults.energyScaleSystErrRawMB, &analysisResults.energyScaleSystErrRawHT1, &analysisResults.energyScaleSystErrRawHT2
	    , EnergyScaleSystErrName.Data(), EnergyScaleSystErrTitle.Data(), histEnergyScaleSystErr, &legendEnergyScaleSystErr
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histEnergyScaleSystErr = showList(&analysisResults.energyScaleSystErr
	    , EnergyScaleSystErrName.Data(), EnergyScaleSystErrTitle.Data(), histEnergyScaleSystErr, &legendEnergyScaleSystErr
	    , analysisSettingsNew.settingsSimMB.drawOptions
	    );
	if (histEnergyScaleSystErr) {
	    histEnergyScaleSystErr->SetMinimum(0.0);
	    histEnergyScaleSystErr->SetMaximum(0.5);
	    histEnergyScaleSystErr->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	calculatePtShiftCorr(analysisResults.spectrumMB, analysisResults.ptShiftCorrMB, analysisResults.ptCenterShiftMB);
	calculatePtShiftCorr(analysisResults.spectrumHT1, analysisResults.ptShiftCorrHT1, analysisResults.ptCenterShiftHT1);
	calculatePtShiftCorr(analysisResults.spectrumHT2, analysisResults.ptShiftCorrHT2, analysisResults.ptCenterShiftHT2);

	TLegend *legendPtShiftCorr = 0;
	TString ptShiftCorrName = analysisSettingsNew.name + "_ptshiftcorr";
	TString ptShiftCorrTitle = analysisSettingsNew.title + ", p_{T} shift correction due to large bin width;p_{T}, GeV/c;f(nominal)/f(true)";
	TH1F *histPtShiftCorr = 0;
	histPtShiftCorr = showLists(&analysisResults.ptShiftCorrMB, &analysisResults.ptShiftCorrHT1, &analysisResults.ptShiftCorrHT2
	    , ptShiftCorrName.Data(), ptShiftCorrTitle.Data(), histPtShiftCorr, &legendPtShiftCorr
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histPtShiftCorr) {
		histPtShiftCorr->SetMinimum(0.6);
		histPtShiftCorr->SetMaximum(1.0);
		histPtShiftCorr->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	TLegend *legendPtShift = 0;
	TString ptShiftName = analysisSettingsNew.name + "_ptshift";
	TString ptShiftTitle = analysisSettingsNew.title + ", p_{T} shift due to large bin width;p_{T}, GeV/c;(true - nominal) / nominal";
	TH1F *histPtShift = 0;
	histPtShift = showLists(&analysisResults.ptCenterShiftMB, &analysisResults.ptCenterShiftHT1, &analysisResults.ptCenterShiftHT2
	    , ptShiftName.Data(), ptShiftTitle.Data(), histPtShift, &legendPtShift
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histPtShift) {
		histPtShift->SetMinimum(-0.1);
		histPtShift->SetMaximum(0.0);
		histPtShift->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

	{
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
	{
	    Float_t cxMB = (analysisSettingsNew.crossectionSimMB != 0) ? (analysisSettingsNew.crossectionSimMBTotalFraction / analysisSettingsNew.crossectionSimMB) : 1.0;
	    Float_t cxHT1 = (analysisSettingsNew.crossectionSimHT1 != 0) ? (analysisSettingsNew.crossectionSimMBTotalFraction / analysisSettingsNew.crossectionSimHT1) : 1.0;
	    Float_t cxHT2 = (analysisSettingsNew.crossectionSimHT2 != 0) ? (analysisSettingsNew.crossectionSimMBTotalFraction / analysisSettingsNew.crossectionSimHT2) : 1.0;
	    divideBinsFunc(analysisResults.spectrumMBembed, &funcOne, analysisResults.crossectionMBembed, cxMB, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT1embed, &funcOne, analysisResults.crossectionHT1embed, cxHT1, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT2embed, &funcOne, analysisResults.crossectionHT2embed, cxHT2, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumembed, &funcOne, analysisResults.crossectionembed, cxMB, 0.0, false, false, false);
	}
	{
	    Float_t cxMB = (analysisSettingsNew.crossectionSimEtaMB != 0) ? (analysisSettingsNew.crossectionSimEtaMBTotalFraction / analysisSettingsNew.crossectionSimEtaMB) : 1.0;
	    Float_t cxHT1 = (analysisSettingsNew.crossectionSimEtaHT1 != 0) ? (analysisSettingsNew.crossectionSimEtaMBTotalFraction / analysisSettingsNew.crossectionSimEtaHT1) : 1.0;
	    Float_t cxHT2 = (analysisSettingsNew.crossectionSimEtaHT2 != 0) ? (analysisSettingsNew.crossectionSimEtaMBTotalFraction / analysisSettingsNew.crossectionSimEtaHT2) : 1.0;
	    divideBinsFunc(analysisResults.spectrumMBembedEta, &funcOne, analysisResults.crossectionMBembedEta, cxMB, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT1embedEta, &funcOne, analysisResults.crossectionHT1embedEta, cxHT1, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumHT2embedEta, &funcOne, analysisResults.crossectionHT2embedEta, cxHT2, 0.0, false, false, false);
	    divideBinsFunc(analysisResults.spectrumembedEta, &funcOne, analysisResults.crossectionembedEta, cxMB, 0.0, false, false, false);
	}

	divideBins(analysisResults.spectrumEta, analysisResults.spectrum, analysisResults.etaToPi0Ratio, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumMBEta, analysisResults.spectrumMB, analysisResults.etaToPi0RatioMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT1Eta, analysisResults.spectrumHT1, analysisResults.etaToPi0RatioHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT2Eta, analysisResults.spectrumHT2, analysisResults.etaToPi0RatioHT2, 1.0, 0.0, false, false);

	TLegend *legendGammaToPi0 = 0;
	TString GammaToPi0Name = analysisSettingsNew.name + "_gammaToPi0";
	TString GammaToPi0Title = analysisSettingsNew.title + ", Gamma to pi0 ratio;p_{T}, GeV/c;R_{#gamma}";
	TH1F *histGammaToPi0 = 0;
	histGammaToPi0 = showLists(&analysisResults.gammaToPi0InclMB, &analysisResults.gammaToPi0InclHT1, &analysisResults.gammaToPi0InclHT2
	    , GammaToPi0Name.Data(), GammaToPi0Title.Data(), histGammaToPi0, &legendGammaToPi0
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	histGammaToPi0 = showLists(&analysisResults.gammaToPi0DecayMB, &analysisResults.gammaToPi0DecayHT1, &analysisResults.gammaToPi0DecayHT2
	    , GammaToPi0Name.Data(), GammaToPi0Title.Data(), histGammaToPi0, &legendGammaToPi0
	    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	    );
	histGammaToPi0 = showLists(&analysisResults.gammaToEtaDecayMB, &analysisResults.gammaToEtaDecayHT1, &analysisResults.gammaToEtaDecayHT2
	    , GammaToPi0Name.Data(), GammaToPi0Title.Data(), histGammaToPi0, &legendGammaToPi0
	    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	    );
	if (histGammaToPi0) {
		histGammaToPi0->SetMinimum(0);
		histGammaToPi0->SetMaximum(1.0);
		histGammaToPi0->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		decayFunc.DrawCopy("same");
		decayFuncMartijn.DrawCopy("same");
	}

	TLegend *legendRGamma = 0;
	TString rGammaName = analysisSettingsNew.name + "_rgamma";
	TString rGammaTitle = analysisSettingsNew.title + ", Double ratio R_{#gamma};p_{T}, GeV/c;R_{#gamma}";
	TH1F *histRGamma = 0;
	histRGamma = showLists(&analysisResults.RGammaMB, &analysisResults.RGammaHT1, &analysisResults.RGammaHT2
	    , rGammaName.Data(), rGammaTitle.Data(), histRGamma, &legendRGamma
	    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	    );
	if (histRGamma) {
		histRGamma->SetMinimum(0);
		histRGamma->SetMaximum(3.0);
		histRGamma->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
	}

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
	if (analysisSettingsNew.showWeightUsed) {
            const Char_t *showpWeightLegend = "Weight";
            TF1 *fWeight = weightUsed.createFunc(true);
	    if (wGood) {
    		TF1 *fWeightNext = w.createFunc(true);
		Float_t c = fWeightNext->Eval(5.0) / fWeight->Eval(5.0);
		delete fWeightNext;
		fWeight->SetParameter(0, c * fWeight->GetParameter(0));
		fWeight->SetParameter(1, c * fWeight->GetParameter(1));
	    }
            analysisSettingsNew.dataPointsInvYield.push_back(TDataPoints(TDrawOptions(TNamed("WEIGHT_used" , showpWeightLegend), TAttLine(kMagenta, kSolid, 1), TAttFill(kMagenta, 0), TAttMarker(kMagenta, kOpenSquare, 1.0), "L", "L", TString(showpWeightLegend)), fWeight));
	}
	if (analysisSettingsNew.showWeightUsed) {
            const Char_t *showpWeightEtaLegend = "Weight #eta";
            TF1 *fWeightEta = weightUsedEta.createFunc(true);
	    if (wGood) {
    		TF1 *fWeightEtaNext = wEta.createFunc(true);
		Float_t c = fWeightEtaNext->Eval(5.0) / fWeightEta->Eval(5.0);
		delete fWeightEtaNext;
		fWeightEta->SetParameter(0, c * fWeightEta->GetParameter(0));
		fWeightEta->SetParameter(1, c * fWeightEta->GetParameter(1));
	    }
            analysisSettingsNew.dataPointsInvYield.push_back(TDataPoints(TDrawOptions(TNamed("WEIGHT_used_eta" , showpWeightEtaLegend), TAttLine(kBlue, kSolid, 1), TAttFill(kBlue, 0), TAttMarker(kBlue, kOpenSquare, 1.0), "L", "L", TString(showpWeightEtaLegend)), fWeightEta));
	}
	if (analysisSettingsNew.showWeightNextIteration && wGood) {
            const Char_t *showpWeightLegend = "Weight, next iteration";
            TF1 *fWeight = w.createFunc(true);
            analysisSettingsNew.dataPointsInvYield.push_back(TDataPoints(TDrawOptions(TNamed("WEIGHT_next" , showpWeightLegend), TAttLine(kMagenta, kDashed, 1), TAttFill(kMagenta, 0), TAttMarker(kMagenta, kOpenSquare, 1.0), "L", "L", TString(showpWeightLegend)), fWeight));
	}
	if (analysisSettingsNew.showWeightNextIteration && wEtaGood) {
            const Char_t *showpWeightLegend = "Weight, next iteration #eta";
            TF1 *fWeight = wEta.createFunc(true);
            analysisSettingsNew.dataPointsInvYield.push_back(TDataPoints(TDrawOptions(TNamed("WEIGHT_next_eta" , showpWeightLegend), TAttLine(kBlue, kDashed, 1), TAttFill(kBlue, 0), TAttMarker(kBlue, kOpenSquare, 1.0), "L", "L", TString(showpWeightLegend)), fWeight));
	}

	TString spectrumname = analysisSettingsNew.name + "_invyield";
	TString spectrumtitle = analysisSettingsNew.title + ", Invariant yield;p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}Y}{d#etadp_{T}}";
	TLegend *legendSpectrum = 0;
	TH1F *histSpectrum = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		histSpectrum = showLists(&analysisResults.spectrumMB, &analysisResults.spectrumHT1, &analysisResults.spectrumHT2
		    , spectrumname.Data(), spectrumtitle.Data(), histSpectrum, &legendSpectrum
		    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
		);
		histSpectrum = showLists(&analysisResults.spectrumMBembed, &analysisResults.spectrumHT1embed, &analysisResults.spectrumHT2embed
		    , spectrumname.Data(), spectrumtitle.Data(), histSpectrum, &legendSpectrum
		    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
		);
	} else {
		{
		TDrawOptions temp = analysisSettingsNew.drawOptions;
		temp.legendTitle = "#pi^{0}";
		histSpectrum = showList(&analysisResults.spectrum
		    , spectrumname.Data(), spectrumtitle.Data(), histSpectrum, &legendSpectrum
		    , temp
		);
		}
		{
		TDrawOptions temp = analysisSettingsNew.drawOptions;
		temp.legendTitle = "#pi^{0} MC";
		histSpectrum = showList(&analysisResults.spectrumembed
	    	    , spectrumname.Data(), spectrumtitle.Data(), histSpectrum, &legendSpectrum
	    	    , temp
		);
		}
	}
	showResultsOthers(histSpectrum, legendSpectrum, analysisSettingsNew.dataPointsInvYield);
	if (histSpectrum) {
		histSpectrum->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histSpectrum->SetMinimum(1e-10);
		histSpectrum->SetMaximum(1e+01);
		gPad->SetLogy();
	}

	TString spectrumGammaname = analysisSettingsNew.name + "_invyieldGamma";
	TString spectrumGammatitle = analysisSettingsNew.title + ", Invariant yield of inclusive #gamma;p_{T}, GeV/c;#frac{1}{2#pip_{T}N} #frac{d^{2}Y}{d#etadp_{T}}";
	TLegend *legendSpectrumGamma = 0;
	TH1F *histSpectrumGamma = 0;
	histSpectrumGamma = showLists(&analysisResults.invYieldGammaMB, &analysisResults.invYieldGammaHT1, &analysisResults.invYieldGammaHT2
		    , spectrumGammaname.Data(), spectrumGammatitle.Data(), histSpectrumGamma, &legendSpectrumGamma
		    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
	);
	showResultsOthers(histSpectrumGamma, legendSpectrumGamma, analysisSettingsNew.dataPointsInvYield);
	if (histSpectrumGamma) {
		histSpectrumGamma->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histSpectrumGamma->SetMinimum(1e-10);
		histSpectrumGamma->SetMaximum(1e+01);
		gPad->SetLogy();
	}

	TString spectrumcname = analysisSettingsNew.name + "_crossection";
	TString spectrumctitle = analysisSettingsNew.title + ", Invariant cross section;p_{T}, GeV/c;E#frac{d^{3}#sigma}{d^{3}p}, mb*GeV^{-2}*c^{3}";
	TLegend *legendSpectrumc = 0;
	TH1F *histSpectrumc = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		histSpectrumc = showLists(&analysisResults.crossectionMB, &analysisResults.crossectionHT1, &analysisResults.crossectionHT2
		    , spectrumcname.Data(), spectrumctitle.Data(), histSpectrumc, &legendSpectrumc
		    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
		);
		histSpectrumc = showLists(&analysisResults.crossectionMBembed, &analysisResults.crossectionHT1embed, &analysisResults.crossectionHT2embed
		    , spectrumcname.Data(), spectrumctitle.Data(), histSpectrumc, &legendSpectrumc
		    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
		);
	} else {
		TDrawOptions temp = analysisSettingsNew.drawOptions;
		temp.legendTitle = "STAR d+Au #pi^{0}";
		histSpectrumc = showList(&analysisResults.crossection
		    , spectrumcname.Data(), spectrumctitle.Data(), histSpectrumc, &legendSpectrumc
		    , temp
		);
	}
	showResultsOthers(histSpectrumc, legendSpectrumc, analysisSettingsNew.dataPointsCrossSection);
	if (histSpectrumc) {
		histSpectrumc->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histSpectrumc->SetMinimum(1e-09);
		histSpectrumc->SetMaximum(1e+03);
		gPad->SetLogy();
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
		    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
		);
	} else {
		TDrawOptions temp = analysisSettingsNew.drawOptions;
		temp.legendTitle = "STAR d+Au #eta / #pi^{0}";
		histEtaToPi0 = showList(&analysisResults.etaToPi0Ratio
		    , etaToPi0name.Data(), etaToPi0title.Data(), histEtaToPi0, &legendEtaToPi0
		    , temp
		);
	}
	showResultsOthers(histEtaToPi0, legendEtaToPi0, analysisSettingsNew.dataPointsEtaToPi);
	if (histEtaToPi0) {
		histEtaToPi0->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histEtaToPi0->SetMinimum(0);
		histEtaToPi0->SetMaximum(1.2);
	}

	TF1 *weightDivSpectrum = 0;
	{
	TString weightDivSpectrumName = analysisSettingsNew.name + "_weightDivSpectrum";
	TWeightCalculator wDivSpectrum(weightDivSpectrumName.Data());
	wDivSpectrum.Fit(analysisResults.spectrum);
	weightDivSpectrum = wDivSpectrum.createFunc(true);

	divideBinsFunc(analysisResults.spectrumMB, weightDivSpectrum, analysisResults.spectrumMBDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT1, weightDivSpectrum, analysisResults.spectrumHT1Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT2, weightDivSpectrum, analysisResults.spectrumHT2Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrum, weightDivSpectrum, analysisResults.spectrumDiv, 1.0, 0.0, false, false, false);

	divideBinsFunc(analysisResults.spectrumMBembed, weightDivSpectrum, analysisResults.spectrumMBembedDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT1embed, weightDivSpectrum, analysisResults.spectrumHT1embedDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT2embed, weightDivSpectrum, analysisResults.spectrumHT2embedDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumembed, weightDivSpectrum, analysisResults.spectrumembedDiv, 1.0, 0.0, false, false, false);
	}

	TF1 *weightDivCrossection = 0;
	{
	TString weightDivCrossectionName = analysisSettingsNew.name + "_weightDivCrossection";
	TWeightCalculator wDivCrossection(weightDivCrossectionName.Data());
	wDivCrossection.Fit(analysisResults.crossection);
	weightDivCrossection = wDivCrossection.createFunc(true);

	divideBinsFunc(analysisResults.crossectionMB, weightDivCrossection, analysisResults.crossectionMBDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT1, weightDivCrossection, analysisResults.crossectionHT1Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT2, weightDivCrossection, analysisResults.crossectionHT2Div, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossection, weightDivCrossection, analysisResults.crossectionDiv, 1.0, 0.0, false, false, false);

	divideBinsFunc(analysisResults.crossectionMBembed, weightDivCrossection, analysisResults.crossectionMBembedDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT1embed, weightDivCrossection, analysisResults.crossectionHT1embedDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT2embed, weightDivCrossection, analysisResults.crossectionHT2embedDiv, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionembed, weightDivCrossection, analysisResults.crossectionembedDiv, 1.0, 0.0, false, false, false);
	}

	TF1 *weightDivSpectrumEta = 0;
	{
	TString weightDivSpectrumEtaName = analysisSettingsNew.name + "_weightDivSpectrumEta";
	TWeightCalculator wDivSpectrumEta(weightDivSpectrumEtaName.Data());
	wDivSpectrumEta.Fit(analysisResults.spectrumEta);
	weightDivSpectrumEta = wDivSpectrumEta.createFunc(true);

	divideBinsFunc(analysisResults.spectrumMBEta, weightDivSpectrumEta, analysisResults.spectrumMBDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT1Eta, weightDivSpectrumEta, analysisResults.spectrumHT1DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT2Eta, weightDivSpectrumEta, analysisResults.spectrumHT2DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumEta, weightDivSpectrumEta, analysisResults.spectrumDivEta, 1.0, 0.0, false, false, false);

	divideBinsFunc(analysisResults.spectrumMBembedEta, weightDivSpectrumEta, analysisResults.spectrumMBembedDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT1embedEta, weightDivSpectrumEta, analysisResults.spectrumHT1embedDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumHT2embedEta, weightDivSpectrumEta, analysisResults.spectrumHT2embedDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.spectrumembedEta, weightDivSpectrumEta, analysisResults.spectrumembedDivEta, 1.0, 0.0, false, false, false);
	}

	TF1 *weightDivCrossectionEta = 0;
	{
	TString weightDivCrossectionEtaName = analysisSettingsNew.name + "_weightDivCrossectionEta";
	TWeightCalculator wDivCrossectionEta(weightDivCrossectionEtaName.Data());
	wDivCrossectionEta.Fit(analysisResults.crossectionEta);
	weightDivCrossectionEta = wDivCrossectionEta.createFunc(true);

	divideBinsFunc(analysisResults.crossectionMBEta, weightDivCrossectionEta, analysisResults.crossectionMBDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT1Eta, weightDivCrossectionEta, analysisResults.crossectionHT1DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT2Eta, weightDivCrossectionEta, analysisResults.crossectionHT2DivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionEta, weightDivCrossectionEta, analysisResults.crossectionDivEta, 1.0, 0.0, false, false, false);

	divideBinsFunc(analysisResults.crossectionMBembedEta, weightDivCrossectionEta, analysisResults.crossectionMBembedDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT1embedEta, weightDivCrossectionEta, analysisResults.crossectionHT1embedDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionHT2embedEta, weightDivCrossectionEta, analysisResults.crossectionHT2embedDivEta, 1.0, 0.0, false, false, false);
	divideBinsFunc(analysisResults.crossectionembedEta, weightDivCrossectionEta, analysisResults.crossectionembedDivEta, 1.0, 0.0, false, false, false);
	}

	data_points_list dataPointsInvYieldDiv = analysisSettingsNew.dataPointsInvYield;
	if (weightDivSpectrum) dataPointsInvYieldDiv /= *weightDivSpectrum;
	data_points_list dataPointsCrossSectionDiv = analysisSettingsNew.dataPointsCrossSection;
	if (weightDivCrossection) dataPointsCrossSectionDiv /= *weightDivCrossection;

	TLegend *legendInvYieldDiv = 0;
	TString InvYieldname = analysisSettingsNew.name + "_invyield_div";
	TString InvYieldtitle = analysisSettingsNew.title + ", Invariant yield divided by fit;p_{T}, GeV/c;Data / fit";
	TH1F *histInvYieldDiv = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		histInvYieldDiv = showLists(&analysisResults.spectrumMBDiv, &analysisResults.spectrumHT1Div, &analysisResults.spectrumHT2Div
		    , InvYieldname.Data(), InvYieldtitle.Data()
		    , histInvYieldDiv, &legendInvYieldDiv
		    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
		);
		histInvYieldDiv = showLists(&analysisResults.spectrumMBembedDiv, &analysisResults.spectrumHT1embedDiv, &analysisResults.spectrumHT2embedDiv
		    , InvYieldname.Data(), InvYieldtitle.Data()
		    , histInvYieldDiv, &legendInvYieldDiv
		    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
		);
	} else {
		{
		    TDrawOptions temp = analysisSettingsNew.settingsMB.drawOptions;
		    temp.legendTitle = "STAR d+Au #pi^{0}";
		    histInvYieldDiv = showList(&analysisResults.spectrumDiv
			, InvYieldname.Data(), InvYieldtitle.Data(), histInvYieldDiv, &legendInvYieldDiv
			, temp
		    );
		}
	}
	showResultsOthers(histInvYieldDiv, legendInvYieldDiv, dataPointsInvYieldDiv);
	if (histInvYieldDiv) {
		histInvYieldDiv->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histInvYieldDiv->SetMinimum(0);
		histInvYieldDiv->SetMaximum(5);
	}

	TLegend *legendCrossSectionDiv = 0;
	TString CrossSectionname = analysisSettingsNew.name + "_crossection_div";
	TString CrossSectiontitle = analysisSettingsNew.title + ", Cross section divided by fit;p_{T}, GeV/c;Data / pQCD";
	TH1F *histCrossSectionDiv = 0;
	if (analysisSettingsNew.showTriggersSeparately) {
		histCrossSectionDiv = showLists(&analysisResults.crossectionMBDiv, &analysisResults.crossectionHT1Div, &analysisResults.crossectionHT2Div
		    , CrossSectionname.Data(), CrossSectiontitle.Data()
		    , histCrossSectionDiv, &legendCrossSectionDiv
		    , analysisSettingsNew.settingsMB.drawOptions, analysisSettingsNew.settingsHT1.drawOptions, analysisSettingsNew.settingsHT2.drawOptions
		);
		histCrossSectionDiv = showLists(&analysisResults.crossectionMBembedDiv, &analysisResults.crossectionHT1embedDiv, &analysisResults.crossectionHT2embedDiv
		    , CrossSectionname.Data(), CrossSectiontitle.Data()
		    , histCrossSectionDiv, &legendCrossSectionDiv
		    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
		);
	} else {
		{
		    TDrawOptions temp = analysisSettingsNew.settingsMB.drawOptions;
		    temp.legendTitle = "#pi^{0}";
		    histCrossSectionDiv = showList(&analysisResults.crossectionDiv
			, CrossSectionname.Data(), CrossSectiontitle.Data(), histCrossSectionDiv, &legendCrossSectionDiv
			, temp
		    );
		}
		{
		    TDrawOptions temp = analysisSettingsNew.settingsSimMB.drawOptions;
		    temp.legendTitle = "#pi^{0} MC";
		    histCrossSectionDiv = showList(&analysisResults.crossectionembedDiv
			, CrossSectionname.Data(), CrossSectiontitle.Data(), histCrossSectionDiv, &legendCrossSectionDiv
			, temp
		    );
		}
	}
	showResultsOthers(histCrossSectionDiv, legendCrossSectionDiv, dataPointsCrossSectionDiv);
	if (histCrossSectionDiv) {
		histCrossSectionDiv->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histCrossSectionDiv->SetMinimum(0);
		histCrossSectionDiv->SetMaximum(5);
	}

	{
	bin_stat_list_type yieldRatioMB, yieldRatioHT1, yieldRatioHT2;
	divideBins(analysisResults.spectrumMBDiv, analysisResults.spectrumMBembedDiv, yieldRatioMB, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT1Div, analysisResults.spectrumHT1embedDiv, yieldRatioHT1, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT2Div, analysisResults.spectrumHT2embedDiv, yieldRatioHT2, 1.0, 0.0, false, false);
	bin_stat_list_type yieldRatioMBEta, yieldRatioHT1Eta, yieldRatioHT2Eta;
	divideBins(analysisResults.spectrumMBDivEta, analysisResults.spectrumMBembedDivEta, yieldRatioMBEta, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT1DivEta, analysisResults.spectrumHT1embedDivEta, yieldRatioHT1Eta, 1.0, 0.0, false, false);
	divideBins(analysisResults.spectrumHT2DivEta, analysisResults.spectrumHT2embedDivEta, yieldRatioHT2Eta, 1.0, 0.0, false, false);
	TLegend *legendYieldRatio = 0;
	TString YieldRationame = analysisSettingsNew.name + "_YieldRatio";
	TString YieldRatiotitle = analysisSettingsNew.title + ", Invariant yield ratio Data/Sim;p_{T}, GeV/c;Data / Sim";
	TH1F *histYieldRatio = 0;
	histYieldRatio = showLists(&yieldRatioMB, &yieldRatioHT1, &yieldRatioHT2
		    , YieldRationame.Data(), YieldRatiotitle.Data()
		    , histYieldRatio, &legendYieldRatio
		    , analysisSettingsNew.settingsSimMB.drawOptions, analysisSettingsNew.settingsSimHT1.drawOptions, analysisSettingsNew.settingsSimHT2.drawOptions
	);
	histYieldRatio = showLists(&yieldRatioMBEta, &yieldRatioHT1Eta, &yieldRatioHT2Eta
		    , YieldRationame.Data(), YieldRatiotitle.Data()
		    , histYieldRatio, &legendYieldRatio
		    , analysisSettingsNew.settingsSimEtaMB.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT1.drawOptionsEta, analysisSettingsNew.settingsSimEtaHT2.drawOptionsEta
	);
	if (histYieldRatio) {
		histYieldRatio->GetXaxis()->SetLimits(0, analysisSettingsNew.settingsHT2.highPt);
		histYieldRatio->SetMinimum(0);
		histYieldRatio->SetMaximum(3);
	}
	}

	if (analysisSettingsNew.saveDataArrays) {
		TString name = analysisSettingsNew.name;
		TString title = analysisSettingsNew.title;
		TString filename = Form(analysisSettingsNew.dataArraysFilenameFormat.Data(), name.Data());
		cout << "Saving data arrays to " << filename << endl;
		ofstream ostr(filename.Data());
		ostr << "#ifndef DATA_ARRAYS_" << name << endl;
		ostr << "#define DATA_ARRAYS_" << name << endl;
		ostr << "// " << title << endl;
		TDatime startTime;
		ostr << "// Generated: " << startTime.AsSQLString() << endl;
		TString arrayName;

		arrayName = "_MBspectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumMB, "pi0 pT spectrum - MinBias", arrayName.Data());
		arrayName = "_HT1spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT1, "pi0 pT spectrum - HighTower-1", arrayName.Data());
		arrayName = "_HT2spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT2, "pi0 pT spectrum - HighTower-2", arrayName.Data());
		arrayName = "_MBcrossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionMB, "pi0 pT cross section - MinBias", arrayName.Data());
		arrayName = "_HT1crossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT1, "pi0 pT cross section - HighTower-1", arrayName.Data());
		arrayName = "_HT2crossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT2, "pi0 pT cross section - HighTower-2", arrayName.Data());
		arrayName = "_MByield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.rawYield, "Raw yield - MinBias", arrayName.Data());
		arrayName = "_HT1yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.rawYield, "Raw yield - HighTower-1", arrayName.Data());
		arrayName = "_HT2yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.rawYield, "Raw yield - HighTower-2", arrayName.Data());
		//arrayName = "_MBbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakBackground, "Background - MinBias", arrayName.Data());
		//arrayName = "_HT1bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakBackground, "Background - HighTower-1", arrayName.Data());
		//arrayName = "_HT2bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakBackground, "Background - HighTower-2", arrayName.Data());
		//arrayName = "_MBsigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakSigBackground, "Signal+Background - MinBias", arrayName.Data());
		//arrayName = "_HT1sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakSigBackground, "Signal+Background - HighTower-1", arrayName.Data());
		//arrayName = "_HT2sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakSigBackground, "Signal+Background - HighTower-2", arrayName.Data());
		arrayName = "_MBeff_"; arrayName += name; outputCArrays(ostr, analysisResults.effMBSmooth, "Correction factor - MinBias", arrayName.Data());
		arrayName = "_HT1eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT1Smooth, "Correction factor - HighTower-1", arrayName.Data());
		arrayName = "_HT2eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT2Smooth, "Correction factor - HighTower-2", arrayName.Data());
		//arrayName = "_MBsim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuMB, "Simulated yield - MinBias", arrayName.Data());
		//arrayName = "_HT1sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuHT1, "Simulated yield - HighTower-1", arrayName.Data());
		//arrayName = "_HT2sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuHT2, "Simulated yield - HighTower-2", arrayName.Data());
		//arrayName = "_MBsimreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimMB.rawYield, "Simulated reco yield - MinBias", arrayName.Data());
		//arrayName = "_HT1simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimHT1.rawYield, "Simulated reco yield - HighTower-1", arrayName.Data());
		//arrayName = "_HT2simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimHT2.rawYield, "Simulated reco yield - HighTower-2", arrayName.Data());

		arrayName = "_EtaMBspectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumMBEta, "eta pT spectrum - MinBias", arrayName.Data());
		arrayName = "_EtaHT1spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT1Eta, "eta pT spectrum - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2spectrum_"; arrayName += name; outputCArrays(ostr, analysisResults.spectrumHT2Eta, "eta pT spectrum - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBcrossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionMBEta, "eta pT cross section - MinBias", arrayName.Data());
		arrayName = "_EtaHT1crossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT1Eta, "eta pT cross section - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2crossection_"; arrayName += name; outputCArrays(ostr, analysisResults.crossectionHT2Eta, "eta pT cross section - HighTower-2", arrayName.Data());
		arrayName = "_EtaMByield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.rawYieldEta, "Raw eta yield - MinBias", arrayName.Data());
		arrayName = "_EtaHT1yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.rawYieldEta, "Raw eta yield - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2yield_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.rawYieldEta, "Raw eta yield - HighTower-2", arrayName.Data());
		//arrayName = "_EtaMBbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakBackgroundEta, "Background eta - MinBias", arrayName.Data());
		//arrayName = "_EtaHT1bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakBackgroundEta, "Background eta - HighTower-1", arrayName.Data());
		//arrayName = "_EtaHT2bg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakBackgroundEta, "Background eta - HighTower-2", arrayName.Data());
		//arrayName = "_EtaMBsigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsMB.peakSigBackgroundEta, "Signal+Background eta - MinBias", arrayName.Data());
		//arrayName = "_EtaHT1sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT1.peakSigBackgroundEta, "Signal+Background eta - HighTower-1", arrayName.Data());
		//arrayName = "_EtaHT2sigbg_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsHT2.peakSigBackgroundEta, "Signal+Background eta - HighTower-2", arrayName.Data());
		arrayName = "_EtaMBeff_"; arrayName += name; outputCArrays(ostr, analysisResults.effMBSmoothEta, "Correction factor eta - MinBias", arrayName.Data());
		arrayName = "_EtaHT1eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT1SmoothEta, "Correction factor eta - HighTower-1", arrayName.Data());
		arrayName = "_EtaHT2eff_"; arrayName += name; outputCArrays(ostr, analysisResults.effHT2SmoothEta, "Correction factor eta - HighTower-2", arrayName.Data());
		//arrayName = "_EtaMBsim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuEtaMB, "Simulated eta yield - MinBias", arrayName.Data());
		//arrayName = "_EtaHT1sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuEtaHT1, "Simulated eta yield - HighTower-1", arrayName.Data());
		//arrayName = "_EtaHT2sim_"; arrayName += name; outputCArrays(ostr, analysisResults.areaSimuEtaHT2, "Simulated eta yield - HighTower-2", arrayName.Data());
		//arrayName = "_EtaMBsimreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimEtaMB.rawYield, "Simulated reco eta yield - MinBias", arrayName.Data());
		//arrayName = "_EtaHT1simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimEtaHT1.rawYield, "Simulated reco eta yield - HighTower-1", arrayName.Data());
		//arrayName = "_EtaHT2simreco_"; arrayName += name; outputCArrays(ostr, analysisResults.resultsSimEtaHT2.rawYield, "Simulated reco eta yield - HighTower-2", arrayName.Data());

		arrayName = "_EtaToPi0MB_"; arrayName += name; outputCArrays(ostr, analysisResults.etaToPi0RatioMB, "eta/pi0 - MinBias", arrayName.Data());
		arrayName = "_EtaToPi0HT1_"; arrayName += name; outputCArrays(ostr, analysisResults.etaToPi0RatioHT1, "eta/pi0 - HighTower-1", arrayName.Data());
		arrayName = "_EtaToPi0HT2_"; arrayName += name; outputCArrays(ostr, analysisResults.etaToPi0RatioHT2, "eta/pi0 - HighTower-2", arrayName.Data());

		ostr << "#endif" << endl;
		ostr << "" << endl;
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
					ostr << "<TD bgcolor=e0e0e0>T " << ((center - width/2.0) >= 0 ? "+" : "-") << " " << TMath::Abs(center - width/2.0) << " &lt; p<SUB>T</SUB> &lt; T " << ((center + width/2.0) >= 0 ? "+" : "-") << " " << TMath::Abs(center + width/2.0) << "</TD>";
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

void saveCanvases(const Char_t *dir, Bool_t saveCR) {
    TIter iter(gROOT->GetListOfCanvases());
    TString indName = dir;
    indName += "/index.html";
    ofstream ofstr(indName.Data());
    ofstr << "<HTML>" << endl;
    ofstr << "<HEAD><TITLE>Analysis plots</TITLE></HEAD>" << endl;
    ofstr << "<BODY>" << endl;
    ofstr << "<H1>Analysis plots</H1>" << endl;
    ofstr << "<FONT size=\"-1\">Generated: " << TDatime().AsSQLString() << "</FONT><BR/>" << endl;
    ofstr << "<BR/>" << endl;
    ofstr << "<TABLE border=\"0\">" << endl;
    Int_t index = 1;
    while (TCanvas *c = dynamic_cast<TCanvas*>(iter())) {
	cout << "Canvas " << index << ": " << c->GetName() << endl;
	TString nameF = c->GetName();
	nameF += ".gif";
	TString nameFS = c->GetName();
	nameFS += ".small";
	nameFS += ".gif";
	TString nameP = c->GetName();
	nameP += ".eps";
	TString nameC = c->GetName();
	nameC += ".C";
	TString nameR = c->GetName();
	nameR += ".root";
	{
	    TString name = dir;
	    name += "/";
	    name += nameP;
	    c->SaveAs(name.Data());
	}
	{
	    TString name = dir;
	    name += "/";
	    name += nameF;
	    c->SaveAs(name.Data());
	}
	{
	    UInt_t width = c->GetWindowWidth();
	    UInt_t height = c->GetWindowHeight();
	    TString nameS = dir;
	    nameS += "/";
	    nameS += nameFS;
	    UInt_t smallWidth = 140;
	    UInt_t smallHeight = (UInt_t)(smallWidth / 1.62);
	    c->SetCanvasSize(smallWidth, smallHeight);
	    c->Resize();
	    gSystem->ProcessEvents();
	    c->SaveAs(nameS.Data());
	    c->SetCanvasSize(width, height);
	}
	if (saveCR) {
	    TString name = dir;
	    name += "/";
	    name += nameC;
	    c->SaveAs(name.Data());
	}
	if (saveCR) {
	    TString name = dir;
	    name += "/";
	    name += nameR;
	    c->SaveAs(name.Data());
	}
	{
	    TString title = c->GetTitle();
	    ofstr << "<TR>";
	    ofstr << "<TD valign=\"middle\" align=\"right\"><A href=\"" << nameF << "\"><IMG src=\"" << nameFS << "\" border=\"0\"/></A></TD>";
	    ofstr << "<TD valign=\"middle\" align=\"left\">";
	    ofstr << title;
	    ofstr << "<BR/>";
	    ofstr << "<FONT size=\"-1\">";
	    ofstr << "<B>" << index << "</B>: ";
	    ofstr << " <A href=\"" << nameF << "\">gif</A> ";
	    ofstr << " <A href=\"" << nameP << "\">eps</A> ";
	    if (saveCR) {
		ofstr << " <A href=\"" << nameC << "\">C</A> ";
	    }
	    if (saveCR) {
		ofstr << " <A href=\"" << nameR << "\">root</A> ";
	    }
	    ofstr << "</FONT>";
	    ofstr << "</TD>";
	    ofstr << "</TR>" << endl;
	}
	index++;
    }
    index--;
    cout << "Finished, total " << index << " canvases." << endl;
    ofstr << "</TABLE>" << endl;
    ofstr << "</BODY>" << endl;
    ofstr << "</HTML>" << endl;
}
