#ifndef MISCFUNCTIONS_H
#define MISCFUNCTIONS_H

#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TColor.h"
#include "TCanvas.h"
#include <vector>

class MiscFunctions{
private:
    TString _my_ClassNameForErrorMessage = "";
    bool _my_PrintSameErrorMessageMoreThanOnce = true;
    std::vector<TString> _my_ErrorMessageLog;

    //TFile * _my_SetStyleHistogramFile;

public:
    MiscFunctions();
    virtual ~MiscFunctions();

    double KeepPhiInRange(double phi, int order);
    void SetClassNameForErrorMessage(TString ClassNameForErrorMessage);
    void DontPrintSameErrorMessageMoreThanOnce();
    void PrintError(TString errorMessage);
    float CalculateResolution(float AvgCos);
    void FillResolutionGraphFromAvgCosProfile(TGraphAsymmErrors * tgae_Resolution, TProfile * tp1_AvgCos);
    void SetStyles(TDirectory * InputDirectory /* This will take a TFile * also */);
    void SetHistStyle(TObject * Histogram);
    void DelimitString(TString InputString, TString Delimiter, std::vector<TString> &DelimitedStrings);
    TString DelimitString(TString InputString, TString Delimiter, unsigned int Column);
};
#endif
