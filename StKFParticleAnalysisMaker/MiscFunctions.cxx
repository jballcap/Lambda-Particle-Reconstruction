#include "MiscFunctions.h"
#include <iostream>

using namespace std;

MiscFunctions::MiscFunctions(){
    //_my_MiscFunctions = new MiscFunctions();
    //_my_MiscFunctions->SetClassNameForErrorMessage("MiscFunctions");
}
//__________________________________________________________________________________________________________________________________________________________//
MiscFunctions::~MiscFunctions(){}
//__________________________________________________________________________________________________________________________________________________________//
double MiscFunctions::KeepPhiInRange(double phi, int order){
    double phiInRange = phi;
    while (phiInRange < -TMath::Pi()/(double)order) phiInRange += 2.*TMath::Pi()/(double)order;
    while (phiInRange >= TMath::Pi()/(double)order) phiInRange -= 2.*TMath::Pi()/(double)order;
    return phiInRange;
} //KeepPhiInRange
//__________________________________________________________________________________________________________________________________________________________//
void MiscFunctions::SetClassNameForErrorMessage(TString ClassNameForErrorMessage){
    _my_ClassNameForErrorMessage = ClassNameForErrorMessage;
} // SetClassNameForErrorMessage
//__________________________________________________________________________________________________________________________________________________________//
void MiscFunctions::DontPrintSameErrorMessageMoreThanOnce(){
    _my_PrintSameErrorMessageMoreThanOnce = false;
} // DontPrintSameErrorMessageMoreThanOnce
//__________________________________________________________________________________________________________________________________________________________//
void MiscFunctions::PrintError(TString errorMessage){
    if (!_my_PrintSameErrorMessageMoreThanOnce){
	for (unsigned int i = 0; i < _my_ErrorMessageLog.size(); i++){
	    if (errorMessage.EqualTo(_my_ErrorMessageLog[i])) return;
	}
	_my_ErrorMessageLog.push_back(errorMessage);
    } // if (Don't print the same error message more than once)
    
    std::cout<<"=============================================================== ERROR =============================================================="<<endl;
    if ( !_my_ClassNameForErrorMessage.EqualTo("") ) std::cout<<"Error in "<<_my_ClassNameForErrorMessage.Data()<<":"<<endl;
    std::cout<<errorMessage.Data()<<endl;
    std::cout<<"===================================================================================================================================="<<endl;
} // PrintError
//__________________________________________________________________________________________________________________________________________________________//
float MiscFunctions::CalculateResolution(float AvgCos){
    TF1 Subf( "Subf_som", "(sqrt(TMath::Pi())/2.)*x*exp(-x*x/2.)*( TMath::BesselI(0,x*x/2.) + TMath::BesselI(1,x*x/2.) ) - [0]", 0, 10.0 );
    float SubCorr = AvgCos;
    float SubRes = sqrt(SubCorr);
    Subf.SetParameters(SubRes, 0.0);
    ROOT::Math::WrappedTF1 wf1(Subf);
    ROOT::Math::BrentRootFinder brf;
    brf.SetFunction( wf1, 0, 10.0 );
    brf.Solve();
    float EqnRoot = brf.Root();//this is the numerical solution of the equation
    float Chi = EqnRoot*sqrt(2.); //chi is really sqrt2*Chi but this is what I plug in for x
    float RealResCorr = (sqrt(TMath::Pi())/2.)*Chi*exp(-Chi*Chi/2.)*( TMath::BesselI(0,Chi*Chi/2.) + TMath::BesselI(1,Chi*Chi/2.) );
    return RealResCorr;
} // CalculateResolution
//__________________________________________________________________________________________________________________________________________________________//
void MiscFunctions::FillResolutionGraphFromAvgCosProfile(TGraphAsymmErrors * tgae_Resolution, TProfile * tp1_AvgCos){
int nDepVarBins = tp1_AvgCos->GetNbinsX();
for(int iDepVarBin = 1; iDepVarBin <= nDepVarBins; iDepVarBin++){
    float AvgCosValue = tp1_AvgCos->GetBinContent(iDepVarBin);
    float AvgCosError = tp1_AvgCos->GetBinError(iDepVarBin);
    float ResolutionValue = CalculateResolution(AvgCosValue);
    tgae_Resolution->SetPoint(iDepVarBin,
       tp1_AvgCos->GetBinCenter(iDepVarBin),
       ResolutionValue);
    tgae_Resolution->SetPointError(iDepVarBin,
       tp1_AvgCos->GetBinWidth(iDepVarBin)/2.,
       tp1_AvgCos->GetBinWidth(iDepVarBin)/2.,
       fabs(CalculateResolution(AvgCosValue - AvgCosError)-ResolutionValue),
       fabs(CalculateResolution(AvgCosValue + AvgCosError)-ResolutionValue));
} // iDepVarBin
} // FillResolutionGraphFromAvgCosProfile

void MiscFunctions::SetStyles(TDirectory * InputDirectory /* This will take a TFile * also */){ // "taken" from tutorials/io/readCode.C
// This function and SetHistStyles stylize histograms, graphs, and canvases to my liking and re-write them (hence why 
// you will see keys appended with ";2" after calling this function.) For reasons unknown to me, and that I no longer 
// care to investigate, some styles will not propagate to the file and must be set with my rootlogon.C file (kept in .root/); 
// likewise, some settings in my rootlogon.C file will not propagate to the plotted histograms. Because of this, the two 
// must work in conjunction. My rootlogon.C is:
//{
//gStyle->SetOptStat(0);
//gStyle->SetPalette(kBlueGreenYellow);
//gStyle->SetTitleFont(82,"T");
//gStyle->SetPadTickX(1); // workes on histograms, but not on canvases!
//gStyle->SetPadTickY(1);
//gStyle->SetPadTopMargin(0.13); // workes on histograms, but not on canvases!
//gStyle->SetPadBottomMargin(0.25);
//gStyle->SetPadLeftMargin(0.17);
//gStyle->SetPadRightMargin(0.25);
//gStyle->SetCanvasDefW(1500);
//gStyle->SetCanvasDefH(900);
//gStyle->SetTitleSize(0.06,"T"); // for title
//}

    TDirectory * InitialDirectory = gDirectory;
    TIter next(InputDirectory->GetListOfKeys());
    TKey * key;
    while( (key = (TKey *)next()) ){
        if( key->IsFolder() ){
           InputDirectory->cd(key->GetName());
           SetStyles(gDirectory);
           InitialDirectory->cd();
        }
        else{
            TObject * Object = key->ReadObj() ;
            if( Object->InheritsFrom("TH1") && !Object->InheritsFrom("TH2") && !Object->InheritsFrom("TH3") ){
                SetHistStyle(Object);
                static_cast<TH1 *>(Object)->SetLineColor(kBlue+1); // This only sets the color as long as (Object) is in scope
                static_cast<TH1 *>(Object)->Write("",TObject::kWriteDelete);
            } // TH1
            if( Object->InheritsFrom("TH2") && !Object->InheritsFrom("TH3") ){
                SetHistStyle(Object);
                static_cast<TH2 *>(Object)->Write("",TObject::kWriteDelete);
            } // TH2
            if( Object->InheritsFrom("TH3") ){
                SetHistStyle(Object);
                static_cast<TH3 *>(Object)->Write("",TObject::kWriteDelete);
            } // TH3
            if( Object->InheritsFrom("TCanvas") ){
                for(const auto&& ObjectInCanvas: *(static_cast<TCanvas *>(Object)->GetListOfPrimitives())){
                    if( ObjectInCanvas->InheritsFrom("TH1") ) SetHistStyle(ObjectInCanvas);
                }
                static_cast<TCanvas *>(Object)->SetTopMargin(0.13);
                static_cast<TCanvas *>(Object)->SetBottomMargin(0.25);
                static_cast<TCanvas *>(Object)->SetLeftMargin(0.17);
                static_cast<TCanvas *>(Object)->SetRightMargin(0.25);
                static_cast<TCanvas *>(Object)->SetGrid(0,0);
                static_cast<TCanvas *>(Object)->SetTickx(1);
                static_cast<TCanvas *>(Object)->SetTicky(1);
                static_cast<TCanvas *>(Object)->Update();
                static_cast<TCanvas *>(Object)->Write("",TObject::kWriteDelete);
            } // TCanvas
        } // anything that's not a directory
    } // loop over keys
} // SetStyles

void MiscFunctions::SetHistStyle(TObject * Histogram){
    static_cast<TH1 *>(Histogram)->GetXaxis()->SetTitleSize(0.06); // TProfiles, TH2s, and TH3s inherit from TH1 :)
    static_cast<TH1 *>(Histogram)->GetYaxis()->SetTitleSize(0.06);
    static_cast<TH1 *>(Histogram)->GetZaxis()->SetTitleSize(0.06);
    static_cast<TH1 *>(Histogram)->GetXaxis()->SetTitleFont(82);
    static_cast<TH1 *>(Histogram)->GetYaxis()->SetTitleFont(82);
    static_cast<TH1 *>(Histogram)->GetZaxis()->SetTitleFont(82);
    static_cast<TH1 *>(Histogram)->GetXaxis()->SetLabelFont(82);
    static_cast<TH1 *>(Histogram)->GetYaxis()->SetLabelFont(82);
    static_cast<TH1 *>(Histogram)->SetLineWidth(2);
}

void MiscFunctions::DelimitString(TString InputString, TString Delimiter, vector<TString> &DelimitedStrings){
    int InputStringLength = InputString.Length(), DelimiterLength = Delimiter.Length();
    int DelimitedStringStart = 0;     // where, in InputString, the deliminated string starts (i.e. the i'th character)
    int DelimitedStringLength = 0;
    bool DelimiterFoundInInputString = false;
    for( int iInputStringChar=0; iInputStringChar<=(InputStringLength-DelimiterLength); ){ // Loop over characters in InputString to look for Delimiters
        if( static_cast<TString>(InputString(iInputStringChar,DelimiterLength)).EqualTo(Delimiter) ){
            DelimitedStrings.push_back(InputString(DelimitedStringStart,DelimitedStringLength));
            DelimiterFoundInInputString = true;                           // So you can store the substring following the /final/ Delimiter
            DelimitedStringStart = iInputStringChar+DelimiterLength;    // The next deliminated string will start after this Delimiter
            DelimitedStringLength = 0;                                    // its length starts at zero, incrementing until we find another Delimiter
            iInputStringChar += DelimiterLength;                          // keep looking, starting with InputString's characters after this Delimiter
        }
        else{ // the "iInputStringChar"'th character is not the beginning of a Delimiter
            iInputStringChar++;         // check the next character
            DelimitedStringLength++;  // Keep adding this until you find a Delimiter
        }
    }  // iInputStringChar ---  loop over characters in InputString to look for Delimiters
    // Store the substring following the /final/ Delimiter:
        if( DelimiterFoundInInputString ) DelimitedStrings.push_back(InputString(DelimitedStringStart,InputStringLength-DelimitedStringStart));
} // DelimitString

TString MiscFunctions::DelimitString(TString InputString, TString Delimiter, unsigned int Column){
    vector<TString> DelimitedStrings;
    DelimitString(InputString,Delimiter,DelimitedStrings);
    if( DelimitedStrings.size()<=Column ) return "";
    return DelimitedStrings[Column];
} // DelimitString