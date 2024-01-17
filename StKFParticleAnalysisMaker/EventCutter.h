// Notes are on the bottom of this file
#ifndef EVENTCUTTER_H
#define EVENTCUTTER_H

#include "MiscFunctions.h"
#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>

class EventCutter{
private:
    TDirectory * _my_CutVarDistHistos;                          // this must be global or else we can't save histograms here
    TDirectory * _my_RefMultDistributionsSeparateCentralities;  // again, this has to be global

    MiscFunctions *_my_MiscFunctions;

    std::vector<TH1I *> _my_h1is_CutVar_nEvts_BforeAllCuts;                 // 1D histograms showing a variable's distribution before and after cuts
    std::vector<TH1I *> _my_h1is_CutVar_nEvts_AfterAllCuts;
    std::vector<TH1I *> _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut;
    std::vector<TH2I *> _my_h2is_VarX_VarY_nEvts_BforeAllCuts;              // 2D histograms showing two variables' distributions before and after cuts
    std::vector<TH2I *> _my_h2is_VarX_VarY_nEvts_AfterAllCuts;
    std::vector<TH2I *> _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut;

    std::vector<std::vector<TH1I *>> _my_h1is_RefMult_nEvts_SeparateCentralities;      // 1D histograms showing RefMult distribution for each cent. bin separately

    bool _my_EventIsGood = true;                                // true until proven otherwise
    bool _my_FinishedCuttingWasCalled = false;                  // make sure the user called the function
    unsigned int _my_NumCutVars;                                // how many cut variables we are using
    std::vector<float> _my_CutVarMinima, _my_CutVarMaxima;      // the bounds of a variable for the cut
    std::vector<bool> _my_CutVarValIsGood;                      // the value is within the bounds
    std::vector<float> _my_CutVarVals;                          // the value of each variable in this event
    std::vector<std::vector<float>> _my_VarValsXY;              // the values of 2D plotting variables in this event

    std::vector<bool> _my_CutWasCalled;     // make sure the "Cut" function was called before we start making histograms for this variable
    std::vector<bool> _my_Cut2DWasCalled;   // the overloaded cut function used for plotting 2D histograms

    bool _my_PreviousRunWasMarkedBad = true;                                                            // to speed things up; most events have the same run number as the previous
    unsigned int _my_PreviousRunNumberForFunctionRunIsMarkedBad = 0, _my_PreviousRunNumberForFunctionChronologicalRunId = 0, // each function should track these separately
                 _my_PreviousChronologicalRunIdForFunctionChronologicalRunId = 0;
    std::vector<std::vector<unsigned int>> _my_ChronologicalBadRunLists, _my_ChronologicalGoodRunLists; // lists of good and bad run numbers for separate productions 
    std::vector<unsigned int>              _my_ChronologicalBadRunList,  _my_ChronologicalGoodRunList;  // which of the above lists corresponds to this analysis's production

    std::vector<std::vector<unsigned int>> _my_RefMultMinimaLimits;
    std::vector<std::vector<float>> _my_CentralityUpperLimits;
    unsigned int _my_NumCentralityBinLists;

    void _my_FindRelevantRunLists(unsigned int RunNumber);

public:
    EventCutter(unsigned int NumCutVars, std::vector<std::vector<float>> CentralityUpperLimits);
    virtual ~EventCutter();
    void InstantiateHists(unsigned int CutVarId, TString CutVarName, TString CutVarTitle, unsigned int NumBins,  float HistMin,  float HistMax);
    void InstantiateHists(unsigned int CutVarId, TString CutVarName, TString CutVarTitle, unsigned int NumBins,  float HistMin,  float HistMax,
                                                 TString VarNameX,   TString VarTitleX,   unsigned int NumBinsX, float HistMinX, float HistMaxX,
                                                 TString VarNameY,   TString VarTitleY,   unsigned int NumBinsY, float HistMinY, float HistMaxY);
    void SetCentralityLimits(unsigned int RefMultCutVarId, TString RefMultHistogramFileName, TString RefMultVariableNameInHistogram, int RefMultMax);
    float GetCentrality(unsigned int RefMult, unsigned int NumCentralityBins);
    void SetCutVarBounds(unsigned int CutVarId, float CutVarMin, float CutVarMax);
    bool RunIsMarkedBad(unsigned int RunNumber);
    int ChronologicalRunId(unsigned int RunNumber);
    void Cut(unsigned int CutVarId, float CutVarVal, float VarValX, float VarValY);
    void Cut(unsigned int CutVarId, float CutVarVal);
    void FinishedCutting();
    void Finish(TFile * MainOutputFile);

    bool EventIsGood = false; // if the user accidentally accesses this value before or without calling FinishedCutting, they will always see a bad event
}; // class EventCutter
#endif

// Author: Joey Adams (adams.1940@osu.edu/joseph.r.adams92@gmail.com)
// Last updated 29-01-2020
//
// This class will help when you make event cuts. Here is how to use it:
//  1. In the header file, create a pointer. e.g.:
//      - EventCutter * MyEventCutter;
//  2. In the implementation file, put an enum with all the variables you will, or even might use as event cuts. You'll use these as inputs 
//     for the subsequent function calls. Put a "size" element at the end since you'll later need to know how many elements are in this. e.g.:
//      - namespace CutVars{ // For MyEventCutter
//            enum{
//              Vz,
//              VpdVz,
//              size // you'll need this when instantiating
//            };
//        }
//  3. Before running over events, (in the "Init" method if you're using a maker) instantiate MyEventCutter, specifying what output file you're 
//     putting histograms into and how many cut variables you will, or even potentially will, use. e.g. 
//      - MyEventCutter = new EventCutter(MyHistogramOutput,CutVars::size);
//  4. Instantiate the histograms, giving the name of the variable (that will be used for the histogram name, so no spaces!), the title of the 
//     variable (that will be used for the histogram titles), and the binning (as usual for histograms; # bins, min, max). e.g.
//      -  MyEventCutter->InstantiateHists(CutVars::Vz,"Vz","v_{z, TPC} (cm)",512,-100,100);
//      -  MyEventCutter->InstantiateHists(CutVars::VpdVz,"VpdVz","v_{z, VPD} (cm)",512,-100,100);
//  5. Set the valid (inclusive) bounds for the variable. If you want to see the variable's distribution after all other cuts, but don't want to 
//     actually cut on the variable, just make these bounds very large. e.g.:
//      -  MyEventCutter->SetCutVarBounds(CutVars::Vz,-10,10);
//      -  MyEventCutter->SetCutVarBounds(CutVars::VpdVz,-100,100) //
//  6. When running over events, after doing your trigger cut (which this class does not do), cut on each of the variables (if there are any you 
//     don't want to cut on, just make the bounds large or simply don't call "Cut" on that variable). e.g.:
//      -  float Vz = PicoEvent->primaryVertex().Z();
//      -  float VpdVz = PicoEvent->vzVpd();
//      -  MyEventCutter->Cut(CutVars::Vz,Vz);
//      -  MyEventCutter->Cut(CutVars::VpdVz,VpdVz); // I don't want
//  7. When finished calling "Cut" on each of your variables, call "Finished Cutting". e.g.:
//      - MyEventCutter->FinishedCutting();
//  8. Access "EventIsGood" to determine whether you should analyze the event. e.g.:
//      - if( !MyEventCutter->EventIsGood ) return kStOK; // do this AFTER calling EventCutter->FinishedCutting()
//  9. After making a RefMult distribution with something like 
//      - MyEventCutter->Cut(CutVars::RefMult,picoEvent->refMult());
//     you can use the function SetCentralityLimits in the Init method (or otherwise before running over evetns) to find the centrality bins; no 
//     need to find them yourself! Give it the standard variable ID (which you use when calling Cut), the name of the file containing your RefMult 
//     distribution histogram, the name of the variable in that histogram (something like "RefMult"; you set this when calling InstantiateHistograms 
//     last time), and an upper limit on RefMult (for pileup) (use -1 if you don't want an upper limit). This will look something like:
//      - MyEventCutter->SetCentralityLimits(CutVars::RefMult,"/star/u/adams92/Desktop/Example.root","RefMult",-1);
//     After calling SetCentralityLimits, you no longer need to call SetCutVarBounds for RefMult (although, as long as that is called before calling 
//     SetCentralityLimits, there's no problem). Now you can call Cut on RefMult as usual, and you will by default reject >80% central events and 
//     pilup.
// 10. You can get the centrality in increments with the function GetCentrality (e.g. 0-5% centrality will by 2.5). You only have to pass it RefMult 
//     and how many centrality bins you are using; for example, the common number of centrality bins used is 9, with centralities 0-5%, 5-10%, 
//     10-20%, ... , 70-80%. This will look something like:
//      - MyEventCutter->GetCentrality(picoEvent->refMult(),9); // Will give 1 of 9 centralities
//     This is currently configured to take values of 8, 9, or 16 ({10,20,30,40,50,60,70,80}, {5,10,20,30,40,50,60,70,80}, 
//     {5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80}); if you want anything else, you can add it to the constructor.
// Other notes: 
//  - Values of -999 are considered to be good
//  - If you want to make 2D histograms of two variables that are related to a single cut (e.g. Vx and Vy for the transverse-vertex cut):
//      1. Call "InstlantiateHists", specifying the names and titles of all variables. Just use the cut variable's value from the enum. e.g.:
//          - MyEventCutter->InstantiateHists(CutVars::Vtrans,"Vtrans","v_{T} (cm)",256,0,50, "Vx","v_{x, TPC} (cm)",128,-20,20, "Vy","v_{y, TPC} (cm)",128,-20,20);
//      2. When calling "Cut", just throw in the x and y variables you want plotted in the 2D histogram. e.g.:
//          - TVector3 PrimaryVertex = PicoEvent->primaryVertex();
//          - float Vtrans = sqrt(pow(PrimaryVertex.X(),2)+pow(PrimaryVertex.Y(),2));
//          - MyEventCutter->Cut(CutVars::Vtrans,Vtrans,PrimaryVertex.X(),PrimaryVertex.Y());