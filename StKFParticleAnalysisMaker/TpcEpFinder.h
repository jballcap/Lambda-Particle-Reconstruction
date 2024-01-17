#ifndef TPCEPFINDER_H
#define TPCEPFINDER_H

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TVector2.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TDirectory.h"

#include "StRoot/StPicoEvent/StPicoTrack.h"

#include "MiscFunctions.h"

#include <vector>
#include <iostream>

class TpcEpFinder{
private:
    MiscFunctions * _my_MiscFunctions;

    TDirectory * TpcEpFinderDirectory; // This needs to be global or else histograms can't be written to it (I don't know why this is)

    unsigned int _my_ePositive, _my_eNegative;  // "e" prefix like "enum"
    TString _my_ChargeName[2];    // used in histogram names and titles

    std::vector<float> _my_ChargedTrackPtBinEdges[2];   // (index: charge ([0] is negative, [1] is positive))
    std::vector<float> _my_ChargedTrackEtaBinEdges[2];
    std::vector<float> _my_ChargedTrackPhiBinEdges[2];

    TH3I * _my_h3is_Pt_Eta_Phi_TrackYield_Charge[2];    // (index: charge ([0] is negative, [1] is positive))
    TH3D * _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[2];
    TH3F * _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[2];
    TH1D * _my_h1is_Pt__WeightHistogramBinEdges[2];
    TH1D * _my_h1is_Eta__WeightHistogramBinEdges[2];
    TH1D * _my_h1is_Phi__WeightHistogramBinEdges[2];

    std::vector<float> _my_SubeventPtMinima, _my_SubeventPtMaxima;    // (index: SubeventUniqueId)
    std::vector<float> _my_SubeventEtaMinima, _my_SubeventEtaMaxima;

    std::vector<unsigned int>               _my_SubeventNumShiftingVariables;                   // (index: SubeventUniqueId)
    std::vector<std::vector<unsigned int>>  _my_SubeventNumDigitsInMaxVariableValueUniqueId;    // (inner index: iShiftingVariable) (outer index: SubeventUniqueId)
    std::vector<std::vector<TH1I *>>        _my_h1is_ShiftingVariable__BinningForSubevent;
    std::vector<TProfile2D *>               _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent;
    std::vector<TProfile2D *>               _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent;
    unsigned int                            _my_NumFourierTerms;      // for shifting

    std::vector<unsigned int>               _my_SubeventNumShiftingVariables_ReadIn;                // (index: SubeventUniqueId)
    std::vector<std::vector<unsigned int>>  _my_SubeventNumDigitsInMaxVariableValueUniqueId_ReadIn; // (inner index: iShiftingVariable) (outer index: SubeventUniqueId)
    std::vector<std::vector<TH1I *>>        _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn;
    std::vector<TProfile2D *>               _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn;
    std::vector<TProfile2D *>               _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn;
    unsigned int                            _my_NumFourierTerms_ReadIn;    // for shifting

    std::vector<TString> _my_SubeventNames;
    std::vector<int>   _my_SubeventOrders;                 // (index: SubeventUniqueId)
    std::vector<TVector2> _my_SubeventRawQVectors;
    std::vector<TVector2> _my_SubeventAccCorrQVectors;
    std::vector<float> _my_SubeventRawEventPlaneAngles;
    std::vector<float> _my_SubeventAccCorrEventPlaneAngles;

    std::vector<int> _my_SubeventNumTracks;
    
    bool _my_FindEventPlaneAnglesHasBeenCalled;
    std::vector<unsigned int> _my_SubeventUniqueIdsPassed; // for FindEventPlaneAngles

public:
    float MinTrackPt = 0.15;    // basic track QA cuts; change these if you want
    float MaxTrackPt = 10.;
    float MinTrackNumHits = 12;   //Standar was 15 modified by Jonathan
    float MinTrackNumHitsPoss = 5;
    float MinTrackNumHitsRatio = 0.52;

    TpcEpFinder(unsigned int NumSubevents);
    virtual ~TpcEpFinder();
    
    void InstantiatePositiveTrackYieldHistograms(std::vector<float> PtBinEdges, std::vector<float> EtaBinEdges, std::vector<float> PhiBinEdges);
    void InstantiateNegativeTrackYieldHistograms(std::vector<float> PtBinEdges, std::vector<float> EtaBinEdges, std::vector<float> PhiBinEdges);
    void SetAcceptanceCorrectingFile(TString AcceptanceCorrectingFileName);
    void FillAcceptanceCorrectionHistograms(TClonesArray * Tracks);

    void SetSubeventPtLimits(int SubeventUniqueId, float PtMin, float PtMax);
    void SetSubeventEtaLimits(int SubeventUniqueId, float EtaMin, float EtaMax);

    void SetUpShiftingHistograms(TString SubeventName, unsigned int SubeventUniqueId, unsigned int Order, std::vector<TString> ShiftingVariablesNames, std::vector<std::vector<float>> ShiftingVariablesBinLimits);
    void FillShiftingHistograms(unsigned int SubeventUniqueId, float Psi, unsigned int Order, std::vector<float> &ShiftingVariableValues);
    void LoadShiftingHistograms(TString ShiftingHistogramsFileName, unsigned int SubeventUniqueId, TString SubeventName, unsigned int Order, std::vector<TString> ShiftingVariablesNames, std::vector<std::vector<float>> ShiftingVariablesBinLimits);

    void FindEventPlaneAngles(TClonesArray * Tracks, std::vector<std::vector<int>> SubeventUniqueIdsAndOrders);
    float RawPsi(unsigned int SubeventUniqueId);
    float AccCorrPsi(unsigned int SubeventUniqueId);
    float ShiftedPsi(unsigned int SubeventUniqueId, float Psi, unsigned int Order, std::vector<float> &ShiftingVariableValues);
    int NumTracks(unsigned int SubeventUniqueId);

    void Finish(TFile * MainOutputFile);
    
}; // class TpcEpFinder

#endif
