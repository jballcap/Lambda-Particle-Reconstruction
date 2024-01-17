#include "TpcEpFinder.h"

//private:

//public:
TpcEpFinder::TpcEpFinder(unsigned int NumSubevents){
    TH1::AddDirectory(false);	std::cout<<"TpcEpFinder is setting TH1::AddDirectory(false)"<<std::endl;

    _my_MiscFunctions = new MiscFunctions();
    _my_MiscFunctions->SetClassNameForErrorMessage("TpcEpFinder");
    _my_MiscFunctions->DontPrintSameErrorMessageMoreThanOnce();
    
    _my_ePositive = 0;
    _my_eNegative = 1;
    _my_ChargeName[_my_ePositive] = "Positive";
    _my_ChargeName[_my_eNegative] = "Negative";

    _my_NumFourierTerms = 12;

    _my_SubeventPtMinima.resize(NumSubevents);
    _my_SubeventPtMaxima.resize(NumSubevents);
    _my_SubeventEtaMinima.resize(NumSubevents);
    _my_SubeventEtaMaxima.resize(NumSubevents);

    _my_SubeventNumShiftingVariables.resize(NumSubevents);
    _my_SubeventNumDigitsInMaxVariableValueUniqueId.resize(NumSubevents);
    _my_h1is_ShiftingVariable__BinningForSubevent.resize(NumSubevents);
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent.resize(NumSubevents);
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent.resize(NumSubevents);

    _my_SubeventNumShiftingVariables_ReadIn.resize(NumSubevents);
    _my_SubeventNumDigitsInMaxVariableValueUniqueId_ReadIn.resize(NumSubevents);
    _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn.resize(NumSubevents);
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn.resize(NumSubevents);
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn.resize(NumSubevents); 

    _my_FindEventPlaneAnglesHasBeenCalled = false;
    _my_SubeventOrders.resize(NumSubevents);
    _my_SubeventRawQVectors.resize(NumSubevents);
    _my_SubeventAccCorrQVectors.resize(NumSubevents);
    _my_SubeventRawEventPlaneAngles.resize(NumSubevents);
    _my_SubeventAccCorrEventPlaneAngles.resize(NumSubevents);
    _my_SubeventNumTracks.resize(NumSubevents);

} // TpcEpFinder

TpcEpFinder::~TpcEpFinder(){
    delete _my_MiscFunctions;

    delete TpcEpFinderDirectory;

    for( int eCharge=0; eCharge<2; eCharge++ ){
        delete _my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge];
        delete _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eCharge];
        delete _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge];
        delete _my_h1is_Pt__WeightHistogramBinEdges[eCharge];
        delete _my_h1is_Eta__WeightHistogramBinEdges[eCharge];
        delete _my_h1is_Phi__WeightHistogramBinEdges[eCharge];
    } // eCharge

    for( unsigned int iSubevent=0; iSubevent<_my_h1is_ShiftingVariable__BinningForSubevent.size(); iSubevent++ ){
        for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_h1is_ShiftingVariable__BinningForSubevent[iSubevent].size(); iShiftingVariable++ ){
            delete _my_h1is_ShiftingVariable__BinningForSubevent[iSubevent][iShiftingVariable];
        } // iShiftingVariable
    } // iSubevent
    for( unsigned int iSubevent=0; iSubevent<_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent.size(); iSubevent++ ){
        delete _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent[iSubevent];
    } // iSubevent
    for( unsigned int iSubevent=0; iSubevent<_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent.size(); iSubevent++ ){
        delete _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent[iSubevent];
    } // iSubevent

    for( unsigned int iSubevent=0; iSubevent<_my_h1is_ShiftingVariable__BinningForSubevent_ReadIn.size(); iSubevent++ ){
        for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[iSubevent].size(); iShiftingVariable++ ){
            delete _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[iSubevent][iShiftingVariable];
        } // iShiftingVariable
    } // iSubevent
    for( unsigned int iSubevent=0; iSubevent<_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn.size(); iSubevent++ ){
        delete _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn[iSubevent];
    } // iSubevent
    for( unsigned int iSubevent=0; iSubevent<_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn.size(); iSubevent++ ){
        delete _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn[iSubevent];
    } // iSubevent

} // ~TpcEpFinder

void TpcEpFinder::InstantiatePositiveTrackYieldHistograms(std::vector<float> PtBinEdges, std::vector<float> EtaBinEdges, std::vector<float> PhiBinEdges){
    _my_ChargedTrackPtBinEdges[_my_ePositive] = PtBinEdges;
    _my_ChargedTrackEtaBinEdges[_my_ePositive] = EtaBinEdges;
    _my_ChargedTrackPhiBinEdges[_my_ePositive] = PhiBinEdges;
    _my_h3is_Pt_Eta_Phi_TrackYield_Charge[_my_ePositive] = new TH3I("Pt_Eta_Phi_Yield_PositiveTracks","Positive track yield in the TPC;p_{T} (GeV/c);#eta;#phi",
                                                           PtBinEdges.size()-1,&PtBinEdges[0], EtaBinEdges.size()-1,&EtaBinEdges[0], PhiBinEdges.size()-1,&PhiBinEdges[0]);
} // InstantiatePositiveTrackYieldHistograms

void TpcEpFinder::InstantiateNegativeTrackYieldHistograms(std::vector<float> PtBinEdges, std::vector<float> EtaBinEdges, std::vector<float> PhiBinEdges){
    _my_ChargedTrackPtBinEdges[_my_eNegative] = PtBinEdges;
    _my_ChargedTrackEtaBinEdges[_my_eNegative] = EtaBinEdges;
    _my_ChargedTrackPhiBinEdges[_my_eNegative] = PhiBinEdges;
    _my_h3is_Pt_Eta_Phi_TrackYield_Charge[_my_eNegative] = new TH3I("Pt_Eta_Phi_Yield_NegativeTracks","Negative track yield in the TPC;p_{T} (GeV/c);#eta;#phi",
                                                           PtBinEdges.size()-1,&PtBinEdges[0], EtaBinEdges.size()-1,&EtaBinEdges[0], PhiBinEdges.size()-1,&PhiBinEdges[0]);
} // InstantiateNegativeTrackYieldHistograms

void TpcEpFinder::SetAcceptanceCorrectingFile(TString AcceptanceCorrectingFileName){
// This function will read in previously created track-yield histograms, and fill "weight" histograms such that the weighted yield will be flat across phi as a function of eta and p_T
    TFile * AcceptanceCorrectingFile = new TFile(AcceptanceCorrectingFileName.Data(),"READ");
    if( !AcceptanceCorrectingFile ){ _my_MiscFunctions->PrintError("Can't find acceptance correction file. returning now"); return;}
    TDirectory * TpcEpFinderDirectory_ReadIn = (TDirectory *)AcceptanceCorrectingFile->Get("TpcEpFinder_CorrectionHistograms");
    if( !TpcEpFinderDirectory_ReadIn ){ _my_MiscFunctions->PrintError("Can't find acceptance correction directory. returning now"); return;}
    TH3I * h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[2];
    for( int i=0; i<2; i++ ) h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[i] = (TH3I *)TpcEpFinderDirectory_ReadIn->Get(Form("Pt_Eta_Phi_Yield_%sTracks",_my_ChargeName[i].Data()));
    if( !h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[0] || !h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[1] ){ _my_MiscFunctions->PrintError("Can't find acceptance correction histogram. returning now"); return;}

    // Get the binning information from the track-yield histogram that you just read in
    int NumPtBins_ReadIn[2], NumEtaBins_ReadIn[2], NumPhiBins_ReadIn[2];
    std::vector<float> ChargedTrackPtBinEdges_ReadIn[2], ChargedTrackEtaBinEdges_ReadIn[2], ChargedTrackPhiBinEdges_ReadIn[2]; // vector indices are bins; array indices are charges ([0] is positve, [1] is negative)
    for( int eCharge=0; eCharge<2; eCharge++ ){
        NumPtBins_ReadIn[eCharge] = h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->GetNbinsX();
        NumEtaBins_ReadIn[eCharge] = h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->GetNbinsY();
        NumPhiBins_ReadIn[eCharge] = h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->GetNbinsZ();

        _my_h1is_Pt__WeightHistogramBinEdges[eCharge] =  h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->ProjectionX();
        _my_h1is_Eta__WeightHistogramBinEdges[eCharge] = h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->ProjectionY();
        _my_h1is_Phi__WeightHistogramBinEdges[eCharge] = h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->ProjectionZ();

        for( int iPtBin=1; iPtBin<=NumPtBins_ReadIn[eCharge]+1; iPtBin++ ) ChargedTrackPtBinEdges_ReadIn[eCharge].push_back(_my_h1is_Pt__WeightHistogramBinEdges[eCharge]->GetBinLowEdge(iPtBin));
        for( int iEtaBin=1; iEtaBin<=NumEtaBins_ReadIn[eCharge]+1; iEtaBin++ ) ChargedTrackEtaBinEdges_ReadIn[eCharge].push_back(_my_h1is_Eta__WeightHistogramBinEdges[eCharge]->GetBinLowEdge(iEtaBin));
        for( int iPhiBin=1; iPhiBin<=NumPhiBins_ReadIn[eCharge]+1; iPhiBin++ ) ChargedTrackPhiBinEdges_ReadIn[eCharge].push_back(_my_h1is_Phi__WeightHistogramBinEdges[eCharge]->GetBinLowEdge(iPhiBin));
    } // eCharge

    // Instantiate the Weight and WeightedYield histograms
    for( int eCharge=0; eCharge<=1; eCharge++ ){
        _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eCharge] = new TH3D(Form("pT_eta_phi_YieldWeights_%sTracks",_my_ChargeName[eCharge].Data()),
    	                                                                Form("%s track efficiency weights;p_{T} (GeV);#eta;#phi",_my_ChargeName[eCharge].Data()),
                                                                        NumPtBins_ReadIn[eCharge],&ChargedTrackPtBinEdges_ReadIn[eCharge][0], 
                                                                        NumEtaBins_ReadIn[eCharge],&ChargedTrackEtaBinEdges_ReadIn[eCharge][0], 
                                                                        NumPhiBins_ReadIn[eCharge],&ChargedTrackPhiBinEdges_ReadIn[eCharge][0]);
        _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge] = new TH3F(Form("pT_eta_phi_WeightedYield_%sTracks",_my_ChargeName[eCharge].Data()),
    	                                                                Form("%s track efficiency-weighted yield;p_{T} (GeV);#eta;#phi",_my_ChargeName[eCharge].Data()),
                                                                        NumPtBins_ReadIn[eCharge],&ChargedTrackPtBinEdges_ReadIn[eCharge][0], 
                                                                        NumEtaBins_ReadIn[eCharge],&ChargedTrackEtaBinEdges_ReadIn[eCharge][0], 
                                                                        NumPhiBins_ReadIn[eCharge],&ChargedTrackPhiBinEdges_ReadIn[eCharge][0]);
    } //eCharge

    // Fill the Weight histograms
    for( int eCharge=0; eCharge<=1; eCharge++ ){
        for( int iPtBin=1; iPtBin<=NumPtBins_ReadIn[eCharge]; iPtBin++ ){
            for( int iEtaBin=1; iEtaBin<=NumEtaBins_ReadIn[eCharge]; iEtaBin++ ){
        	    double NumTracksInEtaSlice = 0.;
        	    for( int iPhiBin=1; iPhiBin<=NumPhiBins_ReadIn[eCharge]; iPhiBin++ ){
        	        NumTracksInEtaSlice += h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->GetBinContent(iPtBin,iEtaBin,iPhiBin);
        	    }// iPhiBin
        	    for( int iPhiBin=1; iPhiBin<=NumPhiBins_ReadIn[eCharge]; iPhiBin++ ){
        	        double Weight = (h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->GetBinContent(iPtBin,iEtaBin,iPhiBin)== 0.)
                                    ?0.:(NumTracksInEtaSlice/((double)NumPhiBins_ReadIn[eCharge]))/h3is_Pt_Eta_Phi_TrackYield_Charge_ReadIn[eCharge]->GetBinContent(iPtBin,iEtaBin,iPhiBin);
        	        _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eCharge]->Fill(
                                    _my_h1is_Pt__WeightHistogramBinEdges[eCharge]->GetBinCenter(iPtBin),
                                    _my_h1is_Eta__WeightHistogramBinEdges[eCharge]->GetBinCenter(iEtaBin),
                                    _my_h1is_Phi__WeightHistogramBinEdges[eCharge]->GetBinCenter(iPhiBin), Weight);
        	    }// iPhiBin
            }// iEtaBin
        } // iPtBin
    } //eCharge

} // SetAcceptanceCorrectingFile

void TpcEpFinder::FillAcceptanceCorrectionHistograms(TClonesArray * Tracks){
    for( int iTrack = 0; iTrack<Tracks->GetEntries(); iTrack++ ){
        StPicoTrack * Track = (StPicoTrack *)(* Tracks)[iTrack];
        if( !Track->isPrimary() ) continue;
        float Pt = Track->gPt(), Eta = Track->gMom().Eta(), Phi = Track->gMom().Phi();
        if( Pt<MinTrackPt || Pt>MaxTrackPt || Track->nHitsPoss()<MinTrackNumHitsPoss || Track->nHits()<MinTrackNumHits || (float)Track->nHitsFit()/(float)Track->nHitsPoss()<MinTrackNumHitsRatio ) continue;
        int eTrackCharge = (Track->charge()>0)?_my_ePositive:_my_eNegative;
        if( _my_h3is_Pt_Eta_Phi_TrackYield_Charge[eTrackCharge] ) _my_h3is_Pt_Eta_Phi_TrackYield_Charge[eTrackCharge]->Fill(Pt,Eta,Phi);
        if( _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eTrackCharge] ){
            float Weight = _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eTrackCharge]->GetBinContent(
                _my_h1is_Pt__WeightHistogramBinEdges[eTrackCharge]->FindBin(Pt), _my_h1is_Eta__WeightHistogramBinEdges[eTrackCharge]->FindBin(Eta), _my_h1is_Phi__WeightHistogramBinEdges[eTrackCharge]->FindBin(Phi));
            _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eTrackCharge]->Fill(Pt,Eta,Phi,Weight);
        }
    } // iTrack

} // FillAcceptanceCorrectionHistograms

void TpcEpFinder::SetSubeventPtLimits(int SubeventUniqueId, float PtMin, float PtMax){
    _my_SubeventPtMinima[SubeventUniqueId] = PtMin;
    _my_SubeventPtMaxima[SubeventUniqueId] = PtMax;
} // SetSubeventPtLimits

void TpcEpFinder::SetSubeventEtaLimits(int SubeventUniqueId, float EtaMin, float EtaMax){
    _my_SubeventEtaMinima[SubeventUniqueId] = EtaMin;
    _my_SubeventEtaMaxima[SubeventUniqueId] = EtaMax;
} // SetSubeventEtaLimits

void TpcEpFinder::SetUpShiftingHistograms(TString SubeventName, unsigned int SubeventUniqueId, unsigned int Order, std::vector<TString> ShiftingVariablesNames, std::vector<std::vector<float>> ShiftingVariablesBinLimits){
    // check uniqueness of SubeventName, variable names?
    for( unsigned int iSubevent=0; iSubevent<_my_SubeventNames.size(); iSubevent++ ){
        if( SubeventName==_my_SubeventNames[iSubevent] ){
            _my_MiscFunctions->PrintError("You already passed SetUpShiftingHistograms a subevent with that name... returning.");
            return;
        }
    }
    _my_SubeventNames.push_back(SubeventName);

    for( unsigned int iShiftingVariable=0; iShiftingVariable<ShiftingVariablesNames.size(); iShiftingVariable++ ){
        for( unsigned int jShiftingVariable=0; jShiftingVariable<iShiftingVariable; jShiftingVariable++ ){
            if( ShiftingVariablesNames[iShiftingVariable]==ShiftingVariablesNames[jShiftingVariable] ){
                _my_MiscFunctions->PrintError("In SetUpShiftingHistograms, ShiftingVariableNames must contain unique variable names for a given subevent. returning");
                return;
            }
        }
    }

    if( ShiftingVariablesNames.size()!=ShiftingVariablesBinLimits.size() ){
        _my_MiscFunctions->PrintError("In SetUpShiftingHistograms, the two vectors passed must be of same dimensionality. returning.");
        return;
    }

    // check for spaces in names
    if( SubeventName.Contains(" ") ){
        _my_MiscFunctions->PrintError("Don't use spaces in SubeventName for SetUpShiftingHistograms. returning");
        return;
    }
    for( TString ShiftingVariableName:ShiftingVariablesNames ){
        if( ShiftingVariableName.Contains(" ") ){
            _my_MiscFunctions->PrintError("Don't use spaces in ShiftingVariablesNames for SetUpShiftingHistograms. returning");
            return;
        }
    }

    _my_SubeventNumShiftingVariables[SubeventUniqueId] = ShiftingVariablesNames.size();

    // When filling the shifting histograms, the set of bins that each variable fall into (in a given event) will be given a single ID.
    // Now we find the maximum value of this variable-value unique ID
    unsigned long long int MaxNumDigitsInVariablesUniqueId = 0;
    for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_SubeventNumShiftingVariables[SubeventUniqueId]; iShiftingVariable++ ){
        unsigned int NumShiftingVariableBins = ShiftingVariablesBinLimits[iShiftingVariable].size() + 1; // +1 is for underflow bin
 
        int NumDigitsInNumShiftingVariableBins = ceil(log10(NumShiftingVariableBins));
        _my_SubeventNumDigitsInMaxVariableValueUniqueId[SubeventUniqueId].push_back(NumDigitsInNumShiftingVariableBins);
        MaxNumDigitsInVariablesUniqueId += NumDigitsInNumShiftingVariableBins;
    } // iShiftingVariable

    // Instantiate binning and shifting histograms. We will read all of these in when shifting.
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent[SubeventUniqueId] = new TProfile2D(Form("ShiftingVariablesUniqueId_iFourier_AvgCosPsi%i_%s",Order,SubeventName.Data()),Form("Subevent %s;Shifting variables unique ID;iFourier;#LTcos(%i #Psi_{%i})#GT",SubeventName.Data(),Order,Order), pow(10,MaxNumDigitsInVariablesUniqueId+1),-0.5,pow(10,MaxNumDigitsInVariablesUniqueId+1)-0.5, _my_NumFourierTerms,0.5,_my_NumFourierTerms+0.5 );
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent[SubeventUniqueId] = new TProfile2D(Form("ShiftingVariablesUniqueId_iFourier_AvgSinPsi%i_%s",Order,SubeventName.Data()),Form("Subevent %s;Shifting variables unique ID;iFourier;#LTsin(%i #Psi_{%i})#GT",SubeventName.Data(),Order,Order), pow(10,MaxNumDigitsInVariablesUniqueId+1),-0.5,pow(10,MaxNumDigitsInVariablesUniqueId+1)-0.5, _my_NumFourierTerms,0.5,_my_NumFourierTerms+0.5 );
    for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_SubeventNumShiftingVariables[SubeventUniqueId]; iShiftingVariable++ ){
        _my_h1is_ShiftingVariable__BinningForSubevent[SubeventUniqueId].push_back( new TH1I(Form("%s__BinningForSubevent%s",ShiftingVariablesNames[iShiftingVariable].Data(),SubeventName.Data()), Form("Binning for shifting subevent %s;%s",SubeventName.Data(),ShiftingVariablesNames[iShiftingVariable].Data()), ShiftingVariablesBinLimits[iShiftingVariable].size()-1,&ShiftingVariablesBinLimits[iShiftingVariable][0] ));
    } // iShiftingVariable

} // SetUpShiftingHistograms

void TpcEpFinder::FillShiftingHistograms(unsigned int SubeventUniqueId, float Psi, unsigned int Order, std::vector<float> &ShiftingVariableValues){
    if( ShiftingVariableValues.size()!=_my_SubeventNumShiftingVariables[SubeventUniqueId] ){
        _my_MiscFunctions->PrintError("In FillShiftingHistograms, you need to pass a vector of the same dimension as when calling SetUpShiftingHistograms. returning.");
        return;
    }

    // Get ShiftingVariablesValuesUniqueId
    unsigned long long int ShiftingVariablesValuesUniqueId = 0;
    unsigned int NumDigitsToOffset = 0;
    for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_SubeventNumShiftingVariables[SubeventUniqueId]; iShiftingVariable++ ){
        ShiftingVariablesValuesUniqueId += _my_h1is_ShiftingVariable__BinningForSubevent[SubeventUniqueId][iShiftingVariable]->FindBin(ShiftingVariableValues[iShiftingVariable])*pow(10,NumDigitsToOffset);
        NumDigitsToOffset+=_my_SubeventNumDigitsInMaxVariableValueUniqueId[SubeventUniqueId][iShiftingVariable];
    } // iShiftingVariable

    for( unsigned int iFourier=1; iFourier<_my_NumFourierTerms; iFourier++ ){
        _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent[SubeventUniqueId]->Fill((float)ShiftingVariablesValuesUniqueId,(float)iFourier,cos(Order*iFourier*Psi)); // explicit type conversion necessary for special case of TProfile::Fill(0,whatever) ambiguity
        _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent[SubeventUniqueId]->Fill((float)ShiftingVariablesValuesUniqueId,(float)iFourier,sin(Order*iFourier*Psi));
    } // iFourier

} // FillShiftingHistograms

void TpcEpFinder::LoadShiftingHistograms(TString ShiftingHistogramsFileName, unsigned int SubeventUniqueId, TString SubeventName, unsigned int Order, std::vector<TString> ShiftingVariablesNames, std::vector<std::vector<float>> ShiftingVariablesBinLimits){
    TFile * ShiftingHistogramFile = new TFile(ShiftingHistogramsFileName.Data(),"READ");
    TDirectory * ShiftingHistogramDirectory = (TDirectory *)ShiftingHistogramFile->Get("TpcEpFinder_CorrectionHistograms");

    if( !ShiftingHistogramFile ){
        _my_MiscFunctions->PrintError("LoadShiftingHistograms can't find the file! returning");
        return;
    }
    if( !ShiftingHistogramDirectory ){
        _my_MiscFunctions->PrintError("LoadShiftingHistograms can't find the directory! returning");
        return;
    }

    for( unsigned int iShiftingVariable=0; iShiftingVariable<ShiftingVariablesNames.size(); iShiftingVariable++ ){
        _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId].push_back( (TH1I *)ShiftingHistogramDirectory->Get(Form("%s__BinningForSubevent%s",ShiftingVariablesNames[iShiftingVariable].Data(),SubeventName.Data())));
        if( !_my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId][iShiftingVariable] ){
            _my_MiscFunctions->PrintError("LoadShiftingHistograms can't find the binning histogram! returning");
            return;
        }

        if( _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId][iShiftingVariable]->GetNbinsX()+1!=(int)ShiftingVariablesBinLimits[iShiftingVariable].size() ){
            _my_MiscFunctions->PrintError("In LoadShiftingHistograms, you are not passing the same variable binning as what's in the read-in binning histograms. returning.");
            return;
        }
        for( int iBin=0; iBin<=_my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId][iShiftingVariable]->GetNbinsX(); iBin++){
            if( _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId][iShiftingVariable]->GetBinLowEdge(iBin+1)!=ShiftingVariablesBinLimits[iShiftingVariable][iBin] ){
                _my_MiscFunctions->PrintError("In LoadShiftingHistograms, you are not passing the same variable binning as what's in the read-in binning histograms. returning.");
                return;
            }
        }
    } // iShiftingVariable


    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn[SubeventUniqueId] = (TProfile2D *)ShiftingHistogramDirectory->Get(Form("ShiftingVariablesUniqueId_iFourier_AvgCosPsi%i_%s",Order,SubeventName.Data()));
    _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn[SubeventUniqueId] = (TProfile2D *)ShiftingHistogramDirectory->Get(Form("ShiftingVariablesUniqueId_iFourier_AvgSinPsi%i_%s",Order,SubeventName.Data()));
    if( !_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn[SubeventUniqueId] || !_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn[SubeventUniqueId] ){
        _my_MiscFunctions->PrintError("LoadShiftingHistograms can't find the shifting histogram! returning");
        return;
    }

    _my_NumFourierTerms_ReadIn = _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn[SubeventUniqueId]->GetNbinsY();

    for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_SubeventNumShiftingVariables[SubeventUniqueId]; iShiftingVariable++ ){
        unsigned int NumShiftingVariableBins = ShiftingVariablesBinLimits[iShiftingVariable].size() + 1; // +1 is for underflow bin
        int NumDigitsInNumShiftingVariableBins = ceil(log10(NumShiftingVariableBins));
 
        _my_SubeventNumDigitsInMaxVariableValueUniqueId_ReadIn[SubeventUniqueId].push_back(NumDigitsInNumShiftingVariableBins);
    } // iShiftingVariable

    _my_SubeventNumShiftingVariables_ReadIn[SubeventUniqueId] = ShiftingVariablesBinLimits.size();
} // LoadShiftingHistograms

void TpcEpFinder::FindEventPlaneAngles(TClonesArray * Tracks, std::vector<std::vector<int>> SubeventUniqueIdsAndOrders){
    // Get the subevent unique IDs and their orders first time around
    if( !_my_FindEventPlaneAnglesHasBeenCalled ){ 
        for( std::vector<int> SubeventUniqueIdAndOrder:SubeventUniqueIdsAndOrders ){
            int SubeventUniqueId = SubeventUniqueIdAndOrder[0];
            _my_SubeventUniqueIdsPassed.push_back(SubeventUniqueId);
            _my_SubeventOrders[SubeventUniqueId] = SubeventUniqueIdAndOrder[1];
        } // SubeventUniqueIdAndOrder
        _my_FindEventPlaneAnglesHasBeenCalled = true;
    } // if !_my_FindEventPlaneAnglesHasBeenCalled

    // Reset the Q vectors
    for( int SubeventUniqueId:_my_SubeventUniqueIdsPassed ){
        _my_SubeventRawQVectors[SubeventUniqueId].Set(0.,0.);
        _my_SubeventAccCorrQVectors[SubeventUniqueId].Set(0.,0.);
        _my_SubeventNumTracks[SubeventUniqueId] = 0;
    } // SubeventUniqueId
    
    for( int iTrack = 0; iTrack<Tracks->GetEntries(); iTrack++ ){
        StPicoTrack * Track = (StPicoTrack *)(* Tracks)[iTrack];
        if( !Track->isPrimary() ) continue;
        float Pt = Track->gPt(), Eta = Track->gMom().Eta(), Phi = Track->gMom().Phi();
        if( Pt<MinTrackPt || Pt>MaxTrackPt || Track->nHitsPoss()<MinTrackNumHitsPoss || Track->nHits()<MinTrackNumHits || (float)Track->nHitsFit()/(float)Track->nHitsPoss()<MinTrackNumHitsRatio ) continue;
        int eTrackCharge = (Track->charge()>0)?_my_ePositive:_my_eNegative;

        for( int SubeventUniqueId:_my_SubeventUniqueIdsPassed ){
            if( Pt<_my_SubeventPtMinima[SubeventUniqueId] || Pt>_my_SubeventPtMaxima[SubeventUniqueId] || Eta<_my_SubeventEtaMinima[SubeventUniqueId] || Eta>_my_SubeventEtaMaxima[SubeventUniqueId] ) continue;

            int Order = _my_SubeventOrders[SubeventUniqueId];
            TVector2 TrackRawQVector(Pt*cos(Order*Phi),Pt*sin(Order*Phi));
            _my_SubeventRawQVectors[SubeventUniqueId] += TrackRawQVector;
            _my_SubeventNumTracks[SubeventUniqueId]++;
            if( _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eTrackCharge] ){
                float Weight = _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eTrackCharge]->GetBinContent(_my_h1is_Pt__WeightHistogramBinEdges[eTrackCharge]->FindBin(Pt), 
                                                                _my_h1is_Eta__WeightHistogramBinEdges[eTrackCharge]->FindBin(Eta), _my_h1is_Phi__WeightHistogramBinEdges[eTrackCharge]->FindBin(Phi));
                TVector2 TrackAccCorrQVector(Weight*Pt*cos(Order*Phi),Weight*Pt*sin(Order*Phi));
                _my_SubeventAccCorrQVectors[SubeventUniqueId] += TrackAccCorrQVector;
            } // if there is TrackYieldWeight histogram
        } // SubeventUniqueId
    } // iTrack

    for( int SubeventUniqueId:_my_SubeventUniqueIdsPassed ){
        int Order = _my_SubeventOrders[SubeventUniqueId];
        _my_SubeventRawEventPlaneAngles[SubeventUniqueId] = _my_MiscFunctions->KeepPhiInRange(_my_SubeventRawQVectors[SubeventUniqueId].Phi()/Order,Order);
        _my_SubeventAccCorrEventPlaneAngles[SubeventUniqueId] = _my_MiscFunctions->KeepPhiInRange(_my_SubeventAccCorrQVectors[SubeventUniqueId].Phi()/Order,Order);
    } // SubeventUniqueId

} // FindEventPlaneAngles

float TpcEpFinder::RawPsi(unsigned int SubeventUniqueId){
    if( SubeventUniqueId>_my_SubeventRawEventPlaneAngles.size() ){
        _my_MiscFunctions->PrintError("RawPsi can't find this subevent. returning 0.");
        return 0.;
    }
    return _my_SubeventRawEventPlaneAngles[SubeventUniqueId];
} // RawPsi

float TpcEpFinder::AccCorrPsi(unsigned int SubeventUniqueId){
    if( SubeventUniqueId>_my_SubeventAccCorrEventPlaneAngles.size() ){
        _my_MiscFunctions->PrintError("AccCorrPsi can't find this subevent. returning 0.");
        return 0.;
    }
    return _my_SubeventAccCorrEventPlaneAngles[SubeventUniqueId];
} // AccCorrPsi

float TpcEpFinder::ShiftedPsi(unsigned int SubeventUniqueId, float Psi, unsigned int Order, std::vector<float> &ShiftingVariableValues){
    if( ShiftingVariableValues.size()!=_my_SubeventNumShiftingVariables_ReadIn[SubeventUniqueId] ){
        _my_MiscFunctions->PrintError("In ShiftedPsi, you are not passing a vector of the correct dimensionality. returning 0.");
       return 0.;
    }

    // Get ShiftingVariablesValuesUniqueId
    unsigned long long int ShiftingVariablesValuesUniqueId = 0;
    unsigned int NumDigitsToOffset = 0;
    for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_SubeventNumShiftingVariables_ReadIn[SubeventUniqueId]; iShiftingVariable++ ){
        if( !_my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId][iShiftingVariable] ){
            _my_MiscFunctions->PrintError("ShiftedPsi can't find the binning histogram. returning 0.");
            return 0.;
        }
        ShiftingVariablesValuesUniqueId += _my_h1is_ShiftingVariable__BinningForSubevent_ReadIn[SubeventUniqueId][iShiftingVariable]->FindBin(ShiftingVariableValues[iShiftingVariable])*pow(10,NumDigitsToOffset);
        NumDigitsToOffset += _my_SubeventNumDigitsInMaxVariableValueUniqueId[SubeventUniqueId][iShiftingVariable];
    } // iShiftingVariable

    float PsiFlattened = Psi;
    for( unsigned int iFourier=1; iFourier<_my_NumFourierTerms_ReadIn; iFourier++ ){
        if( !_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn[SubeventUniqueId] || !_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn[SubeventUniqueId] ){
            _my_MiscFunctions->PrintError("ShiftedPsi can't find shifting histogram. returning 0.");
            return 0.;
        }
        float AvgCosiPsi = _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent_ReadIn[SubeventUniqueId]->GetBinContent((float)ShiftingVariablesValuesUniqueId+1,(float)iFourier,cos(Order*iFourier*Psi)); // explicit type conversion necessary for special case of TProfile::Fill(0,whatever) ambiguity
        float AvgSiniPsi = _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent_ReadIn[SubeventUniqueId]->GetBinContent((float)ShiftingVariablesValuesUniqueId+1,(float)iFourier,sin(Order*iFourier*Psi));
        PsiFlattened += (2./(Order*iFourier))*( AvgCosiPsi*sin(Order*iFourier*Psi) - AvgSiniPsi*cos(Order*iFourier*Psi) );
    } // iFourier
    return _my_MiscFunctions->KeepPhiInRange(PsiFlattened,Order);
} // ShiftedPsi

int TpcEpFinder::NumTracks(unsigned int SubeventUniqueId){
    return _my_SubeventNumTracks[SubeventUniqueId];
} // NumTracks

void TpcEpFinder::Finish(TFile * MainOutputFile){
    MainOutputFile->cd();
    TpcEpFinderDirectory = MainOutputFile->mkdir("TpcEpFinder_CorrectionHistograms");
    TpcEpFinderDirectory->cd();

    // Project 3D yield histograms to 3 2D yield histograms and write
    for( int eCharge=0; eCharge<2; eCharge++ ){
        if( _my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge] && _my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge]->GetEntries()>0 ){
            TH2I * h2i_Pt_Phi_TrackYield_Charge = (TH2I *)_my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge]->Project3D("zx");
                h2i_Pt_Phi_TrackYield_Charge->SetName(Form("Pt_Phi_Yield_%sTracks",_my_ChargeName[eCharge].Data()));
                h2i_Pt_Phi_TrackYield_Charge->SetTitle(Form("%s tracks in the TPC;p_{T} (GeV/c);#phi;Yield",_my_ChargeName[eCharge].Data()));
            TH2I * h2i_Pt_Eta_TrackYield_Charge = (TH2I *)_my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge]->Project3D("yx");
                h2i_Pt_Eta_TrackYield_Charge->SetName(Form("Pt_Eta_Yield_%sTracks",_my_ChargeName[eCharge].Data()));
                h2i_Pt_Eta_TrackYield_Charge->SetTitle(Form("%s tracks in the TPC;p_{T} (GeV/c);#eta;Yield",_my_ChargeName[eCharge].Data()));
            TH2I * h2i_Eta_Phi_TrackYield_Charge = (TH2I *)_my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge]->Project3D("zy");
                h2i_Eta_Phi_TrackYield_Charge->SetName(Form("Eta_Phi_Yield_%sTracks",_my_ChargeName[eCharge].Data()));
                h2i_Eta_Phi_TrackYield_Charge->SetTitle(Form("%s tracks in the TPC;#eta;#phi;Yield",_my_ChargeName[eCharge].Data()));

            _my_h3is_Pt_Eta_Phi_TrackYield_Charge[eCharge]->Write();
            h2i_Pt_Phi_TrackYield_Charge->Write();
            h2i_Pt_Eta_TrackYield_Charge->Write();
            h2i_Eta_Phi_TrackYield_Charge->Write();

            delete h2i_Pt_Phi_TrackYield_Charge;
            delete h2i_Pt_Eta_TrackYield_Charge;
            delete h2i_Eta_Phi_TrackYield_Charge;
        } // if the yield histogram exists and has entries
    } // eCharge

    for( int eCharge=0; eCharge<2; eCharge++ ){
        if( _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge] && _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge]->GetEntries()>0 ){
            TH2F * h2i_Pt_Phi_TrackWeightedYield_Charge = (TH2F *)_my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge]->Project3D("zx");
                h2i_Pt_Phi_TrackWeightedYield_Charge->SetName(Form("Pt_Phi_WeightedYield_%sTracks",_my_ChargeName[eCharge].Data()));
                h2i_Pt_Phi_TrackWeightedYield_Charge->SetTitle(Form("%s tracks in the TPC;p_{T} (GeV/c);#phi;Weighted Yield",_my_ChargeName[eCharge].Data()));
            TH2F * h2i_Pt_Eta_TrackWeightedYield_Charge = (TH2F *)_my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge]->Project3D("yx");
                h2i_Pt_Eta_TrackWeightedYield_Charge->SetName(Form("Pt_Eta_WeightedYield_%sTracks",_my_ChargeName[eCharge].Data()));
                h2i_Pt_Eta_TrackWeightedYield_Charge->SetTitle(Form("%s tracks in the TPC;p_{T} (GeV/c);#eta;Weighted Yield",_my_ChargeName[eCharge].Data()));
            TH2F * h2i_Eta_Phi_TrackWeightedYield_Charge = (TH2F *)_my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge]->Project3D("zy");
                h2i_Eta_Phi_TrackWeightedYield_Charge->SetName(Form("Eta_Phi_WeightedYield_%sTracks",_my_ChargeName[eCharge].Data()));
                h2i_Eta_Phi_TrackWeightedYield_Charge->SetTitle(Form("%s tracks in the TPC;#eta;#phi;Weighted Yield",_my_ChargeName[eCharge].Data()));

            _my_h3is_Pt_Eta_Phi_TrackWeightedYield_Charge[eCharge]->Write();
            h2i_Pt_Phi_TrackWeightedYield_Charge->Write();
            h2i_Pt_Eta_TrackWeightedYield_Charge->Write();
            h2i_Eta_Phi_TrackWeightedYield_Charge->Write();

            delete h2i_Pt_Phi_TrackWeightedYield_Charge;
            delete h2i_Pt_Eta_TrackWeightedYield_Charge;
            delete h2i_Eta_Phi_TrackWeightedYield_Charge;
        } // if the yield histogram exists and has entries
    } // eCharge

    for( int eCharge=0; eCharge<2; eCharge++ ){
        if( _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eCharge] ) _my_h3is_Pt_Eta_Phi_TrackYieldWeight_Charge[eCharge]->Write();
    } // eCharge

    for( unsigned int iSubevent=0; iSubevent<_my_h1is_ShiftingVariable__BinningForSubevent.size(); iSubevent++ ){
        for( unsigned int iShiftingVariable=0; iShiftingVariable<_my_h1is_ShiftingVariable__BinningForSubevent[iSubevent].size(); iShiftingVariable++ ){
            if( _my_h1is_ShiftingVariable__BinningForSubevent[iSubevent][iShiftingVariable] ){ // This is never filled; don't check for entries
                _my_h1is_ShiftingVariable__BinningForSubevent[iSubevent][iShiftingVariable]->Write();
            }
        }
    }
    for( unsigned int iSubevent=0; iSubevent<_my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent.size(); iSubevent++ ){
        if( _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent[iSubevent] && _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent[iSubevent]->GetEntries()>0 ){
            _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgCosPsi_Subevent[iSubevent]->Write();
            _my_tp1s_ShiftingVariablesUniqueId_iFourier_AvgSinPsi_Subevent[iSubevent]->Write();
        }
    }

    MainOutputFile->cd();
} // Finish