#include "StKFParticleAnalysisMaker.h"

#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"

#define LambdaPdgMass 1.11568
#define ProtonPdg          2212
#define PiMinusPdg         -211
#define ProtonPdgMass      0.938272
#define PionPdgMass        0.139570

ClassImp(StKFParticleAnalysisMaker);

namespace TpcSubevents{ // Used for TpcEpFinder, to give each subevent (in this case, only one) it's own unique ID, and to know how many subevents there are (with the size element)
  enum{
    MidAbsEtaPsi1,
    size
  };
}

namespace CutVars{ // For MyEventCutter
  enum{
    VzTpc,
    VtransFromAvg,
    NumPrimaryTracks,
    size
  };
}
//____________________________________________________________________________________________________________________________________________________________________//
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fRefmultCorrUtil(0){
  EventsStarted = 0;
  EventsProcessed = 0;

  memset(mBeg,0,mEnd-mBeg+1);
}
//____________________________________________________________________________________________________________________________________________________________________//
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker(){
  SafeDelete(KFParticleInterface);
  SafeDelete(KFParticlePerformanceInterface);
}
Int_t StKFParticleAnalysisMaker::Init(){
  HistogramOutput = new TFile(Form("%s.histograms.root",mFileNameBase.Data()),"recreate");

  h2d_ComEta_nSigmaProton_ = new TH2D("ComEta_nSigmaProton_",";#eta+|y_{beam}|; nSigmaProton",128,-1,1, 128,-10,10);
  h2d_ComEta_nSigmaPion_ = new TH2D("ComEta_nSigmaPion_",";#eta+|y_{beam}|; nSigmaPion",128,-1,1, 128,-10,10);
  h2d_ComEta_NHitsFit_ = new TH2D("ComEta_NHitsFit_",";#eta+|y_{beam}|; nHitsFit",128,-1,1, 128,0,100);
  h2d_ComEta_NHitsDedx_ = new TH2D("ComEta_NHitsDedx_",";#eta+|y_{beam}|; nHitsDedx",128,-1,1, 128,0,100);
  h2d_ComEta_dEdx_ = new TH2D("ComEta_dEdx_",";#eta+|y_{beam}|; dE/dx",128,-1,1, 128,-10,10);
  h2d_ComEta_DedxError_ = new TH2D("ComEta_DedxError_",";#eta+|y_{beam}|; dE/dx Error",128,-1,1, 128,-100,100);
  h2d_ComEta_gDcaToPvMag_ = new TH2D("ComEta_gDcaToPvMag_",";#eta+|y_{beam}|; gDCA to PV (cm)",128,-1,1, 128,0.,10.);
  h2d_FxtMult3_FxtMult_ = new TH2D("FxtMult3_FxtMult_", "; FxtMult; FxtMult3", 300, 0, 300, 140, 0, 140);
   h1d_Dca_Tracks_ = new TH1D("Dca_Tracks_", "Dca", 100, 0, 4);
  
  
  EpdGeom = new StEpdGeom;
  
  vector<vector<float>> CentralityDefinitionUpperLimits = { {5,10,20,30,40,50,60,70,80} };
  MyEventCutter = new EventCutter(CutVars::size,CentralityDefinitionUpperLimits);
    MyEventCutter->InstantiateHists(CutVars::VzTpc,"VzTpc","v_{z, TPC} (cm)",8192,-300,300);
    MyEventCutter->SetCutVarBounds(CutVars::VzTpc,200.,201.4);
    MyEventCutter->InstantiateHists(CutVars::VtransFromAvg,"VtransFromAvg","#sqrt{(v_{x}-#LTv_{x}#GT)^{2}+(v_{y}-#LTv_{y}#GT)^{2}} (cm)",1024,0,15, "Vx","v_{x, TPC} (cm)",1024,-7,7, "Vy","v_{y, TPC} (cm)",1024,-7,7);
    MyEventCutter->SetCutVarBounds(CutVars::VtransFromAvg,0,1.);
    MyEventCutter->InstantiateHists(CutVars::NumPrimaryTracks,"NumPrimaryTracks","Number of Primary Tracks",501,-0.5,500.5);
    MyEventCutter->SetCutVarBounds(CutVars::NumPrimaryTracks,0,500);

// ------------------------------------------------------------- Load the TPC gain-matching information ------------------------------------------------------------- //
// The TPC is gain-matched on a pT, eta, phi, and charge-dependent basis. One wants to minimize the number of bins to enhance statistics, but maximize the number of 
// bins to capture the non-uniformity. The following variable-sized bins capture the features of the non-uniformity while maintaining enough statistics to ensure
// meaningful gain matching
  MyTpcEpFinder = new TpcEpFinder(TpcSubevents::size);
   //MyTpcEpFinder->SetAcceptanceCorrectingFile("/star/u/adams92/LambdaPolarizationAnalyses/3GeV/haddedFiles/E726875751A13F23129B2B9CE19969A3_846995f8962ceb0a194ce1f48787caf0ee8d77e6/E726875751A13F23129B2B9CE19969A3_846995f8962ceb0a194ce1f48787caf0ee8d77e6.histograms.root"); // Simply comment this out if you want to save some time while testing
   //MyTpcEpFinder->SetAcceptanceCorrectingFile("/star/u/jballcap/LambdaFluctuationAnalyses/3GeV/RCF/haddedFiles/E726875751A13F23129B2B9CE19969A3_846995f8962ceb0a194ce1f48787caf0ee8d77e6/E726875751A13F23129B2B9CE19969A3_846995f8962ceb0a194ce1f48787caf0ee8d77e6.histograms.root"); // Simply comment this out if you want to save some time while testing
    vector<float> TpcEfficiencyHistogramPtBinEdges = {
      0.15, 0.153773, 0.157681, 0.161729, 0.165922, 0.170266, 0.174765, 0.179426, 0.184254, 0.189256, 0.194437, 0.199803, 0.205362, 0.21112, 0.217085, 0.223264, 0.229665, 0.236295, 
      0.243163, 0.250277, 0.257647, 0.26528, 0.273188, 0.281379, 0.289864, 0.298654, 0.307758, 0.31719, 0.326959, 0.337079, 0.347562, 0.358421, 0.369669, 0.381321, 0.393391, 0.405894, 
      0.418845, 0.432261, 0.446158, 0.460553, 0.475465, 0.490912, 0.506912, 0.523487, 0.540656, 0.558441, 0.576864, 0.595948, 0.615716, 0.636193, 0.657405, 0.679377, 0.702138, 0.725715, 
      0.750138, 0.775437, 0.801643, 0.828789, 0.856909, 0.886038, 0.916211, 0.947467, 0.979843, 1.01338, 1.04812, 1.08411, 1.12139, 1.16, 1.2};
    vector<float> TpcEfficiencyHistogramEtaBinEdges; for( float EtaBinEdge=-1.5; EtaBinEdge<=0.; EtaBinEdge+=0.02 )TpcEfficiencyHistogramEtaBinEdges.push_back(EtaBinEdge);
    vector<float> TpcEfficiencyHistogramPhiBinEdges; for( float PhiBinEdge=-TMath::Pi(); PhiBinEdge<=TMath::Pi(); PhiBinEdge+=2.*TMath::Pi()/200. )TpcEfficiencyHistogramPhiBinEdges.push_back(PhiBinEdge);
    MyTpcEpFinder->InstantiatePositiveTrackYieldHistograms(TpcEfficiencyHistogramPtBinEdges,TpcEfficiencyHistogramEtaBinEdges,TpcEfficiencyHistogramPhiBinEdges);
    MyTpcEpFinder->InstantiateNegativeTrackYieldHistograms(TpcEfficiencyHistogramPtBinEdges,TpcEfficiencyHistogramEtaBinEdges,TpcEfficiencyHistogramPhiBinEdges);
    MyTpcEpFinder->SetSubeventPtLimits(TpcSubevents::MidAbsEtaPsi1,0.15,1.2);
    MyTpcEpFinder->SetSubeventEtaLimits(TpcSubevents::MidAbsEtaPsi1,-0.7,-0.4);

  EventTree = new IEvent;
  cout << "in BookTree " << endl;
  FemtoDstFile = new TFile(Form("%s.FemtoDst.root",mFileNameBase.Data()),"recreate");
  Tree = new TTree("LambdaTree","Lambda search");
    Tree->Branch("LambdaEvent","IEvent",&EventTree);   
  cout << "Tree booked! \n";

  TFile * f = GetTFile();
  if( f ){f->cd(); BookVertexPlots();}

  return kStOK;
}
//____________________________________________________________________________________________________________________________________________________________________//
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber){
  return StMaker::InitRun(runumber);
}
//____________________________________________________________________________________________________________________________________________________________________//
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt){
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt 
       << "\tMemory : Total = " << info.fMemTotal 
       << "\tUsed = " << info.fMemUsed
       << "\tFree = " << info.fMemFree
       << "\tSwap Total = " << info.fSwapTotal
       << "\tUsed = " << info.fSwapUsed
       << "\tFree = " << info.fSwapFree << endl;
}
//____________________________________________________________________________________________________________________________________________________________________//
void StKFParticleAnalysisMaker::BookVertexPlots(){
  KFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  KFParticlePerformanceInterface = new StKFParticlePerformanceInterface(KFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
}
//____________________________________________________________________________________________________________________________________________________________________//
Int_t StKFParticleAnalysisMaker::Make(){  
  EventsStarted++;
// -------------------------------------------------------------- Load the PicoDst and make event cuts -------------------------------------------------------------- //
  PicoDst = StPicoDst::instance(); if( !PicoDst ) return kStOK;
  StPicoEvent * PicoEvent =  PicoDst->event();	if( !PicoEvent ) return kStOK;

  if( MyEventCutter->RunIsMarkedBad(PicoEvent->runId()) ) {return kStOK;}

  TVector3 PrimaryVertex = PicoEvent->primaryVertex();
  float VtransFromAvg = sqrt(pow(PrimaryVertex.X()+0.08442,2)+pow(PrimaryVertex.Y()+1.926,2)); // The distance in the transverse plane between this primary vertex and the average
  int NumPrimaryTracks = 0;

  int NumPrimaryTracksNoProtons =0; 
     


for( unsigned int iTrack=0; iTrack<PicoDst->numberOfTracks(); iTrack++ ){
	    StPicoTrack *Track = PicoDst->track(iTrack); if( !Track ) continue;
      double TrackComEta = Track->gMom().PseudoRapidity()+1.045;

    
   // TVector3 Origin = Track->origin() ;
   // double dca = Origin.Mag();

      h1d_Dca_Tracks_ -> Fill(Track->gMom().Mag());
      h2d_ComEta_dEdx_->Fill(TrackComEta,Track->dEdx());
      h2d_ComEta_nSigmaProton_->Fill(TrackComEta,Track->nSigmaProton());
      h2d_ComEta_nSigmaPion_->Fill(TrackComEta,Track->nSigmaPion());
      h2d_ComEta_NHitsFit_->Fill(TrackComEta,Track->nHitsFit());
      h2d_ComEta_NHitsDedx_->Fill(TrackComEta,Track->nHitsDedx());
      h2d_ComEta_DedxError_->Fill(TrackComEta,Track->dEdxError());
      h2d_ComEta_gDcaToPvMag_->Fill(TrackComEta,Track->gMom().Mag());
     

 if( Track->isPrimary()) NumPrimaryTracks++;  //All Primary Tracks
      
      if(!Track->isPrimary()) continue;
	float dca = Track -> gDCA(PicoDst->event()->primaryVertex()).Mag();
	int nHitsFit = Track -> nHitsFit(); 
	int nHitsPoss = Track-> nHitsMax();
	float quality= (float)nHitsFit/(float)nHitsPoss; 

     		if(fabs(dca) > 3.0) continue;
	        if(nHitsFit< 15) continue;
		if(quality< 0.51) continue;
      		if( Track->nSigmaProton()>-3.0) continue;
	 NumPrimaryTracksNoProtons++; 
       //Good Primary tracks except Proton 
      
	} // iTrack

  h2d_FxtMult3_FxtMult_ -> Fill(NumPrimaryTracks, NumPrimaryTracksNoProtons);  //Comparing FxtMult3 vs FxtMult

  MyEventCutter->Cut(CutVars::VzTpc,PrimaryVertex.Z());
  MyEventCutter->Cut(CutVars::VtransFromAvg,VtransFromAvg,PrimaryVertex.X(),PrimaryVertex.Y());
  MyEventCutter->Cut(CutVars::NumPrimaryTracks,NumPrimaryTracks);
  MyEventCutter->FinishedCutting();
  if( !MyEventCutter->EventIsGood ){ return kStOK;} // do this AFTER calling EventCutter->FinishedCutting

// -------------------------------------------------------------- Find event-plane angles from the TPC -------------------------------------------------------------- //
  MyTpcEpFinder->FillAcceptanceCorrectionHistograms(PicoDst->picoArray(StPicoArrays::Track));
  vector<vector<int>> TpcSubeventUniqueIdentifiersAndOrders = { {TpcSubevents::MidAbsEtaPsi1, 1} };
  MyTpcEpFinder->FindEventPlaneAngles(PicoDst->picoArray(StPicoArrays::Track),TpcSubeventUniqueIdentifiersAndOrders);

// -------------------------------------------------------------- Find event-plane angles from the EPD -------------------------------------------------------------- //
  int NumHitsEpdSmallAbsEta = 0, NumHitsEpdLargeAbsEta = 0;
  int nEpdHits = PicoDst->numberOfEpdHits();
  TVector2 QEpdSmallAbsEta, QEpdLargeAbsEta;
  for (int iEpdHit = 0; iEpdHit < nEpdHits; iEpdHit++){
    StPicoEpdHit *epdHit = PicoDst->epdHit(iEpdHit);
    if( epdHit->side()>0 ) continue;
    if( epdHit->nMIP()<0.3 ) continue;
    TVector3 randomPointOnTile = EpdGeom->RandomPointOnTile(epdHit->id());
    TVector3 lineToRandomPoint = randomPointOnTile - PrimaryVertex;
    double eta = lineToRandomPoint.Eta();
    double phi = EpdGeom->TileCenter(epdHit->id()).Phi();
    TVector2 QVector(cos(phi),sin(phi));
    if( eta<-3.3 && eta>-3.9 ){
      QEpdLargeAbsEta -= QVector;
      NumHitsEpdLargeAbsEta++;
    }
    else if( eta<-2.6 && eta>-2.9 ){
      QEpdSmallAbsEta -= QVector;
      NumHitsEpdSmallAbsEta++;
    }
  } // iEpdHit

// ------------------------------------------------------------------ Set things up for KFParticle ------------------------------------------------------------------ //
	int maxGBTrackIndex = -1; //find max global track index
	for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++){
		StPicoTrack *Track = PicoDst->track(iTrack); if ( !Track ) continue;
		if ( Track->id()>maxGBTrackIndex ) maxGBTrackIndex = Track->id();
	} // iTrack
	vector<KFMCTrack> mcTracks(0);
	vector<int> triggeredTracks;
	vector<int> mcIndices(maxGBTrackIndex+1);
	for( unsigned int iIndex=0; iIndex<mcIndices.size(); iIndex++) mcIndices[iIndex] = -1;
	if( maxGBTrackIndex>0 ) KFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
	if ( !KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks) ) return kStOK;
	TrackMap.resize(maxGBTrackIndex+1, -1); //make a map from trackID to track index in global track array
	for( unsigned int iTrack=0; iTrack<PicoDst->numberOfTracks(); iTrack++ ){
		StPicoTrack *Track = PicoDst->track(iTrack); if (!Track) continue;
		int index = Track->id();
		TrackMap[index] = iTrack;
	} // iTrack

	KFParticlePerformanceInterface->SetMCTracks(mcTracks);
	KFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
	KFParticlePerformanceInterface->SetCentralityBin(-1);
	KFParticlePerformanceInterface->SetCentralityWeight(1.);
	KFParticlePerformanceInterface->SetPrintEffFrequency(100000);
	KFParticlePerformanceInterface->PerformanceAnalysis();

// ----------------------------------------------------------------- Use KFParticle to find Lambdas ----------------------------------------------------------------- //
// Note: Here, we use KFParticle to get the track IDs of the daughters, and then we reconstruct the Lambdas ourselves. This may seem strange, since KFParticle has the 
// functionality to get the daughter momenta information, but the daughter information is slightly incorrect (for complicated, but known reasons having to do with 
// using the DCA point as a helix fitting point); one can see this, for example, when looking at the sum of the boosted daughter momenta in the Lambda frame which is 
// not zero as it should be. We therefore use KFParticle to identify Lambdas since it is fast and simple, and then we take that information to extract the more 
// accurate daughter information ourselves.
vector<TLorentzVector> Proton4MomentaAtDaughterDca, Pion4MomentaAtDaughterDca, Lambda4Momenta;
vector<TVector3> LambdaDecayPoints;//, gMomProton, pMomProton, gMomPion, pMomPion;
vector<int> ProtonTrackIndices, PionTrackIndices, ProtonCharges;

	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){ 
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
		if( abs(particle.GetPDG())!=3122 ) continue;                //The Absolute value is used to consider particle and antiparticle
    int ProtonTrackIndex = -99999, PionTrackIndex = -99999;

    TLorentzVector Proton4MomentumAtDaughterDca, Pion4MomentumAtDaughterDca;

    for( int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++ ){ 
    	const int daughterId = particle.DaughterIds()[iDaughter]; 
    	const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
    	const int globalTrackId = daughter.DaughterIds()[0];
    	int trackIndex = TrackMap[globalTrackId];    

      daughter.SetProductionVertex(particle);
    	if( abs(daughter.GetPDG())==2212 ){
        ProtonTrackIndex = trackIndex;
        Proton4MomentumAtDaughterDca.SetXYZM(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), ProtonPdgMass);
      }
    	else if(abs(daughter.GetPDG())==211 ){
        PionTrackIndex = trackIndex;
        Pion4MomentumAtDaughterDca.SetXYZM(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), PionPdgMass);
      }
      else{cout<<"----------------------------------- ERROR: KFP found a Lambda but no daughters --------------------------------"<<endl;}
    }  // iDaughter




    if (ProtonTrackIndex == -99999 || PionTrackIndex == -99999) continue; 
    StPicoTrack * ProtonTrack = PicoDst->track(ProtonTrackIndex);
    StPicoTrack * PionTrack = PicoDst->track(PionTrackIndex);
    if( !ProtonTrack ) continue; if( !PionTrack ) continue;

    TLorentzVector Lambda4Momentum; Lambda4Momentum.SetXYZM(particle.GetPx(), particle.GetPy(), particle.GetPz(), particle.GetMass());
    Proton4MomentumAtDaughterDca.Boost(-Lambda4Momentum.BoostVector());
    Pion4MomentumAtDaughterDca.Boost(-Lambda4Momentum.BoostVector());
    TVector3 LambdaDecayPoint(particle.X(),particle.Y(),particle.Z());

		//---------------------------------------- Store info in vectors --------------------------------------//
    Proton4MomentaAtDaughterDca.push_back(Proton4MomentumAtDaughterDca);
    Pion4MomentaAtDaughterDca.push_back(Pion4MomentumAtDaughterDca);
    Lambda4Momenta.push_back(Lambda4Momentum);
    LambdaDecayPoints.push_back(LambdaDecayPoint);
    ProtonTrackIndices.push_back(ProtonTrackIndex);
    PionTrackIndices.push_back(PionTrackIndex);
    ProtonCharges.push_back(ProtonTrack->charge());
	} // End loop over KFParticles

// -------------------------------------------------------------- Filter Lambdas that share daughters --------------------------------------------------------------- //
// KFParticle will sometimes find multiple Lambdas that share the same daughters. Here, we check for such cases; if we find multiple Lambdas sharing the same daughter,
// we choose the one that has the invariant mass closest to the true Lambda mass of 1.11568 GeV. We call this process "Thunderdome", after the Mad Max movies...
//      
 //

 //                                                                                                                             "Two men enter, one man leaves."
                                                                                                                                 
vector<bool> LambdaHasBeenCheckedForMatchingDaughter, LambdaIsThunderdomeWinner;
for( int iProton=0; iProton<(int)ProtonTrackIndices.size(); iProton++ ){
  LambdaHasBeenCheckedForMatchingDaughter.push_back(false);
  LambdaIsThunderdomeWinner.push_back(false);
} // iProton
for( int iProton=0; iProton<(int)ProtonTrackIndices.size(); iProton++ ){
  if( LambdaHasBeenCheckedForMatchingDaughter[iProton] ) continue;
  vector<float> InvMassesOfLambdasCorrespondingToThisProton = {(float)Lambda4Momenta[iProton].M()};
  vector<int> InfoVectorIndicesOfLambdasCorrespondingToThisProton = {iProton};
  for( int jProton=iProton+1; jProton<(int)ProtonTrackIndices.size(); jProton++ ){
    if( ProtonTrackIndices[jProton]==ProtonTrackIndices[iProton] ){
      InvMassesOfLambdasCorrespondingToThisProton.push_back(Lambda4Momenta[jProton].M());
      InfoVectorIndicesOfLambdasCorrespondingToThisProton.push_back(jProton);
      LambdaHasBeenCheckedForMatchingDaughter[jProton] = true;
    }
  } // jProton
  float BestMassDiff = fabs(Lambda4Momenta[iProton].M()-LambdaPdgMass);
  int InfoVectorIndexCorrespondingToBestMassDiff = iProton;
  for( int iProtonMatch=1; iProtonMatch<(int)InfoVectorIndicesOfLambdasCorrespondingToThisProton.size(); iProtonMatch++ ){
    if( fabs(InvMassesOfLambdasCorrespondingToThisProton[iProtonMatch]-LambdaPdgMass)<BestMassDiff ){
      BestMassDiff = fabs(InvMassesOfLambdasCorrespondingToThisProton[iProtonMatch]-LambdaPdgMass);
      InfoVectorIndexCorrespondingToBestMassDiff = InfoVectorIndicesOfLambdasCorrespondingToThisProton[iProtonMatch];
    }
  } // iProtonMatch
  for( int iProtonMatch=0; iProtonMatch<(int)InfoVectorIndicesOfLambdasCorrespondingToThisProton.size(); iProtonMatch++ ){
    if( InfoVectorIndicesOfLambdasCorrespondingToThisProton[iProtonMatch]==InfoVectorIndexCorrespondingToBestMassDiff ) LambdaIsThunderdomeWinner[InfoVectorIndicesOfLambdasCorrespondingToThisProton[iProtonMatch]] = true;
    else LambdaIsThunderdomeWinner[InfoVectorIndicesOfLambdasCorrespondingToThisProton[iProtonMatch]] = false;
  } // iProtonMatch
} // iProton

for( int iPion=0; iPion<(int)PionTrackIndices.size(); iPion++ ) LambdaHasBeenCheckedForMatchingDaughter[iPion] = false;
for( int iPion=0; iPion<(int)PionTrackIndices.size(); iPion++ ){
  if( LambdaHasBeenCheckedForMatchingDaughter[iPion] ) continue;
  if( !LambdaIsThunderdomeWinner[iPion] ) continue;
  vector<float> InvMassesOfLambdasCorrespondingToThisPion = {(float)Lambda4Momenta[iPion].M()};
  vector<int> InfoVectorIndicesOfLambdasCorrespondingToThisPion = {iPion};
  for( int jPion=iPion+1; jPion<(int)PionTrackIndices.size(); jPion++ ){
    if( PionTrackIndices[jPion]==PionTrackIndices[iPion] && LambdaIsThunderdomeWinner[jPion] ){
      InvMassesOfLambdasCorrespondingToThisPion.push_back(Lambda4Momenta[jPion].M());
      InfoVectorIndicesOfLambdasCorrespondingToThisPion.push_back(jPion);
      LambdaHasBeenCheckedForMatchingDaughter[jPion] = true;
    }
  } // jPion
  float BestMassDiff = fabs(Lambda4Momenta[iPion].M()-LambdaPdgMass);
  int InfoVectorIndexCorrespondingToBestMassDiff = iPion;
  for( int iPionMatch=1; iPionMatch<(int)InfoVectorIndicesOfLambdasCorrespondingToThisPion.size(); iPionMatch++ ){
    if( fabs(InvMassesOfLambdasCorrespondingToThisPion[iPionMatch]-LambdaPdgMass)<BestMassDiff ){
      BestMassDiff = fabs(InvMassesOfLambdasCorrespondingToThisPion[iPionMatch]-LambdaPdgMass);
      InfoVectorIndexCorrespondingToBestMassDiff = InfoVectorIndicesOfLambdasCorrespondingToThisPion[iPionMatch];
    }
  } // iPionMatch
  for( int iPionMatch=0; iPionMatch<(int)InfoVectorIndicesOfLambdasCorrespondingToThisPion.size(); iPionMatch++ ){
    if( InfoVectorIndicesOfLambdasCorrespondingToThisPion[iPionMatch]==InfoVectorIndexCorrespondingToBestMassDiff ) LambdaIsThunderdomeWinner[InfoVectorIndicesOfLambdasCorrespondingToThisPion[iPionMatch]] = true;
    else LambdaIsThunderdomeWinner[InfoVectorIndicesOfLambdasCorrespondingToThisPion[iPionMatch]] = false;
  } // iPionMatch
} // iPion

 



// ---------------------------------------------------------------- Fill the tree and wrap things up ---------------------------------------------------------------- //
  for( int iLambda=0; iLambda<(int)Lambda4Momenta.size(); iLambda++ ){
    if( !LambdaIsThunderdomeWinner[iLambda] ) continue;      //condition for Thunderdome
    ITrack ProtonTree, PionTree;
    ILambda * LambdaTree = new ILambda();
    ProtonTree.SetMomentum(Proton4MomentaAtDaughterDca[iLambda]);
    PionTree.SetMomentum(Pion4MomentaAtDaughterDca[iLambda]);
 //   ProtonTree.SetgMom(gMomProton[iLambda]);
//    ProtonTree.SetpMom(pMomProton[iLambda]);
//    PionTree.SetgMom(gMomPion[iLambda]);
//    PionTree.SetpMom(pMomPion[iLambda]);
    LambdaTree->SetThreeMomentum(Lambda4Momenta[iLambda].Vect());
    LambdaTree->SetDecayPoint(LambdaDecayPoints[iLambda]);
    LambdaTree->SetCalcMass(Lambda4Momenta[iLambda].M());
    if( ProtonCharges[iLambda]==1 ){
    	LambdaTree->SetPosDaughter(ProtonTree);
    	LambdaTree->SetNegDaughter(PionTree);
    	EventTree->AddLambda(LambdaTree);
   }//cout<<"Lambda Found"<<endl;  // if Lambda
    else if( ProtonCharges[iLambda]==-1 ){
    	LambdaTree->SetPosDaughter(PionTree);
    	LambdaTree->SetNegDaughter(ProtonTree);
    	EventTree->AddAntiLambda(LambdaTree);
   //cout<<"antilambda Found"<<endl; 
    } // if AntiLambda
  } // iLambda

  EventTree->SetRunNumber(PicoEvent->runId());
  EventTree->SetRefMult(NumPrimaryTracks);
  EventTree->SetRefMult3(NumPrimaryTracksNoProtons); //Setting RefMult3
  EventTree->SetTriggers(PicoEvent->triggerIds());
  EventTree->SetPrimaryVertex(PrimaryVertex);
  EventTree->SetBfield(PicoEvent->bField());
  EventTree->SetPsi1EpdSmallAbsEta(QEpdSmallAbsEta.Phi());
  EventTree->SetPsi1EpdLargeAbsEta(QEpdLargeAbsEta.Phi());
  EventTree->SetPsi1TpcMidAbsEta(MyTpcEpFinder->AccCorrPsi(TpcSubevents::MidAbsEtaPsi1));
  EventTree->SetNumHitsEpdSmallAbsEta(NumHitsEpdSmallAbsEta);
  EventTree->SetNumHitsEpdLargeAbsEta(NumHitsEpdLargeAbsEta);
  EventTree->SetNumHitsTpcMidAbsEta(MyTpcEpFinder->NumTracks(TpcSubevents::MidAbsEtaPsi1));

  Tree->Fill();
  EventTree->ClearEvent();

  EventsProcessed++;
  return kStOK;
}
//____________________________________________________________________________________________________________________________________________________________________//
Int_t StKFParticleAnalysisMaker::Finish(){
  HistogramOutput->cd();
  MyTpcEpFinder->Finish(HistogramOutput);
  MyEventCutter->Finish(HistogramOutput);
  HistogramOutput->Write();
  HistogramOutput->Close();
  FemtoDstFile->cd();
  Tree->Write();
  FemtoDstFile->Close();

  cout<<endl<<endl<<"Events Started: "<<EventsStarted<<"; Events Fully Processed: "<<EventsProcessed<<endl<<endl<<endl;
  
  return kStOK;
}
