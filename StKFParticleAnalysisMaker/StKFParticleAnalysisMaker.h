#ifndef STAR_StKFParticleAnalysisMaker
#define STAR_StKFParticleAnalysisMaker

#ifndef StMaker_H
#include "StMaker.h"
#endif

#include "TMVA/Reader.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"

#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "StEpdGeom.h"
#include "StPicoEvent/StPicoEpdHit.h"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TH1.h"

#include "TpcEpFinder.h"
#include "EventCutter.h"

#include "IEvent.h"
#include "ITrack.h"
#include "ILambda.h"
#include "TTree.h"


class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class KFParticle;
class StPicoDst;
class StMuDst;

class StKFParticleAnalysisMaker : public StMaker {
 private:
  Char_t mBeg[1];
  Char_t mEnd[1];
  StPicoDst * PicoDst;
  StKFParticleInterface * KFParticleInterface;
  StKFParticlePerformanceInterface * KFParticlePerformanceInterface;
  StRefMultCorr * fRefmultCorrUtil;
  bool fIsPicoAnalysis;
  vector<int> TrackMap;

  StEpdGeom * EpdGeom;
  TTree * Tree;
  IEvent * EventTree;
  TpcEpFinder * MyTpcEpFinder;
  EventCutter * MyEventCutter;
  TFile * FemtoDstFile, * HistogramOutput;
  TString mFileNameBase;

  TH1D * h1d_Dca_Tracks_;
  TH2D * h2d_ComEta_nSigmaProton_;
  TH2D * h2d_ComEta_nSigmaPion_;
  TH2D * h2d_ComEta_NHitsFit_;
  TH2D * h2d_ComEta_NHitsDedx_;
  TH2D * h2d_ComEta_dEdx_;
  TH2D * h2d_ComEta_DedxError_;
  TH2D * h2d_ComEta_gDcaToPvMag_;
  TH2D * h2d_FxtMult3_FxtMult_; 

  ULong_t EventsStarted; // Number of Events read
  ULong_t EventsProcessed; // Number of Events processed and analyze
  
 public: 
  StKFParticleAnalysisMaker(const char *name="KFParticleAnalysis");
  void SetOutputFileNameBase(TString OutputFileName="OutputFileName"){mFileNameBase = OutputFileName;}
  virtual       ~StKFParticleAnalysisMaker();
  virtual Int_t  Init();
  virtual Int_t  InitRun(Int_t runumber);
  void           BookVertexPlots();
  virtual Int_t  Make();
  virtual Int_t  Finish();
  Bool_t         Check();
  void AnalysePicoDst() { fIsPicoAnalysis = true; }
  static void    PrintMem(const Char_t *opt = "");
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StKFParticleAnalysisMaker.h,v 1.0 2017/10/07 11:43:53 mzyzak Exp $ built " __DATE__ " " __TIME__ ; 
    return cvs;
  }
  
  ClassDef(StKFParticleAnalysisMaker,0)
};
#endif
