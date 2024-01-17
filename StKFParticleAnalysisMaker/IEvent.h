// author: Mike Lisa 14 feb 2013 // edited for just lambdas Dec 5 by Isaac Upsal

#ifndef IEVENT_H
#define IEVENT_H

#include "TVector3.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <vector>

class ITrack;
class ILambda;
class IProton;

/*
  Store simple event summary information
*/

class IEvent : public TObject {

 protected:

  Int_t mRunNumber;
  Int_t mChronologicalRunId;
  Int_t mEventID;  // within run
  Int_t mRefMult;
  Int_t mRefMult3; //Excluding protons
  Int_t mCentralityID9;		//centrality ID number from StRefMultCorr RefMult9
  Int_t mCentralityID16;	//centrality ID number from StRefMultCorr RefMult16
  TVector3 mPrimaryVertex;
  Float_t mBfield;
  Float_t mVpdVz;
  Int_t mRefMultPos; // you may want to keep these for the Ach analysis
  Int_t mRefMultNeg;

  std::vector<unsigned int> mTriggers;

  Float_t mPsi1RawEast;  
  Float_t mPsi1RawWest;
  Float_t mPsi1GainMatchedAndv1weightedEast;  
  Float_t mPsi1GainMatchedAndv1weightedWest;

  Float_t mPsi1EpdSmallAbsEta;
  Float_t mPsi1EpdLargeAbsEta;
  Float_t mPsi1TpcMidAbsEta;

  Int_t mNumHitsEpdSmallAbsEta;
  Int_t mNumHitsEpdLargeAbsEta;
  Int_t mNumHitsTpcMidAbsEta;

  TClonesArray* fLambdaCollection;
  TClonesArray* fAntiLambdaCollection;

  void AddLambdaVertex(TClonesArray*, ILambda*);
  ILambda* GetLambdaVertex(TClonesArray*, Int_t index);

 public:

  IEvent();
  ~IEvent();

  void ClearEvent();
  void Init();

  Int_t   NLambda(){return fLambdaCollection->GetEntriesFast();}
  ILambda*    GetLambda(Int_t index);
  void    AddLambda(ILambda* lam);

  Int_t   NAntiLambda(){return fAntiLambdaCollection->GetEntriesFast();}
  ILambda*    GetAntiLambda(Int_t index);
  void    AddAntiLambda(ILambda* lam);


  TVector3 PrimaryVertex() {return mPrimaryVertex;}
  void SetPrimaryVertex(TVector3 pv){mPrimaryVertex = pv;}

  Int_t GetRunNumber(){return mRunNumber;}
  Int_t GetChronologicalRunId(){return mChronologicalRunId;}
  Int_t GetEventID(){return mEventID;}
  Int_t GetRefMult(){return mRefMult;}
  Int_t GetRefMult3(){return mRefMult3;}  //Get RefMult3
  Int_t GetCentralityID9(){return mCentralityID9;}
  Int_t GetCentralityID16(){return mCentralityID16;}
  Float_t Bfield(){return mBfield;}
  Float_t GetVpdVz(){return mVpdVz;}
  Int_t GetRefMultPos(){return mRefMultPos;}
  Int_t GetRefMultNeg(){return mRefMultNeg;}
  std::vector<unsigned int> GetTriggers(){return mTriggers;}

  Float_t GetPsi1EpdSmallAbsEta(){return mPsi1EpdSmallAbsEta;}
  Float_t GetPsi1EpdLargeAbsEta(){return mPsi1EpdLargeAbsEta;}
  Float_t GetPsi1TpcMidAbsEta(){return mPsi1TpcMidAbsEta;}

  Float_t GetNumHitsEpdSmallAbsEta(){return mNumHitsEpdSmallAbsEta;}
  Float_t GetNumHitsEpdLargeAbsEta(){return mNumHitsEpdLargeAbsEta;}
  Float_t GetNumHitsTpcMidAbsEta(){return mNumHitsTpcMidAbsEta;}

  Float_t GetPsi1GainMatchedAndv1weightedEast(){return mPsi1GainMatchedAndv1weightedEast;}
  Float_t GetPsi1GainMatchedAndv1weightedWest(){return mPsi1GainMatchedAndv1weightedWest;}
  Float_t GetPsi1RawEast(){return mPsi1RawEast;}
  Float_t GetPsi1RawWest(){return mPsi1RawWest;}


  void SetRunNumber(Int_t rn){mRunNumber=rn;}
  void SetChronologicalRunId(Int_t rn){mChronologicalRunId=rn;}
  void SetEventID(Int_t en){mEventID=en;}
  void SetRefMult(Int_t rm){mRefMult = rm;}
  void SetRefMult3(Int_t rm){mRefMult3 = rm;} //Set RefMult3
  void SetCentralityID9(Int_t cent){mCentralityID9=cent;}
  void SetCentralityID16(Int_t cent){mCentralityID16=cent;}
  void SetBfield(Float_t f){mBfield=f;}
  void SetVpdVz(Float_t n){mVpdVz=n;}
  void SetRefMultPos(Int_t rm){mRefMultPos = rm;}
  void SetRefMultNeg(Int_t rm){mRefMultNeg = rm;}
  void SetTriggers( const std::vector<unsigned int> & n){mTriggers = n;}

  void SetPsi1RawEast(Float_t p){mPsi1RawEast=p;}
  void SetPsi1RawWest(Float_t p){mPsi1RawWest=p;}
  void SetPsi1GainMatchedAndv1weightedEast(Float_t p){mPsi1GainMatchedAndv1weightedEast=p;}
  void SetPsi1GainMatchedAndv1weightedWest(Float_t p){mPsi1GainMatchedAndv1weightedWest=p;}
  void SetPsi1EpdSmallAbsEta(Float_t p){mPsi1EpdSmallAbsEta=p;}
  void SetPsi1EpdLargeAbsEta(Float_t p){mPsi1EpdLargeAbsEta=p;}
  void SetPsi1TpcMidAbsEta(Float_t p){mPsi1TpcMidAbsEta=p;}

  void SetNumHitsEpdSmallAbsEta(Float_t p){mNumHitsEpdSmallAbsEta=p;}
  void SetNumHitsEpdLargeAbsEta(Float_t p){mNumHitsEpdLargeAbsEta=p;}
  void SetNumHitsTpcMidAbsEta(Float_t p){mNumHitsTpcMidAbsEta=p;}

  ClassDef(IEvent,1)  // my event
};

#endif


