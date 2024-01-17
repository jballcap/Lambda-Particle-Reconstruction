#include "ITrack.h"
#include "TVector3.h"
#include "TLorentzVector.h"

ClassImp(ITrack)

//___________________
ITrack::ITrack(){
  ClearTrack();
}
//_________________________
ITrack::~ITrack(){
  ClearTrack();
}
//_________________________
void ITrack::ClearTrack(){
  mTrackID = -99;
  mMomentum.SetXYZM(-99,-99,-99,-99);
  mgMom.SetXYZ(-99,-99,-99);
  mpMom.SetXYZ(-99,-99,-99);
}
//_________________________
