void analysis(Int_t NumEvents = 10000, TString InputFileName , TString OutputDirectory, TString OutputName)
{
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  gROOT->LoadMacro("lMuDst.C");
  
  lMuDst(-1,InputFileName.Data(),"ry2016,RpicoDst,mysql,kfpAna,quiet,nodefault","pico.root");
    
  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  kfpAnalysis->SetOutputFileNameBase(OutputName);
  
  chain->Init();

//StKFParticleInterface::instance()->CleanLowPVTrackEvents(); //removes out of consideration events, where more than 90% of tracks are considered as secondary.
  

//StKFParticleInterface::instance()->SetSoftTofPidMode(); //Tof is used in case of Tof hit exist
  
  
//  StKFParticleInterface::instance()->SetChiPrimaryCut(10); //selection primary and secondary tracks
 
/* StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1); //cm
  StKFParticleInterface::instance()->SetLCut(3.f); // distance PV - LambdaVerex (cm)
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(10); // Secondary tracks for Lambda > 8


  StKFParticleInterface::instance()->SetChi2Cut2D(1); // Chi^2 < 3 : two tracks p and pi come from the same space point
  StKFParticleInterface::instance()->SetLdLCut2D(10);   // L/dl > 5
*/


/******Trying different cuts to try to get antilambdas******/
//Cuts based on Jakub Presentation 
  //StKFParticleInterface::instance()->SetChiPrimaryCut(12); //selection primary and secondary tracks  
  //StKFParticleInterface::instance()->SetLCut(5.f); // distance PV - LambdaVerex (cm)
  //StKFParticleInterface::instance()->SetChiPrimaryCut2D(8); // Secondary tracks for Lambda > 8
  //StKFParticleInterface::instance()->SetChi2Cut2D(3); // Chi^2 < 3 : two tracks p and pi come from the same space point
  //StKFParticleInterface::instance()->SetLdLCut2D(8);   // L/dl > 5


/************Using cuts from lambda spectra******************/
/*    
    StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(5);
    StKFParticleInterface::instance()->SetLCut(1);
    StKFParticleInterface::instance()->SetChi2Cut2D(20);
    StKFParticleInterface::instance()->SetChiPrimaryCut(5);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
    StKFParticleInterface::instance()->SetLdLCut2D(3);
*/


/*******************BS Cuts****************/
/*
    StKFParticleInterface::instance()->SetChi2Cut2D(5);
    StKFParticleInterface::instance()->SetChiPrimaryCut(10);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
*/
/************** BS Cut 2 *******************/
/*    StKFParticleInterface::instance()->SetChi2Cut2D(5);
    StKFParticleInterface::instance()->SetChiPrimaryCut(12);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
*/

/************** BS Cut 3 *******************/
    StKFParticleInterface::instance()->SetChi2Cut2D(5);
         StKFParticleInterface::instance()->SetChiPrimaryCut(15);
             StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
             


/************** BS Cut 4 *******************/
/*    StKFParticleInterface::instance()->SetChi2Cut2D(5);
    StKFParticleInterface::instance()->SetChiPrimaryCut(5);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
*/

/************** BS Cut 5 *******************/ //This might have been a wrong cut 
//    StKFParticleInterface::instance()->SetChi2Cut2D(10);
//    StKFParticleInterface::instance()->SetChiPrimaryCut(7);
//    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);


/************ BS Cut 6 *********************/
   /* StKFParticleInterface::instance()->SetChi2Cut2D(7);
    StKFParticleInterface::instance()->SetChiPrimaryCut(10);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
    */

/************ BS Cut 7 *********************/
   /* StKFParticleInterface::instance()->SetChi2Cut2D(3);
    StKFParticleInterface::instance()->SetChiPrimaryCut(10);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(5);
*/

/************ BS Cut 8 *********************/ 
/*    StKFParticleInterface::instance()->SetChi2Cut2D(5);
    StKFParticleInterface::instance()->SetChiPrimaryCut(10);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(7);
*/


  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);
  //StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122); //Including antilambdas
  Long64_t nevent = NumEvents;
    StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
    if (! maker) return;
    maker->SetStatus("*",1);
    TChain *tree = maker->chain();
    Long64_t nentries = tree->GetEntries();
    if (nentries <= 0) return;
    nevent = TMath::Min(nevent,nentries);
    cout << nentries << " events in chain " << nevent << " will be read." << endl;
  
  chain->EventLoop(nevent);
#endif
  
}
