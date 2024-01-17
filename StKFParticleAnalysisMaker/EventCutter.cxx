// See the notes at the bottom of "EventCutter.h"
#include "EventCutter.h"

// private:
void EventCutter::_my_FindRelevantRunLists(unsigned int RunNumber){ // Looks at the run number and finds to which run it belongs
    for( unsigned int iChronologicalRunList=0; iChronologicalRunList<_my_ChronologicalBadRunLists.size(); iChronologicalRunList++ ){
        for( unsigned int BadRunNumber:_my_ChronologicalBadRunLists[iChronologicalRunList] ){
            if( RunNumber==BadRunNumber ){
                _my_ChronologicalBadRunList = _my_ChronologicalBadRunLists[iChronologicalRunList];
                _my_ChronologicalGoodRunList = _my_ChronologicalGoodRunLists[iChronologicalRunList];
                return;
            }
        } // BadRunNumber
        for( unsigned int GoodRunNumber:_my_ChronologicalGoodRunLists[iChronologicalRunList] ){
            if( RunNumber==GoodRunNumber ){
                _my_ChronologicalBadRunList = _my_ChronologicalBadRunLists[iChronologicalRunList];
                _my_ChronologicalGoodRunList = _my_ChronologicalGoodRunLists[iChronologicalRunList];
                return;
            }
        } // GoodRunNumber
    } // for all good/bad run lists
    _my_MiscFunctions->PrintError("_my_FindRelevantRunLists: Can't find the run list");
} // _my_FindRelevantRunLists

// public:
EventCutter::EventCutter(unsigned int NumCutVars, std::vector<std::vector<float>> CentralityUpperLimits){
    TH1::AddDirectory(false); std::cout<<"EventCutter is setting TH1::AddDirectory(false)"<<std::endl; // allows this class to write histograms to the output file

    _my_NumCutVars = NumCutVars;

    _my_MiscFunctions = new MiscFunctions();
    _my_MiscFunctions->SetClassNameForErrorMessage("EventCutter");
    _my_MiscFunctions->DontPrintSameErrorMessageMoreThanOnce();
    if( _my_ChronologicalBadRunLists.size()!=_my_ChronologicalBadRunLists.size() ) _my_MiscFunctions->PrintError("_my_FindRelevantRunLists: Missing a run list");

    _my_h1is_CutVar_nEvts_BforeAllCuts.resize(NumCutVars);
    _my_h1is_CutVar_nEvts_AfterAllCuts.resize(NumCutVars);
    _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut.resize(NumCutVars);
    _my_h2is_VarX_VarY_nEvts_BforeAllCuts.resize(NumCutVars);
    _my_h2is_VarX_VarY_nEvts_AfterAllCuts.resize(NumCutVars);
    _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut.resize(NumCutVars);

    _my_CutVarVals.resize(NumCutVars);
    _my_VarValsXY.resize(NumCutVars);
    for( int iCutVar=0; iCutVar<(int)NumCutVars; iCutVar++ ){
        _my_CutVarMinima.push_back(0.);
        _my_CutVarMaxima.push_back(0.);
        _my_CutVarValIsGood.push_back(true);
        _my_CutWasCalled.push_back(false);
        _my_Cut2DWasCalled.push_back(false);

        _my_VarValsXY[iCutVar].resize(2);
    } // iCutVar

    _my_CentralityUpperLimits = CentralityUpperLimits;
    _my_NumCentralityBinLists = _my_CentralityUpperLimits.size();

    /*_my_ChronologicalBadRunLists = { // Note: These lists MUST be ordered from smallest to largest; other functions depend on it!!
        { // production_3p85GeV_fixedTarget_2018
        19151030, 19151035, 19151037, 19151038, 19151040, 19151042, 19151051, 19151057, 19151064, 19152012, 19152015, 19152026, 19152049, 19152050, 
        19152072, 19152077, 19152079, 19152080, 19153005, 19153006, 19153008, 19153030, 19153060, 19154008, 19154009, 19154010, 19154011, 19154043, 
        19154050, 19154062, 19155002, 19155007, 19155012, 19155013, 19155014, 19155015
        } */ // production_3p85GeV_fixedTarget_2018
    _my_ChronologicalBadRunLists = { // Note: These lists MUST be ordered from smallest to largest; other functions depend on it!!
        { // production_3p85GeV_fixedTarget_2018
        19151029, 19151031, 19151034, 19151036, 19151039, 19151041, 19151043,
		19151044, 19151045, 19151046, 19151047, 19151048, 19151049, 19151050,
		19151052, 19151053, 19151054, 19151055, 19151056, 19151066, 19151067,
		19151068, 19151069, 19151070, 19151071, 19151072, 19151082, 19151083,
		19151084, 19152001, 19152002, 19152003, 19152008, 19152009, 19152010,
		19152014, 19152016, 19152021, 19152023, 19152024, 19152025, 19152027,
		19152028, 19152029, 19152030, 19152031, 19152032, 19152033, 19152034,
		19152035, 19152036, 19152037, 19152038, 19152039, 19152040, 19152041,
		19152042, 19152043, 19152044, 19152045, 19152046, 19152048, 19152051,
		19152052, 19152053, 19152054, 19152055, 19152071, 19152073, 19152074,
		19152075, 19152076, 19152078, 19152081, 19153001, 19153002, 19153003,
		19153004, 19153007, 19153009, 19153010, 19153011, 19153012, 19153013,
		19153014, 19153015, 19153016, 19153017, 19153018, 19153019, 19153020,
		19153021, 19153022, 19153023, 19153024, 19153025, 19153027, 19153028,
		19153029, 19153031, 19153032, 19153065, 19154012, 19154013, 19154014,
		19154015, 19154016, 19154017, 19154018, 19154019, 19154020, 19154021,
		19154022, 19154023, 19154024, 19154026, 19154046, 19154051, 19154056
        } // production_3p85GeV_fixedTarget_2018
        ,
        { // production_19GeV_2019
        20093006, 20093008, 20093020, 20093031, 20093032, 20093034, 20092034, 20092043, 20092045, 20091019, 20090027, 20090028, 20090029, 20090030, 
        20090032, 20090033, 20090034, 20090045, 20088007, 20086013, 20086052, 20082023, 20080005, 20080021, 20078008, 20078009, 20078010, 20077032, 
        20076060, 20075010, 20075014, 20075017, 20075022, 20075041, 20073004, 20073009, 20073010, 20073012, 20073014, 20073066, 20072040, 20072042, 
        20072049, 20071032, 20071034, 20070001, 20070036, 20068018, 20068021, 20068024, 20068026, 20066014, 20066067, 20066068, 20066069, 20066077, 
        20065002, 20065011, 20064041, 20063004, 20062007, 20062009, 20062010, 20062011, 20062012, 20061026, 20057042
        } // production_19GeV_2019
        //, Just a reminder for a comma here :)
    }; // _my_ChronologicalBadRunLists

    
  /* _my_ChronologicalGoodRunLists = { // Note: These lists MUST be ordered from smallest to largest; other functions depend on it!!
        { // production_3p85GeV_fixedTarget_2018
        19151029, 19151031, 19151034, 19151036, 19151039, 19151041, 19151043, 19151044, 19151045, 19151046, 19151047, 19151048, 19151049, 19151050, 
        19151052, 19151053, 19151054, 19151055, 19151056, 19151066, 19151067, 19151068, 19151069, 19151070, 19151071, 19151072, 19151082, 19151083, 
        19151084, 19152001, 19152002, 19152003, 19152008, 19152009, 19152010, 19152014, 19152016, 19152021, 19152023, 19152024, 19152025, 19152027, 
        19152028, 19152029, 19152030, 19152031, 19152032, 19152033, 19152034, 19152035, 19152036, 19152037, 19152038, 19152039, 19152040, 19152041, 
        19152042, 19152043, 19152044, 19152045, 19152046, 19152048, 19152051, 19152052, 19152053, 19152054, 19152055, 19152071, 19152073, 19152074, 
        19152075, 19152076, 19152078, 19152081, 19153001, 19153002, 19153003, 19153004, 19153007, 19153009, 19153010, 19153011, 19153012, 19153013, 
        19153014, 19153015, 19153016, 19153017, 19153018, 19153019, 19153020, 19153021, 19153022, 19153023, 19153024, 19153025, 19153027, 19153028, 
        19153029, 19153031, 19153032, 19153033, 19153034, 19153035, 19153036, 19153037, 19153042, 19153043, 19153044, 19153050, 19153051, 19153052, 
        19153053, 19153054, 19153055, 19153056, 19153057, 19153058, 19153059, 19153061, 19153062, 19153063, 19153064, 19153065, 19153066, 19154001, 
        19154002, 19154005, 19154007, 19154012, 19154013, 19154014, 19154015, 19154016, 19154017, 19154018, 19154019, 19154020, 19154021, 19154022, 
        19154023, 19154024, 19154026, 19154027, 19154028, 19154029, 19154030, 19154031, 19154032, 19154036, 19154037, 19154038, 19154039, 19154040, 
        19154041, 19154044, 19154045, 19154046, 19154047, 19154048, 19154049, 19154051, 19154052, 19154053, 19154054, 19154055, 19154056, 19154057, 
        19154058, 19154061, 19154063, 19154064, 19154065, 19154066, 19154067, 19155001, 19155003, 19155004, 19155005, 19155006, 19155008, 19155009, 
        19155010, 19155011, 19155016, 19155017, 19155018, 19155019, 19155020, 19155021, 19155022
        }*/ // production_3p85GeV_fixedTarget_2018   from Joey's code 
    _my_ChronologicalGoodRunLists = { // Note: These lists MUST be ordered from smallest to largest; other functions depend on it!!
        { // production_3p85GeV_fixedTarget_2018
        19153029, 19153031, 19153033, 19153034, 19153035, 19153036, 19153037, 19153042,
        19153043, 19153044, 19153050, 19153051, 19153052, 19153053, 19153054, 19153055,
        19153056, 19153057, 19153058, 19153059, 19153061, 19153062, 19153063, 19153064,
        19153066, 19154001, 19154002, 19154005, 19154007, 19154027, 19154028, 19154029,
        19154030, 19154031, 19154032, 19154036, 19154037, 19154038, 19154039, 19154040,
        19154041, 19154044, 19154045, 19154046, 19154047, 19154048, 19154049, 19154052,
        19154053, 19154054, 19154055, 19154056, 19154057, 19154058, 19154061, 19154063,
        19154064, 19154065, 19154066, 19154067, 19155001, 19155003, 19155004, 19155005,
        19155006, 19155008, 19155009, 19155010, 19155011, 19155016, 19155017, 19155018,
        19155019, 19155020, 19155021, 19155022
        } // production_3p85GeV_fixedTarget_2018   from Good Run Fluctuation Analyses From Cebra
        ,
        { // production_19GeV_2019
        20056032, 20056040, 20057003, 20057004, 20057005, 20057006, 20057007, 20057008, 20057009, 20057011, 20057012, 20057013, 20057014, 20057015, 
        20057016, 20057025, 20057026, 20057027, 20057028, 20057029, 20057032, 20057033, 20057037, 20057038, 20057039, 20057040, 20057043, 20057044, 
        20057046, 20057047, 20057048, 20057049, 20057050, 20058001, 20058002, 20058003, 20058004, 20058005, 20059001, 20059016, 20059017, 20059018, 
        20059060, 20059061, 20059062, 20059073, 20060001, 20060003, 20060004, 20060005, 20060012, 20060017, 20060018, 20060019, 20060020, 20060022, 
        20060023, 20060024, 20060025, 20060026, 20060027, 20060045, 20060046, 20060059, 20060060, 20060061, 20060062, 20060063, 20060064, 20060065, 
        20060066, 20060067, 20060068, 20060069, 20060070, 20061001, 20061003, 20061004, 20061005, 20061006, 20061007, 20061008, 20061009, 20061015, 
        20061016, 20061017, 20061018, 20061021, 20061023, 20061024, 20061028, 20061029, 20061030, 20061031, 20061032, 20061033, 20061034, 20061035, 
        20061036, 20061037, 20061038, 20061039, 20061040, 20061042, 20061043, 20061045, 20062001, 20062002, 20062003, 20062004, 20062005, 20062006, 
        20062013, 20062014, 20062015, 20062017, 20062018, 20062019, 20062036, 20062037, 20062038, 20062041, 20062042, 20062043, 20062044, 20062045, 
        20062046, 20062047, 20062048, 20062049, 20062051, 20062054, 20062055, 20062056, 20062057, 20062059, 20063002, 20063005, 20063007, 20063008, 
        20063009, 20063010, 20063011, 20063012, 20063013, 20063015, 20063017, 20063018, 20063019, 20063020, 20063021, 20063023, 20063028, 20063029, 
        20063032, 20063033, 20063034, 20063035, 20063036, 20063037, 20063039, 20063040, 20063042, 20063043, 20063044, 20063048, 20063050, 20063052, 
        20063053, 20063055, 20063057, 20063059, 20064002, 20064004, 20064005, 20064006, 20064007, 20064008, 20064009, 20064011, 20064012, 20064013, 
        20064014, 20064015, 20064016, 20064017, 20064020, 20064022, 20064025, 20064030, 20064032, 20064033, 20064034, 20064037, 20064038, 20064039, 
        20064040, 20064042, 20064044, 20064045, 20064047, 20064050, 20064051, 20064052, 20065001, 20065004, 20065005, 20065006, 20065007, 20065008, 
        20065009, 20065010, 20065012, 20065015, 20065017, 20065018, 20065019, 20065056, 20065057, 20065058, 20065060, 20065061, 20065062, 20065063, 
        20066001, 20066003, 20066004, 20066005, 20066006, 20066007, 20066008, 20066009, 20066010, 20066012, 20066013, 20066015, 20066017, 20066018, 
        20066019, 20066020, 20066021, 20066023, 20066024, 20066025, 20066026, 20066027, 20066028, 20066070, 20066072, 20066073, 20066074, 20066076, 
        20066078, 20066080, 20066081, 20066082, 20067001, 20067002, 20067003, 20067004, 20067005, 20067007, 20067008, 20067009, 20067010, 20067011, 
        20067012, 20067013, 20067014, 20067015, 20067016, 20067017, 20067018, 20067019, 20067020, 20067022, 20067023, 20067024, 20067028, 20067029, 
        20067030, 20067031, 20067038, 20067039, 20067040, 20067041, 20067044, 20067045, 20067046, 20067047, 20067048, 20067049, 20067051, 20068001, 
        20068002, 20068003, 20068004, 20068005, 20068007, 20068008, 20068009, 20068010, 20068011, 20068012, 20068013, 20068014, 20068017, 20068019, 
        20068020, 20068022, 20068023, 20068025, 20068028, 20068029, 20068030, 20068034, 20068035, 20068036, 20068037, 20068049, 20068050, 20068051, 
        20068052, 20068053, 20068055, 20068056, 20068057, 20068058, 20068059, 20068060, 20068062, 20068063, 20068064, 20068065, 20068066, 20069001, 
        20069002, 20069003, 20069004, 20069005, 20069006, 20069007, 20069008, 20069009, 20069010, 20069011, 20069020, 20069021, 20069022, 20069023, 
        20069024, 20069025, 20069026, 20069027, 20069028, 20069029, 20069030, 20069031, 20069032, 20069033, 20069034, 20069035, 20069036, 20069037, 
        20069038, 20069042, 20069043, 20069044, 20069045, 20069046, 20069050, 20069051, 20069052, 20069053, 20069054, 20069055, 20069057, 20069058, 
        20069059, 20069060, 20069061, 20069062, 20070002, 20070003, 20070004, 20070005, 20070006, 20070007, 20070010, 20070011, 20070012, 20070013, 
        20070014, 20070015, 20070016, 20070017, 20070018, 20070019, 20070020, 20070021, 20070037, 20070038, 20070039, 20070040, 20070041, 20070042, 
        20070043, 20070044, 20070045, 20070047, 20071001, 20071003, 20071004, 20071005, 20071006, 20071007, 20071008, 20071009, 20071010, 20071011, 
        20071012, 20071013, 20071014, 20071015, 20071016, 20071017, 20071018, 20071020, 20071021, 20071022, 20071024, 20071025, 20071026, 20071027, 
        20071029, 20071030, 20071031, 20071036, 20071037, 20071038, 20071040, 20071041, 20071042, 20071043, 20071044, 20071045, 20071046, 20071047, 
        20071048, 20071049, 20071050, 20071051, 20071052, 20071053, 20071054, 20071055, 20071056, 20071057, 20071058, 20071059, 20071061, 20071062, 
        20071063, 20071065, 20072001, 20072002, 20072003, 20072004, 20072005, 20072006, 20072007, 20072008, 20072009, 20072010, 20072011, 20072012, 
        20072013, 20072014, 20072016, 20072017, 20072018, 20072034, 20072035, 20072036, 20072037, 20072038, 20072039, 20072041, 20072045, 20072046, 
        20072047, 20072048, 20072050, 20072051, 20072052, 20072054, 20072055, 20072056, 20073001, 20073002, 20073003, 20073005, 20073006, 20073007, 
        20073008, 20073011, 20073013, 20073015, 20073016, 20073017, 20073018, 20073019, 20073021, 20073022, 20073023, 20073024, 20073025, 20073026, 
        20073027, 20073071, 20073072, 20073074, 20073076, 20074001, 20074002, 20074003, 20074004, 20074005, 20074006, 20074007, 20074008, 20074009, 
        20074010, 20074011, 20074012, 20074014, 20074016, 20074017, 20074018, 20074019, 20074020, 20074021, 20074023, 20074026, 20074027, 20074029, 
        20074030, 20074032, 20074033, 20074034, 20074043, 20074044, 20074045, 20074046, 20075001, 20075002, 20075004, 20075006, 20075007, 20075008, 
        20075009, 20075011, 20075012, 20075013, 20075015, 20075018, 20075019, 20075020, 20075021, 20075023, 20075027, 20075028, 20075031, 20075032, 
        20075033, 20075034, 20075035, 20075036, 20075037, 20075039, 20075040, 20075042, 20075043, 20075044, 20075046, 20075047, 20075048, 20075049, 
        20075050, 20075051, 20075052, 20075054, 20075055, 20075056, 20075057, 20075058, 20075059, 20075060, 20075061, 20075062, 20075063, 20075064, 
        20075065, 20075066, 20075067, 20075068, 20076001, 20076002, 20076003, 20076004, 20076005, 20076006, 20076007, 20076008, 20076009, 20076010, 
        20076011, 20076012, 20076013, 20076014, 20076015, 20076016, 20076017, 20076018, 20076019, 20076021, 20076022, 20076023, 20076025, 20076026, 
        20076027, 20076028, 20076029, 20076030, 20076031, 20076032, 20076033, 20076034, 20076035, 20076036, 20076037, 20076038, 20076039, 20076040, 
        20076042, 20076043, 20076045, 20076046, 20076047, 20076048, 20076049, 20076050, 20076051, 20076052, 20076053, 20076054, 20076055, 20076056, 
        20076057, 20076058, 20076059, 20076061, 20076063, 20077001, 20077002, 20077003, 20077004, 20077005, 20077006, 20077007, 20077008, 20077009, 
        20077010, 20077011, 20077012, 20077013, 20077014, 20077015, 20077016, 20077017, 20077018, 20077019, 20077020, 20077021, 20077034, 20077035, 
        20078001, 20078003, 20078004, 20078005, 20078006, 20078007, 20078011, 20078012, 20078013, 20078014, 20078015, 20078016, 20078017, 20078018, 
        20078019, 20078020, 20078021, 20078022, 20078023, 20078024, 20078026, 20078028, 20078029, 20078030, 20078032, 20078033, 20078034, 20078035, 
        20078038, 20078039, 20078040, 20078041, 20078042, 20078043, 20078044, 20078045, 20078046, 20078047, 20078048, 20078049, 20078051, 20078052, 
        20078053, 20078054, 20078055, 20078056, 20078057, 20078058, 20078059, 20078060, 20078061, 20078062, 20078063, 20078064, 20078065, 20078067, 
        20079002, 20079004, 20079006, 20079007, 20079008, 20079009, 20079010, 20079011, 20079013, 20079014, 20079015, 20079016, 20079017, 20079018, 
        20079019, 20079020, 20079021, 20079022, 20079023, 20079024, 20079027, 20079028, 20079044, 20079045, 20079046, 20080003, 20080004, 20080006, 
        20080007, 20080008, 20080009, 20080010, 20080011, 20080012, 20080013, 20080014, 20080015, 20080016, 20080017, 20080018, 20080019, 20080020, 
        20080022, 20080023, 20080024, 20080025, 20080030, 20080031, 20080032, 20081001, 20081002, 20081003, 20081004, 20081005, 20081006, 20081007, 
        20081008, 20081009, 20081010, 20081012, 20081013, 20081014, 20081015, 20081016, 20081017, 20081018, 20081019, 20081021, 20081025, 20081026, 
        20081027, 20081032, 20081033, 20082001, 20082002, 20082003, 20082004, 20082005, 20082006, 20082007, 20082008, 20082010, 20082011, 20082012, 
        20082013, 20082014, 20082015, 20082016, 20082017, 20082018, 20082019, 20082020, 20082021, 20082024, 20082025, 20082026, 20082027, 20082029, 
        20082030, 20082031, 20082034, 20082035, 20082036, 20082038, 20082041, 20082043, 20082047, 20082048, 20082049, 20082050, 20082051, 20082052, 
        20082053, 20082054, 20082055, 20082056, 20082057, 20082058, 20082059, 20082060, 20082061, 20082062, 20082063, 20082065, 20082066, 20082067, 
        20082068, 20083001, 20083002, 20083003, 20083004, 20083006, 20083009, 20083010, 20083011, 20083014, 20083017, 20083019, 20083020, 20083021, 
        20083022, 20083023, 20083024, 20083025, 20083026, 20083027, 20083029, 20083030, 20083031, 20083032, 20083033, 20083034, 20083072, 20083073, 
        20083074, 20083075, 20083076, 20083077, 20083078, 20083079, 20083080, 20083081, 20084001, 20084002, 20084003, 20084004, 20084005, 20084007, 
        20084008, 20084009, 20084010, 20084011, 20084012, 20084013, 20084014, 20084015, 20084016, 20084017, 20084018, 20084022, 20084023, 20084024, 
        20085001, 20085002, 20085003, 20085004, 20085005, 20085006, 20085007, 20085008, 20085009, 20085010, 20085011, 20085013, 20085016, 20085017, 
        20085018, 20085050, 20086001, 20086002, 20086003, 20086004, 20086005, 20086006, 20086007, 20086008, 20086009, 20086010, 20086011, 20086012, 
        20086014, 20086015, 20086016, 20086017, 20086046, 20086047, 20086048, 20086051, 20086053, 20087001, 20087002, 20087003, 20087004, 20087005, 
        20087006, 20087007, 20087008, 20087009, 20087010, 20087011, 20087012, 20087014, 20087017, 20087018, 20087019, 20087020, 20087021, 20087022, 
        20087023, 20088001, 20088002, 20088003, 20088004, 20088005, 20088006, 20088008, 20088009, 20088010, 20088011, 20088012, 20088013, 20088014, 
        20088015, 20088016, 20088028, 20088029, 20088030, 20088031, 20088032, 20088033, 20088034, 20088036, 20088037, 20088038, 20089002, 20089003, 
        20089004, 20089005, 20089006, 20089007, 20089008, 20089009, 20089010, 20089011, 20089012, 20089013, 20089014, 20089015, 20089016, 20089017, 
        20089018, 20089019, 20089020, 20089026, 20089027, 20089028, 20089029, 20090001, 20090002, 20090003, 20090004, 20090005, 20090006, 20090007, 
        20090008, 20090009, 20090010, 20090011, 20090012, 20090013, 20090014, 20090015, 20090016, 20090017, 20090018, 20090019, 20090020, 20090021, 
        20090022, 20090023, 20090024, 20090025, 20090026, 20090031, 20090035, 20090046, 20090047, 20090048, 20091001, 20091002, 20091003, 20091004, 
        20091005, 20091006, 20091007, 20091008, 20091009, 20091010, 20091011, 20091012, 20091013, 20091014, 20091015, 20091016, 20091017, 20091018, 
        20091020, 20091021, 20091022, 20091023, 20092003, 20092004, 20092005, 20092006, 20092007, 20092008, 20092009, 20092010, 20092011, 20092012, 
        20092013, 20092014, 20092015, 20092016, 20092017, 20092018, 20092019, 20092020, 20092021, 20092022, 20092023, 20092024, 20092025, 20092026, 
        20092027, 20092028, 20092029, 20092030, 20092031, 20092032, 20092033, 20092035, 20092036, 20092037, 20092038, 20092039, 20092044, 20092050, 
        20092051, 20092052, 20092053, 20092054, 20092055, 20092056, 20092057, 20092060, 20092061, 20093001, 20093002, 20093003, 20093005, 20093009, 
        20093010, 20093011, 20093012, 20093016, 20093017, 20093018, 20093024, 20093025, 20093035, 20093036
        } // production_19GeV_2019
        //, Just a reminder for a comma here :)
    }; // _my_ChronologicalBadRunLists
    // if your run (energy, production) list is not here, just add it here! Nothing more to do

} // EventCutter

EventCutter::~EventCutter(){
    for( int iCutVar=0; iCutVar<(int)_my_NumCutVars; iCutVar++ ){
        delete _my_h1is_CutVar_nEvts_BforeAllCuts[iCutVar];
        delete _my_h1is_CutVar_nEvts_AfterAllCuts[iCutVar];
        delete _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut[iCutVar];
        delete _my_h2is_VarX_VarY_nEvts_BforeAllCuts[iCutVar];
        delete _my_h2is_VarX_VarY_nEvts_AfterAllCuts[iCutVar];
        delete _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut[iCutVar];
    } // iCutVar
    delete _my_CutVarDistHistos;
    delete _my_RefMultDistributionsSeparateCentralities;
} // ~EventCutter

void EventCutter::InstantiateHists(unsigned int CutVarId, TString CutVarName, TString CutVarTitle, unsigned int NumBins,  float HistMin,  float HistMax,
                                             TString VarNameX,   TString VarTitleX,   unsigned int NumBinsX, float HistMinX, float HistMaxX,
                                             TString VarNameY,   TString VarTitleY,   unsigned int NumBinsY, float HistMinY, float HistMaxY){
    InstantiateHists(CutVarId, CutVarName, CutVarTitle, NumBins, HistMin, HistMax);
    _my_h2is_VarX_VarY_nEvts_BforeAllCuts[CutVarId] = new TH2I(
        Form("%s_%s_nEvts_BforeAllCuts",VarNameX.Data(),VarNameY.Data()),
        Form("Before any cuts;%s;%s;# Events",VarTitleX.Data(),VarTitleY.Data()),
        NumBinsX,HistMinX,HistMaxX,NumBinsY,HistMinY,HistMaxY);
    _my_h2is_VarX_VarY_nEvts_AfterAllCuts[CutVarId] = new TH2I(
        Form("%s_%s_nEvts_AfterAllCuts",VarNameX.Data(),VarNameY.Data()),
        Form("After all cuts;%s;%s;# Events",VarTitleX.Data(),VarTitleY.Data()),
        NumBinsX,HistMinX,HistMaxX,NumBinsY,HistMinY,HistMaxY);
    _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut[CutVarId] = new TH2I(
        Form("%s_%s_nEvts_AfterAllCutsBut%sCut",VarNameX.Data(),VarNameY.Data(),CutVarName.Data()),
        Form("After all cuts except %s cut;%s;%s;# Events",CutVarTitle.Data(),VarTitleX.Data(),VarTitleY.Data()),
        NumBinsX,HistMinX,HistMaxX,NumBinsY,HistMinY,HistMaxY);
} // InstantiateHists

void EventCutter::InstantiateHists(unsigned int CutVarId, TString CutVarName, TString CutVarTitle, unsigned int NumBins, float HistMin, float HistMax){
    _my_h1is_CutVar_nEvts_BforeAllCuts[CutVarId] = new TH1I(
        Form("%s_nEvts_BforeAllCuts",CutVarName.Data()),
        Form("Before any cuts;%s;# Events",CutVarTitle.Data()),
        NumBins,HistMin,HistMax);
    _my_h1is_CutVar_nEvts_AfterAllCuts[CutVarId] = new TH1I(
        Form("%s_nEvts_AfterAllCuts",CutVarName.Data()),
        Form("After all cuts;%s;# Events",CutVarTitle.Data()),
        NumBins,HistMin,HistMax);
    _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut[CutVarId] = new TH1I(
        Form("%s_nEvts_AfterAllCutsBut%sCut"/* NOTE! The constructor uses this naming convention for reading the RefMult histogram! */,CutVarName.Data(),CutVarName.Data()),
        Form("After all cuts except %s cut;%s;# Events",CutVarTitle.Data(),CutVarTitle.Data()),
        NumBins,HistMin,HistMax);
} // InstantiateHistograms

void EventCutter::SetCentralityLimits(unsigned int RefMultCutVarId, TString RefMultHistogramFileName, TString RefMultVariableNameInHistogram, int RefMultMax /* use -1 for no limit */){
    // This function can take in a RefMult distribution(*) and set _my_CutVarMinima[RefMultCutVarId] to reject >80% central events and set _my_CutVarMaxima[RefMultCutVarId]
    // to reject the user-input RefMultMax (which will likely either be set to reject pileup or set to -1 for "no" upper limit on RefMult)
    // (*) This function looks for a histogram with the RefMult variable's name appended with "_nEvts_AfterAllCutsBut..." (see a few lines below). If you want it to read a 
    // histogram you previously created using your own naming scheme, you'll have to edit this.

    // Load the file containing the RefMult histogram and load the histogram
    TFile * RefMultHistogramFile = new TFile(RefMultHistogramFileName,"READ");
    if( !RefMultHistogramFile ){
        _my_MiscFunctions->PrintError("EventCutter::LoadRefMultHistogram could not find your RefMult histogram file!");
        return;
    }
    TDirectory * CutVariableDistributionHistograms = (TDirectory *)RefMultHistogramFile->Get("EventCutter_CutVariableDistributionHistograms"); // Find the histogram in this directory
    TH1I * h1i_RefMult_nEvts_ReadIn = (TH1I *)CutVariableDistributionHistograms->Get(Form("%s_nEvts_AfterAllCutsBut%sCut",RefMultVariableNameInHistogram.Data(),RefMultVariableNameInHistogram.Data()));
    if( !h1i_RefMult_nEvts_ReadIn ){
        _my_MiscFunctions->PrintError("See below");
        std::cout<<"EventCutter::LoadRefMultHistogram could not find your RefMult histogram! The second argument should be something like \"RefMult\". Your options are:"<<std::endl;
        // look at histograms in file with name containing "_nEvts_AfterAllCutsBut"; if you made a RefMult distribution with this class, you'll have a RefMult distribution with this format
            TIter CutVariableHistogramDirectoryKeyIterator(CutVariableDistributionHistograms->GetListOfKeys());
            TKey * CutVariableHistogramDirectoryKey;
            while( (CutVariableHistogramDirectoryKey = (TKey *)CutVariableHistogramDirectoryKeyIterator()) ){
                if( !CutVariableHistogramDirectoryKey->IsFolder() ){ // doing the next line (CutVariableHistogramDirectoryKey->ReadObj()) on any sub-directories will screw things up
                    TObject * Object = CutVariableHistogramDirectoryKey->ReadObj() ;
                    if( Object->InheritsFrom("TH1") ){
                        TString VariableNameFromHist = _my_MiscFunctions->DelimitString(Object->GetName(),"_nEvts_AfterAllCutsBut",0);
                        if( !VariableNameFromHist.EqualTo("") ) std::cout<<VariableNameFromHist<<std::endl;
                    } // TH1
                } // anything that's not a directory
            } // loop over CutVariableHistogramDirectoryKeys
            std::cout<<std::endl;
            return;
    } // if can't find histogram

    // Check that the RefMult histogram's binning makes sense, and note whether it starts at zero or one
    if( h1i_RefMult_nEvts_ReadIn->GetBinWidth(1)!=1 ){
        _my_MiscFunctions->PrintError("Your RefMult histogram doesn't have bin widths of 1; this is nonsensical and you'll have to find centrality yourself.");
        return;
    }
    int RefMultForFirstBin;
    if( h1i_RefMult_nEvts_ReadIn->FindBin(0)==1) RefMultForFirstBin = 0;
    else if( h1i_RefMult_nEvts_ReadIn->FindBin(1)==1 ) RefMultForFirstBin = 1;
    else{
        _my_MiscFunctions->PrintError("Your RefMult histogram's first bin should contain either 0 or 1; you'll have to find centrality yourself.");
        return;
    }
    int NumRefMultBins = h1i_RefMult_nEvts_ReadIn->GetNbinsX();
    if( h1i_RefMult_nEvts_ReadIn->GetBinContent(NumRefMultBins)!=0 ){
        _my_MiscFunctions->PrintError("Your RefMult histogram's last bin has contents. Consider exanding the RefMult axis.");
    }

    int NumEvtsNoPileup = 0; // Number of entries in your RefMult histogram that aren't pileup
    if( RefMultMax==-1 || RefMultMax>NumRefMultBins ) RefMultMax = h1i_RefMult_nEvts_ReadIn->GetNbinsX()+1; // +1 in case there's anything in the overflow bin
    for( int iRefMultBin=RefMultMax+(-1*RefMultForFirstBin+1); iRefMultBin>0; iRefMultBin-- ) NumEvtsNoPileup += h1i_RefMult_nEvts_ReadIn->GetBinContent(iRefMultBin);

    float RunningSumOfRefMultHistogramBinContents = 0;
    std::vector<unsigned int> CentralityBinWeAreTryingToFind;
    for( unsigned int i=0; i<_my_NumCentralityBinLists; i++ ) CentralityBinWeAreTryingToFind.push_back(0);
    _my_RefMultMinimaLimits.resize(_my_NumCentralityBinLists);
    for( int iRefMultBin=RefMultMax+(-1*RefMultForFirstBin+1); iRefMultBin>0; iRefMultBin-- ){
        RunningSumOfRefMultHistogramBinContents += h1i_RefMult_nEvts_ReadIn->GetBinContent(iRefMultBin);
        float ThisSumsFractionOfTotalEvts = 100.*RunningSumOfRefMultHistogramBinContents/NumEvtsNoPileup;
        float NextSumsFractionOfTotalEvts = 100.*(RunningSumOfRefMultHistogramBinContents+h1i_RefMult_nEvts_ReadIn->GetBinContent(iRefMultBin+1))/NumEvtsNoPileup;
        for( unsigned int iCentralityBinList=0; iCentralityBinList<_my_NumCentralityBinLists; iCentralityBinList++ ){
            if( CentralityBinWeAreTryingToFind[iCentralityBinList]<_my_CentralityUpperLimits[iCentralityBinList].size() ){
                float ThisDistanceInCentralityToNextCentralityBin = fabs(ThisSumsFractionOfTotalEvts-_my_CentralityUpperLimits[iCentralityBinList][CentralityBinWeAreTryingToFind[iCentralityBinList]]);
                float NextDistanceInCentralityToNextCentralityBin = fabs(NextSumsFractionOfTotalEvts-_my_CentralityUpperLimits[iCentralityBinList][CentralityBinWeAreTryingToFind[iCentralityBinList]]);
                if( ThisDistanceInCentralityToNextCentralityBin<NextDistanceInCentralityToNextCentralityBin ){
                    _my_RefMultMinimaLimits[iCentralityBinList].push_back(iRefMultBin+(RefMultForFirstBin-1));
                    CentralityBinWeAreTryingToFind[iCentralityBinList]++;
                } // if this is a Centrality's RefMult Cutoff
            } // if we haven't reached >80% centrality yet
        } // iCentralityBinList
    } // iRefMultBin
    if( _my_RefMultMinimaLimits.size()>0 && _my_RefMultMinimaLimits[0].size()>0){ // just a safeguard
        SetCutVarBounds(RefMultCutVarId,_my_RefMultMinimaLimits[0][_my_RefMultMinimaLimits[0].size()-1],RefMultMax);
    }

    _my_h1is_RefMult_nEvts_SeparateCentralities.resize(_my_NumCentralityBinLists);
    for( unsigned int iCentralityBinList=0; iCentralityBinList<_my_NumCentralityBinLists; iCentralityBinList++ ){
        for( int iCentrality=0; iCentrality<(int)_my_CentralityUpperLimits[iCentralityBinList].size(); iCentrality++ ){
            int CentralityLowerLimit = (iCentrality==0)?0:_my_CentralityUpperLimits[iCentralityBinList][iCentrality-1];
            int CentralityUpperLimit = (int)_my_CentralityUpperLimits[iCentralityBinList][iCentrality];
            _my_h1is_RefMult_nEvts_SeparateCentralities[iCentralityBinList].push_back(new TH1I(
                Form("RefMult_nEvts_%ito%iPercentCentrality_%iCentralityBins",CentralityLowerLimit,CentralityUpperLimit,_my_CentralityUpperLimits[iCentralityBinList].size()),
                Form("%i-%i%% Centrality;;# Events",CentralityLowerLimit,CentralityUpperLimit),
                NumRefMultBins,h1i_RefMult_nEvts_ReadIn->GetBinLowEdge(1),h1i_RefMult_nEvts_ReadIn->GetBinLowEdge(NumRefMultBins+1)));
        } // iCentrality
    } // iCentralityBinList

} // SetCentralityLimits

float EventCutter::GetCentrality(unsigned int RefMult, unsigned int NumCentralityBins){
    // see which of the centrality definitions (_my_CentralityUpperLimits) the user wants from NumCentralityBins
    int ThisCentralityBinList = -1;
    for( unsigned int iCentralityBinList=0; iCentralityBinList<_my_NumCentralityBinLists; iCentralityBinList++ ){
        if( NumCentralityBins==_my_CentralityUpperLimits[iCentralityBinList].size() ){
            ThisCentralityBinList = iCentralityBinList;
            break;
        }
    }
    if( ThisCentralityBinList==-1 ){
       _my_MiscFunctions->PrintError("You passed GetCentrality an invalid second argument (number of centrality bins). See _my_CentralityUpperLimits in this class's constructor");
       return -1;
   }

    for( unsigned int iCentrality=0; iCentrality<_my_CentralityUpperLimits[ThisCentralityBinList].size(); iCentrality++ ){
        if( RefMult>=_my_RefMultMinimaLimits[ThisCentralityBinList][iCentrality] ){
            _my_h1is_RefMult_nEvts_SeparateCentralities[ThisCentralityBinList][iCentrality]->Fill(RefMult); // Should make sure this is filled only once per event (in case user calls this multiple times)
            if( iCentrality==0 ) return _my_CentralityUpperLimits[ThisCentralityBinList][0]/2.; // so, e.g. 0-5% centrality returns 2.5
            else return (_my_CentralityUpperLimits[ThisCentralityBinList][iCentrality]+_my_CentralityUpperLimits[ThisCentralityBinList][iCentrality-1])/2.;
        }
    } // iCentrality
    return 90; // This shouldn't happen
}

void EventCutter::SetCutVarBounds(unsigned int CutVarId, float CutVarMin, float CutVarMax){
    _my_CutVarMinima[CutVarId] = CutVarMin;
    _my_CutVarMaxima[CutVarId] = CutVarMax;
} // SetCutVarBounds

bool EventCutter::RunIsMarkedBad(unsigned int RunNumber){
    if( RunNumber==0 ){ _my_MiscFunctions->PrintError("RunIsMarkedBad: RunNumber is zero; calling this bad."); return true; }
    if( RunNumber==_my_PreviousRunNumberForFunctionRunIsMarkedBad ) return _my_PreviousRunWasMarkedBad;
    _my_PreviousRunNumberForFunctionRunIsMarkedBad = RunNumber; // new run number! Update this
    if( _my_ChronologicalGoodRunList.size()==0 ) _my_FindRelevantRunLists(RunNumber);
    for( unsigned int BadRunNumber:_my_ChronologicalBadRunList ){
        if( RunNumber>BadRunNumber ){
            _my_PreviousRunWasMarkedBad = false;
            return false;
        }
        if( RunNumber==BadRunNumber ){
            _my_PreviousRunWasMarkedBad = true;
            return true;
        }
    } // BadRunNumber
    _my_PreviousRunWasMarkedBad = false; // you need this since a good event number larger than the last element of the "bad" array will not return above
    return false;
} // RunIsMarkedBad

int EventCutter::ChronologicalRunId(unsigned int RunNumber){
    if( RunNumber==_my_PreviousRunNumberForFunctionChronologicalRunId ) return _my_PreviousChronologicalRunIdForFunctionChronologicalRunId;
    _my_PreviousRunNumberForFunctionChronologicalRunId = RunNumber;
    if( _my_ChronologicalGoodRunList.size()==0 ) _my_FindRelevantRunLists(RunNumber);
    for( unsigned int iGoodRunNumber=0; iGoodRunNumber<=_my_ChronologicalGoodRunList.size(); iGoodRunNumber++ ){
        if( RunNumber==_my_ChronologicalGoodRunList[iGoodRunNumber] ){
            int ChronologicalRunId = iGoodRunNumber+1;
            _my_PreviousChronologicalRunIdForFunctionChronologicalRunId = ChronologicalRunId;
            return ChronologicalRunId;
        }
    } // GoodRunNumber
    _my_MiscFunctions->PrintError("ChronologicalRunId: returning -1"); // you probably didn't reject events marked as bad
    return -1; // this shouldn't happen
} // RunIsMarkedBad

void EventCutter::Cut(unsigned int CutVarId, float CutVarVal, float VarValX, float VarValY){
    _my_Cut2DWasCalled[CutVarId] = true;
    _my_h2is_VarX_VarY_nEvts_BforeAllCuts[CutVarId]->Fill(VarValX,VarValY); // don't check that this has been instantiated; easier to trace back w/ seg. fault
    _my_VarValsXY[CutVarId][0] = VarValX;
    _my_VarValsXY[CutVarId][1] = VarValY;
    Cut(CutVarId,CutVarVal);
} // Cut

void EventCutter::Cut(unsigned int CutVarId, float CutVarVal){
    _my_CutWasCalled[CutVarId] = true;
    _my_h1is_CutVar_nEvts_BforeAllCuts[CutVarId]->Fill(CutVarVal);
    _my_CutVarVals[CutVarId] = CutVarVal;
    if ( CutVarVal == -999 ) return; // we'll consider this good; no need to check w.r.t. the bounds
    if ( CutVarVal<_my_CutVarMinima[CutVarId] || CutVarVal>_my_CutVarMaxima[CutVarId] ){
        _my_CutVarValIsGood[CutVarId] = false;
        _my_EventIsGood = false;
    }
} // Cut

void EventCutter::FinishedCutting(){
    EventIsGood = _my_EventIsGood; // the user wants the value of this variable at this point; it gets reset at the end of this function
    for( int iCutVar=0; iCutVar<(int)_my_NumCutVars; iCutVar++ ){
        if( !_my_CutWasCalled[iCutVar] ) continue;
        if( _my_EventIsGood ){
            _my_h1is_CutVar_nEvts_AfterAllCuts[iCutVar]->Fill(_my_CutVarVals[iCutVar]);
            if( _my_Cut2DWasCalled[iCutVar] ) _my_h2is_VarX_VarY_nEvts_AfterAllCuts[iCutVar]->Fill(_my_VarValsXY[iCutVar][0],_my_VarValsXY[iCutVar][1]);
        }
        bool AllOtherCutVarValsAreGood = true;
        for( int jCutVar=0; jCutVar<(int)_my_NumCutVars; jCutVar++ ){ // check if all other cut variables are good
            if( iCutVar==jCutVar ) continue; // we don't care about the current cut (iCutVar)
            if( !_my_CutVarValIsGood[jCutVar] ){
                AllOtherCutVarValsAreGood = false;
                break;
            }
        } // jCutVar
        if( AllOtherCutVarValsAreGood ){
            _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut[iCutVar]->Fill(_my_CutVarVals[iCutVar]);
            if( _my_Cut2DWasCalled[iCutVar] ) _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut[iCutVar]->Fill(_my_VarValsXY[iCutVar][0],_my_VarValsXY[iCutVar][1]);
        }
    } // iCutVar
    _my_FinishedCuttingWasCalled = true;    // make sure the user called this; otherwise there will be problems
    _my_EventIsGood = true;                 // reset for next event; true until proven otherwise
    for( unsigned int iCutVar=0; iCutVar<_my_NumCutVars; iCutVar++) _my_CutVarValIsGood[iCutVar] = true;
} // FinishedCutting

void EventCutter::Finish(TFile * MainOutputFile){

    if( !_my_FinishedCuttingWasCalled ) _my_MiscFunctions->PrintError("You didn't call the function \"FinishedCutting\"! This is bad...");

    MainOutputFile->cd();
    _my_CutVarDistHistos = MainOutputFile->mkdir("EventCutter_CutVariableDistributionHistograms" /* Note! LoadRefMultHistogram depends on this directory's name! */);
    _my_CutVarDistHistos->cd(); // NOTE: The constructor looks for the RefMult histogram in this directory!
    // doing the following in separate loops makes finding histograms easier
    for( unsigned int iCutVar=0; iCutVar<_my_NumCutVars; iCutVar++ ){
        if( _my_h1is_CutVar_nEvts_BforeAllCuts[iCutVar] && _my_h1is_CutVar_nEvts_BforeAllCuts[iCutVar]->GetEntries()>0 ){
            _my_h1is_CutVar_nEvts_BforeAllCuts[iCutVar]->SetLineColor(2);
            _my_h1is_CutVar_nEvts_BforeAllCuts[iCutVar]->Write();
        }
        if( _my_h2is_VarX_VarY_nEvts_BforeAllCuts[iCutVar] && _my_h2is_VarX_VarY_nEvts_BforeAllCuts[iCutVar]->GetEntries()>0 ) _my_h2is_VarX_VarY_nEvts_BforeAllCuts[iCutVar]->Write();
    }
    for( int iCutVar=0; iCutVar<(int)_my_NumCutVars; iCutVar++ ){
        if( _my_h1is_CutVar_nEvts_AfterAllCuts[iCutVar] && _my_h1is_CutVar_nEvts_AfterAllCuts[iCutVar]->GetEntries()>0 ) _my_h1is_CutVar_nEvts_AfterAllCuts[iCutVar]->Write();
        if( _my_h2is_VarX_VarY_nEvts_AfterAllCuts[iCutVar] && _my_h2is_VarX_VarY_nEvts_AfterAllCuts[iCutVar]->GetEntries()>0 ) _my_h2is_VarX_VarY_nEvts_AfterAllCuts[iCutVar]->Write();
    }
    for( int iCutVar=0; iCutVar<(int)_my_NumCutVars; iCutVar++ ){
        if( _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut[iCutVar] && _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut[iCutVar]->GetEntries()>0 ) _my_h1is_CutVar_nEvts_AfterAllCutsButThisVarsCut[iCutVar]->Write();
        if( _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut[iCutVar] && _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut[iCutVar]->GetEntries()>0 ) _my_h2is_VarX_VarY_nEvts_AfterAllCutsButThisVarsCut[iCutVar]->Write();
    }

    _my_RefMultDistributionsSeparateCentralities = _my_CutVarDistHistos->mkdir("CutVariableDistributionHistograms_RefMultDistributionsSeparateCentralies");
    _my_RefMultDistributionsSeparateCentralities->cd();
    for( std::vector<TH1I *> RefMultHistogramsForThisCentralityBinList: _my_h1is_RefMult_nEvts_SeparateCentralities ){
        for( TH1I * RefMultHistogram:RefMultHistogramsForThisCentralityBinList ){
            if( RefMultHistogram && RefMultHistogram->GetEntries()>0 ) RefMultHistogram->Write();
        }
    }

    MainOutputFile->cd();
} // Finish
