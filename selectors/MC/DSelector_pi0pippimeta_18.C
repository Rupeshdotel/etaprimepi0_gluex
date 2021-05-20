#include "DSelector_pi0pippimeta_18.h"
#include "DMCThrown.h"
#include "DKinematicData.h"


void DSelector_pi0pippimeta_18::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "pi0pippimeta_18.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit


	//declare outputfile and qfactortree
	fileout = new TFile("qfactortree.root", "Recreate");
	qfactortree = new TTree("qfactortree", "qfactortree");



	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	//std::deque<Particle_t> MyPhi;
	//MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	//dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	//dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//CUT ON SHOWER QUALITY
	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

	//BEAM ENERGY
	//dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));  // Coherent peak for runs in the range 30000-59999

	//KINEMATICS
	//dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	//dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	//dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	//dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	//dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	
	/*
	qfactortree->Branch("mm2", &mm2, "mm2/D");

	qfactortree->Branch("mm2m", &mm2m, "mm2m/D");
	qfactortree->Branch("dt", &dt, "dt/D");
	qfactortree->Branch("be", &be, "be/D");

	qfactortree->Branch("mpi0", &mpi0, "mpi0/D");
	qfactortree->Branch("meta", &meta, "meta/D");
	qfactortree->Branch("metap", &metap, "metap/D");
	qfactortree->Branch("metappi0", &metappi0, "metappi0/D");

	qfactortree->Branch("mpi0m", &mpi0m, "mpi0m/D");
	qfactortree->Branch("metam", &metam, "metam/D");
	qfactortree->Branch("metapm", &metapm, "metapm/D");
	qfactortree->Branch("metappi0m", &metappi0m, "metappi0m/D");

	qfactortree->Branch("mpi013", &mpi013, "mpi013/D");
	qfactortree->Branch("mpi014", &mpi014, "mpi014/D");
	qfactortree->Branch("mpi023", &mpi023, "mpi023/D");
	qfactortree->Branch("mpi024", &mpi024, "mpi024/D");

	qfactortree->Branch("mpi013m", &mpi013m, "mpi013m/D");
	qfactortree->Branch("mpi014m", &mpi014m, "mpi014m/D");
	qfactortree->Branch("mpi023m", &mpi023m, "mpi023m/D");
	qfactortree->Branch("mpi024m", &mpi024m, "mpi024m/D");

	qfactortree->Branch("mpippimpi0", &mpippimpi0, "mpippimpi0/D");
	qfactortree->Branch("mproeta", &mproeta, "mproeta/D");
	qfactortree->Branch("mpippimpi0m", &mpippimpi0m, "mpippimpi0m/D");

	qfactortree->Branch("mpi0p", &mpi0p, "mpi0p/D");
	qfactortree->Branch("mpi0pm", &mpi0pm, "mpi0pm/D");

	qfactortree->Branch("mant", &mant, "mant/D");
	qfactortree->Branch("mantm", &mantm, "mantm/D");

	qfactortree->Branch("photon1_sq", &photon1_sq, "photon1_sq/D");
	qfactortree->Branch("photon2_sq", &photon2_sq, "photon2_sq/D");
	qfactortree->Branch("photon3_sq", &photon3_sq, "photon3_sq/D");
	qfactortree->Branch("photon4_sq", &photon4_sq, "photon4_sq/D");

	qfactortree->Branch("px_pr", &px_pr, "px_pr/F");
	qfactortree->Branch("px_etapr", &px_etapr, "px_etapr/F");
	qfactortree->Branch("px_pi0", &px_pi0, "px_pi0/F");

	qfactortree->Branch("py_pr", &py_pr, "py_pr/F");
	qfactortree->Branch("py_etapr", &py_etapr, "py_etapr/F");
	qfactortree->Branch("py_pi0", &py_pi0, "py_pi0/F");

	qfactortree->Branch("pz_pr", &pz_pr, "pz_pr/F");
	qfactortree->Branch("pz_etapr", &pz_etapr, "pz_etapr/F");
	qfactortree->Branch("pz_pi0", &pz_pi0, "pz_pi0/F");

	qfactortree->Branch("e_pr", &e_pr, "e_pr/F");
	qfactortree->Branch("e_etapr", &e_etapr, "e_etapr/F");
	qfactortree->Branch("e_pi0", &e_pi0, "e_pi0/F");


	
	qfactortree->Branch("px_beam", &px_beam, "px_beam/F");
	qfactortree->Branch("py_beam", &px_beam, "py_beam/F");
	qfactortree->Branch("pz_beam", &pz_beam, "pz_beam/F");
	qfactortree->Branch("e_beam", &e_beam, "e_beam/F");

	qfactortree->Branch("cos_t", &cos_t, "cos_t/D");
	qfactortree->Branch("costheta_X_cm", &costheta_X_cm, "costheta_X_cm/D");
	qfactortree->Branch("phi_gj", &phi_gj, "phi_gj/D");
	qfactortree->Branch("pol", &pol, "pol/I");

	qfactortree->Branch("event_num", &event_num, "event_num/I");
	qfactortree->Branch("run_num", &run_num, "run_num/I");
		*/

	qfactortree->Branch("mm2m", &mm2m, "mm2m/D");
	qfactortree->Branch("dt", &dt, "dt/D");
	qfactortree->Branch("be", &be, "be/D");
	qfactortree->Branch("mpi013", &mpi013, "mpi013/D");
	qfactortree->Branch("mpi014", &mpi014, "mpi014/D");
	qfactortree->Branch("mpi023", &mpi023, "mpi023/D");
	qfactortree->Branch("mpi024", &mpi024, "mpi024/D");

	
	qfactortree->Branch("metap", &metap, "metap/D");
	qfactortree->Branch("metappi0", &metappi0, "metappi0/D");
	qfactortree->Branch("mpippimpi0", &mpippimpi0, "mpippimpi0/D");
	qfactortree->Branch("mproeta", &mproeta, "mproeta/D");
	qfactortree->Branch("mpi0p", &mpi0p, "mpi0p/D");
	qfactortree->Branch("mant", &mant, "mant/D");
	qfactortree->Branch("photon1_sq", &photon1_sq, "photon1_sq/D");
	qfactortree->Branch("photon2_sq", &photon2_sq, "photon2_sq/D");
	qfactortree->Branch("photon3_sq", &photon3_sq, "photon3_sq/D");
	qfactortree->Branch("photon4_sq", &photon4_sq, "photon4_sq/D");
	qfactortree->Branch("cos_t", &cos_t, "cos_t/D");
	qfactortree->Branch("costheta_X_cm", &costheta_X_cm, "costheta_X_cm/D");
	qfactortree->Branch("phi_gj", &phi_gj, "phi_gj/D");

	qfactortree->Branch("photon1id_thrown", &photon1id_thrown, "photon1id_thrown/I");
	qfactortree->Branch("photon1id_parentid", &photon1id_parentid, "photon1id_parentid/I");

	qfactortree->Branch("photon2id_thrown", &photon2id_thrown, "photon2id_thrown/I");
	qfactortree->Branch("photon3id_thrown", &photon3id_thrown, "photon3id_thrown/I");
	qfactortree->Branch("photon4id_thrown", &photon4id_thrown, "photon4id_thrown/I");

	qfactortree->Branch("protonid_thrown", &protonid_thrown, "protonid_thrown/I");
	qfactortree->Branch("pipid_thrown", &pipid_thrown, "pipid_thrown/I");
	qfactortree->Branch("pimid_thrown", &pimid_thrown, "pimid_thrown/I");

	qfactortree->Branch("beam_thrown_id", &beam_thrown_id, "beam_thrown_id/I");
	qfactortree->Branch("beam_recon_id", &beam_recon_id, "beam_recon_id/I");


	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	//dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_pi0pippimeta_18::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	if(Get_EventNumber()%100000==0) {cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;}
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dAnalysisUtilities.Get_PolarizationAngle(locRunNumber,  dPolarizationAngle); 
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	//Reset_Actions_NewEvent();
	//dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

		/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/

	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);
		beam_thrown_id = dThrownBeam->Get_BeamID();

		//Do stuff with the wrapper here ...
	}



	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously


		double dMinKinFitCL = 1e-03; //throw out bad results within 0.1 percent of resolution
		double locKinFitConLev = dComboWrapper->Get_ConfidenceLevel_KinFit();
		if(locKinFitConLev < dMinKinFitCL){
			
		  dComboWrapper->Set_IsComboCut(true);
		  continue;
		  }

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();

		//Step 2
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();

		//Step 3
		Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();
		Int_t locPhoton4NeutralID = dPhoton4Wrapper->Get_NeutralID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		double BE = locBeamP4.E();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		//Step 2
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		//Step 3
		TLorentzVector locDecayingEtaP4 = dDecayingEtaWrapper->Get_P4();
		TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4();
		TLorentzVector locPhoton4P4 = dPhoton4Wrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locPhoton1P4_Measured = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4_Measured = dPhoton2Wrapper->Get_P4_Measured();
		//Step 2
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		//Step 3
		TLorentzVector locPhoton3P4_Measured = dPhoton3Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton4P4_Measured = dPhoton4Wrapper->Get_P4_Measured();

		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		//Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);

		// calculate accidental subtraction weight based on time difference 
		double locWeight = 0.; // weight to accidentally subtracted histgorams
		bool locAccid = false; // flag to fill separate prompt and accidental histograms for later subtraction

		if(fabs(locDeltaT_RF) < 0.5*4.008) { // prompt signal recieves a weight of 1
		  locWeight = 1.;
		  locAccid = false;
		}
                else { // accidentals recieve a weight of 1/# RF bunches included in TTree (8 in this case)
		  locWeight = -1./8.;
		  locAccid = true;


		  }
		//Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		// Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		// Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		// Int_t locNumOutOfTimeBunchesInTree = 4; //YOU need to specify this number
			//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 

		// Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		// Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		// Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		// Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		// Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		// if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		// 	dComboWrapper->Set_IsComboCut(true); 
		// 	continue; 
		// } 

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		//TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		//locMissingP4_Measured -= locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured + locPiMinusP4_Measured + locPiPlusP4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured;

		TLorentzVector totalP4_Measured =  locBeamP4_Measured + dTargetP4;		
		TLorentzVector locMissingP4_Measured = totalP4_Measured - 
		 (locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured
		  + locPhoton1P4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured);

		TLorentzVector totalP4 =  locBeamP4 + dTargetP4;		
		TLorentzVector locMissingP4 = totalP4 - 
		 (locPiPlusP4 + locPiMinusP4 + locProtonP4
		  + locPhoton1P4 + locPhoton2P4 + locPhoton3P4 + locPhoton4P4);



		TLorentzVector locPhoton12 = locPhoton1P4 + locPhoton2P4;
		TLorentzVector locpi0 = locPhoton1P4 + locPhoton2P4;
		TLorentzVector locPhoton34 = locPhoton3P4 + locPhoton4P4;
		TLorentzVector loceta = locPhoton3P4 + locPhoton4P4;

		TLorentzVector locpeta = locProtonP4 + locPhoton3P4 + locPhoton4P4;

		TLorentzVector locPhoton12m = locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector locPhoton34m = locPhoton3P4_Measured + locPhoton4P4_Measured;

		TLorentzVector locetaprime = locPiPlusP4 + locPiMinusP4	+ locPhoton3P4 + locPhoton4P4;
		TLorentzVector locetaprimem = locPiPlusP4_Measured + locPiMinusP4_Measured	+ locPhoton3P4_Measured + locPhoton4P4_Measured;

		TLorentzVector locetaprimepi0 = locPiPlusP4 + locPiMinusP4	+ locPhoton3P4 + locPhoton4P4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector locetaprimepi0m = locPiPlusP4_Measured + locPiMinusP4_Measured + 
									locPhoton3P4_Measured + locPhoton4P4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;


		TLorentzVector locPhoton13 = locPhoton1P4 + locPhoton3P4;
		TLorentzVector locPhoton14 = locPhoton1P4 + locPhoton4P4;
		TLorentzVector locPhoton23 = locPhoton2P4 + locPhoton3P4;
		TLorentzVector locPhoton24 = locPhoton2P4 + locPhoton4P4;

		TLorentzVector locPhoton13m = locPhoton1P4_Measured + locPhoton3P4_Measured;
		TLorentzVector locPhoton14m = locPhoton1P4_Measured + locPhoton4P4_Measured;
		TLorentzVector locPhoton23m = locPhoton2P4_Measured + locPhoton3P4_Measured;
		TLorentzVector locPhoton24m = locPhoton2P4_Measured + locPhoton4P4_Measured;

		TLorentzVector locpippimpi0 = locPiPlusP4 + locPiMinusP4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector locpippimpi0m = locPiPlusP4_Measured + locPiMinusP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;


		TLorentzVector locpi0p = locProtonP4 + locPhoton1P4 + locPhoton2P4;
		TLorentzVector locpi0pm = locProtonP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

		TLorentzVector loct = locBeamP4 - locetaprimepi0;
		TLorentzVector loctm = locBeamP4_Measured - locetaprimepi0m;

		


		//costheta in GJ frame 
		//redefine centre of mass 4-vectors to 4-vectors in rest frame of etaprimepi0 system
   		TLorentzVector BeamP4_GJ = locBeamP4;
   		TLorentzVector EtaPrimeP4_GJ = locetaprime;
   		TLorentzVector Pi0P4_GJ = locpi0;
   		TLorentzVector EtaPrimePi0P4_GJ = locetaprimepi0;

		TVector3 boostGJ;

  		TVector3 z_GJ;
   		TVector3 z_hat_GJ;

   		TVector3 y_GJ;
   		TVector3 y_hat_GJ;

   		TVector3 x_GJ;
   		TVector3 x_hat_GJ;

   		TVector3 vetaprime;


		// boost to GJ Frame
   		boostGJ = -(locetaprimepi0.Vect())*(1.0/locetaprimepi0.E()); //calculate beta for etaprimepi0 system

		//boost in rest frame of etaprimepi0 system
  		BeamP4_GJ.Boost(boostGJ);
   		EtaPrimeP4_GJ.Boost(boostGJ);
   		Pi0P4_GJ.Boost(boostGJ);
		EtaPrimePi0P4_GJ.Boost(boostGJ);

		z_GJ.SetXYZ(BeamP4_GJ.X(), BeamP4_GJ.Y(), BeamP4_GJ.Z());//z GJ
		z_hat_GJ = z_GJ.Unit();  //beam direction (z-direction) 
		
		y_GJ = locBeamP4.Vect().Cross(locetaprimepi0.Vect());  //y-direction for rest frame of etaprimepi0 system
   		y_hat_GJ = y_GJ.Unit();    

		x_hat_GJ = y_hat_GJ.Cross(z_hat_GJ);//x hat GJ 

		vetaprime.SetXYZ(EtaPrimeP4_GJ.Vect()*x_hat_GJ, EtaPrimeP4_GJ.Vect()*y_hat_GJ, EtaPrimeP4_GJ.Vect()*z_hat_GJ);

		// costheta of resonace in CM frame
		// get cm vector
		TLorentzVector cm_vec = locBeamP4 + dTargetP4;
		
		// define etaprimepi0 4vector  in cm frame
		TLorentzVector locetaprimepi0_cm = locetaprimepi0;

		// get the beta through boostvector and boost in cm frame
		locetaprimepi0_cm.Boost(-cm_vec.BoostVector());

		//get the costhetaXcm
		costheta_X_cm = locetaprimepi0_cm.CosTheta();







		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		//dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		//if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
		//	continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		//if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		//{
			//dHist_BeamEnergy->Fill(locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
			//dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor); // Alternate version with accidental subtraction

			//locUsedSoFar_BeamEnergy.insert(locBeamID);
		//}

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();
		//define vairables for qfactortree

		 /*
		 mm2 = locMissingP4.M2();

		 mpi0 = locPhoton12.M();
		 meta = locPhoton34.M();
		
		 mpippimpi0m = locpippimpi0m.M();
		 mpi0m= locPhoton12m.M();
		 mpi013m = locPhoton13m.M();
		 mpi014m = locPhoton14m.M();
		 mpi023m = locPhoton23m.M();
		 mpi024m = locPhoton24m.M();
		 metam = locPhoton34m.M();
		 metapm = locetaprimem.M();
		 metappi0m = locetaprimepi0m.M();
		 mpi0pm = locpi0pm.M();
		 
		 mantm = (-1)*loctm.M2();
		 
		

		//get costhetat in GJ
   		
		pol = dPolarizationAngle; 

		 px_pr = locProtonP4.Px(); 
		 px_etapr = locetaprime.Px();
		 px_pi0 = locPhoton12.Px();
		


		 py_pr = locProtonP4.Py();
		 py_etapr = locetaprime.Py();
		 py_pi0 = locPhoton12.Py();
		

		 pz_pr = locProtonP4.Pz();
		 pz_etapr = locetaprime.Pz();
		 pz_pi0 = locPhoton12.Pz();
		


		 e_pr = locProtonP4.E();
		 e_etapr = locetaprime.E();
		 e_pi0 = locPhoton12.E();
		

		 px_beam = locBeamP4.Px();
		 py_beam = locBeamP4.Py();
		 pz_beam = locBeamP4.Pz();
		 e_beam = locBeamP4.E();


		chisq = dComboWrapper->Get_ChiSq_KinFit();  //cout <<  chisq << "chisq" << endl;
		ndf = dComboWrapper->Get_NDF_KinFit();     //cout <<  ndf << "ndf" << endl;
		chisq_ndf = chisq/ndf;                     //cout << chisq_ndf << "chisq/ndf" << endl;
		event_num = Get_EventNumber();
		run_num = locRunNumber;
		*/
		
		mm2m = locMissingP4_Measured.M2();
		 
		 mpi013 = locPhoton13.M();
		 mpi014 = locPhoton14.M();
		 mpi023 = locPhoton23.M();
		 mpi024 = locPhoton24.M();
		 
		 metap = locetaprime.M();
		 metappi0 = locetaprimepi0.M();
		 mpippimpi0 = locpippimpi0.M();
		 mpi0p = locpi0p.M();
		 dt = locDeltaT_RF;
		 be = BE;
		 mant = (-1)*loct.M2();
		 mproeta = locpeta.M();
		 photon1_sq = dPhoton1Wrapper->Get_Shower_Quality();
		 photon2_sq = dPhoton2Wrapper->Get_Shower_Quality();
		 photon3_sq = dPhoton3Wrapper->Get_Shower_Quality();
		 photon4_sq = dPhoton4Wrapper->Get_Shower_Quality();
		 cos_t = vetaprime.CosTheta();
		 Double_t pi = 3.14159;
		 phi_gj = vetaprime.Phi() * (180/pi); 

		 //check the match of reconstructed with the thrown info

		 
		 //dThrownWrapper->Set_ArrayIndex(dPhoton1Wrapper->Get_ThrownIndex());
		  photon1id_thrown = dPhoton1Wrapper->Get_ThrownIndex();
		  //dThrownWrapper->Set_ArrayIndex(dPhoton1Wrapper->Get_ThrownIndex());
		  //photon1id_parentid = dPhoton1Wrapper->Get_ParentIndex();
		  
		  
		  photon2id_thrown = dPhoton2Wrapper->Get_ThrownIndex();
		  photon3id_thrown = dPhoton3Wrapper->Get_ThrownIndex();
		  photon4id_thrown = dPhoton4Wrapper->Get_ThrownIndex();
		  protonid_thrown = dProtonWrapper->Get_ThrownIndex();
		  pipid_thrown = dPiPlusWrapper->Get_ThrownIndex();
		  pimid_thrown = dPiMinusWrapper->Get_ThrownIndex();
		  beam_thrown_id = dThrownBeam->Get_BeamID();
		  beam_recon_id = dComboBeamWrapper->Get_BeamID();

		 

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton1NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton2NeutralID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton3NeutralID);
		locUsedThisCombo_MissingMass[Gamma].insert(locPhoton4NeutralID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			//dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
			//dHist_MissingMassSquared->Fill(locMissingMassSquared,locHistAccidWeightFactor); // Alternate version with accidental subtraction


			Bool_t coherentbeam = (be > 8.2) && (be < 8.8);
			Bool_t MissingMassSquaredcut = (mm2m > -0.02) && (mm2m < 0.02);


			Bool_t etaprimemasswindow = (metap > 0.85) && (metap < 1.05);
			Bool_t twindow = (mant > 0.1) && (mant < 0.7);

			Bool_t pi013 = ((mpi013 > 0.10) && (mpi013 < 0.175)); 
			Bool_t pi024 = ((mpi024 > 0.10) && (mpi024 < 0.175));
			Bool_t pi014 = ((mpi014 > 0.10) && (mpi014 < 0.175));
			Bool_t pi023 = ((mpi023 > 0.10) && (mpi023 < 0.175));
			Bool_t omega = (mpippimpi0 > 0.73) && (mpippimpi0 < 0.83);
			Bool_t pi0p_deltap = mpi0p < 1.4;

			/*
			if(MissingMassSquaredcut && coherentbeam && !locAccid &&  etaprimemasswindow && twindow && //selection cuts

			!((pi013 && pi024) || (pi014 && pi023)) && 									  //veto cuts
			!(pi0p_deltap) &&
			!(omega) )
			*/
			
			if(MissingMassSquaredcut && coherentbeam && !locAccid && twindow && etaprimemasswindow)
			{
				

				qfactortree->Fill();
				
			}

			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/

		//FILL FLAT TREE
		//Fill_FlatTree(); //for the active combo
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

		/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/



	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_pi0pippimeta_18::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	fileout->cd(); 
	qfactortree->Write();

	fileout->Close();

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
