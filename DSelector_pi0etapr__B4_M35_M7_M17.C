#include "DSelector_pi0etapr__B4_M35_M7_M17.h"
#include <iostream>
#include <fstream>
using namespace std;
//ofstream myfile ("pi0qfactor.txt");
//ofstream myfile1 ("etaprimeqfactor.data");
//ofstream myfile2 ("t_etaprimepi0.data");
//ofstream myfile3 ("etaqfactor.data");
//ofstream myfile4 ("etaprimepi0qfactor.data");




void DSelector_pi0etapr__B4_M35_M7_M17::Init(TTree *locTree) // ::(called scope) is used to define a func outside a class in cpp
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
  	dOutputFileName = "pi0etapr__B4_M35_M7_M17.root"; //"" for none
	//dOutputTreeFileName = "clonetest.root"; //"" for none
	//dOutputTreeFileName = "custom.root"; //"" for none
	//dOutputTreeName = "customtree"; //"" for none
	//dFlatTreeFileName = "flat_pi0etapr__B4_M35_M7_M17.root"; //output flat tree (one combo per tree entry), "" for none
	//dFlatTreeName = ""; //if blank, default name will be chosen
	

  
	
	
      

	



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

	std::deque<Particle_t> MyEtaPrime;

  	MyEtaPrime.push_back(PiMinus); MyEtaPrime.push_back(PiPlus); MyEtaPrime.push_back(Eta);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//cout << "dCombowrapper " << dComboWrapper << endl;

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.0, 11.6));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	//dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyEtaPrime, 1000, 0.9, 2.4, "CutActionEffect" );

	// Unused energy cut
	//dAnalysisActions.push_back(new DCutAction_Energy_UnusedShowers(dComboWrapper, 0.0));

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	//dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.1, 0.1);

	dHist_MissingEnergy = new TH1I("MissingEnergy", ";Missing Energy (GeV/c^{2})", 600, -2.0, 2.0);


	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	dHist_BeamDeltaT = new TH1F("BeamDeltaT", "; t_{Tagger} - t_{RF} (ns)", 400, -10.0, 10.0);

	//1-D photon pairs

	dHist_allpairs_1D = new TH1I("allpairs_1D", ";M(#gamma#gamma)", 200, 0.0, 1.2);

	dHist_Photons12 = new TH1I("Photons12", ";Inv_MPhotons12 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons13 = new TH1I("Photons13", ";Inv_MPhotons13 (GeV/c^{2})", 200, 0.0, 1.0);

	dHist_Photons1v2_SC = new TH2F("Photons1v2_SC", ";E#gamma_{1} (GeV/c^{2});E#gamma_{2} (GeV/c^{2})", 200, 0.0, 2.0, 200, 0.0, 2.0);
	dHist_Photons1v2 = new TH2F("Photons1v2", ";E#gamma_{1} (GeV/c^{2});E#gamma_{2} (GeV/c^{2})", 200, 0.0, 2.0, 200, 0.0, 2.0);

	dHist_Photons3v4_SC = new TH2F("Photons3v4_SC", ";E#gamma_{3} (GeV/c^{2});E#gamma_{4} (GeV/c^{2})", 200, 0.0, 2.0, 200, 0.0, 2.0);
	dHist_Photons3v4 = new TH2F("Photons3v4", ";E#gamma_{3} (GeV/c^{2});E#gamma_{4} (GeV/c^{2})", 200, 0.0, 2.0, 200, 0.0, 2.0);

	


	dHist_Photons12_exclusiveboth = new TH1F("Photons12_exclusiveboth", ";M(#gamma_{1}#gamma_{2}) (GeV/c^{2}); counts / 0.005 GeV ", 200, 0.0, 1.0);
	dHist_Photons34_exclusiveboth = new TH1F("Photons34_exclusiveboth", ";M(#gamma_{3}#gamma_{4}) (GeV/c^{2}); counts / 0.005 GeV", 200, 0.0, 1.0);

	dHist_Photons12_exclusiveboth_accid = new TH1F("Photons12_exclusiveboth_accid", ";M(#gamma_{1}#gamma_{2}) (GeV/c^{2}); counts / 0.005 GeV ", 200, 0.0, 1.0);
	dHist_Photons34_exclusiveboth_accid = new TH1F("Photons34_exclusiveboth_accid", ";M(#gamma_{3}#gamma_{4}) (GeV/c^{2}); counts / 0.005 GeV", 200, 0.0, 1.0);
	dHist_pippimeta_accid = new TH1F("pippimeta_accid", ";M(#pi^{+}#pi^{-}#eta (GeV/c^{2});", 200, 0.6, 1.6);

	dHist_Photons12_exclusiveboth_prompt = new TH1F("Photons12_exclusiveboth_prompt", ";M(#gamma_{1}#gamma_{2}) (GeV/c^{2}); counts / 0.005 GeV ", 200, 0.0, 1.0);
	dHist_Photons34_exclusiveboth_prompt = new TH1F("Photons34_exclusiveboth_prompt", ";M(#gamma_{3}#gamma_{4}) (GeV/c^{2}); counts / 0.005 GeV", 200, 0.0, 1.0);
	dHist_pippimeta_prompt = new TH1F("pippimeta_prompt", ";M(#pi^{+}#pi^{-}#eta (GeV/c^{2});", 200, 0.6, 1.6);

	
	dHist_Photons12_exclusiveboth_accid->Sumw2();
	dHist_Photons34_exclusiveboth_accid->Sumw2();
	dHist_pippimeta_accid->Sumw2();
	
	dHist_Photons12_exclusiveboth_prompt->Sumw2();
	dHist_Photons34_exclusiveboth_prompt->Sumw2();
	dHist_pippimeta_prompt->Sumw2();


	dHist_Photons12_exclusiveboth_SC = new TH1I("Photons12_exclusiveboth_SC", ";M(#gamma_{1}#gamma_{2}) (GeV/c^{2}); counts / 0.005 GeV ", 200, 0.0, 1.0);
	dHist_Photons34_exclusiveboth_SC = new TH1I("Photons34_exclusiveboth_SC", ";M(#gamma_{3}#gamma_{4}) (GeV/c^{2}); counts / 0.005 GeV", 200, 0.0, 1.0);

	dHist_Photons12_belowthreshold = new TH1I("Photons12_belowthreshold", ";Inv_MPhotons12 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons34_belowthreshold = new TH1I("Photons34_belowthreshold", ";Inv_MPhotons34 (GeV/c^{2})", 200, 0.0, 1.0);

	dHist_Photons12_inclusiveboth = new TH1I("Photons12_inclusiveboth", ";Inv_MPhotons12 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons34_inclusiveboth = new TH1I("Photons34_inclusiveboth", ";Inv_MPhotons34 (GeV/c^{2})", 200, 0.0, 1.0);

	dHist_Photons14 = new TH1I("Photons14", ";Inv_MPhotons14 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons23 = new TH1I("Photons23", ";Inv_MPhotons23 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons24 = new TH1I("Photons24", ";Inv_MPhotons24 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons34 = new TH1I("Photons34", ";Inv_MPhotons34 (GeV/c^{2})", 200, 0.0, 1.0);

	

	// reconstruction of etaprime and other possible combinations

	

	dHist_pippimeta_SC = new TH1I("pippimeta_SC", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2}); counts / 0.0075 GeV", 200, 0.6, 1.6);

	dHist_pippimeta_pi0selectiononly = new TH1I("pippimeta_pi0selectiononly", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);

	dHist_pippim = new TH1I("pippim", ";M(#pi^{+}#pi^{-}) (GeV/c^{2})", 200, 0.0, 1.5);

	dHist_pimpi0 = new TH1I("pimpi0", ";M(#pi^{-}#pi^{0}) (GeV/c^{2})", 200, 0.0, 1.5);

	dHist_pippimetapi0 = new TH1I("pippimetapi0", ";M(#pi^{+}#pi^{-}#eta#pi^{0}) (GeV/c^{2})", 200, 0.8, 3.0);

	dHist_pippimeta_belowthreshold = new TH1I("pippimeta_belowthreshold", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);
	dHist_pippimeta_abovethreshold = new TH1I("pippimeta_abovethreshold", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);

	
	dHist_pippimeta_aboveomega = new TH1I("pippimeta_cutpi0pippim", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);

	dHist_pippimeta = new TH1I("pippimeta", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2});", 200, 0.6, 1.6);
	dHist_pippimeta_cutbr = new TH1I("pippimeta_cutbr", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.6, 1.6);
	dHist_pippimeta_cutomega = new TH1I("pippimeta_cutomega", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.6, 1.6);
	dHist_pippimeta_aboveomega = new TH1I("pippimeta_aboveomega", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.6, 1.6);

	dHist_etaprimepi0 = new TH1I("etaprimepi0", ";M(#eta^{'}#pi^{0})(GeV/c^{2})   ", 100, 1.0, 4.0);
	dHist_etaprimepi0_cutbr = new TH1I("etaprimepi0_cutbr", ";M(#eta^{'}#pi^{0})(GeV/c^{2})  ", 100, 1.0, 4.0);
	dHist_etaprimepi0_cutomega = new TH1I("etaprimepi0_cutomega", ";M(#eta^{'}#pi^{0})(GeV/c^{2})  ", 100, 1.0, 4.0);
	dHist_etaprimepi0_aboveomega = new TH1I("etaprimepi0_aboveomega", ";M(#eta^{'}#pi^{0})(GeV/c^{2})  ", 100, 1.0, 4.0);

	dHist_InvMetapvcostheta = new TH2F("InvMetapvcostheta ", ";M(#eta^{'}#pi^{0})(GeV/c^{2}) ;cos#theta_{GJ} ", 100, 1.0, 4.0, 100, -1.0, 1.0);
	dHist_InvMetapvcostheta_cutbr = new TH2F("InvMetapvcostheta_cutbr ", ";M(#eta^{'}#pi^{0})(GeV/c^{2}) ;cos#theta_{GJ} ", 100, 1.0, 4.0, 100, -1.0, 1.0);
	dHist_InvMetapvcostheta_cutomega = new TH2F("InvMetapvcostheta_cutomega ", ";M(#eta^{'}#pi^{0})(GeV/c^{2}) ;cos#theta_{GJ} ", 100, 1.0, 4.0, 100, -1.0, 1.0);
	dHist_InvMetapvcostheta_aboveomega = new TH2F("InvMetapvcostheta_aboveomega ", ";M(#eta^{'}#pi^{0})(GeV/c^{2}) ;cos#theta_{GJ} ", 100, 1.0, 4.0, 100, -1.0, 1.0);
	
	
	

	dHist_pippimeta_cutomega_SC = new TH1I("pippimeta_cutomega_SC", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.6, 1.6);

	dHist_pippimeta_omega = new TH1I("pippimeta_omega", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.6, 1.6);
	//dHist_pippimeta_cut_pipp_pi0p = new TH1I("pippimeta_cut_pipp_pi0p", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);

	dHist_pippimeta_beloweta = new TH1I("pippimeta_beloweta", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);
	dHist_pippimeta_eta = new TH1I("pippimeta_eta", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);
	dHist_pippimeta_aboveeta = new TH1I("pippimeta_aboveeta", ";M(#pi^{+}#pi^{-}#eta) (GeV/c^{2})", 200, 0.5, 2.0);


	dHist_pippimpi0 = new TH1I("pippimpi0", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.5);

	dHist_pippimpi0_SC = new TH1I("pippimpi0_SC", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.5);


	dHist_pippimpi0_M = new TH1I("pippimpi0_M", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.5);

	dHist_pippimpi0_pi0selectiononly = new TH1I("pippimpi0_pi0selectiononly", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.5);

	dHist_pippimpi0_belowthreshold = new TH1I("pippimpi0_belowthreshold", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.0);
	dHist_pippimp = new TH1I("pippimp", ";M(#pi^{+}#pi^{-}p) (GeV/c^{2})", 200, 1.2, 4.0);
	dHist_pippi0p = new TH1I("pippi0p", ";M(#pi^{+}#pi^{0}p) (GeV/c^{2})", 200, 1.2, 4.0);
	dHist_etapi0 = new TH1I("etapi0", ";M(#eta#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.0);

	dHist_etapi0_belowthreshold = new TH1I("etapi0_belowthreshold", ";M(#eta#pi^{0}) (GeV/c^{2})", 200, 0.0, 2.0);

	dHist_pi0p = new TH1I("pi0p", ";M(#pi^{0}p) (GeV/c^{2})", 200, 0.8, 3.5);
	dHist_pipp = new TH1I("pipp", ";M(#pi^{+}p) (GeV/c^{2})", 200, 0.8, 3.5);
	dHist_pimp = new TH1I("pimp", ";M(#pi^{-}p) (GeV/c^{2})", 200, 0.8, 3.5);

	dHist_etap = new TH1I("etap", ";M(#etap) (GeV/c^{2})", 200, 1.0, 4.5);
	dHist_etap_reverse = new TH1I("etap_reverse", ";M(#etap_reverse) (GeV/c^{2})", 200, 0.0, 2.5);

	dHist_etap_reverse_M = new TH1I("etap_reverse_M", ";M(#etap_reverse) (GeV/c^{2})", 200, 0.0, 2.5);

	
	dHist_eta_cleanetap_alleta = new TH1I("eta_cleanetap_alleta", ";M(#eta_cleanetaprime) (GeV/c^{2})", 200, 0.0, 1.0);



	dHist_pippimetapi0vpippimeta = new TH2F("pippimetapi0vpippimeta", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}#eta)(GeV/c^{2})", 100, 0.0, 4.5, 100, 0.8, 1.6);

	

	dHist_pippimetapi0vpippimeta_aboveomega = new TH2F("pippimetapi0vpippimeta_cutpippimpi0", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}#eta)(GeV/c^{2})", 100, 0.0, 4.5, 100, 0.8, 1.6);

	

	dHist_pippimetapi0vpippimeta_cut_pipp_pi0p = new TH2F("pippimetapi0vpippimeta_cut_pipp_pi0p", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}#eta)(GeV/c^{2})", 100, 0.0, 4.5, 100, 0.8, 1.6);
	
	dHist_pippimetapi0vpippimpi0 = new TH2F("pippimetapi0vpippimpi0", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}#pi^{0})(GeV/c^{2})", 100, 0.0, 4.5, 100, 0.0, 1.5);
	dHist_pippimetapi0vpippimp = new TH2F("pippimetapi0vpippimp", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}p)(GeV/c^{2})", 100, 0.0, 4.5, 100, 1.2, 4.0);
	dHist_pippimetapi0vpippi0p = new TH2F("pippimetapi0vpippi0p", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{0}p)(GeV/c^{2})", 100, 0.0, 4.5, 100, 1.2, 4.0);
	dHist_pippimetapi0vetapi0 = new TH2F("pippimetapi0vetapi0", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#eta#pi^{0})(GeV/c^{2})", 100, 0.0, 4.5, 100, 0.0, 1.8);
	dHist_pippimetapi0vpipp = new TH2F("pippimetapi0vpipp", ";M(#pi^{+}#pi^{-}#eta#pi^{0})(GeV/c^{2});M(#pi^{+}p)(GeV/c^{2})", 100, 0.0, 4.5, 100, 0.0, 3.5);

	dHist_pippimetavpippimpi0 = new TH2F("pippimetavpippimpi0", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2}) counts / 0.015 GeV;M(#pi^{+}#pi^{-}#eta)(GeV/c^{2}) counts / 0.007 GeV", 100, 0.5, 2.0, 100, 0.8, 1.5);
	dHist_pippimetavpippimpi0_cutomega = new TH2F("pippimetavpippimpi0_cutomega", ";M(#pi^{+}#pi^{-}#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}#eta)(GeV/c^{2})", 100, 0.5, 2.0, 100, 0.8, 1.5);
	dHist_pippimetavpippimpi0_omega = new TH2F("pippimetavpippimpi0_omega", ";M(#pi^{+}#pi^{-}#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-}#eta)(GeV/c^{2})", 100, 0.5, 2.0, 100, 0.8, 1.5);


	dHist_pippimetavpi0p = new TH2F("pippimetavpi0p", "; M(#pi^{0}p)(GeV/c^{2});M(#pi^{+}#pi^{-}#eta) (GeV/c^{2}) ", 100, 1.0, 3.0, 100, 0.8, 1.6);
	dHist_pippi0etavpimp = new TH2F("pippi0etavpimp", " ;M(#pi^{-}p)(GeV/c^{2});M(#pi^{+}#pi^{0}#eta) (GeV/c^{2}) ", 100, 1.0, 3.0, 100, 0.8, 2.5);
	dHist_pi0pimetavpipp = new TH2F("pi0pimetavpipp", " ;M(#pi^{+}p)(GeV/c^{2});M(#pi^{0}#pi^{-}#eta) (GeV/c^{2}) ", 100, 1.0, 3.0, 100, 0.8, 2.5);
	dHist_pi0etapvpippim = new TH2F("pi0etapvpippim", " ;M(#pi^{+}#pi^{-})(GeV/c^{2});M(#pi^{0}#etap) (GeV/c^{2}) ", 100, 0.25, 1.0, 100, 1.5, 4.5);

	dHist_pippimetavpippimpi0_SC = new TH2F("pippimetavpippimpi0_SC", ";M(#pi^{+}#pi^{-}#pi^{0}) (GeV/c^{2}) counts / 0.015 GeV;M(#pi^{+}#pi^{-}#eta)(GeV/c^{2}) counts / 0.007 GeV", 100, 0.5, 2.0, 100, 0.8, 1.5);
	

	dHist_pipetavspippi0 = new TH2F("pipetavspippi0", ";M(#pi^{+}#eta)(GeV/c^{2});M(#pi^{+}#pi^{0})(GeV/c^{2})", 100, 0.0, 2.0, 100, 0.0, 2.0);
	dHist_etapi0vspippim = new TH2F("etapi0vspippim", ";M(#eta#pi^{0})(GeV/c^{2});M(#pi^{+}#pi^{-})(GeV/c^{2})", 100, 0.0, 3.0, 100, 0.0, 2.0);

	
	dHist_etaprimepi0_cutomega_SC = new TH1I("etaprimepi0_cutomega_SC", ";M(#eta^{'}#pi^{0})(GeV/c^{2});  counts / 0.035 GeV", 100, 0.0, 4.0);

	//dHist_etaprimepi0_cutomega = new TH1I("etaprimepi0_cutomega", ";M(#eta^{'}#pi^{0}) (GeV/c^{2})", 100, 0.0, 4.0);


	
	

	dHist_InvMetapvcostheta_SC = new TH2F("InvMetapvcostheta_SC ", ";M(#eta^{'}#pi^{0})(GeV/c^{2})  counts / 0.035 GeV;cos#theta_{GJ} counts / 0.02 GeV", 100, 0.5, 4.0, 100, -1.0, 1.0);

	dHist_pi0costheta_GJ = new TH1F("pi0costheta_GJ", ";#pi0cos#theta_{GJ}", 100, -1.0, 1.0);

	dHist_etacostheta_GJ_EtaPrime = new TH1F("etacostheta_GJ_EtaPrime", ";#etacos#theta_{GJ}", 100, -1.1, 1.1);
	dHist_etaphi_GJ_EtaPrime = new TH1F("etaphi_GJ_EtaPrime", ";#eta#phi_{GJ}", 100, -360.0, 360.0);




	//1-D photon pairs measured variables

	dHist_Photons12_M = new TH1I("Photons12_M", ";Inv_MPhotons12 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons13_M = new TH1I("Photons13_M", ";Inv_MPhotons13 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons14_M = new TH1I("Photons14_M", ";Inv_MPhotons14 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons23_M = new TH1I("Photons23_M", ";Inv_MPhotons23 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons24_M = new TH1I("Photons24_M", ";Inv_MPhotons24 (GeV/c^{2})", 200, 0.0, 1.0);
	dHist_Photons34_M = new TH1I("Photons34_M", ";Inv_MPhotons34 (GeV/c^{2})", 200, 0.0, 1.0);



	//2_D photon pairs

	dHist_allpairs_2D = new TH2F("allpairs_2D", ";M(#gamma#gamma)(GeV/c^{2});M(#gamma#gamma)(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);

	dHist_gg12vsgg34 = new TH2F("gg12vsgg34", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});M(#gamma_{3}#gamma_{4})(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	
	dHist_gg12vsgg34_exclusiveboth = new TH2F("gg12vsgg34_exclusiveboth", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});M(#gamma_{3}#gamma_{4})(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);

	


	dHist_gg12vsgg34_inclusiveboth = new TH2F("gg12vsgg34_inclusiveboth", ";InvM_12(GeV/c^{2});InvM_34(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);

	dHist_gg12vsgg13 = new TH2F("gg12vsgg13", ";InvM_12(GeV/c^{2});InvM_13(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg12vsgg14 = new TH2F("gg12vsgg14", ";InvM_12(GeV/c^{2});InvM_14(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg12vsgg23 = new TH2F("gg12vsgg23", ";InvM_12(GeV/c^{2});InvM_23(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg12vsgg24 = new TH2F("gg12vsgg24", ";InvM_12(GeV/c^{2});InvM_24(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg34vsgg13 = new TH2F("gg34vsgg13", ";InvM_34(GeV/c^{2});InvM_13(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg34vsgg14 = new TH2F("gg34vsgg14", ";InvM_34(GeV/c^{2});InvM_14(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg34vsgg23 = new TH2F("gg34vsgg23", ";InvM_34(GeV/c^{2});InvM_23(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg34vsgg24 = new TH2F("gg34vsgg24", ";InvM_34(GeV/c^{2});InvM_24(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);



	
	dHist_gg13vsgg24 = new TH2F("gg13vsgg24", ";M(#gamma_{1}#gamma_{3})(GeV/c^{2});M(#gamma_{2}#gamma_{4})(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg14vsgg23 = new TH2F("gg14vsgg23", ";M(#gamma_{1}#gamma_{4})(GeV/c^{2});M(#gamma_{2}#gamma_{3})(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);

	dHist_gg13vsgg24_excludepi0 = new TH2F("gg13vsgg24_excludepi0", ";InvM_13(GeV/c^{2});InvM_24(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg14vsgg23_excludepi0 = new TH2F("gg14vsgg23_excludepi0", ";InvM_14(GeV/c^{2});InvM_23(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);

	dHist_gg13vsgg24_onlypi0 = new TH2F("gg13vsgg24_onlypi0", ";InvM_13(GeV/c^{2});InvM_24(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);
	dHist_gg14vsgg23_onlypi0 = new TH2F("gg14vsgg23_onlypi0", ";InvM_14(GeV/c^{2});InvM_23(GeV/c^{2})", 100, 0.0, 1.0, 100, 0.0, 1.0);


	

	


	//2_D photon pairs measured variables
	dHist_gg12vsgg34_M = new TH2F("gg12vsgg34_M", ";InvM_12(GeV/c^{2});InvM_34(GeV/c^{2})", 100, 0.1, 0.2, 100, 0.2, 0.8);
 	dHist_gg13vsgg24_M = new TH2F("gg13vsgg24_M", ";InvM_13(GeV/c^{2});InvM_24(GeV/c^{2})", 100, 0.1, 0.2, 100, 0.1, 0.2);
	dHist_gg14vsgg23_M = new TH2F("gg14vsgg23_M", ";InvM_14(GeV/c^{2});InvM_23(GeV/c^{2})", 100, 0.1, 0.2, 100, 0.1, 0.2);

	dHist_ShowerQuality1= new TH1F("ShowerQuality1", ";ShowerQuality", 100, 0.0, 1.0);
	dHist_ShowerQuality2= new TH1F("ShowerQuality2", ";ShowerQuality", 100, 0.0, 1.0);
	dHist_ShowerQuality3= new TH1F("ShowerQuality3", ";ShowerQuality", 100, 0.0, 1.0);
	dHist_ShowerQuality4= new TH1F("ShowerQuality4", ";ShowerQuality", 100, 0.0, 1.0);
	
	dHist_NumNeutrals = new TH1F("NumNeutrals", ";NumNeutrals", 100, 3.0, 10);
	dHist_NumNeutrals_2pi0veto = new TH1F("NumNeutrals_2pi0veto", ";NumNeutrals", 100, 3.0, 10);
	dHist_NumNeutrals_omegaveto_acc = new TH1F("NumNeutrals_omegaveto_acc", ";NumNeutrals", 100, 3.0, 10);
	dHist_NumNeutrals_omegaveto_prompt = new TH1F("NumNeutrals_omegaveto_prompt", ";NumNeutrals", 100, 3.0, 10);

	dHist_Mggvcosthetapi0GJ_prompt = new TH2F("Mggvcosthetapi0GJ_prompt", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});#pi^{0}_{cos#theta}", 100, 0.11, 0.16, 100, -1., 1.);
	dHist_Mggvphipi0GJ_prompt = new TH2F("Mggvphipi0GJ_prompt", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});#pi^{0}_{#phi}", 100, 0.11, 0.16, 100, -180., 180.);
	dHist_MggvBE_pi0_prompt = new TH2F("MggvBE_pi0_prompt", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});E_{beam}", 100, 0.11, 0.16, 100, 5., 12.);


	dHist_MggvcosthetaetaGJ_prompt = new TH2F("MggvcosthetaetaGJ_prompt", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});#eta^{0}_{cos#theta}", 100, 0.42, 0.65, 100, -1., 1.);
	dHist_MggvphietaGJ_prompt = new TH2F("MggvphietaGJ_prompt", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});#eta^{0}_{#phi}", 100, 0.42, 0.65, 100, -180., 180.);
	dHist_MggvBE_eta_prompt = new TH2F("MggvBE_eta_prompt", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});E_{beam}", 100, 0.42, 0.65, 100, 5., 12.);

	dHist_MggvcosthetapippimGJ_prompt = new TH2F("MggvcosthetapippimGJ_prompt", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});#pippim_{cos#theta}", 100, 0.42, 0.65, 100, -1., 1.);
	dHist_MggvphipippimGJ_prompt = new TH2F("MggvphipippimGJ_prompt", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});#pippim_{#phi}", 100, 0.42, 0.65, 100, -180., 180.);

	//dHist_Mggvcosthetapi0GJ_Acc = new TH2F("Mggvcosthetapi0GJ_Acc", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});#pi^{0}_{cos#theta}", 100, 0.11, 0.16, 100, -1., 1.);
	//dHist_Mggvphipi0GJ_Acc = new TH2F("Mggvphipi0GJ_Acc", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});#pi^{0}_{#phi}", 100, 0.11, 0.16, 100, -180., 180.);
	//dHist_MggvBE_pi0_Acc = new TH2F("MggvBE_pi0_Acc", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});E_{beam}", 100, 0.11, 0.16, 100, 5., 12.);
//
//
	//dHist_MggvcosthetaetaGJ_Acc = new TH2F("MggvcosthetaetaGJ_Acc", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});#eta^{0}_{cos#theta}", 100, 0.42, 0.65, 100, -1., 1.);
	//dHist_MggvphietaGJ_Acc = new TH2F("MggvphietaGJ_Acc", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});#eta^{0}_{#phi}", 100, 0.42, 0.65, 100, -180., 180.);
	//dHist_MggvBE_eta_Acc = new TH2F("MggvBE_eta_Acc", ";M(#gamma_{3}#gamma_{4})(GeV/c^{2});E_{beam}", 100, 0.42, 0.65, 100, 5., 12.);



	dHist_Mggvt_etaprimepi0 = new TH2F("Mggvt_etaprimepi0", ";M(#gamma_{1}#gamma_{2})(GeV/c^{2});T_{#eta^{'}#pi^{0}}", 100, 0.11, 0.16, 100, 0., 5.);

	dHist_MetaprimevcosthetaGJ = new TH2F("MetaprimevcosthetaGJ", ";M(#eta^{'})(GeV/c^{2});#eta^{'}_{cos#theta}", 100, 0.9, 1.5, 100, -1., 1.);
	dHist_MetaprimevphiGJ = new TH2F("MetaprimevphiGJ", ";M(#eta^{'})(GeV/c^{2});#eta^{'}_{#phi}", 100, 0.9, 1.5, 100, -180., 180.);
	dHist_MetaprimevBE = new TH2F("MetaprimevBE", ";M(#eta^{'})(GeV/c^{2});E_{beam}", 100, 0.9, 1.5, 100, 5., 12.);
	dHist_Metaprimevt_etaprimepi0 = new TH2F("Metaprimevt_etaprimepi0", ";M(#eta^{'})(GeV/c^{2});T_{#eta^{'}#pi^{0}}", 100, 0.9, 1.0, 100, 0., 5.);

	dHist_t_pi0 = new TH1F("t_pi0", ";t_{#pi^{0}}", 100, 0.0, 8.0);
	dHist_t_EtaPrime = new TH1F("t_EtaPrime", ";t_{#eta^{'}}", 100, 0.0, 5.0);

	dHist_Metaprimepi0vt_pi0 = new TH2F("Metaprimepi0vt_pi0", ";M(#eta^{'}#pi^{0});t_{#pi^{0}}", 100, 0.8, 3.0, 100, 0.0, 8.0);
	dHist_Metaprimepi0vt_EtaPrime = new TH2F("Metaprimepi0vt_EtaPrime", ";M(#eta^{'}#pi^{0});t_{#eta^{'}}", 100, 1.2, 4.0, 100, 0.0, 8.0);

	dHist_gg12VSBeamDeltaT = new TH2F("gg12VSBeamDeltaT" , ";M(gg);deltaT", 100, 0.05, 0.20, 100, -10., 10.);
	dHist_gg34VSBeamDeltaT = new TH2F("gg34VSBeamDeltaT" , ";M(gg);deltaT", 100, 0.1, 1.0, 100, -10., 10.);
	dHist_etaprimeVSBeamDeltaT = new TH2F("etaprimeVSBeamDeltaT" , ";M(#etaprime);deltaT", 100, 0.8, 1.6, 100, -10., 10.);

	dHist_pippimetavspi0 = new TH2F("pippimetavspi0", ";M(pippimeta);M(g1g2)", 100, 0.8, 1.6, 100, 0.1, 0.2);
	







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

	// Add branches to the qfactortree
	qfactortree->Branch("num_combos", &num_combos, "num_combos/I");
	qfactortree->Branch("pi0mass", &pi0mass, "pi0mass/D");
	qfactortree->Branch("etaprimemass", &etaprimemass, "etaprimemass/D");
	qfactortree->Branch("etaprimepi0mass", &etaprimepi0mass, "etaprimemasspi0/D");

	qfactortree->Branch("pi0costhetaGJ", &pi0costhetaGJ, "pi0costhetaGJ/D");
	qfactortree->Branch("pi0phiGJ", &pi0phiGJ, "pi0phiGJ/D");
	qfactortree->Branch("etaprimecosthetaGJ", &etaprimecosthetaGJ, "etaprimecosthetaGJ/D");
	qfactortree->Branch("etaprimephiGJ", &etaprimephiGJ, "etaprimephiGJ/D");

	qfactortree->Branch("combo_number", &combo_number, "combo_number/I");
	qfactortree->Branch("kinfit_CL", &kinfit_CL, "kinfit_CL/D");
	qfactortree->Branch("chisq_ndf", &chisq_ndf, "chisq_ndf/D");
	qfactortree->Branch("event_num", &event_num, "event_num/I");

	qfactortree->Branch("pippimpi0", &pippimpi0, "pippimpi0/D");
	qfactortree->Branch("pipp", &pipp, "pipp/D");
	qfactortree->Branch("pi0p", &pi0p, "pi0p/D");
	qfactortree->Branch("dt", &dt, "dt/D");

	qfactortree->Branch("etamass", &etamass, "etamass/D");
	qfactortree->Branch("BEa", &BEa, "BEa/D");

	qfactortree->Branch("pi0mass13", &pi0mass13, "pi0mass13/D");
	qfactortree->Branch("pi0mass24", &pi0mass24, "pi0mass24/D");
	qfactortree->Branch("pi0mass14", &pi0mass14, "pi0mass14/D");
	qfactortree->Branch("pi0mass23", &pi0mass23, "pi0mass23/D");

	qfactortree->Branch("etapi0mass", &etapi0mass, "etapi0mass/D");

	qfactortree->Branch("num_unusedshowers", &num_unusedshowers, "num_unusedshowers/I");
	
	qfactortree->Branch("mant", &mant, "mant/D");

	qfactortree->Branch("photon1_sq", &photon1_sq, "photon1_sq/D");
	qfactortree->Branch("photon2_sq", &photon2_sq, "photon2_sq/D");
	qfactortree->Branch("photon3_sq", &photon3_sq, "photon3_sq/D");
	qfactortree->Branch("photon4_sq", &photon4_sq, "photon4_sq/D");


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
	//dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	
	

	
	
	



	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("massetaprime");
}



Bool_t DSelector_pi0etapr__B4_M35_M7_M17::Process(Long64_t locEntry)
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
	//TLorentzVector locProductionX4 = Get_X4_Production();

	


	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
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

	//double  locMydouble = 7.;
	
	//for(int loc_i = 0; loc_i < locMydouble; ++loc_i)

	
	
	
	
	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
		
			continue; // Combo has been cut previously

		double dMinKinFitCL = 1e-03; //throw out bad results within 1 percent of resolution
		double locKinFitConLev = dComboWrapper->Get_ConfidenceLevel_KinFit();
		if(locKinFitConLev < dMinKinFitCL){
			
		  dComboWrapper->Set_IsComboCut(true);
		  continue;
		  }


		
		

		//dFlatTreeInterface->Fill_Fundamental<double>("pi0costheta_GJ",-1,loc_i);

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
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4();
		//Step 2
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		//Step 3
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

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector totalP4_Measured =  locBeamP4_Measured + dTargetP4;		
		TLorentzVector locMissingP4_Measured =totalP4_Measured - 
		 (locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured
		  + locPhoton1P4_Measured + locPhoton2P4_Measured + locPhoton3P4_Measured + locPhoton4P4_Measured);

		TLorentzVector locPhoton12 = locPhoton1P4 + locPhoton2P4;
		TLorentzVector locPhoton13 = locPhoton1P4 + locPhoton3P4;
		TLorentzVector locPhoton14 = locPhoton1P4 + locPhoton4P4;
		TLorentzVector locPhoton23 = locPhoton2P4 + locPhoton3P4;
		TLorentzVector locPhoton24 = locPhoton2P4 + locPhoton4P4;
		TLorentzVector locPhoton34 = locPhoton3P4 + locPhoton4P4;

		

		TLorentzVector locetaprime = locPiPlusP4 + locPiMinusP4	+ locPhoton3P4 + locPhoton4P4;
		TLorentzVector locpippi0eta = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 + locPhoton3P4 + locPhoton4P4;
		TLorentzVector locpi0pimeta =  locPhoton1P4 + locPhoton2P4 + locPiMinusP4 + locPhoton3P4 + locPhoton4P4;
		TLorentzVector locpi0etap =  locPhoton1P4 + locPhoton2P4 +  locPhoton3P4 + locPhoton4P4 + locProtonP4;


		TLorentzVector locetapi0 = 	locPhoton1P4 + locPhoton2P4 +  locPhoton3P4 + locPhoton4P4;

		TLorentzVector locetaprimepi0 = locPiPlusP4 + locPiMinusP4	+ locPhoton3P4 + locPhoton4P4 + locPhoton1P4 + locPhoton2P4;

		TLorentzVector locomega = locPiPlusP4 + locPiMinusP4	+ locPhoton1P4 + locPhoton2P4;

		TLorentzVector locomega_M = locPiPlusP4_Measured + locPiMinusP4_Measured + locPhoton1P4_Measured + locPhoton2P4_Measured;

		TLorentzVector locpippimp = locPiPlusP4 + locPiMinusP4	+ locProtonP4;

		TLorentzVector locpippim = locPiPlusP4 + locPiMinusP4;

		

		TLorentzVector locpippi0p = locPiPlusP4 + 	locPhoton1P4 + locPhoton2P4 + locProtonP4;

		TLorentzVector k1 = locBeamP4 + dTargetP4 - locProtonP4 - locPiPlusP4 - locPiMinusP4;

		TLorentzVector k2 = locPhoton1P4 + locPhoton2P4 +  locPhoton3P4 + locPhoton4P4;


		TLorentzVector locdeltap = locProtonP4 + locPhoton1P4 + locPhoton2P4;

		TLorentzVector locdeltapp = locProtonP4 + locPiPlusP4;

		TLorentzVector locdelta0 = locProtonP4 + locPiMinusP4;

		
		TLorentzVector locetap = 	locProtonP4 +  locPhoton3P4 + locPhoton4P4;
		TLorentzVector locetap_reverse = locBeamP4 + dTargetP4 - locetap;


		TLorentzVector locetap_M = 	locProtonP4_Measured +  locPhoton3P4_Measured + locPhoton4P4_Measured;
		TLorentzVector locetap_reverse_M = locBeamP4_Measured + dTargetP4 - locetap_M;
		double locMassetap_reverse_M = locetap_reverse_M.M();

		TLorentzVector locpipeta = locPiPlusP4 + locPhoton3P4 + locPhoton4P4;

		TLorentzVector locpippi0 = locPiPlusP4 + locPhoton1P4 + locPhoton2P4 ;

		TLorentzVector locpimpi0 = locPiMinusP4 + locPhoton1P4 + locPhoton2P4 ;

		TLorentzVector t_Pi0P4 = locBeamP4 - locPhoton1P4 - locPhoton2P4;
		TLorentzVector t_EtaPrimeP4 = locBeamP4 - locetaprime;

		TLorentzVector tP4 = locBeamP4 - locetaprimepi0;

		

		



		TLorentzVector locPhoton12_Measured = locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector locPhoton13_Measured = locPhoton1P4_Measured + locPhoton3P4_Measured;
		TLorentzVector locPhoton14_Measured = locPhoton1P4_Measured + locPhoton4P4_Measured;
		TLorentzVector locPhoton23_Measured = locPhoton2P4_Measured + locPhoton3P4_Measured;
		TLorentzVector locPhoton24_Measured = locPhoton2P4_Measured + locPhoton4P4_Measured;
		TLorentzVector locPhoton34_Measured = locPhoton3P4_Measured + locPhoton4P4_Measured;


		// etaprime pi0 invariant mass and angular distribution in gottfried jackson frame

		/// CoM Boost
		TLorentzVector locEtaPrimeP4=locPiPlusP4 + locPiMinusP4 + locPhoton34 ; //define etaprime 4-vector
		 TLorentzVector locEtaPrimePi0P4=locEtaPrimeP4+locPhoton12; //define etaprime-pi0 4-vector
		TLorentzVector locCoMP4=locBeamP4 + dTargetP4; //define a 4-vector for beam and target together 
		TVector3 boostCoM =-(locCoMP4.Vect())*(1.0/locCoMP4.E()); //calculate beta for centre of mass system 

        //redefine lab 4-vectors to centre of mass 4-vectors
		TLorentzVector locBeamP4_CM=locBeamP4;
		TLorentzVector locEtaPrimeP4_CM=locEtaPrimeP4;
		TLorentzVector locPi0P4_CM=locPhoton12;	
		TLorentzVector locEtaP4_CM=locPhoton34;	 // define eta in CM
		TLorentzVector locpippim_CM = locpippim;
		TLorentzVector locEtaPrimePi0P4_CM=locEtaPrimePi0P4;

		//boost in centre of mass 
		locEtaPrimePi0P4_CM.Boost(boostCoM);
		locBeamP4_CM.Boost(boostCoM);
		locEtaPrimeP4_CM.Boost(boostCoM);
		locPi0P4_CM.Boost(boostCoM);

		locEtaP4_CM.Boost(boostCoM); // boost in CM
		locpippim_CM.Boost(boostCoM);

        //GJ Boost
		TVector3 boostGJ =-(locEtaPrimePi0P4_CM.Vect())*(1.0/locEtaPrimePi0P4_CM.E()); //calculate beta for etaprimepi0 system

		TVector3 boostGJ_EtaPrime =-(locEtaPrimeP4_CM.Vect())*(1.0/locEtaPrimeP4_CM.E()); //calculate beta for etaprime system

        //redefine centre of mass 4-vectors to 4-vectors in rest frame of etaprimepi0 system
		TLorentzVector locEtaPrimePi0P4_GJ=locEtaPrimePi0P4_CM;
		TLorentzVector locEtaPrimeP4_GJ=locEtaPrimeP4_CM;

		TLorentzVector locEtaPrimeP4_GJ_Etaprime=locEtaPrimeP4_CM;

		TLorentzVector locPi0P4_GJ=locPi0P4_CM;
 		TLorentzVector locBeamP4GJ=locBeamP4_CM;

		TLorentzVector locBeamP4GJ_EtaPrime=locBeamP4_CM; //define beam in rest frame of etaprime

		TLorentzVector locEtaP4_GJ_EtaPrime=locEtaP4_CM; // define eta in GJ frame of etaprime

		TLorentzVector locpippim_GJ_EtaPrime=locpippim_CM; // define pippim  in GJ of etaprime

        //boost in rest frame of etaprimepi0 system
		locEtaPrimePi0P4_GJ.Boost(boostGJ);
		locBeamP4GJ.Boost(boostGJ);
		locEtaPrimeP4_GJ.Boost(boostGJ);
		locPi0P4_GJ.Boost(boostGJ);

		//boost in rest frame of etaprime
		locEtaPrimeP4_GJ_Etaprime.Boost(boostGJ_EtaPrime);
		locBeamP4GJ_EtaPrime.Boost(boostGJ_EtaPrime);
		locEtaP4_GJ_EtaPrime.Boost(boostGJ_EtaPrime);
		locpippim_GJ_EtaPrime.Boost(boostGJ_EtaPrime);

		


        TVector3 z_GJ;
		z_GJ.SetXYZ(locBeamP4GJ.X(),locBeamP4GJ.Y(),locBeamP4GJ.Z());//z GJ
		TVector3 z_hat_GJ=z_GJ.Unit();  //beam direction (z-direction) 

		TVector3 z_GJ_EtaPrime;
		z_GJ_EtaPrime.SetXYZ(locBeamP4GJ_EtaPrime.X(),locBeamP4GJ_EtaPrime.Y(),locBeamP4GJ_EtaPrime.Z());//z GJ_EtaPrime
		TVector3 z_hat_GJ_EtaPrime=z_GJ_EtaPrime.Unit();  //beam direction (z-direction) 

		TVector3 y_GJ=locBeamP4_CM.Vect().Cross(locEtaPrimePi0P4_CM.Vect());  //y-direction for rest frame of etaprimepi0 system
		TVector3 y_hat_GJ=y_GJ.Unit();    

		TVector3 y_GJ_EtaPrime=locBeamP4_CM.Vect().Cross(locEtaPrimeP4_CM.Vect());  //y-direction for rest frame of etaprime system
		TVector3 y_hat_GJ_EtaPrime=y_GJ_EtaPrime.Unit();    

		TVector3 x_hat_GJ=y_hat_GJ.Cross(z_hat_GJ);//x hat GJ 

		TVector3 x_hat_GJ_EtaPrime=y_hat_GJ_EtaPrime.Cross(z_hat_GJ_EtaPrime);//x hat GJ_EtaPrime 
		
		TVector3 vetaprime(locEtaPrimeP4_GJ.Vect()*x_hat_GJ, locEtaPrimeP4_GJ.Vect()*y_hat_GJ, // etaprime 3-vector in etaprimepi0 rest frame
		locEtaPrimeP4_GJ.Vect()*z_hat_GJ);


		TVector3 vpi0(locPi0P4_GJ.Vect()*x_hat_GJ, locPi0P4_GJ.Vect()*y_hat_GJ, // pi0 3-vector in etaprimepi0 rest frame
		locPi0P4_GJ.Vect()*z_hat_GJ);

		TVector3 veta(locEtaP4_GJ_EtaPrime.Vect()*x_hat_GJ_EtaPrime, locEtaP4_GJ_EtaPrime.Vect()*y_hat_GJ_EtaPrime, // eta 3-vector in etaprime rest frame
		locEtaP4_GJ_EtaPrime.Vect()*z_hat_GJ_EtaPrime);


		TVector3 vpippim(locpippim_GJ_EtaPrime.Vect()*x_hat_GJ_EtaPrime, locpippim_GJ_EtaPrime.Vect()*y_hat_GJ_EtaPrime, // pippim 3-vector in etaprime rest frame
		locpippim_GJ_EtaPrime.Vect()*z_hat_GJ_EtaPrime);



		double cosThetaEta_GJ_EtaPrime = veta.CosTheta(); 
		double PhiEta_GJ_EtaPrime = veta.Phi()*180./TMath::Pi(); 

		double cosThetapippim_GJ_EtaPrime = vpippim.CosTheta(); 
		double Phipippim_GJ_EtaPrime = vpippim.Phi()*180./TMath::Pi(); 


		double BE_GJ = locBeamP4GJ.E();
		double BE_GJ_EtaPrime = locBeamP4GJ_EtaPrime.E();


		 //cout << "cosThetapippim_GJ_EtaPrime" << cosThetapippim_GJ_EtaPrime << endl;
		 //cout << "Phipippim_GJ_EtaPrime" << Phipippim_GJ_EtaPrime << endl;

		 dHist_etacostheta_GJ_EtaPrime->Fill(cosThetaEta_GJ_EtaPrime);
		 dHist_etaphi_GJ_EtaPrime->Fill(PhiEta_GJ_EtaPrime);





        //double cosTheta = vetaprime.CosTheta();  //cosine theta which is the variable I use (y-axis in my plot)
		double costhetaetaprime_GJ = vetaprime.CosTheta();
		double phietaprime_GJ = vetaprime.Phi()*180./TMath::Pi(); //phi of pi0 in GJ frame

		double locEtaPrimePi0P4mass =  locEtaPrimePi0P4.M(); //inavriant mass of etaprime pi0 system (x-axis in my plot)


		double cosThetapi0_GJ = vpi0.CosTheta(); //cosine theta for pi0 in GJ frame
		double phipi0_GJ = vpi0.Phi()*180./TMath::Pi(); //phi of pi0 in GJ frame

	
    
    double locEtaPrimeP4mass = locEtaPrimeP4.M();

		

		

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		//dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

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

		 

		/******************************************** ACCIDENTAL SUBRACTION INFO *******************************************/
		
		// measured tagger time for combo
		/*TLorentzVector locBeam_X4_Measured = dComboBeamWrapper->Get_X4_Measured(); 

		// measured RF time for combo
		double locRFTime = dComboWrapper->Get_RFTime_Measured(); 

		// time difference between tagger and RF (corrected for production vertex position relative to target center)
		double locBeamDeltaT = locBeam_X4_Measured.T() - (locRFTime + (locBeam_X4_Measured.Z() - dTargetCenter.Z())/29.9792458); 
		dHist_BeamDeltaT->Fill(locBeamDeltaT);*/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		//Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		 Double_t locBeamDeltaT = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		 dHist_BeamDeltaT->Fill(locBeamDeltaT);

		// calculate accidental subtraction weight based on time difference 
		double locWeight = 0.; // weight to accidentally subtracted histgorams
		bool locAccid = false; // flag to fill separate prompt and accidental histograms for later subtraction

		if(fabs(locBeamDeltaT) < 0.5*4.008) { // prompt signal recieves a weight of 1
		  locWeight = 1.;
		  locAccid = false;
		}
                else { // accidentals recieve a weight of 1/# RF bunches included in TTree (4 in this case)
		  locWeight = -1./4.;
		  locAccid = true;


		  }



		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		double locMissingEnergy = locMissingP4_Measured.E();

		
		double locMasspippi0eta = locpippi0eta.M();
		double locMasspi0pimeta = locpi0pimeta.M();
		double locMasspi0etap = locpi0etap.M();

		double locEnergyPhoton1 = locPhoton1P4.E();
		double locEnergyPhoton2 = locPhoton2P4.E();

		double locEnergyPhoton3 = locPhoton3P4.E();
		double locEnergyPhoton4 = locPhoton4P4.E();

		double locMassPhoton12 = locPhoton12.M();
		double locMassPhoton13 = locPhoton13.M();
		double locMassPhoton14 = locPhoton14.M();
		double locMassPhoton23 = locPhoton23.M();
		double locMassPhoton24 = locPhoton24.M();
		double locMassPhoton34 = locPhoton34.M();

		
		double locMassetaprime = locetaprime.M();
		double locMassetapi0 = locetapi0.M();
		double locMassetaprimepi0 = locetaprimepi0.M();
		double locMassomega = locomega.M();

		double locMassomega_M = locomega_M.M();

		double locMasspippimp = locpippimp.M();

		double locMasspippim = locpippim.M();

		double locMasspimpi0_M = locpimpi0.M();

		
		double locMasspippi0p = locpippi0p.M();

		double locMassdeltap = locdeltap.M();
		double locMassdeltapp = locdeltapp.M();
		double locMassdelta0 = locdelta0.M();

		double locMassetap = locetap.M();
		double locMassetap_reverse = locetap_reverse.M();

		

		double locMasspipeta = locpipeta.M();
		double locMasspippi0 = locpippi0.M();  

		double locMassPhoton12_Measured= locPhoton12_Measured.M();
		double locMassPhoton13_Measured= locPhoton13_Measured.M();
		double locMassPhoton14_Measured= locPhoton14_Measured.M();
		double locMassPhoton23_Measured= locPhoton23_Measured.M();
		double locMassPhoton24_Measured= locPhoton24_Measured.M();
		double locMassPhoton34_Measured= locPhoton34_Measured.M();


		double photon1_shower_quality = dPhoton1Wrapper->Get_Shower_Quality();
		dHist_ShowerQuality1->Fill(photon1_shower_quality);
		double photon2_shower_quality = dPhoton2Wrapper->Get_Shower_Quality();
		dHist_ShowerQuality2->Fill(photon2_shower_quality);
		double photon3_shower_quality = dPhoton3Wrapper->Get_Shower_Quality();
		dHist_ShowerQuality3->Fill(photon3_shower_quality);
		double photon4_shower_quality = dPhoton4Wrapper->Get_Shower_Quality();
		dHist_ShowerQuality4->Fill(photon4_shower_quality);

		// get num of neutrals 
		double num_neutrals = Get_NumNeutralHypos();
		dHist_NumNeutrals->Fill(num_neutrals);
		
		double t_pi0 = (-1)*t_Pi0P4.M2();
		double t_EtaPrime = (-1)*t_EtaPrimeP4.M2();

		double t = (-1)*tP4.M2();

		double tr_pi0mass = locMassPhoton12;

		// get extra showers
		double NumUnusedShowers  = dComboWrapper->Get_NumUnusedShowers();

		// Assign values to the previously defined variables for qfactortree
		num_combos =  Get_NumCombos();
		combo_number = loc_i;

		pi0mass = locMassPhoton12;
		etaprimemass = locEtaPrimeP4mass;
		etaprimepi0mass = locEtaPrimePi0P4mass;
		etamass = locMassPhoton34;

		pi0costhetaGJ = cosThetapi0_GJ;
		pi0phiGJ = phipi0_GJ;

		etaprimecosthetaGJ = costhetaetaprime_GJ;
		etaprimephiGJ = phietaprime_GJ;

		pippimpi0 = locMassomega;
		pipp = locMassdeltapp;
		pi0p = locMassdeltap;
		dt = locBeamDeltaT;

		BEa = BE;
		

		kinfit_CL = locKinFitConLev;
		chisq = dComboWrapper->Get_ChiSq_KinFit();  //cout <<  chisq << "chisq" << endl;
		ndf = dComboWrapper->Get_NDF_KinFit();     //cout <<  ndf << "ndf" << endl;
		chisq_ndf = chisq/ndf;                     //cout << chisq_ndf << "chisq/ndf" << endl;
		event_num = Get_EventNumber();

		// variables for 2pi0 , t and extra showers cut experimentation

		etapi0mass = locMassetapi0;

		pi0mass13 = locMassPhoton13;
		pi0mass24 = locMassPhoton24;
		pi0mass14 = locMassPhoton14;
		pi0mass23 = locMassPhoton23;

		num_unusedshowers = NumUnusedShowers;
		mant = t;
		
		photon1_sq = photon1_shower_quality;
		photon2_sq = photon2_shower_quality;
		photon3_sq = photon3_shower_quality;
		photon4_sq = photon4_shower_quality;
		



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
			dHist_MissingMassSquared->Fill(locMissingMassSquared);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);

			dHist_MissingEnergy->Fill(locMissingEnergy);

			

			/******************* Define Bools for event selection *****************/

			Bool_t coherentbeam = (BE > 8.0) && (BE < 9.0);
			Bool_t MissingMassSquaredcut = (locMissingMassSquared > -0.02) && (locMissingMassSquared < 0.02);
			Bool_t FCAL_showerqualitycut = (photon1_shower_quality > 0.5) && (photon2_shower_quality > 0.5)
											 && (photon3_shower_quality > 0.5) && (photon3_shower_quality > 0.5);
			
			Bool_t pi013 = ((locMassPhoton13 > 0.12) && (locMassPhoton13 < 0.15)); 
			Bool_t pi024 = ((locMassPhoton24 > 0.12) && (locMassPhoton24 < 0.15));
			Bool_t pi014 = ((locMassPhoton14 > 0.12) && (locMassPhoton14 < 0.15));
			Bool_t pi023 = ((locMassPhoton23 > 0.12) && (locMassPhoton23 < 0.15));

			Bool_t etapi0masswindow = ((locMassPhoton34 > 0.48) && (locMassPhoton34 < 0.60))
										 && ((locMassPhoton12 > 0.08) && (locMassPhoton12 < 0.18));

			Bool_t etaprimemasswindow = (locEtaPrimeP4mass > 0.85) && (locEtaPrimeP4mass < 1.05);
			Bool_t baryons = (locMassdeltap < 1.35)  || (locMassdeltapp < 1.35);
			Bool_t omega = (locMassomega > 0.75) && (locMassomega < 0.85);
			Bool_t aboveomega = (locMassomega > 0.85);
			
			
			//if  (MissingMassSquaredcut && coherentbeam && FCAL_showerqualitycut && !locAccid)
			if  (MissingMassSquaredcut && coherentbeam  && !locAccid)
			//if  (MissingMassSquaredcut && coherentbeam && FCAL_showerqualitycut) //accidental study
			{

			// Fill the qfactor tree	
			qfactortree->Fill(); 


			//1-D photon pairs kinfitted variables
 
			dHist_Photons12->Fill(locMassPhoton12);
			dHist_Photons13->Fill(locMassPhoton13);
			dHist_Photons14->Fill(locMassPhoton14);
			dHist_Photons23->Fill(locMassPhoton23);
			dHist_Photons24->Fill(locMassPhoton24);
			dHist_Photons34->Fill(locMassPhoton34);

			

			//1-D photon pairs measured variables

			dHist_Photons12_M->Fill(locMassPhoton12_Measured);
			dHist_Photons13_M->Fill(locMassPhoton13_Measured);
			dHist_Photons14_M->Fill(locMassPhoton14_Measured);
			dHist_Photons23_M->Fill(locMassPhoton23_Measured);
			dHist_Photons24_M->Fill(locMassPhoton24_Measured);
			dHist_Photons34_M->Fill(locMassPhoton34_Measured);

			//2-D photon pairs

			dHist_gg12vsgg34->Fill(locMassPhoton12,locMassPhoton34);

			

			dHist_gg12vsgg13->Fill(locMassPhoton12,locMassPhoton13);
			dHist_gg12vsgg14->Fill(locMassPhoton12,locMassPhoton14);

			dHist_gg12vsgg23->Fill(locMassPhoton12,locMassPhoton23);
			dHist_gg12vsgg24->Fill(locMassPhoton12,locMassPhoton24);

			dHist_gg34vsgg13->Fill(locMassPhoton34,locMassPhoton13);
			dHist_gg34vsgg14->Fill(locMassPhoton34,locMassPhoton14);

			dHist_gg34vsgg23->Fill(locMassPhoton34,locMassPhoton23);
			dHist_gg34vsgg24->Fill(locMassPhoton34,locMassPhoton24);

			dHist_allpairs_2D->Fill(locMassPhoton13,locMassPhoton24);
			dHist_allpairs_2D->Fill(locMassPhoton14,locMassPhoton23);

			dHist_gg13vsgg24->Fill(locMassPhoton13,locMassPhoton24);
			dHist_gg14vsgg23->Fill(locMassPhoton14,locMassPhoton23);

			//2-D photon pairs measured variables

			dHist_gg12vsgg34_M->Fill(locMassPhoton12_Measured,locMassPhoton34_Measured);
			dHist_gg13vsgg24_M->Fill(locMassPhoton13_Measured,locMassPhoton24_Measured);
			dHist_gg14vsgg23_M->Fill(locMassPhoton14_Measured,locMassPhoton23_Measured);

			
			// 2pi0 veto 
			//if (!((pi013 && pi024) || (pi014 && pi023)))
			
			
			//{
				//qfactortree->Fill(); 
				
				

			dHist_gg13vsgg24_excludepi0->Fill(locMassPhoton13,locMassPhoton24);
			dHist_gg14vsgg23_excludepi0->Fill(locMassPhoton14,locMassPhoton23);
			dHist_gg12vsgg34_exclusiveboth->Fill(locMassPhoton12,locMassPhoton34);
			dHist_Photons12_exclusiveboth->Fill(locMassPhoton12);
			dHist_Photons34_exclusiveboth->Fill(locMassPhoton34);
			dHist_pippimetavspi0->Fill(locMassetaprime, locMassPhoton12);

			

			dHist_gg12VSBeamDeltaT->Fill(locMassPhoton12, locBeamDeltaT);	
			 dHist_gg34VSBeamDeltaT->Fill(locMassPhoton34, locBeamDeltaT);

			 

			dHist_Photons12_exclusiveboth_prompt->Fill(locMassPhoton12);
					dHist_Photons34_exclusiveboth_prompt->Fill(locMassPhoton34);
					

			
			dHist_NumNeutrals_2pi0veto->Fill(num_neutrals);
			
			
			/*************************  etapi0 mass window  **********************/
				
			//if (etapi0masswindow)

			//{
					
			
			dHist_pi0costheta_GJ->Fill(cosThetapi0_GJ);
			

			

			dHist_Mggvcosthetapi0GJ_prompt->Fill(locMassPhoton12 , cosThetapi0_GJ);
			dHist_Mggvphipi0GJ_prompt->Fill(locMassPhoton12 , phipi0_GJ);
			dHist_MggvBE_pi0_prompt->Fill(locMassPhoton12 , BE);

			dHist_MggvcosthetaetaGJ_prompt->Fill(locMassPhoton34 , cosThetaEta_GJ_EtaPrime);
			dHist_MggvphietaGJ_prompt->Fill(locMassPhoton34 , PhiEta_GJ_EtaPrime);
			dHist_MggvBE_eta_prompt->Fill(locMassPhoton34 , BE);

			dHist_MggvcosthetapippimGJ_prompt->Fill(locMassPhoton34 , cosThetapippim_GJ_EtaPrime);
			dHist_MggvphipippimGJ_prompt->Fill(locMassPhoton34 , Phipippim_GJ_EtaPrime);
			dHist_pippimeta_prompt->Fill(locEtaPrimeP4mass);
		
			
			
			


			dHist_MetaprimevcosthetaGJ->Fill(locMassetaprime , costhetaetaprime_GJ);
			dHist_MetaprimevphiGJ->Fill(locMassetaprime , phietaprime_GJ);
			dHist_MetaprimevBE->Fill(locMassetaprime , BE);
			
			dHist_etaprimeVSBeamDeltaT->Fill(locMassetaprime, locBeamDeltaT);

			dHist_pippimetavpi0p->Fill(locMassdeltap, locMassetaprime);
			dHist_pippi0etavpimp->Fill(locMassdelta0, locMasspippi0eta);
			dHist_pi0pimetavpipp->Fill(locMassdeltapp, locMasspi0pimeta);
			dHist_pi0etapvpippim->Fill(locMasspippim, locMasspi0etap);

			dHist_Photons1v2->Fill(locEnergyPhoton1, locEnergyPhoton2);
			dHist_Photons3v4->Fill(locEnergyPhoton3, locEnergyPhoton4);

			dHist_pippimetapi0->Fill(locMassetaprimepi0);
			dHist_etapi0->Fill(locMassetapi0);
			dHist_pippi0p->Fill(locMasspippi0p);
			dHist_pippimp->Fill(locMasspippimp);
			
			dHist_pippimpi0->Fill(locMassomega);
			
			dHist_pi0p->Fill(locMassdeltap);
			dHist_pipp->Fill(locMassdeltapp);
			dHist_pimp->Fill(locMassdelta0);
			dHist_pippim->Fill(locMasspippim);

			dHist_pimpi0->Fill(locMasspimpi0_M);

			dHist_etap->Fill(locMassetap);
			dHist_etap_reverse->Fill(locMassetap_reverse);

			
			//2-D invariant mass combinations for eta and pi0 selected events

			dHist_pippimetapi0vpippimeta->Fill(locMassetaprimepi0, locMassetaprime);
			dHist_pippimetapi0vetapi0->Fill(locMassetaprimepi0, locMassetapi0);
			dHist_pippimetapi0vpippimpi0->Fill(locMassetaprimepi0, locMassomega);

			

			dHist_pippimetapi0vpippimp->Fill(locMassetaprimepi0, locMasspippimp);
			dHist_pippimetapi0vpippi0p->Fill(locMassetaprimepi0, locMasspippi0p);
			dHist_pippimetapi0vpipp->Fill(locMassetaprimepi0,locMassdeltapp);

			dHist_pipetavspippi0->Fill(locMasspipeta,locMasspippi0);
			dHist_etapi0vspippim->Fill(locMassetapi0,locMasspippim);
			dHist_pippimetavpippimpi0->Fill(locMassomega,locMassetaprime);

			dHist_pippimeta->Fill(locMassetaprime);	
			
			//select etaprime events 
			//if (etaprimemasswindow) 
			//{
				dHist_etaprimepi0->Fill(locEtaPrimePi0P4mass);
			dHist_InvMetapvcostheta->Fill(locEtaPrimePi0P4mass, costhetaetaprime_GJ);

			
			dHist_NumNeutrals_omegaveto_prompt->Fill(num_neutrals);


					//Fill the qfactortree
					 //qfactortree->Fill(); 
			//} //etaprime mass window closes
					
					 
					    
				

				//omega anticut
				if (!aboveomega)
				{
					dHist_pippimetavpippimpi0_omega->Fill( locMassomega, locMassetaprime);
					dHist_pippimeta_omega->Fill(locMassetaprime);

				} //omega anticut ends here

				//plot t_distributions for double regge study			
				
				dHist_t_pi0->Fill(t_pi0);
				dHist_Metaprimepi0vt_pi0->Fill(locMassetaprimepi0, t_pi0);

					
				//} // pi0 and eta mass window ends here
			//} // 2 pi0 veto ends here
			
			/****************** Display 2pi0 events  **********************/

			if (pi013 && pi024) dHist_gg13vsgg24_onlypi0->Fill(locMassPhoton13,locMassPhoton24);

			if (pi014 && pi023) dHist_gg14vsgg23_onlypi0->Fill(locMassPhoton14,locMassPhoton23);
			

			} // MissingMasssquaredcut,coherentbeam,showerqualitycut, coincidence time cut ends
			
		} // end of "if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())"

		

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

			/*
		
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		//dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		//TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		//dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		
		for(int loc_i = 0; loc_i < locMyInt_Flat; ++loc_i)
		{
			//dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			//TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			//dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
			dFlatTreeInterface->Fill_Fundamental<Int_t>("pi0costheta_GJ",cosThetapi0_GJ,loc_i);
			
		} 
		
		*/

			
			
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



void DSelector_pi0etapr__B4_M35_M7_M17::Finalize(void)
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
