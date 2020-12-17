#ifndef DSelector_pi0etapr__B4_M35_M7_M17_h
#define DSelector_pi0etapr__B4_M35_M7_M17_h

#include <iostream>
#include <fstream>
using namespace std;

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TCanvas.h"
#include "TTree.h"

class DSelector_pi0etapr__B4_M35_M7_M17 : public DSelector
{
	public:

		DSelector_pi0etapr__B4_M35_M7_M17(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pi0etapr__B4_M35_M7_M17(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);
	

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dPiPlusWrapper;

		//Step 3
		DParticleComboStep* dStep3Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;
		DNeutralParticleHypothesis* dPhoton4Wrapper;

		
		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:

		TH1I* dHist_MissingMassSquared;

		TH1I* dHist_MissingEnergy;

		TH1I* dHist_BeamEnergy;

		TH1I* dHist_Photons12;
		TH1I* dHist_Photons13;

		TH1I* dHist_k;

		TH1I* dHist_Photons12_exclusive;
		TH1I* dHist_Photons34_exclusive;

		TH1F* dHist_Photons12_exclusiveboth;
		TH1F* dHist_Photons34_exclusiveboth;

		TH1F* dHist_Photons12_exclusiveboth_accid;
		TH1F* dHist_Photons34_exclusiveboth_accid;
		TH1F* dHist_pippimeta_accid;

		TH1F* dHist_Photons12_exclusiveboth_prompt;
		TH1F* dHist_Photons34_exclusiveboth_prompt;
		TH1F* dHist_pippimeta_prompt;


		TH1I* dHist_Photons12_exclusiveboth_SC;
		TH1I* dHist_Photons34_exclusiveboth_SC;

		TH1I* dHist_Photons12_belowthreshold;
		TH1I* dHist_Photons34_belowthreshold;

		

		TH1I* dHist_Photons12_inclusiveboth;
		TH1I* dHist_Photons34_inclusiveboth;

		TH1I* dHist_Photons14;
		TH1I* dHist_Photons23;
		TH1I* dHist_Photons24;
		TH1I* dHist_Photons34;

		TH1I* dHist_pippimeta;
		TH1I* dHist_pippimeta_cutbr;
		TH1I* dHist_pippimeta_cutomega;
		TH1I* dHist_pippimeta_aboveomega;

		TH1I* dHist_etaprimepi0;
		TH1I* dHist_etaprimepi0_cutbr;
		TH1I* dHist_etaprimepi0_cutomega;
		TH1I* dHist_etaprimepi0_aboveomega;

		TH2F* dHist_InvMetapvcostheta;
		TH2F* dHist_InvMetapvcostheta_cutbr;
		TH2F* dHist_InvMetapvcostheta_cutomega;
		TH2F* dHist_InvMetapvcostheta_aboveomega;



		TH1I* dHist_pippimeta_SC;
		TH1I* dHist_pippimeta_cutomega_SC;
		
		TH1I* dHist_pippimeta_pi0selectiononly;

		TH1I* dHist_pippimeta_cutetapi0;
		TH1I* dHist_pippimeta_cutpippimpi0;
		
		TH1I* dHist_pippimeta_omega;
		TH1I* dHist_pippimeta_cut_pipp_pi0p;

		TH1I* dHist_pippimetapi0;

		TH1I* dHist_pippimeta_belowthreshold;
		TH1I* dHist_pippimeta_abovethreshold;

		TH1I* dHist_pippimeta_beloweta;
		TH1I* dHist_pippimeta_aboveeta;
		TH1I* dHist_pippimeta_eta;

		TH1I* dHist_pippimpi0;
		TH1I* dHist_pippimpi0_SC;
		TH1I* dHist_pippimpi0_M;

		TH1I* dHist_pippimpi0_pi0selectiononly;
		TH1I* dHist_pippimpi0_belowthreshold;

		TH1I* dHist_pippimp;
		TH1I* dHist_pippi0p;
		TH1I* dHist_etapi0;

		TH1I* dHist_etapi0_belowthreshold;

		TH1I* dHist_pi0p;
		TH1I* dHist_pipp;

		TH1I* dHist_pimp;

		TH1I* dHist_pippim;
		TH1I* dHist_pimpi0;

		TH1I* dHist_eta_cleanetap;
		TH1I* dHist_eta_cleanetap_alleta;
		TH1I* dHist_etap;
		TH1I* dHist_etap_reverse;
		TH1I* dHist_etap_reverse_M;

		
		TH1I* dHist_etaprimepi0_cutomega_SC;

		TH1F* dHist_pi0costheta_GJ;

		TH2F* dHist_pippimetapi0vpippimeta;

		TH2F* dHist_pippimetavpi0p;
		TH2F* dHist_pippi0etavpimp;
		TH2F* dHist_pi0pimetavpipp;
		TH2F* dHist_pi0etapvpippim;

		TH2F* dHist_pippimetapi0vpippimeta_cutetapi0;

		TH2F* dHist_pippimetavpippimpi0;
		TH2F* dHist_pippimetavpippimpi0_SC;

		
		TH2F* dHist_InvMetapvcostheta_SC;
		
		TH2F* dHist_pippimetavpippimpi0_cutomega;
		TH2F* dHist_pippimetavpippimpi0_omega;

		TH2F* dHist_pippimetapi0vpippimeta_aboveomega;

		TH2F* dHist_pippimetapi0vpippimeta_cutomega;

		TH2F* dHist_pippimetapi0vpippimeta_cut_pipp_pi0p;

		TH2F* dHist_pippimetapi0vpippimpi0;
		TH2F* dHist_pippimetapi0vpippimp;
		TH2F* dHist_pippimetapi0vpippi0p;
		TH2F* dHist_pippimetapi0vetapi0;
		TH2F* dHist_pippimetapi0vpipp;

		TH2F* dHist_pipetavspippi0;
		TH2F* dHist_etapi0vspippim;

		TH2F* dHist_Photons1v2_SC;
		TH2F* dHist_Photons1v2;

		TH2F* dHist_Photons3v4_SC;
		TH2F* dHist_Photons3v4;

		TH1I* dHist_Photons12_M;
		TH1I* dHist_Photons13_M;
		TH1I* dHist_Photons14_M;
		TH1I* dHist_Photons23_M;
		TH1I* dHist_Photons24_M;
		TH1I* dHist_Photons34_M;


		TH2F* dHist_gg12vsgg34;

		TH1I* dHist_allpairs_1D;
		TH2F* dHist_allpairs_2D;

		TH2F* dHist_gg12vsgg34_exclusive;
		TH2F* dHist_gg12vsgg34_exclusiveboth;

		TH2F* dHist_gg12vsgg34_Mcondition;

		TH2F* dHist_gg12vsgg34_inclusiveboth;

		TH2F* dHist_gg13vsgg24;
		TH2F* dHist_gg14vsgg23;

		TH2F* dHist_gg13vsgg24_excludepi0;
		TH2F* dHist_gg14vsgg23_excludepi0;

		TH2F* dHist_gg13vsgg24_onlypi0;
		TH2F* dHist_gg14vsgg23_onlypi0;


		TH2F* dHist_gg12vsgg13;
		TH2F* dHist_gg12vsgg14;

		TH2F* dHist_gg12vsgg23;
		TH2F* dHist_gg12vsgg24;

		TH2F* dHist_gg34vsgg13;
		TH2F* dHist_gg34vsgg14;

		TH2F* dHist_gg34vsgg23;
		TH2F* dHist_gg34vsgg24;

		TH2F* dHist_gg12vsgg34_M;
		TH2F* dHist_gg13vsgg24_M;
		TH2F* dHist_gg14vsgg23_M;

		TH3F* dHist_gg12vsgg133d; 
		TCanvas* c2;


		TH1F* dHist_BeamDeltaT;

		TH1F* dHist_ShowerQuality1;
		TH1F* dHist_ShowerQuality2;
		TH1F* dHist_ShowerQuality3;
		TH1F* dHist_ShowerQuality4;

		TH1F* dHist_NumNeutrals;
		TH1F* dHist_NumNeutrals_2pi0veto;
		TH1F* dHist_NumNeutrals_omegaveto_acc;
		TH1F* dHist_NumNeutrals_omegaveto_prompt;

		TH1F* dHist_etacostheta_GJ_EtaPrime;
		TH1F* dHist_etaphi_GJ_EtaPrime;

		TH2F* dHist_Metaprimepi0vt_pi0;
		TH2F* dHist_Metaprimepi0vt_EtaPrime;

		TH1F* dHist_t_pi0;
		TH1F* dHist_t_EtaPrime;







		TH2F* dHist_Mggvcosthetapi0GJ_prompt;
		TH2F* dHist_Mggvphipi0GJ_prompt;
		TH2F* dHist_MggvBE_pi0_prompt;

		TH2F* dHist_Mggvcosthetapi0GJ_Acc;
		TH2F* dHist_Mggvphipi0GJ_Acc;
		TH2F* dHist_MggvBE_pi0_Acc;

		TH2F* dHist_MggvcosthetaetaGJ_prompt;
		TH2F* dHist_MggvphietaGJ_prompt;
		TH2F* dHist_MggvBE_eta_prompt;

		TH2F* dHist_MggvcosthetapippimGJ_prompt;
		TH2F* dHist_MggvphipippimGJ_prompt;

		TH2F* dHist_MggvcosthetaetaGJ_Acc;
		TH2F* dHist_MggvphietaGJ_Acc;
		TH2F* dHist_MggvBE_eta_Acc;


		TH2F* dHist_Mggvt_etaprimepi0;

		TH2F* dHist_MetaprimevcosthetaGJ;
		TH2F* dHist_MetaprimevphiGJ;
		TH2F* dHist_MetaprimevBE;
		TH2F* dHist_Metaprimevt_etaprimepi0;

		TH2F* dHist_gg12VSBeamDeltaT;
		TH2F* dHist_gg34VSBeamDeltaT;
		TH2F* dHist_etaprimeVSBeamDeltaT;

		TH2F* dHist_pippimetavspi0;

		// extra charged tracks and showers
		TH1F * dHist_NumExS;
		TH1F * dHist_NumExT;  

		//tree stuff (variables for qfactor analysis) goes here
		TFile *fileout; // file for the outputtree
		TTree *qfactortree; //qfactortree

		Int_t num_combos;
		Int_t combo_number;
		Int_t ndf;
		Int_t event_num;

    	Double_t  pi0mass;
		Double_t  etaprimemass;
		Double_t  etaprimepi0mass;
		Double_t etamass; 

		Double_t pi0costhetaGJ;
		Double_t pi0phiGJ;
		Double_t etaprimecosthetaGJ;
		Double_t etaprimephiGJ;

		Double_t kinfit_CL;
		Double_t chisq;
		Double_t chisq_ndf;

		Double_t pippimpi0;
		Double_t pipp;
		Double_t pi0p;
		Double_t dt;

		
		Double_t BEa;
		

		
		



	ClassDef(DSelector_pi0etapr__B4_M35_M7_M17, 0);
};

void DSelector_pi0etapr__B4_M35_M7_M17::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));

	//Step 3
	dStep3Wrapper = dComboWrapper->Get_ParticleComboStep(3);
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep3Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep3Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_pi0etapr__B4_M35_M7_M17_h
