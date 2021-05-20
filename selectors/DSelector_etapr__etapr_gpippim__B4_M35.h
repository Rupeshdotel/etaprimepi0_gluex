#ifndef DSelector_etapr__etapr_gpippim__B4_M35_h
#define DSelector_etapr__etapr_gpippim__B4_M35_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_etapr__etapr_gpippim__B4_M35 : public DSelector
{
	public:

		DSelector_etapr__etapr_gpippim__B4_M35(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_etapr__etapr_gpippim__B4_M35(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		bool dIsMC;

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
		DNeutralParticleHypothesis* dPhotonWrapper;
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;

		//tree stuff (variables for analysis) goes here
		TFile *fileout; // file for the outputtree
		TTree *Tree_etapr_gpippim; //tree

		
		Double_t mm2m;
		
		Double_t be;
		Double_t mant;
		Double_t dt;

		Double_t mpippim;
		Double_t mgpippim;
		Double_t mpipp;
		Double_t mpimp;
		Double_t mpippimp;


		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		//TH1I* dHist_MissingMassSquared;
		//TH1I* dHist_BeamEnergy;

	ClassDef(DSelector_etapr__etapr_gpippim__B4_M35, 0);
};

void DSelector_etapr__etapr_gpippim__B4_M35::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhotonWrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(2));
}

#endif // DSelector_etapr__etapr_gpippim__B4_M35_h
