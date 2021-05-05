#define Selector_gen_events_cxx
// The class definition in Selector_gen_events.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Selector_gen_events.C")
// root> T->Process("Selector_gen_events.C","some options")
// root> T->Process("Selector_gen_events.C+")
//


#include "Selector_gen_events.h"
#include "TLorentzVector.h"
#include <TH2.h>
#include <TStyle.h>

void Selector_gen_events::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();
   //hist = new TH1D("hist","Histogram of pt",100, -2.0, 2.0);

   //declare outputfile and gen_mctree
	fileout = new TFile("gen_mctree.root", "Recreate");
	gen_mctree = new TTree("gen_mctree", "gen_mctree");
   


  


   gen_mctree->Branch("mproton", &mproton, "mproton/D");
   gen_mctree->Branch("metap", &metap, "metap/D");
   gen_mctree->Branch("mpi0", &mpi0, "mpi0/D");
   gen_mctree->Branch("metaprimepi0", &metaprimepi0, "metaprimepi0/D");
   gen_mctree->Branch("cos_theta", &cos_theta, "cos_theta/D");
   gen_mctree->Branch("beam_energy", &beam_energy, "beam_energy/D");
   


}

void Selector_gen_events::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();

}

Bool_t Selector_gen_events::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   TLorentzVector  BeamP4;
   TLorentzVector  ProtonP4;
   TLorentzVector  EtaprimeP4;
   TLorentzVector  Pi0P4;
   TLorentzVector  EtaprimePi0P4;

   TVector3 boostGJ;

   TVector3 z_GJ;
   TVector3 z_hat_GJ;

   TVector3 y_GJ;
   TVector3 y_hat_GJ;

   TVector3 x_GJ;
   TVector3 x_hat_GJ;

   TVector3 vetaprime;


   
   
   // get the required 4 vectors
   BeamP4.SetPxPyPzE(*Px_Beam, *Py_Beam, *Pz_Beam, *E_Beam);

   ProtonP4.SetPxPyPzE(Px_FinalState[0], Py_FinalState[0], Pz_FinalState[0], E_FinalState[0]);

   EtaprimeP4.SetPxPyPzE(Px_FinalState[1],Py_FinalState[1],Pz_FinalState[1],E_FinalState[1]);

   Pi0P4.SetPxPyPzE(Px_FinalState[2],Py_FinalState[2],Pz_FinalState[2],E_FinalState[2]);
    

   //combine 4 vectors
   EtaprimePi0P4 = EtaprimeP4 + Pi0P4;

   // get the required components from 4 vectors
   beam_energy = BeamP4.E();
   mproton = ProtonP4.M();
   metap = EtaprimeP4.M();
   mpi0 = Pi0P4.M();
   metaprimepi0 = EtaprimePi0P4.M();

   // boost to GJ Frame
   boostGJ = -(EtaprimePi0P4.Vect())*(1.0/EtaprimePi0P4.E()); //calculate beta for etaprimepi0 system

      

   //redefine centre of mass 4-vectors to 4-vectors in rest frame of etaprimepi0 system
   TLorentzVector BeamP4_GJ = BeamP4;
   TLorentzVector EtaPrimeP4_GJ = EtaprimeP4;
   TLorentzVector Pi0P4_GJ = Pi0P4;
   TLorentzVector EtaPrimePi0P4_GJ = EtaprimePi0P4;
   
   //boost in rest frame of etaprimepi0 system
   BeamP4_GJ.Boost(boostGJ);
   EtaPrimeP4_GJ.Boost(boostGJ);
   Pi0P4_GJ.Boost(boostGJ);
	EtaPrimePi0P4_GJ.Boost(boostGJ);
	
	
	
	z_GJ.SetXYZ(BeamP4_GJ.X(), BeamP4_GJ.Y(), BeamP4_GJ.Z());//z GJ
	z_hat_GJ = z_GJ.Unit();  //beam direction (z-direction) 

   y_GJ = BeamP4.Vect().Cross(EtaprimePi0P4.Vect());  //y-direction for rest frame of etaprimepi0 system
   y_hat_GJ = y_GJ.Unit();    

   x_hat_GJ = y_hat_GJ.Cross(z_hat_GJ);//x hat GJ 

   vetaprime.SetXYZ(EtaPrimeP4_GJ.Vect()*x_hat_GJ, EtaPrimeP4_GJ.Vect()*y_hat_GJ, EtaPrimeP4_GJ.Vect()*z_hat_GJ);

   //get costhetat in GJ
   cos_theta = vetaprime.CosTheta();  


   // fill the generated MC tree
   gen_mctree->Fill();

   return kTRUE;
}

void Selector_gen_events::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Selector_gen_events::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
//hist->Draw();

//save the tree in outputfile 
fileout->cd(); 
gen_mctree->Write();

fileout->Close();

}