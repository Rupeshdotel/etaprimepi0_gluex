//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  3 15:04:15 2021 by ROOT version 6.22/06
// from TTree kin/Kinematics
// found on file: genamp_gen_resonance.root
//////////////////////////////////////////////////////////

#ifndef Selector_gen_events_h
#define Selector_gen_events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "TH1D.h"

// Headers needed by this particular selector


class Selector_gen_events : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> NumFinalState = {fReader, "NumFinalState"};
   TTreeReaderArray<Float_t> E_FinalState = {fReader, "E_FinalState"};
   TTreeReaderArray<Float_t> Px_FinalState = {fReader, "Px_FinalState"};
   TTreeReaderArray<Float_t> Py_FinalState = {fReader, "Py_FinalState"};
   TTreeReaderArray<Float_t> Pz_FinalState = {fReader, "Pz_FinalState"};
   TTreeReaderValue<Float_t> E_Beam = {fReader, "E_Beam"};
   TTreeReaderValue<Float_t> Px_Beam = {fReader, "Px_Beam"};
   TTreeReaderValue<Float_t> Py_Beam = {fReader, "Py_Beam"};
   TTreeReaderValue<Float_t> Pz_Beam = {fReader, "Pz_Beam"};


   Selector_gen_events(TTree * /*tree*/ =0) { }
   virtual ~Selector_gen_events() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

  

   //declare a tree
   //tree stuff (variables for qfactor analysis) goes here
	TFile *fileout; // file for the outputtree
	TTree *qfactortree; //qfactortree

   //tree variables
   Double_t mproton;
   Double_t metap;
   Double_t mpi0;
   Double_t metaprimepi0;
   Double_t cos_theta;
   Double_t beam_energy;



   ClassDef(Selector_gen_events,0);

};

#endif

#ifdef Selector_gen_events_cxx
void Selector_gen_events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Selector_gen_events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Selector_gen_events_cxx
