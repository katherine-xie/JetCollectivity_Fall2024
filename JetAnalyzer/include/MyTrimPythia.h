//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 30 23:40:55 2024 by ROOT version 6.26/11
// from TTree trackTree/v1
// found on file: MATCH_mc_test.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>

using namespace std;

class MyClass {
public :
   TFile * fFile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::vector<std::string> fileList;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *genJetEta;
   vector<float>   *genJetPt;
   vector<float>   *genJetPhi;
   vector<int>     *genJetChargedMultiplicity;

   vector<vector<int> > *genDau_chg;
   vector<vector<int> > *genDau_pid;
   vector<vector<float> > *genDau_pt;
   vector<vector<float> > *genDau_eta;
   vector<vector<float> > *genDau_phi;

   // List of branches
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetChargedMultiplicity;   //!

   TBranch        *b_genDau_chg;   //!
   TBranch        *b_genDau_pid;   //!
   TBranch        *b_genDau_pt;   //!
   TBranch        *b_genDau_eta;   //!
   TBranch        *b_genDau_phi;   //!


   //MyClass(TTree *tree=0);
   MyClass(std::vector<std::string> _fileList);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int job, std::string fList);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef MyClass_cxx
MyClass::MyClass(std::vector<std::string> _fileList) : fChain(0)
//MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MATCH_mc_test.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("MATCH_mc_test.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("MATCH_mc_test.root:/analyzerOffline");
//      dir->GetObject("trackTree",tree);

//   }
   //Init(tree);
    fileList = _fileList;
}

MyClass::~MyClass()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genJetEta = 0;
   genJetPt = 0;
   genJetPhi = 0;
   genJetChargedMultiplicity = 0;

   genDau_chg = 0;
   genDau_pid = 0;
   genDau_pt = 0;
   genDau_eta = 0;
   genDau_phi = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetChargedMultiplicity", &genJetChargedMultiplicity, &b_genJetChargedMultiplicity);
   
   fChain->SetBranchAddress("genDau_chg", &genDau_chg, &b_genDau_chg);
   fChain->SetBranchAddress("genDau_pid", &genDau_pid, &b_genDau_pid);
   fChain->SetBranchAddress("genDau_pt", &genDau_pt, &b_genDau_pt);
   fChain->SetBranchAddress("genDau_eta", &genDau_eta, &b_genDau_eta);
   fChain->SetBranchAddress("genDau_phi", &genDau_phi, &b_genDau_phi);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
