#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "THelix.h"
#include "TView.h"
#include "TRandom.h"
#include "TAttPad.h"
#include "TMath.h"
#include "TVector3.h"
#include "TView3D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLatex.h"
#include "TLegend.h"

#include <Pythia8/Pythia.h> 

#include <iostream>
#include <string>
#include <vector>

#include "SHARED_LIBRARY/COORDINATE_FUNCTIONS.h"
#include "SHARED_LIBRARY/JETFRAME_FUNCTIONS.h"

Pythia8::Pythia pythia; // Initialize Pythia
Pythia8::ParticleData &particleData = pythia.particleData; // Access the particle data table

// Function to create multiplicity histogram
TCanvas* createNchHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

    TCanvas *c_Nch = new TCanvas("c_Nch", "Canvas for N_{ch} Distribution", 800, 600);
    c_Nch->cd();

    // Creating histogram and setting attributes
    TH1F* hist = new TH1F(legendLabel, "N_{ch} for " + legendLabel, 140, 60, 200);
    hist->SetMarkerColor(colorVal);
    hist->SetMarkerStyle(markerStyle); // Pointer style
    hist->SetMarkerSize(0.5); // Pointer size
    hist->SetLineColor(kBlue); // Error bar color

    // Setup branches for particles
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(*reader, "genDau_pt");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(*reader, "genDau_phi");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(*reader, "genDau_eta");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> pChg(*reader, "genDau_chg");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> pPid(*reader, "genDau_pid");   

    reader->Restart();

    // Event loop
    while (reader->Next()) {

        Int_t nchCounter = 0;

        // Loop through daughter branches (jets)
        for (Int_t i = 0; i < pPt->size(); i++) {
            // Loop through particles within a jet
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {
                // Checking if particle is charged
                if ((*pChg)[i][j] == 0) {continue;}
                nchCounter++;
            }
        }

        if (nchCounter > 0) {hist->Fill(nchCounter);}

    }

    hist->GetYaxis()->SetTitle("N_{ch}^{entries}");
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetXaxis()->SetTitle("N_{ch}");

    hist->Draw("E");
    gPad->SetLogy();
    c_Nch->Write();
    c_Nch->Close();
    delete hist;
    return c_Nch;
}

void nchDistribution() {

    TChain *chain = new TChain("trackTree");
    chain->Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_*.root");

    TObjArray *fileList = chain->GetListOfFiles();

    Int_t countFiles = 0;
    for (Int_t i = 0; i < fileList->GetEntries(); i++) {
        countFiles++;
        TNamed *file = (TNamed*)fileList->At(i);  // Cast to TNamed
        std::cout << i+1 << ") File: " << file->GetTitle() << std::endl;  // Get file name
    }
    std::cout << "Total number of files: " << countFiles << std::endl;
    std::cout << "Total number of entries: " << chain->GetEntries() << std::endl;

    TTreeReader *reader = new TTreeReader(chain);

    // Setup branches for particles
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(*reader, "genDau_pt");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(*reader, "genDau_phi");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(*reader, "genDau_eta");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> pChg(*reader, "genDau_chg");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> pPid(*reader, "genDau_pid");  

    TFile *fout = new TFile("multiplicity.root", "recreate"); // Creating output file

    TCanvas* c_Nch = createNchHist("pp (13 TeV, N_{ch} #geq 60)", kBlack, 21, reader);

    delete c_Nch;
    delete chain;
    delete fout;
}