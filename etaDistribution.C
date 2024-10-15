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

TCanvas* createEtaHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

    TCanvas *c_Eta = new TCanvas("c_Eta", "Canvas for #eta Distribution", 800, 600);
    c_Eta->cd();

    // Creating histogram and setting attributes
    TH1F* hist = new TH1F(legendLabel, "#eta for " + legendLabel, 100, -5, 5);
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
        // Loop through daughter branches (jets)
        for (Int_t i = 0; i < pPt->size(); i++) {
            // Loop through particles within a jet
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {
                // Checking if particle is charged
                if ((*pChg)[i][j] == 0) {continue;}
                hist->Fill((*pEta)[i][j]);
            }
        }
    }

    // ***** NORMALIZATION *****
    Float_t binWidth = 10.0/100.0;
    hist->Scale(1.0/(binWidth));

    hist->Draw("E");
    hist->GetYaxis()->SetTitle("#frac{dN_{entries}}{d#Delta#eta}");
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitle("#eta");

    hist->Draw("E");
    gPad->SetLogy();
    c_Eta->Write();
    c_Eta->Close();

    delete hist;
    return c_Eta;
}

void etaDistribution() {

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

    TFile *fout = new TFile("eta.root", "recreate"); // Creating output file

    TCanvas* c_Eta = createEtaHist("pp (13 TeV, N_{ch} #geq 60)", kBlack, 21, reader);

    delete c_Eta;
    delete chain;
    delete fout;
}