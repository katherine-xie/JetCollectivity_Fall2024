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

TCanvas* createPtHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

    TCanvas *c_Pt = new TCanvas("c_Pt", "Canvas for p_{T} Distribution", 800, 600);
    c_Pt->cd();

    // Creating histogram and setting attributes
    TH1F* hist = new TH1F(legendLabel, "p_{T} for " + legendLabel, 200, 0, 200);
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
                hist->Fill((*pPt)[i][j]);
            }
        }
    }
    
    // ***** NORMALIZATION *****
    Float_t binWidth = 200.0/200.0;
    hist->Scale(1.0/(binWidth));

    hist->GetYaxis()->SetTitle("#frac{dN_{entries}}{d#Deltap_{T}}");
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitle("p_{T}");

    hist->Draw("E");
    gPad->SetLogy();
    c_Pt->Write();
    c_Pt->Close();
    delete hist;
    return c_Pt;
}

void ptDistribution() {

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

    TFile *fout = new TFile("pt.root", "recreate"); // Creating output file

    TCanvas* c_Pt = createPtHist("pp (13 TeV, N_{ch} #geq 60)", kBlack, 21, reader);

    delete c_Pt;
    delete chain;
    delete fout;
}