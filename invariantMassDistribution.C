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

// Function to create invariant mass histogram
TCanvas* createInvariantMassHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

    TCanvas *c_InvariantMass = new TCanvas("c_InvariantMass", "Canvas for Invariant Mass Distribution", 800, 600);
    c_InvariantMass->cd();

    // Creating histogram and setting attributes
    TH1F* hist = new TH1F(legendLabel, "Invariant Mass for " + legendLabel, 410, 40, 450);
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

            // Initializing sum variables (for each jet)
            Float_t sumE = 0;
            Float_t sumPx = 0;
            Float_t sumPy = 0;
            Float_t sumPz = 0;

            // Loop through particles within a jet
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {

                // Checking if particle is charged
                if ((*pChg)[i][j] == 0) {continue;}

                Float_t px = calculatePx((*pPt)[i][j], (*pPhi)[i][j]);
                Float_t py = calculatePy((*pPt)[i][j], (*pPhi)[i][j]);
                Float_t pz = calculatePz((*pPt)[i][j], (*pEta)[i][j]);
                Float_t E = calculateEnergy(px, py, pz, particleData.m0((*pPid)[i][j]));

                sumPx += px;
                sumPy += py;
                sumPz += pz;
                sumE += E;
            }

            Float_t jetMass = sqrt((sumE * sumE) - (sumPx * sumPx) - (sumPy * sumPy) - (sumPz * sumPz));
            hist->Fill(jetMass);
        }
    }

    // ***** NORMALIZATION *****
    Float_t binWidth = 410.0/410.0;
    hist->Scale(1.0/(binWidth));

    hist->Draw("E");
    hist->GetYaxis()->SetTitle("#frac{dN_{entries}}{d#Deltam_{jet}}");
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetXaxis()->SetTitle("m_{jet} [GeV]");
    gPad->SetLogy();

    c_InvariantMass->Write();
    c_InvariantMass->Close();

    delete hist;
    return c_InvariantMass;    
}

void testTChain() {

    TChain *chain = new TChain("trackTree");
    chain->Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_*.root");
    //chain->Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1000.root");    
    //chain->Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1001.root");

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

    TFile *fout = new TFile("testMass.root", "recreate"); // Creating output file

    TCanvas* c_InvariantMass = createInvariantMassHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);

    delete c_InvariantMass;
    delete fout;
}