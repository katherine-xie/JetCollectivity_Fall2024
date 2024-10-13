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

TChain* createValidChain(const char* treeName, TChain* originalChain) {
    
    TChain* validChain = new TChain(treeName);
    
    // Iterate through all files in the original chain
    TObjArray* fileList = originalChain->GetListOfFiles();

    for (int i = 0; i < fileList->GetSize(); i++) {

        TString fileName = fileList->At(i)->GetTitle();
        TFile* file = TFile::Open(fileName, "read");
        
        // Check if the file is valid and can be opened
        if (file && !file->IsZombie()) {

            TTree* currTree = dynamic_cast<TTree*>(file->Get(treeName));

            // Checking if tree is valid (not a nullptr)
            if (currTree) {
                validChain->AddFile(fileName); // Add valid file to the new chain
            } 
            
            else {
                std::cout << "Tree not found in file: " << fileName << std::endl;
            }
            file->Close();

        } 
        else {
            std::cout << "Invalid file: " << fileName << std::endl;
        }
    }
    
    return validChain;
}


// Function to create eta histogram
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

// Function to create jet pT histogram
TCanvas* createJetFrameEtaHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

    TCanvas *c_JetFrameEta = new TCanvas("c_JetFrameEta", "Canvas for #eta* Distribution", 800, 600);
    c_JetFrameEta->cd();

    // Creating histogram and setting attributes
    TH1F* hist = new TH1F(legendLabel, "#eta* for " + legendLabel, 100, -1, 10);
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

    // Counter (for normalization purposes)
    Int_t countSelectedJets = 0;

    reader->Restart();

    // Event loop
    while (reader->Next()) {

        // Loop through daughter branches (jets)
        for (Int_t i = 0; i < pPt->size(); i++) {

            // Vector for each jet with components pT, eta, phi
            TVector3 jet;
            jet.SetPtEtaPhi(calculateJetPt((*pPt)[i], (*pPhi)[i]),  
                            calculateJetEta((*pPt)[i], (*pEta)[i], (*pPhi)[i]),
                            calculateJetPhi((*pPt)[i], (*pPhi)[i])); 

            countSelectedJets++;

            std::vector<TVector3> particlesVec; // Vector to hold the particles 

            // Loop through particles within a jet
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {

                // Checking if particle is charged
                if ((*pChg)[i][j] == 0) {continue;}

                // Initialize vector (for individual particles)
                TVector3 particle;
                    
                particle.SetPtEtaPhi((*pPt)[i][j],
                                     (*pEta)[i][j],
                                     (*pPhi)[i][j]);

                particlesVec.push_back(particle);
            }
            
            for (Int_t j = 0; j < particlesVec.size(); j++) {hist->Fill(etaWRTJet(jet, particlesVec[j]));} // Filling histogram
        }
    }

    // ***** NORMALIZATION *****
    Float_t binWidth = 11.0/100.0;
    hist->Scale(1.0/(countSelectedJets * binWidth));

    hist->Draw("E");
    hist->GetYaxis()->SetTitle("#frac{dN_{jet}^{entries}}{d#Delta#eta*}");
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetXaxis()->SetTitle("#eta*");
    gPad->SetLogy();

    c_JetFrameEta->Write();
    c_JetFrameEta->Close();

    delete hist;
    return c_JetFrameEta;    
}

void etaDistributions() {

    TChain *originalChain = new TChain("trackTree");
    originalChain->Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_*.root");
    //chain->Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1000.root");    
    //chain->Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1001.root");
    
    /*
    // Checking for null pointers
    if (currTree == nullptr) {
        std::cout << "There is an invalid tree (null pointer)" << std::endl;
        return;
    } */

    // Removing invalid tress
    TChain* newChain = createValidChain("trackTree", originalChain);
    delete originalChain;

    TObjArray *fileList = newChain->GetListOfFiles();

    Int_t countFiles = 0;
    for (Int_t i = 0; i < fileList->GetEntries(); i++) {
        countFiles++;
        TNamed *file = (TNamed*)fileList->At(i);  // Cast to TNamed
        std::cout << i+1 << ") File: " << file->GetTitle() << std::endl;  // Get file name
    }

    std::cout << "Total number of files: " << countFiles << std::endl;
    std::cout << "Total number of entries: " << newChain->GetEntries() << std::endl;

    TTreeReader *reader = new TTreeReader(newChain);

    TFile *fout = new TFile("testEtaHist.root", "recreate"); // Creating output file

    TCanvas* c_Eta = createEtaHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);
    TCanvas* c_JetFrameEta = createEtaHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);

    delete c_Eta;
    delete c_JetFrameEta;
    delete fout;
}