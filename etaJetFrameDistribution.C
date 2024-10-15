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

TCanvas* createEtaJetFrameHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

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

void etaJetFrameDistribution() {

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

    TFile *fout = new TFile("etaJetFrame.root", "recreate"); // Creating output file

    TCanvas* c_etaJetFrame = createEtaJetFrameHist("pp (13 TeV, N_{ch} #geq 60)", kBlack, 21, reader);

    delete chain;
    delete fout;
    delete c_etaJetFrame;
}