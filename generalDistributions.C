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

// Global variables
std::string title = "pp (13 TeV), N_{ch} #geq 60";
//std::string fullPathDir = "/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_1050.root"; // Path directory for the source data
//TFile *f = new TFile(fullPathDir.c_str(), "read"); // Opening file */

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
    //gPad->SetLogy();
    c_Nch->Write();
    c_Nch->Close();
    delete hist;
    return c_Nch;
}

// Function to create pT histogram
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
TCanvas* createJetFramePtHist(TString legendLabel, Int_t colorVal, Int_t markerStyle, TTreeReader* reader) {

    TCanvas *c_JetFramePt = new TCanvas("c_JetFramePt", "Canvas for p_{T}* Distribution", 800, 600);
    c_JetFramePt->cd();

    // Creating histogram and setting attributes
    TH1F* hist = new TH1F(legendLabel, "p_{T}* for " + legendLabel, 150, 0, 150);
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

    // Counter for normalization purposes
    Int_t countSelectedJets = 0;

    reader->Restart();

    // Event loop
    while (reader->Next()) {

        // Loop through daughter branches (jets)
        for (Int_t i = 0; i < pPt->size(); i++) {

            countSelectedJets++;

            // Vector for each jet with components pT, eta, phi
            TVector3 jet;
            jet.SetPtEtaPhi(calculateJetPt((*pPt)[i], (*pPhi)[i]),  
                            calculateJetEta((*pPt)[i], (*pEta)[i], (*pPhi)[i]),
                            calculateJetPhi((*pPt)[i], (*pPhi)[i])); 

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
            
            for (Int_t j = 0; j < particlesVec.size(); j++) {hist->Fill(ptWRTJet(jet, particlesVec[j]));} // Filling histogram
        }
    }
    
    // ***** NORMALIZATION *****
    Float_t binWidth = 150.0/150.0;
    hist->Scale(1.0/(countSelectedJets * binWidth));

    hist->Draw("E");
    hist->GetYaxis()->SetTitle("#frac{dN_{jet}^{entries}}{d#Deltap_{T}*}");
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetTitle("p_{T}*");

    gPad->SetLogy();
    c_JetFramePt->Write();
    c_JetFramePt->Close();
    delete hist;
    return c_JetFramePt;    
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

    // Counter (for normalization purposes)
    Int_t countSelectedJets = 0;

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

void generalDistributions() {

    // Chaining files together
    TChain *chain = new TChain("trackTree");
    chain->Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_*.root");

    std::cout << "Num entries in the chain: " << chain->GetEntries() << std::endl;
    
    // Setting up tree
    TTreeReader* reader = new TTreeReader(chain);

    TFile *fout = new TFile("generalHistograms.root", "recreate"); // Creating output file

    TCanvas* c_Nch = createNchHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);
    TCanvas* c_Pt = createPtHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);
    TCanvas* c_Eta = createEtaHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);
    TCanvas* c_JetFramePt = createJetFramePtHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);
    TCanvas* c_JetFrameEta = createJetFrameEtaHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);
    TCanvas* c_InvariantMass = createInvariantMassHist("pp (13 Tev, N_{ch} #geq 60)", kBlack, 21, reader);

    delete fout;
    delete reader;
    delete chain;
}