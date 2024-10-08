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
#include "TLatex.h"
#include "TLegend.h"

#include <iostream>
#include <string>
#include <vector>

#include "SHARED_LIBRARY/COORDINATE_FUNCTIONS.h"
#include "SHARED_LIBRARY/JETFRAME_FUNCTIONS.h"

// Test run for file pp_highMultGen_nChGT60_1065.root
// Global variables
std::string title = "pp (13 TeV), N_{ch} #geq 60";
std::string fullPathDir = "/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_1065.root"; // Path directory for the source data
//Int_t multCutOff = 194; // The pp data is already cut off at N_ch > 60
TFile *f = new TFile(fullPathDir.c_str(), "read"); // Opening file

// Creating TTreeReader object and linking branches
TTreeReader* reader = new TTreeReader("trackTree");

// Setup branches for particles
TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(*reader, "genDau_pt");
TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(*reader, "genDau_phi");
TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(*reader, "genDau_eta");
TTreeReaderValue<std::vector<std::vector<Int_t>>> pChg(*reader, "genDau_chg");

// Function that returns the distribution of the simple signal function 
TH2F createSignalDist_JetFrame(std::vector<Int_t> multiplicityVector, 
                               std::vector<std::vector<Float_t>> jetEtaVals_AllEvents,
                               std::vector<std::vector<Float_t>> jetPhiVals_AllEvents) {
    
    std::string signalTitle = "Normalized Signal Distribution for " + title;

    // Histogram for the signal distribution
    TH2F hSignal("hSignal", signalTitle.c_str(), 100, -5, 5, 100, -TMath::Pi(), 2*TMath::Pi());

    reader->Restart(); // Restarting event loop

    Int_t numSelectedEvents = 0;
    Int_t jetCounter = 0; 
    Int_t numTrigg = 0;

    while (reader->Next()) {

        // Check to see if the event is in the multiplicity bin
        Int_t eventIndex = reader->GetCurrentEntry();
        if (multiplicityVector[eventIndex] == 0) {continue;}
        numSelectedEvents++;
        
        // ***** PARTICLE LOOP ******

        // Particle 1 Loop
        for (Int_t i = 0; i < jetEtaVals_AllEvents[eventIndex].size()-1; i++) {

            // Caluclating eta and phi for particle 1
            Float_t eta1 = jetEtaVals_AllEvents[eventIndex][i];
            //if (fabs(eta1) > etaCutOff) {continue;} // Checking eta range
            if (std::isnan(eta1)) {continue;}
            Float_t phi1 = jetPhiVals_AllEvents[eventIndex][i];
            numTrigg++;
            //std::cout << "Event " << eventIndex << ": Trigger Particle " << i << ", eta: " << eta1 << ", phi: " << phi1 << std::endl;

            // Particle 2 Loop
            for (Int_t j = i + 1; j < jetEtaVals_AllEvents[eventIndex].size(); j++) {

                // Calculating eta and phi for particle 2
                Float_t eta2 = jetEtaVals_AllEvents[eventIndex][j];
                //if (fabs(eta2) > etaCutOff) {continue;} // Checking eta range
                Float_t phi2 = jetPhiVals_AllEvents[eventIndex][j];

                // Calculating delta eta and delta phi
                Float_t deltaEta = fabs(eta2 - eta1);
                Float_t deltaPhi = TMath::ACos(TMath::Cos(phi2-phi1));

                // Filling the histograms multiple times due to symmetries
                hSignal.Fill(deltaEta, deltaPhi);
                hSignal.Fill(deltaEta, -deltaPhi);
                hSignal.Fill(-deltaEta, deltaPhi);
                hSignal.Fill(-deltaEta, -deltaPhi);
                hSignal.Fill(-deltaEta, (2*TMath::Pi()) - deltaPhi);
                hSignal.Fill(deltaEta, (2*TMath::Pi()) - deltaPhi); 
            }
        } 
    
    }
    std::cout << "Number of selected events for the signal distribution: " << numSelectedEvents << std::endl;

    // ***** NORMALIZATION ******
    //hSignal.Scale(1.0/(reader->GetEntries()));
    hSignal.Scale(1.0/numTrigg); 

    // ***** HISTOGRAM CUSTOMIZATION ******
    hSignal.SetXTitle("#Delta#eta*");
    hSignal.SetYTitle("#Delta#phi*");
    hSignal.SetZTitle("S(#Delta#eta*, #Delta#phi*)");
    
    hSignal.SetTitleOffset(1.5, "X");
    hSignal.SetTitleOffset(1.5, "Y");
    hSignal.SetTitleOffset(1.3, "Z");

    hSignal.SetTitleFont(132, "T");
    hSignal.SetTitleFont(132, "XYZ");

    hSignal.SetLabelFont(132, "T");
    hSignal.SetLabelFont(132, "XYZ");

    hSignal.SetStats(0);

    return hSignal;
}



// Function that creates the eta-phi distribution for the pseudoparticle mixing technique
TH2F createEtaPhiDist_JetFrame(std::vector<Int_t> multiplicityVector,
                               std::vector<std::vector<Float_t>> jetEtaVals_AllEvents,
                               std::vector<std::vector<Float_t>> jetPhiVals_AllEvents) {

    // First intialize the eta-phi distribution for all particles
    TH2F etaPhiDist("etaPhiDist", "(#eta*, #phi*) Distribution for all Particles", 50, 0, 7, 50, -TMath::Pi(), TMath::Pi());
     
    reader->Restart(); // Restarting event loop

    while (reader->Next()) {

        // Check to see if the event is in the multiplicity bin
        Int_t currEventIndex = reader->GetCurrentEntry();
        //if (multiplicityVector[currEventIndex] < multCutOff) {continue;}

        for (Int_t i = 0; i < jetEtaVals_AllEvents[currEventIndex].size(); i++) {
            
            if (std::isnan(jetEtaVals_AllEvents[currEventIndex][i])) {continue;}
            //if (fabs(eta) > etaCutOff) {continue;} // Checking eta range
            etaPhiDist.Fill(jetEtaVals_AllEvents[currEventIndex][i], 
                            jetPhiVals_AllEvents[currEventIndex][i]);
        }
    }

    // ***** CUSTOMIZATION *****
    etaPhiDist.SetXTitle("#eta*");
    etaPhiDist.SetYTitle("#phi*");
    etaPhiDist.SetZTitle("Entries");
    etaPhiDist.SetTitleOffset(1.5, "X");
    etaPhiDist.SetTitleOffset(1.5, "Y");
    etaPhiDist.SetTitleOffset(1.2, "Z");

    etaPhiDist.SetTitleFont(132, "T");
    etaPhiDist.SetTitleFont(132, "X");
    etaPhiDist.SetTitleFont(132, "Y");
    etaPhiDist.SetTitleFont(132, "Z");
    
    return etaPhiDist;
}



// Function that returns the distribution of the simple mixed-event background function
TH2F createBackgroundDist_JetFrame(std::vector<Int_t> multiplicityVector, Int_t numMixFactor, TH2F hSignal, TH2F etaPhiDist) {

    std::string backgroundTitle = "Background Distribution for " + title;

    // Histogram for the background distribution
    TH2F hBackground("hBackground", backgroundTitle.c_str(), 100, -5, 5, 100, -TMath::Pi(), 2*TMath::Pi());

    // ***** PSEUDOPARTICLE MIXING *****
    // Calculating the number of pseudoparticles samples:
    // Note: numPseudo needs to first be a double to avoid limits with the int datatype
    int numSigEntries = hSignal.GetEntries();
    double numPseudo = (1 + floor(sqrt(1+(4*2*numMixFactor*(double)numSigEntries))))/2;; // Not sure why floor is added after the sqrt (used in Austin's code)
    numPseudo = (int)numPseudo; // Casting back into int

    std::cout << "Number of entries in signal distribution: " << numSigEntries << std::endl;
    std::cout << "Number of pseudoparticles to select: " << numPseudo << std::endl;

    for (Int_t i = 0; i < numPseudo-1; i++) {
        Double_t selectedEta1, selectedPhi1;
        etaPhiDist.GetRandom2(selectedEta1, selectedPhi1);
        
        for (Int_t j = i+1; j < numPseudo; j++) {
            Double_t selectedEta2, selectedPhi2;
            etaPhiDist.GetRandom2(selectedEta2, selectedPhi2);
            
            // Calculating delta eta and delta phi
            Float_t deltaEta = fabs(selectedEta2 - selectedEta1);
            Float_t deltaPhi = TMath::ACos(TMath::Cos(selectedPhi2-selectedPhi1));

            //std::cout << i << ": (" << deltaEta << ", " << deltaPhi << ")" << std::endl;

            // Filling the histogram multiple times due to symmetries
            hBackground.Fill(deltaEta, deltaPhi);
            hBackground.Fill(deltaEta, -deltaPhi);
            hBackground.Fill(-deltaEta, deltaPhi);
            hBackground.Fill(-deltaEta, -deltaPhi);
            hBackground.Fill(-deltaEta, (2*TMath::Pi()) - deltaPhi);
            hBackground.Fill(deltaEta, (2*TMath::Pi()) - deltaPhi);
        }
    } 

    // ***** HISTOGRAM CUSTOMIZATION *****
    hBackground.SetXTitle("#Delta#eta*");
    hBackground.SetYTitle("#Delta#phi*");
    hBackground.SetZTitle("B(#Delta#eta*, #Delta#phi*)");

    hBackground.SetTitleOffset(1.5, "X");

    hBackground.SetTitleOffset(1.3, "Z");

    hBackground.SetTitleFont(132, "T");
    hBackground.SetTitleFont(132, "XYZ");

    hBackground.SetLabelFont(132, "T");
    hBackground.SetLabelFont(132, "XYZ");

    hBackground.SetStats(0);

    return hBackground; 
} 


void twoParticleCorr_JetFrame() {

    std::vector<Int_t> multVec; // Stores multiplicity values; vector has same size as num. events
    std::vector<std::vector<Float_t>> jetEtaVals_AllEvents; // Stores rotated eta vals in jet frame
    std::vector<std::vector<Float_t>> jetPhiVals_AllEvents; // Stores rotated phi vals in jet frame

    reader->Restart(); // Ensuring event loop starts from beginning

    // ***** EVENT LOOP *****
    while (reader->Next()) {

        // ***** Finding the multiplicities of each event and storing it in a vector *****
        Int_t multiplicity = 0; // Counter for N_ch

        for (Int_t i = 0; i < pPt->size(); i++) {
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {
                if ((*pChg)[i][j] != 0) {multiplicity++;}
            }
        }

        multVec.push_back(multiplicity);
        std::cout << "Event Index " << reader->GetCurrentEntry() << ": " << multiplicity << std::endl; 

        // ***** Calculating the coordinates in the jet frame *****
        std::vector<Float_t> jetEtaVals_SingleEvent;
        std::vector<Float_t> jetPhiVals_SingleEvent;

        for (Int_t i = 0; i < pPt->size(); i++) {

            // Vector for each jet with components pT, eta, phi
            TVector3 jet;
                
            jet.SetPtEtaPhi(calculateJetPt((*pPt)[i], (*pPhi)[i]),  
                            calculateJetEta((*pPt)[i], (*pEta)[i], (*pPhi)[i]),
                            calculateJetPhi((*pPt)[i], (*pPhi)[i])); 

            /* std::cout << "Event " << reader->GetCurrentEntry() << ", Jet Index " << i << ": " << 
            jet.Eta() << ", " << 
            jet.Phi() << ", " << 
            jet.Pt() << std::endl; */

            // Loop through particles within a jet
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {

                if ((*pChg)[i][j] == 0) {continue;}

                TVector3 currParticle;

                currParticle.SetPtEtaPhi((*pPt)[i][j],
                                         (*pEta)[i][j],
                                         (*pPhi)[i][j]);

                Float_t currJetEta = etaWRTJet(jet, currParticle);     
                Float_t currJetPhi = phiWRTJet(jet, currParticle);

                jetEtaVals_SingleEvent.push_back(currJetEta);
                jetPhiVals_SingleEvent.push_back(currJetPhi);

                /* std::cout << "Event " << reader->GetCurrentEntry() << 
                ": Eta = " << currJetEta <<
                ", Phi = " << currJetPhi << std::endl; */
            } 

            /*
            std::cout << "Event " << reader->GetCurrentEntry() << 
            ", Jet " << i <<
            ": Number of Rotated Particles = " << jetEtaVals_SingleEvent.size() << std::endl; */

        }

        jetEtaVals_AllEvents.push_back(jetEtaVals_SingleEvent);
        jetPhiVals_AllEvents.push_back(jetPhiVals_SingleEvent); 

    }

    
    //std::cout << jetEtaVals_AllEvents.size() << std::endl;
    //std::cout << jetPhiVals_AllEvents.size() << std::endl;
    std::cout << "Multiplicity Vector: " << multVec.size() << std::endl;


    TFile *fout = new TFile("twoParticleCorr_pp_highMultGen_nChGT60_1065.root", "recreate"); // Creating output file


    // Creating canvas for the signal histogram
    TCanvas *cSignal = new TCanvas("cSignal", "Canvas for the Signal Distribution", 800, 600);

    TH2F simpleSignalHist = createSignalDist_JetFrame(multVec, jetEtaVals_AllEvents, jetPhiVals_AllEvents);

    cSignal->cd();
    simpleSignalHist.Draw("SURF1");
    cSignal->Write();


    // Testing eta/phi distribution
    TCanvas *cEtaPhi = new TCanvas("cEtaPhi", "Canvas for the Eta-Phi Distribution", 800, 600);
    TH2F etaPhiHist = createEtaPhiDist_JetFrame(multVec, jetEtaVals_AllEvents, jetPhiVals_AllEvents);
    
    cEtaPhi->cd();
    etaPhiHist.Draw("SURF1");
    cEtaPhi->Write();


    // Creating canvas for the background histogram
    TCanvas *cBackground = new TCanvas("cBackground", "Canvas for the Background Distribution", 800, 600);
    TH2F simpleBackgroundHist = createBackgroundDist_JetFrame(multVec, 10, simpleSignalHist, etaPhiHist);

    cBackground->cd();
    simpleBackgroundHist.Draw("SURF1");
    cBackground->Write();

    // Creating canvas for the corrected signal distribution
    TCanvas *cCorrected = new TCanvas("cCorrected", "Canvas for the Corrected Signal Distribution", 1000, 1000);
    TH2F correctedHist = simpleSignalHist;

    correctedHist.Divide(&simpleBackgroundHist);
    std::cout << "B(0,0): " << simpleBackgroundHist.GetBinContent(simpleBackgroundHist.FindBin(0,0)) << std::endl;
    correctedHist.Scale(simpleBackgroundHist.GetBinContent(simpleBackgroundHist.FindBin(0,0)));

    cCorrected->cd();
    cCorrected->SetFillColor(0);
    correctedHist.Draw("SURF1");

    // ***** HISTOGRAM CUSTOMIZATION ***** //
    std::string correctedTitle = "Corrected Signal Distribution for " + title;
    correctedHist.SetTitle(correctedTitle.c_str());
    //correctedHist.SetTitle("");

    correctedHist.GetZaxis()->SetTitle("C(#Delta#eta*, #Delta#phi*)");
    correctedHist.GetZaxis()->SetTitleSize(0.04);
    //correctedHist.GetZaxis()->SetTitleOffset(0.01);

    correctedHist.GetXaxis()->SetTitleOffset(1);
    correctedHist.GetXaxis()->SetTitleSize(0.05);

    correctedHist.GetYaxis()->SetTitleOffset(1);
    correctedHist.GetYaxis()->SetTitleSize(0.05);

    correctedHist.SetTitleOffset(1.1, "Z");
    correctedHist.SetTitleFont(132, "T");
    correctedHist.SetTitleFont(132, "XYZ");
    correctedHist.SetLabelFont(132, "T");
    correctedHist.SetLabelFont(132, "XYZ");

    correctedHist.SetAxisRange(-2, 2, "X");
    correctedHist.SetAxisRange(-2, 5, "Y");
    
    cCorrected->Write(); 

    // Copying the corrected histogram and truncating the z axis (for a better view)
    TCanvas *cTruncated = new TCanvas("cTruncated", "Canvas for the Truncated Corrected Distribution", 1000, 1000);
    cTruncated->cd();

    TH2F truncatedHist = correctedHist;
    truncatedHist.SetMaximum(0.11);
    truncatedHist.Draw("SURF1");
    cTruncated->Write();

    // Creating canvas for the projected delta phi histgram
    TCanvas *cProjection = new TCanvas("cProjection", "Canvas for the Projected Distributions", 800, 600);

    cProjection->cd();
    
    TH2F* correctedCopy = (TH2F*)correctedHist.Clone();
    correctedCopy->SetAxisRange(0, 2, "X");

    TH1D* projectedHist = correctedCopy->ProjectionY("projectedHist", 1, -1);

    projectedHist->SetLineWidth(2);
    projectedHist->SetLineColor(kBlue);
    projectedHist->Draw("HIST");
    projectedHist->GetXaxis()->SetTitleOffset(0.5);
    projectedHist->GetXaxis()->SetTitleFont(132);

    std::string projectionTitle = "Projection of Corrected Distribution for " + title;
    projectedHist->SetTitle(projectionTitle.c_str());

    gPad->SetGrid();
    
    gPad->SetGrid();
    TH1D *projectedSignalHist = simpleSignalHist.ProjectionY("projectedSignalHist", 1, -1);
    projectedSignalHist->SetLineWidth(2);
    projectedSignalHist->SetLineColor(kBlue);
    projectedSignalHist->Draw("HIST L");

    cProjection->cd(2);
    gPad->SetGrid();
    TH1D *projectedBackgroundHist = simpleBackgroundHist.ProjectionY("projectedBackgroundHist", 1, -1);
    projectedBackgroundHist->SetMinimum(0);
    //projectedBackgroundHist->SetMaximum(50e6);
    projectedBackgroundHist->SetLineWidth(2);
    projectedBackgroundHist->SetLineColor(kBlue);
    projectedBackgroundHist->Draw("HIST L SAME"); 

    cProjection->Write(); 
    projectedHist->Write(); 

    delete cSignal;
    delete cEtaPhi;
    delete cBackground;
    delete cCorrected;
    delete cTruncated;
    delete cProjection;  

    f->Close();
    fout->Close();
}