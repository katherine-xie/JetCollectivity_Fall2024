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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <string>
#include <vector>

TStopwatch t;

// //Run locally
// #include "./SHARED_LIB_LOCAL/coordinateTools_SharedLib.h"
// #include "./SHARED_LIB_LOCAL/myFunctions_SharedLib.h"
// R__LOAD_LIBRARY(./SHARED_LIB_LOCAL/myFunctions_SharedLib_C.so);

// #include "./HEADER_FILES/coordinateTools.h"
// #include "./HEADER_FILES/myFunctions.h"

// Run on server:
#include "./SHARED_LIB_SERVER/COORDINATE_FUNCTIONS.h"
#include "./SHARED_LIB_SERVER/JETFRAME_FUNCTIONS.h"
R__LOAD_LIBRARY(./SHARED_LIB_SERVER/COORDINATE_FUNCTIONS_C.so);


// // Global variables
std::string title = "pp (186 #leq N_{ch} #leq 227, 13 TeV, N_{ch})";
TChain chain("trackTree");
TTreeReader reader(&chain);

// Setup branches for particles
TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(reader, "genDau_pt");
TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(reader, "genDau_phi");
TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(reader, "genDau_eta");
TTreeReaderValue<std::vector<std::vector<Int_t>>> pChg(reader, "genDau_chg");
TTreeReaderValue<std::vector<std::vector<Int_t>>> pPid(reader, "genDau_pid");  

TTreeReaderValue<std::vector<Float_t>> jPt(reader, "genJetPt");
TTreeReaderValue<std::vector<Float_t>> jEta(reader, "genJetEta");
TTreeReaderValue<std::vector<Float_t>> jPhi(reader, "genJetPhi");


void initializeChain() {

    // // Local chain
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1000.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1001.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1002.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1003.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1004.root");

    // // Limited Server chain
    // for (Int_t i = 0; i < 4; i++) {   
    //     std::string str = "/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_";
    //     str.append(std::to_string((i+1)));
    //     str.append(".root");
    //     chain.Add(str.c_str());
    // }

    //Server chain
    chain.Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_*.root");

    TObjArray *fileList = chain.GetListOfFiles();

    Int_t countFiles = 0;
    for (Int_t i = 0; i < fileList->GetEntries(); i++) {
        countFiles++;
        TNamed *file = (TNamed*)fileList->At(i);  // Cast to TNamed
        std::cout << i+1 << ") File: " << file->GetTitle() << std::endl;  // Get file name
    }
    std::cout << "Total number of files: " << countFiles << std::endl;
    std::cout << "Total number of entries: " << chain.GetEntries() << std::endl;
}

// Function that returns the distribution of the simple signal function 
TH2F createSignalDist_JetFrame(std::vector<Int_t> multiplicityVector, 
                               std::vector<std::vector<std::vector<Float_t>>> jetEtaVals_AllEvents,
                               std::vector<std::vector<std::vector<Float_t>>> jetPhiVals_AllEvents) {

    std::cout << "------------------ Calculating Signal Distribution ... ------------------" << std::endl;
    
    std::string signalTitle = "Normalized Signal Distribution for " + title;

    // Histogram for the signal distribution
    TH2F hSignal("hSignal", signalTitle.c_str(), 50, -2, 2, 50, -TMath::Pi(), 2*TMath::Pi());

    // // Setup branches for particles
    // TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(reader, "genDau_pt");
    // TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(reader, "genDau_phi");
    // TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(reader, "genDau_eta");
    // TTreeReaderValue<std::vector<std::vector<Int_t>>> pChg(reader, "genDau_chg");
    // TTreeReaderValue<std::vector<std::vector<Int_t>>> pPid(reader, "genDau_pid");  

    Int_t numSelectedEvents = 0;
    Int_t jetCounter = 0; 
    Int_t numTrigg = 0;

    for (Int_t eventIndex = 0; eventIndex < jetEtaVals_AllEvents.size(); eventIndex++) {

        numSelectedEvents++;

        std::cout << "************** Event " << eventIndex << ", Num jets:  " << jPt->size() << " **************" << std::endl;

        // Looping through jets in the jetEtaVals vector
        for (Int_t jet = 0; jet < jetEtaVals_AllEvents[eventIndex].size(); jet++) {
            jetCounter++;

            std::cout << "Jet " << jet << ", Num particles: " << jetEtaVals_AllEvents[eventIndex][jet].size() << std::endl;
            
            // ***** PARTICLE LOOP ******
            // Particle 1 Loop
            for (Int_t i = 0; i < jetEtaVals_AllEvents[eventIndex][jet].size() - 1; i++) {

                // Caluclating eta and phi for particle 1
                Float_t eta1 = jetEtaVals_AllEvents[eventIndex][jet][i];
                Float_t phi1 = jetPhiVals_AllEvents[eventIndex][jet][i];
                numTrigg++;
                //std::cout << "Event " << eventIndex << ": Trigger Particle " << i << ", eta: " << eta1 << ", phi: " << phi1 << std::endl;

                // Particle 2 Loop
                for (Int_t j = i + 1; j < jetEtaVals_AllEvents[eventIndex][jet].size(); j++) {

                    // Calculating eta and phi for particle 2
                    Float_t eta2 = jetEtaVals_AllEvents[eventIndex][jet][j];
                    Float_t phi2 = jetPhiVals_AllEvents[eventIndex][jet][j];
                    // std::cout << "Event " << eventIndex << ": Trigg " << i << ", eta: " << eta1 << ", phi: " << phi1 << std::endl;
                    // std::cout<< "               Assoc " << j << ": eta " << eta2 << ", phi " << phi2 << std::endl;

                    // Calculating delta eta and delta phi
                    Float_t deltaEta = fabs(eta2 - eta1);
                    //std::cout << "delta eta: " << deltaEta << std::endl;

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
    }
    std::cout << "Number of trigger particles: " << numTrigg << std::endl;
    std::cout << "Number of selected events for the signal distribution: " << numSelectedEvents << std::endl;
    std::cout << "Number of selected jets: " << jetCounter << std::endl;

    // ***** NORMALIZATION ******
    //hSignal.Scale(1.0/(reader.GetEntries()));
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

    std::cout << "Done!" << std::endl;

    return hSignal;
}

// Function that creates the eta-phi distribution for the pseudoparticle mixing technique
TH2F createEtaPhiDist_JetFrame(std::vector<Int_t> multiplicityVector,
                               std::vector<std::vector<std::vector<Float_t>>> jetEtaVals_AllEvents,
                               std::vector<std::vector<std::vector<Float_t>>> jetPhiVals_AllEvents) {

    std::cout << "------------------Calculating Eta Phi Distribution ... ------------------" << std::endl;

    // First intialize the eta-phi distribution for all particles
    TH2F etaPhiDist("etaPhiDist", "(#eta*, #phi*) Distribution for all Particles", 50, 0, 7, 50, -TMath::Pi(), TMath::Pi());
    
    // // Setup branches for particles
    // TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(reader, "genDau_pt");
    // TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(reader, "genDau_phi");
    // TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(reader, "genDau_eta");
    // TTreeReaderValue<std::vector<std::vector<Int_t>>> pChg(reader, "genDau_chg");
    // TTreeReaderValue<std::vector<std::vector<Int_t>>> pPid(reader, "genDau_pid");  

    for (Int_t currEventIndex = 0; currEventIndex < jetEtaVals_AllEvents.size(); currEventIndex++) {

        // Jet Loop
        for (Int_t i = 0; i < jetEtaVals_AllEvents[currEventIndex].size(); i++) {
            // Particle Loop
            for (Int_t j = 0; j < jetEtaVals_AllEvents[currEventIndex][i].size(); j++) {
                etaPhiDist.Fill(jetEtaVals_AllEvents[currEventIndex][i][j], 
                                jetPhiVals_AllEvents[currEventIndex][i][j]);
            }
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

    std::cout << "Done!" << std::endl;

    return etaPhiDist;
}


// Function that returns the distribution of the simple mixed-event background function
TH2F createBackgroundDist_JetFrame(std::vector<Int_t> multiplicityVector, Int_t numMixFactor, TH2F hSignal, TH2F etaPhiDist) {

    std::cout << "------------------ Calculating Background Distribution ... ------------------" << std::endl;

    std::string backgroundTitle = "Background Distribution for " + title;

    // Histogram for the background distribution
    TH2F hBackground("hBackground", backgroundTitle.c_str(), 50, -2, 2, 50, -TMath::Pi(), 2*TMath::Pi());

    // ***** PSEUDOPARTICLE MIXING *****
    // Calculating the number of pseudoparticles samples:
    // Note: numPseudo needs to first be a double to avoid limits with the int datatype
    std::cout << "Values of entries directly from the histogram: " << hSignal.GetEntries() << std::endl;
    double numSigEntries = hSignal.GetEntries();
    double numPseudo = (1 + floor(sqrt(1+(4*2*numMixFactor*(double)numSigEntries))))/2;; // Not sure why floor is added after the sqrt (used in Austin's code)
    
    //numSigEntries = (int)numSigEntries;
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
            ///hBackground.Fill(deltaEta, deltaPhi);
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

    //hBackground.SetStats(0);

    std::cout << "Done!" << std::endl;

    return hBackground; 
} 


int twoParticleCorr_MultBin() {

    initializeChain();

    std::vector<Int_t> multVec; // Stores multiplicity values; vector has same size as num. events

    // The following vectors have dimensions of event{jet{particles}}
    std::vector<std::vector<std::vector<Float_t>>> jetEtaVals_AllEvents; // Stores eta* vals in jet frame
    std::vector<std::vector<std::vector<Float_t>>> jetPhiVals_AllEvents; // Stores phi* vals in jet frame

    reader.Restart(); // Ensuring event loop starts from beginning

    // ***** EVENT LOOP *****
    while (reader.Next()) {

        // ***** Finding the multiplicities of each event and storing it in a vector *****
        Int_t multiplicity = 0; // Counter for N_ch
        for (Int_t i = 0; i < pPt->size(); i++) {
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {
                if ((*pChg)[i][j] != 0) {multiplicity++;}
            }
        }
        multVec.push_back(multiplicity);
        //std::cout << "Index " << reader->GetCurrentEntry() << ": " << multiplicity << std::endl; 
    }

    reader.Restart();

    while (reader.Next()) {

        // SELF DEFINED JET FUNCTIONS ARE PROBLEMATIC
        // // ***** Finding the jet pT's and eta (of the whole jet) and storing it in a vector 
        // // This is for the jet selection criteria
        // std::vector<Float_t> ptOfJetVals_SingleEvent;
        // std::vector<Float_t> etaOfJetVals_SingleEvent;
        // // ***** JET LOOP *****
        // for (Int_t i = 0; i < pPt->size(); i++) {
        //    // std::cout << "i: " << i<< std::endl;
        //     // Calculating the pT and etas of the ith jet (NOTE: I am using ALL particles including the non-charged ones)
        //     Float_t currPtOfJet = calculateJetPt((*pPt)[i], (*pPhi)[i]);
        //     Float_t currEtaOfJet = calculateJetEta((*pPt)[i], (*pEta)[i], (*pPhi)[i]); 
        //     //std::cout << "Macro: " << currEtaOfJet << std::endl;
        //     ptOfJetVals_SingleEvent.push_back(currPtOfJet);
        //     etaOfJetVals_SingleEvent.push_back(currEtaOfJet);
        //     // std::cout << "Event " << reader.GetCurrentEntry() << ", Jet " << i << 
        //     // ": Pt of Jet = " << currPtOfJet << 
        //     // ", Eta of Jet = " << currEtaOfJet << std::endl;
        // }

        Int_t currEventIndex = reader.GetCurrentEntry();
        
        //Selecting events
        if (multVec[currEventIndex] == 0) {continue;}
        if (multVec[currEventIndex] < 186 || multVec[currEventIndex] > 227) {continue;}

        std::cout << "Selected Event " << reader.GetCurrentEntry() << ": " << multVec[currEventIndex] << std::endl; 

    
        // ***** Calculating the coordinates in the jet frame *****
        std::vector<std::vector<Float_t>> jetEtaVals_SingleEvent;
        std::vector<std::vector<Float_t>> jetPhiVals_SingleEvent;

        // Jet Loop
        for (Int_t i = 0; i < pPt->size(); i++) {

            //std::cout << "Jet size (before cut)" << i << ": " << (*pPt)[i].size() << std::endl;
            std::vector<Float_t> jetEtaVals_PerJet;
            std::vector<Float_t> jetPhiVals_PerJet;

            // if (std::isnan(ptOfJetVals_SingleEvent[i])) {continue;}
            // if (std::isnan(etaOfJetVals_SingleEvent[i])) {continue;}
            
            if (std::isnan((*jPt)[i])) {continue;}
            if (std::isnan((*jEta)[i])) {continue;}

            // Applying jet selection criteria
            if ((*jPt)[i] <= 550) {continue;}
            if (fabs((*jEta)[i]) >= 1.6) {continue;}

            // std::cout << "Event " << currEventIndex << ", Jet " << i << 
            // ": Pt of Jet = " << (*jPt)[i] << 
            // ", Eta of Jet = " << (*jEta)[i] << std::endl;

            // Vector for each jet with components pT, eta, phi
            TVector3 jet;

            Double_t ptOfJet = (*jPt)[i];
            Double_t etaOfJet = (*jEta)[i];
            Double_t phiOfJet = (*jPhi)[i];

            jet.SetPtEtaPhi(ptOfJet, etaOfJet, phiOfJet); 

            // jet.SetPtEtaPhi(calculateJetPt((*pPt)[i], (*pPhi)[i]),  
            //                 calculateJetEta((*pPt)[i], (*pEta)[i], (*pPhi)[i]),
            //                 calculateJetPhi((*pPt)[i], (*pPhi)[i])); 

            
            // Loop through particles within a jet
            for (Int_t j = 0; j < (*pPt)[i].size(); j++) {

                if ((*pChg)[i][j] == 0) {continue;}

                TVector3 currParticle;

                currParticle.SetPtEtaPhi((*pPt)[i][j],
                                         (*pEta)[i][j],
                                         (*pPhi)[i][j]);

                // Calculating the jT value and applying selection criteria
                Float_t currJetPt = ptWRTJet(jet, currParticle);
                if (std::isnan(currJetPt)) {continue;}
                if (currJetPt <= 0.3 || currJetPt >= 3) {continue;}

                Float_t currJetEta = etaWRTJet(jet, currParticle);     
                Float_t currJetPhi = phiWRTJet(jet, currParticle);
                
                if (std::isnan(currJetEta)) {continue;}
                if (std::isnan(currJetPhi)) {continue;}

                jetEtaVals_PerJet.push_back(currJetEta);
                jetPhiVals_PerJet.push_back(currJetPhi);

                // std::cout << "Event " << reader->GetCurrentEntry() << 
                // ": Eta = " << currJetEta <<
                // ", Phi = " << currJetPhi << 
                // ", jT = " << currJetPt << std::endl;
            }

            // Skip adding the jet if there are no particles that meet the criteria
            if (jetEtaVals_PerJet.size() == 0) {continue;}            
            if (jetPhiVals_PerJet.size() == 0) {continue;}

            jetEtaVals_SingleEvent.push_back(jetEtaVals_PerJet);
            jetPhiVals_SingleEvent.push_back(jetPhiVals_PerJet);
        }

        //if (jetEtaVals_SingleEvent.size() == 0) {continue;}            
        //if (jetPhiVals_SingleEvent.size() == 0) {continue;}
        jetEtaVals_AllEvents.push_back(jetEtaVals_SingleEvent);
        jetPhiVals_AllEvents.push_back(jetPhiVals_SingleEvent);
    }

    // for (Int_t i = 0; i < jetEtaVals_AllEvents.size(); i++) {
    //     std::cout << "******** Event " << i << "********" << std::endl;

    //     for (Int_t j = 0; j < jetEtaVals_AllEvents[i].size(); j++) {
    //         std::cout << "Jet " << j << ", Num particles: " << jetEtaVals_AllEvents[i][j].size() << std::endl;

    //         for (Int_t k = 0; k < jetEtaVals_AllEvents[i][j].size(); k++) {
    //             std::cout << "Eta*: " << jetEtaVals_AllEvents[i][j][k] << ", Phi*: " << jetPhiVals_AllEvents[i][j][k] << std::endl;
    //         }
    //     }
    // }

    TFile *fout = new TFile("Server_AllFiles_FourthBin.root", "recreate"); // Creating output file

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

    // cProjection->cd(2);
    // gPad->SetGrid();
    // TH1D *projectedBackgroundHist = simpleBackgroundHist.ProjectionY("projectedBackgroundHist", 1, -1);
    // projectedBackgroundHist->SetMinimum(0);
    // //projectedBackgroundHist->SetMaximum(50e6);
    // projectedBackgroundHist->SetLineWidth(2);
    // projectedBackgroundHist->SetLineColor(kBlue);
    // projectedBackgroundHist->Draw("HIST L SAME"); 

    cProjection->Write(); 

    delete cSignal;
    delete cEtaPhi;
    delete cBackground;
    delete cCorrected;
    delete cProjection;  

    fout->Close();
    t.Print();
    return 0;
}