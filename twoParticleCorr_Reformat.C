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
#include "TStopwatch.h"

#include <iostream>
#include <string>
#include <vector>

TStopwatch t;

// // Run locally
// #include "./SHARED_LIB_LOCAL/coordinateTools_SharedLib.h"
// #include "./SHARED_LIB_LOCAL/myFunctions_SharedLib.h"
// R__LOAD_LIBRARY(./SHARED_LIB_LOCAL/myFunctions_SharedLib_C.so);

// #include "./HEADER_FILES/coordinateTools.h"
//#include "./HEADER_FILES/myFunctions.h"

// Run on server:
#include "./SHARED_LIB_SERVER/COORDINATE_FUNCTIONS.h"
#include "./SHARED_LIB_SERVER/JETFRAME_FUNCTIONS.h"
R__LOAD_LIBRARY(./SHARED_LIB_SERVER/COORDINATE_FUNCTIONS_C.so);


// // Global variables
std::string title = "pp (13 TeV), Server Files";
TChain chain("trackTree");
TTreeReader reader(&chain);

// Setup branches for particles
TTreeReaderValue<std::vector<std::vector<Float_t>>> pPt(reader, "genDau_pt");
TTreeReaderValue<std::vector<std::vector<Float_t>>> pPhi(reader, "genDau_phi");
TTreeReaderValue<std::vector<std::vector<Float_t>>> pEta(reader, "genDau_eta");
TTreeReaderValue<std::vector<std::vector<int>>> pChg(reader, "genDau_chg");
TTreeReaderValue<std::vector<std::vector<int>>> pPid(reader, "genDau_pid");  

TTreeReaderValue<std::vector<Float_t>> jPt(reader, "genJetPt");
TTreeReaderValue<std::vector<Float_t>> jEta(reader, "genJetEta");
TTreeReaderValue<std::vector<Float_t>> jPhi(reader, "genJetPhi");
TTreeReaderValue<std::vector<int>> jMult(reader, "genJetChargedMultiplicity");

bool isInTrackBin(int currJetMult, int lowBound, int highBound) {
    if (currJetMult >= lowBound && currJetMult < highBound) {return true;} 
    else {return false;}
}

bool isInPtBin(float currPt, float lowBound, float highBound) {
    if (currPt >= lowBound && currPt < highBound) {return true;} 
    else {return false;}
}

// Jet eta < 1.6, Jet pt > 550 GeV
bool jetPass(float jetEta, float jetPt) {
    if (jetEta < 1.6 && jetPt > 550) {return true;}
    else {return false;}
}

void initializeChain() {

    // Local chain (inclusive)
    //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_1.root");
    //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_2.root");
    //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_3.root");
    //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_*.root");


    //Local chain
    //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1000.root");
    //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1001.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1002.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1003.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1004.root");
   //chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_*.root");


    // // Limited Server chain
    // for (int i = 0; i < 100; i++) {   
    //     std::string str = "/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_";
    //     str.append(std::to_string((i+1)));
    //     str.append(".root");
    //     chain.Add(str.c_str());
    // }

    // Server chain
    chain.Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_nChGT60_*.root");

    // Server chain (Inclusive)
    chain.Add("/storage1/users/aab9/Pythia8_CP5_PrivateGen_April27/pp_highMultGen_CP5_inclusive_*.root");

    TObjArray *fileList = chain.GetListOfFiles();

    int countFiles = 0;
    for (int i = 0; i < fileList->GetEntries(); i++) {
        countFiles++;
        TNamed *file = (TNamed*)fileList->At(i);  // Cast to TNamed
        std::cout << i+1 << ") File: " << file->GetTitle() << std::endl;  // Get file name
    }
    std::cout << "Total number of files: " << countFiles << std::endl;
    std::cout << "Total number of entries: " << chain.GetEntries() << std::endl;
}


void twoParticleCorr_Reformat() {

    std::cout << "------------------ Now calculating " << title << " ------------------" << std::endl;

    initializeChain();

    // Jet Multiplicity bin bounds
    const int trackBinLower[10] = {0,  20,  30,  40,  50,  59,  66,  76,  83,   78};
    const int trackBinUpper[10] = {20, 30,  40,  50,  59,  66,  76,  83,  1000, 1000};

    // jT bin bounds
    const float ptBinLower[4] = {0.0, 0.3,  0.5,  1.0};
    const float ptBinUpper[4] = {0.3, 3.0,  3.0,  3.0};

    // Declaring Histograms
    TH1D* hEventPass = new TH1D("hEventPass", "Events Passed Per Track Bin", 10, 0, 10);
    TH2D* hJetPass = new TH1D("hJetPass", "Jets Passed Per Track Bin and pT Bins", 10, 0, 10, 4, 0, 4);
    TH2D* hPairs = new TH2D("hPairs", "Pairs between Track Bins and pT Bins", 10, 0, 10, 4, 0, 4);

    // Histogram Arrays
    TH1D* hBinDist[10];
    TH1D* hJtDist[10][4]; // jT distribution for each track bin and pT bin
    TH2D* hSignal[10][4];
    TH2D* hSignal_Normalized[10][4];
    TH2D* hBackground[10][4];
    TH2D* hEtaPhi[10][4];

    // Counters (for verification)
    int numSelectedJets = 0;

    // Track bin loop
    for (int tBin = 0; tBin < 10; tBin++) {

        // Nch distribution (to later calculate <Nch>)
        std::string binDistName = "hBinDist" + std::to_string(tBin);
        std::string binDistTitle = "Multiplicity Distributon (Track Bin " + std::to_string(tBin) + ")";
        hBinDist[tBin] = new TH1D(binDistName.c_str(), binDistTitle.c_str(), 120, 0, 120);

        // pT bin loop
        for (int pBin = 0; pBin < 4; pBin++) {

            std::string jTBinDistName = "hJtDist" + std::to_string(tBin) + std::to_string(pBin);
            std::string jTBinDistTitle = "jT Distributon (Track Bin " + std::to_string(tBin) + ", pT Bin " + std::to_string(pBin) + ")";
            hJtDist[tBin][pBin] = new TH1D(jTBinDistName.c_str(), jTBinDistTitle.c_str(), 100, 0, 4);

            // Initializing signal hist
            std::string signalName = "hSignal" + std::to_string(tBin) + std::to_string(pBin);
            std::string signalTitle = "Signal Distribution (Track Bin " + std::to_string(tBin) + ", pT Bin " + std::to_string(pBin) + ")";
            hSignal[tBin][pBin] = new TH2D(signalName.c_str(), signalTitle.c_str(), 41, -4, 4, 33, -TMath::Pi(), 2*TMath::Pi());

            // Initializing background hist
            std::string backgroundName = "hBackground" + std::to_string(tBin) + std::to_string(pBin);
            std::string backgroundTitle = "Background Distribution (Track Bin " + std::to_string(tBin) + ", pT Bin " + std::to_string(pBin) + ")";
            hBackground[tBin][pBin] = new TH2D(backgroundName.c_str(), backgroundTitle.c_str(), 41, -4, 4, 33, -TMath::Pi(), 2*TMath::Pi());

            // Initializing eta-phi hist
            std::string etaPhiName = "hEtaPhi" + std::to_string(tBin) + std::to_string(pBin);
            std::string etaPhiTitle = "Eta-Phi Distribution (Track Bin " + std::to_string(tBin) + ", pT Bin " + std::to_string(pBin) + ")";
            hEtaPhi[tBin][pBin] = new TH2D(etaPhiName.c_str(), etaPhiTitle.c_str(), 150, 0, 10, 120, -4, 4);

        }
    }

    reader.Restart(); // Ensuring event loop starts from beginning

    // ***** EVENT LOOP *****
    std::cout << "================== Entering Event Loop ... ==================" << std::endl;

    while (reader.Next()) {

        int currEventIndex = reader.GetCurrentEntry();

            // Jet Loop
            for (int iJet = 0; iJet < pPt->size(); iJet++) {

                int numTracksinJet = (*pPt)[iJet].size(); // Number of particles in the ith jet
                int numTriggCurrJet = 0; // Counter for trigger particle sper jet

                if (!jetPass((*jEta)[iJet], (*jPt)[iJet])) {continue;} // Applying jet criteria 
                numSelectedJets++;

                // Track bin loop 
                for (int tBin = 0; tBin < 10; tBin++) {
                    if (!isInTrackBin((*jMult)[iJet], trackBinLower[tBin], trackBinUpper[tBin])) {continue;} 
                    hJetPass->Fill(tBin);
                    hBinDist[tBin]->Fill((*jMult)[iJet]);
                }

                TVector3 jet;
                double jetPt = (*jPt)[iJet];
                double jetEta = (*jEta)[iJet];
                double jetPhi = (*jPhi)[iJet];

                jet.SetPtEtaPhi(jetPt, jetEta, jetPhi); 

                // Particle loop (to find number of trigger particles) 
                for (int iParticle = 0; iParticle < numTracksinJet; iParticle++) {

                    // Particle selection criteria (lab frame)
                    // pT > 0.3, |eta| < 2.4
                    if ((*pChg)[iJet][iParticle] == 0) {continue;}
                    if((*pPt)[iJet][iParticle] <= 0.3) {continue;}
                    if(fabs((*pEta)[iJet][iParticle]) >= 2.4) continue;

                    numTriggCurrJet++;
                }

                // Particle loop
                for (int iTrack = 0; iTrack < numTracksinJet; iTrack++) {

                    // Particle selection criteria (lab frame)
                    // pT > 0.3, |eta| < 2.4
                    if ((*pChg)[iJet][iTrack] == 0) {continue;}
                    if((*pPt)[iJet][iTrack] <= 0.3) {continue;}
                    if(fabs((*pEta)[iJet][iTrack]) >= 2.4) continue;

                    TVector3 currParticle;

                    currParticle.SetPtEtaPhi((*pPt)[iJet][iTrack],
                                             (*pEta)[iJet][iTrack],
                                             (*pPhi)[iJet][iTrack]);


                    // Calculating jet frame values
                    double jetFramePt = ptWRTJet(jet, currParticle);
                    double jetFrameEta = etaWRTJet(jet, currParticle);     
                    double jetFramePhi = phiWRTJet(jet, currParticle);
                    
                    if (std::isnan(jetFramePt)) {continue;}
                    if (std::isnan(jetFrameEta)) {continue;}
                    if (std::isnan(jetFramePhi)) {continue;}


                    // FILLING ETA PHI HISTOGRAM HERE

                    // Track bin loop
                    for (int tBin = 0; tBin < 10; tBin++) {

                        if (!isInTrackBin((*jMult)[iJet], trackBinLower[tBin], trackBinUpper[tBin])) {continue;}
                        
                        // pT bin loop
                        for (int pBin = 0; pBin < 4; pBin++) {
                            
                            if (!isInPtBin(jetFramePt, ptBinLower[pBin], ptBinUpper[pBin])) {continue;}
                            
                            hJtDist[tBin][pBin]->Fill(jetFramePt);
                            hEtaPhi[tBin][pBin]->Fill(jetFrameEta, jetFramePhi);

                        } // pT bin loop end
                    } // Track bin loop end


                    if (iTrack == numTracksinJet - 1) {continue;} // avoiding dimension errors

                    // FILLING SIGNAL HISTOGRAM HERE
                    for (int jTrack = iTrack + 1; jTrack < numTracksinJet; jTrack++) {

                        // Particle selection criteria (lab frame)
                        // pT > 0.3, |eta| < 2.4
                        if ((*pChg)[iJet][iTrack] == 0) {continue;}
                        if((*pPt)[iJet][iTrack] <= 0.3) {continue;}
                        if(fabs((*pEta)[iJet][iTrack]) >= 2.4) {continue;}
                        
                        TVector3 assocParticle;

                        assocParticle.SetPtEtaPhi((*pPt)[iJet][jTrack],
                                                (*pEta)[iJet][jTrack],
                                                (*pPhi)[iJet][jTrack]);

                        // Calculating jet frame values
                        double jetFramePt_assoc = ptWRTJet(jet, assocParticle);
                        double jetFrameEta_assoc = etaWRTJet(jet, assocParticle);     
                        double jetFramePhi_assoc = phiWRTJet(jet, assocParticle);
                    
                        if (std::isnan(jetFramePt_assoc)) {continue;}
                        if (std::isnan(jetFrameEta_assoc)) {continue;}
                        if (std::isnan(jetFramePhi_assoc)) {continue;}
                        
                        // Track bin
                        for (int tBin = 0; tBin < 10; tBin++) {
                            
                            if (!isInTrackBin((*jMult)[iJet], trackBinLower[tBin], trackBinUpper[tBin])) {continue;}

                            // pT Bin
                            for (int pBin = 0; pBin < 4; pBin++) {
                                    
                                    // Both particles need to be in pT bin
                                    if (isInPtBin(jetFramePt, ptBinLower[pBin], ptBinUpper[pBin]) &&
                                        isInPtBin(jetFramePt_assoc, ptBinLower[pBin], ptBinUpper[pBin])) {
                                        
                                        double deltaEta = fabs((jetFrameEta - jetFrameEta_assoc));
                                        double deltaPhi = fabs(TMath::ACos(TMath::Cos(jetFramePhi - jetFramePhi_assoc)));

                                        //std::cout << "Event " << currEventIndex << ", Track " << jTrack << "/" << numTracksinJet << ", [" << tBin << "][" << pBin << "]: (" << deltaEta << ", " << deltaPhi << ")" << std::endl;

                                        hSignal[tBin][pBin]->Fill(deltaEta, deltaPhi, 1.0/numTriggCurrJet);
                                        hSignal[tBin][pBin]->Fill(-deltaEta, deltaPhi, 1.0/numTriggCurrJet);
                                        hSignal[tBin][pBin]->Fill(deltaEta, -deltaPh, 1.0/numTriggCurrJet);
                                        hSignal[tBin][pBin]->Fill(-deltaEta, -deltaPhi, 1.0/numTriggCurrJet);
                                        hSignal[tBin][pBin]->Fill(deltaEta, 2*TMath::Pi() - deltaPhi, 1.0/numTriggCurrJet);
                                        hSignal[tBin][pBin]->Fill(-deltaEta, 2*TMath::Pi() - deltaPhi, 1.0/numTriggCurrJet);
                                        hPairs->Fill(tBin, pBin);

                                    }
                            } // pT bin loop end
                        } // track bin loop end
                    } // associate particle loop end
                } // Main particle loop end
            } // Jet loop end
    } // Event loop end

    std::cout << "Total selected jets: " << numSelectedJets << std::endl;
    std::cout << "Eta-phi and signal distributions finished!" << std::endl;

    // Creating background distribution
    std::cout << "================== Creating Background Distribution ... ==================" << std::endl;
    int mixFactor = 5;

    for (int tBin = 0; tBin < 10; tBin++) {
        for (int pBin = 0; pBin < 4; pBin++) {

                long int numSigEntries =  hSignal[tBin][pBin]->GetEntries(); // should be 6 times as much (filled 6 times)
                long int numhPairsEntries = hPairs->GetBinContent(hPairs->FindBin(tBin, pBin));


                long int numPseudo = ((1+floor(sqrt(1+(4*2*mixFactor*numhPairsEntries))))/2);


                // double numPseudo = (1 + floor(sqrt(1+(4*2*mixFactor*(double)numSigEntries))))/2;; // Not sure why floor is added after the sqrt (used in Austin's code)
                
                // //numSigEntries = (int)numSigEntries;
                // numPseudo = (int)numPseudo; // Casting back into int

                std::cout << "Track Bin " << tBin << ", jT Bin " << pBin << ": " << std::endl;
                std::cout << "Num Entries in Signal: " << numSigEntries << std::endl;
                std::cout << "Num Entries in hPairs: " << numhPairsEntries << std::endl;
                std::cout << "Num Pseudoparticles: " << numPseudo << std::endl;
                std::cout << "" << std::endl;

                float etaArr[numPseudo];
                float phiArr[numPseudo];
                
                for(int i = 0; i < numPseudo; i++) {

                    gRandom->SetSeed(0);

                    // Making pseudoparticles
                    double eta, phi;
                    hEtaPhi[tBin][pBin]->GetRandom2(eta, phi);

                    etaArr[i] = eta;
                    phiArr[i] = phi;
                }
            
                for(long int i = 0; i < (numPseudo-1); i++){
                    for(long int j = (i+1); j < numPseudo; j++){

                        double deltaEta = fabs((etaArr[i]-etaArr[j]));
                        double deltaPhi = fabs((TMath::ACos(TMath::Cos(phiArr[i]-phiArr[j]))));

                        hBackground[tBin][pBin]->Fill(deltaEta, deltaPhi, 1);
                        hBackground[tBin][pBin]->Fill(-deltaEta, deltaPhi, 1);
                        hBackground[tBin][pBin]->Fill(deltaEta, -deltaPhi, 1);
                        hBackground[tBin][pBin]->Fill(-deltaEta, -deltaPhi,1);
                        hBackground[tBin][pBin]->Fill(deltaEta, 2*TMath::Pi() - deltaPhi, 1);
                        hBackground[tBin][pBin]->Fill(-deltaEta, 2*TMath::Pi() - deltaPhi, 1);

                    }
                }
        }
    }

    TFile *fout = new TFile("fullServerRun_NormalizedSignal.root", "recreate"); // Creating output file

    // Saving histograms to output file
    hJetPass->Write();
    hPairs->Write();

    for (int tBin = 0; tBin < 10; tBin++) {

        hBinDist[tBin]->Write();

        for (int pBin = 0; pBin < 4; pBin++) {
            hEtaPhi[tBin][pBin]->Write();
            hSignal[tBin][pBin]->Write();
            hBackground[tBin][pBin]->Write();
            hJtDist[tBin][pBin]->Write();
        }
    }


    fout->Close(); 
    t.Print();
}