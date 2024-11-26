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
std::string title = "pp (13 TeV), All Server Files ";
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
TTreeReaderValue<std::vector<Int_t>> jMult(reader, "genJetChargedMultiplicity");

bool isInMultBin(Int_t currJetMult, Int_t lowBound, Int_t highBound) {
    if (currJetMult >= lowBound && currJetMult < highBound) {return true;} 
    else {return false;}
}

void initializeChain() {

    // // Local chain (inclusive)
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_1.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_2.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_CP5_inclusive_3.root");

    // //Local chain
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1000.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1001.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1002.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1003.root");
    // chain.Add("/Users/katherinexie/JetCollectivity_Fall2024/Pythia_CP5_SourceData/pp_highMultGen_nChGT60_1004.root");

    // // Limited Server chain
    // for (Int_t i = 0; i < 100; i++) {   
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
TH2D createSignalDist_JetFrame(std::vector<Int_t> multVec, 
                               std::vector<std::vector<std::vector<Double_t>>> jetEtaVals_AllEvents,
                               std::vector<std::vector<std::vector<Double_t>>> jetPhiVals_AllEvents,
                               Int_t arrIndex) {

    std::cout << "------------------ Calculating Signal Distribution ... ------------------" << std::endl;

    std::string objName = "hSignal" + std::to_string(arrIndex);
    std::string signalTitle = "Normalized Signal Distribution for " + title + ", Mult. Bin " + std::to_string(arrIndex);

    // Histogram for the signal distribution
    TH2D hSignal(objName.c_str(), signalTitle.c_str(), 50, -4, 4, 50, -TMath::Pi(), 2*TMath::Pi());

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

        if (multVec[eventIndex] == 0) {continue;}
        numSelectedEvents++;

        // Checking to see if there are any jets in the vector
        if (jetEtaVals_AllEvents[eventIndex].size() == 0) {continue;}            

        //std::cout << "************** Event " << eventIndex << ", Num jets:  " << jetEtaVals_AllEvents[eventIndex].size() << " **************" << std::endl;


        // Looping through jets in the jetEtaVals vector
        for (Int_t jet = 0; jet < jetEtaVals_AllEvents[eventIndex].size(); jet++) {

            jetCounter++;

            //std::cout << "Jet " << jet << ", Num particles: " << jetEtaVals_AllEvents[eventIndex][jet].size() << std::endl;
            
            // ***** PARTICLE LOOP ******
            // Particle 1 Loop
            for (Int_t i = 0; i < jetEtaVals_AllEvents[eventIndex][jet].size() - 1; i++) {

                // Caluclating eta and phi for particle 1
                Double_t eta1 = jetEtaVals_AllEvents[eventIndex][jet][i];
                Double_t phi1 = jetPhiVals_AllEvents[eventIndex][jet][i];
                numTrigg++;
                //std::cout << "Event " << eventIndex << ": Trigger Particle " << i << ", eta: " << eta1 << ", phi: " << phi1 << std::endl;

                // Particle 2 Loop
                for (Int_t j = i + 1; j < jetEtaVals_AllEvents[eventIndex][jet].size(); j++) {

                    // Calculating eta and phi for particle 2
                    Double_t eta2 = jetEtaVals_AllEvents[eventIndex][jet][j];
                    Double_t phi2 = jetPhiVals_AllEvents[eventIndex][jet][j];
                    // std::cout << "Event " << eventIndex << ": Trigg " << i << ", eta: " << eta1 << ", phi: " << phi1 << std::endl;
                    // std::cout<< "               Assoc " << j << ": eta " << eta2 << ", phi " << phi2 << std::endl;

                    // Calculating delta eta and delta phi
                    Double_t deltaEta = fabs(eta2 - eta1);
                    //std::cout << "delta eta: " << deltaEta << std::endl;

                    Double_t deltaPhi = TMath::ACos(TMath::Cos(phi2-phi1));

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
TH2D createEtaPhiDist_JetFrame(std::vector<Int_t> multiplicityVector,
                               std::vector<std::vector<std::vector<Double_t>>> jetEtaVals_AllEvents,
                               std::vector<std::vector<std::vector<Double_t>>> jetPhiVals_AllEvents, 
                               Int_t arrIndex) {

    std::cout << "------------------ Calculating Eta Phi Distribution ... ------------------" << std::endl;

    // First intialize the eta-phi distribution for all particles
    std::string objName = "etaPhiDist" + std::to_string(arrIndex);
    std::string etaPhiTitle = "(#eta*, #phi*) Distribution, Mult. Bin " + std::to_string(arrIndex);

    TH2D etaPhiDist(objName.c_str(), etaPhiTitle.c_str(), 50, 0, 7, 50, -TMath::Pi(), TMath::Pi());
    
    // // Setup branches for particles
    // TTreeReaderValue<std::vector<std::vector<Double_t>>> pPt(reader, "genDau_pt");
    // TTreeReaderValue<std::vector<std::vector<Double_t>>> pPhi(reader, "genDau_phi");
    // TTreeReaderValue<std::vector<std::vector<Double_t>>> pEta(reader, "genDau_eta");
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
TH2D createBackgroundDist_JetFrame(std::vector<Int_t> multiplicityVector, Int_t numMixFactor, TH2D hSignal, TH2D etaPhiDist, Int_t arrIndex) {

    std::cout << "------------------ Calculating Background Distribution ... ------------------" << std::endl;

    std::string objName = "hBackground" + std::to_string(arrIndex);
    std::string backgroundTitle = "Background Distribution for " + title + ", Mult. Bin " + std::to_string(arrIndex);

    // Histogram for the background distribution
    TH2D hBackground(objName.c_str(), backgroundTitle.c_str(), 50, -4, 4, 50, -TMath::Pi(), 2*TMath::Pi());

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
            Double_t deltaEta = fabs(selectedEta2 - selectedEta1);
            Double_t deltaPhi = TMath::ACos(TMath::Cos(selectedPhi2-selectedPhi1));

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

    std::cout << "Number of entries in background: " << hBackground.GetEntries() << std::endl;
    std::cout << "Done!" << std::endl;

    return hBackground; 
} 


int twoParticleCorr_MultBin() {

    std::cout << "------------------ Now calculating " << title << " ------------------" << std::endl;

    initializeChain();

    std::vector<Int_t> multVec; // Stores multiplicity values; vector has same size as num. events

    // The following vectors have dimensions of event{jet{particles}}
    std::vector<std::vector<std::vector<Double_t>>> jetEtaVals_AllEvents[10]; // Stores eta* vals in jet frame
    std::vector<std::vector<std::vector<Double_t>>> jetPhiVals_AllEvents[10]; // Stores phi* vals in jet frame

    // Jet Multiplicity bin bounds
    Int_t jMultLower[10] = {0,  20,  30,  40,  50,  59,  66,  76,  83,   78};
    Int_t jMultUpper[10] = {20, 30,  40,  50,  59,  66,  76,  83,  1000, 1000};

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

    //  ***** EVENT LOOP *****
    reader.Restart();
    while (reader.Next()) {

        // SELF DEFINED JET FUNCTIONS ARE PROBLEMATIC
        // // ***** Finding the jet pT's and eta (of the whole jet) and storing it in a vector 
        // // This is for the jet selection criteria
        // std::vector<Double_t> ptOfJetVals_SingleEvent;
        // std::vector<Double_t> etaOfJetVals_SingleEvent;
        // // ***** JET LOOP *****
        // for (Int_t i = 0; i < pPt->size(); i++) {
        //    // std::cout << "i: " << i<< std::endl;
        //     // Calculating the pT and etas of the ith jet (NOTE: I am using ALL particles including the non-charged ones)
        //     Double_t currPtOfJet = calculateJetPt((*pPt)[i], (*pPhi)[i]);
        //     Double_t currEtaOfJet = calculateJetEta((*pPt)[i], (*pEta)[i], (*pPhi)[i]); 
        //     //std::cout << "Macro: " << currEtaOfJet << std::endl;
        //     ptOfJetVals_SingleEvent.push_back(currPtOfJet);
        //     etaOfJetVals_SingleEvent.push_back(currEtaOfJet);
        //     // std::cout << "Event " << reader.GetCurrentEntry() << ", Jet " << i << 
        //     // ": Pt of Jet = " << currPtOfJet << 
        //     // ", Eta of Jet = " << currEtaOfJet << std::endl;
        // }

        Int_t currEventIndex = reader.GetCurrentEntry();
        
        //Selecting events
        //if (multVec[currEventIndex] == 0) {continue;}
        //if (multVec[currEventIndex] < 186 || multVec[currEventIndex] > 227) {continue;}
        //if (multVec[currEventIndex] < 94) {continue;}

        //std::cout << "Selected Event " << reader.GetCurrentEntry() << ": " << multVec[currEventIndex] << std::endl; 
  
        // Array Loop 
        for (Int_t arrIndex = 0; arrIndex < 10; arrIndex++) {
              
            // ***** Calculating the coordinates in the jet frame for each multipliicity bin*****
            std::vector<std::vector<Double_t>> jetEtaVals_SingleEvent;
            std::vector<std::vector<Double_t>> jetPhiVals_SingleEvent;

            // Jet Loop
            for (Int_t i = 0; i < pPt->size(); i++) {

                //std::cout << "Jet size (before cut)" << i << ": " << (*pPt)[i].size() << std::endl;
                std::vector<Double_t> jetEtaVals_PerJet;
                std::vector<Double_t> jetPhiVals_PerJet;
                    
                if (std::isnan((*jPt)[i])) {continue;}
                if (std::isnan((*jEta)[i])) {continue;}

                // Applying jet selection criteria
                //if ((*jMult)[i] < 100) {continue;}
                if ((*jPt)[i] <= 550) {continue;}
                if (fabs((*jEta)[i]) >= 1.6) {continue;}

                // Applying jet pass for each of the bins
                if (!isInMultBin((*jMult)[i], jMultLower[arrIndex], jMultUpper[arrIndex])) {continue;}

                //std::cout << currEventIndex << ", Array " << arrIndex << ", Jet " << i << ", Jet Mult: " << (*jMult)[i] << std::endl;

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
                    Double_t currJetPt = ptWRTJet(jet, currParticle);
                    if (std::isnan(currJetPt)) {continue;}
                    if (currJetPt <= 0.3 || currJetPt >= 3) {continue;}

                    Double_t currJetEta = etaWRTJet(jet, currParticle);     
                    Double_t currJetPhi = phiWRTJet(jet, currParticle);
                        
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
            jetEtaVals_AllEvents[arrIndex].push_back(jetEtaVals_SingleEvent);
            jetPhiVals_AllEvents[arrIndex].push_back(jetPhiVals_SingleEvent);
        } // Array loop end

        std::cout << "Event " << currEventIndex << " done" << std::endl;
    } // Event Loop End


    // // Testing to see if dimensions match up
    // std::cout << "jMult size (should be total number of events): " << jMult->size() << std::endl; 
    // std::cout << "multVec size (should be total number of events): " << multVec.size() << std::endl; 
    // std::cout << "jetEtaVals_AllEvents size: " << jetEtaVals_AllEvents.size() << std::endl;

    // for (Int_t i = 0; i < jetEtaVals_AllEvents.size(); i++) {
    //     std::cout << "******** Event " << i << "********" << std::endl;

    //     for (Int_t j = 0; j < jetEtaVals_AllEvents[i].size(); j++) {
    //         std::cout << "Jet " << j << ", Num particles: " << jetEtaVals_AllEvents[i][j].size() << std::endl;

    //         for (Int_t k = 0; k < jetEtaVals_AllEvents[i][j].size(); k++) {
    //             std::cout << "Eta*: " << jetEtaVals_AllEvents[i][j][k] << ", Phi*: " << jetPhiVals_AllEvents[i][j][k] << std::endl;
    //         }
    //     }
    // }

    TFile *fout = new TFile("Server_AllFiles.root", "recreate"); // Creating output file

    //Testing if the jetEtaVals array works
    std::cout << "jetEtaVals_AllEvents[0] size: " << jetEtaVals_AllEvents[0].size() << std::endl;
    std::cout << "jetEtaVals_AllEvents[1] size: " << jetEtaVals_AllEvents[1].size() << std::endl;
    std::cout << "jetEtaVals_AllEvents[6] size: " << jetEtaVals_AllEvents[6].size() << std::endl;
    std::cout << "jetEtaVals_AllEvents[9] size: " << jetEtaVals_AllEvents[9].size() << std::endl;

    // for (Int_t event = 0; event < jetEtaVals_AllEvents[8].size(); event++) {
    //     if (jetEtaVals_AllEvents[8][event].size() == 0) {continue;}
    //     std::cout << "Event " << event << ": " << jetEtaVals_AllEvents[8][event].size() << " jets" << std::endl;
    // }

    // for (Int_t j = 0; j < jetEtaVals_AllEvents[7][4194194].size(); j++) {
    //     std::cout << "Num tracks in jet " << j << ": " << jetEtaVals_AllEvents[7][4194194][j].size() << std::endl;
    // }

    // Creating array of histograms
    TH2D hSignalArr[10]; 
    TH2D etaPhiArr[10]; 
    TH2D hBackgroundArr[10]; 
    //TH2D* hCorrectedArr[10]; 

    for (Int_t arrIndex = 0; arrIndex < 10; arrIndex++) {

        std::cout << "Calculating array " << arrIndex << " ... ";

        // Creating canvas for the signal histogram
        //TCanvas *cSignal = new TCanvas("cSignal", "Canvas for the Signal Distribution", 800, 600);
        hSignalArr[arrIndex] = createSignalDist_JetFrame(multVec, jetEtaVals_AllEvents[arrIndex], jetPhiVals_AllEvents[arrIndex], arrIndex);
        hSignalArr[arrIndex].Write();

        // // cSignal->cd();
        // // simpleSignalHist.Draw("SURF1");
        // cSignal->Write();

        // Testing eta/phi distribution
        //TCanvas *cEtaPhi = new TCanvas("cEtaPhi", "Canvas for the Eta-Phi Distribution", 800, 600);
        etaPhiArr[arrIndex] = createEtaPhiDist_JetFrame(multVec, jetEtaVals_AllEvents[arrIndex], jetPhiVals_AllEvents[arrIndex], arrIndex);
        etaPhiArr[arrIndex].Write();
        // cEtaPhi->cd();
        // etaPhiHist.Draw("SURF1");
        // cEtaPhi->Write();

        // Creating canvas for the background histogram
        //TCanvas *cBackground = new TCanvas("cBackground", "Canvas for the Background Distribution", 800, 600);
        hBackgroundArr[arrIndex] = createBackgroundDist_JetFrame(multVec, 5, hSignalArr[arrIndex], etaPhiArr[arrIndex], arrIndex);
        hBackgroundArr[arrIndex].Write();
        // cBackground->cd();
        // simpleBackgroundHist.Draw("SURF1");
        // cBackground->Write();

        // Creating canvas for the corrected signal distribution
        //TCanvas *cCorrected = new TCanvas("cCorrected", "Canvas for the Corrected Signal Distribution", 1000, 1000);
        // //TH2D hCorrectedArr[arrIndex] = simpleSignalHist;
        // hCorrectedArr[arrIndex] = (TH2D*)hSignalArr[arrIndex].Clone();
        // hCorrectedArr[arrIndex]->Divide(&hBackgroundArr[arrIndex]);
        // std::cout << "B(0,0): " << hBackgroundArr[arrIndex].GetBinContent(hBackgroundArr[arrIndex].FindBin(0,0)) << std::endl;
        // hCorrectedArr[arrIndex]->Scale(hBackgroundArr[arrIndex].GetBinContent(hBackgroundArr[arrIndex].FindBin(0,0)));

        // // ***** CORRECTED HISTOGRAM CUSTOMIZATION ***** //
        // std::string correctedName = "hCorrected" + std::to_string(arrIndex);
        // std::string correctedTitle = "Corrected Signal Distribution for " + title + ", Mult. Bin " + std::to_string(arrIndex);
        // hCorrectedArr[arrIndex]->SetNameTitle(correctedName.c_str(), correctedTitle.c_str());

        // hCorrectedArr[arrIndex]->GetZaxis()->SetTitle("C(#Delta#eta*, #Delta#phi*)");
        // hCorrectedArr[arrIndex]->GetZaxis()->SetTitleSize(0.04);
        // hCorrectedArr[arrIndex]->GetZaxis()->SetTitleOffset(0.01);

        // hCorrectedArr[arrIndex]->GetXaxis()->SetTitleOffset(1);
        // hCorrectedArr[arrIndex]->GetXaxis()->SetTitleSize(0.05);

        // hCorrectedArr[arrIndex]->GetYaxis()->SetTitleOffset(1);
        // hCorrectedArr[arrIndex]->GetYaxis()->SetTitleSize(0.05);

        // hCorrectedArr[arrIndex]->SetTitleOffset(1.1, "Z");
        // hCorrectedArr[arrIndex]->SetTitleFont(132, "T");
        // hCorrectedArr[arrIndex]->SetTitleFont(132, "XYZ");
        // hCorrectedArr[arrIndex]->SetLabelFont(132, "T");
        // hCorrectedArr[arrIndex]->SetLabelFont(132, "XYZ");

        // hCorrectedArr[arrIndex]->SetAxisRange(-4, 4, "X");
        // hCorrectedArr[arrIndex]->SetAxisRange(-2, 5, "Y");
        // hCorrectedArr[arrIndex]->Write();

        // Creating canvas for the projected delta phi histgram
        //TCanvas *cProjection = new TCanvas("cProjection", "Canvas for the Projected Distributions", 800, 600);
        //cProjection->cd();
        // TH2D* correctedCopy = (TH2D*)hCorrectedArr[arrIndex]->Clone();
        // correctedCopy->SetAxisRange(0, 2, "X");

        // TH1D* projectedHist = correctedCopy->ProjectionY("projectedHist", 1, -1);
        // projectedHist->SetLineWidth(2);
        // projectedHist->SetLineColor(kBlue);
        // projectedHist->Draw("HIST");
        // projectedHist->GetXaxis()->SetTitleOffset(0.5);
        // projectedHist->GetXaxis()->SetTitleFont(132);

        // std::string projectionTitle = "Projection of Corrected Distribution for " + title;
        // projectedHist->SetTitle(projectionTitle.c_str());
        // gPad->SetGrid();
        // //cProjection->Write();
        // projectedHist->Write();

        // TH1D *projectedSignalHist = simpleSignalHist.ProjectionY("projectedSignalHist", 1, -1);
        // projectedSignalHist->SetLineWidth(2);
        // projectedSignalHist->SetLineColor(kBlue);
        // projectedSignalHist->Draw("HIST");
        // cProjection->cd(2);
        // gPad->SetGrid();

        // TH1D *projectedBackgroundHist = simpleBackgroundHist.ProjectionY("projectedBackgroundHist", 1, -1);
        // projectedBackgroundHist->SetMinimum(0);
        // //projectedBackgroundHist->SetMaximum(50e6);
        // projectedBackgroundHist->SetLineWidth(2);
        // projectedBackgroundHist->SetLineColor(kBlue);
        // //projectedBackgroundHist->Draw("HIST L SAME"); 
        // projectedBackgroundHist->Write();

        std::cout << "Done!" << std::endl;
    }

   

    // delete cSignal;
    // delete cEtaPhi;
    // delete cBackground;
    // delete cCorrected;
    // delete cProjection;  

    fout->Close();
    t.Print();
    return 0;
}