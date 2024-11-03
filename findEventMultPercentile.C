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

// Global variables
std::string title = "pp (N_ch >= 60, 13 TeV)";
TChain chain("trackTree");
TTreeReader reader(&chain);
Float_t percentile = 0.99;

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

void findEventMultPercentile() {

    initializeChain();

    // Finding the multiplicity of particles in each event
    std::vector<Int_t> multVec;

    reader.Restart(); // Ensururing that event loop starts at beginning

    // EVENT MULTIPLICITY
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
    multVec.erase(std::remove(multVec.begin(), multVec.end(), 0), multVec.end());
    sort(multVec.begin(), multVec.end()); // Sorting vector

    for (Int_t i = 0; i < multVec.size(); i++) {
        std::cout << i << ".) " << multVec[i] << std::endl;
    }

    // Finding percentile index (I am rounding up)
    Int_t index = std::ceil(percentile * multVec.size()) - 1; // The - 1 is to adjust for a 0-based indexing system in C++
    
    // Finding number of entries greater or equal to the percentile
    Int_t numEntries = 0;
 
    for (Int_t i = 0; i < multVec.size(); i++) {
        if (multVec[i] < multVec[index]) {continue;}
        numEntries++;
    }
    
    std::cout << "There are a total of " << multVec.size() << " events." << std::endl;
    std::cout << "The 99th percentile for the event N_ch is " << multVec[index] << std::endl;
    std::cout << "There are " << numEntries << " entries greater than the 99th percentile" << std::endl;
}