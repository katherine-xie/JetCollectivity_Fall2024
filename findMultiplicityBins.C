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

Int_t numMultBins = 4; 

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

void findMultiplicityBins() {
    initializeChain();

    reader.Restart(); // Ensuring event loop starts from beginning
    std::vector<Int_t> multVec;

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

    // Find min and max
    auto min = std::min_element(multVec.begin(), multVec.end());
    auto max = std::max_element(multVec.begin(), multVec.end());
    Int_t minMultiplicity = *min;
    Int_t maxMultiplicity = *max;


    // Calculate the bin width
    Int_t range = maxMultiplicity - minMultiplicity;
    Int_t binWidth = (range + 3) / numMultBins; // +3 for rounding up in case of integer division (to accomodate all multiplicities)

    // Define bin boundaries
    std::vector<std::pair<int, int>> bins(4);
    for (Int_t i = 0; i < numMultBins; i++) {
        bins[i].first = minMultiplicity + (i * binWidth);
        bins[i].second = (i == 3) ? maxMultiplicity : (minMultiplicity + (i + 1) * binWidth - 1);
    }

    // Output bin boundaries
    for (int i = 0; i < bins.size(); i++) {
        std::cout << "Bin " << i + 1 << ": [" << bins[i].first << ", " << bins[i].second << "]\n";
    }

    t.Print();

}
