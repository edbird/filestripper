

#include <iostream>
#include <chrono>

// Root headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"



void filesplitter()
{

    // QuickSort Test
    //qst();

    // datalist file:
    // /home/blotsd/NEMO3/Nd150_analysis/DataLists/nd150_61_rot_nd150.lst


    const int NUM_NAMES = 20;
    std::string names[NUM_NAMES];
    names[0] = "Nd150_2eNg_output_split_nd150_101_01.root";
    names[1] = "Nd150_2eNg_output_split_nd150_101_02.root";
    names[2] = "Nd150_2eNg_output_split_nd150_102_01.root";
    names[3] = "Nd150_2eNg_output_split_nd150_102_02.root";
    names[4] = "Nd150_2eNg_output_split_nd150_103_01.root";
    names[5] = "Nd150_2eNg_output_split_nd150_103_02.root";
    names[6] = "Nd150_2eNg_output_split_nd150_104_01.root";
    names[7] = "Nd150_2eNg_output_split_nd150_104_02.root";
    names[8] = "Nd150_2eNg_output_split_nd150_105_01.root";
    names[9] = "Nd150_2eNg_output_split_nd150_105_02.root";
    names[10] = "Nd150_2eNg_output_split_nd150_106_01.root";
    names[11] = "Nd150_2eNg_output_split_nd150_106_02.root";
    names[12] = "Nd150_2eNg_output_split_nd150_107_01.root";
    names[13] = "Nd150_2eNg_output_split_nd150_107_02.root";
    names[14] = "Nd150_2eNg_output_split_nd150_108_01.root";
    names[15] = "Nd150_2eNg_output_split_nd150_108_02.root";
    names[16] = "Nd150_2eNg_output_split_nd150_109_01.root";
    names[17] = "Nd150_2eNg_output_split_nd150_109_02.root";
    names[18] = "Nd150_2eNg_output_split_nd150_110_01.root";
    names[19] = "Nd150_2eNg_output_split_nd150_110_02.root";


    // working at home
    TString finput_dir = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/";
    TFile *finput = new TFile(finput_dir + "Nd150_2eNg_output.root");
    TTree *tinput = (TTree*)finput->Get("Nd150_2eNg/Nd150_2eNg");


    Int_t Event = 0; 
    Int_t Run = 0; 
    Int_t runStatus = 0; 
    Int_t nElectrons = 0; 
    Double_t radonWeight = 0.0; 
    Double_t bi210Weight = 0.0;
    Int_t foilSide = 0;
    Double_t eventTime = 0.0;
    Double_t trueVertexR = 0.0;
    Double_t trueVertexZ = 0.0;
    Double_t trueVertexSector = 0.0;
    Int_t trueVertexLayer = 0;
    Double_t electronEnergy[2];
    Double_t eTrackLength[2];
    Int_t electronSide[2];
    Double_t trackSign[2];
    Double_t electronMeasTime[2];
    Double_t electronDMeasTime[2];
    Int_t electronBlockType[2];
    Double_t internalPull = 0.0;
    Double_t internalProb = 0.0;
    Double_t externalPull = 0.0;
    Double_t externalProb = 0.0;
    Double_t cosee = 0.0;
    Double_t cosee_weight = 0.0;
    Int_t electronPMT[2];
    Int_t electronLDFlag[2];
    Double_t electronLDCorr[2];
    Double_t electronLDCorrErr[2];
    Double_t vertexZ[2];
    Double_t vertexSec[2];
    Double_t vertexR[2];
    Bool_t vertexInHotSpot[2];
    Int_t firstGgHitLayer[2];
    Int_t lastGgHitLayer[2];
    Int_t NAPromptGgHits = 0;
    Int_t NAPromptGgHitsSide[5];
    Double_t NAPromptGgHitsDist2Vertex[5]; //[NAPromptGgHits];
    Double_t NAPromptGgHitsDist2Calo[5]; //[NAPromptGgHits];
    Int_t nGammaClusters = 0;
    Int_t nInCluster[10]; //[nGammaClusters];
    Double_t clusterEnergy[10]; //[nGammaClusters];
    Double_t clusterTimeSpan[10]; //[nGammaClusters];
    Int_t nTotalClusterHits = 0;
    Double_t clusterHitEnergy[20]; //[nTotalClusterHits];
    Int_t clusterHitPMT[20]; //[nTotalClusterHits];
    Int_t clusterHitLDFlag[20]; //[nTotalClusterHits];
    Double_t clusterHitLDCorr[20]; //[nTotalClusterHits];
    Double_t clusterHitLDCorrErr[20]; //[nTotalClusterHits];
    Double_t clusterHitSec[20]; //[nTotalClusterHits];
    Double_t clusterHitZ[20]; //[nTotalClusterHits];

    tinput->SetBranchAddress("Event", &Event);
    tinput->SetBranchAddress("Run", &Run);
    tinput->SetBranchAddress("runStatus", &runStatus);
    tinput->SetBranchAddress("nElectrons", &nElectrons);
    tinput->SetBranchAddress("radonWeight", &radonWeight);
    tinput->SetBranchAddress("bi210Weight", &bi210Weight);
    tinput->SetBranchAddress("foilSide", &foilSide);
    tinput->SetBranchAddress("eventTime", &eventTime);
    tinput->SetBranchAddress("trueVertexR", &trueVertexR);
    tinput->SetBranchAddress("trueVertexZ", &trueVertexZ);
    tinput->SetBranchAddress("trueVertexSector", &trueVertexSector);
    tinput->SetBranchAddress("trueVertexLayer", &trueVertexLayer);
    tinput->SetBranchAddress("electronEnergy", electronEnergy);
    tinput->SetBranchAddress("eTrackLength", eTrackLength);
    tinput->SetBranchAddress("electronSide", electronSide);
    tinput->SetBranchAddress("trackSign", trackSign);
    tinput->SetBranchAddress("electronMeasTime", electronMeasTime);
    tinput->SetBranchAddress("electronDMeasTime", electronDMeasTime);
    tinput->SetBranchAddress("electronBlockType", electronBlockType);
    tinput->SetBranchAddress("internalPull", &internalPull);
    tinput->SetBranchAddress("internalProb", &internalProb);
    tinput->SetBranchAddress("externalPull", &externalPull);
    tinput->SetBranchAddress("externalProb", &externalProb);
    tinput->SetBranchAddress("cosee", &cosee);
    tinput->SetBranchAddress("cosee_weight", &cosee_weight);
    tinput->SetBranchAddress("electronPMT", electronPMT);
    tinput->SetBranchAddress("electronLDFlag", electronLDFlag);
    tinput->SetBranchAddress("electronLDCorr", electronLDCorr);
    tinput->SetBranchAddress("electronLDCorrErr", electronLDCorrErr);
    tinput->SetBranchAddress("vertexZ", vertexZ);
    tinput->SetBranchAddress("vertexSec", vertexSec);
    tinput->SetBranchAddress("vertexR", vertexR);
    tinput->SetBranchAddress("vertexInHotSpot", vertexInHotSpot);
    tinput->SetBranchAddress("firstGgHitLayer", firstGgHitLayer);
    tinput->SetBranchAddress("lastGgHitLayer", lastGgHitLayer);
    tinput->SetBranchAddress("NAPromptGgHits", &NAPromptGgHits);
    tinput->SetBranchAddress("NAPromptGgHitsSide", &NAPromptGgHitsSide);
    tinput->SetBranchAddress("NAPromptGgHitsDist2Vertex", NAPromptGgHitsDist2Vertex);
    tinput->SetBranchAddress("NAPromptGgHitsDist2Calo", NAPromptGgHitsDist2Calo);
    tinput->SetBranchAddress("nGammaClusters", &nGammaClusters);
    tinput->SetBranchAddress("nInCluster", nInCluster);
    tinput->SetBranchAddress("clusterEnergy", clusterEnergy);
    tinput->SetBranchAddress("clusterTimeSpan", clusterTimeSpan);
    tinput->SetBranchAddress("nTotalClusterHits", &nTotalClusterHits);
    tinput->SetBranchAddress("clusterHitEnergy", clusterHitEnergy);
    tinput->SetBranchAddress("clusterHitPMT", clusterHitPMT);
    tinput->SetBranchAddress("clusterHitLDFlag", clusterHitLDFlag);
    tinput->SetBranchAddress("clusterHitLDCorr", clusterHitLDCorr);
    tinput->SetBranchAddress("clusterHitLDCorrErr", clusterHitLDCorrErr);
    tinput->SetBranchAddress("clusterHitSec", clusterHitSec);
    tinput->SetBranchAddress("clusterHitZ", clusterHitZ);

    Long64_t max{tinput->GetEntries()};
    tinput->GetEntry(0);
    Int_t Run_last = Run;

    TFile *foutput = nullptr;
    TTree *toutput = nullptr;

    int output_file_index = 0;

    // create new file

    TString foutput_dir = "/mnt/ramdisk/";
    TString foutput_fname = names[output_file_index];

    ///////////////////////////////////////////////////////////////////////////
    // initialize output file
    ///////////////////////////////////////////////////////////////////////////

    foutput = new TFile(foutput_dir + foutput_fname, "recreate");
    TDirectory *doutput = foutput->mkdir("Nd150_2eNg");
    foutput->cd("Nd150_2eNg");

    toutput = new TTree("Nd150_2eNg", "Nd150_2eNg");

    toutput->Branch("Event", &Event, "Event/I");
    toutput->Branch("Run", &Run, "Run/I");
    toutput->Branch("runStatus", &runStatus, "runStatus/I");
    toutput->Branch("nElectrons", &nElectrons, "nElectrons/I");
    toutput->Branch("radonWeight", &radonWeight, "radonWeight/D");
    toutput->Branch("bi210Weight", &bi210Weight, "bi210Weight/D");
    toutput->Branch("foilSide", &foilSide, "foilSide/I");
    toutput->Branch("eventTime", &eventTime, "eventTime/D");
    toutput->Branch("trueVertexR", &trueVertexR, "trueVertexR/D");
    toutput->Branch("trueVertexZ", &trueVertexZ, "trueVertexZ/D");
    toutput->Branch("trueVertexSector", &trueVertexSector, "trueVertexSector/D");
    toutput->Branch("trueVertexLayer", &trueVertexLayer, "trueVertexLayer/I");
    toutput->Branch("electronEnergy", &electronEnergy, "electronEnergy[2]/D");
    toutput->Branch("eTrackLength", &eTrackLength, "eTrackLength[2]/D");
    toutput->Branch("electronSide", &electronSide, "electronSide[2]/I");
    toutput->Branch("trackSign", &trackSign, "trackSign[2]/D");
    toutput->Branch("electronMeasTime", &electronMeasTime, "electronMeasTime[2]/D");
    toutput->Branch("electronDMeasTime", &electronDMeasTime, "electronDMeasTime[2]/D");
    toutput->Branch("electronBlockType", &electronBlockType, "electronBlockType[2]/I");
    toutput->Branch("internalPull", &internalPull, "internalPull/D");
    toutput->Branch("internalProb", &internalProb, "internalProb/D");
    toutput->Branch("externalPull", &externalPull, "externalPull/D");
    toutput->Branch("externalProb", &externalProb, "externalProb/D");
    toutput->Branch("cosee", &cosee, "cosee/D");
    toutput->Branch("cosee_weight", &cosee_weight, "cosee_weight/D");
    toutput->Branch("electronPMT", &electronPMT, "electronPMT[2]/I");
    toutput->Branch("electronLDFlag", &electronLDFlag, "electronLDFlag[2]/I");
    toutput->Branch("electronLDCorr", &electronLDCorr, "electronLDCorr[2]/D");
    toutput->Branch("electronLDCorrErr", &electronLDCorrErr, "electronLDCorrErr[2]/D");
    toutput->Branch("vertexZ", &vertexZ, "vertexZ[2]/D");
    toutput->Branch("vertexSec", &vertexSec, "vertexSec[2]/D");
    toutput->Branch("vertexR", &vertexR, "vertexR[2]/D");
    toutput->Branch("vertexInHotSpot", &vertexInHotSpot, "vertexInHotSpot[2]/O");
    toutput->Branch("firstGgHitLayer", &firstGgHitLayer, "firstGgHitLayer[2]/I");
    toutput->Branch("lastGgHitLayer", &lastGgHitLayer, "lastGgHitLayer[2]/I");
    toutput->Branch("NAPromptGgHits", &NAPromptGgHits, "NAPromptGgHits/I");
    toutput->Branch("NAPromptGgHitsSide", &NAPromptGgHitsSide, "NAPromptGgHitsSide/I");
    toutput->Branch("NAPromptGgHitsDist2Vertex", &NAPromptGgHitsDist2Vertex, "NAPromptGgHitsDist2Vertex[NAPromptGgHits]/D");
    toutput->Branch("NAPromptGgHitsDist2Calo", &NAPromptGgHitsDist2Calo, "NAPromptGgHitsDist2Calo[NAPromptGgHits]/D");
    toutput->Branch("nGammaClusters", &nGammaClusters, "nGammaClusters/I");
    toutput->Branch("nInCluster", &nInCluster, "nInCluster[nGammaClusters]/I");
    toutput->Branch("clusterEnergy", &clusterEnergy, "clusterEnergy[nGammaClusters]/D");
    toutput->Branch("clusterTimeSpan", &clusterTimeSpan, "clusterTimeSpan[nGammaClusters]/D");
    toutput->Branch("nTotalClusterHits", &nTotalClusterHits, "nTotalClusterHits/I");
    toutput->Branch("clusterHitEnergy", &clusterHitEnergy, "clusterHitEnergy[nTotalClusterHits]/D");
    toutput->Branch("clusterHitPMT", &clusterHitPMT, "clusterHitPMT[nTotalClusterHits]/I");
    toutput->Branch("clusterHitLDFlag", &clusterHitLDFlag, "clusterHitLDFlags[nTotalClusterHits]/I");
    toutput->Branch("clusterHitLDCorr", &clusterHitLDCorr, "clusterHitLDCorr[nTotalClusterHits]/D");
    toutput->Branch("clusterHitLDCorrErr", &clusterHitLDCorrErr, "clusterHitLDCorrErr[nTotalClusterHits]/D");
    toutput->Branch("clusterHitSec", &clusterHitSec, "clusterHitSec[nTotalClusterHits]/D");
    toutput->Branch("clusterHitZ", &clusterHitZ, "clusterHitZ[nTotalClusterHits]/D");

    ++ output_file_index;

    // loop over input
    Long64_t count = 0;
    std::cout << "Processing " << max << " entries" << std::endl;
    for(Long64_t ix_input{0}; ix_input < max; ++ ix_input)
    {
        tinput->GetEntry(ix_input);

        if(Run < Run_last)
        {
            std::cout << "ix_input=" << ix_input << " current file contains " << count << " events" << std::endl;
            count = 0;
            std::cout << "Next file: " << names[output_file_index] << std::endl;
            std::cout << "output_file_index=" << output_file_index << std::endl;

            // close existing file
            if(foutput != nullptr)
            {
                toutput->Write();
                foutput->Close();

                toutput = nullptr;
                foutput = nullptr;
            }

            // create new file

            TString foutput_dir = "/mnt/ramdisk/";
            TString foutput_fname = names[output_file_index];
            ++ output_file_index;

            ///////////////////////////////////////////////////////////////////////////
            // initialize output file
            ///////////////////////////////////////////////////////////////////////////

            foutput = new TFile(foutput_dir + foutput_fname, "recreate");
            TDirectory *doutput = foutput->mkdir("Nd150_2eNg");
            foutput->cd("Nd150_2eNg");

            toutput = new TTree("Nd150_2eNg", "Nd150_2eNg");

            toutput->Branch("Event", &Event, "Event/I");
            toutput->Branch("Run", &Run, "Run/I");
            toutput->Branch("runStatus", &runStatus, "runStatus/I");
            toutput->Branch("nElectrons", &nElectrons, "nElectrons/I");
            toutput->Branch("radonWeight", &radonWeight, "radonWeight/D");
            toutput->Branch("bi210Weight", &bi210Weight, "bi210Weight/D");
            toutput->Branch("foilSide", &foilSide, "foilSide/I");
            toutput->Branch("eventTime", &eventTime, "eventTime/D");
            toutput->Branch("trueVertexR", &trueVertexR, "trueVertexR/D");
            toutput->Branch("trueVertexZ", &trueVertexZ, "trueVertexZ/D");
            toutput->Branch("trueVertexSector", &trueVertexSector, "trueVertexSector/D");
            toutput->Branch("trueVertexLayer", &trueVertexLayer, "trueVertexLayer/I");
            toutput->Branch("electronEnergy", &electronEnergy, "electronEnergy[2]/D");
            toutput->Branch("eTrackLength", &eTrackLength, "eTrackLength[2]/D");
            toutput->Branch("electronSide", &electronSide, "electronSide[2]/I");
            toutput->Branch("trackSign", &trackSign, "trackSign[2]/D");
            toutput->Branch("electronMeasTime", &electronMeasTime, "electronMeasTime[2]/D");
            toutput->Branch("electronDMeasTime", &electronDMeasTime, "electronDMeasTime[2]/D");
            toutput->Branch("electronBlockType", &electronBlockType, "electronBlockType[2]/I");
            toutput->Branch("internalPull", &internalPull, "internalPull/D");
            toutput->Branch("internalProb", &internalProb, "internalProb/D");
            toutput->Branch("externalPull", &externalPull, "externalPull/D");
            toutput->Branch("externalProb", &externalProb, "externalProb/D");
            toutput->Branch("cosee", &cosee, "cosee/D");
            toutput->Branch("cosee_weight", &cosee_weight, "cosee_weight/D");
            toutput->Branch("electronPMT", &electronPMT, "electronPMT[2]/I");
            toutput->Branch("electronLDFlag", &electronLDFlag, "electronLDFlag[2]/I");
            toutput->Branch("electronLDCorr", &electronLDCorr, "electronLDCorr[2]/D");
            toutput->Branch("electronLDCorrErr", &electronLDCorrErr, "electronLDCorrErr[2]/D");
            toutput->Branch("vertexZ", &vertexZ, "vertexZ[2]/D");
            toutput->Branch("vertexSec", &vertexSec, "vertexSec[2]/D");
            toutput->Branch("vertexR", &vertexR, "vertexR[2]/D");
            toutput->Branch("vertexInHotSpot", &vertexInHotSpot, "vertexInHotSpot[2]/O");
            toutput->Branch("firstGgHitLayer", &firstGgHitLayer, "firstGgHitLayer[2]/I");
            toutput->Branch("lastGgHitLayer", &lastGgHitLayer, "lastGgHitLayer[2]/I");
            toutput->Branch("NAPromptGgHits", &NAPromptGgHits, "NAPromptGgHits/I");
            toutput->Branch("NAPromptGgHitsSide", &NAPromptGgHitsSide, "NAPromptGgHitsSide/I");
            toutput->Branch("NAPromptGgHitsDist2Vertex", &NAPromptGgHitsDist2Vertex, "NAPromptGgHitsDist2Vertex[NAPromptGgHits]/D");
            toutput->Branch("NAPromptGgHitsDist2Calo", &NAPromptGgHitsDist2Calo, "NAPromptGgHitsDist2Calo[NAPromptGgHits]/D");
            toutput->Branch("nGammaClusters", &nGammaClusters, "nGammaClusters/I");
            toutput->Branch("nInCluster", &nInCluster, "nInCluster[nGammaClusters]/I");
            toutput->Branch("clusterEnergy", &clusterEnergy, "clusterEnergy[nGammaClusters]/D");
            toutput->Branch("clusterTimeSpan", &clusterTimeSpan, "clusterTimeSpan[nGammaClusters]/D");
            toutput->Branch("nTotalClusterHits", &nTotalClusterHits, "nTotalClusterHits/I");
            toutput->Branch("clusterHitEnergy", &clusterHitEnergy, "clusterHitEnergy[nTotalClusterHits]/D");
            toutput->Branch("clusterHitPMT", &clusterHitPMT, "clusterHitPMT[nTotalClusterHits]/I");
            toutput->Branch("clusterHitLDFlag", &clusterHitLDFlag, "clusterHitLDFlags[nTotalClusterHits]/I");
            toutput->Branch("clusterHitLDCorr", &clusterHitLDCorr, "clusterHitLDCorr[nTotalClusterHits]/D");
            toutput->Branch("clusterHitLDCorrErr", &clusterHitLDCorrErr, "clusterHitLDCorrErr[nTotalClusterHits]/D");
            toutput->Branch("clusterHitSec", &clusterHitSec, "clusterHitSec[nTotalClusterHits]/D");
            toutput->Branch("clusterHitZ", &clusterHitZ, "clusterHitZ[nTotalClusterHits]/D");

        }

        toutput->Fill();
        Run_last = Run;
        ++ count;

    }

    // close existing file
    std::cout << "current file contains " << count << " events" << std::endl;
    count = 0;
    std::cout << "output_file_index=" << output_file_index << std::endl;

    if(foutput != nullptr)
    {
        toutput->Write();
        foutput->Close();

        toutput = nullptr;
        foutput = nullptr;
    }
            
    std::cout << "done" << std::endl;

    finput->Close();

}
