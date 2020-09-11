

#include <iostream>
#include <chrono>

// Root headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"


void GetTrueEnergy(Double_t* trueElectronEnergy, Float_t* Pxntu, Float_t* Pyntu, Float_t* Pzntu)
{

    const Double_t electron_rest_mass = 1.0e-3 * 0.51099895; // GeV
    const Double_t m = electron_rest_mass;
    const Double_t m2 = electron_rest_mass * electron_rest_mass;
    Double_t p0_2 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0]; // GeV
    Double_t electron_energy_0 = std::sqrt(p0_2 + m2) - m; // GeV
    Double_t p1_2 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    Double_t electron_energy_1 = std::sqrt(p1_2 + m2) - m;
    //Double_t electron_energy_0 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0];
    //electron_energy_0 = 1.0e3 * std::sqrt(electron_energy_0);
    //Double_t electron_energy_1 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    //electron_energy_1 = 1.0e3 * std::sqrt(electron_energy_1);
//    if(electron_energy_0 < electron_energy_1) std::swap(electron_energy_0, electron_energy_1);
#if 0
    if(electron_energy_0 < electron_energy_1)
    {
        //std::cout << "swap TRUE" << std::endl;
        ++ swap_count_true;
    }
    else
    {
        //std::cout << "swap FALSE" << std::endl;
        ++ swap_count_false;
    }
#endif

#define DEBUG 0
#if DEBUG
    std::cout << "original files:  trueEnergy(" << 1.0e+03 * electron_energy_0 << ", " << 1.0e+03 * electron_energy_1 << ") [MeV] " << std::endl;
    std::cout << "original files:  recoEnergy(" << 1.0e+03 * Sc[0][8] << ", " << 1.0e+03 * Sc[1][8] << ") [MeV] " << std::endl;
    std::cout << "processed files: recoEnergy(" << electronEnergy[0] << ", " << electronEnergy[1] << ") [MeV] " << std::endl;
    //std::cout << "Electron Energy: " << electronEnergy[0] << ", " << 1.0e3 * electron_energy_0 << std::endl;
    //std::cout << "Electron Energy: " << electronEnergy[1] << ", " << 1.0e3 * electron_energy_1 << std::endl;
    //std::cout << "Ecc: " << Ecc[0] << ", " << Ecc[1] << std::endl;
    //std::cout << "Sc: " << 1.0e+03 * Sc[0][8] << ", " << 1.0e+03 * Sc[1][8] << std::endl;
#endif

    // TODO: can simply set branch addresses to input variables?
    // not sure if multiplication is required
    trueElectronEnergy[0] = 1.0e+03 * electron_energy_0;
    trueElectronEnergy[1] = 1.0e+03 * electron_energy_1;
    
}

void SearchFunction(
    int &return_flag,
    Int_t Run,
    Double_t trueVertexR,
    Double_t trueVertexSector,
    Double_t trueVertexZ,
    Int_t run,
    Float_t Xvntu,
    Float_t Yvntu,
    Float_t Zvntu,
    Long64_t ix,
    Long64_t ix_chain)
{
        
    Double_t rntu = std::sqrt(Xvntu * Xvntu + Yvntu * Yvntu);
    Double_t thetantu = std::atan2(Yvntu, Xvntu);
    Double_t sectorntu = (2.5 / std::atan(1.0)) * thetantu; // / 18.0;

    // match variables
    if(-run == Run)
    {
        // if run matches, check if event is a complete match
        // run will always match

        bool close_match_R = false;
        bool close_match_sector = false;
        bool close_match_Z = false;
        int close_match_count = 0;
        if(std::abs(trueVertexR - rntu) < 5.0e-5)
        //if(trueVertexR == rntu)
        {
            close_match_R = true;
            ++ close_match_count;
        }
        if(std::abs(trueVertexSector - sectorntu) < 1.0e-6)
        //if(trueVertexSector == sectorntu)
        {
            close_match_sector = true;
            ++ close_match_count;
        }
        //if(std::abs(trueVertexZ - Zvntu) < 1.0e-5)
        if(trueVertexZ == Zvntu)
        {
            close_match_Z = true;
            ++ close_match_count;
        }

        // check all
        if(close_match_count == 3)
        {
        /*
            std::cout << std::endl;
            std::cout << "ix_chain=" << ix_chain << " ix=" << ix << " close match for trueVertexR; difference: " << std::abs(trueVertexR - rntu) << std::endl;
            std::cout << "ix_chain=" << ix_chain << " ix=" << ix << " close match for trueVertexSector; difference: " << std::abs(trueVertexSector - sectorntu) << std::endl;
            std::cout << "ix_chain=" << ix_chain << " ix=" << ix << " close match for trueVertexZ; difference: " << std::abs(trueVertexZ - Zvntu) << std::endl;
        */
            //std::cin.get();

            return_flag = 1;
        }
        else
        {
            /*
            if(close_match_count > 0)
            {
                std::cout << std::endl;
            }
            */
            /*
            if(close_match_R == true)
            {
                std::cout << "ix_chain=" << ix_chain << " ix=" << ix << " close match for trueVertexR; difference: " << std::abs(trueVertexR - rntu) << std::endl;
            }
            if(close_match_sector == true)
            {
                std::cout << "ix_chain=" << ix_chain << " ix=" << ix << " close match for trueVertexSector; difference: " << std::abs(trueVertexSector - sectorntu) << std::endl;
            }
            if(close_match_Z == true)
            {
                std::cout << "ix_chain=" << ix_chain << " ix=" << ix << " close match for trueVertexZ; difference: " << std::abs(trueVertexZ - Zvntu) << std::endl;
            }
            */
        }
    }
}

void filestripper()
{

    // QuickSort Test
    //qst();

    // datalist file:
    // /home/blotsd/NEMO3/Nd150_analysis/DataLists/nd150_61_rot_nd150.lst


    const int NUM_NAMES = 20;
    
    // input set A
    // the files which Summer used to create Nd150_2eNg_output.root
    // these are the "small" versions which contain only the minimal
    // information
    std::string names_input_A[NUM_NAMES];
    names_input_A[0] = "foils_nd150_61_rot_nd150_101_01_small.root";
    names_input_A[1] = "foils_nd150_61_rot_nd150_101_02_small.root";
    names_input_A[2] = "foils_nd150_61_rot_nd150_102_01_small.root";
    names_input_A[3] = "foils_nd150_61_rot_nd150_102_02_small.root";
    names_input_A[4] = "foils_nd150_61_rot_nd150_103_01_small.root";
    names_input_A[5] = "foils_nd150_61_rot_nd150_103_02_small.root";
    names_input_A[6] = "foils_nd150_61_rot_nd150_104_01_small.root";
    names_input_A[7] = "foils_nd150_61_rot_nd150_104_02_small.root";
    names_input_A[8] = "foils_nd150_61_rot_nd150_105_01_small.root";
    names_input_A[9] = "foils_nd150_61_rot_nd150_105_02_small.root";
    names_input_A[10] = "foils_nd150_61_rot_nd150_106_01_small.root";
    names_input_A[11] = "foils_nd150_61_rot_nd150_106_02_small.root";
    names_input_A[12] = "foils_nd150_61_rot_nd150_107_01_small.root";
    names_input_A[13] = "foils_nd150_61_rot_nd150_107_02_small.root";
    names_input_A[14] = "foils_nd150_61_rot_nd150_108_01_small.root";
    names_input_A[15] = "foils_nd150_61_rot_nd150_108_02_small.root";
    names_input_A[16] = "foils_nd150_61_rot_nd150_109_01_small.root";
    names_input_A[17] = "foils_nd150_61_rot_nd150_109_02_small.root";
    names_input_A[18] = "foils_nd150_61_rot_nd150_110_01_small.root";
    names_input_A[19] = "foils_nd150_61_rot_nd150_110_02_small.root";

    // input set B
    // Summers Nd150_2eNg_output.root split into individual
    // files
    std::string names_input_B[NUM_NAMES];
    names_input_B[0] = "Nd150_2eNg_output_split_nd150_101_01.root";
    names_input_B[1] = "Nd150_2eNg_output_split_nd150_101_02.root";
    names_input_B[2] = "Nd150_2eNg_output_split_nd150_102_01.root";
    names_input_B[3] = "Nd150_2eNg_output_split_nd150_102_02.root";
    names_input_B[4] = "Nd150_2eNg_output_split_nd150_103_01.root";
    names_input_B[5] = "Nd150_2eNg_output_split_nd150_103_02.root";
    names_input_B[6] = "Nd150_2eNg_output_split_nd150_104_01.root";
    names_input_B[7] = "Nd150_2eNg_output_split_nd150_104_02.root";
    names_input_B[8] = "Nd150_2eNg_output_split_nd150_105_01.root";
    names_input_B[9] = "Nd150_2eNg_output_split_nd150_105_02.root";
    names_input_B[10] = "Nd150_2eNg_output_split_nd150_106_01.root";
    names_input_B[11] = "Nd150_2eNg_output_split_nd150_106_02.root";
    names_input_B[12] = "Nd150_2eNg_output_split_nd150_107_01.root";
    names_input_B[13] = "Nd150_2eNg_output_split_nd150_107_02.root";
    names_input_B[14] = "Nd150_2eNg_output_split_nd150_108_01.root";
    names_input_B[15] = "Nd150_2eNg_output_split_nd150_108_02.root";
    names_input_B[16] = "Nd150_2eNg_output_split_nd150_109_01.root";
    names_input_B[17] = "Nd150_2eNg_output_split_nd150_109_02.root";
    names_input_B[18] = "Nd150_2eNg_output_split_nd150_110_01.root";
    names_input_B[19] = "Nd150_2eNg_output_split_nd150_110_02.root";


    ///////////////////////////////////////////////////////////////////////////
    // variables
    ///////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////
    // input B and output

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

    Double_t trueElectronEnergy[2];


    ///////////////////////////////////////////////////////////////////////////
    // input A

    Int_t Nsc;
    Float_t Sc[2000][12];
    Int_t run = 0;
    Int_t Nvntu = 0;
    Float_t Xvntu = 0.0;
    Float_t Yvntu = 0.0;
    Float_t Zvntu = 0.0;
    Int_t Ntntu = 0;
    Float_t Pxntu[31];
    Float_t Pyntu[31];
    Float_t Pzntu[31];


    ///////////////////////////////////////////////////////////////////////////
    // initialize output file
    ///////////////////////////////////////////////////////////////////////////

    TString foutput_dir = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/output_split";
    TString foutput_fname = "Nd150_2eNg_output_truth_NEWNEW_.root";
    TFile *foutput = new TFile(foutput_dir + foutput_fname, "recreate");
    TDirectory *doutput = foutput->mkdir("Nd150_2eNg");
    foutput->cd("Nd150_2eNg");
    TTree *toutput = new TTree("Nd150_2eNg", "Nd150_2eNg");


    ///////////////////////////////////////////////////////////////////////
    // output tree
    ///////////////////////////////////////////////////////////////////////

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

    toutput->Branch("trueElectronEnergy", &trueElectronEnergy, "trueElectronEnergy[2]/D");



    ///////////////////////////////////////////////////////////////////////////
    // file loop
    ///////////////////////////////////////////////////////////////////////////

    std::cout << "start main loop" << std::endl;
    for(int file_index = 0; file_index < NUM_NAMES; ++ file_index)
    {
        std::cout << "********************************************************************************" << std::endl;
        std::cout << "********************************************************************************" << std::endl;
        std::cout << "file_index=" << file_index << std::endl;

        // time execution of pass through single file
        std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
        
        ///////////////////////////////////////////////////////////////////////
        // initialize files
        ///////////////////////////////////////////////////////////////////////

        TString finput_dir_A = "/mnt/ecb/unix/nemo2/reco/reco_20130206_fix/mc_8.0/rot/";
        TString finput_name_A = names_input_A[file_index];
        TFile *finput_A = new TFile(finput_dir_A + finput_name_A);
        if(!finput_A->IsOpen())
        {
            std::cout << "Error opening file " << finput_name_A << std::endl;
        }
        TTree *tinput_A = (TTree*)finput_A->Get("h10");

        TString finput_dir_B = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/input_split";
        TString finput_name_B = names_input_B[file_index];
        TFile *finput_B = new TFile(finput_dir_B + finput_name_B);
        if(!finput_B->IsOpen())
        {
            std::cout << "Error opening file " << finput_name_B << std::endl;
        }
        TTree *tinput_B = (TTree*)finput_B->Get("Nd150_2eNg/Nd150_2eNg");



        ///////////////////////////////////////////////////////////////////////
        // initialize trees
        ///////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////
        // A tree
        ///////////////////////////////////////////////////////////////////////

        tinput_A->SetBranchAddress("Nsc", &Nsc);
        tinput_A->SetBranchAddress("Sc", Sc);
        tinput_A->SetBranchAddress("run", &run);
        tinput_A->SetBranchAddress("Nvntu", &Nvntu);
        tinput_A->SetBranchAddress("Xvntu", &Xvntu);
        tinput_A->SetBranchAddress("Yvntu", &Yvntu);
        tinput_A->SetBranchAddress("Zvntu", &Zvntu);
        tinput_A->SetBranchAddress("Ntntu", &Ntntu);
        tinput_A->SetBranchAddress("Pxntu", Pxntu);
        tinput_A->SetBranchAddress("Pyntu", Pyntu);
        tinput_A->SetBranchAddress("Pzntu", Pzntu);

        
        ///////////////////////////////////////////////////////////////////////
        // B tree
        ///////////////////////////////////////////////////////////////////////

        tinput_B->SetBranchAddress("Event", &Event);
        tinput_B->SetBranchAddress("Run", &Run);
        tinput_B->SetBranchAddress("runStatus", &runStatus);
        tinput_B->SetBranchAddress("nElectrons", &nElectrons);
        tinput_B->SetBranchAddress("radonWeight", &radonWeight);
        tinput_B->SetBranchAddress("bi210Weight", &bi210Weight);
        tinput_B->SetBranchAddress("foilSide", &foilSide);
        tinput_B->SetBranchAddress("eventTime", &eventTime);
        tinput_B->SetBranchAddress("trueVertexR", &trueVertexR);
        tinput_B->SetBranchAddress("trueVertexZ", &trueVertexZ);
        tinput_B->SetBranchAddress("trueVertexSector", &trueVertexSector);
        tinput_B->SetBranchAddress("trueVertexLayer", &trueVertexLayer);
        tinput_B->SetBranchAddress("electronEnergy", electronEnergy);
        tinput_B->SetBranchAddress("eTrackLength", eTrackLength);
        tinput_B->SetBranchAddress("electronSide", electronSide);
        tinput_B->SetBranchAddress("trackSign", trackSign);
        tinput_B->SetBranchAddress("electronMeasTime", electronMeasTime);
        tinput_B->SetBranchAddress("electronDMeasTime", electronDMeasTime);
        tinput_B->SetBranchAddress("electronBlockType", electronBlockType);
        tinput_B->SetBranchAddress("internalPull", &internalPull);
        tinput_B->SetBranchAddress("internalProb", &internalProb);
        tinput_B->SetBranchAddress("externalPull", &externalPull);
        tinput_B->SetBranchAddress("externalProb", &externalProb);
        tinput_B->SetBranchAddress("cosee", &cosee);
        tinput_B->SetBranchAddress("cosee_weight", &cosee_weight);
        tinput_B->SetBranchAddress("electronPMT", electronPMT);
        tinput_B->SetBranchAddress("electronLDFlag", electronLDFlag);
        tinput_B->SetBranchAddress("electronLDCorr", electronLDCorr);
        tinput_B->SetBranchAddress("electronLDCorrErr", electronLDCorrErr);
        tinput_B->SetBranchAddress("vertexZ", vertexZ);
        tinput_B->SetBranchAddress("vertexSec", vertexSec);
        tinput_B->SetBranchAddress("vertexR", vertexR);
        tinput_B->SetBranchAddress("vertexInHotSpot", vertexInHotSpot);
        tinput_B->SetBranchAddress("firstGgHitLayer", firstGgHitLayer);
        tinput_B->SetBranchAddress("lastGgHitLayer", lastGgHitLayer);
        tinput_B->SetBranchAddress("NAPromptGgHits", &NAPromptGgHits);
        tinput_B->SetBranchAddress("NAPromptGgHitsSide", &NAPromptGgHitsSide);
        tinput_B->SetBranchAddress("NAPromptGgHitsDist2Vertex", NAPromptGgHitsDist2Vertex);
        tinput_B->SetBranchAddress("NAPromptGgHitsDist2Calo", NAPromptGgHitsDist2Calo);
        tinput_B->SetBranchAddress("nGammaClusters", &nGammaClusters);
        tinput_B->SetBranchAddress("nInCluster", nInCluster);
        tinput_B->SetBranchAddress("clusterEnergy", clusterEnergy);
        tinput_B->SetBranchAddress("clusterTimeSpan", clusterTimeSpan);
        tinput_B->SetBranchAddress("nTotalClusterHits", &nTotalClusterHits);
        tinput_B->SetBranchAddress("clusterHitEnergy", clusterHitEnergy);
        tinput_B->SetBranchAddress("clusterHitPMT", clusterHitPMT);
        tinput_B->SetBranchAddress("clusterHitLDFlag", clusterHitLDFlag);
        tinput_B->SetBranchAddress("clusterHitLDCorr", clusterHitLDCorr);
        tinput_B->SetBranchAddress("clusterHitLDCorrErr", clusterHitLDCorrErr);
        tinput_B->SetBranchAddress("clusterHitSec", clusterHitSec);
        tinput_B->SetBranchAddress("clusterHitZ", clusterHitZ);


        // count number of events saved to output file
        Long64_t count{0};
        Long64_t max_B{tinput_B->GetEntries()};

        ///////////////////////////////////////////////////////////////////////////
        // start data processing loop
        ///////////////////////////////////////////////////////////////////////////

        //Long64_t ix_faster = 0;
        for(Long64_t ix{0}; ix < max_B; ++ ix)
        {
            tinput_B->GetEntry(ix);

            //std::cout << "ix=" << ix << " / " << max << std::endl;
            //std::cout << "looking for Run=" << Run << std::endl;

            //std::cout << "Searching, ix=" << ix << std::endl;
        
            // new 2020-01-14
            //Run_current = Run;
            Long64_t max_A{tinput_A->GetEntries()};
            Long64_t ix_A_start, ix_A_end, ix_A_match;
            //std::cout << "ix_faster=" << ix_faster << std::endl;
            //for(Long64_t ix_A{ix_faster}; ix_A < max_A; ++ ix_A)
            for(Long64_t ix_A{0}; ix_A < max_A; ++ ix_A)
            {
                tinput_A->GetEntry(ix_A);
                //std::cout << "ix_A=" << ix_A << " run=" << run << std::endl;
                //std::cin.get();
            
                // find index of first matching run number
                if(-run == Run)
                {
                    ix_A_start = ix_A;
                    //std::cout << "ix_A_start=" << ix_A_start << std::endl;
                    break;
                }
            }
            for(Long64_t ix_A{ix_A_start + 1}; ix_A < max_A; ++ ix_A)
            {
                tinput_A->GetEntry(ix_A);
                //std::cout << "ix_A=" << ix_A << " run=" << run << std::endl;
                //std::cin.get();

                // find index of last matching run number, + 1
                // (this bound is not included in the search)
                if(-run == Run)
                {
                    //ix_A_end = ix_A;
                }
                else
                {
                    ix_A_end = ix_A;
                    //ix_A_faster = ix_A;
                    //std::cout << "ix_A_end=" << ix_A_end << " ix_A_faster=" << ix_A_faster << std::endl;
                    break;
                }
            }

            //std::cout << "ix=" << ix << ", ix_A_start=" << ix_A_start << " ix_A_end=" << ix_A_end << std::endl;
            //std::cin.get();

            // ix_A_start and ix_A_end are now set
            // count the number of matches
            Long64_t multi_match_count = 0;
            for(Long64_t ix_A{ix_A_start}; ix_A < ix_A_end; ++ ix_A)
            {
                tinput_A->GetEntry(ix_A); 
        
                int return_flag = 0;

                /*            
                void SearchFunction(
                    int &return_flag,
                    Int_t Run,
                    Double_t trueVertexR,
                    Double_t trueVertexSector,
                    Double_t trueVertexZ,
                    Int_t run,
                    Float_t Xvntu,
                    Float_t Yvntu,
                    Float_t Zvntu,
                    Long64_t ix,
                    Long64_t ix_A)
                */
                SearchFunction(return_flag, Run, trueVertexR, trueVertexSector, trueVertexZ, run, Xvntu, Yvntu, Zvntu, ix, ix_A);
                if(return_flag == 1)
                {
                    if(multi_match_count == 0)
                    {
                        ix_A_match = ix_A;
                    }
                    multi_match_count ++;
                }
            }

            // check the number of matches
            if(multi_match_count == 1)
            {
                //std::cout << "match for ix=" << ix << " found: ix_A_match=" << ix_A_match << std::endl;
                tinput_A->GetEntry(ix_A_match);

                /*
                Double_t rntu = std::sqrt(Xvntu * Xvntu + Yvntu * Yvntu);
                Double_t thetantu = std::atan2(Yvntu, Xvntu);
                Double_t sectorntu = (2.5 / std::atan(1.0)) * thetantu; // / 18.0;
                std::cout << "ix_A_match=" << ix_A_match << " ix=" << ix << " close match for trueVertexR; difference: " << std::abs(trueVertexR - rntu) << std::endl;
                std::cout << "ix_A_match=" << ix_A_match << " ix=" << ix << " close match for trueVertexSector; difference: " << std::abs(trueVertexSector - sectorntu) << std::endl;
                std::cout << "ix_A_match=" << ix_A_match << " ix=" << ix << " close match for trueVertexZ; difference: " << std::abs(trueVertexZ - Zvntu) << std::endl;
                */

                // write to file
                GetTrueEnergy(trueElectronEnergy, Pxntu, Pyntu, Pzntu);
                ++ count;
                toutput->Fill();
            }
            else if(multi_match_count == 0)
            {
                std::cout << "Error: multi_match_count=0 ! : xi=" << ix << " Run=" << Run <<  std::endl;
                std::cin.get();
            }
            else
            {
                std::cout << "Error: multi_match_count=" << multi_match_count << " ! : ix=" << ix << " Run=" << Run << std::endl;
                std::cin.get();
            }
        }

        std::cout << "found " << count << " valid matches out of " << max_B << " events" << std::endl;


        // close input files
        finput_A->Close();
        finput_B->Close();

        // log time taken per file
        std::chrono::system_clock::time_point current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> runtime_seconds = current_time - start_time;
        std::cout << "Runtime: " << runtime_seconds.count() << " s, " << runtime_seconds.count() / 3600.0 << " h" << std::endl;

        std::cout << std::endl;
    }
    
    foutput->cd();
    toutput->Write();
    foutput->Close();


    std::cout << "done" << std::endl;






    



    // check tchain run numbers in order
    // 2020-09-07: check passed on foils_nd150_61_rot_nd150_1xx_xx.root
    #if 0
    {
        Long64_t max_chain{tchain->GetEntries()};
        tchain->GetEntry(0);
        Long64_t run_last = run;
        for(Long64_t ix_chain{1}; ix_chain < max_chain; ++ ix_chain)
        {
            tchain->GetEntry(ix_chain);

            if(run == run_last)
            {
                // do nothing, ok
            }
            else if(run == run_last - 1)
            {
                // do nothing, ok
            }
            else if(run < run_last - 1)
            {
                std::cout << "run: " << run_last << " -> " << run << std::endl;
            }
            else
            {
                std::cout << "ERROR: run=" << run << " run_last=" << run_last << std::endl;
                std::cin.get();
            }

            run_last = run;
        }
    }
    #endif

    // check tinput run numbers in order TODO
    #if 0
    {
        Long64_t max{tinput->GetEntries()};
        tinput->GetEntry(0);
        Long64_t Run_last = Run;
        for(Long64_t ix{1}; ix < max; ++ ix)
        {
            tinput->GetEntry(ix);

            if(Run == Run_last)
            {
                // do nothing, ok
            }
            else if(Run == Run_last + 1)
            {
                // do nothing, ok
                //std::cout << "Run=" << Run << std::endl;
            }
            else if(Run > Run_last + 1)
            {
                //std::cout << "Run: " << Run_last << " -> " << Run << std::endl;
            }
            else
            {
                std::cout << "ERROR: ix=" << ix << " Run=" << Run << " Run_last=" << Run_last << std::endl;
                //std::cin.get();
            }

            Run_last = Run;
        }
    }

    std::cout << "checks done" << std::endl;
    #endif


}
