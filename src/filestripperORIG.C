

#include <iostream>
#include <chrono>

// Root headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"


std::size_t swap_count_true = 0;
std::size_t swap_count_false = 0;


void GetTrueEnergy(
    Double_t* trueElectronEnergy,
    Float_t* Pxntu, Float_t* Pyntu, Float_t* Pzntu,
    Double_t* electronEnergy,
    Float_t Sc_0_1, Float_t Sc_1_1,
    Float_t Sc_0_2, Float_t Sc_1_2,
    Float_t Sc_0_3, Float_t Sc_1_3,
    Float_t Sc_0_4, Float_t Sc_1_4,
    Float_t Sc_0_8, Float_t Sc_1_8,
    Int_t *electronPMT)
{

    const Double_t electron_rest_mass = 1.0e-3 * 0.51099895; // GeV
    const Double_t m = electron_rest_mass;
    const Double_t m2 = electron_rest_mass * electron_rest_mass;
    Double_t p0_2 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0]; // GeV
    Double_t trueElectronEnergy_0 = std::sqrt(p0_2 + m2) - m; // GeV
    Double_t p1_2 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    Double_t trueElectronEnergy_1 = std::sqrt(p1_2 + m2) - m;
    //Double_t trueElectronEnergy_0 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0];
    //trueElectronEnergy_0 = 1.0e3 * std::sqrt(trueElectronEnergy_0);
    //Double_t trueElectronEnergy_1 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    //trueElectronEnergy_1 = 1.0e3 * std::sqrt(trueElectronEnergy_1);
//    if(trueElectronEnergy_0 < trueElectronEnergy_1) std::swap(trueElectronEnergy_0, trueElectronEnergy_1);

    if(trueElectronEnergy_0 < trueElectronEnergy_1)
    {
        //std::cout << "swap TRUE" << std::endl;
        ++ swap_count_true;
    }
    else
    {
        //std::cout << "swap FALSE" << std::endl;
        ++ swap_count_false;
    }


#if 0
    // constants
    const int nOrientation[] = {0, 34, 73, 85};
    const int nBlock[] = {17, 13, 3, 3};

    // input variables to construct my block number
    int sectorNumber_0 = (int)Sc_0_1;
    int iobtFlag_0 = (int)Sc_0_2;
    int fcll_0 = (int)Sc_0_3;

    // check against this input
    int blockNumber_0 = (int)Sc_0_4;
    
    // input variables to construct my block number
    int sectorNumber_1 = (int)Sc_1_1;
    int iobtFlag_1 = (int)Sc_1_2;
    int fcll_1 = (int)Sc_1_3;

    // check against this input
    int blockNumber_1 = (int)Sc_1_4;

    int nOrientation_0 = nOrientation[iobtFlag_0];
    int nBlock_0 = nBlock[iobtFlag_0];
    int column_0 = (iobtFlag_0 == 0) ? fcll_0 / 2 : fcll_0;

    int nOrientation_1 = nOrientation[iobtFlag_1];
    int nBlock_1 = nBlock[iobtFlag_1];
    int column_1 = (iobtFlag_1 == 0) ? fcll_1 / 2 : fcll_1;
    
    Int_t my_Sc_0_4 = sectorNumber_0 * 97 + nOrientation_0 + column_0 * nBlock_0 + blockNumber_0 + 1;
    Int_t my_Sc_1_4 = sectorNumber_1 * 97 + nOrientation_1 + column_1 * nBlock_1 + blockNumber_1 + 1;
    
    // electronPMT[x] contains the PMT block number as read from processed 
    // input files
    //
    // Sc_0_4 / Sc_1_4 contains the PMT block number as read from the tchain
    //
    // my_Sc_0_4 / my_Sc_1_4 contains the PMT block number as reconstructed
    // from the processing module which generates the <input files> from
    // the tchain
    //
    // so if electronPMT[x] does not match my_Sc_x_4 then this suggests the
    // processing module swapped the values related to the electron
    // which would include electronPMT and electronEnergy
    // but we have no clue what it did with the trueEnergy
    // this variable may never have even been constructed
    // regardless, we have no clue what it did with the true momentum
    // Pxntu etc

    std::cout << "Check block number" << std::endl;
    std::cout << "block numbers for chain: " << Sc_0_4 << ", " << Sc_1_4 << std::endl;
    std::cout << "block numbers for input: " << electronPMT[0] << "," << electronPMT[1] << std::endl;
    std::cout << "my block numbers for input: " << my_Sc_0_4 << "," << my_Sc_1_4 << std::endl;

    
    bool swap = false;
    if(Sc_0_4 == Sc_1_4)
    {
        std::cout << "BLOCK NUMBERS ARE THE SAME" << std::endl;
        std::cin.get();
    }
    else if(Sc_0_4 == my_Sc_1_4)
    {
        if(Sc_1_4 == my_Sc_0_4)
        {
            swap = true;
            std::cout << "SWAP DETECTED USING BLOCK NUMBERS" << std::endl;
            std::cin.get();

            std::swap(Sc_0_4, Sc_1_4); // this is irrelevant futher in function
            std::swap(Sc_0_8, Sc_1_8); // this needs to be swapped
            std::swap(trueElectronEnergy_0, trueElectronEnergy_1); // this also needs to be swapped
        }
    }

    
    if((std::abs(1.0e+03 * Sc_0_8 - electronEnergy[0]) > 1.0e-8) &&
       (std::abs(1.0e+03 * Sc_1_8 - electronEnergy[1]) > 1.0e-8) &&
       (std::abs(1.0e+03 * Sc_0_8 - electronEnergy[1]) < 1.0e-8) &&
       (std::abs(1.0e+03 * Sc_1_8 - electronEnergy[0]) < 1.0e-8))
    {
        // the data in the input chain is swapped relative to Summers values
        // swap all data associated with input chain
        // this includes the reconstructed scintillator hit energy
        // as well as the true momentum

        std::swap(Sc_0_8, Sc_1_8); // does not do anything outside this function
        std::swap(trueElectronEnergy_0, trueElectronEnergy_1); // this does do something outside of this function
        std::cout << "suspicious event fixed - this should now never occur" << std::endl;
    }
#endif

#define DEBUG 0
#if DEBUG
    std::cout << "original files:  trueEnergy(" << 1.0e+03 * trueElectronEnergy_0 << ", " << 1.0e+03 * trueElectronEnergy_1 << ") [MeV] " << std::endl;
    std::cout << "original files:  recoEnergy(" << 1.0e+03 * Sc_0_8 << ", " << 1.0e+03 * Sc_1_8 << ") [MeV] " << std::endl;
    std::cout << "processed files: recoEnergy(" << electronEnergy[0] << ", " << electronEnergy[1] << ") [MeV] " << std::endl;
    //std::cout << "Electron Energy: " << electronEnergy[0] << ", " << 1.0e3 * trueElectronEnergy_0 << std::endl;
    //std::cout << "Electron Energy: " << electronEnergy[1] << ", " << 1.0e3 * trueElectronEnergy_1 << std::endl;
    //std::cout << "Ecc: " << Ecc[0] << ", " << Ecc[1] << std::endl;
    //std::cout << "Sc: " << 1.0e+03 * Sc[0][8] << ", " << 1.0e+03 * Sc[1][8] << std::endl;
#endif

#if 0
    if((std::abs(1.0e+03 * Sc_0_8 - electronEnergy[0]) > 1.0e-3) ||
       (std::abs(1.0e+03 * Sc_1_8 - electronEnergy[1]) > 1.0e-3))
    {
        std::cout << "suspicious event: energy 0: " << 1.0e+03 * Sc_0_8 << ", " << electronEnergy[0] << std::endl;
        std::cout << "suspicious event: energy 1: " << 1.0e+03 * Sc_1_8 << ", " << electronEnergy[1] << std::endl;
    }

    if(
        (
            (trueElectronEnergy_0 < trueElectronEnergy_1) &&
            (electronEnergy[1] < electronEnergy[0])
        ) ||
        (
            (trueElectronEnergy_1 < trueElectronEnergy_0) &&
            (electronEnergy[0] < electronEnergy[1])
        )
    )
    {
        std::cout << "swap detected" << std::endl;
    }
#endif

    std::cout << std::endl;
    // TODO: can simply set branch addresses to input variables?
    // not sure if multiplication is required
    trueElectronEnergy[0] = 1.0e+03 * trueElectronEnergy_0;
    trueElectronEnergy[1] = 1.0e+03 * trueElectronEnergy_1;
   
    std::cin.get();
}
// release_1.0.0/NemoObjects/src/N3CalorimeterHit.cpp
// release_1.0.0/NemoMods/src/NemorInputModule.cpp

void GetTrueEnergy(
    Double_t* trueElectronEnergy,
    Float_t* Pxntu, Float_t* Pyntu, Float_t* Pzntu)
{

    const Double_t electron_rest_mass = 1.0e-3 * 0.51099895; // GeV
    const Double_t m = electron_rest_mass;
    const Double_t m2 = electron_rest_mass * electron_rest_mass;
    Double_t p0_2 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0]; // GeV
    Double_t trueElectronEnergy_0 = std::sqrt(p0_2 + m2) - m; // GeV
    Double_t p1_2 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    Double_t trueElectronEnergy_1 = std::sqrt(p1_2 + m2) - m;
    //Double_t trueElectronEnergy_0 = Pxntu[0] * Pxntu[0] + Pyntu[0] * Pyntu[0] + Pzntu[0] * Pzntu[0];
    //trueElectronEnergy_0 = 1.0e3 * std::sqrt(trueElectronEnergy_0);
    //Double_t trueElectronEnergy_1 = Pxntu[1] * Pxntu[1] + Pyntu[1] * Pyntu[1] + Pzntu[1] * Pzntu[1];
    //trueElectronEnergy_1 = 1.0e3 * std::sqrt(trueElectronEnergy_1);
//    if(trueElectronEnergy_0 < trueElectronEnergy_1) std::swap(trueElectronEnergy_0, trueElectronEnergy_1);

    if(trueElectronEnergy_0 < trueElectronEnergy_1)
    {
        //std::cout << "swap TRUE" << std::endl;
        ++ swap_count_true;
    }
    else
    {
        //std::cout << "swap FALSE" << std::endl;
        ++ swap_count_false;
    }


#if 0
    // constants
    const int nOrientation[] = {0, 34, 73, 85};
    const int nBlock[] = {17, 13, 3, 3};

    // input variables to construct my block number
    int sectorNumber_0 = (int)Sc_0_1;
    int iobtFlag_0 = (int)Sc_0_2;
    int fcll_0 = (int)Sc_0_3;

    // check against this input
    int blockNumber_0 = (int)Sc_0_4;
    
    // input variables to construct my block number
    int sectorNumber_1 = (int)Sc_1_1;
    int iobtFlag_1 = (int)Sc_1_2;
    int fcll_1 = (int)Sc_1_3;

    // check against this input
    int blockNumber_1 = (int)Sc_1_4;

    int nOrientation_0 = nOrientation[iobtFlag_0];
    int nBlock_0 = nBlock[iobtFlag_0];
    int column_0 = (iobtFlag_0 == 0) ? fcll_0 / 2 : fcll_0;

    int nOrientation_1 = nOrientation[iobtFlag_1];
    int nBlock_1 = nBlock[iobtFlag_1];
    int column_1 = (iobtFlag_1 == 0) ? fcll_1 / 2 : fcll_1;
    
    Int_t my_Sc_0_4 = sectorNumber_0 * 97 + nOrientation_0 + column_0 * nBlock_0 + blockNumber_0 + 1;
    Int_t my_Sc_1_4 = sectorNumber_1 * 97 + nOrientation_1 + column_1 * nBlock_1 + blockNumber_1 + 1;
    
    // electronPMT[x] contains the PMT block number as read from processed 
    // input files
    //
    // Sc_0_4 / Sc_1_4 contains the PMT block number as read from the tchain
    //
    // my_Sc_0_4 / my_Sc_1_4 contains the PMT block number as reconstructed
    // from the processing module which generates the <input files> from
    // the tchain
    //
    // so if electronPMT[x] does not match my_Sc_x_4 then this suggests the
    // processing module swapped the values related to the electron
    // which would include electronPMT and electronEnergy
    // but we have no clue what it did with the trueEnergy
    // this variable may never have even been constructed
    // regardless, we have no clue what it did with the true momentum
    // Pxntu etc

    std::cout << "Check block number" << std::endl;
    std::cout << "block numbers for chain: " << Sc_0_4 << ", " << Sc_1_4 << std::endl;
    std::cout << "block numbers for input: " << electronPMT[0] << "," << electronPMT[1] << std::endl;
    std::cout << "my block numbers for input: " << my_Sc_0_4 << "," << my_Sc_1_4 << std::endl;

    
    bool swap = false;
    if(Sc_0_4 == Sc_1_4)
    {
        std::cout << "BLOCK NUMBERS ARE THE SAME" << std::endl;
        std::cin.get();
    }
    else if(Sc_0_4 == my_Sc_1_4)
    {
        if(Sc_1_4 == my_Sc_0_4)
        {
            swap = true;
            std::cout << "SWAP DETECTED USING BLOCK NUMBERS" << std::endl;
            std::cin.get();

            std::swap(Sc_0_4, Sc_1_4); // this is irrelevant futher in function
            std::swap(Sc_0_8, Sc_1_8); // this needs to be swapped
            std::swap(trueElectronEnergy_0, trueElectronEnergy_1); // this also needs to be swapped
        }
    }

    
    if((std::abs(1.0e+03 * Sc_0_8 - electronEnergy[0]) > 1.0e-8) &&
       (std::abs(1.0e+03 * Sc_1_8 - electronEnergy[1]) > 1.0e-8) &&
       (std::abs(1.0e+03 * Sc_0_8 - electronEnergy[1]) < 1.0e-8) &&
       (std::abs(1.0e+03 * Sc_1_8 - electronEnergy[0]) < 1.0e-8))
    {
        // the data in the input chain is swapped relative to Summers values
        // swap all data associated with input chain
        // this includes the reconstructed scintillator hit energy
        // as well as the true momentum

        std::swap(Sc_0_8, Sc_1_8); // does not do anything outside this function
        std::swap(trueElectronEnergy_0, trueElectronEnergy_1); // this does do something outside of this function
        std::cout << "suspicious event fixed - this should now never occur" << std::endl;
    }
#endif

#define DEBUG 0
#if DEBUG
    std::cout << "original files:  trueEnergy(" << 1.0e+03 * trueElectronEnergy_0 << ", " << 1.0e+03 * trueElectronEnergy_1 << ") [MeV] " << std::endl;
    std::cout << "original files:  recoEnergy(" << 1.0e+03 * Sc_0_8 << ", " << 1.0e+03 * Sc_1_8 << ") [MeV] " << std::endl;
    std::cout << "processed files: recoEnergy(" << electronEnergy[0] << ", " << electronEnergy[1] << ") [MeV] " << std::endl;
    //std::cout << "Electron Energy: " << electronEnergy[0] << ", " << 1.0e3 * trueElectronEnergy_0 << std::endl;
    //std::cout << "Electron Energy: " << electronEnergy[1] << ", " << 1.0e3 * trueElectronEnergy_1 << std::endl;
    //std::cout << "Ecc: " << Ecc[0] << ", " << Ecc[1] << std::endl;
    //std::cout << "Sc: " << 1.0e+03 * Sc[0][8] << ", " << 1.0e+03 * Sc[1][8] << std::endl;
#endif

#if 0
    if((std::abs(1.0e+03 * Sc_0_8 - electronEnergy[0]) > 1.0e-3) ||
       (std::abs(1.0e+03 * Sc_1_8 - electronEnergy[1]) > 1.0e-3))
    {
        std::cout << "suspicious event: energy 0: " << 1.0e+03 * Sc_0_8 << ", " << electronEnergy[0] << std::endl;
        std::cout << "suspicious event: energy 1: " << 1.0e+03 * Sc_1_8 << ", " << electronEnergy[1] << std::endl;
    }

    if(
        (
            (trueElectronEnergy_0 < trueElectronEnergy_1) &&
            (electronEnergy[1] < electronEnergy[0])
        ) ||
        (
            (trueElectronEnergy_1 < trueElectronEnergy_0) &&
            (electronEnergy[0] < electronEnergy[1])
        )
    )
    {
        std::cout << "swap detected" << std::endl;
    }
#endif

#if 0
    std::cout << std::endl;
    // TODO: can simply set branch addresses to input variables?
    // not sure if multiplication is required
    trueElectronEnergy[0] = 1.0e+03 * trueElectronEnergy_0;
    trueElectronEnergy[1] = 1.0e+03 * trueElectronEnergy_1;
   
    std::cin.get();
#endif
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

void filestripperORIG()
{

    // QuickSort Test
    //qst();

    // datalist file:
    // /home/blotsd/NEMO3/Nd150_analysis/DataLists/nd150_61_rot_nd150.lst


    const int NUM_NAMES = 20;
    std::string names[NUM_NAMES];
    names[0] = "foils_nd150_61_rot_nd150_101_01.root";
    names[1] = "foils_nd150_61_rot_nd150_101_02.root";
    names[2] = "foils_nd150_61_rot_nd150_102_01.root";
    names[3] = "foils_nd150_61_rot_nd150_102_02.root";
    names[4] = "foils_nd150_61_rot_nd150_103_01.root";
    names[5] = "foils_nd150_61_rot_nd150_103_02.root";
    names[6] = "foils_nd150_61_rot_nd150_104_01.root";
    names[7] = "foils_nd150_61_rot_nd150_104_02.root";
    names[8] = "foils_nd150_61_rot_nd150_105_01.root";
    names[9] = "foils_nd150_61_rot_nd150_105_02.root";
    names[10] = "foils_nd150_61_rot_nd150_106_01.root";
    names[11] = "foils_nd150_61_rot_nd150_106_02.root";
    names[12] = "foils_nd150_61_rot_nd150_107_01.root";
    names[13] = "foils_nd150_61_rot_nd150_107_02.root";
    names[14] = "foils_nd150_61_rot_nd150_108_01.root";
    names[15] = "foils_nd150_61_rot_nd150_108_02.root";
    names[16] = "foils_nd150_61_rot_nd150_109_01.root";
    names[17] = "foils_nd150_61_rot_nd150_109_02.root";
    names[18] = "foils_nd150_61_rot_nd150_110_01.root";
    names[19] = "foils_nd150_61_rot_nd150_110_02.root";
////    names[0] = "foils_nd150_61_rot_nd150_1xx_xx.root";

    TChain *tchain = new TChain("h10", "h10");
    for(int i = 0; i < NUM_NAMES; ++ i)
    {
        std::cout << "add: i=" << i << " name[i]=" << names[i] << std::endl; 
        tchain->Add(names[i].c_str());
    }

//    TString tchain_dir = "/mnt/ramdisk/";
//    TString tchain_fname = "foils_nd150_61_rot_nd150_1xx_xx.root";
//    TChain *tchain = new TChain("h10", "h10");
//    tchain->Add(tchain_dir + tchain_fname);

    std::cout << "All files added in TChain" << std::endl;
    //std::cin.get();

    //TFile *finput = new TFile("/unix/nemo3/users/sblot/Nd150Analysis/newAnalysis/2e/betabeta/data_2e/Nd150_2eNg_output.root");
    //TFile *finput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2b2n_m4/Nd150_2eNg_output_sorted.root");
    //TFile *finput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2b2n_m4/Nd150_2eNg_output.root");
    //TFile *finput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/nd150/nd150_rot_2n2b_m4/Nd150_2eNg_output.root");
                              //unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2b2n_m4/Nd150_2eNg_output.root
    //TTree *tinput = (TTree*)finput->Get("Nd150_2eNg");

    // working at home
    
    //TFile *finput = new TFile("/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/Nd150_2eNg_output.root");
    TFile *finput = new TFile("/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/nd150/nd150_rot_2n2b_m4/Nd150_2eNg_output.root");

    TTree *tinput = (TTree*)finput->Get("Nd150_2eNg/Nd150_2eNg");

    //TFile *foutput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/betabeta/data_2e/Nd150_2eNg_output_truth.root", "recreate");
    //TFile *foutput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2b2n_m4/Nd150_2eNg_output_truth.root", "recreate");
    //TFile *foutput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/nd150/nd150_rot_2n2b_m4/Nd150_2eNg_output_truth.root", "recreate");

    //TFile *foutput = new TFile("/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/Nd150_2eNg_output_truth_NEW.root", "recreate");
    //TString foutput_dir = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/";

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




    
    Int_t Nsc;
    //Int_t Ngg;: Ngg/I                                                  *
    Float_t Sc[2000][12];
    //Float_t Gg;: Gg[Ngg][15]/F                                          *
    /*
    Int_t Nbr_tks;: Nbr_tksc/I                                             *
    Int_t Nbr_pts;: Nbr_ptsc[Nbr_tksc]/B                                   *
    Int_t Ind_points;: Ind_pointsc[Nbr_tksc][200]/b                         *
    Float_t Xcc;: Xcc[Nbr_tksc]/F                                        *
    Float_t Ycc;: Ycc[Nbr_tksc]/F                                        *
    Float_t Zcc;: Zcc[Nbr_tksc]/F                                        *
    Float_t Radcc;: Radcc[Nbr_tksc]/F                                      *
    Float_t Hcc;: Hcc[Nbr_tksc]/F                                        *
    Float_t Ecc;: Ecc[Nbr_tksc]/F                                        *
    */
    //Float_t Ecc[2];
    /*
    Float_t Decc;: Decc[Nbr_tksc]/F                                       *
    Float_t Qcc;: Qcc[Nbr_tksc]/F                                        *
    Float_t Prob_radc;: Prob_radcc[Nbr_tksc]/F                                *
    Float_t Prob_hc;: Prob_hcc[Nbr_tksc]/F                                   *
    Float_t Prob_helix;: Prob_helixc[Nbr_tksc]/F                              *
    Float_t X_foil;: X_foilc[Nbr_tksc]/F                                    *
    Float_t Y_foil;: Y_foilc[Nbr_tksc]/F                                    *
    Float_t Z_foil;: Z_foilc[Nbr_tksc]/F                                    *
    Float_t Cos_dir;: Cos_dirc[Nbr_tksc][3]/F                                *
    Float_t Xcvc;: Xcvc[Nbr_tksc]/F                                       *
    Float_t Ycvc;: Ycvc[Nbr_tksc]/F                                       *
    Float_t Zcvc;: Zcvc[Nbr_tksc]/F                                       *
    Float_t Radcv;: Radcvc[Nbr_tksc]/F                                     *
    Float_t Hcvc;: Hcvc[Nbr_tksc]/F                                       *
    Float_t Thminv;: Thminvc[Nbr_tksc]/F                                    *
    Float_t Thmaxv;: Thmaxvc[Nbr_tksc]/F                                    *
    Float_t X_foilv;: X_foilvc[Nbr_tksc]/F                                   *
    Float_t Y_foilv;: Y_foilvc[Nbr_tksc]/F                                   *
    Float_t Z_foilv;: Z_foilvc[Nbr_tksc]/F                                   *
    Float_t Cos_dirv;: Cos_dirvc[Nbr_tksc][3]/F                               *
    Float_t Qcvc;: Qcvc[Nbr_tksc]/F                                       *
    Float_t Xcsc;: Xcsc[Nbr_tksc]/F                                       *
    Float_t Ycsc;: Ycsc[Nbr_tksc]/F                                       *
    Float_t Zcsc;: Zcsc[Nbr_tksc]/F                                       *
    Float_t Radcs;: Radcsc[Nbr_tksc]/F                                     *
    Float_t Hcsc;: Hcsc[Nbr_tksc]/F                                       *
    Float_t Thmins;: Thminsc[Nbr_tksc]/F                                    *
    Float_t Thmaxs;: Thmaxsc[Nbr_tksc]/F                                    *
    Float_t Cos_dirs;: Cos_dirsc[Nbr_tksc][3]/F                               *
    Float_t Qcsc;: Qcsc[Nbr_tksc]/F                                       *
    Float_t X_scints;: X_scintsc[Nbr_tksc]/F                                  *
    Float_t Y_scints;: Y_scintsc[Nbr_tksc]/F                                  *
    Float_t Z_scints;: Z_scintsc[Nbr_tksc]/F                                  *
    Int_t Myieven;: Myievent/I                                             *
    Int_t Mynbr_tk;: Mynbr_tks/I                                            *
    Float_t X_scinti;: X_scintil[Mynbr_tks]/F                                 *
    Float_t Y_scinti;: Y_scintil[Mynbr_tks]/F                                 *
    Float_t Z_scinti;: Z_scintil[Mynbr_tks]/F                                 *
    Int_t Ind_scinti;: Ind_scintil[Mynbr_tks]/B                             *
    Int_t Myimpac;: Myimpact[Mynbr_tks]/B                                  *
    Float_t Xvert;: Xvert/F                                                *
    Float_t Yvert;: Yvert/F                                                *
    Float_t Zvert;: Zvert/F                                                *
    */
    Int_t run = 0;
//    Int_t date = 0;
//    Int_t time = 0;
///    Int_t evntime = 0;
    /*
    Float_t tau_sc_sav;: tau_sc_save/F                                        *
    */
    Int_t Nvntu = 0;
    Float_t Xvntu = 0.0;
    Float_t Yvntu = 0.0;
    Float_t Zvntu = 0.0;
    /*
    Float_t Tofvnt;: Tofvntu[Nvntu]/F                                       *
    */
    Int_t Ntntu = 0;
    //Float_t *Pxntu = nullptr;
    //Float_t *Pyntu = nullptr;
    //Float_t *Pzntu = nullptr;
    Float_t Pxntu[31];
    Float_t Pyntu[31];
    Float_t Pzntu[31];
    /*
    Float_t Toftnt;: Toftntu[Ntntu]/F                                       *
    Int_t Ivntu;: Ivntu[Ntntu]/b                                         *
    Int_t Ipntu;: Ipntu[Ntntu]/b                                         *
    */
    tchain->SetBranchAddress("Nsc", &Nsc);
    tchain->SetBranchAddress("Sc", Sc);
    tchain->SetBranchAddress("run", &run);
    //tchain->SetBranchAddress("date", &date);
    //tchain->SetBranchAddress("time", &time);
////    tchain->SetBranchAddress("evntime", &evntime);
    //tchain->SetBranchAddress("Ecc", &Ecc);
    tchain->SetBranchAddress("Nvntu", &Nvntu);
    tchain->SetBranchAddress("Xvntu", &Xvntu);
    tchain->SetBranchAddress("Yvntu", &Yvntu);
    tchain->SetBranchAddress("Zvntu", &Zvntu);
    tchain->SetBranchAddress("Ntntu", &Ntntu);
    tchain->SetBranchAddress("Pxntu", Pxntu);
    tchain->SetBranchAddress("Pyntu", Pyntu);
    tchain->SetBranchAddress("Pzntu", Pzntu);


/*
    // check if tchain is in order by event number
    // first check order
    Long64_t okcount{0};
    Long64_t badcount{0};
    Long64_t max_chain{tchain->GetEntries()};
    std::cout << "max_chain=" << max_chain << std::endl;
    tchain->GetEntry(0);
    Int_t run1, run2;
    run1 = -run;
    for(Long64_t ix_chain{1}; ix_chain < max_chain; ++ ix_chain)
    //for(Long64_t ix{0}; ix < max_chain - 1; ++ ix)
    {
        tchain->GetEntry(ix_chain);
        Int_t run2 = -run;

        if(run1 <= run2)
        {
            // ok
            ++ okcount;
        }
        else
        {
            std::cout << run1 << " " << run2 << std::endl;
            ++ badcount;
        }
    
        run1 = run2;
    }
    std::cout << "prog STOP." << std::endl;
    std::cout << "ok: " << okcount << " bad: " << badcount << std::endl;
    std::cin.get();
*/
    /*
    Int_t nElectrons;
    Double_t trueT1;
    Double_t trueT2;
    Double_t el_energy_[2];
    Double_t gen_weight;


    TFile *f_in = new TFile("NewElectronNtuplizerExe_Int_ManDB_output.root");
    TTree *t_in = (TTree*)f_in->Get("NewElectronNtuplizer/NewElectronNtuplizer");


    t_in->SetBranchAddress("nElectrons", &nElectrons);
    t_in->SetBranchAddress("trueT1", &trueT1);
    t_in->SetBranchAddress("trueT2", &trueT2);
    t_in->SetBranchAddress("el_energy_", el_energy_);


    TFile *f_out = new TFile("NewElectronNtuplizerExe_Int_ManDB_output_2e.root", "recreate");
    TDirectory *d_out = f_out->mkdir("NewElectronNtuplizer");
    f_out->cd("NewElectronNtuplizer");
    TTree *t_out = new TTree("NewElectronNtuplizer", "NewElectronNtuplizer");


    t_out->Branch("nElectrons", &nElectrons, "nElectrons/I");
    t_out->Branch("trueT1", &trueT1, "trueT1/D");
    t_out->Branch("trueT2", &trueT2, "trueT2/D");
    t_out->Branch("el_energy_", el_energy_, "el_energy[nElectrons]/D");

    */


    // check tchain run numbers in order
    // 2020-09-07: check passed on foils_nd150_61_rot_nd150_1xx_xx.root
    #if 0
    if(0)
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
                //std::cout << "run: " << run_last << " -> " << run << std::endl;
            }
            else
            {
                std::cout << "ERROR: run=" << run << " run_last=" << run_last << std::endl;
                //std::cin.get();
            }

            run_last = run;
        }
    }

    // check tinput run numbers in order TODO
    if(0)
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

    Long64_t count{0};
    Long64_t max{tinput->GetEntries()};

    //Long64_t max_chain{tchain->GetEntries()};

/*
    Long64_t max_chain{tchain->GetEntries()};
    //Long64_t ix_chain{0};
    Double_t *trueVertexZ_array = new Double_t[max_chain];
    for(Long64_t ix{0}; ix < max_chain; ++ ix)
    {
        tchain->GetEntry(ix);
        trueVertexZ_array[ix] = Zvntu;
    }
    // sort trueVertexZ_array
    QuickSort(trueVertexZ_array, 0, max_chain - 1);
    std::cout << "index built" << std::endl;

    // quick hack check - look for duplicates in trueVertexZ_array
    // first check order
    Long64_t okcount{0};
    Long64_t badcount{0};
    for(Long64_t ix{0}; ix < max_chain - 1; ++ ix)
    {
        if(trueVertexZ_array[ix] <= trueVertexZ_array[ix + 1])
        {
            // ok
            ++ okcount;
        }
        else
        {
            std::cout << trueVertexZ_array[ix] << " " << trueVertexZ_array[ix + 1] << std::endl;
            ++ badcount;
        }
    }
    std::cout << "prog STOP." << std::endl;
    std::cout << "ok: " << okcount << " bad: " << badcount << std::endl;
    std::cin.get();
*/

    //Long64_t ix_chain_start{0};

    // count number of matches to check for double matches
//    Long64_t multi_match_count{0};
    ///Long64_t no_match_start = -1;
    //Long64_t no_match_stop = -1;
    //Long64_t close_match_R_count{0};
    //Long64_t close_match_sector_count{0};
    //Long64_t close_match_Z_count{0};
    //Int_t Run_current;

    tinput->GetEntry(tinput->GetEntries() - 1);
    const Int_t Run_last = Run;
    tinput->GetEntry(0);
    const Int_t Run_first = Run;
    std::cout << "Run_first=" << Run_first << " Run_last=" << Run_last << std::endl;
    // list of file name extensions and start / end run numbers
    // extension    start       end
    // _0           Run_first   2500
    // _1           2500        3000
    // _2           3000        3500
    // _3           3500        4000
    //
    //
    //
    //Int_t Run_min = Run_first; // ... NEW_0
    //Int_t Run_max = 2500;
    
    //Int_t Run_min = 2500; // ... NEW_1
    //Int_t Run_max = 3000;

    //Int_t Run_min = 3000; // ... NEW_2
    //Int_t Run_max = 4000;

    //Int_t Run_min = 4000; // ... NEW_3
    //Int_t Run_max = 4500;

    //Int_t Run_min = 4500; // ... NEW_4
    //Int_t Run_max = 5000;

    //Int_t Run_min = 5000; // ... NEW_5
    //Int_t Run_max = 6000;

    //Int_t Run_min = 6000; // ... NEW_6
    //Int_t Run_max = 7000;

    //Int_t Run_min = 7000; // ... NEW_7
    //Int_t Run_max = 8000;

    //Int_t Run_min = 8000; // ... NEW_8
    //Int_t Run_max = 9000;

    //Int_t Run_min = 9000; // ... NEW_9
    //Int_t Run_max = Run_last;

    // TODO: change output file name

    //TString foutput_dir = "/mnt/ramdisk/";
    //TString foutput_dir = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/nd150/nd150_rot_2n2b_m4/";
    TString foutput_dir = "/mnt/ecb/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2eNg_29Sep2015/nd150/nd150_rot_2n2b_m4/";
    TString foutput_fname = "Nd150_2eNg_output_truth_NEW_all_3.root";

    std::cout << "last chance to check important parameters" << std::endl;
    //std::cout << "Run_min=" << Run_min << " Run_max=" << Run_max << std::endl;
    std::cout << "Output file directory: " << foutput_dir << std::endl;
    std::cout << "Output file name: " << foutput_fname << std::endl;
    //std::cin.get();



    ///////////////////////////////////////////////////////////////////////////
    // initialize output file
    ///////////////////////////////////////////////////////////////////////////

    TFile *foutput = new TFile(foutput_dir + foutput_fname, "recreate");
    TDirectory *doutput = foutput->mkdir("Nd150_2eNg");
    foutput->cd("Nd150_2eNg");
    TTree *toutput = new TTree("Nd150_2eNg", "Nd150_2eNg");

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
    // start main loop
    ///////////////////////////////////////////////////////////////////////////

    //std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();

    Long64_t ix_chain_faster = 0;
    for(Long64_t ix{0}; ix < max; ++ ix)
    {
        tinput->GetEntry(ix);

        if(ix % 10000 == 0)
        {
            std::cout << "ix=" << ix << " / " << max << std::endl;
        }

        /*
        if(Run < Run_min)
        {
            continue;
        }
        else if(Run >= Run_max)
        {
            break;
        }

        if(ix % 1 == 0)
        {
            std::cout << "ix=" << ix << " / " << max << std::endl;
            std::cout << "looking for Run=" << Run << std::endl;

            std::chrono::system_clock::time_point current_time = std::chrono::high_resolution_clock::now();
            //auto runtime_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);
            //auto runtime_hours = std::chrono::duration_cast<std::chrono::hours>(current_time - start_time);
            //std::cout << "Runtime: " << runtime_seconds.count() << " s, " << runtime_hours.count() << " h" << std::endl;
            std::chrono::duration<double> runtime_seconds = current_time - start_time;
            std::cout << "Runtime: " << runtime_seconds.count() << " s, " << runtime_seconds.count() / 3600.0 << " h" << std::endl;
            Int_t total_runs = Run_max - Run_min;
            Int_t completed_runs = Run - Run_min;
            //Int_t todo_runs = total_runs - current_runs;
            double fraction_complete = (double)completed_runs / (double)total_runs;
            std::cout << "Percentage complete: " << 100.0 * fraction_complete << " %" << std::endl;
            try
            {
                double estimated_total_time_seconds = runtime_seconds.count() / fraction_complete;
                double estimated_time_remaining_seconds = estimated_total_time_seconds - runtime_seconds.count();
                double estimated_time_remaining_hours = estimated_time_remaining_seconds / 3600.0;
                std::cout << "Estimated time to completion: " << estimated_time_remaining_seconds << " s, " << estimated_time_remaining_hours << " h" << std::endl;
            }
            catch(...)
            {
                // do nothing
            }
        }
        */

        //std::cout << "Searching, ix=" << ix << std::endl;
//        multi_match_count = 0;

        
        // new 2019-01-14
        //Run_current = Run;
        //std::string expression("run==");
        //expression += std::to_string(-Run);
        //std::cout << expression << std::endl;
        //TTree *tchain_filtered = tchain->CopyTree(expression.c_str());
        //std::cout << tchain_filtered->GetEntries() << std::endl;
        //std::cin.get();
        Long64_t max_chain{tchain->GetEntries()};
        Long64_t ix_chain_start, ix_chain_end, ix_chain_match;
        //std::cout << "ix_chain_faster=" << ix_chain_faster << std::endl;
        for(Long64_t ix_chain{ix_chain_faster}; ix_chain < max_chain; ++ ix_chain)
        //for(Long64_t ix_chain{0}; ix_chain < max_chain; ++ ix_chain)
        {
            tchain->GetEntry(ix_chain);
            //std::cout << "ix_chain=" << ix_chain << " run=" << run << std::endl;
            //std::cin.get();
            
            if(-run == Run)
            {
                ix_chain_start = ix_chain;
                ix_chain_end = ix_chain + 1;
                ix_chain_faster = ix_chain; // 2020-09-11 moved this to here from below
                //std::cout << "ix_chain_start=" << ix_chain_start << std::endl;
                break;
            }
        }
        //std::cout << "setting ix_chain_start=" << ix_chain_start << std::endl;
        for(Long64_t ix_chain{ix_chain_start + 1}; ix_chain < max_chain; ++ ix_chain)
        {
            tchain->GetEntry(ix_chain);
            //std::cout << "ix_chain=" << ix_chain << " run=" << run << std::endl;
            //std::cin.get();

            if(-run == Run)
            {
                //ix_chain_end = ix_chain;
                ix_chain_end = ix_chain + 1; // does not include the final + 1 in search
            }
            else
            {
                //ix_chain_end = ix_chain; // moved above
                //ix_chain_faster = ix_chain;
                //std::cout << "ix_chain_end=" << ix_chain_end << " ix_chain_faster=" << ix_chain_faster << std::endl;
                break;
            }
        }
        //std::cout << "setting ix_chain_end=" << ix_chain_end << std::endl;

        //std::cout << "ix=" << ix << ", ix_chain_start=" << ix_chain_start << " ix_chain_end=" << ix_chain_end << std::endl;
        //std::cin.get();

        // ix_chain_start and ix_chain_end are now set
        Long64_t multi_match_count = 0;
        for(Long64_t ix_chain{ix_chain_start}; ix_chain < ix_chain_end; ++ ix_chain)
        {
            tchain->GetEntry(ix_chain); 
        
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
                Long64_t ix_chain)
            */
            SearchFunction(return_flag, Run, trueVertexR, trueVertexSector, trueVertexZ, run, Xvntu, Yvntu, Zvntu, ix, ix_chain);
            if(return_flag == 1)
            {
                if(multi_match_count == 0)
                {
                    ix_chain_match = ix_chain;
                }
                multi_match_count ++;
            }

            //void GetTrueEnergy(Double_t* trueElectronEnergy, Float_t* Pxntu, Float_t* Pyntu, Float_t* Pzntu)
            //void SearchFunction(int &return_flag, Int_t Run, Doublt_t trueVertexR, Double_t trueVertexSector, Double_t trueVertexZ, Double_t* trueElectronEnergy, Int_t run, Float_t Xvntu, Float_t Yvntu, Float_t Zvntu)
        }

        if(multi_match_count == 1)
        {
            //std::cout << "match for ix=" << ix << " found: ix_chain_match=" << ix_chain_match << std::endl;
            tchain->GetEntry(ix_chain_match);

            /*
            Double_t rntu = std::sqrt(Xvntu * Xvntu + Yvntu * Yvntu);
            Double_t thetantu = std::atan2(Yvntu, Xvntu);
            Double_t sectorntu = (2.5 / std::atan(1.0)) * thetantu; // / 18.0;
            std::cout << "ix_chain_match=" << ix_chain_match << " ix=" << ix << " close match for trueVertexR; difference: " << std::abs(trueVertexR - rntu) << std::endl;
            std::cout << "ix_chain_match=" << ix_chain_match << " ix=" << ix << " close match for trueVertexSector; difference: " << std::abs(trueVertexSector - sectorntu) << std::endl;
            std::cout << "ix_chain_match=" << ix_chain_match << " ix=" << ix << " close match for trueVertexZ; difference: " << std::abs(trueVertexZ - Zvntu) << std::endl;
            */

            // write to file
            #if 0
            GetTrueEnergy(trueElectronEnergy,
                Pxntu, Pyntu, Pzntu,
                electronEnergy,
                Sc[0][1], Sc[1][1],
                Sc[0][2], Sc[1][2],
                Sc[0][3], Sc[1][3],
                Sc[0][4], Sc[1][4],
                Sc[0][8], Sc[1][8],
                electronPMT);
            #endif
            GetTrueEnergy(trueElectronEnergy,
                Pxntu, Pyntu, Pzntu);
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


/*
        // get index of match using binary search
        Long64_t L = 0;
        Long64_t R = max_chain - 1;
        Double_t T = trueVertexZ;
        bool success = false;
        Long64_t return_index = -1;
        while(L <= R)
        {
            Long64_t m = (L + R) / 2;
            if(trueVertexZ_array[m] == T)
            {
                success = true;
                return_index = m;
                break;
            }
            else if(trueVertexZ_array[m] < T)
            {
                L = m + 1;
            }
            else if(trueVertexZ_array[m] > T)
            {
                R = m - 1;
            }
        }
        //success = false;
        Long64_t ix_chain_0 = 0;
        Long64_t ix_chain_1 = 0;
        if(success == true)
        {
            ix_chain_0 = return_index;
            while(trueVertexZ_array[ix_chain_0] == T)
            {
                ix_chain_0 --;
            }
            while(trueVertexZ_array[ix_chain_1] == T)
            {
                ix_chain_1 ++;
            }
            ix_chain_1 = ix_chain_0 + 1;
            // this just makes for loop below work
        }
*/

#if 0
        Long64_t max_chain{tchain->GetEntries()};
        // search for corresponding entry in tchain ttree
        std::cout << "starting search at ix_chain_start=" << ix_chain_start << std::endl;
        for(Long64_t ix_chain{ix_chain_start}; ix_chain < max_chain; ++ ix_chain)
        //for(Long64_t ix_chain{ix_chain_0}; ix_chain < ix_chain_1; ++ ix_chain)
        {
            // set resume search point to be the next event
            //ix_chain_start = ix_chain + 1;
            // not sure if this is the correct place for this statement, so moving
            // it inside the if condition, where I know it will work but will not be
            // optimal
            
            tchain->GetEntry(ix_chain);

            Double_t thetantu = std::atan2(Yvntu, Xvntu);
            Double_t sectorntu = (2.5 / std::atan(1.0)) * thetantu; // / 18.0;
            Double_t rntu = std::sqrt(Xvntu * Xvntu + Yvntu * Yvntu);

            // match variables
            if(-run == Run)
            {
                //ix_chain_start = ix_chain;

                // if run matches, check if event is a complete match

                bool close_match_R = false;
                bool close_match_sector = false;
                bool close_match_Z = false;
                int close_match_count = 0;
                if(std::abs(trueVertexR - rntu) < 1.0e-4)
                //if(trueVertexR == rntu)
                {
                    //std::cout << "close match for R" << std::endl;
                    //close_match = true;
                    close_match_R = true;
                    //++ close_match_R_count;
                    ++ close_match_count;
                }
                if(std::abs(trueVertexSector - sectorntu) < 1.0e-6)
                //if(trueVertexSector == sectorntu)
                {
                    //std::cout << "close match for sector" << std::endl;
                    //close_match = true;
                    close_match_sector = true;
                    //++ close_match_sector_count;
                    ++ close_match_count;
                }
                //if(std::abs(trueVertexZ - Zvntu) < 1.0e-5)
                if(trueVertexZ == Zvntu)
                {
                    //std::cout << "close match for Z" << std::endl;
                    //close_match = true;
                    close_match_Z = true;
                    //++ close_match_Z_count;
                    ++ close_match_count;
                }
                //if(close_match == true)
                if(close_match_count == 3)
                {
                    /***
                    if(no_match_start != -1)
                    {
                        std::cout << "No matching entry for ix range: " << no_match_start << " < ix < " << ix - 1 << std::endl;
                        no_match_start = -1;
                    }
                   ***/
#define DEBUG 0
#if DEBUG
                    std::cout << "ix=" << ix << " -> ix_chain=" << ix_chain << std::endl;
                    std::cout << "delta: r=" << std::abs(trueVertexR - rntu) << ", s=" << std::abs(trueVertexSector - sectorntu) << ", Z=" << std::abs(trueVertexZ - Zvntu) << std::endl;
#endif
#if DEBUG
                    std::cout << "Found corresponding event: ix=" << ix << " ix_chain=" << ix_chain << std::endl;

                    std::cout << "ix=" << ix << " ix_chain=" << ix_chain << std::endl;
                    std::cout << "-run=" << -run << " Run=" << Run << std::endl;
                    //std::cout << "evntime=" << evntime << " Event=" << Event << std::endl;
                    std::cout << "original files:  trueVertex(" << Xvntu << "," << Yvntu << "," << Zvntu << ")" << " [cartesian] " << std::endl;
                    std::cout << "original files:  trueVertex(" << rntu << "," << sectorntu << "," << Zvntu << ")" << " [cylindrical] " << std::endl;
                    std::cout << "processed files: trueVertex(" << trueVertexR << "," << trueVertexSector  << "," << trueVertexZ << ")" << " [cylindrical] " << std::endl;
                    std::cout << "other variables..." << std::endl;
                    std::cout << "Event Time: " << eventTime << ", " << evntime << ", " << time << std::endl;
                    std::cout << "Ntntu=" << Ntntu << std::endl;
                    //std::cout << "Pxntu=" << Pxntu << std::endl;
                    //std::cout << "Pyntu=" << Pyntu << std::endl;
                    //std::cout << "Pzntu=" << Pzntu << std::endl;
                    // Note: think these might be in GeVs
                    std::cout << "Pxntu[0]=" << Pxntu[0] << std::endl;
                    std::cout << "Pyntu[0]=" << Pyntu[0] << std::endl;
                    std::cout << "Pzntu[0]=" << Pzntu[0] << std::endl;
                    std::cout << "Pxntu[1]=" << Pxntu[1] << std::endl;
                    std::cout << "Pyntu[1]=" << Pyntu[1] << std::endl;
                    std::cout << "Pzntu[1]=" << Pzntu[1] << std::endl;
#endif
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
                    if(electron_energy_0 < electron_energy_1) std::swap(electron_energy_0, electron_energy_1);
#if DEBUG
                    std::cout << "original files:  trueEnergy(" << 1.0e+03 * electron_energy_0 << ", " << 1.0e+03 * electron_energy_1 << ") [MeV] " << std::endl;
                    std::cout << "original files:  recoEnergy(" << 1.0e+03 * Sc[0][8] << ", " << 1.0e+03 * Sc[1][8] << ") [MeV] " << std::endl;
                    std::cout << "processed files: recoEnergy(" << electronEnergy[0] << ", " << electronEnergy[1] << ") [MeV] " << std::endl;
                    //std::cout << "Electron Energy: " << electronEnergy[0] << ", " << 1.0e3 * electron_energy_0 << std::endl;
                    //std::cout << "Electron Energy: " << electronEnergy[1] << ", " << 1.0e3 * electron_energy_1 << std::endl;
                    //std::cout << "Ecc: " << Ecc[0] << ", " << Ecc[1] << std::endl;
                    //std::cout << "Sc: " << 1.0e+03 * Sc[0][8] << ", " << 1.0e+03 * Sc[1][8] << std::endl;
#endif

                    //std::cin.get();

                    // TODO: can simply set branch addresses to input variables?
                    // not sure if multiplication is required
                    trueElectronEnergy[0] = 1.0e+03 * electron_energy_0;
                    trueElectronEnergy[1] = 1.0e+03 * electron_energy_1;

                    if(multi_match_count == 0)
                    {
                        std::cout << "match for ix=" << ix << " found: ix_chain=" << ix_chain << std::endl;
                        ix_chain_start = ix_chain;
                        // write to file
                        ++ count;
                        toutput->Fill();
                    }

                    // NOTE: enable these two lines to make code work, but not check for duplicate events
                    ++ ix_chain;
                    break;
                    ++ multi_match_count;
                }
            }
            else // if(run == -Run)
            {
                // run variable does not match
                // if we already found a match, break here
                // we need to start looking for next event in
                // tinput
                if(multi_match_count > 0)
                {
                    // NOTE: this broke my code, because there can be events with repeated run numbers in tinput
                    //ix_chain_start = ix_chain;
                    //std::cout << "ix_chain=" << ix_chain << " ix_chain_start set to " << ix_chain_start << " multi_match_count=" << multi_match_count << std::endl;
                    break;
                }
                // else is no match was found, do nothing?
                else
                {
                    //std::cout << "no match found for ix=" << ix << " ix_chain=" << ix_chain << std::endl;
                }
            }
            

            //++ ix_chain;

            /*
            if(ix_chain % 1000000 == 0)
            {
                std::cout << "ix_chain=" << ix_chain << std::endl;
            }
            */
        }

        if(multi_match_count > 1)
        {
            std::cout << "ix=" << ix << " multiple matches found in data" << std::endl;
        }
#endif


    }
    //if(no_match_start != -1)
    //{
    //    std::cout << "No matching entry for ix range: " << no_match_start << " < ix < " << max - 1 << std::endl;
    //}
    //std::cout << "close_match_R_count=" << close_match_R_count << std::endl;
    //std::cout << "close_match_sector_count=" << close_match_sector_count << std::endl;
    //std::cout << "close_match_Z_count=" << close_match_Z_count << std::endl;
    std::cout << "done" << std::endl;
    std::cout << "swap_count_true=" << swap_count_true << std::endl;
    std::cout << "swap_count_false=" << swap_count_false << std::endl;

    toutput->Write();
    foutput->Close();

    finput->Close();
    //tchain->Close();

    
    std::cout << count << "/" << max << std::endl;;
    

    /*delete [] trueVertexZ_array;*/

}
