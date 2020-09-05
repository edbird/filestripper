

#include <iostream>

// Root headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"



void makechainsmall()
{

    // QuickSort Test
    //qst();

    // datalist file:
    // /home/blotsd/NEMO3/Nd150_analysis/DataLists/nd150_61_rot_nd150.lst

    TString inputdir = "/mnt/ecb/unix/nemo2/reco/reco_20130206_fix/mc_8.0/rot/";
    TString outputdir = "/mnt/ecb/unix/nemo2/reco/reco_20130206_fix/mc_8.0/rot/";

    std::string names[20];
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

    std::string namesout[20];


    Int_t run = 0;

    Int_t Nsc;
    Float_t Sc[2000][12];
     
    Int_t Nvntu = 0;
    Float_t Xvntu[31];
    Float_t Yvntu[31];
    Float_t Zvntu[31];

    Int_t Ntntu = 0;
    Float_t Pxntu[31];
    Float_t Pyntu[31];
    Float_t Pzntu[31];



    TFile *finput = nullptr;
    TFile *foutput = nullptr;
    TTree *tinput = nullptr;
    TTree *toutput = nullptr;
    for(int i = 0; i < 20; ++ i)
    {
        std::cout << "i=" << i << std::endl;

        std::string name = names[i];
        std::string nameout = name.substr(0, name.find(".")) + "_small.root";
        namesout[i] = nameout;

        std::cout << "Reading: " << name << std::endl;
        std::cout << "Writing: " << nameout << std::endl;

        ///////////////////////////////////////////////////////////////////////

        finput = new TFile(inputdir + TString(name.c_str()));
        tinput = (TTree*)finput->Get("h10");

        tinput->SetBranchAddress("run", &run);

        tinput->SetBranchAddress("Nsc", &Nsc);
        tinput->SetBranchAddress("Sc", Sc);
        
        tinput->SetBranchAddress("Nvntu", &Nvntu); // this is always 1
        tinput->SetBranchAddress("Xvntu", Xvntu);
        tinput->SetBranchAddress("Yvntu", Yvntu);
        tinput->SetBranchAddress("Zvntu", Zvntu);

        tinput->SetBranchAddress("Ntntu", &Ntntu); // this is always 2
        tinput->SetBranchAddress("Pxntu", Pxntu);
        tinput->SetBranchAddress("Pyntu", Pyntu);
        tinput->SetBranchAddress("Pzntu", Pzntu);

        ///////////////////////////////////////////////////////////////////////

        foutput = new TFile(outputdir + TString(nameout.c_str()), "RECREATE");
        toutput = new TTree("h10", "h10");

        toutput->Branch("run", &run, "run/I");

        toutput->Branch("Nsc", &Nsc, "NSc/I");
        toutput->Branch("Sc", Sc, "Sc[NSc][12]/F");
        
        toutput->Branch("Nvntu", &Nvntu, "Nvntu/I"); // this is always 1
        toutput->Branch("Xvntu", Xvntu, "Xvntu[Nvntu]/F");
        toutput->Branch("Yvntu", Yvntu, "Yvntu[Nvntu]/F");
        toutput->Branch("Zvntu", Zvntu, "Zvntu[Nvntu]/F");

        toutput->Branch("Ntntu", &Ntntu, "Ntntu/I"); // this is always 2
        toutput->Branch("Pxntu", Pxntu, "Pxntu[Ntntu]/F");
        toutput->Branch("Pyntu", Pyntu, "Pyntu[Ntntu]/F");
        toutput->Branch("Pzntu", Pzntu, "Pzntu[Ntntu]/F");
        
        Long64_t count{tinput->GetEntries()};
        for(Long64_t ix{0}; ix < count; ++ ix)
        {
            tinput->GetEntry(ix);
            toutput->Fill();
        }

        toutput->Write();
        foutput->Close();

        finput->Close();
    }



    std::cout << "done" << std::endl;

}
