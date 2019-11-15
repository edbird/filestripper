

#include <iostream>

// Root headers
#include "TFile.h"
#include "TTree.h"


void filestripper()
{

    TChain *tchain = new TChain("h10", "h10");
    tchain->Add("foils_nd150_61_rot_nd150_101_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_101_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_102_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_102_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_103_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_103_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_104_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_104_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_105_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_105_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_106_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_106_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_107_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_107_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_108_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_108_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_109_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_109_02.root")
    tchain->Add("foils_nd150_61_rot_nd150_110_01.root")
    tchain->Add("foils_nd150_61_rot_nd150_110_02.root")

    TFile *finput = new Tfile("/unix/nemo3/users/sblot/Nd150Analysis/newAnalysis/2e/betabeta/data_2e/Nd150_2eNg_output.root")
    TFile *foutput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/betabeta/data_2e/Nd150_2eNg_output_truth.root");

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


    Long64_t count{0};
    Long64_t max{t_in->GetEntries()};
    for(Long64_t ix{0}; ix < t_in->GetEntries(); ++ ix)
    {

        t_in->GetEntry(ix);
        
        // analysis only valid for 2 electron events
        if(nElectrons != 2) continue;
    
        ++ count;

        t_out->Fill();
    
    }

    t_out->Write();
    f_out->Close();

    f_in->Close();

    std::cout << count << "/" << max << std::endl;;
    */

    finput->Close();
    foutput->Close();

}
