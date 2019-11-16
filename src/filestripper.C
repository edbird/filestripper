

#include <iostream>

// Root headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"


void filestripper()
{

    TChain *tchain = new TChain("h10", "h10");
    tchain->Add("foils_nd150_61_rot_nd150_101_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_101_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_102_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_102_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_103_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_103_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_104_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_104_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_105_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_105_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_106_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_106_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_107_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_107_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_108_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_108_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_109_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_109_02.root");
    tchain->Add("foils_nd150_61_rot_nd150_110_01.root");
    tchain->Add("foils_nd150_61_rot_nd150_110_02.root");

    TFile *finput = new TFile("/unix/nemo3/users/sblot/Nd150Analysis/newAnalysis/2e/betabeta/data_2e/Nd150_2eNg_output.root");
    TTree *tinput = (TTree*)finput->Get("Nd150_2eNg/Nd150_2eNg");

    TFile *foutput = new TFile("/unix/nemo3/users/ebirdsall/Nd150Analysis/newAnalysis/2e/betabeta/data_2e/Nd150_2eNg_output_truth.root", "recreate");
    TDirectory *doutput = foutput->mkdir("Nd150_2eNg");
    foutput->cd("Nd150_2eNg");
    TTree *toutput = new TTree("Nd150_2eNg", "Nd150_2eNg");

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
    Int_t NAPromptGgHitsSide = 0;
    Double_t *NAPromptGgHitsDist2Vertex = NULL; //[NAPromptGgHits];
    Double_t *NAPromptGgHitsDist2Calo = NULL; //[NAPromptGgHits];
    Int_t nGammaClusters = 0;
    Int_t *nInCluster = NULL; //[nGammaClusters];
    Double_t *clusterEnergy = NULL; //[nGammaClusters];
    Double_t *clusterTimeSpan = NULL; //[nGammaClusters];
    Int_t nTotalClusterHits = 0;
    Double_t *clusterHitEnergy = NULL; //[nTotalClusterHits];
    Int_t *clusterHitPMT = NULL; //[nTotalClusterHits];
    Int_t *clusterHitLDFlag = NULL; //[nTotalClusterHits];
    Double_t *clusterHitLDCorr = NULL; //[nTotalClusterHits];
    Double_t *clusterHitLDCorrErr = NULL; //[nTotalClusterHits];
    Double_t *clusterHitSec = NULL; //[nTotalClusterHits];
    Double_t *clusterHitZ = NULL; //[nTotalClusterHits];

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



    /*
    Int_t Nsc;: Nsc/I                                                  *
    Int_t Ngg;: Ngg/I                                                  *
    Float_t Sc;: Sc[Nsc][12]/F                                          *
    Float_t Gg;: Gg[Ngg][15]/F                                          *
    Int_t Nbr_tks;: Nbr_tksc/I                                             *
    Int_t Nbr_pts;: Nbr_ptsc[Nbr_tksc]/B                                   *
    Int_t Ind_points;: Ind_pointsc[Nbr_tksc][200]/b                         *
    Float_t Xcc;: Xcc[Nbr_tksc]/F                                        *
    Float_t Ycc;: Ycc[Nbr_tksc]/F                                        *
    Float_t Zcc;: Zcc[Nbr_tksc]/F                                        *
    Float_t Radcc;: Radcc[Nbr_tksc]/F                                      *
    Float_t Hcc;: Hcc[Nbr_tksc]/F                                        *
    Float_t Ecc;: Ecc[Nbr_tksc]/F                                        *
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
    Int_t date = 0;
    Int_t time = 0;
    Int_t evntime = 0;
    /*
    Float_t tau_sc_sav;: tau_sc_save/F                                        *
    Int_t Nvntu;: Nvntu/I                                                *
    Float_t Xvntu;: Xvntu[Nvntu]/F                                         *
    Float_t Yvntu;: Yvntu[Nvntu]/F                                         *
    Float_t Zvntu;: Zvntu[Nvntu]/F                                         *
    Float_t Tofvnt;: Tofvntu[Nvntu]/F                                       *
    */
    Int_t Ntntu = 0;
    Float_t *Pxntu = nullptr;
    Float_t *Pyntu = nullptr;
    Float_t *Pzntu = nullptr;
    /*
    Float_t Toftnt;: Toftntu[Ntntu]/F                                       *
    Int_t Ivntu;: Ivntu[Ntntu]/b                                         *
    Int_t Ipntu;: Ipntu[Ntntu]/b                                         *
    */
    tchain->SetBranchAddress("run", &run);
    //tchain->SetBranchAddress("date", &date);
    //tchain->SetBranchAddress("time", &time);
    tchain->SetBranchAddress("evntime", &evntime);
    tchain->SetBranchAddress("Ntntu", &Ntntu);
    tchain->SetBranchAddress("Pxntu", &Pxntu);
    tchain->SetBranchAddress("Pyntu", &Pyntu);
    tchain->SetBranchAddress("Pzntu", &Pzntu);

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

    Long64_t count{0};
    Long64_t max{tinput->GetEntries()};
    Long64_t max_chain{tchain->GetEntries()};
    Long64_t ix_chain{0};
    for(Long64_t ix{0}; ix < max; ++ ix)
    {

        tinput->GetEntry(ix);
        
        // search for corresponding entry in tchain ttree
        for(;;)
        {
            tchain->GetEntry(ix_chain);

            // match variables
            if(run == Run)
            {
                if(evntime == Event)
                {
                    std::cout << "Found corresponding event: ix=" << ix << " ix_chain=" << ix_chain << std::endl;
                    ++ ix_chain;
                    break;
                }
            }

            ++ ix_chain;
        }


        ++ count;
        toutput->Fill();
    }

    toutput->Write();
    foutput->Close();

    finput->Close();
    //tchain->Close();

    /*
    std::cout << count << "/" << max << std::endl;;
    */


}
