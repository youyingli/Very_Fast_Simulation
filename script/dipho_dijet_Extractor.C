/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l jetExtractor.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------
struct InputFeature
{

    double dipho_mass;
    double dipho_PToM;

    double leadPho_PToM        ;
    double sublPho_PToM        ;
    double dipho_LeadPhoPt     ;
    double dipho_LeadPhoEta    ;
    double dipho_LeadPhoPhi    ;
    double dipho_LeadPhoE      ;
    double dipho_SubPhoPt      ;
    double dipho_SubPhoEta     ;
    double dipho_SubPhoPhi     ;
    double dipho_SubPhoE       ;
    double dijet_LeadJPt       ;
    double dijet_LeadJEta      ;
    double dijet_LeadJPhi      ;
    double dijet_LeadJE        ;
    double dijet_SubJPt        ;
    double dijet_SubJEta       ;
    double dijet_SubJPhi       ;
    double dijet_SubJE         ;
    double dijet_abs_dEta      ;
    double dijet_Mjj           ;
    double dijet_centrality_gg ;
    double dijet_dphi          ;
    double dijet_minDRJetPho   ;
    double dijet_dipho_dphi_trunc    ;

    double event_pt;
    double event_eta;
    double event_phi;
    double event_mass;
    double event_energy;


};

//------------------------------------------------------------------------------
struct InputFeature_PF
{

    double PF_eta;
    double PF_phi;

    double PF_pt;
    double PF_energy;
    double PF_rel_pt;
    double PF_rel_energy;
    double PF_deltaEta;
    double PF_deltaPhi;
    double PF_charge;

    double PF_isEle;
    double PF_isMu;
    double PF_isChargedHad;
    double PF_isGamma;
    double PF_isNeutralHad;

    double PF_tanhd0;
    double PF_tanhdz;
    double PF_sigmad0;
    double PF_sigmadz;

};

//------------------------------------------------------------------------------

struct pt_sortor {
    bool operator()( const InputFeature_PF &pf1, const InputFeature_PF &pf2 ) {
        return ( pf1.PF_pt > pf2.PF_pt );
    }
};

//------------------------------------------------------------------------------
TLorentzVector UseEtaPhiEMass(double eta, double phi, double energy, double mass)
{
    double p = energy - mass > 0 ? TMath::Sqrt(energy*energy - mass*mass) : 0.0000001;

    double pt = p / TMath::CosH(eta);

    TLorentzVector p4;
    p4.SetPtEtaPhiM(pt, eta, phi, mass);

    return p4;
}

//------------------------------------------------------------------------------

void dipho_dijet_Extractor(const char *inputFile)
{
    gSystem->Load("libDelphes");

    InputFeature dipho_dijet_vars;

    std::vector<double> PF_eta          ; 
    std::vector<double> PF_phi          ; 
    std::vector<double> PF_pt           ; 
    std::vector<double> PF_energy       ; 
    std::vector<double> PF_rel_pt       ; 
    std::vector<double> PF_rel_energy   ; 
    std::vector<double> PF_deltaEta     ; 
    std::vector<double> PF_deltaPhi     ; 
    std::vector<double> PF_charge       ; 
    std::vector<double> PF_isEle        ; 
    std::vector<double> PF_isMu         ; 
    std::vector<double> PF_isChargedHad ; 
    std::vector<double> PF_isGamma      ; 
    std::vector<double> PF_isNeutralHad ; 
    std::vector<double> PF_tanhd0       ; 
    std::vector<double> PF_tanhdz       ; 
    std::vector<double> PF_sigmad0      ; 
    std::vector<double> PF_sigmadz      ; 

    std::vector<double> Jet1_PF_eta          ; 
    std::vector<double> Jet1_PF_phi          ; 
    std::vector<double> Jet1_PF_pt           ; 
    std::vector<double> Jet1_PF_energy       ; 
    std::vector<double> Jet1_PF_rel_pt       ; 
    std::vector<double> Jet1_PF_rel_energy   ; 
    std::vector<double> Jet1_PF_deltaEta     ; 
    std::vector<double> Jet1_PF_deltaPhi     ; 
    std::vector<double> Jet1_PF_charge       ; 
    std::vector<double> Jet1_PF_isEle        ; 
    std::vector<double> Jet1_PF_isMu         ; 
    std::vector<double> Jet1_PF_isChargedHad ; 
    std::vector<double> Jet1_PF_isGamma      ; 
    std::vector<double> Jet1_PF_isNeutralHad ; 
    std::vector<double> Jet1_PF_tanhd0       ; 
    std::vector<double> Jet1_PF_tanhdz       ; 
    std::vector<double> Jet1_PF_sigmad0      ; 
    std::vector<double> Jet1_PF_sigmadz      ; 

    std::vector<double> Jet2_PF_eta          ; 
    std::vector<double> Jet2_PF_phi          ; 
    std::vector<double> Jet2_PF_pt           ; 
    std::vector<double> Jet2_PF_energy       ; 
    std::vector<double> Jet2_PF_rel_pt       ; 
    std::vector<double> Jet2_PF_rel_energy   ; 
    std::vector<double> Jet2_PF_deltaEta     ; 
    std::vector<double> Jet2_PF_deltaPhi     ; 
    std::vector<double> Jet2_PF_charge       ; 
    std::vector<double> Jet2_PF_isEle        ; 
    std::vector<double> Jet2_PF_isMu         ; 
    std::vector<double> Jet2_PF_isChargedHad ; 
    std::vector<double> Jet2_PF_isGamma      ; 
    std::vector<double> Jet2_PF_isNeutralHad ; 
    std::vector<double> Jet2_PF_tanhd0       ; 
    std::vector<double> Jet2_PF_tanhdz       ; 
    std::vector<double> Jet2_PF_sigmad0      ; 
    std::vector<double> Jet2_PF_sigmadz      ; 

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");

    TFile* file = TFile::Open("output_dipho_dijet.root", "recreate");
    TTree* tree = new TTree("dipho_dijet", "");

    TH1F*h1_monitor_njet = new TH1F("h1_monitor_njet", "", 10, 0, 10);
    TH1F*h1_monitor_npf  = new TH1F("h1_monitor_npf", "", 50, 0, 300);

    tree->Branch( "dipho_mass"              , &dipho_dijet_vars.dipho_mass             , "dipho_mass/D"            ); 
    tree->Branch( "dipho_PToM"              , &dipho_dijet_vars.dipho_PToM             , "dipho_PToM/D"            ); 
    tree->Branch( "leadPho_PToM"            , &dipho_dijet_vars.leadPho_PToM           , "leadPho_PToM/D"          ); 
    tree->Branch( "sublPho_PToM"            , &dipho_dijet_vars.sublPho_PToM           , "sublPho_PToM/D"          );
    tree->Branch( "dipho_LeadPhoPt"         , &dipho_dijet_vars.dipho_LeadPhoPt        , "dipho_LeadPhoPt/D"       ); 
    tree->Branch( "dipho_LeadPhoEta"        , &dipho_dijet_vars.dipho_LeadPhoEta       , "dipho_LeadPhoEta/D"      ); 
    tree->Branch( "dipho_LeadPhoPhi"        , &dipho_dijet_vars.dipho_LeadPhoPhi       , "dipho_LeadPhoPhi/D"      ); 
    tree->Branch( "dipho_LeadPhoE"          , &dipho_dijet_vars.dipho_LeadPhoE         , "dipho_LeadPhoE/D"        ); 
    tree->Branch( "dipho_SubPhoPt"          , &dipho_dijet_vars.dipho_SubPhoPt         , "dipho_SubPhoPt/D"        ); 
    tree->Branch( "dipho_SubPhoEta"         , &dipho_dijet_vars.dipho_SubPhoEta        , "dipho_SubPhoEta/D"       ); 
    tree->Branch( "dipho_SubPhoPhi"         , &dipho_dijet_vars.dipho_SubPhoPhi        , "dipho_SubPhoPhi/D"       ); 
    tree->Branch( "dipho_SubPhoE"           , &dipho_dijet_vars.dipho_SubPhoE          , "dipho_SubPhoE/D"         ); 
    tree->Branch( "dijet_LeadJPt"           , &dipho_dijet_vars.dijet_LeadJPt          , "dijet_LeadJPt/D"         ); 
    tree->Branch( "dijet_LeadJEta"          , &dipho_dijet_vars.dijet_LeadJEta         , "dijet_LeadJEta/D"        ); 
    tree->Branch( "dijet_LeadJPhi"          , &dipho_dijet_vars.dijet_LeadJPhi         , "dijet_LeadJPhi/D"        ); 
    tree->Branch( "dijet_LeadJE"            , &dipho_dijet_vars.dijet_LeadJE           , "dijet_LeadJE/D"          ); 
    tree->Branch( "dijet_SubJPt"            , &dipho_dijet_vars.dijet_SubJPt           , "dijet_SubJPt/D"          ); 
    tree->Branch( "dijet_SubJEta"           , &dipho_dijet_vars.dijet_SubJEta          , "dijet_SubJEta/D"         ); 
    tree->Branch( "dijet_SubJPhi"           , &dipho_dijet_vars.dijet_SubJPhi          , "dijet_SubJPhi/D"         ); 
    tree->Branch( "dijet_SubJE"             , &dipho_dijet_vars.dijet_SubJE            , "dijet_SubJE/D"           ); 
    tree->Branch( "dijet_abs_dEta"          , &dipho_dijet_vars.dijet_abs_dEta         , "dijet_abs_dEta/D"        ); 
    tree->Branch( "dijet_Mjj"               , &dipho_dijet_vars.dijet_Mjj              , "dijet_Mjj/D"             ); 
    tree->Branch( "dijet_centrality_gg"     , &dipho_dijet_vars.dijet_centrality_gg    , "dijet_centrality_gg/D"   );
    tree->Branch( "dijet_dphi"              , &dipho_dijet_vars.dijet_dphi             , "dijet_dphi/D"            );
    tree->Branch( "dijet_minDRJetPho"       , &dipho_dijet_vars.dijet_minDRJetPho      , "dijet_minDRJetPho/D"     );
    tree->Branch( "dijet_dipho_dphi_trunc"  , &dipho_dijet_vars.dijet_dipho_dphi_trunc , "dijet_dipho_dphi_trunc/D");
    tree->Branch( "event_pt"      , &dipho_dijet_vars.event_pt     , "event_pt/D"     );
    tree->Branch( "event_eta"     , &dipho_dijet_vars.event_eta    , "event_eta/D"    );
    tree->Branch( "event_phi"     , &dipho_dijet_vars.event_phi    , "event_phi/D"    );
    tree->Branch( "event_mass"    , &dipho_dijet_vars.event_mass   , "event_mass/D"   );
    tree->Branch( "event_energy"  , &dipho_dijet_vars.event_energy , "event_energy/D" );

    tree->Branch( "PF_eta"            , &PF_eta            );
    tree->Branch( "PF_phi"            , &PF_phi            );
    tree->Branch( "PF_pt"             , &PF_pt             );
    tree->Branch( "PF_energy"         , &PF_energy         );
    tree->Branch( "PF_rel_pt"         , &PF_rel_pt         );
    tree->Branch( "PF_rel_energy"     , &PF_rel_energy     );
    tree->Branch( "PF_deltaEta"       , &PF_deltaEta       );
    tree->Branch( "PF_deltaPhi"       , &PF_deltaPhi       );
    tree->Branch( "PF_charge"         , &PF_charge         );
    tree->Branch( "PF_isEle"          , &PF_isEle          );
    tree->Branch( "PF_isMu"           , &PF_isMu           );
    tree->Branch( "PF_isChargedHad"   , &PF_isChargedHad   );
    tree->Branch( "PF_isGamma"        , &PF_isGamma        );
    tree->Branch( "PF_isNeutralHad"   , &PF_isNeutralHad   );
    tree->Branch( "PF_tanhd0"         , &PF_tanhd0         );
    tree->Branch( "PF_tanhdz"         , &PF_tanhdz         );
    tree->Branch( "PF_sigmad0"        , &PF_sigmad0        );
    tree->Branch( "PF_sigmadz"        , &PF_sigmadz        );

    tree->Branch( "Jet1_PF_eta"            , &Jet1_PF_eta            );
    tree->Branch( "Jet1_PF_phi"            , &Jet1_PF_phi            );
    tree->Branch( "Jet1_PF_pt"             , &Jet1_PF_pt             );
    tree->Branch( "Jet1_PF_energy"         , &Jet1_PF_energy         );
    tree->Branch( "Jet1_PF_rel_pt"         , &Jet1_PF_rel_pt         );
    tree->Branch( "Jet1_PF_rel_energy"     , &Jet1_PF_rel_energy     );
    tree->Branch( "Jet1_PF_deltaEta"       , &Jet1_PF_deltaEta       );
    tree->Branch( "Jet1_PF_deltaPhi"       , &Jet1_PF_deltaPhi       );
    tree->Branch( "Jet1_PF_charge"         , &Jet1_PF_charge         );
    tree->Branch( "Jet1_PF_isEle"          , &Jet1_PF_isEle          );
    tree->Branch( "Jet1_PF_isMu"           , &Jet1_PF_isMu           );
    tree->Branch( "Jet1_PF_isChargedHad"   , &Jet1_PF_isChargedHad   );
    tree->Branch( "Jet1_PF_isGamma"        , &Jet1_PF_isGamma        );
    tree->Branch( "Jet1_PF_isNeutralHad"   , &Jet1_PF_isNeutralHad   );
    tree->Branch( "Jet1_PF_tanhd0"         , &Jet1_PF_tanhd0         );
    tree->Branch( "Jet1_PF_tanhdz"         , &Jet1_PF_tanhdz         );
    tree->Branch( "Jet1_PF_sigmad0"        , &Jet1_PF_sigmad0        );
    tree->Branch( "Jet1_PF_sigmadz"        , &Jet1_PF_sigmadz        );

    tree->Branch( "Jet2_PF_eta"            , &Jet2_PF_eta            );
    tree->Branch( "Jet2_PF_phi"            , &Jet2_PF_phi            );
    tree->Branch( "Jet2_PF_pt"             , &Jet2_PF_pt             );
    tree->Branch( "Jet2_PF_energy"         , &Jet2_PF_energy         );
    tree->Branch( "Jet2_PF_rel_pt"         , &Jet2_PF_rel_pt         );
    tree->Branch( "Jet2_PF_rel_energy"     , &Jet2_PF_rel_energy     );
    tree->Branch( "Jet2_PF_deltaEta"       , &Jet2_PF_deltaEta       );
    tree->Branch( "Jet2_PF_deltaPhi"       , &Jet2_PF_deltaPhi       );
    tree->Branch( "Jet2_PF_charge"         , &Jet2_PF_charge         );
    tree->Branch( "Jet2_PF_isEle"          , &Jet2_PF_isEle          );
    tree->Branch( "Jet2_PF_isMu"           , &Jet2_PF_isMu           );
    tree->Branch( "Jet2_PF_isChargedHad"   , &Jet2_PF_isChargedHad   );
    tree->Branch( "Jet2_PF_isGamma"        , &Jet2_PF_isGamma        );
    tree->Branch( "Jet2_PF_isNeutralHad"   , &Jet2_PF_isNeutralHad   );
    tree->Branch( "Jet2_PF_tanhd0"         , &Jet2_PF_tanhd0         );
    tree->Branch( "Jet2_PF_tanhdz"         , &Jet2_PF_tanhdz         );
    tree->Branch( "Jet2_PF_sigmad0"        , &Jet2_PF_sigmad0        );
    tree->Branch( "Jet2_PF_sigmadz"        , &Jet2_PF_sigmadz        );


    // Loop over all events
    for(int entry = 0; entry < treeReader->GetEntries(); ++entry) {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        //---------------------------- Particle Level --------------------------
    	std::vector<InputFeature_PF> PF_cands;

        PF_eta          .clear(); 
        PF_phi          .clear(); 
        PF_pt           .clear(); 
        PF_energy       .clear(); 
        PF_rel_pt       .clear(); 
        PF_rel_energy   .clear(); 
        PF_deltaEta     .clear(); 
        PF_deltaPhi     .clear(); 
        PF_charge       .clear(); 
        PF_isEle        .clear(); 
        PF_isMu         .clear(); 
        PF_isChargedHad .clear(); 
        PF_isGamma      .clear(); 
        PF_isNeutralHad .clear(); 
        PF_tanhd0       .clear(); 
        PF_tanhdz       .clear(); 
        PF_sigmad0      .clear(); 
        PF_sigmadz      .clear(); 

        TLorentzVector event_p4(0.,0.,0.,0.);

        // Electron, Muon, Charged hadron
        for ( int i = 0; i < branchEFlowTrack->GetEntriesFast(); ++i ){
        
            InputFeature_PF PF_cand;

            auto track = (Track*) branchEFlowTrack->At(i);
            auto particle = (GenParticle*) track->Particle.GetObject();

            TLorentzVector track_p4;
            track_p4.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, particle->Mass);

            PF_cand.PF_pt               = track_p4.Pt();
            PF_cand.PF_eta              = track_p4.Eta();
            PF_cand.PF_phi              = track_p4.Phi();
            PF_cand.PF_energy           = track_p4.Energy();
            PF_cand.PF_charge           = track->Charge;
            PF_cand.PF_isEle            = 0.;
            PF_cand.PF_isMu             = 0.;
            PF_cand.PF_isChargedHad     = 0.;
            PF_cand.PF_isGamma          = 0.;
            PF_cand.PF_isNeutralHad     = 0.;
            PF_cand.PF_tanhd0           = TMath::TanH(track->D0);
            PF_cand.PF_tanhdz           = TMath::TanH(track->DZ);
            PF_cand.PF_sigmad0          = track->ErrorD0;
            PF_cand.PF_sigmadz          = track->ErrorDZ;

            if (abs(particle->PID) == 11) PF_cand.PF_isEle = 1.;
            else if (abs(particle->PID) == 13) PF_cand.PF_isMu = 1.;
            else PF_cand.PF_isChargedHad = 1.;

            PF_cands.emplace_back(PF_cand);
            event_p4 += track_p4;
        }

        // Photon
        for ( int i = 0; i < branchEFlowPhoton->GetEntriesFast(); ++i ){

            InputFeature_PF PF_cand;

            auto tower = (Tower*) branchEFlowPhoton->At(i);
            auto pho_p4 = UseEtaPhiEMass(tower->Eta, tower->Phi, tower->E, 0.);

            PF_cand.PF_pt               = pho_p4.Pt();
            PF_cand.PF_eta              = pho_p4.Eta();
            PF_cand.PF_phi              = pho_p4.Phi();
            PF_cand.PF_energy           = pho_p4.Energy();
            PF_cand.PF_charge           = 0.;
            PF_cand.PF_isEle            = 0.;
            PF_cand.PF_isMu             = 0.;
            PF_cand.PF_isChargedHad     = 0.;
            PF_cand.PF_isGamma          = 1.;
            PF_cand.PF_isNeutralHad     = 0.;
            PF_cand.PF_tanhd0           = 0.;
            PF_cand.PF_tanhdz           = 0.;
            PF_cand.PF_sigmad0          = 0.;
            PF_cand.PF_sigmadz          = 0.;

            PF_cands.emplace_back(PF_cand);
            event_p4 += pho_p4;
        }

        // Neutron hadron
        for ( int i = 0; i < branchEFlowNeutralHadron->GetEntriesFast(); ++i ){

            InputFeature_PF PF_cand;
                                                                                         
            auto tower = (Tower*) branchEFlowNeutralHadron->At(i);
            auto particle = (GenParticle*) tower->Particles.At(0);

            auto neutralH_p4 = UseEtaPhiEMass(tower->Eta, tower->Phi, tower->E, particle->Mass);

            PF_cand.PF_pt               = neutralH_p4.Pt();
            PF_cand.PF_eta              = neutralH_p4.Eta();
            PF_cand.PF_phi              = neutralH_p4.Phi();
            PF_cand.PF_energy           = neutralH_p4.Energy();
            PF_cand.PF_charge           = 0.;
            PF_cand.PF_isEle            = 0.;
            PF_cand.PF_isMu             = 0.;
            PF_cand.PF_isChargedHad     = 0.;
            PF_cand.PF_isGamma          = 0.;
            PF_cand.PF_isNeutralHad     = 1.;
            PF_cand.PF_tanhd0           = 0.;
            PF_cand.PF_tanhdz           = 0.;
            PF_cand.PF_sigmad0          = 0.;
            PF_cand.PF_sigmadz          = 0.;

            PF_cands.emplace_back(PF_cand);
            event_p4 += neutralH_p4;
        }

        h1_monitor_npf->Fill( PF_cands.size() );

        dipho_dijet_vars.event_pt     = event_p4.Pt()     ;
        dipho_dijet_vars.event_eta    = event_p4.Eta()    ;
        dipho_dijet_vars.event_phi    = event_p4.Phi()    ;
        dipho_dijet_vars.event_mass   = event_p4.M()      ;
        dipho_dijet_vars.event_energy = event_p4.Energy() ;

        std::sort(PF_cands.begin(), PF_cands.end(), pt_sortor());

        // MET is put in the first element
        PF_eta            .emplace_back( 0.                                                     );
        PF_phi            .emplace_back( ((MissingET*) branchMET->At(0))->Phi                   );
        PF_pt             .emplace_back( ((MissingET*) branchMET->At(0))->MET                   );
        PF_energy         .emplace_back( 0.                                                     );
        PF_rel_pt         .emplace_back( ((MissingET*) branchMET->At(0))->MET/event_p4.Energy() );
        PF_rel_energy     .emplace_back( 0.                                                     );
        PF_deltaEta       .emplace_back( 0.                                                     );
        PF_deltaPhi       .emplace_back( ((MissingET*) branchMET->At(0))->Phi - event_p4.Phi()  );
        PF_charge         .emplace_back( 0.                                                     ); 
        PF_isEle          .emplace_back( 0.                                                     ); 
        PF_isMu           .emplace_back( 0.                                                     ); 
        PF_isChargedHad   .emplace_back( 0.                                                     ); 
        PF_isGamma        .emplace_back( 0.                                                     ); 
        PF_isNeutralHad   .emplace_back( 0.                                                     ); 
        PF_tanhd0         .emplace_back( 0.                                                     ); 
        PF_tanhdz         .emplace_back( 0.                                                     ); 
        PF_sigmad0        .emplace_back( 0.                                                     ); 
        PF_sigmadz        .emplace_back( 0.                                                     ); 


        for (const auto it : PF_cands) {
            PF_eta            .emplace_back( it.PF_eta                        ); 
            PF_phi            .emplace_back( it.PF_phi                        );
            PF_pt             .emplace_back( it.PF_pt                         );
            PF_energy         .emplace_back( it.PF_energy                     );
            PF_rel_pt         .emplace_back( it.PF_pt/event_p4.Energy()       );
            PF_rel_energy     .emplace_back( it.PF_energy/event_p4.Energy()   );
            PF_deltaEta       .emplace_back( it.PF_eta - event_p4.Eta()       );
            PF_deltaPhi       .emplace_back( it.PF_phi - event_p4.Phi()       );
            PF_charge         .emplace_back( it.PF_charge                     );
            PF_isEle          .emplace_back( it.PF_isEle                      );
            PF_isMu           .emplace_back( it.PF_isMu                       );
            PF_isChargedHad   .emplace_back( it.PF_isChargedHad               );
            PF_isGamma        .emplace_back( it.PF_isGamma                    );
            PF_isNeutralHad   .emplace_back( it.PF_isNeutralHad               );
            PF_tanhd0         .emplace_back( it.PF_tanhd0                     ); 
            PF_tanhdz         .emplace_back( it.PF_tanhdz                     ); 
            PF_sigmad0        .emplace_back( it.PF_sigmad0                    ); 
            PF_sigmadz        .emplace_back( it.PF_sigmadz                    ); 
        }

        //---------------------------- Object Level --------------------------
        vector<TLorentzVector> photons_p4;

        // Loop over all photons in events
        for (int i_pho = 0; i_pho < branchPhoton->GetEntriesFast(); ++i_pho) {

            auto photon = (Photon*) branchPhoton->At(i_pho);

            //if( photon->PT < 25. ) continue;
            if( fabs(photon->Eta) > 2.5 ) continue;

            TLorentzVector pho_p4;
            pho_p4.SetPtEtaPhiE(photon->PT, photon->Eta, photon->Phi, photon->E);

            photons_p4.emplace_back(pho_p4);
        }

        if (photons_p4.size() < 2) continue;

        Jet1_PF_eta           .clear();
        Jet1_PF_phi           .clear();
        Jet1_PF_pt            .clear();
        Jet1_PF_energy        .clear();
        Jet1_PF_rel_pt        .clear();
        Jet1_PF_rel_energy    .clear();
        Jet1_PF_deltaEta      .clear();
        Jet1_PF_deltaPhi      .clear();
        Jet1_PF_charge        .clear();
        Jet1_PF_isEle         .clear();
        Jet1_PF_isMu          .clear();
        Jet1_PF_isChargedHad  .clear();
        Jet1_PF_isGamma       .clear();
        Jet1_PF_isNeutralHad  .clear();
        Jet1_PF_tanhd0        .clear();
        Jet1_PF_tanhdz        .clear();
        Jet1_PF_sigmad0       .clear();
        Jet1_PF_sigmadz       .clear();
 
        Jet2_PF_eta           .clear();
        Jet2_PF_phi           .clear();
        Jet2_PF_pt            .clear();
        Jet2_PF_energy        .clear();
        Jet2_PF_rel_pt        .clear();
        Jet2_PF_rel_energy    .clear();
        Jet2_PF_deltaEta      .clear();
        Jet2_PF_deltaPhi      .clear();
        Jet2_PF_charge        .clear();
        Jet2_PF_isEle         .clear();
        Jet2_PF_isMu          .clear();
        Jet2_PF_isChargedHad  .clear();
        Jet2_PF_isGamma       .clear();
        Jet2_PF_isNeutralHad  .clear();
        Jet2_PF_tanhd0        .clear();
        Jet2_PF_tanhdz        .clear();
        Jet2_PF_sigmad0       .clear();
        Jet2_PF_sigmadz       .clear();

        int njet = 0;
        vector<TLorentzVector> jets_p4;

        // Loop over all jets in events
        for (int i_jet = 0; i_jet < branchJet->GetEntriesFast(); ++i_jet) {

            auto jet = (Jet*) branchJet->At(i_jet);

        //    if ( jet->PT < 30. ) continue;
            if ( fabs(jet->Eta) > 4.7 ) continue;
        //    if ( jet->NCharged < 1 ) continue;
            if ( jet->Constituents.GetEntriesFast() < 2 ) continue;

            TLorentzVector jet_p4;
            jet_p4.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);

            if (jet_p4.DeltaR(photons_p4[0]) < 0.4 || jet_p4.DeltaR(photons_p4[1]) < 0.4) continue;

            std::vector<InputFeature_PF> Jet_PF_cands;

            int nCH = 0;
            // Loop over all jet's constituents
            for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j) {
                auto object = jet->Constituents.At(j);

                // Check if the constituent is accessible
                if(object == 0) continue;

                InputFeature_PF Jet_PF_cand;

                Jet_PF_cand.PF_charge        = 0.;
                Jet_PF_cand.PF_isEle         = 0.;
                Jet_PF_cand.PF_isMu          = 0.;
                Jet_PF_cand.PF_isChargedHad  = 0.;
                Jet_PF_cand.PF_isGamma       = 0.;
                Jet_PF_cand.PF_isNeutralHad  = 0.;
                Jet_PF_cand.PF_tanhd0        = 0.; 
                Jet_PF_cand.PF_tanhdz        = 0.; 
                Jet_PF_cand.PF_sigmad0       = 0.; 
                Jet_PF_cand.PF_sigmadz       = 0.; 

                // Charged Particle
                if(object->IsA() == Track::Class()){
                    auto track = (Track*) object;
                    auto particle = (GenParticle*) track->Particle.GetObject();

                    TLorentzVector track_p4;
                    track_p4.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, particle->Mass);

                    Jet_PF_cand.PF_eta = track_p4.Eta();
                    Jet_PF_cand.PF_phi = track_p4.Phi();

                    Jet_PF_cand.PF_pt            = track_p4.Pt()                  ;
                    Jet_PF_cand.PF_energy        = track_p4.E()                   ;
                    Jet_PF_cand.PF_rel_pt        = track_p4.Pt() / jet_p4.Pt()    ;
                    Jet_PF_cand.PF_rel_energy    = track_p4.E() / jet_p4.E()      ;
                    Jet_PF_cand.PF_deltaEta      = track_p4.Eta() - jet_p4.Eta()  ;
                    Jet_PF_cand.PF_deltaPhi      = track_p4.DeltaPhi(jet_p4)      ;
                    Jet_PF_cand.PF_charge        = track->Charge                  ;
                    Jet_PF_cand.PF_tanhd0        = TMath::TanH(track->D0)         ;
                    Jet_PF_cand.PF_tanhdz        = TMath::TanH(track->DZ)         ;
                    Jet_PF_cand.PF_sigmad0       = track->ErrorD0                 ;
                    Jet_PF_cand.PF_sigmadz       = track->ErrorDZ                 ;

                    if (abs(particle->PID) == 11) {
                        // Is Electron
                        Jet_PF_cand.PF_isEle = 1.;
                    } else if (abs(particle->PID) == 13) {
                        // Is Muon
                        Jet_PF_cand.PF_isMu = 1.;
                    } else {
                        // Is Charged Hadron
                        nCH++;
                        Jet_PF_cand.PF_isChargedHad = 1.;
                    }

                }else if(object->IsA() == Tower::Class()){
                    auto tower = (Tower*) object;
                    auto particle = (GenParticle*) tower->Particles.At(0);

                    auto tower_p4 = UseEtaPhiEMass(tower->Eta, tower->Phi, tower->E, particle->Mass);

                    Jet_PF_cand.PF_eta = tower_p4.Eta();
                    Jet_PF_cand.PF_phi = tower_p4.Phi();

                    Jet_PF_cand.PF_pt            = tower_p4.Pt()                 ;
                    Jet_PF_cand.PF_energy        = tower_p4.E()                  ;
                    Jet_PF_cand.PF_rel_pt        = tower_p4.Pt() / jet_p4.Pt()   ;
                    Jet_PF_cand.PF_rel_energy    = tower_p4.E() / jet_p4.E()     ;
                    Jet_PF_cand.PF_deltaEta      = tower_p4.Eta() - jet_p4.Eta() ;
                    Jet_PF_cand.PF_deltaPhi      = tower_p4.DeltaPhi(jet_p4)     ;

                    if (particle->PID == 22) {
                        // Is Photon
                        Jet_PF_cand.PF_isGamma = 1.;
                    } else {
                        // Is Neutral Hadron
                        Jet_PF_cand.PF_isNeutralHad  = 1.;
                    }
                }

                Jet_PF_cands.emplace_back(Jet_PF_cand);

            }//Constituents++

            //if (nCH == 0) continue;

            jets_p4.emplace_back(jet_p4);

            njet++;
            if(njet == 1) {

                std::sort(Jet_PF_cands.begin(), Jet_PF_cands.end(), pt_sortor());

                for (const auto it : Jet_PF_cands) {

                    Jet1_PF_eta            .emplace_back( it.PF_eta           ); 
                    Jet1_PF_phi            .emplace_back( it.PF_phi           );
                    Jet1_PF_pt             .emplace_back( it.PF_pt            );
                    Jet1_PF_energy         .emplace_back( it.PF_energy        );
                    Jet1_PF_rel_pt         .emplace_back( it.PF_rel_pt        );
                    Jet1_PF_rel_energy     .emplace_back( it.PF_rel_energy    );
                    Jet1_PF_deltaEta       .emplace_back( it.PF_deltaEta      );
                    Jet1_PF_deltaPhi       .emplace_back( it.PF_deltaPhi      );
                    Jet1_PF_charge         .emplace_back( it.PF_charge        );
                    Jet1_PF_isEle          .emplace_back( it.PF_isEle         );
                    Jet1_PF_isMu           .emplace_back( it.PF_isMu          );
                    Jet1_PF_isChargedHad   .emplace_back( it.PF_isChargedHad  );
                    Jet1_PF_isGamma        .emplace_back( it.PF_isGamma       );
                    Jet1_PF_isNeutralHad   .emplace_back( it.PF_isNeutralHad  );
                    Jet1_PF_tanhd0         .emplace_back( it.PF_tanhd0        ); 
                    Jet1_PF_tanhdz         .emplace_back( it.PF_tanhdz        ); 
                    Jet1_PF_sigmad0        .emplace_back( it.PF_sigmad0       ); 
                    Jet1_PF_sigmadz        .emplace_back( it.PF_sigmadz       ); 

                }
            } else if (njet == 2) {

                std::sort(Jet_PF_cands.begin(), Jet_PF_cands.end(), pt_sortor());
                                                                                     
                for (const auto it : Jet_PF_cands) {
                                                                                     
                    Jet2_PF_eta            .emplace_back( it.PF_eta           ); 
                    Jet2_PF_phi            .emplace_back( it.PF_phi           );
                    Jet2_PF_pt             .emplace_back( it.PF_pt            );
                    Jet2_PF_energy         .emplace_back( it.PF_energy        );
                    Jet2_PF_rel_pt         .emplace_back( it.PF_rel_pt        );
                    Jet2_PF_rel_energy     .emplace_back( it.PF_rel_energy    );
                    Jet2_PF_deltaEta       .emplace_back( it.PF_deltaEta      );
                    Jet2_PF_deltaPhi       .emplace_back( it.PF_deltaPhi      );
                    Jet2_PF_charge         .emplace_back( it.PF_charge        );
                    Jet2_PF_isEle          .emplace_back( it.PF_isEle         );
                    Jet2_PF_isMu           .emplace_back( it.PF_isMu          );
                    Jet2_PF_isChargedHad   .emplace_back( it.PF_isChargedHad  );
                    Jet2_PF_isGamma        .emplace_back( it.PF_isGamma       );
                    Jet2_PF_isNeutralHad   .emplace_back( it.PF_isNeutralHad  );
                    Jet2_PF_tanhd0         .emplace_back( it.PF_tanhd0        ); 
                    Jet2_PF_tanhdz         .emplace_back( it.PF_tanhdz        ); 
                    Jet2_PF_sigmad0        .emplace_back( it.PF_sigmad0       ); 
                    Jet2_PF_sigmadz        .emplace_back( it.PF_sigmadz       ); 
                                                                                     
                }
            }

        }//Jet++

        h1_monitor_njet->Fill(njet);

        if (njet < 2) continue;

        TLorentzVector pho1_p4, pho2_p4, jet1_p4, jet2_p4, dipho_p4, dijet_p4;
        pho1_p4 = photons_p4[0];
        pho2_p4 = photons_p4[1];
        jet1_p4 = jets_p4[0];
        jet2_p4 = jets_p4[1];
        dipho_p4 = pho1_p4 + pho2_p4;
        dijet_p4 = jet1_p4 + jet2_p4;

        if (pho1_p4.Pt() < 30. || pho2_p4.Pt() < 20. ) continue;
        if (pho1_p4.Pt() / dipho_p4.M() < 1./4 || pho2_p4.Pt() / dipho_p4.M() < 1./5) continue;
        if (dipho_p4.M() < 115. || dipho_p4.M() > 135.) continue;
        if (jet1_p4.Pt() < 30.) continue;
        //if (dijet_p4.M() < 100.) continue;

        dipho_dijet_vars.dipho_mass                = dipho_p4.M();
        dipho_dijet_vars.dipho_PToM                = dipho_p4.Pt()/ dipho_p4.M();
        dipho_dijet_vars.leadPho_PToM              = pho1_p4.Pt() / dipho_p4.M();
        dipho_dijet_vars.sublPho_PToM              = pho2_p4.Pt() / dipho_p4.M();
        dipho_dijet_vars.dipho_LeadPhoPt           = pho1_p4.Pt();
        dipho_dijet_vars.dipho_LeadPhoEta          = pho1_p4.Eta();
        dipho_dijet_vars.dipho_LeadPhoPhi          = pho1_p4.Phi();
        dipho_dijet_vars.dipho_LeadPhoE            = pho1_p4.Energy();
        dipho_dijet_vars.dipho_SubPhoPt            = pho2_p4.Pt();
        dipho_dijet_vars.dipho_SubPhoEta           = pho2_p4.Eta();
        dipho_dijet_vars.dipho_SubPhoPhi           = pho2_p4.Phi();
        dipho_dijet_vars.dipho_SubPhoE             = pho2_p4.Energy();
        dipho_dijet_vars.dijet_LeadJPt             = jet1_p4.Pt();
        dipho_dijet_vars.dijet_LeadJEta            = jet1_p4.Eta();
        dipho_dijet_vars.dijet_LeadJPhi            = jet1_p4.Phi();
        dipho_dijet_vars.dijet_LeadJE              = jet1_p4.Energy();
        dipho_dijet_vars.dijet_SubJPt              = jet2_p4.Pt();
        dipho_dijet_vars.dijet_SubJEta             = jet2_p4.Eta();
        dipho_dijet_vars.dijet_SubJPhi             = jet2_p4.Phi();
        dipho_dijet_vars.dijet_SubJE               = jet2_p4.Energy();
        dipho_dijet_vars.dijet_abs_dEta            = fabs(jet1_p4.Eta() - jet2_p4.Eta()); 
        dipho_dijet_vars.dijet_Mjj                 = dijet_p4.M();
        dipho_dijet_vars.dijet_centrality_gg = TMath::Exp(-4*TMath::Power( fabs( dipho_p4.Eta() - 0.5*(jet1_p4.Eta()+jet2_p4.Eta()) )  / dipho_dijet_vars.dijet_abs_dEta, 2 ) );
        dipho_dijet_vars.dijet_dphi                = fabs( jet1_p4.DeltaPhi(jet2_p4) ); 
        dipho_dijet_vars.dijet_minDRJetPho         = std::min( std::min( jet1_p4.DeltaR(pho1_p4), jet1_p4.DeltaR(pho2_p4) ), 
                                                               std::min( jet2_p4.DeltaR(pho1_p4), jet2_p4.DeltaR(pho2_p4) ) );
        dipho_dijet_vars.dijet_dipho_dphi_trunc    = std::min( fabs(dipho_p4.DeltaPhi(dijet_p4)), 2.9416 );

        tree->Fill();

    }//Event++

    file->Write();
    file->Close();


    delete treeReader;
    delete chain;
}

//------------------------------------------------------------------------------
