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
    double dijet_LeadJPt       ;
    double dijet_LeadJEta      ;
    double dijet_SubJPt        ;
    double dijet_SubJEta       ;
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
    double PF_pt;
    double PF_eta;
    double PF_phi;
    double PF_energy;
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

    int PF_JetIndex;
};

struct InputFeature_Object
{
    double Object_pt;
    double Object_eta;
    double Object_phi;
    double Object_energy;
    double Object_charge;
    double Object_isEle;
    double Object_isMu;
    double Object_isGamma;
    double Object_isGluon;
    double Object_isLight;
    double Object_isCharm;
    double Object_isBottom;
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

void  event_Extractor(const char *inputFile)
{
    gSystem->Load("libDelphes");

    InputFeature dipho_dijet_vars;

    // Object output
    std::vector<double> Object_pt       ; 
    std::vector<double> Object_eta      ; 
    std::vector<double> Object_phi      ; 
    std::vector<double> Object_energy   ; 
    std::vector<double> Object_charge   ; 
    std::vector<double> Object_isEle    ;
    std::vector<double> Object_isMu     ;
    std::vector<double> Object_isGamma  ;
    std::vector<double> Object_isGluon  ;
    std::vector<double> Object_isLight  ;
    std::vector<double> Object_isCharm  ;
    std::vector<double> Object_isBottom ;

    // Particle Flow output
    std::vector<double> PF_eta          ; 
    std::vector<double> PF_phi          ; 
    std::vector<double> PF_pt           ; 
    std::vector<double> PF_energy       ; 
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
    std::vector<int>    PF_JetIndex     ; 

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchPFCandicate = treeReader->UseBranch("ParticleFlowCandidate");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchMET = treeReader->UseBranch("MissingET");

    TFile* file = TFile::Open("output_events.root", "recreate");
    TTree* tree = new TTree("Events", "");

    TH1F*h1_monitor_njet = new TH1F("h1_monitor_njet", "", 10, 0, 10);
    TH1F*h1_monitor_npf  = new TH1F("h1_monitor_npf", "", 50, 0, 300);

    tree->Branch( "dipho_mass"              , &dipho_dijet_vars.dipho_mass             , "dipho_mass/D"            ); 
    tree->Branch( "dipho_PToM"              , &dipho_dijet_vars.dipho_PToM             , "dipho_PToM/D"            ); 
    tree->Branch( "leadPho_PToM"            , &dipho_dijet_vars.leadPho_PToM           , "leadPho_PToM/D"          ); 
    tree->Branch( "sublPho_PToM"            , &dipho_dijet_vars.sublPho_PToM           , "sublPho_PToM/D"          );
    tree->Branch( "dijet_LeadJPt"           , &dipho_dijet_vars.dijet_LeadJPt          , "dijet_LeadJPt/D"         ); 
    tree->Branch( "dijet_LeadJEta"          , &dipho_dijet_vars.dijet_LeadJEta         , "dijet_LeadJEta/D"        ); 
    tree->Branch( "dijet_SubJPt"            , &dipho_dijet_vars.dijet_SubJPt           , "dijet_SubJPt/D"          ); 
    tree->Branch( "dijet_SubJEta"           , &dipho_dijet_vars.dijet_SubJEta          , "dijet_SubJEta/D"         ); 
    tree->Branch( "dijet_abs_dEta"          , &dipho_dijet_vars.dijet_abs_dEta         , "dijet_abs_dEta/D"        ); 
    tree->Branch( "dijet_Mjj"               , &dipho_dijet_vars.dijet_Mjj              , "dijet_Mjj/D"             ); 
    tree->Branch( "dijet_centrality_gg"     , &dipho_dijet_vars.dijet_centrality_gg    , "dijet_centrality_gg/D"   );
    tree->Branch( "dijet_dphi"              , &dipho_dijet_vars.dijet_dphi             , "dijet_dphi/D"            );
    tree->Branch( "dijet_minDRJetPho"       , &dipho_dijet_vars.dijet_minDRJetPho      , "dijet_minDRJetPho/D"     );
    tree->Branch( "dijet_dipho_dphi_trunc"  , &dipho_dijet_vars.dijet_dipho_dphi_trunc , "dijet_dipho_dphi_trunc/D");
    tree->Branch( "event_pt"                , &dipho_dijet_vars.event_pt               , "event_pt/D"              );
    tree->Branch( "event_eta"               , &dipho_dijet_vars.event_eta              , "event_eta/D"             );
    tree->Branch( "event_phi"               , &dipho_dijet_vars.event_phi              , "event_phi/D"             );
    tree->Branch( "event_mass"              , &dipho_dijet_vars.event_mass             , "event_mass/D"            );
    tree->Branch( "event_energy"            , &dipho_dijet_vars.event_energy           , "event_energy/D"          );

    tree->Branch( "Object_pt"         , &Object_pt         ); 
    tree->Branch( "Object_eta"        , &Object_eta        ); 
    tree->Branch( "Object_phi"        , &Object_phi        ); 
    tree->Branch( "Object_energy"     , &Object_energy     ); 
    tree->Branch( "Object_charge"     , &Object_charge     ); 
    tree->Branch( "Object_isEle"      , &Object_isEle      ); 
    tree->Branch( "Object_isMu"       , &Object_isMu       ); 
    tree->Branch( "Object_isGamma"    , &Object_isGamma    ); 
    tree->Branch( "Object_isGluon"    , &Object_isGluon    ); 
    tree->Branch( "Object_isLight"    , &Object_isLight    ); 
    tree->Branch( "Object_isCharm"    , &Object_isCharm    ); 
    tree->Branch( "Object_isBottom"   , &Object_isBottom   ); 

    tree->Branch( "PF_eta"            , &PF_eta            );
    tree->Branch( "PF_phi"            , &PF_phi            );
    tree->Branch( "PF_pt"             , &PF_pt             );
    tree->Branch( "PF_energy"         , &PF_energy         );
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
    tree->Branch( "PF_JetIndex"       , &PF_JetIndex       );

    // Loop over all events
    for(int entry = 0; entry < treeReader->GetEntries(); ++entry) {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        //---------------------------- Particle Level --------------------------
        std::map<ParticleFlowCandidate*, InputFeature_PF> PF_cands_mapping;

        TLorentzVector event_p4(0.,0.,0.,0.);
        for ( int i = 0; i < branchPFCandicate->GetEntriesFast(); ++i ){

            InputFeature_PF PFCandFeature;

            auto PFCand = (ParticleFlowCandidate*) branchPFCandicate->At(i);

            PFCandFeature.PF_pt               = PFCand->PT;
            PFCandFeature.PF_eta              = PFCand->Eta;
            PFCandFeature.PF_phi              = PFCand->Phi;
            PFCandFeature.PF_energy           = PFCand->E;
            PFCandFeature.PF_charge           = PFCand->Charge;
            PFCandFeature.PF_isEle            = (double)( abs(PFCand->PID) == 11 );
            PFCandFeature.PF_isMu             = (double)( abs(PFCand->PID) == 13 );
            PFCandFeature.PF_isChargedHad     = (double)( PFCand->Charge != 0 && abs(PFCand->PID) != 11 && abs(PFCand->PID) != 13 );
            PFCandFeature.PF_isGamma          = (double)( PFCand->PID == 22 );
            PFCandFeature.PF_isNeutralHad     = (double)( PFCand->PID == 0  );
            PFCandFeature.PF_tanhd0           = TMath::TanH(PFCand->D0);
            PFCandFeature.PF_tanhdz           = TMath::TanH(PFCand->DZ);
            PFCandFeature.PF_sigmad0          = PFCand->ErrorD0;
            PFCandFeature.PF_sigmadz          = PFCand->ErrorDZ;

            // Assign PF_JetIndex as -1 first
            PFCandFeature.PF_JetIndex         = -1;

            PF_cands_mapping[PFCand] = PFCandFeature;

            TLorentzVector PFCand_p4;
            PFCand_p4.SetPtEtaPhiE(PFCand->PT, PFCand->Eta, PFCand->Phi, PFCand->E);
            event_p4 += PFCand_p4;
        }

        h1_monitor_npf->Fill( PF_cands_mapping.size() );

        //---------------------------- Object Level --------------------------
        std::vector<InputFeature_Object> Objects;

        // Loop over all photons in events
        std::vector<TLorentzVector> photons_p4;

        for (int i_pho = 0; i_pho < branchPhoton->GetEntriesFast(); ++i_pho) {

            auto photon = (Photon*) branchPhoton->At(i_pho);

            if( photon->PT < 20. ) continue;
            if( fabs(photon->Eta) > 2.5 ) continue;

            TLorentzVector pho_p4;
            pho_p4.SetPtEtaPhiE(photon->PT, photon->Eta, photon->Phi, photon->E);

            photons_p4.emplace_back(pho_p4);

            InputFeature_Object obj_photon_feature;
            obj_photon_feature.Object_pt       = photon->PT;
            obj_photon_feature.Object_eta      = photon->Eta;
            obj_photon_feature.Object_phi      = photon->Phi;
            obj_photon_feature.Object_energy   = photon->E;
            obj_photon_feature.Object_charge   = 0;
            obj_photon_feature.Object_isEle    = 0;
            obj_photon_feature.Object_isMu     = 0;
            obj_photon_feature.Object_isGamma  = 1;
            obj_photon_feature.Object_isGluon  = 0;
            obj_photon_feature.Object_isLight  = 0;
            obj_photon_feature.Object_isCharm  = 0;
            obj_photon_feature.Object_isBottom = 0;

            Objects.emplace_back(obj_photon_feature);

        }

        if (photons_p4.size() < 2) continue;

        // Loop over all jets in events
        int jet_index = 0;
        vector<TLorentzVector> jets_p4;

        for (int i_jet = 0; i_jet < branchJet->GetEntriesFast(); ++i_jet) {

            auto jet = (Jet*) branchJet->At(i_jet);

            if ( jet->PT < 30. ) continue;
            if ( fabs(jet->Eta) > 4.7 ) continue;
            if ( jet->Constituents.GetEntriesFast() < 2 ) continue;

            TLorentzVector jet_p4;
            jet_p4.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);

            if (jet_p4.DeltaR(photons_p4[0]) < 0.4 || jet_p4.DeltaR(photons_p4[1]) < 0.4) continue;
            jets_p4.emplace_back(jet_p4);

            // Link all jet's constituents to the particles in the events
            for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j) {
                auto object = jet->Constituents.At(j);

                if(object == 0) continue;

                if ( object->IsA() == ParticleFlowCandidate::Class() ) {
                    auto PFcand = (ParticleFlowCandidate*) object;
                    PF_cands_mapping[PFcand].PF_JetIndex = jet_index; 
                }
            }

            // Jet Matching to Gen Parton
            // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/mcMatchLayer0/jetMatch_cfi.py#L6
            int pdgid_tmp = -1;                                                                       
            double min_deltaR = 9999.;
            for (int it_genP = 0; it_genP < branchParticle->GetEntriesFast(); ++it_genP) {
                auto genP = (GenParticle*)branchParticle->At(it_genP);

                if (abs(genP->PID) > 5 && abs(genP->PID) != 21 || abs(genP->PID) == 0) continue;
                if (genP->Status != 23) continue;

                TLorentzVector genP_p4;
                genP_p4.SetPtEtaPhiM(genP->PT, genP->Eta, genP->Phi, genP->Mass);

                double deltaR = jet_p4.DeltaR(genP_p4);
                if (deltaR < 0.4 && fabs(genP_p4.Pt() - jet_p4.Pt())/genP_p4.Pt() < 3.) {
                    if (min_deltaR > deltaR) {
                        min_deltaR = deltaR;
                        pdgid_tmp = abs(genP->PID);
                    }
                }
            }

            InputFeature_Object obj_jet_feature;
            obj_jet_feature.Object_pt       = jet_p4.Pt();
            obj_jet_feature.Object_eta      = jet_p4.Eta();
            obj_jet_feature.Object_phi      = jet_p4.Phi();
            obj_jet_feature.Object_energy   = jet_p4.Energy();
            obj_jet_feature.Object_charge   = 0.;
            obj_jet_feature.Object_isEle    = 0.;
            obj_jet_feature.Object_isMu     = 0.;
            obj_jet_feature.Object_isGamma  = 0.;
            obj_jet_feature.Object_isGluon  = (double)(pdgid_tmp == 21);
            obj_jet_feature.Object_isLight  = (double)(pdgid_tmp < 4 && pdgid_tmp > 0);
            obj_jet_feature.Object_isCharm  = (double)(pdgid_tmp == 4);
            obj_jet_feature.Object_isBottom = (double)(pdgid_tmp == 5);

            Objects.emplace_back(obj_jet_feature);

            jet_index++;

        }//Jet++

        h1_monitor_njet->Fill(jet_index);

        if (jet_index < 2) continue;

        TLorentzVector pho1_p4, pho2_p4, jet1_p4, jet2_p4, dipho_p4, dijet_p4;
        pho1_p4 = photons_p4[0];
        pho2_p4 = photons_p4[1];
        jet1_p4 = jets_p4[0];
        jet2_p4 = jets_p4[1];
        dipho_p4 = pho1_p4 + pho2_p4;
        dijet_p4 = jet1_p4 + jet2_p4;

        if (pho1_p4.Pt() < 30.) continue;
        if (pho1_p4.Pt() / dipho_p4.M() < 1./4 || pho2_p4.Pt() / dipho_p4.M() < 1./5) continue;
        if (dipho_p4.M() < 110. || dipho_p4.M() > 140.) continue;

        dipho_dijet_vars.dipho_mass                = dipho_p4.M();
        dipho_dijet_vars.dipho_PToM                = dipho_p4.Pt()/ dipho_p4.M();
        dipho_dijet_vars.leadPho_PToM              = pho1_p4.Pt() / dipho_p4.M();
        dipho_dijet_vars.sublPho_PToM              = pho2_p4.Pt() / dipho_p4.M();
        dipho_dijet_vars.dijet_LeadJPt             = jet1_p4.Pt();
        dipho_dijet_vars.dijet_LeadJEta            = jet1_p4.Eta();
        dipho_dijet_vars.dijet_SubJPt              = jet2_p4.Pt();
        dipho_dijet_vars.dijet_SubJEta             = jet2_p4.Eta();
        dipho_dijet_vars.dijet_abs_dEta            = fabs(jet1_p4.Eta() - jet2_p4.Eta()); 
        dipho_dijet_vars.dijet_Mjj                 = dijet_p4.M();
        dipho_dijet_vars.dijet_centrality_gg       = TMath::Exp(-4*TMath::Power( fabs( dipho_p4.Eta() - 0.5*(jet1_p4.Eta()+jet2_p4.Eta()) )  / dipho_dijet_vars.dijet_abs_dEta, 2 ) );
        dipho_dijet_vars.dijet_dphi                = fabs( jet1_p4.DeltaPhi(jet2_p4) ); 
        dipho_dijet_vars.dijet_minDRJetPho         = std::min( std::min( jet1_p4.DeltaR(pho1_p4), jet1_p4.DeltaR(pho2_p4) ), 
                                                               std::min( jet2_p4.DeltaR(pho1_p4), jet2_p4.DeltaR(pho2_p4) ) );
        dipho_dijet_vars.dijet_dipho_dphi_trunc    = std::min( fabs(dipho_p4.DeltaPhi(dijet_p4)), 2.9416 );

        dipho_dijet_vars.event_pt     = event_p4.Pt()     ;
        dipho_dijet_vars.event_eta    = event_p4.Eta()    ;
        dipho_dijet_vars.event_phi    = event_p4.Phi()    ;
        dipho_dijet_vars.event_mass   = event_p4.M()      ;
        dipho_dijet_vars.event_energy = event_p4.Energy() ;

    	std::vector<InputFeature_PF> PF_cands;
        for (const auto& it : PF_cands_mapping) PF_cands.emplace_back( it.second );
        std::sort(PF_cands.begin(), PF_cands.end(), pt_sortor());


        Object_pt       .clear(); 
        Object_eta      .clear(); 
        Object_phi      .clear(); 
        Object_energy   .clear(); 
        Object_charge   .clear(); 
        Object_isEle    .clear(); 
        Object_isMu     .clear(); 
        Object_isGamma  .clear(); 
        Object_isGluon  .clear(); 
        Object_isLight  .clear(); 
        Object_isCharm  .clear(); 
        Object_isBottom .clear(); 

        PF_eta          .clear(); 
        PF_phi          .clear(); 
        PF_pt           .clear(); 
        PF_energy       .clear(); 
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
        PF_JetIndex     .clear(); 

        // MET is put in the first element
        Object_pt         .emplace_back( ((MissingET*) branchMET->At(0))->MET ); 
        Object_eta        .emplace_back( 0.                                   );
        Object_phi        .emplace_back( ((MissingET*) branchMET->At(0))->Phi );
        Object_energy     .emplace_back( 0.                                   ); 
        Object_charge     .emplace_back( 0.                                   ); 
        Object_isEle      .emplace_back( 0.                                   ); 
        Object_isMu       .emplace_back( 0.                                   ); 
        Object_isGamma    .emplace_back( 0.                                   ); 
        Object_isGluon    .emplace_back( 0.                                   ); 
        Object_isLight    .emplace_back( 0.                                   ); 
        Object_isCharm    .emplace_back( 0.                                   ); 
        Object_isBottom   .emplace_back( 0.                                   ); 

        PF_eta            .emplace_back( 0.                                                     );
        PF_phi            .emplace_back( ((MissingET*) branchMET->At(0))->Phi                   );
        PF_pt             .emplace_back( ((MissingET*) branchMET->At(0))->MET                   );
        PF_energy         .emplace_back( 0.                                                     );
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
        PF_JetIndex       .emplace_back( -1                                                     ); 

        for (const auto it : Objects) {

            Object_pt       .emplace_back( it.Object_pt       ); 
            Object_eta      .emplace_back( it.Object_eta      ); 
            Object_phi      .emplace_back( it.Object_phi      ); 
            Object_energy   .emplace_back( it.Object_energy   ); 
            Object_charge   .emplace_back( it.Object_charge   ); 
            Object_isEle    .emplace_back( it.Object_isEle    ); 
            Object_isMu     .emplace_back( it.Object_isMu     ); 
            Object_isGamma  .emplace_back( it.Object_isGamma  ); 
            Object_isGluon  .emplace_back( it.Object_isGluon  ); 
            Object_isLight  .emplace_back( it.Object_isLight  ); 
            Object_isCharm  .emplace_back( it.Object_isCharm  ); 
            Object_isBottom .emplace_back( it.Object_isBottom ); 

        }

        for (const auto it : PF_cands) {
            PF_eta            .emplace_back( it.PF_eta                        ); 
            PF_phi            .emplace_back( it.PF_phi                        );
            PF_pt             .emplace_back( it.PF_pt                         );
            PF_energy         .emplace_back( it.PF_energy                     );
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
            PF_JetIndex       .emplace_back( it.PF_JetIndex                   ); 
        }

        tree->Fill();

    }//Event++

    file->Write();
    file->Close();


    delete treeReader;
    delete chain;
}

//------------------------------------------------------------------------------
