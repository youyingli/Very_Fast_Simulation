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

    double Jet_PF_eta;
    double Jet_PF_phi;

    double Jet_PF_log_pt;
    double Jet_PF_log_energy;
    double Jet_PF_log_rel_pt;
    double Jet_PF_log_rel_energy;
    double Jet_PF_deltaEta;
    double Jet_PF_deltaPhi;
    double Jet_PF_charge;

    double Jet_PF_isEle;
    double Jet_PF_isMu;
    double Jet_PF_isChargedHad;
    double Jet_PF_isGamma;
    double Jet_PF_isNeutralHad;

    double Jet_PF_tanhd0;
    double Jet_PF_tanhdz;
    double Jet_PF_sigmad0;
    double Jet_PF_sigmadz;

};

//------------------------------------------------------------------------------

struct pt_sortor {
    bool operator()( const InputFeature &pf1, const InputFeature &pf2 ) {
        return ( pf1.Jet_PF_log_pt > pf2.Jet_PF_log_pt );
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

void jetExtractor(const char *inputFile)
{
    gSystem->Load("libDelphes");


    double Jet_pt        ;
    double Jet_eta       ;
    double Jet_phi       ;
    double Jet_energy    ;
    int Jet_flavour      ;
    int Jet_multiplicity ;

    std::vector<InputFeature> PFSet;
    std::vector<double> Jet_PF_eta           ; 
    std::vector<double> Jet_PF_phi           ; 
    std::vector<double> Jet_PF_log_pt        ; 
    std::vector<double> Jet_PF_log_energy    ; 
    std::vector<double> Jet_PF_log_rel_pt    ; 
    std::vector<double> Jet_PF_log_rel_energy; 
    std::vector<double> Jet_PF_deltaEta      ; 
    std::vector<double> Jet_PF_deltaPhi      ; 
    std::vector<double> Jet_PF_charge        ; 
    std::vector<double> Jet_PF_isEle         ; 
    std::vector<double> Jet_PF_isMu          ; 
    std::vector<double> Jet_PF_isChargedHad  ; 
    std::vector<double> Jet_PF_isGamma       ; 
    std::vector<double> Jet_PF_isNeutralHad  ; 
    std::vector<double> Jet_PF_tanhd0        ;
    std::vector<double> Jet_PF_tanhdz        ;
    std::vector<double> Jet_PF_sigmad0       ;
    std::vector<double> Jet_PF_sigmadz       ;

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");


    TFile* file = TFile::Open("output_jet.root", "recreate");
    TTree* tree = new TTree("jetTagging", "");

    tree->Branch( "Jet_pt"            , &Jet_pt,           "Jet_pt/D"             );
    tree->Branch( "Jet_eta"           , &Jet_eta,          "Jet_eta/D"            );
    tree->Branch( "Jet_phi"           , &Jet_phi,          "Jet_phi/D"            );
    tree->Branch( "Jet_energy"        , &Jet_energy,       "Jet_energy/D"         );
    tree->Branch( "Jet_flavour"       , &Jet_flavour,      "Jet_flavour/I"        );
    tree->Branch( "Jet_multiplicity"  , &Jet_multiplicity, "Jet_multiplicity/I"   );

    tree->Branch( "Jet_PF_eta"            , &Jet_PF_eta            );
    tree->Branch( "Jet_PF_phi"            , &Jet_PF_phi            );
    tree->Branch( "Jet_PF_log_pt"         , &Jet_PF_log_pt         );
    tree->Branch( "Jet_PF_log_energy"     , &Jet_PF_log_energy     );
    tree->Branch( "Jet_PF_log_rel_pt"     , &Jet_PF_log_rel_pt     );
    tree->Branch( "Jet_PF_log_rel_energy" , &Jet_PF_log_rel_energy );
    tree->Branch( "Jet_PF_deltaEta"       , &Jet_PF_deltaEta       );
    tree->Branch( "Jet_PF_deltaPhi"       , &Jet_PF_deltaPhi       );
    tree->Branch( "Jet_PF_charge"         , &Jet_PF_charge         );
    tree->Branch( "Jet_PF_isEle"          , &Jet_PF_isEle          );
    tree->Branch( "Jet_PF_isMu"           , &Jet_PF_isMu           );
    tree->Branch( "Jet_PF_isChargedHad"   , &Jet_PF_isChargedHad   );
    tree->Branch( "Jet_PF_isGamma"        , &Jet_PF_isGamma        );
    tree->Branch( "Jet_PF_isNeutralHad"   , &Jet_PF_isNeutralHad   );
    tree->Branch( "Jet_PF_tanhd0"         , &Jet_PF_tanhd0         );
    tree->Branch( "Jet_PF_tanhdz"         , &Jet_PF_tanhdz         );
    tree->Branch( "Jet_PF_sigmad0"        , &Jet_PF_sigmad0        );
    tree->Branch( "Jet_PF_sigmadz"        , &Jet_PF_sigmadz        );


    // Loop over all events
    for(int entry = 0; entry < treeReader->GetEntries(); ++entry) {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);


        bool anomaly = false;
        // Loop over all jets in event
        for(int i = 0; i < branchJet->GetEntriesFast(); ++i) {

            auto jet = (Jet*) branchJet->At(i);

            if ( jet->PT < 30. ) continue;
            if ( fabs(jet->Eta) > 2.4 ) continue;
            if ( jet->NCharged < 1 ) continue;
            if ( jet->Constituents.GetEntriesFast() < 2 ) continue;

            TLorentzVector jet_p4;
            jet_p4.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);


            Jet_pt            = jet_p4.Pt();
            Jet_eta           = jet_p4.Eta();
            Jet_phi           = jet_p4.Phi();
            Jet_energy        = jet_p4.E();
            Jet_multiplicity  = jet->NCharged + jet->NNeutrals;

            int pdgid_tmp = -1;
            double min_deltaR = 9999.;
            for (int it_genP = 0; it_genP < branchParticle->GetEntriesFast(); ++it_genP) {
                auto genP = (GenParticle*)branchParticle->At(it_genP);


                if (abs(genP->PID) > 5 && abs(genP->PID) != 21) continue;
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

            if (pdgid_tmp == -1) continue;
            Jet_flavour = pdgid_tmp;

            PFSet.clear();

            Jet_PF_eta           .clear();
            Jet_PF_phi           .clear();
            Jet_PF_log_pt        .clear();
            Jet_PF_log_energy    .clear();
            Jet_PF_log_rel_pt    .clear();
            Jet_PF_log_rel_energy.clear();
            Jet_PF_deltaEta      .clear();
            Jet_PF_deltaPhi      .clear();
            Jet_PF_charge        .clear();
            Jet_PF_isEle         .clear();
            Jet_PF_isMu          .clear();
            Jet_PF_isChargedHad  .clear();
            Jet_PF_isGamma       .clear();
            Jet_PF_isNeutralHad  .clear();
            Jet_PF_tanhd0        .clear();
            Jet_PF_tanhdz        .clear();
            Jet_PF_sigmad0       .clear();
            Jet_PF_sigmadz       .clear();

            int nCH = 0;
            // Loop over all jet's constituents
            for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j) {
                auto object = jet->Constituents.At(j);

                // Check if the constituent is accessible
                if(object == 0) continue;

                InputFeature feature;

                feature.Jet_PF_charge        = 0.;
                feature.Jet_PF_isEle         = 0.;
                feature.Jet_PF_isMu          = 0.;
                feature.Jet_PF_isChargedHad  = 0.;
                feature.Jet_PF_isGamma       = 0.;
                feature.Jet_PF_isNeutralHad  = 0.;
                feature.Jet_PF_tanhd0        = 0.; 
                feature.Jet_PF_tanhdz        = 0.; 
                feature.Jet_PF_sigmad0       = 0.; 
                feature.Jet_PF_sigmadz       = 0.; 

                // Charged Particle
                if(object->IsA() == Track::Class()){
                    auto track = (Track*) object;
                    auto particle = (GenParticle*) track->Particle.GetObject();

                    TLorentzVector track_p4;
                    track_p4.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, particle->Mass);

                    feature.Jet_PF_eta = track_p4.Eta();
                    feature.Jet_PF_phi = track_p4.Phi();

                    feature.Jet_PF_log_pt            = TMath::Log( track_p4.Pt() );
                    feature.Jet_PF_log_energy        = TMath::Log( track_p4.E() );
                    feature.Jet_PF_log_rel_pt        = TMath::Log( track_p4.Pt() / jet_p4.Pt() );
                    feature.Jet_PF_log_rel_energy    = TMath::Log( track_p4.E() / jet_p4.E() );
                    feature.Jet_PF_deltaEta          = track_p4.Eta() - jet_p4.Eta();
                    feature.Jet_PF_deltaPhi          = track_p4.DeltaPhi(jet_p4);
                    feature.Jet_PF_charge            = track->Charge;
                    feature.Jet_PF_tanhd0            = TMath::TanH(track->D0);
                    feature.Jet_PF_tanhdz            = TMath::TanH(track->DZ);
                    feature.Jet_PF_sigmad0           = track->ErrorD0;
                    feature.Jet_PF_sigmadz           = track->ErrorDZ;

                    if (abs(particle->PID) == 11) {
                        // Is Electron
                        feature.Jet_PF_isEle = 1.;
                    } else if (abs(particle->PID) == 13) {
                        // Is Muon
                        feature.Jet_PF_isMu = 1.;
                    } else {
                        // Is Charged Hadron
                        nCH++;
                        feature.Jet_PF_isChargedHad = 1.;
                    }

                }else if(object->IsA() == Tower::Class()){
                    auto tower = (Tower*) object;
                    auto particle = (GenParticle*) tower->Particles.At(0);

                    auto tower_p4 = UseEtaPhiEMass(tower->Eta, tower->Phi, tower->E, particle->Mass);

                    feature.Jet_PF_eta = tower_p4.Eta();
                    feature.Jet_PF_phi = tower_p4.Phi();

                    feature.Jet_PF_log_pt            = TMath::Log( tower_p4.Pt()               );
                    feature.Jet_PF_log_energy        = TMath::Log( tower_p4.E()                );
                    feature.Jet_PF_log_rel_pt        = TMath::Log( tower_p4.Pt() / jet_p4.Pt() );
                    feature.Jet_PF_log_rel_energy    = TMath::Log( tower_p4.E() / jet_p4.E()   );
                    feature.Jet_PF_deltaEta          = tower_p4.Eta() - jet_p4.Eta();
                    feature.Jet_PF_deltaPhi          = tower_p4.DeltaPhi(jet_p4);


                    if (particle->PID == 22) {
                        // Is Photon
                        feature.Jet_PF_isGamma = 1.;
                    } else {
                        // Is Neutral Hadron
                        feature.Jet_PF_isNeutralHad  = 1.;
                    }
                }
                if(feature.Jet_PF_eta < -1000) anomaly = true;
                PFSet.emplace_back(feature);

            }//Constituents++

            if (nCH == 0 || anomaly) continue;

            std::sort(PFSet.begin(), PFSet.end(), pt_sortor());

            for (const auto it : PFSet) {

                Jet_PF_eta            .emplace_back( it.Jet_PF_eta           ); 
                Jet_PF_phi            .emplace_back( it.Jet_PF_phi           );
                Jet_PF_log_pt         .emplace_back( it.Jet_PF_log_pt        );
                Jet_PF_log_energy     .emplace_back( it.Jet_PF_log_energy    );
                Jet_PF_log_rel_pt     .emplace_back( it.Jet_PF_log_rel_pt    );
                Jet_PF_log_rel_energy .emplace_back( it.Jet_PF_log_rel_energy);
                Jet_PF_deltaEta       .emplace_back( it.Jet_PF_deltaEta      );
                Jet_PF_deltaPhi       .emplace_back( it.Jet_PF_deltaPhi      );
                Jet_PF_charge         .emplace_back( it.Jet_PF_charge        );
                Jet_PF_isEle          .emplace_back( it.Jet_PF_isEle         );
                Jet_PF_isMu           .emplace_back( it.Jet_PF_isMu          );
                Jet_PF_isChargedHad   .emplace_back( it.Jet_PF_isChargedHad  );
                Jet_PF_isGamma        .emplace_back( it.Jet_PF_isGamma       );
                Jet_PF_isNeutralHad   .emplace_back( it.Jet_PF_isNeutralHad  );
                Jet_PF_tanhd0         .emplace_back( it.Jet_PF_tanhd0        ); 
                Jet_PF_tanhdz         .emplace_back( it.Jet_PF_tanhdz        ); 
                Jet_PF_sigmad0        .emplace_back( it.Jet_PF_sigmad0       ); 
                Jet_PF_sigmadz        .emplace_back( it.Jet_PF_sigmadz       ); 

            }

            tree->Fill();

        }//Jet++


    }//Event++

    file->Write();
    file->Close();


    delete treeReader;
    delete chain;
}

//------------------------------------------------------------------------------
