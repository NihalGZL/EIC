#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <iostream>

#include <fastjet/ClusterSequence.hh>

int nentries;


void ReconstructedParticles() {


    const Int_t ndet = 2;
    const Int_t ntype = 2;
    const Int_t max_track = 10000;

    TFile* file = new TFile("pythia8NCDIS_10x100_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_vtxfix_1.0002.eicrecon.tree.edm4eic.root", "READ");
    
    TTree* events = (TTree*)file->Get("events");

    nentries = events->GetEntries();


  // Define the EEC histogram
    
    int nBinsEEC = 100;
    double minEEC = 0.0;
    double maxEEC = 2.0;
    
    TH1F *hist_EEC = new TH1F("EEC", "Energy-Energy Correlation Function", nBinsEEC, minEEC, maxEEC);
  
    TH2F *hist_EEC_R = new TH2F("hist_EEC_R", "hist_EEC_R", 100, 0.0, 2.5, 100, 0, 100);



    // Set branch addresses  momentum


     Float_t pmc[3][max_track];
     
    events->SetBranchAddress("ReconstructedChargedParticles.momentum.x", pmc[0]);
    events->SetBranchAddress("ReconstructedChargedParticles.momentum.y", pmc[1]);
    events->SetBranchAddress("ReconstructedChargedParticles.momentum.z", pmc[2]);


    // Create a vector to hold input particles for FastJet clustering
    
    std::vector<fastjet::PseudoJet> input_particles;

    for (int i = 0; i < nentries; i++) {
        if ((i + 1) % 100 == 0)
            std::cout << "Processing entry == " << i + 1 << " == out of " << nentries << ".\n";

        events->GetEntry(i);

        // Loop through EcalEndcapPClusters and create PseudoJets using ReconstructedChargedParticles energy and momentum
        
       
            
            for (Long64_t j = 0; j < events->GetLeaf("ReconstructedChargedParticles.energy")->GetLen(); j++) {
        
            TLeaf* leaf_Recon_particles_energy = events->GetLeaf("ReconstructedChargedParticles.energy");
            
            double energy = leaf_Recon_particles_energy->GetValue(j);
            
                       

            // Extract momentum components for each particle
            
            Float_t Px = pmc[0][j];
            Float_t Py = pmc[1][j];
            Float_t Pz = pmc[2][j];

            // Create a PseudoJet using momentum components
            
            fastjet::PseudoJet particle(Px, Py, Pz, energy);

            // Add the particle to the input_particles vector
            
            input_particles.push_back(particle);
        }
    }

    // Perform FastJet jet clustering
    
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4); // Adjust the radius parameter as needed
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    std::vector<fastjet::PseudoJet> jets = clust_seq.inclusive_jets(20.0); // Adjust the minimum transverse momentum as needed

    // Loop through the resulting jets and access their properties
    
    for (const fastjet::PseudoJet& jet : jets) {
    
        double jet_px = jet.px();
       
      

        // Calculate EEC for the jet
        
        double EEC = 0.0;
        double totalEnergy = jet.E();
        double deltaR = 0.0;
        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
        

        for (unsigned k1 = 0; k1 < constituents.size(); k1++) {
            double energy1 = constituents[k1].E();
            double pt1 = constituents[k1].pt();
            double eta1 = constituents[k1].eta();
            


            for (unsigned k2 = k1 + 1; k2 < constituents.size(); k2++) {
            
                double energy2 = constituents[k2].E();
                double pt2 = constituents[k2].pt();
                double eta2 = constituents[k2].eta();
                
                 

                // Calculate the angular separation
                
                deltaR = sqrt((eta1 - eta2) * (eta1 - eta2) + (constituents[k1].phi() - constituents[k2].phi()) * (constituents[k1].phi() - constituents[k2].phi()));

                // Calculate the energy-energy correlation term
                EEC += energy1 * energy2 * pow(deltaR, 2);
            }
            
            // Normalize the EEC
            
        EEC /= (totalEnergy * totalEnergy);

        // Fill the EEC histogram
        
        hist_EEC->Fill(EEC);
        
        hist_EEC_R->Fill(deltaR,EEC);
        
        
        }
        
      }  
        
    TFile *eecFile = new TFile("ECC_histograms.root", "recreate");
        
    hist_EEC->Write();
    
    hist_EEC_R->Write();
    
    eecFile->Close();
    
    file->Close();

    }
   
   
      

  int main() {
  
    ReconstructedParticles();
    
    return 0;
}

