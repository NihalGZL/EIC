#include <iostream>
#include <cmath>
#include <algorithm>
#include "fastjet/ClusterSequence.hh"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <vector>
#include "Pythia8/Pythia.h"
#include "TParticle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TProfile.h"


using namespace std;
using namespace fastjet;
using namespace Pythia8;

void print(std::vector<int> const &a) {
    std::cout << "The vector elements are: ";

    for (int i = 0; i < a.size(); i++)
        std::cout << a.at(i) << ' ';
    cout << endl;
}

int main() {

    float min_eta_jet = -4.0;
    float max_eta_jet = 4.0;
    
  //  double pTjetMin = 10.0, pTjetMax = 50.0;
  
      double pTjetMin = 2.0, pTjetMax = 20.0;
    
   
    
    
    // Define the EEC histogram
    
    int nBinsEEC = 100;
    double minEEC = 0.0;
    double maxEEC = 2.0;
    TH1F *hist_EEC = new TH1F("EEC", "Energy-Energy Correlation Function", nBinsEEC, minEEC, maxEEC);
   
   
   
    TH2F *hist_E_vs_eta = new TH2F("Jet_E_vs_eta", "Jet Energy vs. Jet Eta", 100, -2.5, 2.5, 100, 0, 100);

    TCanvas *canvas = new TCanvas("Canvas", "Jet Energy vs. Jet Eta", 800, 600);
    
    //float min_eta_jet, max_eta_jet;
    //double pTjetMin, pTjetMax;
    
    ofstream myFile("Jet_Flavor.csv");
    myFile << "Jet charge multiplicity," << " Jet charge," << " Flavor (1:down, 2:up, 3:charm, 4:beauty)\n";
    TFile *output = new TFile("tutorial5.root", "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    // Define histograms for jet pT and eta
    
   
    TH1F *hist_pT = new TH1F("Jet_pT", "Jet pT Distribution", 100, 0, 30);
    TH1F *hist_eta = new TH1F("Jet_eta", "Jet Eta Distribution", 100, -10, 10);
    
    
    TH1F *hist_mass = new TH1F("hist_mass", "hist_mass", 100, 0, 20);
    TH1F *hist_jet_energy = new TH1F("hist_jet_energy", "hist_jet_energy", 100, 0, 80);

 // TH2F *hist_EEC_R = new TH2F("hist_EEC_R", "hist_EEC_R", 100, 0.02, 0.5, 100, 0, 2.0);  
    
   TProfile *hist_EEC_R = new TProfile("hist_EEC_R","hist_EEC_R",100,0.02,1.0,0,2.0,""); 
    
    
    int id, event, size, no, mother1, mother2, nevents, isfinal, mother_id;
    
    vector<int> motherlist;
    
    double m, px, py, pz, En, Q2min = 10, charge;
    
   
    
    tree->Branch("event", &event, "event/I");
    tree->Branch("size", &size, "size/I");
    tree->Branch("no", &no, "no/I");
    tree->Branch("id", &id, "id/I");
    tree->Branch("m", &m, "m/D");
    tree->Branch("px", &px, "px/D");
    tree->Branch("py", &py, "py/D");
    tree->Branch("pz", &pz, "pz/D");
    tree->Branch("mother1", &mother1, "mother1/I");
    tree->Branch("mother2", &mother2, "mother2/I");
    tree->Branch("isfinal", &isfinal, "isfinal/I");
    tree->Branch("motherlist", &motherlist);
    tree->Branch("En", &En, "En/D");
    tree->Branch("mother_id", &mother_id, "mother_id/I");
    tree->Branch("charge", &charge, "charge/D");
   // double eProton = 920.;
   //double eElectron = 27.5;
   
   // double R = 0.8;
   
      double R = 1.0;

    Pythia8::Pythia pythia;
    pythia.readString("Beams:idA = 11"); // electron 
    pythia.readString("Beams:idB = 2212");
  //  pythia.readString("Beams:eCM = 1300");
    pythia.readString("211:mayDecay = false");
    pythia.readString("-211:mayDecay = false");
    pythia.readString("111:mayDecay = false");
    pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
    pythia.readString("PDF:lepton = off");
    pythia.readString("TimeShower:QEDshowerByL = off");
    pythia.readString("Beams:LHEF = void");
    pythia.readString("HadronLevel:all = on");
    pythia.readString("Beams:eA = 18.");
    pythia.readString("Beams:eB = 275.");
    
    pythia.readString("HardQCD:all = on");  // Turn on hard QCD processes
  //  pythia.readString("WeakSingle:all = on");  // Enable weak decays
    
     

    pythia.init();
    
    JetDefinition jet_def(antikt_algorithm, R);
    vector<vector<PseudoJet>> events;

    vector<double> charge1;
    nevents = 1e7;
    double ka = 0.3;
    for (int i = 0; i < nevents; i++) {
        if (!pythia.next()) continue;
        pythia.settings.parm("Random:setSeed = on");
        pythia.settings.parm("Random:seed = 1000");

        int entries = pythia.event.size();

        std::cout << "Event: " << i << std::endl;
        event = i;
        size = entries;
        int e = 0;
        vector<PseudoJet> parts;
        for (int j = 0; j < entries; j++) {
            id = pythia.event[j].id();
            no = j;
            mother1 = pythia.event[j].mother1();
            mother2 = pythia.event[j].mother2();
            isfinal = pythia.event[j].isFinal();
            charge = pythia.event[j].charge();
            m = pythia.event[j].m();
            px = pythia.event[j].px();
            py = pythia.event[j].py();
            pz = pythia.event[j].pz();
            En = pythia.event[j].e();
            if (isfinal == 1) {
            
            
                PseudoJet part(px, py, pz, En);
                part.set_user_index(j);
                parts.push_back(part);
                
            
            
       
                
           }
        }

         // Anti-kT Algorithm


        ClusterSequence clust(parts, jet_def);
        vector<PseudoJet> jets = sorted_by_pt(clust.inclusive_jets(pTjetMin));

        for (unsigned j = 0; j < jets.size(); j++) {
        
          if (jets[j].pt() < pTjetMin ||
                jets[j].pt() > pTjetMax ||
                jets[j].eta() > max_eta_jet ||
                jets[j].eta() < min_eta_jet) continue;
                
              
                // Calculate jet pT and eta
        double jet_pT = jets[j].pt();
        double jet_eta = jets[j].eta();

        // Fill histograms
        
        hist_pT->Fill(jet_pT);
        hist_eta->Fill(jet_eta);
       
        // Fill the 2D histogram
        
    hist_E_vs_eta->Fill(jet_eta, jets[j].E());  
    
   
                       
   vector<PseudoJet> constituents = jets[j].constituents();
   
   
   // Calculate EEC for the jet
        double EEC = 0.0;
        double totalEnergy = jets[j].E();
        double deltaR = 0.0;

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
        }

        // Normalize the EEC
        EEC /= (totalEnergy * totalEnergy);

        // Fill the EEC histogram
        hist_EEC->Fill(EEC);
        hist_EEC_R->Fill(deltaR,EEC);
            
            
            // Initialize the total energy and momentum of the jet
    double totalE = 0.0;
    double totalPx = 0.0;
    double totalPy = 0.0;
    double totalPz = 0.0;

           for (const PseudoJet& constituent : constituents) {
        totalE += constituent.E();
        totalPx += constituent.px();
        totalPy += constituent.py();
        totalPz += constituent.pz();
    }

           double jetMass = sqrt(totalE * totalE - totalPx * totalPx - totalPy * totalPy - totalPz * totalPz);
              
             hist_mass->Fill(jetMass);
             hist_jet_energy->Fill(totalE);
        
            
            
            double Q_k = 0.;

            int up_count = 0;
            int down_count = 0;
            int charm_count = 0;
            int beauty_count = 0;

            double multiplicity = 0;
            double whole_multiplicity = 0;
            vector<int> mothers1;
            vector<int> mothers2;

            for (unsigned k = 0; k < constituents.size(); k++) {
                double pTJet = jets[j].pt();
                int p_index = constituents[k].user_index();
                double Q_h = pythia.event[p_index].charge();
                double id1 = pythia.event[p_index].id();
               
                if (abs(id1) == 4) { // Charm quark
                    charm_count++;
                } else if (abs(id1) == 5) { // Beauty quark
                    beauty_count++;
                } else if (abs(id1) == 2) { // Up quark
                    up_count++;
                } else if (abs(id1) == 1) { // Down quark
                    down_count++;
                }
                
                if (Q_h > 0) {
                    if (abs(id1) == 211 || abs(id1) == 111) {
                        multiplicity++;
                    } else {
                        continue;
                    }
                }

                whole_multiplicity++;
                double pt_h = constituents[k].pt();
                double z_h = pt_h / pTJet;
                Q_k += pow(z_h, ka) * Q_h;
            }




            int flavor = 0; // Flavor code (1: down, 2: up, 3: charm, 4: beauty)
            if (charm_count > 0) {
                flavor = 3;
            } else if (beauty_count > 0) {
                flavor = 4;
            } else if (up_count > 0) {
                flavor = 2;
            } else if (down_count > 0) {
                flavor = 1;
            }

            myFile << multiplicity << ", " << Q_k << ", " << flavor << "\n";

            tree->Fill();
        }
    }



   TFile *histFile = new TFile("jet_histograms.root", "recreate");
    hist_pT->Write();
    hist_eta->Write();
    hist_mass->Write();
    hist_jet_energy->Write();
  
    histFile->Close();
    
    TFile *eecFile = new TFile("eec_histogram.root", "recreate");
    hist_EEC->Write();
    hist_EEC_R->Write();
    eecFile->Close();
    
    TFile *hist2DFile = new TFile("jet_energy_vs_eta.root", "recreate");
    hist_E_vs_eta->Write();
    hist2DFile->Close();
    
    output->Write();
    output->Close();
    myFile.close();
    return 0;
}

