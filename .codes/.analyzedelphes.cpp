#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string.h>
#include <vector>
#include <filesystem>

#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "ExRootTreeReader.h"
#include "ExRootClasses.h"

// #ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
// #endif

#include "utilities.hpp"


using namespace std;
using json = nlohmann::json;



int main()
{
    // Read the json file for parameters
    json conf = load_config("params.json");

    // define some useful variables
    bool bkg = conf["background"];
    TString base = conf["base_dir"];
    TString final_state = conf["final_state"];
    TString mass[] = {"10", "50", "100"};

    // Delphes files
    TString event_file[] = {
        base + "delphes_sig_10GeV.root",
        base + "delphes_sig_50GeV.root",
        base + "delphes_sig_100GeV.root"
    };

    // Background file
    TString background_file = base + "delphes_bg.root";

    // Observables
    TString obs[] = {
        "dphi_a1_j1",
        "dphi_j1_j2",
        "dphi_a1_j2",
        "dphi_a1_bbar",
        "dphi_a2_bbar",
        "dphi_j1_bbar",
        "cost_bbar"
    };

    // Histogram settings
    gStyle->SetOptStat(0);

    // Scanning through the mass values
    for (int f = 0; f < size(mass); f++)
    {
        // Output
        TString outpath = base + "histograms/m" + mass[f];
        string outpathString(outpath.Data());

        mkdir(outpathString.c_str(), 0774);

        cout << "\033[32m" << "  \xE2\x9C\x93 " << outpath << " created"  << "\033[0m" << endl;

        // Defining histograms
        TString labels = ";#delta #Phi;#frac{1}{#sigma} #frac{d#sigma}{d#Phi}";

        TH1F *h0 =new TH1F(obs[0] + "_m" + mass[f], labels, 20, -3.14159, 3.14159);
        TH1F *h1 =new TH1F(obs[1] + "_m" + mass[f], labels, 20, -3.14159, 3.14159);
        TH1F *h2 =new TH1F(obs[2] + "_m" + mass[f], labels, 20, -3.14159, 3.14159);
        TH1F *h3 =new TH1F(obs[3] + "_m" + mass[f], labels, 20, -3.14159, 3.14159);
        TH1F *h4 =new TH1F(obs[4] + "_m" + mass[f], labels, 20, -3.14159, 3.14159);  
        TH1F *h5 =new TH1F(obs[5] + "_m" + mass[f], labels, 20, -3.14159, 3.14159);
        TH1F *h6 =new TH1F(obs[6] + "_m" + mass[f], ";#cos #theta;#frac{1}{#sigma} #frac{d#sigma}{d#cos#theta}", 20, -1.0, 1.0);

        h0 -> Sumw2();
        h1 -> Sumw2();
        h2 -> Sumw2();
        h3 -> Sumw2();
        h4 -> Sumw2();
        h5 -> Sumw2();
        h6 -> Sumw2();

        // Put histograms in a vector
        TH1F *histograms[] = {h0, h1, h2, h3, h4, h5, h6};

        cout << "\033[32m" << "  \xE2\x9C\x93 " << "Histogram Array created" << "\033[0m" << endl;


        TChain chain("Delphes");
        chain.Add(event_file[f]);

        // Do background file
        if (bkg == true)
        {
            chain.Add(background_file);
        }

        // Objects to class ExRoot
        ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
        Long64_t numberOfEntries = treeReader -> GetEntries();

        // Pointers to branches used in the analysis
        TClonesArray *branchWeight = treeReader->UseBranch("Weight");
        Double_t weight;

        TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
        TClonesArray *branchMET = treeReader->UseBranch("MissingET");
        TClonesArray *branchElectron = treeReader->UseBranch("Electron");
        TClonesArray *branchMuon = treeReader->UseBranch("Muon");
        TClonesArray *branchJet = treeReader->UseBranch("Jet");

        // Now looping over all events
        for (Int_t entry = 0; entry < numberOfEntries; entry++)
        {
            // Loading branch with data
            treeReader -> ReadEntry(entry);

            LHEFWeight *W = (LHEFWeight*) branchWeight->At(0);
            weight = (W->Weight * 1000) / numberOfEntries;
            
            // Dealing with the to Photons
            Photon *gamma1, *gamma2;

            if (branchPhoton->GetEntries() > 1)
            {
                gamma1 = (Photon *) branchPhoton->At(0);
                gamma2 = (Photon *) branchPhoton->At(1);
            }

            // Neutrinos (Taken as the missing transverse energy)
            if (branchMET->GetEntries() > 0)
            {
                MissingET *met = (MissingET *) branchMET->At(0);
            }

            // Jets of which one is a bjet
            Jet *jet1, *jet2, *bjet;

            if (branchJet->GetEntries() > 1)
            {
                vector<Jet*> bJets, lightJets;

                for (Int_t j = 0; j < branchJet->GetEntriesFast(); j++)
                {
                    Jet *jet = (Jet*) branchJet->At(j);

                    if (jet->BTag == 1)
                    {
                        bJets.push_back(jet);
                    } else {
                        lightJets.push_back(jet); // The jets are already ordered in PT here
                    }
                }

                bjet = bJets[0];
                jet1 = lightJets[0];
                jet2 = lightJets[1];
            } else {
                // If there is only one jet, then it has to be a btagged jet
                bjet = (Jet*) branchJet->At(0);
            }

            // Electrons
            if (branchElectron->GetEntries() > 0)
            {
                Electron *electron = (Electron *) branchElectron->At(0);
            }
            
            // Muons
            if (branchMuon->GetEntries() > 0)
            {
                Muon *muon = (Muon *) branchMuon->At(0);
            }

            /**
             * Bosting the bjet into the top quark rest frame
             */

            // Doing the observables


            Double_t observables[8];
            
            if (final_state == "hadronic")
            {
                observables[0] = dphi_photon_jet(gamma1, jet1);
                observables[1] = dphi_jet_jet(jet1, jet2);
                observables[2] = dphi_photon_jet(gamma1, jet2);
                observables[3] = dphi_photon_jet(gamma1, bjet);
                observables[4] = dphi_photon_jet(gamma2, bjet);
                observables[5] = dphi_jet_jet(jet1, bjet);
                observables[6] = TMath::Cos(bjet->P4().Theta());
            }

            // Include else for semileptonic final state above


            // Applying cuts
            TLorentzVector axion = gamma1->P4() + gamma2->P4();

            bool cuts[] = {
                gamma1->PT > 45 && gamma2->PT > 5 && axion.M() > 7 && axion.M() < 13,    // m10
                gamma1->PT > 55 && axion.M() > 57 && axion.M() < 63,                      // m50 (parton level analysis these cuts are for m60)
                gamma1->PT > 75 && axion.M() > 97 && axion.M() < 103                       // m100
            };

            if (gamma1->Eta >= 0.5 && gamma1->Eta <= 5.0 && gamma2->Eta >= 0.5 && gamma2->Eta <= 5.0)
            {
                if (cuts[f])
                {
                    for (Int_t hg = 0; hg < size(histograms); hg++)
                    {
                        histograms[hg] -> Fill(observables[hg], weight);
                    }
                }
            }
        }

        cout << "\033[32m" << "  \xE2\x9C\x93 " << "Done filling histograms" << "\033[0m" << endl;

        // Get bin counts and save the histos in txt file
        Double_t sum = 0;
        string outfilepath;
        ofstream outputFile;

        for (int j = 0; j < size(histograms); j++)
        {
            if (bkg == true) {
                outputFile.open(outpath + "/" + obs[j] + "_m" + mass[f] +  "_background.txt", std::ios::out);
            } else {
                outputFile.open(outpath + "/" + obs[j] + "_m" + mass[f] +  ".txt", std::ios::out);
            }

            outputFile << "bin_center\tx" << endl;

            for (Int_t bin = 1; bin < histograms[j]->GetNbinsX() + 1; bin++)
            {
                outputFile << histograms[j]->GetXaxis()->GetBinCenter(bin);
                outputFile << "\t" << histograms[j]->GetBinContent(bin) << endl;

                sum += histograms[j]->GetBinContent(bin) * weight;
            } 
            outputFile.close();

            // cout << "\033[32m" << "  \xE2\x9C\x93 " << j << " \033[0m";
            cout << "      x-section: " << "\033[1;35m" << sum << "\033[0m" << endl;
            sum = 0;
        }
        
        // delete histograms;
        cout << "\033[32m" << "  \xE2\x9C\x93 " << "Deleting Histograms" << "\033[0m" << endl;
        delete h0, h1, h2, h3, h4, h5, h6;


        cout << "\033[1;36m***********************************************************************\033[0m" << endl;

    }

}