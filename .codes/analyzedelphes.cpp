#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "utilities.hpp"


using namespace std;
using json = nlohmann::json;



int main(int argc, char *argv[]) try
{
    if (argc > 2)
    {
        cerr << "Usage: " << argv[0] << " [params.json]" << endl;
        return 1;
    }

    // Read the json file for parameters
    json conf = load_config(argc == 2 ? argv[1] : "params.json");

    // define some useful variables
    bool bkg = conf["background"];
    filesystem::path base = conf["base_dir"];
    TString final_state = conf["final_state"];
    TString mass[] = {"10", "50", "100"};

    // Delphes files
    vector<string> signal_files = conf["signal_files"];

    TString event_file[] = {
        (base / signal_files[0]).c_str(),
        (base / signal_files[1]).c_str(),
        (base / signal_files[2]).c_str()
    };

    // Background file
    TString background_file = (base / conf["background_file"]).c_str();

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

    // Set this before ROOT initializes its interpreter and loads the dictionary.
    string rootIncludePath = DELPHES_DIR ":" DELPHES_DIR "/external";
    const char *existingIncludePath = gSystem->Getenv("ROOT_INCLUDE_PATH");

    if (existingIncludePath && *existingIncludePath)
    {
        rootIncludePath += string(":") + existingIncludePath;
    }
    gSystem->Setenv("ROOT_INCLUDE_PATH", rootIncludePath.c_str());

    // Histogram settings
    gStyle->SetOptStat(0);

    // Scanning through the mass values
    for (size_t f = 0; f < size(mass); f++)
    {
        bool isBackground = f == size(mass);
        const TString input_file = isBackground ? background_file : event_file[f];
        const TString sampleLabel = isBackground ? TString("background") : TString("m") + mass[f];

        // Output
        TString outpath = (base / "histograms" / sampleLabel.Data()).c_str();
        filesystem::create_directories(outpath.Data());

        cout << "\033[32m" << "  \xE2\x9C\x93 " << outpath << " created"  << "\033[0m" << endl;

        // Defining histograms
        TString labels = ";#delta #Phi;#frac{1}{#sigma} #frac{d#sigma}{d#Phi}";

        TH1F *h0 =new TH1F(obs[0] + "_" + sampleLabel, labels, 20, -3.14159, 3.14159);
        TH1F *h1 =new TH1F(obs[1] + "_" + sampleLabel, labels, 20, -3.14159, 3.14159);
        TH1F *h2 =new TH1F(obs[2] + "_" + sampleLabel, labels, 20, -3.14159, 3.14159);
        TH1F *h3 =new TH1F(obs[3] + "_" + sampleLabel, labels, 20, -3.14159, 3.14159);
        TH1F *h4 =new TH1F(obs[4] + "_" + sampleLabel, labels, 20, -3.14159, 3.14159);
        TH1F *h5 =new TH1F(obs[5] + "_" + sampleLabel, labels, 20, -3.14159, 3.14159);
        TH1F *h6 =new TH1F(obs[6] + "_" + sampleLabel, ";#cos #theta;#frac{1}{#sigma} #frac{d#sigma}{d#cos#theta}", 20, -1.0, 1.0);

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
        chain.Add(input_file);

        // Objects to class ExRoot
        ExRootTreeReader treeReader(&chain);
        Long64_t numberOfEntries = treeReader.GetEntries();
        
        cout << "Input: " << input_file << " (" << numberOfEntries << " entries; "
             << sampleLabel << (isBackground ? ", common selection)" : " selection)") << endl;

        // Pointers to branches used in the analysis
        TClonesArray *branchEvent = treeReader.UseBranch("Event");
        TClonesArray *branchPhoton = treeReader.UseBranch("Photon");
        TClonesArray *branchJet = treeReader.UseBranch("Jet");

        Long64_t missingPhotons = 0, missingJets = 0, selectedEvents = 0;

        // Now looping over all events
        for (Long64_t entry = 0; entry < numberOfEntries; entry++)
        {
            treeReader.ReadEntry(entry);
            auto *event = dynamic_cast<const HepMCEvent *>(branchEvent->At(0));
            Double_t weight = (event->Weight * 1000) / numberOfEntries;
            

            // Photons.
            if (branchPhoton->GetEntriesFast() < 2)
            {
                ++missingPhotons;
                continue;
            }
            Photon *gamma1 = (Photon*) branchPhoton->At(0);
            Photon *gamma2 = (Photon*) branchPhoton->At(1);

            // The hadronic observables require a b-tagged jet and two light jets.
            vector<Jet *> bJets, lightJets;

            for (Int_t j = 0; j < branchJet->GetEntriesFast(); j++)
            {
                Jet *jet = (Jet *) branchJet->At(j);
                if (jet->BTag == 1)
                {
                    bJets.push_back(jet);
                } else {
                    lightJets.push_back(jet); // Already ordered in PT.
                }
            }

            if (bJets.empty() || lightJets.size() < 2)
            {
                ++missingJets;
                continue;
            }
            Jet *bjet = bJets[0], *jet1 = lightJets[0], *jet2 = lightJets[1];

            TLorentzVector top = bjet->P4() + jet1->P4() + jet2->P4();
            TLorentzVector boosted_bjet = transform_bjet(bjet->P4(), top);

            // Now get observables
            Double_t observables[] = {
                calculate_dphi(gamma1->P4(), jet1->P4()),
                calculate_dphi(jet1->P4(), jet2->P4()),
                calculate_dphi(gamma1->P4(), jet2->P4()),
                calculate_dphi(gamma1->P4(), boosted_bjet),
                calculate_dphi(gamma2->P4(), boosted_bjet),
                calculate_dphi(jet1->P4(), boosted_bjet),
                boosted_bjet.CosTheta()
            };

            // Applying cuts
            TLorentzVector axion = gamma1->P4() + gamma2->P4();

            bool cuts[] = {
                gamma1->PT > 5 && gamma2->PT > 5 && axion.M() > 7 && axion.M() < 13,    // m10
                gamma1->PT > 5 && gamma2->PT > 5 && axion.M() > 47 && axion.M() < 53,   // m50
                gamma1->PT > 5 && gamma2->PT > 5 && axion.M() > 97 && axion.M() < 103   // m100
            };

            if (gamma1->Eta >= 0.5 && gamma1->Eta <= 5.0 && gamma2->Eta >= 0.5 && gamma2->Eta <= 5.0)
            {
                if (cuts[f])
                {
                    ++selectedEvents;
                    for (size_t hg = 0; hg < size(histograms); hg++)
                    {
                        histograms[hg] -> Fill(observables[hg], weight);
                    }
                }
            }
        }

        cout << "\033[32m" << "  \xE2\x9C\x93 " << "Done filling histograms" << "\033[0m" << endl;
        cout << "Events: " << numberOfEntries << ", fewer than two photons: " << missingPhotons
             << ", missing required jets after photon selection: " << missingJets
             << ", passed all cuts: " << selectedEvents << endl;
        if (selectedEvents == 0)
        {
            cout << "Warning: no events passed the selection; output histograms are empty." << endl;
        }

        // Get bin counts and save the histos in txt file
        Double_t sum = 0;
        ofstream outputFile;
        outputFile.exceptions(ios::failbit | ios::badbit);

        for (size_t j = 0; j < size(histograms); j++)
        {
            outputFile.open(outpath + "/" + obs[j] + "_" + sampleLabel + ".txt", std::ios::out);

            outputFile << "bin_center\tx" << endl;

            for (Int_t bin = 1; bin < histograms[j]->GetNbinsX() + 1; bin++)
            {
                outputFile << histograms[j]->GetXaxis()->GetBinCenter(bin);
                outputFile << "\t" << histograms[j]->GetBinContent(bin) << endl;

                // Bin contents already include the event weights.
                sum += histograms[j]->GetBinContent(bin);
            } 
            outputFile.close();

            // cout << "\033[32m" << "  \xE2\x9C\x93 " << j << " \033[0m";
            cout << "      " << obs[j] << " integral: " << "\033[1;35m" << sum << "\033[0m"
                 << " (underflow=" << histograms[j]->GetBinContent(0)
                 << ", overflow=" << histograms[j]->GetBinContent(histograms[j]->GetNbinsX() + 1)
                 << ")" << endl;
            sum = 0;
        }
        
        cout << "\033[32m" << "  \xE2\x9C\x93 " << "Deleting Histograms" << "\033[0m" << endl;
        for (auto *histogram : histograms)
        {
            delete histogram;
        }


        cout << "\033[1;36m***********************************************************************\033[0m" << endl;
    }

    return 0;
}
catch (const exception &error)
{
    cerr << "Error: " << error.what() << endl;
    return 1;
}
