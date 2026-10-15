#ifndef UTILITIES_CPP
#define UTILITIES_CPP

#include <iostream>
#include <string.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <vector>

#include"TH1F.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"


using json = nlohmann::json;
using namespace std;

json load_config(const std::string& filename);
Double_t calculate_dphi(TLorentzVector, TLorentzVector);
TLorentzVector transform_bjet(TLorentzVector, const TLorentzVector &);

#endif
