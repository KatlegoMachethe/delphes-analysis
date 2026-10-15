#include "utilities.hpp"
#include <stdexcept>
#include <cmath>


/**
 * load_config - reads configurations pertaining the chi ssqaured calculations.
 * @filename: The name of the .json file to read configurations from.
 * 
 * Return: Returns a json object.
 */
json load_config(const std::string& filename)
{
    std::ifstream f(filename);
    if (!f.is_open())
    {
        throw std::runtime_error("Could not open config file " + filename);
    }

    json config;
    config = json::parse(f);

    return config;
}

/**
 * calculate_dphi - Calculate dphi between two particles.
 * @first_p_v: The Lorentz vector of the first particle.
 * @second_p_v: The Lorentz vector of the second particle.
 * 
 * Return: Returns a Double_t with the dphi
 */
Double_t calculate_dphi(TLorentzVector p1_v, TLorentzVector p2_v)
{
    return p1_v.Phi() - p2_v.Phi();
}

/**
 * Apply exactly the b-jet rotation/boost sequence in analysis/newcode/histo.cpp.
 * The reconstructed top remains unrotated, and only the b-jet is transformed.
 */
TLorentzVector transform_bjet(TLorentzVector bjet, const TLorentzVector &top)
{
    // Rotatate and boost the 
    bjet.RotateZ(-top.Phi());
    bjet.RotateY(-top.Theta());
    bjet.Boost(-top.BoostVector());
    return bjet;
}

/**
 * dphi_photon_jet - Calculate dphi between photons and jets.
 * @photon: Photons properties.
 * @jet: Jet properties.
 * 
 * Return: Returns a Float_t with the dphi
 */
Double_t dphi_photon_jet(Photon *photon, Jet *jet)
{
    return photon->Phi - jet->Phi;
}

/**
 * dphi_jet_jet - Calculate dphi between photons and jets.
 * @first_jet: Photons properties.
 * @second_jet: Jet properties.
 * 
 * Return: Returns a Float_t with the dphi
 */
Double_t dphi_jet_jet(Jet *jet1, Jet *jet2)
{
    return jet1->Phi - jet2->Phi;
}
