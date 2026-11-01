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
    return p1_v.DeltaPhi(p2_v);
}

/**
 * Apply exactly the b-jet rotation/boost sequence in analysis/newcode/histo.cpp.
 * The reconstructed top remains unrotated, and only the b-jet is transformed.
 */
TLorentzVector transform_bjet(TLorentzVector bjet, const TLorentzVector &top)
{
    // Rotatate and boost the
    // bjet.Boost(-top.BoostVector());
    // bjet.RotateZ(-top.Phi());
    // bjet.RotateY(-top.Theta());


    // Ammendments from Karabo Mosala
    // Rotate the top so that it points along +z
    TLorentzVector t = top;
    // t.RotateZ(-top.Phi());
    // t.RotateY(-top.Theta());

    // Apply the same rotation to the b-jet
    // bjet.RotateZ(-top.Phi());
    // bjet.RotateY(-top.Theta());

    // Boost with the ROTATED top velocity (now along +z)
    bjet.Boost(-t.BoostVector());

    return bjet;
}
