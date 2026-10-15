#include "../utilities.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace {
void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

void check_close(double actual, double expected, double tolerance, const char *message)
{
    require(std::abs(actual - expected) <= tolerance, message);
}
} // namespace

int main() try
{
    // Analytic off-axis case: after the two rotations the b-jet momentum is
    // (0,-4,3). Boosting with the ORIGINAL top beta gives the values below.
    // A rotated-top boost would give a different result and fail this test.
    const TLorentzVector bjet(5, 0, 0, 13);
    const TLorentzVector top(3, 4, 0, 13);
    const TLorentzVector savedTop = top;
    const TLorentzVector savedBjet = bjet;
    const TLorentzVector result = transform_bjet_as_reference(bjet, top);
    const TLorentzVector expected(-341.0 / 100.0, -641.0 / 75.0, 3.0, 185.0 / 12.0);
    for (int component = 0; component < 4; ++component)
    {
        check_close(result[component], expected[component], 1e-12,
                    "Reference transformation differs from the analytic result");
        require(top[component] == savedTop[component], "Top four-vector was changed");
        require(bjet[component] == savedBjet[component], "Input b-jet four-vector was changed");
    }

    std::mt19937 generator(20261014);
    std::uniform_real_distribution<double> momentum(-100.0, 100.0);
    for (int i = 0; i < 1000; ++i)
    {
        TLorentzVector daughter, parent;
        daughter.SetXYZM(momentum(generator), momentum(generator), momentum(generator), 4.8);
        parent.SetXYZM(momentum(generator), momentum(generator), momentum(generator), 173.0);
        const auto transformed = transform_bjet_as_reference(daughter, parent);
        check_close(transformed.M2(), daughter.M2(), 1e-8,
                    "Rotation and boost changed the b-jet invariant mass");
        require(std::isfinite(transformed.E()) && transformed.E() > 0,
                "Transformed energy is invalid");
    }

    for (const auto &invalidTop : {TLorentzVector(0, 0, 0, 0),
                                   TLorentzVector(0, 0, 1, 1),
                                   TLorentzVector(2, 0, 0, 1),
                                   TLorentzVector(0, 0, 0, -1),
                                   TLorentzVector(0, 0, 0, std::numeric_limits<double>::quiet_NaN())})
    {
        bool rejected = false;
        try {
            transform_bjet_as_reference(bjet, invalidTop);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid top four-momentum was accepted");
    }
    std::cout << "PASS: analytic reference sequence, unchanged input vectors, "
                 "1000 invariant-mass checks and invalid-top rejection.\n";
    return 0;
}
catch (const std::exception &error)
{
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
}
