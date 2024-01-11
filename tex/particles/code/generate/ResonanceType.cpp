#include "ResonanceType.hpp"

#include <iostream>

// constructor
ResonanceType::ResonanceType(std::string const& name, double mass, int charge,
                             double width)
    : ParticleType{name, mass, charge}, m_width{width} {}

// public methods
double ResonanceType::getWidth() const { return m_width; }

void ResonanceType::print() const {
  std::cout << "Name: " << getName() << '\n'
            << "Mass: " << getMass() << '\n'
            << "Charge: " << getCharge() << '\n'
            << "Width: " << m_width << '\n';
}