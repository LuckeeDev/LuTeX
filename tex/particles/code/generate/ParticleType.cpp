#include "ParticleType.hpp"

#include <iostream>

// constructor

ParticleType::ParticleType(std::string const& name, double mass, int charge)
    : m_name{name}, m_mass{mass}, m_charge{charge} {}

// public methods

std::string ParticleType::getName() const { return m_name; }

double ParticleType::getMass() const { return m_mass; }

int ParticleType::getCharge() const { return m_charge; }

double ParticleType::getWidth() const { return 0.; }

void ParticleType::print() const {
  std::cout << "Name: " << m_name << '\n'
            << "Mass: " << m_mass << '\n'
            << "Charge: " << m_charge << '\n';
}
