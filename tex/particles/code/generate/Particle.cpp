#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "TMath.h"

// init static members

std::vector<std::unique_ptr<ParticleType>> Particle::m_particle_types{};

// momentum constructors

Momentum::Momentum(double x, double y, double z) : x{x}, y{y}, z{z} {}

Momentum::Momentum(PolarVector const& polar)
    : x{polar.r * std::sin(polar.theta) * std::cos(polar.phi)},
      y{polar.r * std::sin(polar.theta) * std::sin(polar.phi)},
      z{polar.r * std::cos(polar.theta)} {}

// momentum functions

PolarVector Momentum::getPolar() const {
  double r = std::sqrt(x * x + y * y + z * z);
  double theta = std::acos(z / r);
  double phi = std::atan(y / x);

  // convert to the correct phi by adjusting the atan(y / x) result based on
  // x and y coordinates
  if (phi > 0 && x < 0) {
    phi += TMath::Pi();
  } else if (phi < 0 && x > 0) {
    phi += TMath::Pi() * 2.;
  } else if (phi < 0 && x < 0) {
    phi += TMath::Pi();
  }

  return {r, theta, phi};
};

// momentum operators

Momentum Momentum::operator+(Momentum const& momentum) const {
  return {x + momentum.x, y + momentum.y, z + momentum.z};
}

double Momentum::operator*(Momentum const& momentum) const {
  return x * momentum.x + y * momentum.y + z * momentum.z;
}

// constructor

Particle::Particle(std::string const& name, Momentum const& momentum)
    : m_momentum{momentum}, m_index{std::nullopt} {
  if (name != "") {
    setIndex(name);
  }
}

// public methods

void Particle::printData() const {
  if (m_index != std::nullopt) {
    std::cout << "Index: " << m_index.value() << '\n'
              << "Name: " << m_particle_types[m_index.value()]->getName()
              << '\n'
              << "Momentum: (" << m_momentum.x << ", " << m_momentum.y << ", "
              << m_momentum.z << ")\n";
  }
}

int Particle::decayToBody(Particle& dau1, Particle& dau2) const {
  if (getMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if (m_index != std::nullopt) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * std::rand() * invnum - 1.0;
      x2 = 2.0 * std::rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = std::sqrt((-2.0 * std::log(w)) / w);
    y1 = x1 * w;

    massMot += m_particle_types[m_index.value()]->getWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf(
        "Decayment cannot be preformed because mass is too low in this "
        "channel\n");
    return 2;
  }

  double pout =
      std::sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = std::rand() * norm;
  double theta = std::rand() * norm * 0.5 - M_PI / 2.;

  dau1.setMomentum(pout * std::sin(theta) * std::cos(phi),
                   pout * std::sin(theta) * std::sin(phi),
                   pout * std::cos(theta));
  dau2.setMomentum(-pout * std::sin(theta) * std::cos(phi),
                   -pout * std::sin(theta) * std::sin(phi),
                   -pout * std::cos(theta));

  double energy =
      std::sqrt(m_momentum.x * m_momentum.x + m_momentum.y * m_momentum.y +
                m_momentum.z * m_momentum.z + massMot * massMot);

  double bx = m_momentum.x / energy;
  double by = m_momentum.y / energy;
  double bz = m_momentum.z / energy;

  dau1.boost(bx, by, bz);
  dau2.boost(bx, by, bz);

  return 0;
}

// getters

std::optional<int> Particle::getIndex() const { return m_index; }

Momentum Particle::getMomentum() const { return m_momentum; }

double Particle::getEnergy() const {
  if (m_index != std::nullopt) {
    return std::sqrt(std::pow(m_particle_types[m_index.value()]->getMass(), 2) +
                     m_momentum * m_momentum);
  } else {
    std::cout
        << "ERROR: This particle has no energy because its index is invalid!"
        << '\n';

    return 0;
  }
}

double Particle::getMass() const {
  if (m_index != std::nullopt) {
    return m_particle_types[m_index.value()]->getMass();
  } else {
    std::cout
        << "ERROR: This particle has no mass because its index is invalid!"
        << '\n';

    return 0;
  }
}

double Particle::getCharge() const {
  if (m_index != std::nullopt) {
    return m_particle_types[m_index.value()]->getCharge();
  } else {
    std::cout
        << "ERROR: This particle has no charge because its index is invalid!"
        << '\n';

    return 0;
  }
}

std::string Particle::getName() const {
  if (m_index != std::nullopt) {
    return m_particle_types[m_index.value()]->getName();
  } else {
    std::cout
        << "ERROR: This particle has no name because its index is invalid!"
        << '\n';

    return "";
  }
}

double Particle::getInvariantMass(Particle const& p) const {
  auto newMomentum = m_momentum + p.getMomentum();

  return std::sqrt(std::pow(getEnergy() + p.getEnergy(), 2) -
                   newMomentum * newMomentum);
}

// setters

void Particle::setIndex(std::string const& name) {
  m_index = mFindParticleIndex(name);
}

void Particle::setIndex(int index) {
  if (index >= 0 && index < static_cast<int>(m_particle_types.size())) {
    m_index = index;
  } else {
    std::cout << "ERROR: The index \"" << index
              << "\" does not refer to a particle type!" << '\n';
  }
}

void Particle::setMomentum(double px, double py, double pz) {
  m_momentum = {px, py, pz};
}

void Particle::setMomentum(Momentum const& momentum) { m_momentum = momentum; }

// static methods

int Particle::countParticleTypes() { return m_particle_types.size(); }

void Particle::addParticleType(std::string const& name, double mass, int charge,
                               double width) {
  auto existing_index = mFindParticleIndex(name);

  if (existing_index == std::nullopt) {
    if (width == 0) {
      m_particle_types.push_back(
          std::unique_ptr<ParticleType>{new ParticleType{name, mass, charge}});
    } else {
      m_particle_types.push_back(std::unique_ptr<ParticleType>{
          new ResonanceType{name, mass, charge, width}});
    }
  } else {
    std::cout << "ERROR: The \"" << name << "\" particle type already exists!"
              << '\n';
  }
}

void Particle::printParticleTypes() {
  auto v_end = m_particle_types.end();

  for (auto it = m_particle_types.begin(); it < v_end; ++it) {
    (*it)->print();

    if (it != v_end - 1) {
      std::cout << '\n';
    }
  }
}

// private methods

std::optional<int> Particle::mFindParticleIndex(std::string const& name) {
  auto v_begin = m_particle_types.begin();
  auto v_end = m_particle_types.end();

  auto it = std::find_if(v_begin, v_end,
                         [&name](std::unique_ptr<ParticleType> const& pt) {
                           return pt->getName() == name;
                         });

  if (it == v_end) {
    return std::nullopt;
  }

  return std::distance(m_particle_types.begin(), it);
}

void Particle::boost(double bx, double by, double bz) {
  double energy = getEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - b2);
  double bp = bx * m_momentum.x + by * m_momentum.y + bz * m_momentum.z;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  m_momentum.x += gamma2 * bp * bx + gamma * bx * energy;
  m_momentum.y += gamma2 * bp * by + gamma * by * energy;
  m_momentum.z += gamma2 * bp * bz + gamma * bz * energy;
}