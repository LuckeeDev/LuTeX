#include <string>

#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

class ParticleType {
 public:
  ParticleType(std::string const&, double, int);
  virtual ~ParticleType() = default;

  std::string getName() const;
  double getMass() const;
  int getCharge() const;
  virtual double getWidth() const;
  virtual void print() const;

 private:
  std::string const m_name;
  double const m_mass;
  int const m_charge;
};

#endif