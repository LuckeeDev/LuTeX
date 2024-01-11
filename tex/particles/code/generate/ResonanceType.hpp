#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "ParticleType.hpp"

class ResonanceType : public ParticleType {
 public:
  ResonanceType(std::string const&, double, int, double);
  ~ResonanceType() override = default;

  double getWidth() const override;
  void print() const override;

 private:
  double const m_width;
};

#endif