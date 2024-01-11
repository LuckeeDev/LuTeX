#include <iostream>
#include <string>
#include <vector>

#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "TBenchmark.h"
#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom.h"

/**
 * Helper function to fill invariant mass histograms with data coming from two
 * particles.
 */
void fillHistograms(Particle const& particle_1, Particle const& particle_2,
                    TH1F* invm_all_h, TH1F* invm_opposite_charge_h,
                    TH1F* invm_same_charge_h, TH1F* invm_pion_kaon_opposite_h,
                    TH1F* invm_pion_kaon_same_h) {
  auto invariant_mass = particle_1.getInvariantMass(particle_2);

  // invariant mass with all particles
  invm_all_h->Fill(invariant_mass);

  // invariant mass with opposite charge particles
  if (particle_2.getCharge() * particle_1.getCharge() < 0) {
    invm_opposite_charge_h->Fill(invariant_mass);
  }

  // invariant mass with same charge particles
  if (particle_2.getCharge() * particle_1.getCharge() > 0) {
    invm_same_charge_h->Fill(invariant_mass);
  }

  // invariant mass with pion+ and kaon- or pion- and kaon+
  if ((particle_1.getName() == "pion+" && particle_2.getName() == "kaon-") ||
      (particle_1.getName() == "kaon-" && particle_2.getName() == "pion+") ||
      (particle_1.getName() == "pion-" && particle_2.getName() == "kaon+") ||
      (particle_1.getName() == "kaon+" && particle_2.getName() == "pion-")) {
    invm_pion_kaon_opposite_h->Fill(invariant_mass);
  }

  // invariant mass with pion+ and kaon+ or piaon- and kaon-
  if ((particle_1.getName() == "pion+" && particle_2.getName() == "kaon+") ||
      (particle_1.getName() == "kaon+" && particle_2.getName() == "pion+") ||
      (particle_1.getName() == "pion-" && particle_2.getName() == "kaon-") ||
      (particle_1.getName() == "kaon-" && particle_2.getName() == "pion-")) {
    invm_pion_kaon_same_h->Fill(invariant_mass);
  }
}

void generate(int n_gen, const char* file_name) {
  gBenchmark->Start("Benchmark");

  R__LOAD_LIBRARY(ParticleType_cpp.so)
  R__LOAD_LIBRARY(ResonanceType_cpp.so)
  R__LOAD_LIBRARY(Particle_cpp.so)

  // mass and width measured in GeV/c^2
  Particle::addParticleType("pion+", 0.13957, 1);
  Particle::addParticleType("pion-", 0.13957, -1);
  Particle::addParticleType("kaon+", 0.49367, 1);
  Particle::addParticleType("kaon-", 0.49367, -1);
  Particle::addParticleType("proton+", 0.93827, 1);
  Particle::addParticleType("proton-", 0.93827, -1);
  Particle::addParticleType("k*", 0.89166, 0, 0.050);

  gRandom->SetSeed();

  TList* histo_list = new TList();

  // particle histograms
  TH1I* particle_types_h =
      new TH1I("particle_types_h", "Particle types", 7, 0, 7);
  histo_list->Add(particle_types_h);  // 0

  TH1F* azimutal_angles_h =
      new TH1F("azimutal_angles_h", "Azimutal angles", 1e3, 0, TMath::Pi());
  histo_list->Add(azimutal_angles_h);  // 1

  TH1F* polar_angles_h =
      new TH1F("polar_angles_h", "Polar angles", 1e3, 0, TMath::Pi() * 2.);
  histo_list->Add(polar_angles_h);  // 2

  TH1F* momentum_h = new TH1F("momentum_h", "Momentum", 1e3, 0, 9);
  histo_list->Add(momentum_h);  // 3

  TH1F* momentum_xy_h = new TH1F("momentum_xy_h", "Momentum xy", 1e3, 0, 9);
  histo_list->Add(momentum_xy_h);  // 4

  TH1F* energy_h = new TH1F("energy_h", "Energy", 1e4, 0, 4);
  histo_list->Add(energy_h);  // 5

  // invariant mass histograms
  TH1F* invm_all_h =
      new TH1F("invm_all_h", "Invariant mass, all particles", 1e4, 0, 9);
  invm_all_h->Sumw2();
  histo_list->Add(invm_all_h);  // 6

  TH1F* invm_opposite_charge_h = new TH1F(
      "invm_opposite_charge_h", "Invariant mass, opposite charge", 1e4, 0, 9);
  invm_opposite_charge_h->Sumw2();
  histo_list->Add(invm_opposite_charge_h);  // 7

  TH1F* invm_same_charge_h =
      new TH1F("invm_same_charge_h", "Invariant mass, same charge", 1e4, 0, 9);
  invm_same_charge_h->Sumw2();
  histo_list->Add(invm_same_charge_h);  // 8

  TH1F* invm_pion_kaon_opposite_h =
      new TH1F("invm_pion_kaon_opposite_h",
               "Invariant mass, pion+ and kaon- or pion- and kaon+", 1e4, 0, 9);
  invm_pion_kaon_opposite_h->Sumw2();
  histo_list->Add(invm_pion_kaon_opposite_h);  // 9

  TH1F* invm_pion_kaon_same_h =
      new TH1F("invm_pion_kaon_same_h",
               "Invariant mass, pion+ and kaon+ or pion- and kaon-", 1e4, 0, 9);
  invm_pion_kaon_same_h->Sumw2();
  histo_list->Add(invm_pion_kaon_same_h);  // 10

  TH1F* invm_decayed_h =
      new TH1F("invm_decayed_h", "Invariant mass, decayed particles from K*",
               1e3, 0.6, 1.2);
  invm_decayed_h->Sumw2();
  histo_list->Add(invm_decayed_h);  // 11

  for (int i{}; i < n_gen; ++i) {
    std::vector<Particle> event_particles(100);

    event_particles.reserve(150);

    for (int j{}; j < 100; ++j) {
      auto r = gRandom->Exp(1);  // GeV
      auto theta = gRandom->Uniform(0, TMath::Pi());
      auto phi = gRandom->Uniform(0, TMath::Pi() * 2.);

      // convert polar to cartesian coordinates
      event_particles[j].setMomentum(Momentum{PolarVector{r, theta, phi}});

      auto x = gRandom->Uniform(0, 1);

      if (x <= 0.4) {
        event_particles[j].setIndex("pion+");
      } else if (x <= 0.8) {
        event_particles[j].setIndex("pion-");
      } else if (x <= 0.85) {
        event_particles[j].setIndex("kaon+");
      } else if (x <= 0.9) {
        event_particles[j].setIndex("kaon-");
      } else if (x <= 0.945) {
        event_particles[j].setIndex("proton+");
      } else if (x <= 0.99) {
        event_particles[j].setIndex("proton-");
      } else {
        event_particles[j].setIndex("k*");

        auto decay_into = gRandom->Uniform(0, 1);

        Particle decay_product_1{};
        Particle decay_product_2{};

        if (decay_into <= 0.5) {
          decay_product_1.setIndex("pion+");
          decay_product_2.setIndex("kaon-");
        } else {
          decay_product_1.setIndex("pion-");
          decay_product_2.setIndex("kaon+");
        }

        event_particles[j].decayToBody(decay_product_1, decay_product_2);

        // fill decay products invariant mass histogram
        auto invariant_mass_products =
            decay_product_1.getInvariantMass(decay_product_2);

        invm_decayed_h->Fill(invariant_mass_products);

        event_particles.push_back(decay_product_1);
        event_particles.push_back(decay_product_2);
      }

      // fill generation histograms
      auto const& new_particle = event_particles[j];

      // type
      particle_types_h->Fill(new_particle.getIndex().value());

      auto momentum = new_particle.getMomentum();
      auto polar_momentum = momentum.getPolar();

      // azimutal angle
      azimutal_angles_h->Fill(polar_momentum.theta);

      // polar angle
      polar_angles_h->Fill(polar_momentum.phi);

      // momentum
      momentum_h->Fill(std::sqrt(momentum * momentum));

      // momentum on xy plane
      momentum_xy_h->Fill(
          std::sqrt(momentum.x * momentum.x + momentum.y * momentum.y));

      // energy
      energy_h->Fill(new_particle.getEnergy());

      // fill invariant mass histograms. This loop improves performance because
      // it avoids unnecessary iterations in the loop after completing the event
      // generation
      if (new_particle.getName() != "k*") {
        for (auto invm_i = j - 1; invm_i >= 0; --invm_i) {
          auto const& invm_particle = event_particles[invm_i];

          if (invm_particle.getName() == "k*") {
            continue;
          }

          fillHistograms(new_particle, invm_particle, invm_all_h,
                         invm_opposite_charge_h, invm_same_charge_h,
                         invm_pion_kaon_opposite_h, invm_pion_kaon_same_h);
        }
      }
    }

    // fill invariant mass histograms with combinations including decayed
    // particles
    auto event_particles_begin = event_particles.begin();
    auto event_particles_end = event_particles.end();

    for (auto it = event_particles_begin + 100; it < event_particles_end;
         ++it) {
      auto const& decayed_particle = *it;

      for (auto invm_it = event_particles_begin; invm_it < it; ++invm_it) {
        auto const& invm_particle = *invm_it;

        if (invm_particle.getName() == "k*") {
          continue;
        }

        fillHistograms(decayed_particle, invm_particle, invm_all_h,
                       invm_opposite_charge_h, invm_same_charge_h,
                       invm_pion_kaon_opposite_h, invm_pion_kaon_same_h);
      }
    }
  }

  TFile* file = new TFile(file_name, "RECREATE");

  histo_list->Write();

  file->Close();

  gBenchmark->Show("Benchmark");
}