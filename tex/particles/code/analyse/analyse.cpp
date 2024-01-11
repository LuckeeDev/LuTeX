#include <array>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"

const char* PARTICLE_NAMES[7]{"pion+",   "pion-",   "kaon+", "kaon-",
                              "proton+", "proton-", "K*"};

// Labels for fit parameters
const char* HEIGHT_LABEL = "Height (p0)";
const char* MEAN_LABEL = "Mean (p1)";
const char* STDDEV_LABEL = "Std. dev. (p2)";

// Histogram limits for gaussian fits
const double H_LOW = 0.7;
const double H_HIGH = 1.1;

// Define fit functions
Double_t gauss(Double_t* xx, Double_t* par) {
  Double_t x = xx[0];
  Double_t val =
      par[0] * TMath::Exp(-(x - par[1]) * (x - par[1]) / (2 * par[2] * par[2]));
  return val;
}

Double_t exp(Double_t* xx, Double_t* par) {
  Double_t x = xx[0];
  Double_t val = par[0] * TMath::Exp(-x / par[1]);
  return val;
}

Double_t uniform(Double_t* xx, Double_t* par) { return par[0]; }

void analyse(const char* file_name) {
  // Set histogram options. Show entries, parameters, errors and chi square/DOF
  gStyle->SetOptFit(001);
  gStyle->SetOptStat("e");

  TFile* file = new TFile(file_name, "READ");

  std::array<TH1*, 12> histo_array;

  // Read histograms from TList and put them inside an array
  auto key_list = file->GetListOfKeys();

  for (int i{}; i < 12; ++i) {
    auto key = (TKey*)key_list->At(i);
    histo_array[i] = (TH1*)file->Get(key->GetName());
  }

  std::array<double, 12> expected_entries{
      1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 5e8, 2.5e8, 2.5e8, 4.46e7, 4.45e7, 1e5};

  // Print expected and real entries
  std::cout << "ENTRIES" << '\n';

  for (int i{}; i < 12; ++i) {
    auto title = histo_array[i]->GetTitle();
    auto entries = histo_array[i]->GetEntries();

    std::cout << title << ": expected " << expected_entries[i] << ", got "
              << entries << '\n';
  }

  auto particle_types_histogram = histo_array[0];

  // Print generated occurrences
  std::cout << "\nPARTICLE OCCURRENCES" << '\n';

  for (int i{}; i < 7; ++i) {
    auto occurrences = particle_types_histogram->GetBinContent(i + 1);

    std::cout << "Particle " << PARTICLE_NAMES[i] << ": " << occurrences
              << " +- " << std::sqrt(occurrences) << '\n';
  }

  // Create first canvas, with particle types, angles and momentum
  TCanvas* particles_canvas = new TCanvas();
  particles_canvas->Divide(2, 2);

  // Add first histogram with particle type occurrences
  particles_canvas->cd(1);
  auto x_axis = particle_types_histogram->GetXaxis();
  x_axis->CenterLabels();
  for (int i{0}; i < 7; ++i) {
    x_axis->SetBinLabel(i + 1, PARTICLE_NAMES[i]);
  }
  particle_types_histogram->SetXTitle("Particle type");
  particle_types_histogram->SetYTitle("Occurrences");
  particle_types_histogram->Draw();

  // Add and fit azimutal angles histogram
  TF1* azimutal_fit = new TF1("azimutal_fit", uniform, 0., TMath::Pi(), 1);
  azimutal_fit->SetParameter(0, 10000);
  azimutal_fit->SetParName(0, HEIGHT_LABEL);

  particles_canvas->cd(2);
  auto azimutal_angles_histogram = histo_array[1];
  azimutal_angles_histogram->SetXTitle("Azimutal angle (rad)");
  azimutal_angles_histogram->SetYTitle("Occurrences");
  azimutal_angles_histogram->Fit(azimutal_fit, "Q");

  std::cout << "\nAZIMUTAL FIT" << '\n'
            << HEIGHT_LABEL << ": " << azimutal_fit->GetParameter(0) << " +- "
            << azimutal_fit->GetParError(0) << '\n'
            << "Chi square/NDF: "
            << azimutal_fit->GetChisquare() / azimutal_fit->GetNDF() << '\n'
            << "Probability: " << azimutal_fit->GetProb() << '\n';

  // Add and fit polar angles histogram
  TF1* polar_fit = new TF1("polar_fit", uniform, 0., TMath::TwoPi(), 1);
  polar_fit->SetParameter(0, 10000);
  polar_fit->SetParName(0, HEIGHT_LABEL);

  particles_canvas->cd(3);
  auto polar_angles_histogram = histo_array[2];
  polar_angles_histogram->SetXTitle("Polar angle (rad)");
  polar_angles_histogram->SetYTitle("Occurrences");
  polar_angles_histogram->Fit(polar_fit, "Q");

  std::cout << "\nPOLAR FIT" << '\n'
            << HEIGHT_LABEL << ": " << polar_fit->GetParameter(0) << " +- "
            << polar_fit->GetParError(0) << '\n'
            << "Chi square/NDF: "
            << polar_fit->GetChisquare() / polar_fit->GetNDF() << '\n'
            << "Probability: " << polar_fit->GetProb() << '\n';

  // Add and fit momentum histogram
  TF1* momentum_fit = new TF1("momentum_fit", exp, 0., 9., 2);
  momentum_fit->SetParameter(0, 90000);
  momentum_fit->SetParName(0, HEIGHT_LABEL);
  momentum_fit->SetParameter(1, 1.);
  momentum_fit->SetParName(1, MEAN_LABEL);

  particles_canvas->cd(4);
  auto momentum_histogram = histo_array[3];
  momentum_histogram->SetXTitle("Momentum (GeV)");
  momentum_histogram->SetYTitle("Occurrences");
  momentum_histogram->Fit(momentum_fit, "Q");

  std::cout << "\nMOMENTUM FIT" << '\n'
            << HEIGHT_LABEL << ": " << momentum_fit->GetParameter(0) << " +- "
            << momentum_fit->GetParError(0) << '\n'
            << MEAN_LABEL << ": " << momentum_fit->GetParameter(1) << " +- "
            << momentum_fit->GetParError(1) << '\n'
            << "Chi square/NDF: "
            << momentum_fit->GetChisquare() / momentum_fit->GetNDF() << '\n'
            << "Probability: " << momentum_fit->GetProb() << '\n';

  // Create second canvas, with invariant mass
  TCanvas* inv_mass_canvas = new TCanvas();
  inv_mass_canvas->Divide(2, 2);
  inv_mass_canvas->cd(1);

  // Add first histogram directly from the generation
  TF1* k_star_fit = new TF1("k_star_fit", gauss, H_LOW, H_HIGH, 3);
  k_star_fit->SetParameter(0, 500);
  k_star_fit->SetParName(0, HEIGHT_LABEL);
  k_star_fit->SetParameter(1, 0.9);
  k_star_fit->SetParName(1, MEAN_LABEL);
  k_star_fit->SetParameter(2, 0.05);
  k_star_fit->SetParName(2, STDDEV_LABEL);
  auto invm_decayed_h = histo_array[11];
  // Rebin in order to have wider bins
  invm_decayed_h->Rebin(5);
  invm_decayed_h->SetAxisRange(H_LOW, H_HIGH);
  invm_decayed_h->SetTitle("Decay products");
  invm_decayed_h->SetXTitle("Invariant mass (GeV)");
  invm_decayed_h->SetYTitle("Occurrences");
  invm_decayed_h->Fit(k_star_fit, "Q");

  // Create second histogram. Find the signal by subtracting same charge
  // particles from opposite charge particles
  auto invm_opposite_charge_h = histo_array[7];
  invm_opposite_charge_h->Rebin(10);
  auto invm_same_charge_h = histo_array[8];
  invm_same_charge_h->Rebin(10);
  TH1F* invm_subtraction_all = new TH1F(*(TH1F*)invm_opposite_charge_h);
  // Cosmetics
  invm_subtraction_all->SetTitle("Opposite charge - same charge");
  invm_subtraction_all->SetName("invm_subtraction_all");
  invm_subtraction_all->SetXTitle("Invariant mass (GeV)");
  invm_subtraction_all->SetYTitle("Occurrences");
  invm_subtraction_all->Add(invm_opposite_charge_h, invm_same_charge_h, 1, -1);
  invm_subtraction_all->SetEntries(invm_opposite_charge_h->GetEntries());
  invm_subtraction_all->SetAxisRange(H_LOW, H_HIGH);

  TF1* invm_all_fit = new TF1("invm_all_fit", gauss, H_LOW, H_HIGH, 3);
  invm_all_fit->SetParameter(0, 7.996);
  invm_all_fit->SetParName(0, HEIGHT_LABEL);
  invm_all_fit->SetParameter(1, 0.8919);
  invm_all_fit->SetParName(1, MEAN_LABEL);
  invm_all_fit->SetParameter(2, 0.04989);
  invm_all_fit->SetParName(2, STDDEV_LABEL);

  inv_mass_canvas->cd(2);
  invm_subtraction_all->Fit(invm_all_fit, "Q");

  std::cout << "\nINVARIANT MASS BETWEEN ALL PARTICLES (OPPOSITE CHARGE - SAME "
               "CHARGE) FIT"
            << '\n'
            << HEIGHT_LABEL << ": " << invm_all_fit->GetParameter(0) << " +- "
            << invm_all_fit->GetParError(0) << '\n'
            << MEAN_LABEL << ": " << invm_all_fit->GetParameter(1) << " +- "
            << invm_all_fit->GetParError(1) << '\n'
            << STDDEV_LABEL << ": " << invm_all_fit->GetParameter(2) << " +- "
            << invm_all_fit->GetParError(2) << '\n'
            << "Chi square/NDF: "
            << invm_all_fit->GetChisquare() / invm_all_fit->GetNDF() << '\n'
            << "Probability: " << invm_all_fit->GetProb() << '\n'
            << "K* mass: " << invm_all_fit->GetParameter(1) << " +- "
            << invm_all_fit->GetParError(1) << '\n'
            << "K* width: " << invm_all_fit->GetParameter(2) << " +- "
            << invm_all_fit->GetParError(2) << '\n';

  // Create the third histogram. Find the signal by subtracting kaon & pion with
  // same charge from kaon & pion with opposite charge
  auto invm_pion_kaon_opposite_h = histo_array[9];
  invm_pion_kaon_opposite_h->Rebin(10);
  auto invm_pion_kaon_same_h = histo_array[10];
  invm_pion_kaon_same_h->Rebin(10);
  TH1F* invm_subtraction_pion_kaon =
      new TH1F(*(TH1F*)invm_pion_kaon_opposite_h);
  // Cosmetics
  invm_subtraction_pion_kaon->SetTitle(
      "Opposite charge - same charge (kaon & pion)");
  invm_subtraction_pion_kaon->SetName("invm_subtraction_pion_kaon");
  invm_subtraction_pion_kaon->SetXTitle("Invariant mass (GeV)");
  invm_subtraction_pion_kaon->SetYTitle("Occurrences");
  invm_subtraction_pion_kaon->Add(invm_pion_kaon_opposite_h,
                                  invm_pion_kaon_same_h, 1, -1);
  invm_subtraction_pion_kaon->SetEntries(
      invm_pion_kaon_opposite_h->GetEntries());
  invm_subtraction_pion_kaon->SetAxisRange(H_LOW, H_HIGH);

  TF1* invm_pion_kaon_fit =
      new TF1("invm_pion_kaon_fit", gauss, H_LOW, H_HIGH, 3);
  invm_pion_kaon_fit->SetParameter(0, 7.996);
  invm_pion_kaon_fit->SetParName(0, HEIGHT_LABEL);
  invm_pion_kaon_fit->SetParameter(1, 0.8919);
  invm_pion_kaon_fit->SetParName(1, MEAN_LABEL);
  invm_pion_kaon_fit->SetParameter(2, 0.04989);
  invm_pion_kaon_fit->SetParName(2, STDDEV_LABEL);

  inv_mass_canvas->cd(3);
  invm_subtraction_pion_kaon->Fit(invm_pion_kaon_fit, "Q");

  std::cout << "\nINVARIANT MASS BETWEEN KAON AND PION (OPPOSITE CHARGE - SAME "
               "CHARGE) FIT"
            << '\n'
            << HEIGHT_LABEL << ": " << invm_pion_kaon_fit->GetParameter(0)
            << " +- " << invm_pion_kaon_fit->GetParError(0) << '\n'
            << MEAN_LABEL << ": " << invm_pion_kaon_fit->GetParameter(1)
            << " +- " << invm_pion_kaon_fit->GetParError(1) << '\n'
            << STDDEV_LABEL << ": " << invm_pion_kaon_fit->GetParameter(2)
            << " +- " << invm_pion_kaon_fit->GetParError(2) << '\n'
            << "Chi square/NDF: "
            << invm_pion_kaon_fit->GetChisquare() / invm_pion_kaon_fit->GetNDF()
            << '\n'
            << "Probability: " << invm_pion_kaon_fit->GetProb() << '\n'
            << "K* mass: " << invm_pion_kaon_fit->GetParameter(1) << " +- "
            << invm_pion_kaon_fit->GetParError(1) << '\n'
            << "K* width: " << invm_pion_kaon_fit->GetParameter(2) << " +- "
            << invm_pion_kaon_fit->GetParError(2) << '\n';
}