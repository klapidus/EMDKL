//
// Created by Kirill Lapidus on 2019-08-19.
//

#ifndef EMDKL_EMDKLANALYSER_H
#define EMDKL_EMDKLANALYSER_H

#include "stdlib.h"
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <map>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>

struct cjet {
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t weight;
  Int_t   ptbin;
  Int_t   origin_index;
  Bool_t  exclude;
  std::vector <fastjet::PseudoJet> particles;
};

class EMDKLAnalyser {

private:

//  TString fname_exp;
//  Int_t sim_sample_size;
//  Int_t exp_sample_size;

//  vector<cjet> *cjets_exp;
//  vector<cjet> *cjets_sim;

  map< TString, vector<vector<Float_t>> >   emd_values_map;
  map< TString, std::vector<cjet> * >       cjets_map;

//  void get_two_non_overlapping_samples(std::vector<cjet> & nsample1, std::vector<cjet> & nsample2, std::vector<cjet> *insample, Int_t desired_size1, Int_t desired_size2 = 0);
  std::vector<cjet> get_custom_size_sample(std::vector<cjet> *sample_big, Int_t const desired_size);
  void fill_emd_values_from_file(TString const & fname, vector< vector<Float_t> > & emd_values, Int_t const jet_sample_size1, Int_t const jet_sample_size2 = 0);
  static vector<Float_t> min_value_sample1_to_sample2(std::vector<cjet> const & sample1, std::vector<cjet> const & sample2, std::vector< std::vector<Float_t> > const & emd_values, vector<Float_t> & weights);
//  vector<Float_t> max_value_sim_to_sim(std::vector<cjet> const & sample1, std::vector<cjet> const & sample2, std::vector< std::vector<Float_t> > const & all_emd_values_sim, vector<Float_t> & weights);
//  vector<Float_t> ave_value_sim_to_sim(std::vector<cjet> const & sample1, std::vector<cjet> const & sample2, std::vector< std::vector<Float_t> > const & all_emd_values_sim, vector<Float_t> & weights);
//
//  TH1F* get_X_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2 = 0, TString const X = "XYY");
//
//  Float_t min_value_exp_to_sim(std::vector< vector<Float_t> > const & emd_values, std::vector<cjet> const & jets_sim_tailored_sample, Int_t const ix1, Float_t & weight);
//  Float_t max_value_exp_to_sim(std::vector< vector<Float_t> > const & emd_values, std::vector<cjet> const & jets_sim_tailored_sample, Int_t const ix1, Float_t & weight);
//  Float_t ave_value_exp_to_sim(std::vector< vector<Float_t> > const & emd_values, std::vector<cjet> const & jets_sim_tailored_sample, Int_t const ix1, Float_t & weight);
//
//  TH1F* get_X_exp_to_sim(Int_t const sim_trial_sample_size, TString const X, TString hname = "hrndmnm");



public:

//  EMDKLAnalyser() : fname_sim(""), fname_exp(""), fname_exp_to_sim("") {};
  EMDKLAnalyser() {};

  void add_sample(TString const sample_name, vector<cjet> *cjets, TString const emd_filename);
  vector<Float_t> calc_sample_average(TString const sample_name);

  //two subsamples that don't overlap
  TH1F* get_min_subs1_to_subs2(TString const sample_name, vector<cjet> *cjets1, vector<cjet> *cjets2, Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2);

//  TH1F* get_min_sim1_to_sim2(TString const sample_name1, TString const sample_name2, Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2);



//  void fill_all_emd_values();
//  vector<Float_t> calc_sample_average();
//  vector<Float_t> sample_emd_values();


//  void add_sample(&jets);



//  TH1F* get_min_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2 = 0);
//  TH1F* get_max_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2 = 0);
//  TH1F* get_ave_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2 = 0);
//
//  TH1F* get_min_exp_to_sim(Int_t const sim_trial_sample_size, TString hname = "hkl_min_exp_to_sim");
//  TH1F* get_max_exp_to_sim(Int_t const sim_trial_sample_size);
//  TH1F* get_ave_exp_to_sim(Int_t const sim_trial_sample_size); //TODO this one might be further cross-checked
//
//  TH1F* get_min_sim1_to_sim2(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2);
//  TH1F* get_min_sim1_to_sim2(Int_t const ntrials, vector<cjet> *cjets1, vector<cjet> *cjets2, Int_t const trial_sample_size1, Int_t const trial_sample_size2);

private:
  void fill_2Dhisto_from_histos(TH1F** h_1D, Int_t const nhistos, TH2F* h_2D);
  void make_profile_from_2Dhisto(TH2F* h_2D, TH1F* h_profile);
//  void set_proper_histo_uncertanties(TH1F* hweight, TH1F *hnoweight);
//
//  Float_t average_vector_value_weighted(std::vector<Float_t> const &values, std::vector<Float_t> const &weights);

};


#endif //EMDKL_EMDKLANALYSER_H
