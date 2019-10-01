//
// Created by Kirill Lapidus on 2019-08-19.
//

#include "EMDKLAnalyser.h"

void EMDKLAnalyser::add_sample(TString const sample_name, vector<cjet> *cjets, TString const emd_filename) {

  vector<vector<Float_t>> emd_values;
  fill_emd_values_from_file( emd_filename, emd_values, cjets->size() );

  emd_values_map[sample_name] = emd_values;
  cjets_map[sample_name] = cjets;

  cout<<"sample name = "<<sample_name.Data()<<endl;
  cout<<"emd_values_size dim 1 = "<<emd_values.size()<<endl;
  cout<<"emd_values_size dim 2 = "<<emd_values[0].size()<<endl;

}

void EMDKLAnalyser::fill_emd_values_from_file(TString const & fname, vector< vector<Float_t> > & emd_values, Int_t const jet_sample_size1, Int_t const jet_sample_size2) {

  cout<<"Extracting all emd values from file"<<endl;
  cout<<"File name: "<<fname.Data()<<endl;

  //same sample (square shape ndim2=ndim1)
  //or two different samples (ndim2!=ndim1)
  Int_t second_dimension_size = jet_sample_size1; //default behaviour, one sample
  if (jet_sample_size2 > 0) second_dimension_size = jet_sample_size2; //rectangular shape, two samples

  FILE *pFile = fopen(fname.Data(), "r");
  for (int i = 0; i < jet_sample_size1; i++) {
//    cout<<"i = "<<i<<endl;
    Float_t tmp_value;
    std::vector<Float_t> emd_values_line;
    for (int j = 0; j < second_dimension_size; j++) {
      if ( j!=(second_dimension_size-1) ) fscanf(pFile, "%f", &tmp_value);
      else fscanf(pFile, "%f\n", &tmp_value);
      emd_values_line.push_back(tmp_value);
    }
    emd_values.push_back(emd_values_line);
    emd_values_line.clear();
  }
  fclose(pFile);

}

vector<Float_t> EMDKLAnalyser::calc_sample_average(TString const sample_name) {

  auto & emd_values = emd_values_map[sample_name];
  auto dim1_size    = emd_values.size(); //rectangular for the moment

  vector<Float_t> average_values;

  for (int i = 0; i < dim1_size; i++) {
    float sum = 0.0;
    int counter = 0;
    for (int j = 0; j < dim1_size; j++) {
      if (i==j) continue;
      sum += emd_values[i][j];
      counter++;
    }
    average_values.push_back( sum/(1.0*counter) );
  }

  return average_values;

}

/*
TH1F* EMDKLAnalyser::get_min_sim1_to_sim2(TString const sample_name1, TString const sample_name2, Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2) {

  Int_t nbins = 0;
  Int_t valuel = 0;
  Int_t valuer = 0;

  nbins = 15;
  valuel = 0;
  valuer = 15;

  TH1F *hemd_X_sim1_to_sim2[ntrials];

  for (int tr = 0; tr < ntrials; tr++) {

    TString hname_X = "hemd_X_sim1_to_sim2_";
    hname_X += tr;
    hemd_X_sim1_to_sim2[tr] = new TH1F(hname_X,hname_X,nbins,valuel,valuer);
    hemd_X_sim1_to_sim2[tr]->Sumw2(kFALSE);

    std::vector<cjet> jets_sim1_trial_sample = get_custom_size_sample(cjets_map[sample_name1], trial_sample_size1);
    std::vector<cjet> jets_sim2_trial_sample = get_custom_size_sample(cjets_map[sample_name2], trial_sample_size2);

    vector<Float_t> weights_emd_X_sim1_to_sim2;
    vector<Float_t> emd_X_sim1_to_sim2;
    //TO DO:
    if ( sample_name1.Data()!=sample_name2.Data() ) cout<<"sim1_to_sim2 not implemented yet for two different samples!"<<endl;
    emd_X_sim1_to_sim2 = min_value_sim_to_sim(jets_sim1_trial_sample, jets_sim2_trial_sample, emd_values_sim, weights_emd_X_sim1_to_sim2);
    for (Int_t j = 0; j < jets_sim1_trial_sample.size(); j++) {
      hemd_X_sim1_to_sim2[tr]->Fill( emd_X_sim1_to_sim2[j], weights_emd_X_sim1_to_sim2[j] ); //this is in the end sum of weigths for each emd bin
    }

  }//ntrials

  //REALLY NEED TO WORK OUT WEIGHTS HERE!

  //normalize the sim-sim distance, needed for the 2D profile
  //divide by the sum of weights???
  for (int tr = 0; tr < ntrials; tr++) {
    auto tnorm_ave = hemd_X_sim1_to_sim2[tr]->Integral();
    hemd_X_sim1_to_sim2[tr]->Scale(1./tnorm_ave);
  }

  auto tmp_ind1 = ( cjets1->at(0) ).origin_index;
  auto tmp_ind2 = ( cjets2->at(0) ).origin_index;

  TString histo_name = "hemdkl_";
  histo_name += "min_";
  histo_name += tmp_ind1;
  histo_name += "_";
  histo_name += tmp_ind2;
  histo_name += "_sim1_to_sim2_2D_profile";

  //2D profile of EMD min sim to sim

  TString tmp_name = "hemd_X_sim1_to_sim2_2D";
  TH2F *hemd_X_sim1_to_sim2_2D = new TH2F(tmp_name,tmp_name,nbins,valuel,valuer,1000,0,1.0);
  TH1F *hemd_X_sim1_to_sim2_2D_profile = new TH1F(histo_name,histo_name,nbins,valuel,valuer);

  fill_2Dhisto_from_histos(hemd_X_sim1_to_sim2, ntrials, hemd_X_sim1_to_sim2_2D);
  make_profile_from_2Dhisto(hemd_X_sim1_to_sim2_2D, hemd_X_sim1_to_sim2_2D_profile);

  for (int tr = 0; tr < ntrials; tr++) {
    delete hemd_X_sim1_to_sim2[tr];
  }
  delete hemd_X_sim1_to_sim2_2D;

  return hemd_X_sim1_to_sim2_2D_profile;
}
 */


TH1F* EMDKLAnalyser::get_min_subs1_to_subs2(TString const sample_name, vector<cjet> *cjets1, vector<cjet> *cjets2, Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2) {

  if ( emd_values_map.find(sample_name) == emd_values_map.end() ) {
    cout<<"error: emd values entry in the map doesn't exist!"<<endl;
  }

  Int_t nbins  = 15;
  Int_t valuel = 0;
  Int_t valuer = 30;

  TH1F *hemd_X_sim1_to_sim2[ntrials];

  for (int tr = 0; tr < ntrials; tr++) {
    TString hname_X = "hemd_X_sim_to_sim_";
    hname_X += tr;
    hemd_X_sim1_to_sim2[tr] = new TH1F(hname_X,hname_X,nbins,valuel,valuer);
    hemd_X_sim1_to_sim2[tr]->Sumw2(kFALSE);
    std::vector<cjet> jets_sim_trial_sample1 = get_custom_size_sample(cjets1, trial_sample_size1);
    std::vector<cjet> jets_sim_trial_sample2 = get_custom_size_sample(cjets2, trial_sample_size2);
    vector<Float_t> weights_emd_X_sim1_to_sim2;
    vector<Float_t> emd_X_sim1_to_sim2;
    //TO DO:
    emd_X_sim1_to_sim2 = min_value_sample1_to_sample2(jets_sim_trial_sample1, jets_sim_trial_sample2, emd_values_map[sample_name], weights_emd_X_sim1_to_sim2);
    for (Int_t j = 0; j < jets_sim_trial_sample1.size(); j++) {
      hemd_X_sim1_to_sim2[tr]->Fill( emd_X_sim1_to_sim2[j], weights_emd_X_sim1_to_sim2[j] ); //this is in the end sum of weigths for each emd bin
    }

  }//ntrials

  //REALLY NEED TO WORK OUT WEIGHTS HERE!

  //normalize the sim-sim distance, needed for the 2D profile
  //divide by the sum of weights???
  for (int tr = 0; tr < ntrials; tr++) {
    auto tnorm_ave = hemd_X_sim1_to_sim2[tr]->Integral();
    hemd_X_sim1_to_sim2[tr]->Scale(1./tnorm_ave);
  }

  auto tmp_ind1 = ( cjets1->at(0) ).origin_index;
  auto tmp_ind2 = ( cjets2->at(0) ).origin_index;

  TString histo_name = "hemdkl_";
  histo_name += sample_name.Data();
  histo_name += "_min_";
  histo_name += tmp_ind1;
  histo_name += "_";
  histo_name += tmp_ind2;
  histo_name += "_sim1_to_sim2_2D_profile";

  //2D profile of EMD min sim to sim

  TString tmp_name = "hemd_X_sim1_to_sim2_2D";
  TH2F *hemd_X_sim1_to_sim2_2D = new TH2F(tmp_name,tmp_name,nbins,valuel,valuer,1000,0,1.0);
  TH1F *hemd_X_sim1_to_sim2_2D_profile = new TH1F(histo_name,histo_name,nbins,valuel,valuer);

  fill_2Dhisto_from_histos(hemd_X_sim1_to_sim2, ntrials, hemd_X_sim1_to_sim2_2D);
  make_profile_from_2Dhisto(hemd_X_sim1_to_sim2_2D, hemd_X_sim1_to_sim2_2D_profile);

  for (int tr = 0; tr < ntrials; tr++) {
    delete hemd_X_sim1_to_sim2[tr];
  }
  delete hemd_X_sim1_to_sim2_2D;

  return hemd_X_sim1_to_sim2_2D_profile;
}


std::vector<cjet> EMDKLAnalyser::get_custom_size_sample(std::vector<cjet> *sample_big, Int_t const desired_size) {

  std::vector<cjet> nsample;
  auto const initial_size = sample_big->size();

  if (desired_size > initial_size) {
    cout << "desired sample size is larger than the pool size!" << endl;
    return nsample;
  }
  else if (desired_size==initial_size) {
    cout<<"desired sample size is the same as the pool size, no need to randomize!"<<endl;
    for (auto j : *sample_big) {
      nsample.push_back(j);
    }
    return nsample;
  }

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> dis(1, initial_size);

  std::vector<int> random_values;
  for (int i = 0; i < desired_size; i++) {
    auto rndvalue = dis(gen);
    rndvalue -= 1; //from 0 to size - 1
    auto iter = std::find(random_values.begin(), random_values.end(), rndvalue);
    if (iter == random_values.end()) {
      random_values.push_back(rndvalue);
    }
    else {
      i--;
    }
  }

  for (auto v : random_values) {
    nsample.push_back( sample_big->at(v) );
  }

//  cout<<"reduced sim sample size = "<<nsample.size()<<endl;

  return nsample;

}






/*
void EMDKLAnalyser::set_sim_emd_fname(TString const fname) {
  fname_sim = fname;
}

void EMDKLAnalyser::set_exp_emd_fname(TString const fname) {
  fname_exp = fname;
}

void EMDKLAnalyser::set_exp_to_sim_emd_fname(TString const fname) {
  fname_exp_to_sim = fname;
}

void EMDKLAnalyser::fill_all_emd_values() {
  if ( !fname_sim.IsNull() ) fill_emd_values_from_file(fname_sim, emd_values_sim, sim_sample_size);
  if ( !fname_sim_pl.IsNull() ) fill_emd_values_from_file(fname_sim_pl, emd_values_sim_pl, sim_sample_size);
  if ( !fname_sim_dl.IsNull() ) fill_emd_values_from_file(fname_sim_dl, emd_values_sim_dl, sim_sample_size);
  if ( !fname_exp.IsNull() ) fill_emd_values_from_file(fname_exp, emd_values_exp, exp_sample_size);
  if ( !fname_exp_to_sim.IsNull() ) fill_emd_values_from_file(fname_exp_to_sim, emd_values_exp_to_sim, exp_sample_size, sim_sample_size);
}

void EMDKLAnalyser::set_sim_sample_size(Int_t const size) {
  sim_sample_size = size;
}

void EMDKLAnalyser::set_exp_sample_size(Int_t const size) {
  exp_sample_size = size;
}

void EMDKLAnalyser::set_exp_cjets(vector<cjet> *cjets) {
  cjets_exp = cjets;
}

void EMDKLAnalyser::set_sim_cjets(vector<cjet> *cjets, int const sample) {
  if (0 == sample)      cjets_sim = cjets;
  else if (1 == sample) cjets_sim1 = cjets;
  else if (2 == sample) cjets_sim2 = cjets;
}

TH1F* EMDKLAnalyser::get_max_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2) {
  return get_X_sim_to_sim(ntrials, trial_sample_size1, trial_sample_size2, "max");
}

TH1F* EMDKLAnalyser::get_min_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2) {
  return get_X_sim_to_sim(ntrials, trial_sample_size1, trial_sample_size2, "min");
}

TH1F* EMDKLAnalyser::get_ave_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2) {
  return get_X_sim_to_sim(ntrials, trial_sample_size1, trial_sample_size2, "ave");
}

TH1F* EMDKLAnalyser::get_X_sim_to_sim(Int_t const ntrials, Int_t const trial_sample_size1, Int_t const trial_sample_size2, TString const X) {

  Int_t nbins = 0;
  Int_t valuel = 0;
  Int_t valuer = 0;

  if  ( 0==strncmp(X.Data(),"min",3) ) {
    nbins = 15;
    valuel = 0;
    valuer = 15;
  }
  else if ( 0==strncmp(X.Data(),"max",3) ) {
    nbins = 40;
    valuel = 240;
    valuer = 400;
  }
  else if ( 0==strncmp(X.Data(),"ave",3) ) {
    nbins = 40;
    valuel = 0;
    valuer = 40;
  }
  else cout<<"Error, unknown quantity asked! "<<X.Data()<<endl;


  TH1F *hemd_X_sim_to_sim[ntrials];

  for (int tr = 0; tr < ntrials; tr++) {
    TString hname_X = "hemd_X_sim_to_sim_";
    hname_X += tr;
    hemd_X_sim_to_sim[tr] = new TH1F(hname_X,hname_X,nbins,valuel,valuer);
    hemd_X_sim_to_sim[tr]->Sumw2(kFALSE);
    std::vector<cjet> jets_sim_trial_sample1;
    std::vector<cjet> jets_sim_trial_sample2;
    get_two_non_overlapping_samples(jets_sim_trial_sample1, jets_sim_trial_sample2, cjets_sim, trial_sample_size1, trial_sample_size2);
    vector<Float_t> weights_emd_X_sim_to_sim;
    vector<Float_t> emd_X_sim_to_sim;
    if ( 0==strncmp(X.Data(),"min",3) ) {
      emd_X_sim_to_sim = min_value_sim_to_sim(jets_sim_trial_sample1, jets_sim_trial_sample2, emd_values_sim, weights_emd_X_sim_to_sim);
    }
    else if ( 0==strncmp(X.Data(),"max",3) ) {
      emd_X_sim_to_sim = max_value_sim_to_sim(jets_sim_trial_sample1, jets_sim_trial_sample2, emd_values_sim, weights_emd_X_sim_to_sim);
    }
    else if ( 0==strncmp(X.Data(),"ave",3) ) {
      emd_X_sim_to_sim = ave_value_sim_to_sim(jets_sim_trial_sample1, jets_sim_trial_sample2, emd_values_sim, weights_emd_X_sim_to_sim);
    }
    else cout<<"value "<<X.Data()<<" not implemented!"<<endl;
    for (Int_t j = 0; j < jets_sim_trial_sample1.size(); j++) {
      hemd_X_sim_to_sim[tr]->Fill( emd_X_sim_to_sim[j], weights_emd_X_sim_to_sim[j] );
    }
  }

  //normalize the sim-sim distance, needed for the 2D profile
  //divide by the sum of weights???
  for (int tr = 0; tr < ntrials; tr++) {
    auto tnorm_ave = hemd_X_sim_to_sim[tr]->Integral();
    hemd_X_sim_to_sim[tr]->Scale(1./tnorm_ave);
  }

  TString histo_name = "hemdkl_";
  histo_name += X;
  histo_name += "_sim_to_sim_2D_profile";

  //2D profile of EMD min sim to sim
  TH2F *hemd_X_sim_to_sim_2D = new TH2F("hemd_X_sim_to_sim_2D","hemd_X_sim_to_sim_2D",nbins,valuel,valuer,1000,0,1.0);
  TH1F *hemd_X_sim_to_sim_2D_profile = new TH1F(histo_name,histo_name,nbins,valuel,valuer);

  fill_2Dhisto_from_histos(hemd_X_sim_to_sim, ntrials, hemd_X_sim_to_sim_2D);
  make_profile_from_2Dhisto(hemd_X_sim_to_sim_2D, hemd_X_sim_to_sim_2D_profile);

  for (int tr = 0; tr < ntrials; tr++) {
    delete hemd_X_sim_to_sim[tr];
  }
  delete hemd_X_sim_to_sim_2D;

  return hemd_X_sim_to_sim_2D_profile;

}


TH1F* EMDKLAnalyser::get_min_sim1_to_sim2(Int_t const ntrials, vector<cjet> *cjets1, vector<cjet> *cjets2, Int_t const trial_sample_size1, Int_t const trial_sample_size2) {

  Int_t nbins = 0;
  Int_t valuel = 0;
  Int_t valuer = 0;

  nbins = 15;
  valuel = 0;
  valuer = 15;

  TH1F *hemd_X_sim1_to_sim2[ntrials];
//  TH1F *hemd_X_sim1_to_sim2_weights[ntrials]; //for each trial, store the weights for 1D histo: one value - one weight

  for (int tr = 0; tr < ntrials; tr++) {

    TString hname_X = "hemd_X_sim1_to_sim2_";
    hname_X += tr;
    hemd_X_sim1_to_sim2[tr] = new TH1F(hname_X,hname_X,nbins,valuel,valuer);
    hemd_X_sim1_to_sim2[tr]->Sumw2(kFALSE);

//    hname_X += "_"; hname_X += "weigths";
//    hemd_X_sim1_to_sim2_weights[tr] = new TH1F(hname_X,hname_X,nbins,valuel,valuer);
//    hemd_X_sim1_to_sim2_weights[tr]->Sumw2(kFALSE);

    std::vector<cjet> jets_sim1_trial_sample = get_custom_size_sample(cjets1, trial_sample_size1);
    std::vector<cjet> jets_sim2_trial_sample = get_custom_size_sample(cjets2, trial_sample_size2);

    vector<Float_t> weights_emd_X_sim1_to_sim2;
    vector<Float_t> emd_X_sim1_to_sim2;
    //TO DO:
    emd_X_sim1_to_sim2 = min_value_sim_to_sim(jets_sim1_trial_sample, jets_sim2_trial_sample, emd_values_sim, weights_emd_X_sim1_to_sim2);
    for (Int_t j = 0; j < jets_sim1_trial_sample.size(); j++) {
      hemd_X_sim1_to_sim2[tr]->Fill( emd_X_sim1_to_sim2[j], weights_emd_X_sim1_to_sim2[j] ); //this is in the end sum of weigths for each emd bin
    }

  }//ntrials

  //REALLY NEED TO WORK OUT WEIGHTS HERE!

  //normalize the sim-sim distance, needed for the 2D profile
  //divide by the sum of weights???
  for (int tr = 0; tr < ntrials; tr++) {
    auto tnorm_ave = hemd_X_sim1_to_sim2[tr]->Integral();
    hemd_X_sim1_to_sim2[tr]->Scale(1./tnorm_ave);
  }

  auto tmp_ind1 = ( cjets1->at(0) ).origin_index;
  auto tmp_ind2 = ( cjets2->at(0) ).origin_index;

  TString histo_name = "hemdkl_";
  histo_name += "min_";
  histo_name += tmp_ind1;
  histo_name += "_";
  histo_name += tmp_ind2;
  histo_name += "_sim1_to_sim2_2D_profile";

  //2D profile of EMD min sim to sim

  TString tmp_name = "hemd_X_sim1_to_sim2_2D";
  TH2F *hemd_X_sim1_to_sim2_2D = new TH2F(tmp_name,tmp_name,nbins,valuel,valuer,1000,0,1.0);
  TH1F *hemd_X_sim1_to_sim2_2D_profile = new TH1F(histo_name,histo_name,nbins,valuel,valuer);

  fill_2Dhisto_from_histos(hemd_X_sim1_to_sim2, ntrials, hemd_X_sim1_to_sim2_2D);
  make_profile_from_2Dhisto(hemd_X_sim1_to_sim2_2D, hemd_X_sim1_to_sim2_2D_profile);

  for (int tr = 0; tr < ntrials; tr++) {
    delete hemd_X_sim1_to_sim2[tr];
  }
  delete hemd_X_sim1_to_sim2_2D;

  return hemd_X_sim1_to_sim2_2D_profile;
}

void EMDKLAnalyser::get_two_non_overlapping_samples(std::vector<cjet> & nsample1, std::vector<cjet> & nsample2, std::vector<cjet> *insample, Int_t desired_size1, Int_t desired_size2) {

  auto initial_size = insample->size();

  if (0==desired_size2) desired_size2 = desired_size1;

  if ( (desired_size1 + desired_size2) > initial_size) {
    cout<<"Warning! Error!"<<endl;
    cout<<"get_two_non_overlapping_samples"<<endl;
    cout<<"desired size1 = "<<desired_size1<<endl;
    cout<<"desired size2 = "<<desired_size2<<endl;
    cout<<"initial size = "<<initial_size<<endl;
  }

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> dis(1, initial_size);

  std::vector<int> random_values1;
  for (int i = 0; i < desired_size1; i++) {
    auto rndvalue = dis(gen) - 1; //from 0 to size - 1
    auto iter = std::find(random_values1.begin(), random_values1.end(), rndvalue);
    if (iter == random_values1.end()) { //only unique indices
      random_values1.push_back(rndvalue);
    }
    else {
      i--;
    }
  }

  std::vector<int> random_values2;
  for (int i = 0; i < desired_size2; i++) {
    auto rndvalue = dis(gen) - 1; //from 0 to size - 1
    //check both arrays:
    auto iter1 = std::find(random_values1.begin(), random_values1.end(), rndvalue);
    auto iter2 = std::find(random_values2.begin(), random_values2.end(), rndvalue);
    if ( ( iter1 == random_values1.end() ) && ( iter2 == random_values2.end() ) ) {
      random_values2.push_back(rndvalue);
    }
    else {
      i--;
    }
  }

  for (int i = 0; i < desired_size1; i++) {
//    nsample1.push_back( insample[ random_values1[i] ] );
    nsample1.push_back( insample->at( random_values1[i] ) );
  }
  for (int i = 0; i < desired_size2; i++) {
//    nsample2.push_back( insample[ random_values2[i] ] );
    nsample2.push_back( insample->at( random_values2[i] ) );

  }

  for (auto v : random_values1) {
    auto iter = std::find(random_values2.begin(), random_values2.end(), v);
    if (iter != random_values2.end()) cout<<"error: index doubling between the samples!"<<endl;
  }


}

std::vector<cjet> EMDKLAnalyser::get_custom_size_sample(std::vector<cjet> *sample_big, Int_t const desired_size) {

  std::vector<cjet> nsample;
  auto const initial_size = sample_big->size();

  if (desired_size > initial_size) {
    cout << "desired sample size is larger than the pool size!" << endl;
    return nsample;
  }
  else if (desired_size==initial_size) {
    cout<<"desired sample size is the same as the pool size, no need to randomize!"<<endl;
    for (auto j : *sample_big) {
      nsample.push_back(j);
    }
    return nsample;
  }

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> dis(1, initial_size);

  std::vector<int> random_values;
  for (int i = 0; i < desired_size; i++) {
    auto rndvalue = dis(gen);
    rndvalue -= 1; //from 0 to size - 1
    auto iter = std::find(random_values.begin(), random_values.end(), rndvalue);
    if (iter == random_values.end()) {
      random_values.push_back(rndvalue);
    }
    else {
      i--;
    }
  }

//  cout<<"size of rndvalues vector = "<<random_values.size()<<endl;
  for (auto v : random_values) {
//    nsample.push_back( sample_big[v] );
    nsample.push_back( sample_big->at(v) );
  }

  cout<<"reduced sim sample size = "<<nsample.size()<<endl;

  return nsample;

}

*/


void EMDKLAnalyser::make_profile_from_2Dhisto(TH2F* h_2D, TH1F* h_profile) {
  Int_t nbinsx = h_2D->GetXaxis()->GetNbins();
  for (Int_t ib = 1; ib <= nbinsx; ib++) {
    auto hproj = h_2D->ProjectionY("tmp_name", ib, ib);
    auto mean = hproj->GetMean();
    auto sigma = hproj->GetStdDev();
    h_profile->SetBinContent(ib, mean);
    h_profile->SetBinError(ib, sigma);
  }
}

void EMDKLAnalyser::fill_2Dhisto_from_histos(TH1F** h_1D, Int_t const nhistos, TH2F* h_2D) {

  for (int tr = 0; tr < nhistos; tr++) {
    Int_t nbinsx = h_1D[tr]->GetXaxis()->GetNbins();
    for (Int_t ib = 1; ib <= nbinsx; ib++) {
      auto bcenter = h_1D[tr]->GetXaxis()->GetBinCenter(ib);
      auto bvalue  = h_1D[tr]->GetBinContent(ib);
      h_2D->Fill(bcenter, bvalue);
    }
  }

}


vector<Float_t> EMDKLAnalyser::min_value_sample1_to_sample2(std::vector<cjet> const & sample1, std::vector<cjet> const & sample2, std::vector< std::vector<Float_t> > const & emd_values, vector<Float_t> & weights) {
  vector<Float_t> min_values;
  Float_t weight = 0.0;
  for (auto j1 : sample1) {
    if (j1.exclude) continue; //skip the jets that are marked for exclusion
    Float_t min_value = 1.e+10;
    for (auto j2 : sample2) {
      if ( j2.exclude ) continue; //skip the jets that are marked for exclusion
      if ( (j1.origin_index < 0) || (j2.origin_index < 0) ) cout<<"origin index error"<<endl;
      if (j1.origin_index == j2.origin_index) continue; //very same jet
      auto emd_value = emd_values[j1.origin_index][j2.origin_index];
      if (emd_value < min_value) {
        min_value = emd_value;
        weight = (j1.weight) * (j2.weight); //?
      }
    }
    min_values.push_back(min_value);
    weights.push_back(weight);
  }
  return min_values;
}

/*

vector<Float_t> EMDKLAnalyser::max_value_sim_to_sim(std::vector<cjet> const & sample1, std::vector<cjet> const & sample2, std::vector< std::vector<Float_t> > const & all_emd_values_sim, vector<Float_t> & weights) {
  vector<Float_t> max_values;
  Float_t weight = 0.0;
  for (auto j1 : sample1) {
    if (j1.exclude) continue; //skip the jets that are marked for exclusion
    Float_t max_value = -1.e+10;
    for (auto j2 : sample2) {
      if ( j2.exclude ) continue; //skip the jets that are marked for exclusion
      if ( (j1.origin_index < 0) || (j2.origin_index < 0) ) cout<<"origin index error"<<endl;
      if (j1.origin_index == j2.origin_index) continue; //very same jet
      auto emd_value = all_emd_values_sim[j1.origin_index][j2.origin_index];
      if (emd_value > max_value) {
        max_value = emd_value;
        weight = (j1.weight) * (j2.weight);
      }
    }
    max_values.push_back(max_value);
    weights.push_back(weight);
  }
  return max_values;
}

vector<Float_t> EMDKLAnalyser::ave_value_sim_to_sim(std::vector<cjet> const & sample1, std::vector<cjet> const & sample2, std::vector< std::vector<Float_t> > const & all_emd_values_sim, vector<Float_t> & weights) {
  vector<Float_t> ave_values;
  for (auto j1 : sample1) {
    if (j1.exclude) continue; //skip the jets that are marked for exclusion
    vector<Float_t> interim_values;
    vector<Float_t> interim_weights;
    for (auto j2 : sample2) {
      if ( j2.exclude ) continue; //skip the jets that are marked for exclusion
      if ( (j1.origin_index < 0) || (j2.origin_index < 0) ) cout<<"origin index error"<<endl;
      if (j1.origin_index == j2.origin_index) continue; //very same jet
      auto emd_value = all_emd_values_sim[j1.origin_index][j2.origin_index];
      interim_values.push_back( emd_value );
      interim_weights.push_back( j2.weight );
    }
    auto tmp_average_value = average_vector_value_weighted(interim_values, interim_weights);
    ave_values.push_back(tmp_average_value);
    weights.push_back( j1.weight );
    interim_values.clear();
    interim_weights.clear();
  }
  return ave_values;
}

Float_t EMDKLAnalyser::min_value_exp_to_sim(std::vector< vector<Float_t> > const & emd_values, std::vector<cjet> const & jets_sim_tailored_sample, Int_t const ix1, Float_t & weight) {
  Float_t min_value = 1.e+10;
  for (auto j : jets_sim_tailored_sample) {
    if (j.origin_index < 0) cout<<"origin index error"<<endl;
    if ( j.exclude ) continue; //skip the jets that are marked for exclusion
    auto emd_value = emd_values[ix1][j.origin_index]; //exp-sim
    if (emd_value < min_value) {
      min_value = emd_value;
      weight = j.weight;
    }
  }
//  cout<<"min_value = "<<min_value<<endl;
  return min_value;
}

Float_t EMDKLAnalyser::max_value_exp_to_sim(std::vector< vector<Float_t> > const & emd_values, std::vector<cjet> const & jets_sim_tailored_sample, Int_t const ix1, Float_t & weight) {
  Float_t max_value = -1.e+5;
  for (auto j : jets_sim_tailored_sample) {
    if (j.origin_index < 0) cout<<"origin index error"<<endl;
    if ( j.exclude ) continue; //skip the jets that are marked for exclusion
    auto emd_value = emd_values[ix1][j.origin_index]; //exp-sim
    if (emd_value > max_value) {
      max_value = emd_value;
      weight = j.weight;
    }
  }
  return max_value;
}

Float_t EMDKLAnalyser::ave_value_exp_to_sim(std::vector< vector<Float_t> > const & emd_values, std::vector<cjet> const & jets_sim_tailored_sample, Int_t const ix1, Float_t & weight) {
  weight = 1.0; //TODO is this correct logic? should think about this
  vector<Float_t> interim_values;
  vector<Float_t> interim_weights;
  for (auto j : jets_sim_tailored_sample) {
    if (j.origin_index < 0) cout<<"origin index error"<<endl;
    if ( j.exclude ) continue; //skip the jets that are marked for exclusion
    auto emd_value = emd_values[ix1][j.origin_index]; //exp-sim
    interim_values.push_back(emd_value);
    interim_weights.push_back(j.weight);
  }
  return average_vector_value_weighted(interim_values, interim_weights);
}


Float_t EMDKLAnalyser::average_vector_value_weighted(std::vector<Float_t> const &values, std::vector<Float_t> const &weights) {
  Float_t sum   = 0.0;
  Float_t sum_w = 0.0;
  if ( values.size() != weights.size() ) cout<<"values-weights length mismatch!"<<endl;
  for (int i = 0; i < values.size(); i++) {
    sum += values[i] * weights[i];
    sum_w += weights[i];
  }
  return sum/sum_w;
}

TH1F* EMDKLAnalyser::get_X_exp_to_sim(Int_t const sim_trial_sample_size, TString const X, TString hname) {

  std::vector<cjet> jets_sim_tailored_sample = get_custom_size_sample( cjets_sim, sim_trial_sample_size ); //same size as exp
  Float_t emd_X_exp_to_sim[exp_sample_size]; //desired (min/max/ave/...) value for each exp jet
  Float_t weight_emd_X_exp_to_sim[exp_sample_size]; //desire (min/max/ave/...) value for each exp jet
  for (int i = 0; i < exp_sample_size; i++) {
    if      ( 0==strncmp(X.Data(),"min",3) ) emd_X_exp_to_sim[i] = min_value_exp_to_sim( emd_values_exp_to_sim, jets_sim_tailored_sample, i, weight_emd_X_exp_to_sim[i] );
    else if ( 0==strncmp(X.Data(),"max",3) ) emd_X_exp_to_sim[i] = max_value_exp_to_sim( emd_values_exp_to_sim, jets_sim_tailored_sample, i, weight_emd_X_exp_to_sim[i] );
    else if ( 0==strncmp(X.Data(),"ave",3) ) emd_X_exp_to_sim[i] = ave_value_exp_to_sim( emd_values_exp_to_sim, jets_sim_tailored_sample, i, weight_emd_X_exp_to_sim[i] );
    else cout<<"Error, unknown quantity asked! "<<X.Data()<<endl;
  }

  TString histo_name = "hemdkl_";
  histo_name += X;
  histo_name += "_exp_to_sim";

  TString histo_name_weight = histo_name;
  histo_name_weight += "_weight";

  Int_t nbins = 0;
  Int_t valuel = 0;
  Int_t valuer = 0;

  if  ( 0==strncmp(X.Data(),"min",3) ) {
    nbins = 15;
    valuel = 0;
    valuer = 15;
  }
  else if ( 0==strncmp(X.Data(),"max",3) ) {
    nbins = 40;
    valuel = 240;
    valuer = 400;
  }
  else if ( 0==strncmp(X.Data(),"ave",3) ) {
    nbins = 40;
    valuel = 0;
    valuer = 40;
  }
  else cout<<"Error, unknown quantity asked! "<<X.Data()<<endl;

  TH1F *hemd_X_exp_to_sim_weight = new TH1F(hname,hname,nbins,valuel,valuer);
  TH1F *hemd_X_exp_to_sim_noweight = new TH1F(histo_name_weight,histo_name_weight,nbins,valuel,valuer);

  for (int i = 0; i < exp_sample_size; i++) {
    hemd_X_exp_to_sim_noweight->Fill( emd_X_exp_to_sim[i] ); //not weighted
    hemd_X_exp_to_sim_weight->Fill( emd_X_exp_to_sim[i], weight_emd_X_exp_to_sim[i] ); //weighted
  }

  set_proper_histo_uncertanties(hemd_X_exp_to_sim_weight, hemd_X_exp_to_sim_noweight);

  delete hemd_X_exp_to_sim_noweight;
  return hemd_X_exp_to_sim_weight;

}

TH1F* EMDKLAnalyser::get_min_exp_to_sim(Int_t const sim_trial_sample_size, TString hname) {
  return get_X_exp_to_sim(sim_trial_sample_size, "min", hname);
}
TH1F* EMDKLAnalyser::get_max_exp_to_sim(Int_t const sim_trial_sample_size) {
  return get_X_exp_to_sim(sim_trial_sample_size, "max");
}
TH1F* EMDKLAnalyser::get_ave_exp_to_sim(Int_t const sim_trial_sample_size) {
  return get_X_exp_to_sim(sim_trial_sample_size, "ave");
}

void EMDKLAnalyser::set_proper_histo_uncertanties(TH1F* hweight, TH1F *hnoweight) {

  auto const nxbins1 = hweight->GetXaxis()->GetNbins();
  auto const nxbins2 = hnoweight->GetXaxis()->GetNbins();
  if (nxbins1!=nxbins2) cout<<"number of bins mismatch!"<<endl;

  for (int i = 1; i <= nxbins1; i++) { //no under- and overflow
    auto binc_weight   = hweight->GetBinContent(i);
    auto binc_noweight = hnoweight->GetBinContent(i);
    auto proper_absolute_error_weight = 0.0;
    if (binc_noweight > 0) {
      proper_absolute_error_weight = binc_weight * sqrt(binc_noweight) / binc_noweight;
    }
    hweight->SetBinError(i, proper_absolute_error_weight);
  }

}

vector<Float_t> EMDKLAnalyser::calc_sample_average() {

  auto dim1_size = emd_values_sim.size(); //rectangular for the moment

  vector<Float_t> average_values;

  for (int i = 0; i < dim1_size; i++) {
    float sum = 0.0;
    int counter = 0;
    for (int j = 0; j < dim1_size; j++) {
      if (i==j) continue;
      sum += emd_values_sim[i][j];
      counter++;
    }
    average_values.push_back( sum/(1.0*counter) );
  }

  return average_values;

}

vector<Float_t> EMDKLAnalyser::sample_emd_values() {

  vector<Float_t> emd_values;

  int dim1_size = emd_values_sim.size();

  for (int i = 0; i < dim1_size; i++) {
    for (int j = 0; j < dim1_size; j++) {
      if (i==j) continue;
      emd_values.push_back( emd_values_sim[i][j] );
    }
  }

  return emd_values;

}
*/


