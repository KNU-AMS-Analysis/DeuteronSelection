
void make_data(char* input_filename, char* output_filename)
{
  TDatime* time = new TDatime();
  cout<<"it start : "<<time->GetDate()<<"/"<<time->GetHour()<<":"<<time->GetMinute()<<":"<<time->GetSecond()<<endl;
  TStopwatch* stopwatch = new TStopwatch();
  stopwatch->Start();
  //TString fname = "/Users/sckang/work/Data/testrun.root";
  TString fname = "/Users/sckang/work/Data/proton.root";
  //TString fname = "/Users/sckang/work/Data/proton_flux.root";

  TFile* f = new TFile(input_filename);
  TTree* tree = (TTree*)f->Get("t");

  //select event cut///////////////////////////////////////////////////////////////////
  bool cut_tof_chi2;  //tof chi2(time) <=10 && tof chi2(coord) <= 10 : true
  bool cut_trk_charge; //0.7<trk_charge<1.5
  bool cut_trk_yhit;  //trk yhit > 5
  bool cut_trk_chi2;  //trk chi2 < 10
  bool cut_tof_charge; //0.5 < tof lower_tof_charge < 3 && 1.25 < tof upper charge
  bool cut_trk_hitL1XY; //hit trk layerXY 1
  bool cut_trk_hitL9XY; //hit trk layerXY 9
  bool cut_tof_cluster; //Ntofcluster(particle)==4 && Ntofcluster in time == 4
  bool cut_tof_beta;
  Int_t trd_vertex3;  //trd vertex == 0
  Int_t ntrdtrack;  //ntrd track == 1
  Int_t ntrtrack;  //ntracker track == 1
  Int_t nacchits;  //nacchits == 0
  Int_t rich_isnaf_d; // rich is naf = 1 not = 0
  Int_t rich_passcut_d;  //rich pass cut == 1
  Int_t trd_mc_cut;
  Int_t trk_mc_cut; // trk particle ID == primary ID
  Int_t tof_mc_cut; // tof particle ID == primary ID
  Int_t rich_mc_cut_d; // rich particle ID == primary ID
  /////////////////////////////////////////////////////////////////////////////////////

  //get value////////////////////////////////////////////////////////////////////////
  Float_t init_momentum_d, beta_d, rich_beta_d, ecal_bdt_d, charge_tof_d, charge_trk_d;
  Float_t trk_edep[9], tof_edep[4];
  Double_t rigidity_max_d, trdk_e2plikelihood_d, trdk_he2plikelihood_d, trdp_e2plikelihood_d, trdp_he2plikelihood_d;
  bool cut_trk_hitL1_d;
  bool cut_trk_hitL9_d;
  Double_t trk_theta_d, trk_phi_d;
  Int_t trk_pattern_X_d, trk_pattern_Y_d, trk_pattern_XY_d, Ntof_cluster_d, Ntof_cluster_intime_d, beta_pattern_d, rich_used_hits_d, rich_Nhits_ring_d, rich_Npmts_ring_d, tof_used_hit_d, tof_sum_hit_d;
  Float_t tof_chi2_time_d, tof_chi2_coord_d, tof_mass_err_d, tof_mass_d, rich_photoelectrons_d, rich_expect_photoelectrons_d, rich_prob_d, rich_track_theta_d, rich_track_phi_d; 
  Float_t tof_cluster_Aedep_d[4], tof_cluster_Dedep_d[4], tof_cluster_hit_time_d[4], tof_cluster_hit_time_err_d[4], tof_cluster_non_correct_time_d[2][4], tof_cluster_correct_time_d[2][4], tof_cluster_Acharge_d[2][4], tof_cluster_DCharge_d[2][3][4];
  Int_t tof_cluster_pattern_d[4], tof_cluster_layer_d[4];
  Float_t charge_tofL_d[4];
  Double_t trk_chi2_max_X_d, trk_chi2_max_Y_d, trk_chi2_inner_X_d, trk_chi2_inner_Y_d;
  Float_t inv_beta_err_d, rich_beta_err_d;
  Double_t inv_rigidity_inner_err_d, inv_rigidity_max_err_d, rigidity_inner_d; 
  Int_t primary_ID, particle_ID;
  Float_t mc_prob[2];
  Float_t ecal_shower_energy_d, ecal_layer_energy_d[18];
  /////////////////////////////////////////////////////////////////////////////////////


  //saved value//////////////////////////////////////////////////////////////////////////
  Float_t beta, mean_trk_edep, mean_tof_edep, rich_beta_naf, rich_beta_aer, init_momentum, ecal_bdt, rich_beta, charge_tof, charge_trk, inv_beta_err, rich_beta_err, inv_beta, rigidity_beta;
  Double_t rigidity_max, rigidity_inner, inv_rigidity_inner_err, inv_rigidity_max_err, inv_rigidity_max, inv_rigidity_inner, trdk_e2plikelihood, trdk_he2plikelihood, trdp_e2plikelihood, trdp_he2plikelihood;
  bool cut_trk_hitL1;
  bool cut_trk_hitL9;
  Int_t rich_passcut, rich_isnaf, trk_pattern_X, trk_pattern_Y, trk_pattern_XY, Ntof_cluster, Ntof_cluster_intime, beta_pattern, tof_cluster_pattern1, tof_cluster_pattern2, tof_cluster_pattern3, tof_cluster_pattern4, tof_cluster_layer1, tof_cluster_layer2, tof_cluster_layer3, tof_cluster_layer4, rich_used_hits, rich_Nhits_ring, rich_Npmts_ring, tof_used_hit, tof_sum_hit, rich_mc_cut;
  Float_t tof_chi2_time, tof_chi2_coord, rich_mass, rich_mass_err, tof_mass_err, tof_mass, tof_cluster_Aedep1, tof_cluster_Aedep2, tof_cluster_Aedep3, tof_cluster_Aedep4, tof_cluster_Dedep1, tof_cluster_Dedep2, tof_cluster_Dedep3, tof_cluster_Dedep4, tof_cluster_hit_time1, tof_cluster_hit_time2, tof_cluster_hit_time3, tof_cluster_hit_time4, tof_cluster_hit_time_err1, tof_cluster_hit_time_err2, tof_cluster_hit_time_err3, tof_cluster_hit_time_err4, tof_cluster_non_correct_time1, tof_cluster_non_correct_time2, tof_cluster_non_correct_time3, tof_cluster_non_correct_time4, tof_cluster_correct_time1, tof_cluster_correct_time2, tof_cluster_correct_time3, tof_cluster_correct_time4, tof_cluster_Acharge1, tof_cluster_Acharge2, tof_cluster_Acharge3, tof_cluster_Acharge4, tof_cluster_Dcharge1, tof_cluster_Dcharge2, tof_cluster_Dcharge3, tof_cluster_Dcharge4, rich_photoelectrons, rich_expect_photoelectrons, rich_prob, rich_track_theta, rich_track_phi;
  Double_t trk_theta, trk_phi;
  Float_t tof_cluster_non_correct_time11, tof_cluster_non_correct_time12, tof_cluster_non_correct_time13, tof_cluster_non_correct_time14;
  Float_t tof_cluster_non_correct_time21, tof_cluster_non_correct_time22, tof_cluster_non_correct_time23, tof_cluster_non_correct_time24;
  Float_t tof_cluster_correct_time11, tof_cluster_correct_time12, tof_cluster_correct_time13, tof_cluster_correct_time14;
  Float_t tof_cluster_correct_time21, tof_cluster_correct_time22, tof_cluster_correct_time23, tof_cluster_correct_time24;
  Float_t tof_cluster_Acharge11, tof_cluster_Acharge12, tof_cluster_Acharge13, tof_cluster_Acharge14;
  Float_t tof_cluster_Acharge21, tof_cluster_Acharge22, tof_cluster_Acharge23, tof_cluster_Acharge24;
  Float_t tof_cluster_Dcharge111, tof_cluster_Dcharge112, tof_cluster_Dcharge113, tof_cluster_Dcharge114;
  Float_t tof_cluster_Dcharge121, tof_cluster_Dcharge122, tof_cluster_Dcharge123, tof_cluster_Dcharge124;
  Float_t tof_cluster_Dcharge131, tof_cluster_Dcharge132, tof_cluster_Dcharge133, tof_cluster_Dcharge134;
  Float_t tof_cluster_Dcharge211, tof_cluster_Dcharge212, tof_cluster_Dcharge213, tof_cluster_Dcharge214;
  Float_t tof_cluster_Dcharge221, tof_cluster_Dcharge222, tof_cluster_Dcharge223, tof_cluster_Dcharge224;
  Float_t tof_cluster_Dcharge231, tof_cluster_Dcharge232, tof_cluster_Dcharge233, tof_cluster_Dcharge234;
  Float_t charge_tofL1, charge_tofL2, charge_tofL3, charge_tofL4;
  Double_t trk_chi2_max_X, trk_chi2_max_Y, trk_chi2_inner_X, trk_chi2_inner_Y;
  Float_t ecal_shower_energy, mean_ecal_energy;
  Float_t diff_rigbeta_tofbeta;
  Double_t diff_trkchi2Y_inner_max;
  ///////////////////////////////////////////////////////////////////////////////////////

  TFile* output = new TFile(output_filename, "recreate");
  //TFile* output = new TFile("data_mc_p0_10000_highcut.root", "recreate");
  //TFile* output = new TFile("data_mc_p1_20_rig_resol.root", "recreate");
  //TFile* output = new TFile("data_mc_p0_10000_rig_resol.root", "recreate");
  //TFile* output = new TFile("data_mc_p_pure.root", "recreate");
  //TFile* output = new TFile("data_mc_p0.root", "recreate");

  TTree* mc_d_tree = new TTree("mc_d_data", "mc_d_data");


  //make tree///////////////////////////////////////////////////////////////////////
  mc_d_tree->Branch("cut_trk_hitL1", &cut_trk_hitL1, "cut_trk_hitL1/O");
  mc_d_tree->Branch("cut_trk_hitL9", &cut_trk_hitL9, "cut_trk_hitL9/O");
  mc_d_tree->Branch("rich_passcut", &rich_passcut, "rich_passcut/I");
  mc_d_tree->Branch("rich_isnaf", &rich_isnaf, "rich_isnaf/I");
  mc_d_tree->Branch("trk_pattern_X", &trk_pattern_X, "trk_pattern_X/I");
  mc_d_tree->Branch("trk_pattern_Y", &trk_pattern_Y, "trk_pattern_Y/I");
  mc_d_tree->Branch("trk_pattern_XY", &trk_pattern_XY, "trk_pattern_XY/I");
  mc_d_tree->Branch("Ntof_cluster", &Ntof_cluster, "Ntof_cluster/I");
  mc_d_tree->Branch("Ntof_cluster_intime", &Ntof_cluster_intime, "Ntof_cluster_intime/I");
  mc_d_tree->Branch("beta_pattern", &beta_pattern, "beta_pattern/I");
  mc_d_tree->Branch("tof_cluster_pattern1", &tof_cluster_pattern1, "tof_cluster_pattern1/I");
  mc_d_tree->Branch("tof_cluster_pattern2", &tof_cluster_pattern2, "tof_cluster_pattern2/I");
  mc_d_tree->Branch("tof_cluster_pattern3", &tof_cluster_pattern3, "tof_cluster_pattern3/I");
  mc_d_tree->Branch("tof_cluster_pattern4", &tof_cluster_pattern4, "tof_cluster_pattern4/I");
  mc_d_tree->Branch("tof_cluster_layer1", &tof_cluster_layer1, "tof_cluster_layer1/I");
  mc_d_tree->Branch("tof_cluster_layer2", &tof_cluster_layer2, "tof_cluster_layer2/I");
  mc_d_tree->Branch("tof_cluster_layer3", &tof_cluster_layer3, "tof_cluster_layer3/I");
  mc_d_tree->Branch("tof_cluster_layer4", &tof_cluster_layer4, "tof_cluster_layer4/I");
  mc_d_tree->Branch("rich_used_hits", &rich_used_hits, "rich_used_hits/I");
  mc_d_tree->Branch("rich_Nhits_ring", &rich_Nhits_ring, "rich_Nhits_ring/I");
  mc_d_tree->Branch("rich_Npmts_ring", &rich_Npmts_ring, "rich_Npmts_ring/I");
  mc_d_tree->Branch("tof_used_hit", &tof_used_hit, "tof_used_hit/I");
  mc_d_tree->Branch("tof_sum_hit", &tof_sum_hit, "tof_sum_hit/I");
  mc_d_tree->Branch("rich_mc_cut", &rich_mc_cut, "rich_mc_cut/I");
  mc_d_tree->Branch("init_momentum", &init_momentum, "init_momentum/F");
  mc_d_tree->Branch("beta", &beta, "beta/F");
  mc_d_tree->Branch("rich_beta", &rich_beta, "rich_beta/F");
  mc_d_tree->Branch("rich_beta_naf", &rich_beta_naf, "rich_beta_naf/F");
  mc_d_tree->Branch("rich_beta_aer", &rich_beta_aer, "rich_beta_aer/F");
  mc_d_tree->Branch("mean_trk_edep", &mean_trk_edep, "mean_trk_edep/F");
  mc_d_tree->Branch("mean_tof_edep", &mean_tof_edep, "mean_tof_edep/F");
  mc_d_tree->Branch("ecal_bdt", &ecal_bdt, "ecal_bdt/F");
  mc_d_tree->Branch("tof_chi2_time", &tof_chi2_time, "tof_chi2_time/F");
  mc_d_tree->Branch("tof_chi2_coord", &tof_chi2_coord, "tof_chi2_coord/F");
  mc_d_tree->Branch("tof_mass", &tof_mass, "tof_mass/F");
  mc_d_tree->Branch("tof_mass_err", &tof_mass_err, "tof_mass_err/F");
  mc_d_tree->Branch("rich_mass", &rich_mass, "rich_mass/F");
  mc_d_tree->Branch("rich_mass_err", &rich_mass_err, "rich_mass_err/F");
  mc_d_tree->Branch("tof_cluster_Aedep1", &tof_cluster_Aedep1, "tof_cluster_Aedep1/F");
  mc_d_tree->Branch("tof_cluster_Aedep2", &tof_cluster_Aedep2, "tof_cluster_Aedep2/F");
  mc_d_tree->Branch("tof_cluster_Aedep3", &tof_cluster_Aedep3, "tof_cluster_Aedep3/F");
  mc_d_tree->Branch("tof_cluster_Aedep4", &tof_cluster_Aedep4, "tof_cluster_Aedep4/F");
  mc_d_tree->Branch("tof_cluster_Dedep1", &tof_cluster_Dedep1, "tof_cluster_Dedep1/F");
  mc_d_tree->Branch("tof_cluster_Dedep2", &tof_cluster_Dedep2, "tof_cluster_Dedep2/F");
  mc_d_tree->Branch("tof_cluster_Dedep3", &tof_cluster_Dedep3, "tof_cluster_Dedep3/F");
  mc_d_tree->Branch("tof_cluster_Dedep4", &tof_cluster_Dedep4, "tof_cluster_Dedep4/F");
  mc_d_tree->Branch("tof_cluster_hit_time1", &tof_cluster_hit_time1, "tof_cluster_hit_time1/F");
  mc_d_tree->Branch("tof_cluster_hit_time2", &tof_cluster_hit_time2, "tof_cluster_hit_time2/F");
  mc_d_tree->Branch("tof_cluster_hit_time3", &tof_cluster_hit_time3, "tof_cluster_hit_time3/F");
  mc_d_tree->Branch("tof_cluster_hit_time4", &tof_cluster_hit_time4, "tof_cluster_hit_time4/F");
  mc_d_tree->Branch("tof_cluster_hit_time_err1", &tof_cluster_hit_time_err1, "tof_cluster_hit_time_err1/F");
  mc_d_tree->Branch("tof_cluster_hit_time_err2", &tof_cluster_hit_time_err2, "tof_cluster_hit_time_err2/F");
  mc_d_tree->Branch("tof_cluster_hit_time_err3", &tof_cluster_hit_time_err3, "tof_cluster_hit_time_err3/F");
  mc_d_tree->Branch("tof_cluster_hit_time_err4", &tof_cluster_hit_time_err4, "tof_cluster_hit_time_err4/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time11", &tof_cluster_non_correct_time11, "tof_cluster_non_correct_time11/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time12", &tof_cluster_non_correct_time12, "tof_cluster_non_correct_time12/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time13", &tof_cluster_non_correct_time13, "tof_cluster_non_correct_time13/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time14", &tof_cluster_non_correct_time14, "tof_cluster_non_correct_time14/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time21", &tof_cluster_non_correct_time21, "tof_cluster_non_correct_time21/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time22", &tof_cluster_non_correct_time22, "tof_cluster_non_correct_time22/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time23", &tof_cluster_non_correct_time23, "tof_cluster_non_correct_time23/F");
  mc_d_tree->Branch("tof_cluster_non_correct_time24", &tof_cluster_non_correct_time24, "tof_cluster_non_correct_time24/F");
  mc_d_tree->Branch("tof_cluster_correct_time11", &tof_cluster_correct_time11, "tof_cluster_correct_time11/F");
  mc_d_tree->Branch("tof_cluster_correct_time12", &tof_cluster_correct_time12, "tof_cluster_correct_time12/F");
  mc_d_tree->Branch("tof_cluster_correct_time13", &tof_cluster_correct_time13, "tof_cluster_correct_time13/F");
  mc_d_tree->Branch("tof_cluster_correct_time14", &tof_cluster_correct_time14, "tof_cluster_correct_time14/F");
  mc_d_tree->Branch("tof_cluster_correct_time21", &tof_cluster_correct_time21, "tof_cluster_correct_time21/F");
  mc_d_tree->Branch("tof_cluster_correct_time22", &tof_cluster_correct_time22, "tof_cluster_correct_time22/F");
  mc_d_tree->Branch("tof_cluster_correct_time23", &tof_cluster_correct_time23, "tof_cluster_correct_time23/F");
  mc_d_tree->Branch("tof_cluster_correct_time24", &tof_cluster_correct_time24, "tof_cluster_correct_time24/F");
  mc_d_tree->Branch("tof_cluster_Acharge11", &tof_cluster_Acharge11, "tof_cluster_Acharge11/F");
  mc_d_tree->Branch("tof_cluster_Acharge12", &tof_cluster_Acharge12, "tof_cluster_Acharge12/F");
  mc_d_tree->Branch("tof_cluster_Acharge13", &tof_cluster_Acharge13, "tof_cluster_Acharge13/F");
  mc_d_tree->Branch("tof_cluster_Acharge14", &tof_cluster_Acharge14, "tof_cluster_Acharge14/F");
  mc_d_tree->Branch("tof_cluster_Acharge21", &tof_cluster_Acharge21, "tof_cluster_Acharge21/F");
  mc_d_tree->Branch("tof_cluster_Acharge22", &tof_cluster_Acharge22, "tof_cluster_Acharge22/F");
  mc_d_tree->Branch("tof_cluster_Acharge23", &tof_cluster_Acharge23, "tof_cluster_Acharge23/F");
  mc_d_tree->Branch("tof_cluster_Acharge24", &tof_cluster_Acharge24, "tof_cluster_Acharge24/F");
  mc_d_tree->Branch("tof_cluster_Dcharge111", &tof_cluster_Dcharge111, "tof_cluster_Dcharge111/F");
  mc_d_tree->Branch("tof_cluster_Dcharge112", &tof_cluster_Dcharge112, "tof_cluster_Dcharge112/F");
  mc_d_tree->Branch("tof_cluster_Dcharge113", &tof_cluster_Dcharge113, "tof_cluster_Dcharge113/F");
  mc_d_tree->Branch("tof_cluster_Dcharge114", &tof_cluster_Dcharge114, "tof_cluster_Dcharge114/F");
  mc_d_tree->Branch("tof_cluster_Dcharge121", &tof_cluster_Dcharge121, "tof_cluster_Dcharge121/F");
  mc_d_tree->Branch("tof_cluster_Dcharge122", &tof_cluster_Dcharge122, "tof_cluster_Dcharge122/F");
  mc_d_tree->Branch("tof_cluster_Dcharge123", &tof_cluster_Dcharge123, "tof_cluster_Dcharge123/F");
  mc_d_tree->Branch("tof_cluster_Dcharge124", &tof_cluster_Dcharge124, "tof_cluster_Dcharge124/F");
  mc_d_tree->Branch("tof_cluster_Dcharge131", &tof_cluster_Dcharge131, "tof_cluster_Dcharge131/F");
  mc_d_tree->Branch("tof_cluster_Dcharge132", &tof_cluster_Dcharge132, "tof_cluster_Dcharge132/F");
  mc_d_tree->Branch("tof_cluster_Dcharge133", &tof_cluster_Dcharge133, "tof_cluster_Dcharge133/F");
  mc_d_tree->Branch("tof_cluster_Dcharge134", &tof_cluster_Dcharge134, "tof_cluster_Dcharge134/F");
  mc_d_tree->Branch("tof_cluster_Dcharge211", &tof_cluster_Dcharge211, "tof_cluster_Dcharge211/F");
  mc_d_tree->Branch("tof_cluster_Dcharge212", &tof_cluster_Dcharge212, "tof_cluster_Dcharge212/F");
  mc_d_tree->Branch("tof_cluster_Dcharge213", &tof_cluster_Dcharge213, "tof_cluster_Dcharge213/F");
  mc_d_tree->Branch("tof_cluster_Dcharge214", &tof_cluster_Dcharge214, "tof_cluster_Dcharge214/F");
  mc_d_tree->Branch("tof_cluster_Dcharge221", &tof_cluster_Dcharge221, "tof_cluster_Dcharge221/F");
  mc_d_tree->Branch("tof_cluster_Dcharge222", &tof_cluster_Dcharge222, "tof_cluster_Dcharge222/F");
  mc_d_tree->Branch("tof_cluster_Dcharge223", &tof_cluster_Dcharge223, "tof_cluster_Dcharge223/F");
  mc_d_tree->Branch("tof_cluster_Dcharge224", &tof_cluster_Dcharge224, "tof_cluster_Dcharge224/F");
  mc_d_tree->Branch("tof_cluster_Dcharge231", &tof_cluster_Dcharge231, "tof_cluster_Dcharge231/F");
  mc_d_tree->Branch("tof_cluster_Dcharge232", &tof_cluster_Dcharge232, "tof_cluster_Dcharge232/F");
  mc_d_tree->Branch("tof_cluster_Dcharge233", &tof_cluster_Dcharge233, "tof_cluster_Dcharge233/F");
  mc_d_tree->Branch("tof_cluster_Dcharge234", &tof_cluster_Dcharge234, "tof_cluster_Dcharge234/F");
  mc_d_tree->Branch("rich_photoelectrons", &rich_photoelectrons, "rich_photoelectrons/F");
  mc_d_tree->Branch("rich_expect_photoelectrons", &rich_expect_photoelectrons, "rich_expect_photoelectrons/F");
  mc_d_tree->Branch("rich_prob", &rich_prob, "rich_prob/F");
  mc_d_tree->Branch("rich_track_theta", &rich_track_theta, "rich_track_theta/F");
  mc_d_tree->Branch("rich_track_phi", &rich_track_phi, "rich_track_phi/F");
  mc_d_tree->Branch("charge_tof", &charge_tof, "charge_tof/F");
  mc_d_tree->Branch("charge_trk", &charge_trk, "charge_trk/F");
  mc_d_tree->Branch("charge_tofL1", &charge_tofL1, "charge_tofL1/F");
  mc_d_tree->Branch("charge_tofL2", &charge_tofL2, "charge_tofL2/F");
  mc_d_tree->Branch("charge_tofL3", &charge_tofL3, "charge_tofL3/F");
  mc_d_tree->Branch("charge_tofL4", &charge_tofL4, "charge_tofL4/F");
  mc_d_tree->Branch("inv_beta_err", &inv_beta_err, "inv_beta_err/F");
  mc_d_tree->Branch("rich_beta_err", &rich_beta_err, "rich_beta_err/F");
  mc_d_tree->Branch("inv_beta", &inv_beta, "inv_beta/F");
  mc_d_tree->Branch("rigidity_beta", &rigidity_beta, "rigidity_beta/F");
  mc_d_tree->Branch("mean_ecal_energy", &mean_ecal_energy, "mean_ecal_energy/F");
  mc_d_tree->Branch("ecal_shower_energy", &ecal_shower_energy, "ecal_shower_energy/F");
  mc_d_tree->Branch("diff_rigbeta_tofbeta", &diff_rigbeta_tofbeta, "diff_rigbeta_tofbeta/F");
  mc_d_tree->Branch("rigidity_max", &rigidity_max, "rigidity_max/D");
  mc_d_tree->Branch("rigidity_inner", &rigidity_inner, "rigidity_inner/D");
  mc_d_tree->Branch("inv_rigidity_max", &inv_rigidity_max, "inv_rigidity_max/D");
  mc_d_tree->Branch("inv_rigidity_inner", &inv_rigidity_inner, "inv_rigidity_inner/D");
  mc_d_tree->Branch("inv_rigidity_max_err", &inv_rigidity_max_err, "inv_rigidity_max_err/D");
  mc_d_tree->Branch("inv_rigidity_inner_err", &inv_rigidity_inner_err, "inv_rigidity_inner_err/D");
  mc_d_tree->Branch("trdk_e2plikelihood", &trdk_e2plikelihood, "trdk_e2plikelihood/D");
  mc_d_tree->Branch("trdk_he2plikelihood", &trdk_he2plikelihood, "trdk_he2plikelihood/D");
  mc_d_tree->Branch("trdp_e2plikelihood", &trdp_e2plikelihood, "trdp_e2plikelihood/D");
  mc_d_tree->Branch("trdp_he2plikelihood", &trdp_he2plikelihood, "trdp_he2plikelihood/D");
  mc_d_tree->Branch("trk_chi2_max_X", &trk_chi2_max_X, "trk_chi2_max_X/D");
  mc_d_tree->Branch("trk_chi2_max_Y", &trk_chi2_max_Y, "trk_chi2_max_Y/D");
  mc_d_tree->Branch("trk_chi2_inner_X", &trk_chi2_inner_X, "trk_chi2_inner_X/D");
  mc_d_tree->Branch("trk_chi2_inner_Y", &trk_chi2_inner_Y, "trk_chi2_inner_Y/D");
  mc_d_tree->Branch("diff_trkchi2Y_inner_max", &diff_trkchi2Y_inner_max, "diff_trk2Y_inner_max/D");
  //////////////////////////////////////////////////////////////////////////////////

  //get tree///////////////////////////////////////////////////////////////////////
  tree->SetBranchAddress("init_momentum", &init_momentum_d);
  tree->SetBranchAddress("beta", &beta_d);
  tree->SetBranchAddress("rich_beta", &rich_beta_d);
  tree->SetBranchAddress("ecal_bdt", &ecal_bdt_d);
  tree->SetBranchAddress("trk_edep", &trk_edep);
  tree->SetBranchAddress("tof_edep", &tof_edep);
  tree->SetBranchAddress("rigidity_max", &rigidity_max_d);
  tree->SetBranchAddress("rigidity_inner", &rigidity_inner_d);
  tree->SetBranchAddress("cut_tof_chi2", &cut_tof_chi2);
  tree->SetBranchAddress("cut_trk_charge", &cut_trk_charge);
  tree->SetBranchAddress("cut_trk_yhit5", &cut_trk_yhit);
  tree->SetBranchAddress("cut_trk_chi2", &cut_trk_chi2);
  tree->SetBranchAddress("cut_tof_charge", &cut_tof_charge);
  tree->SetBranchAddress("cut_trk_hitL1", &cut_trk_hitL1_d);
  tree->SetBranchAddress("cut_trk_hitL9", &cut_trk_hitL9_d);
  tree->SetBranchAddress("cut_tof_cluster", &cut_tof_cluster);
  tree->SetBranchAddress("trd_vertex3", &trd_vertex3);
  tree->SetBranchAddress("ntrdtrack", &ntrdtrack);
  tree->SetBranchAddress("ntrtrack", &ntrtrack);
  tree->SetBranchAddress("nacchits", &nacchits);
  tree->SetBranchAddress("rich_isnaf", &rich_isnaf_d);
  tree->SetBranchAddress("rich_passcut", &rich_passcut_d);
  tree->SetBranchAddress("trdk_e2plikelihood", &trdk_e2plikelihood_d);
  tree->SetBranchAddress("trdk_he2plikelihood", &trdk_he2plikelihood_d);
  tree->SetBranchAddress("trdp_e2plikelihood", &trdp_e2plikelihood_d);
  tree->SetBranchAddress("trdp_he2plikelihood", &trdp_he2plikelihood_d);
  tree->SetBranchAddress("charge_trk", &charge_trk_d);
  tree->SetBranchAddress("charge_tof", &charge_tof_d);
  tree->SetBranchAddress("charge_tofL", &charge_tofL_d);
  tree->SetBranchAddress("trk_theta", &trk_theta_d);
  tree->SetBranchAddress("trk_phi", &trk_phi_d);
  tree->SetBranchAddress("trk_pattern_X", &trk_pattern_X_d);
  tree->SetBranchAddress("trk_pattern_Y", &trk_pattern_Y_d);
  tree->SetBranchAddress("trk_pattern_XY", &trk_pattern_XY_d);
  tree->SetBranchAddress("Ntof_cluster", &Ntof_cluster_d);
  tree->SetBranchAddress("Ntof_cluster_intime", &Ntof_cluster_intime_d);
  tree->SetBranchAddress("beta_pattern", &beta_pattern_d);
  tree->SetBranchAddress("tof_chi2_time", &tof_chi2_time_d);
  tree->SetBranchAddress("tof_chi2_coord", &tof_chi2_coord_d);
  tree->SetBranchAddress("tof_used_hit", &tof_used_hit_d);
  tree->SetBranchAddress("tof_sum_hit", &tof_sum_hit_d);
  tree->SetBranchAddress("tof_cluster_pattern", &tof_cluster_pattern_d);
  tree->SetBranchAddress("tof_cluster_layer", &tof_cluster_layer_d);
  tree->SetBranchAddress("tof_cluster_Aedep", &tof_cluster_Aedep_d);
  tree->SetBranchAddress("tof_cluster_Dedep", &tof_cluster_Dedep_d);
  tree->SetBranchAddress("tof_cluster_hit_time", &tof_cluster_hit_time_d);
  tree->SetBranchAddress("tof_cluster_hit_time_err", &tof_cluster_hit_time_err_d);
  tree->SetBranchAddress("tof_cluster_non_correct_time", &tof_cluster_non_correct_time_d);
  tree->SetBranchAddress("tof_cluster_correct_time", &tof_cluster_correct_time_d);
  tree->SetBranchAddress("tof_cluster_Acharge", &tof_cluster_Acharge_d);
  tree->SetBranchAddress("tof_cluster_DCharge", &tof_cluster_DCharge_d);
  tree->SetBranchAddress("rich_used_hits", &rich_used_hits_d);
  tree->SetBranchAddress("rich_photoelectrons", &rich_photoelectrons_d);
  tree->SetBranchAddress("rich_expect_photoelectrons", &rich_expect_photoelectrons_d);
  tree->SetBranchAddress("rich_prob", &rich_prob_d);
  tree->SetBranchAddress("rich_track_theta", &rich_track_theta_d);
  tree->SetBranchAddress("rich_track_phi", &rich_track_phi_d);
  tree->SetBranchAddress("rich_Nhits_ring", &rich_Nhits_ring_d);
  tree->SetBranchAddress("rich_Npmts_ring", &rich_Npmts_ring_d);
  tree->SetBranchAddress("cut_tof_beta", &cut_tof_beta);
  tree->SetBranchAddress("red_chi2x_max", &trk_chi2_max_X_d);
  tree->SetBranchAddress("red_chi2y_max", &trk_chi2_max_Y_d);
  tree->SetBranchAddress("red_chi2x_inner", &trk_chi2_inner_X_d);
  tree->SetBranchAddress("red_chi2y_inner", &trk_chi2_inner_Y_d);
  tree->SetBranchAddress("inv_beta_err", &inv_beta_err_d);
  tree->SetBranchAddress("rich_beta_err", &rich_beta_err_d);
  tree->SetBranchAddress("inv_rigidity_max_err", &inv_rigidity_max_err_d);
  tree->SetBranchAddress("inv_rigidity_inner_err", &inv_rigidity_inner_err_d);
  tree->SetBranchAddress("particle_ID", &particle_ID);
  tree->SetBranchAddress("primary_ID", &primary_ID);
  tree->SetBranchAddress("mc_prob", &mc_prob);
  tree->SetBranchAddress("trk_mc_cut", &trk_mc_cut);
  tree->SetBranchAddress("tof_mc_cut", &tof_mc_cut);
  tree->SetBranchAddress("rich_mc_cut", &rich_mc_cut_d);
  tree->SetBranchAddress("trd_mc_cut", &trd_mc_cut);
  tree->SetBranchAddress("ecal_Edep", &ecal_layer_energy_d);
  tree->SetBranchAddress("E_ecal", &ecal_shower_energy_d);

  /////////////////////////////////////////////////////////////////////////////////

  TH1F* hEvtCounter = new TH1F("hEvtCounter", "event counter",19,1,20);
  hEvtCounter->GetXaxis()->SetBinLabel(1, "cut tof chi2");
  hEvtCounter->GetXaxis()->SetBinLabel(2, "cut trk charge");
  hEvtCounter->GetXaxis()->SetBinLabel(3, "cut trk yhit");
  hEvtCounter->GetXaxis()->SetBinLabel(4, "cut trk chi2");
  hEvtCounter->GetXaxis()->SetBinLabel(5, "cut tof charge");
  hEvtCounter->GetXaxis()->SetBinLabel(6, "cut tof cluster");
  hEvtCounter->GetXaxis()->SetBinLabel(7, "trd vertex");
  hEvtCounter->GetXaxis()->SetBinLabel(8, "ntrd track");
  hEvtCounter->GetXaxis()->SetBinLabel(9, "ntrk track");
  hEvtCounter->GetXaxis()->SetBinLabel(10, "nacc hits");
  hEvtCounter->GetXaxis()->SetBinLabel(11, "cut tof beta");
  hEvtCounter->GetXaxis()->SetBinLabel(12, "match particle ID");
  hEvtCounter->GetXaxis()->SetBinLabel(13, "rigidity cut");
  hEvtCounter->GetXaxis()->SetBinLabel(14, "ECAL BDT");
  hEvtCounter->GetXaxis()->SetBinLabel(15, "TRDk e2plikelihood");
  hEvtCounter->GetXaxis()->SetBinLabel(16, "TRDk he2plikelhood");
  hEvtCounter->GetXaxis()->SetBinLabel(17, "TRDp e2plikelihood");
  hEvtCounter->GetXaxis()->SetBinLabel(18, "TRDp he2plikelihood");
  hEvtCounter->GetXaxis()->SetBinLabel(19, "chi2 difference");

  Float_t bin[59] = {1.00e+00, 1.15e+00, 1.33e+00, 1.51e+00, 1.71e+00, 1.92e+00, 2.15e+00, 2.40e+00, 2.67e+00, 2.97e+00,
                     3.29e+00, 3.64e+00, 4.02e+00, 4.43e+00, 4.88e+00, 5.37e+00, 5.90e+00, 6.47e+00, 7.09e+00, 7.76e+00,
                     8.48e+00, 9.26e+00, 1.01e+01, 1.10e+01, 1.20e+01, 1.30e+01, 1.41e+01, 1.53e+01, 1.66e+01, 1.80e+01,
                     1.95e+01, 2.11e+01, 2.28e+01, 2.47e+01, 2.67e+01, 2.88e+01, 3.11e+01, 3.35e+01, 3.61e+01, 3.89e+01,
                     4.19e+01, 4.51e+01, 4.85e+01, 5.22e+01, 5.61e+01, 6.03e+01, 6.48e+01, 6.97e+01, 7.49e+01, 8.05e+01,
                     9.30e+01, 1.08e+02, 1.25e+02, 1.47e+02, 1.75e+02, 2.11e+02, 2.59e+02, 3.30e+02, 5.25e+02};

  fstream input;
  input.open("/u/user/sinchul/batch_job/MC_d/chi2_difference_d.txt", ios::in);
  Float_t chi2_diff_cut[31];
  for(Int_t i=0 ; i<31 ; i++)
  {
    input>>chi2_diff_cut[i];
  }

  fstream input1;
  input1.open("/Users/sckang/work/BDT/aachen/rigidity_resolution/rigidity_resolution.txt", ios::in);
  Float_t rigidity_resolution[31];
  for(Int_t i=0 ; i<31 ; i++)
  {
    input1>>rigidity_resolution[31];
  }

  fstream input2;
  input2.open("/Users/sckang/work/BDT/aachen/rigidity_resolution/rigidity_mean.txt", ios::in);
  Float_t rigidity_mean[31];
  for(Int_t i=0 ; i<31 ; i++)
  {
    input2>>rigidity_mean[31];
  }

  ULong64_t nentries = tree->GetEntries();
  ULong64_t x = nentries/100;
  Int_t y=0;

  cout<<"nentries : "<<nentries<<endl;

  //start event loop///////////////////////////////////////////////////////////////////////
  //start event loop///////////////////////////////////////////////////////////////////////
  //start event loop///////////////////////////////////////////////////////////////////////

  for(ULong64_t entry=0 ; entry<nentries ; entry++)
  //for(ULong64_t entry=0 ; entry<100000 ; entry++)
  {
    if(entry%(int)x==0)
    {
      cout<<y<<"%"<<endl;
      y++;
    }

    tree->GetEntry(entry);


    if(!cut_tof_chi2) continue;
    hEvtCounter->Fill(1);
    
    if(!cut_trk_charge) continue;
    hEvtCounter->Fill(2);

    if(!cut_trk_yhit) continue;
    hEvtCounter->Fill(3);

    if(!cut_trk_chi2) continue;
    hEvtCounter->Fill(4);

    if(!cut_tof_charge) continue;
    hEvtCounter->Fill(5);

    //if(!cut_tof_cluster) continue;
    if(Ntof_cluster_d != 4 || Ntof_cluster_intime_d > 2) continue;
    hEvtCounter->Fill(6);

    if(trd_vertex3 != 0) continue;
    hEvtCounter->Fill(7);

    if(ntrdtrack != 1) continue;
    hEvtCounter->Fill(8);

    if(ntrtrack != 1) continue;
    hEvtCounter->Fill(9);

    if(nacchits != 0) continue;
    hEvtCounter->Fill(10);
    if(!cut_tof_beta) continue;
    hEvtCounter->Fill(11);
    
    if(trd_mc_cut == 0 || trk_mc_cut == 0 || tof_mc_cut == 0) continue;
    hEvtCounter->Fill(12);

    if(rigidity_inner_d == -9999) continue;

    if(rigidity_max_d <1 || rigidity_max_d >20) continue;

    Int_t rigidity_cut=1;
    if(rigidity_inner_d < init_momentum_d*0.7 || rigidity_inner_d > init_momentum_d*1.3) rigidity_cut=0;
    //hEvtCounter->Fill(13);

    Int_t chi2_cut=1;
    for(Int_t i=0 ; i<31 ; i++)
    {
      if(rigidity_max_d >= bin[i] && rigidity_max_d < bin[i+1])
      {
        if(chi2_diff_cut[i] < TMath::Abs(trk_chi2_max_Y_d - trk_chi2_inner_Y_d)) chi2_cut=0;
      }
      if(init_momentum_d >= bin[i] && init_momentum_d < bin[i+1])
      {
        //if(rigidity_inner_d < init_momentum_d*0.7 || rigidity_inner_d > init_momentum_d*1.3) rigidity_cut=0;
      }
      if(init_momentum_d < 4.67 && init_momentum_d >= bin[i] && init_momentum_d < bin[i+1])
      {
        //if( TMath::Abs(init_momentum_d - rigidity_inner_d) / ((bin[i]+bin[i+1])/2) > 3*(0.36297*TMath::Power(init_momentum_d,6)-0.514811*TMath::Power(init_momentum_d,5)+0.449582*TMath::Power(init_momentum_d,4)-0.208271*TMath::Power(init_momentum_d,3)+0.0530592*TMath::Power(init_momentum_d,2)-0.00702338*init_momentum_d+0.000377701) ) rigidity_cut=0;
      }
      if(init_momentum_d >= 4.67 && init_momentum_d >= bin[i] && init_momentum_d < bin[i+1])
      {
        //if( TMath::Abs(init_momentum_d - rigidity_inner_d) / ((bin[i]+bin[i+1])/2) > 3*(0.100078*TMath::Power(init_momentum_d,2)+0.000612969*init_momentum_d+9.6167e-5) ) rigidity_cut=0;
      }

    }
    if(rigidity_cut==0) continue;
    hEvtCounter->Fill(13);

    //if(ecal_bdt_d >= -0.98) continue;
    hEvtCounter->Fill(14);

    //if(trdk_e2plikelihood_d <= 0.7) continue;
    hEvtCounter->Fill(15);

    //if(trdk_he2plikelihood_d >= 0.1) continue;
    hEvtCounter->Fill(16);

    //if(trdp_e2plikelihood_d <= 0.7) continue;
    hEvtCounter->Fill(17);

    //if(trdp_he2plikelihood_d >= 0.1) continue;
    hEvtCounter->Fill(18);

    if(chi2_cut==0) continue;
    hEvtCounter->Fill(19);

    beta = beta_d;
    inv_beta = 1/beta_d;
    rigidity_max = rigidity_max_d;
    rigidity_inner = rigidity_inner_d;
    inv_rigidity_max = 1/rigidity_max_d;
    inv_rigidity_inner = 1/rigidity_inner_d;
    inv_beta_err = inv_beta_err_d;
    rich_beta_err = rich_beta_err_d;
    inv_rigidity_max_err = inv_rigidity_max_err_d;
    inv_rigidity_inner_err = inv_rigidity_inner_err_d;
    rich_mc_cut = rich_mc_cut_d;
    rigidity_beta = rigidity_inner_d/TMath::Sqrt(TMath::Power( 0.938272297, 2) + TMath::Power( rigidity_inner_d ,2) ); 
    diff_rigbeta_tofbeta = TMath::Abs(rigidity_beta - beta_d);

    Int_t nhit=0;
    Float_t sum_tof_edep=0;

    for(Int_t i=0 ; i<2 ; i++)
    {
      if(tof_edep[i] == 0) continue;
      nhit++;
      sum_tof_edep += tof_edep[i];
    }

    mean_tof_edep = sum_tof_edep/nhit;


    nhit=0;
    Float_t sum_trk_edep=0;

    for(Int_t i=0 ; i<9 ; i++)
    {
      if(trk_edep[i] == 0) continue;
      nhit++;
      sum_trk_edep += trk_edep[i];
    }

    mean_trk_edep = sum_trk_edep/nhit;

    nhit=0;
    Float_t sum_ecal_edep=0;

    for(Int_t i=0 ; i<18 ; i++)
    {
      if(ecal_layer_energy_d[i]==0) continue;
      nhit++;
      sum_ecal_edep += ecal_layer_energy_d[i];
    }

    mean_ecal_energy = sum_ecal_edep/nhit;
    if(nhit==0) mean_ecal_energy=0;

    ecal_shower_energy = ecal_shower_energy_d;


    rich_beta_naf = -1;
    rich_beta_aer = -1;
    if(rich_passcut==1 && rich_isnaf_d==1) rich_beta_naf = rich_beta_d;
    if(rich_passcut==1 && rich_isnaf_d==0) rich_beta_aer = rich_beta_d;
    rich_beta = rich_beta_d;
    rich_isnaf = rich_isnaf_d;
    rich_passcut = rich_passcut_d;

    cut_trk_hitL1 = cut_trk_hitL1_d;
    cut_trk_hitL9 = cut_trk_hitL9_d;

    init_momentum = init_momentum_d;

    ecal_bdt = ecal_bdt_d;
    trdk_e2plikelihood = trdk_e2plikelihood_d;
    trdk_he2plikelihood = trdk_he2plikelihood_d;
    trdp_e2plikelihood = trdp_e2plikelihood_d;
    trdp_he2plikelihood = trdp_he2plikelihood_d;

    charge_trk = charge_trk_d;
    charge_tof = charge_tof_d;
    trk_theta = trk_theta_d;
    trk_phi = trk_phi_d;
    trk_pattern_X = trk_pattern_X_d;
    trk_pattern_Y = trk_pattern_Y_d;
    trk_pattern_XY = trk_pattern_XY_d;
    Ntof_cluster = Ntof_cluster_d;
    Ntof_cluster_intime = Ntof_cluster_intime_d;
    beta_pattern = beta_pattern_d;
    tof_chi2_time = tof_chi2_time_d;
    tof_chi2_coord = tof_chi2_coord_d;
    tof_mass=-1;
    if(beta_d <1) tof_mass = (charge_trk_d*rigidity_max_d*TMath::Sqrt(1-beta_d*beta_d))/beta_d;//tof_mass_d;
    //tof_mass = tof_mass_d;
    Float_t rigidity_err = inv_rigidity_max_err_d*TMath::Power(rigidity_max_d,2);
    Float_t beta_err = inv_beta_err_d*TMath::Power(beta_d,2);
    tof_mass_err=-1;
    if(beta_d <1) tof_mass_err = TMath::Sqrt( TMath::Power( rigidity_err*TMath::Sqrt(1-beta_d*beta_d)/beta_d , 2) + TMath::Power( rigidity_max_d*beta_err/( beta_d*beta_d*TMath::Sqrt(1-beta_d*beta_d)) ,2) );
    //tof_mass_err = tof_mass_err_d;
    rich_mass=-1;
    rich_mass_err=-1;
    if(rich_beta_d <1 && rich_beta_d != 0) rich_mass = (charge_trk_d*rigidity_max_d*TMath::Sqrt(1-rich_beta_d*rich_beta_d))/rich_beta_d;
    if(rich_beta_d <1 && rich_beta_d != 0) rich_mass_err = TMath::Sqrt(TMath::Power( rigidity_err*TMath::Sqrt(1-rich_beta_d*rich_beta_d)/rich_beta_d ,2) + TMath::Power( rigidity_max_d*rich_beta_err/( rich_beta_d*rich_beta_d*TMath::Sqrt(1-rich_beta_d*rich_beta_d)),2)); 
    tof_used_hit = tof_used_hit_d;
    tof_sum_hit = tof_sum_hit_d;
    tof_cluster_pattern1 = tof_cluster_pattern_d[0];
    tof_cluster_pattern2 = tof_cluster_pattern_d[1];
    tof_cluster_pattern3 = tof_cluster_pattern_d[2];
    tof_cluster_pattern4 = tof_cluster_pattern_d[3];
    tof_cluster_layer1 = tof_cluster_layer_d[0];
    tof_cluster_layer2 = tof_cluster_layer_d[1];
    tof_cluster_layer3 = tof_cluster_layer_d[2];
    tof_cluster_layer4 = tof_cluster_layer_d[3];
    tof_cluster_Aedep1 = tof_cluster_Aedep_d[0];
    tof_cluster_Aedep2 = tof_cluster_Aedep_d[1];
    tof_cluster_Aedep3 = tof_cluster_Aedep_d[2];
    tof_cluster_Aedep4 = tof_cluster_Aedep_d[3];
    tof_cluster_Dedep1 = tof_cluster_Dedep_d[0];
    tof_cluster_Dedep2 = tof_cluster_Dedep_d[1];
    tof_cluster_Dedep3 = tof_cluster_Dedep_d[2];
    tof_cluster_Dedep4 = tof_cluster_Dedep_d[3];
    tof_cluster_hit_time1 = tof_cluster_hit_time_d[0];
    tof_cluster_hit_time2 = tof_cluster_hit_time_d[1];
    tof_cluster_hit_time3 = tof_cluster_hit_time_d[2];
    tof_cluster_hit_time4 = tof_cluster_hit_time_d[3];
    tof_cluster_hit_time_err1 = tof_cluster_hit_time_err_d[0];
    tof_cluster_hit_time_err2 = tof_cluster_hit_time_err_d[1];
    tof_cluster_hit_time_err3 = tof_cluster_hit_time_err_d[2];
    tof_cluster_hit_time_err4 = tof_cluster_hit_time_err_d[3];
    tof_cluster_non_correct_time11 = tof_cluster_non_correct_time_d[0][0];
    tof_cluster_non_correct_time12 = tof_cluster_non_correct_time_d[0][1];
    tof_cluster_non_correct_time13 = tof_cluster_non_correct_time_d[0][2];
    tof_cluster_non_correct_time14 = tof_cluster_non_correct_time_d[0][3];
    tof_cluster_non_correct_time21 = tof_cluster_non_correct_time_d[1][0];
    tof_cluster_non_correct_time22 = tof_cluster_non_correct_time_d[1][1];
    tof_cluster_non_correct_time23 = tof_cluster_non_correct_time_d[1][2];
    tof_cluster_non_correct_time24 = tof_cluster_non_correct_time_d[1][3];
    tof_cluster_correct_time11 = tof_cluster_correct_time_d[0][0];
    tof_cluster_correct_time12 = tof_cluster_correct_time_d[0][1];
    tof_cluster_correct_time13 = tof_cluster_correct_time_d[0][2];
    tof_cluster_correct_time14 = tof_cluster_correct_time_d[0][3];
    tof_cluster_correct_time21 = tof_cluster_correct_time_d[1][0];
    tof_cluster_correct_time22 = tof_cluster_correct_time_d[1][1];
    tof_cluster_correct_time23 = tof_cluster_correct_time_d[1][2];
    tof_cluster_correct_time24 = tof_cluster_correct_time_d[1][3];
    tof_cluster_Acharge11 = tof_cluster_Acharge_d[0][0];
    tof_cluster_Acharge12 = tof_cluster_Acharge_d[0][1];
    tof_cluster_Acharge13 = tof_cluster_Acharge_d[0][2];
    tof_cluster_Acharge14 = tof_cluster_Acharge_d[0][3];
    tof_cluster_Acharge21 = tof_cluster_Acharge_d[1][0];
    tof_cluster_Acharge22 = tof_cluster_Acharge_d[1][1];
    tof_cluster_Acharge23 = tof_cluster_Acharge_d[1][2];
    tof_cluster_Acharge24 = tof_cluster_Acharge_d[1][3];
    tof_cluster_Dcharge111 = tof_cluster_DCharge_d[0][0][0];
    tof_cluster_Dcharge112 = tof_cluster_DCharge_d[0][0][1];
    tof_cluster_Dcharge113 = tof_cluster_DCharge_d[0][0][2];
    tof_cluster_Dcharge114 = tof_cluster_DCharge_d[0][0][3];
    tof_cluster_Dcharge121 = tof_cluster_DCharge_d[0][1][0];
    tof_cluster_Dcharge122 = tof_cluster_DCharge_d[0][1][1];
    tof_cluster_Dcharge123 = tof_cluster_DCharge_d[0][1][2];
    tof_cluster_Dcharge124 = tof_cluster_DCharge_d[0][1][3];
    tof_cluster_Dcharge131 = tof_cluster_DCharge_d[0][2][0];
    tof_cluster_Dcharge132 = tof_cluster_DCharge_d[0][2][1];
    tof_cluster_Dcharge133 = tof_cluster_DCharge_d[0][2][2];
    tof_cluster_Dcharge134 = tof_cluster_DCharge_d[0][2][3];
    tof_cluster_Dcharge211 = tof_cluster_DCharge_d[1][0][0];
    tof_cluster_Dcharge212 = tof_cluster_DCharge_d[1][0][1];
    tof_cluster_Dcharge213 = tof_cluster_DCharge_d[1][0][2];
    tof_cluster_Dcharge214 = tof_cluster_DCharge_d[1][0][3];
    tof_cluster_Dcharge221 = tof_cluster_DCharge_d[1][1][0];
    tof_cluster_Dcharge222 = tof_cluster_DCharge_d[1][1][1];
    tof_cluster_Dcharge223 = tof_cluster_DCharge_d[1][1][2];
    tof_cluster_Dcharge224 = tof_cluster_DCharge_d[1][1][3];
    tof_cluster_Dcharge231 = tof_cluster_DCharge_d[1][2][0];
    tof_cluster_Dcharge232 = tof_cluster_DCharge_d[1][2][1];
    tof_cluster_Dcharge233 = tof_cluster_DCharge_d[1][2][2];
    tof_cluster_Dcharge234 = tof_cluster_DCharge_d[1][2][3];
    rich_used_hits = rich_used_hits_d;
    rich_photoelectrons = rich_photoelectrons_d;
    rich_expect_photoelectrons = rich_expect_photoelectrons_d;
    rich_prob = rich_prob_d;
    rich_track_theta = rich_track_theta_d;
    rich_track_phi = rich_track_phi_d;
    rich_Nhits_ring = rich_Nhits_ring_d;
    rich_Npmts_ring = rich_Npmts_ring_d;
    charge_tofL1 = charge_tofL_d[0];
    charge_tofL2 = charge_tofL_d[1];
    charge_tofL3 = charge_tofL_d[2];
    charge_tofL4 = charge_tofL_d[3];
    trk_chi2_max_X = trk_chi2_max_X_d;
    trk_chi2_max_Y = trk_chi2_max_Y_d;
    trk_chi2_inner_Y = trk_chi2_inner_Y_d;
    trk_chi2_inner_X = trk_chi2_inner_X_d;
    diff_trkchi2Y_inner_max = TMath::Abs( trk_chi2_inner_Y_d - trk_chi2_max_Y_d)/(trk_chi2_inner_Y_d + trk_chi2_max_Y_d);

    mc_d_tree->Fill();
  }

  output->Write();

  stopwatch->Stop();
  stopwatch->Print();

}





