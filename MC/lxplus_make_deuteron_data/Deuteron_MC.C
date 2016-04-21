#include "root.h"
#include "amschain.h"
#include <root_setup.h>
#include <TrTrack.h>
#include <Tofcharge_ihep.h>
#include <bcorr.h>
#include <TrExtAlignDB.h>
#include <TrRecon.h> // for TRCLFFKEY.UseSensolAlign = 0
#include "FileManager.hh"
#include "AMSRootSupport.hh"
#include "AMSRootParticleFactory.hh"
#include <TrdHit.hh>
#include <TrRecHit.h>
#include <AnalysisEvent.hh>
#include <AnalysisParticle.hh>
#include <TrdKCluster.h>

Int_t nlayer;
Float_t qrms;

    HeaderR* header = NULL;
    Level1R* trig = NULL;
    DaqEventR* daq = NULL;
    TrTrackR* trktrack = NULL;
    BetaHR* betaHR = NULL;
    ParticleR* part = NULL;
    AMSEventR* pev = NULL;
    RichRingR* richring = NULL;
    EcalShowerR* ecalshower = NULL;
    TofChargeHR tofchargeH;
    AMSSetupR::RTI rti;
int id_innerspan;
int id_maxspan;
int id_uphalf_span;
int id_dwnhalf_span;
int id_upinhalf_span;
int id_dwninhalf_span;
int hitACCpatt[8]={-1,-1,-1,-1,-1,-1,-1,-1};
int hitACC = -1;
int physBpatt[8] = {-1,-1,-1,-1,-1,-1,-1,-1};

bool CheckBadRun(AMSEventR*);
bool CheckScienceRun(AMSEventR*);
bool GoodHW(DaqEventR*);
bool GoodTrkTrack(TrTrackR*);

typedef struct segment {
  float x;
  float y;
  float z;
} trd_segment;

int main(int argc, char* argv[])
{
  AMSChain* ch = new AMSChain;

  char filelist[128];
  char oputfile[128];
  TString logfilename;

  if(argc == 1)
  {
    cout<<"*** Program is running in single root file mode ***"<<endl;
    //ch->Add("root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6//1363740203.00000001.root");
    //ch->Add("root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6//1411948600.00000001.root");
    //ch->Add("root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6//1364040601.00000001.root");
    ch->Add("root://eosams.cern.ch//eos/ams/MC/AMS02/2014/d.B1030/d.pl1.0_520_GG_Blic/872728226.00000001.root");
    //ch->Add("root://eosams.cern.ch//eos/ams/MC/AMS02/2014/aHe.B1052/ahe3.pl1.21000.112/1208830902.00000001.root");
    //ch->Add("/afs/cern.ch/work/s/sikang/background_study/Deuteron_MC/aachen/link_d_aachen1/664500000.00000001.root");
    //ch->Add("/afs/cern.ch/work/s/sikang/background_study/Deuteron_MC/aachen/664505320.00000001.root");

    TString oputfile1 = "test_Deuteron_MC.root";
    //TString oputfile1 = "test_MC.root";
    strcpy(oputfile, oputfile1);
  }
  else if( argc >= 2 )
  {
    cout<<Form("usage: %s \n",argv[0]);
    strcpy(filelist, argv[1]);
    strcpy(oputfile, argv[2]);

    string file_candidate;
    ifstream file(filelist);
    if(file.is_open())
    {
      while(!file.eof())
      {
        getline( file, file_candidate);
        if(file_candidate.length() == 0) continue;
        ch->Add( file_candidate.c_str() );
        cout<<Form("added '%s' to the ams", file_candidate.c_str()) << endl;
      }
      file.close();
    }
    else
    {
      cout<<Form("Unable to open the filelist %s\n\n exit", filelist) << endl;
      exit(1);
    }
    logfilename = Form("%s_out.txt", oputfile);

  }

  bool cut_tof_chi2;
  bool cut_trk_charge;
  bool cut_trk_yhit5;
  bool cut_trk_yhit4;
  bool cut_trk_chi2;
  bool cut_tof_charge;
  bool cut_trk_hitL1XY;
  bool cut_trk_hitL9XY;
  bool cut_trk_hitL1;
  bool cut_trk_hitL9;
  bool cut_tof_cluster;
  bool cut_tof_beta;

  ULong64_t num = ch->GetEntries();
  ULong64_t x = num/100;
  Int_t y=0;
  Float_t beta;
  Float_t inv_beta_err;
  //Float_t inv_beta_norm_err;
  Double_t rigidity_inner;
  Double_t rigidity_max;
  Double_t inv_rigidity_inner_err;
  Double_t inv_rigidity_max_err;
  Float_t E_ecal;
  Float_t corrected_beta;
  Float_t charge_trk;
  Float_t charge_tof;
  Float_t charge_tof_up;
  Float_t charge_tof_low;
  Float_t charge_tofL[4];
  Float_t livetime;
  Float_t zenith_angle;
  Float_t Rcut[4][2];
  Double_t red_chi2x_inner;
  Double_t red_chi2y_inner;
  Double_t red_chi2x_max;
  Double_t red_chi2y_max;
  Double_t red_chi2x_uphalf;
  Double_t red_chi2y_uphalf;
  Double_t red_chi2y_dwnhalf;
  Double_t red_chi2x_dwnhalf;
  Float_t phystrig;
  Float_t unbtrig;
  Float_t rich_beta;
  Float_t rich_beta_err;
  Int_t rich_isnaf;
  Float_t ecal_bdt;
  Int_t rich_passcut;
  Float_t trk_edep[9];
  Float_t tof_edep[4];
  Float_t init_momentum;
  Int_t particle_ID;
  Int_t primary_ID;
  Int_t trd_mc_cut;
  Int_t tof_mc_cut;
  Int_t trk_mc_cut;
  Int_t rich_mc_cut;

  Float_t mc_prob[2];
  Float_t init_coo[3];
  Float_t init_dir[3];
  Float_t init_mass;
  Float_t init_charge;
  Int_t run_Num;
  Int_t evt_Num;
  Int_t trd_vertex3;
  Int_t trd_vertex2;
  Int_t ntrdtrack;
  Float_t tofcoo[4][3];
  Float_t trkcoo[9][3];
  Float_t richcoo[2][3];
  Float_t ecalcoo[3][3];
  Float_t tofdir[4][2];
  Float_t trkdir[9][2];
  Float_t ecaldir[3][2];
  Float_t richdir[2][2];
  Float_t trk_tof[4][2];
  Float_t trk_trk[9][2];
  Float_t trk_rich[2][2];
  Float_t trk_ecal[3][2];
  Double_t rigidity_up_half;
  Double_t rigidity_dwn_half;
  Int_t ntrtrack;
  Int_t ntrdhsegment;
  Int_t ntrdsegment;
  Float_t ecal_Edep[18];
  Int_t nacchits;
  vector<Float_t> acc_time_diff;
  Float_t acc_time_c1, acc_time_c2;
  Int_t ntofhits;
  vector<Float_t> tof_time;
  Int_t ntrkhits;
  Int_t ntrkxhits;
  Int_t ntrkxyhits;
  Int_t ntrkyhits;
  Int_t trkhit_plane[9];
  vector<Float_t> acc_edep;
  Double_t rigidity_dwn_inhalf;
  Double_t rigidity_up_inhalf;
  trd_segment trdsegment;
  vector<Float_t> evt_segment_coo_x;
  vector<Float_t> evt_segment_coo_y;
  vector<Float_t> evt_segment_coo_z;
  vector<Float_t> part_segment_coo_x;
  vector<Float_t> part_segment_coo_y;
  vector<Float_t> part_segment_coo_z;
  vector<Int_t> evt_ntrd_cluster;
  vector<Int_t> part_ntrd_cluster;
  vector<Int_t> evt_segment_direction;
  vector<Int_t> part_segment_direction;
  Float_t rich_charge2;
  Double_t trdk_e2plikelihood;
  Double_t trdk_he2plikelihood;
  Double_t trdp_e2plikelihood;
  Double_t trdp_he2plikelihood;
  Double_t trk_theta;
  Double_t trk_phi;
  Int_t trk_pattern_X;
  Int_t trk_pattern_Y;
  Int_t trk_pattern_XY;
  Int_t Ntof_cluster;
  Int_t Ntof_cluster_intime;
  Int_t beta_pattern;
  Float_t tof_chi2_time;
  Float_t tof_chi2_coord;
  Float_t tof_mass;
  Float_t tof_mass_err;
  Int_t tof_used_hit;
  Int_t tof_sum_hit;
  Int_t tof_cluster_pattern[4];
  Int_t tof_cluster_layer[4];
  Float_t tof_cluster_Aedep[4];
  Float_t tof_cluster_Dedep[4];
  Float_t tof_cluster_hit_time[4];
  Float_t tof_cluster_hit_time_err[4];
  Float_t tof_cluster_non_correct_time[2][4];
  Float_t tof_cluster_correct_time[2][4];
  Float_t tof_cluster_Acharge[2][4];
  Float_t tof_cluster_Dcharge[2][3][4];
  Int_t rich_used_hits;
  Float_t rich_photoelectrons;
  Float_t rich_expect_photoelectrons;
  Float_t rich_prob;
  Float_t rich_track_theta;
  Float_t rich_track_phi;
  Int_t rich_Nhits_ring;
  Int_t rich_Npmts_ring;
  Int_t rich_isgood;
  Int_t rich_isclean;
  

  TFile* f = new TFile(oputfile, "recreate");
  TTree* t = new TTree("t", "deuteron_selection");
  t->Branch("cut_tof_chi2", &cut_tof_chi2, "cut_tof_chi2/O");
  t->Branch("cut_trk_charge", &cut_trk_charge, "cut_trk_charge/O");
  t->Branch("cut_trk_yhit5", &cut_trk_yhit5, "cut_trk_yhit5/O");
  t->Branch("cut_trk_yhit4", &cut_trk_yhit4, "cut_trk_yhit4/O");
  t->Branch("cut_trk_chi2", &cut_trk_chi2, "cut_trk_chi2/O");
  t->Branch("cut_tof_charge", &cut_tof_charge, "cut_tof_charge/O");
  t->Branch("cut_trk_hitL1XY", &cut_trk_hitL1XY, "cut_trk_hitL1XY/O");
  t->Branch("cut_trk_hitL1", &cut_trk_hitL1, "cut_trk_hitL1/O");
  t->Branch("cut_trk_hitL9XY", &cut_trk_hitL9XY, "cut_trk_hitL9XY/O");
  t->Branch("cut_trk_hitL9", &cut_trk_hitL9, "cut_trk_hitL9/O");
  t->Branch("cut_tof_cluster", &cut_tof_cluster, "cut_tof_cluster/O");
  t->Branch("cut_tof_beta", &cut_tof_beta, "cut_tof_beta/O");
  t->Branch("run_Num", &run_Num, "run_Num/I");
  t->Branch("evt_Num", &evt_Num, "evt_Num/I");
  t->Branch("particle_ID", &particle_ID, "particle_ID/I");
  t->Branch("primary_ID", &primary_ID, "primary_ID/I");
  t->Branch("trk_mc_cut", &trk_mc_cut, "trk_mc_cut/I");
  t->Branch("rich_mc_cut", &rich_mc_cut, "rich_mc_cut/I");
  t->Branch("trd_mc_cut", &trd_mc_cut, "trd_mc_cut/I");
  t->Branch("tof_mc_cut", &tof_mc_cut, "tof_mc_cut/I");
  t->Branch("rich_isnaf", &rich_isnaf, "rich_isnaf/I");
  t->Branch("rich_passcut", &rich_passcut, "rich_passcut/I");
  t->Branch("trd_vertex3", &trd_vertex3, "trd_vertex3/I");
  t->Branch("trd_vertex2", &trd_vertex2, "trd_vertex2/I");
  t->Branch("ntrdtrack", &ntrdtrack, "ntrdtrack/I");
  t->Branch("ntrtrack", &ntrtrack, "ntrtrack/I");
  t->Branch("ntrdsegment", &ntrdsegment, "ntrdsegment/I");
  t->Branch("ntrdhsegment", &ntrdhsegment, "ntrdhsegment/I");
  t->Branch("nacchits", &nacchits, "nacchits/I");
  t->Branch("ntofhits", &ntofhits, "ntofhits/I");
  t->Branch("ntrkhits", &ntrkhits, "ntrkhits/I");
  t->Branch("ntrkxhits", &ntrkxhits, "ntrkxhits/I");
  t->Branch("ntrkyhits", &ntrkyhits, "ntrkyhits/I");
  t->Branch("ntrkxyhits", &ntrkxyhits, "ntrkxyhits/I");
  t->Branch("trkhit_plane", &trkhit_plane, "trkhit_plane[9]/I");
  t->Branch("trk_pattern_X", &trk_pattern_X, "trk_pattern_X/I");
  t->Branch("trk_pattern_Y", &trk_pattern_Y, "trk_pattern_Y/I");
  t->Branch("trk_pattern_XY", &trk_pattern_XY, "trk_pattern_XY/I");
  t->Branch("Ntof_cluster", &Ntof_cluster, "Ntof_cluster/I");
  t->Branch("Ntof_cluster_intime", &Ntof_cluster_intime, "Ntof_cluster_intime/I");
  t->Branch("beta_pattern", &beta_pattern, "beta_pattern/I");
  t->Branch("tof_used_hit", &tof_used_hit, "tof_used_hit/I");
  t->Branch("tof_sum_hit", &tof_sum_hit, "tof_sum_hit/I");
  t->Branch("tof_cluster_pattern", &tof_cluster_pattern, "tof_cluster_pattern[4]/I");
  t->Branch("tof_cluster_layer", &tof_cluster_layer, "tof_cluster_layer[4]/I");
  t->Branch("rich_used_hits", &rich_used_hits, "rich_used_hits/I");
  t->Branch("rich_Nhits_ring", &rich_Nhits_ring, "rich_Nhits_ring/I");
  t->Branch("rich_Npmts_ring", &rich_Npmts_ring, "rich_Npmts_ring/I");
  t->Branch("evt_ntrd_cluster", &evt_ntrd_cluster);
  t->Branch("part_ntrd_cluster", &part_ntrd_cluster);
  t->Branch("evt_segment_direction", &evt_segment_direction);
  t->Branch("part_segment_direction", &part_segment_direction);
  t->Branch("init_momentum", &init_momentum, "init_momentum/F");
  t->Branch("init_coo", &init_coo, "init_coo[3]/F");
  t->Branch("init_dir", &init_dir, "init_dir[3]/F");
  t->Branch("init_mass", &init_mass, "init_mass/F");
  t->Branch("init_charge", &init_charge, "init_charge/F");
  t->Branch("mc_prob", &mc_prob, "mc_prob[2]/F");
  t->Branch("charge_trk", &charge_trk, "charge_trk/F");
  t->Branch("charge_tof", &charge_tof, "charge_tof/F");
  t->Branch("charge_tof_up", &charge_tof_up, "charge_tof_up/F");
  t->Branch("charge_tof_low", &charge_tof_low, "charge_tof_low/F");
  t->Branch("charge_tofL", &charge_tofL, "charge_tofL[4]/F");
  t->Branch("beta", &beta, "beta/F");
  t->Branch("inv_beta_err", &inv_beta_err, "inv_beta_err/F");
  //t->Branch("inv_beta_norm_err", &inv_beta_norm_err, "inv_beta_norm_err/F");
  t->Branch("rich_beta", &rich_beta, "rich_beta/F");
  t->Branch("rich_beta_err", &rich_beta_err, "rich_beta_err/F");
  t->Branch("ecal_bdt", &ecal_bdt, "ecal_bdt/F");
  t->Branch("E_ecal", &E_ecal, "E_ecal/F");
  t->Branch("trk_edep", &trk_edep, "trk_edep[9]/F");
  t->Branch("tof_edep", &tof_edep, "tof_edep[4]/F");
  t->Branch("corrected_beta", &corrected_beta, "corrected_beta/F");
  t->Branch("tofcoo", &tofcoo, "tofcoo[4][3]/F");
  t->Branch("trkcoo", &trkcoo, "trkcoo[9][3]/F");
  t->Branch("richcoo", &richcoo, "richcoo[2][3]/F");
  t->Branch("ecalcoo", &ecalcoo, "ecalcoo[3][3]/F");
  t->Branch("tofdir", &tofdir, "tofdir[4][2]/F");
  t->Branch("trkdir", &trkdir, "trkdir[9][2]/F");
  t->Branch("richdir", &richdir, "richdir[2][2]/F");
  t->Branch("rich_charge2", &rich_charge2, "rich_charge2/F");
  t->Branch("ecaldir", &ecaldir, "ecaldir[3][2]/F");
  t->Branch("trk_tof", &trk_tof, "trk_tof[4][2]/F");
  t->Branch("trk_trk", &trk_trk, "trk_trk[9][2]/F");
  t->Branch("trk_rich", &trk_rich, "trk_rich[3][2]/F");
  t->Branch("trk_ecal", &trk_ecal, "trk_ecal[3][2]/F");
  t->Branch("ecal_Edep", &ecal_Edep, "ecal_Edep[18]/F");
  t->Branch("tof_chi2_time", &tof_chi2_time, "tof_chi2_time/F");
  t->Branch("tof_chi2_coord", &tof_chi2_coord, "tof_chi2_coord/F");
  t->Branch("tof_mass", &tof_mass, "tof_mass/F");
  t->Branch("tof_mass_err", &tof_mass_err, "tof_mass_err/F");
  t->Branch("tof_cluster_Aedep", &tof_cluster_Aedep, "tof_cluster_Aedep[4]/F");
  t->Branch("tof_cluster_Dedep", &tof_cluster_Dedep, "tof_cluster_Dedep[4]/F");
  t->Branch("tof_cluster_hit_time", &tof_cluster_hit_time, "tof_cluster_hit_time[4]/F");
  t->Branch("tof_cluster_hit_time_err", &tof_cluster_hit_time_err, "tof_cluster_hit_time_err[4]/F");
  t->Branch("tof_cluster_non_correct_time", &tof_cluster_non_correct_time, "tof_cluster_non_correct_time[2][4]/F");
  t->Branch("tof_cluster_correct_time", &tof_cluster_correct_time, "tof_cluster_correct_time[2][4]/F");
  t->Branch("tof_cluster_Acharge", &tof_cluster_Acharge, "tof_cluster_Acharge[2][4]/F");
  t->Branch("tof_cluster_DCharge", &tof_cluster_Dcharge, "tof_cluster_Dcharge[2][3][4]/F");
  t->Branch("rich_photoelectrons", &rich_photoelectrons, "rich_photoelectrons/F");
  t->Branch("rich_expect_photoelectrons", &rich_expect_photoelectrons, "rich_expect_photoelectrons/F");
  t->Branch("rich_prob", &rich_prob, "rich_prob/F");
  t->Branch("rich_track_theta", &rich_track_theta, "rich_track_theta/F");
  t->Branch("rich_track_phi", &rich_track_phi, "rich_track_phi/F");
  t->Branch("rich_isgood", &rich_isgood, "rich_isgood/I");
  t->Branch("rich_isclean", &rich_isclean, "rich_isclean/I");
  t->Branch("acc_edep", &acc_edep);
  t->Branch("acc_time_diff", &acc_time_diff);
  t->Branch("tof_time", &tof_time);
  t->Branch("evt_segment_coo_x", &evt_segment_coo_x);
  t->Branch("evt_segment_coo_y", &evt_segment_coo_y);
  t->Branch("evt_segment_coo_z", &evt_segment_coo_z);
  t->Branch("part_segment_coo_x", &part_segment_coo_x);
  t->Branch("part_segment_coo_y", &part_segment_coo_y);
  t->Branch("part_segment_coo_z", &part_segment_coo_z);
  t->Branch("rigidity_inner", &rigidity_inner, "rigidity_inner/D");
  t->Branch("rigidity_max", &rigidity_max, "rigidity_max/D");
  t->Branch("inv_rigidity_inner_err", &inv_rigidity_inner_err, "inv_rigidity_inner_err/D");
  t->Branch("inv_rigidity_max_err", &inv_rigidity_max_err, "inv_rigidity_max_err/D");
  t->Branch("red_chi2x_inner", &red_chi2x_inner, "red_chi2x_inner/D");
  t->Branch("red_chi2y_inner", &red_chi2y_inner, "red_chi2y_inner/D");
  t->Branch("red_chi2x_max", &red_chi2x_max, "red_chi2x_max/D");
  t->Branch("red_chi2y_max", &red_chi2y_max, "red_chi2y_max/D");
  t->Branch("rigidity_up_half", &rigidity_up_half, "rigidity_up_half/D");
  t->Branch("rigidity_up_inhalf", &rigidity_up_inhalf, "rigidity_up_inhalf/D");
  t->Branch("rigidity_dwn_half", &rigidity_dwn_half, "rigidity_dwn_half/D");
  t->Branch("rigidity_dwn_inhalf", &rigidity_dwn_inhalf, "rigidity_dwn_inhalf/D");
  t->Branch("red_chi2x_uphalf", &red_chi2x_uphalf, "red_chi2x_uphalf/D");
  t->Branch("red_chi2y_uphalf", &red_chi2y_uphalf, "red_chi2y_uphalf/D");
  t->Branch("red_chi2x_dwnhalf", &red_chi2x_dwnhalf, "red_chi2x_dwnhalf/D");
  t->Branch("red_chi2y_dwnhalf", &red_chi2y_dwnhalf, "red_chi2y_dwnhalf/D");
  t->Branch("trdk_e2plikelihood", &trdk_e2plikelihood, "trdk_e2plikelihood/D");
  t->Branch("trdk_he2plikelihood", &trdk_he2plikelihood, "trdk_he2plikelihood/D");
  t->Branch("trdp_e2plikelihood", &trdp_e2plikelihood, "trdp_e2plikelihood/D");
  t->Branch("trdp_he2plikelihood", &trdp_he2plikelihood, "trdp_he2plikelihood/D");
  t->Branch("trk_theta", &trk_theta, "trk_theta/D");
  t->Branch("trk_phi", &trk_phi, "trk_phi/D");

  

  TH1D* hEvtCounter = new TH1D("hEvtcounter", "Event Counter", 17, 1,18);
  hEvtCounter->GetXaxis()->SetBinLabel(1, "Total Entries");
  hEvtCounter->GetXaxis()->SetBinLabel(2, "get event");
  hEvtCounter->GetXaxis()->SetBinLabel(3, "pass trigger");
  hEvtCounter->GetXaxis()->SetBinLabel(4, "header");
  hEvtCounter->GetXaxis()->SetBinLabel(5, "gbatch err");
  hEvtCounter->GetXaxis()->SetBinLabel(6, "bad run");
  hEvtCounter->GetXaxis()->SetBinLabel(7, "science run");
  hEvtCounter->GetXaxis()->SetBinLabel(8, "good HW");
  hEvtCounter->GetXaxis()->SetBinLabel(9, "In SAA");
  hEvtCounter->GetXaxis()->SetBinLabel(10, "particle point");
  hEvtCounter->GetXaxis()->SetBinLabel(11, "RTI cut");
  hEvtCounter->GetXaxis()->SetBinLabel(12, "tracker position");
  hEvtCounter->GetXaxis()->SetBinLabel(13, "tracker hit L1L9");
  hEvtCounter->GetXaxis()->SetBinLabel(14, "tof charge");
  hEvtCounter->GetXaxis()->SetBinLabel(15, "tracker charge");
  hEvtCounter->GetXaxis()->SetBinLabel(16, "tof cluster hit");
  hEvtCounter->GetXaxis()->SetBinLabel(17, "chi2 cut");
  
  TH1D* hEvtCounter2 = new TH1D("hEvtcounter2", "Event Counter2", 11, 1, 12);
  hEvtCounter2->GetXaxis()->SetBinLabel(1, "Total Entries");
  hEvtCounter2->GetXaxis()->SetBinLabel(2, "get event");
  hEvtCounter2->GetXaxis()->SetBinLabel(3, "one particle");
  hEvtCounter2->GetXaxis()->SetBinLabel(4, "tof chi2");
  hEvtCounter2->GetXaxis()->SetBinLabel(5, "trk charge");
  hEvtCounter2->GetXaxis()->SetBinLabel(6, "trk yhit");
  hEvtCounter2->GetXaxis()->SetBinLabel(7, "tof cluster");
  hEvtCounter2->GetXaxis()->SetBinLabel(8, "tof cluster-layer");
  hEvtCounter2->GetXaxis()->SetBinLabel(9, "trk chi2");
  hEvtCounter2->GetXaxis()->SetBinLabel(10, "tof charge");
  hEvtCounter2->GetXaxis()->SetBinLabel(11, "filled event");
  


  AMSSetupR::RTI::UseLatest(6);
  TkDBc::UseFinal();


  const bool setAMSRootDefaults = true;

  AMSRootSupport amsRootSupport(AC::MCRun, setAMSRootDefaults);
  Analysis::EventFactory& eventFactory = amsRootSupport.EventFactory();
  Analysis::AMSRootParticleFactory& particleFactory = amsRootSupport.ParticleFactory();

  TH1F* tofcluster_layer = new TH1F("tofcluster", "tofcluster", 4,1,5);
  TH1F* tofcluster[7];
  for(Int_t i=0 ; i<7 ; i++)
  {
    if(i<5) tofcluster[i] = new TH1F(Form("tofcluster_%d",i), Form("tofcluster_%d",i), 15, 0, 15);
    if(i==5) tofcluster[i] = new TH1F("upintime_lowintime", "upintime_lowintime", 15, 0, 15);
    if(i==6) tofcluster[i] = new TH1F("mean_upintime_lowintime", "mean_upintime_lowintime", 15, 0, 15);
  }
  TH1F* tofcluster_part = new TH1F("tofcluster_part", "tofcluster_part", 4, 1, 5);
  TH1F* tofcluster_intime = new TH1F("tofcluster_intime", "tof cluster intime ", 4,1,5);
  for(ULong64_t entry=0 ; entry<num ; entry++)
  //for(ULong64_t entry=0 ; entry<1000 ; entry++)
  {
    if(entry%(int)x==0)
    {
      cout<<y<<"%"<<endl;
      y++;
    }
    
    HeaderR* header = NULL;
    Level1R* trig = NULL;
    DaqEventR* daq = NULL;
    TrTrackR* trktrack = NULL;
    BetaHR* betaHR = NULL;
    ParticleR* part = NULL;
    AMSEventR* pev = NULL;
    RichRingR* richring = NULL;
    EcalShowerR* ecalshower = NULL;
    TofChargeHR tofchargeH;
    AMSSetupR::RTI rti;
    TofClusterHR* tofclusterH = NULL;
    MCEventgR* mcevent = NULL;
    MCTrackR* mctrack = NULL;
    RichMCClusterR* mcrich = NULL;
    TofMCClusterR* mctof = NULL;
    TrdMCClusterR* mctrd = NULL;

    cut_tof_chi2 = true;
    cut_trk_charge = true;
    cut_trk_yhit5 = true;
    cut_trk_yhit4 = true;
    cut_trk_chi2 = true;
    cut_tof_charge = true;
    cut_trk_hitL1XY = true;
    cut_trk_hitL1 = true;
    cut_trk_hitL9XY = true;
    cut_trk_hitL9 = true;
    cut_tof_cluster = true;
    cut_tof_beta = true;

    hEvtCounter->Fill(1);
    hEvtCounter2->Fill(1);
  
    pev = ch->GetEvent(entry);
    if(!pev) continue;

    pev->SetDefaultMCTuningParameters();
    Int_t g4version = pev->G4Version();
    Int_t product_version = pev->ProductionVersion();
    Int_t version = pev->Version();
    //cout<<"g4version : "<<g4version<<"   product version : "<<product_version<<"   version : "<<version<<endl;

    hEvtCounter2->Fill(2);
    hEvtCounter->Fill(2);
    /*
    trig = pev->pLevel1(0);
    if(!trig) continue;

    hEvtCounter->Fill(3);
    
    header = &(pev->fHeader);
    if(!header) continue;

    hEvtCounter->Fill(4);

    bool gbatcherr = pev->Status(30);
    if(gbatcherr) continue;

    hEvtCounter->Fill(5);

    if(CheckBadRun(pev)) continue;

    hEvtCounter->Fill(6);

    if(!CheckScienceRun(pev)) continue;

    hEvtCounter->Fill(7);

    daq = pev->pDaqEvent(0);
    if(!GoodHW(daq)) continue;

    hEvtCounter->Fill(8);


    int time = rti.utime;
    if(pev->IsInSAA(time)) continue;

    hEvtCounter->Fill(9);
*/
    int runNum = pev->Run();
    //if(runNum == 1306219312 || runNum == 1306219522 || runNum == 1306233745 || (runNum >= 1307125541 && runNum <= 1307218054)) continue;
    
    mcevent = pev->pMCEventg(0);
    

    if(mcevent) TrExtAlignDB::SmearExtAlign();

    if(mcevent) TRCLFFKEY.UseSensorAlign = 0;

    run_Num = pev->Run();
    if(!run_Num) continue;
    evt_Num = pev->Event();
    if(!evt_Num) continue;

    int nparticle = 0;
    nparticle = pev->NParticle();
    
    int goodpart_index = 0;
    int nGoodPart=0;
    for(int ipart=0 ; ipart<nparticle ; ipart++)
    {
      part = pev->pParticle(ipart);
      if(!part) continue;

      betaHR = part->pBetaH();
      if(!betaHR) continue;

      trktrack = part->pTrTrack();
      if(!trktrack) continue;

      if(!GoodTrkTrack(trktrack)) continue;

      nGoodPart++;

      goodpart_index = ipart;
    }


    if(nGoodPart==0) continue;
    if(nGoodPart>1) continue;
    hEvtCounter2->Fill(3);

    hEvtCounter->Fill(10);

    part = pev->pParticle(goodpart_index);
    betaHR = part->pBetaH();
    trktrack = part->pTrTrack();
    id_innerspan = trktrack->iTrTrackPar(1,3,20); // iTrTrackPar(1, 0(max), refit, mass, charge, beta, fixrig)
    id_maxspan = trktrack->iTrTrackPar(1,0,20);
    id_uphalf_span = trktrack->iTrTrackPar(1,9999,22);
    id_dwnhalf_span = trktrack->iTrTrackPar(1,999990000,22);
    id_upinhalf_span = trktrack->iTrTrackPar(1,1,22);
    id_dwninhalf_span = trktrack->iTrTrackPar(1,2,22);
    //rigidity_inner = trktrack->GetRigidity(id_innerspan);
    rigidity_max = trktrack->GetRigidity(id_maxspan);
    inv_rigidity_max_err = trktrack->GetErrRinv(id_maxspan);
    rigidity_inner = -9999;
    rigidity_up_half = -9999;
    rigidity_dwn_half = -9999;
    rigidity_up_inhalf = -9999;
    rigidity_dwn_inhalf = -9999;
    if(id_innerspan >= 0) 
    {
      rigidity_inner = trktrack->GetRigidity(id_innerspan);
      inv_rigidity_inner_err = trktrack->GetErrRinv(id_innerspan);
    }
    if(id_uphalf_span >=0) rigidity_up_half = trktrack->GetRigidity(id_uphalf_span);
    if(id_dwnhalf_span >=0) rigidity_dwn_half = trktrack->GetRigidity(id_dwnhalf_span);
    if(id_upinhalf_span >=0) rigidity_up_inhalf = trktrack->GetRigidity(id_upinhalf_span);
    if(id_dwninhalf_span >=0) rigidity_dwn_inhalf = trktrack->GetRigidity(id_dwninhalf_span);
    //if(rigidity_dwn_half==-9999) cout<<"track id : "<<id_dwnhalf_span<<endl;

    Float_t bcor1, bcor2;
    if(!mcevent)
    {
      int bret1 = MagnetVarp::btempcor(bcor1, 0, 1);
      int bret2 = MagnetVarp::btempcor(bcor2, 0, 2);

      if(bret1==0 && bret2==0)
      {
        rigidity_max *= (bcor1+bcor2)/2;
        rigidity_inner *= (bcor1+bcor2)/2;
        rigidity_up_half *= (bcor1+bcor2)/2;
        rigidity_dwn_half *= (bcor1+bcor2)/2;
      }

      else if(bret1 != 0 && bret2 == 0)
      {
        rigidity_max *= bcor2;
        rigidity_inner *= bcor2;
        rigidity_up_half *= bcor2;
        rigidity_dwn_half *= bcor2;
      }
    }

    mcevent = pev->GetPrimaryMC();
    primary_ID = mcevent->Particle;
    init_momentum = mcevent->Momentum;
    init_mass = mcevent->Mass;
    init_charge = mcevent->Charge;
    particle_ID = part->Particle;
    mc_prob[0] = part->Prob[0];
    mc_prob[1] = part->Prob[1];
    trd_mc_cut=0;
    tof_mc_cut=0;
    trk_mc_cut=0;
    rich_mc_cut=0;
    for(Int_t i=0 ; i<pev->nTrdMCCluster() ; i++)
    {
      Int_t a = pev->pTrdMCCluster(i)->ParticleNo;
      if(a == primary_ID) trd_mc_cut=1;
    }

    for(Int_t i=0 ; i<pev->nTofMCCluster() ; i++)
    {
      Int_t a = pev->pTofMCCluster(i)->Particle;
      if(a == primary_ID) tof_mc_cut=1;
    }

    for(Int_t i=0 ; i<pev->nTrMCCluster() ; i++)
    {
      Int_t b = pev->pTrMCCluster(i)->IsPrimary();
      if(b==1) trk_mc_cut=1;
    }


    for(Int_t i=0 ; i<pev->nRichMCCluster() ; i++)
    {
      Int_t a = pev->pRichMCCluster(i)->Id;
      if(a == primary_ID) rich_mc_cut=1;
    }
    
    richring = part->pRichRing();
    rich_passcut = 0;

    if(richring) 
    {
      rich_beta = richring->getBeta();
      rich_beta_err = richring->ErrorBeta;
      rich_isnaf = richring->IsNaF();
      rich_charge2 = richring->getCharge2Estimate(1);
      if(richring->IsGood() || richring->IsClean())
      {
        if(richring->getProb()>0.2 && richring->getUsedHits()>2 && richring->getExpectedPhotoElectrons()>2 && richring->getPhotoElectrons()>(0.5*richring->getExpectedPhotoElectrons())) rich_passcut = 1;
        if(rich_isnaf && rich_beta>1) rich_passcut = 0;
        if(!rich_isnaf && rich_beta>1) rich_passcut = 0;
      }
    }


    ecalshower = part->pEcalShower();
    if(ecalshower) ecal_bdt = ecalshower->GetEcalBDT();
    if(!ecalshower) ecal_bdt = -999;
    if(ecalshower) E_ecal = ecalshower->EnergyE;

    Float_t tof_hit = betaHR->NTofClusterH();
    int rti_good=0;
    /*if(pev->GetRTI(rti)==0)
    {
      for(Int_t j=0 ; j<6 ; j++) if(rti.good&(1<<j)) rti_good++;
      livetime = rti.lf;
      for(Int_t j=0 ; j<4 ; j++)
      {
        for(Int_t k=0 ; k<2 ; k++)
        {
          Rcut[j][k]=rti.cfi[j][k];
        }
      }
      zenith_angle = rti.zenith;
    }
    if(rti_good!=0) continue;
    if(livetime<=0.5) continue;
    if(zenith_angle>=40) continue;
    if(rigidity_inner>0 && rigidity_inner<(1.2*Rcut[0][1])) continue;
    if(rigidity_inner<0 && rigidity_inner>(1.2*Rcut[0][0])) continue;

    hEvtCounter->Fill(11);

    AMSPoint pn1, pn9, pd1, pd9;
    pev->GetRTIdL1L9(0, pn1, pd1, pev->UTime(), 60);
    pev->GetRTIdL1L9(1, pn9, pd9, pev->UTime(), 60);
    if(pd1.y() > 35 || pd9.y() > 45) continue;

    hEvtCounter->Fill(12);
*/
    //if(rigidity<=0.8) continue;

    //betaHR = part->pBetaH();
    beta = betaHR->GetBeta();
    inv_beta_err = betaHR->GetEBetaV();
    //inv_beta_norm_err = betaHR->GetNormEBetaV();
    if(beta<0.3) cut_tof_beta = false; // tof beta cut
    corrected_beta = betaHR->GetBetaC();

    if(betaHR->GetNormChi2T() > 10) cut_tof_chi2 = false; // tof chi2(time) cut
    if(betaHR->GetNormChi2C() > 10) cut_tof_chi2 = false; // tof chi2(coordinate) cut

    hEvtCounter2->Fill(4);




    trktrack = part->pTrTrack();
    if(!(trktrack->TestHitLayerJHasXY(1))) cut_trk_hitL1XY = false; 
    if(!(trktrack->TestHitLayerJ(1))) cut_trk_hitL1 = false; 
    if(!(trktrack->TestHitLayerJHasXY(9))) cut_trk_hitL9XY = false;
    if(!(trktrack->TestHitLayerJ(9))) cut_trk_hitL9 = false;

    hEvtCounter->Fill(13);

    AMSPoint L9;
    L9 = trktrack->GetHitCooLJ(9);
    //if(L9.x()>=33) continue;
    

    
    hEvtCounter->Fill(14);

    if(trktrack->GetInnerQ()<=0.7 || trktrack->GetInnerQ()>=1.5) cut_trk_charge = false;
    hEvtCounter2->Fill(5);
    // full span
    //if(trktrack->GetLayerJQ(1,beta)<=0.6 || trktrack->GetLayerJQ(1,beta)>=1.9) continue;
    //if(trktrack->GetLayerJQ(9,beta)<=0.6 || trktrack->GetLayerJQ(9,beta)>=1.9) continue;

    if(part->pTrTrack()->GetNhitsY() < 5) cut_trk_yhit5 = false;
    if(part->pTrTrack()->GetNhitsY() < 4) cut_trk_yhit4 = false;
    hEvtCounter2->Fill(6);
    // tracker hit y coordinate >= 5 -> pass

    hEvtCounter->Fill(15);
    charge_trk = trktrack->GetQ(beta);
    tofchargeH = betaHR->gTofCharge();
    charge_tof = tofchargeH.GetQ(nlayer,qrms);
    for(Int_t i=0 ; i<4 ; i++)
    {
      charge_tofL[i] = tofchargeH.GetQL(i);
    }

    int ncls[4];
    //if(pev->GetNTofClustersInTime(betaHR, ncls)<4) cut_tof_cluster = false; // 1 cluster each layers
    //if(pev->GetNTofClustersInTime(betaHR, ncls)<3) continue; // 1 cluster each layers
    tofcluster[0]->Fill(ncls[0]);
    tofcluster[1]->Fill(ncls[1]);
    tofcluster[2]->Fill(ncls[2]);
    tofcluster[3]->Fill(ncls[3]);
    tofcluster[4]->Fill(pev->GetNTofClustersInTime(betaHR, ncls));
    tofcluster[5]->Fill(ncls[0]+ncls[2]);
    tofcluster[6]->Fill((ncls[0]+ncls[2])/2);

    //if(ncls[0]<1 || ncls[2]<1) continue;
    //if(ncls[0] + ncls[2] < 3) continue;
    hEvtCounter2->Fill(7);

    //cout<<"up intime : "<<ncls[0]<<"   up offtime : "<<ncls[1]<<"   low intime : "<<ncls[2]<<"   low offtime : "<<ncls[3]<<"   nncls : "<<ncls[4]<<endl;
  /* 
    Int_t times_layer=1;
    Int_t a=0;
    Int_t b[4]={0,0,0,0};
    //for(Int_t i=0 ; i<pev->GetNTofClustersInTime(betaHR, ncls) ; i++)
    for(Int_t i=0 ; i<ncls[0]+ncls[2] ; i++)
    {
      Int_t cluster_layer = pev->pTofClusterH(i)->Layer;
      if(cluster_layer == 0) b[0]=1;
      if(cluster_layer == 1) b[1]=1;
      if(cluster_layer == 2) b[2]=1;
      if(cluster_layer == 3) b[3]=1;
      //tofcluster_layer->Fill(cluster_layer);
      if(cluster_layer == 0) a=1;
      if(cluster_layer != 0) times_layer *= cluster_layer;
    }
    Int_t passed_layer=0;
    for(Int_t i=0 ; i<4 ; i++)
    {
      passed_layer += b[i];
    }
    tofcluster_layer->Fill(passed_layer);
    //if(times_layer % 6 != 0 || a==0) cut_tof_cluster = false; // at least one cluster exists in each layers
    if(passed_layer < 3) continue;//cut_tof_cluster = false;
   */ hEvtCounter2->Fill(8);

    hEvtCounter->Fill(16);

    //double ra = fabs(rigidity);
    //if(1/beta < sqrt(1+1.5*1.5/(ra+0.5)/(ra+0.5))-0.25) continue;

    if((trktrack->GetNormChisqY(id_innerspan)>=10)) cut_trk_chi2 = false;
    hEvtCounter2->Fill(9);

    hEvtCounter->Fill(17);

    //red_chi2x_inner = trktrack->GetNormChisqX(id_innerspan);
    //red_chi2y_inner = trktrack->GetNormChisqY(id_innerspan);
    red_chi2x_max = trktrack->GetNormChisqX(id_maxspan);
    red_chi2y_max = trktrack->GetNormChisqY(id_maxspan);
    red_chi2x_inner = -9999;
    red_chi2y_inner = -9999;
    red_chi2x_uphalf = -9999;
    red_chi2y_uphalf = -9999;
    red_chi2x_dwnhalf = -9999;
    red_chi2y_dwnhalf = -9999;

    if(id_innerspan >= 0)
    {
      red_chi2x_inner = trktrack->GetNormChisqX(id_innerspan);
      red_chi2y_inner = trktrack->GetNormChisqY(id_innerspan);
    }

    if(id_uphalf_span >=0)
    {
      red_chi2x_uphalf = trktrack->GetNormChisqX(id_uphalf_span);
      red_chi2y_uphalf = trktrack->GetNormChisqY(id_uphalf_span);
    }
    if(id_dwnhalf_span >=0)
    {
      red_chi2x_dwnhalf = trktrack->GetNormChisqX(id_dwnhalf_span);
      red_chi2y_dwnhalf = trktrack->GetNormChisqY(id_dwnhalf_span);
    }

    for(Int_t j=0 ; j<4 ; j++)
    {
      tof_edep[j] = 0;
    }

    for(Int_t j=0 ; j < betaHR->NTofClusterH() ; j++)
    {
      tofclusterH = betaHR->pTofClusterH(j);
      if(tofclusterH > 0) tof_edep[tofclusterH->Layer] += tofclusterH->GetEdep();
    }

    for(Int_t j=0 ; j<9 ; j++)
    {
      trk_edep[j] = 0;
    }
    for(Int_t j=0 ; j<trktrack->NTrRecHit(); j++)
    {
      trk_edep[trktrack->pTrRecHit(j)->GetLayerJ()-1] += trktrack->pTrRecHit(j)->GetEdep(2);
    }
    Int_t lay;
    Float_t z;
    Int_t Z;
    Z=betaHR->GetZ(lay,z);

    particleFactory.SetAMSTrTrackR(trktrack);
    particleFactory.SetAMSEcalShowerR(ecalshower);
    particleFactory.SetAMSBetaHR(betaHR);

    if(!amsRootSupport.SwitchToSpecificTrackFitById(id_maxspan)) continue;
    Analysis::Event& event = amsRootSupport.BuildEvent(ch, pev);

    eventFactory.PerformTrdTracking(event);
    eventFactory.PerformTrdVertexFinding(event);
    //eventFactory.CreateParticles(event);

    int productionSteps = 0;
    productionSteps |= Analysis::CreateSplineTrack;
    productionSteps |= Analysis::CreateTrdTrack;
   // productionSteps |= Analysis::FillTrdQt;
    eventFactory.FillParticles(event, productionSteps);
    const Analysis::Particle* particle = event.PrimaryParticle();
    assert(particle);

    trdp_e2plikelihood = particle->CalculateElectronProtonLikelihood();
    trdp_he2plikelihood = particle->CalculateHeliumProtonLikelihood();

    if( particle->LowerTofCharge()<=0.5 || particle->LowerTofCharge()>=2.0) cut_tof_charge = false;
    if( particle->UpperTofCharge() <= 0.5 || particle->UpperTofCharge() >=2.0) cut_tof_charge = false;
    charge_tof_up = particle->UpperTofCharge();
    charge_tof_low = particle->LowerTofCharge();
    hEvtCounter2->Fill(10);


    //const std::vector<Analysis::TrdHit>& TrdHit = trdQtFromTrackerTrack->GetAssignedHits();

    const std::vector<Analysis::TrdVertex>& verticesXZ = event.TrdVerticesXZ();
    const std::vector<Analysis::TrdVertex>& verticesYZ = event.TrdVerticesYZ();

    trd_vertex3=0;
    trd_vertex2=0;

    for(std::vector<Analysis::TrdVertex>::const_iterator xzlter = verticesXZ.begin(); xzlter != verticesXZ.end(); ++xzlter)
    {
      const Analysis::TrdVertex& xzVertex = *xzlter;
      for(std::vector<Analysis::TrdVertex>::const_iterator yzlter = verticesYZ.begin(); yzlter != verticesYZ.end(); ++yzlter)
      {
        const Analysis::TrdVertex& yzVertex = *yzlter;
        if(std::max(xzVertex.NumberOfSegments(), yzVertex.NumberOfSegments() ) >=3) //continue;
        {
          if(fabs(xzVertex.Z() - yzVertex.Z() ) < fabs(xzVertex.ErrorZ() + yzVertex.ErrorZ()))
          {
            trd_vertex3++;
          }
        }

        if(std::max(xzVertex.NumberOfSegments(), yzVertex.NumberOfSegments() ) >=2) //continue;
        {
          if(fabs(xzVertex.Z() - yzVertex.Z() ) < fabs(xzVertex.ErrorZ() + yzVertex.ErrorZ()))
          {
            trd_vertex2++;
          }
        }
      }
    }


    ntrdtrack = -9999;
    ntrdtrack = pev->nTrdTrack();  //number of trd track
//coordinate of each detectors
    for(Int_t i=0 ; i<4 ; i++)
    {
      for(Int_t j=0 ; j<3 ; j++)
      {
        tofcoo[i][j]=-99;
      }
    }

    for(Int_t i=0 ; i<9 ; i++)
    {
      for(Int_t j=0 ; j<3 ; j++)
      {
        trkcoo[i][j]=-99;
      }
    }

    for(Int_t i=0 ; i<4 ; i++)
    {
      for(Int_t j=0 ; j<3 ; j++)
      {
        tofcoo[i][j] = part->TOFCoo[i][j];
      }
    }

    for(Int_t i=0 ; i<9 ; i++)
    {
      for(Int_t j=0 ; j<3 ; j++)
      {
        trkcoo[i][j] = part->TrCoo[i][j];
      }
    }

    for(Int_t i=0 ; i<2 ; i++)
    {
      for(Int_t j=0 ; j<3 ; j++)
      {
        richcoo[i][j] = part->RichCoo[i][j];
      }
    }

    for(Int_t i=0 ; i<3 ; i++)
    {
      for(Int_t j=0 ; j<3 ; j++)
      {
        ecalcoo[i][j] = part->EcalCoo[i][j];
      }
    }


    AMSPoint pnt;
    AMSDir dir;

    //coordinate of interpolated track on each detectors

    for(Int_t i=0 ; i<4 ; i++)
    {
      trktrack->Interpolate(part->TOFCoo[i][2], pnt, dir, id_innerspan);
      for(Int_t j=0 ; j<2 ; j++)
      {
        tofdir[i][j] = dir[j];
        trk_tof[i][j] = pnt[j];
      }
    }

    for(Int_t i=0 ; i<9 ; i++)
    {
      trktrack->Interpolate(part->TrCoo[i][2], pnt, dir, id_innerspan);
      for(Int_t j=0 ; j<2 ; j++)
      {
        trkdir[i][j] = dir[j];
        trk_trk[i][j] = pnt[j];
      }
    }

    for(Int_t i=0 ; i<2 ; i++)
    {
      trktrack->Interpolate(part->RichCoo[i][2], pnt, dir, id_innerspan);
      for(Int_t j=0 ; j<2 ; j++)
      {
        richdir[i][j] = dir[j];
        trk_rich[i][j] = pnt[j];
      }
    }

    for(Int_t i=0 ; i<3 ; i++)
    {
      trktrack->Interpolate(part->EcalCoo[i][2], pnt, dir, id_innerspan);
      for(Int_t j=0 ; j<2 ; j++)
      {
        ecaldir[i][j] = dir[j];
        trk_ecal[i][j] = pnt[j];
      }
    }


    ntrtrack = pev->nTrTrack();

    ntrdsegment = pev->nTrdSegment();
    ntrdhsegment = pev->nTrdHSegment();


    for(Int_t i=0 ; i<18 ; i++) 
    {
      ecal_Edep[i]=0;
    }

    if(ecalshower)
    {
      for(Int_t j=0 ; j<ecalshower->NEcal2DCluster(); j++)
      {
        for(Int_t k=0 ; k<ecalshower->pEcal2DCluster(j)->NEcalCluster(); k++)
        {
          ecal_Edep[ecalshower->pEcal2DCluster(j)->pEcalCluster(k)->Plane] += ecalshower->pEcal2DCluster(j)->pEcalCluster(k)->Edep;
        }
      }
    }


    acc_time_diff.clear();
    tof_time.clear();
    acc_edep.clear();

    nacchits = -9999;
    ntofhits = -9999;

    nacchits = pev->nAntiCluster();
    ntofhits = pev->nTofCluster();

    for(Int_t i=0 ; i<nacchits ; i++)
    {
      for(Int_t j=0 ; j<nacchits ; j++)
      {
        if( i==j || i>j) continue;
        AntiClusterR* pAntiCluster_c1 = pev->pAntiCluster(i);
        AntiClusterR* pAntiCluster_c2 = pev->pAntiCluster(j);
        acc_time_c1 = pAntiCluster_c1->time;
        acc_time_c2 = pAntiCluster_c2->time;
        acc_time_diff.push_back(TMath::Abs(acc_time_c1 - acc_time_c2));
      }
    }

    for(Int_t i=0 ; i<nacchits ; i++)
    {
      AntiClusterR* pAntiCluster = pev->pAntiCluster(i);
      acc_edep.push_back(pAntiCluster->Edep);
    }

    for(Int_t i=0 ; i<ntofhits ; i++)
    {
      TofClusterR* pTofCluster = pev->pTofCluster(i);
      tof_time.push_back(pTofCluster->Time);
    }

    
    ntrkhits = -9999;
    ntrkxhits = -9999;
    ntrkxyhits = -9999;
    ntrkxyhits = -9999;

    ntrkhits = trktrack->GetNhits();
    ntrkxhits = trktrack->GetNhitsX();
    ntrkyhits = trktrack->GetNhitsY();
    ntrkxyhits = trktrack->GetNhitsXY();


    for(Int_t i=0 ; i<9 ; i++)
    {
      trkhit_plane[i] = (trktrack->GetBitPatternJ()>>i)&1;
    }
    //trd_segment trdsegment;
    
    evt_ntrd_cluster.clear();
    part_ntrd_cluster.clear();
    evt_segment_coo_x.clear();
    evt_segment_coo_y.clear();
    evt_segment_coo_z.clear();
    part_segment_coo_x.clear();
    part_segment_coo_y.clear();
    part_segment_coo_z.clear();
    evt_segment_direction.clear();
    part_segment_direction.clear();

    for(Int_t i=0 ; i<pev->NTrdSegment() ; i++)
    {
      if(pev->pTrdSegment(i)) evt_ntrd_cluster.push_back(pev->pTrdSegment(i)->NTrdCluster());
      if(pev->pTrdSegment(i))
      {
        for(Int_t j=0 ; j<pev->pTrdSegment(i)->NTrdCluster() ; j++)
        {
          if(!pev->pTrdSegment(i)->pTrdCluster(j)) continue;
          evt_segment_coo_x.push_back(pev->pTrdSegment(i)->pTrdCluster(j)->Coo[0]);
          evt_segment_coo_y.push_back(pev->pTrdSegment(i)->pTrdCluster(j)->Coo[1]);
          evt_segment_coo_z.push_back(pev->pTrdSegment(i)->pTrdCluster(j)->Coo[2]);
          evt_segment_direction.push_back(pev->pTrdSegment(i)->pTrdCluster(j)->Direction);
          //if(pev->pTrdSegment(i)->pTrdCluster(j)->Coo[1]==0) 
          //cout<<"x : "<<pev->pTrdSegment(i)->pTrdCluster(j)->Coo[0]<<"   y : "<<pev->pTrdSegment(i)->pTrdCluster(j)->Coo[1]<<"   z : "<<pev->pTrdSegment(i)->pTrdCluster(j)->Coo[2]<<endl;
        }
      }
    }
    if(part->pTrdTrack())
    {
      for(Int_t i=0 ; i<part->pTrdTrack()->NTrdSegment() ; i++)
      {
        if(!part->pTrdTrack()->pTrdSegment(i)) continue;

        part_ntrd_cluster.push_back(part->pTrdTrack()->pTrdSegment(i)->NTrdCluster());
        for(Int_t j=0 ; j<part->pTrdTrack()->pTrdSegment(i)->NTrdCluster() ; j++)
        {
        
          //part_segment_coo->emplace_back(part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Coo[0], part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Coo[1], part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Coo[2]);
          if(!part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)) continue;
          part_segment_coo_x.push_back(part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Coo[0]);
          part_segment_coo_y.push_back(part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Coo[1]);
          part_segment_coo_z.push_back(part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Coo[2]);
          part_segment_direction.push_back(part->pTrdTrack()->pTrdSegment(i)->pTrdCluster(j)->Direction);
        }   
      }
    }


    TrdKCluster trdk;
    Double_t likelihood_ratio[3]={-999,-999,-999};
    Double_t likelihood[3];
    Int_t trd_nhits;
    Float_t totalpath, totalamp;
    Int_t trdk_isvalid;
    if(pev->NTrdRawHit() > 0)
    {
      trdk = TrdKCluster(pev, trktrack, id_maxspan);
      trdk_isvalid = trdk.GetLikelihoodRatio_TrTrack(15, likelihood_ratio, likelihood, trd_nhits, totalpath, totalamp, -1, 0);
      if(trdk_isvalid != 0 && trd_nhits != 0)
      {
        trdk_e2plikelihood = likelihood_ratio[0];
        trdk_he2plikelihood = likelihood_ratio[1];
      }
    }


    //calculate number of layers which have cluster

    Int_t b[4] = {0,0,0,0};

    if(betaHR->NTofClusterH()==4) 
    {
      for(Int_t i=0 ; i<betaHR->NTofClusterH() ; i++)
      {
        Int_t tof_layer = betaHR->pTofClusterH(i)->Layer;
        if(tof_layer == 0) b[0]=1;
        if(tof_layer == 1) b[1]=1;
        if(tof_layer == 2) b[2]=1;
        if(tof_layer == 3) b[3]=1;
      }

      Int_t hist_fill=0;
      for(Int_t i=0 ; i<4 ; i++)
      {
        hist_fill +=b[i];
      }
      if(hist_fill != 4) cut_tof_cluster = false;
    }

    if(pev->GetNTofClustersInTime(betaHR, ncls) > 2) cut_tof_cluster = false;

    tofcluster_part->Fill(betaHR->NTofClusterH());

    trk_theta = trktrack->GetTheta(id_maxspan);
    trk_phi = trktrack->GetPhi(id_maxspan);
    trk_pattern_X = trktrack->GetPatternX();
    trk_pattern_Y = trktrack->GetPatternY();
    trk_pattern_XY = trktrack->GetPatternXY();

    Ntof_cluster = betaHR->NTofClusterH();
    Ntof_cluster_intime = pev->GetNTofClustersInTime(betaHR, ncls);
    beta_pattern = betaHR->GetBetaPattern();
    tof_chi2_time = betaHR->GetChi2T();
    tof_chi2_coord = betaHR->GetChi2C();
    tof_mass = betaHR->GetMass();
    tof_mass = betaHR->GetEMass();
    tof_used_hit = betaHR->GetUseHit();
    tof_sum_hit = betaHR->GetSumHit();
    
    for(Int_t i=0 ; i<4 ; i++)
    {
      tof_cluster_pattern[i]=0;
      tof_cluster_layer[i]=0;
      tof_cluster_Aedep[i]=0;
      tof_cluster_Dedep[i]=0;
      tof_cluster_hit_time[i]=0;
      tof_cluster_hit_time_err[i]=0;
      for(Int_t j=0 ; j<2 ; j++)
      {
        tof_cluster_non_correct_time[j][i]=0;
        tof_cluster_correct_time[j][i]=0;
        tof_cluster_Acharge[j][i]=0;
      }
      for(Int_t j=0 ; j<2 ; j++)
      {
        for(Int_t k=0 ; k<3 ; k++)
        {
          tof_cluster_Dcharge[j][k][i]=0;
        }
      }
    }


    for(Int_t i=0 ; i<Ntof_cluster ; i++)
    {
      tof_cluster_pattern[i] = betaHR->pTofClusterH(i)->Pattern;
      tof_cluster_layer[i] = betaHR->pTofClusterH(i)->Layer;
      tof_cluster_Aedep[i] = betaHR->pTofClusterH(i)->AEdep;
      tof_cluster_Dedep[i] = betaHR->pTofClusterH(i)->DEdep;
      tof_cluster_hit_time[i] = betaHR->pTofClusterH(i)->Time;
      tof_cluster_hit_time_err[i] = betaHR->pTofClusterH(i)->ETime;
      for(Int_t j=0 ; j<2 ; j++)
      {
        tof_cluster_non_correct_time[j][i] = betaHR->pTofClusterH(i)->Rtime[j];
        tof_cluster_correct_time[j][i] = betaHR->pTofClusterH(i)->Stime[j];
        tof_cluster_Acharge[j][i] = betaHR->pTofClusterH(i)->AQ2[j];
      }
      for(Int_t j=0 ; j<2 ; j++)
      {
        for(Int_t k=0 ; k<3 ; k++)
        {
          tof_cluster_Dcharge[j][k][i] = betaHR->pTofClusterH(i)->DQ2[j][k];
        }
      }
    }

    rich_used_hits = 0;
    rich_photoelectrons = 0;
    rich_expect_photoelectrons = 0;
    rich_prob = 0;
    rich_track_theta = 0;
    rich_track_phi = 0;
    rich_Nhits_ring = 0;
    rich_Npmts_ring = 0;
    rich_isgood = 0;
    rich_isclean=0;
    
    if(!richring)
    {
      rich_used_hits = 0;
      rich_photoelectrons=0;
      rich_expect_photoelectrons = 0;
      rich_prob=0;
      rich_track_theta=0;
      rich_track_phi=0;
      rich_Nhits_ring=0;
      rich_Npmts_ring=0;
      rich_isgood = 0;
      rich_isclean = 0;
    }

    if(richring)
    {
      rich_used_hits = richring->getUsedHits();
      rich_photoelectrons = richring->getPhotoElectrons();
      rich_expect_photoelectrons = richring->getExpectedPhotoElectrons();
      rich_prob = richring->getProb();
      rich_track_theta = richring->getTrackTheta();
      rich_track_phi = richring->getTrackPhi();
      rich_Nhits_ring = richring->getHits();
      rich_Npmts_ring = richring->getPMTs();
      rich_isgood = richring->IsGood();
      rich_isclean = richring->IsClean();
    }







    //Analysis::Particle* acparticle;

    //trdp_e2plikelihood = acparticle->CalculateElectronProtonLikelihood();
    //trdp_he2plikelihood = acparticle->CalculateHeliumProtonLikelihood();
    


    /*
    for(Int_t i=1 ; i<10 ; i++)
    {
      TrRecHitR* trrechitr = trktrack->GetHitLJ(i);
      if(!trrechitr) continue;
      cout<<"####"<<trrechitr->GetLayerJ();
    }
    cout<<endl;
*/


   
    //if(((trig->PhysBPatt>>2)&1)==1 || ((trig->PhysBPatt>>3)&1)==1) cout<<"Z(bit3 or 4) : "<<Z<<endl;
    //if(((trig->PhysBPatt>>1)&1)==1) cout<<"Z(bit2) : "<<Z<<endl;
    //if(!cut_tof_chi2 || !cut_trk_charge || !cut_trk_yhit || !cut_trk_chi2 || !cut_tof_charge || !cut_tof_cluster) continue;

    hEvtCounter2->Fill(11);

    t->Fill();
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  tofcluster_layer->Draw();

  f->Write();
  f->Close();

  return 0;
}

bool CheckBadRun(AMSEventR *this_evt)
{
  int tag = this_evt->IsBadRun("");
  if(tag==0) return false;
  if(tag==1) return true;
  if(tag==2) return true;
}

bool CheckScienceRun(AMSEventR *this_evt)
{
  header = &(this_evt->fHeader);
  if(!header) return false;
  if( (header->RunType >> 12) != 0xf) return false;
  return true;
}

bool GoodHW(DaqEventR *daq)
{
  bool goodHW = true;
  if(!daq) return false;
  bool errors = false;
  for(int ii=0 ; ii<4 ; ii++) errors |= (bool)( (daq->JINJStatus[ii]>>8) & 0x7F);
  for(int ii=0 ; ii<24 ; ii++) errors |= (bool) ( daq->JError[ii] & 0x7F);
  if(errors) goodHW=false;
  return (goodHW);
}

bool GoodTrkTrack(TrTrackR *this_track)
{
  if(this_track->IsFake()) return false;
  id_maxspan = this_track->iTrTrackPar(1,0,20);
  if(id_maxspan<0) return false;
  return true;
}
