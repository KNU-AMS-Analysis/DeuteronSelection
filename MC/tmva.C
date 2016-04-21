
void tmva()
{
  cout<<"check before running program"<<endl;
  Int_t i;
  cin>>i;
  cout<<"/////////////////////////////start/////////////////////////////"<<endl;
  TDatime* time = new TDatime();
  time->Print();
  TStopwatch* stopwatch = new TStopwatch();
  stopwatch->Start();
  TFile* output;

  //TString fname = "/Users/sckang/work/BDT/aachen/data_d_mc/add_ecal_trd/data_mc_d1_20_rig_resol.root";
  TString fname = "../../Data1/deuteron_B1059.root";
  //TString fname = "/Users/sckang/work/BDT/aachen/data_d_mc/data_mc_d1_20.root";
  //TString fname = "/Users/sckang/work/BDT/aachen/data_d_mc/data_mc_d0.root";
  TFile* input_d = new TFile(fname);

  //TString fname1 = "/Users/sckang/work/BDT/aachen/data_p_mc/add_ecal_trd/data_mc_p1_20_rig_resol.root";
  TString fname1 = "../../Data1/data_p.root";
  //TString fname1 = "/Users/sckang/work/BDT/aachen/data_p_mc/data_mc_p1_20.root";
  //TString fname1 = "/Users/sckang/work/BDT/aachen/data_p_mc/data_mc_p0.root";
  TFile* input_p = new TFile(fname1);

  TTree* sigtree;
  TTree* backtree;

  TMVA::Factory* factory;

  TMVA::Reader* reader;

  Float_t var1, var2, var3, var4, var5, var6, var7, var8, var9, var10,
          var11, var12, var13, var14, var15, var16, var17, var18, var19, var20,
          var21, var22, var23, var24, var25, var26, var27, var28, var29, var30,
          var31, var32, var33, var34, var35, var36, var37, var38, var39, var40,
          var41, var42, var43, var44, var45, var46, var47, var48, var49, var50,
          var51, var52, var53, var54, var55, var56, var57, var58, var59, var60,
          var61, var62, var63, var64, var65, var66, var67, var68, var69, var70,
          var71, var72, var73, var74, var75, var76, var77, var78, var79, var80,
          var81, var82, var83, var84, var85, var86, var87, var88, var89, var90,
          var91, var92, var93, var94, var95, var96, var97, var98, var99, var100,
          var101;

  cout<<"pass here 1"<<endl;

  output = new TFile("tmva_method.root", "recreate");

  //factory = new TMVA::Factory("tmva_method", output, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;G,D:AnalysisType=Classification");
  factory = new TMVA::Factory("tmva_method", output, "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");

  sigtree = (TTree*)input_d->Get("mc_d_data");
  backtree = (TTree*)input_p->Get("mc_p_data");

    cout<<"pass here 2"<<endl;
    
  factory->AddVariable("beta", 'F');
  factory->AddVariable("inv_beta", 'F');
  factory->AddVariable("mean_trk_edep", 'F');
  factory->AddVariable("charge_trk", 'F');
  //factory->AddVariable("trk_pattern_X", 'I');
  //factory->AddVariable("trk_pattern_Y", 'I');
  //factory->AddVariable("Ntof_cluster_intime", 'I');
  factory->AddVariable("rich_Nhits_ring", 'I');
  factory->AddVariable("rich_Npmts_ring", 'I');
  factory->AddVariable("rich_beta_naf", 'F');
  factory->AddVariable("rich_beta_aer", 'F');
  factory->AddVariable("rich_beta_err", 'F');
  //factory->AddVariable("tof_chi2_time", 'F');
  factory->AddVariable("tof_chi2_coord", 'F');
  factory->AddVariable("tof_mass", 'F');
  factory->AddVariable("tof_mass_err", 'F');
  factory->AddVariable("rich_mass", 'F');
  factory->AddVariable("rich_mass_err", 'F');
  //factory->AddVariable("tof_cluster_hit_time1", 'F');
  //factory->AddVariable("tof_cluster_hit_time2", 'F');
  //factory->AddVariable("tof_cluster_hit_time3", 'F');
  //factory->AddVariable("tof_cluster_hit_time4", 'F');
  factory->AddVariable("rich_photoelectrons", 'F');
  factory->AddVariable("rich_prob", 'F');
  factory->AddVariable("rich_track_theta", 'F');
  factory->AddVariable("rich_track_phi", 'F');
  factory->AddVariable("charge_tofL1", 'F');
  factory->AddVariable("charge_tofL2", 'F');
  factory->AddVariable("charge_tofL3", 'F');
  factory->AddVariable("charge_tofL4", 'F');
  //factory->AddVariable("trk_chi2_max_X", 'F');
  //factory->AddVariable("trk_chi2_max_Y", 'F');
  //factory->AddVariable("trk_chi2_inner_X", 'F');
  //factory->AddVariable("trk_chi2_inner_Y", 'F');
  //factory->AddVariable("inv_rigidity_max", 'F');
  //factory->AddVariable("inv_rigidity_inner", 'F');
  factory->AddVariable("inv_rigidity_max_err", 'F');
  //factory->AddVariable("inv_rigidity_inner_err", 'F');
  factory->AddVariable("rigidity_beta", 'F');
  //factory->AddVariable("diff_rigbeta_tofbeta", 'F');
  factory->AddVariable("diff_trkchi2Y_inner_max", 'F');


  factory->AddSpectator("rigidity_max", 'F');
  //factory->AddSpectator("rigidity_inner", 'F');
  factory->AddSpectator("rich_passcut", 'I');
  factory->AddSpectator("rich_isnaf", 'I');

  cout<<"pass here 3"<<endl;

  factory->AddSignalTree(sigtree);
  factory->AddBackgroundTree(backtree);

  factory->BookMethod( TMVA::Types::kBDT, "BDT","!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
  cout<<"pass here@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;

  TMVA::MethodCategory* mcat = 0;

  TString catvar_naf = "beta:inv_beta:mean_trk_edep:charge_trk:rich_Nhits_ring:rich_Npmts_ring:rich_beta_naf:rich_beta_err:tof_chi2_coord:tof_mass:tof_mass_err:rich_mass:rich_mass_err:rich_photoelectrons:rich_prob:rich_track_theta:rich_track_phi:charge_tofL1:charge_tofL2:charge_tofL3:charge_tofL4:inv_rigidity_max_err:rigidity_beta:diff_trkchi2Y_inner_max";

  TString catvar_aer = "beta:inv_beta:mean_trk_edep:charge_trk:rich_Nhits_ring:rich_Npmts_ring:rich_beta_aer:rich_beta_err:tof_chi2_coord:tof_mass:tof_mass_err:rich_mass:rich_mass_err:rich_photoelectrons:rich_prob:rich_track_theta:rich_track_phi:charge_tofL1:charge_tofL2:charge_tofL3:charge_tofL4:inv_rigidity_max_err:rigidity_beta:diff_trkchi2Y_inner_max";

  TString catvar_norich = "beta:inv_beta:mean_trk_edep:charge_trk:tof_chi2_coord:tof_mass:tof_mass_err:charge_tofL1:charge_tofL2:charge_tofL3:charge_tofL4:inv_rigidity_max_err:rigidity_beta:diff_trkchi2Y_inner_max";

  //TString catvar_norich = "beta:mean_trk_edep:mean_tof_edep:charge_trk:trk_pattern_X:trk_pattern_Y:Ntof_cluster_intime:beta_pattern:rich_Nhits_ring:rich_Npmts_ring:tof_used_hit:tof_chi2_time:tof_chi2_coord:tof_mass:tof_cluster_Aedep1:tof_cluster_Aedep2:tof_cluster_Aedep3:tof_cluster_Aedep4:tof_cluster_Dedep1:tof_cluster_Dedep2:tof_cluster_Dedep3:tof_cluster_Dedep4:tof_cluster_hit_time1:tof_cluster_hit_time2:tof_cluster_hit_time3:tof_cluster_hit_time4:tof_cluster_hit_time_err1:tof_cluster_hit_time_err2:tof_cluster_hit_time_err4:tof_cluster_non_correct_time11:tof_cluster_non_correct_time12:tof_cluster_non_correct_time13:tof_cluster_non_correct_time14:tof_cluster_non_correct_time21:tof_cluster_non_correct_time22:tof_cluster_non_correct_time23:tof_cluster_non_correct_time24:tof_cluster_correct_time11:tof_cluster_correct_time12:tof_cluster_correct_time13:tof_cluster_correct_time14:tof_cluster_correct_time21:tof_cluster_correct_time22:tof_cluster_correct_time23:tof_cluster_correct_time24:tof_cluster_Acharge11:tof_cluster_Acharge12:tof_cluster_Acharge13:tof_cluster_Acharge14:tof_cluster_Acharge21:tof_cluster_Acharge22:tof_cluster_Acharge23:tof_cluster_Acharge24:tof_cluster_Dcharge111:tof_cluster_Dcharge112:tof_cluster_Dcharge113:tof_cluster_Dcharge114:tof_cluster_Dcharge121:tof_cluster_Dcharge122:tof_cluster_Dcharge123:tof_cluster_Dcharge124:tof_cluster_Dcharge131:tof_cluster_Dcharge211:tof_cluster_Dcharge212:tof_cluster_Dcharge213:tof_cluster_Dcharge214:tof_cluster_Dcharge221:tof_cluster_Dcharge222:tof_cluster_Dcharge223:tof_cluster_Dcharge224:tof_cluster_Dcharge231:rich_photoelectrons:rich_prob:rich_track_theta:rich_track_phi:charge_tofL1:charge_tofL2:charge_tofL3:charge_tofL4";

  TMVA::MethodBase* bdtcat = factory->BookMethod( TMVA::Types::kCategory, "bdtcat", "");
  //TMVA::IMethod* bdtcat = factory->BookMethod( TMVA::Types::kCategory, "bdtcat", "");
  mcat = dynamic_cast<TMVA::MethodCategory*>(bdtcat);
  mcat->AddMethod("rich_passcut==1 && rich_isnaf == 1", catvar_naf, TMVA::Types::kBDT, "category_bdt_1", "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
  mcat->AddMethod("rich_passcut==1 && rich_isnaf == 0", catvar_aer, TMVA::Types::kBDT, "category_bdt_2", "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
  mcat->AddMethod("rich_passcut==0 ", catvar_norich, TMVA::Types::kBDT, "category_bdt_3", "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
  //mcat->AddMethod("rich_passcut==0", catvar_norich, TMVA::Types::kBDT, "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
  

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  cout<<"pass here#############################"<<endl;

  reader = new TMVA::Reader();
    
  reader->AddVariable("beta", &var1);
  reader->AddVariable("inv_beta", &var44);
  reader->AddVariable("mean_trk_edep", &var2);
  //reader->AddVariable("mean_tof_edep", &var3);
  reader->AddVariable("charge_trk", &var4);
  //reader->AddVariable("trk_pattern_X", &var5);
  //reader->AddVariable("trk_pattern_Y", &var6);
  //reader->AddVariable("Ntof_cluster_intime", &var9);
  reader->AddVariable("rich_Nhits_ring", &var20);
  reader->AddVariable("rich_Npmts_ring", &var21);
  reader->AddVariable("rich_beta_naf", &var23);
  reader->AddVariable("rich_beta_aer", &var24);
  reader->AddVariable("rich_beta_err", &var31);
  //reader->AddVariable("tof_chi2_time", &var25);
  reader->AddVariable("tof_chi2_coord", &var26);
  reader->AddVariable("tof_mass", &var27);
  reader->AddVariable("tof_mass_err", &var28);
  reader->AddVariable("rich_mass", &var29);
  reader->AddVariable("rich_mass_err", &var30);
  //reader->AddVariable("tof_cluster_hit_time1", &var36);
  //reader->AddVariable("tof_cluster_hit_time2", &var37);
  //reader->AddVariable("tof_cluster_hit_time3", &var38);
  //reader->AddVariable("tof_cluster_hit_time4", &var39);
  reader->AddVariable("rich_photoelectrons", &var92);
  reader->AddVariable("rich_prob", &var94);
  reader->AddVariable("rich_track_theta", &var95);
  reader->AddVariable("rich_track_phi", &var96);
  reader->AddVariable("charge_tofL1", &var97);
  reader->AddVariable("charge_tofL2", &var99);
  reader->AddVariable("charge_tofL3", &var98);
  reader->AddVariable("charge_tofL4", &var101);
  //reader->AddVariable("trk_chi2_max_X", &var40);
  //reader->AddVariable("trk_chi2_max_Y", &var41);
  //reader->AddVariable("trk_chi2_inner_X", &var45);
  //reader->AddVariable("trk_chi2_inner_Y", &var46);
  //reader->AddVariable("inv_rigidity_max", &var47);
  //reader->AddVariable("inv_rigidity_inner", &var48);
  reader->AddVariable("inv_rigidity_max_err", &var49);
  //reader->AddVariable("inv_rigidity_inner_err", &var50);
  reader->AddVariable("rigidity_beta", &var50);
  //reader->AddVariable("diff_rigbeta_tofbeta", &var50);
  reader->AddVariable("diff_trkchi2Y_inner_max", &var53);

  reader->AddSpectator("rigidity_max", &var43);
  //reader->AddSpectator("rigidity_inner", &var42);
  reader->AddSpectator("rich_passcut", &var51);
  reader->AddSpectator("rich_isnaf", &var52);

  reader->BookMVA("BDT", "weights/tmva_method_BDT.weights.xml");

  if(!gROOT->IsBatch()) TMVA::TMVAGui("tmva_method.root");

  output->Close();

  cout<<"///////////////////////////////end//////////////////////////////"<<endl;
  stopwatch->Stop();
  time->Print();
  stopwatch->Print();


  TFile* result;
  TTree* testtree;
  TH1F* hsig;
  TH1F* hback;
  THStack* hs;

  result = new TFile("tmva_method.root");
  //result[i] = new TFile(Form("tmva_method_%d.root",i));
  testtree = (TTree*)result->Get("TestTree");
    
  hsig = new TH1F("hsig", "hsig", 200, -1, 1);
  hsig->SetLineColor(2);
  hback = new TH1F("hback", "hback", 200, -1, 1);
  hback->SetLineColor(4);

  testtree->Draw("BDT>>hsig", "classID == 0", "goff");
  testtree->Draw("BDT>>hback", "classID == 1", "goff");
  //cout<<hback->GetNbinsX()<<"   "<<hback[i]->Integral()<<endl;

  hs = new THStack("hs", "BDT;BDT");
  hs->Add(hsig);
  hs->Add(hback);
  

  TCanvas* canvas;

  canvas = new TCanvas("canvas", "BDT", 800, 600);
  canvas->SetLogy();
  canvas->cd();
  hs->Draw("nostack");
  canvas->SaveAs("png/bdt.png");

}
