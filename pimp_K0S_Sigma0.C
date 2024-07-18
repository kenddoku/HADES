#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"


void pimp_K0S_Sigma0(){
  PUtils::SetSeed(0);

  makeDistributionManager()->Enable("strangeness");

  makeStaticData()->AddDecay("Sigma0 --> g + Lambda", "Sigma0", "g, Lambda", 1.);

  TFile * f1 = new TFile("K0Sigma0.root","recreate");

  //*******Creating histograms for every tracked particle's momentum and theta angle*******//--------------

  TH1F * h_mom_p_track1 = new TH1F("h_mom_p_track1","h_mom_p_track1",500,0,4500);
  TH1F * h_theta_p_track1 = new TH1F("h_theta_p_track1","h_theta_p_track1",360,0,180);
  TH1F * h_momCut_p_track1 = new TH1F("h_momCut_p_track1","h_momCut_p_track1",500,0,4500);

  TH1F * h_mom_pim_track2 = new TH1F("h_mom_pim_track2","h_mom_pim_track2",500,0,4500);
  TH1F * h_theta_pim_track2 = new TH1F("h_theta_pim_track2","h_theta_pim_track2",360,0,180);
  TH1F * h_momCut_pim_track2 = new TH1F("h_mom_pimCut_track2","h_momCut_pim_track2",500,0,4500);  

  TH1F * h_mom_pip_track3 = new TH1F("h_mom_pip_track3","h_mom_pip_track3",500,0,4500);
  TH1F * h_theta_pip_track3 = new TH1F("h_theta_pip_track3","h_theta_pip_track3",360,0,180);
  TH1F * h_momCut_pip_track3 = new TH1F("h_momCut_pip_track3","h_momCut_pip_track3",500,0,4500);  

  TH1F * h_mom_pim_track4 = new TH1F("h_mom_pim_track4","h_mom_pim_track4",500,0,4500);
  TH1F * h_theta_pim_track4 = new TH1F("h_theta_pim_track4","h_theta_pim_track4",360,0,180);
  TH1F * h_momCut_pim_track4 = new TH1F("h_momCut_pim_track4","h_momCut_pim_track4",500,0,4500);  

  //*******Creating other histograms*******//--------------------------------------------------------------
  TH1F * h_mass_Lambda = new TH1F("h_mass_Lambda","h_mass_Lambda",100,0,2000);
  TH1F * h_mass_K0 = new TH1F("h_mass_K0","h_mass_K0",600,0,1);

  //*******Defining occuring reactions*******//------------------------------------------------------------
  PReaction reac_1(1.16,"pi-","p","K0S [pi+ pi-] Sigma0 [g Lambda [p pi-]]","pimp_test",1,0,1,1);
  PReaction reac_2(1.16, "pi-", "p", "K0S [pi+ pi-] Sigma0 [g Lambda [n pi0 [g g]]]", "pimp_test", 1,0,1,1);
 // PReaction reac_2(1.16, "pi-", "p", "K0S [pi+ pi-] Sigma0 [Lambda [n pi0 [g g]]", "pimp_test", 1,0,1,1);

 //PReaction my_reaction2(1.16,"pi-","p","K0S [pi+ pi-] Lambda [p pi-]","pimp_test",1,0,1,1);

 //*******Calculating momentum and theta angles of all tracked particles*******//--------------------------
 reac_1.Do("mom_p_track1 = [p,1]->P()*1000.;");
 reac_1.Do("theta_p_track1 = ([p,1]->Theta()*180.)/TMath::Pi();");

 reac_1.Do("mom_pim_track2 = [pi-,2]->P()*1000.;");
 reac_1.Do("theta_pim_track2 = ([pi-,2]->Theta()*180.)/TMath::Pi();");

 reac_1.Do("mom_pip_track3 = [pi+,3]->P()*1000.;");
 reac_1.Do("theta_pip_track3 = ([pi+,3]->Theta()*180.)/TMath::Pi();");

 reac_1.Do("mom_pim_track4 = [pi-,4]->P()*1000.;");
 reac_1.Do("theta_pim_track4 = ([pi-,4]->Theta()*180.)/TMath::Pi();");

 //reac_2.Do("mom_n_track1 = [n,1]->P()*1000.;");

 //*******Calculating masses of parent particles*******//--------------------------------------------------
 reac_1.Do("Lambda_invmass = ([p,1]+[pi-,2])->M();");
 reac_1.Do("K0_invmass = ([pi-,4]+[pi+,3])->M();");
 //my_reaction2.Do("thp = ([p]->Theta() * 180.)/TMath::Pi();");
 //my_reaction2.Do("mom_p = [p]->P()*1000.;");


 //my_reaction2.Do("beam1=[p+p];");
 //my_reaction2.Do("beam2=P3E(0,0,5.3565,6.376);");
 //my_reaction2.Do("mm1=(beam1-([p,1]+[p,2]))->M()*1000;");
 //my_reaction2.Do("mm2=(beam2-([p,1]+[p,2]))->M()*1000;");


 //my_reaction2.Do("invMpippim=([pi-]+[pi+])->M();");
 //my_reaction2.Do("echo $invMpippim;");

//*******Filling histograms*******//-----------------------------------------------------------------------
 reac_1.Do(h_mom_p_track1, "_x=mom_p_track1;");
 reac_1.Do(h_theta_p_track1, "_x=theta_p_track1;");
 reac_1.Do(h_momCut_p_track1, "if theta_p_track1 > 18 && theta_p_track1 < 90; _x=mom_p_track1;");
 
 reac_1.Do(h_mom_pim_track2, "_x=mom_pim_track2;");
 reac_1.Do(h_theta_pim_track2, "_x=theta_pim_track2;");
 reac_1.Do(h_momCut_pim_track2, "if theta_pim_track2 > 18 && theta_pim_track2 < 90; _x=mom_pim_track2;");

 reac_1.Do(h_mom_pip_track3, "_x=mom_pip_track3;");
 reac_1.Do(h_theta_pip_track3, "_x=theta_pip_track3;");
 reac_1.Do(h_momCut_pip_track3, "if theta_pip_track3 > 18 && theta_pip_track3 < 90; _x=mom_pip_track3;");

 reac_1.Do(h_mom_pim_track4, "_x=mom_pim_track4;");
 reac_1.Do(h_theta_pim_track4, "_x=theta_pim_track4;");
 reac_1.Do(h_momCut_pim_track4, "if theta_pim_track4 > 18 && theta_pim_track4 < 90; _x=mom_pim_track4;");

 reac_1.Do(h_mass_K0, "_x=K0_invmass");
 reac_1.Do(h_mass_Lambda, "_x=Lambda_invmass")

 reac_1.Print();
 reac_1.Loop(200000);

 //reac_2.Print();
 //reac_2.Loop(200000);

//  TCanvas *c1 = new TCanvas("theta","theta");
//  c1->Divide(3,1);

 f1->cd();

 h_mom_p_track1->Write();
 h_theta_p_track1->Write();
 h_momCut_p_track1->Write();

 h_mom_pim_track2->Write();
 h_theta_pim_track2->Write();
 h_momCut_pim_track2->Write();

 h_mom_pip_track3->Write();
 h_theta_pip_track3->Write();
 h_momCut_pip_track3->Write();

 h_mom_pim_track4->Write();
 h_theta_pim_track4->Write();
 h_momCut_pim_track4->Write();

 h_mass_K0->Write();
 h_mass_Lambda->Write();

//  h_mass_K0->Write();
//  c1->Write();

 f1->Close();

}



