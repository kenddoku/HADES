

#include "ecal_tests.h"
#include "emcdef.h"

#include <hades.h>
#include <hdst.h>
#include <hparasciifileio.h>
#include <hspectrometer.h>
#include <hparticlecand.h>
#include "hemcneutralcand.h"
#include "hemcneutralcandsim.h"
#include "hemcclustersim.h"
#include "hemccluster.h"
#include "hparticlecandsim.h"
#include "hgeantkine.h"

#include <hgeomvolume.h>
#include <hgeomcompositevolume.h>
#include <hgeomvector.h>
#include "hparticleevtinfo.h"
#include "hparticletracksorter.h"
#include "hphysicsconstants.h"
#include "htool.h"
#include "htime.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2I.h>
#include <TH3I.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TMath.h>

#include <algorithm>
#include <cstdlib>
#include <vector>

#define PI 3.14159265

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229;

Int_t ecal_tests(TString inputlist, TString outfile, TString eventsList, Int_t nev)
{
	HLoop* loop = new HLoop(kTRUE);
         
   	if(inputlist.EndsWith(".list")) { loop->addFilesList(inputlist);}
   	else {loop->addMultFiles(inputlist);}

	//----------------------------------Check if loop was properly initialized--------------------------------------
   	if(!loop->setInput("-*,+HParticleCand,+HEmcNeutralCand,+HEmcCluster"))
		//-------------------------------Reading file structure------------------------------------------------------
		{ // reading file structure
      std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
      std::exit(EXIT_FAILURE);
   	}

   	gHades->setBeamTimeID(HADES::kFeb22);
   	HBeamTime::printBeamTimeInfo();

	//----------------------------------Timer for checking analysis time--------------------------------------------
	TStopwatch timer;
	timer.Reset(); // Reset timer
	timer.Start(); // Start timer (T0)

	//----------------------------------Hades Particle Candidates---------------------------------------------------
	HCategory* fParticleCand = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCandSim");
	if(!fParticleCand) {cout << "No catParticleCand!" << endl;}

	HCategory* fEmcNeutralCand = HCategoryManager::getCategory(catEmcNeutralCand, kTRUE, "catEmcNeutralCandSim");
	if(!fEmcNeutralCand) {cout << "No catEmcNeutralCand!" << endl;}
		
	HCategory* fEmcCluster = HCategoryManager::getCategory(catEmcCluster, 0, "catEmcCluster");
	if(!fEmcCluster) {cout << "No catEmcCluster!" << endl;}

	//----------------------------------Energy losses setup---------------------------------------------------------
	HEnergyLossCorrPar dEdxCorr;
	dEdxCorr.setDefaultPar("feb22");

   //----------------------------------Setting parameters for loop over events-------------------------------------
	//??????????????????????????????????????????????????????????????????????????????????????????????????????????????

	//----------------------------------Specifying output file------------------------------------------------------
	TFile* output_file = new TFile(outfile.Data(), "RECREATE");
	output_file->cd();
	cout << "NEW ROOT TREE " << endl;

   //----------------------------------------------------------------------------------------------
  
	TH2F *hbeta_mom = new TH2F("hbeta_mom","hbeta_mom",1000,-4000.,4000.,700,0.,1.4);
	TH2F *hbeta_mom_p = new TH2F("hbeta_mom_p","hbeta_mom_p",1000,0.,4000.,700,0.,1.4);
	TH2F *hbeta_mom_pip = new TH2F("hbeta_mom_pip","hbeta_mom_pip",1000,0.,4000.,700,0.,1.4);
	TH2F *hbeta_mom_pim = new TH2F("hbeta_mom_pim","hbeta_mom_pim",1000,-4000.,0.,700,0.,1.4);
	TH2F *hbeta_mom_epem = new TH2F("hbeta_mom_epem","hbeta_mom_epem",1000,-4000.,4000.,700,0.,1.4);
	TH2F *hbeta_mom_epem_PID = new TH2F("hbeta_mom_epem_PID","hbeta_mom_epem_PID",1000,-4000.,4000.,700,0.,1.4);
	TH2F *htheta_mom = new TH2F("htheta_mom","htheta_mom (except gammas)",1000,0,2000,360,0,90);
	TH2F *htheta_mom_p = new TH2F("htheta_mom_p","htheta_mom_p",1000,0,2000,360,0,90); // <-- theta vs momentum
	TH2F *htheta_mom_pip = new TH2F("htheta_mom_pip","htheta_mom_pip",1000,0,2000,720,0,90);
	TH2F *htheta_mom_pim = new TH2F("htheta_mom_pim","htheta_mom_pim",1000,0,2000,720,0,90);
	TH2F *htheta_mom_K0_to_pim = new TH2F("htheta_mom_K0_to_pim","htheta_mom_K0_to_pim",1000,0,2000,720,0,90);
	TH2F *htheta_mom_L0_to_pim = new TH2F("htheta_mom_L0_to_pim","htheta_mom_L0_to_pim",1000,0,2000,720,0,90);
	TH2F *htheta_mom_g = new TH2F("htheta_mom_g","htheta_mom_g",1000,0,2000,720,0,90);
   TH2F *hmass_mom = new TH2F("hmass_mom","hmass_mom; p*q [MeV/c]; mass",1000,-2000,4000,1000,0,2000);

   TH1F *hmass = new TH1F("hmass","hmass",1000,-2000,2000);
	TH1F *hinvmass = new TH1F("hinvmass","hinvmass",2000,0,2000); 
	TH1F *hmass_K0 = new TH1F("hmass_K0","hmass_K0",2000,0,2000); 
   TH1F *hmass_L0 = new TH1F("hmass_L0","hmass_L0",2000,0,2000); 
	TH1F *hmass_S0 = new TH1F("hmass_S0","hmass_S0",2000,1140,1260); 
	TH1F *hmass_pi0 = new TH1F("hmass_pi0","hmass_pi0",2000,0,2000);
	TH1F *hmass_n = new TH1F("hmass_n","hmass_n",2000,0,2000); 

   TH1F *hbeta = new TH1F("hbeta","hbeta",100,0,2);
   TH1F *hg_energy = new TH1F("hg_energy","hg_energy",2000,0,2000);

	//----------------------------------Creating polygonal cut for proton identification (mass vs mom)--------------
	TCutG *cutP = new TCutG("cutP",11);
	cutP->SetPoint(0,188.0792,1187.5);
	cutP->SetPoint(1,1559.183,1378.472);
	cutP->SetPoint(2,2203.665,1517.361);
	cutP->SetPoint(3,2677.548,1579.861);
   	cutP->SetPoint(4,2886.057,1392.361);
   	cutP->SetPoint(5,2544.861,1177.083);
	cutP->SetPoint(6,2052.022,888.8889);
	cutP->SetPoint(7,188.0792,593.75);
	cutP->SetPoint(8,86.98397,642.3611);
	cutP->SetPoint(9,200.7161,1190.972);
	cutP->SetPoint(10,188.0792,1187.5);

	//----------------------------------Defining constants and lorentz vecotrs of interacting particles--------------
	const double mp = HPhysicsConstants::mass(14); // <-- p mass (PID = 14)
	const double mpim = HPhysicsConstants::mass(9); // <-- pi- mass (PID = 9)
	
	const double beam_energy = 1168.;
	const double beam_momentum = sqrt(beam_energy*beam_energy - mpim*mpim);

	TLorentzVector proj(0,0,beam_momentum,beam_energy); //beam_momentum <- ped wiazki, beam_energy <- calkowita energia (kinetyczna plus masa)
	TLorentzVector targ(0,0,0,mp);
	TLorentzVector beam(0,0,0,0);
	beam = proj + targ;

	//----------------------------------Creating lorentz vectors for daughter particles (gamma x2, p, pi+, pi-)-----
	vector<TLorentzVector> lv_neutr, lv_neutr1, lv_prot, lv_pip, lv_pim;
	
	string currentFileName;

	//--------------------------------------------------------------------------------------------------------------
	//********************************--START OF: LOOP OVER ALL EVENTS--********************************************	
	Int_t entries = loop->getEntries(); // <-- Number of entries in loop

    for(int i=0; i<entries; i++)
	{
    	if(i % 10000 == 0) {printf("Event nr.: %d\n", i);}
		loop->nextEvent(i); // <-- Get next event. Categories will be cleared before

		HTool::printProgress(i, nev, 1, "Analysis :\n");

		/*
      	TString tempFileName;
		Int_t dayOfYear;
		Int_t hour;
		Int_t min;
		Int_t second;

		if(loop->isNewFile(tempFileName))
		{
			TString type;
			Int_t year;
			Int_t eb;
			currentFileName = tempFileName;
			HTime::splitFileName(HTime::stripFileName(tempFileName), type, year, dayOfYear, hour, min, second, eb);
		}

		Float_t valueMinute = dayOfYear + hour / 24. + min / 24. / 60.;
		Float_t valueSecond = valueMinute + second / 24. / 60. / 60.;
		*/

		lv_neutr.clear();
		lv_prot.clear();
		lv_pip.clear();
		lv_pim.clear();

		/*
		HStart2Hit * fstart = nullptr;
		fstart = (HStart2Hit *) fStart2Hit->getObject(0);
		if(!fstart || fstart->getCorrFlag() == -1) continue;
		fCorrFlag=-1 iTOF
		fCorrFlag=0 only LGAD
		fCorrFlag>=0
		*/
	
		HEventHeader* event_header = NULL;
		if(!(event_header = gHades->getCurrentEvent()->getHeader())) continue;

		Int_t TBit = (Int_t)event_header->getTBit();
		Double_t VertexX = event_header->getVertexReco().getX();
		Double_t VertexY = event_header->getVertexReco().getY();
		Double_t VertexZ = event_header->getVertexReco().getZ();

		if(VertexZ<-200 || VertexZ>0) continue;
		//if(VertexZ<-148. || VertexZ>-97.) continue;

		//         Int_t DF     = (Int_t) event_header->getDownscalingFlag();
      	//         Int_t SeqNum = (Int_t) event_header->getEventSeqNumber();
      	//         Int_t TDec   = (Int_t) event_header->getTriggerDecision();
      
		//--------------------------------------------------------------------------------------------------------------
		//********************************--START OF: LOOP OVER NEUTRAL EVENTS--****************************************
		Int_t nNeutral_ev = fEmcNeutralCand->getEntries();
	
		for(int j=0; j<nNeutral_ev; j++)
		{	
			HEmcNeutralCandSim* neutr_cand = HCategoryManager::getObject(neutr_cand, fEmcNeutralCand, j);

			//Float_t dist  = neutr_cand->getDistanceToEmc();
			Int_t ind = neutr_cand->getEmcClusterIndex();

			HEmcCluster *cl=nullptr;
			cl = HCategoryManager::getObject(cl, fEmcCluster, ind);

			Int_t cl_size = cl->getNCells();
			Int_t sec = cl->getSector();
			Int_t cel = cl->getCell();		
			Int_t pid = neutr_cand->getPID();

			if(cel<33) continue;

			Double_t energy = cl->getEnergy();
			Double_t tof = cl->getTime();
			Double_t beta = neutr_cand->getBeta();
			Float_t theta = cl->getTheta();
			Float_t phi = cl->getPhi();
	
			hbeta->Fill(beta);
			hg_energy->Fill(energy);
			htheta_mom_g->Fill(energy,theta);

			TLorentzVector lvg1, lvg;
			//lvg1.SetXYZM(trackVec.getX(),trackVec.getY(),trackVec.getZ(),0);
				
			lvg1 = *neutr_cand;
			if(energy>100 && pid==1) {lv_neutr.push_back(lvg1);}	    
		}
		//********************************- END OF: LOOP OVER NEUTRAL EVENTS -******************************************
		//--------------------------------------------------------------------------------------------------------------
		
		//--------------------------------------------------------------------------------------------------------------
		//********************************- START OF: LOOP OVER CHARGED EVENTS -****************************************
		if(fParticleCand)
	   	{
			Int_t nPart_ev  = fParticleCand->getEntries();

			for (int j=0; j<nPart_ev; ++j)
			{
				HParticleCandSim* fparticlecand = HCategoryManager::getObject(fparticlecand, fParticleCand, j);
			
				Float_t theta = fparticlecand->getTheta();
				Float_t phi = fparticlecand->getPhi();
				Float_t mom = fparticlecand->getMomentum();
				Float_t beta = fparticlecand->getBeta();
				Float_t tof = fparticlecand->getTof();
				Float_t charge=fparticlecand->getCharge(); 
				Float_t mass= fparticlecand->getMass();
				//Float_t d2meta=fparticlecand->getDistanceToMetaHit();			 
				Int_t sec = fparticlecand->getSector();
				Int_t system = fparticlecand->getSystem();

				int geant_id = fparticlecand->getGeantPID();
				int geant_info1 = fparticlecand->getGeantGeninfo1();
				int parent_id=fparticlecand->getGeantParentPID();
				
				if(fparticlecand->isFlagBit(kIsUsed) && fparticlecand->getGeantParentTrackNum()==0 && fparticlecand->getGeantParentPID()==-1)
				{
					//cout<<fparticlecand->getGeantPID()<<endl;
					//cout<< i << " " << geant_id <<" "<< geant_info1   <<" "<<parent_id<<endl;
					
					hbeta_mom->Fill(charge*mom,beta);
					hmass_mom->Fill(charge*mom,mass);
					htheta_mom->Fill(mom,theta);
					hmass->Fill(charge*mass);

					//---------------------- Getting info about protons (PID = 14) --------------------------------------
					//if(charge>0 &&  cutP->IsInside(charge*mom,mass)) // <-- Polygonal cut (cutP) for protons
					if(fparticlecand->getGeantPID()==14)
					{
						Double_t deltaMom = dEdxCorr.getDeltaMom(14, mom, theta); 
						fparticlecand->setMomentum(mom + deltaMom);                                       
						fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(14));
						hbeta_mom_p->Fill(charge*mom,beta);
						htheta_mom_p->Fill(mom,theta);
			
						//part_pos.push_back(fparticlecand);
						TLorentzVector lv1a;
						lv1a = *fparticlecand;
						lv_prot.push_back(lv1a); 
					}			

					//----------------------- Getting info about p+ (PID = 8) -------------------------------------------
					//if(charge>0 &&  cutPIP->IsInside(charge*mom,mass)) <-- Polygonal cut (cutPIP) for pi+
					if(fparticlecand->getGeantPID()==8)
					{
						Double_t deltaMom = dEdxCorr.getDeltaMom(8, mom, theta); 
						fparticlecand->setMomentum(mom + deltaMom);                                     
						fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(8));
						hbeta_mom_pip->Fill(charge*mom,beta);
						htheta_mom_pip->Fill(mom,theta);	
						
						TLorentzVector lv1b;
						lv1b = *fparticlecand;
						lv_pip.push_back(lv1b);				
					}			

					//----------------------- Getting info about p- (PID = 9) -------------------------------------------
					//if(charge<0 && mass < 550){//PID pi-
					//if(charge<0 && cutPIM->IsInside(charge*mom,mass)) <-- Polygonal cut (cutPIM) for pi-
					if(fparticlecand->getGeantPID()==9)
					{
						Double_t deltaMom = dEdxCorr.getDeltaMom(9, mom, theta); 
						fparticlecand->setMomentum(mom + deltaMom);                                     
						fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(9));
						hbeta_mom_pim->Fill(charge*mom,beta);
						htheta_mom_pim->Fill(mom,theta);

						TLorentzVector lv1c;
						lv1c = *fparticlecand;
						lv_pim.push_back(lv1c);
						
						if(fparticlecand->getGeantParentPID()==16) // <-- pi- from K0S decay
							htheta_mom_K0_to_pim->Fill(mom,theta);
						if(fparticlecand->getGeantParentPID()==18) // <-- pi- from Lambda0 decay
							htheta_mom_L0_to_pim->Fill(mom,theta);	
					}			
		      	}
		   	}
	   	}
		//********************************- END OF: LOOP OVER CHARGED EVENTS -****************************************
		//------------------------------------------------------------------------------------------------------------
		
		"""--------------------------------------------------------------------------------------------------------"""
		"""*******************************- START OF: invMASS and missingMASS CALCULATIONS -***********************"""
		"""********************************- REACTION: pim + p --> K0S [pi- pi+] + SIGMA0 -************************"""
		"""--------------------------------------------------------------------------------------------------------""";

		//------------------------------- Calculating invMass of K0S -------------------------------------------------
		if(lv_pip.size() && lv_pim.size())
		{	
			for (int i=0; i<(int)(lv_pim.size()); i++)
			{
				for (int j=0; j<(int)(lv_pip.size()); j++)
				{	  	  
					TLorentzVector lv1 = lv_pim[i] + lv_pip[j];
					float invM = lv1.M();
					hmass_K0->Fill(invM);

					//---------------------- Calculating missing mass (Sigma0) ------------------------------------------
					if(invM>480 && invM<510) // <-- Check if invMass is equal to K0S'
					{
						TLorentzVector lv2 = beam - lv_pim[i] - lv_pip[j];
						float mm = lv2.M(); // mm - missing mass [MeV]
						float mm2 = lv2.M2(); // mm2 - missing mass squared [MeV^2]
						hmass_S0->Fill(mm);

						for(int k=0; k<(int)(lv_pim.size()); k++)
						{
							//---------------- Calculating invMass of Lambda0 [p pi-] -------------------------------------
							if(i==k || !(lv_prot.size())) continue; // <-- Using only remaining pi- if there exists a proton
							TLorentzVector lv3 = lv_pim[k] + lv_prot[0];
						}					
					}
				}
			}
		}

		//------------------------------- Calculating invMass of Lambda0 [p pi-]------------------------------------
		if(lv_prot.size() && lv_pim.size())
		{
			for(int i=0; i<(int)(lv_prot.size()); i++)
			{
				for(int j=0; j<(int)(lv_pim.size()); j++)
				{
					TLorentzVector lv1 = lv_prot[i] + lv_pim[j];
					float invM = lv1.M();
					hmass_L0->Fill(invM);

					//if(invM > 1110 && invM < 1120)
					//{
					// 	for(int k=0; k<lv_pim.size(); k++)
					// 	{
					// 	if(k!=j) continue;
					// 	TLorentzVector lv1 = beam - lv_prot[i] - lv_pim[j];
					// 	float mm = lv1.M(); // mm - missing mass [MeV]
					// 	float mm2 = lv1.M2(); // mm2 - missing mass squared [MeV^2]
					// 	hinvmass->Fill(mm);
					// 	}
					//}
				}
			}
		}

		//------------------------------- Calculating invMass of pi0 [gamma gamma] -----------------------------------
		if((int)lv_neutr.size()>=2)
		{
			for(int i=0; i<(int)(lv_neutr.size()); i++)
			{
				for(int j=0; j<(int)(lv_neutr.size()); j++)
				{
					if(i==j) continue;
					TLorentzVector lv1 = lv_neutr[i] + lv_neutr[j];
					float invM = lv1.M();
					hmass_pi0->Fill(invM);
				}
			}
		}
		"""*******************************- END OF: invMASS and missingMASS CALCULATIONS -*************************"""
		"""--------------------------------------------------------------------------------------------------------""";
	}
	//********************************--END OF: LOOP OVER ALL EVENTS--**********************************************
	//--------------------------------------------------------------------------------------------------------------

	output_file->cd();

	hbeta_mom->Write();
	hmass_mom->Write();
	hmass->Write();
	hbeta_mom_p->Write();
	hbeta_mom_pip->Write();
	hbeta_mom_pim->Write();
	htheta_mom_g->Write();
	htheta_mom->Write(); // <-- Theta vs Momentum (except for photons)
	htheta_mom_p->Write();
	htheta_mom_pip->Write();
	htheta_mom_pim->Write();
	htheta_mom_K0_to_pim->Write();
	htheta_mom_L0_to_pim->Write();
	hinvmass->Write();
	hmass_K0->Write();
	hmass_L0->Write();
	hmass_S0->Write();
	hmass_pi0->Write();
	hmass_n->Write();

	output_file->Close();
	cout << "writing root tree done" << endl;

	timer.Stop();
	timer.Print();

	return 0;
}
