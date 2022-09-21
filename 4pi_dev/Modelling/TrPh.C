#include <cmath>
#include <iostream>
#include <set>
#include "TrPh.h"
#include "TMatrixDSym.h"
#include "TVector.h"
#include "Hypo2ChPionsPi02Photons.hpp"
#define Mpi 139.57
#define mpiz 134.9766

TrPh::TrPh(TTree* tree) : kfcmd::core::TrPh(tree) {}

TrPh::~TrPh() {}

int callLoop(std::string name_in = "/store11/eakozyrev/4pi/toy_july_2022/raw/2pi2pi0/tr_ph_run000107.root", std::string name_out = "Out_107_5x5.root"){
//int callLoop(std::string name_in = "root://sl10cmd//scan2020/scan2020_tr_ph_fc_e935_v8.root", std::string name_out = "Out_107_5x5.root"){
  TFile *oldfile;
  TTree *oldtree;

  oldfile = TFile::Open(name_in.c_str());
  oldtree = (TTree*)oldfile->Get("tr_ph");
    
  TrPh myTrPh(oldtree);
  myTrPh.Loop(name_out, 1.3);
 
  return 0;
}






double tdedx_low(double x){
    return 0.3*(0.5*(4200.*(1+TMath::Exp(-0.0095*(x-110.))+TMath::Exp(-0.024*(x-110.)))+500)+200.);
}

double tdedx_high(double x){
    return 0.3*(0.5*(4200.*(1+TMath::Exp(-0.0095*(x-190.))+TMath::Exp(-0.024*(x-190.)))+4500)+8500.);
}




Int_t TrPh::Cut(Long64_t entry) {
  return 1;
}





void TrPh::Loop(const std::string& outpath, double magneticField) {
    if (fChain == 0) return;

    const std::set<std::string> sPi0 = {"g1", "g2"};
    const std::set<std::string> sGG = {"g3", "g4"};

    Long64_t nentries = fChain->GetEntriesFast();
   
    TFile newfile(outpath.c_str(),"recreate");
    TTree *newtree = fChain->GetTree()->CloneTree(0);

    TLorentzVector Lpip,Lpim,Lph1,Lph2,Lph3,Lph4;
    double Chi25C,m2g0,m2g;
    int mode;
    double Err_Chi25C;

    int Ngen = nentries;
    int correspond[20];
    newtree->Branch("Ngen",&Ngen);
    newtree->Branch("Chi25C",&Chi25C, "Chi25C/D");
    newtree->Branch("Err_Chi25C",&Err_Chi25C, "Err_Chi25C/D");
    newtree->Branch("Lpim",&Lpim);
    newtree->Branch("Lpip",&Lpip);
    newtree->Branch("Lph1",&Lph1);
    newtree->Branch("Lph2",&Lph2);
    newtree->Branch("Lph3",&Lph3);
    newtree->Branch("Lph4",&Lph4);
    newtree->Branch("m2g0",&m2g0);
    newtree->Branch("m2g",&m2g);
    newtree->Branch("mode",&mode, "mode/I");
    newtree->Branch("correspond",correspond,"correspond[20]/I");
    cout << nentries << endl;
    
    // window for mpi0
    double Cut_mpi0_before = 45; double Cut_mpi0_above = 255;
    double th_ph_min = 0.7;double th_ph_max = 3.1415 - th_ph_min;
    double th_tr_min = 0.8;double th_tr_max = 3.1415 - th_tr_min;
    double ephmin = 30.;

   
    int numbers[100000][4],numbers_alone[100000][2];
    int npi0 = 0;int npi0_alone = 0;
    TLorentzVector PPh[1000];TLorentzVector PTr[1000];

    Long64_t nbytes = 0, nb = 0;
    double N_events = 0;
    std::cout << "SizeOf terr0 = "  << sizeof(terr0[0])/sizeof(terr0[0][0][0]) << std::endl;
    fChain->GetEntry(0)
    Hypo2ChPionsPi02Photons hypo(2.e-3 * emeas, magneticField);
    for(Long64_t jentry=0; jentry<nentries; jentry++){
	
        //Read entry:
    	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
    	nb = fChain->GetEntry(jentry);   nbytes += nb;
 	 
        if(finalstate_id == 3)continue;
    	if(jentry%5000 == 0)  std::cout << (double)jentry/nentries << "  "  << std::endl;
	
	if(fabs(ebeam - emeas)>50)emeas = ebeam;
    	
	
        hypo.setBeamXY(xbeam, ybeam);
        hypo.fixVertexParameter("vtx0", 0, xbeam);
        hypo.fixVertexParameter("vtx0", 1, ybeam);

	//initialize some vars
	Chi25C = 500;
        TLorentzVector Pebeam(0.,0.,0.,2.*emeas);
        m2g = 0;
        m2g0 = 0;
        npi0 = 0;
        npi0_alone = 0;
        bool firsttry = true;

        for(int k = 0; k < nph; k++){
        	PPh[k] = TLorentzVector(phen[k]*TMath::Cos(phphi[k])*TMath::Sin(phth[k]),
				phen[k]*TMath::Sin(phphi[k])*TMath::Sin(phth[k]),
				phen[k]*TMath::Cos(phth[k]),
				phen[k]);
      } 
      // tracks 4-momenta
        for(int k = 0; k < nt; k++){
		PTr[k] = TLorentzVector(tptot[k]*TMath::Cos(tphi[k])*TMath::Sin(tth[k]),
				tptot[k]*TMath::Sin(tphi[k])*TMath::Sin(tth[k]),
				tptot[k]*TMath::Cos(tth[k]),
				sqrt(tptot[k]*tptot[k]+Mpi*Mpi));
        }
        // the list of pi0
        // npi0_alone - the number of pi0 in an event
        for(int i = 0; i < nph-1; i++){
	  for(int j = i + 1; j < nph; j++){
	    if((PPh[i] + PPh[j]). M() < Cut_mpi0_above && (PPh[i] + PPh[j]). M() > Cut_mpi0_before){
	      if(phth[i] < th_ph_min || phth[i] > th_ph_max || phen[i] < ephmin)continue;
	      if(phth[j] < th_ph_min || phth[j] > th_ph_max || phen[j] < ephmin)continue;
	      if((phlxe[i] + phbgo[i]) < 15)continue;
	      if((phlxe[j] + phbgo[j]) < 15)continue;
	    numbers_alone[npi0_alone][0] = i; numbers_alone[npi0_alone][1] = j; npi0_alone++;
	    } 		
	  }
        }

       for(int i=0; i<20; i++){
       	correspond[i] = -1;
       }
       mode = 0;
       double mass = 30;


      //-----------------------------REQUIRE----------------------------------------//
       if(nt < 2 || npi0_alone < 1)continue;         

       for(int i = 0; i < nph-3; i++){
	for(int j = i + 1; j < nph-2; j++){
	  for(int k = j +1; k < nph-1; k++){
	    for(int l = k + 1; l < nph; l++){
	      if(phth[i] < th_ph_min || phth[i] > th_ph_max || phen[i] < ephmin)continue;
	      if(phth[j] < th_ph_min || phth[j] > th_ph_max || phen[j] < ephmin)continue;	      
	      if(phth[k] < th_ph_min || phth[k] > th_ph_max || phen[k] < ephmin)continue;
	      if(phth[l] < th_ph_min || phth[l] > th_ph_max || phen[l] < ephmin)continue;
	      if((phlxe[i] + phbgo[i]) < 15)continue;
	      if((phlxe[j] + phbgo[j]) < 15)continue;
	      if((phlxe[k] + phbgo[k]) < 15)continue;
	      if((phlxe[l] + phbgo[l]) < 15)continue;
	      if((PPh[i] + PPh[j]). M() < Cut_mpi0_above && (PPh[i] + PPh[j]). M() > Cut_mpi0_before
		 && (PPh[k] + PPh[l]). M() < Cut_mpi0_above && (PPh[k] + PPh[l]). M() > Cut_mpi0_before){
		numbers[npi0][0] = i;numbers[npi0][1] = j;numbers[npi0][2] = k;numbers[npi0][3] = l; npi0++;
	      }
	      if((PPh[i] + PPh[k]). M() < Cut_mpi0_above && (PPh[i] + PPh[k]). M() > Cut_mpi0_before
		 && (PPh[j] + PPh[l]). M() < Cut_mpi0_above && (PPh[j] + PPh[l]). M() > Cut_mpi0_before){
		numbers[npi0][0] = i;numbers[npi0][1] = k;numbers[npi0][2] = j;numbers[npi0][3] = l; npi0++;
	      }
	      if((PPh[i] + PPh[l]). M() < Cut_mpi0_above && (PPh[i] + PPh[l]). M() > Cut_mpi0_before
		 && (PPh[k] + PPh[j]). M() < Cut_mpi0_above && (PPh[k] + PPh[j]). M() > Cut_mpi0_before){
		numbers[npi0][0] = i;numbers[npi0][1] = l;numbers[npi0][2] = k;numbers[npi0][3] = j; npi0++;
	      }
	    }
	 }
        }
       }
     
      
     //***************************************** 
      //   2pi0 pi+ pi- hypothesis 
      //*****************************************
      Chi25C = 3000.;
      Err_Chi25C=10;
      for(int one = 0; one < nt-1; one++){
	  	  for(int two = one + 1; two < nt; two++){
		    if(tcharge[one]*tcharge[two] > 0)continue;
		    if(fabs(tz[one]) > 12. || fabs(trho[one]) > 0.1 || tdedx[one] > tdedx_high(tptot[one]) || tdedx[one] < tdedx_low(tptot[one])){continue;}
		    if(fabs(tz[two]) > 12. || fabs(trho[two]) > 0.1 || tdedx[two] > tdedx_high(tptot[two]) || tdedx[two] < tdedx_low(tptot[two])){continue;}
	  	  	  for(int j = 0; j < npi0; j++){
			    int k = numbers[j][0];
			    int l = numbers[j][1];
			    int m = numbers[j][2];
			    int n = numbers[j][3];


			    
			    if(fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).E() - Pebeam.E()) > 500){continue;}
			    if(fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).P()) > 500){continue;}
			    //std::cout << "jentry = " << jentry << std::endl;			    


			    if (!hypo.fillTrack("pi+", one, *this)){
				continue;
			    }
			    if (!hypo.fillTrack("pi-", two, *this)){
				continue;
			    }
			    if (!hypo.fillPhoton("g1", k, *this)){
				continue;
			    }
			    if (!hypo.fillPhoton("g2", l, *this)){
				continue;
			    }
			    if (!hypo.fillPhoton("g3", m, *this)){
				continue;
			    }
			    if (!hypo.fillPhoton("g4", n, *this)){
				continue;
			    }
			   
			    //hypo.updateInitialParams();
			    


			    hypo.optimize();
			    int err=hypo.getErrorCode();
			    //std::cout << "jentry  = " << jentry << "   err = " << err << std::endl;
			    if (err!=0) continue;
 			    if(hypo.getChiSquare() < Chi25C || firsttry == true){
			      //std::cout << "jentry = " << jentry << std::endl;	
			      if(hypo.getFinalMomentum(sGG).M()*1.e3 < Cut_mpi0_before  || hypo.getFinalMomentum(sGG).M()*1.e3 > Cut_mpi0_above)continue;
			      Err_Chi25C = err;
			      firsttry = false;	
			      Chi25C = hypo.getChiSquare();
			      if (Chi25C < 0) std::cout << "Chi25c  = " << Chi25C << std::endl;
			      m2g = hypo.getFinalMomentum(sGG).M()*1.e3;	    	    
			      m2g0 = hypo.getInitialMomentum(sGG).M()*1.e3;

			      Lpip = hypo.getFinalMomentum("pi+")*1.e3;
			      Lpim = hypo.getFinalMomentum("pi-")*1.e3; 
			      Lph1 = hypo.getFinalMomentum("g1")*1.e3;
			      Lph2 = hypo.getFinalMomentum("g2")*1.e3;
			      Lph3 = hypo.getFinalMomentum("g3")*1.e3;
			      Lph4 = hypo.getFinalMomentum("g4")*1.e3;
         
			      
			      correspond[0] = one;
			      correspond[1] = two;  
			      correspond[2] = k;
			      correspond[3] = l;	
			      correspond[4] = m;
			      correspond[5] = n;
			    
			      mode = 1;
			      N_events++;

			    }


		          }
		  }
      }
   
      if(mode  > 0){newtree->Fill();}

   }
    

    cout << "OK1  " << outpath << " " << (double)newtree->GetEntries()/(double)nentries << endl;
    //cout << "N_events = " << N_events << "    " << (double)newtree->GetEntries() << endl;
    newtree->Write();
    delete newtree;
    newfile.Write();
    newfile.Close();

   
    
 
}
