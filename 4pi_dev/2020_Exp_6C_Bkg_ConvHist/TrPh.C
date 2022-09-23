#include <cmath>
#include <iostream>
#include <set>
#include "TrPh.h"
#include "TMatrixDSym.h"
#include "TVector.h"
#include "Hypo2ChPions2pi0.hpp"
#include "Hypo2ChPionsPi02Photons.hpp"
#include "Hypo2ChPions2PhotonsLostPi0.hpp"
#define Mpi 139.57
#define mpiz 134.9766

TrPh::TrPh(TTree* tree) : kfcmd::core::TrPh(tree) {}

TrPh::~TrPh() {}

//int callLoop(std::string name_in = "/store11/eakozyrev/4pi/toy_july_2022/raw/2pi2pi0/tr_ph_run000107.root", std::string name_out = "Out_107_5x5.root"){
int callLoop(std::string name_in = "root://sl10cmd//scan2020/scan2020_tr_ph_fc_e935_v8.root", std::string name_out = "Out_scan2020_tr_ph_fc_e935_v8.root"){
  TFile *oldfile;
  TTree *oldtree;

  oldfile = TFile::Open(name_in.c_str());
  oldtree = (TTree*)oldfile->Get("tr_ph");
    
  TrPh myTrPh(oldtree);
  myTrPh.Loop(name_out, 1.0);
 
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

    const std::set<std::string> s_g12 = {"g1", "g2"};
    const std::set<std::string> s_g34 = {"g3", "g4"};

    Long64_t nentries = fChain->GetEntriesFast();
   
    TFile newfile(outpath.c_str(),"recreate");
    TTree *newtree = fChain->GetTree()->CloneTree(0);
    TTree *tree_err_5C = new TTree("tree_err_5C", "tree_err_5C");
    TTree *tree_err_6C = new TTree("tree_err_6C", "tree_err_6C");
    TTree *tree_err_miss = new TTree("tree_err_miss", "tree_err_miss");

    TLorentzVector Lpip,Lpim,Lph1,Lph2,Lph3,Lph4;
    double Chi25C, Chi26C, m2g0, m2g;
    double m2pirest, E_poln, P_poln;
    int mode1, mode2;
    int Err_Chi25C, Err_Chi26C, Err_Chi2Miss;

    double mass, misspi0mass, m2g_miss, m2g0_miss;
    double theta_miss_pi, mom_miss_pi;
    double Xi2misspi0;
    TLorentzVector Ltr1, Ltr2, Lgamma1, Lgamma2;

    int Ngen = nentries;
    int correspond[20];
    newtree->Branch("Ngen",&Ngen);
    newtree->Branch("Chi25C",&Chi25C, "Chi25C/D");
    newtree->Branch("Chi26C",&Chi26C);
    newtree->Branch("Lpim",&Lpim);
    newtree->Branch("Lpip",&Lpip);
    newtree->Branch("Lph1",&Lph1);
    newtree->Branch("Lph2",&Lph2);
    newtree->Branch("Lph3",&Lph3);
    newtree->Branch("Lph4",&Lph4);
    newtree->Branch("E_poln",&E_poln);
    newtree->Branch("P_poln",&P_poln);
    newtree->Branch("m2g0",&m2g0);
    newtree->Branch("m2g",&m2g);
    newtree->Branch("mode1",&mode1, "mode1/I");
    newtree->Branch("mode2",&mode2, "mode2/I");
    newtree->Branch("correspond",correspond,"correspond[20]/I");
    newtree->Branch("m2pirest",&m2pirest);

    newtree->Branch("m2g0_miss",&m2g0_miss);
    newtree->Branch("m2g_miss",&m2g_miss);
    newtree->Branch("Ltr1",&Ltr1);
    newtree->Branch("Ltr2",&Ltr2);
    newtree->Branch("Lgamma1",&Lgamma1);
    newtree->Branch("Lgamma2",&Lgamma2);
    newtree->Branch("Xi2misspi0",&Xi2misspi0);
    newtree->Branch("mom_miss_pi",&mom_miss_pi);
    newtree->Branch("theta_miss_pi",&theta_miss_pi);


    tree_err_5C->Branch("Err_Chi25C",&Err_Chi25C, "Err_Chi25C/I");
    tree_err_6C->Branch("Err_Chi26C",&Err_Chi26C, "Err_Chi26C/I");
    tree_err_miss->Branch("Err_Chi2Miss",&Err_Chi2Miss, "Err_Chi2Miss/I");
    cout << nentries << endl;
     
    //variables to draw "fit convergence histo"
    int flag_5CCut, flag_err5C_0, flag_err5C_1, flag_err5C_2;
    int flag_6CCut, flag_err6C_0, flag_err6C_1, flag_err6C_2;
    int flag_MissCut, flag_errMiss_0, flag_errMiss_1, flag_errMiss_2;

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
    fChain->GetEntry(0);
    Hypo2ChPionsPi02Photons hypo_5C(2.e-3 * emeas, magneticField);
    Hypo2ChPions2pi0 hypo_6C(2.e-3 * emeas, magneticField);
    Hypo2ChPions2PhotonsLostPi0 hypo_miss(2.e-3 * emeas, magneticField);
    int runnum_switch = 90250;
    int flag_m_switch = 0;
    for(Long64_t jentry=2400000; jentry<2600000; jentry++){

        //Read entry:
    	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
    	nb = fChain->GetEntry(jentry);   nbytes += nb;
 	
        if(finalstate_id == 3)continue;
    	if(jentry%50000 == 0)  std::cout << (double)jentry/nentries << "  "  << std::endl;

	if(fabs(ebeam - emeas)>50)emeas = ebeam;
	hypo_5C.getParticle("origin")->fixParameter(3, 2.e-3*emeas);
	hypo_6C.getParticle("origin")->fixParameter(3, 2.e-3*emeas);
	hypo_miss.getParticle("origin")->fixParameter(3, 2.e-3*emeas);

	hypo_5C.setBeamXY(xbeam, ybeam);
	hypo_5C.fixVertexParameter("vtx0", 0, xbeam);
        hypo_5C.fixVertexParameter("vtx0", 1, ybeam);
	
	if (runnum >= runnum_switch && flag_m_switch == 0){
		flag_m_switch = 1;
		hypo_5C.addConstant("#m-field", 1.3);
		hypo_6C.addConstant("#m-field", 1.3);
		hypo_miss.addConstant("#m-field", 1.3);
	}

	flag_5CCut = 0;
	flag_6CCut = 0;
	flag_MissCut = 0;
        flag_err5C_0 = 0;
	flag_err5C_1 = 0;
	flag_err5C_2 = 0;
        flag_err6C_0 = 0;
	flag_err6C_1 = 0;
	flag_err6C_2 = 0;
        flag_errMiss_0 = 0;
	flag_errMiss_1 = 0;
	flag_errMiss_2 = 0;
	Err_Chi25C = -10;
	Err_Chi26C = -10;
	Err_Chi2Miss = -10;


	
	//
	//std::cout << "M field = " << hypo_5C.getConstant("#m-field") << std::endl;
       
	
	//std::cout << "ENERGY = " <<  hypo_5C.getParticle("origin")->getBeginParameters()(3) << std::endl;
	
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
       mode1 = 0;
       mode2 = 0;
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
			    

			    if (!hypo_5C.fillTrack("pi+", one, *this))continue;
			    if (!hypo_5C.fillTrack("pi-", two, *this))continue;
			    if (!hypo_5C.fillPhoton("g1", k, *this))continue;
			    if (!hypo_5C.fillPhoton("g2", l, *this))continue;
			    if (!hypo_5C.fillPhoton("g3", m, *this))continue;
			    if (!hypo_5C.fillPhoton("g4", n, *this))continue;
			   
			    //hypo.updateInitialParams();
			    


			    hypo_5C.optimize();
			    //std::cout << "ENERGY = " <<  hypo_5C.getParticle("origin")->getFinalParameters()(3) << std::endl;
			    int err = hypo_5C.getErrorCode();
			    flag_5CCut = 1;
			    if (err == 0){
				flag_err5C_0 = 1;
			    }
                            if (err == 1){
				flag_err5C_1 = 1;
			    }
			    if (err == 2){
			 	flag_err5C_2 = 1;
			    }
			    
			    
                            
			    if (err!=0) continue;
 			    if(hypo_5C.getChiSquare() < Chi25C || firsttry == true){
			      //std::cout << "jentry = " << jentry << std::endl;	
			      if(hypo_5C.getFinalMomentum(s_g34).M()*1.e3 < Cut_mpi0_before  || hypo_5C.getFinalMomentum(s_g34).M()*1.e3 > Cut_mpi0_above)continue;
			      firsttry = false;	
			      Chi25C = hypo_5C.getChiSquare();
			      
			      m2g = hypo_5C.getFinalMomentum(s_g34).M()*1.e3;	    	    
			      m2g0 = hypo_5C.getInitialMomentum(s_g34).M()*1.e3;

			      //Initialize 6C Hypothesis and do 6C fit
			      hypo_6C.setBeamXY(xbeam, ybeam);
        		      hypo_6C.fixVertexParameter("vtx0", 0, xbeam);
      			      hypo_6C.fixVertexParameter("vtx0", 1, ybeam);
			      if (!hypo_6C.fillTrack("pi+", one, *this))continue;
			      if (!hypo_6C.fillTrack("pi-", two, *this))continue;
			      if (!hypo_6C.fillPhoton("g1", k, *this))continue;
			      if (!hypo_6C.fillPhoton("g2", l, *this))continue;
			      if (!hypo_6C.fillPhoton("g3", m, *this))continue;
			      if (!hypo_6C.fillPhoton("g4", n, *this))continue;


			      hypo_6C.optimize();
			      Chi26C = hypo_6C.getChiSquare();

			      flag_6CCut = 1;
			      err = hypo_6C.getErrorCode();
			      if (err == 0){
				flag_err6C_0 = 1;
			      }
                              if (err == 1){
				flag_err6C_1 = 1;
			      }
			      if (err == 2){
			 	flag_err6C_2 = 1;
			      }
			     
			      Lpip = hypo_6C.getFinalMomentum("pi+")*1.e3;
			      Lpim = hypo_6C.getFinalMomentum("pi-")*1.e3; 
			      Lph1 = hypo_6C.getFinalMomentum("g1")*1.e3;
			      Lph2 = hypo_6C.getFinalMomentum("g2")*1.e3;
			      Lph3 = hypo_6C.getFinalMomentum("g3")*1.e3;
			      Lph4 = hypo_6C.getFinalMomentum("g4")*1.e3;
		      
			      m2pirest = (Pebeam - Lpip - Lpim).M();
			      correspond[0] = one;
			      correspond[1] = two;  
			      correspond[2] = k;
			      correspond[3] = l;	
			      correspond[4] = m;
			      correspond[5] = n;
			      E_poln = fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).E() - Pebeam.E());
			      P_poln = fabs((PTr[one] + PTr[two] + PPh[k] + PPh[l]+ PPh[m]+ PPh[n]).P());
			      mode1 = 1;
			     
			      N_events++; 

			    }


		          }
		  }
      }
      if (flag_5CCut == 1){
      	if (flag_err5C_0 == 1){
		Err_Chi25C = 0;
      	}
      	else if (flag_err5C_1 == 1){
		Err_Chi25C = 1;
      	}
      	else if (flag_err5C_2 == 1){
		Err_Chi25C = 2;
      	}
        tree_err_5C->Fill();
      }

      if (flag_6CCut == 1){
      	if (flag_err6C_0 == 1){
		Err_Chi26C = 0;
      	}
      	else if (flag_err6C_1 == 1){
		Err_Chi26C = 1;
      	}
      	else if (flag_err6C_2 == 1){
		Err_Chi26C = 2;
      	}
	tree_err_6C->Fill();
      }

      

      // end end end end end end end end end end 
      //      2pi0 pi+ pi- hypothesis
      // end end end end end end end end end end





      //***************************************** 
      //   1pi0 pi+ pi- hypothesis 
      //*****************************************
      mass = 300.*ebeam/500.;
      Xi2misspi0 = 1000.;
      for(int one = 0; one < nt-1; one++){
	  	  for(int two = one + 1; two < nt; two++){
			  if(tcharge[one]*tcharge[two] > 0){continue;}
			  if(nt > 4)continue;
			  if(nph > 8)continue;
	  	  	  if(fabs(tz[one]) > 12. || fabs(trho[one]) > 0.1 || tdedx[one] > tdedx_high(tptot[one]) || tdedx[one] < tdedx_low(tptot[one])){continue;}
	  	  	  if(fabs(tz[two]) > 12. || fabs(trho[two]) > 0.1 || tdedx[two] > tdedx_high(tptot[two]) || tdedx[two] < tdedx_low(tptot[two])){continue;} 
	  	  	  for(int l = 0; l < npi0_alone; l++){
				  int m = numbers_alone[l][0];
	  	  	  	  int n = numbers_alone[l][1];
				  if(fabs((Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n]).M() - Mpi) < mass){
					mass = fabs((Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n]).M() - Mpi);
					misspi0mass = (Pebeam - PTr[one] - PTr[two] - PPh[m] - PPh[n]).M();
				      	
					hypo_miss.setBeamXY(xbeam, ybeam);
        				hypo_miss.fixVertexParameter("vtx0", 0, xbeam);
       					hypo_miss.fixVertexParameter("vtx0", 1, ybeam);

			      		if (!hypo_miss.fillTrack("pi+", one, *this))continue;
			      		if (!hypo_miss.fillTrack("pi-", two, *this))continue;
			      		if (!hypo_miss.fillPhoton("g1", m, *this))continue;
			      		if (!hypo_miss.fillPhoton("g2", n, *this))continue;


					hypo_miss.optimize();
				        int err = hypo_miss.getErrorCode();

					flag_MissCut = 1;
					if (err == 0){
						flag_errMiss_0 = 1;
			      		}
                              		if (err == 1){
						flag_errMiss_1 = 1;
			      		}
			      		if (err == 2){
			 			flag_errMiss_2 = 1;
			      		}
			     
					if (err != 0)continue;
					if (hypo_miss.getChiSquare() < Xi2misspi0){
						if((hypo_miss.getFinalMomentum(s_g12)).M()*1.e3 < Cut_mpi0_before || (hypo_miss.getFinalMomentum(s_g12)).M()*1.e3   > Cut_mpi0_above)continue;
						if((hypo_miss.getFinalMomentum(s_g12)).P()*1.e3 < 150. || (hypo_miss.getFinalMomentum(s_g12)).P()*1.e3 > 650.)continue;
						
						m2g0_miss = (hypo_miss.getInitialMomentum(s_g12)).M()*1.e3;
						m2g_miss = (hypo_miss.getFinalMomentum(s_g12)).M()*1.e3;
						Xi2misspi0 = hypo_miss.getChiSquare();

					  	correspond[11] = m;
					  	correspond[12] = n;
					  	correspond[13] = one;
					  	correspond[14] = two;
					  	Ltr1 = hypo_miss.getFinalMomentum("pi+")*1.e3;
					  	Ltr2 = hypo_miss.getFinalMomentum("pi-")*1.e3;
					  	Lgamma1 = hypo_miss.getFinalMomentum("g1")*1.e3;
					  	Lgamma2 = hypo_miss.getFinalMomentum("g2")*1.e3;
					  	mode2 = 1;
					  	theta_miss_pi = (Pebeam - Ltr1 - Ltr2 - Lgamma1 - Lgamma2).Theta();
					  	mom_miss_pi = (Pebeam - Ltr1 - Ltr2 - Lgamma1 - Lgamma2).P();
					}
				  }
			  }
		  }
      }
      

      if (flag_MissCut == 1){
      	if (flag_errMiss_0 == 1){
		Err_Chi2Miss = 0;
      	}
      	else if (flag_errMiss_1 == 1){
		Err_Chi2Miss = 1;
      	}
      	else if (flag_errMiss_2 == 1){
		Err_Chi2Miss = 2;
      	}
	tree_err_miss->Fill();
      }
       // end end end end end end end end end end end end
       //    1 pi0 pi+ pi- hypothesis
       // end end end end end end end end end end end end
 


      if((mode1 + mode2) > 0){newtree->Fill();}
   }


    cout << "OK1  " << outpath << " " << (double)newtree->GetEntries()/(double)nentries << endl;
    //cout << "N_events = " << N_events << "    " << (double)newtree->GetEntries() << endl;
    newtree->Write();
    tree_err_5C->Write();
    tree_err_6C->Write();
    tree_err_miss->Write();
    delete newtree;
    delete tree_err_5C;
    delete tree_err_6C;
    delete tree_err_miss;
    newfile.Write();
    newfile.Close();




   
    
 
}
