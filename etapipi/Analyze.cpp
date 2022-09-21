TFile my_file("Out.root", "READ");


TTree* mytree = (TTree*)my_file.Get("tr_ph;5");
int Nentries = mytree->GetEntries();
printf("N = %d\n", Nentries);

mytree->SetBranchAddress("Lpip", &Pip);
mytree->SetBranchAddress("Lph1", &Pph1);
mytree->SetBranchAddress("m2g", &m2g);
mytree->SetBranchAddress("m2g0", &m2g0);
mytree->SetBranchAddress("Chi25C", &Chi25C);

mytree->SetBranchAddress("correspond", correspond);

int counter_over = 0;
int counter_all = 0;

TH1D hist_pPl("pPl", "pPl", 100, 0, 2);
TH1D hist_pPh("pPh", "pPh", 100, 0, 2);
TH1D hist_m2g("M2g", "M2g", 100, 0.03, 0.3);
TH1D hist_m2g0("M2g0", "M2g0", 100, 0.03, 0.3);
TH1D hist_Chi25C("Chi25C", "Chi25C", 1000, 0, 1000);

for(int i=0; i<Nentries; i++){
	mytree->GetEntry(i);
	if (Chi25C > 20){
		//continue;
	}
	hist_pPl.Fill(Pip->P()/1000.);
	hist_pPh.Fill(Pph1->P()/1000.);
	hist_m2g.Fill(m2g/1000.);
	hist_m2g0.Fill(m2g0/1000.);
	hist_Chi25C.Fill(Chi25C);

	if (Pip->E() > 2000 and correspond[0] != -1){
		counter_over++;	
	}

}
printf("counter overflow = %d\n", counter_over);

TCanvas c1("c1", "", 1000, 1000); 
hist_pPl.Draw();

TCanvas c2("c2", "", 1000, 1000);
hist_pPh.Draw();

TCanvas c3("c3", "", 1000, 1000);
hist_m2g.SetLineColor(kRed);
hist_m2g0.Draw();
hist_m2g.SetLineWidth(3);
hist_m2g0.SetLineWidth(3);
hist_m2g.Draw("same");
TLegend legend1(0, 1, 0.1, 0.9);
legend1.AddEntry("M2g", "M2g", "lep");
legend1.AddEntry("M2g0", "M2g 0", "lep");
legend1.Draw();
gPad->SetLogy();


TCanvas c4("c4", "", 1000, 1000);
hist_Chi25C.Draw();
hist_Chi25C.SetLineWidth(3);
 TFile f_res("Results.root", "RECREATE");
 c3.Write();
 f_res.Close();

