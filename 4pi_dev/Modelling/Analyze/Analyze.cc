{


TFile my_file("Out_107_5x5.root", "READ");
TLorentzVector* Pip;
TLorentzVector* Pim;
TLorentzVector* Pph1;
TLorentzVector* Pph2;
double m2g, m2g0;
int chi2_err;
double Chi25C;
int correspond[20] = {0};

TTree* mytree = (TTree*)my_file.Get("tr_ph");
int Nentries = mytree->GetEntries();
printf("N = %d\n", Nentries);

mytree->SetBranchAddress("Lpip", &Pip);
mytree->SetBranchAddress("Lph1", &Pph1);
mytree->SetBranchAddress("m2g", &m2g);
mytree->SetBranchAddress("m2g0", &m2g0);
mytree->SetBranchAddress("Chi25C", &Chi25C);
mytree->SetBranchAddress("Err_Chi25C", &chi2_err);

mytree->SetBranchAddress("correspond", correspond);

int counter_over = 0;
int counter_all = 0;

TH1D hist_pPl("pPl", "pPl", 100, 0, 2);
TH1D hist_pPh("pPh", "pPh", 100, 0, 2);
TH1D hist_m2g("M2g", "M2g", 200, 0.03, 0.3);
TH1D hist_m2g0("M2g0", "M2g0", 200, 0.03, 0.3);
TH1D hist_Chi25C("Chi25C", "Chi25C", 2000, -1000, 1000);
TH1I hist_Chi25CErr("Chi25C_err", "Chi25C_err", 3, 0, 3);

for(int i=0; i<Nentries; i++){
	mytree->GetEntry(i);
	if (Chi25C > 0){
		hist_pPl.Fill(Pip->P()/1000.);
		hist_pPh.Fill(Pph1->P()/1000.);
		hist_m2g.Fill(m2g/1000.);
		hist_m2g0.Fill(m2g0/1000.);
	}
	hist_Chi25C.Fill(Chi25C);
	hist_Chi25CErr.Fill(chi2_err);
	if (Pip->E() > 2000 and correspond[0] != -1){
		counter_over++;	
	}

}
printf("counter overflow = %d\n", counter_over);

TCanvas c1("c1", "", 1000, 1000); 
hist_pPl.Draw();

TCanvas c2("c2", "", 1000, 1000);
hist_pPh.Draw();

TCanvas c3("1", "1", 1000, 1000);
hist_m2g.SetLineColor(kRed);
hist_m2g.GetXaxis()->SetRangeUser(0.11, 0.16);
hist_m2g.SetLineWidth(3);
hist_m2g0.SetLineWidth(3);
hist_m2g.SetTitle("0 < Chi25C < 50");
hist_m2g.Draw();
hist_m2g0.Draw("same");
gPad->SetGrid();
//gPad->SetLogy();

TLegend legend1(0, 1, 0.1, 0.9);
legend1.AddEntry("M2g", "M2g", "lep");
legend1.AddEntry("M2g0", "M2g0", "lep");
legend1.Draw();
c3.Draw();
TCanvas c4("c4", "", 1000, 1000);
hist_Chi25C.SetLineWidth(3);
hist_Chi25C.Draw();

TCanvas c5("c5", "", 1000, 1000);
hist_Chi25CErr.SetLineWidth(3);
hist_Chi25CErr.Draw();


 TFile f_res("Results.root", "RECREATE");
 c3.Write();
 f_res.Close();
//TH1D hist_test("test", "test", 100, -1, 1);
//for(int i=0; i<








}
