{


TFile my_file("Out_scan2020_tr_ph_fc_e935_v8.root", "READ");
TLorentzVector* Pip;
TLorentzVector* Pim;
TLorentzVector* Pph1;
TLorentzVector* Pph2;
double m2g, m2g0;
double m2g_miss, m2g0_miss;
int chi2_err, mode1, mode2;
double Chi25C, Chi2Miss;
int Err_Chi25C;
int correspond[20] = {0};



TTree *mytree = (TTree*)my_file.Get("tr_ph;12");



TTree* tree_err_5C = (TTree*)my_file.Get("tree_err_5C");
int Nentries = mytree->GetEntriesFast();
printf("N = %d\n", Nentries);

mytree->SetBranchAddress("Lpip", &Pip);
mytree->SetBranchAddress("Lph1", &Pph1);
mytree->SetBranchAddress("m2g", &m2g);
mytree->SetBranchAddress("m2g0", &m2g0);
mytree->SetBranchAddress("m2g_miss", &m2g_miss);
mytree->SetBranchAddress("m2g0_miss", &m2g0_miss);
mytree->SetBranchAddress("Chi25C", &Chi25C);
mytree->SetBranchAddress("Xi2misspi0", &Chi2Miss);
mytree->SetBranchAddress("mode1", &mode1);
mytree->SetBranchAddress("mode2", &mode2);
mytree->SetBranchAddress("correspond", correspond);

tree_err_5C->SetBranchAddress("Err_Chi25C", &Err_Chi25C);
tree_err_6C->SetBranchAddress("Err_Chi26C", &Err_Chi25C);

int counter_over = 0;
int counter_all = 0;

TH1D hist_pPl("pPl", "pPl", 100, 0, 2);
TH1D hist_pPh("pPh", "pPh", 100, 0, 2);
TH1D hist_m2g("M2g", "M2g", 200, 0.03, 0.3);
TH1D hist_m2g0("M2g0", "M2g0", 200, 0.03, 0.3);
TH1D hist_m2g_miss("M2g_miss", "M2g_miss", 200, 0.03, 0.3);
TH1D hist_m2g0_miss("M2g0_miss", "M2g0_miss", 200, 0.03, 0.3);
TH1D hist_Chi25C("Chi25C", "Chi25C", 2000, -1000, 1000);
TH1D hist_Chi26C("Chi26C", "Chi26C", 2000, -1000, 1000);
TH1D hist_Chi2Miss("Chi2Miss", "Chi2Miss", 2000, -1000, 1000);
TH1D hist_err_5C("Err_5C", "Err_5C", 3, 0, 3);
TH1D hist_err_6C("Err_5C", "Err_5C", 3, 0, 3);

for(int i=0; i<Nentries; i++){
	mytree->GetEntry(i);
	
	if (Chi25C < 50){
		hist_pPl.Fill(Pip->P()/1000.);
		hist_pPh.Fill(Pph1->P()/1000.);
		hist_m2g.Fill(m2g/1000.);
		hist_m2g0.Fill(m2g0/1000.);
	}
	if (mode1 == 1){
		hist_Chi25C.Fill(Chi25C);
		hist_Chi26C.Fill(Chi26C);
	}
	if (mode2 == 1){
		hist_Chi2Miss.Fill(Chi2Miss);	
	}

	if (mode2 == 1 && Chi2Miss < 10){
		hist_m2g_miss.Fill(m2g_miss/1000.);
		hist_m2g0_miss.Fill(m2g0_miss/1000.);
	}
	if (Pip->E() > 2000 and correspond[0] != -1){
		counter_over++;	
	}

}

double ineff = 0;

for(int i=0; i<tree_err_5C->GetEntries(); i++){
	tree_err_5C->GetEntry(i);
	if (Err_Chi25C>0)ineff++;
	hist_err_5C.Fill(Err_Chi25C);
}
ineff /=  tree_err_5C->GetEntries();
std::cout << "Inefficiency = " << ineff << std::endl;


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
hist_Chi2Miss.SetLineWidth(3);
hist_Chi2Miss.Draw();

TCanvas c6("2", "2", 1000, 1000);
hist_m2g_miss.SetLineColor(kRed);
hist_m2g_miss.GetXaxis()->SetRangeUser(0.11, 0.16);
hist_m2g_miss.SetLineWidth(3);
hist_m2g0_miss.SetLineWidth(3);
hist_m2g_miss.SetTitle("0 < Chi2Miss < 10");
hist_m2g_miss.Draw();
hist_m2g0_miss.Draw("same");
gPad->SetGrid();
//gPad->SetLogy();

TLegend legend2(0, 1, 0.1, 0.9);
legend2.AddEntry("M2g_miss", "M2g_miss", "lep");
legend2.AddEntry("M2g0_miss", "M2g0_miss", "lep");
legend2.Draw();

TCanvas c7("c7", "", 1000, 1000);
hist_err_5C.SetLineWidth(3);
hist_err_5C.Draw();

 TFile f_res("Results.root", "RECREATE");
 c3.Write();
 f_res.Close();
//TH1D hist_test("test", "test", 100, -1, 1);
//for(int i=0; i<








}
