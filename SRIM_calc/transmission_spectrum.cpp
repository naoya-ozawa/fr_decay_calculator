#include <bits/stdc++.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLatex.h>
#include <TRint.h>
using namespace std;

int main(int argc, char** argv){
	TRint rootapp("app",&argc,argv);
	TCanvas *c1 = new TCanvas();
//	c1->Divide(2,1);

	fstream file;
	string word, filename;

//	filename = "./TRANSMIT_Be_Window_100.txt";
	filename = "./TRANSMIT_Be_Window_10000.txt";

	file.open(filename.c_str());

	TH1D *spectrum = new TH1D("spectrum","Energy Spectrum of Transmitted Ions",1500,90.,105.);
	spectrum->GetXaxis()->SetTitle("Energy (MeV)");
	spectrum->GetYaxis()->SetTitle("Ions / 0.01 MeV");

	int dataflag = 0;
	int datacol = 0;

	while (file >> word){
//		cout << word << endl;
		if (word == "T") dataflag = 1;
		if (dataflag == 1){
			++datacol;
			if (datacol%10 == 1){
				string head = word;
			}else if (datacol%10 == 2){
				int Nion = stoi(word);
			}else if (datacol%10 == 3){
				int A = stoi(word);
			}else if (datacol%10 == 4){
				double energy = stod(word)*TMath::Power(10.,-6);
				spectrum->Fill(energy);
			}else if (datacol%10 == 5){
				double depth = stod(word);
			}else if (datacol%10 == 6){
				double LposY = stod(word);
			}else if (datacol%10 == 7){
				double LposZ = stod(word);
			}else if (datacol%10 == 8){
				double cosx = stod(word);
			}else if (datacol%10 == 9){
				double cosy = stod(word);
			}else{
				double cosz = stod(word);
			}
		}
	}

	cout << datacol/10 << "ions analyzed" << endl;

//	c1->cd(1);

	spectrum->Draw();

	TF1 *fitspect = new TF1("fitspect","gaus(0)",90.,105.);
	fitspect->SetParameters(10000.,100.,1.);
	spectrum->Fit("fitspect","M");

	TLatex *l_spect = new TLatex();
	l_spect->SetTextAlign(12);
	l_spect->SetTextSize(0.05);
	l_spect->DrawLatex(91.,250,Form("<KE> = %g #pm %g MeV",fitspect->GetParameter(1),fitspect->GetParError(1)));
	l_spect->DrawLatex(91.,230,Form("#Delta E = %g #pm %g MeV",fitspect->GetParameter(2),fitspect->GetParError(2)));

//	c1->cd(2);

	// The emittance can be calculated using cosx = v_x/v_tot, etc.
	//
	//
	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;
}
