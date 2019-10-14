#include <iostream>
#include <TCanvas.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TMath.h>
#include <TRint.h>
using namespace std;

double parent(double *x,double *par){

	double flux = par[0]; // particles per second as incoming beam
	double lifetime = par[1]; // seconds
	double time = x[0]; // seconds

	return flux * lifetime * ( 1.0 - TMath::Exp(-time/lifetime));
}

double daughter(double *x,double *par){

	double flux = par[0]; // particles per second as incoming beam
	double lifetime = par[1]; // seconds
	double alphaParentFlux = par[2]; // particles per second as incoming beam
	double alphaParentLifetime = par[3]; // seconds
	double alphaParentBR = par[4]; // alpha-branching ratio
	double betaParentFlux = par[5]; // particles per second as incoming beam
	double betaParentLifetime = par[6]; // seconds
	double betaParentBR = par[7]; // alpha-branching ratio
	double time = x[0]; // seconds

	double coef_daughter = flux * lifetime;
	double coef_alphaParent = (alphaParentBR*alphaParentFlux*lifetime*alphaParentLifetime) / (alphaParentLifetime - lifetime);
	double coef_betaParent = ((1.0-betaParentBR)*betaParentFlux*lifetime*betaParentLifetime) / (betaParentLifetime - lifetime);
	
	double exp_daughter = 1.0 - TMath::Exp(-time/lifetime);
	double exp_alphaParent = TMath::Exp(-time/alphaParentLifetime) - TMath::Exp(-time/lifetime);
	double exp_betaParent = TMath::Exp(-time/betaParentLifetime) - TMath::Exp(-time/lifetime);

	return coef_daughter*exp_daughter + coef_alphaParent*exp_alphaParent + coef_betaParent*exp_betaParent;
}

double ni_ratio(double T, double E_wf, double E_ip, double sw){
	double k = 8.61733 * TMath::Power(10.,-5); // eV/K
	return sw * TMath::Exp( (E_wf - E_ip) / (k*T) );
}

double ionization(double T, double E_wf, double E_ip, double sw){
	return ni_ratio(T,E_wf,E_ip,sw) / (1.0 + ni_ratio(T,E_wf,E_ip,sw));
}

int main (int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	TCanvas *c1 = new TCanvas();

	double seconds = 1.0;
	double minutes = 60.0;
	double hours = 60.0 * minutes;
	double days = 24.0 * hours;
	double years = 365.0 * days;

	double secondaryFlux = 1.0*TMath::Power(10.,7); // Particles per second
	double T = 800.0; // K
	double timelimit = 10.0*minutes;

	double E_wf_Mo = 4.6; // eV
	double E_wf_Au = 5.1; // eV

	double E_ip_Fr = 4.07; // eV

	double R_210fr = 0.632;
	double br_210fr = 0.71;
	double t_210fr = 3.18*minutes;
	double e_210fr = 6.5; // MeV

	TF1 *N_210fr = new TF1("N_210Fr",parent,0.,timelimit,2);
	N_210fr->SetParameters(secondaryFlux*R_210fr*ionization(T,E_wf_Au,E_ip_Fr,0.5),t_210fr);
	TGraph *g_210fr = new TGraph(N_210fr);
	g_210fr->SetMarkerColor(kRed);
	g_210fr->SetLineColor(kRed);

	double R_211fr = 0.034;
	double br_211fr = 0.80;
	double t_211fr = 3.10*minutes;
	double e_211fr = 6.7; // MeV

	TF1 *N_211fr = new TF1("N_211Fr",parent,0.,timelimit,2);
	N_211fr->SetParameters(secondaryFlux*R_211fr*ionization(T,E_wf_Au,E_ip_Fr,0.5),t_211fr);
	TGraph *g_211fr = new TGraph(N_211fr);
	g_211fr->SetMarkerColor(kBlue);
	g_211fr->SetLineColor(kBlue);

	double R_209fr = 0.262;
	double br_209fr = 0.89;
	double t_209fr = 50.0*seconds;
	double e_209fr = 6.9; // MeV

	TF1 *N_209fr = new TF1("N_209Fr",parent,0.,timelimit,2);
	N_209fr->SetParameters(secondaryFlux*R_209fr*ionization(T,E_wf_Au,E_ip_Fr,0.5),t_209fr);
	TGraph *g_209fr = new TGraph(N_209fr);
	g_209fr->SetMarkerColor(kOrange);
	g_209fr->SetLineColor(kOrange);

	
	TMultiGraph *N = new TMultiGraph("N","Number of Ions at MCP Surface; Time Elapsed (s); Ions");
	N->Add(g_210fr);
	N->Add(g_211fr);
	N->Add(g_209fr);

	N->Draw("AL");

	c1->BuildLegend();
	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;

//	double R_208rn = 0.001;
//	double br_208rn = ;

//	double R_211rn = 0.003;
//	double br_211rn = 0.274;
//	double t_211rn = 14.6*hours;
//	double e_211rn = 6.0; // MeV

//	double R_208at = 0.015;
//	double br_208at;

//	double R_207at = 0.049;
//	double br_207at;
	
//	double R_206at = 0.006;


//	double br_211rn = 0.274;
//	double brerr_211rn = 0.017;
//	double t_211rn = 14.6*hours;
//	double terr_211rn = 0.2*hours;

//	double br_211at = 0.4180;
//	double brerr_211at = 0.0008;
//	double t_211at = 7.214*hours;
//	double terr_211at = 0.007*hours;

//	double br_211po = 1.0;
//	double brerr_211po = 0.0;
//	double t_211po = 0.516*seconds;
//	double terr_211po = 0.003*seconds;

//	double br_207at = 0.086;
//	double brerr_207at = 0.01;
//	double t_207at = 1.80*hours;
//	double terr_207at = 0.04*hours;

//	double br_207po = 0.00021;
//	double brerr_207po = 0.00002;
//	double t_207po = 5.80*hours;
//	double terr_207po = 0.02*hours;

//	double br_207bi = 0.0;
//	double brerr_207bi = 0.0;
//	double t_207bi = 31.55*years;
//	double terr_207bi = 0.05*years;

//	double br_207pb = 0.0;
//	double brerr_207pb = 0.0;
//	double t_207pb = 10000.*years;
//	double terr_207pb = 0.0;


}
