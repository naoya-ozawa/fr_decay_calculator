#include <iostream>
#include <TCanvas.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TMath.h>
#include <TLatex.h>
#include <TRint.h>
using namespace std;

double parent(double *x,double *par){

	double flux = par[0]; // particles per second as incoming beam
	double lifetime = par[1]; // seconds
	double n_init = par[2]; // initial number N_Fr(t=0)
//	double irrad_time = par[3]; // beam irradiation time (from start) t_0
	double time = x[0]; // seconds

//	if (time < irrad_time){
		return flux * lifetime  +  ( n_init - flux*lifetime) * TMath::Exp(-time/lifetime);
//	}else{
//		return flux*lifetime*(1.0-TMath::Exp(-irrad_time/lifetime)) * TMath::Exp(-(time-irrad_time)/lifetime);
//	}
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
	double irrad_time = par[8]; // beam irradiation time
	double time = x[0]; // seconds

	double coef_daughter = flux * lifetime;
	double dt_alphaParent = 1./((1./lifetime)-(1./alphaParentLifetime));
	double coef_alphaParent = alphaParentFlux*alphaParentBR*dt_alphaParent;
	double dt_betaParent = 1./((1./lifetime)-(1./betaParentLifetime));
	double coef_betaParent = betaParentFlux*(1.0-betaParentBR)*dt_betaParent;

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
	c1->Divide(2,2);

	double seconds = 1.0;
	double minutes = 60.0;
	double hours = 60.0 * minutes;
	double days = 24.0 * hours;
	double years = 365.0 * days;

	double beam_current = 4000.; // enA of 18-O-6+
	double primaryFlux = beam_current*TMath::Power(10.,-9)/6./(1.6*TMath::Power(10.,-19)); // Particles per second: for all the Fr isotopes total
	double T = 1100.0; // K
	double timelimit = 30.*minutes;

	double E_wf_Mo = 4.6; // eV
	double E_wf_Au = 5.1; // eV

	double E_ip_Fr = 4.07; // eV
	double E_ip_Rn = 10.74850; // eV
	double E_ip_At = 9.31751; // eV

	double diff_fac = 0.5; // probability that the Fr diffuses to the surface of the target
	//double escape_eff = ?; // probability that the Fr escapes from the target within its lifetime --> dependent on temperature
	double si_eff = ionization(T,E_wf_Au,E_ip_Fr,0.5); // Surface ionizatin efficiency
	double trans_eff = 0.3; // transportation efficiency from the target to the MCP surface
	double att_eff = 0.37; // Open-area-ratio of the MCP-IN surface
	double det_eff = 0.005; // detection efficiency of the SSD is around 0.5%
	double dir_prob = 0.5; // direction probability of the alpha emission = 1/2

	c1->cd(1);
	TLatex *ls_ratio = new TLatex();
	ls_ratio->SetTextAlign(12);
	ls_ratio->SetTextSize(0.05);
	ls_ratio->DrawLatex(0.1,0.9,Form("Surface ionization probabilities at %g #circC",T-273.));
	ls_ratio->DrawLatex(0.2,0.8,Form("Fr: %g %%",100.*ionization(T,E_wf_Au,E_ip_Fr,0.5)));
	ls_ratio->DrawLatex(0.2,0.7,Form("Rn: %g %%",100.*ionization(T,E_wf_Au,E_ip_Rn,1.0)));
	ls_ratio->DrawLatex(0.2,0.6,Form("At: %g %%",100.*ionization(T,E_wf_Au,E_ip_At,0.25)));
	ls_ratio->DrawLatex(0.1,0.5,"We assume that only Fr is extracted from the target.");


	c1->cd(2);
	double beam_energy = 6.94; // MeV/u --> 112 MeV injected to target
	double R_210fr = 4.4e-06; // normalized flux from stancari2006
	double R_209fr = 3.1e-06; // normalized flux for 114 MeV
	double R_211fr = 1.2e-06; // normalized flux from stancari2006
	double R_208fr = 0.3e-06; // normalized flux from stancari2006


	double flux_208fr = primaryFlux*R_208fr*diff_fac*si_eff*trans_eff;
	double br_208fr = 0.89;
	double t_208fr = 59.1*seconds/TMath::Log(2.);
	double e_208fr = 6.641; // MeV
	TF1 *N_208fr = new TF1("{}^{208}Fr",parent,0.,timelimit,3);
	N_208fr->SetParameters(flux_208fr,t_208fr,0.0);
	TGraph *g_N208fr = new TGraph(N_208fr);
	g_N208fr->SetMarkerColor(3);
	g_N208fr->SetLineColor(3);

	double flux_209fr = primaryFlux*R_209fr*diff_fac*si_eff*trans_eff;
	double br_209fr = 0.89;
	double t_209fr = 50.0*seconds/TMath::Log(2.);
	double e_209fr = 6.646; // MeV
	TF1 *N_209fr = new TF1("{}^{209}Fr",parent,0.,timelimit,3);
	N_209fr->SetParameters(flux_209fr,t_209fr,0.0);
	TGraph *g_N209fr = new TGraph(N_209fr);
	g_N209fr->SetMarkerColor(4);
	g_N209fr->SetLineColor(4);

	double flux_210fr = primaryFlux*R_210fr*diff_fac*si_eff*trans_eff;
	double br_210fr = 0.71;
	double t_210fr = 3.18*minutes/TMath::Log(2.);
	double e_210fr = 6.545; // MeV
	TF1 *N_210fr = new TF1("{}^{210}Fr",parent,0.,timelimit,3);
	N_210fr->SetParameters(flux_210fr,t_210fr,0.0);
	TGraph *g_N210fr = new TGraph(N_210fr);
	g_N210fr->SetMarkerColor(2);
	g_N210fr->SetLineColor(2);

	double flux_211fr = primaryFlux*R_211fr*diff_fac*si_eff*trans_eff;
	double br_211fr = 0.80;
	double t_211fr = 3.10*minutes/TMath::Log(2.);
	double e_211fr = 6.537; // MeV
	TF1 *N_211fr = new TF1("{}^{211}Fr",parent,0.,timelimit,3);
	N_211fr->SetParameters(flux_211fr,t_211fr,0.0);
	TGraph *g_N211fr = new TGraph(N_211fr);
	g_N211fr->SetMarkerColor(5);
	g_N211fr->SetLineColor(5);


	TMultiGraph *N = new TMultiGraph("N","Ions on the MCP Surface / in the MCP; Time Elapsed (s); Ions");
	N->Add(g_N208fr);
	N->Add(g_N209fr);
	N->Add(g_N210fr);
	N->Add(g_N211fr);
	N->Draw("AL");

	c1->cd(3);

	TGraph *g_alpha208fr = new TGraph(g_N208fr->GetN());
	for (int i=0; i<g_N208fr->GetN(); ++i) {
		g_alpha208fr->GetX()[i] = g_N208fr->GetX()[i];
		g_alpha208fr->GetY()[i] = g_N208fr->GetY()[i]*att_eff*det_eff*dir_prob*br_208fr/t_208fr;
	}
	g_alpha208fr->SetTitle(Form("{}^{208}Fr: %g MeV",e_208fr));
	g_alpha208fr->SetMarkerColor(3);
	g_alpha208fr->SetLineColor(3);
	g_alpha208fr->SetLineWidth(2);

	TGraph *g_alpha209fr = new TGraph(g_N209fr->GetN());
	for (int i=0; i<g_N209fr->GetN(); ++i) {
		g_alpha209fr->GetX()[i] = g_N209fr->GetX()[i];
		g_alpha209fr->GetY()[i] = g_N209fr->GetY()[i]*att_eff*det_eff*dir_prob*br_209fr/t_209fr;
	}
	g_alpha209fr->SetTitle(Form("{}^{209}Fr: %g MeV",e_209fr));
	g_alpha209fr->SetMarkerColor(4);
	g_alpha209fr->SetLineColor(4);
	g_alpha209fr->SetLineWidth(2);

	TGraph *g_alpha210fr = new TGraph(g_N210fr->GetN());
	for (int i=0; i<g_N210fr->GetN(); ++i) {
		g_alpha210fr->GetX()[i] = g_N210fr->GetX()[i];
		g_alpha210fr->GetY()[i] = g_N210fr->GetY()[i]*att_eff*det_eff*dir_prob*br_210fr/t_210fr;
	}
	g_alpha210fr->SetTitle(Form("{}^{210}Fr: %g MeV",e_210fr));
	g_alpha210fr->SetMarkerColor(2);
	g_alpha210fr->SetLineColor(2);
	g_alpha210fr->SetLineWidth(2);

	TGraph *g_alpha211fr = new TGraph(g_N211fr->GetN());
	for (int i=0; i<g_N211fr->GetN(); ++i) {
		g_alpha211fr->GetX()[i] = g_N211fr->GetX()[i];
		g_alpha211fr->GetY()[i] = g_N211fr->GetY()[i]*att_eff*det_eff*dir_prob*br_211fr/t_211fr;
	}
	g_alpha211fr->SetTitle(Form("{}^{211}Fr: %g MeV",e_211fr));
	g_alpha211fr->SetMarkerColor(5);
	g_alpha211fr->SetLineColor(5);
	g_alpha211fr->SetLineWidth(2);

	TMultiGraph *alpha = new TMultiGraph("alpha","Flux of #alpha Particles Detected at the SSD; Time Elapsed (s); #alpha Particles (/s)");
	alpha->Add(g_alpha208fr);
	alpha->Add(g_alpha209fr);
	alpha->Add(g_alpha210fr);
	alpha->Add(g_alpha211fr);
	alpha->Draw("AL");



	c1->cd(4);

	double flux_net = flux_208fr+flux_209fr+flux_210fr+flux_211fr;
	double alpharate_208fr = g_alpha208fr->Eval(timelimit);
	double alpharate_209fr = g_alpha209fr->Eval(timelimit);
	double alpharate_210fr = g_alpha210fr->Eval(timelimit);
	double alpharate_211fr = g_alpha211fr->Eval(timelimit);
	double alpharate_net = alpharate_208fr+alpharate_209fr+alpharate_210fr+alpharate_211fr;

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
	l.DrawLatex(0.1,0.9,Form("After %g seconds of beam irradiation (%g enA):",timelimit,beam_current));
	l.DrawLatex(0.2,0.8,Form("{}^{208}Fr: %g/s production, #alpha detection %g/s (%g MeV)",flux_208fr,alpharate_208fr,e_208fr));
	l.DrawLatex(0.2,0.7,Form("{}^{209}Fr: %g/s production, #alpha detection %g/s (%g MeV)",flux_209fr,alpharate_209fr,e_209fr));
	l.DrawLatex(0.2,0.6,Form("{}^{210}Fr: %g/s production, #alpha detection %g/s (%g MeV)",flux_210fr,alpharate_210fr,e_210fr));
	l.DrawLatex(0.2,0.5,Form("{}^{211}Fr: %g/s production, #alpha detection %g/s (%g MeV)",flux_211fr,alpharate_211fr,e_211fr));
	l.DrawLatex(0.1,0.3,Form("TOTAL: %g/s production, #alpha detection %g/s",flux_net,alpharate_net));


	c1->cd(2)->BuildLegend();
	c1->cd(3)->BuildLegend();
	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;

}
