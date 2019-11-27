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
	double irrad_time = par[3]; // beam irradiation time (from start) t_0
	double time = x[0]; // seconds

	if (time < irrad_time){
		return flux * lifetime  +  ( n_init - flux*lifetime) * TMath::Exp(-time/lifetime);
	}else{
		return (flux * lifetime  +  ( n_init - flux*lifetime) * TMath::Exp(-irrad_time/lifetime)) * TMath::Exp(-(time-irrad_time)/lifetime);
	}
}

double daughter_rn(double *x,double *par){

	double lifetime = par[0]; // seconds
	double n_init_fr = par[1]; // N_Fr(t=0)
	double n_init_rn = par[2]; // N_Rn(t=0)
	double frFlux = par[3]; // particles per second as incoming beam
	double frLifetime = par[4]; // seconds
	double frBR = par[5]; // alpha-branching ratio
	double irrad_time = par[6]; // beam irradiation time
	double time = x[0]; // seconds

	double DeltaTauFr = 1.0/((1.0/lifetime)-(1.0/frLifetime));

	double constant_term = (1.-frBR)*frFlux*lifetime;
	double frdecay_coef = (1.-frBR)*((n_init_fr/frLifetime)-frFlux)*DeltaTauFr;
	double rndecay_coef = n_init_rn - constant_term - frdecay_coef;

	if (time < irrad_time){
		return constant_term + frdecay_coef*TMath::Exp(-time/frLifetime) + rndecay_coef*TMath::Exp(-time/lifetime);
	}else{
		double n_fr_at_t0 = frFlux*frLifetime+(n_init_fr-frFlux*frLifetime)*TMath::Exp(-irrad_time/frLifetime);
		double n_rn_at_t0 = constant_term + frdecay_coef*TMath::Exp(-irrad_time/frLifetime) + rndecay_coef*TMath::Exp(-irrad_time/lifetime);
		double fr_coef = (1.-frBR)*(DeltaTauFr/frLifetime)*n_fr_at_t0;
		double rn_coef = n_rn_at_t0 - fr_coef;
		return fr_coef*TMath::Exp(-(time-irrad_time)/frLifetime) + rn_coef*TMath::Exp(-(time-irrad_time)/lifetime);
	}
}

double daughter_at(double *x,double *par){

	double lifetime = par[0]; // seconds
	double n_init_fr = par[1]; // N_Fr(t=0)
	double n_init_at = par[2]; // N_At(t=0)
	double frFlux = par[3]; // particles per second as incoming beam
	double frLifetime = par[4]; // seconds
	double frBR = par[5]; // alpha-branching ratio
	double irrad_time = par[6]; // beam irradiation time
	double time = x[0]; // seconds

	double DeltaTauFr = 1.0/((1.0/lifetime)-(1.0/frLifetime));

	double constant_term = frBR*frFlux*lifetime;
	double frdecay_coef = frBR*((n_init_fr/frLifetime)-frFlux)*DeltaTauFr;
	double atdecay_coef = n_init_at - constant_term - frdecay_coef;

	if (time < irrad_time){
		return constant_term + frdecay_coef*TMath::Exp(-time/frLifetime) + atdecay_coef*TMath::Exp(-time/lifetime);
	}else{
		double n_fr_at_t0 = frFlux*frLifetime+(n_init_fr-frFlux*frLifetime)*TMath::Exp(-irrad_time/frLifetime);
		double n_at_at_t0 = constant_term + frdecay_coef*TMath::Exp(-irrad_time/frLifetime) + atdecay_coef*TMath::Exp(-irrad_time/lifetime);
		double fr_coef = frBR*(DeltaTauFr/frLifetime)*n_fr_at_t0;
		double at_coef = n_at_at_t0 - fr_coef;
		return fr_coef*TMath::Exp(-(time-irrad_time)/frLifetime) + at_coef*TMath::Exp(-(time-irrad_time)/lifetime);
	}
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
	double T = 1170.0; // K
	double irradiation_time = 15.*minutes;
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
	TF1 *N_208fr = new TF1("{}^{208}Fr",parent,0.,timelimit,4);
	N_208fr->SetParameters(flux_208fr,t_208fr,0.0,irradiation_time);
	TGraph *g_N208fr = new TGraph(N_208fr);
	g_N208fr->SetMarkerColor(3);
	g_N208fr->SetLineColor(3);

	double flux_209fr = primaryFlux*R_209fr*diff_fac*si_eff*trans_eff;
	double br_209fr = 0.89;
	double t_209fr = 50.0*seconds/TMath::Log(2.);
	double e_209fr = 6.646; // MeV
	TF1 *N_209fr = new TF1("{}^{209}Fr",parent,0.,timelimit,4);
	N_209fr->SetParameters(flux_209fr,t_209fr,0.0,irradiation_time);
	TGraph *g_N209fr = new TGraph(N_209fr);
	g_N209fr->SetMarkerColor(4);
	g_N209fr->SetLineColor(4);

	double flux_210fr = primaryFlux*R_210fr*diff_fac*si_eff*trans_eff;
	double br_210fr = 0.71;
	double t_210fr = 3.18*minutes/TMath::Log(2.);
	double e_210fr = 6.545; // MeV
	TF1 *N_210fr = new TF1("{}^{210}Fr",parent,0.,timelimit,4);
	N_210fr->SetParameters(flux_210fr,t_210fr,0.0,irradiation_time);
	TGraph *g_N210fr = new TGraph(N_210fr);
	g_N210fr->SetMarkerColor(2);
	g_N210fr->SetLineColor(2);

	double flux_211fr = primaryFlux*R_211fr*diff_fac*si_eff*trans_eff;
	double br_211fr = 0.80;
	double t_211fr = 3.10*minutes/TMath::Log(2.);
	double e_211fr = 6.537; // MeV
	TF1 *N_211fr = new TF1("{}^{211}Fr",parent,0.,timelimit,4);
	N_211fr->SetParameters(flux_211fr,t_211fr,0.0,irradiation_time);
	TGraph *g_N211fr = new TGraph(N_211fr);
	g_N211fr->SetMarkerColor(5);
	g_N211fr->SetLineColor(5);

	double br_208rn = 0.62;
	double t_208rn = 24.35*minutes/TMath::Log(2.);
	double e_208rn = 6.1401; // MeV
	TF1 *N_208rn = new TF1("{}^{208}Rn",daughter_rn,0.,timelimit,8);
	N_208rn->SetParameters(t_208rn,0.0,0.0,flux_208fr,t_208fr,br_208fr,irradiation_time);
	TGraph *g_N208rn = new TGraph(N_208rn);
	g_N208rn->SetMarkerColor(3);
	g_N208rn->SetLineColor(3);
	g_N208rn->SetLineStyle(6);

	double br_209rn = 0.17;
	double t_209rn = 28.5*minutes/TMath::Log(2.);
	double e_209rn = 6.039; // MeV
	TF1 *N_209rn = new TF1("{}^{209}Rn",daughter_rn,0.,timelimit,8);
	N_209rn->SetParameters(t_209rn,0.0,0.0,flux_209fr,t_209fr,br_209fr,irradiation_time);
	TGraph *g_N209rn = new TGraph(N_209rn);
	g_N209rn->SetMarkerColor(4);
	g_N209rn->SetLineColor(4);
	g_N209rn->SetLineStyle(6);

	double br_210rn = 0.96;
	double t_210rn = 2.4*hours/TMath::Log(2.);
	double e_210rn = 6.041; // MeV
	TF1 *N_210rn = new TF1("{}^{210}Rn",daughter_rn,0.,timelimit,8);
	N_210rn->SetParameters(t_210rn,0.0,0.0,flux_210fr,t_210fr,br_210fr,irradiation_time);
	TGraph *g_N210rn = new TGraph(N_210rn);
	g_N210rn->SetMarkerColor(2);
	g_N210rn->SetLineColor(2);
	g_N210rn->SetLineStyle(6);

	double br_211rn = 0.27;
	double t_211rn = 14.6*hours/TMath::Log(2.);
	double e_211rn = 5.7839; // MeV
	TF1 *N_211rn = new TF1("{}^{211}Rn",daughter_rn,0.,timelimit,8);
	N_211rn->SetParameters(t_211rn,0.0,0.0,flux_211fr,t_211fr,br_211fr,irradiation_time);
	TGraph *g_N211rn = new TGraph(N_211rn);
	g_N211rn->SetMarkerColor(5);
	g_N211rn->SetLineColor(5);
	g_N211rn->SetLineStyle(6);

	double br_204at = 0.04;
	double t_204at = 9.2*minutes/TMath::Log(2.);
	double e_204at = 5.9503; // MeV
	TF1 *N_204at = new TF1("{}^{204}At",daughter_at,0.,timelimit,8);
	N_204at->SetParameters(t_204at,0.0,0.0,flux_208fr,t_208fr,br_208fr,irradiation_time);
	TGraph *g_N204at = new TGraph(N_204at);
	g_N204at->SetMarkerColor(3);
	g_N204at->SetLineColor(3);
	g_N204at->SetLineStyle(2);

	double br_205at = 0.10;
	double t_205at = 26.2*minutes/TMath::Log(2.);
	double e_205at = 5.902; // MeV
	TF1 *N_205at = new TF1("{}^{205}At",daughter_at,0.,timelimit,8);
	N_205at->SetParameters(t_205at,0.0,0.0,flux_209fr,t_209fr,br_209fr,irradiation_time);
	TGraph *g_N205at = new TGraph(N_205at);
	g_N205at->SetMarkerColor(4);
	g_N205at->SetLineColor(4);
	g_N205at->SetLineStyle(2);

	double br_206at = 0.0;
	double t_206at = 30.6*minutes/TMath::Log(2.);
	double e_206at = 0.0; // MeV
	TF1 *N_206at = new TF1("{}^{206}At",daughter_at,0.,timelimit,8);
	N_206at->SetParameters(t_206at,0.0,0.0,flux_210fr,t_210fr,br_210fr,irradiation_time);
	TGraph *g_N206at = new TGraph(N_206at);
	g_N206at->SetMarkerColor(2);
	g_N206at->SetLineColor(2);
	g_N206at->SetLineStyle(2);

	double br_207at = 0.09;
	double t_207at = 1.80*hours/TMath::Log(2.);
	double e_207at = 5.758; // MeV
	TF1 *N_207at = new TF1("{}^{207}At",daughter_at,0.,timelimit,8);
	N_207at->SetParameters(t_207at,0.0,0.0,flux_211fr,t_211fr,br_211fr,irradiation_time);
	TGraph *g_N207at = new TGraph(N_207at);
	g_N207at->SetMarkerColor(5);
	g_N207at->SetLineColor(5);
	g_N207at->SetLineStyle(2);


	TMultiGraph *N = new TMultiGraph("N","Ions on the MCP Surface / in the MCP; Time Elapsed (s); Ions");
	N->Add(g_N208fr);
	N->Add(g_N209fr);
	N->Add(g_N210fr);
	N->Add(g_N211fr);
	N->Add(g_N208rn);
	N->Add(g_N209rn);
	N->Add(g_N210rn);
	N->Add(g_N211rn);
	N->Add(g_N204at);
	N->Add(g_N205at);
	N->Add(g_N206at);
	N->Add(g_N207at);
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

	TGraph *g_alpha208rn = new TGraph(g_N208rn->GetN());
	for (int i=0; i<g_N208rn->GetN(); ++i) {
		g_alpha208rn->GetX()[i] = g_N208rn->GetX()[i];
		g_alpha208rn->GetY()[i] = g_N208rn->GetY()[i]*att_eff*det_eff*dir_prob*br_208rn/t_208rn;
	}
	g_alpha208rn->SetTitle(Form("{}^{208}Rn: %g MeV",e_208rn));
	g_alpha208rn->SetMarkerColor(3);
	g_alpha208rn->SetLineColor(3);
	g_alpha208rn->SetLineWidth(2);
	g_alpha208rn->SetLineStyle(6);

	TGraph *g_alpha209rn = new TGraph(g_N209rn->GetN());
	for (int i=0; i<g_N209rn->GetN(); ++i) {
		g_alpha209rn->GetX()[i] = g_N209rn->GetX()[i];
		g_alpha209rn->GetY()[i] = g_N209rn->GetY()[i]*att_eff*det_eff*dir_prob*br_209rn/t_209rn;
	}
	g_alpha209rn->SetTitle(Form("{}^{209}Rn: %g MeV",e_209rn));
	g_alpha209rn->SetMarkerColor(4);
	g_alpha209rn->SetLineColor(4);
	g_alpha209rn->SetLineWidth(2);
	g_alpha209rn->SetLineStyle(6);

	TGraph *g_alpha210rn = new TGraph(g_N210rn->GetN());
	for (int i=0; i<g_N210rn->GetN(); ++i) {
		g_alpha210rn->GetX()[i] = g_N210rn->GetX()[i];
		g_alpha210rn->GetY()[i] = g_N210rn->GetY()[i]*att_eff*det_eff*dir_prob*br_210rn/t_210rn;
	}
	g_alpha210rn->SetTitle(Form("{}^{210}Rn: %g MeV",e_210rn));
	g_alpha210rn->SetMarkerColor(2);
	g_alpha210rn->SetLineColor(2);
	g_alpha210rn->SetLineWidth(2);
	g_alpha210rn->SetLineStyle(6);

	TGraph *g_alpha211rn = new TGraph(g_N211rn->GetN());
	for (int i=0; i<g_N211rn->GetN(); ++i) {
		g_alpha211rn->GetX()[i] = g_N211rn->GetX()[i];
		g_alpha211rn->GetY()[i] = g_N211rn->GetY()[i]*att_eff*det_eff*dir_prob*br_211rn/t_211rn;
	}
	g_alpha211rn->SetTitle(Form("{}^{211}Rn: %g MeV",e_211rn));
	g_alpha211rn->SetMarkerColor(5);
	g_alpha211rn->SetLineColor(5);
	g_alpha211rn->SetLineWidth(2);
	g_alpha211rn->SetLineStyle(6);

	TGraph *g_alpha204at = new TGraph(g_N204at->GetN());
	for (int i=0; i<g_N204at->GetN(); ++i) {
		g_alpha204at->GetX()[i] = g_N204at->GetX()[i];
		g_alpha204at->GetY()[i] = g_N204at->GetY()[i]*att_eff*det_eff*dir_prob*br_204at/t_204at;
	}
	g_alpha204at->SetTitle(Form("{}^{204}At: %g MeV",e_204at));
	g_alpha204at->SetMarkerColor(3);
	g_alpha204at->SetLineColor(3);
	g_alpha204at->SetLineWidth(2);
	g_alpha204at->SetLineStyle(2);

	TGraph *g_alpha205at = new TGraph(g_N205at->GetN());
	for (int i=0; i<g_N205at->GetN(); ++i) {
		g_alpha205at->GetX()[i] = g_N205at->GetX()[i];
		g_alpha205at->GetY()[i] = g_N205at->GetY()[i]*att_eff*det_eff*dir_prob*br_205at/t_205at;
	}
	g_alpha205at->SetTitle(Form("{}^{205}At: %g MeV",e_205at));
	g_alpha205at->SetMarkerColor(4);
	g_alpha205at->SetLineColor(4);
	g_alpha205at->SetLineWidth(2);
	g_alpha205at->SetLineStyle(2);

	TGraph *g_alpha206at = new TGraph(g_N206at->GetN());
	for (int i=0; i<g_N206at->GetN(); ++i) {
		g_alpha206at->GetX()[i] = g_N206at->GetX()[i];
		g_alpha206at->GetY()[i] = g_N206at->GetY()[i]*att_eff*det_eff*dir_prob*br_206at/t_206at;
	}
	g_alpha206at->SetTitle(Form("{}^{206}At: %g MeV",e_206at));
	g_alpha206at->SetMarkerColor(2);
	g_alpha206at->SetLineColor(2);
	g_alpha206at->SetLineWidth(2);
	g_alpha206at->SetLineStyle(2);

	TGraph *g_alpha207at = new TGraph(g_N207at->GetN());
	for (int i=0; i<g_N207at->GetN(); ++i) {
		g_alpha207at->GetX()[i] = g_N207at->GetX()[i];
		g_alpha207at->GetY()[i] = g_N207at->GetY()[i]*att_eff*det_eff*dir_prob*br_207at/t_207at;
	}
	g_alpha207at->SetTitle(Form("{}^{207}At: %g MeV",e_207at));
	g_alpha207at->SetMarkerColor(5);
	g_alpha207at->SetLineColor(5);
	g_alpha207at->SetLineWidth(2);
	g_alpha207at->SetLineStyle(2);


	TMultiGraph *alpha = new TMultiGraph("alpha","Flux of #alpha Particles Detected at the SSD; Time Elapsed (s); #alpha Particles (/s)");
	alpha->Add(g_alpha208fr);
	alpha->Add(g_alpha209fr);
	alpha->Add(g_alpha210fr);
	alpha->Add(g_alpha211fr);
	alpha->Add(g_alpha208rn);
	alpha->Add(g_alpha209rn);
	alpha->Add(g_alpha210rn);
	alpha->Add(g_alpha211rn);
	alpha->Add(g_alpha204at);
	alpha->Add(g_alpha205at);
	alpha->Add(g_alpha206at);
	alpha->Add(g_alpha207at);
	alpha->Draw("AL");



	c1->cd(4);

	double flux_net = flux_208fr+flux_209fr+flux_210fr+flux_211fr;
	double alpharate_208fr = g_alpha208fr->Eval(irradiation_time);
	double alpharate_209fr = g_alpha209fr->Eval(irradiation_time);
	double alpharate_210fr = g_alpha210fr->Eval(irradiation_time);
	double alpharate_211fr = g_alpha211fr->Eval(irradiation_time);
	double alpharate_net = alpharate_208fr+alpharate_209fr+alpharate_210fr+alpharate_211fr;

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
	l.DrawLatex(0.1,0.9,Form("After %g seconds of beam irradiation (%g enA):",irradiation_time,beam_current));
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
