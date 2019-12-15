#include <iostream>
#include <TCanvas.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>
#include <TRint.h>
using namespace std;

// Database for the lifetimes of each species/isotope
double tau_AZ(int A, const char* Z){
	double seconds = 1.0;
	double minutes = 60.0 * seconds;
	double hours = 60.0 * minutes;
	double days = 24.0 * hours;
	double years = 365.0 * days;

	double lifetime = -9999.;
	double lifetime_err = -9999.;

	if (Z == "Fr"){
		if (A == 208){
			lifetime = 59.1 * seconds;
			lifetime_err = 0.3 * seconds;
		}else if (A == 209){
			lifetime = 50.5 * seconds;
			lifetime_err = 0.7 * seconds;
		}else if (A == 210){
			lifetime = 3.18 * minutes;
			lifetime_err = 0.06 * minutes;
		}else if (A == 211){
			lifetime = 3.10 * minutes;
			lifetime_err = 0.02 * minutes;
		}
	}else if (Z == "Rn"){
		if (A == 208){
			lifetime = 24.35 * minutes;
			lifetime_err = 0.14 * minutes;
		}else if (A == 209){
			lifetime = 28.8 * minutes;
			lifetime_err = 0.9 * minutes;
		}else if (A == 210){
			lifetime = 2.4 * hours;
			lifetime_err = 0.1 * hours;
		}else if (A == 211){
			lifetime = 14.6 * hours;
			lifetime_err = 0.2 * hours;
		}
	}else if (Z == "At"){
		if (A == 204){
			lifetime = 9.22 * minutes;
			lifetime_err = 0.13 * minutes;
		}else if (A == 205){
			lifetime = 26.9 * minutes;
			lifetime_err = 0.8 * minutes;
		}else if (A == 206){
			lifetime = 30.6 * minutes;
			lifetime_err = 0.8 * minutes;
		}else if (A == 207){
			lifetime = 1.80 * hours;
			lifetime_err = 0.04 * hours;
		}else if (A == 209){
			lifetime = 5.41 * hours;
			lifetime_err = 0.05 * hours;
		}else if (A == 211){
			lifetime = 7.214 * hours;
			lifetime_err = 0.007 * hours;
		}
	}else if (Z == "Po"){
		if (A == 206){
			lifetime = 8.8 * days;
			lifetime_err = 0.1 * days;
		}else if (A == 211){
			lifetime = 0.516 * seconds;
			lifetime_err = 0.003 * seconds;
		}
	}

	// Half-life --> lifetime
	lifetime /= TMath::Log(2.0);
	lifetime_err /= TMath::Log(2.0);

	return lifetime;
}

// Database for the alpha-branching ratio of each species/isotope
double b_AZ(int A, const char* Z){
	double b = -9999.;
	double b_err = -9999.;

	if (Z == "Fr"){
		if (A == 208){
			b = 89.;
			b_err = 3.;
		}else if (A == 209){
			b = 89.;
			b_err = 3.;
		}else if (A == 210){
			b = 71.;
			b_err = 4.;
		}else if (A == 211){
			b = 87.;
			b_err = 3.;
		}
	}else if (Z == "Rn"){
		if (A == 208){
			b = 62.;
			b_err = 7.;
		}else if (A == 209){
			b = 17.;
			b_err = 2.;
		}else if (A == 210){
			b = 96.;
			b_err = 1.;
		}else if (A == 211){
			b = 27.4;
			b_err = 1.7;
		}
	}else if (Z == "At"){
		if (A == 204){
			b = 3.91;
			b_err = 0.16;
		}else if (A == 205){
			b = 10.;
			b_err = 2.;
		}else if (A == 206){
			b = 0.90;
			b_err = 0.08;
		}else if (A == 207){
			b = 8.6;
			b_err = 1.0;
		}else if (A == 209){
			b = 4.1;
			b_err = 0.5;
		}else if (A == 211){
			b = 41.80;
			b_err = 0.08;
		}
	}else if (Z == "Po"){
		if (A == 206){
			b = 5.45;
			b = 0.05;
		}else if (A == 211){
			b = 100.0;
			b_err = 0.0;
		}
	}

	// % --> ratio
	b /= 100.;
	b_err /= 100.;

	return b;
}

// Dataset for N (t = t_{beamON}) for each species/isotope
double N_on(int A, const char* Z){
	double number = -9999.;
	if (Z == "Fr"){
		if (A == 208){
			number = 0.0;
		}else if (A == 209){
			number = 0.0;
		}else if (A == 210){
			number = 0.0;
		}else if (A == 211){
			number = 0.0;
		}
	}else if (Z == "Rn"){
		if (A == 208){
			number = 0.0;
		}else if (A == 209){
			number = 0.0;
		}else if (A == 210){
			number = 0.0;
		}else if (A == 211){
			number = 0.0;
		}
	}else if (Z == "At"){
		if (A == 204){
			number = 0.0;
		}else if (A == 205){
			number = 0.0;
		}else if (A == 206){
			number = 0.0;
		}else if (A == 207){
			number = 0.0;
		}else if (A == 209){
			number = 0.0;
		}else if (A == 211){
			number = 0.0;
		}
	}else if (Z == "Po"){
		if (A == 206){
			number = 0.0;
		}else if (A == 211){
			number = 0.0;
		}
	}
	return number;
}

// Calculate the primary beam energy (MeV/u) after transmitting through the Be window
double Be_degraded(double beam_energy){
	// Fitted parameters from the dependence_plotter.cpp
	double C_E0 = -8.76; // MeV
	double C_E1 = 828.91; // MeV*MeV_incident
	double C_Be = 0.68; // MeV/um_Be
	double C_He = 0.0875; // MeV/mm_He

	// Be foil data
	double T_Be_05 = 9.75; // um; Haba Gr MEF#101147032 #16
	double T_Be_06 = 9.70; // um; Haba Gr MEF#101147032 #15
	double T_Be = T_Be_05+T_Be_06; // um; for 11/22 expt
	double T_He = 4.2 + 17.5 + 6.7; // mm; for 11/22 expt

	double E_loss = C_E0 + (C_E1/(beam_energy*18.)) + C_Be*T_Be + C_He*T_He; // MeV
//	cout << "Degraded energy: " << E*18.-E_loss << " MeV" << endl;
	return beam_energy - (E_loss/18.);
}


// Estimate the normalized production based on Stancari2006
double normprod(double incident_energy, int A){
	// normalized production data from Stancari2006
	double energy[10] = {82.,86.,90.,94.,98.,102.,106.,110.,114.,118.}; // MeV incident
	double fr208[10] = {0.0,0.0,0.0,1.8,2.1e+01,3.1e+03,4.7e+04,2.3e+05,5.9e+05,1.8e+06}; // Hz for j = 1e+12 particles/s
	double fr209[10] = {2.3e+01,7.6e+01,5.7e+03,1.2e+05,6.0e+05,1.5e+06,2.4e+06,2.9e+06,3.2e+06,3.2e+06};
	double fr210[10] = {2.2e+04,3.8e+05,1.5e+06,2.9e+06,3.9e+06,4.3e+06,4.4e+06,4.4e+06,4.4e+06,4.4e+06};
	double fr211[10] = {4.3e+05,8.9e+05,1.1e+06,1.2e+06,1.2e+06,1.2e+06,1.2e+06,1.2e+06,1.2e+06,1.2e+06};
	for (int i=0; i<10; ++i){
		fr208[i] /= TMath::Power(10.,12);
		fr209[i] /= TMath::Power(10.,12);
		fr210[i] /= TMath::Power(10.,12);
		fr211[i] /= TMath::Power(10.,12);
	}

	double E = incident_energy*18.; // MeV/u --> MeV, incident energy

	if (A == 208){
		TGraph *g208 = new TGraph(10,energy,fr208);
		TSpline3 *s208 = new TSpline3("Spline Fit for 208",g208);
		return s208->Eval(E);
	}else if (A==209){
		TGraph *g209 = new TGraph(10,energy,fr209);
		TSpline3 *s209 = new TSpline3("Spline Fit for 209",g209);
		return s209->Eval(E);
	}else if (A==210){
		TGraph *g210 = new TGraph(10,energy,fr210);
		TSpline3 *s210 = new TSpline3("Spline Fit for 210",g210);
		return s210->Eval(E);
	}else if (A==211){
		TGraph *g211 = new TGraph(10,energy,fr211);
		TSpline3 *s211 = new TSpline3("Spline Fit for 211",g211);
		return s211->Eval(E);
	}else{
		TCanvas *c_normprod = new TCanvas();
		TGraph *g208 = new TGraph(10,energy,fr208);
		TSpline3 *s208 = new TSpline3("Spline Fit for 208",g208);
		g208->SetTitle("{}^{208}Fr");
		g208->SetLineColor(3);
		g208->SetLineWidth(2);
		s208->SetLineColor(3);
		s208->SetLineWidth(2);
		TGraph *g209 = new TGraph(10,energy,fr209);
		TSpline3 *s209 = new TSpline3("Spline Fit for 209",g209);
		g209->SetTitle("{}^{209}Fr");
		g209->SetLineColor(4);
		g209->SetLineWidth(2);
		s209->SetLineColor(4);
		s209->SetLineWidth(2);
		TGraph *g210 = new TGraph(10,energy,fr210);
		TSpline3 *s210 = new TSpline3("Spline Fit for 210",g210);
		g210->SetTitle("{}^{210}Fr");
		g210->SetLineColor(2);
		g210->SetLineWidth(2);
		s210->SetLineColor(2);
		s210->SetLineWidth(2);
		TGraph *g211 = new TGraph(10,energy,fr211);
		TSpline3 *s211 = new TSpline3("Spline Fit for 211",g211);
		g211->SetTitle("{}^{211}Fr");
		g211->SetLineColor(5);
		g211->SetLineWidth(2);
		s211->SetLineColor(5);
		s211->SetLineWidth(2);
		TMultiGraph *mg_normprod = new TMultiGraph();
		mg_normprod->SetTitle("Normalized Production Based on [Stancari2006];Incident {}^{18}O beam energy (MeV);Normalized production P/j");
		mg_normprod->Add(g208);
		mg_normprod->Add(g209);
		mg_normprod->Add(g210);
		mg_normprod->Add(g211);
		mg_normprod->Draw("AP*");
		c_normprod->BuildLegend();
		s208->Draw("SAME");
		s209->Draw("SAME");
		s210->Draw("SAME");
		s211->Draw("SAME");
		TLine *l_208 = new TLine(81.,s208->Eval(E),119.,s208->Eval(E));
		l_208->SetLineColor(3);
		l_208->SetLineStyle(2);
		l_208->SetLineWidth(2);
		l_208->Draw();
		TLine *l_209 = new TLine(81.,s209->Eval(E),119.,s209->Eval(E));
		l_209->SetLineColor(4);
		l_209->SetLineStyle(2);
		l_209->SetLineWidth(2);
		l_209->Draw();
		TLine *l_210 = new TLine(81.,s210->Eval(E),119.,s210->Eval(E));
		l_210->SetLineColor(2);
		l_210->SetLineStyle(2);
		l_210->SetLineWidth(2);
		l_210->Draw();
		TLine *l_211 = new TLine(81.,s211->Eval(E),119.,s211->Eval(E));
		l_211->SetLineColor(5);
		l_211->SetLineStyle(2);
		l_211->SetLineWidth(2);
		l_211->Draw();
		TLine *l_eval = new TLine(E,1.*TMath::Power(10.,-7),E,4.9*TMath::Power(10.,-6));
		l_eval->SetLineColor(2);
		l_eval->SetLineStyle(2);
		l_eval->Draw();
		return -9999.;
	}
}

// Surface ionization probability
double ni_ratio(double T, double E_wf, double E_ip, double sw){
	double k = 8.61733 * TMath::Power(10.,-5); // eV/K
	return sw * TMath::Exp( (E_wf - E_ip) / (k*T) );
}

double ionization(double T, double E_wf, double E_ip, double sw){
	return ni_ratio(T,E_wf,E_ip,sw) / (1.0 + ni_ratio(T,E_wf,E_ip,sw));
}


// Diffusion coefficient calculation
double Au_selfdiffusion(double T){
  double Tm = 1338.0; // melting point (K)
  double D_0_1 = 0.025; // cm^2 s^-1
  double Q_1 = 1.7; // eV
  double D_0_2 = 0.83; // cm^2 s^-1
  double Q_2 = 2.2; // eV

  if (T >= Tm){ // Lu2006
    double D_Tm = 2.50*TMath::Power(10.0,-9); // m^2/s
//    double M = 197.0; // g/mol
//    double gamma_lv = 1.21; // J/m^2
//    double r = 0.82*TMath::Power(10.0,-10); // m
    double p = 0.13;
    double q = -0.12;
    return D_Tm*TMath::Power(10.,4) * TMath::Power(T/Tm,3./2.) / ( (1.+p-(p*T/Tm)) * TMath::Power(1.-q+(q*T/Tm),2./3.) );
  }else{ // Neumann1986
    double kB = 8.6*TMath::Power(10.,-5); // eV/K
    return D_0_1*TMath::Exp(-Q_1/(kB*T)) + D_0_2*TMath::Exp(-Q_2/(kB*T)); // cm^2 s^-1
  }
}

double D_m(double T, int A){
  double m_Au = 197.0; // u
  // Schoen1958
  return TMath::Sqrt(m_Au/double(A))*Au_selfdiffusion(T);
}

// production depth calculation
double range(double peak_energy,double incident_energy){
	// based on fit of the energy_dependence.cpp in the fr_production_calculator/SRIM_calculations
	// incident_energy: MeV/u
	return -2.6 + 0.34 * incident_energy*18. - 0.31 * peak_energy; // um
}

double depth_Fr(double incident_energy, int isotope){
	// energy: MeV/u (injection energy into Au target)
	if (isotope == 208){
		double peak = 115.; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 209){
		double peak = 100.; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 210){
		double peak = 90.; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 211){
		double peak = 82.; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else{
		return -9999.;
	}
}

// escape efficiency calculation
double y_Fr(double tau,double D,double d){
	double R_Au = 6.0 * 0.1; // mm --> cm
	double alpha = tau*D/(d*d);
	double x = R_Au/d;
	return 0.5 * (1. - TMath::Cos(TMath::ATan(x))) * TMath::Sqrt(alpha) * TMath::TanH(1./TMath::Sqrt(alpha));
}

// ε_{extraction} = ε_{transportation} ε_{desorption} ε_{ionization} ε_{escape}
double extraction_eff(double incident_energy, double temp, int A){

	double E_wf_Mo = 4.6; // eV
	double E_wf_Au = 5.1; // eV
	double E_ip_Fr = 4.07; // eV
	double E_ip_Rn = 10.74850; // eV
	double E_ip_At = 9.31751; // eV
	double E_ip_Po = 8.418070; // eV
	double sw_Fr = 0.5; // statistical weight factor for Fr
	double des_eff = 1.0; // desorption efficiency from the Au surface before decay (assumption)
	double trans_eff = 1.0; // transportation efficiency from the target to the MCP surface (assumption)

	double escape = y_Fr(tau_AZ(A,"Fr"),D_m(temp,A),depth_Fr(incident_energy,A))*ionization(temp,E_wf_Au,E_ip_Fr,sw_Fr)*des_eff*trans_eff;

	return escape;
}


// Function to calculate f_{^AFr} (pps) for each isotope A, given the O beam energy
// based on the Stancari2006 normalized flux
// Includes extraction efficiency calculation
double calc_flux(int A, double beam_energy, double beam_current, double temp){
	double E_deg = Be_degraded(beam_energy); // MeV/u, beam energy --> incidnet energy

	double conv_factor = TMath::Power(10.,-9)/(6.0*1.6*TMath::Power(10.,-19)); // enA --> pps
	double init_flux = beam_current * conv_factor; // enA --> pps
//	cout << "Primary beam flux: " << init_flux << endl;

	double prod_flux = -9999.;
	bool selection = (A > 207)&&(A < 212);
	if (selection){
		prod_flux = extraction_eff(E_deg,temp,A)*normprod(E_deg,A)*init_flux;
	}

	return prod_flux;
}


// N(t) for each species/isotopes
double FranciumA(double *x,double *par){
	int A = int(par[0]); // the species of interest (208/209/210/211)
	double E = par[1]; // the primary beam energy (MeV/u)
	double j = par[2]; // the primary beam current (enA)
	double T = par[3]; // temperature of Au (K)
	double irrad_time = par[4]; // beam irradiation time (from start) t_0
	double n_init = par[5]; // initial number N_AFr(t=0)
	double t = x[0]; // seconds

	double f_AFr = calc_flux(A,E,j,T); // Fr transported (pps) to the MCP
	double tau_AFr = tau_AZ(A,"Fr"); // seconds (lifetime)

	double tildeQ_A = f_AFr * tau_AFr;
	double tildeQ_AFr = n_init - tildeQ_A;

	if (t < irrad_time){
		return tildeQ_A + tildeQ_AFr * TMath::Exp(-t/tau_AFr);
	}else{
		double n_beamoff = tildeQ_A + tildeQ_AFr*TMath::Exp(-irrad_time/tau_AFr);
		double hatQ_AFr = n_beamoff;
		return hatQ_AFr * TMath::Exp(-(t-irrad_time)/tau_AFr);
	}
}

double AstatineAalpha(double *x,double *par){
	int A = par[0]; // the species of interest (204/205/206/207)
	double E = par[1]; // the primary beam energy (MeV/u)
	double j = par[2]; // the primary beam current (enA)
	double T = par[3]; // temperature of Au (K)
	double irrad_time = par[4]; // beam irradiation time (from start) t_0
	double n_init = par[5]; // initial number N_AFr(t=0)
	double n_init_at = par[6]; // initial number of At N_AAt(t=0)
	double t = x[0]; // seconds

	double f_Ap4Fr = calc_flux(A+4,E,j,T); // Fr transported (pps) to the MCP
	double tau_Ap4Fr = tau_AZ(A+4,"Fr"); // seconds (lifetime)
	double tau_AAt = tau_AZ(A,"At"); // seconds (lifetime)
	double Dt_AtFr = 1. / ((1./tau_AAt) - (1./tau_Ap4Fr)); // \Delta tau
	double b_Ap4Fr = b_AZ(A+4,"Fr"); // alpha-branching ratio

	double tildeQ_Ap4 = f_Ap4Fr * tau_Ap4Fr;
	double tildeQ_Ap4Fr = n_init - tildeQ_Ap4;
	double tildeP_A = b_Ap4Fr*tildeQ_Ap4*(tau_AAt/tau_Ap4Fr);
	double tildeP_Ap4Fr = b_Ap4Fr*tildeQ_Ap4Fr*(Dt_AtFr/tau_Ap4Fr);
	double tildeP_AAt = n_init_at - tildeP_A - tildeP_Ap4Fr;

	if (t < irrad_time){
		return tildeP_A + tildeP_Ap4Fr*TMath::Exp(-t/tau_Ap4Fr) + tildeP_AAt*TMath::Exp(-t/tau_AAt);
	}else{
		double n_beamoff = tildeP_A + tildeP_Ap4Fr*TMath::Exp(-irrad_time/tau_Ap4Fr) + tildeP_AAt*TMath::Exp(-irrad_time/tau_AAt);
		double hatQ_Ap4Fr = tildeQ_Ap4 + tildeQ_Ap4Fr*TMath::Exp(-irrad_time/tau_Ap4Fr);
		double hatP_Ap4Fr = b_Ap4Fr*hatQ_Ap4Fr*(Dt_AtFr/tau_Ap4Fr);
		double hatP_AAt = n_beamoff - hatP_Ap4Fr;
		return hatP_Ap4Fr*TMath::Exp(-(t-irrad_time)/tau_Ap4Fr) + hatP_AAt*TMath::Exp(-(t-irrad_time)/tau_AAt);
	}
}

double RadonA(double *x, double *par){
	int A = par[0]; // the species of interest (208/209/210/211)
	double E = par[1]; // the primary beam energy (MeV/u)
	double j = par[2]; // the primary beam current (enA)
	double T = par[3]; // temperature of Au (K)
	double irrad_time = par[4]; // beam irradiation time (from start) t_0
	double n_init = par[5]; // initial number N_AFr(t=0)
	double n_init_rn = par[6]; // initial number N_ARn(t=0)
	double t = x[0]; // seconds

	double f_AFr = calc_flux(A,E,j,T); // Fr transported (pps) to the MCP
	double tau_AFr = tau_AZ(A,"Fr"); // seconds (lifetime)
	double tau_ARn = tau_AZ(A,"Rn"); // seconds (lifetime)
	double Dt_RnFr = 1. / ((1./tau_ARn) - (1./tau_AFr)); // \Delta tau
	double b_AFr = b_AZ(A,"Fr"); // alpha-branching ratio

	double tildeQ_A = f_AFr * tau_AFr;
	double tildeQ_AFr = n_init - tildeQ_A;
	double tildeN_A = (1.-b_AFr)*tildeQ_A*(tau_ARn/tau_AFr);
	double tildeN_AFr = (1.-b_AFr)*tildeQ_AFr*(Dt_RnFr/tau_AFr);
	double tildeN_ARn = n_init_rn - tildeN_A - tildeN_AFr;

	if (t < irrad_time){
		return tildeN_A + tildeN_AFr*TMath::Exp(-t/tau_AFr) + tildeN_ARn*TMath::Exp(-t/tau_ARn);
	}else{
		double n_beamoff = tildeN_A + tildeN_AFr*TMath::Exp(-irrad_time/tau_AFr) + tildeN_ARn*TMath::Exp(-irrad_time/tau_ARn);
		double hatQ_AFr = tildeQ_A + tildeQ_AFr*TMath::Exp(-irrad_time/tau_AFr);
		double hatN_AFr = (1.-b_AFr)*hatQ_AFr*(Dt_RnFr/tau_AFr);
		double hatN_ARn = n_beamoff - hatN_AFr;
		return hatN_AFr*TMath::Exp(-(t-irrad_time)/tau_AFr) + hatN_ARn*TMath::Exp(-(t-irrad_time)/tau_ARn);
	}
}

double AstatineAbeta(double *x, double *par){
	int A = par[0]; // the species of interest (209/211)
	double E = par[1]; // the primary beam energy (MeV/u)
	double j = par[2]; // the primary beam current (enA)
	double T = par[3]; // temperature of Au (K)
	double irrad_time = par[4]; // beam irradiation time (from start) t_0
	double n_init = par[5]; // initial number N_AFr(t=0)
	double n_init_rn = par[6]; // initial number of Rn N_ARn(t=0)
	double n_init_at = par[7]; // initial number of At N_AAt(t=0)
	double t = x[0]; // seconds

	double f_AFr = calc_flux(A,E,j,T); // Fr transported (pps) to the MCP
	double tau_AFr = tau_AZ(A,"Fr"); // seconds (lifetime)
	double tau_ARn = tau_AZ(A,"Rn"); // seconds (lifetime)
	double tau_AAt = tau_AZ(A,"At"); // seconds (lifetime)
	double Dt_RnFr = 1. / ((1./tau_ARn) - (1./tau_AFr)); // \Delta tau
	double Dt_AtFr = 1. / ((1./tau_AAt) - (1./tau_AFr)); // \Delta tau
	double Dt_AtRn = 1. / ((1./tau_AAt) - (1./tau_ARn)); // \Delta tau
	double b_AFr = b_AZ(A,"Fr"); // alpha-branching ratio
	double b_ARn = b_AZ(A,"Rn"); // alpha-branching ratio

	double tildeQ_A = f_AFr * tau_AFr;
	double tildeQ_AFr = n_init - tildeQ_A;
	double tildeN_A = (1.-b_AFr)*tildeQ_A*(tau_ARn/tau_AFr);
	double tildeN_AFr = (1.-b_AFr)*tildeQ_AFr*(Dt_RnFr/tau_AFr);
	double tildeN_ARn = n_init_rn - tildeN_A - tildeN_AFr;
	double tildeM_A = (1.-b_ARn)*tildeN_A*(tau_AAt/tau_ARn);
	double tildeM_AFr = (1.-b_ARn)*tildeN_AFr*(Dt_AtFr/tau_ARn);
	double tildeM_ARn = (1.-b_ARn)*tildeN_ARn*(Dt_AtRn/tau_ARn);
	double tildeM_AAt = n_init_at - tildeM_A - tildeM_AFr - tildeM_ARn;

	if (t < irrad_time){
		return tildeM_A + tildeM_AFr*TMath::Exp(-t/tau_AFr) + tildeM_ARn*TMath::Exp(-t/tau_ARn) + tildeM_AAt*TMath::Exp(-t/tau_AAt);
	}else{
		double n_beamoff = tildeM_A + tildeM_AFr*TMath::Exp(-irrad_time/tau_AFr) + tildeM_ARn*TMath::Exp(-irrad_time/tau_ARn) + tildeM_AAt*TMath::Exp(-irrad_time/tau_AAt);
		double hatQ_AFr = tildeQ_A + tildeQ_AFr*TMath::Exp(-irrad_time/tau_AFr);
		double hatN_AFr = (1.-b_AFr)*hatQ_AFr*(Dt_RnFr/tau_AFr);
		double hatN_ARn = tildeN_A + tildeN_AFr*TMath::Exp(-irrad_time/tau_AFr) + tildeN_ARn*TMath::Exp(-irrad_time/tau_ARn) - hatN_AFr;
		double hatM_AFr = (1.-b_ARn)*hatN_AFr*(Dt_AtFr/tau_ARn);
		double hatM_ARn = (1.-b_ARn)*hatN_ARn*(Dt_AtRn/tau_ARn);
		double hatM_AAt = n_beamoff - hatM_AFr - hatM_ARn;
		return hatM_AFr*TMath::Exp(-(t-irrad_time)/tau_AFr) + hatM_ARn*TMath::Exp(-(t-irrad_time)/tau_ARn) + hatM_AAt*TMath::Exp(-(t-irrad_time)/tau_AAt);
	}
}

double Polonium206(double *x, double *par){
	int A = 206; // the species of interest (206)
	double E = par[0]; // the primary beam energy (MeV/u)
	double j = par[1]; // the primary beam current (enA)
	double T = par[2]; // temperature of Au (K)
	double irrad_time = par[3]; // beam irradiation time (from start) t_0
	double n_init = par[4]; // initial number N_A+4Fr(t=0)
	double n_init_rn = par[5]; // initial number of Rn N_A+4Rn(t=0)
	double n_init_at = par[6]; // initial number of At N_AAt(t=0)
	double n_init_po = par[7]; // initial number of At N_APo(t=0)
	double t = x[0]; // seconds

	double f_Ap4Fr = calc_flux(A+4,E,j,T); // Fr transported (pps) to the MCP
	double tau_Ap4Fr = tau_AZ(A+4,"Fr"); // seconds (lifetime)
	double tau_Ap4Rn = tau_AZ(A+4,"Rn"); // seconds (lifetime)
	double tau_AAt = tau_AZ(A,"At"); // seconds (lifetime)
	double tau_APo = tau_AZ(A,"Po"); // seconds (lifetime)
	double Dt_RnFr = 1. / ((1./tau_Ap4Rn) - (1./tau_Ap4Fr)); // \Delta tau
	double Dt_AtFr = 1. / ((1./tau_AAt) - (1./tau_Ap4Fr)); // \Delta tau
	double Dt_AtRn = 1. / ((1./tau_AAt) - (1./tau_Ap4Rn)); // \Delta tau
	double Dt_PoFr = 1. / ((1./tau_APo) - (1./tau_Ap4Fr)); // \Delta tau
	double Dt_PoRn = 1. / ((1./tau_APo) - (1./tau_Ap4Rn)); // \Delta tau
	double Dt_PoAt = 1. / ((1./tau_APo) - (1./tau_AAt)); // \Delta tau
	double b_Ap4Fr = b_AZ(A+4,"Fr"); // alpha-branching ratio
	double b_Ap4Rn = b_AZ(A+4,"Rn"); // alpha-branching ratio
	double b_AAt = b_AZ(A,"At"); // alpha-branching ratio

	double tildeQ_Ap4 = f_Ap4Fr * tau_Ap4Fr;
	double tildeQ_Ap4Fr = n_init - tildeQ_Ap4;
	double tildeP_A = b_Ap4Fr*tildeQ_Ap4*(tau_AAt/tau_Ap4Fr);
	double tildeP_Ap4Fr = b_Ap4Fr*tildeQ_Ap4Fr*(Dt_AtFr/tau_Ap4Fr);
	double tildeP_AAt = n_init_at - tildeP_A - tildeP_Ap4Fr;
	double tildeN_Ap4 = (1.-b_Ap4Fr)*tildeQ_Ap4*(tau_Ap4Rn/tau_Ap4Fr);
	double tildeN_Ap4Fr = (1.-b_Ap4Fr)*tildeQ_Ap4Fr*(Dt_RnFr/tau_Ap4Fr);
	double tildeN_Ap4Rn = n_init_rn - tildeN_Ap4 - tildeN_Ap4Fr;
	double tildeL_A = b_Ap4Rn*tildeN_Ap4*(tau_APo/tau_Ap4Rn) + (1.-b_AAt)*tildeP_A*(tau_APo/tau_AAt);
	double tildeL_Ap4Fr = b_Ap4Rn*tildeN_Ap4Fr*(Dt_PoFr/tau_Ap4Rn) + (1.-b_AAt)*tildeP_Ap4Fr*(Dt_PoFr/tau_AAt);
	double tildeL_Ap4Rn = b_Ap4Rn*tildeN_Ap4Rn*(Dt_PoRn/tau_Ap4Rn);
	double tildeL_AAt = (1.-b_AAt)*tildeP_AAt*(Dt_PoAt/tau_AAt);
	double tildeL_APo = n_init_po - tildeL_A - tildeL_Ap4Fr - tildeL_Ap4Rn - tildeL_AAt;

	if (t < irrad_time){
		return tildeL_A + tildeL_Ap4Fr*TMath::Exp(-t/tau_Ap4Fr) + tildeL_Ap4Rn*TMath::Exp(-t/tau_Ap4Rn) + tildeL_AAt*TMath::Exp(-t/tau_AAt) + tildeL_APo*TMath::Exp(-t/tau_APo);
	}else{
		double n_beamoff = tildeL_A + tildeL_Ap4Fr*TMath::Exp(-irrad_time/tau_Ap4Fr) + tildeL_Ap4Rn*TMath::Exp(-irrad_time/tau_Ap4Rn) + tildeL_AAt*TMath::Exp(-irrad_time/tau_AAt) + tildeL_APo*TMath::Exp(-irrad_time/tau_APo);
		double hatQ_Ap4Fr = tildeQ_Ap4 + tildeQ_Ap4Fr*TMath::Exp(-irrad_time/tau_Ap4Fr);
		double hatP_Ap4Fr = b_Ap4Fr*hatQ_Ap4Fr*(Dt_AtFr/tau_Ap4Fr);
		double hatP_AAt = tildeP_A + tildeP_Ap4Fr*TMath::Exp(-irrad_time/tau_Ap4Fr) + tildeP_AAt*TMath::Exp(-irrad_time/tau_AAt) - hatP_Ap4Fr;
		double hatN_Ap4Fr = (1.-b_Ap4Fr)*hatQ_Ap4Fr*(Dt_RnFr/tau_Ap4Fr);
		double hatN_Ap4Rn = tildeN_Ap4 + tildeN_Ap4Fr*TMath::Exp(-irrad_time/tau_Ap4Fr) + tildeN_Ap4Rn*TMath::Exp(-irrad_time/tau_Ap4Rn) - hatN_Ap4Fr;
		double hatL_Ap4Fr = b_Ap4Rn*hatN_Ap4Fr*(Dt_PoFr/tau_Ap4Rn) + (1.-b_AAt)*hatP_Ap4Fr*(Dt_PoFr/tau_AAt);
		double hatL_Ap4Rn = b_Ap4Rn*hatN_Ap4Rn*(Dt_PoRn/tau_AAt);
		double hatL_AAt = (1.-b_AAt)*hatP_AAt*(Dt_PoAt/tau_AAt);
		double hatL_APo = n_beamoff - hatL_Ap4Fr - hatL_Ap4Rn - hatL_AAt;
		return hatL_Ap4Fr*TMath::Exp(-(t-irrad_time)/tau_Ap4Fr) + hatL_Ap4Rn*TMath::Exp(-(t-irrad_time)/tau_Ap4Rn) + hatL_AAt*TMath::Exp(-(t-irrad_time)/tau_AAt) + hatL_APo*TMath::Exp(-(t-irrad_time)/tau_APo);
	}
}

double Polonium211(double *x, double *par){
	int A = 211; // the species of interest (211)
	double E = par[0]; // the primary beam energy (MeV/u)
	double j = par[1]; // the primary beam current (enA)
	double T = par[2]; // temperature of Au (K)
	double irrad_time = par[3]; // beam irradiation time (from start) t_0
	double n_init = par[4]; // initial number N_AFr(t=0)
	double n_init_rn = par[5]; // initial number of Rn N_ARn(t=0)
	double n_init_at = par[6]; // initial number of At N_AAt(t=0)
	double n_init_po = par[7]; // initial number of Po N_APo(t=0)
	double t = x[0]; // seconds

	double f_AFr = calc_flux(A,E,j,T); // Fr transported (pps) to the MCP
	double tau_AFr = tau_AZ(A,"Fr"); // seconds (lifetime)
	double tau_ARn = tau_AZ(A,"Rn"); // seconds (lifetime)
	double tau_AAt = tau_AZ(A,"At"); // seconds (lifetime)
	double tau_APo = tau_AZ(A,"Po"); // seconds (lifetime)
	double Dt_RnFr = 1. / ((1./tau_ARn) - (1./tau_AFr)); // \Delta tau
	double Dt_AtFr = 1. / ((1./tau_AAt) - (1./tau_AFr)); // \Delta tau
	double Dt_AtRn = 1. / ((1./tau_AAt) - (1./tau_ARn)); // \Delta tau
	double Dt_PoFr = 1. / ((1./tau_APo) - (1./tau_AFr)); // \Delta tau
	double Dt_PoRn = 1. / ((1./tau_APo) - (1./tau_ARn)); // \Delta tau
	double Dt_PoAt = 1. / ((1./tau_APo) - (1./tau_AAt)); // \Delta tau
	double b_AFr = b_AZ(A,"Fr"); // alpha-branching ratio
	double b_ARn = b_AZ(A,"Rn"); // alpha-branching ratio
	double b_AAt = b_AZ(A,"At"); // alpha-branching ratio

	double tildeQ_A = f_AFr * tau_AFr;
	double tildeQ_AFr = n_init - tildeQ_A;
	double tildeN_A = (1.-b_AFr)*tildeQ_A*(tau_ARn/tau_AFr);
	double tildeN_AFr = (1.-b_AFr)*tildeQ_AFr*(Dt_RnFr/tau_AFr);
	double tildeN_ARn = n_init_rn - tildeN_A - tildeN_AFr;
	double tildeM_A = (1.-b_ARn)*tildeN_A*(tau_AAt/tau_ARn);
	double tildeM_AFr = (1.-b_ARn)*tildeN_AFr*(Dt_AtFr/tau_ARn);
	double tildeM_ARn = (1.-b_ARn)*tildeN_ARn*(Dt_AtRn/tau_ARn);
	double tildeM_AAt = n_init_at - tildeM_A - tildeM_AFr - tildeM_ARn;
	double tildeK_A = (1.-b_AAt)*tildeM_A*(tau_APo/tau_AAt);
	double tildeK_AFr = (1.-b_AAt)*tildeM_AFr*(Dt_PoFr/tau_AAt);
	double tildeK_ARn = (1.-b_AAt)*tildeM_ARn*(Dt_PoRn/tau_AAt);
	double tildeK_AAt = (1.-b_AAt)*tildeM_AAt*(Dt_PoAt/tau_AAt);
	double tildeK_APo = n_init_po - tildeK_A - tildeK_AFr - tildeK_ARn - tildeK_AAt;

	if (t < irrad_time){
		return tildeK_A + tildeK_AFr*TMath::Exp(-t/tau_AFr) + tildeK_ARn*TMath::Exp(-t/tau_ARn) + tildeK_AAt*TMath::Exp(-t/tau_AAt) + tildeK_APo*TMath::Exp(-t/tau_APo);
	}else{
		double n_beamoff = tildeK_A + tildeK_AFr*TMath::Exp(-irrad_time/tau_AFr) + tildeK_ARn*TMath::Exp(-irrad_time/tau_ARn) + tildeK_AAt*TMath::Exp(-irrad_time/tau_AAt) + tildeK_APo*TMath::Exp(-irrad_time/tau_APo);
		double hatQ_AFr = tildeQ_A + tildeQ_AFr*TMath::Exp(-irrad_time/tau_AFr);
		double hatN_AFr = (1.-b_AFr)*hatQ_AFr*(Dt_RnFr/tau_AFr);
		double hatN_ARn = tildeN_A + tildeN_AFr*TMath::Exp(-irrad_time/tau_AFr) + tildeN_ARn*TMath::Exp(-irrad_time/tau_ARn) - hatN_AFr;
		double hatM_AFr = (1.-b_ARn)*hatN_AFr*(Dt_AtFr/tau_ARn);
		double hatM_ARn = (1.-b_ARn)*hatN_ARn*(Dt_AtRn/tau_ARn);
		double hatM_AAt = tildeM_A + tildeM_AFr*TMath::Exp(-irrad_time/tau_AFr) + tildeM_ARn*TMath::Exp(-irrad_time/tau_ARn) + tildeM_AAt*TMath::Exp(-irrad_time/tau_AAt) - hatM_AFr - hatM_ARn;
		double hatK_AFr = (1.-b_AAt)*hatM_AFr*(Dt_PoFr/tau_AAt);
		double hatK_ARn = (1.-b_AAt)*hatM_ARn*(Dt_PoRn/tau_AAt);
		double hatK_AAt = (1.-b_AAt)*hatM_AAt*(Dt_PoAt/tau_AAt);
		double hatK_APo = n_beamoff - hatK_AFr - hatK_ARn - hatK_AAt;
		return hatK_AFr*TMath::Exp(-(t-irrad_time)/tau_AFr) + hatK_ARn*TMath::Exp(-(t-irrad_time)/tau_ARn) + hatK_AAt*TMath::Exp(-(t-irrad_time)/tau_AAt) + hatK_APo*TMath::Exp(-(t-irrad_time)/tau_APo);
	}
}

double detection_eff(int A, const char* Z){
	double br = b_AZ(A,Z);

	double att_eff = 0.37; // Open-area-ratio of the MCP-IN surface
	double geo_eff = 0.0042/2.; // detection efficiency of the SSD

	return br*att_eff*geo_eff;
}

int main (int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	double seconds = 1.0;
	double minutes = 60.0;
	double hours = 60.0 * minutes;
	double days = 24.0 * hours;
	double years = 365.0 * days;

	double beam_current = 4000.; // enA of 18-O-6+
	double primaryFlux = beam_current*TMath::Power(10.,-9)/(6.*1.6*TMath::Power(10.,-19)); // Particles per second: for 18-O
	double T = 900. + 273.; // K
	double irradiation_time = 15.*minutes;
	double timelimit = 1.*hours;

	double beam_energy = 6.94; // MeV/u --> 125 MeV injected to Be window

//	// Extraction efficiencies: from the Au target to the MCP
//	double ext_208 = extraction_eff(beam_energy,T,208);
//	double ext_209 = extraction_eff(beam_energy,T,209);
//	double ext_210 = extraction_eff(beam_energy,T,210);
//	double ext_211 = extraction_eff(beam_energy,T,211);

	double draw_normprod = normprod(Be_degraded(beam_energy),0);

	TCanvas *c1 = new TCanvas();
	c1->Divide(2,2);

	c1->cd(1);
	TLatex *ls_ratio = new TLatex();
	ls_ratio->SetTextAlign(12);
	ls_ratio->SetTextSize(0.05);
//	ls_ratio->DrawLatex(0.1,0.9,Form("Surface ionization probabilities at %g #circC",T-273.));
//	ls_ratio->DrawLatex(0.2,0.8,Form("Fr: %g %%",100.*ionization(T,E_wf_Au,E_ip_Fr,1./2.)));
//	ls_ratio->DrawLatex(0.2,0.7,Form("Rn: %g %%",100.*ionization(T,E_wf_Au,E_ip_Rn,4./1.)));
//	ls_ratio->DrawLatex(0.2,0.6,Form("At: %g %%",100.*ionization(T,E_wf_Au,E_ip_At,5./4.)));
//	ls_ratio->DrawLatex(0.2,0.5,Form("Po: %g %%",100.*ionization(T,E_wf_Au,E_ip_Po,4./5.)));
	ls_ratio->DrawLatex(0.1,0.9,Form("E_{{}^{18}O} = %3.2f#rightarrow%3.2f MeV, I_{{}^{18}O}: %3.2f e#muA (%3.2f#times10^{12} pps), T_{Au} = %3.0f #circC",beam_energy*18.,Be_degraded(beam_energy)*18.,beam_current*TMath::Power(10.,-3),primaryFlux*TMath::Power(10.,-12),T-273.));
	ls_ratio->DrawLatex(0.1,0.8,"f_{{}^{A}Fr} = #varepsilon_{transportation} #varepsilon_{desorption} #varepsilon_{ionization} #varepsilon_{escape} #frac{P_{{}^{A}Fr}}{j} I_{{}^{18}O}");
	ls_ratio->DrawLatex(0.15,0.7,Form("P_{{}^{208/209/210/211}Fr}: %3.1f / %3.1f / %3.1f / %3.1f [#times10^{6}/s]",normprod(Be_degraded(beam_energy),208)*primaryFlux*TMath::Power(10.,-6),normprod(Be_degraded(beam_energy),209)*primaryFlux*TMath::Power(10.,-6),normprod(Be_degraded(beam_energy),210)*primaryFlux*TMath::Power(10.,-6),normprod(Be_degraded(beam_energy),211)*primaryFlux*TMath::Power(10.,-6)));
	ls_ratio->DrawLatex(0.15,0.6,Form("#varepsilon_{escape}({}^{208/209/210/211}Fr): %3.1f / %3.1f / %3.1f / %3.1f %%",100.*y_Fr(tau_AZ(208,"Fr"),D_m(T,208),depth_Fr(Be_degraded(beam_energy),208)),100.*y_Fr(tau_AZ(209,"Fr"),D_m(T,209),depth_Fr(Be_degraded(beam_energy),209)),100.*y_Fr(tau_AZ(210,"Fr"),D_m(T,210),depth_Fr(Be_degraded(beam_energy),210)),100.*y_Fr(tau_AZ(211,"Fr"),D_m(T,211),depth_Fr(Be_degraded(beam_energy),211))));
	ls_ratio->DrawLatex(0.15,0.5,Form("#varepsilon_{ionization}: %3.1f%%, #varepsilon_{desorption}: %3.1f%% (?), #varepsilon_{transportation}: %3.1f%% (?)",100.*ionization(T,5.1,4.07,1./2.),100.,100.));
	ls_ratio->DrawLatex(0.15,0.4,Form("f_{{}^{208/209/210/211}Fr}: %3.1f / %3.1f / %3.1f / %3.1f [#times10^{6}/s]",calc_flux(208,beam_energy,beam_current,T)*TMath::Power(10.,-6),calc_flux(209,beam_energy,beam_current,T)*TMath::Power(10.,-6),calc_flux(210,beam_energy,beam_current,T)*TMath::Power(10.,-6),calc_flux(211,beam_energy,beam_current,T)*TMath::Power(10.,-6)));
	ls_ratio->DrawLatex(0.1,0.3,"Number of ions at MCP surface");
	ls_ratio->DrawLatex(0.1,0.2,"N_{{}^{A}Fr}(t) = f_{{}^{A}Fr}#tau_{{}^{A}Fr} + #left[ N_{{}^{A}Fr}(t_{ON}) - f_{{}^{A}Fr}#tau_{{}^{A}Fr} #right] e^{- #frac{t - t_{ON}}{#tau_{{}^{A}Fr}}} #rightarrow N_{{}^{A}Fr}(t_{OFF}) e^{- #frac{t - t_{OFF}}{#tau_{{}^{A}Fr}}}");
	ls_ratio->DrawLatex(0.1,0.1,"**We assume that only Fr is extracted from the target.");

	c1->cd(2);

	TF1 *N_208fr = new TF1("{}^{208}Fr",FranciumA,0.,timelimit,6);
	N_208fr->SetParameters(208,beam_energy,beam_current,T,irradiation_time,N_on(208,"Fr"));
	TGraph *g_N208fr = new TGraph(N_208fr);
	g_N208fr->SetMarkerColor(3);
	g_N208fr->SetLineColor(3);

	TF1 *N_209fr = new TF1("{}^{209}Fr",FranciumA,0.,timelimit,6);
	N_209fr->SetParameters(209,beam_energy,beam_current,T,irradiation_time,N_on(209,"Fr"));
	TGraph *g_N209fr = new TGraph(N_209fr);
	g_N209fr->SetMarkerColor(4);
	g_N209fr->SetLineColor(4);

	TF1 *N_210fr = new TF1("{}^{210}Fr",FranciumA,0.,timelimit,6);
	N_210fr->SetParameters(210,beam_energy,beam_current,T,irradiation_time,N_on(210,"Fr"));
	TGraph *g_N210fr = new TGraph(N_210fr);
	g_N210fr->SetMarkerColor(2);
	g_N210fr->SetLineColor(2);

	TF1 *N_211fr = new TF1("{}^{211}Fr",FranciumA,0.,timelimit,6);
	N_211fr->SetParameters(211,beam_energy,beam_current,T,irradiation_time,N_on(211,"Fr"));
	TGraph *g_N211fr = new TGraph(N_211fr);
	g_N211fr->SetMarkerColor(5);
	g_N211fr->SetLineColor(5);

	TF1 *N_208rn = new TF1("{}^{208}Rn",RadonA,0.,timelimit,7);
	N_208rn->SetParameters(208,beam_energy,beam_current,T,irradiation_time,N_on(208,"Fr"),N_on(208,"Rn"));
	TGraph *g_N208rn = new TGraph(N_208rn);
	g_N208rn->SetMarkerColor(3);
	g_N208rn->SetLineColor(3);
	g_N208rn->SetLineStyle(6);

	TF1 *N_209rn = new TF1("{}^{209}Rn",RadonA,0.,timelimit,7);
	N_209rn->SetParameters(209,beam_energy,beam_current,T,irradiation_time,N_on(209,"Fr"),N_on(209,"Rn"));
	TGraph *g_N209rn = new TGraph(N_209rn);
	g_N209rn->SetMarkerColor(4);
	g_N209rn->SetLineColor(4);
	g_N209rn->SetLineStyle(6);

	TF1 *N_210rn = new TF1("{}^{210}Rn",RadonA,0.,timelimit,7);
	N_210rn->SetParameters(210,beam_energy,beam_current,T,irradiation_time,N_on(210,"Fr"),N_on(210,"Rn"));
	TGraph *g_N210rn = new TGraph(N_210rn);
	g_N210rn->SetMarkerColor(2);
	g_N210rn->SetLineColor(2);
	g_N210rn->SetLineStyle(6);

	TF1 *N_211rn = new TF1("{}^{211}Rn",RadonA,0.,timelimit,7);
	N_211rn->SetParameters(211,beam_energy,beam_current,T,irradiation_time,N_on(211,"Fr"),N_on(211,"Rn"));
	TGraph *g_N211rn = new TGraph(N_211rn);
	g_N211rn->SetMarkerColor(5);
	g_N211rn->SetLineColor(5);
	g_N211rn->SetLineStyle(6);

	TF1 *N_204at = new TF1("{}^{204}At",AstatineAalpha,0.,timelimit,7);
	N_204at->SetParameters(204,beam_energy,beam_current,T,irradiation_time,N_on(208,"Fr"),N_on(204,"At"));
	TGraph *g_N204at = new TGraph(N_204at);
	g_N204at->SetMarkerColor(3);
	g_N204at->SetLineColor(3);
	g_N204at->SetLineStyle(2);

	TF1 *N_205at = new TF1("{}^{205}At",AstatineAalpha,0.,timelimit,7);
	N_205at->SetParameters(205,beam_energy,beam_current,T,irradiation_time,N_on(209,"Fr"),N_on(205,"At"));
	TGraph *g_N205at = new TGraph(N_205at);
	g_N205at->SetMarkerColor(4);
	g_N205at->SetLineColor(4);
	g_N205at->SetLineStyle(2);

	TF1 *N_206at = new TF1("{}^{206}At",AstatineAalpha,0.,timelimit,7);
	N_206at->SetParameters(206,beam_energy,beam_current,T,irradiation_time,N_on(210,"Fr"),N_on(206,"At"));
	TGraph *g_N206at = new TGraph(N_206at);
	g_N206at->SetMarkerColor(2);
	g_N206at->SetLineColor(2);
	g_N206at->SetLineStyle(2);

	TF1 *N_207at = new TF1("{}^{207}At",AstatineAalpha,0.,timelimit,7);
	N_207at->SetParameters(207,beam_energy,beam_current,T,irradiation_time,N_on(211,"Fr"),N_on(207,"At"));
	TGraph *g_N207at = new TGraph(N_207at);
	g_N207at->SetMarkerColor(5);
	g_N207at->SetLineColor(5);
	g_N207at->SetLineStyle(2);

	TF1 *N_209at = new TF1("{}^{209}At",AstatineAbeta,0.,timelimit,8);
	N_209at->SetParameters(209,beam_energy,beam_current,T,irradiation_time,N_on(209,"Fr"),N_on(209,"Rn"),N_on(209,"At"));
	TGraph *g_N209at = new TGraph(N_209at);
	g_N209at->SetMarkerColor(4);
	g_N209at->SetLineColor(4);
	g_N209at->SetLineStyle(4);

	TF1 *N_211at = new TF1("{}^{211}At",AstatineAbeta,0.,timelimit,8);
	N_211at->SetParameters(211,beam_energy,beam_current,T,irradiation_time,N_on(211,"Fr"),N_on(211,"Rn"),N_on(211,"At"));
	TGraph *g_N211at = new TGraph(N_211at);
	g_N211at->SetMarkerColor(5);
	g_N211at->SetLineColor(5);
	g_N211at->SetLineStyle(4);

	TF1 *N_206po = new TF1("{}^{206}Po",Polonium206,0.,timelimit,8);
	N_206po->SetParameters(beam_energy,beam_current,T,irradiation_time,N_on(210,"Fr"),N_on(210,"Rn"),N_on(206,"At"),N_on(206,"Po"));
	TGraph *g_N206po = new TGraph(N_206po);
	g_N206po->SetMarkerColor(2);
	g_N206po->SetLineColor(2);
	g_N206po->SetLineStyle(9);

	TF1 *N_211po = new TF1("{}^{211}Po",Polonium211,0.,timelimit,8);
	N_211po->SetParameters(beam_energy,beam_current,T,irradiation_time,N_on(211,"Fr"),N_on(211,"Rn"),N_on(211,"At"),N_on(211,"Po"));
	TGraph *g_N211po = new TGraph(N_211po);
	g_N211po->SetMarkerColor(5);
	g_N211po->SetLineColor(5);
	g_N206po->SetLineStyle(9);


	TMultiGraph *N = new TMultiGraph("N","Ions at the MCP;Time Elapsed (s);Ions N(t)");
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
	N->Add(g_N209at);
	N->Add(g_N211at);
	N->Add(g_N206po);
	N->Add(g_N211po);
	N->Draw("AL");

	cout << "Beam Energy = " << beam_energy*18. << "-->" << Be_degraded(beam_energy)*18. << " MeV, Beam Current = " << beam_current*TMath::Power(10.,-9)/(6.*1.6*TMath::Power(10.,-19)) << " pps" << endl;

	TLine *limit_208 = new TLine(irradiation_time-120.,calc_flux(208,beam_energy,beam_current,T)*tau_AZ(208,"Fr"),irradiation_time+120.,calc_flux(208,beam_energy,beam_current,T)*tau_AZ(208,"Fr"));
	limit_208->SetLineColor(3);
	limit_208->SetLineWidth(2);
	limit_208->SetLineStyle(2);
//	limit_208->Draw();
	cout << "f_{208Fr} = " << calc_flux(208,beam_energy,beam_current,T) << " pps --> ";
	cout << "f*tau(208Fr): " << calc_flux(208,beam_energy,beam_current,T)*tau_AZ(208,"Fr") << " ions" << endl;

	TLine *limit_209 = new TLine(irradiation_time-120.,calc_flux(209,beam_energy,beam_current,T)*tau_AZ(209,"Fr"),irradiation_time-120.,calc_flux(209,beam_energy,beam_current,T)*tau_AZ(209,"Fr"));
	limit_209->SetLineColor(4);
	limit_209->SetLineWidth(2);
	limit_209->SetLineStyle(2);
//	limit_209->Draw();
	cout << "f_{209Fr} = " << calc_flux(209,beam_energy,beam_current,T) << " pps --> ";
	cout << "f*tau(209Fr): " << calc_flux(209,beam_energy,beam_current,T)*tau_AZ(209,"Fr") << " ions" << endl;

	TLine *limit_210 = new TLine(irradiation_time-120.,calc_flux(210,beam_energy,beam_current,T)*tau_AZ(210,"Fr"),irradiation_time-120.,calc_flux(210,beam_energy,beam_current,T)*tau_AZ(210,"Fr"));
	limit_210->SetLineColor(2);
	limit_210->SetLineWidth(2);
	limit_210->SetLineStyle(2);
//	limit_210->Draw();
	cout << "f_{210Fr} = " << calc_flux(210,beam_energy,beam_current,T) << " pps --> ";
	cout << "f*tau(210Fr): " << calc_flux(210,beam_energy,beam_current,T)*tau_AZ(210,"Fr") << " ions" << endl;

	TLine *limit_211 = new TLine(irradiation_time-120.,calc_flux(211,beam_energy,beam_current,T)*tau_AZ(211,"Fr"),irradiation_time-120.,calc_flux(211,beam_energy,beam_current,T)*tau_AZ(211,"Fr"));
	limit_211->SetLineColor(5);
	limit_211->SetLineWidth(2);
	limit_211->SetLineStyle(2);
//	limit_211->Draw();
	cout << "f_{211Fr} = " << calc_flux(211,beam_energy,beam_current,T) << " pps --> ";
	cout << "f*tau(211Fr): " << calc_flux(211,beam_energy,beam_current,T)*tau_AZ(211,"Fr") << " ions" << endl;



	c1->cd(1);
//	ls_ratio->DrawLatex(0.1,0.6,"Escape efficiencies #varepsilon_{escape}:");
//	ls_ratio->DrawLatex(0.15,0.5,Form("{}^{208}Fr: %3.2f%%, {}^{209}Fr: %3.2f%%, {}^{210}Fr: %3.2f%%, {}^{211}Fr: %3.2f%%",100.*escape_eff_208,100.*escape_eff_209,100.*escape_eff_210,100.*escape_eff_211));


	c1->cd(3);

	double e_208fr = 6.641; // MeV
	TGraph *g_alpha208fr = new TGraph(g_N208fr->GetN());
	for (int i=0; i<g_N208fr->GetN(); ++i) {
		g_alpha208fr->GetX()[i] = g_N208fr->GetX()[i];
		g_alpha208fr->GetY()[i] = g_N208fr->GetY()[i]*detection_eff(208,"Fr")/tau_AZ(208,"Fr");
	}
	g_alpha208fr->SetTitle(Form("{}^{208}Fr: %g MeV",e_208fr));
	g_alpha208fr->SetMarkerColor(3);
	g_alpha208fr->SetLineColor(3);
	g_alpha208fr->SetLineWidth(2);

	double e_209fr = 6.646; // MeV
	TGraph *g_alpha209fr = new TGraph(g_N209fr->GetN());
	for (int i=0; i<g_N209fr->GetN(); ++i) {
		g_alpha209fr->GetX()[i] = g_N209fr->GetX()[i];
		g_alpha209fr->GetY()[i] = g_N209fr->GetY()[i]*detection_eff(209,"Fr")/tau_AZ(209,"Fr");
	}
	g_alpha209fr->SetTitle(Form("{}^{209}Fr: %g MeV",e_209fr));
	g_alpha209fr->SetMarkerColor(4);
	g_alpha209fr->SetLineColor(4);
	g_alpha209fr->SetLineWidth(2);

	double e_210fr = 6.545; // MeV
	TGraph *g_alpha210fr = new TGraph(g_N210fr->GetN());
	for (int i=0; i<g_N210fr->GetN(); ++i) {
		g_alpha210fr->GetX()[i] = g_N210fr->GetX()[i];
		g_alpha210fr->GetY()[i] = g_N210fr->GetY()[i]*detection_eff(210,"Fr")/tau_AZ(210,"Fr");
	}
	g_alpha210fr->SetTitle(Form("{}^{210}Fr: %g MeV",e_210fr));
	g_alpha210fr->SetMarkerColor(2);
	g_alpha210fr->SetLineColor(2);
	g_alpha210fr->SetLineWidth(2);

	double e_211fr = 6.537; // MeV
	TGraph *g_alpha211fr = new TGraph(g_N211fr->GetN());
	for (int i=0; i<g_N211fr->GetN(); ++i) {
		g_alpha211fr->GetX()[i] = g_N211fr->GetX()[i];
		g_alpha211fr->GetY()[i] = g_N211fr->GetY()[i]*detection_eff(211,"Fr")/tau_AZ(211,"Fr");
	}
	g_alpha211fr->SetTitle(Form("{}^{211}Fr: %g MeV",e_211fr));
	g_alpha211fr->SetMarkerColor(5);
	g_alpha211fr->SetLineColor(5);
	g_alpha211fr->SetLineWidth(2);

	double e_208rn = 6.1401; // MeV
	TGraph *g_alpha208rn = new TGraph(g_N208rn->GetN());
	for (int i=0; i<g_N208rn->GetN(); ++i) {
		g_alpha208rn->GetX()[i] = g_N208rn->GetX()[i];
		g_alpha208rn->GetY()[i] = g_N208rn->GetY()[i]*detection_eff(208,"Rn")/tau_AZ(208,"Rn");
	}
	g_alpha208rn->SetTitle(Form("{}^{208}Rn: %g MeV",e_208rn));
	g_alpha208rn->SetMarkerColor(3);
	g_alpha208rn->SetLineColor(3);
	g_alpha208rn->SetLineWidth(2);
	g_alpha208rn->SetLineStyle(6);

	double e_209rn = 6.039; // MeV
	TGraph *g_alpha209rn = new TGraph(g_N209rn->GetN());
	for (int i=0; i<g_N209rn->GetN(); ++i) {
		g_alpha209rn->GetX()[i] = g_N209rn->GetX()[i];
		g_alpha209rn->GetY()[i] = g_N209rn->GetY()[i]*detection_eff(209,"Rn")/tau_AZ(209,"Rn");
	}
	g_alpha209rn->SetTitle(Form("{}^{209}Rn: %g MeV",e_209rn));
	g_alpha209rn->SetMarkerColor(4);
	g_alpha209rn->SetLineColor(4);
	g_alpha209rn->SetLineWidth(2);
	g_alpha209rn->SetLineStyle(6);

	double e_210rn = 6.041; // MeV
	TGraph *g_alpha210rn = new TGraph(g_N210rn->GetN());
	for (int i=0; i<g_N210rn->GetN(); ++i) {
		g_alpha210rn->GetX()[i] = g_N210rn->GetX()[i];
		g_alpha210rn->GetY()[i] = g_N210rn->GetY()[i]*detection_eff(210,"Rn")/tau_AZ(210,"Rn");
	}
	g_alpha210rn->SetTitle(Form("{}^{210}Rn: %g MeV",e_210rn));
	g_alpha210rn->SetMarkerColor(2);
	g_alpha210rn->SetLineColor(2);
	g_alpha210rn->SetLineWidth(2);
	g_alpha210rn->SetLineStyle(6);

	double e_211rn = 5.7839; // MeV
	TGraph *g_alpha211rn = new TGraph(g_N211rn->GetN());
	for (int i=0; i<g_N211rn->GetN(); ++i) {
		g_alpha211rn->GetX()[i] = g_N211rn->GetX()[i];
		g_alpha211rn->GetY()[i] = g_N211rn->GetY()[i]*detection_eff(211,"Rn")/tau_AZ(211,"Rn");
	}
	g_alpha211rn->SetTitle(Form("{}^{211}Rn: %g MeV",e_211rn));
	g_alpha211rn->SetMarkerColor(5);
	g_alpha211rn->SetLineColor(5);
	g_alpha211rn->SetLineWidth(2);
	g_alpha211rn->SetLineStyle(6);

	double e_204at = 5.9503; // MeV
	TGraph *g_alpha204at = new TGraph(g_N204at->GetN());
	for (int i=0; i<g_N204at->GetN(); ++i) {
		g_alpha204at->GetX()[i] = g_N204at->GetX()[i];
		g_alpha204at->GetY()[i] = g_N204at->GetY()[i]*detection_eff(204,"At")/tau_AZ(204,"At");
	}
	g_alpha204at->SetTitle(Form("{}^{204}At: %g MeV",e_204at));
	g_alpha204at->SetMarkerColor(3);
	g_alpha204at->SetLineColor(3);
	g_alpha204at->SetLineWidth(2);
	g_alpha204at->SetLineStyle(2);

	double e_205at = 5.902; // MeV
	TGraph *g_alpha205at = new TGraph(g_N205at->GetN());
	for (int i=0; i<g_N205at->GetN(); ++i) {
		g_alpha205at->GetX()[i] = g_N205at->GetX()[i];
		g_alpha205at->GetY()[i] = g_N205at->GetY()[i]*detection_eff(205,"At")/tau_AZ(205,"At");
	}
	g_alpha205at->SetTitle(Form("{}^{205}At: %g MeV",e_205at));
	g_alpha205at->SetMarkerColor(4);
	g_alpha205at->SetLineColor(4);
	g_alpha205at->SetLineWidth(2);
	g_alpha205at->SetLineStyle(2);

	double e_206at = 5.7026; // MeV
	TGraph *g_alpha206at = new TGraph(g_N206at->GetN());
	for (int i=0; i<g_N206at->GetN(); ++i) {
		g_alpha206at->GetX()[i] = g_N206at->GetX()[i];
		g_alpha206at->GetY()[i] = g_N206at->GetY()[i]*detection_eff(206,"At")/tau_AZ(206,"At");
	}
	g_alpha206at->SetTitle(Form("{}^{206}At: %g MeV",e_206at));
	g_alpha206at->SetMarkerColor(2);
	g_alpha206at->SetLineColor(2);
	g_alpha206at->SetLineWidth(2);
	g_alpha206at->SetLineStyle(2);

	double e_207at = 5.758; // MeV
	TGraph *g_alpha207at = new TGraph(g_N207at->GetN());
	for (int i=0; i<g_N207at->GetN(); ++i) {
		g_alpha207at->GetX()[i] = g_N207at->GetX()[i];
		g_alpha207at->GetY()[i] = g_N207at->GetY()[i]*detection_eff(207,"At")/tau_AZ(207,"At");
	}
	g_alpha207at->SetTitle(Form("{}^{207}At: %g MeV",e_207at));
	g_alpha207at->SetMarkerColor(5);
	g_alpha207at->SetLineColor(5);
	g_alpha207at->SetLineWidth(2);
	g_alpha207at->SetLineStyle(2);

	double e_209at = 5.647; // MeV
	TGraph *g_alpha209at = new TGraph(g_N209at->GetN());
	for (int i=0; i<g_N209at->GetN(); ++i) {
		g_alpha209at->GetX()[i] = g_N209at->GetX()[i];
		g_alpha209at->GetY()[i] = g_N209at->GetY()[i]*detection_eff(209,"At")/tau_AZ(209,"At");
	}
	g_alpha209at->SetTitle(Form("{}^{209}At: %g MeV",e_209at));
	g_alpha209at->SetMarkerColor(4);
	g_alpha209at->SetLineColor(4);
	g_alpha209at->SetLineWidth(2);
	g_alpha209at->SetLineStyle(4);

	double e_211at = 5.8695; // MeV
	TGraph *g_alpha211at = new TGraph(g_N211at->GetN());
	for (int i=0; i<g_N211at->GetN(); ++i) {
		g_alpha211at->GetX()[i] = g_N211at->GetX()[i];
		g_alpha211at->GetY()[i] = g_N211at->GetY()[i]*detection_eff(211,"At")/tau_AZ(211,"At");
	}
	g_alpha211at->SetTitle(Form("{}^{211}At: %g MeV",e_211at));
	g_alpha211at->SetMarkerColor(5);
	g_alpha211at->SetLineColor(5);
	g_alpha211at->SetLineWidth(2);
	g_alpha211at->SetLineStyle(4);

	double e_206po = 5.2237; // MeV
	TGraph *g_alpha206po = new TGraph(g_N206po->GetN());
	for (int i=0; i<g_N206po->GetN(); ++i) {
		g_alpha206po->GetX()[i] = g_N206po->GetX()[i];
		g_alpha206po->GetY()[i] = g_N206po->GetY()[i]*detection_eff(206,"Po")/tau_AZ(206,"Po");
	}
	g_alpha206po->SetTitle(Form("{}^{206}Po: %g MeV",e_206po));
	g_alpha206po->SetMarkerColor(2);
	g_alpha206po->SetLineColor(2);
	g_alpha206po->SetLineWidth(2);
	g_alpha206po->SetLineStyle(9);

	double e_211po = 7.4503; // MeV
	TGraph *g_alpha211po = new TGraph(g_N211po->GetN());
	for (int i=0; i<g_N211po->GetN(); ++i) {
		g_alpha211po->GetX()[i] = g_N211po->GetX()[i];
		g_alpha211po->GetY()[i] = g_N211po->GetY()[i]*detection_eff(211,"Po")/tau_AZ(211,"Po");
	}
	g_alpha211po->SetTitle(Form("{}^{211}Po: %g MeV",e_211po));
	g_alpha211po->SetMarkerColor(5);
	g_alpha211po->SetLineColor(5);
	g_alpha211po->SetLineWidth(5);
	g_alpha211po->SetLineStyle(9);

	TMultiGraph *alpha = new TMultiGraph("alpha","Flux of #alpha Particles Detected at the SSD; Time Elapsed (s); Detected #alpha Particles Y (/s)");
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
	alpha->Add(g_alpha209at);
	alpha->Add(g_alpha211at);
	alpha->Add(g_alpha206po);
	alpha->Add(g_alpha211po);
	alpha->Draw("AL");

	c1->cd(4);

//	double flux_net = flux_208fr+flux_209fr+flux_210fr+flux_211fr;
	double alpharate_208fr = g_alpha208fr->Eval(irradiation_time);
	double alpharate_209fr = g_alpha209fr->Eval(irradiation_time);
	double alpharate_210fr = g_alpha210fr->Eval(irradiation_time);
	double alpharate_211fr = g_alpha211fr->Eval(irradiation_time);
//	double alpharate_net = alpharate_208fr+alpharate_209fr+alpharate_210fr+alpharate_211fr;

	TLatex l;
	l.SetTextAlign(12);
	l.SetTextSize(0.05);
	l.DrawLatex(0.1,0.9,"Detected signal Y_{{}^{A}Fr}(t) = #varepsilon_{signal} #varepsilon_{SSD} (1-#varepsilon_{OAR}) #frac{b_{{}^{A}Fr}}{#tau_{{}^{A}Fr}} N_{{}^{A}Fr}(t)");
	l.DrawLatex(0.1,0.8,Form("#varepsilon_{signal} = %3.1f%% (?), #varepsilon_{SSD} = %3.2f%%, #varepsilon_{OAR} = %3.1f%%",100.,0.2080,63.));
	l.DrawLatex(0.1,0.7,Form("After %3.2f seconds of beam irradiation:",irradiation_time));
	l.DrawLatex(0.2,0.6,Form("{}^{208}Fr: #alpha detection %3.2f pps (%3.2f MeV)",alpharate_208fr,e_208fr));
	l.DrawLatex(0.2,0.5,Form("{}^{209}Fr: #alpha detection %3.2f pps (%3.2f MeV)",alpharate_209fr,e_209fr));
	l.DrawLatex(0.2,0.4,Form("{}^{210}Fr: #alpha detection %3.2f pps (%3.2f MeV)",alpharate_210fr,e_210fr));
	l.DrawLatex(0.2,0.3,Form("{}^{211}Fr: #alpha detection %3.2f pps (%3.2f MeV)",alpharate_211fr,e_211fr));



	c1->cd(2)->BuildLegend();
	c1->cd(3)->BuildLegend();
	c1->Update();
	c1->Modified();

	rootapp.Run();

	return 0;

}
