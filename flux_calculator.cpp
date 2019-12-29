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
//		l_208->Draw();
		TLine *l_209 = new TLine(81.,s209->Eval(E),119.,s209->Eval(E));
		l_209->SetLineColor(4);
		l_209->SetLineStyle(2);
		l_209->SetLineWidth(2);
//		l_209->Draw();
		TLine *l_210 = new TLine(81.,s210->Eval(E),119.,s210->Eval(E));
		l_210->SetLineColor(2);
		l_210->SetLineStyle(2);
		l_210->SetLineWidth(2);
//		l_210->Draw();
		TLine *l_211 = new TLine(81.,s211->Eval(E),119.,s211->Eval(E));
		l_211->SetLineColor(5);
		l_211->SetLineStyle(2);
		l_211->SetLineWidth(2);
//		l_211->Draw();
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
	double d = -2.6 + 0.34 * incident_energy*18. - 0.31 * peak_energy; // um
	if (d < 0.0){
		return -9999.;
	}else{
		return d;
	}
}

double depth_Fr(double incident_energy, int isotope){
	// energy: MeV/u (injection energy into Au target)
	if (isotope == 208){
		double peak = 113.; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 209){
		double peak = 102.; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 210){
		double peak = 91.; // MeV based on stancari2006 plot
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

double escape_time(double incident_energy, double temp, int A){
	// based on Fujioka1981
	double t_ave = depth_Fr(incident_energy,A)*depth_Fr(incident_energy,A)/(3.*D_m(temp,A));
	if (depth_Fr(incident_energy,A) < 0.0){
		return 0.0;
	}else{
		return t_ave;
	}
}

void draw_extraction_eff(double temp){
	int npts = 50;
	double emin = 80.; // MeV (primary)
	double emax = 130.; // MeV (primary)
	double energy[npts];
	double eff208[npts];
	double eff209[npts];
	double eff210[npts];
	double eff211[npts];
	double esc208[npts];
	double esc209[npts];
	double esc210[npts];
	double esc211[npts];

	for (int i=0; i<npts; ++i){
		energy[i] = emin + double(i)*(emax-emin)/double(npts);
		eff208[i] = 100.*extraction_eff(Be_degraded(energy[i]/18.),temp,208);
		eff209[i] = 100.*extraction_eff(Be_degraded(energy[i]/18.),temp,209);
		eff210[i] = 100.*extraction_eff(Be_degraded(energy[i]/18.),temp,210);
		eff211[i] = 100.*extraction_eff(Be_degraded(energy[i]/18.),temp,211);
		esc208[i] = escape_time(Be_degraded(energy[i]/18.),temp,208);
		esc209[i] = escape_time(Be_degraded(energy[i]/18.),temp,209);
		esc210[i] = escape_time(Be_degraded(energy[i]/18.),temp,210);
		esc211[i] = escape_time(Be_degraded(energy[i]/18.),temp,211);
	}

	TGraph *g_eff208 = new TGraph(npts,energy,eff208);
	g_eff208->SetTitle("{}^{208}Fr");
	g_eff208->SetLineColor(3);
	g_eff208->SetLineWidth(2);
	TGraph *g_eff209 = new TGraph(npts,energy,eff209);
	g_eff209->SetTitle("{}^{209}Fr");
	g_eff209->SetLineColor(4);
	g_eff209->SetLineWidth(2);
	TGraph *g_eff210 = new TGraph(npts,energy,eff210);
	g_eff210->SetTitle("{}^{210}Fr");
	g_eff210->SetLineColor(2);
	g_eff210->SetLineWidth(2);
	TGraph *g_eff211 = new TGraph(npts,energy,eff211);
	g_eff211->SetTitle("{}^{211}Fr");
	g_eff211->SetLineColor(5);
	g_eff211->SetLineWidth(2);

	TMultiGraph *mg_eff = new TMultiGraph();
	mg_eff->SetTitle(Form("Escape Probability Before Decay at T_{Au} = %3.0f#circC;Primary Beam Energy E_{{}^{18}O} (MeV);Probability #varepsilon_{extraction} = #varepsilon_{transportation} #varepsilon_{desorption} #varepsilon_{ionization} #varepsilon_{escape} (%%)",temp-273.));
	mg_eff->Add(g_eff208);
	mg_eff->Add(g_eff209);
	mg_eff->Add(g_eff210);
	mg_eff->Add(g_eff211);

	TCanvas *c_extraction_efficiency = new TCanvas();
	c_extraction_efficiency->Divide(2,1);
	c_extraction_efficiency->cd(1);
	mg_eff->Draw("ALP");
	c_extraction_efficiency->cd(1)->BuildLegend();

	TGraph *g_esc208 = new TGraph(npts,energy,esc208);
	g_esc208->SetTitle("{}^{208}Fr");
	g_esc208->SetLineColor(3);
	g_esc208->SetLineWidth(2);
	TGraph *g_esc209 = new TGraph(npts,energy,esc209);
	g_esc209->SetTitle("{}^{209}Fr");
	g_esc209->SetLineColor(4);
	g_esc209->SetLineWidth(2);
	TGraph *g_esc210 = new TGraph(npts,energy,esc210);
	g_esc210->SetTitle("{}^{210}Fr");
	g_esc210->SetLineColor(2);
	g_esc210->SetLineWidth(2);
	TGraph *g_esc211 = new TGraph(npts,energy,esc211);
	g_esc211->SetTitle("{}^{211}Fr");
	g_esc211->SetLineColor(5);
	g_esc211->SetLineWidth(2);

	TMultiGraph *mg_esc = new TMultiGraph();
	mg_esc->SetTitle(Form("Average Escape Time at T_{Au} = %3.0f#circC;Primary Beam Energy E_{{}^{18}O} (MeV);Escape Time <t> = #frac{d^{2}}{3D} (s)",temp-273.));
	mg_esc->Add(g_esc208);
	mg_esc->Add(g_esc209);
	mg_esc->Add(g_esc210);
	mg_esc->Add(g_esc211);

	c_extraction_efficiency->cd(2);
	mg_esc->Draw("ALP");
  c_extraction_efficiency->cd(2)->SetLogy();
	c_extraction_efficiency->cd(2)->BuildLegend();

	return;
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

double database(int timestamp, const char* dataset){
  double I = -9999.; // enA
  double T = -9999.; // K
  if (timestamp == 122112){
    I = 1000.;
    T = 326.;
  }else if (timestamp == 141651){
    I = 3500.;
    T = 687.;
  }else if (timestamp == 145015){
    I = 3500.;
    T = 405.;
  }else if (timestamp == 175332){
    I = 2007.;
    T = 737.;
  }else if (timestamp == 183718){
    I = 2007.;
    T = 879.;
  }else if (timestamp == 190330){
    I = 2007.;
    T = 960.;
  }else if (timestamp == 204011){
    I = 3500.;
    T = 670.;
  }else if (timestamp == 205046){
    I = 4040.;
    T = 670.;
  }

  if (dataset == "beam_current"){
    return I;
  }else if (dataset == "T"){
    return T+273.;
  }
}

int main(int argc, char** argv){

  int timestamp = 205046;

	TRint rootapp("app",&argc,argv);

	double seconds = 1.0;
	double minutes = 60.0;
	double hours = 60.0 * minutes;
	double days = 24.0 * hours;
	double years = 365.0 * days;

	double beam_current = database(timestamp,"beam_current"); // enA of 18-O-6+
	double primaryFlux = beam_current*TMath::Power(10.,-9)/(6.*1.6*TMath::Power(10.,-19)); // Particles per second: for 18-O
	double T = database(timestamp,"T"); // K

//	double t_on = 30.0*seconds;
//	double irradiation_time = 15.*minutes;
//	double t_off = t_on+irradiation_time;
//	double timelimit = 20.*minutes;
	double beam_energy = 6.94; // MeV/u --> 125 MeV injected to Be window

//	// Extraction efficiencies: from the Au target to the MCP
//	double ext_208 = extraction_eff(beam_energy,T,208);
//	double ext_209 = extraction_eff(beam_energy,T,209);
//	double ext_210 = extraction_eff(beam_energy,T,210);
//	double ext_211 = extraction_eff(beam_energy,T,211);

	double draw_normprod = normprod(Be_degraded(beam_energy),0);
	draw_extraction_eff(T);

	TCanvas *c1 = new TCanvas();

  TLatex *ls_ratio = new TLatex();
	ls_ratio->SetTextAlign(12);
	ls_ratio->SetTextSize(0.03);
  ls_ratio->DrawLatex(0.1,0.9,Form("For Dataset 20191122-%d:",timestamp));
	ls_ratio->DrawLatex(0.1,0.8,Form("E_{{}^{18}O} = %3.2f#rightarrow%3.2f MeV, I_{{}^{18}O}: %3.2f e#muA (%3.2f#times10^{12} pps), T_{Au} = %3.0f #circC",beam_energy*18.,Be_degraded(beam_energy)*18.,beam_current*TMath::Power(10.,-3),primaryFlux*TMath::Power(10.,-12),T-273.));
	ls_ratio->DrawLatex(0.1,0.7,"f_{{}^{A}Fr} = #varepsilon_{transportation} #varepsilon_{desorption} #varepsilon_{ionization} #varepsilon_{escape} #frac{P_{{}^{A}Fr}}{j} I_{{}^{18}O}");
//	ls_ratio->DrawLatex(0.15,0.7,Form("P_{{}^{208/209/210/211}Fr}: %3.1f / %3.1f / %3.1f / %3.1f [#times10^{6}/s]",normprod(Be_degraded(beam_energy),208)*primaryFlux*TMath::Power(10.,-6),normprod(Be_degraded(beam_energy),209)*primaryFlux*TMath::Power(10.,-6),normprod(Be_degraded(beam_energy),210)*primaryFlux*TMath::Power(10.,-6),normprod(Be_degraded(beam_energy),211)*primaryFlux*TMath::Power(10.,-6)));
	ls_ratio->DrawLatex(0.15,0.5,Form("#varepsilon_{escape}({}^{208/209/210/211}Fr): %3.1f / %3.1f / %3.1f / %3.1f %%, #varepsilon_{ionization}: %3.1f%%",100.*y_Fr(tau_AZ(208,"Fr"),D_m(T,208),depth_Fr(Be_degraded(beam_energy),208)),100.*y_Fr(tau_AZ(209,"Fr"),D_m(T,209),depth_Fr(Be_degraded(beam_energy),209)),100.*y_Fr(tau_AZ(210,"Fr"),D_m(T,210),depth_Fr(Be_degraded(beam_energy),210)),100.*y_Fr(tau_AZ(211,"Fr"),D_m(T,211),depth_Fr(Be_degraded(beam_energy),211)),100.*ionization(T,5.1,4.07,1./2.)));
//	ls_ratio->DrawLatex(0.15,0.5,Form("#varepsilon_{ionization}: %3.1f%%, #varepsilon_{desorption}: %3.1f%% (?), #varepsilon_{transportation}: %3.1f%% (?)",100.*ionization(T,5.1,4.07,1./2.),100.,100.));
//	ls_ratio->DrawLatex(0.15,0.4,Form("f_{{}^{208/209/210/211}Fr}: %3.1f / %3.1f / %3.1f / %3.1f [#times10^{6}/s]",calc_flux(208,beam_energy,beam_current,T)*TMath::Power(10.,-6),calc_flux(209,beam_energy,beam_current,T)*TMath::Power(10.,-6),calc_flux(210,beam_energy,beam_current,T)*TMath::Power(10.,-6),calc_flux(211,beam_energy,beam_current,T)*TMath::Power(10.,-6)));
//	ls_ratio->DrawLatex(0.1,0.3,"Number of ions at MCP surface");
//	ls_ratio->DrawLatex(0.1,0.2,"N_{{}^{A}Fr}(t) = f_{{}^{A}Fr}#tau_{{}^{A}Fr} + #left[ N_{{}^{A}Fr}(t_{ON}) - f_{{}^{A}Fr}#tau_{{}^{A}Fr} #right] e^{- #frac{t - t_{ON}}{#tau_{{}^{A}Fr}}} #rightarrow N_{{}^{A}Fr}(t_{OFF}) e^{- #frac{t - t_{OFF}}{#tau_{{}^{A}Fr}}}");
//	ls_ratio->DrawLatex(0.1,0.4,"**We assume that only Fr is extracted from the target.");


//  ls_ratio->DrawLatex(0.1,0.9,Form("Calculated Yield for dataset %d:",timestamp));
  double p_208 = normprod(Be_degraded(beam_energy),208)*primaryFlux; // pps
  double p_209 = normprod(Be_degraded(beam_energy),209)*primaryFlux; // pps
  double p_210 = normprod(Be_degraded(beam_energy),210)*primaryFlux; // pps
  double p_211 = normprod(Be_degraded(beam_energy),211)*primaryFlux; // pps
  double f_208 = calc_flux(208,beam_energy,beam_current,T); // pps
  double f_209 = calc_flux(209,beam_energy,beam_current,T); // pps
  double f_210 = calc_flux(210,beam_energy,beam_current,T); // pps
  double f_211 = calc_flux(211,beam_energy,beam_current,T); // pps
  double b_208 = b_AZ(208,"Fr");
  double b_209 = b_AZ(209,"Fr");
  double b_210 = b_AZ(210,"Fr");
  double b_211 = b_AZ(211,"Fr");
  double sum_Pb = p_208*b_208 + p_209*b_209 + p_210*b_210 + p_211*b_211;
  ls_ratio->DrawLatex(0.15,0.6,Form("#sum_{A = 208}^{211} P_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",sum_Pb));
  double sum_fb = f_208*b_208 + f_209*b_209 + f_210*b_210 + f_211*b_211;
  ls_ratio->DrawLatex(0.15,0.4,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 100%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",sum_fb));
  ls_ratio->DrawLatex(0.15,0.3,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 80%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",0.8*sum_fb));
  ls_ratio->DrawLatex(0.15,0.2,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 60%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",0.6*sum_fb));
  ls_ratio->DrawLatex(0.15,0.1,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 40%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",0.4*sum_fb));


  rootapp.Run();

  return 0;
}
