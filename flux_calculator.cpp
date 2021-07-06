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
double b_AZ(int A, const char* Z, const char* data = "0"){
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

	if (data == "error"){
		return b_err;
	}else{
		return b;
	}
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
  }else{
	return -9999.;
  }
}

int main(int argc, char** argv){

  int timestamp = 190330;

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
//  ls_ratio->DrawLatex(0.15,0.4,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 100%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",sum_fb));
//  ls_ratio->DrawLatex(0.15,0.3,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 80%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",0.8*sum_fb));
//  ls_ratio->DrawLatex(0.15,0.2,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 60%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",0.6*sum_fb));
//  ls_ratio->DrawLatex(0.15,0.1,Form("#varepsilon_{transportation} #varepsilon_{desorption} = 40%%: #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f pps",0.4*sum_fb));

	double p_210_flux = p_210/TMath::Sqrt(primaryFlux);
	double p_210_energyh = normprod(Be_degraded(beam_energy*1.03),210)*primaryFlux - p_210;
	double p_210_energyl = normprod(Be_degraded(beam_energy*0.97),210)*primaryFlux - p_210;
	double p_210_cs = 0.2 * p_210;

	double varesc = y_Fr(tau_AZ(210,"Fr"),D_m(T,210),depth_Fr(Be_degraded(beam_energy),210));
	double varesc_energyh = y_Fr(tau_AZ(210,"Fr"),D_m(T,210),depth_Fr(Be_degraded(beam_energy*0.97),210)) - varesc;
	double varesc_energyl = varesc - y_Fr(tau_AZ(210,"Fr"),D_m(T,210),depth_Fr(Be_degraded(beam_energy*1.03),210));
	double varesc_Th = y_Fr(tau_AZ(210,"Fr"),D_m(T+200.,210),depth_Fr(Be_degraded(beam_energy),210)) - varesc;
	double varesc_Tl = varesc - y_Fr(tau_AZ(210,"Fr"),D_m(T-200.,210),depth_Fr(Be_degraded(beam_energy),210));
	double varion = ionization(T,5.1,4.07,1./2.);
	double varion_Th = ionization(T-200.,5.1,4.07,1./2.) - varion;
	double varion_Tl = varion - ionization(T+200.,5.1,4.07,1./2.);
	double f_208_flux = calc_flux(208,beam_energy,beam_current+beam_current/TMath::Sqrt(primaryFlux),T) - f_208;
	double f_208_energyh = calc_flux(208,beam_energy*1.03,beam_current,T) - f_208;
	double f_208_energyl = f_208 - calc_flux(208,beam_energy*0.97,beam_current,T);
	double f_208_Th = calc_flux(208,beam_energy,beam_current,T+200.) - f_208;
	double f_208_Tl = f_208 - calc_flux(208,beam_energy,beam_current,T-200.);
	double f_208_cs = 0.2 * f_208;
	double b_208_err = b_AZ(208,"Fr","error");
	double f_209_flux = calc_flux(209,beam_energy,beam_current+beam_current/TMath::Sqrt(primaryFlux),T) - f_209;
	double f_209_energyh = calc_flux(209,beam_energy*0.97,beam_current,T) - f_209;
	double f_209_energyl = f_209 - calc_flux(209,beam_energy*1.03,beam_current,T);
	double f_209_Th = calc_flux(209,beam_energy,beam_current,T+200.) - f_209;
	double f_209_Tl = f_209 - calc_flux(209,beam_energy,beam_current,T-200.);
	double f_209_cs = 0.2 * f_209;
	double b_209_err = b_AZ(209,"Fr","error");
	double f_210_flux = calc_flux(210,beam_energy,beam_current+beam_current/TMath::Sqrt(primaryFlux),T) - f_210;
	double f_210_energyh = calc_flux(210,beam_energy*0.97,beam_current,T) - f_210;
	double f_210_energyl = f_210 - calc_flux(210,beam_energy*1.03,beam_current,T);
	double f_210_Th = calc_flux(210,beam_energy,beam_current,T+200.) - f_210;
	double f_210_Tl = f_210 - calc_flux(210,beam_energy,beam_current,T-200.);
	double f_210_cs = 0.2 * f_210;
	double b_210_err = b_AZ(210,"Fr","error");
	double f_211_flux = calc_flux(211,beam_energy,beam_current+beam_current/TMath::Sqrt(primaryFlux),T) - f_211;
	double f_211_energyh = calc_flux(211,beam_energy*0.97,beam_current,T) - f_211;
	double f_211_energyl = f_211 - calc_flux(211,beam_energy*1.03,beam_current,T);
	double f_211_Th = calc_flux(211,beam_energy,beam_current,T+200.) - f_211;
	double f_211_Tl = f_211 - calc_flux(211,beam_energy,beam_current,T-200.);
	double f_211_cs = 0.2 * f_211;
	double b_211_err = b_AZ(211,"Fr","error");

	double sum_fb = f_208*b_AZ(208,"Fr") + f_209*b_AZ(209,"Fr") + f_210*b_AZ(210,"Fr") + f_211*b_AZ(211,"Fr");
	double f_208_errh = TMath::Sqrt(f_208_flux*f_208_flux + f_208_energyh*f_208_energyh + f_208_Th*f_208_Th + f_208_cs*f_208_cs);
	double f_208_errl = TMath::Sqrt(f_208_flux*f_208_flux + f_208_energyl*f_208_energyl + f_208_Tl*f_208_Tl + f_208_cs*f_208_cs);
	double fb_208_errh = TMath::Sqrt((b_AZ(208,"Fr")*f_208_errh)*(b_AZ(208,"Fr")*f_208_errh) + (f_208*b_208_err)*(f_208*b_208_err));
	double fb_208_errl = TMath::Sqrt((b_AZ(208,"Fr")*f_208_errl)*(b_AZ(208,"Fr")*f_208_errl) + (f_208*b_208_err)*(f_208*b_208_err));
	double f_209_errh = TMath::Sqrt(f_209_flux*f_209_flux + f_209_energyh*f_209_energyh + f_209_Th*f_209_Th + f_209_cs*f_209_cs);
	double f_209_errl = TMath::Sqrt(f_209_flux*f_209_flux + f_209_energyl*f_209_energyl + f_209_Tl*f_209_Tl + f_209_cs*f_209_cs);
	double fb_209_errh = TMath::Sqrt((b_AZ(209,"Fr")*f_209_errh)*(b_AZ(209,"Fr")*f_209_errh) + (f_209*b_209_err)*(f_209*b_209_err));
	double fb_209_errl = TMath::Sqrt((b_AZ(209,"Fr")*f_209_errl)*(b_AZ(209,"Fr")*f_209_errl) + (f_209*b_209_err)*(f_209*b_209_err));
	double f_210_errh = TMath::Sqrt(f_210_flux*f_210_flux + f_210_energyh*f_210_energyh + f_210_Th*f_210_Th + f_210_cs*f_210_cs);
	double f_210_errl = TMath::Sqrt(f_210_flux*f_210_flux + f_210_energyl*f_210_energyl + f_210_Tl*f_210_Tl + f_210_cs*f_210_cs);
	double fb_210_errh = TMath::Sqrt((b_AZ(210,"Fr")*f_210_errh)*(b_AZ(210,"Fr")*f_210_errh) + (f_210*b_210_err)*(f_210*b_210_err));
	double fb_210_errl = TMath::Sqrt((b_AZ(210,"Fr")*f_210_errl)*(b_AZ(210,"Fr")*f_210_errl) + (f_210*b_210_err)*(f_210*b_210_err));
	double f_211_errh = TMath::Sqrt(f_211_flux*f_211_flux + f_211_energyh*f_211_energyh + f_211_Th*f_211_Th + f_211_cs*f_211_cs);
	double f_211_errl = TMath::Sqrt(f_211_flux*f_211_flux + f_211_energyl*f_211_energyl + f_211_Tl*f_211_Tl + f_211_cs*f_211_cs);
	double fb_211_errh = TMath::Sqrt((b_AZ(211,"Fr")*f_211_errh)*(b_AZ(211,"Fr")*f_211_errh) + (f_211*b_211_err)*(f_211*b_211_err));
	double fb_211_errl = TMath::Sqrt((b_AZ(211,"Fr")*f_211_errl)*(b_AZ(211,"Fr")*f_211_errl) + (f_211*b_211_err)*(f_211*b_211_err));
	double sum_fb_errh = TMath::Sqrt(fb_208_errh*fb_208_errh + fb_209_errh*fb_209_errh + fb_210_errh*fb_210_errh + fb_211_errh*fb_211_errh);
	double sum_fb_errl = TMath::Sqrt(fb_208_errl*fb_208_errl + fb_209_errl*fb_209_errl + fb_210_errl*fb_210_errl + fb_211_errl*fb_211_errl);

	double r_210 = f_210/sum_fb;
	double r_210_errh = TMath::Sqrt((f_210_errh/sum_fb)*(f_210_errh/sum_fb) + (sum_fb_errl/sum_fb)*(sum_fb_errl/sum_fb)*r_210*r_210);
	double r_210_errl = TMath::Sqrt((f_210_errl/sum_fb)*(f_210_errl/sum_fb) + (sum_fb_errh/sum_fb)*(sum_fb_errh/sum_fb)*r_210*r_210);

	TCanvas *c2 = new TCanvas("c2");

	TLatex *l_210 = new TLatex();
	l_210->SetTextAlign(12);
	l_210->SetTextSize(0.02);
	l_210->DrawLatex(0.05,0.95,"Formulation/Condition:");
	l_210->DrawLatex(0.1,0.90,"f_{{}^{210}Fr} = #varepsilon_{transportation} #varepsilon_{desorption} (T) #varepsilon_{ionization} (T) #varepsilon_{escape} (E_{{}^{18}O}, T) P_{{}^{210}Fr} (E_{{}^{18}O}, I_{{}^{18}O})");
	l_210->DrawLatex(0.1,0.85,Form("E_{{}^{18}O} = %3.1f#rightarrow%3.1f (#pm3%% = +%3.1f, -%3.1f) MeV, I_{{}^{18}O} = %3.3f (#pm#frac{I}{#sqrt{j}} = %3.3f) e#muA, T_{Au} = %3.0f (#pm 200) #circC",beam_energy*18.,Be_degraded(beam_energy)*18.,Be_degraded(beam_energy*1.3)*18.-Be_degraded(beam_energy)*18.,Be_degraded(beam_energy)*18.-Be_degraded(beam_energy*0.97)*18.,beam_current*TMath::Power(10.,-3),beam_current*TMath::Power(10.,-3)/TMath::Sqrt(primaryFlux),T-273.));

	l_210->DrawLatex(0.05,0.80,"Model Calculation:");
	l_210->DrawLatex(0.1,0.75,Form("P_{{}^{210}Fr} = %3.1f pps, #DeltaP_{flux} = #pm%3.1f pps, #DeltaP_{energy} = (+%3.1f, -%3.1f) pps, #DeltaP_{#sigma} = #pm%3.1f pps",p_210,p_210_flux,p_210_energyh,p_210_energyl,p_210_cs));
	l_210->DrawLatex(0.1,0.70,Form("#varepsilon_{escape} = %3.3f, #Delta#varepsilon_{escape, energy} = (+%3.3f, -%3.3f), #Delta#varepsilon_{escape, T} = (+%3.3f, -%3.3f)",varesc,varesc_energyh,varesc_energyl,varesc_Th,varesc_Tl));
	l_210->DrawLatex(0.1,0.65,Form("#varepsilon_{ionization} = %3.3f, #Delta#varepsilon_{ionization, T} = (+%3.3f, -%3.3f)",varion,varion_Th,varion_Tl));
	l_210->DrawLatex(0.1,0.60,Form("f_{{}^{208}Fr} = %3.1f/s, #Deltaf_{flux} = #pm%3.1f/s, #Deltaf_{energy} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{T} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{#sigma} = #pm%3.1f/s, b_{{}^{208}Fr} = %.3f#pm%.3f",f_208,f_208_flux,f_208_energyh,f_208_energyl,f_208_Th,f_208_Tl,f_208_cs,b_AZ(208,"Fr"),b_208_err));
	l_210->DrawLatex(0.1,0.55,Form("f_{{}^{209}Fr} = %3.1f/s, #Deltaf_{flux} = #pm%3.1f/s, #Deltaf_{energy} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{T} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{#sigma} = #pm%3.1f/s, b_{{}^{209}Fr} = %.3f#pm%.3f",f_209,f_209_flux,f_209_energyh,f_209_energyl,f_209_Th,f_209_Tl,f_209_cs,b_AZ(209,"Fr"),b_209_err));
	l_210->DrawLatex(0.1,0.50,Form("f_{{}^{210}Fr} = %3.1f/s, #Deltaf_{flux} = #pm%3.1f/s, #Deltaf_{energy} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{T} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{#sigma} = #pm%3.1f/s, b_{{}^{210}Fr} = %.3f#pm%.3f",f_210,f_210_flux,f_210_energyh,f_210_energyl,f_210_Th,f_210_Tl,f_210_cs,b_AZ(210,"Fr"),b_210_err));
	l_210->DrawLatex(0.1,0.45,Form("f_{{}^{211}Fr} = %3.1f/s, #Deltaf_{flux} = #pm%3.1f/s, #Deltaf_{energy} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{T} = ^{+%3.1f}_{-%3.1f}/s, #Deltaf_{#sigma} = #pm%3.1f/s, b_{{}^{211}Fr} = %.3f#pm%.3f",f_211,f_211_flux,f_211_energyh,f_211_energyl,f_211_Th,f_211_Tl,f_211_cs,b_AZ(211,"Fr"),b_211_err));
	l_210->DrawLatex(0.1,0.40,Form("fb = #sum_{A = 208}^{211} f_{{}^{A}Fr} b_{{}^{A}Fr} = %3.1f ^{+%3.1f}_{-%3.1f} /s",sum_fb,sum_fb_errh,sum_fb_errl));
	l_210->DrawLatex(0.1,0.35,Form("R_{210} = #frac{f_{{}^{210}Fr}}{fb} = %3.3f ^{+%3.3f}_{-%3.3f}",r_210,r_210_errh,r_210_errl));

	// for 190330
	l_210->DrawLatex(0.05,0.3,"Experimental Value of Data 190330:");
	double fb_exp = 2048020.;
	double fb_stat = 31373.5;
	double fb_tsys = 44217.1;
	double fb_ssdsys = 1058927.1;
	double fb_dtsys = 18001.8;
	double fb_exp_sysh = TMath::Sqrt(fb_tsys*fb_tsys + fb_dtsys*fb_dtsys);
	double fb_exp_sysl = TMath::Sqrt(fb_ssdsys*fb_ssdsys);

	double f210_exp = fb_exp*r_210;
	double f210_stat = fb_stat*r_210;
	double f210_sysh = TMath::Sqrt((fb_exp_sysh*r_210)*(fb_exp_sysh*r_210) + (f210_exp*r_210_errh)*(f210_exp*r_210_errh));
	double f210_sysl = TMath::Sqrt((fb_exp_sysl*r_210)*(fb_exp_sysl*r_210) + (f210_exp*r_210_errl)*(f210_exp*r_210_errl));

	l_210->DrawLatex(0.1,0.25,Form("f_{{}^{210}Fr, exp} = %3.1f#pm%3.1f_{stat} pps, #Deltaf_{{}^{210}Fr, sys} = ^{+%3.1f}_{-%3.1f} pps",f210_exp,f210_stat,f210_sysh,f210_sysl));

	double eff_tdie = f210_exp/p_210;
	double eff_tdie_stat = f210_stat/p_210;
	double p_210_sysh = TMath::Sqrt(p_210_flux*p_210_flux + p_210_energyh*p_210_energyh + p_210_cs*p_210_cs);
	double p_210_sysl = TMath::Sqrt(p_210_flux*p_210_flux + p_210_energyl*p_210_energyl + p_210_cs*p_210_cs);
	double eff_tdie_sysh = TMath::Sqrt((f210_sysh/p_210)*(f210_sysh/p_210) + (p_210_sysl/p_210)*(p_210_sysl/p_210)*eff_tdie*eff_tdie);
	double eff_tdie_sysl = TMath::Sqrt((f210_sysl/p_210)*(f210_sysl/p_210) + (p_210_sysh/p_210)*(p_210_sysh/p_210)*eff_tdie*eff_tdie);
	l_210->DrawLatex(0.1,0.20,Form("#varepsilon_{transportation} #varepsilon_{desorption} (T) #varepsilon_{ionization} (T) #varepsilon_{escape} (E_{{}^{18}O}, T) = %3.3f #pm %3.3f (stat.) ^{+%3.3f}_{-%3.3f} (syst.)",eff_tdie,eff_tdie_stat,eff_tdie_sysh,eff_tdie_sysl));

	double eff_tdi = eff_tdie/varesc;
	double eff_tdi_stat = eff_tdie_stat/varesc;
	double varesc_sysh = TMath::Sqrt(varesc_energyh*varesc_energyh + varesc_Th*varesc_Th);
	double varesc_sysl = TMath::Sqrt(varesc_energyl*varesc_energyl + varesc_Tl*varesc_Tl);
	double eff_tdi_sysh = TMath::Sqrt((eff_tdie_sysh/varesc)*(eff_tdie_sysh/varesc) + (varesc_sysl/varesc)*(varesc_sysl/varesc)*eff_tdi*eff_tdi);
	double eff_tdi_sysl = TMath::Sqrt((eff_tdie_sysl/varesc)*(eff_tdie_sysl/varesc) + (varesc_sysh/varesc)*(varesc_sysh/varesc)*eff_tdi*eff_tdi);
	l_210->DrawLatex(0.1,0.15,Form("#varepsilon_{transportation} #varepsilon_{desorption} (T) #varepsilon_{ionization} = %3.3f #pm %3.3f (stat.) ^{+%3.3f}_{-%3.3f} (syst.)",eff_tdi,eff_tdi_stat,eff_tdi_sysh,eff_tdi_sysl));

	double eff_td = eff_tdi/varion;
	double eff_td_stat = eff_tdi_stat/varion;
	double eff_td_sysh = TMath::Sqrt((eff_tdi_sysh/varion)*(eff_tdi_sysh/varion) + (varion_Tl/varion)*(varion_Tl/varion)*eff_td*eff_td);
	double eff_td_sysl = TMath::Sqrt((eff_tdi_sysl/varion)*(eff_tdi_sysl/varion) + (varion_Th/varion)*(varion_Th/varion)*eff_td*eff_td);
	l_210->DrawLatex(0.1,0.1,Form("#varepsilon_{transportation} #varepsilon_{desorption} (T) = %3.3f #pm %3.3f (stat.) ^{+%3.3f}_{-%3.3f} (syst.)",eff_td,eff_td_stat,eff_td_sysh,eff_td_sysl));

	l_210->DrawLatex(0.1,0.05,Form("#varepsilon_{transportation} > %3.3f %%",(eff_td-eff_td_sysl)*100.));


	TCanvas *c_R = new TCanvas("c_R");

	// value
	l_210->DrawLatex(0.1,0.90,Form("R_{210} = %3.3f",r_210));

	// b_err
	double fb_berr = TMath::Sqrt((f_208*b_208_err)*(f_208*b_208_err) + (f_209*b_209_err)*(f_209*b_209_err) + (f_210*b_210_err)*(f_210*b_210_err) + (f_211*b_211_err)*(f_211*b_211_err));
	double r_berr = f_210 * fb_berr / (sum_fb*sum_fb);
	l_210->DrawLatex(0.15,0.85,Form("#Delta R_{210} (#Delta b_{210}) = #pm%3.3f", r_berr));

	// flux
	double fb_fluxerr = TMath::Sqrt((f_208_flux*b_208)*(f_208_flux*b_208) + (f_209_flux*b_209)*(f_209_flux*b_209) + (f_210_flux*b_210)*(f_210_flux*b_210) + (f_211_flux*b_211)*(f_211_flux*b_211));
	double r_fluxerr = TMath::Sqrt((f_210_flux/sum_fb)*(f_210_flux/sum_fb) + (fb_fluxerr/sum_fb)*(fb_fluxerr/sum_fb)*r_210*r_210);
	l_210->DrawLatex(0.15,0.80,Form("#Delta R_{210} (#Delta I_{{}^{18}O}) = #pm%3.3f", r_fluxerr));

	// energy
	double fb_energyerrh = TMath::Sqrt((f_208_energyh*b_208)*(f_208_energyh*b_208) + (f_209_energyh*b_209)*(f_209_energyh*b_209) + (f_210_energyh*b_210)*(f_210_energyh*b_210) + (f_211_energyh*b_211)*(f_211_energyh*b_211));
	double fb_energyerrl = TMath::Sqrt((f_208_energyl*b_208)*(f_208_energyl*b_208) + (f_209_energyl*b_209)*(f_209_energyl*b_209) + (f_210_energyl*b_210)*(f_210_energyl*b_210) + (f_211_energyl*b_211)*(f_211_energyl*b_211));
	double r_energyerrh = TMath::Sqrt((f_210_energyh/sum_fb)*(f_210_energyh/sum_fb) + (fb_energyerrl/sum_fb)*(fb_energyerrl/sum_fb)*r_210*r_210);
	double r_energyerrl = TMath::Sqrt((f_210_energyl/sum_fb)*(f_210_energyl/sum_fb) + (fb_energyerrh/sum_fb)*(fb_energyerrh/sum_fb)*r_210*r_210);
	l_210->DrawLatex(0.15,0.75,Form("#Delta R_{210} (#Delta E_{{}^{18}O}) = +%3.3f, -%3.3f", r_energyerrh, r_energyerrl));

	// temp
	double fb_temperrh = TMath::Sqrt((f_208_Th*b_208)*(f_208_Th*b_208) + (f_209_Th*b_209)*(f_209_Th*b_209) + (f_210_Th*b_210)*(f_210_Th*b_210) + (f_211_Th*b_211)*(f_211_Th*b_211));
	double fb_temperrl = TMath::Sqrt((f_208_Tl*b_208)*(f_208_Tl*b_208) + (f_209_Tl*b_209)*(f_209_Tl*b_209) + (f_210_Tl*b_210)*(f_210_Tl*b_210) + (f_211_Tl*b_211)*(f_211_Tl*b_211));
	double r_temperrh = TMath::Sqrt((f_210_Th/sum_fb)*(f_210_Th/sum_fb) + (fb_temperrl/sum_fb)*(fb_temperrl/sum_fb)*r_210*r_210);
	double r_temperrl = TMath::Sqrt((f_210_Tl/sum_fb)*(f_210_Tl/sum_fb) + (fb_temperrh/sum_fb)*(fb_temperrh/sum_fb)*r_210*r_210);
	l_210->DrawLatex(0.15,0.70,Form("#Delta R_{210} (#Delta T_{Au}) = +%3.3f, -%3.3f", r_temperrh, r_temperrl));

	// cs
	double fb_cserr = TMath::Sqrt((f_208_cs*b_208)*(f_208_cs*b_208) + (f_209_cs*b_209)*(f_209_cs*b_209) + (f_210_cs*b_210)*(f_210_cs*b_210) + (f_211_cs*b_211)*(f_211_cs*b_211));
	double r_cserr = TMath::Sqrt((f_210_cs/sum_fb)*(f_210_cs/sum_fb) + (fb_cserr/sum_fb)*(fb_cserr/sum_fb)*r_210*r_210);
	l_210->DrawLatex(0.15,0.65,Form("#Delta R_{210} (#Delta #sigma_{{}^{210}Fr}) = #pm%3.3f", r_cserr));


	// exp-data 190330
	double f210_tsys = fb_tsys*r_210;
	double f210_ssdsys = fb_ssdsys*r_210;
	double f210_dtsys = fb_dtsys*r_210;
	double f210_berr = fb_exp*r_berr;
	double f210_fluxerr = fb_exp*r_fluxerr;
	double f210_energyerrh = fb_exp*r_energyerrh;
	double f210_energyerrl = fb_exp*r_energyerrl;
	double f210_temperrh = fb_exp*r_temperrh;
	double f210_temperrl = fb_exp*r_temperrl;
	double f210_cserr = fb_exp*r_cserr;

	l_210->DrawLatex(0.1,0.55,Form("f_{210, exp} = %3.3f pps",f210_exp));
	l_210->DrawLatex(0.15,0.50,Form("#Delta f_{210, exp} = #pm%3.3f (stat) #pm%3.3f (t) -%3.3f (#varepsilon_{SSD}) +%3.3f (DT)",f210_stat,f210_tsys,f210_ssdsys,f210_dtsys));
	l_210->DrawLatex(0.20,0.45,Form("#pm%3.3f (b_{A}) #pm%3.3f (I) ^{+%3.3f}_{-%3.3f} (E) ^{+%3.3f}_{-%3.3f} (T) #pm%3.3f (#sigma)",f210_berr,f210_fluxerr,f210_energyerrh,f210_energyerrl,f210_temperrh,f210_temperrl,f210_cserr));


	double eff_tdie_tsys = f210_tsys/p_210;
	double eff_tdie_ssdsys = f210_ssdsys/p_210;
	double eff_tdie_dtsys = f210_dtsys/p_210;
	double eff_tdie_berr = f210_berr/p_210;
	double eff_tdie_fluxerr = TMath::Sqrt((f210_fluxerr/p_210)*(f210_fluxerr/p_210) + (p_210_flux/p_210)*(p_210_flux/p_210)*eff_tdie*eff_tdie);
	double eff_tdie_energyerrh = TMath::Sqrt((f210_energyerrh/p_210)*(f210_energyerrh/p_210) + (p_210_energyl/p_210)*(p_210_energyl/p_210)*eff_tdie*eff_tdie);
	double eff_tdie_energyerrl = TMath::Sqrt((f210_energyerrl/p_210)*(f210_energyerrl/p_210) + (p_210_energyh/p_210)*(p_210_energyh/p_210)*eff_tdie*eff_tdie);
	double eff_tdie_temperrh = f210_temperrh/p_210;
	double eff_tdie_temperrl = f210_temperrl/p_210;
	double eff_tdie_cserr = TMath::Sqrt((f210_cserr/p_210)*(f210_cserr/p_210) + 0.2*0.2*eff_tdie*eff_tdie);

	l_210->DrawLatex(0.1,0.40,Form("From P_{210} = %3.3f #pm %3.3f pps (#sigma),",p_210,p_210_cs));
	l_210->DrawLatex(0.1,0.35,Form("#frac{f_{210, exp}}{P_{210}} = %3.3f #pm%3.3f (stat) #pm%3.3f (t) -%3.3f (#varepsilon_{SSD}) +%3.3f (DT)",eff_tdie,eff_tdie_stat,eff_tdie_tsys,eff_tdie_ssdsys,eff_tdie_dtsys));
	l_210->DrawLatex(0.15,0.30,Form("#pm%3.3f (b_{A}) #pm%3.3f (I) ^{+%3.3f}_{-%3.3f} (E) ^{+%3.3f}_{-%3.3f} (T) #pm%3.3f (#sigma)",eff_tdie_berr,eff_tdie_fluxerr,eff_tdie_energyerrh,eff_tdie_energyerrl,eff_tdie_temperrh,eff_tdie_temperrl,eff_tdie_cserr));

	l_210->DrawLatex(0.1,0.20,Form("#varepsilon_{ionization} = 1, assume #varepsilon_{desorption} = 1, #varepsilon_{escape} = %3.3f ^{+%3.3f}_{-%3.3f} (E) ^{+%3.3f}_{-%3.3f} (T)",varesc,varesc_energyh,varesc_energyl,varesc_Th,varesc_Tl));

	double eff_t = eff_tdie/varesc;
	double eff_t_stat = eff_tdie_stat/varesc;
	double eff_t_tsys = eff_tdie_tsys/varesc;
	double eff_t_ssdsys = eff_tdie_ssdsys/varesc;
	double eff_t_dtsys = eff_tdie_dtsys/varesc;
	double eff_t_berr = eff_tdie_berr/varesc;
	double eff_t_fluxerr = eff_tdie_fluxerr/varesc;
	double eff_t_energyerrh = TMath::Sqrt((eff_tdie_energyerrh/varesc)*(eff_tdie_energyerrh/varesc) + (varesc_energyl/varesc)*(varesc_energyl/varesc)*eff_t*eff_t);
	double eff_t_energyerrl = TMath::Sqrt((eff_tdie_energyerrl/varesc)*(eff_tdie_energyerrl/varesc) + (varesc_energyh/varesc)*(varesc_energyh/varesc)*eff_t*eff_t);
	double eff_t_temperrh = TMath::Sqrt((eff_tdie_temperrh/varesc)*(eff_tdie_temperrh/varesc) + (varesc_Tl/varesc)*(varesc_Tl/varesc)*eff_t*eff_t);
	double eff_t_temperrl = TMath::Sqrt((eff_tdie_temperrl/varesc)*(eff_tdie_temperrl/varesc) + (varesc_Th/varesc)*(varesc_Th/varesc)*eff_t*eff_t);
	double eff_t_cserr = eff_tdie_cserr/varesc;

	l_210->DrawLatex(0.1,0.15,Form("#varepsilon_{transportation} = %3.3f #pm%3.3f (stat) #pm%3.3f (t) -%3.3f (#varepsilon_{SSD}) +%3.3f (DT)",eff_t,eff_t_stat,eff_t_tsys,eff_t_ssdsys,eff_t_dtsys));
	l_210->DrawLatex(0.15,0.10,Form("#pm%3.3f (b_{A}) #pm%3.3f (I) ^{+%3.3f}_{-%3.3f} (E) ^{+%3.3f}_{-%3.3f} (T) #pm%3.3f (#sigma)",eff_t_berr,eff_t_fluxerr,eff_t_energyerrh,eff_t_energyerrl,eff_t_temperrh,eff_t_temperrl,eff_t_cserr));

  rootapp.Run();

  return 0;
}
