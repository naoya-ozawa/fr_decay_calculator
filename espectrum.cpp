#include <iostream>
#include <TCanvas.h>
#include <TPie.h>
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

// Calculate the primary beam energy (MeV/u) after transmitting through the Be window
double Be_degraded(double beam_energy, double T_u, double T_d, double T_He, bool plot=true){
	// Fitted parameters from
	// https://github.com/naoya-ozawa/Be_window_calculations/blob/master/SRIM_calculations/dependence_plotter.cpp
	double C_E0 = -9.; // MeV
	double C_E1 = 800.; // MeV*MeV_incident
	double C_Be = 0.68; // MeV/um_Be
	double C_He = 0.09; // MeV/mm_He
	double T_Be = T_u+T_d;

	if (plot == true){
		TCanvas *c_bedegraded = new TCanvas("c_bedegraded");
		TF1 *f_bedegraded = new TF1("E_{loss}","[0] + [1]/(18.*x)",6.0,8.0);
		f_bedegraded->SetParameter(0,C_E0+C_Be*T_Be+C_He*T_He);
		f_bedegraded->SetParameter(1,C_E1);
		TGraph *g_bedegraded = new TGraph(f_bedegraded);
		g_bedegraded->SetTitle("Energy loss at Be window;{}^{18}O beam energy (MeV/u);Total energy loss (MeV)");
		g_bedegraded->Draw("APL");
		c_bedegraded->cd()->SetGridx();
		c_bedegraded->cd()->SetGridy();
		TLatex *l_bedegraded = new TLatex();
		l_bedegraded->DrawLatex(7.5,TMath::MaxElement(g_bedegraded->GetN(),g_bedegraded->GetY()),Form("T_{Be:Up} = %3.1f #mum",T_u));
		l_bedegraded->DrawLatex(7.5,TMath::MaxElement(g_bedegraded->GetN(),g_bedegraded->GetY())-0.2,Form("T_{He} = %3.1f mm",T_He));
		l_bedegraded->DrawLatex(7.5,TMath::MaxElement(g_bedegraded->GetN(),g_bedegraded->GetY())-0.4,Form("T_{Be:Down} = %3.1f #mum",T_d));
	}

	double E_loss = C_E0 + (C_E1/(beam_energy*18.)) + C_Be*T_Be + C_He*T_He; // MeV
//	cout << "Degraded energy: " << E*18.-E_loss << " MeV" << endl;
	return beam_energy - (E_loss/18.);
}

// Estimate the normalized production based on Stancari2006
double normprod(double incident_energy, int A, bool plot=true){
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
	if (plot == true){
		TCanvas *c_normprod = new TCanvas("c_normprod");
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
	}
	double E = incident_energy*18.; // MeV/u --> MeV, incident energy
	if (plot == true){
		TLatex *poverj = new TLatex();
		TLine *l_208 = new TLine(81.,s208->Eval(E),119.,s208->Eval(E));
		l_208->SetLineColor(3);
		l_208->SetLineStyle(2);
		l_208->SetLineWidth(2);
		l_208->Draw();
		poverj->DrawLatex(82.,s208->Eval(E),Form("%2.1e",s208->Eval(E)));
		TLine *l_209 = new TLine(81.,s209->Eval(E),119.,s209->Eval(E));
		l_209->SetLineColor(4);
		l_209->SetLineStyle(2);
		l_209->SetLineWidth(2);
		l_209->Draw();
		poverj->DrawLatex(82.,s209->Eval(E),Form("%2.1e",s209->Eval(E)));
		TLine *l_210 = new TLine(81.,s210->Eval(E),119.,s210->Eval(E));
		l_210->SetLineColor(2);
		l_210->SetLineStyle(2);
		l_210->SetLineWidth(2);
		l_210->Draw();
		poverj->DrawLatex(82.,s210->Eval(E),Form("%2.1e",s210->Eval(E)));
		TLine *l_211 = new TLine(81.,s211->Eval(E),119.,s211->Eval(E));
		l_211->SetLineColor(5);
		l_211->SetLineStyle(2);
		l_211->SetLineWidth(2);
		l_211->Draw();
		poverj->DrawLatex(82.,s211->Eval(E),Form("%2.1e",s211->Eval(E)));
		TLine *l_eval = new TLine(E,1.*TMath::Power(10.,-7),E,4.9*TMath::Power(10.,-6));
		l_eval->SetLineColor(2);
		l_eval->SetLineStyle(2);
		l_eval->SetLineWidth(5);
		l_eval->Draw();
		poverj->DrawLatex(E+1.,2.5e-6,Form("%3.1f MeV",E));
	}

	if (A == 208){
		return s208->Eval(E);
	}else if (A==209){
		return s209->Eval(E);
	}else if (A==210){
		return s210->Eval(E);
	}else if (A==211){
		return s211->Eval(E);
	}else{
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

double depth_Fr(double incident_energy, int isotope, bool plot=true){
	// energy: MeV/u (injection energy into Au target)
	// use the derivative of normalized production

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


	// obtain energy dependence of dEdX
	if (plot == true){
		TCanvas *c_Fr_dedx = new TCanvas("c_Fr_dedx");
	}

	double energy_dedx[13] = {50.,55.,60.,65.,70.,80.,90.,100.,110.,120.,130.,140.,150.}; // MeV
	double dedx[13] = {2.03,1.97,1.90,1.84,1.79,1.69,1.60,1.52,1.46,1.39,1.34,1.29,1.24}; // MeV/(mg/cm2)
	TGraph *g_dedx = new TGraph(13,energy_dedx,dedx);
	g_dedx->SetTitle("Energy loss <dE/dX> of {}^{18}O in Au;Beam Energy E (MeV);<dE/dX> (MeV/(mg/cm^{2}))");
	TF1 *f_dedx = new TF1("f_dedx","[0]+[1]*TMath::Power(x,[2])",50.,150.);
	f_dedx->SetParameters(1.,1.,-0.3);
	if (plot == true){
		g_dedx->Fit("f_dedx");
	}else{
		g_dedx->Fit("f_dedx","Q");
	}
	if (plot == true){
		g_dedx->Draw("AP*");
		f_dedx->Draw("SAME");
	}

	// plot cross section

	TGraph *g208 = new TGraph(10,energy,fr208);
	TSpline3 *s208 = new TSpline3("Spline Fit for 208",g208);
	TGraph *g209 = new TGraph(10,energy,fr209);
	TSpline3 *s209 = new TSpline3("Spline Fit for 209",g209);
	TGraph *g210 = new TGraph(10,energy,fr210);
	TSpline3 *s210 = new TSpline3("Spline Fit for 210",g210);
	TGraph *g211 = new TGraph(10,energy,fr211);
	TSpline3 *s211 = new TSpline3("Spline Fit for 211",g211);

	TGraph *dg208 = new TGraph();
	TGraph *dg209 = new TGraph();
	TGraph *dg210 = new TGraph();
	TGraph *dg211 = new TGraph();
	for (int i=0; i<60; ++i){
		double e = 85. + double(i)*0.5; // MeV
		double mass = 196.96655 * 1000.; // mg/mol
		double avogadro = 6.0 * TMath::Power(10.,23); // /mol
		double transcoef = TMath::Power(10.,27); // mb/cm2
		dg208->SetPoint(i,e,transcoef*f_dedx->Eval(e)*mass*s208->Derivative(e)/avogadro);
		dg209->SetPoint(i,e,transcoef*f_dedx->Eval(e)*mass*s209->Derivative(e)/avogadro);
		dg210->SetPoint(i,e,transcoef*f_dedx->Eval(e)*mass*s210->Derivative(e)/avogadro);
		dg211->SetPoint(i,e,transcoef*f_dedx->Eval(e)*mass*s211->Derivative(e)/avogadro);
	}

	TSpline3 *ds208 = new TSpline3("Spline fit for 208",dg208);
	dg208->SetTitle("{}^{208}Fr");
	dg208->SetLineColor(3);
	dg208->SetLineWidth(2);
	ds208->SetLineColor(3);
	ds208->SetLineWidth(2);
	TSpline3 *ds209 = new TSpline3("Spline fit for 209",dg209);
	dg209->SetTitle("{}^{209}Fr");
	dg209->SetLineColor(4);
	dg209->SetLineWidth(2);
	ds209->SetLineColor(4);
	ds209->SetLineWidth(2);
	TSpline3 *ds210 = new TSpline3("Spline fit for 210",dg210);
	dg210->SetTitle("{}^{210}Fr");
	dg210->SetLineColor(2);
	dg210->SetLineWidth(2);
	ds210->SetLineColor(2);
	ds210->SetLineWidth(2);
	TSpline3 *ds211 = new TSpline3("Spline fit for 211",dg211);
	dg211->SetTitle("{}^{211}Fr");
	dg211->SetLineColor(5);
	dg211->SetLineWidth(2);
	ds211->SetLineColor(5);
	ds211->SetLineWidth(2);
	if (plot == true){
		TCanvas *c_depth_Fr = new TCanvas("c_depth_Fr");
		TMultiGraph *mg_depth_Fr = new TMultiGraph();
		mg_depth_Fr->SetTitle("Evaporation Cross Section Based on [Stancari2006];Incident {}^{18}O beam energy (MeV);Fusion-evaporation cross section #sigma_{EvR} (mb)");
		mg_depth_Fr->Add(dg208);
		mg_depth_Fr->Add(dg209);
		mg_depth_Fr->Add(dg210);
		mg_depth_Fr->Add(dg211);
		mg_depth_Fr->Draw("AP");
		c_depth_Fr->BuildLegend();
		ds208->Draw("SAME");
		ds209->Draw("SAME");
		ds210->Draw("SAME");
		ds211->Draw("SAME");
	}

	// find peak of cross sections
	double max208 = 0., max209 = 0., max210 = 0., max211 = 0.;
	double peak208 = 0., peak209 = 0., peak210 = 0., peak211 = 0.;
	for (int i=0; i<60; ++i){
		double e = 85. + double(i)*0.5; // MeV
		if (ds208->Eval(e) > max208){
			max208 = ds208->Eval(e);
			peak208 = e;
		}
		if (ds209->Eval(e) > max209){
			max209 = ds209->Eval(e);
			peak209 = e;
		}
		if (ds210->Eval(e) > max210){
			max210 = ds210->Eval(e);
			peak210 = e;
		}
		if (ds211->Eval(e) > max211){
			max211 = ds211->Eval(e);
			peak211 = e;
		}
	}
	if (plot == true){
		TLine *l_peak208 = new TLine(peak208,1.,peak208,200.);
		l_peak208->SetLineColor(3);
		l_peak208->SetLineWidth(3);
		l_peak208->SetLineStyle(2);
		l_peak208->Draw();
		TLine *l_peak209 = new TLine(peak209,1.,peak209,200.);
		l_peak209->SetLineColor(4);
		l_peak209->SetLineWidth(3);
		l_peak209->SetLineStyle(2);
		l_peak209->Draw();
		TLine *l_peak210 = new TLine(peak210,1.,peak210,200.);
		l_peak210->SetLineColor(2);
		l_peak210->SetLineWidth(3);
		l_peak210->SetLineStyle(2);
		l_peak210->Draw();
		TLine *l_peak211 = new TLine(peak211,1.,peak211,200.);
		l_peak211->SetLineColor(5);
		l_peak211->SetLineWidth(3);
		l_peak211->SetLineStyle(2);
		l_peak211->Draw();
	}

	if (isotope == 208){
		double peak = peak208; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 209){
		double peak = peak209; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 210){
		double peak = peak210; // MeV based on stancari2006 plot
		return range(peak,incident_energy)*TMath::Power(10.,-4); // cm
	}else if (isotope == 211){
		double peak = peak211; // MeV based on stancari2006 plot
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

// desorption efficiency calculation
double des_eff_Fr(double tau0, double E_ads, double T){

	double kB = 8.6*TMath::Power(10.,-5); // eV/K
	double res_time = tau0 * TMath::Exp(E_ads/(kB*T)); // sec

	// escape probability before decay, following Fujioka1981
	return tau_AZ(210,"Fr")/(tau_AZ(210,"Fr")+res_time);
}

// ε_{extraction} = ε_{transportation} ε_{desorption} ε_{ionization} ε_{escape}
double extraction_eff(double incident_energy, double temp, int A, bool plot=true){

	const double E_wf_Mo = 4.6; // eV
	const double E_wf_Au = 5.1; // eV
	const double E_ip_Fr = 4.07; // eV
	const double E_ip_Rn = 10.74850; // eV
	const double E_ip_At = 9.31751; // eV
	const double E_ip_Po = 8.418070; // eV
	const double sw_Fr = 0.5; // statistical weight factor for Fr
	const double tau_0_FrAu = 1.9*TMath::Power(10.,-13); // s; Cs-Re from Delhuille2002
	const double E_ads_FrAu = 2.0; // eV; Cs-Re from Delhuille2002

	double des_eff = des_eff_Fr(tau_0_FrAu,E_ads_FrAu,temp); // desorption efficiency from the Au surface before decay
	double trans_eff = 1.0; // transportation efficiency from the target to the MCP surface (assumption)

	double escape = y_Fr(tau_AZ(A,"Fr"),D_m(temp,A),depth_Fr(incident_energy,A,plot))*ionization(temp,E_wf_Au,E_ip_Fr,sw_Fr)*des_eff*trans_eff;

	return escape;
}

double escape_time(double incident_energy, double tau0, double E_ads, double T, int A, bool plot=true){
	// based on Fujioka1981 & Delhuille2002
	double kB = 8.6*TMath::Power(10.,-5); // eV/K

	if (plot == true){
		// plot temperature-dependence and energy-dependence
		TCanvas *c_escape_time = new TCanvas("c_escape_time");
		c_escape_time->Divide(1,2);

		c_escape_time->cd(1);
		TGraph *g_et_tempdep = new TGraph();
		g_et_tempdep->SetTitle(Form("Temperature dependence of escape time at E_{{}^{18}O} = %3.1f MeV;Au temperature (K);#tau = #tau_{esc}+#tau_{des} (s)",incident_energy*18.));
		for (int i=0; i<200; ++i){
			double temp = 800. + 273. + double(i);
			double time = depth_Fr(incident_energy,A,false)*depth_Fr(incident_energy,A,false)/(3.*D_m(temp,A)) + tau0 * TMath::Exp(E_ads/(kB*temp));
			g_et_tempdep->SetPoint(i,temp,time);
		}
		g_et_tempdep->Draw("ALP");


		c_escape_time->cd(2);
		TGraph *g_et_enerdep = new TGraph();
		g_et_enerdep->SetTitle(Form("Energy dependence of escape time at T_{Au} = %3.1f degC;{}^{18}O energy (MeV);#tau = #tau_{esc}+#tau_{des} (s)",T-273));
		for (int i=0; i<100; ++i){
			double ener = 6. + double(i)/50.;
			double time = depth_Fr(ener,A,false)*depth_Fr(ener,A,false)/(3.*D_m(T,A)) + tau0 * TMath::Exp(E_ads/(kB*T));
			g_et_enerdep->SetPoint(i,ener,time);
		}
		g_et_enerdep->Draw("ALP");
	}

	double t_ave = depth_Fr(incident_energy,A,false)*depth_Fr(incident_energy,A,false)/(3.*D_m(T,A));
	double tau_D = tau0 * TMath::Exp(E_ads/(kB*T));
	if (depth_Fr(incident_energy,A,false) < 0.0){
		return 0.0;
	}else{
		return t_ave+tau_D;
	}
}

double eua_to_pps(double eua){
	double charge = 6.0;
	double e = 1.6 * TMath::Power(10.,-19); // C
	double u_to_0 = TMath::Power(10.,-6); // A/uA
	return eua * u_to_0 / (charge * e);
}

// calculate the Fr flux (steady-state) in pps
// f = ε_{extraction} * normprod * 18-O current
double Fr_flux(int isotope, double energy_18O, double current_18O, double Au_temperature, double T_Be_upstream, double T_Be_downstream, double T_Helium, bool plot=true){

	const double tau_0_FrAu = 1.9*TMath::Power(10.,-13); // s; Cs-Re from Delhuille2002
	const double E_ads_FrAu = 2.0; // eV; Cs-Re from Delhuille2002

	double E = Be_degraded(energy_18O,T_Be_upstream,T_Be_downstream,T_Helium,plot);
	double ee = extraction_eff(E,Au_temperature,isotope,plot);
	double pa = normprod(E,isotope,plot);
	double J = eua_to_pps(current_18O);

	cout << "==================" << endl;
	cout << "For " << isotope << "-Fr:" << endl;
	cout << "ε = " << ee*100. << "%" << endl;
	cout << "P = " << pa*J << " pps" << endl;
	cout << "Escape time = " << escape_time(E,tau_0_FrAu,E_ads_FrAu,Au_temperature,isotope,plot) << " s" << endl;
	cout << "Flux = " << ee*pa*J << " pps" << endl;
	return ee * pa * J;
}

void prod_chart(double *values, int *colors, const char **lbls, double E_0, double J, double T_Au){
	TCanvas *c_prod_chart = new TCanvas("c_prod_chart");
	const char* ttl = Form("#splitline{Fr isotopes extracted from Au target}{#left(E_{0} = %3.1f MeV/u, J = %3.1f e#muA, T_{Au} = %3.1f #circC#right)}",E_0,J,T_Au);
	TPie *pie_prod = new TPie("pie_prod",ttl,4,values,colors,lbls);
	pie_prod->SetRadius(.2);
	pie_prod->SetLabelsOffset(.01);
	pie_prod->SetLabelFormat("#splitline{%val (%perc)}{%txt}");
	pie_prod->SetValueFormat("%.1e");
	pie_prod->SetAngle3D(0);
	pie_prod->SetAngularOffset(90);
	pie_prod->Draw("nol");
}

int main (int argc, char** argv){

	TRint rootapp("app",&argc,argv);

	const double energy_18O = 7.0; // MeV/u
	const double current_18O = 6.0; // euA
	const double Au_temperature = 700.+273.; // K

//	// NP1712-AVF52-02
//	const double T_Be_upstream = 9.75; // um // MEF#101435899 #1
//	const double T_Be_downstream = 9.70; // um // MEF#101435899 #3
//	const double T_Helium = 4.2+17.5+6.7; // mm // for NP1217-AVF52-03

	// // NP1712-AVF52-03
	// const double T_Be_upstream = 10.0; // um // MEF#101435899 #1
	// const double T_Be_downstream = 10.2; // um // MEF#101435899 #3
	// const double T_Helium = 7.0; // mm // for NP1217-AVF52-03

	// NP2012-AVF72-01
	const double T_Be_upstream = 11.2; // um
	const double T_Be_downstream = 12.4; // um
	const double T_Helium = 7.0; // mm

	cout << "T_Au = " << Au_temperature-273. << " degC" << endl;
	cout << "E_0 = " << energy_18O << " MeV/u --> Be Window --> " << Be_degraded(energy_18O,T_Be_upstream,T_Be_downstream,T_Helium,false)*18. << " MeV" << endl;
	cout << "J_18-O = " << eua_to_pps(current_18O) << " pps (" << current_18O << " euA)" << endl;

	double values[4];
	values[0] = Fr_flux(208,energy_18O,current_18O,Au_temperature,T_Be_upstream,T_Be_downstream,T_Helium,false);
	values[1] = Fr_flux(209,energy_18O,current_18O,Au_temperature,T_Be_upstream,T_Be_downstream,T_Helium,false);
	values[2] = Fr_flux(210,energy_18O,current_18O,Au_temperature,T_Be_upstream,T_Be_downstream,T_Helium,true);
	values[3] = Fr_flux(211,energy_18O,current_18O,Au_temperature,T_Be_upstream,T_Be_downstream,T_Helium,false);
	int colors[4] = {3,4,2,5};
	const char *lbls[4] = {"{}^{208}Fr","{}^{209}Fr","{}^{210}Fr","{}^{211}Fr"};
	prod_chart(values,colors,lbls,energy_18O,current_18O,Au_temperature-273.);

	rootapp.Run();

	return 0;
}
