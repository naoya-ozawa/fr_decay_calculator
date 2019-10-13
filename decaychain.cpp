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

int main (int argc, char** argv){

	double seconds = 1.0;
	double minutes = 60.0;
	double hours = 60.0 * minutes;
	double days = 24.0 * hours;
	double years = 365.0 * days;

	double R_211fr = 0.19;

	double br_211fr = 0.80;
	double brerr_211fr = 0.20;
	double t_211fr = 3.10*minutes;
	double terr_211fr = 0.02*minutes;

	double br_211rn = 0.274;
	double brerr_211rn = 0.017;
	double t_211rn = 14.6*hours;
	double terr_211rn = 0.2*hours;

	double br_211at = 0.4180;
	double brerr_211at = 0.0008;
	double t_211at = 7.214*hours;
	double terr_211at = 0.007*hours;

	double br_211po = 1.0;
	double brerr_211po = 0.0;
	double t_211po = 0.516*seconds;
	double terr_211po = 0.003*seconds;

	double br_207at = 0.086;
	double brerr_207at = 0.01;
	double t_207at = 1.80*hours;
	double terr_207at = 0.04*hours;

	double br_207po = 0.00021;
	double brerr_207po = 0.00002;
	double t_207po = 5.80*hours;
	double terr_207po = 0.02*hours;

	double br_207bi = 0.0;
	double brerr_207bi = 0.0;
	double t_207bi = 31.55*years;
	double terr_207bi = 0.05*years;

	double br_207pb = 0.0;
	double brerr_207pb = 0.0;
	double t_207pb = 10000.*years;
	double terr_207pb = 0.0;
