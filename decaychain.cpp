double parent(double flux,double lifetime,double time){
	return flux * lifetime * ( 1.0 - TMath::Exp(-time/lifetime));
}

double daughter(double flux,double lifetime,double alphaParentFlux,double alphaParentLifetime,double alphaParantBR,double betaParentFlux,double betaParentLifetime,double betaParentBR,double time){
	double coef_daughter = flux * lifetime;
	double coef_alphaParent = (alphaParentBR*alphaParentFlux*lifetime*alphaParentLifetime) / (alphaParentLifetime - lifetime);
	double coef_betaParent = ((1.0-betaParentBR)*betaParentFlux*lifetime*betaParentLifetime) / (betaParentLifetime - lifetime);
	
	double exp_daughter = 1.0 - TMath::Exp(-time/lifetime);
	double exp_alphaParent = TMath::Exp(-time/alphaParentLifetime) - TMath::Exp(-time/lifetime);
	double exp_betaParent = TMath::Exp(-time/betaParentLifetime) - TMath::Exp(-time/lifetime);

	return coef_daughter*exp_daughter + coef_alphaParent*exp_alphaParent + coef_betaParent*exp_betaParent;
}
