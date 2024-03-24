//==================================================================================================
//
// Computation of the thermal and k-SZ effect using explicit integration collision term.
// Using a precomputed kernel over the angular directions. 
// TODO: Update this!!!!
// The cluster is assumed to be isothermal and moving at a given speed betac with direction muc 
// relative to the line of sight. The integrals are carried out in the cluster frame.
//
//==================================================================================================
//
// Author: Elizabeth Lee
// Based on work by Jens Chluba (CITA, University of Toronto)
//
// first implementation: October 2018
// last modification: March 2021
//
//==================================================================================================

#include "SZ_multipoleKernel.h"
#include <cmath>
#include "global_functions.h"


//==================================================================================================
// Multipole Kernel
//==================================================================================================
MultipoleKernel::MultipoleKernel()
    : MultipoleKernel(0, 0.0, 1.0, 1.0e-4) {}

MultipoleKernel::MultipoleKernel(int l_i, double s_i, double eta_i)
    : MultipoleKernel(l_i, s_i, eta_i, 1.0e-4) {}

MultipoleKernel::MultipoleKernel(int l_i, double s_i, double eta_i, double Int_eps_i){
    Update_l(l_i);
    Update_s(s_i);
    eta = eta_i;
    Int_eps = Int_eps_i;
    electron_anisotropy = false;
    gamma0 = sqrt(1+eta*eta);
    s_low = log((gamma0-eta)/(gamma0+eta));
    s_high = log((gamma0+eta)/(gamma0-eta));
    Lpart = mus = mu = r = 0;
    P.resize(l+3);
    h0.resize(l+1);
    h2.resize(l+1);
    h4.resize(l+1);
    K1.resize(l+1);
    K2.resize(l+1);
}

vector<double> MultipoleKernel::s_limits(){
    vector<double> lims;
    lims.resize(2);
    lims[0] = s_low;
    lims[1] = s_high;
    return lims;
}

void MultipoleKernel::Update_s(double s_i){
    s = s_i;
    t = exp(s);
}

void MultipoleKernel::Update_l(int l_i){
    l = l_i;
    if (l < 0){
        l = 0;
        print_error("l is out of bounds, it must be >=0. l is now set to 0.");
    }
}

void MultipoleKernel::Calculate_integral_variables(){
    //Standard definitions
    //==============================================================================================
    // mu == cosine of angle of gamma  and beta 
    // mup== cosine of angle of gamma' and beta
    // mus== cosine of angle of gamma  and gamma'
    //==============================================================================================   
    double beta0 = eta/gamma0; 
    double zeta=t/(pow(gamma0-eta*mu,2));
    double alpha_sc=1.0-mus;

    //d2sigma/dmu/dmup as defined equation 2 CNSN
    double dsig = 3.0/8.0/PI/gamma0*gamma0*(1.0-zeta*alpha_sc*(1.0-0.5*zeta*alpha_sc));
    double dphidt = t/sqrt(beta0*beta0*t*t*(-1+mu*mu)*(-1+mus*mus)-pow((t-1)+beta0*mu*(1-t*mus),2));
    double Plmu = 0;
    if (!electron_anisotropy) {
        for (int k=0; k<l+1; k++){
            Plmu += pow(-1,k)*Binomial(l,k)*Binomial(l+k,k)*pow((1-mus)/2,k);
        }
    }
    else {
        for (int k=0; k<l+1; k++){
            Plmu += pow(-1,k)*Binomial(l,k)*Binomial(l+k,k)*pow((1-mu)/2,k);
        }
    }
    r = Plmu*dphidt*dsig;
}

double MultipoleKernel::sigl_Boltzmann_Compton(double int_mu){
    mu = int_mu;
    Calculate_integral_variables();
    return r;
}

double MultipoleKernel::mu_Int(double int_mus){
    mus = int_mus;
    double mus_cr = (2*eta*eta*t-(1-t)*(1-t))/(2*eta*eta*t);
    double limitvar = eta*t*sqrt(2*t*(1-mus*mus)*(mus_cr-mus));
    double a=((1 - t)*gamma0*(1- t*mus)-limitvar)/(eta*(1+t*t-2*t*mus));
    double b=((1 - t)*gamma0*(1- t*mus)+limitvar)/(eta*(1+t*t-2*t*mus));
    double epsrel=Int_eps*0.8, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sigl_Boltzmann_Compton(int_var);});
}

double MultipoleKernel::mus_Int(){
    double mus_cr = (2*eta*eta*t-(1-t)*(1-t))/(2*eta*eta*t);
    double a=-1.0;
    double b=mus_cr;
    double epsrel=Int_eps*0.9, epsabs=1.0e-300;

    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mu_Int(int_var);});
}

void MultipoleKernel::Calculate_formula_variables(){
    Lpart = 0.5*(fabs(s)-2.0*asinh(eta));
    for (int m = 0; m < l+3; m++){
        double ppart = pow(eta,2*m+1)/pow(1+eta*eta,2.5);
        double tpart = fabs(pow(t-1,2*m+1))/pow(4*t,m-2)/pow(1+t,5.0);
        P[m] = ppart-tpart;
    }

    h0[0] = -(2*P[2]+5*P[1])/15.0;
    h2[0] = -(P[2]+P[1])/3.0;
    h4[0] = Lpart+P[2]+2*P[1]+P[0];

    if (l >= 1){
        h0[1] = P[2]/5.0;
        h2[1] = -(3*Lpart+4*P[2]+7*P[1]+3*P[0])/3.0;
        h4[1] = (3*Lpart+P[3]+5*P[2]+7*P[1]+3*P[0])/2.0;

        if (l >= 2){
            h0[2] = Lpart+(23*P[2]+35*P[1]+15*P[0])/15.0;
            h2[2] = (-5*h0[2]-P[3])/2.0;
            h4[2] = (-3*h2[2]-P[4]-P[3])/4.0;
            
            if (l >= 3){
                for (int m = 3; m<l+1; m++){
                    h0[m] = ((2*m+1)*h0[m-1]-pow(-1,m)*P[m])/(2*(m-2));
                    h2[m] = (-5*h0[m]-pow(-1,m)*P[m+1])/(2*(m-1));
                    h4[m] = (15*h0[m]-pow(-1,m)*(2*m-5)*P[m+1]-pow(-1,m)*2*(m-1)*P[m+2])/(4*m*(m-1));
                }
            }
        }
    }

    for (int m = 0; m<l+1; m++){
        //Calculating K1
        double invpref = eta*gamma0*pow(4*t,m);
        double ksum = 0;
        for (int k = 0; k < m+1; k++){
            double pref = Binomial(m,k)*pow(-1,k)*pow(t-1,2*(m-k))/(2*k+1);
            ksum += pref*(pow(t+1,2*k+1)-(pow(gamma0*fabs(t-1)/eta,2*k+1)));
        }
        K1[m] = ksum/invpref;

        //Calculating K2
        K2[m] = 5*(1+eta*eta)*(1+t)*(1+t)*h0[m];
        K2[m] += -((5+2*eta*eta)*(1+t*t)+(22+16*eta*eta)*t)*h2[m];
        K2[m] += 4*(3+2*eta*eta)*t*h4[m];
    }
}

double MultipoleKernel::Calculate_electron_multipoles(){
    if (s<s_low|| s>s_high) { return 0.0; }

    Lpart = 0.5*(fabs(s)-2.0*asinh(eta));
    double beta0 = eta/gamma0;
    vector<double> K(l+5,0.0);
    for (int k = 0; k< l+5; k++){
        int m = k-3;
        double part1 = (pow(t,m)+1)*(pow(1+beta0,m)-pow(1-beta0,m));
        double part2 = (pow(t,m)-1)*(pow(1+beta0,m)+pow(1-beta0,m));
        double signt = t<1 ? -1.0 : 1.0;
        K[k] = (part1-signt*part2)/(2*pow(1-beta0,m));
    }

    vector<double> X0(l+1,0.0), X1(l+1,0.0), X2(l+1,0.0), X3(l+1,0.0), X4(l+1,0.0);
    for (int m = 0; m < l+1; m++){
        for (int k = 0; k < m+1; k++){
            double factor = Binomial(m,k)*pow(-1,k);
            X0[m] += factor*K[k+4]/(k+1.0);
            X1[m] += (k==0) ? -2.0*factor*Lpart : factor*K[k+3]/k;
            X2[m] += (k==1) ? -2.0*factor*Lpart : factor*K[k+2]/(k-1.0);
            X3[m] += (k==2) ? -2.0*factor*Lpart : factor*K[k+1]/(k-2.0);
            X4[m] += (k==3) ? -2.0*factor*Lpart : factor*K[k]/(k-3.0);
        }
        double pref = pow((gamma0-eta)/eta,m+1)/pow(2.0,m);
        X0[m] *= pref;
        X1[m] *= pref/(gamma0-eta);
        X2[m] *= pref/pow(gamma0-eta,2);
        X3[m] *= pref/pow(gamma0-eta,3);
        X4[m] *= pref/pow(gamma0-eta,4);
    }

    vector<double> Pn(l+1,0.0);
    for (int m = 0; m<l+1; m++){
        Pn[m] = 3*t*t*X4[m]-6*gamma0*t*(1+t)*X3[m]+((3+2*eta*eta)*(1+t*t)+12*gamma0*gamma0*t)*X2[m];
        Pn[m] += -2*gamma0*(3+2*eta*eta)*(1+t)*X1[m]+(3+4*eta*eta+4*eta*eta*eta*eta)*X0[m];
    }

    double prefactor = 3*t/(32*gamma0*pow(eta,5));
    double lsum = 0;
    for (int m = 0; m<l+1; m++){
        double pref = Binomial(l,m)*Binomial(l+m,m);
        lsum += pref*Pn[m];
    }
    return prefactor*lsum;
}

double MultipoleKernel::Calculate_integrated(){
    if (s<s_low|| s>s_high) { return 0.0; }
    return mus_Int()*t;
}

double MultipoleKernel::Calculate_formula(){
    if (!electron_anisotropy) {
        if (s<s_low|| s>s_high) { return 0.0; }
        Calculate_formula_variables();
        double prefactor = 3.0/(32.0*pow(eta,6));
        double lsum = 0;
        for (int n = 0; n<l+1; n++){
            double pref = Binomial(l,n)*Binomial(l+n,n);
            lsum += pref*(4*pow(eta,6)*K1[n]+(1+t)*K2[n]/t/pow(eta,2*n));
        }
        return prefactor*lsum*t;
    }
    return Calculate_electron_multipoles();
}

double MultipoleKernel::Calculate_stable(){
    if (s<s_low|| s>s_high) { return 0.0; }
    switch (l) {
    case 0:
        if (electron_anisotropy && eta >= 0.00151){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.00111){return Calculate_formula();}
        return Calculate_integrated();
    case 1:
        if (electron_anisotropy && eta >= 0.00599){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.0121){return Calculate_formula();}
        return Calculate_integrated();
    case 2:
        if (electron_anisotropy && eta >= 0.0161){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.0416){return Calculate_formula();}
        return Calculate_integrated();
    case 3:
        if (electron_anisotropy && eta >= 0.0331){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.0936){return Calculate_formula();}
        return Calculate_integrated();
    case 4:
        if (electron_anisotropy && eta >= 0.0538){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.1620){return Calculate_formula();}
        return Calculate_integrated();
    case 5:
        if (electron_anisotropy && eta >= 0.0840){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.244){return Calculate_formula();}
        return Calculate_integrated();
    case 6:
        if (electron_anisotropy && eta >= 0.126){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.333){return Calculate_formula();}
        return Calculate_integrated();
    case 7:
        if (electron_anisotropy && eta >= 0.165){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.419){return Calculate_formula();}
        return Calculate_integrated();
    case 8:
        if (electron_anisotropy && eta >= 0.219){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.508){return Calculate_formula();}
        return Calculate_integrated();
    case 9:
        if (electron_anisotropy && eta >= 0.301){return Calculate_formula();}
        if (!electron_anisotropy && eta >= 0.632){return Calculate_formula();}
        return Calculate_integrated();
    default:
        return Calculate_integrated();
    }
}

//==================================================================================================
// Integration Class
//==================================================================================================
IntegralKernel::IntegralKernel()
    : IntegralKernel(0.0, 0.0, 0.0, 1.0e-4) {}

IntegralKernel::IntegralKernel(double x_i, double betac_i, double muc_i, double eps_Int_i){
    Int_eps = eps_Int_i;
    betac = betac_i;
    muc = muc_i;
    x = x_i;
    eta=s=xp=0.0;
    l = 0;
    run_mode = "";
    etaDistribution = plainDistribution;
    MK = MultipoleKernel(l, s, eta, Int_eps);
    xfac = 1.0;
    fixed_eta=false;
    electron_anisotropy = true;
}

IntegralKernel::IntegralKernel(double x_i, Parameters fp)
    : IntegralKernel(fp.calc.xfac*x_i, fp.betac, fp.calc.mucc, fp.relative_accuracy) {
        xfac = fp.calc.xfac;
    }

void IntegralKernel::Update_x(double x_i){
    x = xfac*x_i;
}

void IntegralKernel::Calculate_shared_variables(){
    xp = x*exp(s);
    MK = MultipoleKernel(l, s, eta, Int_eps);
    MK.electron_anisotropy = electron_anisotropy;
}

double IntegralKernel::Calculate_kernel(int l_i){
    MK.Update_l(l_i);
    return MK.Calculate_stable();
}

double IntegralKernel::Calculate_monopole(int l_i){
    double Sx = xk_dk_nPl(0,x);
    double Sp= xk_dk_nPl(0,xp);
    double dist = etaDistribution(eta);
    double F = Calculate_kernel(l_i);
    return dist*F*(Sp-Sx);
}

double IntegralKernel::Calculate_dipole(int l_i){
    double Gx = xk_dk_nPl(1,x);
    double Gp = xk_dk_nPl(1,xp);
    double dist = etaDistribution(eta);
    double F = Calculate_kernel(l_i);
    return dist*F*(Gp-Gx);
}

double IntegralKernel::Calculate_quadrupole(int l_i){
    double Qx = xk_dk_nPl(2,x);
    double Qp = xk_dk_nPl(2,xp);
    double dist = etaDistribution(eta);
    double F = Calculate_kernel(l_i);

    return dist*F*(Qp-Qx);
}

double IntegralKernel::Calculate_monopole_correction(int l_i){
    double Gx = xk_dk_nPl(1,x);
    double Gp = xk_dk_nPl(1,xp);
    double Qx = xk_dk_nPl(2,x);
    double Qp = xk_dk_nPl(2,xp);
    double dist = etaDistribution(eta);
    double F = Calculate_kernel(l_i);
    
    return dist*F*((Qp-Qx)+3*(Gp-Gx));
}

double IntegralKernel::sig_Boltzmann_Compton(double int_eta){
    eta = int_eta;
    Calculate_shared_variables();

    double r = 0.0; 

    if(run_mode=="monopole"){
        r = Calculate_monopole(l);
    } 
    else if(run_mode=="dipole"){ 
        r = Calculate_dipole(l);
    }
    else if(run_mode=="quadrupole"){
        r = Calculate_quadrupole(l);
    }
    else if(run_mode=="monopole_corr"){ 
        r = Calculate_monopole_correction(l);
    }
    else if(run_mode=="kernel"){
        double dist = etaDistribution(eta);
        double F = Calculate_kernel(l);
        r = dist*F;
    }
    else {
        r = Calculate_monopole(l);
    }
    return eta*eta*r;
}

double IntegralKernel::eta_Int(double int_s){
    s = int_s;
    if (s==0){ s = 1e-15; }
    if (fixed_eta){ return sig_Boltzmann_Compton(eta)/eta/eta; }

    double a=sinh(fabs(s)/2.0), b = 30.0;//lim=30.0, b=lim*(1.0+0.5*lim*0.05);
    double epsrel=Int_eps, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sig_Boltzmann_Compton(int_var);});
}

double IntegralKernel::s_Int(){
    double a=-2.0, b=2.0;
    double epsrel=Int_eps, epsabs=1.0e-300;
    double integral = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->eta_Int(int_var);});
    return integral;
}

double IntegralKernel::compute_kernel(int l_i, double s_i, electronDistribution eDistribution, bool e_anis){
    electron_anisotropy = e_anis;
    run_mode = "kernel";
    l = l_i;
    etaDistribution = eDistribution;
    electron_anisotropy = false;
    return eta_Int(s_i);
}

double IntegralKernel::compute_distortion(string mode, electronDistribution eDistribution, int l_i, bool e_anis){
    electron_anisotropy = e_anis;
    run_mode="monopole";
    if(mode=="monopole" || mode=="dipole" || mode=="quadrupole" || mode=="monopole_corr"){
        run_mode = mode;
    }
    etaDistribution = eDistribution;
    l=l_i;
    double result = s_Int();
    electron_anisotropy = false;
    return result;
}

double IntegralKernel::compute_distortion_fixed_eta(string mode, double eta_i, int l_i, bool e_anis){
    fixed_eta = true;
    eta = eta_i;
    double result = compute_distortion(mode, etaDistribution, l_i, e_anis);
    fixed_eta = false;
    return result;
}

//==================================================================================================
//
// 3D integration carried out using Patterson scheme
//
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
//
// mode == "monopole"      --> only scattering of monopole without second order kinematic corr
// mode == "dipole"        --> only scattering of dipole     (first order kinematic correction)
// mode == "quadrupole"    --> only scattering of quadrupole (second order kinematic correction)
// mode == "monopole_corr" --> only scattering of second order kinematic correction to monopole
// mode == "all"           --> all terms added
// mode == "kin"           --> only kinematic terms
//
//==================================================================================================

void compute_SZ_distortion_kernel(vector<double> &Dn, Parameters &fp, bool DI, 
                                        std::function<double(double)> eDistribution, int l, bool e_anis){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    for(int k = 0; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(fp.rare.RunMode, eDistribution, l, e_anis);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

void compute_averaged_kernel(vector<double> &Dn, Parameters &fp, std::function<double(double)> eDistribution, int l_i, bool e_anis){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(0.1, fp);
    for(int k = 0; k < fp.gridpoints; k++){
        Dn[k] = szDistortion.compute_kernel(l_i, fp.kernel.srange[k], eDistribution, e_anis);
    }
}

void compute_SZ_distortion_fixed_eta(vector<double> &Dn, Parameters &fp, bool DI, double eta, int l_i, bool e_anis){
    Dn.resize(fp.gridpoints);
    IntegralKernel szDistortion = IntegralKernel(fp.xcmb[0], fp);
    for(int k = 0; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion_fixed_eta(fp.rare.RunMode, eta, l_i, e_anis);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}

//==================================================================================================
//
// 5D integration carried out using Patterson scheme
// eps_Int : relative accuracy for numerical integration (lower than 10^-6 is hard to achieve)
// Calculates the 5D distortion for an arbitrary 
//
//==================================================================================================

nonThermal5D::nonThermal5D(){
    Int_eps=x=mu=mup=muep=eta=phi=phip=dx=xp=dsig=0.0;
    exp_mx=dex=exp_mxp=dexp=nx=nxp=0.0;
    xfac = 1.0;
}

nonThermal5D::nonThermal5D(double x_i, double eps_Int_i){
    Int_eps = eps_Int_i;
    x = x_i;
    mu=mup=muep=eta=phi=phip=dx=xp=dsig=exp_mx=dex=exp_mxp=dexp=nx=nxp=0.0;
    xfac = 1.0;
    etamuphiDist = plainFullDistribution;
}

nonThermal5D::nonThermal5D(double x_i, Parameters fp)
    : nonThermal5D(fp.calc.xfac*x_i, fp.relative_accuracy) {
        xfac = fp.calc.xfac;
    }

void nonThermal5D::Update_x(double x_i){
    x = xfac*x_i;
}

void nonThermal5D::Calculate_shared_variables(){
    //Standard definitions
    //==============================================================================================
    // mu == cosine of angle of gamma  and beta 
    // mup== cosine of angle of gamma' and beta
    //==============================================================================================    
    double beta=eta/sqrt(1.0+eta*eta);
    double gamma2=1.0+eta*eta;
    //
    //muep=mu*mup+cos(phi-phip)*sqrt( (1.0-mu*mu)*(1.0-mup*mup) );
    mu = muep*mup+cos(phi-phip)*sqrt((1.0-muep*muep)*(1.0-mup*mup));
    //
    double kappa =1.0-beta*mu;
    double kappap=1.0-beta*muep;
    double zeta=1.0/gamma2/kappa/kappap;
    //
    double alpha_sc=1.0-mup;

    //The variables based on x
    dx=x*beta*(muep-mu)/kappap;
    xp=x+dx;

    //d2sigma/dmu/dmup as defined equation 2 CNSN
    dsig = 3.0/8.0/PI*kappa/pow(kappap, 2)/gamma2*(1.0-zeta*alpha_sc*(1.0-0.5*zeta*alpha_sc));
}

void nonThermal5D::Calculate_HigherOrder_Variables(){
    exp_mx=exp(-x);
    dex=one_minus_exp_mx(x, exp_mx); 
    exp_mxp=exp(-xp);
    dexp=one_minus_exp_mx(xp, exp_mxp);
    nx=exp_mx/dex;
    nxp=exp_mxp/dexp;
}

double nonThermal5D::Calculate_monopole(){
    Calculate_HigherOrder_Variables();
    return (nxp-nx);
}

double nonThermal5D::sig_Boltzmann_Compton(double int_phip){
    phip = int_phip;
    Calculate_shared_variables();
    
    double r = 0.0; 

    r = Calculate_monopole();
    double dist = etamuphiDist(eta, muep, phip)/4.0/PI;
    return dist*dsig*r;
}

double nonThermal5D::phip_Int(double int_phi){
    phi = int_phi;
    double a=0.0, b=TWOPI;
    double epsrel=Int_eps*0.6, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->sig_Boltzmann_Compton(int_var);});;
}

double nonThermal5D::phi_Int(double int_mup){
    mup = int_mup;
    double a=0.0, b=TWOPI;
    double epsrel=Int_eps*0.7, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->phip_Int(int_var);});
}

double nonThermal5D::mup_Int(double mu_int){
    muep = mu_int;
    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.8, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->phi_Int(int_var);});
}

double nonThermal5D::mue_Int(double eta_int){
    eta = eta_int;
    double a=-1.0, b=1.0;
    double epsrel=Int_eps*0.9, epsabs=1.0e-300;
    double integral = Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mup_Int(int_var);});
    return integral;
}

double nonThermal5D::eta_Int(){
    double a=0.0, lim=30.0, b=lim*(1.0+0.5*lim*0.05);
    double epsrel=Int_eps, epsabs=1.0e-300;
    
    return Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, [this](double int_var) { return this->mue_Int(int_var);});
}

double nonThermal5D::compute_distortion(fullElectronDistribution empDistribution){
    etamuphiDist = empDistribution;
    return eta_Int();
}

void compute_SZ_distortion_5DnonThermal(vector<double> &Dn, Parameters &fp, bool DI, 
                                        std::function<double(double, double, double)> empDistribution){
                                            Dn.resize(fp.gridpoints);
    nonThermal5D szDistortion = nonThermal5D(fp.xcmb[0], fp);
    for(int k = 0; k < fp.gridpoints; k++){
        szDistortion.Update_x(fp.xcmb[k]);
        Dn[k] = fp.Dtau*szDistortion.compute_distortion(empDistribution);
        if (DI) { Dn[k] *= pow(fp.xcmb[k],3.0)*fp.rare.Dn_DI_conversion(); }
    }
}