#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <nlopt.hpp>

#include "iofunctions.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif 

#ifndef QUASI1D_H_INCLUDED
#define QUASI1D_H_INCLUDED

//Efficient way of computing binomial coefficients
double binomialCoefficients(int n, int k) {
   double C[k+1];
   memset(C, 0, sizeof(C));
   C[0] = 1;
   for (int i = 1; i <= n; i++) {
      for (int j = std::min(i, k); j > 0; j--)
         C[j] = C[j] + C[j-1];
   }
   return C[k];
}

//Class Quasi1d
class Dumbbells2 {       
  public:         
    int chainsize;
    int nsize; //size of the discretization (number of species).
    double bp; //pressure
    double beta; //beta
    double r0; //corona size
    double density;
    double amin;
    double corrlength1;
    double corrlength2;
    double l0;
    std::vector<double> limvaluesample;
    double fmin, fmax;

    double A2; //A coming from the eigenvalue equation.
    Eigen::Matrix<double, Eigen::Dynamic, 1> xvalues; //vector of molar fractions.
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> am;
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> M;
    //Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> distances;
    Eigen::Matrix<double,Eigen::Dynamic, 1> distances;

    //To compute eigenvalues
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> P;

    //Constructor
    Dumbbells2(InputData *indat){
        chainsize = indat->chainsize;
        nsize = indat->nsize;
        beta = 1.0;

        //Update fmin and fmax
        fmin = dfmin(indat->bp);
        fmax = 0.52394;
       
        //Internal stuff
        limvaluesample = std::vector<double>(3,0.0);
        am.resize(nsize, nsize);
        M.resize(nsize, nsize);
        //distances.resize(chainsize*chainsize,1);
        distances.resize(4*chainsize,1);
        UpdateMatrixA();
        amin=am.minCoeff();
        std::cout << "Mininum contact distance = " << amin << std::endl;
        
    }

    double sample(int s, int i, int j){
        if (s==0)  return cos(2.0*(phi(i)-phi(j)));
        else if (s==1) return cos(2.0*phi(i))*cos(2.0*phi(j));
        else if (s==2) return sin(2.0*phi(i)) * sin(2.0*phi(j));
        return 0.0;
    }

    //Set pressure of the system and update xvalues and Lmean.
    int SetPressure(double bp0){
        //Create matrix M
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                M(i,j) = omega(bp0, i, j);
            }
        }
        //Solve eigenvalue equation
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es;
        es.compute(M);
        xvalues = es.eigenvectors().col(nsize-1);
        //Make it symmetric in case it was wrongly calculated
        auto xvaluesold = xvalues;
        for(int i = 0; i<nsize; i++){
            xvalues(i) = (xvaluesold(i)+xvaluesold(nsize-1-i))/2.0;
        }
        //std::cout << "eig = " << es.eigenvalues() << std::endl;
        A2 = 1.0/es.eigenvalues()(nsize-1);
        //Prepare values
        double sum = 0;
        for(int i = 0; i<nsize; i++){
            xvalues(i) = xvalues(i)*xvalues(i);
            sum += xvalues(i);
        }
        xvalues /= sum;
        
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sorteigvals = es.eigenvalues().cwiseAbs();

        corrlength1 = 1.0/log( sorteigvals(nsize-1)/sorteigvals(nsize-2) );
        corrlength2 = 1.0/( log(sorteigvals(nsize-1)) - log(sorteigvals(nsize-2)) );
        l0 = sorteigvals(nsize-1);

        //Update internal bp
        this->bp = bp0;
        this->density = Density(bp);
        limvaluesample[0] = GetLimValueSample(0); 
        limvaluesample[1] = GetLimValueSample(1);
        limvaluesample[2] = GetLimValueSample(2);

        //std::cout << "xvalues = " << xvalues << std::endl;
        //std::cout << "Aval = " << Aval << std::endl;

        return 0;
    }
    
    //If instead of pressure the user wishes to set the density, use this function.
    void SetDensity(double rho){
        std::cout << rho << std::endl;
        //Bolzano method to compute the pressure associated with that density
        double a = 0.000001;
        double b=750;
  
        double c = a;
        while ((b-a) >= 0.00000001){
 
            // Find middle point
            c = (a+b)/2.0;
           
            // Decide the side to repeat the steps
            if ((Density(c)-rho)*(Density(a)-rho) < 0)   b = c;
            else  a = c;
        }

        if (abs(c-750.0)<=0.0001){
            std::cout << "WARNING: required density set too high, numerical errors might occur." << std::endl;
        }

        //Set the corresponding pressure
        this->density = rho;
        SetPressure(c);
        
    }

    // Compute density of the system at a certain bp
    double Density(double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double nsum = 0;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                nsum += sqrt(xvalues(i)*xvalues(j))*omegad(bp,i,j)*A2;
            }
        }
        return -beta/nsum;
    }

    double ComputeZ(double bp0){
        if (bp0 != bp)  SetPressure(bp0); //Check the pressure is correct
        double Zx = 0.0;
        
        double temp = 0.0;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                Zx = Zx + sqrt(xvalues(i)*xvalues(j))*exp(-bp*(am(i,j)-1.0))*am(i,j);
            }
        }
        Zx = A2*Zx + 1.0;

        return Zx;
    }

    double PosMaxAngle(double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        Eigen::MatrixXf::Index   maxIndex;
        xvalues.maxCoeff(&maxIndex);
        return fabs(phi(maxIndex));
    }

    double HeightMaxAngle(double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        return  xvalues.maxCoeff();
    }

    double OrderParameter(double bp0){
        if (bp0 != bp)  SetPressure(bp0);

        double nsum = 0.0;
        for(int j=0; j<nsize; j++){
            nsum += xvalues(j)*cos(2.0*phi(j));
        }

        return nsum;
    }

    double GetLimValueSample(int s){
        double lim = 0.0;
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                if (ix!=jx){
                    lim += 2.0*xvalues(ix)*xvalues(jx)*sample(s,ix,jx);
                }
                else{
                    lim += xvalues(ix)*xvalues(jx)*sample(s,ix,jx);
                } 
            }
        }
        return lim;
    }
    
    //Calculamos el determinante para los polos
    double dtm(const std::vector<double>& x)
    {
    using cplx = std::complex<double>;

    // s = - (x0 + i x1)
    const cplx s = cplx{x[0], x[1]};

    // z = s + bp
    const cplx z = s + this->bp;

    // Guard against division by (near) zero if you want robustness
    // (optional — Python would just produce inf/nan)
    // if (std::abs(z) < 1e-15) return std::numeric_limits<double>::infinity();

    // omega(z) = exp(-a*z)/z, with ELEMENT-WISE exp (NumPy behavior)
    const Eigen::ArrayXXcd omega =
        (-this->am.array().cast<cplx>() * z).exp() / z;

    //std::cout << omega << "\n" <<std::endl;
    // P = I - A2 * omega   (A2 is scalar)
    const int n = static_cast<int>(this->am.rows());
    Eigen::MatrixXcd P = Eigen::MatrixXcd::Identity(n, n) - (exp(bp)*this->A2 * omega.matrix());
    //std::cout << P << std::endl;

    // dtm(x) = |det(P)|
    return std::abs(P.determinant());
    }

    static double dtm_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* data)
    {
        (void)grad; // Nelder–Mead ignores gradients
        // Cast back to your object
        Dumbbells2* self = static_cast<Dumbbells2*>(data);
        return self->dtm(x);
    }

    std::complex<double> minimize(std::vector<double> x0, double *funcvalue){
        nlopt::opt opt(nlopt::LN_NELDERMEAD, 2);

        // Optional bounds (SciPy allows unbounded; NLopt works either way)
        if (fabs(x0[1])>1e-5){
            opt.set_upper_bounds({0.85*x0[0], 1.15*x0[1]});
            opt.set_lower_bounds({1.15*x0[0], 0.85*x0[1]});
        }
        else {
            opt.set_upper_bounds({0.90*x0[0], HUGE_VAL});
            opt.set_lower_bounds({1.10*x0[0], -0.1});

        }

        opt.set_min_objective(Dumbbells2::dtm_wrapper, this);

        // Similar stopping criteria you’d set in SciPy options
        opt.set_xtol_rel(1e-20);     // like xatol-ish (relative)
        opt.set_ftol_rel(1e-29);    // like fatol-ish (relative)
        opt.set_maxeval(5000);

        // Initial guess (like x0 in SciPy)
        //std::vector<double> x = {-1.0, 5.0};
        double funcmin = 0.0;

        try {
            nlopt::result r = opt.optimize(x0, funcmin);
            //std::cout << "Result code: " << r << "\n";
            //std::cout << "x* = [" << x0[0] << ", " << x0[1] << "]\n";
            //std::cout << "f(x*) = " << funcmin << "\n";
            //std::cout << "------------------------" << "\n";
        } catch (const std::exception& e) {
            std::cerr << "NLopt failed: " << e.what() << "\n";
        }
        *funcvalue=funcmin;
        return std::complex<double>(x0[0],x0[1]);
    }
    //Calculamos Gx de forma analitica cuando x es pequeño
    std::vector<double> GxAna(double x, double bp0, int nmax){
        if (bp0 != bp)  SetPressure(bp0);
        std::vector<double> gtotal = std::vector<double>(8,0.0);

        //Optimizacion: vemos en que capa de proximos vecinos cae x y ajustamos el valor de nmax si es necesario
        int nshell = int(x/amin);
        if (nshell == 0) return gtotal;
        else if (nshell < nmax)  nmax = nshell; //keep lowest n value because they are going to yield the same result

        double gtemp = 0.0;
        double gmean = 0.0;
        double gorient = 0.0;
        double gcos = 0.0;
        double gsin = 0.0;
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                //Choose potential
                gtemp = Gijx(x,bp,nmax,ix,jx);
                if (ix!=jx){
                    gmean += 2.0*xvalues(ix)*xvalues(jx)*gtemp;
                    gorient += 2.0*xvalues(ix)*xvalues(jx)*gtemp*sample(0,ix,jx);
                    gcos += 2.0*xvalues(ix)*xvalues(jx)*gtemp*sample(1,ix,jx);
                    gsin += 2.0*xvalues(ix)*xvalues(jx)*gtemp*sample(2,ix,jx);
                }
                else{
                    gmean += xvalues(ix)*xvalues(jx)*gtemp;
                    gorient += xvalues(ix)*xvalues(jx)*gtemp*sample(0,ix,jx);
                    gcos += xvalues(ix)*xvalues(jx)*gtemp*sample(1,ix,jx);
                    gsin += xvalues(ix)*xvalues(jx)*gtemp*sample(2,ix,jx);
                } 
                //Get partials G(r)
                if(ix==(nsize-1)/2 && jx==(nsize-1)/2){
                    gtotal[4] = gtemp -1.0;
                }
                if(ix==(nsize-1)/3 && jx==2*(nsize-1)/3){
                    gtotal[5] = gtemp -1.0;
                }
                if(ix==(nsize-1)/3 && jx==(nsize-1)/3){
                    gtotal[6] = gtemp-1.0;
                }
                if(ix==(nsize-1)/3 && jx==(nsize-1)/2){
                    gtotal[7] = gtemp-1.0;
                }
            }
        }
        gtotal[0] = gmean-1.0;
        //Divide all orientation functions by g(z)
        gtotal[1] = gorient/gmean - limvaluesample[0];
        gtotal[2] = gcos/gmean - limvaluesample[1];
        gtotal[3] = gsin/gmean - limvaluesample[2];

        return gtotal;
    }

    double Gijx(double x, double bp0, int nmax, int i, int j){
        if (bp0 != bp)  SetPressure(bp0);

        if (int(x/am(i,j)) == 0)  return 0;

        double gij = 0.0;
        for (int ni=0; ni<nmax; ni++){
            gij += qij(x,bp,ni+1,i,j);
        }
        gij /=(density*sqrt(xvalues(i)*xvalues(j)));
        return gij;
    }

    //Calculamos Gx como la transformada inversa de Laplace de Gs
    std::vector<double> GxInv(double t, double bp0, int ntr, int meuler, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::vector<std::complex<double>> Gstemp = Gs<std::complex<double>>(xx,bp0,Qtt, Qttm);
        std::complex<double> suma = Gstemp[0]/2.0;
        std::complex<double> suma1 = Gstemp[1]/2.0;
        std::complex<double> suma2 = Gstemp[2]/2.0;
        std::complex<double> suma3 = Gstemp[3]/2.0;
        std::complex<double> suma4 = Gstemp[4]/2.0;
        std::complex<double> suma5 = Gstemp[5]/2.0;
        std::complex<double> suma6 = Gstemp[6]/2.0;
        std::complex<double> suma7 = Gstemp[7]/2.0;

        std::complex<double> dummy;
        for (int i=1;i<ntr+1;i++){
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx, i*hh),bp0,Qtt, Qttm);
            dummy = std::complex<double>(pow(-1,i),0);
            suma += dummy * Gstemp[0];
            suma1 += dummy * Gstemp[1];
            suma2 += dummy * Gstemp[2];
            suma3 += dummy * Gstemp[3];
            suma4 += dummy * Gstemp[4];
            suma5 += dummy * Gstemp[5];
            suma6 += dummy * Gstemp[6];
            suma7 += dummy * Gstemp[7];
        }

        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su1(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su2(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su3(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su4(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su5(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su6(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su7(meuler+2,std::complex<double>(0,0));
        su[0] = suma;
        su1[0] = suma1;
        su2[0] = suma2;
        su3[0] = suma3;
        su4[0] = suma4;
        su5[0] = suma5;
        su6[0] = suma6;
        su7[0] = suma7;

        //#avgsu
        //std::complex<double> avgsu(0,0);
        std::complex<double> avgsu(0,0);
        std::complex<double> avgsu1(0,0);
        std::complex<double> avgsu2(0,0);
        std::complex<double> avgsu3(0,0);
        std::complex<double> avgsu4(0,0);
        std::complex<double> avgsu5(0,0);
        std::complex<double> avgsu6(0,0);
        std::complex<double> avgsu7(0,0);
        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            dummy = pow(-1,nn);
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx,nn*hh),bp0,Qtt, Qttm);
            su[i+1] = su[i] + dummy * Gstemp[0];
            su1[i+1] = su1[i] + dummy * Gstemp[1];
            su2[i+1] = su2[i] + dummy * Gstemp[2];
            su3[i+1] = su3[i] + dummy * Gstemp[3];
            su4[i+1] = su4[i] + dummy * Gstemp[4];
            su5[i+1] = su5[i] + dummy * Gstemp[5];
            su6[i+1] = su6[i] + dummy * Gstemp[6];
            su7[i+1] = su7[i] + dummy * Gstemp[7];
            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            //avgsu += bn*su[i];
            avgsu += bn*su[i+1];
            avgsu1 += bn*su1[i+1];
            avgsu2 += bn*su2[i+1];
            avgsu3 += bn*su3[i+1];
            avgsu4 += bn*su4[i+1];
            avgsu5 += bn*su5[i+1];
            avgsu6 += bn*su6[i+1];
            avgsu7 += bn*su7[i+1];

        }
        double pp = pow(2.0,meuler);
        //std::complex<double> fun = uu*avgsu/(pp);
        std::vector<double> fun = std::vector<double>(8,0.0);
        fun[0] = real(uu*avgsu/(pp)) - 1.0;
        fun[1] = real(uu*avgsu1/(pp))/(fun[0]+1.0) - limvaluesample[0];
        fun[2] = real(uu*avgsu2/(pp))/(fun[0]+1.0) - limvaluesample[1];
        fun[3] = real(uu*avgsu3/(pp))/(fun[0]+1.0) - limvaluesample[2];
        fun[4] = real(uu*avgsu4/(pp))-1.0;
        fun[5] = real(uu*avgsu5/(pp))-1.0;
        fun[6] = real(uu*avgsu6/(pp))-1.0;
        fun[7] = real(uu*avgsu7/(pp))-1.0;

        return fun;
        //return 0.0;

    }

    template<class T>
    std::vector<T> Gs(T s, double bp0, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){

        std::vector<T> gtotal = std::vector<T>(8,0.0);

        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                Qtt(i,j) = A2*omega<T>(s+bp0,i,j)*exp(-s);
                Qttm(i,j) = -Qtt(i,j);
                if(i==j) Qttm(i,j) = 1.0+Qttm(i,j);
            }
        }

        //Invert matrix and all that
        Qtt = Qtt*(Qttm.inverse());
        T gmean = 0;
        T gorient  = 0;
        double gmeanRE = 0;
        double gmeanIM = 0;
        double gmeanREorient = 0;
        double gmeanIMorient = 0;
        double gREcos = 0;
        double gIMcos = 0;
        double gREsin = 0;
        double gIMsin = 0;

        double density = Density(bp);
        T Gijs = 0.0;
        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                //std::cout << s << "-> " << 1.0/s << std::endl;
                Gijs = sqrt(xvalues(i)*xvalues(j))*Qtt(i,j);
                if (i!=j){
                    gmeanRE += std::real(2.0*Gijs);
                    gmeanIM += std::imag(2.0*Gijs);
                    gmeanREorient += std::real(2.0*Gijs)*sample(0,i,j);
                    gmeanIMorient += std::imag(2.0*Gijs)*sample(0,i,j);
                    gREcos += std::real(2.0*Gijs)*sample(1,i,j);
                    gIMcos += std::imag(2.0*Gijs)*sample(1,i,j);
                    gREsin += std::real(2.0*Gijs)*sample(2,i,j);
                    gIMsin += std::imag(2.0*Gijs)*sample(2,i,j);
                }
                else{
                    gmeanRE += std::real(Gijs);
                    gmeanIM += std::imag(Gijs);
                    gmeanREorient += std::real(Gijs)*sample(0,i,j);
                    gmeanIMorient += std::imag(Gijs)*sample(0,i,j);
                    gREcos += std::real(Gijs)*sample(1,i,j);
                    gIMcos += std::imag(Gijs)*sample(1,i,j);
                    gREsin += std::real(Gijs)*sample(2,i,j);
                    gIMsin += std::imag(Gijs)*sample(2,i,j);
                }
                //Get partials G(s)
                if(i==(nsize-1)/2 && j==(nsize-1)/2){
                    gtotal[4] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                }
                if(i==(nsize-1)/3 && j==2*(nsize-1)/3){      
                    gtotal[5] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                }
                if(i==(nsize-1)/3 && j==(nsize-1)/3){                  
                    gtotal[6] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                }
                if(i==(nsize-1)/3 && j==(nsize-1)/2){                  
                    gtotal[7] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                }
            }
        }
        if constexpr (std::is_same_v<T, std::complex<double>>){
            gtotal[0] = std::complex<double>(gmeanRE/density, gmeanIM/density);
            gtotal[1] = std::complex<double>(gmeanREorient/density, gmeanIMorient/density);
            gtotal[2] = std::complex<double>(gREcos/density, gIMcos/density);
            gtotal[3] = std::complex<double>(gREsin/density, gIMsin/density);
        }
        else{
            gtotal[0] = gmeanRE/density;
            gtotal[1] = gmeanREorient/density;
            gtotal[2] = gREcos/density;
            gtotal[3] = gREsin/density;
        }
        return gtotal;
    }

    
    double GijxInv(double t, double bp0, int ntr, int meuler, int ii, int jj, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::complex<double> suma = Gijs(xx,bp0, ii, jj, Qtt, Qttm)/2.0;

        std::complex<double> argument(0,0);
        for (int i=1;i<ntr+1;i++){
            suma = suma + std::complex<double>(pow(-1,i),0) * Gijs(std::complex<double>(xx, i*hh),bp0, ii, jj, Qtt, Qttm);
        }

        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        su[0] = suma;

        //#avgsu
        std::complex<double> avgsu(0,0);
        std::complex<double> avgsu1(0,0);
        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            su[i+1] = su[i] + pow(-1,nn) * Gijs(std::complex<double>(xx,nn*hh),bp0, ii, jj, Qtt, Qttm);
            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            avgsu += bn*su[i];
            avgsu1 += bn*su[i+1];

        }
        double pp = pow(2.0,meuler);
        std::complex<double> fun = uu*avgsu/(pp);
        std::complex<double> fun1 = uu*avgsu1/(pp);
        return real(fun1);
        //return 0.0;

    }

    template<class T>
    T Gijs(T s,  double bp0, int ii, int jj, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        if (bp0 != bp && bp0 >1e-10)  SetPressure(bp0);
        std::complex<double> gtotal;

        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                Qtt(i,j) = A2*omega<T>(s+bp0,i,j)*exp(-s);
                Qttm(i,j) = -Qtt(i,j);
                if(i==j) Qttm(i,j) = 1.0+Qttm(i,j);
            }
        }

        //Invert matrix and all that
        Qtt = Qtt*(Qttm.inverse());

        gtotal += Qtt(ii,jj)/sqrt(xvalues(ii)*xvalues(jj));
        return gtotal.real()/density;
    }

    double ProbDensity(double x, int n, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        
        double gtotal = 0.0;

        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                if (ix!=jx){
                    gtotal += 2.0*sqrt(xvalues(jx)*xvalues(ix))*qij(x, bp0, n, ix, jx);
                }
                else{
                    gtotal += sqrt(xvalues(jx)*xvalues(ix))*qij(x, bp0, n, ix, jx);
                } 
            }

        }
        return gtotal;
    }
  
    /*double a(int i, int j){ //Real system. Comment out for the Toy model
        double f1 = phi(i);
        double f2= phi(j);
        double temp;
        double dist = 0;
        double khinit = -((double)chainsize-1.0)/2.0;
        double kk, hh;
        for(int k=0; k<chainsize; k++){
            for(int h=0; h<chainsize; h++){
                kk = khinit+k;
                hh = khinit+h;
                temp = (kk*cos(f1) + hh*cos(f2))*(kk*cos(f1) + hh*cos(f2));
                if(temp >= 1.0){
                    distances(k,h) = 0.0;
                }
                else{
                    distances(k,h) = kk*sin(f1) + hh*sin(f2) + sqrt(1.0 - temp);
                }
            }
        }
        return distances.maxCoeff(); // /(double)chainsize;
    }
    */

    double contact(double i, double j, double kk, double hh){
        double f1 = phi(i);
        double f2 = phi(j);
        double temp = (kk*cos(f1) + hh*cos(f2))*(kk*cos(f1) + hh*cos(f2));
        if(temp >= 1.0)
            return 0.0;
        else
            return kk*sin(f1) + hh*sin(f2) + sqrt(1.0 - temp);
    }

    double a(int i, int j){
        double khinit = -(chainsize-1.0)/2.0;
        int idx = 0;
        double klist[] = {-(chainsize-1.0)/2.0, (chainsize-1.0)/2.0};
        for (double kk : klist){
            for (int h=0; h<chainsize; h++){
                distances(idx,0) = contact(i,j,kk,khinit+h);
                idx = idx+1;
                distances(idx,0) = contact(i,j,khinit+h,kk);
                idx = idx+1;
            }
        }
        return distances.maxCoeff();
    }
      
    template<class T>
    T omega(T s, int i, int j){
        return exp(-(am(i,j)-1.0)*s)/s;
    }

    template<class T>
    T omegad(T s, int i, int j){
        return -exp(-(am(i,j)-1.0)*s)*(1.0+am(i,j)*s)/(s*s);
    }

    double dfmin(double bp0){
        double k=0.99;
        if (bp0>k*1000 && bp0<2000.0) return 0.40;
        else if (bp0>k*2000.0 && bp0<3000.0) return 0.42;
        else if (bp0>k*3000.0 && bp0<4000.0) return 0.43;
        else if (bp0>k*4000.0 && bp0<5000.0) return 0.44;
        else if (bp0>k*5000.0 && bp0<6000.0) return 0.45;
        else if (bp0>k*6000.0 && bp0<7000.0) return 0.46;
        else if (bp0>k*7000.0 && bp0<8000.0) return 0.465;
        else if (bp0>k*8000.0 && bp0<10000.0) return 0.47;
        else if (bp0>k*10000.0 && bp0<20000.0) return 0.48;
        else if (bp0>k*20000.0 && bp0<30000.0) return 0.49;
        else if (bp0>k*30000.0 && bp0<100000.0) return 0.50;
        else if (bp0>k*100000) return 0.511882;
        else return 0.0;

    }

    double phi(int i){
        if (nsize==2){
            return -M_PI/6.0 + (double)i*M_PI/(double)(nsize-1.0)/3;
        }
         return -M_PI/2.0 + (double)i*M_PI/(double)(nsize-1.0); //Esta es la línea buena

       /* double angle;
        if (i<nsize/2){
            angle = -fmax + (double)i*(fmax-fmin)/((nsize/2-1.0));
        }
        else{
            angle = fmin + ((double)i-nsize/2)*(fmax-fmin)/((nsize/2-1.0));

        }
        return angle;
        */
    }

    int UpdateMatrixA(){
        double temp;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                //std::cout << i << ", " << j << " -> " << a(i,j) << std::endl;
                am(i,j) = a(i,j);
            }
        }
        return 0;
    }  

    double f(int n, double x, double b){
        if(x-b>0){
            return exp(-bp*x)*pow(x-b,n-1)/tgamma(n); //tgamma(n+1)=factorial(n)
        }
        else{
            return 0.0;
        }
    }

    double qij(double x, double bp0, int n, int i, int j){
      if (bp0 != bp)  SetPressure(bp0);
        if (n>1){
            std::vector<int> kvec(n+1,0); 
            kvec[0]=i;
            kvec[n]=j;

            int p = 1; //Used to increment all of the indicies correctly, at the end of each loop.
            double stotal = 0.0;
            while (kvec[n]==j) {

                stotal += qikj(x,bp0,kvec,n);
                //Update indices correctly
                kvec[1]++;
                while(kvec[p]==nsize && p<n) {
                    kvec[p]=0;
                    kvec[++p]++; //increase p by 1, and increase the next (p+1)th index
                    if(kvec.at(p)!=nsize)
                        p=1;
                }
            }
         return pow(A2*exp(bp), n)*stotal;
        }
        else {
            std::vector<int> kvec{i,j};
            return pow(A2*exp(bp), n)*qikj(x,bp0,kvec,1);
        }
      return 0.0;

    }

    double qikj(double x, double bp0, std::vector<int> kvec, int n){
        if (bp0 != bp)  SetPressure(bp0);
        if (kvec.size() != n+1){
            std::cout << "Wrong size of kvec! Something went wrong with the input parameters \n";
            return 0.0;
        }
        double sa = x;
        //for (int h=0; h<kvec.size()-1; h++){
        //    sa -= am(kvec.at(h),kvec.at(h+1));
        //}

        int h=0;
        while (sa >=0.0 && h<kvec.size()-1){
            sa -= am(kvec.at(h),kvec.at(h+1));
            h++;
        }

        if (sa >= 0){
            return exp(-bp0*x)*pow(sa, n-1.0)/tgamma(n); //tgamma(n+1)=factorial(n)
        }
        else{
            return 0.0;
        }
    }

    
};

#endif
