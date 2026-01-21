#include <iostream>
#include <string> //getline
#include <iomanip>
#include <fstream>
#include <stdio.h> 
#include <unistd.h>
#include <cmath>
#include <pthread.h>
#include <omp.h>
#include <complex>
#include <vector>


#include "Dumbbells.hpp"
#include "iofunctions.hpp"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>

int main() 
{
    InputData indata;
    if(ReadInputFile("input.txt", &indata)){
        std::cout << "Exiting...\n" << std::endl;
        std::cout << "\n Press enter to exit \n" << std::endl;
        std::cin.get();
        return 1;
    }
    if(CheckData(&indata)){
        std::cout << "Exiting...\n" << std::endl;
        std::cout << "\n Press enter to exit \n" << std::endl;
        std::cin.get();
        return 1;
    }

    //Discretization parameter
    int nsize = indata.nsize;

    //nmax for shell calculation
    int nmax = 3;

    //Lo primero es ver si el usuario ha dado presion o temperatura
    bool bPressure = true;
    Dumbbells2 f = Dumbbells2(&indata);
    //std::cout << "Contact distance matrix\n"<< std::endl;
    //std::cout << f.am << std::endl;
    if (indata.density>0.0){
        f.SetDensity(indata.density);
        std::cout << std::fixed << std::setprecision(9) << "The longitudinal pressure of the system has been set to " << f.bp << std::endl;
        std::cout << "Therefore compressiblity factor  is Zx = " << f.bp/indata.density << std::endl;
        std::cout << "The order parameter is " << f.OrderParameter(f.bp) << std::endl;
        indata.bp = f.bp;
        bPressure = false;
    }
    else{
        f.SetPressure(indata.bp);
        double density = f.Density(indata.bp);
        std::cout << "The density of the system has been set to " << density << std::endl;
        std::cout << "The compressiblity factor is " << f.bp/density << std::endl;
        std::cout << "The order parameter is " << f.OrderParameter(f.bp) << std::endl;
        bPressure = true;
    }
    std::ofstream ofile, ofile2;

    
    if (indata.comValue == 2){      //Compute equation of state
        std::vector<double> pvalues (indata.npoints, 0.0); //vector storing pressure values in the interval
        std::vector<double> dvalues (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues
        std::vector<double> svalues (indata.npoints, 0.0); //vector storing order parameter for the corresponding pressures in pvalues
        std::vector<double> Zxvalues (indata.npoints, 0.0); //vector storing pressure values in the interval
        std::vector<double> corrvalues (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues
        std::vector<double> corrvalues2 (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues
        std::vector<double> l0values (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues
        std::vector<double> heightangle (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues
        std::vector<double> posangle (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues

        double kmax = 0.0;

        //Compute step size of the pressure interval
        double pstep = pow(indata.xmax/indata.xmin,1.0/(double)indata.npoints); //multiplicando

        //Compute density for each of the pressure values in the interval
        double xtemp = indata.xmin;
        for(int pidx=0;pidx<indata.npoints;pidx++){
            pvalues[pidx] = xtemp;
            f.SetPressure(pvalues[pidx]);
            dvalues[pidx] = f.density;
            svalues[pidx] = f.OrderParameter(pvalues[pidx]);
            Zxvalues[pidx] = f.ComputeZ(pvalues[pidx]);
            corrvalues[pidx] = f.corrlength1;
            corrvalues2[pidx] = f.corrlength2;
            posangle[pidx] = f.PosMaxAngle(pvalues[pidx]);
            heightangle[pidx] = f.HeightMaxAngle(pvalues[pidx]);
            l0values[pidx] = f.l0;
            //update value
            xtemp = xtemp*pstep;

            ProgressBar(double(pidx)/(double(indata.npoints)-1),70,&kmax);

        }
        std::cout << std::endl;
        std::cout << "Done!\n" << std::endl;

        //Write results
        std::string fname;
        fname = "output_EoS_nsize_"+transformDouble(f.nsize,1,0)+".dat";
        ofile.open(fname);
        ofile << "#Parameters of the system: chainsize = " << indata.chainsize << ", temperature = " << 1.0/f.beta << ", nsize = " << nsize << std::endl;
        ofile << "bp density Zx, S, CorrLength l0/l1, CorrLength l0/l2, l0, PosOfMaxAngle, HeightMaxAngle" << std::endl;

        for(int i=0;i<pvalues.size();i++){
                ofile << std::fixed << std::setprecision(13)<< pvalues[i] << "\t" << dvalues[i] << "\t" << Zxvalues[i] <<  "\t" << svalues[i]<<  "\t" << corrvalues[i]<<  "\t" << corrvalues2[i]<<"\t"<< l0values[i] <<"\t"<< posangle[i] << "\t" << heightangle[i]/(M_PI/(double)(nsize-1.0)) <<std::endl;                   
        }
        ofile.close();

    }
    else if (indata.comValue==3){   //Compute Gx

        //vectors to store results
        double xi;
        double rstep = (indata.xmax-indata.xmin)/(indata.npoints-1);
        std::vector<double> rvalues(indata.npoints,0.0);
        std::vector<double> Fvalues(indata.npoints,0.0);

        std::vector<double> Grtemp;
        std::vector<double> g2values(indata.npoints,0.0);
        std::vector<double> gvalues00(indata.npoints,0.0);
        std::vector<double> g2values00(indata.npoints,0.0);
        std::vector<double> gvaluespm(indata.npoints,0.0);
        std::vector<double> g2valuespm(indata.npoints,0.0);
        std::vector<double> gvaluespp(indata.npoints,0.0);
        std::vector<double> g2valuespp(indata.npoints,0.0);
    
        double kmax = 0.0;

        //Limit for analitical calculation, which does not depend on the potential, only on the hard-core size which is always the same
        double xlimmax = (nmax+1)*f.amin;
        std::cout << "Computing..." << std::endl;
        #pragma omp parallel for num_threads(8) schedule(dynamic) private(xi, Grtemp)
        for(int k=0;k<indata.npoints;k++){
            Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtt;
            Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qttm;
            Qtt.resize(f.nsize,f.nsize);
            Qttm.resize(f.nsize,f.nsize);
            xi = indata.xmin + k*rstep;
            rvalues[k] = xi;
            if(xi < xlimmax){
                Grtemp = f.GxAna(xi, f.bp, nmax);
                Fvalues[k] = Grtemp[0];
                g2values[k] = Grtemp[1];
                gvalues00[k] = Grtemp[2];
                g2values00[k] = Grtemp[3];
                gvaluespm[k] = Grtemp[4];
                g2valuespm[k] = Grtemp[5];
                gvaluespp[k] = Grtemp[6];
                g2valuespp[k] = Grtemp[7];
  
            }
            else{
                Grtemp = f.GxInv(xi, f.bp, 300,40,Qtt, Qttm);
                Fvalues[k] = Grtemp[0];
                g2values[k] = Grtemp[1];  
                gvalues00[k] = Grtemp[2];
                g2values00[k] = Grtemp[3];
                gvaluespm[k] = Grtemp[4];
                g2valuespm[k] = Grtemp[5];
                gvaluespp[k] = Grtemp[6];
                g2valuespp[k] = Grtemp[7];
            } 

            //Progress bar
            #pragma omp critical
            ProgressBar(double(k)/(double(indata.npoints)-1),70,&kmax);

        }
        std::cout << std::endl;
        std::cout << "Done!\n" << std::endl;

        //Write results to file
        double S = f.OrderParameter(f.bp);
        std::cout << "******************\nWriting results to file " << "output_Gx.dat" << std::endl;
        std::ofstream ofile;
        std::string fname;
        if(bPressure){
            fname = "output_Gx_bp_"+transformDouble(f.bp,1,2)+"_m_"+transformDouble(f.nsize,1,0)+".dat";
        }
        else{
            fname = "output_Gx_d_"+transformDouble(indata.density,1,2)+"_m_"+transformDouble(f.nsize,1,0)+".dat";
        }
        ofile.open(fname);
        ofile << "#Parameters of the system: chainsize = " << indata.chainsize << ", temperature = " << 1.0/f.beta << ", nsize = " << nsize << ", density = " << f.density <<std::endl;
        ofile << "x, g(x), g_2(x), g_cos(x), g_sin(), g_00(x), g_+-(x), g_++(x), g_+0(x)" << std::endl;
        for(int k=0;k<(int)indata.npoints;k++){
            ofile << std::fixed<< std::setprecision(15) << rvalues[k] << "\t" << Fvalues[k] << "\t" << g2values[k]<< "\t" << gvalues00[k]<< "\t" << g2values00[k]<< "\t"<< gvaluespm[k]<< "\t" << g2valuespm[k] << "\t"<< gvaluespp[k]<< "\t" << g2valuespp[k] << "\n";
        }
        ofile.close();
        std::cout << "******************" << std::endl;
    }
    else if (indata.comValue==1){   //Compute density profile
    //Export density profile
        std::string fname;
        if(bPressure){
            fname = "output_DensityProfile_bp_"+transformDouble(f.bp,1,2)+"_m_"+transformDouble(f.nsize,1,0)+".dat";
        }
        else{
            fname = "output_DensityProfile_d_"+transformDouble(f.density,1,2)+"_m_"+transformDouble(f.nsize,1,0)+".dat";
        }
        ofile.open(fname);
        ofile << "#Parameters of the system: chainsize = " << indata.chainsize << ", temperature = " << 1.0/f.beta << ", pressure = " << f.bp << ", density = " << f.density << ", nsize = " << nsize << std::endl;
        ofile << "l0 = " << std::fixed << std::setprecision(15) << f.l0*(M_PI/(double)(nsize-1.0)) << std::endl;
        ofile << "corrlength_1 = " << std::fixed << std::setprecision(15) << f.corrlength1 << std::endl;
        ofile << "corrlength_2 = " << std::fixed << std::setprecision(15) << f.corrlength2 << std::endl;
        ofile << "y, phi^2, phi" << std::endl;
        std::cout << "Limits are fmin = " << f.fmin << " and fmax = " << f.fmax << std::endl;
        double gap = M_PI;// f.fmax - f.fmin; //M_PI;
        double dnsize = (double)nsize;
        for(int i=0;i<f.xvalues.size();i++){
            ofile << f.phi(i) << "\t" << f.xvalues[i]/(gap/(dnsize-1.0)) << "\t" << sqrt(f.xvalues[i]/(M_PI/(double)(nsize-1.0))) << std::endl;           
        }
        ofile.close();
    }
    else if (indata.comValue == 4){ //Compute the poles
        //General el vector de presiones
        std::vector<double> ps;
        ps.reserve(static_cast<size_t>(indata.npoints));
        const double r = std::pow(indata.xmax / indata.xmin, 1.0 / (indata.npoints - 1));
        for (int i = 0; i < indata.npoints; ++i) {
            ps.push_back(indata.xmin * std::pow(r, i));
        }
        ps.back() = indata.xmax;

        std::vector<std::vector<double>> points = { //m=109 todos los polos que creo que son importantes
            //{-4.95335, -9.64621e-10}, {-3.93956, 1.91701}// p=0.05
            {-0.0495858, 3.28066e-23}, {-0.0344255, 5.90309} //p=30
            //{-0.0236107, -2.09829e-23}, {-0.0210445, 5.99932} //p=40
            //{-9.53549e-05, -8.49571e-24}, {-0.00359955, 6.16453} //p=100
        };


        std::vector<double> funcvalues = {0,0};

        std::vector<std::vector<double>> pointsold = points;
        std::vector<std::vector<double>> seed = points;
        int progress = 0;
        double bpold = f.bp;
        double bp = f.bp;

        std::cout << f.dtm(points.at(0)) << std::endl;
        std::cout << f.dtm(points.at(1)) << std::endl;


        //Abrimos el archivo
        std::string fname;
        if(bPressure){
            fname = "output_Poles_m_"+transformDouble(f.nsize,1,0)+".dat";
        }
        else{
            fname = "output_Poles_m_"+transformDouble(f.nsize,1,0)+".dat";
        }
        ofile.open(fname);
        ofile << "#Parameters of the system: chainsize = " << indata.chainsize << ", temperature = " << 1.0/f.beta << ", pressure = " << f.bp << ", density = " << f.density << ", nsize = " << nsize << std::endl;
        ofile << "y, phi^2, phi" << std::endl;
            // Itera presiones

        // 0) seed = point
        //std::vector<double> seed = {std::real(point), std::imag(point)};
        for (double bpc : ps) {
            //std::cout << bpc << "\t";
            f.SetPressure(bpc);

            // Guarda resultados z para esta presión (uno por punto)
            std::vector<std::complex<double>> zs;
            zs.reserve(points.size());

            // Itera puntos
            for (int i=0; i<points.size();i++) {

                if (progress>1){
                    seed[i][0] = pointsold[i][0] + (points[i][0] - pointsold[i][0]) * (bpc - bpold) / (bp - bpold);
                    seed[i][1] = pointsold[i][1] + (points[i][1] - pointsold[i][1]) * (bpc - bpold) / (bp - bpold);
                }

                // 2) z = f.minimize(seed)
                std::complex<double> z = f.minimize(seed[i], &funcvalues[i]);

                // 3) seed = (real(z), imag(z))  -> seed = z
                pointsold[i] = points[i];
                points[i] = {std::real(z), std::imag(z)};
                zs.push_back(z);
            }
        bpold = bp;
        bp=bpc;

        ofile << bpc << " " << f.density;
        for (int i=0;i<points.size();i++) {
            ofile << ' ' << std::real(zs[i]) << ' ' << std::imag(zs[i]);
            std::cout << bpc << ' ' << std::real(zs[i]) << ' ' << std::imag(zs[i]) << ' ' << funcvalues[i];
        }
        ofile << '\n';
        std::cout << "\n";
        progress++;
        }

    }
   
    
    return 0;
    
}
