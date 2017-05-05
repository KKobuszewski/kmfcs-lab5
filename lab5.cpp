#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <cmath>


#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


using namespace std::literals::complex_literals;


// Ev0 [eV]
// hbar^2/2m0 = hbar^2c^2/2m0c^2 = 197.3269 / 0.51099 * 2

#define HBAR2_OVER_2M0 ((double) 197.3269*197.3269 / (-0.51099 * 2))  // minus because of holes!

#define A  ((double) 5.56 )  // GaAs lattice constant in Angstroems

// k  = M_PI / A

#define Ev0  ((double) 0.8)  // band edge energy

// Luttinger constants
#define LUTTINGER_GAMMA_1 ((double) 6.98)
#define LUTTINGER_GAMMA_2 ((double) 2.06)
#define LUTTINGER_GAMMA_3 ((double) 2.93)



inline double Pk(double kx, double ky, double kz)
{
    return HBAR2_OVER_2M0 * LUTTINGER_GAMMA_1 * ( kx*kx + ky*ky + kz*kz ) - Ev0;
}

inline double Qk(double kx, double ky, double kz)
{
    return HBAR2_OVER_2M0 * LUTTINGER_GAMMA_2 * ( kx*kx + ky*ky - 2.0 * kz*kz );
}

inline std::complex<double> Rk(double kx, double ky, double kz)
{
    return -sqrt(3.0) * HBAR2_OVER_2M0 * LUTTINGER_GAMMA_2 * ( kx*kx - ky*ky )  +  1i * 2.0*sqrt(3) * HBAR2_OVER_2M0 * LUTTINGER_GAMMA_3 * kx*ky;
}

inline std::complex<double> Sk(double kx, double ky, double kz)
{
    return  2.0*sqrt(3) * HBAR2_OVER_2M0 * LUTTINGER_GAMMA_3 * (kx - 1i *ky)*kz;
}



void create_hamiltonian(Eigen::MatrixXcd &H, const double kx, const double ky, const double kz)
{
    H(0,0) =  Pk(kx,ky,kz) + Qk(kx,ky,kz);
    H(1,0) = -Sk(kx,ky,kz);
    H(2,0) =  Rk(kx,ky,kz);
    H(3,0) =  0.0;
    
    H(1,1) =  Pk(kx,ky,kz) - Qk(kx,ky,kz);
    H(2,1) =  0.0;
    H(3,1) =  Rk(kx,ky,kz);
    
    H(2,2) =  Pk(kx,ky,kz) - Qk(kx,ky,kz);
    H(3,2) =  Sk(kx,ky,kz);
    
    H(3,3) =  Pk(kx,ky,kz) + Qk(kx,ky,kz);
    
    H(0,1) = std::conj(H(1,0));
    H(0,2) = std::conj(H(2,0));
    H(0,3) = std::conj(H(3,0));
    
    H(1,2) = std::conj(H(2,1));
    H(1,3) = std::conj(H(3,1));
    
    H(2,3) = std::conj(H(3,2));
}


void iterate_over_path(Eigen::MatrixXcd &H, Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>& eigensolver, std::ofstream &file,
                       const double kx_start, const double ky_start, const double kz_start,
                       const double kx_end, const double ky_end, const double kz_end,
                       const unsigned num_points = 100)
{
    double kx,ky,kz; 
    double dkx, dky, dkz;
    dkx = 2.0*M_PI * (kx_end - kx_start)/((double) num_points)/A;
    dky = 2.0*M_PI * (ky_end - ky_start)/((double) num_points)/A;
    dkz = 2.0*M_PI * (kz_end - kz_start)/((double) num_points)/A;
    
    kx = kx_start;
    ky = ky_start;
    kz = kz_start;
    
    for (unsigned ix = 0; ix <= num_points; ix++)
    {
        kx += dkx;
        ky += dky;
        kz += dkz;
        
        create_hamiltonian(H,kx,ky,kz);
        
        eigensolver.compute(H);
        
        std::cout << std::setprecision(2) << std::fixed;
        std::cout << "(" << kx << "," << ky << "," << kz << ")" << "\t";
        std::cout << eigensolver.eigenvalues().transpose() << std::endl;
        file      << std::setprecision(15) << std::fixed;
        file      << kx << "\t" << ky << "\t" << kz << "\t";
        file      << eigensolver.eigenvalues().transpose()/1e04 << std::endl;
        
    }
}



int main(int argc, char* argv[])
{
    
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
    Eigen::MatrixXcd H(4,4);
    //Eigen::MatrixXcd S(2,2);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver;
    
    
    
    std::ofstream file;
    file.open("k-p.dat");
    
    
    iterate_over_path(H,eigensolver,file,0.0,0.0,0.0,0.0,0.0,0.2);
    
    
    // path M-Gamma
    std::cout << "M-Gamma" <<std::endl;
    
    // path Gamma-K
    std::cout << "Gamma-K" <<std::endl;
    
    
    // path K-M
    std::cout << "K-M" <<std::endl;
    
    
    
    
    file.close();
    
    return EXIT_SUCCESS;
}