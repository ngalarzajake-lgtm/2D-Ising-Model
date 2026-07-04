#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <fstream>
#include <string>




class ising {
public: 
int N;
int spins;

std::vector<std::vector<int>> spin = std::vector<std::vector<int>>(N, std::vector<int>(N));
std::uniform_int_distribution<int> random_site{0, N - 1};
std::bernoulli_distribution coin{0.5};
std::uniform_real_distribution<double> unif_dis{0.0, 1.0};
std::random_device rd;
std::mt19937 rng;

 ising(int lattice_size)
        : N(lattice_size),
          spins(lattice_size * lattice_size),
          spin(lattice_size, std::vector<int>(lattice_size)),
          random_site(0, lattice_size - 1),
          coin(0.5),
          unif_dis(0.0, 1.0),
          rng(rd())

    {}



    void initial_lattice_spin(int i, int j){
        for (j = 0; j<N; j++){
            for (i = 0; i<N; i++) {
                if (coin(rng) == 1){
                spin[i][j] = 1;
                }
                else{
                spin[i][j] = -1;
                } 
            }            
        }

    }

    double total_energy_i (int i, int j){
        double E=0.0;
        double dE;
        for (j=0; j<N; j++){
            for (i=0; i<N; i++){
                int up = spin[i][(j+1)%N];
                int down = spin[i][(j-1+N)%N];
                int right = spin[(i+1)%N][j];
                int left = spin[(i-1+N)%N][j];
                dE = -(spin[i][j] * (up+down+right+left))/2;
                E=E+dE;
                
            }
        }
        return E;
    }


    bool metropolis_algo (int i, int j, double T){
        i = random_site(rng);
        j = random_site (rng);
        spin[i][j]*=-1;
        double dE;
        int up = spin[i][(j+1)%N];
        int down = spin[i][(j-1+N)%N];
        int right = spin[(i+1)%N][j];
        int left = spin[(i-1+N)%N][j];
        dE = -2*(spin[i][j] * (up+down+right+left));

        if (dE <= 0) {;
            return true;
        }
        else {
            double probability = std::exp(-dE / T); //insert different temperatures here
            if (unif_dis(rng) < probability) {
                ;
                return true;
            }
            else {
                spin[i][j] *= -1;
                return false;
            }
        }
    }
        
    double calculate_total_magnetization(){
        double M=0;
        for(int j=0; j<N; j++){
            for(int i=0; i<N; i++){
                M = M+spin[i][j];
            }
        }
        return M;
    }

  
        

};


    void simulation(int N){
        ising model(N); 
        
        std::ofstream data("ising_data_N" + std::to_string(N) + ".csv");   
        
        double Tmax=4;
        double Tmin=1; 
        double step_add= 0.01;

        int sweep = N*N;
        int equilibriation_sweep = 100000;
        int sample_sweep = 300000;
        int total_runs = (equilibriation_sweep + sample_sweep) * sweep;
        int run = 0;    
        

        //acceptance ratio
        double accepted = 0.0;
        double rejected = 0.0;

        //raw data
        double E = 0.0;
        double M = 0.0;

        //expectation values
    
    
        double E_2 = 0.0;
        double M_2 = 0.0;

    


        
    // Each column stores one temperature, and each row stores one measured quantity.
    int temperature_points = static_cast<int>(std::round((Tmax - Tmin) / step_add)) + 1;
    std::vector<std::vector<double>> statistics = std::vector<std::vector<double>>(8, std::vector<double>(temperature_points));


        for (int temp_index = 0; temp_index < temperature_points; temp_index++){
            double alpha = Tmin + temp_index * step_add;
            model.initial_lattice_spin(0, 0);
            int samples = 0;
            int last_percent = -1;
            for (int r = 0; r<total_runs; r++) {
                int current_sweep = (r+1) / sweep;
                    bool decision = model.metropolis_algo(0 , 0, alpha);
                    if (decision) {
                        accepted = accepted + 1;
                    }
                    else {
                        rejected = rejected + 1;
                    }

                    if ((r+1) % sweep == 0) {
                        if (r >= equilibriation_sweep*sweep) {
                            double Ef = model.total_energy_i(0, 0);
                            double Mf = model.calculate_total_magnetization();
                            E = E + Ef;
                            E_2 = E_2 + Ef * Ef;

                            M = M + std::abs(Mf);
                            M_2 = M_2 + Mf *Mf;

                            samples++;
                        
                        }
                        
                    }
                            
                            int percent = static_cast<int>(100.0 * current_sweep / (equilibriation_sweep + sample_sweep));

                            if (percent != last_percent || current_sweep == (equilibriation_sweep + sample_sweep)) {
                            std::cout << "\rAlpha: " << alpha << " | " << percent << "% / 100% ";
                            std::cout.flush();
                            last_percent = percent;
                            }
            
            }

            double norm = 1.0/(samples);
            std::cout << "\n";
            statistics[0][run] = E*norm; //<E>
            statistics[1][run] = E_2*norm; //<E^2>
            statistics[2][run] = M*norm; // <M>_abs
            statistics[3][run] = M_2*norm; // <M^2>
            statistics[4][run] = accepted/(rejected+accepted);
            statistics[5][run] = alpha;
            statistics[6][run] = ((E_2*norm)-(E*norm*E*norm))/(alpha*alpha*N*N); //heat capacity
            statistics[7][run] = ((M_2*norm)-(M*norm*M*norm))/(alpha*N*N); // magnetic susceptibility 

            accepted=0;
            rejected=0;
            run = run + 1;
            E=0.0;
            E_2=0.0;
            M=0.0;
            M_2=0.0;  

            
        }
        
        data << "alpha,E_avg,E_per_spin,M_avg_abs,M_per_spin,Cv_per_spin,chi_per_spin,acceptance_ratio\n";

        std::cout << "---------------------------------------------------------------------------------------\n";

        for (int d = 0; d < run; d++) {

            double alpha = statistics[5][d];
            double E_avg = statistics[0][d];
            double M_avg_abs = statistics[2][d];
            double Cv = statistics[6][d];
            double chi = statistics[7][d];
            double acceptance_ratio = statistics[4][d];

            data << alpha << ","
                << E_avg << ","
                << E_avg / (N * N) << ","
                << M_avg_abs << ","
                << M_avg_abs / (N * N) << ","
                << Cv << ","
                << chi << ","
                << acceptance_ratio << "\n";

            std::cout << "Temperature: " << alpha
                    << "  <E>: " << E_avg
                    << "  <E>/spin: " << E_avg / (N * N)
                    << "  <M>: " << M_avg_abs
                    << "  <M>/spin: " << M_avg_abs / (N * N)
                    << "  Cv: " << Cv
                    << "  chi: " << chi
                    << "  Acceptance: " << acceptance_ratio
                    << "\n";
        }

        data.close();

        std::cout << "---------------------------------------------------------------------------------------\n";

    
    }
int main() {
    simulation(20);
    simulation(40);
    return 0;
}



        
       
        



