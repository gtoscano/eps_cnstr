#include <iostream>
#include <random>
#include <vector>
#include <float.h>

/*
 * A structure for individuals
 */

extern int ngenes;
extern int ngenes_selector;
extern int L;
extern std::string uuid;
extern std::vector<int> ubounds;
extern int pollutant_id;
extern double Theta;
extern double eta_c;
extern double eta_m;

struct IND{
	std::vector<int> chrom_selector;
	std::vector<double> chrom;
	double constr;
    double fitness;

    IND();
	IND(const IND& ind );
	IND& operator=(const IND& ind);
	//int decode(int idx);
    //void encode(std::vector<int> &x_vars, std::vector<int> &bin_str);
    double get_chrom(int idx);

    void insert_solution(std::vector<double> &my_final_x, std::vector<int> &my_final_selector);
	friend std::ostream& operator<<(std::ostream &os, const IND& ind) {
		for (auto& i : ind.chrom) {
			os << i;
		}
		os << std::endl;
		return os;
	}

	void xover(IND &ind);
	void mutation();
};



//The genetic algorithm
class SGA {
	int pop_size = 100;
	int ngenerations = 100;
	std::vector<IND> pop;
	std::vector<IND> pop_new;
	IND best_ind;

public:
    bool activate_local_search=false;
	SGA();
	void evaluate();
	void evaluate_best_solution();
	void selection();
    void local_search(); 
	void xover_and_mutation();
	void elitism();
	void new_generation();
	void run();
    void get_best_individual(std::vector<double> &best_individual);
};
