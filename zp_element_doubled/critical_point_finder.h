#ifndef critical_point_finder_h
#define critical_point_finder_h

#include "mscomplex.h"
namespace msc2d
{
	class CPFinder
	{
	public:
		CPFinder(MSComplex2D& _msc);
		~CPFinder();
		void findCriticalPoints();
		bool JacbiCor(double * pMatrix, int nDim, double *pdblVects, double *pdbEigenValues, double dbEps, int nJt);
	private:
		CriticalPointType getPointType(double eig_value1,double eig_value2)const;
		MSComplex2D &msc;
	};

	class reconstructor
	{
	public:
		reconstructor(char*, char*);
		~reconstructor(void);
		void set_dimension(int, int);
		void set_step(double, double);
		void initialize();
		void read_data();
		void reconstruct();

		bool with_header;
	private:
		char* original_data_path;
		char* output_data_path;
		//point_evaluator *pe;
		//global_functions *gf;
		int L;
		int W;
		double stepx;
		double stepy;
		int nx;
		int ny;

		std::vector<double> x_ordinates;
		std::vector<double> y_ordinates;

		std::vector<std::vector<std::vector<double> > > data;
		//	std::vector<std::vector<double> > reconstructed_piece;


	};
		
}
#endif