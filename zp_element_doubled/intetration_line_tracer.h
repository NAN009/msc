#ifndef intetration_line_tracer_h_
#define intetration_line_tracer_h_
#include <utility>
#include <vector>
#include <iostream>
#include "mscomplex.h"
namespace msc2d
{
	class CriticalPoint;
	class ILTracer
	{
		class WEdge
		{
		public:
			//存放有所有最大值的序数
			std::vector<std::pair<size_t, size_t>> max_ranges;			
			std::vector<std::pair<size_t, size_t>> min_ranges;
		};
	public:
		ILTracer(MSComplex2D&);
		~ILTracer();
		bool traceIntegrationLine();
		bool traceAscendingPath();
		int Round(double r);
	private:
		bool createWEdge();
		

		pair<int, int> getGradDirection(pair<int, int> xy);
		std::vector<WEdge> wedge_vec;
		MSComplex2D& msc;
	};
}
#endif