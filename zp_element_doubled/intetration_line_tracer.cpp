#include "intetration_line_tracer.h"
#include "critical_point_finder.h"
#include <iostream>
#include <vector>
//#include "vtkReader.h"
#include "mscomplex.h"
#include<memory>
const int vr1_size = 201;
//shared_ptr<vtkReader> vr; 
//vtkReader vr1;


using namespace std;
namespace msc2d
{
	ILTracer::ILTracer(MSComplex2D& _msc) :msc(_msc){}
	ILTracer::~ILTracer(){}

	/*bool ILTracer::createWEdge()
	{
		wedge_vec.clear();
		wedge_vec.resize(vr1_size);
		for (size_t vid = 0; vid < vr1_size; ++vid)
		{
			WEdge& we = wedge_vec[vid];
			const vertHandleArray &adj_vertices=
		}
	}*/
	
	bool ILTracer::traceAscendingPath()
	{
		cout << "Trace ascending path" << endl;
		for (vector<CriticalPoint>::iterator it = msc.cp_vec.begin(); it != msc.cp_vec.end(); ++it)
		{
			if (it->type == SADDLE)
			{
				CriticalPoint& sad = *it;
				
				int prev_vid = -1, curr_x = sad.xy_local.first,curr_y=sad.xy_local.second,curr_vid=sad.meshIndex;
					
					msc.il_vec.push_back(IntegrationLine());
					IntegrationLine &il = msc.il_vec[msc.il_vec.size() - 1];
					PATH& mesh_path = il.path;
					mesh_path.push_back(make_pair(curr_x,curr_y));

					int k1x = floor(sad.eig_vector1.first);
					int k2x = floor(msc.cp_vec[(curr_x+k1x)*vr1_size+curr_y+1].dif.first);
					int k3x = floor(msc.cp_vec[(curr_x + k2x)*vr1_size + curr_y + 1].dif.first);
					//int next_x = curr_x + floor((sad.eig_vector1.first + 2 * msc.cp_vec[(curr_x + k1x)*vr1_size + curr_y + 1].dif.first + 2 * msc.cp_vec[(curr_x + k2x)*vr1_size + curr_y + 1].dif.first + msc.cp_vec[(curr_x + k3x)*vr1_size + curr_y + 2].dif.first) / 3);
					int k1y = floor(sad.eig_vector1.second);
					int k2y = floor(msc.cp_vec[curr_x + 1 + (curr_y + k1y)*vr1_size].dif.second);
					int k3y = floor(msc.cp_vec[curr_x + 1 + (curr_y + k2y)*vr1_size].dif.second);

					int next_x = curr_x + floor((sad.eig_vector1.first + 2 * (msc.cp_vec[(curr_x + k1x)*vr1_size + curr_y + 1].dif.first) + 2 * (msc.cp_vec[(curr_x + k2x)*vr1_size + curr_y + 1].dif.first) + msc.cp_vec[(curr_x + k3x)*vr1_size + curr_y + 2].dif.first) / 3);
					int next_y = curr_y + floor((sad.eig_vector1.second + 2 * (msc.cp_vec[curr_x + 1 + (curr_y + k1y)*vr1_size].dif.second) + 2 * (msc.cp_vec[curr_x + 1 + (curr_y + k2y)*vr1_size].dif.second) + msc.cp_vec[curr_x + 2 + (curr_y + k3y)*vr1_size].dif.second) / 3);

					curr_vid = next_x*vr1_size + next_y;
					while (msc.cp_vec[curr_vid].type != MAXIMAL)
					{
						if (msc.cp_vec[curr_vid].type == SADDLE)
						{
							continue;
						}
						pair<int, int> tmp_xy;
						tmp_xy = getGradDirection(msc.cp_vec[curr_vid].xy_local);
						mesh_path.push_back(tmp_xy);                  
						
						prev_vid = curr_vid;
						curr_vid = tmp_xy.first*vr1_size + tmp_xy.second;
					}
					il.startIndex = it->xy_local;
					il.endIndex = mesh_path[mesh_path.size()-1];
					
				}
			}

		
		return true;
	}
	
	pair<int,int> ILTracer::getGradDirection(pair<int,int> xy)
	{
		//四阶龙格库塔法 h=2,待检验，floor向下取整
		int x = xy.first, y = xy.second;
		int k1x = floor(msc.cp_vec[x*vr1_size+y].dif.first);
		int k2x=floor(msc.cp_vec[(x + k1x)*vr1_size + y + 1].dif.first);
		int k3x = floor(msc.cp_vec[(x + k2x)*vr1_size + y + 1].dif.first);
		int next_X = x + floor((msc.cp_vec[x*vr1_size + y].dif.first + 2 * msc.cp_vec[(x + k1x)*vr1_size + y + 1].dif.first + 2 * msc.cp_vec[(x + k2x)*vr1_size + y + 1].dif.first + msc.cp_vec[(x + k3x)*vr1_size + y + 2].dif.first) / 3);
		int k1y = floor(msc.cp_vec[y*vr1_size + x].dif.second);
		int k2y = floor(msc.cp_vec[(y + k1y)*vr1_size + x + 1].dif.second);
		int k3y = floor(msc.cp_vec[(y + k2y)*vr1_size + x + 1].dif.second);
		int next_Y = y + floor((msc.cp_vec[y*vr1_size + x].dif.second + 2 * msc.cp_vec[(y + k1y)*vr1_size + x + 1].dif.second + 2 * msc.cp_vec[(y + k2y)*vr1_size + x + 1].dif.second + msc.cp_vec[(y + k3y)*vr1_size + x +2].dif.second) / 3);

		return make_pair(next_X, next_Y); 

	}
}

