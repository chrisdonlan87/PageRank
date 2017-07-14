#include <iostream>
#include <fstream>
#include<sstream>
#include <list>
#include <array>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;
double round(double number,int digits){
	double round_shift = pow(10.,(double)digits);
	return floor(number*round_shift + 0.5)/round_shift;
}
template<typename Out>
void split(const std::string &s, char delim, Out result){
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while(std::getline(ss,item,delim)){
		*(result++) = item;
	}
}
double errorCheck(string input_file_string,vector<double> * page_rank){
	if (input_file_string == "small.txt"){
		ifstream soutput("small-output.txt");
		int u = 0;double max_error = 0.;
		string l;
		while(getline(soutput,l)) {
			if(fabs(stod(l) - (*page_rank)[u]) > max_error)
				max_error = fabs(stod(l)-(*page_rank)[u]);
			u++;
		}
		soutput.close();
		return max_error;
	}
	return 0.;
}
void currentRank(std::vector<int>*rank2Vertex,std::vector<int>*vertex2Rank,std::vector<double>*page_rank_vector){
#pragma omp parallel for
	for (int i = 0; i < (*rank2Vertex).size(); i++)
		(*rank2Vertex)[i] = i;
	std::sort(rank2Vertex->begin(),rank2Vertex->end(),[page_rank_vector](int i,int j){
		return ((*page_rank_vector)[i] > (*page_rank_vector)[j]);
	});
#pragma omp parallel for
	for (int i = 0; i < rank2Vertex->size();i++){
		(*vertex2Rank)[(*rank2Vertex)[i]] = i;
	}
}

void runPageRankParallel(double convergence_req, double alpha, vector<list<int>>*Gf,vector<list<int>>*Gr,vector<double>*page_rank){
	unsigned int max_u = page_rank->size();
	vector<double> prior_page_rank(page_rank->size(), 1. / (double)page_rank->size());

	double absolute_diff = 1.0;
	double checksum = 0.;
	int iterations = 0;
	double pr_value,pr_l,pr_r;
	int u,v;
	int x = 1;
	while (absolute_diff > convergence_req) {
		for (u = 0; u < max_u; u++) {
			pr_l = (1. - alpha) / (max_u);
			pr_r = 0.;
			for (int v:(*Gr)[u])
				pr_r += prior_page_rank[v] / (double) (*Gf)[v].size();
			(*page_rank)[u] = pr_l + alpha*pr_r;
		}
		absolute_diff = 0.;
#pragma omp parallel for
		for (u = 0; u < max_u; u++) {
			if (fabs((*page_rank)[u] - prior_page_rank[u]) > absolute_diff)
				absolute_diff = fabs((*page_rank)[u] - prior_page_rank[u]);
			prior_page_rank[u] = (*page_rank)[u];
		}
		if (iterations % 10 == 0)
			cout << iterations << " " << absolute_diff << endl;
		iterations++;
	}
}
void runPageRank(double convergence_req, double alpha, vector<list<int>>*Gf,vector<list<int>>*Gr,vector<double>*page_rank){
	unsigned int max_u = page_rank->size();
	vector<double> prior_page_rank(page_rank->size(), 1. / (double)page_rank->size());
	for (int i = 0; i < max_u; i++)
		(*page_rank)[i] = 1./(double)max_u;

	int gf_size = Gf->size();
	int gr_size = Gr->size();
	for(int i = 0; i < gf_size;i++){
		for ( int v : (*Gf)[i]){
			int k = v;
		}
		for (int v : (*Gr)[i])
			int x = v;
	}


	double absolute_diff = 1.0;
	int iterations = 0;
	double pr_value,pr_l,pr_r;
	int u,v;

	while (absolute_diff > convergence_req) {
		int x = 1;
		for (u = 0; u < max_u; u++) {
			if (u == 800)
				x = 2;
			pr_l = (1. - alpha) / (max_u);
			pr_r = 0.;
			for (int v:(*Gr)[u])
				pr_r += prior_page_rank[v] / (double) (*Gf)[v].size();
			(*page_rank)[u] = pr_l + alpha*pr_r;
		}
		absolute_diff = 0.;
		for (u = 0; u < max_u; u++) {
			if (fabs((*page_rank)[u] - prior_page_rank[u]) > absolute_diff)
				absolute_diff = fabs((*page_rank)[u] - prior_page_rank[u]);
			prior_page_rank[u] = (*page_rank)[u];
		}
		//cout << iterations << " " << absolute_diff << endl;
		iterations++;
	}
}

void vary_alpha(double low,double high,double step,double convergence_requirement,vector<list<int>>*Gi,vector<list<int>>*Gr){
	for (double i = low; i < high; i+= step){
		vector<double>page_rank(Gi->size());
		runPageRank(convergence_requirement,i,Gi,Gr,&page_rank);

		string alpha_string;
		stringstream astream;
		astream << (int)(i*100.);
		alpha_string = "results" + astream.str() + ".txt";
		ofstream alpha_file(alpha_string);
		for (int u = 0; u < Gi->size(); u++)
			alpha_file<<setprecision(10) << page_rank[u] << endl;
		alpha_file.close();
	}
}
list<int> my_binary_search(int min_,int max_,std::function<bool(int,int,list<int>*)>*subproblem,bool*search_success,bool verbose){
	list<int> neighbors_added;
	int pos = min_;
	int prev_pos = min_;
	int upper,lower;
	int step = 1;
	bool success = false;
	bool started = false;
	while(pos < max_){
		if(!started) {
			step = 1;
			started=true;
		} else step = step*2;

		prev_pos = pos;
		pos = min(pos+step,max_);
		if (verbose)
			std::cout<<"binary search phase 1 position = "<<pos<<std::endl;
		list<int>neighbors;
		success = (*subproblem)(min_,pos,&neighbors);
		if(success){
			if (verbose)
				std::cout<<"binary search phase 1 successful"<<std::endl;
			upper = pos;
			lower = prev_pos;
			break;
		}
	}
	if(!success){
		if (verbose)
			std::cout<<"binary search phase 1 failed"<<std::endl;
		*search_success = false;
		list<int> neighbors;
		return neighbors;
	}

	int split_point;
	list<int>neighbors;
	while (upper - lower > 1){
		split_point = ((upper-lower) / 2)+lower;
		if (verbose)
			std::cout<<"binary search phase 2 split point = "<<split_point<<std::endl;
		neighbors.clear();
		success = (*subproblem)(min_,split_point,&neighbors);
		if (success)
			upper = split_point;
		else
			lower = split_point;
	}
	*search_success = true;
	return neighbors;
};
void ex2_binSearch(int min_,int max_,std::function<bool(int,list<int>*)>*subproblem,bool*search_success,list<int> * neighbors_added,bool verbose){
	int pos = min_;
	int prev_pos = min_;
	int upper,lower;
	int step = 1;
	bool success = false;
	bool started = false;
	while(pos < max_){
		if(!started) {
			step = 1;
			started=true;
		} else step = step*2;

		prev_pos = pos;
		pos = min(pos+step,max_);
		if (verbose)
			std::cout<<"binary search phase 1 position = "<<pos<<std::endl;
		list<int>neighbors;
		success = (*subproblem)(pos,&neighbors);
		if(success){
			if (verbose)
				std::cout<<"binary search phase 1 successful"<<std::endl;
			upper = pos;
			lower = prev_pos;
			break;
		}
	}
	if(!success){
		if (verbose)
			std::cout<<"binary search phase 1 failed"<<std::endl;
		*search_success = false;
		return;
	}

	int split_point;
	while (upper - lower > 1){
		split_point = ((upper-lower) / 2)+lower;
		if (verbose)
			std::cout<<"binary search phase 2 split point = "<<split_point<<std::endl;
		neighbors_added->clear();
		success = (*subproblem)(split_point,neighbors_added);
		if (success)
			upper = split_point;
		else
			lower = split_point;
	}
	*search_success = true;
	return;
};

list<int> my_DP(int min,int max,int bin_search_max,int step,function<bool(int,int,list<int>*)>*subproblem,function<double(list<int>)>*cost,vector<int>*orig_vertex2Rank){
	list<int> L;
	double T= numeric_limits<double>::max();
	double c;
	for(int i = 0; i < max-min; i+=step){
		bool success;
		list<int> neighbors = my_binary_search(i+min,bin_search_max,subproblem,&success,false);
		if (success) {
			// min of cost of binary search starting from i, T[i-1]
			c = (*cost)(neighbors);
			if (c < T) {
				T = c;
				copy(neighbors.begin(), neighbors.end(), L.begin());
				ofstream ex1("ex1_results.txt",ios_base::app);
				ex1 << "cost: " << c << endl;
				cout << "cost: " << c << endl;
				ex1 << "vertex,rank" << endl;
				cout << "vertex,rank" << endl;
				for (int j : neighbors) {
					ex1 << j << "," << (*orig_vertex2Rank)[j] << endl;
					cout << j << "," << (*orig_vertex2Rank)[j] << endl;
				}
				cout<<endl << endl;
				ex1.close();
			}
		}
	}
	return L;
}


list<int> extraCredit1(vector<list<int>>*Gi,vector<list<int>>*Gr,double alpha,double convergence_req){
	vector<list<int>> Gi_1(Gi->size()+1); copy(Gi->begin(),Gi->end(),Gi_1.begin());
	vector<list<int>> Gf_1(Gi->size()+1);
	vector<list<int>> Gr_1(Gr->size()+1); copy(Gr->begin(),Gr->end(),Gr_1.begin());
	// node index!
	int node = Gi->size();
	Gi_1[Gi->size()].push_back((int) Gi->size());
	Gr_1[Gr->size()].push_back((int) Gi->size());
	vector<double> page_rank(Gi_1.size());
	vector<int>  rank2Vertex(Gi_1.size());
	vector<int>  vertex2Rank(Gi_1.size());
	vector<int> orig_vertex2Rank(Gi_1.size());
	vector<int> orig_rank2Vertex(Gi_1.size());
	list<int> neighbors_added;
	int top_one_percent = (int)((double)Gi_1.size() * 0.01);
	// page rank subproblem
	function<bool(int,int,list<int>*)> PageRank_SubProblem = [&page_rank,node,&Gi_1,&Gf_1,&Gr_1,top_one_percent,&vertex2Rank,&orig_vertex2Rank,&rank2Vertex,&orig_rank2Vertex,alpha,convergence_req]
			(int ranks_lower,int ranks_upper,list<int>*neighbors) {
		copy(Gi_1.begin(),Gi_1.end(),Gf_1.begin());
		Gr_1.~vector();
		new(&Gr_1) vector<list<int>>(Gi_1.size());
#pragma omp parallel for
		for(int u = 0; u < Gi_1.size(); u++){
			for(int v : Gi_1[u])
				Gr_1[v].push_back(u);
		}
		neighbors->clear();
#pragma omp parallel for
		for(int i = ranks_lower; i < ranks_upper + 1; i++){
			int e = orig_rank2Vertex[orig_rank2Vertex.size()-1 - i];
			if (e != node){
				Gf_1[e].push_back(node);
				Gr_1[node].push_back(e);

				neighbors->push_front(e);
			}
		}
		runPageRank(convergence_req,alpha,&Gf_1,&Gr_1,&page_rank);
		currentRank(&rank2Vertex,&vertex2Rank,&page_rank);
		int node_rank = vertex2Rank[node];
		return (vertex2Rank[node] <= top_one_percent);
	};
	function<double(list<int>)> Cost_SubProblem = [&orig_vertex2Rank](list<int> neighbors){
		double cost = 0.; int rank;
		for(int v : neighbors) {
			rank = orig_vertex2Rank[v];
			cost += pow(876000 - rank, 2.);
		}
		return cost;
	};



	runPageRank(convergence_req,alpha,&Gi_1,&Gr_1,&page_rank);
	currentRank(&rank2Vertex,&vertex2Rank,&page_rank);
	copy(rank2Vertex.begin(),rank2Vertex.end(),orig_rank2Vertex.begin()); // so as to accurately add websites
	copy(vertex2Rank.begin(),vertex2Rank.end(),orig_vertex2Rank.begin()); // so as to accurately calculate the cost

	bool search_success;
	neighbors_added = my_DP(0,(int)Gi_1.size(),(int)Gi_1.size(),10,&PageRank_SubProblem,&Cost_SubProblem,&orig_vertex2Rank);
}
void extraCredit2(vector<list<int>>*Gi,vector<list<int>>*Gr,double alpha,double convergence_req,list<int>*final_neighbors){
	unsigned long original_size = Gi->size()+1;
	int top_one_percent = (int) ((double)Gi->size() *0.01);
	vector<list<int>> Gi_1(Gi->size()+1);
	vector<list<int>> Gir_1(Gi->size()+1);
	for(int i = 0; i < Gi->size(); i++) {
		for (int u : (*Gi)[i])
			Gi_1[i].push_back(u);
		for (int u : (*Gr)[i])
			Gir_1[i].push_back(u);
	}
	int x = 0;

	// node index!
	int node = (int) Gi->size();
	Gi_1[Gi->size()].push_back(node);
	Gir_1[Gi->size()].push_back(node);
	list<int>fake_sites;
	//region PageRank_SubProblem
	function<bool(int,list<int>*)> PageRank_SubProblem = [node,&Gir_1,&Gi_1,top_one_percent,alpha,convergence_req,&fake_sites]

			(int fake_site_count,list<int>*neighbors) {

		vector<double> page_rank1(Gi_1.size());
		vector<int>  rank2Vertex1(Gi_1.size());
		vector<int>  vertex2Rank1(Gi_1.size());

		vector<list<int>> Gf_1(Gi_1.size()); copy(Gi_1.begin(),Gi_1.end(),Gf_1.begin());
		vector<list<int>> Gr_1(Gir_1.size()); copy(Gir_1.begin(),Gir_1.end(),Gr_1.begin());
		int xx = 0;
		for (int i = 0; i < Gf_1.size(); i++)
			for(int u : Gf_1[i])
				int v = u;
		xx = 1;
		for (int i = 0; i < Gr_1.size(); i++)
			for(int u : Gr_1[i])
				int v = u;
		xx = 2;

		int sz = (int) Gi_1.size();
		fake_sites.clear();
		for(int i = 0; i < fake_site_count; i++){
			fake_sites.push_back(sz+i);
			list<int>fedge_list;fedge_list.push_back(sz+i);Gf_1.push_back(fedge_list);
			list<int>redge_list;redge_list.push_back(sz+i);Gr_1.push_back(redge_list);

			int sz2 = Gf_1.size();
			vertex2Rank1.push_back(sz+i);
			rank2Vertex1.push_back(sz+i);
			page_rank1.push_back(1./(double)(sz+i));
		}
		neighbors->clear();
		for(int u : fake_sites){
			Gf_1[u].push_back(node);
			Gr_1[node].push_back(u);
			neighbors->push_front(u);
		}

		runPageRankParallel(convergence_req,alpha,&Gf_1,&Gr_1,&page_rank1);
		currentRank(&rank2Vertex1,&vertex2Rank1,&page_rank1);
		int node_rank = vertex2Rank1[node];
		return (vertex2Rank1[node] <= top_one_percent);
	};
	//endregion
	//region costSubProblem
	function<double(list<int>)> costSubProblem = [&original_size](list<int> neighbors){
		double cost = 0.; int rank;
		for(int v : neighbors) {
			cost += 1000;
		}
		return cost;
	};
	//endregion
	bool search_success;
	cout << "starting bin search" << endl;
	ex2_binSearch(0,1000,&PageRank_SubProblem,&search_success,final_neighbors,true);
	double c = costSubProblem(*final_neighbors);

	ofstream ex2_out("ex2_results.txt");
	ex2_out << "Cost: " << c << endl;
	for (int i : *final_neighbors)
		ex2_out << i << endl;
	ex2_out.close();
}
int main(int argc, char *argv[]) {
	double alpha = atof(argv[2]);
	string input_file_str = argv[1];
	fstream ifile(input_file_str);
	double convergence_requirement = 1e-10;
	string l, ustr;
	int u,v; unsigned long max_u = 0;
	//region Get max_u
	while (getline(ifile, l)) {
		std::stringstream ls(l);
		getline(ls, ustr, ':');
		u = std::stoi(ustr);
		if (u > max_u) max_u = (unsigned long) u;
		while (getline(ls,ustr,',')){
			u = std::stoi(ustr);
			if (u > max_u) max_u = (unsigned long) u;
		}
	}
	max_u += 1; // Add one for basing at zero!
	ifile.clear();
	ifile.seekg(0,ios::beg);
	//endregion
	vector<double> page_rank(max_u, 1. / max_u);
	//region Read file into Gi, Gr, create self edges
	vector<list<int>> Gi(max_u),Gr(max_u);
	while(getline(ifile,l)){
		stringstream ls(l);
		getline(ls,ustr,':');
		u = std::stoi(ustr);
		while(getline(ls,ustr,',')){
			Gi[u].push_back(stoi(ustr));
			Gr[stoi(ustr)].push_back(u);
		}
	}
	ifile.close();
	for (u = 0; u < max_u; u++){
		Gi[u].push_back(u);
		Gr[u].push_back(u);
	}
	//endregion
//	cout<<"starting first page rank calculation"<<endl;
//	time_t timer;
//	time(&timer);
	runPageRank(convergence_requirement,alpha,&Gi,&Gr,&page_rank);
//	double deltat = difftime(time(NULL),timer);
//	cout<<"completed first page rank calculation, time: " << deltat<<endl;


//	double max_error = errorCheck(input_file_str,&page_rank);
	//vary_alpha(0.75,0.96,0.01,convergence_requirement,&Gi,&Gr);
//	list<int> top1 = extraCredit1(&Gi,&Gr,alpha,convergence_requirement);
//	list<int> top2;
//	cout<<"starting extra credit 2"<<endl;
//	extraCredit2(&Gi,&Gr,alpha,convergence_requirement,&top2);
//	cout<< "Number of fake sites: " << top2.size() << endl;

	// Output results
	for (u = 0; u < Gi.size();u++)
		cout << setprecision(11) << page_rank[u] << endl;
    return 0;
}
