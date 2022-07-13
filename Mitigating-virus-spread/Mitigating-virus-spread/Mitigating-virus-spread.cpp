// Mitigating-virus-spread.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include <fstream>
#include <time.h>
#include <sstream>
#include <string>
#include <Windows.h>
#include <iostream>
#include <tchar.h>
#include <strsafe.h>
#include "Header.h"



using namespace std;

#define file_network     "Data/SocialNetwork1.dat"					//name of the file containing the edges between nodes (format: nodeId nodeId)
#define file_community   "Data/SocialNetwork1_Community.dat"			//name of the file containing communities (format: nodeId communityId)
#define file_seed        "Data/SocialNetwork1_InitialInfected_"		//prefix of the filename containing initial infected set 
#define MAX_THREADS 30											//the number of SIR repeatation 
#define node_count 200											//the number of nodes in the network
#define seed_set_count 2										// the numnber of initial seed sets
#define alpha 3
#define q1 1
#define q2 1
int Population_size=50;
int g_max=50;
#define BUF_SIZE 255


class Node_list;

float **community_avg_weight,*original_infected_nodes;
int community_count=0,edge_count=0;
Node_list *communities;
int **community_links;
float **neighbours;



double combination(int n,int r){
	double a=1,b=1,x=n;
	if(r>n/2) r=n-r;
	while(x>(n-r))
		a*=x--;
	x=r;
	while(x>1)
		b*=x--;
	return a/b;
}
float calcuate_rho(int i,int j,Node_list *communities,float **limited_community_links,int community_count){
	float rho=0;
	int N=communities[i].size();
	int n=ceil(limited_community_links[0][i]);
	float alp=(float)alpha/communities[j].size();
	int rij=limited_community_links[i][j];
	if(n>0)
	{
		for(int f=1;f<=min(n,rij);f++)
		{
			float PF=0;
			for (int l=f;l<=min(n,rij);l++)
			{
				float AL=combination(n,l)*combination(N-n,rij-l)/combination(N,rij);
				float BLF=combination(l,f)*powf(alp,f)*powf(1-alp,l-f);
				PF+=AL*BLF;
			}
			rho+=(1-powf(1-community_avg_weight[j][3],f))*PF;
		}
	}
	return rho;
}
float calcuate_Pjtr(int i,int j,Node_list *communities,float **limited_community_links,int community_count){
	float p=0;
	int N=communities[i].size();
	int n=ceil(limited_community_links[0][i]);
	int rij=alpha;
	float wi=community_avg_weight[i][3];
	if(n>0)
	{
		for(int f=1;f<=min(n,rij);f++)
		{
				float AL=combination(n,f)*combination(N-n,rij-f)/combination(N,rij);
				p+=(1-powf(1-community_avg_weight[i][3],f))*AL;
		}
	}
	return p;
}
void Determine_passengers(Node_list *communities,float **community_links,int community_count,Node_list *passengers,float **temp_neighbours,float **comm_avg_w){
	Node_list *Guests = (Node_list *)malloc(sizeof(Node_list)*(community_count+1));
	Node_list *Passengs = (Node_list *)malloc(sizeof(Node_list)*(community_count+1));
	for (int i = 0; i <=community_count; i++)
	{
		Guests[i]=Node_list(node_count);
		Passengs[i]=Node_list(node_count);
	}
		
	int *Populations = (int *)malloc(sizeof(int)*(community_count + 1));
	Populations[0]=0;
	for(int i=1;i<=community_count;i++)
		Populations[i]=communities[i].size()+Populations[i-1];

	for (int i=1;i<=community_count;i++)
		for(int j=1;j<=community_count;j++)
			if(i!=j)
			{
				int cnt=0;
				while(cnt<community_links[i][j])
				{
					int rnd=(1+rand()%communities[i].size())+Populations[i-1];
					if(!Guests[j].contains(rnd))
					{
						cnt++;
						Guests[j].insert(rnd);
						if(!Passengs[i].contains(rnd))
							Passengs[i].insert(rnd);
						if(!passengers->contains(rnd))
							passengers->insert(rnd);
						
					}
				}
			}
	
	Node_list tmp_chance=Node_list(node_count);

	for(int j=1;j<=community_count;j++)
	{
		for (int i=1;i<=Guests[j].size();i++)
		{
			int g=Guests[j].get(i);
			int cnt=0,tmp_neigh_count=alpha;
			if(communities[j].size()-Passengs[j].size()<alpha)
				tmp_neigh_count=communities[j].size()-Passengs[j].size();
			tmp_chance.reset();
			for (int i2=1;i2<=communities[j].size();i2++)
			{
				int nd=communities[j].get(i2);
				if(!passengers->contains(nd))
				{
					tmp_chance.insert(nd);
				}
			}
			
			while((cnt<tmp_neigh_count)&&(tmp_chance.size()>0))
			{
				int rnd=(rand()%tmp_chance.size())+1;
				int node=(int)tmp_chance.get(rnd);
				if((temp_neighbours[g][node]==0)&&(!passengers->contains(node)))
				{
					cnt++;
					temp_neighbours[g][node]=comm_avg_w[j][3];
					temp_neighbours[node][g]=comm_avg_w[j][3];
					tmp_chance.remove(node);
				}
			}
		}
	}
	
	tmp_chance.release();
	for (int i = 0; i <=community_count; i++)
	{
		Guests[i].release();
		Passengs[i].release();
	}
	delete Guests;
	delete Passengs;
	delete []Populations;
}
void Determin_Potential_seeds(Node_list Seeds,float **neighbours,float **temp_neighbours,Pseed_list *Potential,Node_list Infetced)
{
	Node_list first_order=Node_list(node_count);
	for(int i=1;i<=Seeds.size();i++)
	{
		int seed=Seeds.get(i);
		for (int j=1;j<=node_count;j++)
			if(((neighbours[seed][j]>0))&&(!first_order.contains(j))&&(!Seeds.contains(j))&&(!Infetced.contains(j)))
			{
				first_order.insert(j);
			}
	}

	for(int i=1;i<=first_order.size();i++)
	{
		int j=first_order.get(i);
		float np=1;
		for(int i2=1;i2<=Seeds.size();i2++)
		{
			int seed=Seeds.get(i2);
			if((neighbours[seed][j]>0)||(temp_neighbours[seed][j]>0))
				np=np*(1-neighbours[seed][j]);

		}
		np=1-np;
		Pseed Ps=Pseed(j,np);
		Potential->insert(Ps);
	}
	
	first_order.release();
}
int VIKOR(float **objectives, int size){
	int objective_count=2;
	float *sum = (float*)malloc(sizeof(float)*(objective_count));
	
	//**************************************************************************calcualte weights*****
	float **p = (float **)malloc(sizeof(float)*(size));
	for (int i = 0; i < size; i++)
		p[i] = (float*)malloc(sizeof(float)*(objective_count));
	for (int j = 0; j < objective_count; j++){
		sum[j] = 0;
		for (int i = 0; i <size; i++)
			sum[j] += objectives[i][j];
	}
	for (int i = 0; i < size; i++)
		for (int j = 0; j <objective_count; j++)
			p[i][j] = (float)objectives[i][j] / sum[j];
	//*******calcualte W*****
	float *E = (float*)malloc(sizeof(float)*(objective_count));
	float *D = (float*)malloc(sizeof(float)*(objective_count));
	float *W = (float*)malloc(sizeof(float)*(objective_count));
	float sum_D = 0;
	for (int j = 0; j <objective_count; j++){
		E[j] = 0;
		for (int i = 0; i <size; i++)
		if (p[i][j] != 0)
			E[j] += p[i][j] * logf(p[i][j]);
		E[j] = -((float)1 / logf(size))*E[j];
		D[j] = 1 - E[j];
		sum_D += D[j];
	}
	for (int j = 0; j <objective_count; j++)
		W[j] = (D[j] / sum_D);
	W[0]=(q1*W[0])/(q1*W[0]+q2*W[1]);
	W[1]=(q2*W[1])/(q1*W[0]+q2*W[1]);

	//**********************calculate Min and Max***********************************************************
	float *fPositive = (float*)malloc(sizeof(float)*(objective_count));
	float *fNegative = (float*)malloc(sizeof(float)*(objective_count));
	for (int j = 0; j <objective_count; j++){
		float min = objectives[0][j];
		float max = objectives[0][j];
		int locmin = 0, licmax = 0;
		for (int i = 1; i <size; i++)
		{
			if (objectives[i][j] < min)
			{
				min = objectives[i][j];
				locmin = i;
			}
			if (objectives[i][j]>max)
			{
				max = objectives[i][j];
				licmax = i;
			}
		}
		if(j==1)
		{
			fPositive[j]=min;
			fNegative[j]=max;
		}
		else
		{
			fPositive[j]=max;
			fNegative[j]=min;
		}
	}
	//******************************************************************calculate S and R and Q***********
	float *S = (float*)malloc(sizeof(float)*(size));
	float *R = (float*)malloc(sizeof(float)*(size));
	
	float Sb=-10000,Ss=10000,Rb=-10000,Rs=10000;
	
	for (int i = 0; i < size; i++)
	{
		S[i]=0;
		R[i]=0;
		for (int j = 0; j <objective_count; j++)
		{
			float tmp=W[j]*(float)(fPositive[j]-objectives[i][j])/(fPositive[j]-fNegative[j]);
			S[i]+=tmp;
			if(tmp>R[i])
				R[i]=tmp;
		}
		if(S[i]<Ss)
			Ss=S[i];
		if(S[i]>Sb)
			Sb=S[i];
		if(R[i]<Rs)
			Rs=R[i];
		if(R[i]>Rb)
			Rb=R[i];
	}

	float **Q =(float **)malloc(sizeof(float)*(size));
	for (int i = 0; i < size; i++)
		Q[i] = (float*)malloc(sizeof(float)*(objective_count));
	for (int i = 0; i < size; i++)
	{
		Q[i][0]=i;
		Q[i][1]=0.5*(float)(S[i]-Ss)/(Sb-Ss)+(1-0.5)*(float)(R[i]-Rs)/(Rb-Rs);
	}
		
	//************************************************************************find Best******************
	for (int i = 0; i<size; i++)
	{
		float min = Q[i][1];
		float loc = i;
		for (int j=i+1;j<size;j++)
		if (Q[j][1] <min)
		{
			min= Q[j][1];
			loc = j;
		}
		for (int j=0;j<2;j++)
		{
			float tmp=Q[(int)loc][j];
			Q[(int)loc][j]=Q[i][j];
			Q[i][j]=tmp;
		}
	}
	int A1=(int)Q[0][0],A2=(int)Q[1][0];

	for (int i = 0; i < size; i++)
	{
		delete []p[i] ;	
	}
	delete []Q[0];
	delete []Q[1];
	delete []Q;
	delete []sum;
	delete []fPositive;
	delete []fNegative;
	delete []S;
	delete []R;
	delete []E;
	delete []D;
	delete []W;
	return A1;

}
void fitness(float **neighbours,float ** limited_neighbours,Pseed_list Potential,Node_list *communities,int **community_links,int community_count,float **limited_community_links,Node_list Infetced,float *X,Edge_list CE,Edge_list CL,float *MaxC,float *f1,float *f2,float *f3,int key){
	float *One_hop = (float *)malloc(sizeof(float)*(node_count+1)); 
	float *two_hop = (float *)malloc(sizeof(float)*(node_count+1));
	float *potential_prob = (float *)malloc(sizeof(float)*(node_count+1));
	float **rho = (float **)malloc(sizeof(float)*(community_count+1));
	float **Pjtr = (float **)malloc(sizeof(float)*(community_count+1));
	for (int i = 0; i <= community_count; i++)
	{
		rho[i] = (float *)malloc(sizeof(float)*(community_count+1));
		Pjtr[i] = (float *)malloc(sizeof(float)*(community_count+1));
	}
		

	for (int i=1;i<=node_count;i++)
		potential_prob[i]=0;
	for (int i=1;i<=Potential.size();i++)
	{
		Pseed ps=Potential.get(i);
		potential_prob[ps.id]=ps.prob;
	}

	for (int i=1;i<=node_count;i++)
		for (int j=1;j<=node_count;j++)
			limited_neighbours[i][j]=neighbours[i][j];
	for (int i=1;i<=community_count;i++)
		for (int j=1;j<=community_count;j++)
			limited_community_links[i][j]=community_links[i][j];
	int cnt=0;
	for(int i=1;i<=CE.size();i++)
	{
		Edge e=CE.get(i);
		limited_neighbours[e.from][e.to]=X[cnt];
		limited_neighbours[e.to][e.from]=X[cnt];
		cnt++;
	}
	for(int i=1;i<=CL.size();i++)
	{
		Edge e=CL.get(i);
		limited_community_links[e.from][e.to]=X[cnt];
		limited_community_links[e.to][e.from]=X[cnt];
		cnt++;
	}


	for (int i=1;i<=community_count;i++)
	{
		limited_community_links[i][0]=1;
		for (int j=1;j<=community_count;j++)
			if(i!=j) limited_community_links[i][0]*=(1-(float)limited_community_links[i][j]/communities[i].size());
		limited_community_links[i][0]=1-limited_community_links[i][0];
	}
			

	for(int i=1;i<=community_count;i++)
		for (int j=1;j<=community_count;j++)
		{
			if(i==j)
			{
					rho[i][j]=0;
					Pjtr[i][j]=0;
			}
			else
			{
				rho[i][j]=calcuate_rho(i,j,communities,limited_community_links,community_count);	
				Pjtr[j][i]=calcuate_Pjtr(i,j,communities,limited_community_links,community_count)*(float)limited_community_links[j][i]/communities[j].size();	
			}
		}
	for(int j=1;j<=community_count;j++)
	{
		Pjtr[j][0]=1;
		for (int i=1;i<=community_count;i++)
			Pjtr[j][0]*=(1-Pjtr[j][i]);
		Pjtr[j][0]=1-Pjtr[j][0];
	}

	for (int i=1;i<=node_count;i++)
	{
		One_hop[i]=1;
		two_hop[i]=1;
	}
	for (int j=1;j<=node_count;j++)
	{
		if((potential_prob[j]==0)&&(!Infetced.contains(j)))
		{
			int tmp_com=neighbours[j][0];
			for(int i=1;i<=community_count;i++)
				One_hop[j]*=(1-rho[i][tmp_com]);
			for (int i=1;i<=node_count;i++)
			{
				One_hop[j]*=1-limited_neighbours[j][i]*potential_prob[i]*(1-(float)limited_community_links[tmp_com][0]);
			}
			One_hop[j]=(1-One_hop[j])*(1-(float)limited_community_links[tmp_com][0]);
			One_hop[j]+=Pjtr[tmp_com][0];
		}
	}
	for (int j=1;j<=node_count;j++)
	{
		if((potential_prob[j]==0)&&(!Infetced.contains(j)))
		{
			two_hop[j]*=(1-One_hop[j]);
			for (int i=1;i<=node_count;i++)
				if((potential_prob[i]==0)&&(!Infetced.contains(i)))
					{
						two_hop[j]*=(1-One_hop[i]*limited_neighbours[i][j]);
					}
			two_hop[j]=1-two_hop[j];
		}
	}
	if(key==1)
		for (int j=1;j<=node_count;j++)
			original_infected_nodes[j]=two_hop[j];
	else
	{
		*f1=0;
		float counter=0;
		for (int j=1;j<=node_count;j++)
			if((potential_prob[j]==0)&&(!Infetced.contains(j))&&(original_infected_nodes[j]>0))
			{
				*f1=*f1+(float)(original_infected_nodes[j]-two_hop[j]) /original_infected_nodes[j];
				counter++;
			}
		if(counter>0)
			*f1=(float)*f1/counter;
		*f2=0;
		*f3=0;

		cnt=0;
		counter=0;
		for(int i=1;i<=CE.size();i++)
		{
			Edge e=CE.get(i);
			*f2=*f2+(float)(MaxC[cnt]-X[cnt]) /MaxC[cnt];
			cnt++;
			counter++;
		}
		
		float counter3=0;
		for(int i=1;i<=CL.size();i++)
		{
			Edge e=CL.get(i);
			*f3=*f3+((float)(MaxC[cnt]-X[cnt]) /MaxC[cnt])*alpha*2;
			counter3+=MaxC[cnt]*alpha*2;
			cnt++;
		}
		
		if((*f2>0)&&(counter>0))
			*f2=(float)*f2/counter;
		if((*f3>0)&&(counter3>0))
			*f3=(float)*f3/counter3;
		*f2=*f2+*f3;
		
	}
	
	

	for (int i = 0; i <= community_count; i++)
	{
		delete[] rho[i];
		delete[] Pjtr[i];
	}
	delete[] potential_prob;
	delete[] rho;	
	delete[] Pjtr;
	delete[] One_hop;
	delete[] two_hop;

}
void Set_Edges_Capasity(float **neighbours,float ** limited_neighbours,Pseed_list Potential,Node_list *communities,int **community_links,int community_count,float **limited_community_links,Node_list Infetced,float *cost){
	
	Edge_list Candidate_edges=Edge_list(edge_count);
	Edge_list Candidate_links=Edge_list(edge_count);
	Node_list Potential_nodes=Node_list(node_count);
	Node_list Potential_communities=Node_list(node_count);
	
	for (int i=1;i<=community_count;i++)
	{
		community_links[i][0]=0;
		limited_community_links[0][i]=0;
	}
		
	for (int i=1;i<=Potential.size();i++)
	{
		Pseed ps=Potential.get(i);
		Potential_nodes.insert(ps.id);
		for (int i2=1;i2<=community_count;i2++)
			if(communities[i2].contains(ps.id))
			{
				community_links[i2][0]+=1;
				limited_community_links[0][i2]+=ps.prob;
				if (!Potential_communities.contains(i2))
				{
					Potential_communities.insert(i2);
				}
			}
	}

	for(int i=1;i<=Potential.size();i++)
	{
		Pseed ps=Potential.get(i);
		for (int j=1;j<=node_count;j++)
			if((!Infetced.contains(j))&&(neighbours[ps.id][j]>0.001)&&((!Potential_nodes.contains(j))||(ps.id<j)))
			{
				Edge Ed=Edge(ps.id,j,neighbours[ps.id][j]);
				Candidate_edges.insert(Ed);
			}
	}
	
	for (int i2=1;i2<=Potential_communities.size();i2++)
	{
		int i=Potential_communities.get(i2);
		for (int j=1;j<=community_count;j++)
			if((community_links[i][j]>0)&&((!Potential_communities.contains(j))||(i<j)))
			{
				Edge Ed=Edge(i,j,community_links[i][j]);
				Candidate_links.insert(Ed);
			}
	}
	
	float inertiai=0.9,inertiaf=0.4,c1i=2.5,c1f=0.5,c2i=0.5,c2f=2.5;
	int solution_size=Candidate_edges.size()+Candidate_links.size();
	float *Max_capasity = (float *)malloc(sizeof(float)*(solution_size));
	int cnt=0;
	
	for(int i=1;i<=Candidate_edges.size();i++)
	{
		Max_capasity[cnt++]=Candidate_edges.get(i).w;
	}
	
	for(int i=1;i<=Candidate_links.size();i++)
	{
		Max_capasity[cnt++]=Candidate_links.get(i).w;
	}

	float **X = (float **)malloc(sizeof(float)*(Population_size));
	float **V = (float **)malloc(sizeof(float)*(Population_size));
	float **P_Best = (float **)malloc(sizeof(float)*(Population_size));
	float **P_Best_fitness = (float **)malloc(sizeof(float)*(Population_size));
	float *G_Best = (float *)malloc(sizeof(float)*(solution_size));
	float *G_Best_fitness = (float *)malloc(sizeof(float)*(2));
	for (int i = 0; i < Population_size; i++)
	{
		X[i] = (float *)malloc(sizeof(float)*(solution_size));
		V[i] = (float *)malloc(sizeof(float)*(solution_size));
		P_Best[i] = (float *)malloc(sizeof(float)*(solution_size));
		P_Best_fitness[i] = (float *)malloc(sizeof(float)*(2));
	}
	float **tmp_compare = (float **)malloc(sizeof(float)*(2));
	tmp_compare[0] = (float *)malloc(sizeof(float)*(2));
	tmp_compare[1]= (float *)malloc(sizeof(float)*(2));
	//*********************************initial population******************	
	for (int i = 0; i <Population_size; i++)
	{
		P_Best_fitness[i][0]=0;
		P_Best_fitness[i][1]=0;

		for (int j = 0; j <solution_size; j++)
			V[i][j]=0;
	}
	G_Best_fitness[0]=0;
	G_Best_fitness[1]=0;
	
	for (int i = 0; i <Population_size; i++)
	{
		for (int j = 0; j < Candidate_edges.size(); j++)
		{
			X[i][j]=Max_capasity[j];
			if((float) (1 + rand() % (1000))/1000<0.5)
			{
				float rnd = (float) (rand() % ((int)(Max_capasity[j]*1000)))/1000;
				X[i][j]=rnd;
			}
		}
		for (int j = Candidate_edges.size(); j < solution_size; j++)
		{
			X[i][j]=Max_capasity[j];
			if((float) (1 + rand() % (1000))/1000<0.5)
			{
				int rnd = (rand() % ((int)(Max_capasity[j])));
				X[i][j]=rnd;
			}
		}
	}
	
	float f1,f2,f3;
	fitness(neighbours,limited_neighbours,Potential,communities,community_links,community_count,limited_community_links,Infetced,Max_capasity,Candidate_edges,Candidate_links,Max_capasity,&f1,&f2,&f3,1);
	for (int i = 0; i <Population_size; i++)
	{
		fitness(neighbours,limited_neighbours,Potential,communities,community_links,community_count,limited_community_links,Infetced,X[i],Candidate_edges,Candidate_links,Max_capasity,&f1,&f2,&f3,0);
		P_Best_fitness[i][0]=f1;
		P_Best_fitness[i][1]=f2;
		for (int j = 0; j <solution_size; j++)
			P_Best[i][j]=X[i][j];
	}
	int bst_solution=VIKOR(P_Best_fitness,Population_size);
	G_Best_fitness[0]=P_Best_fitness[bst_solution][0];
	G_Best_fitness[1]=P_Best_fitness[bst_solution][1];
	for (int j = 0; j <solution_size; j++)
		G_Best[j]=X[bst_solution][j];
	int g=0;
	int counter_progress=5;
	while((g<g_max)&&(counter_progress>0))
	{
		bool Update_P_Bests=false;
		for (int i = 0; i <Population_size; i++)
		{
			
			for (int j = 0; j <solution_size; j++)
			{
				float r1=(float) (rand() % (1000+1))/1000;
				float r2=(float) (rand() % (1000+1))/1000;
				float inertia=(inertiai-inertiaf)*((float)(g_max-g)/g_max)+inertiaf;
				float c1=(c1f-c1i)*((float)g/g_max)+c1i;
				float c2=(c2f-c2i)*((float)g/g_max)+c2i;
				V[i][j]=inertia*V[i][j]+c1*r1*(P_Best[i][j]-X[i][j])+c2*r2*(G_Best[j]-X[i][j]);
				if(j<Candidate_edges.size())
					X[i][j]=X[i][j]+V[i][j];
				else
					X[i][j]=(int)(X[i][j]+V[i][j]);
				if(X[i][j]>Max_capasity[j])
						X[i][j]=Max_capasity[j];
				if(X[i][j]<0)
					X[i][j]=0;
				
			}
			
			fitness(neighbours,limited_neighbours,Potential,communities,community_links,community_count,limited_community_links,Infetced,X[i],Candidate_edges,Candidate_links,Max_capasity,&f1,&f2,&f3,0);
			tmp_compare[0][0]=P_Best_fitness[i][0];
			tmp_compare[0][1]=P_Best_fitness[i][1];
			tmp_compare[1][0]=f1;
			tmp_compare[1][1]=f2;
			int better_solution=VIKOR(tmp_compare,2);
			if(better_solution==1)
			{
				for (int j = 0; j <solution_size; j++)
					P_Best[i][j]=X[i][j];
				P_Best_fitness[i][0]=f1;
				P_Best_fitness[i][1]=f2;
				Update_P_Bests=true;
			}
		}
		bool progress=false;
		if(Update_P_Bests)
		{
			bst_solution=VIKOR(P_Best_fitness,Population_size);

			if((G_Best_fitness[0]<P_Best_fitness[bst_solution][0])||(G_Best_fitness[1]>P_Best_fitness[bst_solution][1]))
				progress=true;

			G_Best_fitness[0]=P_Best_fitness[bst_solution][0];
			G_Best_fitness[1]=P_Best_fitness[bst_solution][1];
			for (int j = 0; j <solution_size; j++)
				G_Best[j]=X[bst_solution][j];
		}
		if (progress)
			counter_progress=5;
		else
			counter_progress--;
		//cout<<"g="<<g<<" Best F1="<<G_Best_fitness[0]<<"  F2="<<G_Best_fitness[1]<<"  "<<counter_progress<<endl;
		g=g+1;		
	}
	for (int i=1;i<=node_count;i++)
		for (int j=1;j<=node_count;j++)
			limited_neighbours[i][j]=neighbours[i][j];
	for (int i=1;i<=community_count;i++)
		for (int j=1;j<=community_count;j++)
			limited_community_links[i][j]=community_links[i][j];


	cnt=0;
	float tmpf=0;
	for(int i=1;i<=Candidate_edges.size();i++)
	{
		Edge e=Candidate_edges.get(i);
		limited_neighbours[e.from][e.to]=G_Best[cnt];
		limited_neighbours[e.to][e.from]=G_Best[cnt];
		tmpf+=(float)(Max_capasity[cnt]-G_Best[cnt]);
		cnt++;
	}
	*cost+=tmpf;
	//cost_in_timestamps[timestamp][0]+=tmpf;

	float tmpf2=0;
	for(int i=1;i<=Candidate_links.size();i++)
	{
		Edge e=Candidate_links.get(i);
		limited_community_links[e.from][e.to]=G_Best[cnt];
		limited_community_links[e.to][e.from]=G_Best[cnt];
		tmpf2+=((float)(Max_capasity[cnt]-G_Best[cnt]))*alpha*2;
		cnt++;
	}
	*cost+=tmpf2;
	//cost_in_timestamps[timestamp][1]+=tmpf2;

	Candidate_edges.release();
	Candidate_links.release();
	Potential_nodes.release();
	Potential_communities.release();
	delete[] Max_capasity;
	for (int i = 0; i < Population_size; i++)
	{
		delete[] X[i];
		delete[] V[i];
		delete[] P_Best[i];
		delete[] P_Best_fitness[i];
	}
	
	delete[] X ;
	delete[] V ;
	delete[] P_Best;
	delete[] P_Best_fitness;
	delete[] G_Best ;
	delete[] G_Best_fitness ;
	delete[] tmp_compare[0];
	delete[] tmp_compare[1];
	delete[] tmp_compare;
	
}
DWORD WINAPI Spreading_simulation_thread(LPVOID lpParam)
{
	srand(time(nullptr) + GetCurrentThreadId());// srand(GetCurrentThreadId());// srand(time(NULL));
	HANDLE hStdout;
	PThreadData pDataArray;
	TCHAR msgBuf[BUF_SIZE];
	size_t cchStringSize;
	hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	pDataArray = (PThreadData)lpParam;
	//*******************************************************************
	
	pDataArray->SIR[0].reset();
	pDataArray->SIR[1].reset();
	int timestamp=0;
	//********************************Seeds Reading*************************************************
	
	for(int i=1;i<=pDataArray->Seeds.size();i++)
	{
		int sd=pDataArray->Seeds.get(i);
		pDataArray->SIR[1].insert(sd);
	}
			
	pDataArray->Potential_Seeds.reset();
			
	Determin_Potential_seeds(pDataArray->Seeds,neighbours,pDataArray->temp_neighbours,&pDataArray->Potential_Seeds,pDataArray->SIR[1]);
			
	pDataArray->Seeds.reset();
	float tmpcost=0;		
	Set_Edges_Capasity(neighbours,pDataArray->limit_neighbours,pDataArray->Potential_Seeds,communities,community_links,community_count,pDataArray->limit_community_links,pDataArray->SIR[1],&tmpcost);
	pDataArray->cost=pDataArray->cost+tmpcost;		
	for(int i=1;i<=pDataArray->Potential_Seeds.size();i++)
	{
		Pseed ps=pDataArray->Potential_Seeds.get(i);
		float rnd =(float) (1 + rand() % (1000))/1000;
		if(rnd<=ps.prob)
		{
			pDataArray->SIR[0].insert(ps.id);
			pDataArray->Seeds.insert(ps.id);
		}
	}

	for (int j = 1; j <= node_count; j++)
		if ((!pDataArray->SIR[0].contains(j))&&(!pDataArray->SIR[1].contains(j)))
			pDataArray->SIR[2].set_cell(j,1);
				
	while (pDataArray->SIR[0].size() > 0)
	{
		timestamp++;
				
		pDataArray->passengers.reset();
		for (int i = 1; i <= node_count; i++)
			for (int j = 1; j <= node_count; j++)
				pDataArray->temp_neighbours[i][j]=0;

		Determine_passengers(communities,pDataArray->limit_community_links,community_count,&pDataArray->passengers,pDataArray->temp_neighbours,community_avg_weight);
		
		for(int i=1;i<=node_count;i++)
			for (int j=1;j<=node_count;j++)
				if((!pDataArray->passengers.contains(i))&&(!pDataArray->passengers.contains(j)))
					pDataArray->temp_neighbours[i][j]=pDataArray->limit_neighbours[i][j];

		//+++++++++++++++++++++++ first step++++++++++++++++++++++++++++++++++++++
		pDataArray->temp.reset();
		for(int i=1;i<=pDataArray->SIR[0].size();i++)
		{
			int s=pDataArray->SIR[0].get(i);
			pDataArray->temp.insert(s);
		}
				
		for (int i=1;i<=pDataArray->temp.size();i++)
		{
			int s=pDataArray->temp.get(i);
			for (int sn=1;sn<=node_count;sn++) 
			{
				if ((pDataArray->temp_neighbours[s][sn]>0)&&(pDataArray->SIR[2].get(sn)==1))
				{
					float rnd = (float) (1 + rand() % (1000))/1000;
					if (rnd <= pDataArray->temp_neighbours[s][sn])
					{
						pDataArray->SIR[0].insert(sn);

						pDataArray->SIR[2].set_cell(sn,0);
								
						pDataArray->Seeds.insert(sn);
					}
				}
			}
					
			pDataArray->SIR[0].remove(s);
			pDataArray->SIR[1].insert(s);
		}
		
		//+++++++++++++++++++++++ Second step++++++++++++++++++++++++++++++++++++++
		pDataArray->temp.reset();
		for(int i=1;i<=pDataArray->SIR[0].size();i++)
		{
			int s=pDataArray->SIR[0].get(i);
			pDataArray->temp.insert(s);
		}

		for (int i=1;i<=pDataArray->temp.size();i++)
		{
			int s=pDataArray->temp.get(i);
			for (int sn=1;sn<=node_count;sn++) 
			{
				if ((pDataArray->limit_neighbours[s][sn]>0)&&(pDataArray->SIR[2].get(sn)==1))
				{
					float rnd = (float) (1 + rand() % (1000))/1000;
					if (rnd <= pDataArray->limit_neighbours[s][sn])
					{

						pDataArray->SIR[1].insert(sn);

						pDataArray->SIR[2].set_cell(sn,0);
								
						pDataArray->Seeds.insert(sn);

					}
				}
			}
			pDataArray->SIR[0].remove(s);
			pDataArray->SIR[1].insert(s);
		}
		
		//++++++++++++++++++++++set values for the next timestamp++++++++++++++++++++++++++++++
		pDataArray->Potential_Seeds.reset();
		Determin_Potential_seeds(pDataArray->Seeds,pDataArray->limit_neighbours,pDataArray->temp_neighbours,&pDataArray->Potential_Seeds,pDataArray->SIR[1]);

		pDataArray->Seeds.reset();
		for (int i = 1; i <= node_count; i++)
			for (int j = 1; j <= node_count; j++)
			{
				pDataArray->temp_neighbours[i][j]=0;
				pDataArray->limit_neighbours[i][j]=neighbours[i][j];
			}
		for (int i = 1; i <= community_count; i++)
			for (int j = 1; j <= community_count; j++)
				pDataArray->limit_community_links[i][j]=community_links[i][j];
		float tmpcost=0;
		Set_Edges_Capasity(neighbours,pDataArray->limit_neighbours,pDataArray->Potential_Seeds,communities,community_links,community_count,pDataArray->limit_community_links,pDataArray->SIR[1],&tmpcost);
		pDataArray->cost=pDataArray->cost+tmpcost;
		for(int i=1;i<=pDataArray->Potential_Seeds.size();i++)
		{
			Pseed ps=pDataArray->Potential_Seeds.get(i);
			float rnd =(float) (1 + rand() % (1000))/1000;
			if(rnd<=ps.prob)
			{
				pDataArray->SIR[0].insert(ps.id);
				
				pDataArray->SIR[2].set_cell(ps.id,0);

				pDataArray->Seeds.insert(ps.id);
			}
		}
		
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}
	//*****************************************************
	pDataArray->Result+=pDataArray->SIR[1].size();
	pDataArray->status = 1;//done
	return pDataArray->Result;
}

int _tmain(int argc, _TCHAR* argv[])
{
	PThreadData pDataArray[MAX_THREADS+1];
	DWORD   dwThreadIdArray[MAX_THREADS+1];
	HANDLE  hThreadArray[MAX_THREADS+1];
	//***********************************************
	clock_t Run_time = clock();
	int total_result=0;
	float cost_f2=0;
	int number_of_simulations = MAX_THREADS,inf_community;
	ofstream Out("output.dat");


	for(inf_community=1;inf_community<=seed_set_count;inf_community++)
	{
		cout<<endl<<"Simulation for random initial infected set "<<inf_community<<" of "<<seed_set_count<<":   Started"<<endl;

		
		//************************read communities from file*********************************************
		
		int a,b;
		ifstream fp1_community(file_community);
		while (fp1_community)
		{
			fp1_community>>a;
			if(fp1_community)
			{
				fp1_community>>b;
				if(b>community_count)
					community_count=b;
			}
		}
		fp1_community.close();

		communities=(Node_list*)malloc(sizeof(Node_list)*(community_count+ 1));
		for (int i=0;i<=community_count;i++)
			communities[i]=Node_list(node_count);
		for (int i=0;i<=community_count;i++)
			for (int j=0;j<=node_count;j++)
				communities[i].set_cell(j,0);
		
		
		fp1_community=ifstream(file_community);
		while (fp1_community)
		{
			fp1_community>>a;
			if(fp1_community)
			{
				fp1_community>>b;
				communities[b].insert(a);
			}
		}
		fp1_community.close();

		//***********************read neighbours and determine community weights****************
		community_links = (int **)malloc(sizeof(int)*(community_count + 1));
		for (int i = 1; i <= community_count; i++)
			community_links[i] = (int *)malloc(sizeof(int)*(community_count + 1));
		for (int i = 1; i <= community_count; i++)
			for (int j = 1; j <= community_count; j++)
				community_links[i][j]=0;
		
		community_avg_weight = (float **)malloc(sizeof(float)*(community_count + 1));
		for (int i=0;i<=community_count;i++)
			community_avg_weight[i] = (float *)malloc(sizeof(float)*(3+1));
		for (int i=0;i<=community_count;i++)
			for(int j=0;j<=3;j++)
				community_avg_weight[i][j]=0;

		neighbours = (float **)malloc(sizeof(float)*(node_count + 1));
		for (int i = 1; i <= node_count; i++)
			neighbours[i] = (float *)malloc(sizeof(float)*(node_count + 1));
		for (int i = 1; i <= node_count; i++)
			for (int j = 1; j <= node_count; j++)
				neighbours[i][j]=0;
		
		
		ifstream fp1_network(file_network);
		float w,max=-1,min=1,average=0,average_linksW=0;
		int link_count=0;
		while(fp1_network)
		{
			fp1_network>>a;
			if(fp1_network)
			{
				fp1_network>>b;
				fp1_network>>w;
				if(a<b)
				{
					int i=0,j=0;
					for (int t=1;t<=community_count;t++)
					{
						if(communities[t].contains(a))
							i=t;
						if(communities[t].contains(b))
							j=t;
					}
					if(i==j)
					{
						neighbours[a][b]=w;
						neighbours[b][a]=w;
						if(w>max)
							max=w;
						if(w<min)
							min=w;
						average+=w;
						edge_count+=1;
						community_links[i][j]+=1;
						community_avg_weight[i][1]+=w;
						community_avg_weight[i][2]++;
					}
					else
					{
						average_linksW+=w;
						link_count+=1;
						community_links[i][j]+=1;
						community_links[j][i]+=1;
					}
				}
			}
		}
		fp1_network.close();
		for(int i=1;i<=community_count;i++)
			community_avg_weight[i][3]=(float)community_avg_weight[i][1]/community_avg_weight[i][2];
		for(int i=1;i<=community_count;i++)
			for (int j=1;j<=communities[i].size();j++)
			{
				int nde=communities[i].get(j);
				neighbours[nde][0]=i;
			}
		original_infected_nodes=(float *)malloc(sizeof(float)*(node_count + 1));
		srand(time(0));
		

		
		for(int sim=1;sim<=number_of_simulations;sim++)
		{
			ifstream fp1_seed(file_seed+std::to_string(inf_community)+".dat");
			Node_list Seeds = Node_list(node_count);
			int seed;
			while(fp1_seed)
			{
				fp1_seed>>seed;
				if(fp1_seed)
				{
					Seeds.insert(seed);
				}
			}
			
			// Allocate memory for thread data.
			pDataArray[sim] = (PThreadData)HeapAlloc(GetProcessHeap(), HEAP_ZERO_MEMORY,sizeof(ThreadData));
			pDataArray[sim]->id=sim;
			pDataArray[sim]->Result=0;
			pDataArray[sim]->cost=0;
			pDataArray[sim]->status=0;
			pDataArray[sim]->print=0;

			pDataArray[sim]->Seeds=Node_list(node_count);
			for(int i=1;i<=Seeds.size();i++)
				pDataArray[sim]->Seeds.insert(Seeds.get(i));

			for (int i = 0; i <3; i++)
			{
				pDataArray[sim]->SIR[i] = Node_list(node_count);
			}

			pDataArray[sim]->temp =Node_list(node_count);
			pDataArray[sim]->passengers =Node_list(node_count);
			pDataArray[sim]->Potential_Seeds =Pseed_list(node_count);

			pDataArray[sim]->temp_neighbours = (float **)malloc(sizeof(float)*(node_count + 1));
			pDataArray[sim]->limit_neighbours = (float **)malloc(sizeof(float)*(node_count + 1));
			for (int i = 1; i <= node_count; i++)
			{
				pDataArray[sim]->temp_neighbours[i] = (float *)malloc(sizeof(float)*(node_count + 1));
				pDataArray[sim]->limit_neighbours[i] = (float *)malloc(sizeof(float)*(node_count + 1));
			}
			for (int i = 1; i <= node_count; i++)
				for (int j = 1; j <= node_count; j++)
				{
					pDataArray[sim]->temp_neighbours[i][j]=0;
					pDataArray[sim]->limit_neighbours[i][j]=neighbours[i][j];
				}
			pDataArray[sim]->limit_community_links = (float **)malloc(sizeof(float)*(community_count + 1));
			for (int i = 0; i <= community_count; i++)
				pDataArray[sim]->limit_community_links[i] = (float *)malloc(sizeof(float)*(community_count + 1));
			for (int i = 1; i <= community_count; i++)
				for (int j = 1; j <= community_count; j++)
					pDataArray[sim]->limit_community_links[i][j]=community_links[i][j];
		

			hThreadArray[sim] = CreateThread(
				NULL,						
				0,							
				Spreading_simulation_thread,
				pDataArray[sim],			
				0,							
				&dwThreadIdArray[sim]);		
			
			Seeds.release();
			cout<<"\nThread "<<pDataArray[sim]->id<< "  of "<<MAX_THREADS<<": Started"<<endl;
		}
		
		
		bool flag=true;
		while(flag)
		{
			flag=false;
			for(int sim=1;sim<=number_of_simulations;sim++)
				if(pDataArray[sim]->status==0)
					flag=true;
				else
					if(pDataArray[sim]->print==0)
					{
						pDataArray[sim]->print=1;
						cout<<"\nThread "<<pDataArray[sim]->id<< "  of "<<MAX_THREADS<<": Done"<<endl;
					}
		}
		DWORD dw=WaitForMultipleObjects(MAX_THREADS, hThreadArray, TRUE, INFINITE);

		float tmp_infected=0,tmp_cost=0;
		for (int i = 1; i<=MAX_THREADS; i++)
		{
			CloseHandle(hThreadArray[i]);
			if (pDataArray[i] != NULL)
			{
				tmp_infected+=pDataArray[i]->Result;
				tmp_cost+=pDataArray[i]->cost;
				HeapFree(GetProcessHeap(), 0, pDataArray[i]);
				pDataArray[i] = NULL;    // Ensure address is not reused.
			}
		}
		total_result+=tmp_infected;
		cost_f2+=tmp_cost;
		Out<<(float)tmp_infected/MAX_THREADS<<"\t"<<(float)tmp_cost/MAX_THREADS<<endl;
		cout<<endl<<"Simulation for random initial infected set "<<inf_community<<" of "<<seed_set_count<<":   Done"<<endl;
	}
	Out.close();
	cout<<"\n----------------------Results--------------------------------\n";
	cout<<"total infected number->"<<(float)total_result/(number_of_simulations*seed_set_count)<<endl;
	cout<<"total Cost->"<<(float)cost_f2/(number_of_simulations*seed_set_count)<<endl;
	cout << "\ntime for executing=" << ((float)Run_time) / CLOCKS_PER_SEC << "  second\n";
	system("pause");
}

