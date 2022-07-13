
class base{
public:
	int count;
	int size();
	void reset();
};
int base::size()
{
	return count;
}
void base::reset()
{
	count=0;
}

class Pseed{
public:	
	int id;
	float prob;
	Pseed(int i,float p){id=i;prob=p;}
};
class Pseed_list:public base{
	Pseed *elements;
public:
	Pseed_list(int);
	Pseed get(int);
	void insert(Pseed);
};
Pseed_list::Pseed_list(int nd_count)
{
	count=0;
	elements=(Pseed *)malloc(sizeof(Pseed)*(nd_count + 1));
}
Pseed Pseed_list::get(int index)
{
	Pseed Ps=Pseed(elements[index].id,elements[index].prob);
	return Ps;
}
void Pseed_list::insert(Pseed Ps)
{
	count++;
	elements[count].id=Ps.id;
	elements[count].prob=Ps.prob;
}


class Edge{
public:	
	int from;
	int to;
	float w;
	Edge(int i,int j,float p){from=i;to=j;w=p;}
};
class Edge_list:public base
{
	Edge *elements;
public:
	Edge_list(int);
	Edge get(int);
	void insert(Edge);
	void release();
};
Edge_list::Edge_list(int ed_count)
{
	count=0;
	elements=(Edge *)malloc(sizeof(Edge)*(ed_count + 1));
}
Edge Edge_list::get(int index)
{
	Edge Ed=Edge(elements[index].from,elements[index].to,elements[index].w);
	return Ed;
}
void Edge_list::insert(Edge Ed)
{
	count++;
	elements[count].from=Ed.from;
	elements[count].to=Ed.to;
	elements[count].w=Ed.w;
}
void Edge_list::release()
{
	delete []elements;
}

class Node_list:public base
{
	int *elements;
public:
	Node_list(int);
	int get(int);
	void insert(int);
	void release();
	int contains(int);
	void set_cell(int,int);
	void remove(int);
};
Node_list::Node_list(int nd_count)
{
	count=0;
	elements=(int *)malloc(sizeof(int)*(nd_count + 1));
}
int Node_list::get(int index)
{
	int Nd=elements[index];
	return Nd;
}
void Node_list::insert(int Nd)
{
	count++;
	elements[count]=Nd;

}
void Node_list::release()
{
	delete []elements;
}
int Node_list::contains(int x)
{
	for (int i=1;i<=count;i++)
		if(elements[i]==x)
			return 1;
	return 0;
}
void Node_list::set_cell(int index,int value)
{
	elements[index]=value;
}
void Node_list::remove(int x)
{
	int loc=0;
	for (int i=1;(i<=count)&&(!loc);i++)
		if(elements[i]==x)
			loc=i;
	for (int i=loc;i<count;i++)
		elements[i]=elements[i+1];
	count--;
}

typedef struct ThreadData {
	int id;
	Node_list Seeds;
	Node_list SIR[3];
	Node_list temp;
	Node_list passengers;
	float **temp_neighbours;
	float **limit_neighbours;
	Pseed_list Potential_Seeds;
    float **limit_community_links;
	float Result;
	float cost;
	int status;
	int print;
} ThreadData, *PThreadData;
