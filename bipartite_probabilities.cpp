////////////////////////////////////////////////////////////////////////////////////
//This file is part of the implementaiton of lpa algorithms for bipartite graphs  //
//The program calculates the relation matrix, which is used during the propagation//
////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<set>
#include<stdlib.h>
#include<map>

using namespace std;

class BipartiteProbabilities
{

    vector<vector<int> > gseed;
    map<int,int> seed_nodes;
    set<int> labeling_nodes;
    string graph, seeds,labeling, output;
    int total_labels, label_number, groups;

public:
    BipartiteProbabilities(string graph, string seeds, string labeling, string output, int total_labels, int label_number, int groups)
	{
    	this->graph=graph;
    	this->seeds=seeds;
    	this->labeling=labeling;
    	this->output=output;
    	this->total_labels=total_labels;
    	this->label_number=label_number;
    	this->groups=groups;
	}

    ~BipartiteProbabilities(){}

    //getting nodes which should be labaled after label propagation is done

    void readNodestoLabel()
	{
    	int id;
        ifstream ls(labeling.c_str());
        while(ls>>id)
        {
        	labeling_nodes.insert(id);
        }
        ls.close();
	}

    //getting seed nodes with their corresponding labels

    void readSeedNodes()
    {
    	int label;
    	int id;
    	ifstream seedStream(seeds.c_str());
		while(seedStream>>id)
		{
			for(int i=0; i<total_labels; i++)
			{
				seedStream>>label;
				if(i==label_number)
				{
					seed_nodes[id]=label;
				}
			}
		}
		seedStream.close();
    }

    //creating the graph using seed nodes to compute relation matrix

    void createSeedGraph()
    {
    	int nodes=0;
    	ifstream gr(graph.c_str());
    	int id1, id2;
    	while(gr>>id1)
    	{
    		if(id1>nodes) nodes=id1;
    	}
    	nodes++;
    	gr.close();
    	gseed.resize(nodes);
    	for(int i=0; i<nodes; i++)
    	{
    		gseed[i].resize(groups);
    	}
    	gr.open(graph.c_str());
    	while(gr>>id1)
    	{
    		gr>>id2;
    		if(labeling_nodes.find(id1)==labeling_nodes.end())
    		{
    			if(seed_nodes.find(id2)!=seed_nodes.end())
    			{
    				gseed[id1][seed_nodes[id2]]++;
    			}
    		}
    		else
    		{
    			if(seed_nodes.find(id1)!=seed_nodes.end())
    			{
    				gseed[id2][seed_nodes[id1]]++;
    			}
    		}
    	}
    }

    //calculation and printing of relation matrix

    void getProbabilities()
    {
    	readNodestoLabel();
    	readSeedNodes();
    	createSeedGraph();
    	vector<vector<double> > probabilities;
    	probabilities.resize(groups);
    	for(int i=0; i<groups; i++)
    		probabilities[i].resize(groups);
    	for(int i=0; i<gseed.size(); i++)
    	{
    		if(gseed[i].size()>1)
    		{
    			for(int j=0; j<groups; j++)
    			{
    				for(int k=0; k<groups; k++)
    				{
    					if(j==k)
    					{
    						probabilities[j][k]+=gseed[i][j]*(gseed[i][j]-1);
    					}
    					else
    						probabilities[j][k]+=gseed[i][j]*gseed[i][k];
    				}
    			}
    		}
    	}

    	double sum;
    	for(int i=0; i<groups; i++)
    	{
    		sum=0;
    		for(int j=0; j<groups; j++)
    			sum+=probabilities[j][i];
    		if(sum>0)
    		{
    			for(int j=0; j<groups; j++)
					probabilities[j][i]/=sum;
    		}
    	}

		ofstream result(output.c_str());
		for(int i=0; i<groups; i++)
		{
			for(int j=0; j<groups; j++)
			{
				result<<probabilities[i][j]<<"\t";
			}
			result<<endl;
		}
    }

};


int main(int argc, char * argv[])
{
    cout<<"started\n";
    if(argc<8)
    {
        cout<<"Four paths to files are needed:\n";
        cout<<"graph, seed nodes, nodes wanted to label and output file\n";
        cout<<"Also total number of labels, label number which is wanted to propagate and number of groups is needed\n";
        return 0;
    }
    BipartiteProbabilities biProb(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
    biProb.getProbabilities();
    return 0;
}
