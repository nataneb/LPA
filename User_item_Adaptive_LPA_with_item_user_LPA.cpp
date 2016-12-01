////////////////////////////////////////////////////////////////////////
// This file is part of Grappa, a system for scaling irregular
// applications on commodity clusters.

// Copyright (C) 2010-2014 University of Washington and Battelle
// Memorial Institute. University of Washington authorizes use of this
// Grappa software.

// Grappa is free software: you can redistribute it and/or modify it
// under the terms of the Affero General Public License as published
// by Affero, Inc., either version 1 of the License, or (at your
// option) any later version.

// Grappa is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Affero General Public License for more details.

// You should have received a copy of the Affero General Public
// License along with this program. If not, you may obtain one from
// http://www.affero.org/oagl.html.
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// Use the GraphLab API to implement Connected Components
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<stdint.h>
#include <Grappa.hpp>
#include<vector>
#include<set>
#include "graphlab.hpp"
#include <GlobalHashMap.hpp>

DEFINE_bool( metrics, false, "Dump metrics");

DEFINE_int32(scale, 10, "Log2 number of vertices.");
DEFINE_int32(edgefactor, 16, "Average number of edges per vertex.");

DEFINE_int32(trials, 1, "Number of timed trials to run and average over.");

DEFINE_int32(labels_total, 1, "Total number of different labels (gender,age,etc.)");
DEFINE_int32(label_number, 0, "Number of label (in seed file) used in label propagation");
DEFINE_int32(clamping, 1, "Number of sequence of iterations after which clamping happens");

DEFINE_string(path, " ", "Path to graph source file.");
DEFINE_string(format, "tsv", "Format of graph source file.");
DEFINE_string(seedsPath, "", "Path to the seeds");
DEFINE_string(outputPath, "", "Path and file name to the result.");
DEFINE_string(transitionMatrixPath, "", "Path and file name to the transition matrix for Bipartite Adaptive LPA");
DEFINE_string(labelingVerticesPath,"", "Path and file name to the vertices we want to label");

GRAPPA_DEFINE_METRIC(SimpleMetric<double>, init_time, 0);
GRAPPA_DEFINE_METRIC(SimpleMetric<double>, tuple_time, 0);
GRAPPA_DEFINE_METRIC(SimpleMetric<double>, construction_time, 0);
GRAPPA_DEFINE_METRIC(SummarizingMetric<double>, total_time, 0);

const int Unknown = -1;
const int number_of_groups = 2; //number of ifferent classes

double transition[number_of_groups][number_of_groups];

// Graph vertex data
struct VertexData : public GraphlabVertexData<VertexData> {
  int label;
  VertexID id;
  int vid;
  int iter;
  bool check_seed;
  bool labeling;
};

//Structure for propagation

struct PropagationInformation {
  double labelsDistribution[number_of_groups];
  double nadj;
  int vId;
  bool labeling;

  explicit PropagationInformation(int label, int _nadj, int id, bool labeling) {
    CHECK(_nadj >= 0) << "We expect all nodes to have at least one neighbor";
    CHECK(label < number_of_groups) << "Label should be less than supported number of groups.";
    for (int i =0 ; i < number_of_groups; i++) {
       labelsDistribution[i] = 0.0;
    }

    nadj = _nadj;
    this->labeling=labeling;
    vId=id;
    if (label != Unknown)
      labelsDistribution[label] = 1.0;
  }

  explicit PropagationInformation() {
    for (int i =0 ; i < number_of_groups; i++) {
       labelsDistribution[i] = 0.0;
    }
  }

  explicit PropagationInformation(int _nadj, int id, bool labeling) {
     CHECK(_nadj >= 0) << "We expect all nodes to have at least one neighbor";
    vId=id;
    nadj = _nadj;
    // Init the array
    for (int i =0 ; i < number_of_groups; i++) {
       labelsDistribution[i] = 0.0;
    }
    this->labeling=labeling;
  }

  void setID(int id)
  {
	  vId=id;
  }

  //Normalizing an array

  void normalize()
  {
	  double sum =0;
	  for(int i=0; i<number_of_groups; i++)
		  sum+=labelsDistribution[i];
	  if(sum==0) return ;
	  for(int i=0; i<number_of_groups; i++)
		  labelsDistribution[i]/=sum;
  }

  //Setting a probability for a specific class

  void addValue(int index, double value) {
    CHECK( index > -1 && index < number_of_groups);
    labelsDistribution[index] = value;
  }

  //Overloading of the operator based on the type of the node, used during Gather and Scatter phases

  PropagationInformation& operator+=(const PropagationInformation& other) {
    CHECK(other.nadj > 0) << "Neighbor should have atleast has another neighbor.";
    
    double sum=0;
    if(!(other.labeling))
        {
    		for (int i = 0; i < number_of_groups; i++) {
    			labelsDistribution[i] += (other.labelsDistribution[i])/other.nadj;
    		}
        }
        else
        {
        	for(int i=0; i<number_of_groups; i++)
        	{
        		sum=0;
        		for(int j=0; j<number_of_groups; j++)
        		{
        			sum+=other.labelsDistribution[j]*transition[i][j];
        		}
        		sum/=(other.nadj);
        		labelsDistribution[i] +=sum;
        	}
        }	

    return *this;
  }
};

using G = Graph<VertexData,Empty>;

//Defining Gather, Apply and Scatter phases for label propagation

struct LabelPropagation : public GraphlabVertexProgram<G, PropagationInformation> {
  bool do_scatter;
  PropagationInformation changelabel; //Defining the information for propagation for each vertex

  LabelPropagation(Vertex& v) {
    do_scatter = false;
    changelabel = PropagationInformation(v.nadj, v->vid, v->labeling);
  }

  bool gather_edges(const Vertex& v) const { return true; }

  Gather gather(const Vertex& v, Edge& e) const {

    return PropagationInformation(v->label, v.nadj, v->vid, v->labeling);
  }

  /* Getting the most probable label for a node and making a decision to activate() a node for next iteration
   activate() for a node means that it does not meet the criterion we use for stopping,
   namely none of the nodes change a label using Seed-Clamping strategy*/
  void apply(Vertex& v,  Gather& total) {
	 if(v->check_seed){
		 	 for(int i=0; i<number_of_groups; i++)
		 	 {
		 		 total.labelsDistribution[i]=0.0;
		 	 }
		 	 changelabel.nadj=v.nadj;
			 changelabel.labelsDistribution[v->label]=1.0;
		 return;
	 }

    double maxCount = 0.0;
    for(int i=0; i<number_of_groups; i++)
    {
    	changelabel.labelsDistribution[i]=total.labelsDistribution[i];
    	total.labelsDistribution[i]=0;
    }
    changelabel.normalize();
    int maxLabel = v->label;
    changelabel.vId=v->vid;

    for (int i = 0; i < number_of_groups; i++) {

        if (changelabel.labelsDistribution[i] > maxCount) {
          maxCount = changelabel.labelsDistribution[i];
          maxLabel = i;
        }
    }
		if (maxLabel == Unknown||maxCount==0)
		{
		  //do_scatter = false;
		  v->activate();
		}
		else if (maxLabel != v->label)
		{
			v->label = maxLabel;
			v->iter=0;
			v->activate();
		}
  }

  //Decides if Scatter should occur from the node or not

  bool scatter_edges(const Vertex& v) {
	  if(v->label==Unknown) return false;
	  return true;
  }

  Gather scatter(const Edge& e, Vertex& target) const {
    return changelabel;
  }
};

// loading seed nodes with corresponing labels from the file

std::map<int64_t, int> loadSeed(std::string path)
{
  std::ifstream file(path);
  LOG(INFO)<<path<<" loading";
  std::map<int64_t , int> seeds;
  int id;
  int label;
  while (file>>id)
  {
      CHECK(id >= 0);
      for(int i=0; i<FLAGS_labels_total; i++)
      {
    	  file>>label;
    	  if(i==FLAGS_label_number)
    	  	  seeds[id] = label;
      }
  }
  LOG(INFO) << "Loaded Seed file. Count: " << seeds.size();

  return seeds;
}


//loading transition matrix from file

void loadTransitionMatrix(std::string path)
{
	std::ifstream file(path);
	for(int i=0; i<number_of_groups; i++)
	{
		for(int j=0; j<number_of_groups; j++)
		{
			file>>transition[i][j];
		}
	}
}

//loading the vertices (not intermediaries) for which the algorithm should infer labels

std::set<int> loadVerticesForLabeling(std::string path)
{
	std::set<int> vertices;
	int id;
	std::ifstream file(path);
	while(file>>id)
	{
		CHECK(id>=0)<<"Vertex Id can't be less than zero";
		vertices.insert(id);
	}
	LOG(INFO) << "Loaded labeling file. Count: " << vertices.size();
	return vertices;
}

int main(int argc, char* argv[]) {

  init(&argc, &argv);
  run([]{

    double t;

    TupleGraph tg;

    GRAPPA_TIME_REGION(tuple_time) {
      if (FLAGS_path.empty()) {
        LOG(INFO) << "We need to have a path to a graph.";
      } else {
        LOG(INFO) << "loading " << FLAGS_path;
        tg = TupleGraph::Load(FLAGS_path, FLAGS_format);
      }
    }

    LOG(INFO) << tuple_time;
    LOG(INFO) << "constructing graph";
    t = walltime();

    auto g = G::Undirected(tg);
    LOG(INFO) << "Finished loading graph";
    construction_time = walltime()-t;
    LOG(INFO) << construction_time;

    for (int i = 0; i < FLAGS_trials; i++) {
      if (FLAGS_trials > 1) LOG(INFO) << "trial " << i;
      // Load seed file and vertices which need labels
      LOG(INFO) << "start loading seed file";
      {
      symmetric_static std::map<int64_t, int> symmetric_seeds;
      symmetric_static std::set<int> verticesForLabeling;
          on_all_cores( []{
              symmetric_seeds = loadSeed(FLAGS_seedsPath);
              verticesForLabeling=loadVerticesForLabeling(FLAGS_labelingVerticesPath);
              loadTransitionMatrix(FLAGS_transitionMatrixPath);
            });

       forall(g, [](VertexID i, G::Vertex& v){
                    v->vid= i;
                    v->iter=0;
                  });
        LOG(INFO) << "Init the labels";
        LOG(INFO) << "Size of the graph is: " << g->nv;
        forall(g, [](VertexID i, G::Vertex& v){

          //setting the labels for seed nodes to use them during propagation
          if (symmetric_seeds.find(i) == symmetric_seeds.end())
          {
            v->label = -1.0;
            v->check_seed=false;
            v->activate();
          } else {
            v->label = symmetric_seeds[i];
            v->check_seed=true;
            v->labeling=true;
            v->activate();
          }

          //setting types for verices
          if(verticesForLabeling.find(i)!=verticesForLabeling.end())
		  {
        	v->labeling=true;
		  }
          else
          {
        	v->labeling=false;
          }
        });
      }

      GRAPPA_TIME_REGION(total_time) {
        LOG(INFO) << "Init is complete";
        NaiveGraphlabEngine<G,LabelPropagation>::run_sync(g);
      }
      LOG(INFO)<<"Complete";
    }

    LOG(INFO) << total_time;

    // for easy parallel writes, we make a separate file per core
    {
      symmetric_static std::ofstream myFile;
      int pid = getpid();
      on_all_cores( [pid] {
          std::ostringstream oss;
          oss <<FLAGS_outputPath<<"-"<<pid<<"-" << mycore();
          new (&myFile) std::ofstream(oss.str());
          if (!myFile.is_open()) exit(1);
        });
      forall(g, [](VertexID i, G::Vertex& v){
          myFile << i << " " << v->label << "\n";
        });
      on_all_cores( [] {
          myFile.close();
        });
    }

    if (FLAGS_metrics) Metrics::merge_and_print();
    else { std::cerr << total_time << "\n"<< iteration_time << "\n";

    }
    Metrics::merge_and_dump_to_file();
  });
  finalize();
}
