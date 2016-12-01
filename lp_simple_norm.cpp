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
#include "graphlab.hpp"
#include<math.h>
#include <GlobalHashMap.hpp>

DEFINE_bool( metrics, false, "Dump metrics");

DEFINE_int32(scale, 10, "Log2 number of vertices.");
DEFINE_int32(edgefactor, 16, "Average number of edges per vertex.");

DEFINE_int32(trials, 1, "Number of timed trials to run and average over.");

DEFINE_int32(labels_total, 1, "Total number of different labels (gender,age,etc.)");
DEFINE_int32(label_number, 0, "Number of label (in seed file) used in label propagation");


DEFINE_string(path, " ", "Path to graph source file.");
DEFINE_string(format, "tsv", "Format of graph source file.");
DEFINE_string(seedsPath, "", "Path to the seeds");
DEFINE_string(outputPath, "", "Path and file name to the result.");

GRAPPA_DEFINE_METRIC(SimpleMetric<double>, init_time, 0);
GRAPPA_DEFINE_METRIC(SimpleMetric<double>, tuple_time, 0);
GRAPPA_DEFINE_METRIC(SimpleMetric<double>, construction_time, 0);
GRAPPA_DEFINE_METRIC(SummarizingMetric<double>, total_time, 0);

const int Unknown = -1;
const int number_of_groups = 2; //number of different classes
double threshold=0.01;

// Graph vertex data
struct CCData : public GraphlabVertexData<CCData> {
  int label;
  VertexID id;
  int vid;
  int iter;
  bool labeled=false;
  bool check_seed;
};

//Structure for propagation

struct label_counter {
  double label_count[number_of_groups];
  double nadj;
  int vId;

  explicit label_counter(int label, int _nadj, int id) {
    CHECK(_nadj >= 0) << "We expect all nodes to have at least one neighbor";
    CHECK(label < number_of_groups) <<label<<" "<<id<< " Label should be less than supported number of groups.";
    for (int i =0 ; i < number_of_groups; i++) {
       label_count[i] = 0.0;
    }

    nadj = _nadj;
    vId=id;
    if (label != Unknown)
      label_count[label] = 1.0;
  }

  explicit label_counter() {
    for (int i =0 ; i < number_of_groups; i++) {
       label_count[i] = 0.0;
    }
  }

  explicit label_counter(int _nadj, int id) {
     CHECK(_nadj >= 0) << "We expect all nodes to have at least one neighbor";
    vId=id;
    nadj = _nadj;
    // Init the array
    for (int i =0 ; i < number_of_groups; i++) {
       label_count[i] = 0.0;
    }
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
		  sum+=label_count[i];
	  if(sum==0) return ;
	  for(int i=0; i<number_of_groups; i++)
		  label_count[i]/=sum;
  }

  //Setting a probability for a specific class

  void addValue(int index, double value) {
    CHECK( index > -1 && index < number_of_groups);
    label_count[index] = value;
  }

  //Overloading of the operator used during Gather and Scatter phases

  label_counter& operator+=(const label_counter& other) {

    CHECK(other.nadj > 0) << "Neighbor should have at least one other neighbor.";


    for (int i = 0; i < number_of_groups; i++) {
      label_count[i] += (other.label_count[i])/other.nadj;
    }

    return *this;
  }

};

using G = Graph<CCData,Empty>;

//Defining Gather, Apply and Scatter phases for label propagation

struct LabelPropagation : public GraphlabVertexProgram<G, label_counter> {
  bool do_scatter;
  label_counter changelabel; //Defining the distribution over labels for each vertex

  LabelPropagation(Vertex& v) {
    do_scatter = false;
    changelabel = label_counter(v.nadj, v->vid);
  }

  bool gather_edges(const Vertex& v) const { return true; }

  Gather gather(const Vertex& v, Edge& e) const {
    return label_counter(v->label, v.nadj, v->vid);
  }

  /* Getting the most probable label for a node and making a decision to activate() a node for next iteration
     activate() for a node means that it does not meet the criterion we use for stopping, namely the threshold for L2 norm*/
  void apply(Vertex& v,  Gather& total) {
      
	 if(v->check_seed){
		 	 for(int i=0; i<number_of_groups; i++){
		 		 total.label_count[i]=0.0;
		 	 }
		 	 changelabel.nadj=v.nadj;
			 changelabel.label_count[v->label]=1.0;
		 return;
	 }
    double maxCount = 0.0;
    label_counter tmp;
    for(int i=0; i<number_of_groups; i++){
    	tmp.label_count[i]=changelabel.label_count[i];
    	changelabel.label_count[i]=total.label_count[i];
    	total.label_count[i]=0;
    }
    changelabel.normalize();
    double norm=0;
    for(int i=0; i<number_of_groups; i++){
    	norm+=(tmp.label_count[i]-changelabel.label_count[i])*(tmp.label_count[i]-changelabel.label_count[i]);
    }
    if(sqrt(norm)>=threshold) v->activate();
    int maxLabel = v->label;
    changelabel.vId=v->vid;

    for (int i = 0; i < number_of_groups; i++) {
        if (changelabel.label_count[i] > maxCount) {
          maxCount = changelabel.label_count[i];
          maxLabel = i;
        }
    }
	if (maxLabel != v->label){
		if(maxCount>0){
				v->label = maxLabel;
		}
	}
  }

  //Decides if Scatter should happen from the node or not
  bool scatter_edges(const Vertex& v) {
	  if(v->label==Unknown) return false;
	  return true;
  }

  Gather scatter(const Edge& e, Vertex& target) const {
    return changelabel;
  }
};

// loading seed nodes with corresponing label from the file

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
      // Load seed file
      LOG(INFO) << "start loading seed file";
      {
      symmetric_static std::map<int64_t, int> symmetric_seeds;
          on_all_cores( []{
              symmetric_seeds = loadSeed(FLAGS_seedsPath);
            });

       forall(g, [](VertexID i, G::Vertex& v){
                    v->vid= i;
                    v->iter=0;
                  });
        LOG(INFO) << "Init the labels";
        LOG(INFO) << "Size of the graph is: " << g->nv;

        //setting the labels for seed nodes to use them during propagation

        forall(g, [](VertexID i, G::Vertex& v){
          if (symmetric_seeds.find(i) == symmetric_seeds.end())
          {
            v->label = -1.0;
            v->check_seed=false;
            v->iter=0;
            v->activate();
          } else {
            v->label = symmetric_seeds[i];
            v->check_seed=true;
            v->activate();
            v->iter=0;
          }
        });
      }

      GRAPPA_TIME_REGION(total_time) {
        LOG(INFO) << "Init is complete";
        //starting propagation
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
 
