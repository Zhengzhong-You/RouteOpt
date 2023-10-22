
#include "CVRP.hpp"
#include "VRPTW.hpp"

using namespace std;

void CVRP::setResourceInBucketGraph() {
  max_main_resource = cap;
  meet_point_resource_in_bi_dir = max_main_resource / 2;
  main_resource_across_arcs_in_forward_sense.resize(dim);
  for (auto &vertex : main_resource_across_arcs_in_forward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
#ifdef SYMMETRY_PROHIBIT
	  main_resource_across_arcs_in_forward_sense[i][j] = info_vertex[i][3];
#else
	  main_resource_across_arcs_in_forward_sense[i][j] = (info_vertex[i][3] + info_vertex[j][3]) / 2;
#endif
	}
  }
#ifdef SYMMETRY_PROHIBIT
  main_resource_across_arcs_in_backward_sense.resize(dim);
  for (auto &vertex : main_resource_across_arcs_in_backward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  main_resource_across_arcs_in_backward_sense[i][j] = info_vertex[j][3];
	}
  }
#endif
  lb4_vertex.resize(dim, 0);
  ub4_vertex.resize(dim, max_main_resource);
}

void CVRP::initializeBucketGraph() {
  label_array_in_forward_sense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	label_array_in_forward_sense[i] = new VecLabel[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  label_array_in_forward_sense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
	}
  }

  if_exist_extra_labels_in_forward_sense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	if_exist_extra_labels_in_forward_sense[i] = new VecLabel[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  if_exist_extra_labels_in_forward_sense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
	}
  }
  rc2_till_this_bin_in_forward_sense = new double *[dim];
  for (int i = 0; i < dim; ++i) {
	rc2_till_this_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
  }
  rc2_bin_in_forward_sense = new double *[dim];
  for (int i = 0; i < dim; ++i) {
	rc2_bin_in_forward_sense[i] = new double[num_buckets_per_vertex];
  }
#ifdef SYMMETRY_PROHIBIT
  label_array_in_backward_sense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	label_array_in_backward_sense[i] = new VecLabel[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  label_array_in_backward_sense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
	}
  }

  if_exist_extra_labels_in_backward_sense = new VecLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	if_exist_extra_labels_in_backward_sense[i] = new VecLabel[num_buckets_per_vertex];
	for (int j = 0; j < num_buckets_per_vertex; ++j) {
	  if_exist_extra_labels_in_backward_sense[i][j].first.resize(LABEL_LIMIT_PER_BIN);
	}
  }
  rc2_till_this_bin_in_backward_sense = new double *[dim];
  for (int i = 0; i < dim; ++i) {
	rc2_till_this_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
  }
  rc2_bin_in_backward_sense = new double *[dim];
  for (int i = 0; i < dim; ++i) {
	rc2_bin_in_backward_sense[i] = new double[num_buckets_per_vertex];
  }
#endif
}
