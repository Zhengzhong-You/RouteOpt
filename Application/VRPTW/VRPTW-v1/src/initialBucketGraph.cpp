
#include "CVRP.hpp"
#include "VRPTW.hpp"

using namespace std;

void CVRP::setResourceInBucketGraph() {
  resource.first_res = roundAndConvertResLong(cap);
  meet_point_resource_in_bi_dir = double(resource.first_res) / 2;
  resource_across_arcs_in_forward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_forward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
#ifdef SYMMETRY_PROHIBIT
	  resource_across_arcs_in_forward_sense[i][j].first_res = roundAndConvertResLong(info_vertex[i][3]);
#else
	  resource_across_arcs_in_forward_sense[i][j].first_res =
		  roundAndConvertResLong((info_vertex[i][3] + info_vertex[j][3]) / 2);
#endif
	}
  }
#ifdef SYMMETRY_PROHIBIT
  resource_across_arcs_in_backward_sense.resize(dim);
  for (auto &vertex : resource_across_arcs_in_backward_sense) vertex.resize(dim);
  for (int i = 0; i < dim; ++i) {
	for (int j = 0; j < dim; ++j) {
	  resource_across_arcs_in_backward_sense[i][j].first_res = roundAndConvertResLong(info_vertex[j][3]);
	}
  }
#endif
  lb4_vertex.resize(dim, {});
  ub4_vertex.resize(dim, {roundAndConvertResLong(cap)});
}

void CVRP::initializeBucketGraph() {
  label_array_in_forward_sense = new ListLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	label_array_in_forward_sense[i] = new ListLabel[num_buckets_per_vertex];
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
  label_array_in_backward_sense = new ListLabel *[dim];
  for (int i = 0; i < dim; ++i) {
	label_array_in_backward_sense[i] = new ListLabel[num_buckets_per_vertex];
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
