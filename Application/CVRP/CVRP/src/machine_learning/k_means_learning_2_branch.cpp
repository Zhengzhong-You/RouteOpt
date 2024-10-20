#include "machine_learning.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <numeric>

using namespace std;

void kMeans(const vector<vector<double>> &data,
			int k,
			vector<int> &labels,
			vector<vector<double>> &centroids,
			int seed = ML_RANDOM_SEED);
void kMeansPlusPlus(const vector<vector<double>> &data, int k, vector<vector<double>> &centroids, int seed);
double calculateWCSS(const vector<vector<double>> &data,
					 const vector<int> &labels,
					 const vector<vector<double>> &centroids);
int findOptimalCluster(const vector<vector<double>> &data, int maxK, int seed = ML_RANDOM_SEED);

double calculateDistance(const vector<double> &a, const vector<double> &b) {
  return sqrt(pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
}

void kMeansPlusPlus(const vector<vector<double>> &data, int k, vector<vector<double>> &centroids, int seed) {
  int n = (int)data.size();
  centroids.resize(k);

  mt19937 gen(seed);
  uniform_int_distribution<> dis(0, n - 1);

  centroids[0] = data[dis(gen)];

  vector<double> minDistances(n, numeric_limits<double>::max());

  for (int i = 1; i < k; ++i) {
	double totalDist = 0.0;

	for (int j = 0; j < n; ++j) {
	  double dist = calculateDistance(data[j], centroids[i - 1]);
	  if (dist < minDistances[j]) {
		minDistances[j] = dist;
	  }
	  totalDist += minDistances[j];
	}

	double target = uniform_real_distribution<>(0, totalDist)(gen);
	for (int j = 0; j < n; ++j) {
	  if ((target -= minDistances[j]) <= 0) {
		centroids[i] = data[j];
		break;
	  }
	}
  }
}

void kMeans(const vector<vector<double>> &data,
			int k,
			vector<int> &labels,
			vector<vector<double>> &centroids,
			int seed) {
  int n = (int)data.size();
  const int dim = 3;
  labels.assign(n, -1);
  centroids.resize(k, vector<double>(dim, 0.0));

  kMeansPlusPlus(data, k, centroids, seed);

  bool changed = true;
  while (changed) {
	changed = false;

	for (int i = 0; i < n; ++i) {
	  double minDist = numeric_limits<double>::max();
	  int bestCluster = -1;

	  for (int j = 0; j < k; ++j) {
		double dist = calculateDistance(data[i], centroids[j]);
		if (dist < minDist) {
		  minDist = dist;
		  bestCluster = j;
		}
	  }

	  if (labels[i] != bestCluster) {
		labels[i] = bestCluster;
		changed = true;
	  }
	}

	vector<vector<double>> newCentroids(k, vector<double>(dim, 0.0));
	vector<int> counts(k, 0);

	for (int i = 0; i < n; ++i) {
	  int cluster = labels[i];
	  for (int j = 1; j <= 2; ++j) {
		newCentroids[cluster][j] += data[i][j];
	  }
	  counts[cluster]++;
	}

	for (int j = 0; j < k; ++j) {
	  if (counts[j] > 0) {
		for (int d = 1; d <= 2; ++d) {
		  newCentroids[j][d] /= counts[j];
		}
	  }
	}

	centroids = newCentroids;
  }
}

double calculateWCSS(const vector<vector<double>> &data,
					 const vector<int> &labels,
					 const vector<vector<double>> &centroids) {
  double wcss = 0.0;
  for (int i = 0; i < data.size(); ++i) {
	int cluster = labels[i];
	wcss += pow(calculateDistance(data[i], centroids[cluster]), 2);
  }
  return wcss;
}

int findOptimalCluster(const vector<vector<double>> &data, int maxK, int seed) {
  vector<double> wcssValues;

  for (int k = 1; k <= maxK; ++k) {
	vector<int> labels;
	vector<vector<double>> centroids;
	kMeans(data, k, labels, centroids, seed);
	double wcss = calculateWCSS(data, labels, centroids);
	wcssValues.emplace_back(wcss);
  }

  int optimalCluster = 1;
  double maxDelta = 0.0;

  for (int k = 1; k < (int)wcssValues.size() - 1; ++k) {
	double delta = wcssValues[k - 1] - wcssValues[k];
	if (delta > maxDelta) {
	  maxDelta = delta;
	  optimalCluster = k + 1;
	}
  }

  return optimalCluster;
}

void MachineLearning::calculateClusteringCoefficient() {
  vector<double> k_vec(MAX_NUM_RUN_OPTIMAL_K);
  for (int i = 0; i < MAX_NUM_RUN_OPTIMAL_K; ++i) {
	k_vec[i] = findOptimalCluster(cvrp->info_vertex, MAX_NUM_CLUSTER, ML_RANDOM_SEED);
  }
  int optimalCluster = (int)std::round(accumulate(k_vec.begin(), k_vec.end(), 0.0) / MAX_NUM_RUN_OPTIMAL_K);
  verbose_call(cout << "optimalCluster: " << optimalCluster << endl;)

  vector<int> labels;
  vector<vector<double>> centroids;
  kMeans(cvrp->info_vertex, optimalCluster, labels, centroids, ML_RANDOM_SEED);

  unordered_map<int, vector<int>> cluster_map;
  for (int i = 0; i < cvrp->info_vertex.size(); ++i) {
	cluster_map[labels[i]].emplace_back(i);
  }

  double sum_dis_in_cluster = 0;
  int cnt = 0;
  for (auto &pr : cluster_map) {
	auto &vec = pr.second;
	double dis_in = 0;
	if (vec.size() <= 1) continue;
	++cnt;
	for (int i = 0; i < vec.size(); ++i) {
	  for (int j = i + 1; j < vec.size(); ++j) {
		dis_in += calculateDistance(cvrp->info_vertex[vec[i]], cvrp->info_vertex[vec[j]]);
	  }
	}
	dis_in /= double(vec.size() * int(vec.size() - 1)) / 2;
	sum_dis_in_cluster += dis_in;
  }
  sum_dis_in_cluster /= cnt * max_edge_cost;

  cluster_coeff = 1 / sum_dis_in_cluster;
}
