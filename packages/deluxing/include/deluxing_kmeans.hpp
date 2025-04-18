/*
 * Copyright (c) 2025 Yu Yang & Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

/*
 * @file deluxing_kmeans.hpp
 * @brief Header file for the k-means clustering algorithm.
 *
 * This header defines the k-means clustering algorithm, which is used to partition a std::set of data points
 * into k clusters. The algorithm is generic across the data type and the dimension of each data item.
 */

#ifndef ROUTE_OPT_DELUXING_KMEANS_HPP
#define ROUTE_OPT_DELUXING_KMEANS_HPP


#include <cassert>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>


namespace RouteOpt::DeLuxing::KMeans {
    /*
    These functions are all private implementation details and shouldn't be referenced outside of this
    file.
    */
    namespace details {
        /*
        Calculate the square of the distance between two points.
        */
        template
        <typename T>
        T distance_squared(const std::vector<T> &point_a, const std::vector<T> &point_b) {
            T d_squared = T();
            for (typename std::vector<T>::size_type i = 0; i < point_a.size(); ++i) {
                auto delta = point_a[i] - point_b[i];
                d_squared += delta * delta;
            }
            return d_squared;
        }

        template
        <typename T>
        T distance(const std::vector<T> &point_a, const std::vector<T> &point_b) {
            return std::sqrt(distance_squared(point_a, point_b));
        }


        /*
        Calculate the smallest distance between each of the data points and any of the input means.
        */
        template
        <typename T>
        void closest_distance(
            const std::vector<std::vector<T> > &means, const std::vector<std::vector<T> > &data,
            std::vector<T> &distances) {
            int it = 0;
            for (auto &d: data) {
                T closest = distance_squared(d, means[0]);
                for (auto &m: means) {
                    T distance = distance_squared(d, m);
                    if (distance < closest)
                        closest = distance;
                }
                distances[it++];
            }
        }

        /*
        This is an alternate initialization method based on the [kmeans++](https://en.wikipedia.org/wiki/K-means%2B%2B)
        initialization algorithm.
        */
        template
        <typename T>
        void random_plusplus(const std::vector<std::vector<T> > &data, std::vector<std::vector<T> > &means,
                             uint32_t k, uint64_t seed) {
            assert(k > 0);
            assert(data.size() > 0);
            using input_size_t = typename std::vector<T>::size_type;
            // Using a very simple PRBS generator, parameters selected according to
            // https://en.wikipedia.org/wiki/Linear_congruential_generator#Parameters_in_common_use
            std::linear_congruential_engine<uint64_t, 6364136223846793005, 1442695040888963407, UINT64_MAX>
                    rand_engine(seed);

            // Select first mean at random from the std::set
            {
                std::uniform_int_distribution<input_size_t> uniform_generator(0, data.size() - 1);
                means.push_back(data[uniform_generator(rand_engine)]);
            }

            std::vector<T> distances(data.size());
            for (uint32_t count = 1; count < k; ++count) {
                // Calculate the distance to the closest mean for each data point
                details::closest_distance(means, data, distances);
                // Pick a random point weighted by the distance from existing means
                // TODO: This might convert floating point weights to ints, distorting the distribution for small weights
#if !defined(_MSC_VER) || _MSC_VER >= 1900
                std::discrete_distribution<input_size_t> generator(distances.begin(), distances.end());
#else  // MSVC++ older than 14.0
		input_size_t i = 0;
		std::discrete_distribution<input_size_t> generator(distances.size(), 0.0, 0.0, [&distances, &i](double) { return distances[i++]; });
#endif
                means.push_back(data[generator(rand_engine)]);
            }
        }

        /*
        Calculate the index of the mean a particular data point is closest to (euclidean distance)
        */
        template
        <typename T>
        uint32_t closest_mean(const std::vector<T> &point, const std::vector<std::vector<T> > &means) {
            assert(!means.empty());
            T smallest_distance = distance_squared(point, means[0]);
            typename std::vector<T>::size_type index = 0;
            T distance;
            for (size_t i = 1; i < means.size(); ++i) {
                distance = distance_squared(point, means[i]);
                if (distance < smallest_distance) {
                    smallest_distance = distance;
                    index = i;
                }
            }
            return index;
        }

        /*
        Calculate the index of the mean each data point is closest to (euclidean distance).
        */
        template
        <typename T>
        void calculate_clusters(
            const std::vector<std::vector<T> > &data,
            std::vector<uint32_t> &clusters, const std::vector<std::vector<T> > &means) {
            int it = 0;
            for (auto &point: data) {
                clusters[it++] = closest_mean(point, means);
            }
        }

        /*
        Calculate means based on data points and their cluster assignments.
        */
        template
        <typename T>
        void calculate_means(const std::vector<std::vector<T> > &data,
                             std::vector<std::vector<T> > &means,
                             const std::vector<uint32_t> &clusters,
                             const std::vector<std::vector<T> > &old_means,
                             uint32_t k) {
            for (auto &it: means) std::fill(it.begin(), it.end(), 0.);
            std::vector<T> count(k, T());
            for (size_t i = 0; i < std::min(clusters.size(), data.size()); ++i) {
                auto &mean = means[clusters[i]];
                count[clusters[i]] += 1;
                for (size_t j = 0; j < std::min(data[i].size(), mean.size()); ++j) {
                    mean[j] += data[i][j];
                }
            }
            for (size_t i = 0; i < k; ++i) {
                if (count[i] == 0) {
                    means[i] = old_means[i];
                } else {
                    for (size_t j = 0; j < means[i].size(); ++j) {
                        means[i][j] /= count[i];
                    }
                }
            }
        }

        template
        <typename T>
        std::vector<T> deltas(
            const std::vector<std::vector<T> > &old_means, const std::vector<std::vector<T> > &means) {
            std::vector<T> distances;
            distances.reserve(means.size());
            assert(old_means.size() == means.size());
            for (size_t i = 0; i < means.size(); ++i) {
                distances.push_back(distance(means[i], old_means[i]));
            }
            return distances;
        }

        template
        <typename T>
        bool deltas_below_limit(const std::vector<T> &deltas, T min_delta) {
            for (T d: deltas) {
                if (d > min_delta) {
                    return false;
                }
            }
            return true;
        }

        /*
        Calculate the smallest distance between each of the data points and any of the input means.
        */
        template
        <typename T>
        std::vector<T> closest_distance_parallel(
            const std::vector<std::vector<T> > &means, const std::vector<std::vector<T> > &data) {
            std::vector<T> distances(data.size(), T());
#pragma omp parallel for
            for (size_t i = 0; i < data.size(); ++i) {
                T closest = distance_squared(data[i], means[0]);
                for (const auto &m: means) {
                    T distance = distance_squared(data[i], m);
                    if (distance < closest)
                        closest = distance;
                }
                distances[i] = closest;
            }
            return distances;
        }

        /*
        This is an alternate initialization method based on the [kmeans++](https://en.wikipedia.org/wiki/K-means%2B%2B)
        initialization algorithm.
        */
        template
        <
            typename T>


        void random_plusplus_parallel(const std::vector<std::vector<T> > &data, std::vector<std::vector<T> > &means,
                                      uint32_t k, uint64_t seed) {
            assert(k > 0);
            assert(data.size() > 0);
            using input_size_t = typename std::vector<T>::size_type;
            // Using a very simple PRBS generator, parameters selected according to
            // https://en.wikipedia.org/wiki/Linear_congruential_generator#Parameters_in_common_use
            std::linear_congruential_engine<uint64_t, 6364136223846793005, 1442695040888963407, UINT64_MAX>
                    rand_engine(seed);

            // Select first mean at random from the std::set
            {
                std::uniform_int_distribution<input_size_t> uniform_generator(0, data.size() - 1);
                means.push_back(data[uniform_generator(rand_engine)]);
            }

            for (uint32_t count = 1; count < k; ++count) {
                // Calculate the distance to the closest mean for each data point
                auto distances = details::closest_distance_parallel(means, data);
                // Pick a random point weighted by the distance from existing means
                // TODO: This might convert floating point weights to ints, distorting the distribution for small weights
#if !defined(_MSC_VER) || _MSC_VER >= 1900
                std::discrete_distribution<input_size_t> generator(distances.begin(), distances.end());
#else  // MSVC++ older than 14.0
		input_size_t i = 0;
		std::discrete_distribution<input_size_t> generator(distances.size(), 0.0, 0.0, [&distances, &i](double) { return distances[i++]; });
#endif
                means.push_back(data[generator(rand_engine)]);
            }
        }

        /*
        Calculate the index of the mean each data point is closest to (euclidean distance).
        */
        template
        <
            typename T>


        void calculate_clusters_parallel(
            const std::vector<std::vector<T> > &data,
            std::vector<uint32_t> &clusters,
            const std::vector<std::vector<T> > &means) {
#pragma omp parallel for
            for (size_t i = 0; i < data.size(); ++i) {
                clusters[i] = closest_mean(data[i], means);
            }
        }
    } // namespace details

    template
    <typename T>
    class clustering_parameters {
    public:
        explicit clustering_parameters(uint32_t k) : _k(k),
                                                     _has_max_iter(false), _max_iter(),
                                                     _has_min_delta(false), _min_delta(),
                                                     _has_rand_seed(false), _rand_seed() {
        }

        void set_max_iteration(uint64_t max_iter) {
            _max_iter = max_iter;
            _has_max_iter = true;
        }

        void set_min_delta(T min_delta) {
            _min_delta = min_delta;
            _has_min_delta = true;
        }

        void set_random_seed(uint64_t rand_seed) {
            _rand_seed = rand_seed;
            _has_rand_seed = true;
        }

        bool has_max_iteration() const { return _has_max_iter; }
        bool has_min_delta() const { return _has_min_delta; }
        bool has_random_seed() const { return _has_rand_seed; }

        uint32_t get_k() const { return _k; };
        uint64_t get_max_iteration() const { return _max_iter; }
        T get_min_delta() const { return _min_delta; }
        uint64_t get_random_seed() const { return _rand_seed; }

    private:
        uint32_t _k;
        bool _has_max_iter;
        uint64_t _max_iter;
        bool _has_min_delta;
        T _min_delta;
        bool _has_rand_seed;
        uint64_t _rand_seed;
    };

    /*
    clustering_parameters is the configuration used for running the kmeans_lloyd algorithm.

    It requires a k value for initialization, and can subsequently be configured with your choice
    of optional parameters, including:
    * Maximum iteration count; the algorithm will terminate if it reaches this iteration count
      before converging on a solution. The results returned are the means and cluster assignments
      calculated in the last iteration before termination.
    * Minimum delta; the algorithm will terminate if the change in position of all means is
      smaller than the specified distance.
    * Random seed; if present, this will be used in place of `std::random_device` for kmeans++
      initialization. This can be used to ensure reproducible/deterministic behavior.
    */

    /*
    Implementation of k-means generic across the data type and the dimension of each data item. Expects
    the data to be a std::vector of fixed-size arrays. Generic parameters are the type of the base data (T)
    and the dimensionality of each data point (N). All points must have the same dimensionality.

    e.g. points of the form (X, Y, Z) would be N = 3.

    Takes a `clustering_parameters` struct for algorithm configuration. See the comments for the
    `clustering_parameters` struct for more information about the configuration values and how they
    affect the algorithm.

    Returns a std::tuple containing:
      0: A std::vector holding the means for each cluster from 0 to k-1.
      1: A std::vector containing the cluster number (0 to k-1) for each corresponding element of the input
         data std::vector.

    Implementation details:
    This implementation of k-means uses [Lloyd's Algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
    with the [kmeans++](https://en.wikipedia.org/wiki/K-means%2B%2B)
    used for initializing the means.

    */
    template
    <typename T>
    void kmeans_lloyd(const std::vector<std::vector<T> > &data,
                      std::vector<std::vector<T> > &means,
                      std::vector<uint32_t> &clusters,
                      const clustering_parameters<T> &parameters) {
        static_assert(std::is_arithmetic<T>::value && std::is_signed<T>::value,
                      "kmeans_lloyd requires the template parameter T to be a signed arithmetic type (e.g. float, double, int)")
                ;
        assert(parameters.get_k() > 0); // k must be greater than zero
        assert(data.size() >= parameters.get_k()); // there must be at least k data points
        std::random_device rand_device;
        uint64_t seed = parameters.has_random_seed() ? parameters.get_random_seed() : rand_device();
        details::random_plusplus(data, means, parameters.get_k(), seed);

        std::vector<std::vector<T> > old_means;
        std::vector<std::vector<T> > old_old_means;
        // Calculate new means until convergence is reached or we hit the maximum iteration count
        uint64_t count = 0;
        do {
            details::calculate_clusters(data, clusters, means);
            old_old_means = old_means;
            old_means = means;
            details::calculate_means(data, means, clusters, old_means, parameters.get_k());
            ++count;
        } while (means != old_means && means != old_old_means
                 && !(parameters.has_max_iteration() && count == parameters.get_max_iteration())
                 && !(parameters.has_min_delta() && details::deltas_below_limit(
                          details::deltas(old_means, means), parameters.get_min_delta())));
    }

    /*
    This overload exists to support legacy code which uses this signature of the kmeans_lloyd function.
    Any code still using this signature should move to the version of this function that uses a
    `clustering_parameters` struct for configuration.
    */
    template
    <typename T>
    void kmeans_lloyd(
        const std::vector<std::vector<T> > &data,
        std::vector<std::vector<T> > &means, std::vector<uint32_t> &clusters,
        uint32_t k,
        uint64_t max_iter, T min_delta) {
        clustering_parameters<T> parameters(k);
        if (max_iter != 0) {
            parameters.set_max_iteration(max_iter);
        }
        if (min_delta != 0) {
            parameters.set_min_delta(min_delta);
        }
        kmeans_lloyd(data, means, clusters, parameters);
    }

    /*
    Implementation of k-means generic across the data type and the dimension of each data item. Expects
    the data to be a std::vector of fixed-size arrays. Generic parameters are the type of the base data (T)
    and the dimensionality of each data point (N). All points must have the same dimensionality.

    e.g. points of the form (X, Y, Z) would be N = 3.

    Returns a std::tuple containing:
      0: A std::vector holding the means for each cluster from 0 to k-1.
      1: A std::vector containing the cluster number (0 to k-1) for each corresponding element of the input
         data std::vector.

    Implementation details:
    This implementation of k-means uses [Lloyd's Algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
    with the [kmeans++](https://en.wikipedia.org/wiki/K-means%2B%2B)
    used for initializing the means.
    */
    template
    <typename T>
    void kmeans_lloyd_parallel(
        const std::vector<std::vector<T> > &data,
        std::vector<std::vector<T> > &means, std::vector<uint32_t> &clusters,
        const clustering_parameters<T> &parameters) {
        static_assert(std::is_arithmetic<T>::value && std::is_signed<T>::value,
                      "kmeans_lloyd requires the template parameter T to be a signed arithmetic type (e.g. float, double, int)")
                ;
        assert(parameters.get_k() > 0); // k must be greater than zero
        assert(data.size() >= parameters.get_k()); // there must be at least k data points
        std::random_device rand_device;
        uint64_t seed = parameters.has_random_seed() ? parameters.get_random_seed() : rand_device();
        details::random_plusplus_parallel(data, means, parameters.get_k(), seed);

        std::vector<std::vector<T> > old_means;
        std::vector<std::vector<T> > old_old_means;
        // Calculate new means until convergence is reached or we hit the maximum iteration count
        uint64_t count = 0;
        do {
            details::calculate_clusters_parallel(data, clusters, means);
            old_old_means = old_means;
            old_means = means;
            details::calculate_means(data, means, clusters, old_means, parameters.get_k());
            ++count;
        } while ((means != old_means && means != old_old_means)
                 && !(parameters.has_max_iteration() && count == parameters.get_max_iteration())
                 && !(parameters.has_min_delta() && details::deltas_below_limit(
                          details::deltas(old_means, means), parameters.get_min_delta())));
    }

    /*
    This overload exists to support legacy code which uses this signature of the kmeans_lloyd function.
    Any code still using this signature should move to the version of this function that uses a
    `clustering_parameters` struct for configuration.
    */
    template
    <typename T>
    void kmeans_lloyd_parallel(
        const std::vector<std::vector<T> > &data,
        std::vector<std::vector<T> > &means, std::vector<uint32_t> &clusters,
        uint32_t k,
        uint64_t max_iter, T min_delta) {
        clustering_parameters<T> parameters(k);
        if (max_iter != 0) {
            parameters.set_max_iteration(max_iter);
        }
        if (min_delta != 0) {
            parameters.set_min_delta(min_delta);
        }
        kmeans_lloyd_parallel(data, means, clusters, parameters);
    }
}


#endif /* KMEANS_H */
