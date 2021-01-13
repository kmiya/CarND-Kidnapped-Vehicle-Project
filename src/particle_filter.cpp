/**
 * particle_filter.cpp
 *
 * Updated on: Jan 12, 2021
 * Author: Kazuki Miyahara
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(const double x, const double y, const double theta, const double std[]) {
  /**
   * Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * Add random Gaussian noise to each particle.
   */
  std::random_device seed_gen;
  random_engine_ = std::default_random_engine(seed_gen());
  std::normal_distribution<> dist_x(x, std[0]);
  std::normal_distribution<> dist_y(y, std[1]);
  std::normal_distribution<> dist_yaw(theta, std[2]);

  num_particles_ = 10;
  particles_.reserve(num_particles_);
  for (size_t i = 0; i < num_particles_; ++i) {
    Particle p;
    p.id = i;
    p.x = dist_x(random_engine_);
    p.y = dist_y(random_engine_);
    p.theta = dist_yaw(random_engine_);
    p.weight = 1.;
    particles_.emplace_back(p);
  }
  is_initialized_ = true;
}

void ParticleFilter::prediction(const double delta_t, const double std_pos[],
                                const double velocity, const double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   */
  std::normal_distribution<> dist_x(0., std_pos[0]);
  std::normal_distribution<> dist_y(0., std_pos[1]);
  std::normal_distribution<> dist_yaw(0., std_pos[2]);
  for (auto &p : particles_) {
    // update using motion model
    if (std::abs(yaw_rate) > eps_) {
      const double d_theta = p.theta + yaw_rate * delta_t;
      p.x += velocity * (sin(d_theta) - sin(p.theta)) / yaw_rate;
      p.y += velocity * (cos(p.theta) - cos(d_theta)) / yaw_rate;
      p.theta += yaw_rate * delta_t;
    } else {  // prevent division by zero
      p.x += velocity * delta_t * cos(p.theta);
      p.y += velocity * delta_t * sin(p.theta);
    }

    // add noise
    p.x += dist_x(random_engine_);
    p.y += dist_y(random_engine_);
    p.theta += dist_yaw(random_engine_);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> &predicted,
                                     const std::unordered_map<int, Map::single_landmark_s> &map_landmarks) {
  /**
   * Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   */
  for (auto &pred : predicted) {
    int nearest_obj_id = -1;
    double min_distance = std::numeric_limits<double>::max();
    for (const auto &obj : map_landmarks) {
      const double distance = dist(pred.x, pred.y, obj.second.x_f, obj.second.y_f);
      if (distance >= min_distance)
        continue;
      min_distance = distance;
      nearest_obj_id = obj.second.id_i;
    }
    pred.id = nearest_obj_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, const double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a multi-variate Gaussian
   *   distribution.
   */
  for (auto &p : particles_) {
    vector<LandmarkObs> predicted;
    predicted.reserve(observations.size());
    for (const auto &obs : observations) {
      LandmarkObs pred{};
      pred.id = -1;
      // convert to map coordinates
      pred.x = p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);
      pred.y = p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);
      predicted.emplace_back(pred);
    }
    // filter landmarks using sensor_range
    std::unordered_map<int, Map::single_landmark_s> marks;
    for (const auto &mark : map_landmarks.landmark_list) {
      const double distance = dist(p.x, p.y, mark.x_f, mark.y_f);
      if (distance <= sensor_range)
        marks.emplace(mark.id_i, mark);
    }
    dataAssociation(predicted, marks);
    p.weight = 1.;
    for (const auto &pred : predicted) {
      const double map_x = map_landmarks.landmark_map.at(pred.id).x_f;
      const double map_y = map_landmarks.landmark_map.at(pred.id).y_f;
      const double std_x = std_landmark[0];
      const double std_y = std_landmark[1];
      p.weight *= multiv_prob(std_x, std_y, pred.x, pred.y, map_x, map_y);
    }
  }
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                                   double mu_x, double mu_y) const {
  // calculate normalization term
  const double gauss_norm = 1. / (2. * M_PI * sig_x * sig_y);

  // calculate exponent
  const double exponent = (pow(x_obs - mu_x, 2) / (2. * pow(sig_x, 2)))
      + (pow(y_obs - mu_y, 2) / (2. * pow(sig_y, 2)));

  // calculate weight using normalization terms and exponent
  return gauss_norm * exp(-exponent);
}

void ParticleFilter::resample() {
  /**
   * Resample particles with replacement with probability proportional
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::vector<Particle> new_particles(num_particles_);
  vector<double> weights(num_particles_, 0.);
  for (size_t i = 0; i < num_particles_; ++i) {
    weights[i] = particles_[i].weight;
  }
  std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
  for (auto &np : new_particles) {
    const auto resample_idx = dist(random_engine_);
    np = particles_[resample_idx];
  }
  particles_ = new_particles;
}

void ParticleFilter::SetAssociations(Particle &particle,
                                     const vector<int> &associations,
                                     const vector<double> &sense_x,
                                     const vector<double> &sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(const Particle &best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(const Particle &best, const string &coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}