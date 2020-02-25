/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::default_random_engine;
using std::normal_distribution;
using std::uniform_real_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  default_random_engine gen;
  
  // set the standard deviations for x, y, and theta
  const double std_x = std[0];
  const double std_y = std[1];
  const double std_theta = std[2];
  
  //set the normal distributions for x, y, and theta
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);
  
  //initialize particles and set weight to 1
  Particle particle;
  particle.weight = 1.0;
  
  for (int i = 0; i < num_particles; ++i){
   particle.id = i; 
   particle.x = dist_x(gen); 
   particle.y = dist_y(gen); 
   particle.theta = dist_theta(gen); 
   particles.push_back(particle); 
   weights.push_back(particle.weight); 
  }

  is_initialized = true; 
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  default_random_engine gen; 
  
  double predict_x, predict_y, predict_theta; 
  
  for (int i = 0; i < num_particles; ++i){
    if (yaw_rate != 0){
      predict_x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)); 
      predict_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)); 
      predict_theta = particles[i].theta + yaw_rate * delta_t; 
    }
    else if (yaw_rate == 0){
      predict_x = particles[i].x + (velocity * delta_t) * (cos(particles[i].theta)); 
      predict_y = particles[i].y + (velocity * delta_t) * (sin(particles[i].theta)); 
      predict_theta = particles[i].theta;
    }
    //normal distribution for particle prediction x, y, and theta
    normal_distribution<double> pred_dist_x(predict_x, std_pos[0]); 
    normal_distribution<double> pred_dist_y(predict_y, std_pos[1]); 
    normal_distribution<double> pred_dist_theta(predict_theta, std_pos[2]); 
    
    // add gaussian noise to particle prediction x, y, and theta
    particles[i].x = pred_dist_x(gen); 
    particles[i].y = pred_dist_y(gen); 
    particles[i].theta = pred_dist_theta(gen); 
  }
 
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  //initialize id of predicted landmark
  int predicted_landmark_id;
  //initialize x, y, and current distance
  double delta_x = 0.0, delta_y = 0.0, dist = 0.0;
  
  //loop through sensor measurements and set the closest landmark to be far away to allow for looping through all landmarks and saving the closest
  for (unsigned int i = 0; i < observations.size(); ++i){ 
    double closest_obs_dist = 10000; 
    predicted_landmark_id = 0;
    
    //loop through landmarks to find the closest landmark by finding the range
    for (unsigned int j = 0; j < predicted.size(); ++j){
      delta_x = observations[i].x - predicted[j].x; 
      delta_y = observations[i].y - predicted[j].y; 
      dist = sqrt(pow(delta_x, 2.0) + pow(delta_y, 2.0)); 
      
      //if the range is less than the previous the id is saved
      if (dist < closest_obs_dist){
        closest_obs_dist = dist; 
        predicted_landmark_id = predicted[j].id; 
      }
  }
  // pair the closest landmark with sensor measurements   
  observations[i].id = predicted_landmark_id; 
    
 }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  double x_part, y_part, x_obs, y_obs, theta; 
  LandmarkObs t_observation; 
  LandmarkObs current_landmark; 
  double sig_x, sig_y, x_tobs, y_tobs, mu_x, mu_y, weight, final_weight; 
  unsigned int tobs_landmark_id, pred_landmark_id; 
  sig_x = std_landmark[0]; 
  sig_y = std_landmark[1]; 
  
  for (int i = 0; i < num_particles; ++i){
    x_part = particles[i].x; 
    y_part = particles[i].y; 
    theta = particles[i].theta; 
    
    //coordinate transform from vehicle to map coordinates 
    vector<LandmarkObs> transformed_observations; 
    for (unsigned int j = 0; j < observations.size(); ++j){
      x_obs = observations[j].x; 
      y_obs = observations[j].y;
      
      t_observation.id = observations[j].id; 
      t_observation.x = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs); 
      t_observation.y = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs); 
      transformed_observations.push_back(t_observation); 
    }
    
    //retain landmarks within range of the current particle
    vector<LandmarkObs> predicted_landmarks; 
    for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); ++k){
      current_landmark.x = map_landmarks.landmark_list[k].x_f; 
      current_landmark.y = map_landmarks.landmark_list[k].y_f; 
      current_landmark.id = map_landmarks.landmark_list[k].id_i; 
      double landmark_dist = dist(current_landmark.x, current_landmark.y, x_part, y_part); 
      if (landmark_dist <= sensor_range){
        predicted_landmarks.push_back(current_landmark); 
      }
    }
    
    //find the predicted landmarks that are closest to the transformed observations
    dataAssociation(predicted_landmarks, transformed_observations);
    
    //use multivariate gaussian distribution to find the particle weight
    final_weight = 1.0; 
    
    for (unsigned int m = 0; m < transformed_observations.size(); ++m){ 
      x_tobs = transformed_observations[m].x; 
      y_tobs = transformed_observations[m].y; 
      tobs_landmark_id = transformed_observations[m].id; 
      
      for (unsigned int n = 0; n < predicted_landmarks.size(); ++n){
        mu_x = predicted_landmarks[n].x; 
        mu_y = predicted_landmarks[n].y; 
        pred_landmark_id = predicted_landmarks[n].id; 
        
        if (tobs_landmark_id == pred_landmark_id){
          weight = multiv_prob(sig_x, sig_y, x_tobs, y_tobs, mu_x, mu_y); 
          if (weight != 0){
            final_weight *= weight; 
          }
        }
      }
    }
    particles[i].weight = final_weight; 
    weights[i] = particles[i].weight; 
  }   
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::discrete_distribution<> dist(weights.begin(), weights.end());
  std::vector <Particle> resampled_particles(particles.size());
  std::default_random_engine gen;
  
  for (int i = 0; i < num_particles; ++i){
    resampled_particles[i] = particles[dist(gen)];
  }
  particles = std::move(resampled_particles); 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}