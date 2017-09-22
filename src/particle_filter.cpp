/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

 #include <random>
 #include <algorithm>
 #include <iostream>
 #include <numeric>
 #include <math.h> 
 #include <iostream>
 #include <sstream>
 #include <string>
 #include <iterator>
 
 #include "particle_filter.h"
 
 using namespace std;
 
 void ParticleFilter::init(double x, double y, double theta, double std[]) {
	 // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	 //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	 // Add random Gaussian noise to each particle.
	 // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	 
	 default_random_engine random_eng;
 
	 num_particles = 42;
 
	 // standard deviations for x, y, and theta
	 double std_x, std_y, std_theta;
	 std_x = std[0];
	 std_y = std[1];
	 std_theta = std[2];
 
	 normal_distribution<double> dist_x(x, std_x);
	 normal_distribution<double> dist_y(y, std_y);
	 normal_distribution<double> dist_theta(theta, std_theta);
 
	 for(int i = 0; i < num_particles; ++i){
 
		 Particle a_particle;
 
		 a_particle.id = i;
 
		 // sample from normal distrubtions using random engine
		 a_particle.x = dist_x(random_eng);
		 a_particle.y = dist_y(random_eng);
		 a_particle.theta = dist_theta(random_eng);
		 a_particle.weight = 1.0;
 
		 particles.push_back(a_particle);
 
		 weights.push_back(1.0);
	 }
 
	 is_initialized = true;
 }
 
 void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	 // TODO: Add measurements to each particle and add random Gaussian noise.
	 // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	 //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	 //  http://www.cplusplus.com/reference/random/default_random_engine/
 
	 default_random_engine random_eng;
 
	 // standard deviations for x, y, and theta
	 double std_x, std_y, std_theta;
	 
	 std_x = std_pos[0];
	 std_y = std_pos[1];
	 std_theta = std_pos[2];
 
	 // new values for x, y, and theta
	 double new_x, new_y, new_theta;
 
 
	 for(int i = 0; i < num_particles; ++i){
		 
		 if(fabs(yaw_rate) <= 0.00001){
			 new_x = particles[i].x + (velocity * delta_t * cos(particles[i].theta));
			 new_y = particles[i].y + (velocity * delta_t * sin(particles[i].theta));
			 new_theta = particles[i].theta;
		 }else{
			 /*
			 The equations for updating x, y and the yaw angle when the yaw rate is not equal to zero
			 x​f​​ = x​0​​+​​θ​˙​​​​v​​[sin(θ​0​​+​θ​˙​​(dt))−sin(θ​0​​)]
			 y​f ​​= y​0​​+​​θ​˙​​​​v​​[cos(θ​0​​)−cos(θ​0​​+​θ​˙​​(dt))]
			 θ​f​​ = θ​0 ​​+​ θ​˙​​(dt)
			 */
			 double tyd = particles[i].theta + (yaw_rate * delta_t);
			 
			 new_x = particles[i].x + (velocity/yaw_rate) * (sin(tyd) - sin(particles[i].theta));
			 new_y = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(tyd));
			 new_theta = tyd;
		 }
 
		 //add random Gaussian noise to account for sensor uncertainty
		 normal_distribution<double> dist_x(new_x, std_x);
		 normal_distribution<double> dist_y(new_y, std_y);
		 normal_distribution<double> dist_theta(new_theta, std_theta);
 
		 particles[i].x = dist_x(random_eng);
		 particles[i].y = dist_y(random_eng);
		 particles[i].theta = dist_theta(random_eng);
 
 
	 }//END: for(int i = 0; i < num_particles; i++)
 
 }
 
 void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	 // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	 //   observed measurement to this particular landmark.
	 // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	 //   implement this method and use it as a helper during the updateWeights phase.
	double min_dist;
	 
	for(int i = 0; i < observations.size(); i++){

		// init minimum distance to the maximum possible
		min_dist = numeric_limits<double>::max();

		for(unsigned int j = 0; j < predicted.size(); j++){
			
			// grab current prediction
			LandmarkObs predicted_obs = predicted[j];

			// calc distance between observation and predicted observation
			double distance = dist(observations[i].x, observations[i].y, predicted_obs.x, predicted_obs.y);

			// find the predicted landmark nearest the current observed landmark
			if(distance < min_dist){
				min_dist = distance;
				observations[i].id = predicted[j].id;
			}

		}//END: for(unsigned int j = 0; j < predicted.size(); j++)

	}//END: for(int i = 0; i < observations.size(); i++)
 }
 
 void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		 const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	 // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	 //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	 // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	 //   according to the MAP'S coordinate system. You will need to transform between the two systems.
	 //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	 //   The following is a good resource for the theory:
	 //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	 //   and the following is a good resource for the actual equation to implement (look at equation 
	 //   3.33
	 //   http://planning.cs.uiuc.edu/node99.html
 
	 // standard deivations for x, y, theta
	 double std_x = std_landmark[0];
	 double std_y = std_landmark[0];
 
	 double weight_c1 = (1/(2 * M_PI * std_x * std_y));

	 for(int i = 0; i < num_particles; ++i){
		 
		 double obs_x, obs_y;
 
		 double p_x = particles[i].x;
		 double p_y = particles[i].y;
		 double p_theta = particles[i].theta;
 
		 // transform car cordinates to map cordinates
		 vector<LandmarkObs> transformed_obs;
 
		 for(int j = 0; j < observations.size(); j++){
			 obs_x = observations[j].x;
			 obs_y = observations[j].y;
 
			 LandmarkObs transformed_lm_obs;
 
			 transformed_lm_obs.id = observations[j].id;
			 transformed_lm_obs.x = p_x + (cos(p_theta) * obs_x) - (sin(p_theta) * obs_y);
			 transformed_lm_obs.y = p_y + (sin(p_theta) * obs_x) + (cos(p_theta) * obs_y);
 
			 transformed_obs.push_back(transformed_lm_obs);
		 }
 
		 // find landmarks that are within sensor range/car location
		 vector<LandmarkObs> predicted_observations;
 
		 for(int j = 0; j < map_landmarks.landmark_list.size(); ++j){

			double distance = dist(p_x,p_y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
			 
			if(distance < sensor_range){
				predicted_observations.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});
			}
		 }
 
		 // calc nearest neighbour by performing landmark or data association 
		 dataAssociation(predicted_observations,transformed_obs);

		 vector<int> associations;
		 vector<double> sense_x;
		 vector<double> sense_y;
 
		 // calculate the particle final weight by the product of
		 // each measurement's Multivariate-Gaussian probability.
 
		 particles[i].weight = 1.0;
 
		 double particle_weight = 1.0;
 
		 for(int j=0; j < transformed_obs.size(); ++j){
 
			 for(int k = 0; k < predicted_observations.size(); k++){
			 
				 if(predicted_observations[k].id == transformed_obs[j].id){
					 double obs_x, obs_y, pred_x, pred_y;
					 
					 obs_x = transformed_obs[j].x;
					 obs_y = transformed_obs[j].y;
 
					 pred_x = predicted_observations[k].x;
					 pred_y = predicted_observations[k].y;
 
					 double c2 = pow(obs_x - pred_x, 2) / pow(std_x, 2);
					 double c3 = pow(obs_y - pred_y, 2) / pow(std_y, 2);
					 double weight = weight_c1 * exp(-0.5 * (c2 + c3));
		   
					 if (weight < 0.0001){
						 weight = 0.0001;
					 }
					 particle_weight *= weight;
 
					 //break once multi variate gaussian weight calculated for predicted landmark found 
					 break;
				 }
			 }
			
			 associations.push_back(transformed_obs[j].id);
			 sense_x.push_back(transformed_obs[j].x);
			 sense_y.push_back(transformed_obs[j].y);

		 }//END: for(int j=0; j < transformed_obs.size(); j++)
 
		 // update particle weight with new weight
		 particles[i].weight = particle_weight;

		 particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
		 weights[i] = particles[i].weight;
 
	 }//END: for(int i = 0; i < num_particles; i++)
 }
 
 void ParticleFilter::resample() {
	 // TODO: Resample particles with replacement with probability proportional to their weight. 
	 // NOTE: You may find std::discrete_distribution helpful here.
	 //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
 
	 default_random_engine random_eng;
 
	 discrete_distribution<int> dist_weights(weights.begin(), weights.end());

	 vector<Particle> resampled_particles;
 
	 for(int i = 0; i < num_particles; i++){
		 resampled_particles.push_back(particles[dist_weights(random_eng)]);
	 }
 
	 particles = resampled_particles;
 
 }
 
 Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
 {
	 //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	 // associations: The landmark id that goes along with each listed association
	 // sense_x: the associations x mapping already converted to world coordinates
	 // sense_y: the associations y mapping already converted to world coordinates
 
	 //Clear the previous associations
	 particle.associations.clear();
	 particle.sense_x.clear();
	 particle.sense_y.clear();
 
	 particle.associations= associations;
	  particle.sense_x = sense_x;
	  particle.sense_y = sense_y;
 
	  return particle;
 }
 
 string ParticleFilter::getAssociations(Particle best)
 {
	 vector<int> v = best.associations;
	 stringstream ss;
	 copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	 string s = ss.str();
	 s = s.substr(0, s.length()-1);  // get rid of the trailing space
	 return s;
 }
 string ParticleFilter::getSenseX(Particle best)
 {
	 vector<double> v = best.sense_x;
	 stringstream ss;
	 copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	 string s = ss.str();
	 s = s.substr(0, s.length()-1);  // get rid of the trailing space
	 return s;
 }
 string ParticleFilter::getSenseY(Particle best)
 {
	 vector<double> v = best.sense_y;
	 stringstream ss;
	 copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	 string s = ss.str();
	 s = s.substr(0, s.length()-1);  // get rid of the trailing space
	 return s;
 }
 