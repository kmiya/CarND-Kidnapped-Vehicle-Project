/**
 * map.h
 *
 * Created on: Dec 12, 2016
 * Author: mufferm
 */

#ifndef MAP_H_
#define MAP_H_

#include <vector>
#include <unordered_map>

class Map {
 public:  
  struct single_landmark_s {
    int id_i ; // Landmark ID
    float x_f; // Landmark x-position in the map (global coordinates)
    float y_f; // Landmark y-position in the map (global coordinates)
  };

  std::vector<single_landmark_s> landmark_list; // List of landmarks in the map
  std::unordered_map<int, single_landmark_s> landmark_map;
};

#endif  // MAP_H_
