//
//  nis_tracker.cpp
//  UnscentedKF
//
//  tracks nis values
//

#include "nis_tracker.h"

NisTracker::NisTracker() {
    history_ = vector<float>();
}

void NisTracker::add(float nis_rating) {
    history_.push_back(nis_rating);
}

float NisTracker::average() {
    return accumulate(history_.begin(), history_.end(), 0.0) / history_.size();
}
