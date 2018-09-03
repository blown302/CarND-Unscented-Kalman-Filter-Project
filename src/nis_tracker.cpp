//
//  nis_tracker.cpp
//  UnscentedKF
//
//  tracks nis values
//

#include "nis_tracker.h"

NisTracker::NisTracker(double low_limit, double high_limit) {
    low_limit_ = low_limit;
    high_limit_ = high_limit;
    violation_count_ = 0;
    history_ = vector<double>();
}

void NisTracker::add(double nis_rating) {
    if (nis_rating < low_limit_ || nis_rating > high_limit_) {
        ++violation_count_;
    }

    history_.push_back(nis_rating);
}

double NisTracker::average() {
    return accumulate(history_.begin(), history_.end(), 0.0) / history_.size();
}

int NisTracker::get_violations() {
    return violation_count_;
}
