//
//  nis_tracker.h
//  UnscentedKF
//
//  Tracks NIS values for parameter tuning.
//

#ifndef NIS_TRACKER_H
#define NIS_TRACKER_H

#include <vector>
#include <numeric>
using namespace std;

class NisTracker {
public:
    NisTracker();
    void add(float nis_rating);
    float average();
private:
    vector<float> history_;
};



#endif /* NIS_TRACKER_H */
