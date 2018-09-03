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
    NisTracker(double low_limit, double high_limit);
    void add(double nis_rating);
    double average();
    int get_violations();
private:
    int violation_count_;
    double high_limit_, low_limit_;
    vector<double> history_;
};

#endif /* NIS_TRACKER_H */
