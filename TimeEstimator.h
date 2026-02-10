#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>

/**
 * Class for estimating remaining time of a long-running process
 */
class TimeEstimator {
public:
    /**
     * Constructor - initializes the estimator with total work units
     * @param total_work Total number of work units to be processed
     */
    TimeEstimator(size_t total_work)
        : total_work_(total_work), work_done_(0) {
        // Record the starting time point
        start_time_ = std::chrono::steady_clock::now();
    }

    /**
     * Update the progress with current completed work units
     * @param work_done Number of work units completed so far
     */
    void update(size_t work_done) {
        work_done_ = work_done;
        auto now = std::chrono::steady_clock::now();
        // Calculate elapsed time in milliseconds
        elapsed_time_ = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_);
    }

    /**
     * Calculate estimated remaining time
     * @return Remaining time in milliseconds
     */
    double get_remaining_time() const {
        if (work_done_ == 0) return 0.0;  // Avoid division by zero

        // Calculate processing speed (work units per millisecond)
        double speed = static_cast<double>(work_done_) / elapsed_time_.count();
       // double speed = work_done_ / elapsed_time_.count(); //Wrong !!! Careful
        // Remaining time = remaining work / speed
       
        return (total_work_ - work_done_) / speed;

      
    }

    /**
     * Display current progress and estimated remaining time
     */
    void print_progress() const {
        // Calculate completion percentage
        double progress = 100.0 * work_done_ / total_work_;
        // Get remaining time in seconds
        double remaining_sec = get_remaining_time() / 1000.0;

        // Format output with 1 decimal place
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "Progress: " << progress << "% ";
        std::cout << "Remaining: " << remaining_sec << "s\r";
        std::cout.flush();  // Use \r to update in place
    }
    size_t get_total_work() const { return total_work_; }
    size_t get_work_done() const { return work_done_; }
    double get_elapsed_ms() const {
        return elapsed_time_.count();
    }
   

private:
    size_t total_work_;      // Total work units to process
    size_t work_done_;       // Work units completed so far
    std::chrono::steady_clock::time_point start_time_;  // Process start time
    std::chrono::milliseconds elapsed_time_;  // Time elapsed since start
};
