#pragma once

#include <iostream>
#include <chrono>

class Timer
{
  private:
    using time_type = decltype(std::chrono::system_clock::now());
    std::string timer_;
    bool verbose_destructor_value_{false};
    std::vector<time_type> tbeg_, tend_;
    double dtmin_, dtmax_, dtmean_, dtsigma_;

    static time_type get_time() 
    {
      return std::chrono::system_clock::now();
    }
    static void time_stamp(std::vector<time_type> &time)
    {
      time.emplace_back(std::move(get_time()));
    }
    static double compute_dt(const time_type &t0, const time_type &t1)
    {
      return std::chrono::duration_cast<std::chrono::duration<double>>(t1-t0).count();
    }


  public:
    struct verbose_destructor {};
    Timer(std::string timer) : timer_(std::move(timer)) {}
    Timer(std::string timer, verbose_destructor) : 
      timer_(std::move(timer)), verbose_destructor_value_{true} {}
    void tbeg()
    {
      time_stamp(tbeg_);
    }
    void tend()
    {
      time_stamp(tend_);
    }
    void finalize()
    {
      assert(tbeg_.size() == tend_.size());
      const auto n = tbeg_.size();
      dtmin_   =  1e10;
      dtmax_   = -1e10;
      dtmean_  = 0;
      dtsigma_ = 0;
      for (size_t i = 0; i < n; ++i)
      {
        const auto dt = compute_dt(tbeg_[i], tend_[i]);
        dtmin_    = std::min(dtmin_, dt);
        dtmax_    = std::max(dtmax_, dt);
        dtmean_  += dt;
        dtsigma_ += dt*dt;
      }
      dtmean_  /= n;
      dtsigma_ /= n;
      dtsigma_  = std::sqrt(dtsigma_ - dtmean_*dtmean_);
    }
    template<typename F>
      void print(F f)  const
      {
        std::cout << " * " << timer_ << ": ";
        std::cout << "dt= " << dtmean_ << " +/- " << dtsigma_ << " sec ";
        std::cout << "(dtmin= " << dtmin_ << ", dtmax= " << dtmax_  << ")";
        f(dtmean_, dtsigma_);
        std::cout << "\n";
      }
    void print() const  
    {
      print([](double dt, double ddt) { });
    }
    double dtmean() const
    {
      return dtmean_;
    }
    double dtsigma() const
    {
      return dtsigma_;
    }
    ~Timer()
    {
      if (verbose_destructor_value_)
      {
        finalize();
        print();
      }
    }
};
