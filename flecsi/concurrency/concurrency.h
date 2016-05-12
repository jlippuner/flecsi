/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_concurrency_h
#define flecsi_concurrency_h

#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <functional>
#include <thread>

/*!
 * \file concurrency.h
 * \authors nickm
 * \date Initial file creation: Oct 09, 2015
 */

namespace flecsi
{

  class virtual_semaphore{
  public:
    
    virtual_semaphore(int count=0)
    : count_(count),
    max_count_(0){
      done_ = false;
    }
    
    virtual_semaphore(int count, int max_count)
    : count_(count),
    max_count_(max_count){
      done_ = false;
    }
    
    ~virtual_semaphore(){}
    
    bool acquire(){
      if(done_){
        return false;
      }

      lock_.lock();
      
      while(count_ <= 0){
        cond_.wait(lock_);
      }
      
      --count_;
      lock_.unlock();
      
      return true;
    }
    
    bool try_acquire(){
      if(done_){
        return false;
      }

      lock_.lock();

      if(count_ > 0){
        --count_;
        lock_.unlock();
        return true;
      }

      lock_.unlock();
      return false;
    }
    
    void release(){
      lock_.lock();

      if(max_count_ == 0 || count_ < max_count_){
        ++count_;
      }

      cond_.notify_one();
      lock_.unlock();
    }

    void interrupt(){
      done_ = true;
      cond_.notify_all();
    }
    
    virtual_semaphore& operator=(const virtual_semaphore&) = delete;
    
    virtual_semaphore(const virtual_semaphore&) = delete;

  private:
    std::unique_lock<std::mutex> lock_;
    std::condition_variable cond_;

    int count_;
    int max_count_;
    atomic_bool done_;
  };

  class thread_pool{
  public:
    using function_type = std::function<void(void)>;

    thread_pool(){
      done_ = false;
    }

    void run_(){
      for(;;){
        sem_.acquire();
        mutex_.lock();

        if(done_){
          mutex_.unlock();
          return;
        }

        auto f = queue_.front();
        queue_.pop_front();
        mutex_.unlock();
        f();
      }
    }

    void queue(function_type f){
      mutex_.lock();
      queue_.emplace_back(f);
      mutex_.unlock();
      sem_.release();
    }

    void start(){
      auto t = new thread(&thread_pool::run_, this);
      threads_.push_back(t);
    }

    void join(){
      for(auto t : threads_){
        t->join();
      }
    }

    void stop(){
      done_ = true;
    }

  private:
    std::mutex mutex_;
    std::deque<function_type> queue_;
    std::vector<std::thread*> threads_;
    virtual_semaphore sem_;
    std::atomic_bool done_;    
  };

} // namespace flecsi

#endif // flecsi_concurrency_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
