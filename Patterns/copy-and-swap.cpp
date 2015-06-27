#include <iostream>

struct MyDouble
{
  double x_;
  MyDouble(double x) : x_(x) { std::cout << "   MyDouble:: ctor\n"; }
  ~MyDouble() { std::cout << "   MyDouble:: dtor\n"; }

  MyDouble(const MyDouble &rhs) : x_(rhs.x_) { std::cout << "   MyDouble:: copy ctor \n"; }
  MyDouble& operator=(const MyDouble& rhs)
  {
    std::cout << "   MyDouble:: copy assign \n";
    x_ = rhs.x_;
    return *this;
  }
#if 1  /* set to 0 to disable move semantics */
  MyDouble(MyDouble &&rhs) : x_(std::move(rhs.x_)) {std::cout << "   MyDouble:: move ctor \n";} 
  MyDouble& operator=(MyDouble&& rhs)
  {
    std::cout << "   MyDouble:: move assign \n";
    x_ = std::move(rhs.x_);
    rhs.x_ = 0;
    return *this;
  }
#endif
  friend std::ostream& operator << (std::ostream &os, const MyDouble &x)
  {
    return os << x.x_;
  }
};

template<typename T>
struct CopyAndSwap
{
  T x_;

  CopyAndSwap(const double x) : x_(x) { std::cout << "CopyAndSwap:: ctor \n"; }
  CopyAndSwap(const CopyAndSwap& rhs) : x_(rhs.x_) { std::cout << "CopyAndSwap:: copy ctor \n"; }
  CopyAndSwap(CopyAndSwap&& rhs) : x_(std::move(rhs.x_)) { std::cout << "CopyAndSwap:: move ctor \n"; }
  ~CopyAndSwap() { std::cout << "CopyAndSwap:: dtor \n";}

  CopyAndSwap& operator=(CopyAndSwap rhs)
  {
    std::cout << "CopyAndSwap:: assign \n";
    swap(*this,rhs);
    return *this;
  }
  friend void swap(CopyAndSwap& first, CopyAndSwap &second)
  {
    std::cout << "CopyAndSwap:: swap \n";
    std::swap(first.x_, second.x_);
  }
  friend std::ostream& operator << (std::ostream &os, const CopyAndSwap &x)
  {
    return os << x.x_;
  }
};

template<typename T>
struct CopyMove
{
  T x_;

  CopyMove(const double x) : x_(x) { std::cout << "CopyMove:: ctor \n"; }
  CopyMove(const CopyMove&rhs): x_(rhs.x_) { std::cout << "CopyMove:: copy ctor \n"; }
  explicit CopyMove(CopyMove&& rhs)  : x_(std::move(rhs.x_)) { std::cout << "CopyMove:: move ctor \n"; }
  ~CopyMove() { std::cout << "CopyMove:: dtor \n";}

  CopyMove& operator=(CopyMove&& rhs)
  {
    std::cout << "CopyMove:: move assignment\n";
    x_ = std::move(rhs.x_);
    return *this;
  }
  CopyMove& operator=(const CopyMove& rhs)
  {
    std::cout << "CopyMove:: copy assignment\n";
    x_ = rhs.x_;
    return *this;
  }
  friend std::ostream& operator << (std::ostream &os, const CopyMove &x)
  {
    return os << x.x_;
  }
};

template<typename T>
void test()
{
  std::cout << " >construct< \n";
  T x{1};
  T y{2};
  std::cout << "---------------output: x= " << x << " y= " << y << std::endl;
  std::cout << " >copy< \n";
  y = x;
  std::cout << "---------------output: x= " << x << " y= " << y << std::endl;
  std::cout << " >move< \n";
  y = std::move(x);
  std::cout << "---------------output: x= " << x << " y= " << y << std::endl;
  std::cout << " >destroy< \n";
}

int main()
{
#if 0
  using base_type = double;
#else
  using base_type = MyDouble;
#endif
  
  std::cout << "\n";
  std::cout << " ~~~~~~~~~~~~ \n";
  std::cout << "\n";

  test<CopyAndSwap<base_type>>();

  std::cout << "\n";
  std::cout << " ~~~~~~~~~~~~ \n";
  std::cout << "\n";

  test<CopyMove<base_type>>();


  std::cout << "\n";
  std::cout << " ~~~~~~~~~~~~ \n";
  std::cout << "\n";
}
