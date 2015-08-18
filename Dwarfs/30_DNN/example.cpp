#include <iostream>
#include <memory>
#include <cmath>

#include "parse_arguments.h"

using namespace std;

struct Unit
{
  using value_type = double;
  value_type value_{0}, grad_{0};
};

struct UnaryGate
{
  Unit u0_, utop_;
  virtual Unit forward(Unit) = 0;
  virtual void backward() = 0;
};
struct BinaryGate 
{
  Unit u0_,u1_,utop_;
  virtual Unit forward (Unit, Unit) = 0;
  virtual void backward() = 0;
};

struct Multiply : BinaryGate
{
  virtual Unit forward(Unit u0, Unit u1) override 
  {
    u0_   = std::move(u0);
    u1_   = std::move(u1);
    utop_ = Unit{u0_.value_ * u1_.value_, 0};
    return utop_;
  }
  virtual void backward() override
  {
    u0_.grad_ += u1_.value_ * utop_.grad_;
    u1_.grad_ += u0_.value_ * utop_.grad_;
  }
};

struct Add : BinaryGate
{
  virtual Unit forward(Unit u0, Unit u1) override
  {
    u0_   = std::move(u0);
    u1_   = std::move(u1);
    utop_ = Unit{u0_.value_ + u1_.value_,0};
    return utop_;
  }
  virtual void backward() override
  {
    u0_.grad_ += 1 * utop_.grad_;
    u1_.grad_ += 1 * utop_.grad_;
  }
};

struct Sigmoid : UnaryGate
{
  static auto sig(Unit::value_type x) 
  {
    return 1/(1 + std::exp(-x));
  }
  static auto grad_sig(Unit::value_type x)
  {
    auto s = sig(x);
    return s*(s-1);
  }
  virtual Unit forward(Unit u0) override
  {
    u0_   = std::move(u0);
    utop_ = Unit{sig(u0_.value_), 0};
    return utop_;
  }
  virtual void backward() override
  {
    u0_.grad_ += grad_sig(u0_.value_) * utop_.grad_;
  }
};

auto make_unit(Unit::value_type value, Unit::value_type gradient = 0)
{
  return Unit{value, gradient};
}


int main(int argc, char * argv[])
{
  using namespace std;

  auto a = Unit{ 1,0};
  auto b = Unit{ 2,0};
  auto c = Unit{-3,0};
  auto x = Unit{-1,0};
  auto y = Unit{ 3,0};

  auto mulg0 = Multiply{};
  auto mulg1 = Multiply{};
  auto addg0 = Add{};
  auto addg1 = Add{};
  auto sg0   = Sigmoid{};

  auto forwardNeuron = [&]()
  {
    auto ax      = mulg0.forward(a,x);       // a*x 
    auto by      = mulg1.forward(b,y);       // b*y;
    auto axpby   = addg0.forward(ax,by);     // a*x + b*y
    auto axpbypc = addg1.forward(axpby,c);   // a*x + b*y + c
    auto s       = sg0  .forward(axpbypc);   // sig(a*x + b*y + c);
    return s;
  };

  auto s = forwardNeuron();

  cout << "circuit output: "  << s.value_ << endl;

  return 0;
}
