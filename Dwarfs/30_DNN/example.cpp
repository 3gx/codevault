#include <iostream>
#include <memory>
#include <cmath>

#include "parse_arguments.h"

using namespace std;

struct Unit
{
  using value_type = double;
  value_type value_{0}, grad_{0};
  Unit(value_type value, value_type grad) : value_(value), grad_(grad) {}
};

static auto make_unit(Unit::value_type value, Unit::value_type gradient = 0)
{
  return std::make_unique<Unit>(value, gradient);
}

struct UnaryGate
{
  Unit* u0_;
  unique_ptr<Unit> utop_;
  virtual Unit* forward(Unit*) = 0;
  virtual Unit* forward(Unit &u0)
  {
    return forward(&u0);
  }
  virtual void backward() = 0;
};
struct BinaryGate 
{
  Unit *u0_,*u1_;
  unique_ptr<Unit> utop_;
  virtual Unit* forward (Unit*, Unit*) = 0;
  virtual Unit* forward (Unit &u0, Unit &u1)  {    return forward(&u0, &u1); }
  virtual Unit* forward (Unit &u0, Unit *u1)  {    return forward(&u0,  u1); }
  virtual Unit* forward (Unit *u0, Unit &u1)  {    return forward( u0, &u1); }
  virtual void backward() = 0;
};

struct Multiply : BinaryGate
{
  virtual Unit* forward(Unit* u0, Unit* u1) override 
  {
    u0_   = u0;
    u1_   = u1;
    utop_ = make_unit(u0_->value_ * u1_->value_, 0);
    return utop_.get();
  }
  using BinaryGate::forward;
  virtual void backward() override
  {
    u0_->grad_ += u1_->value_ * utop_->grad_;
    u1_->grad_ += u0_->value_ * utop_->grad_;
  }
};

struct Add : BinaryGate
{
  virtual Unit* forward(Unit* u0, Unit* u1) override
  {
    u0_   = u0;
    u1_   = u1;
    utop_ = make_unit(u0_->value_ + u1_->value_,0);
    return utop_.get();
  }
  using BinaryGate::forward;
  virtual void backward() override
  {
    u0_->grad_ += 1 * utop_->grad_;
    u1_->grad_ += 1 * utop_->grad_;
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
    return s*(1-s);
  }
  virtual Unit* forward(Unit* u0) override
  {
    u0_   = u0;
    utop_ = make_unit(sig(u0_->value_), 0);
    return utop_.get();
  }
  virtual void backward() override
  {
    u0_->grad_ += grad_sig(u0_->value_) * utop_->grad_;
  }
};



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
    auto s       = sg0 .forward(axpbypc);   // sig(a*x + b*y + c);
    return s;
  };

  auto s = forwardNeuron();
  cout << "circuit output: "  << s->value_ << endl;

  s->grad_ = 1.0;

  sg0  .backward();
  addg1.backward();
  addg0.backward();
  mulg1.backward();
  mulg0.backward();

  auto step_size = 0.01;
  a.value_ += step_size * a.grad_;
  b.value_ += step_size * b.grad_;
  c.value_ += step_size * c.grad_;
  x.value_ += step_size * x.grad_;
  y.value_ += step_size * y.grad_;

  auto s1 = forwardNeuron();
  cout << "circuit output after one backprop: "  << s1->value_ << endl;


  return 0;
}
