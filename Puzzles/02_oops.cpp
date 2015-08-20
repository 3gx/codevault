#include <iostream>

using namespace std;

class Base
{
  public:
    virtual void func()
    {
      cout << "Base's func() now runnign \n";
    };
};

class Derived : public Base
{
  public:
    void func()
    {
      Base:func();
      cout << "Derived's func() now runnign()\n";
    }
};

int main()
{
  Base *bp = new Derived;
  bp->func();
}
