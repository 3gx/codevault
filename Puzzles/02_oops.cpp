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
#if 1
      Base:func();  
#else  /* this is not Base::func() but Base: is a label */
      Base:
           func();
#endif
      cout << "Derived's func() now runnign()\n";
    }
};

int main()
{
  Base *bp = new Derived;
  bp->func();
}
