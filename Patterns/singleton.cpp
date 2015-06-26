#include <iostream>
using namespace std;

class Singleton
{
  string value;
  private:
    Singleton()
    {
      value = "activated";
      cerr << "Singleton activated\n";
    }
    ~Singleton()
    {
      cerr << "Singleton destroyed\n";
    }
  public:
    static Singleton& getInstance()
    {
      static Singleton s;
      return s;
    }
    string getValue() const 
    {
      return getInstance().value;
    };
    void setValue(string str) const
    {
      getInstance().value = str;
    }
};

void print()
{
  decltype(auto) s = Singleton::getInstance();

  cerr << s.getValue() << endl;
}

int main()
{
  decltype(auto) s = Singleton::getInstance();

  print();

  s.setValue("test");
  
  print();

}
