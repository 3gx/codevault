#include <iostream>
#include <string>
#include <cassert>
#include "parse_arguments.h"

template<size_t ID>
class WidgetT
{
  private:
    std::string widget_name_;
  public:
    WidgetT(const std::string& widget_name) : widget_name_{widget_name} 
    {
//      std::cout << "    Widget \"" << widget_name_ << "\" ctor. \n";
    }
    ~WidgetT()  
    {
//      std::cout << "    Widget \"" << widget_name_ << "\" dtor. \n";
    }

    friend std::unique_ptr<WidgetT> makeSet(const WidgetT& a, const WidgetT& b, const WidgetT &c)
    {
      return std::make_unique<WidgetT>("set: "+a.widget_name_+" with "+b.widget_name_+" with "+c.widget_name_);
    }
    std::string name() const { return widget_name_; }
};

template<size_t ID>
class Factory
{
  private:
    std::string factory_name_;
    using Widget = WidgetT<ID>;
  public:
    Factory(const std::string& factory_name_) : factory_name_{factory_name_}
    {
      static bool firstInstance = true;
      if (!firstInstance)
        throw std::runtime_error(std::string("\n\tFATAL: Instance of a Factory \"")+ factory_name_ + "\"" + ", with an internal unique-id <" + std::to_string(ID) + ">, is already in use!\n");
      firstInstance = false;
 //     std::cout << "Factory \"" << factory_name_ << "\" ctor. \n";
    }
    ~Factory()  
    {
  //    std::cout << "Factory \"" << factory_name_ << "\" dtor. \n";
    }
    std::unique_ptr<Widget> makeWidget(const std::string &name) const { return std::make_unique<Widget>(factory_name_+"_"+name); }
    const std::string& factoryName() const { return factory_name_; }
};

#define buildFactory(name) Factory<__COUNTER__>(name)

int main(int argc, char *argv[])
{
  using namespace std;
#if 0
  using namespace parse_arguments;
  auto nfactory = 1;
  auto params = pack(argc, argv,
      param("number of factories", nfactory, "n","")
      );
  cerr << params.parse_all();
#endif

#if 1
  auto factory1 = Factory<1>("silverware");
  auto factory2 = Factory<2>("bronzeware");
#else
  auto factory1 = buildFactory("silverware");
  auto factory2 = buildFactory("bronzeware");
#endif

  auto widget1A = factory1.makeWidget("spoon");
  auto widget1B = factory1.makeWidget("fork");
  auto widget1C = factory1.makeWidget("knife");
  auto set1     = makeSet(*widget1A, *widget1B, *widget1C);
  cerr << set1->name() << endl;

 
  auto widget2A = factory2.makeWidget("spoon");
  auto widget2B = factory2.makeWidget("fork");
  auto widget2C = factory2.makeWidget("knife");
  auto set2     = makeSet(*widget2A, *widget2B, *widget2C);
  cerr << set2->name() << endl;
  
 
#if 0
  auto set3 = makeSet(*widget1A, *widget1B, *widget2C);
  cerr << set3->name() << endl;
#endif
}
