#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include "parse_arguments.h"

template<size_t ID>
class WidgetT
{
  private:
    std::string widget_name_;
  public:
    WidgetT(const std::string& widget_name) : widget_name_{widget_name} 
    {
    }
    ~WidgetT()  
    {
    }

    friend std::unique_ptr<WidgetT> makeCutlerySet(const WidgetT& a, const WidgetT& b, const WidgetT &c)
    {
      return std::make_unique<WidgetT>("cutlery set: "+a.widget_name_+" with "+b.widget_name_+" with "+c.widget_name_);
    }
    std::string name() const { return widget_name_; }
};

template<size_t ID>
class Factory
{
  private:
    const std::string factory_name_;
    using Widget = WidgetT<ID>;

    enum {ACQUIRE=1, RELEASE=0};
    void raii(bool acquire)
    {
      static std::string original_name = factory_name_;
      static bool firstInstance = true;
      if (acquire)
      {
        if (!firstInstance)
          throw std::runtime_error(std::string("\n  FATAL: Instance of a Factory  with an unique-id <" + std::to_string(ID) + "> and name \""+original_name+"\", can't be used with a Factory \""+factory_name_+"\".\n"));
        original_name = factory_name_;
        firstInstance = false;
      }
      else
      {
        firstInstance = true;
        original_name = "";
      }
    }
  public:
    Factory(const std::string& factory_name)  : factory_name_(factory_name)
    {
      raii(ACQUIRE);
    }
    ~Factory()  
    {
      raii(RELEASE);
    }
    std::unique_ptr<Widget> makeWidget(const std::string &name) const { return std::make_unique<Widget>(name+"_with_"+factory_name_); }
    const std::string& factoryName() const { return factory_name_; }
};

constexpr size_t hash2(size_t x) { return       ((x >> 16) ^ x);             }
constexpr size_t hash1(size_t x) { return hash2(((x >> 16) ^ x) * 0x45d9f3b);}
constexpr size_t hash (size_t x) { return hash1(((x >> 16) ^ x) * 0x45d9f3b);}
template<size_t ID>
static std::unique_ptr<Factory<hash(ID)>> createFactoryT(const std::string& name)
{
  return std::make_unique<Factory<hash(ID)>>(name);
}
#define createFactory(name) createFactoryT<__COUNTER__>(name)

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

  auto factory1 = createFactory("circles");
  auto widget1A = factory1->makeWidget("spoon");
  auto widget1B = factory1->makeWidget("fork");
  auto widget1C = factory1->makeWidget("knife");
  auto set1     = makeCutlerySet(*widget1A, *widget1B, *widget1C);
  cerr << set1->name() << endl;


#if 1
  std::vector<std::string> types{"flowers","cars","buildings"};
  for (auto type : types)
  {
    auto factory = createFactory(type);
    auto widgetA = factory->makeWidget("spoon");
    auto widgetB = factory->makeWidget("fork");
    auto widgetC = factory->makeWidget("knife");
    auto set     = makeCutlerySet(*widgetA, *widgetB, *widgetC);
    cerr << set->name() << endl;
  }
#endif

  auto factory2 = createFactory("squares");
  auto widget2A = factory2->makeWidget("spoon");
  auto widget2B = factory2->makeWidget("fork");
  auto widget2C = factory2->makeWidget("knife");
  auto set2     = makeCutlerySet(*widget2A, *widget2B, *widget2C);
  cerr << set2->name() << endl;
  
 
#if 0
  auto set3 = makeCutlerySet(*widget1A, *widget1B, *widget2C);
  cerr << set3->name() << endl;
#endif
}
