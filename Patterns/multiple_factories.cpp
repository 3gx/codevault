#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include "parse_arguments.h"

template<size_t ID_>
class WidgetT
{
  private:
    std::string widget_name_;
  public:
    enum {ID = ID_};
    WidgetT(const std::string& widget_name) : widget_name_{widget_name} 
    {
    }
    ~WidgetT()  
    {
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

template<size_t ID>
static std::unique_ptr<Factory<ID>> createFactoryT(const std::string& name)
{
  return std::make_unique<Factory<ID>>(name);
}

template<size_t x> class HASH
{
  enum {y = ((x >> 16) ^ x) * 0x45d9f3b, z = ((y >> 16) ^ y) * 0x45d9f3b};
  public:
  enum {value = (z >> 16) ^ z};
};
#define createFactory(name) createFactoryT<HASH<__COUNTER__>::value>(name)

template<typename Widget>
std::unique_ptr<Widget> makeCutlerySet(const Widget& a, const Widget& b, const Widget& c)
{
  return std::make_unique<Widget>("cutlery set: "+a.name()+" with "+b.name() +" with "+c.name());
}

template<typename Widget1, typename Widget2, typename Widget3>
auto makeCutlerySet(const Widget1& a, const Widget2& b, const Widget3& c)
{
  using namespace std;
  constexpr auto newID = Widget1::ID + Widget2::ID + Widget3::ID;
  using NewWidget = WidgetT<HASH<newID>::value>;
  const auto set = tie(a,b,c);
  if (drand48() < 0.3)
    return std::make_unique<NewWidget>("cutlery mix: "+get<0>(set).name()+" with "+get<1>(set).name() +" with "+get<2>(set).name());
  else if (drand48() < 0.6)
    return std::make_unique<NewWidget>("cutlery mix: "+get<1>(set).name()+" with "+get<2>(set).name() +" with "+get<0>(set).name());
  else
    return std::make_unique<NewWidget>("cutlery mix: "+get<2>(set).name()+" with "+get<0>(set).name() +" with "+get<1>(set).name());
}

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

  auto factory2 = createFactory("squares");
  auto widget2A = factory2->makeWidget("spoon");
  auto widget2B = factory2->makeWidget("fork");
  auto widget2C = factory2->makeWidget("knife");
  auto set2     = makeCutlerySet(*widget2A, *widget2B, *widget2C);
  cerr << set2->name() << endl;

  cerr << " ---- \n";

  std::vector<std::string> types{"flowers","cars","buildings"};
  for (auto type : types)
  {
    auto factory = createFactory(type);
    auto widgetA = factory->makeWidget("spoon");
    auto widgetB = factory->makeWidget("fork");
    auto widgetC = factory->makeWidget("knife");
    auto mix = makeCutlerySet(*widgetA, *widget1B, *widget2C);
    cerr << mix->name() << endl;
  }

  cerr << " ---- \n";
  
 
#if 0
  auto set3 = makeCutlerySet(*widget1A, *widget1B, *widget2C);
  cerr << set3->name() << endl;
#endif
}
