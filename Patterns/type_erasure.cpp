#include <iostream>
#include <string>
#include <typeinfo>
#include <typeindex>
#include <memory>
#include <string>
#include <vector>

using namespace std;

#if 0
class object_t
{
  public:
    virtual ~object_t() {}
    virtual void draw(ostream&, size_t) const = 0;
};

using document_t = vector<shared_ptr<object_t>>;

void draw(const document_t& x, ostream& out, size_t position)
{
  out << string(position, ' ') << "<document>" << endl;
  for (const auto& e : x )
    e->draw(out, position + 2);
  out << string(position, ' ') << "</document>" << endl;
}

class my_class_t : public object_t
{
  public:
    void draw(ostream &out, size_t position) const
    {
      out << string(position, ' ') << "my_class_t" << endl;
    }
};

int main(int argc, char * argv[])
{
  document_t document;
  document.emplace_back(new my_class_t{});
  draw(document,cout, 1);
  return 0;
}
#else
class object_t
{
  public:
    template<typename T> 
      object_t(T x) : self_(make_shared<model<T>>(move(x))) {}
    void draw(ostream& out, size_t position) const
    {
      self_->draw_(out,position);
    }

  private:
    struct concept_t
    {
      virtual ~concept_t() = default;
      virtual void draw_(ostream&, size_t) const = 0;
    };

    template<typename T>
    struct model : concept_t
    {
      model(T x) : data_(move(x)) {}
      void draw_(ostream& out, size_t position) const
      {
        data_.draw(out,position);
      }
      T data_;
    };

    shared_ptr<const concept_t> self_;
};

using document_t = vector<object_t>;

void draw(const document_t& x, ostream& out, size_t position)
{
  out << string(position, ' ') << "<document>" << endl;
  for (const auto& e : x )
    e.draw(out, position + 2);
  out << string(position, ' ') << "</document>" << endl;
}

class my_class_t
{
  public:
    void draw(ostream& out, size_t position) const
    {
      out << string(position, ' ') << "my_class_t" << endl;
    }
};

int main()
{
  document_t document;
  document.emplace_back(my_class_t{});
  draw(document,cout,0);
}
#endif
