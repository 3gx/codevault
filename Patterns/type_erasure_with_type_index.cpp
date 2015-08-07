#include <iostream>
#include <string>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>
#include <tuple>
#include <memory>
#include <string>
#include <vector>

using namespace std;

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

class my_class_t
{
  public:
    void draw(ostream& out, size_t position) const
    {
      out << string(position, ' ') << "my_class_t" << endl;
    }
    static string get_type() { return "my_class_t"; }
};

class book_t
{
  public:
  string name;
    void draw(ostream& out, size_t position) const
    {
      out << string(position, ' ') << "book_t{" << name << "}" << endl;
    }
    static string get_type() { return "book_t"; }
};
class name_t
{
  public:
  string name;
    void draw(ostream& out, size_t position) const
    {
      out << string(position, ' ') << "name_t{" << name << "}" << endl;
    }
    static string get_type() { return "name_t"; }
};


using document_t = vector<tuple<object_t,type_index>>;
using document_types = unordered_map<type_index,string>;
void draw(const document_t& x, ostream& out, size_t position)
{
  out << string(position, ' ') << "<document>" << endl;
  for (const auto& e : x )
    get<0>(e).draw(out, position + 2);
  out << string(position, ' ') << "</document>" << endl;
  out << string(position, ' ') << "<types>" << endl;
  for (const auto& e : x )
  {
    auto&& type = get<1>(e);
  }
  out << string(position, ' ') << "</types>" << endl;
}

template<typename T>
tuple<T,type_index> make_doc(T&& doc)
{
  return make_tuple(std::forward<T>(doc),type_index(typeid(T)));
}
int main()
{
  document_t document;
  document.emplace_back(make_doc(my_class_t{}));
  document.emplace_back(make_doc(book_t{"book1"}));
  document.emplace_back(make_doc(name_t{"name1"}));
  document.emplace_back(make_doc(book_t{"book2"}));
  document.emplace_back(make_doc(name_t{"name2"}));
  document.emplace_back(make_doc(my_class_t{}));
  draw(document,cout,0);
}
