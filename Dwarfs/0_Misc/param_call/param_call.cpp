#include <tuple>
#include <utility>
#include <iostream>
using namespace std;

#ifndef CXX14
#if __cplusplus >= 201202L
#define CXX14
#endif
#endif /* CXX14 */

#ifndef CXX14
template<size_t... I>
struct index_sequence {};

template<size_t N, size_t... S>
struct make_index_sequence_impl : make_index_sequence_impl<N-1,N-1,S...> {};

template<size_t... S>
struct make_index_sequence_impl<0,S...>
{
    using type = index_sequence<S...>;
};

template<size_t I>
using make_index_sequence = typename make_index_sequence_impl<I>::type;

template<size_t I, typename Tuple>
using tuple_element_t = typename tuple_element<I,Tuple>::type;
#endif



template<typename Functor, typename Tuple, size_t... I>
void packCallImpl(Functor&& func, Tuple&& t, index_sequence<I...>)
{
  forward_as_tuple(func(I,forward<tuple_element_t<I,Tuple>>(get<I>(t)))...);
}

template<typename Functor, typename... Args>
void packCall(Functor&& func, Args&&... args)
{
  packCallImpl(func, forward_as_tuple(std::forward<Args>(args)...),make_index_sequence<sizeof...(Args)>());
}

struct Foo
{
  template<typename T>
    int operator()(int i, T&& x)
    {
      cout << "i= " << i << ": value= " << x << endl;
      return 0;
    }
};

int main()
{
  cout << "-----------\n functor\n";
  packCall(Foo{}, 1.0,2.0,"std::string", -3);
#ifdef CXX14
  cout << "----------\n generic lambda \n";
  packCall(
      [](auto i, auto x) {cout << "i= " << i << ": value= " << x << endl; return 0;},
  1.0,2.0,"std::string", -3);
#endif
}
