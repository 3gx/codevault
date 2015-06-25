#pragma once

template<typename T>
class make_range_iterator
{
  T _begin, _end, _skip;
  public:
    class iterator : std::iterator<std::bidirectional_iterator_tag, T>
    {
      private:
        T __value, __skip;
      public:
        iterator(T value, T skip) : __value(value), __skip(skip) {}
        operator T() const  { return __value; }
        operator T& ()      { return __value; }
        T operator*() const { return __value; }
        bool operator!=(const iterator& it) const { return __value != it.__value; }
        const iterator& operator++()  { __value += __skip; return *this; }
    };

    make_range_iterator(T begin, T end) : _begin(begin), _end(end), _skip(1) {}
    make_range_iterator(T end) : _begin(0), _end(end), _skip(1) {}
    iterator begin() { return iterator{_begin,_skip}; } 
    iterator end  () { return iterator{_end,  _skip}; }
}; 

template<size_t BEG, size_t END> //, size_t SKIP=1>
class make_range_iteratorT
{
  enum {SKIP = 1};
  using T = size_t;
  public:
    class iterator : std::iterator<std::bidirectional_iterator_tag, T>
    {
      private:
        T __value;
      public:
        iterator(T value) : __value(value) {}
        operator T() const  { return __value; }
        operator T& ()      { return __value; }
        T operator*() const { return __value; }
        bool operator!=(const iterator& it) const { return __value != it.__value; }
        const iterator& operator++()  { __value += SKIP; return *this; }
    };

    iterator begin() { return iterator{BEG}; } 
    iterator end  () { return iterator{END}; }
}; 
