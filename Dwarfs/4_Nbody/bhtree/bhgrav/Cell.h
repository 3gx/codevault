#ifndef __CELL_H__
#define __CELL_H__

struct Cell
{
  enum CellType {ROOT};
  typedef std::vector<Cell> Vector;
  private:
  int _addr;
  int _id;
  int _np;
  int iPad;

  public:
  Cell() : _addr(EMPTY), _id(EMPTY), _np(0) {}
  Cell(const CellType type) : _addr(EMPTY+1), _id(EMPTY), _np(0) {}
  Cell(const int addr, const int id, const int np = 0) : _addr(addr), _id(id), _np(np) {setTouched();}
  bool operator==(const Cell &v) const {return _addr == v._addr && _id == v._id;}
  bool isClean() const { return _addr == EMPTY && _id == EMPTY;}
  bool isEmpty() const { return _addr == EMPTY;}
  bool isLeaf () const { return _addr < EMPTY;}
  bool isNode () const { return _addr > EMPTY;}

  int operator++(const int) {_np++; return _np;}
  int operator--(const int) {_np--; return _np;}
  int np() const {return _np;}

  void set_addr(const int addr) {_addr = addr;}
  int addr() const {return _addr;}
  int id() const {return _id & 0x7FFFFFFF;}
  bool    isTouched() const {return _id & 0x80000000;}
  void   setTouched()       { _id |= 0x80000000;}
  void unsetTouched()       { _id = id();}

  int leafIdx() const 
  { 
    assert(_addr < EMPTY); 
    return reverse_int(_addr); 
  }
};

#endif /* __CELL_H__ */
