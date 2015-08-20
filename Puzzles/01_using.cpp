#include <iostream>

int main()
{
#ifndef PART2
  using namespace std;
#else
  using std::cout;
  using std::endl;
#endif

  cout << "Hello, puzzling world!" << endl;

  int cout;
  cout = 10;
  std::cout << "cout = " << cout << endl;
  std::cout << "cout << 5 = " << (cout << 5)  << endl;

  return 0;
}
