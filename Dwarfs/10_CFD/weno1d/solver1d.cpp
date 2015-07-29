#include <iostream>
#include "parse_arguments.h"
#include "timer.h"

int main(int argc, char*argv[])
{
  using namespace std;
  using namespace parse_arguments;
  auto ncell = 128;

  auto params = pack(argc, argv,
      param("cell count", ncell, "n", "ncell")
      );

  cerr << params.parse_all();

  return 0;
}
