#include <iostream>

#define CATCH_CONFIG_RUNNER
#include "catch/catch.hpp"

int Factorial( int number ) {
//    return number <= 1 ? number : Factorial( number - 1 ) * number;  // fail
return number <= 1 ? 1      : Factorial( number - 1 ) * number;  // pass
}

TEST_CASE( "Factorial of 0 is 1 (fail)", "[single-file]" ) {
    REQUIRE( Factorial(0) == 1 );
}

int main( int argc, char* argv[] ) {

#if _MSC_VER
    std::cout << "Hey!! youyouyou" << std::endl;

	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	void *testWhetherMemoryLeakDetectionWorks = malloc(1);
#endif

  // global setup...
  int result = Catch::Session().run( argc, argv );
  // global clean-up...

  return result;
}