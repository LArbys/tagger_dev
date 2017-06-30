#include <iostream>

int main( int nargs, char** argv ) {
  
  std::cout << "COMPARE TRACK CHARGE" << std::endl;

  // What we're doing.

  // 1) Load information using DataCoordinator
  // 2) Load truth end points for cosmics
  // 3) use the end points to reco the tracks (either linearcharge or astar tracker)
  // 4) use the reco track points to calculate the flash hypothesis chi-2
  // 5) make plots comparing the relative values on the planes -- is the distribution good enough?

  // 6) make reco end points. repeat exercise.  But we can now compare flash-matching for good and bad tracks
  
  return 0;
}
