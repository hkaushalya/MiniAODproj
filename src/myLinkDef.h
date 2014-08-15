#include <vector>
#include <string>
#include <utility>
#ifdef __MAKECINT__
#include <TLorentzVector.h>
#pragma link C++ class std::vector< std::pair<unsigned, unsigned> >+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<vector<int> >+;
#pragma link C++ class std::vector<std::vector<std::vector<int> > >;
#endif
