#include <vector>
#include <string>
#ifdef __MAKECINT__
#include <TLorentzVector.h>
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<vector<int> >+;
#pragma link C++ class std::vector<std::vector<std::vector<int> > >;
#endif
