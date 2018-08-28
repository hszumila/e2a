// Autogenerated Class (Header File)
// Author : G.Gavalian
// Date   : Thu Jul 16 01:32:08 EDT 2009
//

#ifndef __TUserParameter__
#define __TUserParameter__
#include <iostream>
#include <TROOT.h>
#include <TVector3.h>
#include <TObject.h>
#include <TString.h>
#include <TNamed.h>

class TUserParameter : public TNamed {

private:

  TString  fParName;
  TString  fParValue;

public:

TUserParameter ();
 TUserParameter (const char *pname, const char *pval);
~TUserParameter ();

 void      Set(const char *pname, const char *pval)
 {fParName = pname; fParValue = pval; SetName(pname);}
 TString   GetValue(){return fParValue;};

ClassDef(TUserParameter,1)


};
#endif
