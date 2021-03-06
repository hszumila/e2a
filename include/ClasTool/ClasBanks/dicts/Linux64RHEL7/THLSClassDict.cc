//
// File generated by /u/apps/root/5.34.36/root/bin/rootcint at Fri Feb 24 17:38:22 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL7dITHLSClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "THLSClassDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOTShadow {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOTShadow
// END OF SHADOWS

namespace ROOTDict {
   void THLSClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_THLSClass(void *p = 0);
   static void *newArray_THLSClass(Long_t size, void *p);
   static void delete_THLSClass(void *p);
   static void deleteArray_THLSClass(void *p);
   static void destruct_THLSClass(void *p);
   static void streamer_THLSClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::THLSClass*)
   {
      ::THLSClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::THLSClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("THLSClass", ::THLSClass::Class_Version(), "./THLSClass.h", 23,
                  typeid(::THLSClass), ::ROOT::DefineBehavior(ptr, ptr),
                  &::THLSClass::Dictionary, isa_proxy, 0,
                  sizeof(::THLSClass) );
      instance.SetNew(&new_THLSClass);
      instance.SetNewArray(&newArray_THLSClass);
      instance.SetDelete(&delete_THLSClass);
      instance.SetDeleteArray(&deleteArray_THLSClass);
      instance.SetDestructor(&destruct_THLSClass);
      instance.SetStreamerFunc(&streamer_THLSClass);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::THLSClass*)
   {
      return GenerateInitInstanceLocal((::THLSClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::THLSClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr THLSClass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *THLSClass::Class_Name()
{
   return "THLSClass";
}

//______________________________________________________________________________
const char *THLSClass::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::THLSClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int THLSClass::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::THLSClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void THLSClass::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::THLSClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *THLSClass::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::THLSClass*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void THLSClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class THLSClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> TenMhz;
      R__b >> OTR1;
      R__b >> OTR2;
      R__b >> SLM;
      R__b >> LVL1R;
      R__b >> LRC;
      R__b >> LRA;
      R__b >> Fci;
      R__b >> Pmt1;
      R__b >> Pmt2;
      R__b >> Pmt3;
      R__b >> Pmt4;
      R__b >> Res1;
      R__b >> Res2;
      R__b >> HelAcc;
      R__b >> HLSAcc;
      R__b.CheckByteCount(R__s, R__c, THLSClass::IsA());
   } else {
      R__c = R__b.WriteVersion(THLSClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << TenMhz;
      R__b << OTR1;
      R__b << OTR2;
      R__b << SLM;
      R__b << LVL1R;
      R__b << LRC;
      R__b << LRA;
      R__b << Fci;
      R__b << Pmt1;
      R__b << Pmt2;
      R__b << Pmt3;
      R__b << Pmt4;
      R__b << Res1;
      R__b << Res2;
      R__b << HelAcc;
      R__b << HLSAcc;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void THLSClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class THLSClass.
      TClass *R__cl = ::THLSClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "TenMhz", &TenMhz);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "OTR1", &OTR1);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "OTR2", &OTR2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "SLM", &SLM);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "LVL1R", &LVL1R);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "LRC", &LRC);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "LRA", &LRA);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Fci", &Fci);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Pmt1", &Pmt1);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Pmt2", &Pmt2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Pmt3", &Pmt3);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Pmt4", &Pmt4);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Res1", &Res1);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Res2", &Res2);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "HelAcc", &HelAcc);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "HLSAcc", &HLSAcc);
      TObject::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrappers around operator new
   static void *new_THLSClass(void *p) {
      return  p ? new(p) ::THLSClass : new ::THLSClass;
   }
   static void *newArray_THLSClass(Long_t nElements, void *p) {
      return p ? new(p) ::THLSClass[nElements] : new ::THLSClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_THLSClass(void *p) {
      delete ((::THLSClass*)p);
   }
   static void deleteArray_THLSClass(void *p) {
      delete [] ((::THLSClass*)p);
   }
   static void destruct_THLSClass(void *p) {
      typedef ::THLSClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_THLSClass(TBuffer &buf, void *obj) {
      ((::THLSClass*)obj)->::THLSClass::Streamer(buf);
   }
} // end of namespace ROOTDict for class ::THLSClass

/********************************************************
* dicts/Linux64RHEL7/THLSClassDict.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableTHLSClassDict();

extern "C" void G__set_cpp_environmentTHLSClassDict() {
  G__cpp_reset_tagtableTHLSClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTHLSClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* THLSClass */
static int G__THLSClassDict_184_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   THLSClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new THLSClass[n];
     } else {
       p = new((void*) gvp) THLSClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new THLSClass;
     } else {
       p = new((void*) gvp) THLSClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   THLSClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new THLSClass((THLSClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) THLSClass((THLSClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetTenMhz());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetOTR1());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetOTR2());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetSLM());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetLVL1R());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetLRC());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetLRA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetFci());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetPmt1());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetPmt2());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetPmt3());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetPmt4());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetRes1());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetRes2());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetHelAcc());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((THLSClass*) G__getstructoffset())->GetHLSAcc());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((THLSClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) THLSClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) THLSClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) THLSClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      THLSClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((THLSClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_28(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) THLSClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_29(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) THLSClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_30(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) THLSClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__THLSClassDict_184_0_31(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) THLSClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__THLSClassDict_184_0_32(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   THLSClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new THLSClass(*(THLSClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef THLSClass G__TTHLSClass;
static int G__THLSClassDict_184_0_33(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (THLSClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((THLSClass*) (soff+(sizeof(THLSClass)*i)))->~G__TTHLSClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (THLSClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((THLSClass*) (soff))->~G__TTHLSClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__THLSClassDict_184_0_34(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   THLSClass* dest = (THLSClass*) G__getstructoffset();
   *dest = *(THLSClass*) libp->para[0].ref;
   const THLSClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* THLSClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTHLSClassDict {
 public:
  G__Sizep2memfuncTHLSClassDict(): p(&G__Sizep2memfuncTHLSClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTHLSClassDict::*p)();
};

size_t G__get_sizep2memfuncTHLSClassDict()
{
  G__Sizep2memfuncTHLSClassDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceTHLSClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass))) {
     THLSClass *G__Lderived;
     G__Lderived=(THLSClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass),G__get_linked_tagnum(&G__THLSClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTHLSClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__THLSClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__THLSClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THLSClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__THLSClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THLSClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__THLSClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__THLSClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THLSClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__THLSClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__THLSClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* THLSClass */
static void G__setup_memvarTHLSClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass));
   { THLSClass *p; p=(THLSClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->TenMhz)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"TenMhz=",0,"10 Mhz Clock.");
   G__memvar_setup((void*)((long)(&p->OTR1)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"OTR1=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->OTR2)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"OTR2=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->SLM)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"SLM=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->LVL1R)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"LVL1R=",0,"Level 1 trigger rate");
   G__memvar_setup((void*)((long)(&p->LRC)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"LRC=",0,"Left Right Miller Coincidences");
   G__memvar_setup((void*)((long)(&p->LRA)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"LRA=",0,"Left Right Accidentals");
   G__memvar_setup((void*)((long)(&p->Fci)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Fci=",0,"Faraday Cup Current Amplitude");
   G__memvar_setup((void*)((long)(&p->Pmt1)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Pmt1=",0,"Luminosity monitor on beam pipe 1");
   G__memvar_setup((void*)((long)(&p->Pmt2)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Pmt2=",0,"Luminosity monitor on beam pipe 2");
   G__memvar_setup((void*)((long)(&p->Pmt3)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Pmt3=",0,"Luminosity monitor on beam pipe 3");
   G__memvar_setup((void*)((long)(&p->Pmt4)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Pmt4=",0,"Luminosity monitor on beam pipe 4");
   G__memvar_setup((void*)((long)(&p->Res1)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Res1=",0,"Reserved 1");
   G__memvar_setup((void*)((long)(&p->Res2)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"Res2=",0,"Reserved 2");
   G__memvar_setup((void*)((long)(&p->HelAcc)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"HelAcc=",0,"Helicity states accumulating counter.");
   G__memvar_setup((void*)((long)(&p->HLSAcc)-(long)(p)),104,0,0,-1,G__defined_typename("UInt_t"),-1,1,"HLSAcc=",0,"HLS Banks accumulatiing counter.");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__THLSClassDictLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTHLSClassDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTHLSClass(void) {
   /* THLSClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass));
   G__memfunc_setup("THLSClass",817,G__THLSClassDict_184_0_1, 105, G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("THLSClass",817,G__THLSClassDict_184_0_2, 105, G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass), -1, 0, 1, 1, 1, 0, "U 'THLSClass' - 0 - TmpHLS", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetTenMhz",886,G__THLSClassDict_184_0_3, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Sync Signal", (void*) NULL, 0);
   G__memfunc_setup("GetOTR1",582,G__THLSClassDict_184_0_4, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Current Amplitude 2C22", (void*) NULL, 0);
   G__memfunc_setup("GetOTR2",583,G__THLSClassDict_184_0_5, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "X Position Amplitude on BPM 2C22", (void*) NULL, 0);
   G__memfunc_setup("GetSLM",524,G__THLSClassDict_184_0_6, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Y Position Amplitude on BPM 2C22", (void*) NULL, 0);
   G__memfunc_setup("GetLVL1R",657,G__THLSClassDict_184_0_7, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Current Amplitude 2C24", (void*) NULL, 0);
   G__memfunc_setup("GetLRC",513,G__THLSClassDict_184_0_8, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "X Position Amplitude on BPM 2C24", (void*) NULL, 0);
   G__memfunc_setup("GetLRA",511,G__THLSClassDict_184_0_9, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Y Position Amplitude on BPM 2C24", (void*) NULL, 0);
   G__memfunc_setup("GetFci",562,G__THLSClassDict_184_0_10, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Faraday Cup Current Amplitude", (void*) NULL, 0);
   G__memfunc_setup("GetPmt1",642,G__THLSClassDict_184_0_11, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Luminosity monitor on beam pipe 1", (void*) NULL, 0);
   G__memfunc_setup("GetPmt2",643,G__THLSClassDict_184_0_12, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Luminosity monitor on beam pipe 2", (void*) NULL, 0);
   G__memfunc_setup("GetPmt3",644,G__THLSClassDict_184_0_13, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Luminosity monitor on beam pipe 3", (void*) NULL, 0);
   G__memfunc_setup("GetPmt4",645,G__THLSClassDict_184_0_14, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Luminosity monitor on beam pipe 4", (void*) NULL, 0);
   G__memfunc_setup("GetRes1",635,G__THLSClassDict_184_0_15, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Reserved 1", (void*) NULL, 0);
   G__memfunc_setup("GetRes2",636,G__THLSClassDict_184_0_16, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", "Reserved 2", (void*) NULL, 0);
   G__memfunc_setup("GetHelAcc",832,G__THLSClassDict_184_0_17, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetHLSAcc",782,G__THLSClassDict_184_0_18, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__THLSClassDict_184_0_19, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__THLSClassDict_184_0_20, 85, G__get_linked_tagnum(&G__THLSClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&THLSClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__THLSClassDict_184_0_21, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&THLSClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__THLSClassDict_184_0_22, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&THLSClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__THLSClassDict_184_0_23, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&THLSClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__THLSClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__THLSClassDict_184_0_27, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__THLSClassDict_184_0_28, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&THLSClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__THLSClassDict_184_0_29, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&THLSClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__THLSClassDict_184_0_30, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&THLSClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__THLSClassDict_184_0_31, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&THLSClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("THLSClass", 817, G__THLSClassDict_184_0_32, (int) ('i'), G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass), -1, 0, 1, 1, 1, 0, "u 'THLSClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~THLSClass", 943, G__THLSClassDict_184_0_33, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__THLSClassDict_184_0_34, (int) ('u'), G__get_linked_tagnum(&G__THLSClassDictLN_THLSClass), -1, 1, 1, 1, 1, 0, "u 'THLSClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTHLSClassDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalTHLSClassDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {
}

static void G__cpp_setup_func18() {
}

static void G__cpp_setup_func19() {
}

static void G__cpp_setup_func20() {
}

static void G__cpp_setup_func21() {
}

static void G__cpp_setup_func22() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTHLSClassDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
  G__cpp_setup_func18();
  G__cpp_setup_func19();
  G__cpp_setup_func20();
  G__cpp_setup_func21();
  G__cpp_setup_func22();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__THLSClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__THLSClassDictLN_THLSClass = { "THLSClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTHLSClassDict() {
  G__THLSClassDictLN_TClass.tagnum = -1 ;
  G__THLSClassDictLN_TBuffer.tagnum = -1 ;
  G__THLSClassDictLN_TMemberInspector.tagnum = -1 ;
  G__THLSClassDictLN_TObject.tagnum = -1 ;
  G__THLSClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__THLSClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__THLSClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__THLSClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__THLSClassDictLN_THLSClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTHLSClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__THLSClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__THLSClassDictLN_THLSClass),sizeof(THLSClass),-1,62720,"Scaler Bank HLS for Helicity related signals",G__setup_memvarTHLSClass,G__setup_memfuncTHLSClass);
}
extern "C" void G__cpp_setupTHLSClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTHLSClassDict()");
  G__set_cpp_environmentTHLSClassDict();
  G__cpp_setup_tagtableTHLSClassDict();

  G__cpp_setup_inheritanceTHLSClassDict();

  G__cpp_setup_typetableTHLSClassDict();

  G__cpp_setup_memvarTHLSClassDict();

  G__cpp_setup_memfuncTHLSClassDict();
  G__cpp_setup_globalTHLSClassDict();
  G__cpp_setup_funcTHLSClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTHLSClassDict();
  return;
}
class G__cpp_setup_initTHLSClassDict {
  public:
    G__cpp_setup_initTHLSClassDict() { G__add_setup_func("THLSClassDict",(G__incsetup)(&G__cpp_setupTHLSClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTHLSClassDict() { G__remove_setup_func("THLSClassDict"); }
};
G__cpp_setup_initTHLSClassDict G__cpp_setup_initializerTHLSClassDict;

