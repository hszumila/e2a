//
// File generated by /u/apps/root/5.34.13/root/bin/rootcint at Fri Jan 20 19:15:32 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL6dITICHBClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TICHBClassDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
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

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void TICHBClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TICHBClass(void *p = 0);
   static void *newArray_TICHBClass(Long_t size, void *p);
   static void delete_TICHBClass(void *p);
   static void deleteArray_TICHBClass(void *p);
   static void destruct_TICHBClass(void *p);
   static void streamer_TICHBClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TICHBClass*)
   {
      ::TICHBClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TICHBClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TICHBClass", ::TICHBClass::Class_Version(), "./TICHBClass.h", 13,
                  typeid(::TICHBClass), DefineBehavior(ptr, ptr),
                  &::TICHBClass::Dictionary, isa_proxy, 0,
                  sizeof(::TICHBClass) );
      instance.SetNew(&new_TICHBClass);
      instance.SetNewArray(&newArray_TICHBClass);
      instance.SetDelete(&delete_TICHBClass);
      instance.SetDeleteArray(&deleteArray_TICHBClass);
      instance.SetDestructor(&destruct_TICHBClass);
      instance.SetStreamerFunc(&streamer_TICHBClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TICHBClass*)
   {
      return GenerateInitInstanceLocal((::TICHBClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TICHBClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TICHBClass::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TICHBClass::Class_Name()
{
   return "TICHBClass";
}

//______________________________________________________________________________
const char *TICHBClass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TICHBClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TICHBClass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TICHBClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TICHBClass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TICHBClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TICHBClass::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TICHBClass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TICHBClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class TICHBClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Eclust;
      R__b >> Eclmax;
      R__b >> Tclust;
      R__b >> T_next;
      R__b >> xclust;
      R__b >> yclust;
      R__b >> zclust;
      R__b >> xclmax;
      R__b >> yclmax;
      R__b >> M2_x;
      R__b >> M2_y;
      R__b >> M3_x;
      R__b >> M3_y;
      R__b >> ncryst;
      R__b.CheckByteCount(R__s, R__c, TICHBClass::IsA());
   } else {
      R__c = R__b.WriteVersion(TICHBClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Eclust;
      R__b << Eclmax;
      R__b << Tclust;
      R__b << T_next;
      R__b << xclust;
      R__b << yclust;
      R__b << zclust;
      R__b << xclmax;
      R__b << yclmax;
      R__b << M2_x;
      R__b << M2_y;
      R__b << M3_x;
      R__b << M3_y;
      R__b << ncryst;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TICHBClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TICHBClass.
      TClass *R__cl = ::TICHBClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Eclust", &Eclust);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Eclmax", &Eclmax);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Tclust", &Tclust);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "T_next", &T_next);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "xclust", &xclust);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "yclust", &yclust);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "zclust", &zclust);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "xclmax", &xclmax);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "yclmax", &yclmax);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "M2_x", &M2_x);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "M2_y", &M2_y);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "M3_x", &M3_x);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "M3_y", &M3_y);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "ncryst", &ncryst);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TICHBClass(void *p) {
      return  p ? new(p) ::TICHBClass : new ::TICHBClass;
   }
   static void *newArray_TICHBClass(Long_t nElements, void *p) {
      return p ? new(p) ::TICHBClass[nElements] : new ::TICHBClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_TICHBClass(void *p) {
      delete ((::TICHBClass*)p);
   }
   static void deleteArray_TICHBClass(void *p) {
      delete [] ((::TICHBClass*)p);
   }
   static void destruct_TICHBClass(void *p) {
      typedef ::TICHBClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TICHBClass(TBuffer &buf, void *obj) {
      ((::TICHBClass*)obj)->::TICHBClass::Streamer(buf);
   }
} // end of namespace ROOT for class ::TICHBClass

/********************************************************
* dicts/Linux64RHEL6/TICHBClassDict.cc
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

extern "C" void G__cpp_reset_tagtableTICHBClassDict();

extern "C" void G__set_cpp_environmentTICHBClassDict() {
  G__cpp_reset_tagtableTICHBClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTICHBClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TICHBClass */
static int G__TICHBClassDict_183_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TICHBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TICHBClass[n];
     } else {
       p = new((void*) gvp) TICHBClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TICHBClass;
     } else {
       p = new((void*) gvp) TICHBClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TICHBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TICHBClass((TICHBClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TICHBClass((TICHBClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TICHBClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TICHBClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TICHBClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TICHBClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TICHBClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TICHBClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TICHBClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TICHBClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TICHBClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TICHBClassDict_183_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TICHBClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TICHBClassDict_183_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TICHBClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TICHBClass(*(TICHBClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TICHBClass G__TTICHBClass;
static int G__TICHBClassDict_183_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (TICHBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TICHBClass*) (soff+(sizeof(TICHBClass)*i)))->~G__TTICHBClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TICHBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TICHBClass*) (soff))->~G__TTICHBClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TICHBClassDict_183_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TICHBClass* dest = (TICHBClass*) G__getstructoffset();
   *dest = *(TICHBClass*) libp->para[0].ref;
   const TICHBClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TICHBClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTICHBClassDict {
 public:
  G__Sizep2memfuncTICHBClassDict(): p(&G__Sizep2memfuncTICHBClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTICHBClassDict::*p)();
};

size_t G__get_sizep2memfuncTICHBClassDict()
{
  G__Sizep2memfuncTICHBClassDict a;
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
extern "C" void G__cpp_setup_inheritanceTICHBClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass))) {
     TICHBClass *G__Lderived;
     G__Lderived=(TICHBClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass),G__get_linked_tagnum(&G__TICHBClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTICHBClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TICHBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TICHBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TICHBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TICHBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TICHBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TICHBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TICHBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TICHBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TICHBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TICHBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TICHBClass */
static void G__setup_memvarTICHBClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass));
   { TICHBClass *p; p=(TICHBClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->Eclust)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Eclust=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->Eclmax)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Eclmax=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->Tclust)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Tclust=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->T_next)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"T_next=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->xclust)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"xclust=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->yclust)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"yclust=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->zclust)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"zclust=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->xclmax)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"xclmax=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->yclmax)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"yclmax=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->M2_x)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"M2_x=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->M2_y)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"M2_y=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->M3_x)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"M3_x=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->M3_y)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"M3_y=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->ncryst)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"ncryst=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TICHBClassDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTICHBClassDict() {
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
static void G__setup_memfuncTICHBClass(void) {
   /* TICHBClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass));
   G__memfunc_setup("TICHBClass",864,G__TICHBClassDict_183_0_1, 105, G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TICHBClass",864,G__TICHBClassDict_183_0_2, 105, G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass), -1, 0, 1, 1, 1, 0, "U 'TICHBClass' - 0 - Tmp", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__TICHBClassDict_183_0_3, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TICHBClassDict_183_0_4, 85, G__get_linked_tagnum(&G__TICHBClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TICHBClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TICHBClassDict_183_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TICHBClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TICHBClassDict_183_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TICHBClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TICHBClassDict_183_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TICHBClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TICHBClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TICHBClassDict_183_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TICHBClassDict_183_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TICHBClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TICHBClassDict_183_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TICHBClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TICHBClassDict_183_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TICHBClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TICHBClassDict_183_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TICHBClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TICHBClass", 864, G__TICHBClassDict_183_0_16, (int) ('i'), G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass), -1, 0, 1, 1, 1, 0, "u 'TICHBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TICHBClass", 990, G__TICHBClassDict_183_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TICHBClassDict_183_0_18, (int) ('u'), G__get_linked_tagnum(&G__TICHBClassDictLN_TICHBClass), -1, 1, 1, 1, 1, 0, "u 'TICHBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTICHBClassDict() {
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
extern "C" void G__cpp_setup_globalTICHBClassDict() {
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

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTICHBClassDict() {
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
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__TICHBClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TICHBClassDictLN_TICHBClass = { "TICHBClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTICHBClassDict() {
  G__TICHBClassDictLN_TClass.tagnum = -1 ;
  G__TICHBClassDictLN_TBuffer.tagnum = -1 ;
  G__TICHBClassDictLN_TMemberInspector.tagnum = -1 ;
  G__TICHBClassDictLN_TObject.tagnum = -1 ;
  G__TICHBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TICHBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TICHBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TICHBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TICHBClassDictLN_TICHBClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTICHBClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TICHBClassDictLN_TICHBClass),sizeof(TICHBClass),-1,62720,(char*)NULL,G__setup_memvarTICHBClass,G__setup_memfuncTICHBClass);
}
extern "C" void G__cpp_setupTICHBClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTICHBClassDict()");
  G__set_cpp_environmentTICHBClassDict();
  G__cpp_setup_tagtableTICHBClassDict();

  G__cpp_setup_inheritanceTICHBClassDict();

  G__cpp_setup_typetableTICHBClassDict();

  G__cpp_setup_memvarTICHBClassDict();

  G__cpp_setup_memfuncTICHBClassDict();
  G__cpp_setup_globalTICHBClassDict();
  G__cpp_setup_funcTICHBClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTICHBClassDict();
  return;
}
class G__cpp_setup_initTICHBClassDict {
  public:
    G__cpp_setup_initTICHBClassDict() { G__add_setup_func("TICHBClassDict",(G__incsetup)(&G__cpp_setupTICHBClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTICHBClassDict() { G__remove_setup_func("TICHBClassDict"); }
};
G__cpp_setup_initTICHBClassDict G__cpp_setup_initializerTICHBClassDict;

