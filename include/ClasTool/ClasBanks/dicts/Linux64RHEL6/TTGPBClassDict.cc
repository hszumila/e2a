//
// File generated by /u/apps/root/5.34.13/root/bin/rootcint at Fri Jan 20 19:15:47 2017

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dictsdILinux64RHEL6dITTGPBClassDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TTGPBClassDict.h"

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
   void TTGPBClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TTGPBClass(void *p = 0);
   static void *newArray_TTGPBClass(Long_t size, void *p);
   static void delete_TTGPBClass(void *p);
   static void deleteArray_TTGPBClass(void *p);
   static void destruct_TTGPBClass(void *p);
   static void streamer_TTGPBClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TTGPBClass*)
   {
      ::TTGPBClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TTGPBClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TTGPBClass", ::TTGPBClass::Class_Version(), "./TTGPBClass.h", 16,
                  typeid(::TTGPBClass), DefineBehavior(ptr, ptr),
                  &::TTGPBClass::Dictionary, isa_proxy, 0,
                  sizeof(::TTGPBClass) );
      instance.SetNew(&new_TTGPBClass);
      instance.SetNewArray(&newArray_TTGPBClass);
      instance.SetDelete(&delete_TTGPBClass);
      instance.SetDeleteArray(&deleteArray_TTGPBClass);
      instance.SetDestructor(&destruct_TTGPBClass);
      instance.SetStreamerFunc(&streamer_TTGPBClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TTGPBClass*)
   {
      return GenerateInitInstanceLocal((::TTGPBClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TTGPBClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TTGPBClass::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TTGPBClass::Class_Name()
{
   return "TTGPBClass";
}

//______________________________________________________________________________
const char *TTGPBClass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TTGPBClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TTGPBClass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TTGPBClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TTGPBClass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TTGPBClass*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TTGPBClass::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TTGPBClass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TTGPBClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class TTGPBClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> pointer;
      R__b >> Time;
      R__b >> Energy;
      R__b >> dt;
      R__b.CheckByteCount(R__s, R__c, TTGPBClass::IsA());
   } else {
      R__c = R__b.WriteVersion(TTGPBClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << pointer;
      R__b << Time;
      R__b << Energy;
      R__b << dt;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TTGPBClass::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TTGPBClass.
      TClass *R__cl = ::TTGPBClass::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "pointer", &pointer);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Time", &Time);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "Energy", &Energy);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "dt", &dt);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TTGPBClass(void *p) {
      return  p ? new(p) ::TTGPBClass : new ::TTGPBClass;
   }
   static void *newArray_TTGPBClass(Long_t nElements, void *p) {
      return p ? new(p) ::TTGPBClass[nElements] : new ::TTGPBClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_TTGPBClass(void *p) {
      delete ((::TTGPBClass*)p);
   }
   static void deleteArray_TTGPBClass(void *p) {
      delete [] ((::TTGPBClass*)p);
   }
   static void destruct_TTGPBClass(void *p) {
      typedef ::TTGPBClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TTGPBClass(TBuffer &buf, void *obj) {
      ((::TTGPBClass*)obj)->::TTGPBClass::Streamer(buf);
   }
} // end of namespace ROOT for class ::TTGPBClass

/********************************************************
* dicts/Linux64RHEL6/TTGPBClassDict.cc
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

extern "C" void G__cpp_reset_tagtableTTGPBClassDict();

extern "C" void G__set_cpp_environmentTTGPBClassDict() {
  G__cpp_reset_tagtableTTGPBClassDict();
}
#include <new>
extern "C" int G__cpp_dllrevTTGPBClassDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TTGPBClass */
static int G__TTGPBClassDict_183_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TTGPBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TTGPBClass[n];
     } else {
       p = new((void*) gvp) TTGPBClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TTGPBClass;
     } else {
       p = new((void*) gvp) TTGPBClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TTGPBClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TTGPBClass((TTGPBClass*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) TTGPBClass((TTGPBClass*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((TTGPBClass*) G__getstructoffset())->GetPointer());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TTGPBClass*) G__getstructoffset())->GetTime());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TTGPBClass*) G__getstructoffset())->GetEnergy());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((TTGPBClass*) G__getstructoffset())->GetDt());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TTGPBClass*) G__getstructoffset())->Print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TTGPBClass::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TTGPBClass::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TTGPBClass::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TTGPBClass::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TTGPBClass*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TTGPBClass::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TTGPBClass::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TTGPBClass::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TTGPBClassDict_183_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TTGPBClass::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TTGPBClassDict_183_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TTGPBClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TTGPBClass(*(TTGPBClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TTGPBClass G__TTTGPBClass;
static int G__TTGPBClassDict_183_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (TTGPBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TTGPBClass*) (soff+(sizeof(TTGPBClass)*i)))->~G__TTTGPBClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TTGPBClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TTGPBClass*) (soff))->~G__TTTGPBClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TTGPBClassDict_183_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TTGPBClass* dest = (TTGPBClass*) G__getstructoffset();
   *dest = *(TTGPBClass*) libp->para[0].ref;
   const TTGPBClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TTGPBClass */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTTGPBClassDict {
 public:
  G__Sizep2memfuncTTGPBClassDict(): p(&G__Sizep2memfuncTTGPBClassDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTTGPBClassDict::*p)();
};

size_t G__get_sizep2memfuncTTGPBClassDict()
{
  G__Sizep2memfuncTTGPBClassDict a;
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
extern "C" void G__cpp_setup_inheritanceTTGPBClassDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass))) {
     TTGPBClass *G__Lderived;
     G__Lderived=(TTGPBClass*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass),G__get_linked_tagnum(&G__TTGPBClassDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTTGPBClassDict() {

   /* Setting up typedef entry */
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("Float_t",102,-1,0,-1);
   G__setnewtype(-1,"Float 4 bytes (float)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TTGPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TTGPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TTGPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TTGPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TTGPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TTGPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TTGPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TTGPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TTGPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TTGPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TTGPBClass */
static void G__setup_memvarTTGPBClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass));
   { TTGPBClass *p; p=(TTGPBClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->pointer)-(long)(p)),105,0,0,-1,G__defined_typename("Int_t"),-1,1,"pointer=",0,"pointer to TAGR");
   G__memvar_setup((void*)((long)(&p->Time)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Time=",0,"starttime_TAG at interaction point(ns)");
   G__memvar_setup((void*)((long)(&p->Energy)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"Energy=",0,"photon energy (GeV)");
   G__memvar_setup((void*)((long)(&p->dt)-(long)(p)),102,0,0,-1,G__defined_typename("Float_t"),-1,1,"dt=",0,"starttime_ST - starttime_TAG (ns) if no starttime_ST, dt = -starttime_TAG");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TTGPBClassDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTTGPBClassDict() {
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
static void G__setup_memfuncTTGPBClass(void) {
   /* TTGPBClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass));
   G__memfunc_setup("TTGPBClass",887,G__TTGPBClassDict_183_0_1, 105, G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TTGPBClass",887,G__TTGPBClassDict_183_0_2, 105, G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass), -1, 0, 1, 1, 1, 0, "U 'TTGPBClass' - 0 - TmpSTPB", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetPointer",1025,G__TTGPBClassDict_183_0_3, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetTime",687,G__TTGPBClassDict_183_0_4, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetEnergy",906,G__TTGPBClassDict_183_0_5, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetDt",472,G__TTGPBClassDict_183_0_6, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Print",525,G__TTGPBClassDict_183_0_7, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TTGPBClassDict_183_0_8, 85, G__get_linked_tagnum(&G__TTGPBClassDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TTGPBClass::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TTGPBClassDict_183_0_9, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TTGPBClass::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TTGPBClassDict_183_0_10, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TTGPBClass::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TTGPBClassDict_183_0_11, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TTGPBClass::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TTGPBClassDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TTGPBClassDict_183_0_15, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TTGPBClassDict_183_0_16, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TTGPBClass::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TTGPBClassDict_183_0_17, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TTGPBClass::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TTGPBClassDict_183_0_18, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TTGPBClass::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TTGPBClassDict_183_0_19, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TTGPBClass::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TTGPBClass", 887, G__TTGPBClassDict_183_0_20, (int) ('i'), G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass), -1, 0, 1, 1, 1, 0, "u 'TTGPBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TTGPBClass", 1013, G__TTGPBClassDict_183_0_21, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TTGPBClassDict_183_0_22, (int) ('u'), G__get_linked_tagnum(&G__TTGPBClassDictLN_TTGPBClass), -1, 1, 1, 1, 1, 0, "u 'TTGPBClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTTGPBClassDict() {
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
extern "C" void G__cpp_setup_globalTTGPBClassDict() {
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

extern "C" void G__cpp_setup_funcTTGPBClassDict() {
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
G__linked_taginfo G__TTGPBClassDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TTGPBClassDictLN_TTGPBClass = { "TTGPBClass" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTTGPBClassDict() {
  G__TTGPBClassDictLN_TClass.tagnum = -1 ;
  G__TTGPBClassDictLN_TBuffer.tagnum = -1 ;
  G__TTGPBClassDictLN_TMemberInspector.tagnum = -1 ;
  G__TTGPBClassDictLN_TObject.tagnum = -1 ;
  G__TTGPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TTGPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TTGPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TTGPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TTGPBClassDictLN_TTGPBClass.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTTGPBClassDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TTGPBClassDictLN_TTGPBClass),sizeof(TTGPBClass),-1,62720,"Class for accessing the TGPB bank: Tagger",G__setup_memvarTTGPBClass,G__setup_memfuncTTGPBClass);
}
extern "C" void G__cpp_setupTTGPBClassDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTTGPBClassDict()");
  G__set_cpp_environmentTTGPBClassDict();
  G__cpp_setup_tagtableTTGPBClassDict();

  G__cpp_setup_inheritanceTTGPBClassDict();

  G__cpp_setup_typetableTTGPBClassDict();

  G__cpp_setup_memvarTTGPBClassDict();

  G__cpp_setup_memfuncTTGPBClassDict();
  G__cpp_setup_globalTTGPBClassDict();
  G__cpp_setup_funcTTGPBClassDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTTGPBClassDict();
  return;
}
class G__cpp_setup_initTTGPBClassDict {
  public:
    G__cpp_setup_initTTGPBClassDict() { G__add_setup_func("TTGPBClassDict",(G__incsetup)(&G__cpp_setupTTGPBClassDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTTGPBClassDict() { G__remove_setup_func("TTGPBClassDict"); }
};
G__cpp_setup_initTTGPBClassDict G__cpp_setup_initializerTTGPBClassDict;

