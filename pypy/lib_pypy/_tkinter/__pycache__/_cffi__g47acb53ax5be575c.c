
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>   /* XXX for ssize_t on some platforms */

/* this block of #ifs should be kept exactly identical between
   c/_cffi_backend.c, cffi/vengine_cpy.py, cffi/vengine_gen.py */
#if defined(_MSC_VER)
# include <malloc.h>   /* for alloca() */
# if _MSC_VER < 1600   /* MSVC < 2010 */
   typedef __int8 int8_t;
   typedef __int16 int16_t;
   typedef __int32 int32_t;
   typedef __int64 int64_t;
   typedef unsigned __int8 uint8_t;
   typedef unsigned __int16 uint16_t;
   typedef unsigned __int32 uint32_t;
   typedef unsigned __int64 uint64_t;
   typedef __int8 int_least8_t;
   typedef __int16 int_least16_t;
   typedef __int32 int_least32_t;
   typedef __int64 int_least64_t;
   typedef unsigned __int8 uint_least8_t;
   typedef unsigned __int16 uint_least16_t;
   typedef unsigned __int32 uint_least32_t;
   typedef unsigned __int64 uint_least64_t;
   typedef __int8 int_fast8_t;
   typedef __int16 int_fast16_t;
   typedef __int32 int_fast32_t;
   typedef __int64 int_fast64_t;
   typedef unsigned __int8 uint_fast8_t;
   typedef unsigned __int16 uint_fast16_t;
   typedef unsigned __int32 uint_fast32_t;
   typedef unsigned __int64 uint_fast64_t;
   typedef __int64 intmax_t;
   typedef unsigned __int64 uintmax_t;
# else
#  include <stdint.h>
# endif
# if _MSC_VER < 1800   /* MSVC < 2013 */
   typedef unsigned char _Bool;
# endif
#else
# include <stdint.h>
# if (defined (__SVR4) && defined (__sun)) || defined(_AIX)
#  include <alloca.h>
# endif
#endif


#include <tcl.h>
#include <tk.h>

char *get_tk_version() { return TK_VERSION; }
char *get_tcl_version() { return TCL_VERSION; }

Tcl_Command _cffi_f_Tcl_CreateCommand(Tcl_Interp * x0, char const * x1, int(* x2)(void *, Tcl_Interp *, int, char const * *), void * x3, void(* x4)(void *))
{
  return Tcl_CreateCommand(x0, x1, x2, x3, x4);
}

Tcl_Interp * _cffi_f_Tcl_CreateInterp(void)
{
  return Tcl_CreateInterp();
}

void _cffi_f_Tcl_DecrRefCount(Tcl_Obj * x0)
{
  Tcl_DecrRefCount(x0);
}

int _cffi_f_Tcl_DeleteCommand(Tcl_Interp * x0, char const * x1)
{
  return Tcl_DeleteCommand(x0, x1);
}

void _cffi_f_Tcl_DeleteInterp(Tcl_Interp * x0)
{
  Tcl_DeleteInterp(x0);
}

int _cffi_f_Tcl_DoOneEvent(int x0)
{
  return Tcl_DoOneEvent(x0);
}

int _cffi_f_Tcl_Eval(Tcl_Interp * x0, char const * x1)
{
  return Tcl_Eval(x0, x1);
}

int _cffi_f_Tcl_EvalFile(Tcl_Interp * x0, char const * x1)
{
  return Tcl_EvalFile(x0, x1);
}

int _cffi_f_Tcl_EvalObjv(Tcl_Interp * x0, int x1, Tcl_Obj * * x2, int x3)
{
  return Tcl_EvalObjv(x0, x1, x2, x3);
}

int _cffi_f_Tcl_ExprBoolean(Tcl_Interp * x0, char const * x1, int * x2)
{
  return Tcl_ExprBoolean(x0, x1, x2);
}

int _cffi_f_Tcl_ExprDouble(Tcl_Interp * x0, char const * x1, double * x2)
{
  return Tcl_ExprDouble(x0, x1, x2);
}

int _cffi_f_Tcl_ExprLong(Tcl_Interp * x0, char const * x1, long * x2)
{
  return Tcl_ExprLong(x0, x1, x2);
}

int _cffi_f_Tcl_ExprString(Tcl_Interp * x0, char const * x1)
{
  return Tcl_ExprString(x0, x1);
}

void _cffi_f_Tcl_FindExecutable(char * x0)
{
  Tcl_FindExecutable(x0);
}

void _cffi_f_Tcl_Free(char * x0)
{
  Tcl_Free(x0);
}

int _cffi_f_Tcl_GetBoolean(Tcl_Interp * x0, char const * x1, int * x2)
{
  return Tcl_GetBoolean(x0, x1, x2);
}

unsigned char * _cffi_f_Tcl_GetByteArrayFromObj(Tcl_Obj * x0, int * x1)
{
  return Tcl_GetByteArrayFromObj(x0, x1);
}

int _cffi_f_Tcl_GetCharLength(Tcl_Obj * x0)
{
  return Tcl_GetCharLength(x0);
}

Tcl_ThreadId _cffi_f_Tcl_GetCurrentThread(void)
{
  return Tcl_GetCurrentThread();
}

int _cffi_f_Tcl_GetDouble(Tcl_Interp * x0, char const * x1, double * x2)
{
  return Tcl_GetDouble(x0, x1, x2);
}

int _cffi_f_Tcl_GetInt(Tcl_Interp * x0, char const * x1, int * x2)
{
  return Tcl_GetInt(x0, x1, x2);
}

Tcl_Obj * _cffi_f_Tcl_GetObjResult(Tcl_Interp * x0)
{
  return Tcl_GetObjResult(x0);
}

Tcl_ObjType const * _cffi_f_Tcl_GetObjType(char const * x0)
{
  return Tcl_GetObjType(x0);
}

char * _cffi_f_Tcl_GetString(Tcl_Obj * x0)
{
  return Tcl_GetString(x0);
}

char * _cffi_f_Tcl_GetStringFromObj(Tcl_Obj * x0, int * x1)
{
  return Tcl_GetStringFromObj(x0, x1);
}

char const * _cffi_f_Tcl_GetStringResult(Tcl_Interp * x0)
{
  return Tcl_GetStringResult(x0);
}

unsigned short * _cffi_f_Tcl_GetUnicode(Tcl_Obj * x0)
{
  return Tcl_GetUnicode(x0);
}

char const * _cffi_f_Tcl_GetVar(Tcl_Interp * x0, char const * x1, int x2)
{
  return Tcl_GetVar(x0, x1, x2);
}

Tcl_Obj * _cffi_f_Tcl_GetVar2Ex(Tcl_Interp * x0, char const * x1, char const * x2, int x3)
{
  return Tcl_GetVar2Ex(x0, x1, x2, x3);
}

void _cffi_f_Tcl_IncrRefCount(Tcl_Obj * x0)
{
  Tcl_IncrRefCount(x0);
}

int _cffi_f_Tcl_Init(Tcl_Interp * x0)
{
  return Tcl_Init(x0);
}

int _cffi_f_Tcl_ListObjGetElements(Tcl_Interp * x0, Tcl_Obj * x1, int * x2, Tcl_Obj * * * x3)
{
  return Tcl_ListObjGetElements(x0, x1, x2, x3);
}

int _cffi_f_Tcl_ListObjIndex(Tcl_Interp * x0, Tcl_Obj * x1, int x2, Tcl_Obj * * x3)
{
  return Tcl_ListObjIndex(x0, x1, x2, x3);
}

int _cffi_f_Tcl_ListObjLength(Tcl_Interp * x0, Tcl_Obj * x1, int * x2)
{
  return Tcl_ListObjLength(x0, x1, x2);
}

Tcl_Obj * _cffi_f_Tcl_NewBooleanObj(int x0)
{
  return Tcl_NewBooleanObj(x0);
}

Tcl_Obj * _cffi_f_Tcl_NewByteArrayObj(unsigned char * x0, int x1)
{
  return Tcl_NewByteArrayObj(x0, x1);
}

Tcl_Obj * _cffi_f_Tcl_NewDoubleObj(double x0)
{
  return Tcl_NewDoubleObj(x0);
}

Tcl_Obj * _cffi_f_Tcl_NewListObj(int x0, Tcl_Obj * * x1)
{
  return Tcl_NewListObj(x0, x1);
}

Tcl_Obj * _cffi_f_Tcl_NewLongObj(long x0)
{
  return Tcl_NewLongObj(x0);
}

Tcl_Obj * _cffi_f_Tcl_NewStringObj(char const * x0, int x1)
{
  return Tcl_NewStringObj(x0, x1);
}

Tcl_Obj * _cffi_f_Tcl_NewUnicodeObj(unsigned short const * x0, int x1)
{
  return Tcl_NewUnicodeObj(x0, x1);
}

void _cffi_f_Tcl_SetObjResult(Tcl_Interp * x0, Tcl_Obj * x1)
{
  Tcl_SetObjResult(x0, x1);
}

char const * _cffi_f_Tcl_SetVar(Tcl_Interp * x0, char const * x1, char const * x2, int x3)
{
  return Tcl_SetVar(x0, x1, x2, x3);
}

char const * _cffi_f_Tcl_SetVar2(Tcl_Interp * x0, char const * x1, char const * x2, char const * x3, int x4)
{
  return Tcl_SetVar2(x0, x1, x2, x3, x4);
}

Tcl_Obj * _cffi_f_Tcl_SetVar2Ex(Tcl_Interp * x0, char const * x1, char const * x2, Tcl_Obj * x3, int x4)
{
  return Tcl_SetVar2Ex(x0, x1, x2, x3, x4);
}

int _cffi_f_Tcl_SplitList(Tcl_Interp * x0, char * x1, int * x2, char const * * * x3)
{
  return Tcl_SplitList(x0, x1, x2, x3);
}

int _cffi_f_Tcl_UnsetVar2(Tcl_Interp * x0, char const * x1, char const * x2, int x3)
{
  return Tcl_UnsetVar2(x0, x1, x2, x3);
}

int _cffi_f_Tk_GetNumMainWindows(void)
{
  return Tk_GetNumMainWindows();
}

int _cffi_f_Tk_Init(Tcl_Interp * x0)
{
  return Tk_Init(x0);
}

char * _cffi_f_get_tcl_version(void)
{
  return get_tcl_version();
}

char * _cffi_f_get_tk_version(void)
{
  return get_tk_version();
}

int _cffi_const_TCL_DONT_WAIT(long long *out_value)
{
  *out_value = (long long)(TCL_DONT_WAIT);
  return (TCL_DONT_WAIT) <= 0;
}

int _cffi_const_TCL_ERROR(long long *out_value)
{
  *out_value = (long long)(TCL_ERROR);
  return (TCL_ERROR) <= 0;
}

int _cffi_const_TCL_EVAL_DIRECT(long long *out_value)
{
  *out_value = (long long)(TCL_EVAL_DIRECT);
  return (TCL_EVAL_DIRECT) <= 0;
}

int _cffi_const_TCL_EVAL_GLOBAL(long long *out_value)
{
  *out_value = (long long)(TCL_EVAL_GLOBAL);
  return (TCL_EVAL_GLOBAL) <= 0;
}

int _cffi_const_TCL_EXCEPTION(long long *out_value)
{
  *out_value = (long long)(TCL_EXCEPTION);
  return (TCL_EXCEPTION) <= 0;
}

int _cffi_const_TCL_GLOBAL_ONLY(long long *out_value)
{
  *out_value = (long long)(TCL_GLOBAL_ONLY);
  return (TCL_GLOBAL_ONLY) <= 0;
}

int _cffi_const_TCL_LEAVE_ERR_MSG(long long *out_value)
{
  *out_value = (long long)(TCL_LEAVE_ERR_MSG);
  return (TCL_LEAVE_ERR_MSG) <= 0;
}

int _cffi_const_TCL_OK(long long *out_value)
{
  *out_value = (long long)(TCL_OK);
  return (TCL_OK) <= 0;
}

int _cffi_const_TCL_READABLE(long long *out_value)
{
  *out_value = (long long)(TCL_READABLE);
  return (TCL_READABLE) <= 0;
}

int _cffi_const_TCL_WRITABLE(long long *out_value)
{
  *out_value = (long long)(TCL_WRITABLE);
  return (TCL_WRITABLE) <= 0;
}

static void _cffi_check_struct_Tcl_Obj(struct Tcl_Obj *p)
{
  /* only to generate compile-time warnings or errors */
  (void)p;
  { char * *tmp = &p->bytes; (void)tmp; }
  (void)((p->length) << 1);
  { Tcl_ObjType * *tmp = &p->typePtr; (void)tmp; }
  /* cannot generate 'union $1' in field 'internalRep': unknown type name */
}
intptr_t _cffi_layout_struct_Tcl_Obj(intptr_t i)
{
  struct _cffi_aligncheck { char x; struct Tcl_Obj y; };
  static intptr_t nums[] = {
    sizeof(struct Tcl_Obj),
    offsetof(struct _cffi_aligncheck, y),
    offsetof(struct Tcl_Obj, bytes),
    sizeof(((struct Tcl_Obj *)0)->bytes),
    offsetof(struct Tcl_Obj, length),
    sizeof(((struct Tcl_Obj *)0)->length),
    offsetof(struct Tcl_Obj, typePtr),
    sizeof(((struct Tcl_Obj *)0)->typePtr),
    offsetof(struct Tcl_Obj, internalRep),
    sizeof(((struct Tcl_Obj *)0)->internalRep),
    -1
  };
  return nums[i];
  /* the next line is not executed, but compiled */
  _cffi_check_struct_Tcl_Obj(0);
}

static void _cffi_check_struct_Tcl_ObjType(struct Tcl_ObjType *p)
{
  /* only to generate compile-time warnings or errors */
  (void)p;
  { char * *tmp = &p->name; (void)tmp; }
}
intptr_t _cffi_layout_struct_Tcl_ObjType(intptr_t i)
{
  struct _cffi_aligncheck { char x; struct Tcl_ObjType y; };
  static intptr_t nums[] = {
    sizeof(struct Tcl_ObjType),
    offsetof(struct _cffi_aligncheck, y),
    offsetof(struct Tcl_ObjType, name),
    sizeof(((struct Tcl_ObjType *)0)->name),
    -1
  };
  return nums[i];
  /* the next line is not executed, but compiled */
  _cffi_check_struct_Tcl_ObjType(0);
}

void init_cffi__g47acb53ax5be575c(void) { }

