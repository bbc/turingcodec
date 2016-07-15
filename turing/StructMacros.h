/*
Copyright (C) 2016 British Broadcasting Corporation, Parabola Research
and Queen Mary University of London.

This file is part of the Turing codec.

The Turing codec is free software; you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as
published by the Free Software Foundation.

The Turing codec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Commercial support and intellectual property rights for
the Turing codec are also available under a proprietary license.
For more information, contact us at info @ turingcodec.org.
 */

// Convenience macros for defining structures containing a small number of integer parameters.

#ifndef INCLUDED_StructMacros_h
#define INCLUDED_StructMacros_h

#pragma once

#define DEFINE_STRUCT_ARITY_0(name) \
    struct name \
    { \
        static const int arity = 0;\
        static const char *argName(int n) \
        { \
            assert(!"should never be called"); \
            return 0; \
        } \
        int &argValue (int n) \
        { \
            assert(!"should never be called"); \
            static int dummy; \
            return dummy; \
        } \
    };

#define DEFINE_STRUCT_ARITY_1(name, arg0) \
    struct name \
    { \
        int arg0; \
        name(int arg0=0) \
        : \
          arg0(arg0) \
          { \
          } \
          static const int arity = 1;\
          static const char *argName(int n) \
          { \
              const char *names[] = { #arg0 }; \
              return names[n]; \
          } \
          int &argValue (int n) \
          { \
              return this->arg0; \
          } \
    };

#define DEFINE_STRUCT_ARITY_2(name, arg0, arg1) \
    struct name \
    { \
        int arg0; \
        int arg1; \
        name(int arg0=0, int arg1=0) \
        : \
          arg0(arg0), \
          arg1(arg1) \
          { \
          } \
          static const int arity = 2;\
          static const char *argName(int n) \
          { \
              const char *names[] = { #arg0, #arg1 }; \
              return names[n]; \
          } \
          int &argValue (int n) \
          { \
              if (n==1) return this->arg1; \
              return this->arg0; \
          } \
    };

#define DEFINE_STRUCT_ARITY_3(name, arg0, arg1, arg2) \
    struct name \
    { \
        int arg0; \
        int arg1; \
        int arg2; \
        name(int arg0=0, int arg1=0, int arg2=0) \
        : \
          arg0(arg0), \
          arg1(arg1), \
          arg2(arg2) \
          { \
          } \
          static const int arity = 3;\
          static const char *argName(int n) \
          { \
              const char *names[] = { #arg0, #arg1, #arg2 }; \
              return names[n]; \
          } \
          int &argValue (int n) \
          { \
              if (n==2) return this->arg2; \
              if (n==1) return this->arg1; \
              return this->arg0; \
          } \
          bool operator==(name const &rhs) const \
          { \
              return this->arg0 == rhs.arg0 && this->arg1 == rhs.arg1 && this->arg2 == rhs.arg2; \
          } \
    };

#define DEFINE_STRUCT_ARITY_4(name, arg0, arg1, arg2, arg3) \
    struct name \
    { \
        int arg0; \
        int arg1; \
        int arg2; \
        int arg3; \
        name(int arg0=0, int arg1=0, int arg2=0, int arg3=0) \
        : \
          arg0(arg0), \
          arg1(arg1), \
          arg2(arg2), \
          arg3(arg3) \
          { \
          } \
          static const int arity = 4;\
          static const char *argName(int n) \
          { \
              const char *names[] = { #arg0, #arg1, #arg2, #arg3 }; \
              return names[n]; \
          } \
          int &argValue (int n) \
          { \
              if (n==3) return this->arg3; \
              if (n==2) return this->arg2; \
              if (n==1) return this->arg1; \
              return this->arg0; \
          } \
          bool operator==(const name &other) const \
          { \
              return \
                      this->arg3 == other.arg3 && \
                      this->arg2 == other.arg2 && \
                      this->arg1 == other.arg1 && \
                      this->arg0 == other.arg0; \
          } \
    };

#define DEFINE_STRUCT_ARITY_5(name, arg0, arg1, arg2, arg3, arg4) \
    struct name \
    { \
        int arg0; \
        int arg1; \
        int arg2; \
        int arg3; \
        int arg4; \
        name(int arg0=0, int arg1=0, int arg2=0, int arg3=0, int arg4=0) \
        : \
          arg0(arg0), \
          arg1(arg1), \
          arg2(arg2), \
          arg3(arg3), \
          arg4(arg4) \
          { \
          } \
          static const int arity = 7;\
          static const char *argName(int n) \
          { \
              const char *names[] = { #arg0, #arg1, #arg2, #arg3, #arg4 }; \
              return names[n]; \
          } \
          int &argValue (int n) \
          { \
              if (n==4) return this->arg4; \
              if (n==3) return this->arg3; \
              if (n==2) return this->arg2; \
              if (n==1) return this->arg1; \
              return this->arg0; \
          } \
          bool operator==(const name &other) const \
          { \
              return \
                      this->arg4 == other.arg4 && \
                      this->arg3 == other.arg3 && \
                      this->arg2 == other.arg2 && \
                      this->arg1 == other.arg1 && \
                      this->arg0 == other.arg0; \
          } \
    };

#define DEFINE_STRUCT_ARITY_7(name, arg0, arg1, arg2, arg3, arg4, arg5, arg6) \
    struct name \
    { \
        int arg0; \
        int arg1; \
        int arg2; \
        int arg3; \
        int arg4; \
        int arg5; \
        int arg6; \
        name(int arg0=0, int arg1=0, int arg2=0, int arg3=0, int arg4=0, int arg5=0, int arg6=0) \
        : \
          arg0(arg0), \
          arg1(arg1), \
          arg2(arg2), \
          arg3(arg3), \
          arg4(arg4), \
          arg5(arg5), \
          arg6(arg6) \
          { \
          } \
          static const int arity = 7;\
          static const char *argName(int n) \
          { \
              const char *names[] = { #arg0, #arg1, #arg2, #arg3, #arg4, #arg5, #arg6 }; \
              return names[n]; \
          } \
          int &argValue (int n) \
          { \
              if (n==6) return this->arg6; \
              if (n==5) return this->arg5; \
              if (n==4) return this->arg4; \
              if (n==3) return this->arg3; \
              if (n==2) return this->arg2; \
              if (n==1) return this->arg1; \
              return this->arg0; \
          } \
    };

#endif
