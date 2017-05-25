/*
Copyright (C) 2016 British Broadcasting Corporation, Parabola Research
and Queen Mary University of London.

This file is part of the Turing codec.

The Turing codec is free software; you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as
published by the Free Software Foundation.

The Turing codec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

Commercial support and intellectual property rights for
the Turing codec are also available under a proprietary license.
For more information, contact us at info @ turingcodec.org.
 */

 // Context models' state.
 // Review: consider renaming to "StateContextModels"

#ifndef INCLUDED_Cabac_h
#define INCLUDED_Cabac_h

#pragma once

#include "ContextModel.h"
#include "Global.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstdint>


// Metafunction with function "static int get(int initType)"
template <class Tag> struct ContextOffset;

// Metafunction with static "int value"
template <class Tag> struct ContextTotalSize;

// Metafunction with static "int value"
template <class Tag> struct ContextMaxSize;


#define TERNARY_MAX(A, B) ((A)>(B)?(A):(B))

#define CONTEXT_OFFSETS(name, o0, o1, o2, o3) \
template <> struct ContextOffset<name> \
{ \
    static int get(int initType) \
    { \
            static const int value[4] = { o0, o1, o2, o3 }; \
            return value[initType]; \
    } \
}; \
template <> struct ContextTotalSize<name> \
{ \
    static const int value = o3; \
}; \
template <> struct ContextMaxSize<name> \
{ \
    static const int value = TERNARY_MAX(o1 - o0, TERNARY_MAX(o2 - o1, o3 - o2)); \
};



// Metafunction to tell compiler which context to use for a given syntax value V.
// Default behaviour: pass type through as most elements have their own contexts, e.g. split_cu_flag.
// Some, e.g. sao_merge_*_flag, share a context and their specialisms are below.
template <class V>
struct Context
{
    typedef V Type;
};


// sao()

// shared context
struct sao_merge_X_flag { };
template <> struct Context<sao_merge_left_flag> { typedef sao_merge_X_flag Type; };
template <> struct Context<sao_merge_up_flag> { typedef sao_merge_X_flag Type; };
CONTEXT_OFFSETS(sao_merge_X_flag, 0, 1, 2, 3);


// shared context
struct sao_type_idx_X { };
template <> struct Context<sao_type_idx_luma> { typedef sao_type_idx_X Type; };
template <> struct Context<sao_type_idx_chroma> { typedef sao_type_idx_X Type; };
CONTEXT_OFFSETS(sao_type_idx_X, 0, 1, 2, 3);


// coding_quadtree()
CONTEXT_OFFSETS(split_cu_flag, 0, 3, 6, 9);


// coding_unit()
CONTEXT_OFFSETS(cu_transquant_bypass_flag, 0, 1, 2, 3);
CONTEXT_OFFSETS(cu_skip_flag, 0, 0, 3, 6);
CONTEXT_OFFSETS(cu_qp_delta_abs, 0, 2, 4, 6);
CONTEXT_OFFSETS(cu_chroma_qp_offset_flag, 0, 1, 2, 3);
CONTEXT_OFFSETS(cu_chroma_qp_offset_idx, 0, 1, 2, 3);
CONTEXT_OFFSETS(pred_mode_flag, 0, 0, 1, 2);
CONTEXT_OFFSETS(part_mode, 0, 1, 5, 9);


// prediction_unit()
CONTEXT_OFFSETS(prev_intra_luma_pred_flag, 0, 1, 2, 3);
CONTEXT_OFFSETS(intra_chroma_pred_mode, 0, 1, 2, 3);
CONTEXT_OFFSETS(merge_flag, 0, 0, 1, 2);
CONTEXT_OFFSETS(merge_idx, 0, 0, 1, 2);
CONTEXT_OFFSETS(inter_pred_idc, 0, 0, 5, 10);
struct ref_idx_lX { };
template <> struct Context<ref_idx_l0> { typedef ref_idx_lX Type; };
template <> struct Context<ref_idx_l1> { typedef ref_idx_lX Type; };
CONTEXT_OFFSETS(ref_idx_lX, 0, 0, 2, 4);
struct abs_mvd_greaterX_flag { };
CONTEXT_OFFSETS(abs_mvd_greater0_flag, 0, 0, 1, 2);
CONTEXT_OFFSETS(abs_mvd_greater1_flag, 2, 2, 3, 4);
template <> struct ContextTotalSize<abs_mvd_greaterX_flag> { static const int value = 4; };
struct mvp_lX_flag { };
template <> struct Context<mvp_l0_flag> { typedef mvp_lX_flag Type; };
template <> struct Context<mvp_l1_flag> { typedef mvp_lX_flag Type; };
CONTEXT_OFFSETS(mvp_lX_flag, 0, 0, 1, 2);

// transform_tree()
CONTEXT_OFFSETS(rqt_root_cbf, 0, 0, 1, 2);
CONTEXT_OFFSETS(split_transform_flag, 0, 3, 6, 9);
CONTEXT_OFFSETS(cbf_luma, 0, 2, 4, 6);

struct cbf_cX { };
template <> struct Context<cbf_cb> { typedef cbf_cX Type; };
template <> struct Context<cbf_cr> { typedef cbf_cX Type; };
CONTEXT_OFFSETS(cbf_cX, 0, 4, 8, 12);

// residual_coding()

// review how transform_skip_flag is handled: a lot of complexity here just so that initialisation table is same order as in the text
struct transform_skip_flag_X { };
struct transform_skip_flag_Y { };
struct transform_skip_flag_C { };
template <> struct Context<transform_skip_flag_0> { typedef transform_skip_flag_Y Type; };
template <> struct Context<transform_skip_flag_1> { typedef transform_skip_flag_C Type; };
template <> struct Context<transform_skip_flag_2> { typedef transform_skip_flag_C Type; };
CONTEXT_OFFSETS(transform_skip_flag_Y, 0, 1, 2, 3);
CONTEXT_OFFSETS(transform_skip_flag_C, 3, 4, 5, 6);
template <> struct ContextTotalSize<transform_skip_flag_X> { static const int value = 6; };
CONTEXT_OFFSETS(explicit_rdpcm_flag, 0, 0, 2, 4);
CONTEXT_OFFSETS(explicit_rdpcm_dir_flag, 0, 0, 2, 4);
CONTEXT_OFFSETS(last_sig_coeff_x_prefix, 0, 18, 36, 54);
CONTEXT_OFFSETS(last_sig_coeff_y_prefix, 0, 18, 36, 54);
CONTEXT_OFFSETS(coded_sub_block_flag, 0, 4, 8, 12);
CONTEXT_OFFSETS(sig_coeff_flag, 0, 44, 88, 132);
CONTEXT_OFFSETS(coeff_abs_level_greater1_flag, 0, 24, 48, 72);
CONTEXT_OFFSETS(coeff_abs_level_greater2_flag, 0, 6, 12, 18);


// cross_comp_pred()
CONTEXT_OFFSETS(log2_res_scale_abs_plus1, 0, 8, 16, 24);
CONTEXT_OFFSETS(res_scale_sign_flag, 0, 2, 4, 6);


template <class Tag>
struct ContextInit
{
    static int initValue[ContextTotalSize<Tag>::value];
};

template<> struct ContextInit<abs_mvd_greater0_flag> : ContextInit<abs_mvd_greaterX_flag> { };
template<> struct ContextInit<abs_mvd_greater1_flag> : ContextInit<abs_mvd_greaterX_flag> { };

template<> struct ContextInit<transform_skip_flag_Y> : ContextInit<transform_skip_flag_X> { };
template<> struct ContextInit<transform_skip_flag_C> : ContextInit<transform_skip_flag_X> { };


// review: duplication of TypeName<> code here

static bool trimb2(std::string &s)
{
    size_t pos = s.find_first_of(" :<");
    if (pos == std::string::npos) return false;
    if (s[pos] == '<')
    {
        s.erase(pos);
        return false;
    }
    pos += (s[pos] == ':') ? 2 : 1;
    s.erase(0, pos);
    return true;
}

static std::string trimb(std::string name)
{
    while (trimb2(name));
    while (name.find_first_of("0123456789") == 0) name = name.substr(1);
    return name;
}

template <class T>
struct TypeName2
{
    static const std::string value;
};

template <class T>
const std::string TypeName2<T>::value(trimb(typeid(T).name()));


template <class Tag, int size>
struct ContextsBase2
{
    ContextsBase2() {}

    ContextsBase2(int qp, int initType)
    {
        for (int ctxIdx = ContextOffset<Tag>::get(initType); ctxIdx < ContextOffset<Tag>::get(initType + 1); ++ctxIdx)
        {
            int const ctxInc = ctxIdx - ContextOffset<Tag>::get(initType);
            this->cm[ctxInc] = ContextModel(qp, ContextInit<Tag>::initValue[ctxIdx]);
        }
    }

    // review: generic accept() (visitor pattern) might make this struct cleaner
    void print(std::ostream &o, int initType)
    {
        o << "\t" << TypeName2<Tag>::value << ": ";
        for (int ctxIdx = ContextOffset<Tag>::get(initType); ctxIdx < ContextOffset<Tag>::get(initType + 1); ++ctxIdx)
        {
            const int ctxInc = ctxIdx - ContextOffset<Tag>::get(initType);
            o << (int)this->cm[ctxInc].state << ", ";
        }
        o << "\n";
    }

    static int offset(int initType)
    {
        return ContextOffset<Tag>::get(initType);
    }

    std::array<ContextModel, size> cm;
    static const int n = size;

    void checkSameAs(const ContextsBase2 &other, int initType) const
    {
        for (int ctxIdx = ContextOffset<Tag>::get(initType); ctxIdx < ContextOffset<Tag>::get(initType + 1); ++ctxIdx)
        {
            int const ctxInc = ctxIdx - ContextOffset<Tag>::get(initType);
            ASSERT(this->cm[ctxInc] == other.cm[ctxInc]);
        }
    }
};


// Specialisation for zero size
template <class Tag>
struct ContextsBase2<Tag, 0>
{
    ContextsBase2() {}
    ContextsBase2(int qp, int initType) {}

    bool operator==(const ContextsBase2 &) const { return true; }
};


template <class Tag>
struct ContextsBase :
    ContextsBase2<Tag, ContextMaxSize<Tag>::value>
{
    typedef Tag Type;
    typedef ContextOffset<Tag> Offset;

    ContextsBase()
    {
    }
    ContextsBase(int qp, int initType) :
        ContextsBase2<Tag, ContextMaxSize<Tag>::value>(qp, initType)
    {
    }
};


template <template <class> class BaseTemplate>
struct HevcContextsResidual :
    BaseTemplate<rqt_root_cbf>,
    BaseTemplate<split_transform_flag>,
    BaseTemplate<cbf_luma>,
    BaseTemplate<cbf_cX>,
    BaseTemplate<transform_skip_flag_Y>,
    BaseTemplate<transform_skip_flag_C>,
    BaseTemplate<explicit_rdpcm_flag>,
    BaseTemplate<explicit_rdpcm_dir_flag>,
    BaseTemplate<last_sig_coeff_x_prefix>,
    BaseTemplate<last_sig_coeff_y_prefix>,
    BaseTemplate<coded_sub_block_flag>,
    BaseTemplate<sig_coeff_flag>,
    BaseTemplate<coeff_abs_level_greater1_flag>,
    BaseTemplate<coeff_abs_level_greater2_flag>,
    BaseTemplate<log2_res_scale_abs_plus1>,
    BaseTemplate<res_scale_sign_flag>
{
    template <class Visitor>
    void accept(Visitor &v)
    {
        v.visitContextsBase(static_cast<BaseTemplate<rqt_root_cbf>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<split_transform_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cbf_luma>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cbf_cX>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<transform_skip_flag_Y>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<transform_skip_flag_C>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<explicit_rdpcm_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<explicit_rdpcm_dir_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<last_sig_coeff_x_prefix>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<last_sig_coeff_y_prefix>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<coded_sub_block_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<sig_coeff_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<coeff_abs_level_greater1_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<coeff_abs_level_greater2_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<log2_res_scale_abs_plus1>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<res_scale_sign_flag>*>(this));
    }
};


template <template <class> class BaseTemplate>
struct HevcContexts :
    BaseTemplate<sao_merge_X_flag>,
    BaseTemplate<sao_type_idx_X>,
    BaseTemplate<split_cu_flag>,
    BaseTemplate<cu_transquant_bypass_flag>,
    BaseTemplate<cu_skip_flag>,
    BaseTemplate<cu_qp_delta_abs>,
    BaseTemplate<cu_chroma_qp_offset_flag>,
    BaseTemplate<cu_chroma_qp_offset_idx>,
    BaseTemplate<pred_mode_flag>,
    BaseTemplate<part_mode>,
    BaseTemplate<prev_intra_luma_pred_flag>,
    BaseTemplate<intra_chroma_pred_mode>,
    BaseTemplate<merge_flag>,
    BaseTemplate<merge_idx>,
    BaseTemplate<inter_pred_idc>,
    BaseTemplate<ref_idx_lX>,
    BaseTemplate<abs_mvd_greater0_flag>,
    BaseTemplate<abs_mvd_greater1_flag>,
    BaseTemplate<mvp_lX_flag>,
    HevcContextsResidual<BaseTemplate>
{
    template <class Visitor>
    void accept(Visitor &v)
    {
        v.visitContextsBase(static_cast<BaseTemplate<sao_merge_X_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<sao_type_idx_X>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<split_cu_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cu_transquant_bypass_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cu_skip_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cu_qp_delta_abs>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cu_chroma_qp_offset_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<cu_chroma_qp_offset_idx>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<pred_mode_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<part_mode>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<prev_intra_luma_pred_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<intra_chroma_pred_mode>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<merge_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<merge_idx>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<inter_pred_idc>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<ref_idx_lX>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<abs_mvd_greater0_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<abs_mvd_greater1_flag>*>(this));
        v.visitContextsBase(static_cast<BaseTemplate<mvp_lX_flag>*>(this));
        this->HevcContextsResidual<BaseTemplate>::accept(v);
    }
};


struct Initialize
{
    Initialize(int qp, int initType) : qp(qp), initType(initType) { }
    template <class T>
    void visitContextsBase(T *t)
    {
        *t = T(this->qp, this->initType);
    }
    int qp;
    int initType;
};


template <class C>
struct CheckSame
{
    CheckSame(C const &other, int initType) :
        other(other),
        initType(initType)
    {
    }

    template <class T>
    void visitContextsBase(T *t)
    {
        T const &other = this->other;

        t->checkSameAs(other, this->initType);
    }

    C const &other;
    int initType;
};


struct Contexts :
    HevcContexts<ContextsBase>,
    ValueHolder<StatCoeff>
{
    void initialize(int qp, int initType)
    {
        Initialize initialize(qp, initType);
        this->accept<Initialize>(initialize);
        this->initType = initType;
    }

    template <class Tag>
    ContextModel &get(int ctxInc)
    {
        auto size = this->ContextsBase<Tag>::cm.size();
        return this->ContextsBase<Tag>::cm[ctxInc];
    }

    int initType; // review: this increases size of contexts' copy - remove?

    void checkSameAs(const Contexts &other)
    {
        CheckSame<Contexts> checkSame(other, this->initType);
        this->accept(checkSame);
    }
};


struct ContextPrinter
{
    ContextPrinter(std::ostream &o, int initType) : o(o), initType(initType) { }
    template <class T>
    void visitContextsBase(T *t)
    {
        t->print(this->o, this->initType);
    }
    int initType;
    std::ostream &o;
};


static std::ostream& operator<<(std::ostream& o, Contexts& contexts)
{
    ContextPrinter contextPrinter(o, contexts.initType);
    contexts.accept(contextPrinter);
    return o;
}


static inline uint8_t rangeTabLPS(int pStateIdx, int qRangeIdx)
{
    static std::array<std::array<const uint8_t, 4>, 64> table = { {
        { 128, 176, 208, 240 },
        { 128, 167, 197, 227 },
        { 128, 158, 187, 216 },
        { 123, 150, 178, 205 },
        { 116, 142, 169, 195 },
        { 111, 135, 160, 185 },
        { 105, 128, 152, 175 },
        { 100, 122, 144, 166 },
        { 95, 116, 137, 158 },
        { 90, 110, 130, 150 },
        { 85, 104, 123, 142 },
        { 81, 99, 117, 135 },
        { 77, 94, 111, 128 },
        { 73, 89, 105, 122 },
        { 69, 85, 100, 116 },
        { 66, 80, 95, 110 },
        { 62, 76, 90, 104 },
        { 59, 72, 86, 99 },
        { 56, 69, 81, 94 },
        { 53, 65, 77, 89 },
        { 51, 62, 73, 85 },
        { 48, 59, 69, 80 },
        { 46, 56, 66, 76 },
        { 43, 53, 63, 72 },
        { 41, 50, 59, 69 },
        { 39, 48, 56, 65 },
        { 37, 45, 54, 62 },
        { 35, 43, 51, 59 },
        { 33, 41, 48, 56 },
        { 32, 39, 46, 53 },
        { 30, 37, 43, 50 },
        { 29, 35, 41, 48 },
        { 27, 33, 39, 45 },
        { 26, 31, 37, 43 },
        { 24, 30, 35, 41 },
        { 23, 28, 33, 39 },
        { 22, 27, 32, 37 },
        { 21, 26, 30, 35 },
        { 20, 24, 29, 33 },
        { 19, 23, 27, 31 },
        { 18, 22, 26, 30 },
        { 17, 21, 25, 28 },
        { 16, 20, 23, 27 },
        { 15, 19, 22, 25 },
        { 14, 18, 21, 24 },
        { 14, 17, 20, 23 },
        { 13, 16, 19, 22 },
        { 12, 15, 18, 21 },
        { 12, 14, 17, 20 },
        { 11, 14, 16, 19 },
        { 11, 13, 15, 18 },
        { 10, 12, 15, 17 },
        { 10, 12, 14, 16 },
        { 9, 11, 13, 15 },
        { 9, 11, 12, 14 },
        { 8, 10, 12, 14 },
        { 8, 9, 11, 13 },
        { 7, 9, 11, 12 },
        { 7, 9, 10, 12 },
        { 7, 8, 10, 11 },
        { 6, 8, 9, 11 },
        { 6, 7, 9, 10 },
        { 6, 7, 8, 9 },
        { 2, 2, 2, 2 } } };

    return table[pStateIdx][qRangeIdx];
}

#endif
