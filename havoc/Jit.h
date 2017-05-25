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

#pragma once

#define NOMINMAX

#include "xbyak.h"
#include "havoc.h"
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <iostream>
#include <algorithm>


namespace Jit {

struct Buffer
{
    Buffer(size_t n, havoc_instruction_set isa = HAVOC_NONE)
        :
        buffer(n),
        isa(isa)
    {
        this->update(&this->buffer[0]);
        Xbyak::CodeArray::protect(this->p, this->nRemainingBytes(), true);
    }

    void update(uint8_t *p)
    {
        this->p = Xbyak::CodeArray::getAlignedAddress(p);
        assert(this->nRemainingBytes() >= 0);
    }

    void increment(size_t n)
    {
        this->update(this->p + n);
    }

    intptr_t nRemainingBytes() const
    {
        return &this->buffer[0] + this->buffer.size() - this->p;
    }

    uint8_t *pointer() const
    {
        return this->p;
    }

    havoc_instruction_set const isa;

private:
    std::vector<uint8_t> buffer;
    uint8_t *p;
};


template <typename FunctionType>
struct CountArguments;


template<typename R, typename ...Args>
struct CountArguments<R(Args...)>
{
    static const size_t value = sizeof...(Args);
};


struct Function :
    Xbyak::CodeGenerator
{
    int nArguments;

    Function(Buffer *buffer, int nArguments) :
        buffer(buffer),
        Xbyak::CodeGenerator(buffer ? buffer->nRemainingBytes() : 4096, buffer ? buffer->pointer() : 0),
        nArguments(nArguments)
    {
    }

    int pass = 1;

    virtual void assemble() = 0;
    virtual void data() { };

    void build()
    {
        assert(pass == 1);

        this->assemble();
        this->data();

        if (this->getSize() == 0) return;

        if (this->debug)
        {
            std::cerr << "function has "
                << this->nArguments << " arguments, "
                << this->variables << " variables and uses the first "
                << this->mmRegisters << " XMM registers and requies "
                << this->stackSize << " stack bytes\n";
        }

        this->reset();
        this->pass = 2;

        this->prologue(this->variables, this->mmRegisters, this->stackSize);
        this->variables = 0;
        this->mmRegisters = 0;
        this->stackSize = 0;
        this->assemble();
        this->epilogue(this->variables, this->mmRegisters, this->stackSize);
        this->data();

        this->buffer->increment(this->getSize());
    }

    void buildSinglePass(int variables, int mmRegisters, int stackSize)
    {
        this->pass = 0; // indicate special mode
        this->prologue(variables, mmRegisters, stackSize);
        this->variables = 0;
        this->mmRegisters = 0;
        this->stackSize = 0;
        this->assemble();
        this->epilogue(variables, mmRegisters, stackSize);
        this->data();
        this->buffer->increment(this->getSize());
    }

    template <typename F>
    operator F *()
    {
        if (this->getSize() == 0) return 0;
        assert(CountArguments<F>::value == this->nArguments);
        return reinterpret_cast<F *>(this->getCode());
    }

protected:
    int frameSize;
    int stackOffset = 0;
    int stackSize = 0;

#ifdef WIN32
    int const mmRegistersVolatile = 8;
#else
    int const mmRegistersVolatile = 16;
#endif


    void prologue(int variables, int mmRegisters, int stack)
    {
        assert(mmRegisters <= 16);
        //db({ 0xcc });
        //nop();
        //nop();

        this->stackOffset = 0x8; // rbp
#ifdef WIN32
        this->stackOffset += 0x20;  // shadow space
#endif

        this->frameSize = 0;

        //std::cout << "nArguments " << nArguments << "\n";
        //std::cout << "variables " << variables << "\n";
        //std::cout << "stack " << stack << "\n";

        int registers = nArguments + variables;
#ifdef WIN32
        for (int i = 7; i < registers; ++i)
#else
        for (int i = 9; i < registers; ++i)
#endif
        {
            push(reg64(i));
            this->frameSize += 8;
        }

        if (mmRegisters > mmRegistersVolatile || stack)
        {
            if (this->frameSize & 8) this->stackOffset += 8;
            if (mmRegisters > mmRegistersVolatile) this->stackOffset += (mmRegisters - mmRegistersVolatile) * 16;
            this->stackOffset += stack;
            this->frameSize += this->stackOffset;

            sub(rsp, this->stackOffset);
#ifdef WIN32
            // chkstk implementation
            for (int d = this->stackOffset - 4096; d >= 0; d -= 4096)
            {
                mov(ptr[rsp + d], rsp);
            }
#endif
        }

#ifdef WIN32
        // Store xmm6 and xmm7 in shadow space
        if (mmRegisters > 6) vmovaps(ptr[rsp + this->frameSize + 8], xmm6);
        if (mmRegisters > 7) vmovaps(ptr[rsp + this->frameSize + 24], xmm7);
#endif

        // Store callee-saved XMM registers in allocated stack
        for (int i = mmRegistersVolatile; i < mmRegisters; ++i)
        {
            vmovaps(ptr[rsp + stack + (i - 6) * 0x10], regXmm(i));
        }

#ifdef WIN32
        for (int i = 4; i < nArguments; ++i)
        {
            mov(arg64(i), ptr[rsp + this->frameSize + 8 * (i + 1)]);
        }
#else
        for (int i = 6; i < nArguments; ++i)
        {
            mov(arg64(i), ptr[rsp + this->frameSize + 8 * (i - 6 + 1)]);
        }
#endif
    }

    void epilogue(int variables, int mmRegisters, int stack)
    {
#ifdef WIN32
        if (mmRegisters > 6) vmovaps(xmm6, ptr[rsp + this->frameSize + 8]);
        if (mmRegisters > 7) vmovaps(xmm7, ptr[rsp + this->frameSize + 24]);
#endif
        for (int i = mmRegistersVolatile; i < mmRegisters; ++i)
        {
            vmovaps(regXmm(i), ptr[rsp + stack + (i - 6) * 0x10]);
        }

        if (mmRegisters > mmRegistersVolatile || stack)
        {
            add(rsp, this->stackOffset);
        }

        int registers = nArguments + variables;

#ifdef WIN32
        for (int i = registers - 1; i >= 7; --i)
#else
        for (int i = registers - 1; i >= 9; --i)
#endif
        {
            pop(reg64(i));
        }
        ret();
    }

    //virtual void appendix() { }

    Xbyak::Reg64 const &reg64(int i)
    {
        if (pass) this->variables = std::max(this->variables, i - nArguments + 1);
#ifdef WIN32
        Xbyak::Reg64 const *lookup[15] = { &rcx, &rdx, &r8, &r9, &r10, &r11, &rax, &rdi, &rsi, &rbx, &rbp, &r12, &r13, &r14, &r15 };
#else
        //  The first six integer or pointer arguments are passed in registers RDI, RSI, RDX, RCX, R8, and R9.
        // Registers RBP, RBX, and R12-R15 are callee-save registers; all others must be saved by the caller if it wishes to preserve their values.
        Xbyak::Reg64 const *lookup[15] = { &rdi, &rsi, &rdx, &rcx, &r8, &r9, &rax, &r10, &r11, &rbx, &rbp, &r12, &r13, &r14, &r15 };
#endif
        return *lookup[i];
    }

    Xbyak::Reg64 const &arg64(int i)
    {
        assert(i < nArguments);
        return reg64(i);
    }

    Xbyak::Xmm const &regXmm(int i)
    {
        if (pass) this->mmRegisters = std::max(this->mmRegisters, i + 1);
        Xbyak::Xmm const *lookup[] = { &xmm0, &xmm1, &xmm2, &xmm3, &xmm4, &xmm5, &xmm6, &xmm7, &xmm8, &xmm9, &xmm10, &xmm11, &xmm12, &xmm13, &xmm14, &xmm15 };
        return *lookup[i];
    }

    template <class T> void db(std::initializer_list<T> args, int times = 1)
    {
        for (int i = 0; i < times; ++i)
            for (auto arg : args)
                Xbyak::CodeGenerator::db(arg);
    }

    template <class T> void dw(std::initializer_list<T> args, int times = 1)
    {
        for (int i = 0; i < times; ++i)
            for (auto arg : args)
                Xbyak::CodeGenerator::dw(arg);
    }

    template <class T> void dd(std::initializer_list<T> args, int times = 1)
    {
        for (int i = 0; i < times; ++i)
            for (auto arg : args)
                Xbyak::CodeGenerator::dd(arg);
    }

    template <class T> void dq(std::initializer_list<T> args, int times = 1)
    {
        for (int i = 0; i < times; ++i)
            for (auto arg : args)
                Xbyak::CodeGenerator::dq(arg);
    }

    void breakpoint()
    {
        db({ 0xcc });
    }

    bool debug = false;

    havoc_instruction_set isa() const
    {
        return this->buffer->isa;
    }

private:
    Buffer *buffer;
    int increment = 0;
    int variables = 0;
    int mmRegisters = 0;
};

inline int order(int a, int b, int c, int d)
{
    return ((a << 6) | (b << 4) | (c << 2) | d);
}


}

//
//namespace Jit2 {
//
//template <class T>
//struct Function;
//
//struct Buffer
//{
//	Buffer(havoc_instruction_set isa, size_t size) { }
//};
//
//
//template <typename T>
//struct Sad
//{
//	Sad(int bitDepth=8*sizeof(T), int width=0, int height=0) { }
//	typedef int Type(T *p, intptr_t stride);
//	Type *make(Buffer &buffer)
//	{
//	}
//};
//
//
//struct Test
//{
//	Test()
//	{
//		Buffer buffer(HAVOC_AVX2, 2000);
//		Sad<uint16_t>::Type *function = Sad<uint16_t>(10, 32, 32).make(buffer);
//	}
//};
//
//static Test test;
//
//}
