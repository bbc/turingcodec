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

#ifndef INCLUDED_Snake_h
#define INCLUDED_Snake_h

#pragma once

#include <array>
#include <vector>
#include <type_traits>
#include <cstring>


#if defined(_DEBUG)
// When SNAKE_DEBUG is set, each entry has an associated diagonal position that can be used to check
// the location of the snake is as expected.
# define SNAKE_DEBUG
#endif


#ifdef SNAKE_DEBUG
# ifdef WIN32
#  define SNAKE_ASSERT(a) if (!(a)) __debugbreak();
# else
#  define SNAKE_ASSERT(a) assert(a)
# endif
#else
# define SNAKE_ASSERT(a)
#endif


struct Backup{};
struct Restore{};
struct Check{};


// A snake describes state along a the advancing edge of a wavefront. Each entry in a snake array represents
// one point along the wavefront. As the wavefront advances over the picture, entries  are overwritten by
// data points to their bottom right. T is the stored entry type.
// Review: coherency of x0, y0 and origin when snake is offset from (0, 0)
template <class T>
struct Snake
{
    // A pointer into a snake buffer.
    struct Pointer
    {
        typedef T Type;

        Pointer() { }

#ifdef SNAKE_DEBUG
        T *begin;
        T *end;
        int log2Resolution;
        int x0; //< scaled X offset
        int y0; //< scaled Y offset
        int *originPosition;
        int positionOffset;

        void checkPosition(int x, int y, int log2Resolution = 0) const
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            x >>= log2Resolution;
            y >>= log2Resolution;
            T &entry = this->origin[x - y - 1];
            SNAKE_ASSERT(&entry >= this->begin);
            SNAKE_ASSERT(&entry < this->end);
            int &entryPosition = this->originPosition[x - y - 1];
            if (entryPosition != x + y)
            {
                // simultaneous equation to obtain actual coordinates
                // ax - ay - 1 = x - y - 1
                // ax + ay = entryPosition
                // 2*ax - 1 = x - y - 1 + entryPosition
                // 2*ax = x - y + entryPosition
                const int ax = (x - y + entryPosition) / 2;
                const int ay = entryPosition - ax;
                const int round = (1 << log2Resolution) - 1;
                std::cout
                << "expected position (" << (x << log2Resolution) + round << ", " << (y << log2Resolution) + round << ")"
                << " but snake actually at (" << (ax << log2Resolution) + round << ", " << (ay << log2Resolution) + round << ")\n";
                SNAKE_ASSERT(false);
            }
        }

        void forcePosition(int x, int y, int log2Resolution = 0)
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            x >>= log2Resolution;
            y >>= log2Resolution;
            T &entry = this->origin[x - y - 1];
            SNAKE_ASSERT(&entry >= this->begin);
            SNAKE_ASSERT(&entry < this->end);
            int &entryPosition = this->originPosition[x - y - 1];
            entryPosition = x + y;
        }

        template <class Rectangle>
        void checkPosition(When when, const Rectangle &rectangular) const
        {
            const auto x = xPositionOf(rectangular);
            const auto y = yPositionOf(rectangular);
            const auto w = widthOf(rectangular);
            const auto h = heightOf(rectangular);

            const int width = w >> this->log2Resolution;
            const int height = h >> this->log2Resolution;

            if (when == after)
            {
                const int x0 = x >> this->log2Resolution;
                const int y0 = y >> this->log2Resolution;
                for (int x = 0; x < width; ++x)
                {
                    this->checkPosition((x0 + x) << this->log2Resolution, (y0 + height - 1) << this->log2Resolution, this->log2Resolution);
                }

                for (int y = height - 2; y >= 0; --y)
                {
                    this->checkPosition((x0 + width - 1) << this->log2Resolution, (y0 + y) << this->log2Resolution, this->log2Resolution);
                }
            }
            else if (when == before)
            {
                const int x0 = x >> this->log2Resolution;
                const int y0 = y >> this->log2Resolution;
                for (int x = -1; x < width; ++x)
                {
                    this->checkPosition((x0 + x) << this->log2Resolution, (y0 - 1) << this->log2Resolution, this->log2Resolution);
                }

                for (int y = height - 1; y >= 0; --y)
                {
                    this->checkPosition((x0 - 1) << this->log2Resolution, (y0 + y) << this->log2Resolution, this->log2Resolution);
                }
            }
            else if (when == corner)
            {
                const int x0 = x >> this->log2Resolution;
                const int y0 = y >> this->log2Resolution;
                this->checkPosition((x0 - 1) << this->log2Resolution, (y0 - 1) << this->log2Resolution, this->log2Resolution);
            }
        }

        template <class Rectangle>
        void forcePosition(When when, const Rectangle &rectangular)
        {
            const auto x = xPositionOf(rectangular);
            const auto y = yPositionOf(rectangular);
            const auto w = widthOf(rectangular);
            const auto h = heightOf(rectangular);

            const int width = w >> this->log2Resolution;
            const int height = h >> this->log2Resolution;

            if (when == after)
            {
                const int x0 = x >> this->log2Resolution;
                const int y0 = y >> this->log2Resolution;
                for (int x = 0; x < width; ++x)
                {
                    this->forcePosition((x0 + x) << this->log2Resolution, (y0 + height - 1) << this->log2Resolution, this->log2Resolution);
                }

                for (int y = height - 2; y >= 0; --y)
                {
                    this->forcePosition((x0 + width - 1) << this->log2Resolution, (y0 + y) << this->log2Resolution, this->log2Resolution);
                }
            }
            else if (when == before)
            {
                const int x0 = x >> this->log2Resolution;
                const int y0 = y >> this->log2Resolution;
                for (int x = -1; x < width; ++x)
                {
                    this->forcePosition((x0 + x) << this->log2Resolution, (y0 - 1) << this->log2Resolution, this->log2Resolution);
                }

                for (int y = height - 1; y >= 0; --y)
                {
                    this->forcePosition((x0 - 1) << this->log2Resolution, (y0 + y) << this->log2Resolution, this->log2Resolution);
                }
            }
            else if (when == corner)
            {
                const int x0 = x >> this->log2Resolution;
                const int y0 = y >> this->log2Resolution;
                this->forcePosition((x0 - 1) << this->log2Resolution, (y0 - 1) << this->log2Resolution, this->log2Resolution);
            }
        }
#else
        template <class Rectangle>
        void checkPosition(When when, const Rectangle &rectangular) const { }

        template <class Rectangle>
        void forcePosition(When when, const Rectangle &rectangular) const { }
#endif

        T *origin;

        template <class Position>
        const T &get(int x, int y, int log2Resolution) const
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            return this->origin[(x >> log2Resolution) - (y >> log2Resolution) - 1];
        }

        T &operator()(int x, int y)
        {
            return this->origin[x - y - 1];
        }

        const T &operator()(int x, int y) const
        {
            return this->origin[x - y - 1];
        }

        const T &entry(int x, int y, int log2Resolution = 0) const
        {
#ifdef SNAKE_DEBUG
            checkPosition(x, y, log2Resolution);
#endif
            x >>= log2Resolution;
            y >>= log2Resolution;
            return this->origin[x - y - 1];
        }

        const T &at(int x, int y, int log2Resolution = 0) const
        {
            return this->entry(x, y, log2Resolution);
        }

        void commit(const T &t, int x, int y, int log2Resolution = 0)
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            x >>= log2Resolution;
            y >>= log2Resolution;
            T &entry = this->origin[x - y - 1];
            SNAKE_ASSERT(&entry >= this->begin);
            SNAKE_ASSERT(&entry < this->end);
            static_cast<T &>(entry) = t;
#ifdef SNAKE_DEBUG
            this->originPosition[x - y - 1] = x + y;
#endif
        }

        template <class Rectangular, class UnaryFunction>
        void foreach(const Rectangular &rectangular, int log2Resolution, int padX, int padY, UnaryFunction f)
        {
            const auto x0 = xPositionOf(rectangular);
            const auto y0 = yPositionOf(rectangular);
            const auto width = widthOf(rectangular);
            const auto height = heightOf(rectangular);

            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            SNAKE_ASSERT((x0 & ((1 << log2Resolution) - 1)) == 0);
            SNAKE_ASSERT((y0 & ((1 << log2Resolution) - 1)) == 0);
            SNAKE_ASSERT((width & ((1 << log2Resolution) - 1)) == 0);
            SNAKE_ASSERT((height & ((1 << log2Resolution) - 1)) == 0);

            const int resolution = 1 << log2Resolution;

            padX <<= log2Resolution;
            padY <<= log2Resolution;

            for (int x = -resolution; x < width + padX; x += resolution)
            {
                T &entry = this->origin[((x0 + x) >> log2Resolution) - ((y0 - resolution) >> log2Resolution) - 1];
                f(entry);
            }

            for (int y = 0; y < height + padY; y += resolution)
            {
                T &entry = this->origin[((x0 - resolution) >> log2Resolution) - ((y0 + y) >> log2Resolution) - 1];
                f(entry);
            }
        }

        template <class Rectangular>
        void copyBlockFrom(const Pointer &other, When when, const Rectangular &rectangular, int log2Resolution = 0, int padX = 0, int padY = 0)
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            SNAKE_ASSERT(log2Resolution == other.log2Resolution);

            const auto x = xPositionOf(rectangular);
            const auto y = yPositionOf(rectangular);
            const auto w = widthOf(rectangular);
            const auto h = heightOf(rectangular);

            const int width = w >> log2Resolution;
            const int height = h >> log2Resolution;

#ifdef SNAKE_DEBUG

            other.checkPosition(when, rectangular);

            for (int i = -1 - height - padY; i < width + padX; ++i)
            {
                auto p1 = &this->origin[(x >> log2Resolution) - (y >> log2Resolution) + i];
                SNAKE_ASSERT(p1 >= this->begin);
                SNAKE_ASSERT(p1 < this->end);

                auto p2 = &other.origin[(x >> log2Resolution) - (y >> log2Resolution) + i];
                SNAKE_ASSERT(p2 >= other.begin);
                SNAKE_ASSERT(p2 < other.end);

                *p1 = *p2;
                auto pos1 = &this->originPosition[(x >> log2Resolution) - (y >> log2Resolution) + i];
                auto pos2 = &other.originPosition[(x >> log2Resolution) - (y >> log2Resolution) + i];

                *pos1 = *pos2;
            }
            this->checkPosition(when, rectangular);
#else
            const int begin = ((x - y) >> log2Resolution) - 1 - height - padY;
            const int n = width + padX + 1 + height + padY;
            memcpy(&this->origin[begin], &other.origin[begin], n * sizeof(T));
#endif
        }

        void copy(const Pointer &other, int width, int height, int log2Resolution = 0, int padX = 0, int padY = 0)
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            SNAKE_ASSERT(log2Resolution == other.log2Resolution);

            width >>= this->log2Resolution;
            height >>= this->log2Resolution;

            for (int i = 1 - height - padY; i < width + padX; ++i)
            {
                auto p = &this->origin[i - 1];
                SNAKE_ASSERT(p >= this->begin);
                SNAKE_ASSERT(p < this->end);
                *p = other.origin[i - 1];
            }
        }

        template <class Rectangular>
        void print(std::ostream &o, Rectangular const &rectangular, int log2Resolution = 0, int padX = 0, int padY = 0, bool singleLine=false)
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            const auto x = xPositionOf(rectangular);
            const auto y = yPositionOf(rectangular);
            const auto w = widthOf(rectangular);
            const auto h = heightOf(rectangular);

            const int width = w >> log2Resolution;
            const int height = h >> log2Resolution;

            for (int i = -1 - height - padY; i < width + padX; ++i)
            {
#ifdef SNAKE_DEBUG
                if (!singleLine)
                {
                    auto &entryPosition = this->originPosition[(x >> log2Resolution) - (y >> log2Resolution) + i];
                    o << entryPosition << ": ";
                }
#endif
                auto &entry = this->origin[(x >> log2Resolution) - (y >> log2Resolution) + i];

                printSnakeEntry(o, entry);

                o << (singleLine ? ", " : "\n");
            }

            if (singleLine) o << "\n";
        }

        // Writes a rectangle's worth of value to the snake.
        template <class Rectangle>
        void commitRectangle(const Rectangle &rectangular, const T &value, int log2Resolution, bool force=false)
        {
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            const auto x0 = xPositionOf(rectangular) >> log2Resolution;
            const auto y0 = yPositionOf(rectangular) >> log2Resolution;
            const auto width = widthOf(rectangular) >> log2Resolution;
            const auto height = heightOf(rectangular) >> log2Resolution;

#ifdef SNAKE_DEBUG
            for (int x = -1; x < width && !force; ++x)
            {
                this->checkPosition((x0 + x) << this->log2Resolution, (y0 - 1) << log2Resolution, log2Resolution);
            }

            for (int y = height - 1; y >= 0 && !force; --y)
            {
                this->checkPosition((x0 - 1) << this->log2Resolution, (y0 + y) << log2Resolution, log2Resolution);
            }
#endif

            for (int x = 0; x < width; ++x)
            {
                this->commit(value, (x0 + x) << log2Resolution, (y0 + height - 1) << log2Resolution, log2Resolution);
            }

            for (int y = height - 2; y >= 0; --y)
            {
                this->commit(value, (x0 + width - 1) << log2Resolution, (y0 + y) << log2Resolution, log2Resolution);
            }
        }

        Pointer offset(int x0, int y0, int log2Resolution = 0)
        {
            Pointer pointer = *this;
#ifdef SNAKE_DEBUG
            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            SNAKE_ASSERT(x0 == x0 >> log2Resolution << log2Resolution);
            SNAKE_ASSERT(y0 == y0 >> log2Resolution << log2Resolution);
            pointer.x0 += x0 >> log2Resolution;
            pointer.y0 += y0 >> log2Resolution;
            pointer.positionOffset += (x0 - y0) >> log2Resolution;
#endif
            pointer.origin += (x0 - y0) >> log2Resolution;
            return pointer;
        }

        // direction = 'R' for right, 'D' for down
        template <char direction>
        void pad(int x0, int y0, int n, int log2Resolution = 0)
        {
            Pointer pointer = *this;
            pointer.origin += (x0 - y0) >> log2Resolution;
            n >>= log2Resolution;
            auto const value = *pointer.origin;
            if (direction == 'R')
            {
                while (n--)
                    *(++pointer.origin) = value;
            }
            else
            {
                while (n--)
                    *(--pointer.origin) = value;
            }
        }

        // Copys new values from a 2-D array to the snake. The values are copied from the bottom and right
        // edges of the shape described by rectangular
        // TwoDimensionalAccess must have member "const T &operator()(int x, int y)" or equivalent.
        template <class Rectangular, class TwoDimensionalAccess>
        void copyFrom2D(const Rectangular &rectangular, const TwoDimensionalAccess &src, int log2Resolution = 0)
        {
            //std::cout << "copyFrom2D " << rectangular << " " << log2Resolution << "\n";
            this->checkPosition(before, rectangular);

            int x0 = xPositionOf(rectangular);
            int y0 = yPositionOf(rectangular);
            int width = widthOf(rectangular);
            int height = widthOf(rectangular);

            SNAKE_ASSERT(log2Resolution == this->log2Resolution);
            SNAKE_ASSERT(x0 == x0 >> log2Resolution << log2Resolution);
            SNAKE_ASSERT(y0 == y0 >> log2Resolution << log2Resolution);
            SNAKE_ASSERT(width == width >> log2Resolution << log2Resolution);
            SNAKE_ASSERT(height == height >> log2Resolution << log2Resolution);

            x0 >>= log2Resolution;
            y0 >>= log2Resolution;
            width >>= log2Resolution;
            height >>= log2Resolution;

            {
                // copy the bottom edge of the block
                const int y = y0 + height - 1;
#ifdef SNAKE_DEBUG
                for (int x = x0; x < x0 + width; ++x)
                {
                    this->commit(src(x, y), x << log2Resolution, y << log2Resolution, log2Resolution);
                }
#else
                memcpy(&this->origin[x0 - y - 1], &src(x0, y), width * sizeof(T));
#endif
            }

            {
                // copy the left edge of the block
                const int x = x0 + width - 1;
                for (int y = y0 + height - 1; y >= y0; --y)
                {
                    this->commit(src(x, y), x << log2Resolution, y << log2Resolution, log2Resolution);
                }
            }
        }
    };

    // Structure allowing convenient access to snake data whilst working on a rectangular area of
    // picture. Via Cursor, callers can read and write the current rectangular's value. Values above and
    // left, and to the above-left of the rectangular are available to read.
    struct Cursor :
        Pointer
        {
            template <class Rectangular>
            void relocate(Pointer &snake, const Rectangular &rectangular, int log2Resolution = 0, bool force = false)
            {
                if (!force) snake.checkPosition(corner, rectangular);
                static_cast<Pointer &>(*this) = snake.offset(0, 0, log2Resolution);
                this->x0 = xPositionOf(rectangular);
                this->y0 = yPositionOf(rectangular);
#ifdef SNAKE_DEBUG
                this->width = widthOf(rectangular);
                this->height =heightOf(rectangular);
#endif
                SNAKE_ASSERT(log2Resolution == this->log2Resolution);
                SNAKE_ASSERT((xPositionOf(rectangular) & ((1 << this->log2Resolution) - 1)) == 0);
                SNAKE_ASSERT((yPositionOf(rectangular) & ((1 << this->log2Resolution) - 1)) == 0);
                SNAKE_ASSERT((widthOf(rectangular) & ((1 << this->log2Resolution) - 1)) == 0);
                SNAKE_ASSERT((heightOf(rectangular) & ((1 << this->log2Resolution) - 1)) == 0);
                this->p = &snake.origin[(xPositionOf(rectangular) - yPositionOf(rectangular)) >> log2Resolution];
            }

            template <class Position>
            const T &offset(int dx, int dy, int log2Resolution = 0) const
            {
                SNAKE_ASSERT(log2Resolution == this->log2Resolution);
                dx >>= log2Resolution;
                dy >>= log2Resolution;

                if (std::is_same<Position, Left>::value)
                {
                    SNAKE_ASSERT(dx == -1);
                    return this->p[-1 - dy - 1];
                }

                if (std::is_same<Position, Up>::value)
                {
                    SNAKE_ASSERT(dy == -1);
                    return this->p[dx];
                }

                if (std::is_same<Position, Current>::value)
                {
                    SNAKE_ASSERT(dx == 0);
                    SNAKE_ASSERT(dy == 0);
                    return this->value;
                }

                if (dx < 0)
                {
                    SNAKE_ASSERT(dx == -1 || dy >= (this->height >> log2Resolution));
                    SNAKE_ASSERT(dy >= -1);
                    return this->p[dx - dy - 1];
                }
                else if (dy < 0)
                {
                    SNAKE_ASSERT(dx >= -1);
                    SNAKE_ASSERT(dy == -1 || dx >= (this->width >> log2Resolution));
                    return this->p[dx - dy - 1];
                }
                else
                {
                    return this->value;
                }
            }

            T &current(int dx, int dy, int log2Resolution = 0)
            {
                SNAKE_ASSERT(log2Resolution == this->log2Resolution);
                return this->value;
            }

            template <class Rectangle>
            void commit(const Rectangle &rectangular, int log2Resolution, bool force=false)
            {
                this->Pointer::commitRectangle(rectangular, this->value, log2Resolution, force);
            }

            T *p; //< review: could this become a SNAKE_DEBUG parameter?
            T value;
            int x0, y0; //< review: could these become SNAKE_DEBUG parameters?
#ifdef SNAKE_DEBUG
            int width, height;
#endif
        };

    // Represents a backup that can restore original state of the snake after processing
    // a rectangular area of picture.
    template <int maxWidth, int maxHeight>
    struct Backup
    {
        T buffer[maxWidth + maxHeight];
#ifdef SNAKE_DEBUG
        int bufferPosition[maxWidth + maxHeight];
        int x, y, width, height;
#endif

        template <class What, class Rectangle>
        void perform(Pointer &snake, const Rectangle &rectangular, int log2Resolution)
        {
            int x = xPositionOf(rectangular);
            int y = yPositionOf(rectangular);
            int width = widthOf(rectangular);
            int height = heightOf(rectangular);

            T *p = &snake.origin[((x - (y + height - (1 << log2Resolution))) >> log2Resolution) - 1];
#ifdef SNAKE_DEBUG
            int *pPosition = &snake.originPosition[((x - (y + height - (1 << log2Resolution))) >> log2Resolution) - 1];

            if (std::is_same<What, ::Backup>::value)
            {
                this->x = x;
                this->y = y;
                this->width = width;
                this->height = height;
            }
            else
            {
                SNAKE_ASSERT(this->x == x);
                SNAKE_ASSERT(this->y == y);
                SNAKE_ASSERT(this->width == width);
                SNAKE_ASSERT(this->height == height);
            }
#endif

            SNAKE_ASSERT(log2Resolution == snake.log2Resolution);
            x >>= log2Resolution;
            y >>= log2Resolution;
            width >>= log2Resolution;
            height >>= log2Resolution;
            SNAKE_ASSERT(width <= maxWidth);
            SNAKE_ASSERT(height <= maxHeight);

            for (int i = 0; i < width + height - 1; ++i)
            {
                if (std::is_same<What, ::Backup>::value)
                {
                    this->buffer[i] = p[i];
#ifdef SNAKE_DEBUG
                    this->bufferPosition[i] = pPosition[i];
#endif
                }
                if (std::is_same<What, ::Restore>::value)
                {
                    p[i] = this->buffer[i];
#ifdef SNAKE_DEBUG
                    pPosition[i] = this->bufferPosition[i];
#endif
                }
                if (std::is_same<What, ::Check>::value)
                {
                    SNAKE_ASSERT(this->buffer[i] == p[i]);
                    SNAKE_ASSERT(this->bufferPosition[i] == pPosition[i]);
                }
            }
        }
    };

    template<int padX, int padY>
    struct Vector :
        Pointer
        {
            template <class Rectangular>
            void resize(const Rectangular &rectangular, int log2Resolution)
            {
                const auto width = widthOf(rectangular);
                const auto height = heightOf(rectangular);
                const auto x = xPositionOf(rectangular);
                const auto y = yPositionOf(rectangular);

                const int h = (height >> log2Resolution) + padY;
                const int w = (width >> log2Resolution) + padX;

                this->buffer.resize(h + 1 + w, T());
                this->origin = &buffer[h + 1];
                this->origin -= (x >> log2Resolution);
                this->origin += (y >> log2Resolution);

#ifdef SNAKE_DEBUG
                this->bufferPosition.resize(h + 1 + w);
                this->originPosition = &bufferPosition[h + 1];
                this->originPosition -= (x >> log2Resolution);
                this->originPosition += (y >> log2Resolution);
                this->positionOffset = (x - y) >> log2Resolution;

                const int resolution = 1 << log2Resolution;
                SNAKE_ASSERT(width % resolution == 0);
                SNAKE_ASSERT(height % resolution == 0);
                SNAKE_ASSERT(x % resolution == 0);
                SNAKE_ASSERT(y % resolution == 0);

                this->begin = &this->buffer.front();
                this->end = this->begin + this->buffer.size();
                this->log2Resolution = log2Resolution;
                this->x0 = x >> log2Resolution;
                this->y0 = y >> log2Resolution;

                for (int dx = -1; dx < w; ++dx)
                {
                    int dy = -1;
                    this->originPosition[this->x0 + dx - this->y0 - dy - 1] = this->x0 + dx + this->y0 + dy;
                }
                for (int dy = 0; dy < h; ++dy)
                {
                    int dx = -1;
                    this->originPosition[this->x0 + dx - this->y0 - dy - 1] = this->x0 + dx + this->y0 + dy;
                }
#endif
            }

        private:
            std::vector<T> buffer;
#ifdef SNAKE_DEBUG
            std::vector<int> bufferPosition;
#endif
        };

    template <int maxPositionsX, int maxPositionsY, int padX, int padY>
    struct Array :
        Pointer
        {
            template <class Rectangular>
            void resize(const Rectangular &rectangular, int log2Resolution)
            {
                const auto width = widthOf(rectangular);
                const auto height = heightOf(rectangular);
                const auto x = xPositionOf(rectangular);
                const auto y = yPositionOf(rectangular);

                const int w = (width >> log2Resolution) + padX;
                const int h = (height >> log2Resolution) + padY;
                this->origin = &buffer[0];
                this->origin += h + 1 + (y >> log2Resolution) - (x >> log2Resolution);
                this->origin = reinterpret_cast<T *>((reinterpret_cast<intptr_t>(this->origin) + 31) & ~31); // alignment of origin

#ifdef SNAKE_DEBUG
                this->originPosition = &bufferPosition[0];
                this->originPosition += h + 1 + (y >> log2Resolution) - (x >> log2Resolution);
                this->positionOffset = (x - y) >> log2Resolution;

                const int resolution = 1 << log2Resolution;
                SNAKE_ASSERT(width % resolution == 0);
                SNAKE_ASSERT(height % resolution == 0);
                SNAKE_ASSERT(x % resolution == 0);
                SNAKE_ASSERT(y % resolution == 0);

                this->begin = &buffer[0];
                this->end = this->begin + this->buffer.size();
                this->log2Resolution = log2Resolution;
                this->x0 = x >> log2Resolution;
                this->y0 = y >> log2Resolution;

                for (int dx = -1; dx < w; ++dx)
                {
                    int dy = -1;
                    this->originPosition[this->x0 + dx - this->y0 - dy - 1] = this->x0 + dx + this->y0 + dy;
                }
                for (int dy = 0; dy < h; ++dy)
                {
                    int dx = -1;
                    this->originPosition[this->x0 + dx - this->y0 - dy - 1] = this->x0 + dx + this->y0 + dy;
                }
#endif
            }

        private:
            alignas(32) std::array<T, padY + maxPositionsY + 1 + maxPositionsX + padX + 31> buffer;
#ifdef SNAKE_DEBUG
            std::array<int, padY + maxPositionsY + 1 + maxPositionsX + padX> bufferPosition;
#endif
        };
};

#endif
