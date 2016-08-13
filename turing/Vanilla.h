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

// Operations (usually ones that affect state) common to read, write, encode and decode.

#ifndef INCLUDED_Vanilla_h
#define INCLUDED_Vanilla_h

#pragma once


#include "Global.h"
#include "Cabac.h"
#include "GlobalState.h"


struct PictureBegin { };

template <class F> struct Vanilla;

template <>
struct Vanilla<PictureBegin>
{
    template <class H> static void go(PictureBegin, H &h)
    {
        StatePictures *statePictures = h;
        statePictures->sliceHeaderDone(h);
    }
};


template <class H>
void sliceHeaderDone2(H &h)
{
    if (!h[dependent_slice_segment_flag()])
    {
        h[SliceAddrRs()] = h[slice_segment_address()];
    }

    if (h[first_slice_segment_in_pic_flag()])
    {
        std::vector<int> colWidth(h[num_tile_columns_minus1()]+1);
        if( h[uniform_spacing_flag()] )
        {
            for( int i = 0; i  <=  h[num_tile_columns_minus1()]; i++ )
                colWidth[ i ] = ( ( i + 1 ) * h[PicWidthInCtbsY()] ) / ( h[num_tile_columns_minus1()] + 1 ) -
                ( i * h[PicWidthInCtbsY()] ) / ( h[num_tile_columns_minus1()] + 1 );
        }
        else
        {
            colWidth[ h[num_tile_columns_minus1()] ] = h[PicWidthInCtbsY()];
            for( int i = 0; i < h[num_tile_columns_minus1()]; i++ ) {
                colWidth[ i ] = h[column_width_minus1(i)] + 1;
                colWidth[ h[num_tile_columns_minus1()] ]  -=  colWidth[ i ];
            }
        }

        std::vector<int> rowHeight(h[num_tile_rows_minus1()]+1);
        if( h[uniform_spacing_flag()] )
        {
            for( int i = 0; i  <=  h[num_tile_rows_minus1()]; i++ )
                rowHeight[ i ] = ( ( i + 1 ) * h[PicHeightInCtbsY()] ) / ( h[num_tile_rows_minus1()] + 1 ) -
                ( i * h[PicHeightInCtbsY()] ) / ( h[num_tile_rows_minus1()] + 1 );
        }
        else
        {
            rowHeight[ h[num_tile_rows_minus1()] ] = h[PicHeightInCtbsY()];
            for( int i = 0; i < h[num_tile_rows_minus1()]; i++ ) {
                rowHeight[ i ] = h[row_height_minus1(i)] + 1;
                rowHeight[ h[num_tile_rows_minus1()] ]  -=  rowHeight[ i ];
            }
        }

        h[ContainerOf<colBd>()].resize(h[num_tile_columns_minus1()] + 2);
        for(int i = 0; i <= h[num_tile_columns_minus1()]; i++ )
        {
            h[colBd( i + 1)] = h[colBd(i)] + colWidth[ i ];
        }

        h[ContainerOf<rowBd>()].resize(h[num_tile_rows_minus1()] + 2);
        for(int i = 0; i <= h[num_tile_rows_minus1()]; i++ )
        {
            h[rowBd( i + 1 )] = h[rowBd(i)] + rowHeight[ i ];
        }

        h[ContainerOf<CtbAddrRsToTs>()].resize(h[PicSizeInCtbsY()] + 2);
        for( int ctbAddrRS = 0; ctbAddrRS < h[PicSizeInCtbsY()]; ctbAddrRS++ )
        {
            const int tbX = ctbAddrRS %  h[PicWidthInCtbsY()];
            const int tbY = ctbAddrRS /  h[PicWidthInCtbsY()];
            int tileX;
            for( int i = 0; i <= h[num_tile_columns_minus1()]; i++ )
            {
                if( tbX >= h[colBd(i)] )
                {
                    tileX = i;
                }
            }
            int tileY;
            for( int j = 0; j <= h[num_tile_rows_minus1()]; j++ )
            {
                if( tbY >=  h[rowBd(j)])
                {
                    tileY = j;
                }
            }
            h[CtbAddrRsToTs(ctbAddrRS)] = 0;
            for( int i = 0; i < tileX; i++ )
            {
                h[CtbAddrRsToTs(ctbAddrRS)] += rowHeight[ tileY ] * colWidth[ i ];
            }
            for( int j = 0; j < tileY; j++ )
            {
                h[CtbAddrRsToTs(ctbAddrRS)] += h[PicWidthInCtbsY()] * rowHeight[ j ];
            }
            h[CtbAddrRsToTs(ctbAddrRS)] += ( tbY - h[rowBd( tileY)] ) * colWidth[ tileX ] + tbX - h[colBd(tileX) ];
        }

        h[ContainerOf<CtbAddrTsToRs>()].resize(h[PicSizeInCtbsY()] + 1);
        for( int ctbAddrRS = 0; ctbAddrRS < h[PicSizeInCtbsY()]; ctbAddrRS++ )
        {
            h[CtbAddrTsToRs(h[CtbAddrRsToTs(ctbAddrRS )])] = ctbAddrRS;
        }
        h[CtbAddrTsToRs(h[PicSizeInCtbsY()])] = h[PicSizeInCtbsY()];

        h[ContainerOf<TileId>()].resize(h[PicSizeInCtbsY()] + 2);
        h[TileId(0)] = -2;
        for( int j = 0, tIdx = 0; j <= h[num_tile_rows_minus1()]; j++ )
        {
            for( int i = 0; i <= h[num_tile_columns_minus1()]; i++, tIdx++ )
            {
                for( int y = h[rowBd( j )]; y < h[rowBd( j + 1) ]; y++ )
                {
                    for( int x = h[colBd( i )]; x < h[colBd( i + 1) ]; x++ )
                    {
                        h[TileId(h[CtbAddrRsToTs( y* h[PicWidthInCtbsY()]+ x ) ] ) ] = tIdx;
                    }
                }
            }
        }
        h[TileId(h[PicSizeInCtbsY()])] = -1;
    }

}


static const int tablesDs = 0;
static const int tablesWpp = 1;


template <>
struct Syntax<ContextsInitialize>
{
    template <class H> static void go(const ContextsInitialize &e, H &h)
    {
        Contexts& contexts = *static_cast<Contexts *>(h);
        contexts.initialize(h[SliceQpY()], h[initType()]);
        h[StatCoeff(0)] = 0;
        h[StatCoeff(1)] = 0;
        h[StatCoeff(2)] = 0;
        h[StatCoeff(3)] = 0;
    }
};

template <>
struct Syntax<ContextsSave>
{
    template <class H> static void go(const ContextsSave &e, H &h)
    {
        Contexts *contexts = h;
        StateSlice *stateSlice = h;

        stateSlice->savedContexts[e.i] = *contexts;
    }
};

template <>
struct Syntax<ContextsRestore>
{
    template <class H> static void go(const ContextsRestore &e, H &h)
    {
        Contexts *contexts = h;
        StateSlice *stateSlice = h;

        *contexts = stateSlice->savedContexts[e.i];
    }
};


// Local state that precomputes and caches the constant value scanIdx
struct ResidualCodingState :
    ValueCache<scanIdx>
    {
        template <class H>
        ResidualCodingState(H &h):
            ValueCache<scanIdx>(h)
            {
            }
    };


template <class H>
void preCtu(H &h)
{
    static_cast<AvailabilityCtu *>(h)->init(h);

    const bool codingTreeUnitIsFirstInSlice = h[SliceAddrRs()] == h[CtbAddrInRs()];
    const bool codingTreeUnitIsFirstInTile = h[CtbAddrInTs()] == 0 || h[TileId(h[CtbAddrInTs()])] != h[TileId(h[CtbAddrInTs()] - 1)];
    const bool codingTreeUnitIsFirstInRow = h[CtbAddrInRs()] % h[PicWidthInCtbsY()] == 0;

    const bool resetQp =
            codingTreeUnitIsFirstInSlice ||
            codingTreeUnitIsFirstInTile ||
            (codingTreeUnitIsFirstInRow && h[entropy_coding_sync_enabled_flag()] == 1);

    if (resetQp)
    {
        h[QpY()] = h[SliceQpY()];
    }

    const bool initializeCabac =
            h[CtbAddrInRs()] == h[slice_segment_address()] ||
            codingTreeUnitIsFirstInTile ||
            (codingTreeUnitIsFirstInRow && h[entropy_coding_sync_enabled_flag()] == 1);

    if (initializeCabac)
    {
        h(CabacRestart());

        if (codingTreeUnitIsFirstInTile)
        {
            h(ContextsInitialize());

            // start of tile: everything that went before is now unavailable
            StateSpatial &stateSpatial = *static_cast<StateSpatial *>(h);

            Turing::Rectangle rectangle;
            rectangle.x0 = 0;
            rectangle.y0 = 0;
            rectangle.width = h[PicWidthInCtbsY()];
            rectangle.height = h[PicHeightInCtbsY()];

            stateSpatial.snakeSaoCtuData.reset(rectangle, 0, 0, 0);

            rectangle.width <<= h[CtbLog2SizeY()];
            rectangle.height <<= h[CtbLog2SizeY()];

            Neighbourhood *neighbourhood = h;
            neighbourhood->snakeMerge.reset(rectangle, h[MinCbLog2SizeY()] - 1, 1, 1);
            neighbourhood->snake.reset(rectangle, h[MinCbLog2SizeY()] - 1, 1, 1);
        }
        else if (h[entropy_coding_sync_enabled_flag()] == 1 && codingTreeUnitIsFirstInRow)
        {
            const int x0 = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()];
            assert(x0 == 0);
            const int y0 = (h[CtbAddrInRs()] / h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()];
            const int xNbT = x0 + h[CtbSizeY()];
            const int yNbT = y0 - h[CtbSizeY()];
            const bool availableFlagT = h[availableX(x0, y0, xNbT, yNbT)];
            if (availableFlagT)
            {
                // synchronize WPP
                h(ContextsRestore(tablesWpp));
            }
            else
            {
                h(ContextsInitialize());
            }
        }
        else if (h[CtbAddrInRs()] == h[slice_segment_address()] && h[dependent_slice_segment_flag()])
        {
            // synchronize DS
            h(ContextsRestore(tablesDs));
        }
        else
        {
            h(ContextsInitialize());
        }

        if (h[CtbAddrInRs()] == h[slice_segment_address()] && !h[dependent_slice_segment_flag()])
        {
            // start of slice: everything that went before is now unavailable
            const int leftSlice = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]) << h[CtbLog2SizeY()];

            int left = 0;
            for (int i = 0; i<h[num_tile_columns_minus1()]; ++i)
            {
                const int colX = h[colBd(i + 1)] << h[CtbLog2SizeY()];
                if (colX > leftSlice) break;
                left = colX;
            }

            // review: duplicated
            StateSpatial &stateSpatial = *static_cast<StateSpatial *>(h);

            Turing::Rectangle rectangle;
            rectangle.x0 = 0;
            rectangle.y0 = 0;
            rectangle.width = h[PicWidthInCtbsY()];
            rectangle.height = h[PicHeightInCtbsY()];

            stateSpatial.snakeSaoCtuData.reset(rectangle, 0, 0, 0);

            rectangle.width <<= h[CtbLog2SizeY()];
            rectangle.height <<= h[CtbLog2SizeY()];

            Neighbourhood *neighbourhood = h;
            neighbourhood->snakeMerge.reset(rectangle, h[MinCbLog2SizeY()] - 1, 1, 1);
            neighbourhood->snake.reset(rectangle, h[MinCbLog2SizeY()] - 1, 1, 1);
        }
    }

    const int rx = (h[CtbAddrInRs()] % h[PicWidthInCtbsY()]);
    const int ry = (h[CtbAddrInRs()] / h[PicWidthInCtbsY()]);

    StateSpatial *stateSpatial = h;
    SaoCtuData saoCtuData = SaoCtuData();
    stateSpatial->snakeSaoCtuData.commitRectangle(Turing::Rectangle{ rx, ry, 1, 1 }, saoCtuData, 0);
}


template <class H>
void postCtu(H &h)
{
    if (h[entropy_coding_sync_enabled_flag()] == 1 && h[CtbAddrInRs()] % h[PicWidthInCtbsY()] == 1)
    {
        h(ContextsSave(tablesWpp));
    }
}

#endif
