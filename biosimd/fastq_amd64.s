// Copyright 2019 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

        DATA ·Mask0f0f<>+0x00(SB)/8, $0x0f0f0f0f0f0f0f0f
        DATA ·Mask0f0f<>+0x08(SB)/8, $0x0f0f0f0f0f0f0f0f
        GLOBL ·Mask0f0f<>(SB), 24, $16
        // NOPTR = 16, RODATA = 8

TEXT ·fillFastqRecordBodyFromNibblesSSSE3Asm(SB),4,$0-40
        // Requires nBase >= 32.
        // We perform the following three operations in sequence:
        // 1. Render all but the last 0-32 bytes of both the base and qual
        //    lines, 32 bytes at a time, with the usual PackedNibbleLookup
        //    strategy over a full input-vector at a time.  We gain a small
        //    speed boost over the usual pair of PackedNibbleLookup calls since
        //    we only have to read each input-vector once.
        // 2. Render last 31-32 bytes of both the base and qual lines, usually
        //    overlapping the previous pair of writes, AND POSSIBLY WRITING 1
        //    CHARACTER PAST THE END on both lines.  By doing the latter, we
        //    remove the need for a "if srcOdd == 1" code branch, and it
        //    doesn't actually force the function to have any more
        //    preconditions, because...
        // 3. Write line-breaks for base and qual lines, and "+" second line.
        //    Since dst is validated to have enough space to fit these, the
        //    potential 1-byte overwrite in step 2 is guaranteed to be
        //    harmless.

        // The first two steps mirror packedNibbleLookupOddSSSE3Asm in
        // base/simd/simd_amd64.s.

        // In the main loop:
        //   AX = dst[] (bases) write pointer
        // * BX = quals[] write pointer
        //   CX = src[] read pointer
        // * DX = position of final (usually unaligned) src load
        // * SI = parity of nBase
        //   DI = end of bases[] (where the linebreak will be added)
        //   R8 = end of quals[]
        // Registers marked with '*' are first used as temporaries during
        // initialization.
        // (Previously, I avoided CX since there used to be some weird build
        // setting under which it couldn't be used, and it's not like I ever
        // needed all 14 regular registers.  But I no longer see any mention of
        // such a condition in the documentation, and in the meantime I've
        // learned that instruction encodings can be shorter for
        // AX/BX/CX/DX/SI/DI than for R8..R15.  So I'm now structuring my
        // assembly code so that more of the action involves the first six
        // registers.
        // Incidentally, there are a few instructions, like variable-shift,
        // which *require* the CX register as an argument; a bit of extra
        // planning is required when they're present.  But we aren't using
        // those as of this writing.)
        MOVQ    dst+0(FP), AX
        MOVQ    nBase+32(FP), SI
        MOVQ    src+8(FP), CX
        LEAQ    0(AX)(SI*1), DI
        LEAQ    3(DI)(SI*1), R8

        // X0 = baseTable
        // X1 = qualTable
        // X2 = Mask0f0f
        MOVQ    baseTablePtr+16(FP), BX
        MOVQ    qualTablePtr+24(FP), DX
        MOVOU   (BX), X0
        MOVOU   (DX), X1
        MOVOU   ·Mask0f0f<>(SB), X2

        // Final src load is from offset (nBase - 31) >> 1.
        LEAQ    -31(SI), BX
        ANDQ    $1, SI
        SHRQ    $1, BX
        LEAQ    0(CX)(BX*1), DX

        LEAQ    3(DI), BX

fillFastqRecordBodyFromNibblesSSSE3Loop:
        // Load next 16 bytes from src.
        MOVOU   (CX), X4
        // Split into high and low nibbles.
        MOVO    X4, X3
        PSRLQ   $4, X4
        PAND    X2, X3
        PAND    X2, X4
        // X3 now stores low nibbles of src, and X4 the high nibbles.

        // Parallel-lookup of baseTable.
        MOVO    X0, X5
        MOVO    X0, X6
        PSHUFB  X3, X5
        PSHUFB  X4, X6
        // Stitch results together and flush to dst[].
        MOVO    X5, X7
        PUNPCKLBW       X6, X5
        PUNPCKHBW       X6, X7
        MOVOU   X5, (AX)
        MOVOU   X7, 16(AX)
        // Parallel-lookup of qualTable.
        MOVO    X1, X5
        MOVO    X1, X6
        PSHUFB  X3, X5
        PSHUFB  X4, X6
        // Stitch results together and flush to quals[].
        MOVO    X5, X7
        PUNPCKLBW       X6, X5
        PUNPCKHBW       X6, X7
        MOVOU   X5, (BX)
        MOVOU   X7, 16(BX)

        // Advance to next vector.
        ADDQ    $16, CX
        ADDQ    $32, AX
        ADDQ    $32, BX
        CMPQ    DX, CX
        JG      fillFastqRecordBodyFromNibblesSSSE3Loop

        // Final usually-unaligned base/qual read and write.
        MOVOU   (DX), X4

        // Same inner-loop logic as above, but we're now free to clobber
        // X0/X1/X2 the last time we reference them.  (We could simplify this
        // function by reusing the loop above, but unfortunately this wipes out
        // ~30% of the speed improvement on the single-threaded benchmark in my
        // testing.)
        MOVO    X4, X3
        PSRLQ   $4, X4
        PAND    X2, X3
        PAND    X2, X4
        MOVO    X0, X5
        PSHUFB  X3, X5
        PSHUFB  X4, X0

        MOVO    X5, X7
        PUNPCKLBW       X0, X5
        PUNPCKHBW       X0, X7

        // Final bases[] write position is <end of bases[]> - 32 + <parity>.
        // = DI - 32 + SI
        // This potentially writes a character past the end, but that's okay
        // because we unconditionally overwrite it with a linebreak.
        LEAQ    -32(DI)(SI*1), AX
        MOVOU   X5, (AX)
        MOVOU   X7, 16(AX)

        MOVO    X1, X5
        PSHUFB  X3, X5
        PSHUFB  X4, X1
        MOVO    X5, X7
        PUNPCKLBW       X1, X5
        PUNPCKHBW       X1, X7

        LEAQ    -32(R8)(SI*1), BX
        MOVOU   X5, (BX)
        MOVOU   X7, 16(BX)

        // Write "\n+\n" to end of bases[].
        MOVW    $0x2b0a, (DI)
        ADDQ    $2, DI
        MOVB    $0x0a, (DI)

        // Write \n to end of quals[].
        MOVB    $0x0a, (R8)
        RET
