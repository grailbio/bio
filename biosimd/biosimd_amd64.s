// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

        DATA ·Mask0f0f<>+0x00(SB)/8, $0x0f0f0f0f0f0f0f0f
        DATA ·Mask0f0f<>+0x08(SB)/8, $0x0f0f0f0f0f0f0f0f
        GLOBL ·Mask0f0f<>(SB), 24, $16
        // NOPTR = 16, RODATA = 8

        DATA ·GatherOddLow<>+0x00(SB)/8, $0x0f0d0b0907050301
        DATA ·GatherOddLow<>+0x08(SB)/8, $0xffffffffffffffff
        GLOBL ·GatherOddLow<>(SB), 24, $16
        DATA ·GatherOddHigh<>+0x00(SB)/8, $0xffffffffffffffff
        DATA ·GatherOddHigh<>+0x08(SB)/8, $0x0f0d0b0907050301
        GLOBL ·GatherOddHigh<>(SB), 24, $16

        DATA ·Reverse16<>+0x00(SB)/8, $0x08090a0b0c0d0e0f
        DATA ·Reverse16<>+0x08(SB)/8, $0x0001020304050607
        GLOBL ·Reverse16<>(SB), 24, $16
        DATA ·Reverse16Minus16<>+0x00(SB)/8, $0xf8f9fafbfcfdfeff
        DATA ·Reverse16Minus16<>+0x08(SB)/8, $0xf0f1f2f3f4f5f6f7
        GLOBL ·Reverse16Minus16<>(SB), 24, $16
        DATA ·ReverseCompLookup<>+0x00(SB)/8, $0x0e060a020c040800
        DATA ·ReverseCompLookup<>+0x08(SB)/8, $0x0f070b030d050901
        GLOBL ·ReverseCompLookup<>(SB), 24, $16

// This was forked from github.com/willf/bitset .
// Some form of AVX2/AVX-512 detection will probably be added later.
TEXT ·hasSSE42Asm(SB),4,$0-1
        MOVQ    $1, AX
        CPUID
        SHRQ    $23, CX
        ANDQ    $1, CX
        MOVB    CX, ret+0(FP)
        RET

TEXT ·unpackSeqSSE2Asm(SB),4,$0-24
        // Based on packedNibbleLookupSSSE3Asm() in base/simd/simd_amd64.s.
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nDstVec+16(FP), CX

        MOVOU   ·Mask0f0f<>(SB), X0

        // AX = pointer to last relevant word of src[].
        // (note that 8 src bytes -> 16 dst bytes)
        LEAQ    -8(DI)(CX*8), AX
        CMPQ    AX, DI
        JLE     unpackSeqSSE2Final

unpackSeqSSE2Loop:
        MOVOU   (DI), X1
        // Isolate high and low nibbles.
        MOVOU   X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        // Use unpacklo/unpackhi to stitch results together.
        // Even bytes (0, 2, 4, ...) are in X1/X3, odd in X2.
        MOVOU   X1, X3
        PUNPCKLBW       X2, X1
        PUNPCKHBW       X2, X3
        MOVOU   X1, (R8)
        MOVOU   X3, 16(R8)
        ADDQ    $16, DI
        ADDQ    $32, R8
        CMPQ    AX, DI
        JG      unpackSeqSSE2Loop
unpackSeqSSE2Final:
        JL      unpackSeqSSE2Finish

        // Necessary to write one more vector.  We skip unpackhi, but must
        // execute the rest of the loop body.
        MOVOU   (DI), X1
        // Isolate high and low nibbles.
        MOVOU   X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        PUNPCKLBW       X2, X1
        MOVOU   X1, (R8)

unpackSeqSSE2Finish:
        RET

TEXT ·unpackSeqOddSSE2Asm(SB),4,$0-24
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nSrcFullByte+16(FP), CX

        MOVOU   ·Mask0f0f<>(SB), X0

        // set AX to 32 bytes before end of dst[].
        // change CX to 16 bytes before end of src[].
        SUBQ    $16, CX
        LEAQ    0(R8)(CX*2), AX
        ADDQ    DI, CX

unpackSeqOddSSE2Loop:
        MOVOU   (DI), X1
        // Isolate high and low nibbles, then parallel-lookup.
        MOVOU   X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        // Use unpacklo/unpackhi to stitch results together.
        // Even bytes (0, 2, 4, ...) are in X1/X3, odd in X2.
        MOVOU   X1, X3
        PUNPCKLBW       X2, X1
        PUNPCKHBW       X2, X3
        MOVOU   X1, (R8)
        MOVOU   X3, 16(R8)
        ADDQ    $16, DI
        ADDQ    $32, R8
        CMPQ    CX, DI
        JG      unpackSeqOddSSE2Loop

        // Final usually-unaligned read and write.
        MOVOU   (CX), X1
        MOVOU   X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        MOVOU   X1, X3
        PUNPCKLBW       X2, X1
        PUNPCKHBW       X2, X3
        MOVOU   X1, (AX)
        MOVOU   X3, 16(AX)

unpackSeqOddSSE2Finish:
        RET

TEXT ·packSeqSSE41Asm(SB),4,$0-24
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nSrcVec+16(FP), CX

        MOVOU   ·GatherOddLow<>(SB), X0
        MOVOU   ·GatherOddHigh<>(SB), X1

        // AX = pointer to last relevant word of dst[].
        // (note that 16 src bytes -> 8 dst bytes)
        LEAQ    -8(R8)(CX*8), AX
        CMPQ    AX, R8
        JLE     packSeqSSE41Final

packSeqSSE41Loop:
        MOVOU   (DI), X2
        MOVOU   16(DI), X3
        MOVOU   X2, X4
        MOVOU   X3, X5
        PSLLQ   $12, X2
        PSLLQ   $12, X3
        POR     X4, X2
        POR     X5, X3
        // If all bytes of src[] were <16, the odd positions of X2/X3 now
        // contain the values of interest.  Gather them.
        PSHUFB  X0, X2
        PSHUFB  X1, X3
        POR     X3, X2
        MOVOU   X2, (R8)
        ADDQ    $32, DI
        ADDQ    $16, R8
        CMPQ    AX, R8
        JG      packSeqSSE41Loop
packSeqSSE41Final:
        JL      packSeqSSE41Finish

        // Necessary to write one more word.
        MOVOU   (DI), X2
        MOVOU   X2, X4
        PSLLQ   $12, X2
        POR     X4, X2
        PSHUFB  X0, X2
        PEXTRQ  $0, X2, (R8)

packSeqSSE41Finish:
        RET

TEXT ·packSeqOddSSSE3Asm(SB),4,$0-24
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nDstFullByte+16(FP), CX

        MOVOU   ·GatherOddLow<>(SB), X0
        MOVOU   ·GatherOddHigh<>(SB), X1

        // Set AX to 32 bytes before end of src[], and change CX to 16 bytes
        // before end of dst[].
        SUBQ    $16, CX
        LEAQ    0(DI)(CX*2), AX
        ADDQ    R8, CX

packSeqOddSSSE3Loop:
        MOVOU   (DI), X2
        MOVOU   16(DI), X3
        MOVOU   X2, X4
        MOVOU   X3, X5
        PSLLQ   $12, X2
        PSLLQ   $12, X3
        POR     X4, X2
        POR     X5, X3
        // If all bytes of src[] were <16, the odd positions of X2/X3 now
        // contain the values of interest.  Gather them.
        PSHUFB  X0, X2
        PSHUFB  X1, X3
        POR     X3, X2
        MOVOU   X2, (R8)
        ADDQ    $32, DI
        ADDQ    $16, R8
        CMPQ    AX, DI
        JG      packSeqOddSSSE3Loop

        // Final usually-unaligned read and write.
        MOVOU   (AX), X2
        MOVOU   16(AX), X3
        MOVOU   X2, X4
        MOVOU   X3, X5
        PSLLQ   $12, X2
        PSLLQ   $12, X3
        POR     X4, X2
        POR     X5, X3
        PSHUFB  X0, X2
        PSHUFB  X1, X3
        POR     X3, X2
        MOVOU   X2, (CX)
        RET

TEXT ·reverseCompInplaceTinySSSE3Asm(SB),4,$0-16
        // Critical to avoid single-byte-at-a-time table lookup whenever
        // possible.
        // (Could delete this function and force caller to use the non-inplace
        // version; only difference is one extra stack push/pop.)
        // (todo: benchmark this against base/simd reverse functions)
        MOVQ    seq8+0(FP), SI
        MOVD    nByte+8(FP), X2

        MOVOU   ·Reverse16Minus16<>(SB), X0
        MOVOU   ·ReverseCompLookup<>(SB), X1
        PXOR    X3, X3
        PSHUFB  X3, X2
        // all bytes of X2 are now equal to nByte
        PADDB   X0, X2
        // now X2 is {nByte-1, nByte-2, ...}

        MOVOU   (SI), X3
        PSHUFB  X2, X3
        PSHUFB  X3, X1
        MOVOU   X1, (SI)
        RET

TEXT ·reverseCompInplaceSSSE3Asm(SB),4,$0-16
        // This is only called with nByte > 16.  So we can safely divide this
        // into two cases:
        // 1. (nByte+15) % 32 in {0..15}.  Execute (nByte+15)/32 normal
        //    iterations and exit.  Last two writes usually overlap.
        // 2. (nByte+15) % 32 in {16..31}.  Execute (nByte-17)/32 normal
        //    iterations.  Then we have between 33 and 48 central bytes left;
        //    handle them by processing *three* vectors at once at the end.
        MOVQ    seq8+0(FP), SI
        MOVQ    nByte+8(FP), AX

        // DI iterates backwards from the end of seq8[].
        LEAQ    -16(SI)(AX*1), DI

        MOVOU   ·Reverse16<>(SB), X0
        MOVOU   ·ReverseCompLookup<>(SB), X1
        SUBQ    $1, AX
        SHRQ    $1, AX
        MOVQ    AX, BX
        ANDQ    $8, BX
        // BX is now 0 when we don't need to process 3 vectors at the end, and
        // 8 when we do.
        LEAQ    0(AX)(BX*2), CX
        // CX is now (nByte+31)/2 when we don't need to process 3 vectors at
        // the end, and (nByte-1)/2 when we do.
        LEAQ    -24(SI)(CX*1), AX
        // AX can now be used for the loop termination check:
        //   if nByte == 17, CX == 24, so AX == &(seq8[0]).
        //   if nByte == 32, CX == 31, so AX == &(seq8[7]).
        //   if nByte == 33, CX == 16, so AX == &(seq8[-8]).
        //   if nByte == 48, CX == 23, so AX == &(seq8[-1]).
        CMPQ    AX, SI
        JL      reverseCompInplaceSSSE3LastThree

reverseCompInplaceSSSE3Loop:
        MOVOU   (SI), X2
        MOVOU   (DI), X3
        PSHUFB  X0, X2
        PSHUFB  X0, X3
        MOVOU   X1, X4
        MOVOU   X1, X5
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        MOVOU   X5, (SI)
        MOVOU   X4, (DI)
        ADDQ    $16, SI
        SUBQ    $16, DI
        CMPQ    AX, SI
        JGE     reverseCompInplaceSSSE3Loop

        TESTQ   BX, BX
        JNE     reverseCompInplaceSSSE3Ret
reverseCompInplaceSSSE3LastThree:
        MOVOU   (SI), X2
        MOVOU   16(SI), X3
        MOVOU   (DI), X4
        PSHUFB  X0, X2
        PSHUFB  X0, X3
        PSHUFB  X0, X4
        MOVOU   X1, X5
        MOVOU   X1, X6
        MOVOU   X1, X7
        PSHUFB  X2, X5
        PSHUFB  X3, X6
        PSHUFB  X4, X7
        MOVOU   X7, (SI)
        MOVOU   X6, -16(DI)
        MOVOU   X5, (DI)

reverseCompInplaceSSSE3Ret:
        RET

TEXT ·reverseCompTinySSSE3Asm(SB),4,$0-24
        MOVQ    dst+0(FP), DI
        MOVQ    src+8(FP), SI
        MOVD    nByte+16(FP), X2

        MOVOU   ·Reverse16Minus16<>(SB), X0
        MOVOU   ·ReverseCompLookup<>(SB), X1
        PXOR    X3, X3
        PSHUFB  X3, X2
        // all bytes of X2 are now equal to nByte
        PADDB   X0, X2
        // now X2 is {nByte-1, nByte-2, ...}

        MOVOU   (SI), X3
        PSHUFB  X2, X3
        PSHUFB  X3, X1
        MOVOU   X1, (DI)
        RET

TEXT ·reverseCompSSSE3Asm(SB),4,$0-24
        // This is only called with nByte >= 16.  Fortunately, this doesn't
        // have the same complications re: potentially clobbering data we need
        // to keep that the in-place function must deal with.
        MOVQ    dst+0(FP), DI
        MOVQ    src+8(FP), BX
        MOVQ    nByte+16(FP), AX

        // SI iterates backwards from the end of src[].
        LEAQ    -16(BX)(AX*1), SI
        // May as well save start of final dst[] vector.
        LEAQ    -16(DI)(AX*1), CX

        MOVOU   ·Reverse16<>(SB), X0
        MOVOU   ·ReverseCompLookup<>(SB), X1

reverseCompSSSE3Loop:
        MOVOU   (SI), X2
        PSHUFB  X0, X2
        MOVOU   X1, X3
        PSHUFB  X2, X3
        MOVOU   X3, (DI)
        SUBQ    $16, SI
        ADDQ    $16, DI
        CMPQ    CX, DI
        JG      reverseCompSSSE3Loop

        MOVOU   (BX), X2
        PSHUFB  X0, X2
        MOVOU   X1, X3
        PSHUFB  X2, X3
        MOVOU   X3, (CX)
        RET
