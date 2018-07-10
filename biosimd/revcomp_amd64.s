// Copyright 2018 GRAIL, Inc.  All rights reserved.
// Use of this source code is governed by the Apache-2.0
// license that can be found in the LICENSE file.

// +build amd64,!appengine

        DATA ·Mask0f0f<>+0x00(SB)/8, $0x0f0f0f0f0f0f0f0f
        DATA ·Mask0f0f<>+0x08(SB)/8, $0x0f0f0f0f0f0f0f0f
        GLOBL ·Mask0f0f<>(SB), 24, $16
        // NOPTR = 16, RODATA = 8

        DATA ·Reverse8<>+0x00(SB)/8, $0x08090a0b0c0d0e0f
        DATA ·Reverse8<>+0x08(SB)/8, $0x0001020304050607
        GLOBL ·Reverse8<>(SB), 24, $16
        DATA ·Reverse8Minus16<>+0x00(SB)/8, $0xf8f9fafbfcfdfeff
        DATA ·Reverse8Minus16<>+0x08(SB)/8, $0xf0f1f2f3f4f5f6f7
        GLOBL ·Reverse8Minus16<>(SB), 24, $16

        DATA ·ReverseComp8NoTTable<>+0x00(SB)/8, $0x434e4e4e474e544e
        DATA ·ReverseComp8NoTTable<>+0x08(SB)/8, $0x4e4e4e4e4e4e4e4e
        GLOBL ·ReverseComp8NoTTable<>(SB), 24, $16
        DATA ·Capitalizer<>+0x00(SB)/8, $0xdfdfdfdfdfdfdfdf
        DATA ·Capitalizer<>+0x08(SB)/8, $0xdfdfdfdfdfdfdfdf
        GLOBL ·Capitalizer<>(SB), 24, $16
        DATA ·All64<>+0x00(SB)/8, $0x4040404040404040
        DATA ·All64<>+0x08(SB)/8, $0x4040404040404040
        GLOBL ·All64<>(SB), 24, $16
        DATA ·AllT<>+0x00(SB)/8, $0x5454545454545454
        DATA ·AllT<>+0x08(SB)/8, $0x5454545454545454
        GLOBL ·AllT<>(SB), 24, $16

TEXT ·reverseCompInplaceTinyLookupSSSE3Asm(SB),4,$0-24
        // Critical to avoid single-byte-at-a-time table lookup whenever
        // possible.
        // (Could delete this function and force caller to use the non-inplace
        // version; only difference is one extra stack push/pop.)
        // (todo: benchmark this against base/simd reverse functions)
        MOVQ    seq8+0(FP), SI
        MOVQ    tablePtr+8(FP), BX
        MOVD    nByte+16(FP), X2

        MOVOU   ·Reverse8Minus16<>(SB), X0
        MOVOU   (BX), X1
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

TEXT ·reverseCompInplaceLookupSSSE3Asm(SB),4,$0-24
        // This is only called with nByte > 16.  So we can safely divide this
        // into two cases:
        // 1. (nByte+15) % 32 in {0..15}.  Execute (nByte+15)/32 normal
        //    iterations and exit.  Last two writes usually overlap.
        // 2. (nByte+15) % 32 in {16..31}.  Execute (nByte-17)/32 normal
        //    iterations.  Then we have between 33 and 48 central bytes left;
        //    handle them by processing *three* vectors at once at the end.
        MOVQ    main+0(FP), SI
        MOVQ    tablePtr+8(FP), BX
        MOVQ    nByte+16(FP), AX

        // DI iterates backwards from the end of main[].
        LEAQ    -16(SI)(AX*1), DI

        MOVOU   ·Reverse8<>(SB), X0
        MOVOU   (BX), X1
        SUBQ    $1, AX
        SHRQ    $1, AX
        MOVQ    AX, BX
        ANDQ    $8, BX
        // BX is now 0 when we don't need to process 3 vectors at the end, and
        // 8 when we do.
        LEAQ    0(AX)(BX*2), R9
        // R9 is now (nByte+31)/2 when we don't need to process 3 vectors at
        // the end, and (nByte-1)/2 when we do.
        LEAQ    -24(SI)(R9*1), AX
        // AX can now be used for the loop termination check:
        //   if nByte == 17, R9 == 24, so AX == &(seq8[0]).
        //   if nByte == 32, R9 == 31, so AX == &(seq8[7]).
        //   if nByte == 33, R9 == 16, so AX == &(seq8[-8]).
        //   if nByte == 48, R9 == 23, so AX == &(seq8[-1]).
        CMPQ    AX, SI
        JL      reverseCompInplaceLookupSSSE3LastThree

reverseCompInplaceLookupSSSE3Loop:
        MOVOU   (SI), X2
        MOVOU   (DI), X3
        PSHUFB  X0, X2
        PSHUFB  X0, X3
        MOVO    X1, X4
        MOVO    X1, X5
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        MOVOU   X5, (SI)
        MOVOU   X4, (DI)
        ADDQ    $16, SI
        SUBQ    $16, DI
        CMPQ    AX, SI
        JGE     reverseCompInplaceLookupSSSE3Loop

        TESTQ   BX, BX
        JNE     reverseCompInplaceLookupSSSE3Ret
reverseCompInplaceLookupSSSE3LastThree:
        MOVOU   (SI), X2
        MOVOU   16(SI), X3
        MOVOU   (DI), X4
        PSHUFB  X0, X2
        PSHUFB  X0, X3
        PSHUFB  X0, X4
        MOVO    X1, X5
        MOVO    X1, X6
        PSHUFB  X4, X1
        PSHUFB  X2, X5
        PSHUFB  X3, X6
        MOVOU   X1, (SI)
        MOVOU   X6, -16(DI)
        MOVOU   X5, (DI)

reverseCompInplaceLookupSSSE3Ret:
        RET

TEXT ·reverseCompTinyLookupSSSE3Asm(SB),4,$0-32
        MOVQ    dst+0(FP), DI
        MOVQ    src+8(FP), SI
        MOVQ    tablePtr+16(FP), BX
        MOVD    nByte+24(FP), X2

        MOVOU   ·Reverse8Minus16<>(SB), X0
        MOVOU   (BX), X1
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

TEXT ·reverseCompLookupSSSE3Asm(SB),4,$0-32
        // This is only called with nByte >= 16.  Fortunately, this doesn't
        // have the same complications re: potentially clobbering data we need
        // to keep that the in-place function must deal with.
        MOVQ    dst+0(FP), DI
        MOVQ    src+8(FP), BX
        MOVQ    tablePtr+16(FP), R8
        MOVQ    nByte+24(FP), AX

        // SI iterates backwards from the end of src[].
        LEAQ    -16(BX)(AX*1), SI
        // May as well save start of final dst[] vector.
        LEAQ    -16(DI)(AX*1), R9

        MOVOU   ·Reverse8<>(SB), X0
        MOVOU   (R8), X1

reverseCompLookupSSSE3Loop:
        MOVOU   (SI), X2
        PSHUFB  X0, X2
        MOVO    X1, X3
        PSHUFB  X2, X3
        MOVOU   X3, (DI)
        SUBQ    $16, SI
        ADDQ    $16, DI
        CMPQ    R9, DI
        JG      reverseCompLookupSSSE3Loop

        MOVOU   (BX), X2
        PSHUFB  X0, X2
        PSHUFB  X2, X1
        MOVOU   X1, (R9)
        RET

TEXT ·reverseComp8InplaceSSSE3Asm(SB),4,$0-16
        // (Iteration logic is identical to reverseCompInplaceLookupSSSE3Asm.)
        // This is only called with nByte > 16.  So we can safely divide this
        // into two cases:
        // 1. (nByte+15) % 32 in {0..15}.  Execute (nByte+15)/32 normal
        //    iterations and exit.  Last two writes usually overlap.
        // 2. (nByte+15) % 32 in {16..31}.  Execute (nByte-17)/32 normal
        //    iterations.  Then we have between 33 and 48 central bytes left;
        //    handle them by processing *three* vectors at once at the end.
        MOVQ    main+0(FP), SI
        MOVQ    nByte+8(FP), AX

        // DI iterates backwards from the end of main[].
        LEAQ    -16(SI)(AX*1), DI

        MOVOU   ·Reverse8<>(SB), X0
        // Can't perform nibble-lookup since we need to map e.g. 'S' to 'N'.
        // Instead, we do the following:
        // 1. Capitalize all letters by masking with 0xdf.  (May also affect
        //    non-letter bytes, but only in inconsequential ways.)  Removing
        //    this step would speed things up by a tiny bit, but if you
        //    really want to squeeze out every last bit of speed, you almost
        //    certainly shouldn't be computing on the ASCII representation
        //    anyway; this function is about drop-in convenience.
        // 2. Parallel-equality check for 'T'.
        // 3. Check for high 4 bits == 0x04.
        // 4. Apply mask, perform non-T reverse-complement.
        // 5. Merge in 'T' results and save.
        MOVOU   ·ReverseComp8NoTTable<>(SB), X1
        MOVOU   ·Capitalizer<>(SB), X7
        MOVOU   ·All64<>(SB), X8
        // 'N' ^ 'A' happens to be 0x0f, which we already need for other
        // reasons!
        MOVOU   ·Mask0f0f<>(SB), X9
        MOVOU   ·AllT<>(SB), X10
        SUBQ    $1, AX
        SHRQ    $1, AX
        MOVQ    AX, BX
        ANDQ    $8, BX
        // BX is now 0 when we don't need to process 3 vectors at the end, and
        // 8 when we do.
        LEAQ    0(AX)(BX*2), R9
        // R9 is now (nByte+31)/2 when we don't need to process 3 vectors at
        // the end, and (nByte-1)/2 when we do.
        LEAQ    -24(SI)(R9*1), AX
        // AX can now be used for the loop termination check:
        //   if nByte == 17, R9 == 24, so AX == &(seq8[0]).
        //   if nByte == 32, R9 == 31, so AX == &(seq8[7]).
        //   if nByte == 33, R9 == 16, so AX == &(seq8[-8]).
        //   if nByte == 48, R9 == 23, so AX == &(seq8[-1]).
        CMPQ    AX, SI
        JL      reverseComp8InplaceSSSE3LastThree

reverseComp8InplaceSSSE3Loop:
        MOVOU   (SI), X2
        MOVOU   (DI), X3
        // Reverse.
        PSHUFB  X0, X2
        PSHUFB  X0, X3
        // Capitalize.
        PAND    X7, X2
        PAND    X7, X3
        // Check for 'T's.
        MOVO    X10, X11
        MOVO    X10, X12
        PCMPEQB X2, X11
        PCMPEQB X3, X12
        PAND    X9, X11
        PAND    X9, X12
        // Check for high bits == 0x40.
        MOVO    X9, X13
        MOVO    X9, X14
        PANDN   X2, X13
        PANDN   X3, X14
        PCMPEQB X8, X13
        PCMPEQB X8, X14
        // X13/X14 now describe which non-T bases are valid; zero out the rest.
        PAND    X13, X2
        PAND    X14, X3
        MOVO    X1, X4
        MOVO    X1, X5
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        PXOR    X11, X4
        PXOR    X12, X5
        MOVOU   X5, (SI)
        MOVOU   X4, (DI)
        ADDQ    $16, SI
        SUBQ    $16, DI
        CMPQ    AX, SI
        JGE     reverseComp8InplaceSSSE3Loop

        TESTQ   BX, BX
        JNE     reverseComp8InplaceSSSE3Ret
reverseComp8InplaceSSSE3LastThree:
        MOVOU   (SI), X2
        MOVOU   16(SI), X3
        MOVOU   (DI), X4
        PSHUFB  X0, X2
        PSHUFB  X0, X3
        PSHUFB  X0, X4
        PAND    X7, X2
        PAND    X7, X3
        PAND    X7, X4

        MOVO    X10, X11
        MOVO    X10, X12
        // Don't need X10 any more; use it as the third 'T' tracker.
        PCMPEQB X2, X11
        PCMPEQB X3, X12
        PCMPEQB X4, X10
        PAND    X9, X11
        PAND    X9, X12
        PAND    X9, X10

        MOVO    X9, X13
        MOVO    X9, X14
        PANDN   X2, X13
        PANDN   X3, X14
        PANDN   X4, X9
        PCMPEQB X8, X13
        PCMPEQB X8, X14
        PCMPEQB X8, X9
        // X13/X14/X9 now describe which non-T bases are valid.
        PAND    X13, X2
        PAND    X14, X3
        PAND    X9, X4
        MOVO    X1, X5
        MOVO    X1, X6
        PSHUFB  X4, X1
        PSHUFB  X2, X5
        PSHUFB  X3, X6
        PXOR    X10, X1
        PXOR    X11, X5
        PXOR    X12, X6
        MOVOU   X1, (SI)
        MOVOU   X6, -16(DI)
        MOVOU   X5, (DI)

reverseComp8InplaceSSSE3Ret:
        RET
