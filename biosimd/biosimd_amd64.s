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

// This was forked from github.com/willf/bitset .
// Some form of AVX2/AVX-512 detection will probably be added later.
TEXT ·hasSSE42Asm(SB),4,$0-1
        MOVQ    $1, AX
        CPUID
        SHRQ    $23, CX
        ANDQ    $1, CX
        MOVB    CX, ret+0(FP)
        RET

TEXT ·unpackSeqSSE2Asm(SB),4,$0-32
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

TEXT ·unpackSeqOddSSE2Asm(SB),4,$0-32
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

TEXT ·packSeqSSE41Asm(SB),4,$0-32
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

TEXT ·packSeqOddSSSE3Asm(SB),4,$0-32
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
