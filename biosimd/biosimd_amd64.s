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

        DATA ·CleanNoTTable<>+0x00(SB)/8, $0x474e4e4e434e414e
        DATA ·CleanNoTTable<>+0x08(SB)/8, $0x4e4e4e4e4e4e4e4e
        GLOBL ·CleanNoTTable<>(SB), 24, $16
        DATA ·Capitalizer<>+0x00(SB)/8, $0xdfdfdfdfdfdfdfdf
        DATA ·Capitalizer<>+0x08(SB)/8, $0xdfdfdfdfdfdfdfdf
        GLOBL ·Capitalizer<>(SB), 24, $16
        DATA ·All64<>+0x00(SB)/8, $0x4040404040404040
        DATA ·All64<>+0x08(SB)/8, $0x4040404040404040
        GLOBL ·All64<>(SB), 24, $16
        DATA ·AllT<>+0x00(SB)/8, $0x5454545454545454
        DATA ·AllT<>+0x08(SB)/8, $0x5454545454545454
        GLOBL ·AllT<>(SB), 24, $16
        DATA ·AllNXorT<>+0x00(SB)/8, $0x1a1a1a1a1a1a1a1a
        DATA ·AllNXorT<>+0x08(SB)/8, $0x1a1a1a1a1a1a1a1a
        GLOBL ·AllNXorT<>(SB), 24, $16

        DATA ·ASCIIToSeq8NoTTable<>+0x00(SB)/8, $0x040f0f0f020f010f
        DATA ·ASCIIToSeq8NoTTable<>+0x08(SB)/8, $0x0f0f0f0f0f0f0f0f
        GLOBL ·ASCIIToSeq8NoTTable<>(SB), 24, $16
        DATA ·ASCIIToSeq8NXorT<>+0x00(SB)/8, $0x0707070707070707
        DATA ·ASCIIToSeq8NXorT<>+0x08(SB)/8, $0x0707070707070707
        GLOBL ·ASCIIToSeq8NXorT<>(SB), 24, $16

        DATA ·Maskd0d0<>+0x00(SB)/8, $0xd0d0d0d0d0d0d0d0
        DATA ·Maskd0d0<>+0x08(SB)/8, $0xd0d0d0d0d0d0d0d0
        GLOBL ·Maskd0d0<>(SB), 24, $16
        DATA ·MaskNoTTable<>+0x00(SB)/8, $0xff000000ff00ff00
        DATA ·MaskNoTTable<>+0x08(SB)/8, $0x0000000000000000
        GLOBL ·MaskNoTTable<>(SB), 24, $16
        DATA ·AllN<>+0x00(SB)/8, $0x4e4e4e4e4e4e4e4e
        DATA ·AllN<>+0x08(SB)/8, $0x4e4e4e4e4e4e4e4e
        GLOBL ·AllN<>(SB), 24, $16


// This was forked from github.com/willf/bitset .
// Some form of AVX2/AVX-512 detection will probably be added later.
TEXT ·hasSSE42Asm(SB),4,$0-1
        MOVQ    $1, AX
        CPUID
        // CPUID function fills CX register.
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
        MOVQ	nSrcByte+16(FP), R9

        MOVOU   ·Mask0f0f<>(SB), X0

        // AX = pointer to last relevant word of src[].
        // (note that 8 src bytes -> 16 dst bytes)
        LEAQ    -8(DI)(R9*1), AX
        CMPQ    AX, DI
        JLE     unpackSeqSSE2Final

unpackSeqSSE2Loop:
        MOVOU   (DI), X1
        // Isolate high and low nibbles.
        MOVO    X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        // Use unpacklo/unpackhi to stitch results together.
        // Even bytes (0, 2, 4, ...) are in X1/X3, odd in X2.
        MOVO    X1, X3
        PUNPCKLBW       X2, X1
        PUNPCKHBW       X2, X3
        MOVOU   X1, (R8)
        MOVOU   X3, 16(R8)
        ADDQ    $16, DI
        ADDQ    $32, R8
        CMPQ    AX, DI
        JG      unpackSeqSSE2Loop
unpackSeqSSE2Final:
        // Necessary to write one more vector.  We skip unpackhi, but must
        // execute the rest of the loop body.
        MOVOU   (DI), X1
        // Isolate high and low nibbles.
        MOVO    X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        PUNPCKLBW       X2, X1
        MOVOU   X1, (R8)
        RET

TEXT ·unpackSeqOddSSE2Asm(SB),4,$0-24
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nSrcFullByte+16(FP), R9

        MOVOU   ·Mask0f0f<>(SB), X0

        // set AX to 32 bytes before end of dst[].
        // change R9 to 16 bytes before end of src[].
        SUBQ    $16, R9
        LEAQ    0(R8)(R9*2), AX
        ADDQ    DI, R9

unpackSeqOddSSE2Loop:
        MOVOU   (DI), X1
        // Isolate high and low nibbles, then parallel-lookup.
        MOVO    X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        // Use unpacklo/unpackhi to stitch results together.
        // Even bytes (0, 2, 4, ...) are in X1/X3, odd in X2.
        MOVO    X1, X3
        PUNPCKLBW       X2, X1
        PUNPCKHBW       X2, X3
        MOVOU   X1, (R8)
        MOVOU   X3, 16(R8)
        ADDQ    $16, DI
        ADDQ    $32, R8
        CMPQ    R9, DI
        JG      unpackSeqOddSSE2Loop

        // Final usually-unaligned read and write.
        MOVOU   (R9), X1
        MOVO    X1, X2
        PSRLQ   $4, X1
        PAND    X0, X2
        PAND    X0, X1
        MOVO    X1, X3
        PUNPCKLBW       X2, X1
        PUNPCKHBW       X2, X3
        MOVOU   X1, (AX)
        MOVOU   X3, 16(AX)
        RET

TEXT ·packSeqSSE41Asm(SB),4,$0-24
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nSrcByte+16(FP), R9

        MOVOU   ·GatherOddLow<>(SB), X0
        MOVOU   ·GatherOddHigh<>(SB), X1

        // AX = pointer to last relevant word of src[].
        // (note that 16 src bytes -> 8 dst bytes)
        LEAQ    -16(DI)(R9*1), AX
        CMPQ    AX, DI
        JLE     packSeqSSE41Final

packSeqSSE41Loop:
        MOVOU   (DI), X2
        MOVOU   16(DI), X3
        MOVO    X2, X4
        MOVO    X3, X5
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
        JG      packSeqSSE41Loop
packSeqSSE41Final:
        // Necessary to write one more word.
        MOVOU   (DI), X2
        MOVO    X2, X4
        PSLLQ   $12, X2
        POR     X4, X2
        PSHUFB  X0, X2
        PEXTRQ  $0, X2, (R8)
        RET

TEXT ·packSeqOddSSSE3Asm(SB),4,$0-24
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	nDstFullByte+16(FP), R9

        MOVOU   ·GatherOddLow<>(SB), X0
        MOVOU   ·GatherOddHigh<>(SB), X1

        // Set AX to 32 bytes before end of src[], and change R9 to 16 bytes
        // before end of dst[].
        SUBQ    $16, R9
        LEAQ    0(DI)(R9*2), AX
        ADDQ    R8, R9

packSeqOddSSSE3Loop:
        MOVOU   (DI), X2
        MOVOU   16(DI), X3
        MOVO    X2, X4
        MOVO    X3, X5
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
        MOVO    X2, X4
        MOVO    X3, X5
        PSLLQ   $12, X2
        PSLLQ   $12, X3
        POR     X4, X2
        POR     X5, X3
        PSHUFB  X0, X2
        PSHUFB  X1, X3
        POR     X3, X2
        MOVOU   X2, (R9)
        RET

TEXT ·unpackAndReplaceSeqSSSE3Asm(SB),4,$0-32
        // Identical to packedNibbleLookupSSSE3Asm, except with even/odd
        // swapped.
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	tablePtr+16(FP), SI
        MOVQ	nSrcByte+24(FP), R9

        MOVOU   (SI), X0
        MOVOU   ·Mask0f0f<>(SB), X1

        // AX = pointer to last relevant word of src[].
        // (note that 8 src bytes -> 16 dst bytes)
        LEAQ    -8(DI)(R9*1), AX
        CMPQ    AX, DI
        JLE     unpackAndReplaceSeqSSSE3Final

unpackAndReplaceSeqSSSE3Loop:
        MOVOU   (DI), X3
        MOVO    X0, X4
        MOVO    X0, X5
        // Isolate high and low nibbles, then parallel-lookup.
        MOVO    X3, X2
        PSRLQ   $4, X3
        PAND    X1, X2
        PAND    X1, X3
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        // Use unpacklo/unpackhi to stitch results together.
        // Odd bytes (1, 3, 5, ...) are in X4, even in X3/X5.
        MOVO    X5, X3
        PUNPCKLBW       X4, X5
        PUNPCKHBW       X4, X3
        MOVOU   X5, (R8)
        MOVOU   X3, 16(R8)
        ADDQ    $16, DI
        ADDQ    $32, R8
        CMPQ    AX, DI
        JG      unpackAndReplaceSeqSSSE3Loop
unpackAndReplaceSeqSSSE3Final:
        // Necessary to write one more vector.  We skip unpackhi, but must
        // execute the rest of the loop body.
        MOVOU   (DI), X3
        MOVO    X0, X4
        MOVO    X0, X5
        MOVO    X3, X2
        PSRLQ   $4, X3
        PAND    X1, X2
        PAND    X1, X3
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        PUNPCKLBW       X4, X5
        MOVOU   X5, (R8)
        RET

TEXT ·unpackAndReplaceSeqOddSSSE3Asm(SB),4,$0-32
        // Identical to packedNibbleLookupOddSSSE3Asm, except with even/odd
        // swapped.
        // DI = pointer to current src[] element.
        // R8 = pointer to current dst[] element.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ	tablePtr+16(FP), SI
        MOVQ	nSrcFullByte+24(FP), R9

        MOVOU   (SI), X0
        MOVOU   ·Mask0f0f<>(SB), X1

        // set AX to 32 bytes before end of dst[].
        // change R9 to 16 bytes before end of src[].
        SUBQ    $16, R9
        LEAQ    0(R8)(R9*2), AX
        ADDQ    DI, R9

unpackAndReplaceSeqOddSSSE3Loop:
        MOVOU   (DI), X3
        MOVO    X0, X4
        MOVO    X0, X5
        // Isolate high and low nibbles, then parallel-lookup.
        MOVO    X3, X2
        PSRLQ   $4, X3
        PAND    X1, X2
        PAND    X1, X3
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        // Use unpacklo/unpackhi to stitch results together.
        // Odd bytes (1, 3, 5, ...) are in X4, even in X3/X5.
        MOVO    X5, X3
        PUNPCKLBW       X4, X5
        PUNPCKHBW       X4, X3
        MOVOU   X5, (R8)
        MOVOU   X3, 16(R8)
        ADDQ    $16, DI
        ADDQ    $32, R8
        CMPQ    R9, DI
        JG      unpackAndReplaceSeqOddSSSE3Loop

        // Final usually-unaligned read and write.
        MOVOU   (R9), X3
        MOVO    X0, X4
        MOVO    X0, X5
        MOVO    X3, X2
        PSRLQ   $4, X3
        PAND    X1, X2
        PAND    X1, X3
        PSHUFB  X2, X4
        PSHUFB  X3, X5
        MOVO    X5, X3
        PUNPCKLBW       X4, X5
        PUNPCKHBW       X4, X3
        MOVOU   X5, (AX)
        MOVOU   X3, 16(AX)
        RET

TEXT ·cleanASCIISeqInplaceSSSE3Asm(SB),4,$0-16
        MOVQ    ascii8+0(FP), SI
        MOVQ    nByte+8(FP), AX

        MOVOU   ·Capitalizer<>(SB), X0
        MOVOU   ·AllT<>(SB), X1
        MOVOU   ·Mask0f0f<>(SB), X2
        MOVOU   ·All64<>(SB), X3
        MOVOU   ·CleanNoTTable<>(SB), X4
        MOVOU   ·AllNXorT<>(SB), X5

        // set DI to 32 bytes before end of ascii8[].
        LEAQ    -32(SI)(AX*1), DI

        CMPQ    DI, SI
        JLE     cleanASCIISeqInplaceSSSE3Finish

cleanASCIISeqInplaceSSSE3Loop:
        MOVOU   (SI), X6
        // Capitalize.
        PAND    X0, X6
        // Check for 'T's, save positions to X7.
        MOVO    X1, X7
        PCMPEQB X6, X7
        PAND    X5, X7
        // Check for high bits == 0x40.
        MOVO    X2, X8
        PANDN   X6, X8
        PCMPEQB X3, X8
        // X8 now describes which non-T bases are potentially valid.
        PAND    X8, X6
        // Set everything other than A/C/G to N.
        MOVO    X4, X9
        PSHUFB  X6, X9
        // Xor the Ts back in, save result.
        PXOR    X7, X9
        MOVOU   X9, (SI)
        ADDQ    $16, SI
        CMPQ    DI, SI
        JG      cleanASCIISeqInplaceSSSE3Loop

cleanASCIISeqInplaceSSSE3Finish:
        // These loads usually overlap, so they must both occur before the
        // first write-back.
        ADDQ    $16, DI
        MOVOU   (SI), X6
        MOVOU   (DI), X10
        // Capitalize.
        PAND    X0, X6
        PAND    X0, X10
        // Check for 'T's, save positions to X7/X1.
        MOVO    X1, X7
        PCMPEQB X6, X7
        PCMPEQB X10, X1
        PAND    X5, X7
        PAND    X5, X1
        // Check for high bits == 0x40.
        MOVO    X2, X8
        PANDN   X6, X8
        PANDN   X10, X2
        PCMPEQB X3, X8
        PCMPEQB X3, X2
        // X8/X2 now describe which non-T bases are potentially valid.
        PAND    X8, X6
        PAND    X2, X10
        // Set everything other than A/C/G to N.
        MOVO    X4, X9
        PSHUFB  X6, X9
        PSHUFB  X10, X4
        // Xor the Ts back in, save result.
        PXOR    X7, X9
        PXOR    X1, X4
        MOVOU   X9, (SI)
        MOVOU   X4, (DI)
        RET

TEXT ·cleanASCIISeqNoCapitalizeInplaceSSSE3Asm(SB),4,$0-16
        // A bit different from cleanASCIISeqNoCapitalizeInplaceSSSE3Asm due to
        // how we preserve capitalization.
        MOVQ    ascii8+0(FP), SI
        MOVQ    nByte+8(FP), AX

        MOVOU   ·Capitalizer<>(SB), X0
        MOVOU   ·AllT<>(SB), X1
        MOVOU   ·Maskd0d0<>(SB), X2
        MOVOU   ·All64<>(SB), X3
        MOVOU   ·MaskNoTTable<>(SB), X4
        MOVOU   ·AllN<>(SB), X5

        // set DI to 32 bytes before end of ascii8[].
        LEAQ    -32(SI)(AX*1), DI

        CMPQ    DI, SI
        JLE     cleanASCIISeqNoCapitalizeInplaceSSSE3Finish

cleanASCIISeqNoCapitalizeInplaceSSSE3Loop:
        MOVOU   (SI), X6
        // Save off original bytes to X10, before we perform ACGTacgt check on
        // the capitalized version.
        MOVO    X6, X10
        // Capitalize.
        PAND    X0, X6
        // Check for 'T's, save positions to X7.
        MOVO    X1, X7
        PCMPEQB X6, X7
        // Check for high bits == 0x40/0x60.
        MOVO    X6, X8
        PAND    X2, X8
        PCMPEQB X3, X8
        // X8 now describes which non-T bases have potentially valid high bits.
        PAND    X8, X6
        // Create a mask that marks the positions of the original A/C/G
        // bytes.
        MOVO    X4, X9
        PSHUFB  X6, X9
        // Merge the Ts back in, apply mask to original bytes.
        PXOR    X7, X9
        PAND    X9, X10
        // Fill the rest of the bytes with Ns.
        PANDN   X5, X9
        PXOR    X9, X10
        // Save result, advance to next 16 bytes.
        MOVOU   X10, (SI)
        ADDQ    $16, SI
        CMPQ    DI, SI
        JG      cleanASCIISeqNoCapitalizeInplaceSSSE3Loop

cleanASCIISeqNoCapitalizeInplaceSSSE3Finish:
        // These loads usually overlap, so they must both occur before the
        // first write-back.
        ADDQ    $16, DI
        MOVOU   (SI), X6
        MOVOU   (DI), X11
        // Save original bytes to X10/X11.
        MOVO    X6, X10
        // Save capitalized bytes to X6/X0.
        PAND    X0, X6
        PAND    X11, X0
        // Check for 'T's, save positions to X7/X1.
        MOVO    X1, X7
        PCMPEQB X6, X7
        PCMPEQB X0, X1
        // Check for high bits == 0x40/0x60.
        MOVO    X6, X8
        PAND    X2, X8
        PAND    X0, X2
        PCMPEQB X3, X8
        PCMPEQB X3, X2
        // X8/X2 now describe which non-T bases have potentially valid high
        // bits.
        PAND    X8, X6
        PAND    X2, X0
        // Create masks that mark the positions of the original A/C/G bytes.
        MOVO    X4, X9
        PSHUFB  X6, X9
        PSHUFB  X0, X4
        // Merge the Ts back in, apply mask to original bytes.
        PXOR    X7, X9
        PXOR    X1, X4
        PAND    X9, X10
        PAND    X4, X11
        // Fill the rest of the bytes with Ns.
        PANDN   X5, X9
        PANDN   X5, X4
        PXOR    X9, X10
        PXOR    X4, X11
        // Save results.
        MOVOU   X10, (SI)
        MOVOU   X11, (DI)
        RET

TEXT ·isNonACGTPresentSSE41Asm(SB),4,$0-32
        MOVQ    ascii8+0(FP), SI
        MOVQ    nonTTablePtr+8(FP), AX
        MOVQ    nByte+16(FP), R9

        MOVOU   ·AllT<>(SB), X0
        MOVOU   ·Mask0f0f<>(SB), X1
        MOVOU   ·All64<>(SB), X2
        MOVOU   (AX), X3

        // set DI to 16 bytes before end of ascii8[].
        LEAQ    -16(SI)(R9*1), DI

isNonACGTPresentSSE41Loop:
        MOVOU   (SI), X4
        // Check for 'T's, save positions to X5.
        MOVO    X0, X5
        PCMPEQB X4, X5
        // Check for high bits == 0x40.
        MOVO    X1, X6
        PANDN   X4, X6
        PCMPEQB X2, X6
        // X6 now describes which non-T bases are potentially valid.
        PAND    X6, X4
        MOVOU   X3, X7
        PSHUFB  X4, X7
        // Check for invalid chars, treating Ts as valid.
        PTEST   X7, X5
        JNC     isNonACGTPresentSSE41Found

        ADDQ    $16, SI
        CMPQ    DI, SI
        JG      isNonACGTPresentSSE41Loop

        // Final usually-unaligned read.
        MOVOU   (DI), X4
        // Check for 'T's, save positions to X0.
        PCMPEQB X4, X0
        // Check for high bits == 0x40.
        PANDN   X4, X1
        PCMPEQB X2, X1
        // X1 now describes which non-T bases are potentially valid.
        PAND    X1, X4
        PSHUFB  X4, X3
        // Check for invalid chars, treating Ts as valid.
        PTEST   X3, X0
        JNC     isNonACGTPresentSSE41Found
        MOVQ    $0, ret+24(FP)
        RET
isNonACGTPresentSSE41Found:
        MOVQ    $1, ret+24(FP)
        RET

TEXT ·asciiToSeq8SSSE3Asm(SB),4,$0-24
        // Main loop is similar to that of cleanASCIISeqInplaceSSSE3Asm.
        MOVQ    dst+0(FP), R8
        MOVQ    src+8(FP), DI
        MOVQ    nByte+16(FP), R9

        MOVOU   ·Capitalizer<>(SB), X0
        MOVOU   ·AllT<>(SB), X1
        MOVOU   ·Mask0f0f<>(SB), X2
        MOVOU   ·All64<>(SB), X3
        MOVOU   ·ASCIIToSeq8NoTTable<>(SB), X4
        MOVOU   ·ASCIIToSeq8NXorT<>(SB), X5

        // set AX to 16 bytes before end of src[].
        // change R9 to 16 bytes before end of dst[].
        SUBQ    $16, R9
        LEAQ    0(DI)(R9*1), AX
        ADDQ    R8, R9

asciiToSeq8SSSE3Loop:
        MOVOU   (DI), X6
        // Capitalize.
        PAND    X0, X6
        // Check for 'T's, save positions to X7.
        MOVO    X1, X7
        PCMPEQB X6, X7
        PAND    X5, X7
        // Check for high bits == 0x40.
        MOVO    X2, X8
        PANDN   X6, X8
        PCMPEQB X3, X8
        // X8 now describes which non-T bases are potentially valid.
        PAND    X8, X6
        // Set everything other than A/C/G to N.
        MOVO    X4, X9
        PSHUFB  X6, X9
        // Xor the Ts back in, save result.
        PXOR    X7, X9
        MOVOU   X9, (R8)
        ADDQ    $16, DI
        ADDQ    $16, R8
        CMPQ    AX, DI
        JG      asciiToSeq8SSSE3Loop

        // Final usually-unaligned read and write.
        MOVOU   (AX), X6
        PAND    X0, X6
        PCMPEQB X6, X1
        PAND    X5, X1
        PANDN   X6, X2
        PCMPEQB X3, X2
        PAND    X2, X6
        PSHUFB  X6, X4
        PXOR    X1, X4
        MOVOU   X4, (R9)
        RET
