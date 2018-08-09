// Copyright 2018 GRAIL, Inc. All rights reserved.
// Use of this source code is governed by the Apache 2.0
// license that can be found in the LICENSE file.

// Package pam implements PAM reader and writer. PAM is a more compact and
// faster alternative to BAM.
//
// Most people, however, will want to use the bamprovider
// (https://godoc.org/github.com/grailbio/bio/encoding/bamprovider) read PAM
// instead.  The bamprovider works for both BAM and PAM files transparently.
//
// REAMDE.md (https://github.com/grailbio/bio/blob/master/encoding/pam/README.md) contains
// More detailed information about the PAM file format.
package pam
