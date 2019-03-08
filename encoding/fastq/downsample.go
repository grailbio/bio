package fastq

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"io"
	"math/rand"
	"os"

	"github.com/grailbio/base/errors"
	"github.com/grailbio/base/file"
	"github.com/klauspost/compress/gzip"
)

const linesPerRead = 4

// FileHandle captures state needed for reading a fastq file.  Thread
// compatible.
type fileHandle struct {
	path string
	f    file.File
	r    io.ReadCloser // gzip or plaintext reader.
	errp *errors.Once  // accumulates any error encountered by the file.
}

func (fh *fileHandle) seek(ctx context.Context, off int64) {
	n, err := fh.f.Reader(ctx).Seek(off, io.SeekStart)
	if err != nil {
		fh.errp.Set(err)
		return
	}
	if n != off {
		panic(fh)
	}
}

func (fh *fileHandle) reader(ctx context.Context) io.ReadCloser {
	if fh.f == nil || fh.r != nil {
		panic("multiple calls to reader")
	}
	r := fh.f.Reader(ctx)
	// TODO(saito) Don't hardcode gzip, instead autodetect on magic prefix
	// bytes. That would suffice for fastq.
	var err error
	if fh.r, err = gzip.NewReader(r); err != nil {
		fh.errp.Set(err)
		fh.r = file.NewError(err)
	}
	return fh.r
}

func (fh *fileHandle) close(ctx context.Context) {
	if fh.r != nil {
		if err := fh.r.Close(); err != nil {
			fh.errp.Set(errors.E(err, "gzip close", fh.path))
		}
	}
	if fh.f != nil {
		if err := fh.f.Close(ctx); err != nil {
			fh.errp.Set(errors.E(err, "close", fh.path))
		}
	}
}

func newFileHandle(ctx context.Context, path string, errp *errors.Once) *fileHandle {
	fh := &fileHandle{path: path, errp: errp}
	var err error
	fh.f, err = file.Open(ctx, path)
	errp.Set(err)
	return fh
}

func doDownsample(ctx context.Context, rate float64, fh1, fh2 *fileHandle, r1Out, r2Out io.Writer, errp *errors.Once) {
	random := rand.New(rand.NewSource(0))
	r1Scanner := bufio.NewScanner(fh1.reader(ctx))
	r2Scanner := bufio.NewScanner(fh2.reader(ctx))
	for {
		r1, r1Err := scanRead(r1Scanner)
		if r1Err != nil && r1Err != io.EOF {
			errp.Set(errors.E(r1Err, "read", fh1.path))
			return
		}
		r2, r2Err := scanRead(r2Scanner)
		if r2Err != nil && r2Err != io.EOF {
			errp.Set(errors.E(r2Err, "read", fh2.path))
			return
		}
		if r1Err == io.EOF && r2Err == io.EOF {
			// Both readers ended after the same number of reads, as expected.
			break
		} else if r1Err == io.EOF {
			errp.Set(errors.E("more reads in R2 input than in R1 input", fh1.path, fh2.path))
			return
		} else if r2Err == io.EOF {
			errp.Set(errors.E("more reads in R1 input than in R2 input", fh1.path, fh2.path))
			return
		}
		if random.Float64() < rate {
			if _, err := r1Out.Write(r1); err != nil {
				errp.Set(errors.E(err, "write R1"))
				return
			}
			if _, err := r2Out.Write(r2); err != nil {
				errp.Set(errors.E(err, "write R2"))
				return
			}
		}
	}
}

// Downsample writes read pairs from the two files to r1Out and r2Out. Read
// pairs will be randomly selected for inclusion in the output at the given
// sampling rate.
func Downsample(ctx context.Context, rate float64, r1Path, r2Path string, r1Out, r2Out io.Writer) error {
	if rate < 0.0 {
		return errors.New("rate must be >= 0.0")
	}
	if rate > 1.0 {
		fmt.Fprintln(os.Stderr, "warning: Downsample called with input rate > 1.0, interpreting it as rate = 1.0")
		rate = 1.0
		// TODO(kshashidhar): Just copy input to output in this case (better fix the asymmetry in current API first)
	}
	e := errors.Once{}
	fh1 := newFileHandle(ctx, r1Path, &e)
	fh2 := newFileHandle(ctx, r2Path, &e)
	if e.Err() == nil {
		doDownsample(ctx, rate, fh1, fh2, r1Out, r2Out, &e)
	}
	return e.Err()
}

// DownsampleToCount writes read pairs from the two files to r1Out and
// r2Out. Read pairs will be randomly selected for inclusion in the output to
// downsample to the given count, approximately.
func DownsampleToCount(ctx context.Context, count int64, r1Path, r2Path string, r1Out, r2Out io.Writer) (err error) {
	if count <= 0 {
		return errors.E("count must be >= 1")
	}
	e := errors.Once{}
	defer func() { err = e.Err() }()

	fh1 := newFileHandle(ctx, r1Path, &e)
	defer fh1.close(ctx)
	fh2 := newFileHandle(ctx, r2Path, &e)
	defer fh2.close(ctx)
	if e.Err() != nil {
		return
	}

	// Translate the sample count to sampling rate by counting the # of reads in
	// the the first 100KiB of the file, then estimating the # of bytes per read.
	f1Stat, err := fh1.f.Stat(ctx)
	if err != nil {
		e.Set(err)
		return
	}
	prefixLen := int64(100 << 10)
	if prefixLen > f1Stat.Size() {
		prefixLen = f1Stat.Size()
	}
	r1PrefixReader, err := gzip.NewReader(&io.LimitedReader{R: fh1.f.Reader(ctx), N: prefixLen})
	if err != nil {
		e.Set(err)
		return
	}
	// Count the # of reads in the first 100KiB
	scanner := bufio.NewScanner(r1PrefixReader)
	nLine := 0
	for scanner.Scan() {
		nLine++
	}
	var (
		nRead = nLine / linesPerRead
		rate  = 1.0
	)
	if nRead > 0 {
		_ = r1PrefixReader.Close() // this should fail with a checksum error.
		// Approximate the # bytes per read, compressed
		approxCompressedReadLen := float64(prefixLen) / float64(nRead)
		rate = float64(count) * approxCompressedReadLen / float64(f1Stat.Size())
	}
	fh1.seek(ctx, 0)
	doDownsample(ctx, rate, fh1, fh2, r1Out, r2Out, &e)
	return
}

func scanRead(scanner *bufio.Scanner) ([]byte, error) {
	var buffer bytes.Buffer
	for i := 0; i < linesPerRead; i++ {
		if !scanner.Scan() {
			if i == 0 && scanner.Err() == nil {
				// Reached end of input.
				return nil, io.EOF
			}
			// Something went wrong.
			if scanner.Err() != nil {
				return nil, scanner.Err()
			}
			return nil, fmt.Errorf("too few lines in FASTQ record: want %d, got %d", linesPerRead, i)
		}
		buffer.WriteString(scanner.Text())
		buffer.WriteString("\n")
	}
	return buffer.Bytes(), nil
}
