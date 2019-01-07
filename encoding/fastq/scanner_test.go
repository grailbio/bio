package fastq

import (
	"bytes"
	"testing"
)

const fq = `@NB500956:89:HW2FHBGX2:1:11101:25648:1069 1:N:0:ATCACG
ATACAGGCCTGANCCACTGTGCCCAGNCTANNTNATTANTGAANANAGAATNGTTNTAAATANANNNNNTNTNNNC
+
AAAAAEEEEEEE#EEAEEEEEEEEEE#EEE##E#EEEE#EEEE#E#EEEEE#EEE#EEEAEE#A#####E#E###E
@NB500956:89:HW2FHBGX2:1:11101:13871:1070 1:N:0:ATCACG
CTCAACTCTGAGNCAGACAGAAATACNTTTNNTNTGAGTTACANCNTTCTTTTTCNACATATNCNNNNNTNGNNNT
+
AAAAAEEEEEEE#EEEEEEEEEEEEE#EEE##E#EEEEEEEEE#E#EEEEEEEEE#EAEEEE#A#####E#A###E
@NB500956:89:HW2FHBGX2:1:11101:9975:1070 1:N:0:ATCACG
GAGTAACCACGTNCCCATGGCCACAGNTGANNGNGTCACACCTNANCCGGGAGAGNCAATCCNGNNNNNGNANNNC
+
AAAAAEEEEEEE#EEEEEEEEEAEEE#EEA##E#EEEEEEEE<#E#<EEEEEEEE#<EEEA/#/#####A#E###A
@NB500956:89:HW2FHBGX2:1:11101:20247:1070 1:N:0:ATCACG
GATCGGAAGAGCNCACGTCTGAACTCNAGTNNCNTCCCGATCTNGNATGCCGTCTNCTGCTTNANNNNNANANNNG
+
AAAAAEEEEEEE#EEEEEEEEEEEEE#AEE##E#A////6AE<#E#EEEEEEEEA#A/EE/E#E#####/#E###E
@NB500956:89:HW2FHBGX2:1:11101:17754:1070 1:N:0:ATCACG
CAAGCAACTTACNTTACTTTAGGCTGNAAANNGNCTGCCTGAANTNCCTGCTCACNAATCCCNCNNNNNCNTNNNT
+
AAAAAEEEEEEE#EEAEEEEEEEEEE#EEE##E#EEEEEEEEE#E#EEEEEEEEE#EAEAEA#/#####E#A###E
@NB500956:89:HW2FHBGX2:1:11101:26223:1070 1:N:0:ATCACG
TCAATTTCAGAACTTTTTATTGGTCTNTTCNNGNATTCATCTTNTNCCTGGTTTANTCTTGGNANNNNNTNTNNNT
+
AAAAAEEEEEEEEEEEEEEEEEEEEE#EEA##E#EEEEEEEEE#E#<EAEEEEEE#EEEEEE#E#####E#E###E
`

func stringScanner(s string) *Scanner {
	return NewScanner(bytes.NewReader([]byte(s)), All)
}

func scanErr(s string) error {
	scan := stringScanner(s)
	var r Read
	for scan.Scan(&r) {
	}
	return scan.Err()
}

func TestFASTQ(t *testing.T) {
	s := stringScanner(fq)
	var r Read
	if !s.Scan(&r) {
		t.Fatal(s.Err())
	}
	expect := Read{
		ID:   "@NB500956:89:HW2FHBGX2:1:11101:25648:1069 1:N:0:ATCACG",
		Seq:  "ATACAGGCCTGANCCACTGTGCCCAGNCTANNTNATTANTGAANANAGAATNGTTNTAAATANANNNNNTNTNNNC",
		Unk:  "+",
		Qual: "AAAAAEEEEEEE#EEAEEEEEEEEEE#EEE##E#EEEE#EEEE#E#EEEEE#EEE#EEEAEE#A#####E#E###E",
	}
	if got, want := r, expect; got != want {
		t.Errorf("got %v, want %v", got, want)
	}
	var n int
	for s.Scan(&r) {
		n++
	}
	if got, want := n, 5; got != want {
		t.Errorf("got %v, want %v", got, want)
	}
	if err := s.Err(); err != nil {
		t.Errorf("unexpected error %v", err)
	}
}

func TestBadFASTQ(t *testing.T) {
	if got, want := scanErr("12312#"), ErrInvalid; got != want {
		t.Errorf("got %v, want %v", got, want)
	}
	if got, want := scanErr("@1234\n123"), ErrShort; got != want {
		t.Errorf("got %v, want %v", got, want)
	}
}

func TestWriter(t *testing.T) {
	var (
		s = stringScanner(fq)
		b = new(bytes.Buffer)
		w = NewWriter(b)
		r Read
	)
	for s.Scan(&r) {
		if err := w.Write(&r); err != nil {
			t.Fatal(err)
		}
	}
	if err := s.Err(); err != nil {
		t.Fatal(err)
	}
	if got, want := b.String(), fq; got != want {
		t.Errorf("got %v, want %v", got, want)
	}
}
