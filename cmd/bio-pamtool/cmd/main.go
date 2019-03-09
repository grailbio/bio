package cmd

import (
	"fmt"
	"log"
	"strings"

	"github.com/grailbio/base/cmdutil"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/bio/encoding/converter"
	"github.com/grailbio/bio/encoding/pam"
	"v.io/x/lib/cmdline"
)

func newCmdView() *cmdline.Command {
	cmd := &cmdline.Command{
		Name:     "view",
		Short:    "View PAM file metadata",
		ArgsName: "path",
	}
	flags := viewFlags{
		bamIndex:   cmd.Flags.String("index", "", "Input BAM index filename. By default set to input bampath + .bai"),
		headerOnly: cmd.Flags.Bool("header", false, "Print only the header in SAM format"),
		withHeader: cmd.Flags.Bool("with-header", false, "Print header before body"),
		regions: cmd.Flags.String("regions", "", `A comma-separated list of regions to show.
Format of each region is either 'chr:begin-end' or 'chr0:pos0:seq0-chr1:pos1:seq1'.

The first format is the same as samtool's. [begin,end] is a 1-based, closed
interval. For example, 'chr1:123-456' will show reads on chr1, at starting
alignment positions in range [123, 456].

The second format specifies the (chromosome, position, sequence) range as a
0-based, half-open interval. The sequence is a 0-based index that disambiguates
when multiple reads are aligned at the same (chromosome, position).  For
example, 'chr1:123:0-chr3:456:10'. An empty 'chr' part means unmapped reads,
e.g., ':0:1000-:0:2000' will show 1000th to 2000th (0-based) unmapped reads.`),
		filter: cmd.Flags.String("filter", "", filterHelp),
	}
	cmd.Runner = cmdutil.RunnerFunc(func(env *cmdline.Env, argv []string) error {
		if len(argv) != 1 {
			return fmt.Errorf("view takes one pathname argument, but got %v", argv)
		}
		return view(flags, argv[0])
	})
	return cmd
}

func newCmdFlagstat() *cmdline.Command {
	cmd := &cmdline.Command{
		Name:     "flagstat",
		Short:    "Show stats of either a PAM or a BAM file. This command is a clone of 'samtools flagstat'.",
		ArgsName: "path",
	}
	cmd.Runner = cmdutil.RunnerFunc(func(env *cmdline.Env, argv []string) error {
		if len(argv) != 1 {
			return fmt.Errorf("flagstat takes one pathname argument, but got %v", argv)
		}
		return flagstat(argv[0])
	})
	return cmd
}

func newCmdConvert() *cmdline.Command {
	cmd := &cmdline.Command{
		Name:     "convert",
		Short:    "Convert between BAM and PAM",
		ArgsName: "srcpath destpath",
	}
	baiFlag := cmd.Flags.String("index", "", "Input BAM index filename. By default, set to input bampath + .bai")
	bytesPerShardFlag := cmd.Flags.Int64("bytes-per-shard", 4<<30, "A goal size of a PAM file shard")
	bytesPerBlockFlag := cmd.Flags.Int("bytes-per-block", 8<<20, "A goal size of a PAM recordio block")
	formatFlag := cmd.Flags.String("format", "", `
Output file format. Value is either \"bam\" or \"pam\".
If empty, the format is guessed from the input file
(if the input is bam, output is pam and vice versa).`)
	transformersFlag := cmd.Flags.String("transformers", "", `Comma-separated list of transformers to apply during PAM generation.
For example, "-transform=zstd 20".`)
	cmd.Runner = cmdutil.RunnerFunc(func(env *cmdline.Env, argv []string) error {
		if len(argv) != 2 {
			return fmt.Errorf("convert takes srcpath destpath, but found %v", argv)
		}
		srcPath := argv[0]
		destPath := argv[1]
		destFormat := bamprovider.Unknown
		if *formatFlag != "" {
			destFormat = bamprovider.ParseFileType(*formatFlag)
			if destFormat == bamprovider.Unknown {
				return fmt.Errorf("unknown output format \"%s\"", *formatFlag)
			}
		} else {
			switch bamprovider.GuessFileType(srcPath) {
			case bamprovider.BAM:
				destFormat = bamprovider.PAM
			case bamprovider.PAM:
				destFormat = bamprovider.BAM
			}
		}
		switch destFormat {
		case bamprovider.PAM:
			transformers := []string{}
			if *transformersFlag != "" {
				transformers = strings.Split(*transformersFlag, ",")
			}
			return converter.ConvertToPAM(pam.WriteOpts{
				MaxBufSize:   *bytesPerBlockFlag,
				Transformers: transformers,
			}, destPath, srcPath, *baiFlag, *bytesPerShardFlag)
		case bamprovider.BAM:
			p := bamprovider.NewProvider(srcPath, bamprovider.ProviderOpts{Index: *baiFlag})
			err := converter.ConvertToBAM(destPath, p)
			if e := p.Close(); e != nil && err == nil {
				err = e
			}
			return err
		default:
			return fmt.Errorf("cannot determine the output format for conversion from %s to %s",
				srcPath, destPath)
		}
	})
	return cmd
}

func newCmdChecksum() *cmdline.Command {
	cmd := &cmdline.Command{
		Name: "checksum",
		Short: `Compute a checksum of a BAM or PAM file.
The checksum is a JSON string describing the summary of various attributes of the reads`,
		ArgsName: "path",
	}
	opts := checksumOpts{}
	cmd.Flags.StringVar(&opts.baiPath, "index", "", "Input BAM index filename. By default, set to input BAM filename + .bai")
	cmd.Flags.BoolVar(&opts.name, "name", false, "Checksum the name field")
	cmd.Flags.BoolVar(&opts.tempLen, "templen", false, "Checksum the templen field")
	cmd.Flags.BoolVar(&opts.seq, "seq", false, "Checksum the seq field")
	cmd.Flags.BoolVar(&opts.cigar, "cigar", false, "Checksum the cigar field")
	cmd.Flags.BoolVar(&opts.aux, "aux", false, "Checksum the aux field")
	cmd.Flags.BoolVar(&opts.mapQ, "mapq", false, "Checksum the mapq field")
	cmd.Flags.BoolVar(&opts.matePos, "matePos", false, "Checksum the mateRef and matePos fields")
	cmd.Flags.BoolVar(&opts.qual, "qual", false, "Checksum the qual fields")
	cmd.Flags.BoolVar(&opts.all, "all", false, "Checksum the all the fields")
	cmd.Runner = cmdutil.RunnerFunc(func(env *cmdline.Env, argv []string) error {
		if len(argv) != 1 {
			return fmt.Errorf("verify takes a path, but found %v", argv)
		}
		return checksum(argv[0], opts)
	})
	return cmd
}

func Run() {
	log.SetFlags(log.Ldate | log.Ltime | log.Lmicroseconds | log.Lshortfile)
	cmdline.HideGlobalFlagsExcept()
	cmdline.Main(
		&cmdline.Command{
			Name:     "bio-pamtool",
			Short:    "Tools for working with PAM format files",
			LookPath: false,
			Children: []*cmdline.Command{
				newCmdConvert(),
				newCmdFlagstat(),
				newCmdView(),
				newCmdChecksum(),
			},
		})
}
