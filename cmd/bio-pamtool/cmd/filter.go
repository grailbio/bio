package cmd

// Functions for parsing and evaluating the --filter expressions.
// Syntax is very similar to sambamba's:
//  https://github.com/biod/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax.
//
// TODO(saito): Matching tag, seq, cigar, and qual.
import (
	"bytes"
	"errors"
	"fmt"
	"go/ast"
	"go/parser"
	"go/token"
	"log"
	"regexp"
	"strconv"

	"github.com/grailbio/hts/sam"
)

const filterHelp = `Filter expression defines a boolean condition on a single record.

EXAMPLES:
   map_quality >= 60 && sequence_length < 150
   (paired && first_of_pair) || unmapped
   re(ref_name, "^name:[0-9]+$")

SYNTAX:

  Expressions are parsed using the Go parser. The operator precedence rules
  follow Go's.

  expr = intliteral | stringliteral
       re(expr, regexp) |  // Partial regex match.
       binary_op | equality_op
       logical_op
       (expr) |
       symbol

  # Args to a binary op can be integers, strings.
  # The two args must be of the same type.
  binary_op = expr > expr | expr >= expr | expr < expr | expr <= expr

  # Args to an equality op can be integers, strings, or bools.
  # The two args must be of the same type.
  equality_op = expr == expr |
        expr != expr

  logical_op = expr && expr |
        expr || expr |
        !expr

  // The following expressions extract a field value from a record.
  symbol = string_field | int_field | boolean_flag

  string_field = ref_name |  // sam.Record.Ref.Name()
       mate_ref_name |       // sam.Record.MateRef.Name()
       rec_name              // sam.Record.Name

  int_field = ref_id |  // sam.Record.Ref.ID()
       position |       // sam.Record.Pos
       mate_ref_id |    // sam.Record.MateRef.ID()
       mate_position |  // sam.Record.MatePos
       sequence_length |  // sam.Record.Seq.Length
       mapping_quality |  // sam.Record.MapQ
       template_length    // sam.Record.TempLen

  // The following expressions extracts values from sam.Record.Flags.
  booolean_flag = paired | proper_pair | unmapped | mate_is_unmapped |
       is_reverse_strand | mate_is_reverse_strand |
       first_of_pair | second_of_pair |  // R1 or R2
       secondary_alignment| failed_quality_control| duplicate | supplementary |
       chimeric

  Note: flag 'chimeric' is shorthand for (paired && !unmapped && !made_unmapped && (ref_id != mate_ref_id))

  intliteral is 0, 1, 0x10, etc.
  stringliteral is "foo", "文字", etc. It supports all golang string escape sequences.

`

// nodeType defines the type of filterExpr node.
type nodeType int

const (
	nodeInvalid  nodeType = iota
	nodeIntConst          // integer literal
	nodeStrConst          // string literal
	nodeNOT               // !
	nodeLAND              // &&
	nodeLOR               // ||
	nodeEQL               // ==
	nodeNEQ               // !=
	nodeGEQ               // >=
	nodeLEQ               // <=
	nodeLSS               // <
	nodeGTR               // >
	nodeRegex             // regex match

	// Field extractors.
	nodeRecName   // sam.Record.Name
	nodeRefName   // sam.Record.Ref.Name()
	nodeRefID     // sam.Record.Ref.ID()
	nodePos       // sam.Record.Pos
	nodeSeqLength // sam.Record.Seq.Length
	nodeMateRefName
	nodeMateRefID
	nodeMatePos
	nodeMapq    // sam.Record.MapQ
	nodeTempLen // sam.Record.TempLen

	// Predicates on sam.Record.Flags bits.
	nodePaired
	nodeProperPair
	nodeUnmapped
	nodeMateUnmapped
	nodeReverse
	nodeMateReverse
	nodeRead1
	nodeRead2
	nodeSecondary
	nodeQCFail
	nodeDuplicate
	nodeSupplementary
	nodeChimeric
)

type valueType int

const (
	valueTypeInt valueType = iota
	valueTypeStr
	valueTypeBool
	valueTypeRegex
)

// Result of evaluating an exprNode.
type exprValue struct {
	vtype     valueType
	intValue  int64
	strValue  string
	boolValue bool
}

// AST node.
type filterExpr struct {
	ntype    nodeType
	vtype    valueType
	x, y     *filterExpr    // Used by unary or binary ops.
	intConst int64          // set if ntype==nodeIntConst
	strConst string         // set if ntype==nodeStrConst
	regexp   *regexp.Regexp // set if ntype==nodeRegexpConst
}

type exprParser struct {
	err error
}

func doassert(cond bool, expr exprValue) {
	if !cond {
		log.Panicf("Broken expr: %+v", expr)
	}
}

func boolValue(v bool) exprValue {
	return exprValue{vtype: valueTypeBool, boolValue: v}
}

func intValue(v int64) exprValue {
	return exprValue{vtype: valueTypeInt, intValue: v}
}

func strValue(v string) exprValue {
	return exprValue{vtype: valueTypeStr, strValue: v}
}

func (expr *filterExpr) evaluate(rec *sam.Record) exprValue {
	switch expr.ntype {
	case nodeIntConst:
		return intValue(expr.intConst)
	case nodeStrConst:
		return strValue(expr.strConst)
	case nodeRegex:
		x := expr.x.evaluate(rec)
		doassert(x.vtype == valueTypeStr, x)
		return boolValue(expr.regexp.MatchString(x.strValue))
	case nodeRecName:
		return strValue(rec.Name)
	case nodeRefName:
		return strValue(rec.Ref.Name())
	case nodeRefID:
		return intValue(int64(rec.Ref.ID()))
	case nodePos:
		return intValue(int64(rec.Pos))
	case nodeSeqLength:
		return intValue(int64(rec.Len()))
	case nodeMateRefName:
		return strValue(rec.MateRef.Name())
	case nodeMateRefID:
		return intValue(int64(rec.MateRef.ID()))
	case nodeMatePos:
		return intValue(int64(rec.MatePos))
	case nodeMapq:
		return intValue(int64(rec.MapQ))
	case nodeTempLen:
		return intValue(int64(rec.TempLen))
	case nodeReverse:
		return boolValue((rec.Flags & sam.Reverse) != 0)
	case nodeMateReverse:
		return boolValue((rec.Flags & sam.MateReverse) != 0)
	case nodePaired:
		return boolValue((rec.Flags & sam.Paired) != 0)
	case nodeProperPair:
		return boolValue((rec.Flags & sam.ProperPair) != 0)
	case nodeUnmapped:
		return boolValue((rec.Flags & sam.Unmapped) != 0)
	case nodeMateUnmapped:
		return boolValue((rec.Flags & sam.MateUnmapped) != 0)
	case nodeRead1:
		return boolValue((rec.Flags & sam.Read1) != 0)
	case nodeRead2:
		return boolValue((rec.Flags & sam.Read2) != 0)
	case nodeSecondary:
		return boolValue((rec.Flags & sam.Secondary) != 0)
	case nodeQCFail:
		return boolValue((rec.Flags & sam.QCFail) != 0)
	case nodeDuplicate:
		return boolValue((rec.Flags & sam.Duplicate) != 0)
	case nodeSupplementary:
		return boolValue((rec.Flags & sam.Supplementary) != 0)
	case nodeChimeric:
		return boolValue((rec.Flags&sam.Paired) != 0 &&
			(rec.Flags&sam.Unmapped) == 0 &&
			(rec.Flags&sam.MateUnmapped) == 0 &&
			(rec.Ref.ID() != rec.MateRef.ID()))
	case nodeNOT:
		x := expr.x.evaluate(rec)
		doassert(x.vtype == valueTypeBool, x)
		return boolValue(!x.boolValue)
	case nodeLAND, nodeLOR:
		x := expr.x.evaluate(rec)
		y := expr.y.evaluate(rec)
		doassert(x.vtype == valueTypeBool, x)
		doassert(y.vtype == valueTypeBool, y)
		if expr.ntype == nodeLAND {
			return boolValue(x.boolValue && y.boolValue)
		}
		return boolValue(x.boolValue || y.boolValue)
	case nodeGEQ, nodeLEQ, nodeLSS, nodeGTR, nodeEQL, nodeNEQ:
		x := expr.x.evaluate(rec)
		y := expr.y.evaluate(rec)
		switch x.vtype {
		case valueTypeInt:
			switch expr.ntype {
			case nodeGEQ:
				return boolValue(x.intValue >= y.intValue)
			case nodeLEQ:
				return boolValue(x.intValue <= y.intValue)
			case nodeLSS:
				return boolValue(x.intValue < y.intValue)
			case nodeGTR:
				return boolValue(x.intValue > y.intValue)
			case nodeEQL:
				return boolValue(x.intValue == y.intValue)
			case nodeNEQ:
				return boolValue(x.intValue != y.intValue)
			}
		case valueTypeStr:
			switch expr.ntype {
			case nodeGEQ:
				return boolValue(x.strValue >= y.strValue)
			case nodeLEQ:
				return boolValue(x.strValue <= y.strValue)
			case nodeLSS:
				return boolValue(x.strValue < y.strValue)
			case nodeGTR:
				return boolValue(x.strValue > y.strValue)
			case nodeEQL:
				return boolValue(x.strValue == y.strValue)
			case nodeNEQ:
				return boolValue(x.strValue != y.strValue)
			}
		default:
			log.Panicf("Illegal type for op %v: %v, %v", expr.ntype, x, y)
		}
	}
	log.Panicf("Unknown expr: %+v", expr)
	return exprValue{}
}

func (p *exprParser) setError(err error) {
	if err != nil && p.err == nil {
		p.err = err
	}
}

func (p *exprParser) doassert(cond bool, message string, node interface{}) {
	if !cond {
		p.err = errors.New(message + ":" + astDebugString(node))
	}
}

func (p *exprParser) parse(node interface{}) *filterExpr {
	switch e := node.(type) {
	case *ast.ParenExpr:
		return p.parse(e.X)
	case *ast.CallExpr:
		fun, ok := e.Fun.(*ast.Ident)
		if !ok {
			p.setError(fmt.Errorf("expect ident, got %v", astDebugString(e.Fun)))
			return nil
		}
		if fun.Name == "re" {
			if len(e.Args) != 2 {
				p.setError(fmt.Errorf("expect two args for re(), but found %v",
					astDebugString(node)))
				return nil
			}
			x := p.parse(e.Args[0])
			y := p.parse(e.Args[1])
			if p.err != nil {
				return nil
			}
			p.doassert(x.vtype == valueTypeStr && y.vtype == valueTypeStr, "Operang for re() must be string", node)
			re, err := regexp.Compile(y.strConst)
			p.setError(err)
			return &filterExpr{
				ntype:  nodeRegex,
				vtype:  valueTypeBool,
				x:      x,
				regexp: re,
			}
		}
	case *ast.UnaryExpr:
		if e.Op == token.NOT {
			x := p.parse(e.X)
			p.doassert(x.vtype == valueTypeBool, "Operang for ! must be bool", node)
			if p.err != nil {
				return nil
			}
			return &filterExpr{
				ntype: nodeNOT,
				vtype: valueTypeBool,
				x:     x,
			}
		}
	case *ast.BinaryExpr:
		ntype := nodeInvalid
		x, y := p.parse(e.X), p.parse(e.Y)
		if p.err != nil {
			return nil
		}
		switch e.Op {
		case token.LAND:
			p.doassert(x.vtype == valueTypeBool && y.vtype == valueTypeBool, "Operands must be boolean", node)
			ntype = nodeLAND
		case token.LOR:
			p.doassert(x.vtype == valueTypeBool && y.vtype == valueTypeBool, "Operands must be boolean", node)
			ntype = nodeLOR
		case token.NEQ, token.EQL:
			if x.vtype == valueTypeStr {
				p.doassert(x.vtype == y.vtype || y.vtype == valueTypeRegex, "Operands must of the same type", node)
			} else {
				p.doassert(x.vtype == y.vtype, "Operands must of the same type", node)
			}
			ntype = nodeEQL
			if e.Op == token.NEQ {
				ntype = nodeNEQ
			}
		case token.GEQ, token.LEQ, token.LSS, token.GTR:
			p.doassert(x.vtype == y.vtype &&
				(x.vtype == valueTypeInt || x.vtype == valueTypeStr),
				"Wrong operand type", node)
			ntype = nodeGEQ
			if e.Op == token.LEQ {
				ntype = nodeLEQ
			} else if e.Op == token.LSS {
				ntype = nodeLSS
			} else if e.Op == token.GTR {
				ntype = nodeGTR
			}
		default:
			p.setError(fmt.Errorf("unknown binary op: %v", astDebugString(node)))
			return &filterExpr{}
		}
		return &filterExpr{
			ntype: ntype,
			vtype: valueTypeBool,
			x:     x,
			y:     y,
		}
	case *ast.BasicLit:
		switch e.Kind {
		case token.STRING:
			v, err := strconv.Unquote(e.Value)
			p.setError(err)
			return &filterExpr{
				ntype:    nodeStrConst,
				vtype:    valueTypeStr,
				strConst: v,
			}
		case token.INT:
			v, err := strconv.ParseInt(e.Value, 0, 64)
			if err != nil {
				p.err = err
				v = -1
			}
			return &filterExpr{
				ntype:    nodeIntConst,
				vtype:    valueTypeInt,
				intConst: v,
			}
		}
	case *ast.Ident:
		switch e.Name {
		case "ref_name":
			return &filterExpr{ntype: nodeRefName, vtype: valueTypeStr}
		case "mate_ref_name":
			return &filterExpr{ntype: nodeMateRefName, vtype: valueTypeStr}
		case "rec_name":
			return &filterExpr{ntype: nodeRecName, vtype: valueTypeStr}
		case "ref_id":
			return &filterExpr{ntype: nodeRefID, vtype: valueTypeInt}
		case "position":
			return &filterExpr{ntype: nodePos, vtype: valueTypeInt}
		case "mate_ref_id":
			return &filterExpr{ntype: nodeMateRefID, vtype: valueTypeInt}
		case "mate_position":
			return &filterExpr{ntype: nodeMatePos, vtype: valueTypeInt}
		case "sequence_length":
			return &filterExpr{ntype: nodeSeqLength, vtype: valueTypeInt}
		case "mapping_quality":
			return &filterExpr{ntype: nodeMapq, vtype: valueTypeInt}
		case "template_length":
			return &filterExpr{ntype: nodeTempLen, vtype: valueTypeInt}
		case "paired":
			return &filterExpr{ntype: nodePaired, vtype: valueTypeBool}
		case "proper_pair":
			return &filterExpr{ntype: nodeProperPair, vtype: valueTypeBool}
		case "unmapped":
			return &filterExpr{ntype: nodeUnmapped, vtype: valueTypeBool}
		case "mate_is_unmapped":
			return &filterExpr{ntype: nodeMateUnmapped, vtype: valueTypeBool}
		case "is_reverse_strand":
			return &filterExpr{ntype: nodeReverse, vtype: valueTypeBool}
		case "mate_is_reverse_strand":
			return &filterExpr{ntype: nodeMateReverse, vtype: valueTypeBool}
		case "first_of_pair":
			return &filterExpr{ntype: nodeRead1, vtype: valueTypeBool}
		case "second_of_pair":
			return &filterExpr{ntype: nodeRead2, vtype: valueTypeBool}
		case "secondary_alignment":
			return &filterExpr{ntype: nodeSecondary, vtype: valueTypeBool}
		case "failed_quality_control":
			return &filterExpr{ntype: nodeQCFail, vtype: valueTypeBool}
		case "duplicate":
			return &filterExpr{ntype: nodeDuplicate, vtype: valueTypeBool}
		case "supplementary":
			return &filterExpr{ntype: nodeSupplementary, vtype: valueTypeBool}
		case "chimeric":
			return &filterExpr{ntype: nodeChimeric, vtype: valueTypeBool}
		}
	}
	p.setError(fmt.Errorf("unknown expr type %v", astDebugString(node)))
	return nil
}

// Pretty-print a golang AST object.
func astDebugString(node interface{}) string {
	out := bytes.Buffer{}
	fset := token.NewFileSet()
	if err := ast.Fprint(&out, fset, node, nil); err != nil {
		panic(err)
	}
	return out.String()
}

// Parse a filter expression.
func parseFilterExpr(str string) (*filterExpr, error) {
	expr, err := parser.ParseExpr(str)
	if err != nil {
		return nil, err
	}
	p := exprParser{}
	node := p.parse(expr)
	if p.err != nil {
		return nil, p.err
	}
	if node.vtype != valueTypeBool {
		return nil, fmt.Errorf("not a boolean expression: %v", astDebugString(expr))
	}
	return node, nil
}

// Given a parsed filter expression and a sam record, check if the record
// matches the expression condition.
func evaluateFilterExpr(expr *filterExpr, rec *sam.Record) bool {
	val := expr.evaluate(rec)
	doassert(val.vtype == valueTypeBool, val)
	return val.boolValue
}
