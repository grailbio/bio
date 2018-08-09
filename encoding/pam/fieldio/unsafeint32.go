// Code generated from " ../../../../base/gtl/generate.py --prefix=unsafe -DELEM=int32 --package=fieldio --output=unsafeint32.go ../../../../base/gtl/unsafe.go.tpl ". DO NOT EDIT.
package fieldio

import (
	"reflect"
	"unsafe"
)

// unsafeint32sToBytes casts []int32 to []byte without reallocating.
func unsafeint32sToBytes(src []int32) (d []byte) { // nolint: deadcode
	if len(src) == 0 {
		return nil
	}
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dh := (*reflect.SliceHeader)(unsafe.Pointer(&d))
	const elemSize = int(unsafe.Sizeof(src[0]))
	dh.Data = sh.Data
	dh.Len = sh.Len * elemSize
	dh.Cap = sh.Cap * elemSize
	return d
}

// unsafeBytesToint32s casts []byte to []int32 without reallocating.
func unsafeBytesToint32s(src []byte) (d []int32) { // nolint: deadcode
	if len(src) == 0 {
		return nil
	}
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&src))
	dh := (*reflect.SliceHeader)(unsafe.Pointer(&d))
	const elemSize = int(unsafe.Sizeof(d[0]))
	dh.Data = sh.Data
	dh.Len = sh.Len / elemSize
	dh.Cap = sh.Cap / elemSize
	return d
}
