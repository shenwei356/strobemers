//Package strobemers is a Go implementation of the https://github.com/ksahlin/strobemers.

package strobemers

import (
	"fmt"
	"strings"
	"testing"
)

func TestRandStrobeOrder2(t *testing.T) {
	_s := "ACGATCTGGTACCTAG"
	s := []byte(_s)

	n := 2
	l := 3
	wMin := 3
	wMax := 5
	rs, err := NewRandStrobes(&s, n, l, wMin, wMax)
	if err != nil {
		t.Error(err)
	}

	var h uint64
	var ok bool
	var i1, i2 int
	for {
		h, ok = rs.Next()
		if !ok {
			break
		}

		i1, i2, _ = rs.Index()
		fmt.Printf("%s len:%d\n", _s, len(_s))
		fmt.Printf("%s%s i1:%d\n", strings.Repeat(" ", i1), _s[i1:i1+l], i1)
		fmt.Printf("%s%s i2:%d\n", strings.Repeat(" ", i2), _s[i2:i2+l], i2)
		fmt.Printf("%s%d\n", strings.Repeat(" ", len(_s)+1), h)
		fmt.Println()
	}
}

func TestRandStrobeOrder3(t *testing.T) {
	_s := "ACGATCTGGTACCTAG"
	s := []byte(_s)

	n := 3
	l := 3
	wMin := 3
	wMax := 5
	rs, err := NewRandStrobes(&s, n, l, wMin, wMax)
	if err != nil {
		t.Error(err)
	}

	var h uint64
	var ok bool
	var i1, i2, i3 int
	for {
		h, ok = rs.Next()
		if !ok {
			break
		}

		i1, i2, i3 = rs.Index()
		fmt.Printf("%s len:%d\n", _s, len(_s))
		fmt.Printf("%s%s i1:%d\n", strings.Repeat(" ", i1), _s[i1:i1+l], i1)
		fmt.Printf("%s%s i2:%d\n", strings.Repeat(" ", i2), _s[i2:i2+l], i2)
		fmt.Printf("%s%s i2:%d\n", strings.Repeat(" ", i3), _s[i3:i3+l], i3)
		fmt.Printf("%s%d\n", strings.Repeat(" ", len(_s)+1), h)
		fmt.Println()
	}
}
