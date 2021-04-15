package strobemers

import (
	"fmt"
	"strings"
	"testing"
)

func TestRandStrobesOrder2(t *testing.T) {
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
	var ps []int
	var i1, i2 int
	for {
		h, ok = rs.Next()
		if !ok {
			break
		}

		if !debug {
			continue
		}

		ps = rs.Indexes()
		i1, i2 = ps[0], ps[1]
		fmt.Printf("%s len:%d\n", _s, len(_s))
		fmt.Printf("%s%s i1:%d\n", strings.Repeat(" ", i1), _s[i1:i1+l], i1)
		fmt.Printf("%s%s i2:%d\n", strings.Repeat(" ", i2), _s[i2:i2+l], i2)
		fmt.Printf("%s%d\n", strings.Repeat(" ", len(_s)+1), h)
		fmt.Println()
	}
}

func TestRandStrobesOrder3(t *testing.T) {
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
	var ps []int
	var i1, i2, i3 int
	for {
		h, ok = rs.Next()
		if !ok {
			break
		}

		if !debug {
			continue
		}

		ps = rs.Indexes()
		i1, i2, i3 = ps[0], ps[1], ps[2]
		fmt.Printf("%s len:%d\n", _s, len(_s))
		fmt.Printf("%s%s i1:%d\n", strings.Repeat(" ", i1), _s[i1:i1+l], i1)
		fmt.Printf("%s%s i2:%d\n", strings.Repeat(" ", i2), _s[i2:i2+l], i2)
		fmt.Printf("%s%s i3:%d\n", strings.Repeat(" ", i3), _s[i3:i3+l], i3)
		fmt.Printf("%s%d\n", strings.Repeat(" ", len(_s)+1), h)
		fmt.Println()
	}
}
