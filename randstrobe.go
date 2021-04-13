//Package strobemers is a Go implementation of the https://github.com/ksahlin/strobemers.

package strobemers

import "github.com/will-rowe/nthash"

const (
	prime uint64 = 997
)

// RandStrobes is a iterator for randstrobes
type RandStrobes struct {
	seq *[]byte // DNA sequence

	N    int // strobemer order
	L    int // strobes length
	Wmin int // minimum window offset
	Wmax int // maximum window offset

	Idx int // current index

	hasher *nthash.NTHi // nthash
}

func NewRandStrobes(seq *[]byte, n int, l int, wMin int, wMax int) (*RandStrobes, error) {

	return nil, nil
}

// Next returns the next hash value of randstrobe
func (rs *RandStrobes) Next() (uint64, bool) {

	return 0, true
}
