//Package strobemers is a Go implementation of the https://github.com/ksahlin/strobemers.

package strobemers

import (
	"fmt"

	"github.com/will-rowe/nthash"
)

// defaultPrimeNumber is the prime number in minimizing h(m)+h(mj) mod q.
// In this package, we use (h(m)+h(mj)) & q, where q = roundup(q) - 1
var defaultPrimeNumber uint64 = (1 << 20) - 1

// RandStrobes is a iterator for randstrobes
type RandStrobes struct {
	seq *[]byte // DNA sequence

	n    int // strobemer order
	l    int // strobes length
	wMin int // minimum window offset
	wMax int // maximum window offset

	idx, idx2, idx3     int    // indexes of m1, m2, m3
	hash1, hash2, hash3 uint64 // hash value of m1, m2, m3

	hasher *nthash.NTHi // nthash
	hashes []uint64     // precomputed ntHash values of l-mers

	endHash int // position of the last l-mer
	endIdx  int // position of the last m1

	wStart, wEnd, w2Start, w2End int // window start and end

	prime uint64

	// shrink the last searching window for positions near the end of sequence.
	shrinkWindow bool

	// tmp variable
	i    int
	hash uint64
}

// ErrOrderNotSupported means a big strobemer order is not supported.
var ErrOrderNotSupported = fmt.Errorf("strobemers: strobemer order not supported")

// NewRandStrobes creates a RandStrobes iterator.
// Parameters:
//     n    - strobemer order
//     l    - strobes length
//     wMin - minimum window offset, wMin > 0
//     wMax - maximum window offset, wMin <= wMax.
func NewRandStrobes(seq *[]byte, n int, l int, wMin int, wMax int) (*RandStrobes, error) {
	if seq == nil || len(*seq) == 0 {
		return nil, fmt.Errorf("strobemers: invalid DNA sequence")
	}
	if n < 2 {
		return nil, fmt.Errorf("strobemers: strobemer order too small")
	}
	if n > 3 {
		return nil, ErrOrderNotSupported
	}
	if len(*seq) < (n-1)*(wMax+1) {
		return nil, fmt.Errorf("strobemers: sequence too short")
	}
	if l < 1 {
		return nil, fmt.Errorf("strobemers: strobe length too small")
	}
	if !(wMin > 0 && wMax > 0 && wMin <= wMax) {
		return nil, fmt.Errorf("strobemers: window offset should be > 0, and wMin <= wMax")
	}

	rs := &RandStrobes{
		seq:  seq,
		n:    n,
		l:    l,
		wMin: wMin,
		wMax: wMax,

		endHash: len(*seq) - l,           // position of the last l-mer
		endIdx:  len(*seq) - l - (n-1)*l, // position of the last m1

		shrinkWindow: true,

		prime: defaultPrimeNumber,
	}

	return rs, rs.computeHashes()
}

// SetPrime sets the prime number (q) in minimizing h(m)+h(mj) mod q.
// In this package, we use (h(m)+h(mj)) & q, where q = roundup(q) - 1
func (rs *RandStrobes) SetPrime(p uint64) {
	rs.prime = roundup64(p) - 1
}

// SetWindowShrink decides whether shrink the search window at positions
// near the end of the sequence. Default is true.
func (rs *RandStrobes) SetWindowShrink(shrink bool) {
	rs.shrinkWindow = shrink
}

func (rs *RandStrobes) computeHashes() error {
	hasher, err := nthash.NewHasher(rs.seq, uint(rs.l))
	if err != nil {
		return err
	}
	rs.hasher = hasher

	hashes := make([]uint64, len(*rs.seq)-rs.l+1)
	var hash uint64
	var ok bool
	var i int
	for {
		hash, ok = hasher.Next(true)
		if !ok {
			break
		}
		hashes[i] = hash
		i++
	}

	if i != len(*rs.seq)-rs.l+1 {
		return fmt.Errorf("strobemers: incomplete hash values")
	}

	rs.hashes = hashes
	return nil
}

// Index returns the current index (0-based) of strobemers
func (rs *RandStrobes) Index() int {
	return rs.idx - 1
}

// Indexes returns current indexes (0-based) of strobes
func (rs *RandStrobes) Indexes() []int {
	return []int{rs.idx - 1, rs.idx2, rs.idx3}
}

// Next returns the next hash value of randstrobe
func (rs *RandStrobes) Next() (uint64, bool) {
	switch rs.n {
	case 2:
		return rs.nextOrder2()
	case 3:
		return rs.nextOrder3()
	default:
	}

	return 0, false
}

func (rs *RandStrobes) nextOrder2() (uint64, bool) {
	if rs.idx > rs.endIdx {
		return 0, false
	}

	rs.wStart = rs.idx + rs.wMin
	rs.wEnd = rs.idx + rs.wMax

	// for positions near the end of the sequence, shrink the window size from the right
	if rs.wEnd > rs.endHash {
		if !rs.shrinkWindow {
			return 0, false
		}
		rs.wEnd = rs.endHash
	}

	// fmt.Printf("i:%d, window (%d-%d)\n", rs.idx, rs.wStart, rs.wEnd)

	rs.hash1 = rs.hashes[rs.idx]
	rs.hash2 = rs.prime + 356 // initializate with a big number
	for rs.i = rs.wStart; rs.i <= rs.wEnd; rs.i++ {
		rs.hash = (rs.hash1 + rs.hashes[rs.i]) & rs.prime
		if rs.hash < rs.hash2 {
			rs.idx2 = rs.i
			rs.hash2 = rs.hash
		}
	}
	rs.hash2 = rs.hash1 - rs.hashes[rs.idx2]

	rs.idx++
	return rs.hash2, true
}

func (rs *RandStrobes) nextOrder3() (uint64, bool) {
	if rs.idx > rs.endIdx {
		return 0, false
	}

	rs.w2Start = rs.idx + rs.wMax + rs.wMin
	rs.w2End = rs.idx + rs.wMax<<1
	if rs.w2Start > rs.endHash {
		return 0, false
	}
	// for positions near the end of the sequence, shrink the last window size from the right
	if rs.w2End > rs.endHash {
		if !rs.shrinkWindow {
			return 0, false
		}
		rs.w2End = rs.endHash
	}

	rs.wStart = rs.idx + rs.wMin
	rs.wEnd = rs.idx + rs.wMax

	// fmt.Printf("i:%d, window (%d-%d)\n", rs.idx, rs.wStart, rs.wEnd)
	// fmt.Printf("i:%d, window2 (%d-%d)\n", rs.idx, rs.w2Start, rs.w2End)

	rs.hash1 = rs.hashes[rs.idx]
	rs.hash2 = rs.prime + 356
	for rs.i = rs.wStart; rs.i <= rs.wEnd; rs.i++ {
		rs.hash = (rs.hash1 + rs.hashes[rs.i]) & rs.prime
		if rs.hash < rs.hash2 {
			rs.idx2 = rs.i
			rs.hash2 = rs.hash
		}
	}
	rs.hash2 = rs.hash1 - rs.hashes[rs.idx2]

	rs.hash3 = rs.prime + 356
	for rs.i = rs.w2Start; rs.i <= rs.w2End; rs.i++ {
		rs.hash = (rs.hash2 + rs.hashes[rs.i]) & rs.prime
		if rs.hash < rs.hash3 {
			rs.idx3 = rs.i
			rs.hash3 = rs.hash
		}
	}
	rs.hash3 = rs.hash2 + rs.hashes[rs.idx3]<<1

	rs.idx++
	return rs.hash3, true
}
