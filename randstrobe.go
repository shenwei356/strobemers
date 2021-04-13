//Package strobemers is a Go implementation of the https://github.com/ksahlin/strobemers.

package strobemers

import (
	"fmt"

	"github.com/will-rowe/nthash"
)

var Prime uint64 = 997

// RandStrobes is a iterator for randstrobes
type RandStrobes struct {
	seq *[]byte // DNA sequence

	n    int // strobemer order
	l    int // strobes length
	wMin int // minimum window offset
	wMax int // maximum window offset
	k    int // k = N * L

	idx, idx2, idx3     int    // indexes of m1, m2, m3
	hash1, hash2, hash3 uint64 // hash value of m1, m2, m3

	hasher *nthash.NTHi // nthash
	hashes []uint64     // precomputed ntHash values of l-mers

	endHash int
	endIdx  int

	wStart, wEnd, w2Start, w2End int // window start and end

	// tmp variable
	i    int
	hash uint64
}

// ErrOrderNotSupported means a big strobemer order is not supported.
var ErrOrderNotSupported = fmt.Errorf("strobemers: strobemer order not supported")

// NewRandStrobes creates a RandStrobes iterator.
// Parameters:
//   n    - strobemer order
//   l    - strobes length
//   wMin - minimum window offset, wMin > 0
//   wMax - maximum window offset, wMin <= wMax.
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

		endHash: len(*seq) - l,
		endIdx:  len(*seq) - l - (n-1)*l,
		k:       n * l,
	}

	rs.computeHashes()

	return rs, nil
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
		return fmt.Errorf("strobemers: something wrong?")
	}

	rs.hashes = hashes
	return nil
}

// Index returns current indexes (0-based) of strobes
func (rs *RandStrobes) Index() (int, int, int) {
	return rs.idx - 1, rs.idx2, rs.idx3
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
	// fmt.Println()
	// fmt.Printf("i: %d, endIdx: %d, endHash: %d\n", rs.idx, rs.endIdx, rs.endHash)
	if rs.idx > rs.endIdx {
		return 0, false
	}

	// window region
	if rs.idx+rs.wMax <= rs.endHash { // middle of sequence
		// fmt.Println("m1 is in middle of sequence")
		rs.wStart = rs.idx + rs.wMin
		rs.wEnd = rs.idx + rs.wMax
	} else { // near the end of sequence
		// fmt.Println("m1 is near the end of sequence")
		rs.wStart = rs.endHash - (rs.wMax - rs.wMin) // move the window to left
		if rs.wStart <= rs.idx {
			rs.wStart = rs.idx + 1 // make sure wMin > 0
		}
		rs.wEnd = rs.idx + rs.wMax
		if rs.wEnd > rs.endHash {
			rs.wEnd = rs.endHash
		}
	}
	// fmt.Printf("i:%d, window (%d-%d)\n", rs.idx, rs.wStart, rs.wEnd)

	rs.hash1 = rs.hashes[rs.idx]
	rs.hash2 = Prime // initializate with a big number
	for rs.i = rs.wStart; rs.i <= rs.wEnd; rs.i++ {
		rs.hash = (rs.hash1 - rs.hashes[rs.i]) % Prime
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
	if rs.idx > rs.endIdx { // end of sequence: endHash- 2*l
		return 0, false
	}

	// window region
	if rs.idx+rs.wMax<<1 <= rs.endHash { // middle of sequence
		rs.w2Start = rs.idx + rs.wMax + rs.wMin
		rs.w2End = rs.idx + rs.wMax + rs.wMax
	} else { // near the end of sequence
		// fmt.Println("m3 is near the end of sequence")
		rs.w2Start = rs.endHash - (rs.wMax - rs.wMin) // move the window to left
		if rs.w2Start <= rs.idx+1 {
			rs.w2Start = rs.idx + 2
		}
		rs.w2End = rs.idx + rs.wMax + rs.wMax
		if rs.w2End > rs.endHash {
			rs.w2End = rs.endHash
		}
	}
	// fmt.Printf("i:%d, window2 (%d-%d)\n", rs.idx, rs.w2Start, rs.w2End)

	if rs.idx+rs.wMax <= rs.endHash { // middle of sequence
		rs.wStart = rs.idx + rs.wMin
		rs.wEnd = rs.idx + rs.wMax
	} else {
		// fmt.Println("m2 is near the end of sequence")
		rs.wStart = rs.endHash - rs.l - (rs.wMax-rs.wMin)<<1 // move the window to left
		if rs.wStart <= rs.idx {
			rs.wStart = rs.idx + 1
		}
		rs.wEnd = rs.idx + rs.wMax
		if rs.wEnd > rs.endHash {
			rs.wEnd = rs.endHash
		}
	}

	// fmt.Printf("i:%d, window (%d-%d)\n", rs.idx, rs.wStart, rs.wEnd)

	rs.hash1 = rs.hashes[rs.idx]
	rs.hash2 = Prime
	for rs.i = rs.wStart; rs.i <= rs.wEnd; rs.i++ {
		rs.hash = (rs.hash1 + rs.hashes[rs.i]) % Prime
		if rs.hash < rs.hash2 {
			rs.idx2 = rs.i
			rs.hash2 = rs.hash
		}
	}
	rs.hash2 = rs.hash1 - rs.hashes[rs.idx2]

	rs.hash3 = Prime
	for rs.i = rs.w2Start; rs.i <= rs.w2End; rs.i++ {
		rs.hash = (rs.hash2 + rs.hashes[rs.i]) % Prime
		if rs.hash < rs.hash3 {
			rs.idx3 = rs.i
			rs.hash3 = rs.hash
		}
	}
	rs.hash3 = rs.hash2 + rs.hashes[rs.idx3]<<1

	rs.idx++
	return rs.hash3, true
}
