package strobemers

import "math"

// RandStrobes is a iterator for randstrobes
type RandStrobes struct {
	seq *[]byte // DNA sequence

	n    int // strobemer order
	l    int // strobes length
	wMin int // minimum window offset
	wMax int // maximum window offset

	idx, idx2, idx3     int    // indexes of m1, m2, m3
	hash1, hash2, hash3 uint64 // hash value of m1, m2, m3

	hashes []uint64 // precomputed ntHash values of l-mers

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

// NewRandStrobes creates a RandStrobes iterator.
// Parameters:
//     n    - strobemer order
//     l    - strobes length
//     wMin - minimum window offset, wMin > 0
//     wMax - maximum window offset, wMin <= wMax.
func NewRandStrobes(seq *[]byte, n int, l int, wMin int, wMax int) (*RandStrobes, error) {
	if seq == nil || len(*seq) == 0 {
		return nil, ErrInvalidSequence
	}
	if n < 2 {
		return nil, ErrInvalidOrder
	}
	if n > 3 {
		return nil, ErrOrderNotSupported
	}
	if len(*seq) < (n-1)*(wMax+1) {
		return nil, ErrSequenceTooShort
	}
	if l < 1 {
		return nil, ErrStrobeLengthTooSmall
	}
	if !(wMin > 0 && wMax > 0 && wMin <= wMax) {
		return nil, ErrInvalidWindowOffsets
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

	var err error
	rs.hashes, err = computeHashes(seq, l)

	return rs, err
}

// SetPrime sets the prime number (q) in minimizing h(m)+h(mj) mod q.
// In this package, we use (h(m)+h(mj)) & q, where q = roundup(q) - 1.
// The value should not be too small, at least 256.
func (rs *RandStrobes) SetPrime(q uint64) {
	if q < 256 {
		q = 256
	}
	rs.prime = roundup64(q) - 1
}

// SetWindowShrink decides whether shrink the search window at positions
// near the end of the sequence. Default is true.
func (rs *RandStrobes) SetWindowShrink(shrink bool) {
	rs.shrinkWindow = shrink
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
	rs.hash2 = math.MaxUint64
	for rs.i = rs.wStart; rs.i <= rs.wEnd; rs.i++ {
		rs.hash = (rs.hash1 + rs.hashes[rs.i]) & rs.prime
		if rs.hash < rs.hash2 {
			rs.idx2 = rs.i
			rs.hash2 = rs.hash
		}
	}
	rs.hash2 = rs.hash1/2 + rs.hashes[rs.idx2]/3

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
	rs.hash2 = math.MaxUint64
	for rs.i = rs.wStart; rs.i <= rs.wEnd; rs.i++ {
		rs.hash = (rs.hash1 + rs.hashes[rs.i]) & rs.prime
		if rs.hash < rs.hash2 {
			rs.idx2 = rs.i
			rs.hash2 = rs.hash
		}
	}
	rs.hash2 = rs.hash1/3 + rs.hashes[rs.idx2]/4

	rs.hash3 = math.MaxUint64
	for rs.i = rs.w2Start; rs.i <= rs.w2End; rs.i++ {
		rs.hash = (rs.hash2 + rs.hashes[rs.i]) & rs.prime
		if rs.hash < rs.hash3 {
			rs.idx3 = rs.i
			rs.hash3 = rs.hash
		}
	}
	rs.hash3 = rs.hash2 + rs.hashes[rs.idx3]/5

	rs.idx++
	return rs.hash3, true
}
