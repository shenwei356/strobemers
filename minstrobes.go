package strobemers

import (
	"math"
)

// MinStrobes is a iterator for MinStrobes
type MinStrobes struct {
	seq *[]byte // DNA sequence

	n    int // strobemer order
	l    int // strobes length
	wMin int // minimum window offset
	wMax int // maximum window offset

	idx, idx2, idx3     int    // indexes of m1, m2, m3
	hash1, hash2, hash3 uint64 // hash value of m1, m2, m3

	hashes []uint64 // precomputed ntHash values of l-mers

	minlocs   []int    // locations of min hash
	minhashes []uint64 // minhashes of window [i-w,i]

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

// NewMinStrobes creates a MinStrobes iterator.
// Parametems:
//     n    - strobemer order
//     l    - strobes length
//     wMin - minimum window offset, wMin > 0
//     wMax - maximum window offset, wMin <= wMax.
func NewMinStrobes(seq *[]byte, n int, l int, wMin int, wMax int) (*MinStrobes, error) {
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

	ms := &MinStrobes{
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
	ms.hashes, err = computeHashes(seq, l)
	if err != nil {
		return nil, err
	}

	ms.minlocs, ms.minhashes = computeMinHashes(ms.hashes, wMax-wMin+1)

	return ms, err
}

// SetPrime sets the prime number (q) in minimizing h(m)+h(mj) mod q.
// In this package, we use (h(m)+h(mj)) & q, where q = roundup(q) - 1.
// The value should not be too small, at least 256.
func (ms *MinStrobes) SetPrime(q uint64) {
	if q < 256 {
		q = 256
	}
	ms.prime = roundup64(q) - 1
}

// SetWindowShrink decides whether shrink the search window at positions
// near the end of the sequence. Default is true.
func (ms *MinStrobes) SetWindowShrink(shrink bool) {
	ms.shrinkWindow = shrink
}

// Index returns the current index (0-based) of strobemers
func (ms *MinStrobes) Index() int {
	return ms.idx - 1
}

// Indexes returns current indexes (0-based) of strobes
func (ms *MinStrobes) Indexes() []int {
	return []int{ms.idx - 1, ms.idx2, ms.idx3}
}

// Next returns the next hash value of randstrobe
func (ms *MinStrobes) Next() (uint64, bool) {
	switch ms.n {
	case 2:
		return ms.nextOrder2()
	case 3:
		return ms.nextOrder3()
	default:
	}

	return 0, false
}

func (ms *MinStrobes) nextOrder2() (uint64, bool) {
	if ms.idx > ms.endIdx {
		return 0, false
	}

	ms.wStart = ms.idx + ms.wMin
	ms.wEnd = ms.idx + ms.wMax

	// for positions near the end of the sequence, shrink the window size from the right
	if ms.wEnd > ms.endHash {
		if !ms.shrinkWindow {
			return 0, false
		}
		ms.wEnd = ms.endHash

		// fmt.Printf("i:%d, window (%d-%d)\n", ms.idx, ms.wStart, ms.wEnd)

		ms.hash1 = ms.hashes[ms.idx]
		ms.hash2 = math.MaxUint64
		for ms.i = ms.wStart; ms.i <= ms.wEnd; ms.i++ {
			ms.hash = ms.hashes[ms.i]
			if ms.hash < ms.hash2 {
				ms.idx2 = ms.i
				ms.hash2 = ms.hash
			}
		}
		// For 1) asymmetry, 2) avoid value overflow
		ms.hash2 = ms.hash1/2 + ms.hashes[ms.idx2]/3
	} else { // use precomputed min hashes
		ms.hash1 = ms.hashes[ms.idx]
		ms.idx2 = ms.minlocs[ms.wEnd]
		ms.hash2 = ms.hash1/2 + ms.minhashes[ms.wEnd]/3
	}

	ms.idx++
	return ms.hash2, true
}

func (ms *MinStrobes) nextOrder3() (uint64, bool) {
	if ms.idx > ms.endIdx {
		return 0, false
	}

	ms.w2Start = ms.idx + ms.wMax + ms.wMin
	ms.w2End = ms.idx + ms.wMax<<1
	if ms.w2Start > ms.endHash {
		return 0, false
	}

	ms.wStart = ms.idx + ms.wMin
	ms.wEnd = ms.idx + ms.wMax

	// use precomputed min hashes
	ms.hash1 = ms.hashes[ms.idx]
	ms.idx2 = ms.minlocs[ms.wEnd]
	ms.hash2 = ms.hash1/3 + ms.minhashes[ms.wEnd]/4

	// for positions near the end of the sequence, shrink the last window size from the right
	if ms.w2End > ms.endHash {
		if !ms.shrinkWindow {
			return 0, false
		}
		ms.w2End = ms.endHash

		ms.hash3 = math.MaxUint64
		for ms.i = ms.w2Start; ms.i <= ms.w2End; ms.i++ {
			ms.hash = (ms.hash2 + ms.hashes[ms.i]) & ms.prime
			if ms.hash < ms.hash3 {
				ms.idx3 = ms.i
				ms.hash3 = ms.hash
			}
		}
		ms.hash3 = ms.hash2 + ms.hashes[ms.idx3]/5
	} else {
		ms.idx3 = ms.minlocs[ms.w2End]
		ms.hash3 = ms.hash2 + ms.minhashes[ms.w2End]/5
	}

	// fmt.Printf("i:%d, window (%d-%d)\n", ms.idx, ms.wStart, ms.wEnd)
	// fmt.Printf("i:%d, window2 (%d-%d)\n", ms.idx, ms.w2Start, ms.w2End)

	ms.idx++
	return ms.hash3, true
}
