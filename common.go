package strobemers

import (
	"fmt"
	"sort"

	"github.com/will-rowe/nthash"
)

// defaultPrimeNumber is the prime number in minimizing h(m)+h(mj) mod q.
// In this package, we use (h(m)+h(mj)) & q, where q = roundup(q) - 1
var defaultPrimeNumber uint64 = (1 << 20) - 1

// ------------------------------------------------------------------------
// errors

// ErrOrderNotSupported means a big strobemer order is not supported.
var ErrOrderNotSupported = fmt.Errorf("strobemers: strobemer order not supported")

// ErrInvalidOrder means
var ErrInvalidOrder = fmt.Errorf("strobemers: strobemer order too small")

// ErrInvalidSequence means the given sequence is invalid
var ErrInvalidSequence = fmt.Errorf("strobemers: invalid DNA sequence")

// ErrSequenceTooShort means the sequence is too short
var ErrSequenceTooShort = fmt.Errorf("strobemers: sequence too short")

// ErrStrobeLengthTooSmall means the strobe length is too small
var ErrStrobeLengthTooSmall = fmt.Errorf("strobemers: strobe length too small")

// ErrInvalidWindowOffsets means invalid window offsets
var ErrInvalidWindowOffsets = fmt.Errorf("strobemers: window offset should be > 0, and wMin <= wMax")

// ErrIncompleteHashValues means incomplete hash values
var ErrIncompleteHashValues = fmt.Errorf("strobemers: incomplete hash values")

var ErrPrimeNumberTooSmall = fmt.Errorf("strobemers: the primer number is too small")

// ------------------------------------------------------------------------

func computeHashes(sequence *[]byte, k int) ([]uint64, error) {
	hasher, err := nthash.NewHasher(sequence, uint(k))
	if err != nil {
		return nil, err
	}

	hashes := make([]uint64, len(*sequence)-k+1)
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

	if i != len(*sequence)-k+1 {
		return nil, ErrIncompleteHashValues
	}

	return hashes, nil
}

func computeMinHashes(hashes []uint64, w int) ([]int, []uint64) {
	locs := make([]int, len(hashes))
	if w == 1 {
		for i := range hashes {
			locs[i] = i
		}
		return locs, hashes
	}

	minHashes := make([]uint64, len(hashes))

	var hash uint64
	var i, idxMw, b, e, t int
	var i2v IdxValue
	var flag bool

	buf := make([]IdxValue, 0, w)
	end := len(hashes)
	r := w - 1 // last position in the buffer

	for idx := 0; idx < end; idx++ { // idx is end position of a window
		hash = hashes[idx]

		if idx < r { // front of w
			buf = append(buf, IdxValue{Idx: idx, Val: hash}) // add current hash to buf
			continue
		}

		if idx == r { // position w
			buf = append(buf, IdxValue{Idx: idx, Val: hash}) // add current hash to buf
			sort.Sort(idxValues(buf))

			i2v = buf[0] // the smallest one
			locs[idx] = i2v.Idx
			minHashes[idx] = i2v.Val
			continue
		}

		// find min k-mer

		// remove k-mer not in this window.
		// have to check position/index one by one
		idxMw = idx - w
		for i, i2v = range buf {
			if i2v.Idx == idxMw {
				if i < r { // not the last element
					copy(buf[i:r], buf[i+1:])
				} // happen to be at the end
				buf = buf[:r]
				break
			}
		}

		// add new k-mer
		flag = false
		// using binary search, faster han linear search
		b, e = 0, r-1
		for {
			t = b + (e-b)/2
			if hash < buf[t].Val {
				e = t - 1 // end search here
				if e <= b {
					flag = true
					i = b
					break
				}
			} else {
				b = t + 1 // start here
				if b >= r {
					flag = false
					break
				}
				if b >= e {
					flag = true
					i = e // right here
					break
				}
			}
		}
		if !flag { // it's the biggest one, append to the end
			buf = append(buf, IdxValue{idx, hash})
		} else {
			if hash >= buf[i].Val { // have to check again
				i++
			}
			buf = append(buf, blankI2V) // append one element
			copy(buf[i+1:], buf[i:r])   // move right
			buf[i] = IdxValue{idx, hash}
		}

		i2v = buf[0] // the smallest one
		locs[idx] = i2v.Idx
		minHashes[idx] = i2v.Val
	}

	return locs, minHashes
}

type IdxValue struct {
	Idx int    // index
	Val uint64 // hash
}

var blankI2V = IdxValue{0, 0}

type idxValues []IdxValue

func (l idxValues) Len() int               { return len(l) }
func (l idxValues) Less(i int, j int) bool { return l[i].Val < l[j].Val }
func (l idxValues) Swap(i int, j int)      { l[i], l[j] = l[j], l[i] }
