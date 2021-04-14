//Package strobemers is a Go implementation of the https://github.com/ksahlin/strobemers.

package strobemers

import (
	"fmt"
	"math/rand"
	"strings"
	"testing"

	"github.com/shenwei356/util/bytesize"
	"github.com/will-rowe/nthash"
	"github.com/zeebo/xxh3"
)

var seqs [][]byte

var bit2base = [4]byte{'A', 'C', 'G', 'T'}

func init() {
	rand.Seed(11)

	sizes := []int{1 << 10} //, 1 << 20} //, 10 << 20}
	seqs = make([][]byte, len(sizes))
	for i, size := range sizes {
		sequence := make([]byte, size)
		for j := 0; j < size; j++ {
			sequence[j] = bit2base[rand.Intn(4)]
		}
		seqs[i] = sequence
	}
}

var _hash uint64
var _k int = 30
var _n2 int = 2
var _l2 int = 15
var _n3 int = 3
var _l3 int = 10
var _w_min int = 15
var _w_max int = 25

func BenchmarkNTHash(b *testing.B) {
	for i := range seqs {
		size := len(seqs[i])
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			for j := 0; j < b.N; j++ {
				var hash uint64
				var ok bool
				var hasher *nthash.NTHi
				var err error

				hasher, err = nthash.NewHasher(&seqs[i], uint(_k))
				if err != nil {
					b.Errorf("fail to create ntHasher iterator. seq length: %d", size)
				}

				for {
					hash, ok = hasher.Next(true)
					if !ok {
						break
					}

					_hash = hash
				}
			}
		})
	}
}

func BenchmarkKmers(b *testing.B) {
	for i := range seqs {
		size := len(seqs[i])
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			for j := 0; j < b.N; j++ {
				var hash, hashrc uint64
				var end int
				var seq []byte
				var rc []byte
				var _i, _j int

				rc = make([]byte, _k)
				seq = seqs[i]
				end = len(seq) - _k + 1

				for i := 0; i < end; i++ {
					hash = xxh3.Hash(seq[i : i+_k])

					// complementary sequence
					for _i = 0; _i < _k; _i++ {
						rc[_i] = cbases[seq[i+_i]]
					}
					// reverse
					for _i, _j = 0, _k-1; _i < _j; _i, _j = _i+1, _j-1 {
						rc[_i], rc[_j] = rc[_j], rc[_i]
					}
					hashrc = xxh3.Hash(rc)

					// canonical kmer
					if hash < hashrc {
						_hash = hash
					} else {
						_hash = hashrc
					}
				}
			}
		})
	}
}

func BenchmarkRandStrobesOrder2(b *testing.B) {
	for i := range seqs {
		size := len(seqs[i])
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			for j := 0; j < b.N; j++ {
				var hash uint64
				var ok bool
				var rs *RandStrobes
				var err error

				rs, err = NewRandStrobes(&seqs[i], _n2, _l2, _w_min, _w_max)
				if err != nil {
					b.Errorf("fail to create RandStrobes. seq length: %d", size)
				}

				for {
					hash, ok = rs.Next()
					if !ok {
						break
					}

					_hash = hash
				}
			}
		})
	}
}

func BenchmarkRandStrobesOrder3(b *testing.B) {
	for i := range seqs {
		size := len(seqs[i])
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			for j := 0; j < b.N; j++ {
				var hash uint64
				var ok bool
				var rs *RandStrobes
				var err error

				rs, err = NewRandStrobes(&seqs[i], _n3, _l3, _w_min, _w_max)
				if err != nil {
					b.Errorf("fail to create RandStrobes. seq length: %d", size)
				}

				for {
					hash, ok = rs.Next()
					if !ok {
						break
					}

					_hash = hash
				}
			}
		})
	}
}

var debug = false

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
		fmt.Printf("%s%s i2:%d\n", strings.Repeat(" ", i3), _s[i3:i3+l], i3)
		fmt.Printf("%s%d\n", strings.Repeat(" ", len(_s)+1), h)
		fmt.Println()
	}
}
