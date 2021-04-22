package main

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"

	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/strobemers"
	"github.com/will-rowe/nthash"
)

func main() {
	args := os.Args
	if len(args) != 3 {
		checkError(fmt.Errorf("usage: %s query.fasta ref.fasta", os.Args[0]))
	}

	fileQuery, fileRef := args[1], args[2]

	q, _ := filepathTrimExtension(filepath.Base(fileQuery))
	r, _ := filepathTrimExtension(filepath.Base(fileRef))

	type Test struct {
		n    int
		l    int
		wMin int
		wMax int
	}

	tests := []Test{
		{n: 2, l: 10, wMin: 12, wMax: 12},
		{n: 3, l: 7, wMin: 9, wMax: 9},

		{n: 2, l: 10, wMin: 12, wMax: 16},
		{n: 3, l: 7, wMin: 9, wMax: 13},
	}

	seqsQ := readSeqs(fileQuery)
	seqsR := readSeqs(fileRef)

	var kmersQ, kmersR map[uint64]interface{}
	var rstrobesSQ, rstrobesSR, rstrobesQ, rstrobesR map[uint64]interface{}
	var mstrobesSQ, mstrobesSR, mstrobesQ, mstrobesR map[uint64]interface{}
	var kmersInter, rstrobesSInter, rstrobesInter, mstrobesSInter, mstrobesInter int

	fmt.Printf("query\tref\tmethod\tnQuery\tnRef\tnCommon\tqCov\n")
	runtime.GOMAXPROCS(10)
	for _, t := range tests {
		var wg sync.WaitGroup
		wg.Add(10)

		// kmers
		go func() {
			kmersQ = list2map(computeKmers(seqsQ, t.n*t.l))
			wg.Done()
		}()
		go func() {
			kmersR = list2map(computeKmers(seqsR, t.n*t.l))
			wg.Done()
		}()

		// randstrobes
		go func() {
			rstrobesSQ = list2map(computeRandStrobes(seqsQ, t.n, t.l, t.wMin, t.wMax, true))
			wg.Done()
		}()
		go func() {
			rstrobesSR = list2map(computeRandStrobes(seqsR, t.n, t.l, t.wMin, t.wMax, true))
			wg.Done()
		}()
		go func() {
			rstrobesQ = list2map(computeRandStrobes(seqsQ, t.n, t.l, t.wMin, t.wMax, false))
			wg.Done()
		}()
		go func() {
			rstrobesR = list2map(computeRandStrobes(seqsR, t.n, t.l, t.wMin, t.wMax, false))
			wg.Done()
		}()

		// minstrobes
		go func() {
			mstrobesSQ = list2map(computeMinStrobes(seqsQ, t.n, t.l, t.wMin, t.wMax, true))
			wg.Done()
		}()
		go func() {
			mstrobesSR = list2map(computeMinStrobes(seqsR, t.n, t.l, t.wMin, t.wMax, true))
			wg.Done()
		}()
		go func() {
			mstrobesQ = list2map(computeMinStrobes(seqsQ, t.n, t.l, t.wMin, t.wMax, false))
			wg.Done()
		}()
		go func() {
			mstrobesR = list2map(computeMinStrobes(seqsR, t.n, t.l, t.wMin, t.wMax, false))
			wg.Done()
		}()
		wg.Wait()

		// intersection

		var wg2 sync.WaitGroup
		wg2.Add(5)
		go func() {
			kmersInter = intersection(kmersQ, kmersR)
			wg2.Done()
		}()

		go func() {
			rstrobesSInter = intersection(rstrobesSQ, rstrobesSR)
			wg2.Done()
		}()
		go func() {
			rstrobesInter = intersection(rstrobesQ, rstrobesR)
			wg2.Done()
		}()

		go func() {
			mstrobesSInter = intersection(mstrobesSQ, mstrobesSR)
			wg2.Done()
		}()
		go func() {
			mstrobesInter = intersection(mstrobesQ, mstrobesR)
			wg2.Done()
		}()

		wg2.Wait()

		// kmers
		fmt.Printf("%s\t%s\tKmer(%d)\t%d\t%d\t%d\t%.2f\n",
			q, r, t.n*t.l, len(kmersQ), len(kmersR),
			kmersInter, float64(kmersInter)/float64(len(kmersQ))*100)

		// minstrobes
		fmt.Printf("%s\t%s\tMinStrobes(%d,%d,%d,%d,shrink)\t%d\t%d\t%d\t%.2f\n",
			q, r, t.n, t.l, t.wMin, t.wMax, len(mstrobesSQ), len(mstrobesSR),
			mstrobesSInter, float64(mstrobesSInter)/float64(len(mstrobesSQ))*100)
		fmt.Printf("%s\t%s\tMinStrobes(%d,%d,%d,%d)\t%d\t%d\t%d\t%.2f\n",
			q, r, t.n, t.l, t.wMin, t.wMax, len(mstrobesQ), len(mstrobesR),
			mstrobesInter, float64(mstrobesInter)/float64(len(mstrobesQ))*100)

		// randstrobes
		fmt.Printf("%s\t%s\tRankStrobes(%d,%d,%d,%d,shrink)\t%d\t%d\t%d\t%.2f\n",
			q, r, t.n, t.l, t.wMin, t.wMax, len(rstrobesSQ), len(rstrobesSR),
			rstrobesSInter, float64(rstrobesSInter)/float64(len(rstrobesSQ))*100)
		fmt.Printf("%s\t%s\tRankStrobes(%d,%d,%d,%d)\t%d\t%d\t%d\t%.2f\n",
			q, r, t.n, t.l, t.wMin, t.wMax, len(rstrobesQ), len(rstrobesR),
			rstrobesInter, float64(rstrobesInter)/float64(len(rstrobesQ))*100)

		fmt.Printf(" \t \t \t \t \t \t \n")
	}

}

func checkError(e error) {
	if e != nil {
		fmt.Fprintf(os.Stderr, "%s\n", e)
		os.Exit(0)
	}
}

func readSeqs(file string) [][]byte {
	reader, err := fastx.NewDefaultReader(file)
	checkError(err)

	sequences := make([][]byte, 0, 8)

	var record *fastx.Record
	for {
		record, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}

		sequences = append(sequences, record.Seq.Seq)
	}

	return sequences
}

func computeKmers(sequences [][]byte, k int) []uint64 {
	hashes := make([]uint64, 0, 1024)

	var hash uint64
	var ok bool
	var hasher *nthash.NTHi
	var err error
	for _, _seq := range sequences {
		hasher, err = nthash.NewHasher(&_seq, uint(k))
		checkError(err)

		for {
			hash, ok = hasher.Next(true)
			if !ok {
				break
			}

			hashes = append(hashes, hash)
		}
	}

	return hashes
}

func computeRandStrobes(sequences [][]byte, n int, l int, wMin int, wMax int, shrink bool) []uint64 {
	hashes := make([]uint64, 0, 1024)

	var hash uint64
	var ok bool
	var rs *strobemers.RandStrobes
	var err error

	for _, _seq := range sequences {
		rs, err = strobemers.NewRandStrobes(&_seq, n, l, wMin, wMax)
		checkError(err)

		rs.SetWindowShrink(shrink)
		for {
			hash, ok = rs.Next()
			if !ok {
				break
			}

			hashes = append(hashes, hash)
		}
	}

	return hashes
}

func computeMinStrobes(sequences [][]byte, n int, l int, wMin int, wMax int, shrink bool) []uint64 {
	hashes := make([]uint64, 0, 1024)

	var hash uint64
	var ok bool
	var rs *strobemers.MinStrobes
	var err error

	for _, _seq := range sequences {
		rs, err = strobemers.NewMinStrobes(&_seq, n, l, wMin, wMax)
		checkError(err)

		rs.SetWindowShrink(shrink)
		for {
			hash, ok = rs.Next()
			if !ok {
				break
			}

			hashes = append(hashes, hash)
		}
	}

	return hashes
}

func list2map(data []uint64) map[uint64]interface{} {
	m := make(map[uint64]interface{}, len(data))
	for _, k := range data {
		m[k] = struct{}{}
	}
	return m
}

func intersection(m1, m2 map[uint64]interface{}) int {
	n := 0
	var ok bool
	for k := range m1 {
		if _, ok = m2[k]; ok {
			n++
		}
	}
	return n
}

func filepathTrimExtension(file string) (string, string) {
	gz := strings.HasSuffix(file, ".gz") || strings.HasSuffix(file, ".GZ")
	if gz {
		file = file[0 : len(file)-3]
	}
	extension := filepath.Ext(file)
	name := file[0 : len(file)-len(extension)]
	if gz {
		extension += ".gz"
	}
	return name, extension
}
