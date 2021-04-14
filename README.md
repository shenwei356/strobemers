# Strobemers in Go

[![GoDoc](https://godoc.org/github.com/shenwei356/strobemers?status.svg)](https://godoc.org/github.com/shenwei356/strobemers).
[![Go Report Card](https://goreportcard.com/badge/github.com/shenwei356/strobemers)](https://goreportcard.com/report/github.com/shenwei356/strobemers)

## Introduction

This is a Go implementation of the [strobemers](https://github.com/ksahlin/strobemers) (only RandStrobes),
with some [differences](#differences).

The implementation of `Randstrobes` has a not-bad performance (2~3X slower) compared to regular k-mer,
while it's 10~20X slower than [ntHash](https://github.com/will-rowe/nthash/). 
see [benchmark](#benchmark).


## Installation

    go get github.com/shenwei356/strobemers

## Quick Start

We followed the code style of [ntHash](https://github.com/will-rowe/nthash/).

```go
n := 2
l := 3
w_min := 3
w_max := 5
rs, err := NewRandStrobes(seq, n, l, w_min, w_max)
checkError(err)

var hash uint64
var ok bool
var i int  // 0-based index
var positions []int // 0-based indexes of all strobes

rs.SetWindowShrink(true)
for {
    hash, ok = rs.Next()
    if !ok {
        break
    }

    i = rs.Index()
    positions = rs.Indexes()
}

```

## Differences

Differences compared to the original implementation

see [discussion](https://github.com/ksahlin/strobemers/issues/2)

item                |orginal                |this                              |comment
:-------------------|:----------------------|:---------------------------------|:----------------------------------------
window range        |`w_min < w_max`        |`w_min <= w_max`                  |allowing a fixed position
shrinking window    |all `w_min` and `w_max`|optional shrinking last `w_max`   |see figures below
number of strobemers|`len(seq)-n*l+1`       |`len(seq)-n*l+1-(n-1)*l`          |window shrinked
number of strobemers|                       |`len(seq)-n*l+1-(n-1)*(l+w_min-1)`|window not shrinked
choice of min hash  |`(h(m)+h(mj))%q`       |`(h(m)+h(mj))&q`                  |`&` is faster than `%` for `q=pow(2,n)-1`

<img src="illustration_randstrobes_order2.jpg" width="750" />

<img src="illustration_randstrobes_order3.jpg" width="750" />

## Benchmark

    $ go test . -bench=Benchmark* -benchmem \
        | grep Bench \
        | perl -pe 's/\s\s+/\t/g' \
        | csvtk cut -Ht -f 1,3-5 \
        | csvtk add-header -t -n test,time,memory,allocs \
        | csvtk pretty -t -r

                                    test           time      memory        allocs
    -------------------------------------   ------------   ---------   -----------
               BenchmarkNTHash/1.00_KB-16     8410 ns/op     48 B/op   1 allocs/op
                BenchmarkKmers/1.00_KB-16    53865 ns/op     32 B/op   1 allocs/op
    BenchmarkRandStrobesOrder2/1.00_KB-16    92350 ns/op   8432 B/op   3 allocs/op
    BenchmarkRandStrobesOrder3/1.00_KB-16   148057 ns/op   8432 B/op   3 allocs/op


## Similar Projects

- [strobemer_cpptest](https://github.com/BGI-Qingdao/strobemer_cpptest)

## References

- [ntHash](http://dx.doi.org/10.1093/bioinformatics/btw397)
- [strobemers](https://doi.org/10.1101/2021.01.28.428549)