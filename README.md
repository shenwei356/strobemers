# Strobemers in Go

[![GoDoc](https://godoc.org/github.com/shenwei356/strobemers?status.svg)](https://godoc.org/github.com/shenwei356/strobemers)
[![Go Report Card](https://goreportcard.com/badge/github.com/shenwei356/strobemers)](https://goreportcard.com/report/github.com/shenwei356/strobemers)

## Introduction

This is a Go implementation of the [strobemers](https://github.com/ksahlin/strobemers).

The hashing functions is [ntHash](https://github.com/will-rowe/nthash/).


## Differences compared to the original implementation

see [discussion](https://github.com/ksahlin/strobemers/issues/2)

item                |orginal                |this                              |comment
:-------------------|:----------------------|:---------------------------------|:------------------------
window range        |`w_min < w_max`        |`w_min <= w_max`                  |allowing a fixed position
shrinking window    |all `w_min` and `w_max`|optional shrinking last `w_max`   |see figures below
number of strobemers|`len(seq)-n*l+1`       |`len(seq)-n*l+1-(n-1)*l`          |window shrinked
number of strobemers|                       |`len(seq)-n*l+1-(n-1)*(l+w_min-1)`|window not shrinked

<img src="illustration_randstrobes_order2.jpg" width="750"  />

<img src="illustration_randstrobes_order3.jpg" width="750" />


## Installation

    go get github.com/shenwei356/strobemers

## Usage


```go
n := 2
l := 3
w_min := 3
w_max := 5
rs, err := NewRandStrobes(seq, n, l, w_min, w_max)
checkError(err)

var hash uint64
var ok bool
var i1, i2, i3 int // 0-based indexes m1, m2, m3
for {
    hash, ok = rs.Next()
    if !ok {
        break
    }

    i1, i2, i3 = rs.Index()
    ...
}

```

## Similar Projects

- [strobemer_cpptest](https://github.com/BGI-Qingdao/strobemer_cpptest)