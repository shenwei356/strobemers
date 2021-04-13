# Strobemers in Go

[![GoDoc](https://godoc.org/github.com/shenwei356/strobemers?status.svg)](https://godoc.org/github.com/shenwei356/strobemers)
[![Go Report Card](https://goreportcard.com/badge/github.com/shenwei356/strobemers)](https://goreportcard.com/report/github.com/shenwei356/strobemers)

## Introduction

This is a Go implementation of the [strobemers](https://github.com/ksahlin/strobemers).

The hashing functions is [ntHash](https://github.com/will-rowe/nthash/).

## Installation

    go get github.com/shenwei356/strobemers

## Usage


```go
rs, err := NewRandStrobes(seq, 3, 6, 10, 20)
checkError(err)

var h uint64
var ok bool
var i1, i2, i3 int
for {
    h, ok = rs.Next()
    if !ok {
        break
    }

    i1, i2, i3 = rs.Index()
    ...
}

```