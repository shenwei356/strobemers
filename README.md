# Strobemers in Go

[![GoDoc](https://godoc.org/github.com/shenwei356/strobemers?status.svg)](https://godoc.org/github.com/shenwei356/strobemers)
[![Go Report Card](https://goreportcard.com/badge/github.com/shenwei356/strobemers)](https://goreportcard.com/report/github.com/shenwei356/unikmer)


This is a Go implementation of the [strobemers](https://github.com/ksahlin/strobemers).

## Installation

    go get github.com/shenwei356/strobemers

## Usage


```go
rs, err := NewRandStrobes(seq, 3, 6, 10, 20)
checkError(err)

var h uint64
var ok bool
for {
    h, ok = rs.Next()
    if !ok {
        break
    }

    
}

```