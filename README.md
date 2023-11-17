# 2O
This is a companion code for the paper "Even better sums of squares over quintic and cyclotomic fields" by Vitezslav Kala and Pavlo Yatsyna.

The code, written in Magma, is in the file 2O.magma. Part of the code is Robinson's algorithm (as described by McKee and Smyth in [MS] and implemented in Magma by Gary Greaves at https://github.com/grwgrvs/equiangular17).

There's an output file for some of the computations. Additionally, we include all the degree five polynomials with a house bounded by 2+sqrt{6}, which were computed by James McKee. Originally, there were 9 834 226 polynomials, with a range of traces between -13 and 13. Shifting all the roots to be positive and omitting repetition leaves us with 5 587 445 polynomials in total. Here is the breakdown of the numbers by trace:

6 polynomials of trace 12

94 polynomials of trace 13

609 polynomials of trace 14

2 705 polynomials of trace 15

9 709 polynomials of trace 16

29 018 polynomials of trace 17

73 661 polynomials of trace 18

161 959 polynomials of trace 19

311 643 polynomials of trace 20

529 442 polynomials of trace 21

796 871 polynomials of trace 22

1 065 105 polynomials of trace 23

1 265 959 polynomials of trace 24

1 340 664 polynomials of trace 25

[Note that some of these computations were previously associated with "On Kitaoka's conjecture and lifting problem for universal quadratic forms" (https://arxiv.org/abs/2110.06260)]﻿
