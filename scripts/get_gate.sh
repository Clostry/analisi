#!/bin/bash
do_the_thing() {
  ./gate_batch ../acquisizioni/20151030/cf_1070_nopoli_63h.root ../acquisizioni/20151030/Na_2070_th25.root $1 $2 >> gate_newdata-$2
}
export -f do_the_thing
parallel -v -a lista3.txt -a lista2.txt do_the_thing
