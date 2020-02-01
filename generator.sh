#!/bin/bash

grain=10
grid=128

for i in {1..10}
do
  fname="${grain}grain_${grid}grid${i}"
  seeds_fromRandom -N $grain --grid $grid $grid $grid > "${fname}.seeds"
  geom_fromVoronoiTessellation --grid $grid $grid $grid "${fname}.seeds"
  geom_check "${fname}.geom"
done
