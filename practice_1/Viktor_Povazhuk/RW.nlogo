globals [
  reached-border
  reach-border-ticks
]

turtles-own [
  distance-to-zero
  cur-angle-change
  prev-angle-change
  cur-angle
  prev-angle
  prev-xcor
  prev-ycor
  cur-step-size
]

to setup
  clear-all

  ask patches [
    set pcolor white
  ]

  create-turtles num [
    set cur-angle-change (random 360) / 360
    set prev-angle-change cur-angle-change
    pen-down
  ]

  set reached-border false
  set reach-border-ticks 0

  ask turtles [
    set cur-angle heading
    set prev-angle heading
  ]

  reset-ticks
end

to go

  if all? patches [pcolor = yellow] [stop]

  ask turtles [
    set prev-xcor xcor
    set prev-ycor ycor
  ]

  if rw-type = "cd" [
  ask turtles [
    set heading 360 / num-directions * (random num-directions)
    set cur-step-size random-float step-size
    ;wait(0.1)
  ]
  ]

  if rw-type = "levy" [
  ask turtles [
    set heading random 360
    (ifelse
      dist-type = "weibull" [
        set cur-step-size random-weibull 0.5 step-size 0 5
      ]
      dist-type = "log-normal" [
        set cur-step-size random-lognormal 0.5 step-size
      ]
      dist-type = "cauchy" [
        set cur-step-size random-cauchy 0 step-size
      ])
  ]
  ]

  if rw-type = "deterministic" [
  ask turtles [
    set prev-angle-change cur-angle-change
    set cur-angle-change r_param * cur-angle-change * (1 - cur-angle-change)
    set heading heading + cur-angle-change * 360
    set cur-step-size step-size
    ;wait(0.1)
  ]
  ]

;  ask turtles [
;    if patch-at-heading-and-distance heading cur-step-size = nobody and not reached-border[
;      set reached-border true
;      set reach-border-ticks ticks
;    ]
;  ]

  if any? turtles with [((-1 * cur-step-size <= xcor - (max-pxcor + 0.5) and xcor - (max-pxcor + 0.5) <= 0) or (0 <= xcor + (max-pxcor + 0.5) and xcor + (max-pxcor + 0.5) <= cur-step-size) or (-1 * cur-step-size <= ycor - (max-pycor + 0.5) and ycor - (max-pycor + 0.5) <= 0) or (0 <= ycor + (max-pycor + 0.5) and ycor + (max-pycor + 0.5) <= cur-step-size)) and not reached-border]
    [
      set reached-border true
      set reach-border-ticks ticks
  ]

  ask turtles [
    forward cur-step-size
  ]

  ask turtles [
    set distance-to-zero distancexy 0 0
  ]

  ask turtles [
    if [pcolor] of patch-here != "yellow" [
      ask patch-here [set pcolor yellow]
    ]
  ]

  ask turtles [
    set prev-angle cur-angle
    set cur-angle heading
  ]

  if ticks != 0 [plot-angles]

  tick
end

to plot-angles
  set-current-plot "Angles"
  set-current-plot-pen "pen"

  ask turtles [
    if rw-type = "levy" or rw-type = "cd" [plotxy cur-angle prev-angle]
    if rw-type = "deterministic" [plotxy cur-angle-change prev-angle-change]
  ]

  update-plots
end

to-report get-reach-border-time
  report reach-border-ticks
end

to-report get-cover-time
  report ifelse-value all? patches [pcolor = yellow]
  [ticks]
  [0]
end

to-report get-displacement
  report sqrt (sum (map [ i -> i * i ] (reduce [[a b] -> (map + a b)] [list (cur-step-size * dx) (cur-step-size * dy)] of turtles)))
end

to-report get-msd
  report sqrt (sum (map [vec -> sum (map [cor -> cor ^ 2] vec)] [list xcor ycor] of turtles) / count turtles)
end

;; additional probability distributions for NetLogo

to-report random-lognormal [mn sd]
  let lvar ln (1 + (sd ^ 2) / (mn ^ 2))
  let lmn ln mn - (lvar / 2)
  let lsd sqrt lvar
  report exp (lmn + lsd * random-normal 0 1)
end


; This binomial algorithm from
; Devroye. L. 1960. Generating the maximum of independent identically
; distributed random variables. Computers and Mathematics with
; Applications 6, 305-315.
; Based on code from
; https://stackoverflow.com/questions/23561551/a-efficient-binomial-random-number-generator-code-in-java#23574723
to-report random-binomial [n p]
  ;; check this, but ...
  ;; if p > 1 [ print "WARNING: random-binomial p > 1" ]
  ;; p < 0 [ print "WARNING: random-binomial p < 0" ]
  ;; ... forgive it!
  if p = 1 [ report n ]
  if p = 0 [ report 0 ]
  let ln-q ln (1 - p)
  let x 0
  let s 0
  ; also need to avoid x = n
  while [x < n] [
    set s s + ln (random-float 1) / (n - x)
    if s < ln-q [
      report x
    ]
    set x x + 1
  ]
  report x
end


;; negative binomial implemented as a Poisson draw parameterised
;; with a Gamma distributed variate per
;; https://en.wikipedia.org/wiki/Negative_binomial_distribution#Gamma%E2%80%93Poisson_mixture
to-report random-negative-binomial [r p]
  report random-poisson random-gamma r ((1 - p) / p)
end

;; wrapper to get a negative binomial variate with mean m and variance mean ratio vmr
;; the standard nbin parameters are r and p
to-report random-negative-binomial-with-mean-vmr [m vmr]
  let r m / (vmr - 1)
  let p 1 - 1 / vmr
  report random-negative-binomial r p
end



to-report random-multinomial-int [n frequencies]
  report reduce [ [a p] ->
    lput (ifelse-value (p > 0) [random-binomial (n - sum a) p] [0]) a
  ] fput [] conditional-probabilities frequencies
end

;; for a provided list of relative frequencies [f_i]
;; returns the conditional probabilities [f_i / sum_i..n f_i]
to-report conditional-probabilities [frequencies]
  report (map [ [f s] -> ifelse-value (s > 0) [f / s] [0] ]
                frequencies
                cumulative-remainder frequencies)
end



to-report random-multinomial [n frequencies simple-frequencies?]
  let conditional-probs frequencies
  if simple-frequencies? [
    set conditional-probs conditional-probabilities-with-forcing frequencies
  ]
  report reduce [ [a p] ->
    lput (ifelse-value (p > 0) [random-binomial (n - sum a) p] [0]) a
  ] fput [] conditional-probs
end

;; for a provided list of relative frequencies [f_i]
;; returns the conditional probabilities [f_i / sum_i..n f_i]
to-report conditional-probabilities-with-forcing [frequencies]
  let c-p (map [ [f s] -> ifelse-value (s > 0) [f / s] [0] ]
                frequencies
                cumulative-remainder frequencies)
  report replace-item (last-positive c-p) c-p 1
end


to-report last-positive [L]
  report last-position true (map [x -> x > 0] L)
end

to-report last-position [x L]
  report length L - position x reverse L - 1
end


;; for a list [x_i] returns [sum_i..n x_i]
;; e.g. 1 2 3 4 5 6 --> 21 20 18 15 11 6
;; same as reverse dists-cumulative-sum reverse L
to-report cumulative-remainder [L]
  report but-last reduce [ [a b] -> lput (last a - b) a ] fput (list sum L) L
end

to-report dists-cumulative-sum [L]
  report but-first reduce [ [a b] -> lput (last a + b) a ] fput [0] L
end


to-report random-cauchy [locn scl]
  report locn + scl * (tan (-90 + random-float 180))
end

;; https://stackoverflow.com/questions/33611708/random-number-generator-with-generalized-pareto-distribution-and-weilbull-distri
to-report random-weibull [shp scl lower-limit upper-limit]
  let result upper-limit
  while [result >= upper-limit or result < lower-limit] [
    set result scl * (-1 * ln (random-float 1)) ^ (1 / shp)
  ]
  report result
end

to-report random-gamma-with-mean-sd [mn sd]
  let g-scale mn / sd / sd
  let g-shape mn * g-scale
  report random-gamma g-shape g-scale
end

to-report population-standard-deviation [x]
  report sqrt population-variance x
end

to-report population-variance [x]
  let m mean x
  report mean map [xi -> (xi - m) ^ 2] x
end

;; The MIT License (MIT)
;;
;; Copyright (c) 2021 David O'Sullivan
;;
;; Permission is hereby granted, free of charge, to any person
;; obtaining a copy of this software and associated documentation
;; files (the "Software"), to deal in the Software without restriction,
;; including without limitation the rights to use, copy, modify, merge,
;; publish, distribute, sublicense, and/or sell copies of the Software,
;; and to  permit persons to whom the Software is furnished to do so,
;; subject to the following conditions:
;;
;; The above copyright notice and this permission notice shall be included
;; in all copies or substantial portions of the Software.
;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
;; OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
;; THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
;; FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
;; DEALINGS IN THE SOFTWARE.
@#$#@#$#@
GRAPHICS-WINDOW
278
30
591
344
-1
-1
5.0
1
10
1
1
1
0
1
1
1
-30
30
-30
30
0
0
1
ticks
30.0

SLIDER
35
260
207
293
num
num
0
400
100.0
1
1
NIL
HORIZONTAL

CHOOSER
44
169
192
214
rw-type
rw-type
"cd" "levy" "deterministic"
1

SLIDER
36
317
208
350
step-size
step-size
0
2
1.05
0.01
1
NIL
HORIZONTAL

SLIDER
439
497
611
530
r_param
r_param
0
4
3.8
0.1
1
NIL
HORIZONTAL

BUTTON
48
29
121
62
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
148
30
211
63
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
147
78
213
111
30
repeat 30 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
770
10
1133
291
MSD
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot get-msd"
"pen-1" 1.0 0 -2674135 true "" "plot step-size * sqrt ticks"

PLOT
1182
10
1568
292
Displacement
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot get-displacement"

MONITOR
770
315
852
360
cover-time
get-cover-time
17
1
11

MONITOR
769
384
901
429
reach-border-time
get-reach-border-time
17
1
11

PLOT
1170
321
1594
587
Angles
previous
current
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"pen" 1.0 2 -16777216 true "" ""

SLIDER
44
494
216
527
num-directions
num-directions
1
6
5.0
1
1
NIL
HORIZONTAL

CHOOSER
255
489
393
534
dist-type
dist-type
"cauchy" "log-normal" "weibull"
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="cd-run-end" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>get-cover-time</metric>
    <metric>get-reach-border-time</metric>
    <steppedValueSet variable="num" first="100" step="100" last="300"/>
    <steppedValueSet variable="step-size" first="0.5" step="0.5" last="2"/>
    <steppedValueSet variable="num-directions" first="3" step="1" last="6"/>
  </experiment>
  <experiment name="cd-run-step-2" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>get-msd</metric>
    <metric>step-size * sqrt ticks</metric>
    <steppedValueSet variable="num" first="100" step="100" last="300"/>
    <steppedValueSet variable="step-size" first="1" step="0.5" last="2"/>
    <steppedValueSet variable="num-directions" first="3" step="1" last="6"/>
  </experiment>
  <experiment name="determ-run-end" repetitions="6" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="6000"/>
    <metric>get-cover-time</metric>
    <metric>get-reach-border-time</metric>
    <enumeratedValueSet variable="num">
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="step-size" first="0.5" step="0.5" last="2"/>
    <enumeratedValueSet variable="r">
      <value value="0.5"/>
      <value value="3"/>
      <value value="3.6"/>
      <value value="3.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-pxcor">
      <value value="10"/>
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="determ-run-step-msd" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="6000"/>
    <metric>get-msd</metric>
    <metric>step-size * sqrt ticks</metric>
    <enumeratedValueSet variable="num">
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="step-size" first="0.5" step="0.5" last="2"/>
    <enumeratedValueSet variable="r">
      <value value="0.5"/>
      <value value="3"/>
      <value value="3.6"/>
      <value value="3.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-pxcor">
      <value value="10"/>
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="determ-run-step-disp" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="6000"/>
    <metric>get-displacement</metric>
    <enumeratedValueSet variable="num">
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="step-size" first="0.5" step="0.5" last="2"/>
    <enumeratedValueSet variable="r">
      <value value="0.5"/>
      <value value="3"/>
      <value value="3.6"/>
      <value value="3.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-pxcor">
      <value value="10"/>
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="levy-end" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000000"/>
    <metric>get-cover-time</metric>
    <metric>get-reach-border-time</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.5"/>
      <value value="1"/>
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num">
      <value value="50"/>
      <value value="100"/>
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rw-type">
      <value value="&quot;levy&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dist-type">
      <value value="&quot;log-normal&quot;"/>
      <value value="&quot;weibull&quot;"/>
      <value value="&quot;cauchy&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="levy-msd" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3000000"/>
    <metric>get-msd</metric>
    <metric>step-size * sqrt ticks</metric>
    <enumeratedValueSet variable="step-size">
      <value value="0.5"/>
      <value value="1"/>
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num">
      <value value="50"/>
      <value value="100"/>
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dist-type">
      <value value="&quot;log-normal&quot;"/>
      <value value="&quot;weibull&quot;"/>
      <value value="&quot;cauchy&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
