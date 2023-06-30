# Pendulum

This is a program I used to generate the motion of a 3D pendulum, observing the path from a birds-eye view.

## Process

Initially we can derive the equations of motion of the pendulum using the second Euler-Lagrange equation. 
By referring to each position of the particle as defined by its zenith, azimuth, and the (fixed) length
of the string to which the particle is attached, we get two differential equations, one of which we can
reduce to a first order differential equation, and the other which remains a second order differential equation.
By then solving both numerically with sufficient precision we then get a relatively accurate approximation.

After calculating this, we have two different ways of plotting. The first is to just plot a pixel for each
position the particle is in at each point, and the second is to draw a line between each consecutive point.
Honestly there's not too much difference between the two from what I've seen, so I've generally opted for 
former, given that it's a bit easier to work with when applying additional operations.

To get decently nice colours out of this, we can get a notion of density from incrementing at each pixel the
number of times a particle runs over it. We then run a convolution over this data to determine the colour of 
the pixel according to some predefined palette, to develop the image a bit from just being black and white.

To calculate the palettes I've tried a few different methods. Initially I went for using monotonic cubic
interpolation, but this ended up giving weird results, so in the end I just opted for linear interpolation
using RGB (tried using HSV for a bit, but that also gave quite weird results due to the angle system).

The colour palettes are mostly derived from common ggplot data visualisation palettes.
