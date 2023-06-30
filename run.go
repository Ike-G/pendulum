package main

import (
	"fmt"
	"math"
	"math/rand"
	"image"
	"image/color"
	"image/png"
	"os"
	"time"
)

const G = float64(9.8)
const LENGTH = 10000000 // Becomes a solid circle at 10 million
const H = 0.00009765625
const L = float64(10) // H and L should be roughly proportional
const WIDTH = 5000
const HEIGHT = 5000
const NAME = "p72.png"
const TIMELIMIT = 60 * 1e9
const LINES = false
const CONVWIDTH = 4
const MAXCONV = (2 * CONVWIDTH + 1) * (2 * CONVWIDTH + 1)
const REVERSE = true
const SCHEME = 3
// 0 - mako
// 1 - rocket
// 2 - mako2 (dark endpoint)
// 3 - rocket2 (dark endpoint)
// 4 - kovesi pink
// 5 - kovesi blue
const ACTUALDENSITY = false

const RXL = 0 // Default 0
const RYL = 0 // Default 0
const RD = 1000 // Default 1000
const RXR = RXL + RD
const RYR = RYL + RD
const RESTLOWERX = L * (float64(RXL) / 500 - 1)
const RESTUPPERX = L * (float64(RXR) / 500 - 1)
const RESTDX = RESTUPPERX - RESTLOWERX
const RESTLOWERY = L * (float64(RYL) / 500 - 1)
const RESTUPPERY = L * (float64(RYR) / 500 - 1)
const RESTDY = RESTUPPERY - RESTLOWERY

type State struct {
	u, v, udash float64 // v zenith, u azimuth
	c, c2 float64 // Initialising constant and its square
}

func polarToCoord(zenith, azimuth float64) (float64, float64) {
	ct := L * math.Cos(azimuth)
	return ct * math.Sin(zenith), ct * math.Cos(zenith)
}

func coordToPixel(x, y float64) (int, int) {
	return int(float64(WIDTH) * (x - RESTLOWERX) / RESTDX), int(float64(HEIGHT) * (y - RESTLOWERY) / RESTDY)
}

func (s *State) updateState() {
	cosu := math.Cos(s.u)
	if cosu == 0 { panic("Cos(u) == 0") }

	vdash := s.c / math.Pow(cosu, 2)

	ut := s.u + H * s.udash
	// fmt.Printf("du %f\n", H * s.udash)
	for math.Cos(ut) == 0 {
		ut += H * s.udash
	}

	udd := G * cosu / L - s.c2 * math.Sin(s.u) / math.Pow(cosu, 3)
	if udd == 0 {
		fmt.Println("0 azimuth acceleration for azimuth %f", s.u)
		if s.udash == 0 { panic("Oh no")}
	}
	udasht := s.udash + H * udd
	for udasht == 0 {
		udasht += H * udd
	}

	s.v += H * vdash
	s.u = ut
	s.udash = udasht
	// fmt.Printf("New State: %+v\n", State{s.u, s.v, s.udash, s.c, s.c2})
}

func initialise(zenith, zenithspeed, azimuth, azimuthspeed float64) State {
	// Illegal values: Cos(azimuth) == 0
	// Reasonable values: zenith in (-pi, pi), azimuth in (-pi / 2, pi / 2),
	//                    zenithspeed, azimuthspeed in (-pi, pi)
	// zenithspeed == 0 gives a constant zenith, as expected
	// This will either be maintained, or ought never to occur.
	c := zenithspeed * math.Pow(math.Cos(azimuth), 2)
	c2 := math.Pow(c,2)
	return State{azimuth, zenith, azimuthspeed, c, c2}
}

type HSV struct {
	h, s, v float64 // h is an angle in [0,360), s and v are in [0,1]
}

func HSVtoRGB(z HSV) color.RGBA {
	if z.s <= 0. {
		return color.RGBA{uint8(z.v*255), uint8(z.v*255), uint8(z.v*255), 0xff}
	}
	a := z.h
	if a >= 360. {
		a = 0.
	}
	a /= 60.
	b := a - float64(int32(a))
	c := z.v * (1. - z.s)
	d := z.v * (1. - (z.s * b))
	e := z.v * (1. - (z.s * (1. - b)))
	switch int32(a) {
	case 0: return color.RGBA{uint8(z.v*255), uint8(e*255), uint8(c*255), 0xff}
	case 1: return color.RGBA{uint8(d*255), uint8(z.v*255), uint8(c*255), 0xff}
	case 2: return color.RGBA{uint8(c*255), uint8(z.v*255), uint8(e*255), 0xff}
	case 3: return color.RGBA{uint8(c*255), uint8(d*255), uint8(z.v*255), 0xff}
	case 4: return color.RGBA{uint8(e*255), uint8(c*255), uint8(z.v*255), 0xff}
	case 5: return color.RGBA{uint8(z.v*255), uint8(c*255), uint8(d*255), 0xff}
	default: return color.RGBA{uint8(z.v*255), uint8(c*255), uint8(d*255), 0xff}
	}
}

func interpolateHSV(c1, c2 HSV, fraction float64) HSV {
	return HSV{ (1. - fraction) * c1.h + fraction * c2.h,
				(1. - fraction) * c1.s + fraction * c2.s,
				(1. - fraction) * c1.v + fraction * c2.v }
}

func interpolateRGBA(c1, c2 color.RGBA, fraction float64) color.RGBA {
	return color.RGBA{ uint8((1. - fraction) * float64(c1.R) + fraction * float64(c2.R)),
		uint8((1. - fraction) * float64(c1.G) + fraction * float64(c2.G)),
		uint8((1. - fraction) * float64(c1.B) + fraction * float64(c2.B)),
		uint8((1. - fraction) * float64(c1.A) + fraction * float64(c2.A)) }
}

// func interpolate(xs []float64, ys []HSV, samples int) []color.RGBA {
// 	palette := make([]color.RGBA, samples)
// 	i := 0
// 	for j := 0; j < len(xs)-1; j++ {
// 		cap := int(float64(samples) * xs[j+1])
// 		for ; i < cap; i++ {
// 			palette[i] = HSVtoRGB(interpolateHSV(ys[j], ys[j+1], (float64(i) / float64(samples) - xs[j]) / (xs[j+1]-xs[j])))
// 		}
// 	}
// 	return palette
// }

func interpolate(xs []float64, ys []color.RGBA, samples int) []color.RGBA {
	palette := make([]color.RGBA, samples)
	i := 0
	for j := 0; j < len(xs)-1; j++ {
		cap := int(float64(samples) * xs[j+1])
		for ; i < cap ; i++ {
			palette[i] = interpolateRGBA(ys[j], ys[j+1], (float64(i) / float64(samples) - xs[j]) / (xs[j+1] - xs[j]))
		}
	}
	return palette
}

func printPalette(filename string, palette []color.RGBA) {
	img := image.NewRGBA(image.Rectangle{image.Point{0, 0}, image.Point{len(palette), 256}})
	red := color.RGBA{255, 0, 0, 0xff}
	green := color.RGBA{0, 255, 0, 0xff}
	blue := color.RGBA{0, 0, 255, 0xff}
	for x := 0; x < len(palette); x++ {
		for y := 0; y < 256; y++ {
			img.Set(x, y, palette[x])
		}
		img.Set(x, int(palette[x].R), red)
		img.Set(x, int(palette[x].G), blue)
		img.Set(x, int(palette[x].B), green)
	}
	f, _ := os.Create("palettes/" + filename + ".png")
	png.Encode(f, img)
}

func min(a, b int) int {
	if a < b { return a }
	return b
}

func main() {
	quintpoints := []float64{0.0, 0.25, 0.5, 0.75, 1.0}
	darkendpoints := []float64{0.0, 0.2, 0.35, 0.55, 0.75, 1.0}
	quadpoints := []float64{0.0, 0.3, 0.6, 1.0}
	_ = quadpoints
	tripoints := []float64{0.0, 0.5, 1.0}
	_ = tripoints
	makocolours := []color.RGBA{
		color.RGBA{221, 244, 228, 0xff},
		color.RGBA{76, 195, 173, 0xff},
		color.RGBA{53, 123, 162, 0xff},
		color.RGBA{61, 50, 102, 0xff},
		color.RGBA{12, 3, 4, 0xff}}

	makocolours2 := []color.RGBA{
		color.RGBA{3, 5, 26, 0xff},
		color.RGBA{76, 195, 173, 0xff},
		color.RGBA{53, 123, 162, 0xff},
		color.RGBA{61, 50, 102, 0xff},
		color.RGBA{12, 3, 4, 0xff}}
	// printPalette("mako2", interpolate(quintpoints, makocolours, 1000))
	// makocolours := []HSV{
	// 	HSV{138., 0.09, 0.96},
	// 	HSV{179., 0.67, 0.68},
	// 	HSV{236, 0.53, 0.54},
	// 	HSV{353., 0.67, 0.04} }
	// _, _ = makopoints, makocolours
	// printPalette("mako", interpolate(makopoints, makocolours, 1000))
	rocketcolours2 := []color.RGBA{
		color.RGBA{4, 5, 27, 0xff},
		color.RGBA{250, 232, 216, 0xff},
		color.RGBA{243, 132, 93, 0xff},
		color.RGBA{203, 30, 78, 0xff},
		color.RGBA{91, 30, 81, 0xff},
		color.RGBA{3, 5, 26, 0xff}}
	rocketcolours := []color.RGBA{
		color.RGBA{250, 232, 216, 0xff},
		color.RGBA{243, 132, 93, 0xff},
		color.RGBA{203, 30, 78, 0xff},
		color.RGBA{91, 30, 81, 0xff},
		color.RGBA{3, 5, 26, 0xff}}

	kovesipoints := []float64{0.0, 0.25, 0.55, 0.85, 1.0}
	kovesipink := []color.RGBA{
		color.RGBA{254, 231, 254, 0xff},
		color.RGBA{249, 96, 253, 0xff},
		color.RGBA{4, 49, 247, 0xff},
		color.RGBA{0, 11, 91, 0xff},
		color.RGBA{0, 6, 32, 0xff}}
	kovesiblue := []color.RGBA{
		color.RGBA{175, 254, 246, 0xff},
		color.RGBA{49, 156, 253, 0xff},
		color.RGBA{31, 47, 240, 0xff},
		color.RGBA{0, 9, 81, 0xff},
		color.RGBA{0, 4, 21, 0xff}}

	var palette []color.RGBA
	if !LINES {
		switch SCHEME {
		case 0:
			palette = interpolate(quintpoints, makocolours, MAXCONV+1)
		case 1:
			palette = interpolate(quintpoints, rocketcolours, MAXCONV+1)
		case 2:
			palette = interpolate(quintpoints, makocolours2, MAXCONV+1)
		case 3:
			palette = interpolate(darkendpoints, rocketcolours2, MAXCONV+1)
		case 4:
			palette = interpolate(kovesipoints, kovesipink, MAXCONV+1)
		case 5:
			palette = interpolate(kovesipoints, kovesiblue, MAXCONV+1)
		}
		// printPalette("purpleToGreen", palette)
	}
	RESTART:
	img := image.NewRGBA(image.Rectangle{image.Point{0, 0}, image.Point{WIDTH, HEIGHT}})
	var density, bitmap [WIDTH+2*CONVWIDTH][HEIGHT+2*CONVWIDTH]int
	var actualdensity [WIDTH][HEIGHT]int
	// cyan := color.RGBA{100, 200, 200, 0xff}
	for x := 0; x < WIDTH; x++ {
		for y := 0; y < HEIGHT; y++ {
			img.Set(x, y, color.White)
		}
	}

	state := initialise((rand.Float64() - 0.5) * 2 * math.Pi,
				 		(rand.Float64() - 0.5) * math.Pi / 2,
						(rand.Float64() - 0.5) * math.Pi,
						(rand.Float64() - 0.5) * math.Pi / 2)
	start := time.Now()
	currx, curry := coordToPixel(polarToCoord(state.v, state.u))
	var prevx, prevy int
	// x, y := polarToCoord(state.v, state.u)
	for i := 0; i <= LENGTH; i++ {
		state.updateState()
		prevx, prevy = currx, curry
		currx, curry = coordToPixel(polarToCoord(state.v, state.u))
		if 0 <= currx && currx < WIDTH && 0 <= curry && curry < HEIGHT {
			if LINES {
				Bresenham(img, prevx, prevy, currx, curry, color.Black)
			} else {
				bitmap[currx+CONVWIDTH][curry+CONVWIDTH] = 1
				// img.Set(prevx, prevy, color.Black)
				for x := -CONVWIDTH; x <= CONVWIDTH; x++ {
					for y := -CONVWIDTH; y <= CONVWIDTH; y++ {
						density[currx+x+CONVWIDTH][curry+y+CONVWIDTH]++
					}
				}
			}
		}
		if time.Since(start) > TIMELIMIT {
			fmt.Println("Took too long, restarting")
			goto RESTART
		}
	}

	if ACTUALDENSITY {
		for x := 0; x < WIDTH; x++ {
			for y := 0; y < HEIGHT; y++ {
				for cx := -CONVWIDTH; cx <= CONVWIDTH; cx++ {
					for cy := -CONVWIDTH; cy <= CONVWIDTH; cy++ {
						actualdensity[x][y] += bitmap[x+cx+CONVWIDTH][y+cy+CONVWIDTH]
					}
				}
			}
		}
	}
	var idx int
	if !LINES {
		for x := 0; x < WIDTH; x++ {
			for y := 0; y < HEIGHT; y++ {
				// fmt.Println(x, y)
				if ACTUALDENSITY {
					idx = actualdensity[x][y]
				} else {
					idx = min(density[x+CONVWIDTH][y+CONVWIDTH], MAXCONV)
				}
				if REVERSE {
					idx = MAXCONV - idx
				}
				img.Set(x, y, palette[idx])
			}
		}
	}

	f, _ := os.Create(NAME)
	png.Encode(f, img)
}

// Not my code - taken from https://github.com/StephaneBunel/bresenham
type Plotter interface { Set(x, y int, c color.Color) }

func Bresenham(p Plotter, x1, y1, x2, y2 int, col color.Color) {
	var dx, dy, e, slope int

	// Because drawing p1 -> p2 is equivalent to draw p2 -> p1,
	// I sort points in x-axis order to handle only half of possible cases.
	if x1 > x2 {
		x1, y1, x2, y2 = x2, y2, x1, y1
	}

	dx, dy = x2-x1, y2-y1
	// Because point is x-axis ordered, dx cannot be negative
	if dy < 0 {
		dy = -dy
	}

	switch {

	// Is line a point ?
	case x1 == x2 && y1 == y2:
		p.Set(x1, y1, col)

	// Is line an horizontal ?
	case y1 == y2:
		for ; dx != 0; dx-- {
			p.Set(x1, y1, col)
			x1++
		}
		p.Set(x1, y1, col)

	// Is line a vertical ?
	case x1 == x2:
		if y1 > y2 {
			y1, y2 = y2, y1
		}
		for ; dy != 0; dy-- {
			p.Set(x1, y1, col)
			y1++
		}
		p.Set(x1, y1, col)

	// Is line a diagonal ?
	case dx == dy:
		if y1 < y2 {
			for ; dx != 0; dx-- {
				p.Set(x1, y1, col)
				x1++
				y1++
			}
		} else {
			for ; dx != 0; dx-- {
				p.Set(x1, y1, col)
				x1++
				y1--
			}
		}
		p.Set(x1, y1, col)

	// wider than high ?
	case dx > dy:
		if y1 < y2 {
			// BresenhamDxXRYD(img, x1, y1, x2, y2, col)
			dy, e, slope = 2*dy, dx, 2*dx
			for ; dx != 0; dx-- {
				p.Set(x1, y1, col)
				x1++
				e -= dy
				if e < 0 {
					y1++
					e += slope
				}
			}
		} else {
			// BresenhamDxXRYU(img, x1, y1, x2, y2, col)
			dy, e, slope = 2*dy, dx, 2*dx
			for ; dx != 0; dx-- {
				p.Set(x1, y1, col)
				x1++
				e -= dy
				if e < 0 {
					y1--
					e += slope
				}
			}
		}
		p.Set(x2, y2, col)

	// higher than wide.
	default:
		if y1 < y2 {
			// BresenhamDyXRYD(img, x1, y1, x2, y2, col)
			dx, e, slope = 2*dx, dy, 2*dy
			for ; dy != 0; dy-- {
				p.Set(x1, y1, col)
				y1++
				e -= dx
				if e < 0 {
					x1++
					e += slope
				}
			}
		} else {
			// BresenhamDyXRYU(img, x1, y1, x2, y2, col)
			dx, e, slope = 2*dx, dy, 2*dy
			for ; dy != 0; dy-- {
				p.Set(x1, y1, col)
				y1--
				e -= dx
				if e < 0 {
					x1++
					e += slope
				}
			}
		}
		p.Set(x2, y2, col)
	}
}

// End
// func interpolate(xs []float64, ys []color.RGBA, samples int) []color.RGBA {
// 	if len(xs) != len(ys) { panic("ill-defined function to interpolate") }
// 	n := len(xs)
// 	delta := make([][4]float64, len(xs)-1)
// 	dxs := make([]float64, len(xs)-1)
// 	dys := make([][4]float64, len(xs)-1)
// 	for k := 0; k < n-1; k++ {
// 		if xs[k] >= xs[k+1] {
// 			panic("xs not sorted")
// 		}
// 		dxs[k] = xs[k+1]-xs[k]
// 		dys[k][0] = float64(ys[k+1].R - ys[k].R)
// 		dys[k][1] = float64(ys[k+1].G - ys[k].G)
// 		dys[k][2] = float64(ys[k+1].B - ys[k].B)
// 		dys[k][3] = float64(ys[k+1].A - ys[k].A)
// 		for c := 0; c < 4; c++ {
// 			delta[k][c] = dys[k][c] / dxs[k]
// 		}
// 	}
// 	m := make([][4]float64, len(xs))
// 	for c := 0; c < 4; c++ {
// 		m[0][c] = delta[0][c]
// 		m[n-1][c] = delta[n-2][c]
// 	}
// 	for k := 1; k < n-1; k++ {
// 		common := dxs[k-1] + dxs[k]
// 		for c := 0; c < 4; c++ {
// 			if delta[k-1][c] * delta[k][c] <= 0 {
// 				m[k][c] = 0
// 			} else {
// 				m[k][c] = 3*common / ((common + dxs[k])/delta[k-1][c] + (common+dxs[k-1])/delta[k][c])
// 			}
// 		}
// 	}
//
// 	xi := make([][4]float64, n-1)
// 	zeta := make([][4]float64, n-1)
// 	for k := 0; k < n-1; k++ {
// 		for c := 0; c < 4; c++ {
// 			xi[k][c] = (3*delta[k][c]-2*m[k][c]-m[k+1][c])/dxs[k]
// 			zeta[k][c] = (m[k][c]+m[k+1][c]-2*delta[k][c])/math.Pow(dxs[k],2)
// 		}
// 	}
//
// 	dest := make([]color.RGBA, samples)
// 	i := 0
// 	var d, d2, d3 float64
// 	var cap int
// 	// h is the highest value such that xs[h] < i/samples
// 	// if xs[0] >= i/samples, set j = 0, otherwise j = h
// 	//
// 	// Goal: go through every i from 0 to samples-1
// 	// At any point of time operate with the most recent point in xs (diff should always be >=0)
// 	for j := 0; j < len(xs)-1; j++ {
// 		cap = int(float64(samples) * xs[j+1])
// 		for ; i < cap; i++ {
// 			d = float64(i) / float64(samples) - xs[j]
// 			d2 = d*d
// 			d3 = d2*d
// 			dest[i].R = uint8(float64(ys[j].R) + m[j][0]*d + xi[j][0]*d2 + zeta[j][0]*d3)
// 			dest[i].G = uint8(float64(ys[j].G) + m[j][1]*d + xi[j][1]*d2 + zeta[j][1]*d3)
// 			dest[i].B = uint8(float64(ys[j].B) + m[j][2]*d + xi[j][2]*d2 + zeta[j][2]*d3)
// 			dest[i].A = uint8(float64(ys[j].A) + m[j][3]*d + xi[j][3]*d2 + zeta[j][3]*d3)
// 		}
// 	}
// 	for ; i < samples; i++ {
// 		d = float64(i) / float64(samples) - xs[n-1]
// 		d2 = d*d
// 		d3 = d2*d
// 		dest[i].R = uint8(float64(ys[n-1].R) + m[n-2][0]*d + xi[n-2][0]*d2 + zeta[n-2][0]*d3)
// 		dest[i].G = uint8(float64(ys[n-1].G) + m[n-2][1]*d + xi[n-2][1]*d2 + zeta[n-2][1]*d3)
// 		dest[i].B = uint8(float64(ys[n-1].B) + m[n-2][2]*d + xi[n-2][2]*d2 + zeta[n-2][2]*d3)
// 		dest[i].A = uint8(float64(ys[n-1].A) + m[n-2][3]*d + xi[n-2][3]*d2 + zeta[n-2][3]*d3)
// 	}
// 	return dest
// }
