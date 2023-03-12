package strconv

import (
	"math"
	"math/bits"
)

const fnParseFloat = "ParseFloat"

// commonPrefixLenIgnoreCase returns the length of the common
// prefix of s and prefix, with the character case of s ignored.
// The prefix argument must be all lower-case.
func commonPrefixLenIgnoreCase(s, prefix string) int {
	n := len(prefix)
	if n > len(s) {
		n = len(s)
	}
	for i := 0; i < n; i++ {
		c := s[i]
		if 'A' <= c && c <= 'Z' {
			c += 'a' - 'A'
		}
		if c != prefix[i] {
			return i
		}
	}
	return n
}

// special returns the floating-point value for the special,
// possibly signed floating-point representations inf, infinity,
// and NaN. The result is ok if a prefix of s contains one
// of these representations and n is the length of that prefix.
// The character case is ignored.
func special(s string) (f float64, n int, ok bool) {
	if len(s) == 0 {
		return 0, 0, false
	}
	sign := 1
	nsign := 0
	switch s[0] {
	case '+', '-':
		if s[0] == '-' {
			sign = -1
		}
		nsign = 1
		s = s[1:]
		fallthrough
	case 'i', 'I':
		n := commonPrefixLenIgnoreCase(s, "infinity")
		// Anything longer than "inf" is ok, but if we
		// don't have "infinity", only consume "inf".
		if 3 < n && n < 8 {
			n = 3
		}
		if n == 3 || n == 8 {
			return math.Inf(sign), nsign + n, true
		}
	case 'n', 'N':
		if commonPrefixLenIgnoreCase(s, "nan") == 3 {
			return math.NaN(), 3, true
		}
	}
	return 0, 0, false
}

// readFloat reads a decimal or hexadecimal mantissa and exponent from a float
// string representation in s; the number may be followed by other characters.
// readFloat reports the number of bytes consumed (i), and whether the number
// is valid (ok).
func readFloat(s string) (mantissa uint64, exp int, neg, trunc, hex bool, i int, ok bool) {
	underscores := false

	// optional sign
	if i >= len(s) {
		return
	}
	switch {
	case s[i] == '+':
		i++
	case s[i] == '-':
		neg = true
		i++
	}

	// digits
	base := uint64(10)
	maxMantDigits := 19 // 10^19 fits in uint64
	expChar := byte('e')
	if i+2 < len(s) && s[i] == '0' && lower(s[i+1]) == 'x' {
		base = 16
		maxMantDigits = 16 // 16^16 fits in uint64
		i += 2
		expChar = 'p'
		hex = true
	}
	sawdot := false
	sawdigits := false
	nd := 0
	ndMant := 0
	dp := 0
loop:
	for ; i < len(s); i++ {
		switch c := s[i]; true {
		case c == '_':
			underscores = true
			continue

		case c == '.':
			if sawdot {
				break loop
			}
			sawdot = true
			dp = nd
			continue

		case '0' <= c && c <= '9':
			sawdigits = true
			if c == '0' && nd == 0 { // ignore leading zeros
				dp--
				continue
			}
			nd++
			if ndMant < maxMantDigits {
				mantissa *= base
				mantissa += uint64(c - '0')
				ndMant++
			} else if c != '0' {
				trunc = true
			}
			continue

		case base == 16 && 'a' <= lower(c) && lower(c) <= 'f':
			sawdigits = true
			nd++
			if ndMant < maxMantDigits {
				mantissa *= 16
				mantissa += uint64(lower(c) - 'a' + 10)
				ndMant++
			} else {
				trunc = true
			}
			continue
		}
		break
	}
	if !sawdigits {
		return
	}
	if !sawdot {
		dp = nd
	}

	if base == 16 {
		dp *= 4
		ndMant *= 4
	}

	// optional exponent moves decimal point.
	// if we read a very large, very long number,
	// just be sure to move the decimal point by
	// a lot (say, 100000).  it doesn't matter if it's
	// not the exact number.
	if i < len(s) && lower(s[i]) == expChar {
		i++
		if i >= len(s) {
			return
		}
		esign := 1
		if s[i] == '+' {
			i++
		} else if s[i] == '-' {
			i++
			esign = -1
		}
		if i >= len(s) || s[i] < '0' || s[i] > '9' {
			return
		}
		e := 0
		for ; i < len(s) && ('0' <= s[i] && s[i] <= '9' || s[i] == '_'); i++ {
			if s[i] == '_' {
				underscores = true
				continue
			}
			if e < 10000 {
				e = e*10 + int(s[i]) - '0'
			}
		}
		dp += e * esign
	} else if base == 16 {
		// Must have exponent.
		return
	}

	if mantissa != 0 {
		exp = dp - ndMant
	}

	if underscores && !underscoreOK(s[:i]) {
		return
	}

	ok = true
	return
}

// TODO: move elsewhere?
type floatInfo struct {
	mantbits uint
	expbits  uint
	bias     int
}

var float32info = floatInfo{23, 8, -127}
var float64info = floatInfo{52, 11, -1023}

var optimize = true // set to false to force slow-path conversions for testing

// atofHex converts the hex floating-point string s
// to a rounded float32 or float64 value (depending on flt==&float32info or flt==&float64info)
// and returns it as a float64.
// The string s has already been parsed into a mantissa, exponent, and sign (neg==true for negative).
// If trunc is true, trailing non-zero bits have been omitted from the mantissa.
func atofHex(s string, flt *floatInfo, mantissa uint64, exp int, neg, trunc bool) (float64, error) {
	maxExp := 1<<flt.expbits + flt.bias - 2
	minExp := flt.bias + 1
	exp += int(flt.mantbits) // mantissa now implicitly divided by 2^mantbits.

	// Shift mantissa and exponent to bring representation into float range.
	// Eventually we want a mantissa with a leading 1-bit followed by mantbits other bits.
	// For rounding, we need two more, where the bottom bit represents
	// whether that bit or any later bit was non-zero.
	// (If the mantissa has already lost non-zero bits, trunc is true,
	// and we OR in a 1 below after shifting left appropriately.)
	for mantissa != 0 && mantissa>>(flt.mantbits+2) == 0 {
		mantissa <<= 1
		exp--
	}
	if trunc {
		mantissa |= 1
	}
	for mantissa>>(1+flt.mantbits+2) != 0 {
		mantissa = mantissa>>1 | mantissa&1
		exp++
	}

	// If exponent is too negative,
	// denormalize in hopes of making it representable.
	// (The -2 is for the rounding bits.)
	for mantissa > 1 && exp < minExp-2 {
		mantissa = mantissa>>1 | mantissa&1
		exp++
	}

	// Round using two bottom bits.
	round := mantissa & 3
	mantissa >>= 2
	round |= mantissa & 1 // round to even (round up if mantissa is odd)
	exp += 2
	if round == 3 {
		mantissa++
		if mantissa == 1<<(1+flt.mantbits) {
			mantissa >>= 1
			exp++
		}
	}

	if mantissa>>flt.mantbits == 0 { // Denormal or zero.
		exp = flt.bias
	}
	var err error
	if exp > maxExp { // infinity and range error
		mantissa = 1 << flt.mantbits
		exp = maxExp + 1
		err = rangeError(fnParseFloat, s)
	}

	bits := mantissa & (1<<flt.mantbits - 1)
	bits |= uint64((exp-flt.bias)&(1<<flt.expbits-1)) << flt.mantbits
	if neg {
		bits |= 1 << flt.mantbits << flt.expbits
	}
	if flt == &float32info {
		return float64(math.Float32frombits(uint32(bits))), err
	}
	return math.Float64frombits(bits), err
}

// Exact powers of 10.
var float64pow10 = []float64{
	1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,
	1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
	1e20, 1e21, 1e22,
}
var float32pow10 = []float32{1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10}

// If possible to convert decimal representation to 64-bit float f exactly,
// entirely in floating-point math, do so, avoiding the expense of decimalToFloatBits.
// Three common cases:
//
//	value is exact integer
//	value is exact integer * exact power of ten
//	value is exact integer / exact power of ten
//
// These all produce potentially inexact but correctly rounded answers.
func atof64exact(mantissa uint64, exp int, neg bool) (f float64, ok bool) {
	if mantissa>>float64info.mantbits != 0 {
		return
	}
	f = float64(mantissa)
	if neg {
		f = -f
	}
	switch {
	case exp == 0:
		// an integer.
		return f, true
	// Exact integers are <= 10^15.
	// Exact powers of ten are <= 10^22.
	case exp > 0 && exp <= 15+22: // int * 10^k
		// If exponent is big but number of digits is not,
		// can move a few zeros into the integer part.
		if exp > 22 {
			f *= float64pow10[exp-22]
			exp = 22
		}
		if f > 1e15 || f < -1e15 {
			// the exponent was really too large.
			return
		}
		return f * float64pow10[exp], true
	case exp < 0 && exp >= -22: // int / 10^k
		return f / float64pow10[-exp], true
	}
	return
}

// If possible to compute mantissa*10^exp to 32-bit float f exactly,
// entirely in floating-point math, do so, avoiding the machinery above.
func atof32exact(mantissa uint64, exp int, neg bool) (f float32, ok bool) {
	if mantissa>>float32info.mantbits != 0 {
		return
	}
	f = float32(mantissa)
	if neg {
		f = -f
	}
	switch {
	case exp == 0:
		return f, true
	// Exact integers are <= 10^7.
	// Exact powers of ten are <= 10^10.
	case exp > 0 && exp <= 7+10: // int * 10^k
		// If exponent is big but number of digits is not,
		// can move a few zeros into the integer part.
		if exp > 10 {
			f *= float32pow10[exp-10]
			exp = 10
		}
		if f > 1e7 || f < -1e7 {
			// the exponent was really too large.
			return
		}
		return f * float32pow10[exp], true
	case exp < 0 && exp >= -10: // int / 10^k
		return f / float32pow10[-exp], true
	}
	return
}

func eiselLemire64(man uint64, exp10 int, neg bool) (f float64, ok bool) {
	// The terse comments in this function body refer to sections of the
	// https://nigeltao.github.io/blog/2020/eisel-lemire.html blog post.

	// Exp10 Range.
	if man == 0 {
		if neg {
			f = math.Float64frombits(0x8000000000000000) // Negative zero.
		}
		return f, true
	}
	if exp10 < detailedPowersOfTenMinExp10 || detailedPowersOfTenMaxExp10 < exp10 {
		return 0, false
	}

	// Normalization.
	clz := bits.LeadingZeros64(man)
	man <<= uint(clz)
	const float64ExponentBias = 1023
	retExp2 := uint64(217706*exp10>>16+64+float64ExponentBias) - uint64(clz)

	// Multiplication.
	xHi, xLo := bits.Mul64(man, detailedPowersOfTen[exp10-detailedPowersOfTenMinExp10][1])

	// Wider Approximation.
	if xHi&0x1FF == 0x1FF && xLo+man < man {
		yHi, yLo := bits.Mul64(man, detailedPowersOfTen[exp10-detailedPowersOfTenMinExp10][0])
		mergedHi, mergedLo := xHi, xLo+yHi
		if mergedLo < xLo {
			mergedHi++
		}
		if mergedHi&0x1FF == 0x1FF && mergedLo+1 == 0 && yLo+man < man {
			return 0, false
		}
		xHi, xLo = mergedHi, mergedLo
	}

	// Shifting to 54 Bits.
	msb := xHi >> 63
	retMantissa := xHi >> (msb + 9)
	retExp2 -= 1 ^ msb

	// Half-way Ambiguity.
	if xLo == 0 && xHi&0x1FF == 0 && retMantissa&3 == 1 {
		return 0, false
	}

	// From 54 to 53 Bits.
	retMantissa += retMantissa & 1
	retMantissa >>= 1
	if retMantissa>>53 > 0 {
		retMantissa >>= 1
		retExp2 += 1
	}
	// retExp2 is a uint64. Zero or underflow means that we're in subnormal
	// float64 space. 0x7FF or above means that we're in Inf/NaN float64 space.
	//
	// The if block is equivalent to (but has fewer branches than):
	//   if retExp2 <= 0 || retExp2 >= 0x7FF { etc }
	if retExp2-1 >= 0x7FF-1 {
		return 0, false
	}
	retBits := retExp2<<52 | retMantissa&0x000FFFFFFFFFFFFF
	if neg {
		retBits |= 0x8000000000000000
	}
	return math.Float64frombits(retBits), true
}

// detailedPowersOfTen{Min,Max}Exp10 is the power of 10 represented by the
// first and last rows of detailedPowersOfTen. Both bounds are inclusive.
const (
	detailedPowersOfTenMinExp10 = -348
	detailedPowersOfTenMaxExp10 = +347
)

func eiselLemire32(man uint64, exp10 int, neg bool) (f float32, ok bool) {
	// The terse comments in this function body refer to sections of the
	// https://nigeltao.github.io/blog/2020/eisel-lemire.html blog post.
	//
	// That blog post discusses the float64 flavor (11 exponent bits with a
	// -1023 bias, 52 mantissa bits) of the algorithm, but the same approach
	// applies to the float32 flavor (8 exponent bits with a -127 bias, 23
	// mantissa bits). The computation here happens with 64-bit values (e.g.
	// man, xHi, retMantissa) before finally converting to a 32-bit float.

	// Exp10 Range.
	if man == 0 {
		if neg {
			f = math.Float32frombits(0x80000000) // Negative zero.
		}
		return f, true
	}
	if exp10 < detailedPowersOfTenMinExp10 || detailedPowersOfTenMaxExp10 < exp10 {
		return 0, false
	}

	// Normalization.
	clz := bits.LeadingZeros64(man)
	man <<= uint(clz)
	const float32ExponentBias = 127
	retExp2 := uint64(217706*exp10>>16+64+float32ExponentBias) - uint64(clz)

	// Multiplication.
	xHi, xLo := bits.Mul64(man, detailedPowersOfTen[exp10-detailedPowersOfTenMinExp10][1])

	// Wider Approximation.
	if xHi&0x3FFFFFFFFF == 0x3FFFFFFFFF && xLo+man < man {
		yHi, yLo := bits.Mul64(man, detailedPowersOfTen[exp10-detailedPowersOfTenMinExp10][0])
		mergedHi, mergedLo := xHi, xLo+yHi
		if mergedLo < xLo {
			mergedHi++
		}
		if mergedHi&0x3FFFFFFFFF == 0x3FFFFFFFFF && mergedLo+1 == 0 && yLo+man < man {
			return 0, false
		}
		xHi, xLo = mergedHi, mergedLo
	}

	// Shifting to 54 Bits (and for float32, it's shifting to 25 bits).
	msb := xHi >> 63
	retMantissa := xHi >> (msb + 38)
	retExp2 -= 1 ^ msb

	// Half-way Ambiguity.
	if xLo == 0 && xHi&0x3FFFFFFFFF == 0 && retMantissa&3 == 1 {
		return 0, false
	}

	// From 54 to 53 Bits (and for float32, it's from 25 to 24 bits).
	retMantissa += retMantissa & 1
	retMantissa >>= 1
	if retMantissa>>24 > 0 {
		retMantissa >>= 1
		retExp2 += 1
	}
	// retExp2 is a uint64. Zero or underflow means that we're in subnormal
	// float32 space. 0xFF or above means that we're in Inf/NaN float32 space.
	//
	// The if block is equivalent to (but has fewer branches than):
	//   if retExp2 <= 0 || retExp2 >= 0xFF { etc }
	if retExp2-1 >= 0xFF-1 {
		return 0, false
	}
	retBits := retExp2<<23 | retMantissa&0x007FFFFF
	if neg {
		retBits |= 0x80000000
	}
	return math.Float32frombits(uint32(retBits)), true
}

type decimal struct {
	d     [800]byte // digits, big-endian representation
	nd    int       // number of digits used
	dp    int       // decimal point
	neg   bool      // negative flag
	trunc bool      // discarded nonzero digits beyond d[:nd]
}

func (b *decimal) set(s string) (ok bool) {
	i := 0
	b.neg = false
	b.trunc = false

	// optional sign
	if i >= len(s) {
		return
	}
	switch {
	case s[i] == '+':
		i++
	case s[i] == '-':
		b.neg = true
		i++
	}

	// digits
	sawdot := false
	sawdigits := false
	for ; i < len(s); i++ {
		switch {
		case s[i] == '_':
			// readFloat already checked underscores
			continue
		case s[i] == '.':
			if sawdot {
				return
			}
			sawdot = true
			b.dp = b.nd
			continue

		case '0' <= s[i] && s[i] <= '9':
			sawdigits = true
			if s[i] == '0' && b.nd == 0 { // ignore leading zeros
				b.dp--
				continue
			}
			if b.nd < len(b.d) {
				b.d[b.nd] = s[i]
				b.nd++
			} else if s[i] != '0' {
				b.trunc = true
			}
			continue
		}
		break
	}
	if !sawdigits {
		return
	}
	if !sawdot {
		b.dp = b.nd
	}

	// optional exponent moves decimal point.
	// if we read a very large, very long number,
	// just be sure to move the decimal point by
	// a lot (say, 100000).  it doesn't matter if it's
	// not the exact number.
	if i < len(s) && lower(s[i]) == 'e' {
		i++
		if i >= len(s) {
			return
		}
		esign := 1
		if s[i] == '+' {
			i++
		} else if s[i] == '-' {
			i++
			esign = -1
		}
		if i >= len(s) || s[i] < '0' || s[i] > '9' {
			return
		}
		e := 0
		for ; i < len(s) && ('0' <= s[i] && s[i] <= '9' || s[i] == '_'); i++ {
			if s[i] == '_' {
				// readFloat already checked underscores
				continue
			}
			if e < 10000 {
				e = e*10 + int(s[i]) - '0'
			}
		}
		b.dp += e * esign
	}

	if i != len(s) {
		return
	}

	ok = true
	return
}

type leftCheat struct {
	delta  int    // number of new digits
	cutoff string // minus one digit if original < a.
}

var leftcheats = []leftCheat{
	// Leading digits of 1/2^i = 5^i.
	// 5^23 is not an exact 64-bit floating point number,
	// so have to use bc for the math.
	// Go up to 60 to be large enough for 32bit and 64bit platforms.
	/*
		seq 60 | sed 's/^/5^/' | bc |
		awk 'BEGIN{ print "\t{ 0, \"\" }," }
		{
			log2 = log(2)/log(10)
			printf("\t{ %d, \"%s\" },\t// * %d\n",
				int(log2*NR+1), $0, 2**NR)
		}'
	*/
	{0, ""},
	{1, "5"},                                           // * 2
	{1, "25"},                                          // * 4
	{1, "125"},                                         // * 8
	{2, "625"},                                         // * 16
	{2, "3125"},                                        // * 32
	{2, "15625"},                                       // * 64
	{3, "78125"},                                       // * 128
	{3, "390625"},                                      // * 256
	{3, "1953125"},                                     // * 512
	{4, "9765625"},                                     // * 1024
	{4, "48828125"},                                    // * 2048
	{4, "244140625"},                                   // * 4096
	{4, "1220703125"},                                  // * 8192
	{5, "6103515625"},                                  // * 16384
	{5, "30517578125"},                                 // * 32768
	{5, "152587890625"},                                // * 65536
	{6, "762939453125"},                                // * 131072
	{6, "3814697265625"},                               // * 262144
	{6, "19073486328125"},                              // * 524288
	{7, "95367431640625"},                              // * 1048576
	{7, "476837158203125"},                             // * 2097152
	{7, "2384185791015625"},                            // * 4194304
	{7, "11920928955078125"},                           // * 8388608
	{8, "59604644775390625"},                           // * 16777216
	{8, "298023223876953125"},                          // * 33554432
	{8, "1490116119384765625"},                         // * 67108864
	{9, "7450580596923828125"},                         // * 134217728
	{9, "37252902984619140625"},                        // * 268435456
	{9, "186264514923095703125"},                       // * 536870912
	{10, "931322574615478515625"},                      // * 1073741824
	{10, "4656612873077392578125"},                     // * 2147483648
	{10, "23283064365386962890625"},                    // * 4294967296
	{10, "116415321826934814453125"},                   // * 8589934592
	{11, "582076609134674072265625"},                   // * 17179869184
	{11, "2910383045673370361328125"},                  // * 34359738368
	{11, "14551915228366851806640625"},                 // * 68719476736
	{12, "72759576141834259033203125"},                 // * 137438953472
	{12, "363797880709171295166015625"},                // * 274877906944
	{12, "1818989403545856475830078125"},               // * 549755813888
	{13, "9094947017729282379150390625"},               // * 1099511627776
	{13, "45474735088646411895751953125"},              // * 2199023255552
	{13, "227373675443232059478759765625"},             // * 4398046511104
	{13, "1136868377216160297393798828125"},            // * 8796093022208
	{14, "5684341886080801486968994140625"},            // * 17592186044416
	{14, "28421709430404007434844970703125"},           // * 35184372088832
	{14, "142108547152020037174224853515625"},          // * 70368744177664
	{15, "710542735760100185871124267578125"},          // * 140737488355328
	{15, "3552713678800500929355621337890625"},         // * 281474976710656
	{15, "17763568394002504646778106689453125"},        // * 562949953421312
	{16, "88817841970012523233890533447265625"},        // * 1125899906842624
	{16, "444089209850062616169452667236328125"},       // * 2251799813685248
	{16, "2220446049250313080847263336181640625"},      // * 4503599627370496
	{16, "11102230246251565404236316680908203125"},     // * 9007199254740992
	{17, "55511151231257827021181583404541015625"},     // * 18014398509481984
	{17, "277555756156289135105907917022705078125"},    // * 36028797018963968
	{17, "1387778780781445675529539585113525390625"},   // * 72057594037927936
	{18, "6938893903907228377647697925567626953125"},   // * 144115188075855872
	{18, "34694469519536141888238489627838134765625"},  // * 288230376151711744
	{18, "173472347597680709441192448139190673828125"}, // * 576460752303423488
	{19, "867361737988403547205962240695953369140625"}, // * 1152921504606846976
}

// Is the leading prefix of b lexicographically less than s?
func prefixIsLessThan(b []byte, s string) bool {
	for i := 0; i < len(s); i++ {
		if i >= len(b) {
			return true
		}
		if b[i] != s[i] {
			return b[i] < s[i]
		}
	}
	return false
}

// Binary shift left (* 2) by k bits.  k <= maxShift to avoid overflow.
func leftShift(a *decimal, k uint) {
	delta := leftcheats[k].delta
	if prefixIsLessThan(a.d[0:a.nd], leftcheats[k].cutoff) {
		delta--
	}

	r := a.nd         // read index
	w := a.nd + delta // write index

	// Pick up a digit, put down a digit.
	var n uint
	for r--; r >= 0; r-- {
		n += (uint(a.d[r]) - '0') << k
		quo := n / 10
		rem := n - 10*quo
		w--
		if w < len(a.d) {
			a.d[w] = byte(rem + '0')
		} else if rem != 0 {
			a.trunc = true
		}
		n = quo
	}

	// Put down extra digits.
	for n > 0 {
		quo := n / 10
		rem := n - 10*quo
		w--
		if w < len(a.d) {
			a.d[w] = byte(rem + '0')
		} else if rem != 0 {
			a.trunc = true
		}
		n = quo
	}

	a.nd += delta
	if a.nd >= len(a.d) {
		a.nd = len(a.d)
	}
	a.dp += delta
	trim(a)
}

// Binary shift right (/ 2) by k bits.  k <= maxShift to avoid overflow.
func rightShift(a *decimal, k uint) {
	r := 0 // read pointer
	w := 0 // write pointer

	// Pick up enough leading digits to cover first shift.
	var n uint
	for ; n>>k == 0; r++ {
		if r >= a.nd {
			if n == 0 {
				// a == 0; shouldn't get here, but handle anyway.
				a.nd = 0
				return
			}
			for n>>k == 0 {
				n = n * 10
				r++
			}
			break
		}
		c := uint(a.d[r])
		n = n*10 + c - '0'
	}
	a.dp -= r - 1

	var mask uint = (1 << k) - 1

	// Pick up a digit, put down a digit.
	for ; r < a.nd; r++ {
		c := uint(a.d[r])
		dig := n >> k
		n &= mask
		a.d[w] = byte(dig + '0')
		w++
		n = n*10 + c - '0'
	}

	// Put down extra digits.
	for n > 0 {
		dig := n >> k
		n &= mask
		if w < len(a.d) {
			a.d[w] = byte(dig + '0')
			w++
		} else if dig > 0 {
			a.trunc = true
		}
		n = n * 10
	}

	a.nd = w
	trim(a)
}

// trim trailing zeros from number.
// (They are meaningless; the decimal point is tracked
// independent of the number of digits.)
func trim(a *decimal) {
	for a.nd > 0 && a.d[a.nd-1] == '0' {
		a.nd--
	}
	if a.nd == 0 {
		a.dp = 0
	}
}

const uintSize = 32 << (^uint(0) >> 63)
const maxShift = uintSize - 4

// Binary shift left (k > 0) or right (k < 0).
func (a *decimal) Shift(k int) {
	switch {
	case a.nd == 0:
		// nothing to do: a == 0
	case k > 0:
		for k > maxShift {
			leftShift(a, maxShift)
			k -= maxShift
		}
		leftShift(a, uint(k))
	case k < 0:
		for k < -maxShift {
			rightShift(a, maxShift)
			k += maxShift
		}
		rightShift(a, uint(-k))
	}
}

// decimal power of ten to binary power of two.
var powtab = []int{1, 3, 6, 9, 13, 16, 19, 23, 26}

func (d *decimal) floatBits(flt *floatInfo) (b uint64, overflow bool) {
	var exp int
	var mant uint64

	// Zero is always a special case.
	if d.nd == 0 {
		mant = 0
		exp = flt.bias
		goto out
	}

	// Obvious overflow/underflow.
	// These bounds are for 64-bit floats.
	// Will have to change if we want to support 80-bit floats in the future.
	if d.dp > 310 {
		goto overflow
	}
	if d.dp < -330 {
		// zero
		mant = 0
		exp = flt.bias
		goto out
	}

	// Scale by powers of two until in range [0.5, 1.0)
	exp = 0
	for d.dp > 0 {
		var n int
		if d.dp >= len(powtab) {
			n = 27
		} else {
			n = powtab[d.dp]
		}
		d.Shift(-n)
		exp += n
	}
	for d.dp < 0 || d.dp == 0 && d.d[0] < '5' {
		var n int
		if -d.dp >= len(powtab) {
			n = 27
		} else {
			n = powtab[-d.dp]
		}
		d.Shift(n)
		exp -= n
	}

	// Our range is [0.5,1) but floating point range is [1,2).
	exp--

	// Minimum representable exponent is flt.bias+1.
	// If the exponent is smaller, move it up and
	// adjust d accordingly.
	if exp < flt.bias+1 {
		n := flt.bias + 1 - exp
		d.Shift(-n)
		exp += n
	}

	if exp-flt.bias >= 1<<flt.expbits-1 {
		goto overflow
	}

	// Extract 1+flt.mantbits bits.
	d.Shift(int(1 + flt.mantbits))
	mant = d.RoundedInteger()

	// Rounding might have added a bit; shift down.
	if mant == 2<<flt.mantbits {
		mant >>= 1
		exp++
		if exp-flt.bias >= 1<<flt.expbits-1 {
			goto overflow
		}
	}

	// Denormalized?
	if mant&(1<<flt.mantbits) == 0 {
		exp = flt.bias
	}
	goto out

overflow:
	// ±Inf
	mant = 0
	exp = 1<<flt.expbits - 1 + flt.bias
	overflow = true

out:
	// Assemble bits.
	bits := mant & (uint64(1)<<flt.mantbits - 1)
	bits |= uint64((exp-flt.bias)&(1<<flt.expbits-1)) << flt.mantbits
	if d.neg {
		bits |= 1 << flt.mantbits << flt.expbits
	}
	return bits, overflow
}

// Extract integer part, rounded appropriately.
// No guarantees about overflow.
func (a *decimal) RoundedInteger() uint64 {
	if a.dp > 20 {
		return 0xFFFFFFFFFFFFFFFF
	}
	var i int
	n := uint64(0)
	for i = 0; i < a.dp && i < a.nd; i++ {
		n = n*10 + uint64(a.d[i]-'0')
	}
	for ; i < a.dp; i++ {
		n *= 10
	}
	if shouldRoundUp(a, a.dp) {
		n++
	}
	return n
}

// If we chop a at nd digits, should we round up?
func shouldRoundUp(a *decimal, nd int) bool {
	if nd < 0 || nd >= a.nd {
		return false
	}
	if a.d[nd] == '5' && nd+1 == a.nd { // exactly halfway - round to even
		// if we truncated, a little higher than what's recorded - always round up
		if a.trunc {
			return true
		}
		return nd > 0 && (a.d[nd-1]-'0')%2 != 0
	}
	// not halfway - digit tells all
	return a.d[nd] >= '5'
}

func atof32(s string) (f float32, n int, err error) {
	if val, n, ok := special(s); ok {
		return float32(val), n, nil
	}

	mantissa, exp, neg, trunc, hex, n, ok := readFloat(s)
	if !ok {
		return 0, n, syntaxError(fnParseFloat, s)
	}

	if hex {
		f, err := atofHex(s[:n], &float32info, mantissa, exp, neg, trunc)
		return float32(f), n, err
	}

	if optimize {
		// Try pure floating-point arithmetic conversion, and if that fails,
		// the Eisel-Lemire algorithm.
		if !trunc {
			if f, ok := atof32exact(mantissa, exp, neg); ok {
				return f, n, nil
			}
		}
		f, ok := eiselLemire32(mantissa, exp, neg)
		if ok {
			if !trunc {
				return f, n, nil
			}
			// Even if the mantissa was truncated, we may
			// have found the correct result. Confirm by
			// converting the upper mantissa bound.
			fUp, ok := eiselLemire32(mantissa+1, exp, neg)
			if ok && f == fUp {
				return f, n, nil
			}
		}
	}

	// Slow fallback.
	var d decimal
	if !d.set(s[:n]) {
		return 0, n, syntaxError(fnParseFloat, s)
	}
	b, ovf := d.floatBits(&float32info)
	f = math.Float32frombits(uint32(b))
	if ovf {
		err = rangeError(fnParseFloat, s)
	}
	return f, n, err
}

func atof64(s string) (f float64, n int, err error) {
	if val, n, ok := special(s); ok {
		return val, n, nil
	}

	mantissa, exp, neg, trunc, hex, n, ok := readFloat(s)
	if !ok {
		return 0, n, syntaxError(fnParseFloat, s)
	}

	if hex {
		f, err := atofHex(s[:n], &float64info, mantissa, exp, neg, trunc)
		return f, n, err
	}

	if optimize {
		// Try pure floating-point arithmetic conversion, and if that fails,
		// the Eisel-Lemire algorithm.
		if !trunc {
			if f, ok := atof64exact(mantissa, exp, neg); ok {
				return f, n, nil
			}
		}
		f, ok := eiselLemire64(mantissa, exp, neg)
		if ok {
			if !trunc {
				return f, n, nil
			}
			// Even if the mantissa was truncated, we may
			// have found the correct result. Confirm by
			// converting the upper mantissa bound.
			fUp, ok := eiselLemire64(mantissa+1, exp, neg)
			if ok && f == fUp {
				return f, n, nil
			}
		}
	}

	// Slow fallback.
	var d decimal
	if !d.set(s[:n]) {
		return 0, n, syntaxError(fnParseFloat, s)
	}
	b, ovf := d.floatBits(&float64info)
	f = math.Float64frombits(b)
	if ovf {
		err = rangeError(fnParseFloat, s)
	}
	return f, n, err
}

// ParseFloat converts the string s to a floating-point number
// with the precision specified by bitSize: 32 for float32, or 64 for float64.
// When bitSize=32, the result still has type float64, but it will be
// convertible to float32 without changing its value.
//
// ParseFloat accepts decimal and hexadecimal floating-point numbers
// as defined by the Go syntax for [floating-point literals].
// If s is well-formed and near a valid floating-point number,
// ParseFloat returns the nearest floating-point number rounded
// using IEEE754 unbiased rounding.
// (Parsing a hexadecimal floating-point value only rounds when
// there are more bits in the hexadecimal representation than
// will fit in the mantissa.)
//
// The errors that ParseFloat returns have concrete type *NumError
// and include err.Num = s.
//
// If s is not syntactically well-formed, ParseFloat returns err.Err = ErrSyntax.
//
// If s is syntactically well-formed but is more than 1/2 ULP
// away from the largest floating point number of the given size,
// ParseFloat returns f = ±Inf, err.Err = ErrRange.
//
// ParseFloat recognizes the string "NaN", and the (possibly signed) strings "Inf" and "Infinity"
// as their respective special floating point values. It ignores case when matching.
//
// [floating-point literals]: https://go.dev/ref/spec#Floating-point_literals
func ParseFloat(s string, bitSize int) (float64, error) {
	f, n, err := parseFloatPrefix(s, bitSize)
	if n != len(s) && (err == nil || err.(*NumError).Err != ErrSyntax) {
		return 0, syntaxError(fnParseFloat, s)
	}
	return f, err
}

func parseFloatPrefix(s string, bitSize int) (float64, int, error) {
	if bitSize == 32 {
		f, n, err := atof32(s)
		return float64(f), n, err
	}
	return atof64(s)
}

// detailedPowersOfTen contains 128-bit mantissa approximations (rounded down)
// to the powers of 10. For example:
//
//   - 1e43 ≈ (0xE596B7B0_C643C719                   * (2 ** 79))
//   - 1e43 = (0xE596B7B0_C643C719_6D9CCD05_D0000000 * (2 ** 15))
//
// The mantissas are explicitly listed. The exponents are implied by a linear
// expression with slope 217706.0/65536.0 ≈ log(10)/log(2).
//
// The table was generated by
// https://github.com/google/wuffs/blob/ba3818cb6b473a2ed0b38ecfc07dbbd3a97e8ae7/script/print-mpb-powers-of-10.go
var detailedPowersOfTen = [...][2]uint64{
	{0x1732C869CD60E453, 0xFA8FD5A0081C0288}, // 1e-348
	{0x0E7FBD42205C8EB4, 0x9C99E58405118195}, // 1e-347
	{0x521FAC92A873B261, 0xC3C05EE50655E1FA}, // 1e-346
	{0xE6A797B752909EF9, 0xF4B0769E47EB5A78}, // 1e-345
	{0x9028BED2939A635C, 0x98EE4A22ECF3188B}, // 1e-344
	{0x7432EE873880FC33, 0xBF29DCABA82FDEAE}, // 1e-343
	{0x113FAA2906A13B3F, 0xEEF453D6923BD65A}, // 1e-342
	{0x4AC7CA59A424C507, 0x9558B4661B6565F8}, // 1e-341
	{0x5D79BCF00D2DF649, 0xBAAEE17FA23EBF76}, // 1e-340
	{0xF4D82C2C107973DC, 0xE95A99DF8ACE6F53}, // 1e-339
	{0x79071B9B8A4BE869, 0x91D8A02BB6C10594}, // 1e-338
	{0x9748E2826CDEE284, 0xB64EC836A47146F9}, // 1e-337
	{0xFD1B1B2308169B25, 0xE3E27A444D8D98B7}, // 1e-336
	{0xFE30F0F5E50E20F7, 0x8E6D8C6AB0787F72}, // 1e-335
	{0xBDBD2D335E51A935, 0xB208EF855C969F4F}, // 1e-334
	{0xAD2C788035E61382, 0xDE8B2B66B3BC4723}, // 1e-333
	{0x4C3BCB5021AFCC31, 0x8B16FB203055AC76}, // 1e-332
	{0xDF4ABE242A1BBF3D, 0xADDCB9E83C6B1793}, // 1e-331
	{0xD71D6DAD34A2AF0D, 0xD953E8624B85DD78}, // 1e-330
	{0x8672648C40E5AD68, 0x87D4713D6F33AA6B}, // 1e-329
	{0x680EFDAF511F18C2, 0xA9C98D8CCB009506}, // 1e-328
	{0x0212BD1B2566DEF2, 0xD43BF0EFFDC0BA48}, // 1e-327
	{0x014BB630F7604B57, 0x84A57695FE98746D}, // 1e-326
	{0x419EA3BD35385E2D, 0xA5CED43B7E3E9188}, // 1e-325
	{0x52064CAC828675B9, 0xCF42894A5DCE35EA}, // 1e-324
	{0x7343EFEBD1940993, 0x818995CE7AA0E1B2}, // 1e-323
	{0x1014EBE6C5F90BF8, 0xA1EBFB4219491A1F}, // 1e-322
	{0xD41A26E077774EF6, 0xCA66FA129F9B60A6}, // 1e-321
	{0x8920B098955522B4, 0xFD00B897478238D0}, // 1e-320
	{0x55B46E5F5D5535B0, 0x9E20735E8CB16382}, // 1e-319
	{0xEB2189F734AA831D, 0xC5A890362FDDBC62}, // 1e-318
	{0xA5E9EC7501D523E4, 0xF712B443BBD52B7B}, // 1e-317
	{0x47B233C92125366E, 0x9A6BB0AA55653B2D}, // 1e-316
	{0x999EC0BB696E840A, 0xC1069CD4EABE89F8}, // 1e-315
	{0xC00670EA43CA250D, 0xF148440A256E2C76}, // 1e-314
	{0x380406926A5E5728, 0x96CD2A865764DBCA}, // 1e-313
	{0xC605083704F5ECF2, 0xBC807527ED3E12BC}, // 1e-312
	{0xF7864A44C633682E, 0xEBA09271E88D976B}, // 1e-311
	{0x7AB3EE6AFBE0211D, 0x93445B8731587EA3}, // 1e-310
	{0x5960EA05BAD82964, 0xB8157268FDAE9E4C}, // 1e-309
	{0x6FB92487298E33BD, 0xE61ACF033D1A45DF}, // 1e-308
	{0xA5D3B6D479F8E056, 0x8FD0C16206306BAB}, // 1e-307
	{0x8F48A4899877186C, 0xB3C4F1BA87BC8696}, // 1e-306
	{0x331ACDABFE94DE87, 0xE0B62E2929ABA83C}, // 1e-305
	{0x9FF0C08B7F1D0B14, 0x8C71DCD9BA0B4925}, // 1e-304
	{0x07ECF0AE5EE44DD9, 0xAF8E5410288E1B6F}, // 1e-303
	{0xC9E82CD9F69D6150, 0xDB71E91432B1A24A}, // 1e-302
	{0xBE311C083A225CD2, 0x892731AC9FAF056E}, // 1e-301
	{0x6DBD630A48AAF406, 0xAB70FE17C79AC6CA}, // 1e-300
	{0x092CBBCCDAD5B108, 0xD64D3D9DB981787D}, // 1e-299
	{0x25BBF56008C58EA5, 0x85F0468293F0EB4E}, // 1e-298
	{0xAF2AF2B80AF6F24E, 0xA76C582338ED2621}, // 1e-297
	{0x1AF5AF660DB4AEE1, 0xD1476E2C07286FAA}, // 1e-296
	{0x50D98D9FC890ED4D, 0x82CCA4DB847945CA}, // 1e-295
	{0xE50FF107BAB528A0, 0xA37FCE126597973C}, // 1e-294
	{0x1E53ED49A96272C8, 0xCC5FC196FEFD7D0C}, // 1e-293
	{0x25E8E89C13BB0F7A, 0xFF77B1FCBEBCDC4F}, // 1e-292
	{0x77B191618C54E9AC, 0x9FAACF3DF73609B1}, // 1e-291
	{0xD59DF5B9EF6A2417, 0xC795830D75038C1D}, // 1e-290
	{0x4B0573286B44AD1D, 0xF97AE3D0D2446F25}, // 1e-289
	{0x4EE367F9430AEC32, 0x9BECCE62836AC577}, // 1e-288
	{0x229C41F793CDA73F, 0xC2E801FB244576D5}, // 1e-287
	{0x6B43527578C1110F, 0xF3A20279ED56D48A}, // 1e-286
	{0x830A13896B78AAA9, 0x9845418C345644D6}, // 1e-285
	{0x23CC986BC656D553, 0xBE5691EF416BD60C}, // 1e-284
	{0x2CBFBE86B7EC8AA8, 0xEDEC366B11C6CB8F}, // 1e-283
	{0x7BF7D71432F3D6A9, 0x94B3A202EB1C3F39}, // 1e-282
	{0xDAF5CCD93FB0CC53, 0xB9E08A83A5E34F07}, // 1e-281
	{0xD1B3400F8F9CFF68, 0xE858AD248F5C22C9}, // 1e-280
	{0x23100809B9C21FA1, 0x91376C36D99995BE}, // 1e-279
	{0xABD40A0C2832A78A, 0xB58547448FFFFB2D}, // 1e-278
	{0x16C90C8F323F516C, 0xE2E69915B3FFF9F9}, // 1e-277
	{0xAE3DA7D97F6792E3, 0x8DD01FAD907FFC3B}, // 1e-276
	{0x99CD11CFDF41779C, 0xB1442798F49FFB4A}, // 1e-275
	{0x40405643D711D583, 0xDD95317F31C7FA1D}, // 1e-274
	{0x482835EA666B2572, 0x8A7D3EEF7F1CFC52}, // 1e-273
	{0xDA3243650005EECF, 0xAD1C8EAB5EE43B66}, // 1e-272
	{0x90BED43E40076A82, 0xD863B256369D4A40}, // 1e-271
	{0x5A7744A6E804A291, 0x873E4F75E2224E68}, // 1e-270
	{0x711515D0A205CB36, 0xA90DE3535AAAE202}, // 1e-269
	{0x0D5A5B44CA873E03, 0xD3515C2831559A83}, // 1e-268
	{0xE858790AFE9486C2, 0x8412D9991ED58091}, // 1e-267
	{0x626E974DBE39A872, 0xA5178FFF668AE0B6}, // 1e-266
	{0xFB0A3D212DC8128F, 0xCE5D73FF402D98E3}, // 1e-265
	{0x7CE66634BC9D0B99, 0x80FA687F881C7F8E}, // 1e-264
	{0x1C1FFFC1EBC44E80, 0xA139029F6A239F72}, // 1e-263
	{0xA327FFB266B56220, 0xC987434744AC874E}, // 1e-262
	{0x4BF1FF9F0062BAA8, 0xFBE9141915D7A922}, // 1e-261
	{0x6F773FC3603DB4A9, 0x9D71AC8FADA6C9B5}, // 1e-260
	{0xCB550FB4384D21D3, 0xC4CE17B399107C22}, // 1e-259
	{0x7E2A53A146606A48, 0xF6019DA07F549B2B}, // 1e-258
	{0x2EDA7444CBFC426D, 0x99C102844F94E0FB}, // 1e-257
	{0xFA911155FEFB5308, 0xC0314325637A1939}, // 1e-256
	{0x793555AB7EBA27CA, 0xF03D93EEBC589F88}, // 1e-255
	{0x4BC1558B2F3458DE, 0x96267C7535B763B5}, // 1e-254
	{0x9EB1AAEDFB016F16, 0xBBB01B9283253CA2}, // 1e-253
	{0x465E15A979C1CADC, 0xEA9C227723EE8BCB}, // 1e-252
	{0x0BFACD89EC191EC9, 0x92A1958A7675175F}, // 1e-251
	{0xCEF980EC671F667B, 0xB749FAED14125D36}, // 1e-250
	{0x82B7E12780E7401A, 0xE51C79A85916F484}, // 1e-249
	{0xD1B2ECB8B0908810, 0x8F31CC0937AE58D2}, // 1e-248
	{0x861FA7E6DCB4AA15, 0xB2FE3F0B8599EF07}, // 1e-247
	{0x67A791E093E1D49A, 0xDFBDCECE67006AC9}, // 1e-246
	{0xE0C8BB2C5C6D24E0, 0x8BD6A141006042BD}, // 1e-245
	{0x58FAE9F773886E18, 0xAECC49914078536D}, // 1e-244
	{0xAF39A475506A899E, 0xDA7F5BF590966848}, // 1e-243
	{0x6D8406C952429603, 0x888F99797A5E012D}, // 1e-242
	{0xC8E5087BA6D33B83, 0xAAB37FD7D8F58178}, // 1e-241
	{0xFB1E4A9A90880A64, 0xD5605FCDCF32E1D6}, // 1e-240
	{0x5CF2EEA09A55067F, 0x855C3BE0A17FCD26}, // 1e-239
	{0xF42FAA48C0EA481E, 0xA6B34AD8C9DFC06F}, // 1e-238
	{0xF13B94DAF124DA26, 0xD0601D8EFC57B08B}, // 1e-237
	{0x76C53D08D6B70858, 0x823C12795DB6CE57}, // 1e-236
	{0x54768C4B0C64CA6E, 0xA2CB1717B52481ED}, // 1e-235
	{0xA9942F5DCF7DFD09, 0xCB7DDCDDA26DA268}, // 1e-234
	{0xD3F93B35435D7C4C, 0xFE5D54150B090B02}, // 1e-233
	{0xC47BC5014A1A6DAF, 0x9EFA548D26E5A6E1}, // 1e-232
	{0x359AB6419CA1091B, 0xC6B8E9B0709F109A}, // 1e-231
	{0xC30163D203C94B62, 0xF867241C8CC6D4C0}, // 1e-230
	{0x79E0DE63425DCF1D, 0x9B407691D7FC44F8}, // 1e-229
	{0x985915FC12F542E4, 0xC21094364DFB5636}, // 1e-228
	{0x3E6F5B7B17B2939D, 0xF294B943E17A2BC4}, // 1e-227
	{0xA705992CEECF9C42, 0x979CF3CA6CEC5B5A}, // 1e-226
	{0x50C6FF782A838353, 0xBD8430BD08277231}, // 1e-225
	{0xA4F8BF5635246428, 0xECE53CEC4A314EBD}, // 1e-224
	{0x871B7795E136BE99, 0x940F4613AE5ED136}, // 1e-223
	{0x28E2557B59846E3F, 0xB913179899F68584}, // 1e-222
	{0x331AEADA2FE589CF, 0xE757DD7EC07426E5}, // 1e-221
	{0x3FF0D2C85DEF7621, 0x9096EA6F3848984F}, // 1e-220
	{0x0FED077A756B53A9, 0xB4BCA50B065ABE63}, // 1e-219
	{0xD3E8495912C62894, 0xE1EBCE4DC7F16DFB}, // 1e-218
	{0x64712DD7ABBBD95C, 0x8D3360F09CF6E4BD}, // 1e-217
	{0xBD8D794D96AACFB3, 0xB080392CC4349DEC}, // 1e-216
	{0xECF0D7A0FC5583A0, 0xDCA04777F541C567}, // 1e-215
	{0xF41686C49DB57244, 0x89E42CAAF9491B60}, // 1e-214
	{0x311C2875C522CED5, 0xAC5D37D5B79B6239}, // 1e-213
	{0x7D633293366B828B, 0xD77485CB25823AC7}, // 1e-212
	{0xAE5DFF9C02033197, 0x86A8D39EF77164BC}, // 1e-211
	{0xD9F57F830283FDFC, 0xA8530886B54DBDEB}, // 1e-210
	{0xD072DF63C324FD7B, 0xD267CAA862A12D66}, // 1e-209
	{0x4247CB9E59F71E6D, 0x8380DEA93DA4BC60}, // 1e-208
	{0x52D9BE85F074E608, 0xA46116538D0DEB78}, // 1e-207
	{0x67902E276C921F8B, 0xCD795BE870516656}, // 1e-206
	{0x00BA1CD8A3DB53B6, 0x806BD9714632DFF6}, // 1e-205
	{0x80E8A40ECCD228A4, 0xA086CFCD97BF97F3}, // 1e-204
	{0x6122CD128006B2CD, 0xC8A883C0FDAF7DF0}, // 1e-203
	{0x796B805720085F81, 0xFAD2A4B13D1B5D6C}, // 1e-202
	{0xCBE3303674053BB0, 0x9CC3A6EEC6311A63}, // 1e-201
	{0xBEDBFC4411068A9C, 0xC3F490AA77BD60FC}, // 1e-200
	{0xEE92FB5515482D44, 0xF4F1B4D515ACB93B}, // 1e-199
	{0x751BDD152D4D1C4A, 0x991711052D8BF3C5}, // 1e-198
	{0xD262D45A78A0635D, 0xBF5CD54678EEF0B6}, // 1e-197
	{0x86FB897116C87C34, 0xEF340A98172AACE4}, // 1e-196
	{0xD45D35E6AE3D4DA0, 0x9580869F0E7AAC0E}, // 1e-195
	{0x8974836059CCA109, 0xBAE0A846D2195712}, // 1e-194
	{0x2BD1A438703FC94B, 0xE998D258869FACD7}, // 1e-193
	{0x7B6306A34627DDCF, 0x91FF83775423CC06}, // 1e-192
	{0x1A3BC84C17B1D542, 0xB67F6455292CBF08}, // 1e-191
	{0x20CABA5F1D9E4A93, 0xE41F3D6A7377EECA}, // 1e-190
	{0x547EB47B7282EE9C, 0x8E938662882AF53E}, // 1e-189
	{0xE99E619A4F23AA43, 0xB23867FB2A35B28D}, // 1e-188
	{0x6405FA00E2EC94D4, 0xDEC681F9F4C31F31}, // 1e-187
	{0xDE83BC408DD3DD04, 0x8B3C113C38F9F37E}, // 1e-186
	{0x9624AB50B148D445, 0xAE0B158B4738705E}, // 1e-185
	{0x3BADD624DD9B0957, 0xD98DDAEE19068C76}, // 1e-184
	{0xE54CA5D70A80E5D6, 0x87F8A8D4CFA417C9}, // 1e-183
	{0x5E9FCF4CCD211F4C, 0xA9F6D30A038D1DBC}, // 1e-182
	{0x7647C3200069671F, 0xD47487CC8470652B}, // 1e-181
	{0x29ECD9F40041E073, 0x84C8D4DFD2C63F3B}, // 1e-180
	{0xF468107100525890, 0xA5FB0A17C777CF09}, // 1e-179
	{0x7182148D4066EEB4, 0xCF79CC9DB955C2CC}, // 1e-178
	{0xC6F14CD848405530, 0x81AC1FE293D599BF}, // 1e-177
	{0xB8ADA00E5A506A7C, 0xA21727DB38CB002F}, // 1e-176
	{0xA6D90811F0E4851C, 0xCA9CF1D206FDC03B}, // 1e-175
	{0x908F4A166D1DA663, 0xFD442E4688BD304A}, // 1e-174
	{0x9A598E4E043287FE, 0x9E4A9CEC15763E2E}, // 1e-173
	{0x40EFF1E1853F29FD, 0xC5DD44271AD3CDBA}, // 1e-172
	{0xD12BEE59E68EF47C, 0xF7549530E188C128}, // 1e-171
	{0x82BB74F8301958CE, 0x9A94DD3E8CF578B9}, // 1e-170
	{0xE36A52363C1FAF01, 0xC13A148E3032D6E7}, // 1e-169
	{0xDC44E6C3CB279AC1, 0xF18899B1BC3F8CA1}, // 1e-168
	{0x29AB103A5EF8C0B9, 0x96F5600F15A7B7E5}, // 1e-167
	{0x7415D448F6B6F0E7, 0xBCB2B812DB11A5DE}, // 1e-166
	{0x111B495B3464AD21, 0xEBDF661791D60F56}, // 1e-165
	{0xCAB10DD900BEEC34, 0x936B9FCEBB25C995}, // 1e-164
	{0x3D5D514F40EEA742, 0xB84687C269EF3BFB}, // 1e-163
	{0x0CB4A5A3112A5112, 0xE65829B3046B0AFA}, // 1e-162
	{0x47F0E785EABA72AB, 0x8FF71A0FE2C2E6DC}, // 1e-161
	{0x59ED216765690F56, 0xB3F4E093DB73A093}, // 1e-160
	{0x306869C13EC3532C, 0xE0F218B8D25088B8}, // 1e-159
	{0x1E414218C73A13FB, 0x8C974F7383725573}, // 1e-158
	{0xE5D1929EF90898FA, 0xAFBD2350644EEACF}, // 1e-157
	{0xDF45F746B74ABF39, 0xDBAC6C247D62A583}, // 1e-156
	{0x6B8BBA8C328EB783, 0x894BC396CE5DA772}, // 1e-155
	{0x066EA92F3F326564, 0xAB9EB47C81F5114F}, // 1e-154
	{0xC80A537B0EFEFEBD, 0xD686619BA27255A2}, // 1e-153
	{0xBD06742CE95F5F36, 0x8613FD0145877585}, // 1e-152
	{0x2C48113823B73704, 0xA798FC4196E952E7}, // 1e-151
	{0xF75A15862CA504C5, 0xD17F3B51FCA3A7A0}, // 1e-150
	{0x9A984D73DBE722FB, 0x82EF85133DE648C4}, // 1e-149
	{0xC13E60D0D2E0EBBA, 0xA3AB66580D5FDAF5}, // 1e-148
	{0x318DF905079926A8, 0xCC963FEE10B7D1B3}, // 1e-147
	{0xFDF17746497F7052, 0xFFBBCFE994E5C61F}, // 1e-146
	{0xFEB6EA8BEDEFA633, 0x9FD561F1FD0F9BD3}, // 1e-145
	{0xFE64A52EE96B8FC0, 0xC7CABA6E7C5382C8}, // 1e-144
	{0x3DFDCE7AA3C673B0, 0xF9BD690A1B68637B}, // 1e-143
	{0x06BEA10CA65C084E, 0x9C1661A651213E2D}, // 1e-142
	{0x486E494FCFF30A62, 0xC31BFA0FE5698DB8}, // 1e-141
	{0x5A89DBA3C3EFCCFA, 0xF3E2F893DEC3F126}, // 1e-140
	{0xF89629465A75E01C, 0x986DDB5C6B3A76B7}, // 1e-139
	{0xF6BBB397F1135823, 0xBE89523386091465}, // 1e-138
	{0x746AA07DED582E2C, 0xEE2BA6C0678B597F}, // 1e-137
	{0xA8C2A44EB4571CDC, 0x94DB483840B717EF}, // 1e-136
	{0x92F34D62616CE413, 0xBA121A4650E4DDEB}, // 1e-135
	{0x77B020BAF9C81D17, 0xE896A0D7E51E1566}, // 1e-134
	{0x0ACE1474DC1D122E, 0x915E2486EF32CD60}, // 1e-133
	{0x0D819992132456BA, 0xB5B5ADA8AAFF80B8}, // 1e-132
	{0x10E1FFF697ED6C69, 0xE3231912D5BF60E6}, // 1e-131
	{0xCA8D3FFA1EF463C1, 0x8DF5EFABC5979C8F}, // 1e-130
	{0xBD308FF8A6B17CB2, 0xB1736B96B6FD83B3}, // 1e-129
	{0xAC7CB3F6D05DDBDE, 0xDDD0467C64BCE4A0}, // 1e-128
	{0x6BCDF07A423AA96B, 0x8AA22C0DBEF60EE4}, // 1e-127
	{0x86C16C98D2C953C6, 0xAD4AB7112EB3929D}, // 1e-126
	{0xE871C7BF077BA8B7, 0xD89D64D57A607744}, // 1e-125
	{0x11471CD764AD4972, 0x87625F056C7C4A8B}, // 1e-124
	{0xD598E40D3DD89BCF, 0xA93AF6C6C79B5D2D}, // 1e-123
	{0x4AFF1D108D4EC2C3, 0xD389B47879823479}, // 1e-122
	{0xCEDF722A585139BA, 0x843610CB4BF160CB}, // 1e-121
	{0xC2974EB4EE658828, 0xA54394FE1EEDB8FE}, // 1e-120
	{0x733D226229FEEA32, 0xCE947A3DA6A9273E}, // 1e-119
	{0x0806357D5A3F525F, 0x811CCC668829B887}, // 1e-118
	{0xCA07C2DCB0CF26F7, 0xA163FF802A3426A8}, // 1e-117
	{0xFC89B393DD02F0B5, 0xC9BCFF6034C13052}, // 1e-116
	{0xBBAC2078D443ACE2, 0xFC2C3F3841F17C67}, // 1e-115
	{0xD54B944B84AA4C0D, 0x9D9BA7832936EDC0}, // 1e-114
	{0x0A9E795E65D4DF11, 0xC5029163F384A931}, // 1e-113
	{0x4D4617B5FF4A16D5, 0xF64335BCF065D37D}, // 1e-112
	{0x504BCED1BF8E4E45, 0x99EA0196163FA42E}, // 1e-111
	{0xE45EC2862F71E1D6, 0xC06481FB9BCF8D39}, // 1e-110
	{0x5D767327BB4E5A4C, 0xF07DA27A82C37088}, // 1e-109
	{0x3A6A07F8D510F86F, 0x964E858C91BA2655}, // 1e-108
	{0x890489F70A55368B, 0xBBE226EFB628AFEA}, // 1e-107
	{0x2B45AC74CCEA842E, 0xEADAB0ABA3B2DBE5}, // 1e-106
	{0x3B0B8BC90012929D, 0x92C8AE6B464FC96F}, // 1e-105
	{0x09CE6EBB40173744, 0xB77ADA0617E3BBCB}, // 1e-104
	{0xCC420A6A101D0515, 0xE55990879DDCAABD}, // 1e-103
	{0x9FA946824A12232D, 0x8F57FA54C2A9EAB6}, // 1e-102
	{0x47939822DC96ABF9, 0xB32DF8E9F3546564}, // 1e-101
	{0x59787E2B93BC56F7, 0xDFF9772470297EBD}, // 1e-100
	{0x57EB4EDB3C55B65A, 0x8BFBEA76C619EF36}, // 1e-99
	{0xEDE622920B6B23F1, 0xAEFAE51477A06B03}, // 1e-98
	{0xE95FAB368E45ECED, 0xDAB99E59958885C4}, // 1e-97
	{0x11DBCB0218EBB414, 0x88B402F7FD75539B}, // 1e-96
	{0xD652BDC29F26A119, 0xAAE103B5FCD2A881}, // 1e-95
	{0x4BE76D3346F0495F, 0xD59944A37C0752A2}, // 1e-94
	{0x6F70A4400C562DDB, 0x857FCAE62D8493A5}, // 1e-93
	{0xCB4CCD500F6BB952, 0xA6DFBD9FB8E5B88E}, // 1e-92
	{0x7E2000A41346A7A7, 0xD097AD07A71F26B2}, // 1e-91
	{0x8ED400668C0C28C8, 0x825ECC24C873782F}, // 1e-90
	{0x728900802F0F32FA, 0xA2F67F2DFA90563B}, // 1e-89
	{0x4F2B40A03AD2FFB9, 0xCBB41EF979346BCA}, // 1e-88
	{0xE2F610C84987BFA8, 0xFEA126B7D78186BC}, // 1e-87
	{0x0DD9CA7D2DF4D7C9, 0x9F24B832E6B0F436}, // 1e-86
	{0x91503D1C79720DBB, 0xC6EDE63FA05D3143}, // 1e-85
	{0x75A44C6397CE912A, 0xF8A95FCF88747D94}, // 1e-84
	{0xC986AFBE3EE11ABA, 0x9B69DBE1B548CE7C}, // 1e-83
	{0xFBE85BADCE996168, 0xC24452DA229B021B}, // 1e-82
	{0xFAE27299423FB9C3, 0xF2D56790AB41C2A2}, // 1e-81
	{0xDCCD879FC967D41A, 0x97C560BA6B0919A5}, // 1e-80
	{0x5400E987BBC1C920, 0xBDB6B8E905CB600F}, // 1e-79
	{0x290123E9AAB23B68, 0xED246723473E3813}, // 1e-78
	{0xF9A0B6720AAF6521, 0x9436C0760C86E30B}, // 1e-77
	{0xF808E40E8D5B3E69, 0xB94470938FA89BCE}, // 1e-76
	{0xB60B1D1230B20E04, 0xE7958CB87392C2C2}, // 1e-75
	{0xB1C6F22B5E6F48C2, 0x90BD77F3483BB9B9}, // 1e-74
	{0x1E38AEB6360B1AF3, 0xB4ECD5F01A4AA828}, // 1e-73
	{0x25C6DA63C38DE1B0, 0xE2280B6C20DD5232}, // 1e-72
	{0x579C487E5A38AD0E, 0x8D590723948A535F}, // 1e-71
	{0x2D835A9DF0C6D851, 0xB0AF48EC79ACE837}, // 1e-70
	{0xF8E431456CF88E65, 0xDCDB1B2798182244}, // 1e-69
	{0x1B8E9ECB641B58FF, 0x8A08F0F8BF0F156B}, // 1e-68
	{0xE272467E3D222F3F, 0xAC8B2D36EED2DAC5}, // 1e-67
	{0x5B0ED81DCC6ABB0F, 0xD7ADF884AA879177}, // 1e-66
	{0x98E947129FC2B4E9, 0x86CCBB52EA94BAEA}, // 1e-65
	{0x3F2398D747B36224, 0xA87FEA27A539E9A5}, // 1e-64
	{0x8EEC7F0D19A03AAD, 0xD29FE4B18E88640E}, // 1e-63
	{0x1953CF68300424AC, 0x83A3EEEEF9153E89}, // 1e-62
	{0x5FA8C3423C052DD7, 0xA48CEAAAB75A8E2B}, // 1e-61
	{0x3792F412CB06794D, 0xCDB02555653131B6}, // 1e-60
	{0xE2BBD88BBEE40BD0, 0x808E17555F3EBF11}, // 1e-59
	{0x5B6ACEAEAE9D0EC4, 0xA0B19D2AB70E6ED6}, // 1e-58
	{0xF245825A5A445275, 0xC8DE047564D20A8B}, // 1e-57
	{0xEED6E2F0F0D56712, 0xFB158592BE068D2E}, // 1e-56
	{0x55464DD69685606B, 0x9CED737BB6C4183D}, // 1e-55
	{0xAA97E14C3C26B886, 0xC428D05AA4751E4C}, // 1e-54
	{0xD53DD99F4B3066A8, 0xF53304714D9265DF}, // 1e-53
	{0xE546A8038EFE4029, 0x993FE2C6D07B7FAB}, // 1e-52
	{0xDE98520472BDD033, 0xBF8FDB78849A5F96}, // 1e-51
	{0x963E66858F6D4440, 0xEF73D256A5C0F77C}, // 1e-50
	{0xDDE7001379A44AA8, 0x95A8637627989AAD}, // 1e-49
	{0x5560C018580D5D52, 0xBB127C53B17EC159}, // 1e-48
	{0xAAB8F01E6E10B4A6, 0xE9D71B689DDE71AF}, // 1e-47
	{0xCAB3961304CA70E8, 0x9226712162AB070D}, // 1e-46
	{0x3D607B97C5FD0D22, 0xB6B00D69BB55C8D1}, // 1e-45
	{0x8CB89A7DB77C506A, 0xE45C10C42A2B3B05}, // 1e-44
	{0x77F3608E92ADB242, 0x8EB98A7A9A5B04E3}, // 1e-43
	{0x55F038B237591ED3, 0xB267ED1940F1C61C}, // 1e-42
	{0x6B6C46DEC52F6688, 0xDF01E85F912E37A3}, // 1e-41
	{0x2323AC4B3B3DA015, 0x8B61313BBABCE2C6}, // 1e-40
	{0xABEC975E0A0D081A, 0xAE397D8AA96C1B77}, // 1e-39
	{0x96E7BD358C904A21, 0xD9C7DCED53C72255}, // 1e-38
	{0x7E50D64177DA2E54, 0x881CEA14545C7575}, // 1e-37
	{0xDDE50BD1D5D0B9E9, 0xAA242499697392D2}, // 1e-36
	{0x955E4EC64B44E864, 0xD4AD2DBFC3D07787}, // 1e-35
	{0xBD5AF13BEF0B113E, 0x84EC3C97DA624AB4}, // 1e-34
	{0xECB1AD8AEACDD58E, 0xA6274BBDD0FADD61}, // 1e-33
	{0x67DE18EDA5814AF2, 0xCFB11EAD453994BA}, // 1e-32
	{0x80EACF948770CED7, 0x81CEB32C4B43FCF4}, // 1e-31
	{0xA1258379A94D028D, 0xA2425FF75E14FC31}, // 1e-30
	{0x096EE45813A04330, 0xCAD2F7F5359A3B3E}, // 1e-29
	{0x8BCA9D6E188853FC, 0xFD87B5F28300CA0D}, // 1e-28
	{0x775EA264CF55347D, 0x9E74D1B791E07E48}, // 1e-27
	{0x95364AFE032A819D, 0xC612062576589DDA}, // 1e-26
	{0x3A83DDBD83F52204, 0xF79687AED3EEC551}, // 1e-25
	{0xC4926A9672793542, 0x9ABE14CD44753B52}, // 1e-24
	{0x75B7053C0F178293, 0xC16D9A0095928A27}, // 1e-23
	{0x5324C68B12DD6338, 0xF1C90080BAF72CB1}, // 1e-22
	{0xD3F6FC16EBCA5E03, 0x971DA05074DA7BEE}, // 1e-21
	{0x88F4BB1CA6BCF584, 0xBCE5086492111AEA}, // 1e-20
	{0x2B31E9E3D06C32E5, 0xEC1E4A7DB69561A5}, // 1e-19
	{0x3AFF322E62439FCF, 0x9392EE8E921D5D07}, // 1e-18
	{0x09BEFEB9FAD487C2, 0xB877AA3236A4B449}, // 1e-17
	{0x4C2EBE687989A9B3, 0xE69594BEC44DE15B}, // 1e-16
	{0x0F9D37014BF60A10, 0x901D7CF73AB0ACD9}, // 1e-15
	{0x538484C19EF38C94, 0xB424DC35095CD80F}, // 1e-14
	{0x2865A5F206B06FB9, 0xE12E13424BB40E13}, // 1e-13
	{0xF93F87B7442E45D3, 0x8CBCCC096F5088CB}, // 1e-12
	{0xF78F69A51539D748, 0xAFEBFF0BCB24AAFE}, // 1e-11
	{0xB573440E5A884D1B, 0xDBE6FECEBDEDD5BE}, // 1e-10
	{0x31680A88F8953030, 0x89705F4136B4A597}, // 1e-9
	{0xFDC20D2B36BA7C3D, 0xABCC77118461CEFC}, // 1e-8
	{0x3D32907604691B4C, 0xD6BF94D5E57A42BC}, // 1e-7
	{0xA63F9A49C2C1B10F, 0x8637BD05AF6C69B5}, // 1e-6
	{0x0FCF80DC33721D53, 0xA7C5AC471B478423}, // 1e-5
	{0xD3C36113404EA4A8, 0xD1B71758E219652B}, // 1e-4
	{0x645A1CAC083126E9, 0x83126E978D4FDF3B}, // 1e-3
	{0x3D70A3D70A3D70A3, 0xA3D70A3D70A3D70A}, // 1e-2
	{0xCCCCCCCCCCCCCCCC, 0xCCCCCCCCCCCCCCCC}, // 1e-1
	{0x0000000000000000, 0x8000000000000000}, // 1e0
	{0x0000000000000000, 0xA000000000000000}, // 1e1
	{0x0000000000000000, 0xC800000000000000}, // 1e2
	{0x0000000000000000, 0xFA00000000000000}, // 1e3
	{0x0000000000000000, 0x9C40000000000000}, // 1e4
	{0x0000000000000000, 0xC350000000000000}, // 1e5
	{0x0000000000000000, 0xF424000000000000}, // 1e6
	{0x0000000000000000, 0x9896800000000000}, // 1e7
	{0x0000000000000000, 0xBEBC200000000000}, // 1e8
	{0x0000000000000000, 0xEE6B280000000000}, // 1e9
	{0x0000000000000000, 0x9502F90000000000}, // 1e10
	{0x0000000000000000, 0xBA43B74000000000}, // 1e11
	{0x0000000000000000, 0xE8D4A51000000000}, // 1e12
	{0x0000000000000000, 0x9184E72A00000000}, // 1e13
	{0x0000000000000000, 0xB5E620F480000000}, // 1e14
	{0x0000000000000000, 0xE35FA931A0000000}, // 1e15
	{0x0000000000000000, 0x8E1BC9BF04000000}, // 1e16
	{0x0000000000000000, 0xB1A2BC2EC5000000}, // 1e17
	{0x0000000000000000, 0xDE0B6B3A76400000}, // 1e18
	{0x0000000000000000, 0x8AC7230489E80000}, // 1e19
	{0x0000000000000000, 0xAD78EBC5AC620000}, // 1e20
	{0x0000000000000000, 0xD8D726B7177A8000}, // 1e21
	{0x0000000000000000, 0x878678326EAC9000}, // 1e22
	{0x0000000000000000, 0xA968163F0A57B400}, // 1e23
	{0x0000000000000000, 0xD3C21BCECCEDA100}, // 1e24
	{0x0000000000000000, 0x84595161401484A0}, // 1e25
	{0x0000000000000000, 0xA56FA5B99019A5C8}, // 1e26
	{0x0000000000000000, 0xCECB8F27F4200F3A}, // 1e27
	{0x4000000000000000, 0x813F3978F8940984}, // 1e28
	{0x5000000000000000, 0xA18F07D736B90BE5}, // 1e29
	{0xA400000000000000, 0xC9F2C9CD04674EDE}, // 1e30
	{0x4D00000000000000, 0xFC6F7C4045812296}, // 1e31
	{0xF020000000000000, 0x9DC5ADA82B70B59D}, // 1e32
	{0x6C28000000000000, 0xC5371912364CE305}, // 1e33
	{0xC732000000000000, 0xF684DF56C3E01BC6}, // 1e34
	{0x3C7F400000000000, 0x9A130B963A6C115C}, // 1e35
	{0x4B9F100000000000, 0xC097CE7BC90715B3}, // 1e36
	{0x1E86D40000000000, 0xF0BDC21ABB48DB20}, // 1e37
	{0x1314448000000000, 0x96769950B50D88F4}, // 1e38
	{0x17D955A000000000, 0xBC143FA4E250EB31}, // 1e39
	{0x5DCFAB0800000000, 0xEB194F8E1AE525FD}, // 1e40
	{0x5AA1CAE500000000, 0x92EFD1B8D0CF37BE}, // 1e41
	{0xF14A3D9E40000000, 0xB7ABC627050305AD}, // 1e42
	{0x6D9CCD05D0000000, 0xE596B7B0C643C719}, // 1e43
	{0xE4820023A2000000, 0x8F7E32CE7BEA5C6F}, // 1e44
	{0xDDA2802C8A800000, 0xB35DBF821AE4F38B}, // 1e45
	{0xD50B2037AD200000, 0xE0352F62A19E306E}, // 1e46
	{0x4526F422CC340000, 0x8C213D9DA502DE45}, // 1e47
	{0x9670B12B7F410000, 0xAF298D050E4395D6}, // 1e48
	{0x3C0CDD765F114000, 0xDAF3F04651D47B4C}, // 1e49
	{0xA5880A69FB6AC800, 0x88D8762BF324CD0F}, // 1e50
	{0x8EEA0D047A457A00, 0xAB0E93B6EFEE0053}, // 1e51
	{0x72A4904598D6D880, 0xD5D238A4ABE98068}, // 1e52
	{0x47A6DA2B7F864750, 0x85A36366EB71F041}, // 1e53
	{0x999090B65F67D924, 0xA70C3C40A64E6C51}, // 1e54
	{0xFFF4B4E3F741CF6D, 0xD0CF4B50CFE20765}, // 1e55
	{0xBFF8F10E7A8921A4, 0x82818F1281ED449F}, // 1e56
	{0xAFF72D52192B6A0D, 0xA321F2D7226895C7}, // 1e57
	{0x9BF4F8A69F764490, 0xCBEA6F8CEB02BB39}, // 1e58
	{0x02F236D04753D5B4, 0xFEE50B7025C36A08}, // 1e59
	{0x01D762422C946590, 0x9F4F2726179A2245}, // 1e60
	{0x424D3AD2B7B97EF5, 0xC722F0EF9D80AAD6}, // 1e61
	{0xD2E0898765A7DEB2, 0xF8EBAD2B84E0D58B}, // 1e62
	{0x63CC55F49F88EB2F, 0x9B934C3B330C8577}, // 1e63
	{0x3CBF6B71C76B25FB, 0xC2781F49FFCFA6D5}, // 1e64
	{0x8BEF464E3945EF7A, 0xF316271C7FC3908A}, // 1e65
	{0x97758BF0E3CBB5AC, 0x97EDD871CFDA3A56}, // 1e66
	{0x3D52EEED1CBEA317, 0xBDE94E8E43D0C8EC}, // 1e67
	{0x4CA7AAA863EE4BDD, 0xED63A231D4C4FB27}, // 1e68
	{0x8FE8CAA93E74EF6A, 0x945E455F24FB1CF8}, // 1e69
	{0xB3E2FD538E122B44, 0xB975D6B6EE39E436}, // 1e70
	{0x60DBBCA87196B616, 0xE7D34C64A9C85D44}, // 1e71
	{0xBC8955E946FE31CD, 0x90E40FBEEA1D3A4A}, // 1e72
	{0x6BABAB6398BDBE41, 0xB51D13AEA4A488DD}, // 1e73
	{0xC696963C7EED2DD1, 0xE264589A4DCDAB14}, // 1e74
	{0xFC1E1DE5CF543CA2, 0x8D7EB76070A08AEC}, // 1e75
	{0x3B25A55F43294BCB, 0xB0DE65388CC8ADA8}, // 1e76
	{0x49EF0EB713F39EBE, 0xDD15FE86AFFAD912}, // 1e77
	{0x6E3569326C784337, 0x8A2DBF142DFCC7AB}, // 1e78
	{0x49C2C37F07965404, 0xACB92ED9397BF996}, // 1e79
	{0xDC33745EC97BE906, 0xD7E77A8F87DAF7FB}, // 1e80
	{0x69A028BB3DED71A3, 0x86F0AC99B4E8DAFD}, // 1e81
	{0xC40832EA0D68CE0C, 0xA8ACD7C0222311BC}, // 1e82
	{0xF50A3FA490C30190, 0xD2D80DB02AABD62B}, // 1e83
	{0x792667C6DA79E0FA, 0x83C7088E1AAB65DB}, // 1e84
	{0x577001B891185938, 0xA4B8CAB1A1563F52}, // 1e85
	{0xED4C0226B55E6F86, 0xCDE6FD5E09ABCF26}, // 1e86
	{0x544F8158315B05B4, 0x80B05E5AC60B6178}, // 1e87
	{0x696361AE3DB1C721, 0xA0DC75F1778E39D6}, // 1e88
	{0x03BC3A19CD1E38E9, 0xC913936DD571C84C}, // 1e89
	{0x04AB48A04065C723, 0xFB5878494ACE3A5F}, // 1e90
	{0x62EB0D64283F9C76, 0x9D174B2DCEC0E47B}, // 1e91
	{0x3BA5D0BD324F8394, 0xC45D1DF942711D9A}, // 1e92
	{0xCA8F44EC7EE36479, 0xF5746577930D6500}, // 1e93
	{0x7E998B13CF4E1ECB, 0x9968BF6ABBE85F20}, // 1e94
	{0x9E3FEDD8C321A67E, 0xBFC2EF456AE276E8}, // 1e95
	{0xC5CFE94EF3EA101E, 0xEFB3AB16C59B14A2}, // 1e96
	{0xBBA1F1D158724A12, 0x95D04AEE3B80ECE5}, // 1e97
	{0x2A8A6E45AE8EDC97, 0xBB445DA9CA61281F}, // 1e98
	{0xF52D09D71A3293BD, 0xEA1575143CF97226}, // 1e99
	{0x593C2626705F9C56, 0x924D692CA61BE758}, // 1e100
	{0x6F8B2FB00C77836C, 0xB6E0C377CFA2E12E}, // 1e101
	{0x0B6DFB9C0F956447, 0xE498F455C38B997A}, // 1e102
	{0x4724BD4189BD5EAC, 0x8EDF98B59A373FEC}, // 1e103
	{0x58EDEC91EC2CB657, 0xB2977EE300C50FE7}, // 1e104
	{0x2F2967B66737E3ED, 0xDF3D5E9BC0F653E1}, // 1e105
	{0xBD79E0D20082EE74, 0x8B865B215899F46C}, // 1e106
	{0xECD8590680A3AA11, 0xAE67F1E9AEC07187}, // 1e107
	{0xE80E6F4820CC9495, 0xDA01EE641A708DE9}, // 1e108
	{0x3109058D147FDCDD, 0x884134FE908658B2}, // 1e109
	{0xBD4B46F0599FD415, 0xAA51823E34A7EEDE}, // 1e110
	{0x6C9E18AC7007C91A, 0xD4E5E2CDC1D1EA96}, // 1e111
	{0x03E2CF6BC604DDB0, 0x850FADC09923329E}, // 1e112
	{0x84DB8346B786151C, 0xA6539930BF6BFF45}, // 1e113
	{0xE612641865679A63, 0xCFE87F7CEF46FF16}, // 1e114
	{0x4FCB7E8F3F60C07E, 0x81F14FAE158C5F6E}, // 1e115
	{0xE3BE5E330F38F09D, 0xA26DA3999AEF7749}, // 1e116
	{0x5CADF5BFD3072CC5, 0xCB090C8001AB551C}, // 1e117
	{0x73D9732FC7C8F7F6, 0xFDCB4FA002162A63}, // 1e118
	{0x2867E7FDDCDD9AFA, 0x9E9F11C4014DDA7E}, // 1e119
	{0xB281E1FD541501B8, 0xC646D63501A1511D}, // 1e120
	{0x1F225A7CA91A4226, 0xF7D88BC24209A565}, // 1e121
	{0x3375788DE9B06958, 0x9AE757596946075F}, // 1e122
	{0x0052D6B1641C83AE, 0xC1A12D2FC3978937}, // 1e123
	{0xC0678C5DBD23A49A, 0xF209787BB47D6B84}, // 1e124
	{0xF840B7BA963646E0, 0x9745EB4D50CE6332}, // 1e125
	{0xB650E5A93BC3D898, 0xBD176620A501FBFF}, // 1e126
	{0xA3E51F138AB4CEBE, 0xEC5D3FA8CE427AFF}, // 1e127
	{0xC66F336C36B10137, 0x93BA47C980E98CDF}, // 1e128
	{0xB80B0047445D4184, 0xB8A8D9BBE123F017}, // 1e129
	{0xA60DC059157491E5, 0xE6D3102AD96CEC1D}, // 1e130
	{0x87C89837AD68DB2F, 0x9043EA1AC7E41392}, // 1e131
	{0x29BABE4598C311FB, 0xB454E4A179DD1877}, // 1e132
	{0xF4296DD6FEF3D67A, 0xE16A1DC9D8545E94}, // 1e133
	{0x1899E4A65F58660C, 0x8CE2529E2734BB1D}, // 1e134
	{0x5EC05DCFF72E7F8F, 0xB01AE745B101E9E4}, // 1e135
	{0x76707543F4FA1F73, 0xDC21A1171D42645D}, // 1e136
	{0x6A06494A791C53A8, 0x899504AE72497EBA}, // 1e137
	{0x0487DB9D17636892, 0xABFA45DA0EDBDE69}, // 1e138
	{0x45A9D2845D3C42B6, 0xD6F8D7509292D603}, // 1e139
	{0x0B8A2392BA45A9B2, 0x865B86925B9BC5C2}, // 1e140
	{0x8E6CAC7768D7141E, 0xA7F26836F282B732}, // 1e141
	{0x3207D795430CD926, 0xD1EF0244AF2364FF}, // 1e142
	{0x7F44E6BD49E807B8, 0x8335616AED761F1F}, // 1e143
	{0x5F16206C9C6209A6, 0xA402B9C5A8D3A6E7}, // 1e144
	{0x36DBA887C37A8C0F, 0xCD036837130890A1}, // 1e145
	{0xC2494954DA2C9789, 0x802221226BE55A64}, // 1e146
	{0xF2DB9BAA10B7BD6C, 0xA02AA96B06DEB0FD}, // 1e147
	{0x6F92829494E5ACC7, 0xC83553C5C8965D3D}, // 1e148
	{0xCB772339BA1F17F9, 0xFA42A8B73ABBF48C}, // 1e149
	{0xFF2A760414536EFB, 0x9C69A97284B578D7}, // 1e150
	{0xFEF5138519684ABA, 0xC38413CF25E2D70D}, // 1e151
	{0x7EB258665FC25D69, 0xF46518C2EF5B8CD1}, // 1e152
	{0xEF2F773FFBD97A61, 0x98BF2F79D5993802}, // 1e153
	{0xAAFB550FFACFD8FA, 0xBEEEFB584AFF8603}, // 1e154
	{0x95BA2A53F983CF38, 0xEEAABA2E5DBF6784}, // 1e155
	{0xDD945A747BF26183, 0x952AB45CFA97A0B2}, // 1e156
	{0x94F971119AEEF9E4, 0xBA756174393D88DF}, // 1e157
	{0x7A37CD5601AAB85D, 0xE912B9D1478CEB17}, // 1e158
	{0xAC62E055C10AB33A, 0x91ABB422CCB812EE}, // 1e159
	{0x577B986B314D6009, 0xB616A12B7FE617AA}, // 1e160
	{0xED5A7E85FDA0B80B, 0xE39C49765FDF9D94}, // 1e161
	{0x14588F13BE847307, 0x8E41ADE9FBEBC27D}, // 1e162
	{0x596EB2D8AE258FC8, 0xB1D219647AE6B31C}, // 1e163
	{0x6FCA5F8ED9AEF3BB, 0xDE469FBD99A05FE3}, // 1e164
	{0x25DE7BB9480D5854, 0x8AEC23D680043BEE}, // 1e165
	{0xAF561AA79A10AE6A, 0xADA72CCC20054AE9}, // 1e166
	{0x1B2BA1518094DA04, 0xD910F7FF28069DA4}, // 1e167
	{0x90FB44D2F05D0842, 0x87AA9AFF79042286}, // 1e168
	{0x353A1607AC744A53, 0xA99541BF57452B28}, // 1e169
	{0x42889B8997915CE8, 0xD3FA922F2D1675F2}, // 1e170
	{0x69956135FEBADA11, 0x847C9B5D7C2E09B7}, // 1e171
	{0x43FAB9837E699095, 0xA59BC234DB398C25}, // 1e172
	{0x94F967E45E03F4BB, 0xCF02B2C21207EF2E}, // 1e173
	{0x1D1BE0EEBAC278F5, 0x8161AFB94B44F57D}, // 1e174
	{0x6462D92A69731732, 0xA1BA1BA79E1632DC}, // 1e175
	{0x7D7B8F7503CFDCFE, 0xCA28A291859BBF93}, // 1e176
	{0x5CDA735244C3D43E, 0xFCB2CB35E702AF78}, // 1e177
	{0x3A0888136AFA64A7, 0x9DEFBF01B061ADAB}, // 1e178
	{0x088AAA1845B8FDD0, 0xC56BAEC21C7A1916}, // 1e179
	{0x8AAD549E57273D45, 0xF6C69A72A3989F5B}, // 1e180
	{0x36AC54E2F678864B, 0x9A3C2087A63F6399}, // 1e181
	{0x84576A1BB416A7DD, 0xC0CB28A98FCF3C7F}, // 1e182
	{0x656D44A2A11C51D5, 0xF0FDF2D3F3C30B9F}, // 1e183
	{0x9F644AE5A4B1B325, 0x969EB7C47859E743}, // 1e184
	{0x873D5D9F0DDE1FEE, 0xBC4665B596706114}, // 1e185
	{0xA90CB506D155A7EA, 0xEB57FF22FC0C7959}, // 1e186
	{0x09A7F12442D588F2, 0x9316FF75DD87CBD8}, // 1e187
	{0x0C11ED6D538AEB2F, 0xB7DCBF5354E9BECE}, // 1e188
	{0x8F1668C8A86DA5FA, 0xE5D3EF282A242E81}, // 1e189
	{0xF96E017D694487BC, 0x8FA475791A569D10}, // 1e190
	{0x37C981DCC395A9AC, 0xB38D92D760EC4455}, // 1e191
	{0x85BBE253F47B1417, 0xE070F78D3927556A}, // 1e192
	{0x93956D7478CCEC8E, 0x8C469AB843B89562}, // 1e193
	{0x387AC8D1970027B2, 0xAF58416654A6BABB}, // 1e194
	{0x06997B05FCC0319E, 0xDB2E51BFE9D0696A}, // 1e195
	{0x441FECE3BDF81F03, 0x88FCF317F22241E2}, // 1e196
	{0xD527E81CAD7626C3, 0xAB3C2FDDEEAAD25A}, // 1e197
	{0x8A71E223D8D3B074, 0xD60B3BD56A5586F1}, // 1e198
	{0xF6872D5667844E49, 0x85C7056562757456}, // 1e199
	{0xB428F8AC016561DB, 0xA738C6BEBB12D16C}, // 1e200
	{0xE13336D701BEBA52, 0xD106F86E69D785C7}, // 1e201
	{0xECC0024661173473, 0x82A45B450226B39C}, // 1e202
	{0x27F002D7F95D0190, 0xA34D721642B06084}, // 1e203
	{0x31EC038DF7B441F4, 0xCC20CE9BD35C78A5}, // 1e204
	{0x7E67047175A15271, 0xFF290242C83396CE}, // 1e205
	{0x0F0062C6E984D386, 0x9F79A169BD203E41}, // 1e206
	{0x52C07B78A3E60868, 0xC75809C42C684DD1}, // 1e207
	{0xA7709A56CCDF8A82, 0xF92E0C3537826145}, // 1e208
	{0x88A66076400BB691, 0x9BBCC7A142B17CCB}, // 1e209
	{0x6ACFF893D00EA435, 0xC2ABF989935DDBFE}, // 1e210
	{0x0583F6B8C4124D43, 0xF356F7EBF83552FE}, // 1e211
	{0xC3727A337A8B704A, 0x98165AF37B2153DE}, // 1e212
	{0x744F18C0592E4C5C, 0xBE1BF1B059E9A8D6}, // 1e213
	{0x1162DEF06F79DF73, 0xEDA2EE1C7064130C}, // 1e214
	{0x8ADDCB5645AC2BA8, 0x9485D4D1C63E8BE7}, // 1e215
	{0x6D953E2BD7173692, 0xB9A74A0637CE2EE1}, // 1e216
	{0xC8FA8DB6CCDD0437, 0xE8111C87C5C1BA99}, // 1e217
	{0x1D9C9892400A22A2, 0x910AB1D4DB9914A0}, // 1e218
	{0x2503BEB6D00CAB4B, 0xB54D5E4A127F59C8}, // 1e219
	{0x2E44AE64840FD61D, 0xE2A0B5DC971F303A}, // 1e220
	{0x5CEAECFED289E5D2, 0x8DA471A9DE737E24}, // 1e221
	{0x7425A83E872C5F47, 0xB10D8E1456105DAD}, // 1e222
	{0xD12F124E28F77719, 0xDD50F1996B947518}, // 1e223
	{0x82BD6B70D99AAA6F, 0x8A5296FFE33CC92F}, // 1e224
	{0x636CC64D1001550B, 0xACE73CBFDC0BFB7B}, // 1e225
	{0x3C47F7E05401AA4E, 0xD8210BEFD30EFA5A}, // 1e226
	{0x65ACFAEC34810A71, 0x8714A775E3E95C78}, // 1e227
	{0x7F1839A741A14D0D, 0xA8D9D1535CE3B396}, // 1e228
	{0x1EDE48111209A050, 0xD31045A8341CA07C}, // 1e229
	{0x934AED0AAB460432, 0x83EA2B892091E44D}, // 1e230
	{0xF81DA84D5617853F, 0xA4E4B66B68B65D60}, // 1e231
	{0x36251260AB9D668E, 0xCE1DE40642E3F4B9}, // 1e232
	{0xC1D72B7C6B426019, 0x80D2AE83E9CE78F3}, // 1e233
	{0xB24CF65B8612F81F, 0xA1075A24E4421730}, // 1e234
	{0xDEE033F26797B627, 0xC94930AE1D529CFC}, // 1e235
	{0x169840EF017DA3B1, 0xFB9B7CD9A4A7443C}, // 1e236
	{0x8E1F289560EE864E, 0x9D412E0806E88AA5}, // 1e237
	{0xF1A6F2BAB92A27E2, 0xC491798A08A2AD4E}, // 1e238
	{0xAE10AF696774B1DB, 0xF5B5D7EC8ACB58A2}, // 1e239
	{0xACCA6DA1E0A8EF29, 0x9991A6F3D6BF1765}, // 1e240
	{0x17FD090A58D32AF3, 0xBFF610B0CC6EDD3F}, // 1e241
	{0xDDFC4B4CEF07F5B0, 0xEFF394DCFF8A948E}, // 1e242
	{0x4ABDAF101564F98E, 0x95F83D0A1FB69CD9}, // 1e243
	{0x9D6D1AD41ABE37F1, 0xBB764C4CA7A4440F}, // 1e244
	{0x84C86189216DC5ED, 0xEA53DF5FD18D5513}, // 1e245
	{0x32FD3CF5B4E49BB4, 0x92746B9BE2F8552C}, // 1e246
	{0x3FBC8C33221DC2A1, 0xB7118682DBB66A77}, // 1e247
	{0x0FABAF3FEAA5334A, 0xE4D5E82392A40515}, // 1e248
	{0x29CB4D87F2A7400E, 0x8F05B1163BA6832D}, // 1e249
	{0x743E20E9EF511012, 0xB2C71D5BCA9023F8}, // 1e250
	{0x914DA9246B255416, 0xDF78E4B2BD342CF6}, // 1e251
	{0x1AD089B6C2F7548E, 0x8BAB8EEFB6409C1A}, // 1e252
	{0xA184AC2473B529B1, 0xAE9672ABA3D0C320}, // 1e253
	{0xC9E5D72D90A2741E, 0xDA3C0F568CC4F3E8}, // 1e254
	{0x7E2FA67C7A658892, 0x8865899617FB1871}, // 1e255
	{0xDDBB901B98FEEAB7, 0xAA7EEBFB9DF9DE8D}, // 1e256
	{0x552A74227F3EA565, 0xD51EA6FA85785631}, // 1e257
	{0xD53A88958F87275F, 0x8533285C936B35DE}, // 1e258
	{0x8A892ABAF368F137, 0xA67FF273B8460356}, // 1e259
	{0x2D2B7569B0432D85, 0xD01FEF10A657842C}, // 1e260
	{0x9C3B29620E29FC73, 0x8213F56A67F6B29B}, // 1e261
	{0x8349F3BA91B47B8F, 0xA298F2C501F45F42}, // 1e262
	{0x241C70A936219A73, 0xCB3F2F7642717713}, // 1e263
	{0xED238CD383AA0110, 0xFE0EFB53D30DD4D7}, // 1e264
	{0xF4363804324A40AA, 0x9EC95D1463E8A506}, // 1e265
	{0xB143C6053EDCD0D5, 0xC67BB4597CE2CE48}, // 1e266
	{0xDD94B7868E94050A, 0xF81AA16FDC1B81DA}, // 1e267
	{0xCA7CF2B4191C8326, 0x9B10A4E5E9913128}, // 1e268
	{0xFD1C2F611F63A3F0, 0xC1D4CE1F63F57D72}, // 1e269
	{0xBC633B39673C8CEC, 0xF24A01A73CF2DCCF}, // 1e270
	{0xD5BE0503E085D813, 0x976E41088617CA01}, // 1e271
	{0x4B2D8644D8A74E18, 0xBD49D14AA79DBC82}, // 1e272
	{0xDDF8E7D60ED1219E, 0xEC9C459D51852BA2}, // 1e273
	{0xCABB90E5C942B503, 0x93E1AB8252F33B45}, // 1e274
	{0x3D6A751F3B936243, 0xB8DA1662E7B00A17}, // 1e275
	{0x0CC512670A783AD4, 0xE7109BFBA19C0C9D}, // 1e276
	{0x27FB2B80668B24C5, 0x906A617D450187E2}, // 1e277
	{0xB1F9F660802DEDF6, 0xB484F9DC9641E9DA}, // 1e278
	{0x5E7873F8A0396973, 0xE1A63853BBD26451}, // 1e279
	{0xDB0B487B6423E1E8, 0x8D07E33455637EB2}, // 1e280
	{0x91CE1A9A3D2CDA62, 0xB049DC016ABC5E5F}, // 1e281
	{0x7641A140CC7810FB, 0xDC5C5301C56B75F7}, // 1e282
	{0xA9E904C87FCB0A9D, 0x89B9B3E11B6329BA}, // 1e283
	{0x546345FA9FBDCD44, 0xAC2820D9623BF429}, // 1e284
	{0xA97C177947AD4095, 0xD732290FBACAF133}, // 1e285
	{0x49ED8EABCCCC485D, 0x867F59A9D4BED6C0}, // 1e286
	{0x5C68F256BFFF5A74, 0xA81F301449EE8C70}, // 1e287
	{0x73832EEC6FFF3111, 0xD226FC195C6A2F8C}, // 1e288
	{0xC831FD53C5FF7EAB, 0x83585D8FD9C25DB7}, // 1e289
	{0xBA3E7CA8B77F5E55, 0xA42E74F3D032F525}, // 1e290
	{0x28CE1BD2E55F35EB, 0xCD3A1230C43FB26F}, // 1e291
	{0x7980D163CF5B81B3, 0x80444B5E7AA7CF85}, // 1e292
	{0xD7E105BCC332621F, 0xA0555E361951C366}, // 1e293
	{0x8DD9472BF3FEFAA7, 0xC86AB5C39FA63440}, // 1e294
	{0xB14F98F6F0FEB951, 0xFA856334878FC150}, // 1e295
	{0x6ED1BF9A569F33D3, 0x9C935E00D4B9D8D2}, // 1e296
	{0x0A862F80EC4700C8, 0xC3B8358109E84F07}, // 1e297
	{0xCD27BB612758C0FA, 0xF4A642E14C6262C8}, // 1e298
	{0x8038D51CB897789C, 0x98E7E9CCCFBD7DBD}, // 1e299
	{0xE0470A63E6BD56C3, 0xBF21E44003ACDD2C}, // 1e300
	{0x1858CCFCE06CAC74, 0xEEEA5D5004981478}, // 1e301
	{0x0F37801E0C43EBC8, 0x95527A5202DF0CCB}, // 1e302
	{0xD30560258F54E6BA, 0xBAA718E68396CFFD}, // 1e303
	{0x47C6B82EF32A2069, 0xE950DF20247C83FD}, // 1e304
	{0x4CDC331D57FA5441, 0x91D28B7416CDD27E}, // 1e305
	{0xE0133FE4ADF8E952, 0xB6472E511C81471D}, // 1e306
	{0x58180FDDD97723A6, 0xE3D8F9E563A198E5}, // 1e307
	{0x570F09EAA7EA7648, 0x8E679C2F5E44FF8F}, // 1e308
	{0x2CD2CC6551E513DA, 0xB201833B35D63F73}, // 1e309
	{0xF8077F7EA65E58D1, 0xDE81E40A034BCF4F}, // 1e310
	{0xFB04AFAF27FAF782, 0x8B112E86420F6191}, // 1e311
	{0x79C5DB9AF1F9B563, 0xADD57A27D29339F6}, // 1e312
	{0x18375281AE7822BC, 0xD94AD8B1C7380874}, // 1e313
	{0x8F2293910D0B15B5, 0x87CEC76F1C830548}, // 1e314
	{0xB2EB3875504DDB22, 0xA9C2794AE3A3C69A}, // 1e315
	{0x5FA60692A46151EB, 0xD433179D9C8CB841}, // 1e316
	{0xDBC7C41BA6BCD333, 0x849FEEC281D7F328}, // 1e317
	{0x12B9B522906C0800, 0xA5C7EA73224DEFF3}, // 1e318
	{0xD768226B34870A00, 0xCF39E50FEAE16BEF}, // 1e319
	{0xE6A1158300D46640, 0x81842F29F2CCE375}, // 1e320
	{0x60495AE3C1097FD0, 0xA1E53AF46F801C53}, // 1e321
	{0x385BB19CB14BDFC4, 0xCA5E89B18B602368}, // 1e322
	{0x46729E03DD9ED7B5, 0xFCF62C1DEE382C42}, // 1e323
	{0x6C07A2C26A8346D1, 0x9E19DB92B4E31BA9}, // 1e324
	{0xC7098B7305241885, 0xC5A05277621BE293}, // 1e325
	{0xB8CBEE4FC66D1EA7, 0xF70867153AA2DB38}, // 1e326
	{0x737F74F1DC043328, 0x9A65406D44A5C903}, // 1e327
	{0x505F522E53053FF2, 0xC0FE908895CF3B44}, // 1e328
	{0x647726B9E7C68FEF, 0xF13E34AABB430A15}, // 1e329
	{0x5ECA783430DC19F5, 0x96C6E0EAB509E64D}, // 1e330
	{0xB67D16413D132072, 0xBC789925624C5FE0}, // 1e331
	{0xE41C5BD18C57E88F, 0xEB96BF6EBADF77D8}, // 1e332
	{0x8E91B962F7B6F159, 0x933E37A534CBAAE7}, // 1e333
	{0x723627BBB5A4ADB0, 0xB80DC58E81FE95A1}, // 1e334
	{0xCEC3B1AAA30DD91C, 0xE61136F2227E3B09}, // 1e335
	{0x213A4F0AA5E8A7B1, 0x8FCAC257558EE4E6}, // 1e336
	{0xA988E2CD4F62D19D, 0xB3BD72ED2AF29E1F}, // 1e337
	{0x93EB1B80A33B8605, 0xE0ACCFA875AF45A7}, // 1e338
	{0xBC72F130660533C3, 0x8C6C01C9498D8B88}, // 1e339
	{0xEB8FAD7C7F8680B4, 0xAF87023B9BF0EE6A}, // 1e340
	{0xA67398DB9F6820E1, 0xDB68C2CA82ED2A05}, // 1e341
	{0x88083F8943A1148C, 0x892179BE91D43A43}, // 1e342
	{0x6A0A4F6B948959B0, 0xAB69D82E364948D4}, // 1e343
	{0x848CE34679ABB01C, 0xD6444E39C3DB9B09}, // 1e344
	{0xF2D80E0C0C0B4E11, 0x85EAB0E41A6940E5}, // 1e345
	{0x6F8E118F0F0E2195, 0xA7655D1D2103911F}, // 1e346
	{0x4B7195F2D2D1A9FB, 0xD13EB46469447567}, // 1e347
}
