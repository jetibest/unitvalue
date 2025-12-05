package unitvalue

import (
	"fmt"
	"strconv"
	"strings"
	"math"
	"unicode"
	"golang.org/x/exp/utf8string"
)

// note: UnitValue is also a UnitDefinition, in which case uv.Offset != nil, at that time, the uv.Value is only used for conversion (in combination with uv.Offset)
// so while only a UnitValue struct is defined, UnitDefinition is basically a UnitValue where Offset != nil

var (
	// without dimensions
	Dimless = NewUnitDefinition(1, 0, UnitDimensions{})
	Percent = NewUnitDefinition(0.01, 0, UnitDimensions{})
	PartsPerCentMille = NewUnitDefinition(0.00001, 0, UnitDimensions{})
	PartsPerMillion = NewUnitDefinition(0.000001, 0, UnitDimensions{})
	PartsPerBillion = NewUnitDefinition(0.000000001, 0, UnitDimensions{})
	PartsPerTrillion = NewUnitDefinition(0.000000000001, 0, UnitDimensions{})
	PartsPerQuadrillion = NewUnitDefinition(0.000000000000001, 0, UnitDimensions{})
	
	// scientific base units
	Metre = NewUnitDefinition(1, 0, UnitDimensions{"m": 1})
	Length = Metre
	Meter = Metre
	Kilogram = NewUnitDefinition(1, 0, UnitDimensions{"kg": 1})
	Mass = Kilogram
	Second = NewUnitDefinition(1, 0, UnitDimensions{"s": 1})
	Time = Second
	Ampere = NewUnitDefinition(1, 0, UnitDimensions{"A": 1})
	Current = Ampere
	AbsoluteKelvin = NewUnitDefinition(1, 0, UnitDimensions{"K": 1})
	Temperature = AbsoluteKelvin
	Mole = NewUnitDefinition(1, 0, UnitDimensions{"mol": 1})
	Candela = NewUnitDefinition(1, 0, UnitDimensions{"cd": 1})
	LuminousIntensity = Candela
	
	// common units derived from base units
	Gram = NewUnitDefinition(0.001, 0, UnitDimensions{"kg": 1})
	Coulomb = NewUnitDefinition(1, 0, UnitDimensions{"s": 1, "A": 1})
	Hertz = NewUnitDefinition(1, 0, UnitDimensions{"s": -1})
	Newton = NewUnitDefinition(1, 0, UnitDimensions{"kg": 1, "m": 1, "s": -2})
	Litre = NewUnitDefinition(0.001, 0, UnitDimensions{"m": 3})
	Liter = Litre
	Hour = NewUnitDefinition(3600, 0, UnitDimensions{"s": 1})
	Voltage = NewUnitDefinition(1, 0, UnitDimensions{"kg": 1, "m": 2, "s": -3, "A": -1})
	Joule = NewUnitDefinition(1, 0, UnitDimensions{"kg": 1, "m": 2, "s": -2})
	Watt = NewUnitDefinition(1, 0, UnitDimensions{"kg": 1, "m": 2, "s": -3})
	WattHour = NewUnitDefinition(3600, 0, UnitDimensions{"kg": 1, "m": 2, "s": -2})
	VoltAmpere = NewUnitDefinition(1, 0, UnitDimensions{"kg": 1, "m": 2, "s": -3})
	AbsoluteCelsius = NewUnitDefinition(1, 273.15, UnitDimensions{"K": 1})
	AbsoluteFahrenheit = NewUnitDefinition(5.0/9.0, 459.67, UnitDimensions{"K": 1})
	Mile = NewUnitDefinition(1609.344, 0, UnitDimensions{"m": 1})
	Inch = NewUnitDefinition(0.0254, 0, UnitDimensions{"m": 1})
	Foot = NewUnitDefinition(0.3048, 0, UnitDimensions{"m": 1})
	Yard = NewUnitDefinition(0.9144, 0, UnitDimensions{"m": 1})
	
	// custom units (not in Base SI)
	Byte = NewUnitDefinition(1, 0, UnitDimensions{"B": 1})
	Bit = NewUnitDefinition(1.0/8.0, 0, UnitDimensions{"B": 1})
	Euro = NewUnitDefinition(1, 0, UnitDimensions{"EUR": 1}) // because currencies are very volatile, and there is no base currency to compare against as it depends on your perspective (although USD is common), this should probably be defined outside of the library as a custom UnitDefinition
	Radian = NewUnitDefinition(1, 0, UnitDimensions{"m": 1-1}) // in base SI it should be m/m, which is unitless
	Steradian = NewUnitDefinition(1, 0, UnitDimensions{"m": 2-2}) // in base SI it should be m2/m2, which is unitless
	Lumen = NewUnitDefinition(1, 0, UnitDimensions{"cd": 1, "m": 2-2})
	Lux = NewUnitDefinition(1, 0, UnitDimensions{"cd": 1, "m": 2-2 - 2})
	Degree = NewUnitDefinition(math.Pi / 180, 0, UnitDimensions{"m": 1-1})
	RelativeKelvin = NewUnitDefinition(1, 0, UnitDimensions{"ΔK": 1})
	RelativeCelsius = RelativeKelvin
	RelativeFahrenheit = NewUnitDefinition(1, 0, UnitDimensions{"ΔF": 1})
)

var Units = map[string]*UnitValue{
	"": Dimless,
	"%": Percent,
	"pcm": PartsPerCentMille,
	"ppm": PartsPerMillion,
	"ppb": PartsPerBillion,
	"ppt": PartsPerTrillion,
	"ppq": PartsPerQuadrillion,
	
	// NOTE: these units below, should be logarithmic in their mathematics, not just dimless, but should be based off of dB, but with the given logarithmic-unit properties
	// factor of power, energy, or mass/time(-3) is typically factor=10, otherwise it's factor=20 (for amplitude quantities)
	"dB": Dimless, // {base(refUnit): '', ref(refLevel): 1, refFactor: 20, factor: 10, type: 'generic'}
	"dBV": Dimless, // {refUnit: 'V', refLevel: 1, refFactor: 20, factor: 20, type: 'voltage'}
	"dBm": Dimless, // {refUnit: 'W', refLevel: 1e-3, refFactor: 10, factor: 10, type: 'power'}
	"dBW": Dimless, // {refUnit: 'W', refLevel: 1, refFactor: 10, factor: 10, type: 'power'}
	"dBA": Dimless, // {refUnit: 'Pa', refLevel: 20e-6, refFactor: 20, factor: 10, type: 'sound'}
	"dBuV": Dimless, // {refUnit: 'V', refLevel: 1e-6, refFactor: 20, factor: 20, type: 'voltage'}
	"pH": Dimless, // {refUnit: 'mol/L', refLevel: 1, refFactor: -1, factor: -1, type: 'acidity'}
	
	"CO2": Dimless,
	"CO2e": Dimless,
	"CO₂": Dimless,
	"CO₂e": Dimless,
	
	"m": Metre,
	"kg": Kilogram,
	"s": Second,
	"A": Ampere,
	"K": AbsoluteKelvin,
	"mol": Mole,
	"cd": Candela,
	
	"g": Gram,
	"C": Coulomb,
	"Hz": Hertz,
	"N": Newton,
	"l": Litre, // lowercase 'l' is the official spelling
	"L": Litre, // but L is equally valid for clarity due to confusion with '1'
	"h": Hour,
	"V": Voltage,
	"J": Joule,
	"W": Watt,
	"Wh": WattHour,
	"VA": VoltAmpere, // Apparent Power
	"var": VoltAmpere, // Volt-Ampere Reactive, or Reactive Power (lowercase 'var' is the official spelling)
	"VAr": VoltAmpere,
	"ΔK": RelativeKelvin,
	"°C": AbsoluteCelsius,
	"℃": AbsoluteCelsius,
	"Δ°C": RelativeCelsius,
	"Δ℃": RelativeCelsius,
	
	"in": Inch,
	"ft": Foot,
	"yd": Yard,
	"mi": Mile,
	"mph": Mile.Divide(Hour).AsDefinition(),
	
	"B": Byte,
	"Bps": Byte.Divide(Second).AsDefinition(),
	"bps": Bit.Divide(Second).AsDefinition(),
	"bit": Bit,
	"EUR": Euro,
	"rad": Radian,
	"sr": Steradian,
	"lm": Lumen,
	"lx": Lux,
	"deg": Degree,
	"°": Degree,
}

type UnitDimensions = map[string]float64

type UnitValue struct {
	Value float64
	Offset *float64 // defaults to nil, if not 0, then the unit is a unit definition, is only used during parsing or formatting, will be immediately set to zero to make Value remain consistent as a true value in base units
	Unit string // useful when doing abs(-10 °C) which would otherwise not result in 10 °C
	Dimensions UnitDimensions
}

func NewUnitDefinition(slope float64, bias float64, dims UnitDimensions) *UnitValue {
	
	return &UnitValue{Value: slope, Offset: &bias, Dimensions: dims}
}

func FormatUnitDimensions(dims UnitDimensions) string {
	
	// return a unit representation in the highest known order
	s := ""
	
	// first handle the positive ones
	for dim, exp := range dims {
		
		if exp > 0 {
			
			if len(s) != 0 {
				s += "*"
			}
			
			if exp == 1 {
				s += fmt.Sprintf("%s", dim)
			} else {
				s += fmt.Sprintf("%s^%.16g", dim, exp)
			}
		}
	}
	// then the negative ones
	for dim, exp := range dims {
		
		if exp < 0 {
			
			if len(s) != 0 {
				s += "/"
				exp = exp * -1
			}
			
			if exp == 1 {
				s += fmt.Sprintf("%s", dim)
			} else {
				s += fmt.Sprintf("%s^%.16g", dim, exp)
			}
		}
	}
	// ignore if equal to 0
	
	return s
}

// expects a float parseable value, with some optional unit behind it
func Parse(s string) (*UnitValue, error) {
	
	var v float64
	var u string
	n, err := fmt.Sscanf(s, "%f%s", &v, &u)
	if err != nil {
		if n == 1 {
			// the unit is empty
			u = ""
		} else {
			return nil, err
		}
	}
	u = strings.TrimSpace(u)
	
	return ParseUnitWithValue(v, u)
}

func ParseUnitWithValue(v float64, s string) (*UnitValue, error) {
	
	uv, err := ParseUnit(s)
	if err != nil {
		return nil, err
	}
	return uv.From(v), nil
}

func ParseUnit(s string) (*UnitValue, error) {
	
	// s => "kg * m^2 / s^-3", try to parse multiple units in one "combo"
	var combined_uv *UnitValue = nil
	var identifier strings.Builder
	operator := '?'
	for _, c := range s {
		
		if c == ' ' || c == '*' || c == '/' {
			
			if identifier.Len() != 0 {
				
				uv, err := ParseUnit(identifier.String())
				if err != nil {
					return nil, err
				}
				if combined_uv == nil {
					combined_uv = uv
				} else {
					if operator == '/' {
						combined_uv = combined_uv.Divide(uv)
					} else {
						combined_uv = combined_uv.Multiply(uv)
					}
				}
				identifier.Reset()
				
				if c == '*' {
					operator = '*'
				} else if c == '/' {
					operator = '/'
				} else if operator == '?' {
					operator = '*'
				}
			}
			
		} else {
			identifier.WriteRune(c)
		}
	}
	if combined_uv != nil {
		
		if identifier.Len() != 0 {
			uv, err := ParseUnit(identifier.String())
			if err != nil {
				return nil, err
			}
			if operator == '/' {
				combined_uv = combined_uv.Divide(uv)
			} else {
				combined_uv = combined_uv.Multiply(uv)
			}
		}
		
		// store original Unit string (after one operation where the base units change, this Unit will be lost)
		combined_uv.Unit = s
		
		return combined_uv.AsDefinition(), nil
	}
	
	// quick look-up map of base units, custom units, and common combinations (for performance)
	if mappedUnit, ok := Units[s]; ok {
		uv := mappedUnit.Copy()
		uv.Unit = s
		return uv, nil
	}
	
	utf8str := utf8string.NewString(s)
	rc := utf8str.RuneCount()
	
	// handle most common prefixes (both Mi for Byte and M for SI)
	var factor float64 = 1
	prefix := ""
	
	if rc > 1 {
		
		// NOTE: checking utf8next because if we have "m^2" or "m2" then "m" should not be considered a prefix because there is no letter following it, and thus must be interpreted as meter
		// however, if mm^2, then the first m is milli, and the second is the meter
		
		if rc > 2 {
			utf8prefix := utf8str.Slice(0, 2)
			utf8next := utf8str.At(2)
			
			if unicode.IsLetter(rune(utf8next)) {
				
				switch utf8prefix {
					case "Ki":
						factor = 1024
					case "Mi":
						factor = 1024 * 1024
					case "Gi":
						factor = 1024 * 1024 * 1024
					case "Ti":
						factor = 1024 * 1024 * 1024 * 1024
					case "Pi":
						factor = 1024 * 1024 * 1024 * 1024 * 1024
					case "Ei":
						factor = 1024 * 1024 * 1024 * 1024 * 1024 * 1024
					case "da":
						factor = 1e1
					case "µ":
						factor = 1e-6
				}
				
				// update values
				if factor != 1 {
					prefix = utf8prefix
					s = s[len(prefix):]
					utf8str = utf8string.NewString(s)
					rc = utf8str.RuneCount()
				}
			}
		}
		// if not yet trimmed a prefix
		if factor == 1 {
			
			utf8prefix := utf8str.Slice(0, 1)
			utf8next := utf8str.At(1)
			
			if unicode.IsLetter(rune(utf8next)) {
				
				switch utf8prefix {
					case "n":
						factor = 1e-9
					case "u":
						factor = 1e-6
					case "μ":
						factor = 1e-6
					case "m":
						factor = 1e-3
					case "c":
						factor = 1e-2
					case "d":
						factor = 1e-1
					case "h":
						factor = 1e2
					case "k":
						factor = 1e3
					case "M":
						factor = 1e6
					case "G":
						factor = 1e9
					case "T":
						factor = 1e12
					case "P":
						factor = 1e15
					case "E":
						factor = 1e18
				}
				
				if factor != 1 {
					
					prefix = utf8prefix
					s = s[len(prefix):]
					utf8str = utf8string.NewString(s)
					rc = utf8str.RuneCount()
				}
			}
		}
	}
	
	if factor != 1 {
		v, err := ParseUnit(s)
		if err != nil {
			return nil, err
		}
		v.Unit = prefix + v.UnitString()
		return v.Multiply(Dimless.New(factor)).AsDefinition(), nil
	}
	
	// power is only allowed if no further / or * exist
	// check for ^x suffix
	i := strings.LastIndex(s, "^")
	if i != -1 {
		
		u, err := ParseUnit(s[0:i])
		if err != nil {
			return nil, err
		}
		i += 1
		power, err := strconv.ParseFloat(s[i:], 64)
		if err != nil {
			return nil, err
		}
		u.Unit = fmt.Sprintf("%s^%.16g", u.Unit, power)
		res, err := u.Power(power)
		if err != nil {
			return nil, err
		}
		return res.AsDefinition(), nil
		
	} else if strings.HasSuffix(s, "2") || strings.HasSuffix(s, "²") { // e.g. "m2"
		
		u, err := ParseUnit(strings.TrimSuffix(strings.TrimSuffix(s, "²"), "2"))
		if err != nil {
			return nil, err
		}
		i += 1
		u.Unit = fmt.Sprintf("%s^2", u.Unit)
		res, err := u.Power(2.0)
		if err != nil {
			return nil, err
		}
		return res.AsDefinition(), nil
		
	} else if strings.HasSuffix(s, "3") || strings.HasSuffix(s, "³") { // e.g. "m3"
		
		u, err := ParseUnit(strings.TrimSuffix(strings.TrimSuffix(s, "³"), "3"))
		if err != nil {
			return nil, err
		}
		i += 1
		u.Unit = fmt.Sprintf("%s^3", u.Unit)
		res, err := u.Power(3.0)
		if err != nil {
			return nil, err
		}
		return res.AsDefinition(), nil
	}
	
	// try to turn 'p' into '/' as it is common to write b/s as bps
	parts := strings.Split(s, "p")
	if len(parts) == 2 {
		return ParseUnit(parts[0] + "/" + parts[1])
	}
	
	// all possibilities exhausted, failed to parse unit
	return nil, fmt.Errorf("unsupported unit (%s)", s)
}

// copies UnitValue and ensures Offset is defined as a non-nil value, which makes it a UnitDefinition
func (uv *UnitValue) AsDefinition() *UnitValue {
	
	newDims := UnitDimensions{}
	
	for k, v := range uv.Dimensions {
		
		if v == 0 {
			continue
		}
		
		newDims[k] = v
	}
	
	newOffset := 0.0
	if uv.Offset != nil {
		newOffset = *uv.Offset
	}
	
	return &UnitValue{
		Value: uv.Value,
		Offset: &newOffset,
		Unit: uv.Unit,
		Dimensions: newDims,
	}
}

// true if any non-zero dimensions exist
func (uv *UnitValue) HasUnitDimensions() bool {
	
	c := 0
	for _, exp := range uv.Dimensions {
		
		if exp == 0 {
			continue
		}
		
		c += 1
	}
	
	return c > 0
}

// From is only supposed to be used by the const declared units! where value is valid for the given unit as its context
func (uv *UnitValue) From(value float64) *UnitValue {
	
	newDims := UnitDimensions{}
	
	for k, v := range uv.Dimensions {
		
		if v == 0 {
			continue
		}
		
		newDims[k] = v
	}
	
	offset := 0.0
	if uv.Offset != nil {
		offset = *uv.Offset
	}
	
	return &UnitValue{
		Value: uv.Value * (value + offset),
		Unit: uv.Unit,
		Dimensions: newDims,
	}
}

// note: value here is the real underlying value expressed in its base units!
func (uv *UnitValue) New(value float64) *UnitValue {
	
	newDims := UnitDimensions{}
	
	for k, v := range uv.Dimensions {
		
		if v == 0 {
			continue
		}
		
		newDims[k] = v
	}
	
	return &UnitValue{
		Value: value,
		Unit: uv.Unit,
		Dimensions: newDims,
	}
}

func (uv *UnitValue) Copy() *UnitValue {
	
	newDims := UnitDimensions{}
	
	for k, v := range uv.Dimensions {
		
		if v == 0 {
			continue
		}
		
		newDims[k] = v
	}
	
	var newOffsetPtr *float64
	if uv.Offset != nil {
		newOffset := *uv.Offset
		newOffsetPtr = &newOffset
	}
	
	return &UnitValue{
		Value: uv.Value,
		Offset: newOffsetPtr,
		Unit: uv.Unit,
		Dimensions: newDims,
	}
}

func (uv *UnitValue) UnitString() string {
	
	if uv.Unit != "" {
		
		return uv.Unit
	}
	
	return uv.FormatUnit()
}

// express the UnitValue in target unit string, and return the value that would only make sense in terms of the given target unit string
func (uv *UnitValue) GetValueInUnitString(target string) (float64, error) {
	
	target_unit, err := ParseUnit(target)
	if err != nil {
		return 0.0, err
	}
	
	return uv.GetValueInUnit(target_unit)
}

// express the UnitValue in target UnitDefinition, and return the value that would only make sense in terms of the given target unit definition
func (uv *UnitValue) GetValueInUnit(target_unit *UnitValue) (float64, error) {
	
	if target_unit.Offset == nil {
		// this is not a UnitDefinition but a UnitValue, which is not allowed (invalid usage!)
		return 0.0, fmt.Errorf("Conversion error: target_unit is not a UnitDefinition (%+v)", target_unit)
	}
	
	if ! uv.MatchesUnit(target_unit) {
		return 0.0, fmt.Errorf("Conversion error: dimensions do not match (%+v vs %+v)", uv, target_unit)
	}
	
	// for example: let's say target_unit is kWh, and the value is 1000 instead of 1, then we must divide by 1000 to get the Wh into kWh
	// or say, the value is 293.15 K, and we want to convert to °C, then we end up with 20.0 °C
	// or if 20 °C and we want to convert to K, then 293.15 / 1 - 0 = 293.15 K
	
	return uv.Value / target_unit.Value - *target_unit.Offset, nil
}

// format the unit dimensions (raw)
func (uv *UnitValue) FormatUnit() string {
	
	return FormatUnitDimensions(uv.Dimensions)
}

// format the UnitValue or UnitDefinition, but use the original Unit if stored in the struct (overriding the default raw base units format)
func (uv *UnitValue) Format() string {
	
	if uv.Unit != "" {
		
		// don't print the value at all, if it's a UnitDefinition (so the Value is just used as a conversion factor, and not actually representing a real value)
		if uv.Offset != nil {
			return uv.Unit
		}
		
		cv, err := uv.GetValueInUnitString(uv.Unit)
		if err == nil {
			return fmt.Sprintf("%.16g %s", cv, uv.Unit)
		}
	}
	
	return uv.String()
}

// string representation of the UnitValue or UnitDefinition
// note: do not use String(), but use Format() instead:
// because uv.Unit depends on GetValueInUnitString, which can also print the value with .String(), we cannot use uv.Unit in String(), and thus Format() should be preferred
func (uv *UnitValue) String() string {
	
	s := ""
	
	// don't print the value at all, if it's a UnitDefinition (so the Value is just used as a conversion factor, and not actually representing a real value)
	if uv.Offset == nil {
		s += fmt.Sprintf("%.16g", uv.Value)
	}
	
	ustr := uv.FormatUnit()
	if len(ustr) > 0 {
		
		if len(s) > 0 {
		
			s += " "
		}
		s += ustr
	}
	
	return s
}

// check if dimensions match, however, any zero-power dimensions (^0) are ignored
func (uv *UnitValue) MatchesUnit(x *UnitValue) bool {
	
	for dim, uv_exp := range uv.Dimensions {
		
		if uv_exp == 0 {
			continue
		}
		
		x_exp, ok := x.Dimensions[dim]
		if ! ok || uv_exp != x_exp {
			return false
		}
	}
	for dim, x_exp := range x.Dimensions {
		
		if x_exp == 0 {
			continue
		}
		
		uv_exp, ok := uv.Dimensions[dim]
		if ! ok || x_exp != uv_exp {
			return false
		}
	}
	
	return true
}

// do UnitValue to the power of a number, returns an error if UnitValue is a UnitDefinition unless the offset is zero in which case it is allowed
func (uv *UnitValue) Power(x float64) (*UnitValue, error) {
	
	// if exponent x is a logarithmic unit, then we just apply the regular power, without considering x's unit, just using the number value
	// if base is a logarithmic unit, then we convertToLinear(base) before power(), and then do convertFromLinear() after the power is done, and reapplying the unit of the base as it were
	// similarly, with sqrt, abs, log10, etc. we convertToLinear, apply the function, and then convertFromLinear and reapply the unit (except for log10, we don't do convertFromLinear, and don't apply any unit anymore)
	
	newDims := UnitDimensions{}
	
	for dim, uv_exp := range uv.Dimensions {
		
		if uv_exp == 0 {
			continue
		}
		
		newDims[dim] = uv_exp * x
	}
	
	// any UnitValue with non-zero offset is a UnitDefinition with the formula f(x) = Value*(x + Offset), and does not yet have a real value!
	if uv.Offset != nil {
		// remember, °C^2 is not a valid unit, there is no such thing.. however, we'll allow any ^2 on unit definitions with a zero-offset
		
		if *uv.Offset != 0 {
			return nil, fmt.Errorf("unsupported power (%.16g) for unit (%s)", x, uv.UnitString())
		} else {
			// (K)^2 = K^2, the uv.Value is not the actual Value, this is a UnitDefinition, not a UnitValue!
			
			return &UnitValue{
				Value: uv.Value,
				Dimensions: newDims,
			}, nil
		}
	}
	
	return &UnitValue{
		Value: math.Pow(uv.Value, x),
		Dimensions: newDims,
	}, nil
}

// shared function for multiplication and division depending on op
// applies Offset when a UnitDefinition is multiplied by a Dimless number for the first time
// tries to preserve Unit on either UnitValue if present
func multiply_unitvalues(x *UnitValue, y *UnitValue, op float64) *UnitValue {
	
	// when x and y are logarithmic units, that's an error.. but if one of them is, it's just a regular multiplier factor, where the number value is multiplied, and the unit remains intact
	// however, when dividing x/y, and both are logarithmic, we actually return convertToLinear(x)/convertToLinear(y), as a dimless value
	
	newDims := UnitDimensions{}
	
	x_dims := 0
	y_dims := 0
	changed := 0
	
	for dim, x_exp := range x.Dimensions {
		
		if x_exp == 0 {
			continue
		}
		
		x_dims += 1
		
		changed = 1
		
		newDims[dim] = x_exp
	}
	
	for dim, y_exp := range y.Dimensions {
		
		if y_exp == 0 {
			continue
		}
		
		y_dims += 1
		
		if changed == 1 {
			changed = 2
		}
		
		e, ok := newDims[dim]
		if ! ok {
			e = 0
		}
		newDims[dim] = e + op * y_exp
	}
	
	// only use Unit for explicitly user-input defined unit formats
	// preserve unit and offset if possible
	
	// x is dimless, so we can preserve the unit of y
	newUnit := y.Unit
	if changed == 1 {
		// y is dimless, so we can preserve the unit of x
		newUnit = x.Unit
	} else if changed == 2 {
		// neither x nor y is dimless
		newUnit = ""
	}
	
	x_val := x.Value
	y_val := y.Value
	
	// apply offset when converting to a number, when converting UnitDefinition for the first time
	if x_dims != 0 && y_dims == 0 {
		if x.Offset != nil {
			y_val += *x.Offset
		}
	} else if x_dims == 0 && y_dims != 0 {
		if y.Offset != nil {
			x_val += *y.Offset
		}
	}
	
	return &UnitValue{
		Value: x_val * math.Pow(y_val, op), // op=1 for multiplication, op=-1 for division
		Unit: newUnit,
		Dimensions: newDims,
	}
}

// shared function for adding or subtracting unitvalues based on op
// cannot add/subtract UnitDefinition, this will throw an error
// tries to preserve Unit if present
// unit dimensions must be exactly equal, with one exception, Δ is allowed as a prefix on a base-unit as a relative unit, in which case e.g. K + ΔK = K
func add_unitvalues(x *UnitValue, y *UnitValue, op float64) (*UnitValue, error) {
	
	// if x and y are both log units
	// if x.type is generic, then convert x to specific y's unit
	// if y.type is generic, then convert y to specific x's unit
	// if either x or y is type generic OR base and factor of both x and y are both equal, then we can succeed, otherwise incompatible logarithmic units
	// the resulting unitvalue is convertFromLinear(convertToLinear(x) + convertToLinear(y), unit(x))
	// where convertToLinear(x) = x.ref * pow(10, x.number / x.factor)
	// where convertFromLinear(x, unit) = x.factor * log10(linearValue / x.ref)
	
	if x.Offset != nil || y.Offset != nil {
		opstr := "+"
		if op < 0 {
			opstr = "-"
		}
		return nil, fmt.Errorf("unit definition cannot be added or subtracted (%+v %s %+v)", x.Format(), opstr, y.Format())
	}
	
	// we cannot mix unit, choose one
	newUnit := x.Unit
	if newUnit == "" {
		newUnit = y.Unit
	}
	
	newDims := UnitDimensions{}
	
	for dim, x_exp := range x.Dimensions {
		
		if x_exp == 0 {
			continue
		}
		
		y_exp, ok := y.Dimensions[dim]
		if ! ok || x_exp != y_exp {
			
			// if dim is ΔK, then it may be equivalent to y.Dimensions[K], in which case the new unit must be absolute
			if strings.HasPrefix(dim, "Δ") {
				
				tmpdim := strings.TrimPrefix(dim, "Δ")
				y_exp, ok = y.Dimensions[tmpdim]
				if ok && x_exp == y_exp {
					
					// it's impossible to preserve unit here (because what if there is a delta in x, and another delta in y?)
					newUnit = ""
					
					newDims[tmpdim] = y_exp
					
					continue
				}
			}
			
			// or it's the other way around
			y_exp, ok = y.Dimensions["Δ" + dim]
			if ok && x_exp == y_exp {
				
				// it's impossible to preserve unit here (because what if there is a delta in x, and another delta in y?)
				newUnit = ""
				
				newDims[dim] = x_exp
				
				continue
			}
			
			return nil, fmt.Errorf("dimensions mismatch for [%+v].AddOrSub([%+v])", x, y)
		}
		
		newDims[dim] = x_exp
	}
	for dim, y_exp := range y.Dimensions {
		
		if y_exp == 0 {
			continue
		}
		
		x_exp, ok := x.Dimensions[dim]
		if ! ok || y_exp != x_exp {
			
			// if dim is ΔK, then it may be equivalent to y.Dimensions[K], in which case the new unit must be absolute
			if strings.HasPrefix(dim, "Δ") {
				
				tmpdim := strings.TrimPrefix(dim, "Δ")
				x_exp, ok = x.Dimensions[tmpdim]
				if ok && y_exp == x_exp {
					
					// it's impossible to preserve unit here (because what if there is a delta in uv, and another delta in x?)
					newUnit = ""
					
					newDims[tmpdim] = y_exp
					
					continue
				}
			}
			
			// or it's the other way around
			x_exp, ok = x.Dimensions["Δ" + dim]
			if ok && x_exp == y_exp {
				
				// it's impossible to preserve unit here (because what if there is a delta in x, and another delta in y?)
				newUnit = ""
				
				newDims[dim] = y_exp
				
				continue
			}
			
			return nil, fmt.Errorf("dimensions mismatch for [%+v].AddOrSub([%+v])", x, y)
		}
		
		newDims[dim] = y_exp
	}
	
	return &UnitValue{
		Value: x.Value + op * y.Value,
		Unit: newUnit,
		Dimensions: newDims,
	}, nil
}

func (uv *UnitValue) Multiply(x *UnitValue) *UnitValue {
	
	return multiply_unitvalues(uv, x, 1)
}

func (uv *UnitValue) Divide(x *UnitValue) *UnitValue {
	
	return multiply_unitvalues(uv, x, -1)
}

func (uv *UnitValue) Add(x *UnitValue) (*UnitValue, error) {
	
	return add_unitvalues(uv, x, 1)
}

func (uv *UnitValue) Subtract(x *UnitValue) (*UnitValue, error) {
	
	return add_unitvalues(uv, x, -1)
}
