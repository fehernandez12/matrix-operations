package matrices

import (
	"errors"
	"fmt"
	"math"
)

const w_app_id = "A3327A-EPTEETJEAA"

type stack []float64

func (s *stack) isEmpty() bool {
	return len(*s) == 0
}

func (s *stack) push(n float64) {
	*s = append(*s, n)
}

func (s *stack) pop() (float64, bool) {
	if s.isEmpty() {
		return 0, false
	}
	i := len(*s) - 1
	n := (*s)[i]
	*s = (*s)[:i]
	return n, true
}

func (s *stack) ToSlice() []float64 {
	return *s
}

func elimination(matrix [][]float64, index []int) {
	n := len(index)
	c := make([]float64, n)
	for i := 0; i < n; i++ {
		index[i] = i
	}

	for i := 0; i < n; i++ {
		c1 := 0.0
		for j := 0; j < n; j++ {
			c0 := math.Abs(matrix[i][j])
			if c0 > c1 {
				c1 = c0
			}
		}
		c[i] = c1
	}

	k := 0
	for j := 0; j < n-1; j++ {
		pivot1 := 0.0
		for i := j; i < n; i++ {
			pivot0 := math.Abs(matrix[index[i]][j])
			pivot0 /= c[index[i]]
			if pivot0 > pivot1 {
				pivot1 = pivot0
				k = i
			}
		}
		itmp := index[j]
		index[j] = index[k]
		index[k] = itmp
		for i := j + 1; i < n; i++ {
			pj := matrix[index[i]][j] / matrix[index[j]][j]
			matrix[index[i]][j] = pj
			for l := j + 1; l < n; l++ {
				matrix[index[i]][l] -= pj * matrix[index[j]][l]
			}
		}
	}
}

func Invert(matrix [][]float64) [][]float64 {
	n := len(matrix)
	x := make([][]float64, n)
	for i := range x {
		x[i] = make([]float64, n)
	}
	b := make([][]float64, n)
	for i := range x {
		b[i] = make([]float64, n)
	}
	index := make([]int, n)
	for i := 0; i < n; i++ {
		b[i][i] = 1
	}
	elimination(matrix, index)
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			for k := 0; k < n; k++ {
				b[index[j]][k] -= matrix[index[j]][i] * b[index[i]][k]
			}
		}
	}
	for i := 0; i < n; i++ {
		x[n-1][i] = b[index[n-1]][i] / matrix[index[n-1]][n-1]
		for j := n - 2; j >= 0; j-- {
			x[j][i] = b[index[j]][i]
			for k := j + 1; k < n; k++ {
				x[j][i] -= matrix[index[j]][k] * x[k][i]
			}
			x[j][i] /= matrix[index[j]][j]
		}
	}
	return x
}

func Transpose[T any](matrix [][]T) [][]T {
	m := len(matrix)
	n := len(matrix[0])
	transposed := make([][]T, n)
	for i := range transposed {
		transposed[i] = make([]T, m)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			transposed[i][j] = matrix[j][i]
		}
	}
	return transposed
}

func Negate(matrix [][]float64) [][]float64 {
	negated := make([][]float64, len(matrix))
	for i := range negated {
		negated[i] = make([]float64, len(matrix[0]))
	}
	for i := 0; i < len(matrix); i++ {
		for j := 0; j < len(matrix[0]); j++ {
			negated[i][j] = matrix[i][j] * -1
		}
	}
	return negated
}

func Multiply(a [][]float64, b [][]float64) ([][]float64, error) {
	if len(a[0]) != len(b) {
		err := errors.New("las matrices no se pueden multiplicar")
		return nil, err
	}
	result := make([][]float64, len(a))
	for i := range result {
		result[i] = make([]float64, len(b[0]))
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b[0]); j++ {
			result[i][j] = 0
			for k := 0; k < len(a[0]); k++ {
				result[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return result, nil
}

func Add(a [][]float64, b [][]float64) ([][]float64, error) {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		err := errors.New("las matrices no se pueden sumar")
		return nil, err
	}
	result := make([][]float64, len(a))
	for i := range result {
		result[i] = make([]float64, len(a[0]))
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a[0]); j++ {
			result[i][j] = a[i][j] + b[i][j]
		}
	}
	return result, nil
}

func Substract(a [][]float64, b [][]float64) ([][]float64, error) {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		err := errors.New("las matrices no se pueden restar")
		return nil, err
	}
	negB := Negate(b)
	return Add(a, negB)
}

func Determinant(mat [][]float64) (float64, error) {
	if len(mat) != len(mat[0]) {
		return 0.0, errors.New("determinant can only be performed on square matrices")
	}
	if len(mat) == 1 {
		return (mat[0][0]), nil
	}
	if len(mat) == 2 {
		return (mat[0][0] * mat[1][1]) - (mat[0][1] * mat[1][0]), nil
	}

	s := 0.0 // accumulator
	for i := 0; i < len(mat[0]); i++ {
		sm := subMat(mat[1:][:], i)
		z, err := Determinant(sm)

		if err == nil {
			if i%2 != 0 {
				s -= mat[0][i] * z
			} else {
				s += mat[0][i] * z
			}
		}
	}
	return s, nil
}

func subMat(mat [][]float64, p int) [][]float64 {
	stacks := make([]stack, len(mat))
	for n := range mat {
		stacks[n] = stack{}
		for j := range mat[n] {
			if j != p {
				stacks[n].push(mat[n][j])
			}
		}
	}
	out := make([][]float64, len(mat))
	for k := range stacks {
		out[k] = stacks[k].ToSlice()
	}
	return out
}

func Print[T any](matrix [][]T) {
	for _, v := range matrix {
		fmt.Printf("%v\n", v)
	}
}

func Equals[T comparable](a [][]T, b [][]T) bool {
	if len(a) != len(b) || len(a[0]) != len(b[0]) {
		return false
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a[0]); j++ {
			if a[i][j] != b[i][j] {
				return false
			}
		}
	}
	return true
}

func MultiplyByScalar(matrix [][]float64, scalar float64) [][]float64 {
	result := make([][]float64, len(matrix))
	for i := range result {
		result[i] = make([]float64, len(matrix[0]))
	}
	for i := 0; i < len(result); i++ {
		for j := 0; j < len(result[0]); j++ {
			result[i][j] = matrix[i][j] * float64(scalar)
		}
	}
	return result
}

func GaussianElimination(a0 [][]float64, b0 []float64) ([]float64, error) {
	m := len(b0)
	a := make([][]float64, m)
	for i, ai := range a0 {
		row := make([]float64, m+1)
		copy(row, ai)
		row[m] = b0[i]
		a[i] = row
	}
	for k := range a {
		iMax := 0
		max := -1.
		for i := k; i < m; i++ {
			row := a[i]
			// compute scale factor s = max abs in row
			s := -1.
			for j := k; j < m; j++ {
				x := math.Abs(row[j])
				if x > s {
					s = x
				}
			}
			// scale the abs used to pick the pivot.
			if abs := math.Abs(row[k]) / s; abs > max {
				iMax = i
				max = abs
			}
		}
		if a[iMax][k] == 0 {
			return nil, errors.New("singular")
		}
		a[k], a[iMax] = a[iMax], a[k]
		for i := k + 1; i < m; i++ {
			for j := k + 1; j <= m; j++ {
				a[i][j] -= a[k][j] * (a[i][k] / a[k][k])
			}
			a[i][k] = 0
		}
	}
	x := make([]float64, m)
	for i := m - 1; i >= 0; i-- {
		x[i] = a[i][m]
		for j := i + 1; j < m; j++ {
			x[i] -= a[i][j] * x[j]
		}
		x[i] /= a[i][i]
	}
	return x, nil
}

func GaussianJordanElimination(a0 [][]float64, b0 []float64) ([]float64, error) {
	m := len(b0)
	a := make([][]float64, m)
	for i, ai := range a0 {
		row := make([]float64, m+1)
		copy(row, ai)
		row[m] = b0[i]
		a[i] = row
	}
	for k := range a {
		iMax := 0
		max := -1.
		for i := k; i < m; i++ {
			row := a[i]
			// compute scale factor s = max abs in row
			s := -1.
			for j := k; j < m; j++ {
				x := math.Abs(row[j])
				if x > s {
					s = x
				}
			}
			// scale the abs used to pick the pivot.
			if abs := math.Abs(row[k]) / s; abs > max {
				iMax = i
				max = abs
			}
		}
		if a[iMax][k] == 0 {
			return nil, errors.New("singular")
		}
		a[k], a[iMax] = a[iMax], a[k]
		for i := k + 1; i < m; i++ {
			for j := k + 1; j <= m; j++ {
				a[i][j] -= a[k][j] * (a[i][k] / a[k][k])
			}
			a[i][k] = 0
		}
	}
	x := make([]float64, m)
	for i := m - 1; i >= 0; i-- {
		x[i] = a[i][m]
		for j := i + 1; j < m; j++ {
			x[i] -= a[i][j] * x[j]
		}
		x[i] /= a[i][i]
	}
	return x, nil
}

func VectorSum(x []float64, y []float64) []float64 {
	if len(x) != len(y) {
		return nil
	}
	z := make([]float64, len(x))
	for i := range x {
		z[i] = x[i] + y[i]
	}
	return z
}

// VectorSub returns a new vector with each element of x minus the corresponding element of y.
func VectorSub(x []float64, y []float64) []float64 {
	if len(x) != len(y) {
		return nil
	}
	z := make([]float64, len(x))
	for i := range x {
		z[i] = x[i] - y[i]
	}
	return z
}

// VectorScalar returns a new vector with each element multiplied by scalar.
func VectorScalar(x []float64, scalar float64) []float64 {
	z := make([]float64, len(x))
	for i := range x {
		z[i] = x[i] * scalar
	}
	return z
}

// VectorDot returns the dot product of x and y.
func VectorDot(x []float64, y []float64) float64 {
	if len(x) != len(y) {
		return 0
	}
	sum := 0.
	for i := range x {
		sum += x[i] * y[i]
	}
	return sum
}

// VectorNorm returns the norm of x.
func VectorNorm(x []float64) float64 {
	return math.Sqrt(VectorDot(x, x))
}

// VectorNormalize returns a new vector with the same direction as x but with norm 1.
func VectorNormalize(x []float64) []float64 {
	return VectorScalar(x, 1/VectorNorm(x))
}

// VectorCross returns the cross product of x and y.
func VectorCross(x []float64, y []float64) []float64 {
	if len(x) != 3 || len(y) != 3 {
		return nil
	}
	return []float64{
		x[1]*y[2] - x[2]*y[1],
		x[2]*y[0] - x[0]*y[2],
		x[0]*y[1] - x[1]*y[0],
	}
}

// VectorAngleDeg returns the angle between x and y in degrees.
func VectorAngleDeg(x []float64, y []float64) float64 {
	return math.Acos(VectorDot(x, y)/(VectorNorm(x)*VectorNorm(y))) * 180 / math.Pi
}

// VectorAngleRad returns the angle between x and y in radians.
func VectorAngleRad(x []float64, y []float64) float64 {
	return math.Acos(VectorDot(x, y) / (VectorNorm(x) * VectorNorm(y)))
}

// VectorMultiply returns a new vector with each element multiplied by the corresponding element in y.
func VectorMultiply(x []float64, y []float64) []float64 {
	if len(x) != len(y) {
		return nil
	}
	z := make([]float64, len(x))
	for i := range x {
		z[i] = x[i] * y[i]
	}
	return z
}
