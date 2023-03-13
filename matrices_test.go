package matrices

import "testing"

func TestMultiply(t *testing.T) {
	a := [][]float64{
		{1, 2, 3, 4},
		{5, 6, 7, 8},
		{9, 10, 11, 12},
		{13, 14, 15, 16},
	}
	b := [][]float64{
		{16, 15, 14, 13},
		{12, 11, 10, 9},
		{8, 7, 6, 5},
		{4, 3, 2, 1},
	}
	expected := [][]float64{
		{17, 17, 17, 17},
		{17, 17, 17, 17},
		{17, 17, 17, 17},
		{17, 17, 17, 17},
	}
	result, err := Add(a, b)

	if err != nil {
		t.Errorf("%v", err)
	}

	if !Equals(result, expected) {
		t.Errorf("Error adding the matrices.\nExpected:\n%v\nGot:\n%v\n", expected, result)
	}
}
