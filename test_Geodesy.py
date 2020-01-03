import pytest
from hypothesis import given, strategies as st
import Geodesy


class TestAngle(object):
    def test_null_instance_has_zero_degrees(self):
        a = Geodesy.Angle()
        assert a.degrees == 0

    @given(decimal_degrees=st.floats(allow_nan=False, allow_infinity=False))
    def test_one_parameter_instance_has_non_zero_degrees(self, decimal_degrees):
        a = Geodesy.Angle(decimal_degrees)
        assert decimal_degrees == a.degrees

    @given(decimal_degrees=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_two_parameter_instance_converts_positive_parameter_from_minutes_to_decimal_degrees(self,
                                                                                                decimal_degrees,
                                                                                                decimal_minutes):
        a = Geodesy.Angle(decimal_degrees, decimal_minutes)
        assert decimal_degrees + decimal_minutes / 60.0 == a.degrees

    @given(decimal_degrees=st.floats(max_value=-0.000000000001, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_two_parameter_instance_converts_negative_parameter_from_minutes_to_decimal_degrees(self,
                                                                                                decimal_degrees,
                                                                                                decimal_minutes):
        a = Geodesy.Angle(decimal_degrees, decimal_minutes)
        assert decimal_degrees - decimal_minutes / 60.0 == a.degrees

    @given(decimal_degrees=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_seconds=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_three_parameter_instance_converts_positive_parameter_from_minutes_to_decimal_degrees(self,
                                                                                                  decimal_degrees,
                                                                                                  decimal_minutes,
                                                                                                  decimal_seconds):
        a = Geodesy.Angle(decimal_degrees, decimal_minutes, decimal_seconds)
        assert decimal_degrees + decimal_minutes / 60.0 + decimal_seconds / 3600.0 == pytest.approx(a.degrees)

    @given(decimal_degrees=st.floats(max_value=-0.000000000001, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_seconds=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_three_parameter_instance_converts_negative_parameter_from_minutes_to_decimal_degrees(self,
                                                                                                  decimal_degrees,
                                                                                                  decimal_minutes,
                                                                                                  decimal_seconds):
        a = Geodesy.Angle(decimal_degrees, decimal_minutes, decimal_seconds)
        assert decimal_degrees - decimal_minutes / 60.0 - decimal_seconds / 3600.0 == pytest.approx(a.degrees)

    @given(decimal_degrees=st.floats(allow_nan=False, allow_infinity=False))
    def test_one_named_parameter_instance_has_non_zero_degrees(self, decimal_degrees):
        a = Geodesy.Angle(degrees=decimal_degrees)
        assert decimal_degrees == a.degrees

    @given(decimal_degrees=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_two_named_parameter_instance_converts_positive_parameter_from_minutes_to_decimal_degrees(
            self, decimal_degrees, decimal_minutes):
        a = Geodesy.Angle(degrees=decimal_degrees, minutes=decimal_minutes)
        assert decimal_degrees + decimal_minutes / 60.0 == a.degrees

    @given(decimal_degrees=st.floats(max_value=-0.000000000001, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_two_named_parameter_instance_converts_negative_parameter_from_minutes_to_decimal_degrees(
            self, decimal_degrees, decimal_minutes):
        a = Geodesy.Angle(degrees=decimal_degrees, minutes=decimal_minutes)
        assert decimal_degrees - decimal_minutes / 60.0 == a.degrees

    @given(decimal_degrees=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_seconds=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_three_named_parameter_instance_converts_positive_parameter_from_minutes_to_decimal_degrees(
            self, decimal_degrees, decimal_minutes, decimal_seconds):
        a = Geodesy.Angle(degrees=decimal_degrees, minutes=decimal_minutes, seconds=decimal_seconds)
        assert decimal_degrees + decimal_minutes / 60.0 + decimal_seconds / 3600.0 == pytest.approx(a.degrees)

    @given(decimal_degrees=st.floats(max_value=-0.000000000001, allow_nan=False, allow_infinity=False),
           decimal_minutes=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False),
           decimal_seconds=st.floats(min_value=0.0, allow_nan=False, allow_infinity=False))
    def test_three_named_parameter_instance_converts_negative_parameter_from_minutes_to_decimal_degrees(
            self, decimal_degrees, decimal_minutes, decimal_seconds):
        a = Geodesy.Angle(degrees=decimal_degrees, minutes=decimal_minutes, seconds=decimal_seconds)
        assert decimal_degrees - decimal_minutes / 60.0 - decimal_seconds / 3600.0 == pytest.approx(a.degrees)
