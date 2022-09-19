from subfunctions import *
import unittest

class TestTauDcMotor(unittest.TestCase):
    # add default_motor as a class attribute
    def setUp(self):
        self.default_motor = MARVIN_DICT['rover']['wheel_assembly']['motor']

    # test that the function returns the correct value
    def test_tau_dcmotor_accuracy(self):
        test_cases = [
            {
                'omega': np.array([-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]),
                'motor': self.default_motor,
                'expected': np.array([170, 170, 170, 147.631578947, 125.263157895, 102.894736842, 80.5263157895, 58.1578947368, 35.7894736842, 13.4210526316, 0, 0, 0, 0, 0])
            },
            {
                'omega': -1,
                'motor': self.default_motor,
                'expected': np.array([170])
            },
            {
                'omega': 0,
                'motor': self.default_motor,
                'expected': np.array([170])
            },
            {
                'omega': 6,
                'motor': self.default_motor,
                'expected': np.array([0])
            },
            {
                'omega': [2.5, 2],
                'motor': self.default_motor,
                'expected': np.array([58.1578947368, 80.5263157895])
            }
        ]
        # run subtests for each input
        for test_case in test_cases:
            omega, motor, expected = test_case['omega'], test_case['motor'], test_case['expected']
            with self.subTest(omega=omega, motor=motor):
                # calculate the actual output
                actual = tau_dcmotor(omega, motor)
                # check that the actual output matches the expected output to within 1e-6
                self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    # test that the function raises an error as designed when the input is invalid
    def test_tau_dcmotor_input_checking(self):
        # test that the function raises a ValueError when the inputs are the wrong type
        test_cases = [
            {
                'omega': 'string',
                'motor': self.default_motor,
                'expected': TypeError
            },
            {
                'omega': np.array([-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]),
                'motor': 'string',
                'expected': TypeError
            },
            {
                'omega': np.array([-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]),
                'motor': {'incorrect': 'value'},
                'expected': ValueError
            }

        ]
        # run subtests for each input
        for test_case in test_cases:
            omega, motor, expected = test_case['omega'], test_case['motor'], test_case['expected']
            with self.subTest(omega=omega, motor=motor, expected=expected):
                # check that the function raises a ValueError
                self.assertRaises(expected, tau_dcmotor, omega, motor)


class TestGetGearRatio(unittest.TestCase):
    # add default_speed_reducer as a class attribute
    def setUp(self):
        self.default_speed_reducer = MARVIN_DICT['rover']['wheel_assembly']['speed_reducer']

    # test that the function returns the correct value
    def test_get_gear_ratio_accuracy(self):
        test_cases = [
            {
                'speed_reducer': self.default_speed_reducer,
                'expected': 3.062500
            },
            {
                'speed_reducer': {
                    'type': 'reverted',
                    'diam_pinion': 0.03,
                    'diam_gear': 0.07,
                    'mass': 1.5,
                },
                'expected': 5.444444444
            }
        ]
        # run subtests for each input
        for test_case in test_cases:
            speed_reducer, expected = test_case['speed_reducer'], test_case['expected']
            with self.subTest(speed_reducer=speed_reducer):
                # calculate the actual output
                actual = get_gear_ratio(speed_reducer)
                # check that the actual output matches the expected output to within 1e-6
                self.assertTrue(np.allclose(actual, expected, atol=1e-6))


class TestFDrive(unittest.TestCase):
    def setUp(self):
        self.default_rover = MARVIN_DICT['rover']

    # test that the function returns the correct value
    def test_F_drive_accuracy(self):
        test_cases = [
            {
                'rover': self.default_rover,
                'omega': np.array([1]),
                'expected': np.array([1278.72807])
            }
        ]
        # run subtests for each input
        for test_case in test_cases:
            rover, omega, expected = test_case['rover'], test_case['omega'], test_case['expected']
            with self.subTest(rover=rover, omega=omega):
                # calculate the actual output
                actual = F_drive(omega, rover)
                # check that the actual output matches the expected output to within 1e-6
                self.assertTrue(np.allclose(actual, expected, atol=1e-6))


if __name__ == '__main__':
    unittest.main()