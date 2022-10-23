from graphing_scripts.analysis_combined_terrain import MARVIN_DICT
from subfunctions import *
import define_rovers
import unittest

# class TestWierd(unittest.TestCase):
#     # this is an example of how to do a unit test
#     def test_wierd(self):

#         function_output = wierd()
#         expected_output = 7

#         assert function_output == expected_output


# class TestGetMass(unittest.TestCase):
#     # this is an example of how to do a unit test
#     def test_get_mass(self):
#         rover = MARVIN_DICT['rover']

#         function_output = get_mass(rover)
#         expected_output = 9

#         assert function_output == expected_output
    
MARVIN_DICT = define_rovers.rover_1_dict

class TestTauDcMotor(unittest.TestCase):
    # add default_motor as a class attribute
    def setUp(self):
        self.default_motor = MARVIN_DICT['rover']['wheel_assembly']['motor']

    def test_from_class_slides(self):
        omega = np.array([0, 0.5, 1, 2, 3, 3.8])
        motor = self.default_motor
        expected = np.array([170, 147.6316, 125.2632, 80.5263, 35.7895, 0])
        actual = tau_dcmotor(omega, motor)
        self.assertTrue(np.allclose(actual, expected, atol=1e-4))

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
                'omega': np.array([2.5, 2]),
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
                self.assertTrue(np.allclose(actual, expected, atol=1e-5))


    # test that the function raises an error as designed when the input is invalid
    def test_tau_dcmotor_input_checking(self):
        # test that the function raises a ValueError when the inputs are the wrong type
        test_cases = [
            {
                'omega': 'string',
                'motor': self.default_motor,
                'expected': Exception
            },
            {
                'omega': np.array([-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]),
                'motor': 'string',
                'expected': Exception
            },
            {
                'omega': np.array([-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]),
                'motor': {'incorrect': 'value'},
                'expected': Exception
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
                'omega': np.array([0.00, 0.50, 1.00, 2.00, 3.00, 3.80]),
                'expected': np.array([10412.50, 9042.4342, 7672.3684, 4932.2368, 2192.1053, 0.00])  # Based on Class Slides
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


class TestGetMass(unittest.TestCase):
    def setUp(self):
        self.default_rover = MARVIN_DICT['rover']

    # test that the function returns the correct value
    def test_get_mass_accuracy(self):
        test_cases = [
            {
                'rover': self.default_rover,
                'expected': 869.0
            },
            # {
            #     'rover': "string",
            #     'expected': Exception
            # }
        ]
        # run subtests for each input
        for test_case in test_cases:
            rover, expected = test_case['rover'], test_case['expected']
            with self.subTest(rover=rover):
                # calculate the actual output
                actual = get_mass(rover)
                # check that the actual output matches the expected output to within 1e-6
                self.assertTrue(np.allclose(actual, expected, atol=1e-6))


class TestFNet(unittest.TestCase):
    def setUp(self):
        self.default_rover = MARVIN_DICT['rover']
        self.default_planet = MARVIN_DICT['planet']

    def test_0_case(self):
        rover = self.default_rover
        planet = self.default_planet
        omega = np.array([0.00])
        terrain_angle = np.array([0.00])
        Crr = 0.01
        expected = np.array([10412.50])

        actual = F_net(omega, terrain_angle, rover, planet, Crr)

        print(f'expected: {list(expected)}')
        print(f'actual: {list(actual)}')

        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    # test that the function returns the correct value
    def test_F_net_accuracy(self):
        rover = self.default_rover
        planet = self.default_planet
        omega = np.array([0.00, 0.50, 1.00, 2.00, 3.00, 3.80])
        terrain_angle = np.array([-5.0, 0, 5.0, 10.0, 20.0, 30.0])
        Crr = 0.1

        expected = np.array([10694.2466, 8720.9744, 7068.5839, 4052.5310, 782.6910, -1896.2983])  # Based on Class Slides
        actual = F_net(omega, terrain_angle, rover, planet, Crr)
        # check that the actual output matches the expected output to within 1e-6
        # print the actual and expected output
        # print(f'testing F_net')
        # print(f'actual: {list(actual)}')
        # print(f'expected: {list(expected)}')
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))


class TestFGravity(unittest.TestCase):
    def setUp(self) -> None:
        self.default_rover = MARVIN_DICT['rover']
        self.default_planet = MARVIN_DICT['planet']

    def test_F_gravity_accuracy_1(self):
        terrain_angle = np.array([-5.0, 0, 5.0, 10.0, 20.0, 30.0])
        rover = self.default_rover
        planet =self.default_planet
        expected = np.array([281.746626465, 0, -281.746626465, -561.34899098, -1105.64167693, -1616.34])

        # calculate the actual output
        actual = F_gravity(terrain_angle, rover, planet)
        # check that the actual output matches the expected output to within 1e-6
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    def test_0_is_0(self):
        terrain_angle = np.array([0.00])
        rover = self.default_rover
        planet = self.default_planet
        expected = np.array([0.00])

        # calculate the actual output
        actual = F_gravity(terrain_angle, rover, planet)
        # check that the actual output matches the expected output to within 1e-6
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))


class TestFRolling(unittest.TestCase):
    def setUp(self) -> None:
        self.default_rover = MARVIN_DICT['rover']
        self.default_planet = MARVIN_DICT['planet']

    def test_0_is_0(self):
        omega = np.array([0.00])
        terrain_angle = np.array([0.00])
        rover = self.default_rover
        planet = self.default_planet
        Crr = 0.01
        expected = np.array([0.00])

        actual = F_rolling(omega, terrain_angle, rover, planet, Crr)
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    def test_F_rolling_accuracy_2(self):
        test_cases = [
            {
                'omega': np.array([0, 0, 1, 1, 3.8, 3.8]),
                'terrain_angle': np.array([0, 30, 0, 30, 0, 30]),
                'rover': self.default_rover,
                'planet': self.default_planet,
                'Crr': 0.1,
                'expected': np.array([0, 0, -323.2679903, -279.95829183, -323.268, -279.958300231])
            },
            {
                'omega': np.array([0, 0, 1, -1, 3.8, -3.8]),
                'terrain_angle': np.array([-10, -20, -30, -40, -50, -60]),
                'rover': self.default_rover,
                'planet': self.default_planet,
                'Crr': 0.1,
                'expected': np.array([0, 0, -279.95829183, 247.637647608, -207.792665008, 161.634])
            }
        ]
        # run subtests for each input
        for test_case in test_cases:
            omega, terrain_angle, rover, planet, Crr = test_case['omega'], test_case['terrain_angle'], test_case['rover'], test_case['planet'], test_case['Crr']
            expected = test_case['expected']
            with self.subTest(omega=omega, terrain_angle=terrain_angle, rover=rover, planet=planet, Crr=Crr, expected=expected):
                # calculate the actual output
                actual = F_rolling(omega, terrain_angle, rover, planet, Crr)

                # print the actual and expected output
                print(f'testing F_rolling')
                print(f'actual: {list(actual)}')
                print(f'expected: {list(expected)}')

                # check that the actual output matches the expected output to within 1e-6
                self.assertTrue(np.allclose(actual, expected, atol=1e-6))


if __name__ == '__main__':
    unittest.main()