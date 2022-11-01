from re import L
from subfunctions import *
import define_rovers
import define_experiment
import unittest
import copy

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

class TestTauDcMotor(unittest.TestCase):
    # add default_motor as a class attribute
    def setUp(self):
        self.default_motor = define_rovers.rover1()['wheel_assembly']['motor']

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
        self.default_speed_reducer = define_rovers.rover1()['wheel_assembly']['speed_reducer']

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
        self.default_rover = define_rovers.rover1()

    # test that the function returns the correct value

    def test_F_drive_accuracy(self):
        test_cases = [
            {
                'rover': self.default_rover,
                'omega': np.array([0.00, 0.50, 1.00, 2.00, 3.00, 3.80]),
                # Based on Class Slides
                'expected': np.array([10412.50, 9042.4342, 7672.3684, 4932.2368, 2192.1053, 0.00])
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
        self.default_rover = define_rovers.rover1()

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
        self.default_rover = define_rovers.rover1()
        self.default_planet = define_rovers.planet1()

    def test_0_case(self):
        rover = self.default_rover
        planet = self.default_planet
        omega = np.array([0.00])
        terrain_angle = np.array([0.00])
        Crr = 0.01
        expected = np.array([10412.50])

        actual = F_net(omega, terrain_angle, rover, planet, Crr)

        # print(f'expected: {list(expected)}')
        # print(f'actual: {list(actual)}')

        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    # test that the function returns the correct value
    def test_F_net_accuracy(self):
        rover = self.default_rover
        planet = self.default_planet
        omega = np.array([0.00, 0.50, 1.00, 2.00, 3.00, 3.80])
        terrain_angle = np.array([-5.0, 0, 5.0, 10.0, 20.0, 30.0])
        Crr = 0.1

        expected = np.array([10694.2466, 8720.9744, 7068.5839,
                            4052.5310, 782.6910, -1896.2983])  # Based on Class Slides
        actual = F_net(omega, terrain_angle, rover, planet, Crr)
        # check that the actual output matches the expected output to within 1e-6
        # print the actual and expected output
        # print(f'testing F_net')
        # print(f'actual: {list(actual)}')
        # print(f'expected: {list(expected)}')
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))


class TestFGravity(unittest.TestCase):
    def setUp(self) -> None:
        self.default_rover = define_rovers.rover1()
        self.default_planet = define_rovers.planet1()

    def test_F_gravity_accuracy_1(self):
        terrain_angle = np.array([-5.0, 0, 5.0, 10.0, 20.0, 30.0])
        rover = self.default_rover
        planet = self.default_planet
        expected = np.array(
            [281.746626465, 0, -281.746626465, -561.34899098, -1105.64167693, -1616.34])

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
        self.default_rover = define_rovers.rover1()
        self.default_planet = define_rovers.planet1()

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
            omega, terrain_angle, rover, planet, Crr = test_case['omega'], test_case[
                'terrain_angle'], test_case['rover'], test_case['planet'], test_case['Crr']
            expected = test_case['expected']
            with self.subTest(omega=omega, terrain_angle=terrain_angle, rover=rover, planet=planet, Crr=Crr, expected=expected):
                # calculate the actual output
                actual = F_rolling(omega, terrain_angle, rover, planet, Crr)

                # print the actual and expected output
                # print(f'testing F_rolling')
                # print(f'actual: {list(actual)}')
                # print(f'expected: {list(expected)}')

                # check that the actual output matches the expected output to within 1e-6
                self.assertTrue(np.allclose(actual, expected, atol=1e-6))


class TestMotorW(unittest.TestCase):
    def setUp(self) -> None:
        self.default_rover = define_rovers.rover1()

    def test_from_slide_1(self):
        v = np.array([0.1])
        expected = np.array([1.02])
        actual = motorW(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-2))

    def test_from_slide_2(self):
        v = np.array([0.3])
        expected = np.array([3.06])
        actual = motorW(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-2))

    def test_nparray(self):
        v = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        expected = np.array([0, 10.20833333, 20.41666667, 30.625, 40.83333333,
                            51.04166667, 61.25, 71.45833333, 81.66666667, 91.875, 102.08333333])
        actual = motorW(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    def test_int(self):
        v = 1
        expected = 10.20833333
        actual = motorW(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    def test_float(self):
        v = 1.0
        expected = 10.20833333
        actual = motorW(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-6))

    def test_string_exception(self):
        v = '1'
        with self.assertRaises(TypeError):
            motorW(v, self.default_rover)

    def test_list_exception(self):
        v = [1]
        with self.assertRaises(TypeError):
            motorW(v, self.default_rover)

    def test_invalid_rover_exception(self):
        v = 1
        rover = copy.deepcopy(self.default_rover)
        del rover['wheel_assembly']['speed_reducer']
        with self.assertRaises(KeyError):
            motorW(v, rover)

class TestRoverDynamics(unittest.TestCase):
    def setUp(self) -> None :
        self.default_rover = define_rovers.rover1()
        self.default_planet = define_rovers.planet1()
        self.default_experiment, self.default_endevent = define_experiment.experiment1()
        self.default_experiment2, self.default_endevent2 = define_experiment.experiment2()

    def test_from_slide1(self):
        t = 0
        y = np.array([.33, 0])
        expected = np.array([.253, .33])
        actual = rover_dynamics(t, y, self.default_rover, self.default_planet, self.default_experiment)
        self.assertTrue(np.allclose(actual, expected, atol=1e-2))

    def test_from_slide2(self):
        t = 20
        y = np.array([.25, 500])
        expected = np.array([2.86, .25])
        actual = rover_dynamics(t, y, self.default_rover, self.default_planet, self.default_experiment)
        self.assertTrue(np.allclose(actual, expected, atol= 1e-2))

    def test_float(self):
        t = 20.0
        y = np.array([.25, 500.0])
        expected = np.array([2.86, .25])
        actual = rover_dynamics(t, y, self.default_rover, self.default_planet, self.default_experiment)
        self.assertTrue(np.allclose(actual, expected, atol = 1e-2))

    def test_string_exception(self):
        t = 'string'
        y = np.array([0, 0])
        with self.assertRaises(Exception):
            val = rover_dynamics(t, y, self.default_rover, self.default_planet, self.default_experiment)
            print(val)

    def test_dict_exception(self):
        t = 10
        y = np.array([10, 11])
        rover_test = 10
        with self.assertRaises(Exception):
            rover_dynamics(t, y, rover_test, self.default_planet, self.default_experiment)

    def test_array_exception(self):
        t = 10
        y = 10
        with self.assertRaises(Exception):
            rover_dynamics(t, y, self.default_rover, self.default_planet, self.default_experiment)

    def test_y_size(self):
        t = 10
        y = np.array([1, 2, 3])
        with self.assertRaises(Exception):
            rover_dynamics(t, y, self.default_rover, self.default_planet, self.default_experiment)


class TestMechPower(unittest.TestCase):
    def setUp(self) -> None:
        self.default_rover = define_rovers.rover1()

    def test_from_slide_1(self):
        v = np.array([0.05])
        expected = np.array([75.1])
        actual = mechpower(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-1))

    def test_from_slide_2(self):
        v = np.array([0.25])
        expected = np.array([142])
        actual = mechpower(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e0))

    def test_nparray(self):
        v = np.array([0.05, 0.25])
        expected = np.array([75.1, 142])
        actual = mechpower(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1))

    def test_float(self):
        v = 0.05
        expected = 75.1
        actual = mechpower(v, self.default_rover)
        self.assertTrue(np.allclose(actual, expected, atol=1e-1))

    def test_string_exception(self):
        v = '1'
        with self.assertRaises(TypeError):
            mechpower(v, self.default_rover)

    def test_list_exception(self):
        v = [1]
        with self.assertRaises(TypeError):
            mechpower(v, self.default_rover)

    def test_invalid_rover_exception(self):
        v = np.array([1.0])
        rover = copy.deepcopy(self.default_rover)
        del rover['wheel_assembly']['speed_reducer']
        with self.assertRaises(KeyError):
            mechpower(v, rover)
    

class Testbattenergy(unittest.TestCase):
    def setUp(self) -> None:
        self.default_rover = define_rovers.rover1()

    def test_case_from_slide_1(self):
        t = np.array([0, 1, 2, 3, 4, 5, 6])
        v = np.array([0.33, 0.32, 0.33, 0.2, 0.2, 0.25, 0.28])
        expected = 6.8e3
        actual = battenergy(t, v, self.default_rover)
        # allowed error is big here because I think the test case given to us is not very accurate
        self.assertAlmostEqual(actual, expected, delta=0.5e3)
    
    def test_string_exception(self):
        v = '1'
        t = '1'
        with self.assertRaises(TypeError):
            battenergy(t, v, self.default_rover)
    
    def test_list_exception(self):
        v = [1]
        t = [1]
        with self.assertRaises(TypeError):
            battenergy(t, v, self.default_rover)
    
    def test_invalid_rover_exception(self):
        v = np.array([1])
        t = np.array([1])
        rover = copy.deepcopy(self.default_rover)
        del rover['wheel_assembly']['motor']
        with self.assertRaises(KeyError):
            battenergy(t, v, rover)


if __name__ == '__main__':
    unittest.main()
