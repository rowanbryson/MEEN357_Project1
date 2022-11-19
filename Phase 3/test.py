import timeit

time = timeit.timeit('updated_drag_model.__get_coefficients_for_drag_model()',setup = 'import updated_drag_model ',  number=1)

print(time)