from timeit import timeit

print("Testing NumPy speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_numpy',
             stmt='compute_numpy.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 0.015417880000313744

print("Testing Python speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_py',
             stmt='compute_py.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 25.86546850400191

print("Testing Cython speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_cy',
             stmt='compute_cy.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 20.847004840004956

print("Testing Cython (typed) speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_typed',
             stmt='compute_typed.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 8.942294986001798

print("Testing Cython (typed+memview) speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_memview',
             stmt='compute_memview.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 0.015932380003505386

print("Testing Cython (typed+memview+unsafe) speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_unsafe',
             stmt='compute_unsafe.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 0.008792597996944096

print("Testing Cython (typed+memview+unsafe+contiguous) speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_contiguous',
             stmt='compute_contiguous.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 0.007633947003341746

print("Testing Cython (typed+memview+unsafe+contiguous+infer) speed:", end=" ")
print(timeit(setup='import input_data as d; import compute_infer_types',
             stmt='compute_infer_types.compute(d.array_1, d.array_2, d.a, d.b, d.c)',
             number=1))
# Seconds: 0.008812448002572637
